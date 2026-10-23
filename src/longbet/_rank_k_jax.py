# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Cross-sectional JAX BCF + conjugate-block sampler for Rank-K Exposure-Basis BCF.

Implements the per-mode cross-sectional generative model (Section 3.2 of the
companion paper):

    dy_{ik} = mu_k(x_i, pi_hat_i, df_ik^(m)) + theta_k df_ik^(m)
              + lambda_i' f_k + zeta_{ik} + z_i tau_k(x_i) + eps_{ik},

where:
  * ``mu_k`` and ``tau_k`` are ``bartz`` BART forests (with ``change_step``,
    ``change_internal_step``, and ``regrow_step`` matching ``longbet._step``),
  * ``theta_k`` is the mode-k CUPED forecast slope with truncated-Normal prior
    on ``[theta_min, theta_max]``,
  * ``f_k ~ N(0, s_f^2 I_q)`` are post-anchor factor shocks multiplying the
    pre-anchor residual loadings ``lambda_i`` (from pre-adoption history only),
  * ``zeta_k ~ N(0, s_zeta^2 C(ell_k))`` is a Matérn-3/2 spatial random field
    over great-circle distances (or cluster random effects when only
    ``cluster_ids`` are supplied),
  * ``eps_{ik} ~ N(0, sigma_k^2 / w_{ik})`` is heteroskedastic noise with known
    exposure weights ``w_{ik}``.

Under ``baseline_cells="untreated"`` (the default), all baseline components
(``mu_k, theta_k, f_k, zeta_{k,I_0}, ell_k, s_zeta^2, sigma_k^2``) are updated
strictly from never-treated units ``I_0``, ``zeta_{k,I_1}`` is kriged from
``zeta_{k,I_0}`` immediately before the treatment forest step, and ``tau_k`` is
updated on ``I_1``.

With the BART forests (``use_forests=True``, the only setting used by
:func:`longbet.fit_rank_k_longbet`) a sweep is a partially collapsed,
approximate Gibbs step, so the sampler is not an exact posterior sampler: the
``mu_k`` and ``tau_k`` steps use diagonal working likelihoods with noise
variances ``sigma_k^2 + s_zeta^2 / 2`` and ``sigma_k^2 + s_zeta^2`` in place of
the spatially correlated likelihood, and the treated residual holds
``zeta_{k,I_1}`` at its kriging mean ``A_g E[zeta_{k,I_0} | control residuals]``
(the deviation of the kriging draw from that mean enters only the draw of the
scalar treatment intercept). ``use_forests=False`` replaces both forests by
constants (the linear-Gaussian verification limit, used only in the tests,
which compare it with :func:`closed_form_linear_gaussian_posterior`).
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from pathlib import Path
import tempfile
from typing import Any

_DEFAULT_JAX_CACHE = str(Path(tempfile.gettempdir()) / "longbet_jax_xla_cache")
os.environ.setdefault("JAX_COMPILATION_CACHE_DIR", _DEFAULT_JAX_CACHE)
os.environ.setdefault("JAX_PERSISTENT_CACHE_MIN_COMPILE_TIME_SECS", "0")

import equinox as eqx
import jax
import jax.numpy as jnp
from jax import lax, random
import jax.scipy.linalg
import jax.scipy.special
import numpy as np

try:
    Path(_DEFAULT_JAX_CACHE).mkdir(parents=True, exist_ok=True)
    jax.config.update("jax_compilation_cache_dir", _DEFAULT_JAX_CACHE)
    jax.config.update("jax_persistent_cache_min_compile_time_secs", 0)
except Exception:
    pass

from bartz._jaxext import field
from bartz.mcmcstep import Forest, State, Wishart, init, make_p_nonterminal
from bartz.mcmcstep import step as bartz_step
from bartz.mcmcstep._axes import CHAIN_AXIS
from bartz.mcmcstep._step import apply_moves_to_leaf_indices
from bartz.prepcovars import UniqueQuantileBinner

from longbet._change_internal_move import change_internal_step
from longbet._change_move import change_step
from longbet._forest_cache import current_forest_fit, refresh_prec_tree
from longbet._regrow_move import regrow_step
from longbet._ridge import compute_active_leaf_stats
from longbet._scales import compute_treatment_scale_attrs
from longbet._x64 import enable_x64


EARTH_RADIUS_KM = 6371.0088
DEFAULT_LENGTHSCALE_GRID_KM = (100.0, 200.0, 400.0, 800.0, 1600.0)
# Variance-ratio grid r = s_zeta^2 / sigma_k^2 used for collapsed marginal updates
DEFAULT_RATIO_GRID = (
    0.0, 0.02, 0.05, 0.10, 0.20, 0.35, 0.55, 0.85, 1.30, 2.00, 3.20, 5.00, 8.00, 13.00, 22.00, 36.00
)


def great_circle_distance_km(coords_a: np.ndarray, coords_b: np.ndarray) -> np.ndarray:
    """Pairwise great-circle (haversine) distance in km between (lat, lon) degrees."""
    a = np.deg2rad(np.asarray(coords_a, dtype=np.float64))
    b = np.deg2rad(np.asarray(coords_b, dtype=np.float64))
    lat1 = a[:, 0:1]
    lon1 = a[:, 1:2]
    lat2 = b[:, 0:1].T
    lon2 = b[:, 1:2].T
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    h = np.sin(0.5 * dlat) ** 2 + np.cos(lat1) * np.cos(lat2) * (np.sin(0.5 * dlon) ** 2)
    h = np.clip(h, 0.0, 1.0)
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(h))


def matern32_corr(dist_km: np.ndarray, ell_km: float) -> np.ndarray:
    """Matérn-3/2 correlation matrix at range ``ell_km``."""
    scaled = math.sqrt(3.0) * np.asarray(dist_km, dtype=np.float64) / max(float(ell_km), 1e-6)
    return (1.0 + scaled) * np.exp(-scaled)


class SpatialKernelCache(eqx.Module):
    """Precomputed spectral and kriging operators across the lengthscale grid."""

    enabled: bool = eqx.field(static=True)
    n_ell: int = eqx.field(static=True)
    n_ratio: int = eqx.field(static=True)
    ell_grid: jax.Array
    ratio_grid: jax.Array
    ratio_log_prior: jax.Array
    eig_vecs_00: jax.Array
    eig_vals_00: jax.Array
    krig_mat_10: jax.Array
    schur_chol_11: jax.Array
    sqrt_w0: jax.Array
    inv_sqrt_w0: jax.Array


@enable_x64()
def _kriging_operator(C00: np.ndarray, C10: np.ndarray) -> np.ndarray:
    """Return the kriging operator ``A = C10 C00^{-1}``, robust to a non-PD ``C00``.

    The Matérn-3/2 correlation is positive definite in the plane but not
    guaranteed to be on the sphere with great-circle distances, so for long
    ranges and spread-out (e.g. global) unit sets the jittered control matrix
    ``C00`` can be indefinite. The Cholesky solve is used whenever it succeeds,
    so those results are unchanged; otherwise ``C00`` is inverted through its
    symmetric eigen-decomposition with the eigenvalues floored at ``1e-6`` times
    the largest one.

    Args:
      C00: ``(n0, n0)`` correlation matrix of the control units, jittered by
        ``1e-6 I``.
      C10: ``(n1, n0)`` correlation between treated and control units.

    Returns:
      The ``(n1, n0)`` kriging matrix.
    """
    try:
        L00 = np.linalg.cholesky(C00)
    except np.linalg.LinAlgError:
        vals, vecs = np.linalg.eigh(0.5 * (C00 + C00.T))
        vals = np.maximum(vals, 1e-6 * float(np.max(vals)))
        return ((C10 @ vecs) / vals[None, :]) @ vecs.T
    # Solve C00 X = C01 -> X = C00^{-1} C10^T
    invC00_C01 = jax.scipy.linalg.cho_solve((L00, True), C10.T)
    return np.asarray(invC00_C01, dtype=np.float64).T


@enable_x64()
def build_spatial_kernel_cache(
    *,
    ctrl_idx: np.ndarray,
    trt_idx: np.ndarray,
    w_unit: np.ndarray,
    coords: np.ndarray | None = None,
    cluster_ids: np.ndarray | None = None,
    ell_grid_km: tuple[float, ...] = DEFAULT_LENGTHSCALE_GRID_KM,
    ratio_grid: tuple[float, ...] = DEFAULT_RATIO_GRID,
    use_spatial_gp: bool = True,
) -> SpatialKernelCache:
    """Precompute eigendecompositions and kriging matrices on ``I_0`` and ``I_1``."""
    n0 = int(ctrl_idx.size)
    n1 = int(trt_idx.size)
    w0 = np.maximum(np.asarray(w_unit[ctrl_idx], dtype=np.float64), 1e-8)
    sqrt_w0 = np.sqrt(w0)
    inv_sqrt_w0 = 1.0 / sqrt_w0

    if not use_spatial_gp or (coords is None and cluster_ids is None):
        return SpatialKernelCache(
            enabled=False,
            n_ell=1,
            n_ratio=1,
            ell_grid=jnp.asarray([1.0], dtype=jnp.float64),
            ratio_grid=jnp.asarray([0.0], dtype=jnp.float64),
            ratio_log_prior=jnp.asarray([0.0], dtype=jnp.float64),
            eig_vecs_00=jnp.eye(n0, dtype=jnp.float64)[None, :, :],
            eig_vals_00=jnp.zeros((1, n0), dtype=jnp.float64),
            krig_mat_10=jnp.zeros((1, n1, n0), dtype=jnp.float64),
            schur_chol_11=jnp.zeros((1, n1, n1), dtype=jnp.float64),
            sqrt_w0=jnp.asarray(sqrt_w0, dtype=jnp.float64),
            inv_sqrt_w0=jnp.asarray(inv_sqrt_w0, dtype=jnp.float64),
        )

    corr_list: list[np.ndarray] = []
    ell_vals: list[float] = []
    if coords is not None:
        c_arr = np.asarray(coords, dtype=np.float64)
        dist_all = great_circle_distance_km(c_arr, c_arr)
        for ell in ell_grid_km:
            corr_list.append(matern32_corr(dist_all, float(ell)))
            ell_vals.append(float(ell))
    else:
        cid = np.asarray(cluster_ids)
        corr_mat = (cid[:, None] == cid[None, :]).astype(np.float64)
        corr_list.append(corr_mat)
        ell_vals.append(1.0)

    eig_vecs_list: list[np.ndarray] = []
    eig_vals_list: list[np.ndarray] = []
    krig_list: list[np.ndarray] = []
    schur_list: list[np.ndarray] = []

    for C in corr_list:
        C00 = C[np.ix_(ctrl_idx, ctrl_idx)] + 1e-6 * np.eye(n0, dtype=np.float64)
        C10 = C[np.ix_(trt_idx, ctrl_idx)]
        C11 = C[np.ix_(trt_idx, trt_idx)] + 1e-6 * np.eye(n1, dtype=np.float64)

        # Symmetrized weighted covariance: W0^{1/2} C00 W0^{1/2}
        C00_w = (sqrt_w0[:, None] * C00) * sqrt_w0[None, :]
        C00_w = 0.5 * (C00_w + C00_w.T)
        d_g, U_g = np.linalg.eigh(C00_w)
        d_g = np.maximum(d_g, 0.0)

        # Kriging operator A_g = C10 C00^{-1} and conditional Schur complement
        A_g = _kriging_operator(C00, C10)
        S_g = C11 - A_g @ C10.T
        S_g = 0.5 * (S_g + S_g.T)
        vals_s, vecs_s = np.linalg.eigh(S_g)
        vals_s = np.maximum(vals_s, 1e-8)
        S_g_pd = (vecs_s * vals_s[None, :]) @ vecs_s.T
        S_g_pd = 0.5 * (S_g_pd + S_g_pd.T)
        L_schur = np.linalg.cholesky(S_g_pd)

        eig_vecs_list.append(U_g)
        eig_vals_list.append(d_g)
        krig_list.append(A_g)
        schur_list.append(L_schur)

    r_arr = np.asarray(ratio_grid, dtype=np.float64)
    # Discrete prior on the variance-ratio grid r = s_zeta^2 / sigma_k^2 (not a
    # penalized-complexity prior): grid point r gets weight proportional to
    # (1 - rho)^2 = (1 + r)^(-2), with rho = r / (1 + r) the spatial variance share,
    # i.e. the Beta(1, 3) density of rho evaluated at the grid points and renormalized.
    # The weight is largest at rho = 0 (no spatial effect). The lengthscale grid has a
    # discrete uniform prior.
    rho_arr = r_arr / (1.0 + r_arr)
    ratio_log_prior = 2.0 * np.log(np.maximum(1.0 - rho_arr, 1e-12))
    ratio_log_prior = ratio_log_prior - scipy_logsumexp_np(ratio_log_prior)

    return SpatialKernelCache(
        enabled=True,
        n_ell=len(ell_vals),
        n_ratio=int(r_arr.size),
        ell_grid=jnp.asarray(ell_vals, dtype=jnp.float64),
        ratio_grid=jnp.asarray(r_arr, dtype=jnp.float64),
        ratio_log_prior=jnp.asarray(ratio_log_prior, dtype=jnp.float64),
        eig_vecs_00=jnp.asarray(np.stack(eig_vecs_list, axis=0), dtype=jnp.float64),
        eig_vals_00=jnp.asarray(np.stack(eig_vals_list, axis=0), dtype=jnp.float64),
        krig_mat_10=jnp.asarray(np.stack(krig_list, axis=0), dtype=jnp.float64),
        schur_chol_11=jnp.asarray(np.stack(schur_list, axis=0), dtype=jnp.float64),
        sqrt_w0=jnp.asarray(sqrt_w0, dtype=jnp.float64),
        inv_sqrt_w0=jnp.asarray(inv_sqrt_w0, dtype=jnp.float64),
    )


def scipy_logsumexp_np(x: np.ndarray) -> float:
    m = float(np.max(x))
    return m + float(np.log(np.sum(np.exp(x - m))))


@enable_x64()
def _sample_truncated_normal_64(
    key: jax.Array,
    mean: jax.Array,
    sd: jax.Array,
    low: float,
    high: float,
) -> jax.Array:
    """Exact inverse-CDF draw from N(mean, sd^2) truncated to [low, high] in float64."""
    f64 = jnp.float64
    mean64 = jnp.asarray(mean, dtype=f64)
    sd64 = jnp.maximum(jnp.asarray(sd, dtype=f64), 1e-12)
    a = (jnp.asarray(low, dtype=f64) - mean64) / sd64
    b = (jnp.asarray(high, dtype=f64) - mean64) / sd64
    u = random.uniform(key, (), dtype=jnp.float32, minval=1e-6, maxval=1.0 - 1e-6).astype(f64)

    # Normal CDF branch when [a, b] has appreciable mass; fallback to clipped mode in extreme tails
    cdf_a = jax.scipy.special.ndtr(a)
    cdf_b = jax.scipy.special.ndtr(b)
    mass = cdf_b - cdf_a
    p = jnp.clip(cdf_a + u * jnp.maximum(mass, 1e-15), 1e-15, 1.0 - 1e-15)
    z_cdf = jax.scipy.special.ndtri(p)
    z_clip = jnp.clip(z_cdf, a, b)
    # If both a, b are far in the same tail and mass underflows, fall back to boundary closest to mean
    z_safe = jnp.where(mass > 1e-12, z_clip, jnp.clip(0.0, a, b))
    return mean64 + sd64 * z_safe


@enable_x64()
def _sample_inv_gamma_64(key: jax.Array, shape_a: jax.Array, rate_b: jax.Array) -> jax.Array:
    """Exact Inverse-Gamma(shape_a, rate_b) draw with division in float64."""
    rate = jnp.asarray(rate_b, dtype=jnp.float64)
    g = random.gamma(key, jnp.asarray(shape_a, dtype=jnp.float32), dtype=jnp.float32)
    return rate / jnp.maximum(jnp.asarray(g, dtype=jnp.float64), 1e-30)


@enable_x64()
def _slice_sample_scalar_64(
    key: jax.Array,
    x0: jax.Array,
    log_density: Any,
    lower: jax.Array,
    upper: jax.Array,
    width: jax.Array,
    max_steps: int = 32,
) -> jax.Array:
    """One stepping-out / shrinkage slice-sampling update (Neal 2003) of a scalar in float64.

    ``log_density`` must return ``-inf`` outside the support, which is further
    restricted to ``[lower, upper]`` (either may be infinite). The stepping-out
    loops are bounded by ``max_steps`` and the shrinkage loop by ``4 * max_steps``,
    so the update is jittable; if shrinkage exhausts its budget the current value
    is kept (a valid, if lazy, Markov transition).
    """
    f64 = jnp.float64
    k_u, k_pos, k_shrink = random.split(key, 3)
    x0 = jnp.asarray(x0, dtype=f64)
    lower = jnp.asarray(lower, dtype=f64)
    upper = jnp.asarray(upper, dtype=f64)
    width = jnp.maximum(jnp.asarray(width, dtype=f64), 1e-12)
    log_u = jnp.log(random.uniform(k_u, (), dtype=jnp.float32, minval=1e-12, maxval=1.0).astype(f64))
    log_y = log_density(x0) + log_u
    u_pos = random.uniform(k_pos, (), dtype=jnp.float32).astype(f64)
    left = jnp.maximum(x0 - u_pos * width, lower)
    right = jnp.minimum(left + width, upper)

    def _cond_left(carry):
        pos, i = carry
        return (i < max_steps) & (pos > lower) & (log_density(pos) > log_y)

    def _body_left(carry):
        pos, i = carry
        return jnp.maximum(pos - width, lower), i + 1

    left, _ = lax.while_loop(_cond_left, _body_left, (left, jnp.int32(0)))

    def _cond_right(carry):
        pos, i = carry
        return (i < max_steps) & (pos < upper) & (log_density(pos) > log_y)

    def _body_right(carry):
        pos, i = carry
        return jnp.minimum(pos + width, upper), i + 1

    right, _ = lax.while_loop(_cond_right, _body_right, (right, jnp.int32(0)))

    def _cond_shrink(carry):
        _, _, _, i, done, _ = carry
        return (~done) & (i < 4 * max_steps)

    def _body_shrink(carry):
        x, lo, hi, i, done, k = carry
        k, k1 = random.split(k)
        x1 = lo + random.uniform(k1, (), dtype=jnp.float32).astype(f64) * (hi - lo)
        accept = log_density(x1) >= log_y
        lo_new = jnp.where(accept | (x1 >= x0), lo, x1)
        hi_new = jnp.where(accept | (x1 < x0), hi, x1)
        return jnp.where(accept, x1, x), lo_new, hi_new, i + 1, done | accept, k

    x_new, _, _, _, _, _ = lax.while_loop(
        _cond_shrink, _body_shrink, (x0, left, right, jnp.int32(0), jnp.bool_(False), k_shrink)
    )
    return x_new


def _log_prior_tau_64(
    x: jax.Array,
    family: int,
    prior_mean: jax.Array,
    prior_sd: jax.Array,
    p1: jax.Array,
    p2: jax.Array,
) -> jax.Array:
    """Log prior density (up to a constant) of the ATT / mode-1 intercept under ``tau_prior_family``.

    * ``0`` Normal(prior_mean, prior_sd^2);
    * ``1`` half-Student-t: ``x >= 0``, ``df = p1``, scale ``p2 * prior_sd``
      (the YTLR tier: ``T^+_{2.5}(0, (0.5 sigma_0)^2)``);
    * ``2`` Gamma spike: ``x > 0``, shape ``p1``, scale ``p2 * prior_sd``
      (``Gamma(0.40, 1.60 sigma_0)``, unbounded density at 0).
    """
    f64 = jnp.float64
    x = jnp.asarray(x, dtype=f64)
    sd = jnp.maximum(jnp.asarray(prior_sd, dtype=f64), 1e-12)
    if family == 1:
        df = jnp.asarray(p1, dtype=f64)
        scale = jnp.maximum(jnp.asarray(p2, dtype=f64) * sd, 1e-12)
        val = -0.5 * (df + 1.0) * jnp.log1p(jnp.square(x / scale) / df)
        return jnp.where(x >= 0.0, val, -jnp.inf)
    if family == 2:
        shape = jnp.asarray(p1, dtype=f64)
        scale = jnp.maximum(jnp.asarray(p2, dtype=f64) * sd, 1e-12)
        val = (shape - 1.0) * jnp.log(jnp.maximum(x, 1e-300)) - x / scale
        return jnp.where(x > 0.0, val, -jnp.inf)
    return -0.5 * jnp.square((x - jnp.asarray(prior_mean, dtype=f64)) / sd)


class RankKModeState(State):
    """Cross-sectional MCMC state for one mode ``k`` of Rank-K BCF."""

    forest_nu: Forest = field()
    resid_nu: jax.Array = field(chains=CHAIN_AXIS, data=-1)
    prec_scale_nu: jax.Array | None = field(chains=CHAIN_AXIS, data=-1)
    inv_sdev_scale_nu: jax.Array | None = field(chains=CHAIN_AXIS, data=-1)
    inv_sdev_unit_nu: jax.Array = field(chains=CHAIN_AXIS)
    resid_unit_nu: jax.Array = field()
    leaf_prior_cov_inv_nu: jax.Array = field()
    resid_eff_scale_nu: Any = field(chains=CHAIN_AXIS)
    resid_inexact_integral_nu: Any = field(chains=CHAIN_AXIS)

    # Cached forest fits and parametric blocks (per chain)
    mu_fit: jax.Array = field(chains=CHAIN_AXIS, data=-1)
    nu_fit: jax.Array = field(chains=CHAIN_AXIS, data=-1)
    intercept: jax.Array = field(chains=CHAIN_AXIS)
    tau_intercept: jax.Array = field(chains=CHAIN_AXIS)
    theta_f: jax.Array = field(chains=CHAIN_AXIS)
    f_loadings: jax.Array = field(chains=CHAIN_AXIS)
    s_f2: jax.Array = field(chains=CHAIN_AXIS)
    zeta: jax.Array = field(chains=CHAIN_AXIS, data=-1)
    s_zeta2: jax.Array = field(chains=CHAIN_AXIS)
    ell_idx: jax.Array = field(chains=CHAIN_AXIS)
    sigma2: jax.Array = field(chains=CHAIN_AXIS)

    # Shared data arrays across chains for a fixed forecast draw m
    z_vec: jax.Array = field(data=-1)
    w_unit: jax.Array = field(data=-1)
    df_vec: jax.Array = field(data=-1)
    df_point: jax.Array = field(data=-1)
    lambda_mat: jax.Array = field(data=-1)
    ctrl_mask: jax.Array = field(data=-1)
    trt_mask: jax.Array = field(data=-1)
    ctrl_idx: jax.Array = field()
    trt_idx: jax.Array = field()

    # Scalar prior hyperparameters (traced JAX arrays so varying theta_f_clip/tau_clip never recompiles XLA)
    theta_min: jax.Array = field()
    theta_max: jax.Array = field()
    theta_prior_mean: jax.Array = field()
    theta_prior_sd: jax.Array = field()
    intercept_prior_sd: jax.Array = field()
    sigma_prior_a: jax.Array = field()
    sigma_prior_b: jax.Array = field()
    sf_prior_a: jax.Array = field()
    sf_prior_b: jax.Array = field()
    tau_prior_mean: jax.Array = field()
    tau_prior_sd: jax.Array = field()
    tau_min: jax.Array = field()
    tau_max: jax.Array = field()
    # ATT weights a_i = L_{G_i} / |O_trt| on treated units (0 elsewhere). When
    # ``att_centering`` is set, the mode-1 treatment forest is re-centered with
    # these weights instead of the likelihood weights, so that the scalar
    # intercept equals the ATT and any prior / truncation on it acts on the ATT.
    w_att: jax.Array = field(data=-1)
    # Shape parameters of the non-Gaussian ATT prior families
    # (``tau_prior_family``: 1 = half-Student-t (df, scale multiplier),
    # 2 = Gamma spike (shape, scale multiplier)); unused for the Normal family.
    tau_prior_p1: jax.Array = field()
    tau_prior_p2: jax.Array = field()
    # Weight of the spatial variance in the working noise variance of the tau_k step.
    nu_working_zeta_weight: jax.Array = field()
    # Half-Cauchy scale (squared, data units, per leaf) of the learned leaf prior
    # variance of the treatment forest (``learn_sigma_nu``).
    sigma_nu_prior_scale2: jax.Array = field()
    # Treated-arm excess residual variance s_eta^2 (``treated_extra_variance``; 0 otherwise).
    s_eta2: jax.Array = field(chains=CHAIN_AXIS)
    # Parameter-expansion auxiliary of the half-Cauchy leaf-scale prior (``learn_sigma_nu``).
    xi_nu: jax.Array = field(chains=CHAIN_AXIS)

    # Static configuration flags (structural / control-flow only)
    N_units: int = field(static=True, default=0)
    N_ctrl: int = field(static=True, default=0)
    N_trt: int = field(static=True, default=0)
    # Units with a positive precision weight (all of them on complete panels). A unit
    # with ``w_unit == 0`` -- a masked panel's unit whose window does not support the
    # mode -- adds nothing to the weighted sums of squares, so it must not count in
    # the shape parameter of the inverse-gamma variance updates either.
    N_ctrl_dof: int = field(static=True, default=0)
    N_trt_dof: int = field(static=True, default=0)
    N_all_dof: int = field(static=True, default=0)
    q_factors: int = field(static=True, default=0)
    has_forecast: bool = field(static=True, default=True)
    use_forests: bool = field(static=True, default=True)
    use_factor_loadings: bool = field(static=True, default=True)
    use_spatial_gp: bool = field(static=True, default=True)
    baseline_untreated_only: bool = field(static=True, default=True)
    # Static gates of the extensions; all False / 0 reproduces the paper's sampler
    # bit-for-bit because the gated code is not traced at all.
    att_centering: bool = field(static=True, default=False)
    tau_prior_family: int = field(static=True, default=0)
    treated_extra_variance: bool = field(static=True, default=False)
    learn_sigma_nu: bool = field(static=True, default=False)


def _load(r: jax.Array, resid_unit: jax.Array) -> jax.Array:
    return (r / resid_unit).astype(jnp.float32)


@enable_x64(False)
def init_rank_k_mode(
    *,
    X_binned: jax.Array,
    max_split_mu: jax.Array,
    max_split_nu: jax.Array,
    dy_k: jax.Array,
    df_k: jax.Array,
    z_vec: jax.Array,
    w_unit: jax.Array,
    lambda_mat: jax.Array,
    ctrl_idx: np.ndarray,
    trt_idx: np.ndarray,
    num_trees_mu: int = 75,
    num_trees_nu: int = 35,
    sigma_mu_total: float = 0.70,
    sigma_nu_total: float = 0.25,
    has_forecast: bool = True,
    use_forests: bool = True,
    use_factor_loadings: bool = True,
    use_spatial_gp: bool = True,
    baseline_untreated_only: bool = True,
    theta_f_clip: tuple[float, float] = (0.0, 1.5),
    theta_prior_mean: float = 1.0,
    theta_prior_sd: float = 0.5,
    sigma_prior_a: float = 1.5,
    sigma_prior_b: float = 0.5,
    sf_prior_a: float = 2.0,
    sf_prior_b: float = 0.2,
    intercept_prior_sd: float = 10.0,
    tau_prior_mean: float = 0.0,
    tau_prior_sd: float | None = None,
    tau_clip: tuple[float, float] = (-float("inf"), float("inf")),
    att_unit_weights: np.ndarray | None = None,
    tau_prior_family: str = "normal",
    tau_prior_params: tuple[float, float] | None = None,
    alpha_nu: float = 0.25,
    beta_nu: float = 3.0,
    nu_working_zeta_weight: float = 1.0,
    learn_sigma_nu: bool = False,
    sigma_nu_prior_scale: float = 0.5,
    treated_extra_variance: bool = False,
) -> RankKModeState:
    """Initialize a single-chain :class:`RankKModeState`.

    The keyword defaults of the last block (``att_unit_weights=None``,
    ``tau_prior_family="normal"``, ``alpha_nu=0.25``, ``beta_nu=3.0``,
    ``nu_working_zeta_weight=1.0``, ``learn_sigma_nu=False``,
    ``treated_extra_variance=False``) reproduce the paper's sampler exactly.

    ``att_centering`` (re-centering the treatment forest with the ATT weights
    ``att_unit_weights`` so that ``tau_intercept`` is the ATT) is switched on
    automatically, as a static flag, whenever a prior or a truncation acts on the
    intercept: ``tau_prior_sd`` is given, ``tau_clip`` is finite on either side,
    or ``tau_prior_family`` is not ``"normal"``.
    """
    family_code = {"normal": 0, "halfstudentt": 1, "gamma_spike": 2}.get(str(tau_prior_family))
    if family_code is None:
        raise ValueError(
            f"tau_prior_family must be 'normal', 'halfstudentt' or 'gamma_spike'; got {tau_prior_family!r}"
        )
    if family_code != 0 and tau_prior_sd is None:
        raise ValueError(f"tau_prior_family={tau_prior_family!r} needs tau_prior_sd (its scale)")
    if tau_prior_params is None:
        tau_prior_params = {0: (0.0, 0.0), 1: (2.5, 0.5), 2: (0.40, 1.60)}[family_code]
    att_centering = bool(
        tau_prior_sd is not None
        or math.isfinite(float(tau_clip[0]))
        or math.isfinite(float(tau_clip[1]))
        or family_code != 0
    )
    dy32 = jnp.asarray(dy_k, dtype=jnp.float32)
    df32 = jnp.asarray(df_k, dtype=jnp.float32)
    z32 = jnp.asarray(z_vec, dtype=jnp.float32)
    w32 = jnp.asarray(w_unit, dtype=jnp.float32)
    lam32 = jnp.asarray(lambda_mat, dtype=jnp.float32)

    N = int(dy32.shape[0])
    n0 = int(ctrl_idx.size)
    n1 = int(trt_idx.size)
    q = int(lam32.shape[1])
    w_np = np.asarray(w_unit, dtype=np.float64)
    n0_dof = int(np.sum(w_np[np.asarray(ctrl_idx, dtype=np.int64)] > 0.0))
    n1_dof = int(np.sum(w_np[np.asarray(trt_idx, dtype=np.int64)] > 0.0))
    n_all_dof = int(np.sum(w_np > 0.0))

    ctrl_mask = z32 == 0.0
    trt_mask = z32 == 1.0
    base_mask = ctrl_mask if baseline_untreated_only else jnp.ones(N, dtype=jnp.bool_)

    if att_unit_weights is None:
        w_att_np = np.where(np.asarray(z_vec, dtype=np.float64) == 1.0, np.asarray(w_unit, dtype=np.float64), 0.0)
    else:
        w_att_np = np.where(
            np.asarray(z_vec, dtype=np.float64) == 1.0, np.asarray(att_unit_weights, dtype=np.float64), 0.0
        )
    w_att_np = w_att_np / max(float(np.sum(w_att_np)), 1e-12)
    w_att32 = jnp.asarray(w_att_np, dtype=jnp.float32)

    p_nonterminal_mu = make_p_nonterminal(6, 0.95, 2.0)
    p_nonterminal_nu = make_p_nonterminal(6, float(alpha_nu), float(beta_nu))

    sigma_mu = float(sigma_mu_total) / math.sqrt(max(int(num_trees_mu), 1))
    sigma_nu = float(sigma_nu_total) / math.sqrt(max(int(num_trees_nu), 1))
    leaf_prior_cov_inv_mu = jnp.array(1.0 / (sigma_mu ** 2), dtype=jnp.float32)
    leaf_prior_cov_inv_nu = jnp.array(1.0 / (sigma_nu ** 2), dtype=jnp.float32)

    filter_splitless_mu = int(jnp.sum(max_split_mu == 0))
    filter_splitless_nu = int(jnp.sum(max_split_nu == 0))

    # Initial CUPED slope and residual on never-treated controls (strictly pre/control-only)
    c_idx = jnp.asarray(ctrl_idx, dtype=jnp.int32)
    if has_forecast:
        df_c = df32[c_idx] - jnp.mean(df32[c_idx])
        dy_c = dy32[c_idx] - jnp.mean(dy32[c_idx])
        denom = jnp.sum(w32[c_idx] * jnp.square(df_c)) + 1e-4
        num = jnp.sum(w32[c_idx] * df_c * dy_c)
        theta_init = jnp.clip(num / denom, float(theta_f_clip[0]), float(theta_f_clip[1])).astype(jnp.float32)
    else:
        theta_init = jnp.float32(0.0)

    resid_init = (dy32 - theta_init * df32).astype(jnp.float32)
    sigma2_init = jnp.maximum(jnp.var(resid_init[c_idx]), 1e-3).astype(jnp.float32)
    # Deterministic unit-scale reference vector in [-2, 2] so bartz sets
    # resid_unit = 4.0 and leaf_unit = 1.0 identically regardless of dy_k or N_0/N.
    y_ref = jnp.linspace(-2.0, 2.0, N, dtype=jnp.float32)

    # Observation error_scale = 1 / sqrt(w_unit)
    sqrt_w = jnp.sqrt(jnp.maximum(w32, 1e-8))
    mask_mu = base_mask
    err_scale_mu = jnp.where(mask_mu, jnp.reciprocal(sqrt_w), 1.0)
    err_scale_nu = jnp.where(trt_mask, jnp.reciprocal(sqrt_w), 1.0)

    def _err_cov() -> Wishart:
        return Wishart(nu=jnp.array(3.0), rate=jnp.array(3.0), value=jnp.array(1.0))

    state_mu = init(
        X=jnp.copy(X_binned),
        y=jnp.copy(y_ref),
        outcome_type="continuous",
        offset=0.0,
        max_split=jnp.copy(max_split_mu),
        num_trees=int(num_trees_mu),
        p_nonterminal=p_nonterminal_mu,
        leaf_prior_cov_inv=leaf_prior_cov_inv_mu,
        error_cov_inv=_err_cov(),
        error_scale=err_scale_mu,
        missing=~mask_mu,
        min_points_per_leaf=5,
        filter_splitless_vars=filter_splitless_mu,
        num_chains=None,
    )
    state_nu = init(
        X=jnp.copy(X_binned),
        y=jnp.copy(y_ref),
        outcome_type="continuous",
        offset=0.0,
        max_split=jnp.copy(max_split_nu),
        num_trees=int(num_trees_nu),
        p_nonterminal=p_nonterminal_nu,
        leaf_prior_cov_inv=leaf_prior_cov_inv_nu,
        error_cov_inv=_err_cov(),
        error_scale=err_scale_nu,
        missing=~trt_mask,
        min_points_per_leaf=5,
        filter_splitless_vars=filter_splitless_nu,
        num_chains=None,
    )

    def _freeze_var(st: State) -> State:
        frozen = eqx.tree_at(
            lambda w: (w.nu, w.rate), st.error_cov_inv, (None, None), is_leaf=lambda x: x is None
        )
        return eqx.tree_at(lambda s: s.error_cov_inv, st, frozen)

    state_mu = _freeze_var(state_mu)
    state_nu = _freeze_var(state_nu)
    eff_tau_prior_sd = (
        float(tau_prior_sd)
        if tau_prior_sd is not None
        else (10.0 if use_forests else float(sigma_nu_total))
    )
    # Half-Cauchy prior scale of the per-leaf SD of the treatment forest: a total
    # scale ``sigma_nu_prior_scale`` spread over ``num_trees_nu`` trees (data units).
    sigma_nu_prior_leaf = float(sigma_nu_prior_scale) / math.sqrt(max(int(num_trees_nu), 1))

    return RankKModeState(
        _chain_anchor=state_mu._chain_anchor,
        X=state_mu.X,
        y=state_mu.y,
        z=state_mu.z,
        binary_indices=state_mu.binary_indices,
        resid=resid_init,
        resid_unit=state_mu.resid_unit,
        resid_eff_scale=state_mu.resid_eff_scale,
        resid_inexact_integral=state_mu.resid_inexact_integral,
        error_cov_inv=state_mu.error_cov_inv,
        error_scale=state_mu.error_scale,
        prec_scale=state_mu.prec_scale,
        inv_sdev_scale=state_mu.inv_sdev_scale,
        inv_sdev_unit=state_mu.inv_sdev_unit,
        n_non_missing=state_mu.n_non_missing,
        sum_diag_prec_scale=state_mu.sum_diag_prec_scale,
        forest=state_mu.forest,
        config=state_mu.config,
        forest_nu=state_nu.forest,
        resid_nu=state_nu.resid,
        prec_scale_nu=state_nu.prec_scale,
        inv_sdev_scale_nu=state_nu.inv_sdev_scale,
        inv_sdev_unit_nu=state_nu.inv_sdev_unit,
        resid_unit_nu=state_nu.resid_unit,
        leaf_prior_cov_inv_nu=state_nu.forest.leaf_prior_cov_inv,
        resid_eff_scale_nu=state_nu.resid_eff_scale,
        resid_inexact_integral_nu=state_nu.resid_inexact_integral,
        mu_fit=jnp.zeros(N, dtype=jnp.float32),
        nu_fit=jnp.zeros(N, dtype=jnp.float32),
        intercept=jnp.float32(0.0),
        tau_intercept=jnp.float32(0.0),
        theta_f=theta_init,
        f_loadings=jnp.zeros(q, dtype=jnp.float32),
        s_f2=jnp.float32(0.1),
        zeta=jnp.zeros(N, dtype=jnp.float32),
        s_zeta2=jnp.float32(0.25 * sigma2_init),
        ell_idx=jnp.int32(min(2, 0)),
        sigma2=sigma2_init,
        z_vec=z32,
        w_unit=w32,
        df_vec=df32,
        df_point=df32,
        lambda_mat=lam32,
        ctrl_mask=ctrl_mask,
        trt_mask=trt_mask,
        ctrl_idx=jnp.asarray(ctrl_idx, dtype=jnp.int32),
        trt_idx=jnp.asarray(trt_idx, dtype=jnp.int32),
        N_units=N,
        N_ctrl=n0,
        N_trt=n1,
        N_ctrl_dof=n0_dof,
        N_trt_dof=n1_dof,
        N_all_dof=n_all_dof,
        q_factors=q,
        has_forecast=bool(has_forecast),
        use_forests=bool(use_forests),
        use_factor_loadings=bool(use_factor_loadings and q > 0),
        use_spatial_gp=bool(use_spatial_gp),
        baseline_untreated_only=bool(baseline_untreated_only),
        theta_min=jnp.float32(theta_f_clip[0]),
        theta_max=jnp.float32(theta_f_clip[1]),
        theta_prior_mean=jnp.float32(theta_prior_mean),
        theta_prior_sd=jnp.float32(theta_prior_sd),
        sigma_prior_a=jnp.float32(sigma_prior_a),
        sigma_prior_b=jnp.float32(sigma_prior_b),
        sf_prior_a=jnp.float32(sf_prior_a),
        sf_prior_b=jnp.float32(sf_prior_b),
        intercept_prior_sd=jnp.float32(intercept_prior_sd),
        tau_prior_mean=jnp.float32(tau_prior_mean),
        tau_prior_sd=jnp.float32(eff_tau_prior_sd),
        tau_min=jnp.float32(tau_clip[0]),
        tau_max=jnp.float32(tau_clip[1]),
        w_att=w_att32,
        tau_prior_p1=jnp.float32(tau_prior_params[0]),
        tau_prior_p2=jnp.float32(tau_prior_params[1]),
        nu_working_zeta_weight=jnp.float32(nu_working_zeta_weight),
        sigma_nu_prior_scale2=jnp.float32(sigma_nu_prior_leaf**2),
        s_eta2=jnp.float32(0.0),
        xi_nu=jnp.float32(1.0),
        att_centering=att_centering,
        tau_prior_family=int(family_code),
        treated_extra_variance=bool(treated_extra_variance),
        learn_sigma_nu=bool(learn_sigma_nu),
    )


@enable_x64()
def _subspace_mu_zeta_transfer(
    key: jax.Array,
    forest_mu: Forest,
    mu_fit: jax.Array,
    intercept: jax.Array,
    R: jax.Array,
    w_unit: jax.Array,
    ctrl_mask: jax.Array,
    sigma2: jax.Array,
    intercept_prior_sd: jax.Array,
) -> tuple[Forest, jax.Array, jax.Array, jax.Array]:
    """Exact conjugate subspace move shifting ``mu_k`` leaf means vs. ``intercept``.

    Along ``delta in R``, shifting every active leaf of ``mu_k`` by ``delta / m_mu``
    and ``intercept`` by ``-delta`` leaves the sum ``mu_k + intercept`` invariant in
    the likelihood and refreshes the leaf mean prior against the intercept prior in
    one exact Gaussian draw.
    """
    f64 = jnp.float64
    ids_mu = apply_moves_to_leaf_indices(
        forest_mu.leaf_indices, forest_mu.to_prune, forest_mu.move_node
    )
    _, _, is_leaf_mu = compute_active_leaf_stats(forest_mu.split_tree, forest_mu.leaf_tree)
    m_mu = jnp.float64(forest_mu.leaf_tree.shape[0])
    leaf_unit = jnp.asarray(forest_mu.leaf_unit, dtype=f64)
    leaf_prec = jnp.asarray(forest_mu.leaf_prior_cov_inv, dtype=f64) / jnp.maximum(jnp.square(leaf_unit), 1e-12)

    # Active leaves in data units: leaf_data = leaf_unit * leaf_tree
    leaf_data = jnp.where(is_leaf_mu, leaf_unit * jnp.asarray(forest_mu.leaf_tree, dtype=f64), 0.0)
    n_active = jnp.sum(is_leaf_mu.astype(f64))
    sum_leaves = jnp.sum(leaf_data)

    # Let leaf_data_new = leaf_data + (delta / m_mu) * is_leaf_mu, intercept_new = intercept - delta
    # Prior quadratic form in delta:
    # leaf prior: 0.5 * leaf_prec * sum_{t,l in active} (leaf_data_{tl} + delta / m_mu)^2
    # intercept prior: 0.5 / s0^2 * (intercept - delta)^2
    s0_inv2 = 1.0 / jnp.maximum(jnp.square(jnp.asarray(intercept_prior_sd, dtype=f64)), 1e-12)
    prec_delta = leaf_prec * (n_active / (m_mu ** 2)) + s0_inv2
    rhs_delta = -leaf_prec * (sum_leaves / m_mu) + s0_inv2 * jnp.asarray(intercept, dtype=f64)

    # Exact induced change in mu_fit on each unit (equals delta if all m_mu trees have active leaves)
    active_per_unit = jnp.sum(
        jax.vmap(lambda act_t, idx: act_t[idx].astype(f64))(is_leaf_mu, ids_mu),
        axis=0,
    )
    h_unit = active_per_unit / m_mu - 1.0  # net shift in (mu_fit + intercept)
    w_ctrl = jnp.where(ctrl_mask, jnp.asarray(w_unit, dtype=f64), 0.0)
    sig2_64 = jnp.maximum(jnp.asarray(sigma2, dtype=f64), 1e-12)
    prec_delta = prec_delta + jnp.sum(w_ctrl * jnp.square(h_unit)) / sig2_64
    rhs_delta = rhs_delta + jnp.sum(w_ctrl * h_unit * jnp.asarray(R, dtype=f64)) / sig2_64

    sd_delta = 1.0 / jnp.sqrt(jnp.maximum(prec_delta, 1e-12))
    mean_delta = rhs_delta / jnp.maximum(prec_delta, 1e-12)
    delta = mean_delta + sd_delta * random.normal(key, (), dtype=jnp.float32).astype(f64)

    leaf_tree_shift = jnp.where(
        is_leaf_mu,
        (delta / (m_mu * jnp.maximum(leaf_unit, 1e-12))).astype(jnp.float32),
        0.0,
    )
    forest_mu_new = eqx.tree_at(
        lambda f: f.leaf_tree, forest_mu, (forest_mu.leaf_tree + leaf_tree_shift).astype(jnp.float32)
    )
    mu_fit_new = current_forest_fit(forest_mu_new)
    intercept_new = (jnp.asarray(intercept, dtype=f64) - delta).astype(jnp.float32)
    R_new = (R + (mu_fit - mu_fit_new) + (intercept - intercept_new)).astype(jnp.float32)
    return forest_mu_new, mu_fit_new, intercept_new, R_new


@enable_x64()
def _orthogonalize_mu_fit(
    mu_raw: jax.Array,
    state: RankKModeState,
    base_mask: jax.Array,
) -> jax.Array:
    """Project ``mu_raw`` orthogonal to the linear fixed-effect columns ``(1, df_point)`` on ``base_mask``.

    Enforces Frisch-Waugh-Lovell orthogonality between the nonparametric
    prognostic forest ``mu_k`` and the linear point-forecast baseline
    ``intercept + theta_k * df_point`` on ``base_mask`` without projecting out
    between-draw forecast perturbations ``df_vec - df_point``.
    """
    use_orth = state.has_forecast
    if not use_orth:
        return mu_raw
    f64 = jnp.float64
    w_orth = jnp.where(base_mask, jnp.asarray(state.w_unit, dtype=f64), 0.0)
    cols = [
        jnp.ones(state.N_units, dtype=f64),
        jnp.asarray(state.df_point, dtype=f64),
    ]
    B = jnp.column_stack(cols)
    Bw = B * w_orth[:, None]
    gram = B.T @ Bw + 1e-6 * jnp.eye(B.shape[1], dtype=f64)
    rhs = Bw.T @ jnp.asarray(mu_raw, dtype=f64)
    coef = jnp.linalg.solve(gram, rhs)
    return (jnp.asarray(mu_raw, dtype=f64) - B @ coef).astype(jnp.float32)


@enable_x64()
def _update_conjugate_blocks(
    key: jax.Array,
    state: RankKModeState,
    R: jax.Array,
    mu_fit: jax.Array,
    forest_mu: Forest,
    spatial_cache: SpatialKernelCache,
) -> tuple[
    Forest,
    jax.Array,
    jax.Array,
    jax.Array,
    jax.Array,
    jax.Array,
    jax.Array,
    jax.Array,
    jax.Array,
    jax.Array,
    jax.Array,
]:
    """Update ``(intercept, theta_k, f_k, s_f^2, ell_k, s_zeta^2, sigma_k^2, zeta_k)`` and ``R``.

    Implements a Partially Collapsed Gibbs (PCG; Van Dyk & Park 2008) update:
    ``zeta_{k,I_0}`` is marginalized out of the ``(intercept, theta_k, f_k)`` updates
    using the precomputed spectral decomposition of ``W_0^{1/2} C_{00}(ell_g) W_0^{1/2}``,
    and then ``(ell_k, s_zeta^2, sigma_k^2, zeta_{k,I_0}, zeta_{k,I_1})`` are drawn
    jointly from their conditional distribution given ``(intercept, theta_k, f_k, mu_k)``.
    With forests (``state.use_forests``) the returned field and residual hold
    ``zeta_{k,I_1}`` at its kriging mean ``A_g E[zeta_{k,I_0} | control residuals]``
    and the returned ``delta_zeta1`` (kriging draw minus that mean) enters only the
    treatment-intercept draw of :func:`rank_k_mode_single_step`, one of the
    approximations of the forest path (see the module docstring); without forests the
    kriging draw itself enters the residual and ``delta_zeta1`` is zero.
    """
    f64 = jnp.float64
    k_sub, k_int, k_th, k_f, k_sf, k_grid, k_sig, k_z0, k_z1 = random.split(key, 9)

    w_all = jnp.asarray(state.w_unit, dtype=f64)
    base_mask = state.ctrl_mask if state.baseline_untreated_only else jnp.ones_like(state.ctrl_mask)
    w_base = jnp.where(base_mask, w_all, 0.0)
    sig2 = jnp.maximum(jnp.asarray(state.sigma2, dtype=f64), 1e-10)
    ctrl_idx = state.ctrl_idx
    trt_idx = state.trt_idx

    # 2a. Subspace transfer between mu_k leaves and intercept
    intercept = state.intercept
    if state.use_forests:
        forest_mu, mu_fit, intercept, R = _subspace_mu_zeta_transfer(
            k_sub,
            forest_mu,
            mu_fit,
            intercept,
            R,
            state.w_unit,
            base_mask,
            state.sigma2,
            state.intercept_prior_sd,
        )

    # Remove old zeta_k and anchor the control-only baseline updates on df_point so
    # between-draw forecast perturbations (df_vec - df_point) are not absorbed by
    # (intercept, theta_k, f_k, zeta_0) on I_0.
    zeta_old64 = jnp.asarray(state.zeta, dtype=f64)
    df_diff64 = (
        jnp.asarray(state.df_vec, dtype=f64) - jnp.asarray(state.df_point, dtype=f64)
        if state.has_forecast
        else jnp.zeros(state.N_units, dtype=f64)
    )
    R_no_zeta = (
        jnp.asarray(R, dtype=f64)
        + zeta_old64
        + jnp.asarray(state.theta_f, dtype=f64) * df_diff64
    )

    use_collapsed_omega = bool(
        state.use_spatial_gp and spatial_cache.enabled and state.baseline_untreated_only
    )
    if use_collapsed_omega:
        g_curr = jnp.clip(state.ell_idx, 0, spatial_cache.n_ell - 1)
        r_curr = jnp.maximum(jnp.asarray(state.s_zeta2, dtype=f64) / sig2, 0.0)
        U_curr = spatial_cache.eig_vecs_00[g_curr]  # (N0, N0)
        d_curr = spatial_cache.eig_vals_00[g_curr]  # (N0,)
        w_eff = 1.0 / (1.0 + r_curr * d_curr)       # (N0,)

        def _proj_base(v_full: jax.Array) -> jax.Array:
            v0 = v_full[ctrl_idx]
            if v0.ndim == 1:
                return U_curr.T @ (spatial_cache.sqrt_w0 * v0)
            return U_curr.T @ (spatial_cache.sqrt_w0[:, None] * v0)

    else:
        w_eff = w_base

        def _proj_base(v_full: jax.Array) -> jax.Array:
            return v_full

    # 2b. Intercept update on baseline cells (marginalized over zeta_{k, I_0})
    r_int = R_no_zeta + jnp.asarray(intercept, dtype=f64)
    p_ones = _proj_base(jnp.ones(state.N_units, dtype=f64))
    p_r_int = _proj_base(r_int)
    int_prior_var = jnp.maximum(jnp.square(jnp.asarray(state.intercept_prior_sd, dtype=f64)), 1e-12)
    prec_int = jnp.sum(w_eff * jnp.square(p_ones)) / sig2 + 1.0 / int_prior_var
    rhs_int = jnp.sum(w_eff * p_ones * p_r_int) / sig2
    int_new64 = rhs_int / prec_int + random.normal(k_int, (), dtype=jnp.float32).astype(f64) / jnp.sqrt(prec_int)
    R_no_zeta = r_int - int_new64
    intercept_new = int_new64.astype(jnp.float32)

    # 2c. Forecast slope theta_k update (truncated Normal on [theta_min, theta_max], marginalized over zeta_{k, I_0})
    if state.has_forecast:
        df64 = jnp.asarray(state.df_point, dtype=f64)
        r_th = R_no_zeta + jnp.asarray(state.theta_f, dtype=f64) * df64
        p_df = _proj_base(df64)
        p_r_th = _proj_base(r_th)
        th_prior_var = jnp.maximum(jnp.square(jnp.asarray(state.theta_prior_sd, dtype=f64)), 1e-12)
        prec_th = jnp.sum(w_eff * jnp.square(p_df)) / sig2 + 1.0 / th_prior_var
        rhs_th = (
            jnp.sum(w_eff * p_df * p_r_th) / sig2
            + jnp.asarray(state.theta_prior_mean, dtype=f64) / th_prior_var
        )
        mean_th = rhs_th / prec_th
        sd_th = 1.0 / jnp.sqrt(prec_th)
        th_new64 = _sample_truncated_normal_64(k_th, mean_th, sd_th, state.theta_min, state.theta_max)
        R_no_zeta = r_th - th_new64 * df64
        theta_new = th_new64.astype(jnp.float32)
    else:
        th_new64 = jnp.float64(0.0)
        theta_new = jnp.float32(0.0)

    # 2d. Pre-anchor factor loading shocks f_k and variance s_f^2 (marginalized over zeta_{k, I_0})
    if state.use_factor_loadings:
        lam64 = jnp.asarray(state.lambda_mat, dtype=f64)
        f_old64 = jnp.asarray(state.f_loadings, dtype=f64)
        r_f = R_no_zeta + lam64 @ f_old64
        sf2_64 = jnp.maximum(jnp.asarray(state.s_f2, dtype=f64), 1e-8)
        p_lam = _proj_base(lam64)  # (N0, q) or (N, q)
        p_r_f = _proj_base(r_f)
        p_lam_w = p_lam * w_eff[:, None]
        Q_f = (p_lam.T @ p_lam_w) / sig2 + jnp.eye(state.q_factors, dtype=f64) / sf2_64
        Q_f = 0.5 * (Q_f + Q_f.T)
        b_f = (p_lam_w.T @ p_r_f) / sig2
        L_f = jnp.linalg.cholesky(Q_f)
        z_f = random.normal(k_f, (state.q_factors,), dtype=jnp.float32).astype(f64)
        f_new64 = (
            jax.scipy.linalg.cho_solve((L_f, True), b_f)
            + jax.scipy.linalg.solve_triangular(L_f.T, z_f, lower=False)
        )
        R_no_zeta = r_f - lam64 @ f_new64
        sf2_new64 = _sample_inv_gamma_64(
            k_sf,
            jnp.asarray(state.sf_prior_a, dtype=f64) + 0.5 * state.q_factors,
            jnp.asarray(state.sf_prior_b, dtype=f64) + 0.5 * jnp.sum(jnp.square(f_new64)),
        )
        f_new = f_new64.astype(jnp.float32)
        s_f2_new = jnp.clip(sf2_new64, 1e-6, 100.0).astype(jnp.float32)
    else:
        f_new = state.f_loadings
        s_f2_new = state.s_f2

    # 2e. Spatial / cluster random field zeta_k and variance sigma_k^2 (collapsed over zeta_{k,I_0})
    r_z = R_no_zeta
    r0 = r_z[ctrl_idx]

    if state.use_spatial_gp and spatial_cache.enabled:
        # Weighted control residual u_0 = W_0^{1/2} r_0
        u0 = spatial_cache.sqrt_w0 * r0  # (N0,)
        # Project onto eigenvectors for every lengthscale g: v_g = U_g^T u_0, shape (G, N0)
        v_g = jnp.einsum("gij,i->gj", spatial_cache.eig_vecs_00, u0)
        v_g2 = jnp.square(v_g)  # (G, N0)
        d_g = spatial_cache.eig_vals_00  # (G, N0)

        one_plus_r = 1.0 + spatial_cache.ratio_grid[None, :, None]  # (1, H, 1)
        denom_gh = (1.0 + spatial_cache.ratio_grid[None, :, None] * d_g[:, None, :]) / one_plus_r  # (G, H, N0)
        logdet_gh = jnp.sum(jnp.log(jnp.maximum(denom_gh, 1e-12)), axis=-1)         # (G, H)
        quad_gh = jnp.sum(v_g2[:, None, :] / denom_gh, axis=-1)                     # (G, H)

        a_post = jnp.asarray(state.sigma_prior_a, dtype=f64) + 0.5 * state.N_ctrl_dof
        b_post_gh = jnp.asarray(state.sigma_prior_b, dtype=f64) + 0.5 * quad_gh
        log_post_gh = (
            spatial_cache.ratio_log_prior[None, :]
            - 0.5 * logdet_gh
            - a_post * jnp.log(jnp.maximum(b_post_gh, 1e-12))
        )
        flat_logits = log_post_gh.reshape(-1)
        flat_idx = random.categorical(k_grid, flat_logits)
        g_idx = (flat_idx // spatial_cache.n_ratio).astype(jnp.int32)
        h_idx = (flat_idx % spatial_cache.n_ratio).astype(jnp.int32)

        r_sel = spatial_cache.ratio_grid[h_idx]
        d_sel = d_g[g_idx]        # (N0,)
        v_sel = v_g[g_idx]        # (N0,)
        U_sel = spatial_cache.eig_vecs_00[g_idx]  # (N0, N0)
        if state.use_forests:
            b_sel = b_post_gh[g_idx, h_idx]
            sig2_tot_new64 = _sample_inv_gamma_64(k_sig, a_post, b_sel)
            sig2_new64 = jnp.maximum(sig2_tot_new64 / (1.0 + r_sel), 1e-8)
            s_zeta2_new64 = (r_sel / (1.0 + r_sel)) * sig2_tot_new64
        else:
            b_raw_sel = jnp.asarray(state.sigma_prior_b, dtype=f64) + 0.5 * jnp.sum(
                v_g2[g_idx] / (1.0 + r_sel * d_sel)
            )
            sig2_new64 = jnp.maximum(_sample_inv_gamma_64(k_sig, a_post, b_raw_sel), 1e-8)
            s_zeta2_new64 = r_sel * sig2_new64

        # Draw zeta_{k, I_0} | r_0, g_idx, r_sel, sig2_new64 in the spectral basis:
        shrink = (r_sel * d_sel) / (1.0 + r_sel * d_sel)
        alpha_mean = shrink * v_sel
        alpha_sd = jnp.sqrt(jnp.maximum(sig2_new64 * shrink, 0.0))
        alpha_draw = alpha_mean + alpha_sd * random.normal(k_z0, (state.N_ctrl,), dtype=jnp.float32).astype(f64)
        zeta0_mean = spatial_cache.inv_sqrt_w0 * (U_sel @ alpha_mean)
        zeta0_draw = spatial_cache.inv_sqrt_w0 * (U_sel @ alpha_draw)

        # Draw zeta_{k, I_1} | zeta_{k, I_0}, g_idx, s_zeta2_new64 via exact kriging:
        A_sel = spatial_cache.krig_mat_10[g_idx]      # (N1, N0)
        L11_sel = spatial_cache.schur_chol_11[g_idx]  # (N1, N1)
        zeta1_mean_clean = A_sel @ zeta0_mean
        zeta1_mean = A_sel @ zeta0_draw
        zeta1_draw = zeta1_mean + jnp.sqrt(jnp.maximum(s_zeta2_new64, 0.0)) * (
            L11_sel @ random.normal(k_z1, (state.N_trt,), dtype=jnp.float32).astype(f64)
        )
        delta_zeta1 = zeta1_draw - zeta1_mean_clean

        if state.use_forests:
            # Approximation of the forest path: the residual keeps zeta_{k,I_1} at its
            # kriging mean; delta_zeta1 enters only the treatment-intercept draw.
            zeta_new64 = (
                jnp.zeros(state.N_units, dtype=f64)
                .at[ctrl_idx].set(zeta0_draw)
                .at[trt_idx].set(zeta1_mean_clean)
            )
        else:
            zeta_new64 = (
                jnp.zeros(state.N_units, dtype=f64)
                .at[ctrl_idx].set(zeta0_draw)
                .at[trt_idx].set(zeta1_draw)
            )
            delta_zeta1 = jnp.zeros(state.N_trt, dtype=f64)
        R64 = r_z - zeta_new64 - th_new64 * df_diff64
        ell_idx_new = g_idx
    else:
        zeta_new64 = jnp.zeros(state.N_units, dtype=f64)
        delta_zeta1 = jnp.zeros(state.N_trt, dtype=f64)
        if state.baseline_untreated_only:
            n_base = jnp.float64(state.N_ctrl_dof)  # positive-weight never-treated units (all of them on complete panels)
            ssr_base = jnp.sum(w_base * jnp.square(r_z))
            sig2_new64 = _sample_inv_gamma_64(
                k_sig,
                jnp.asarray(state.sigma_prior_a, dtype=f64) + 0.5 * n_base,
                jnp.asarray(state.sigma_prior_b, dtype=f64) + 0.5 * ssr_base,
            )
        else:
            # When baseline_cells == "all", sigma_k^2 is updated once at the end of the sweep after tau_k
            sig2_new64 = sig2
        R64 = r_z - th_new64 * df_diff64
        s_zeta2_new64 = jnp.float64(0.0)
        ell_idx_new = jnp.int32(0)

    return (
        forest_mu,
        mu_fit,
        intercept_new,
        theta_new,
        f_new,
        s_f2_new,
        zeta_new64.astype(jnp.float32),
        s_zeta2_new64.astype(jnp.float32),
        ell_idx_new,
        jnp.clip(sig2_new64, 1e-8, 1e6).astype(jnp.float32),
        R64.astype(jnp.float32),
        delta_zeta1.astype(jnp.float32),
    )


@enable_x64()
@jax.jit
@enable_x64(False)
def rank_k_mode_single_step(
    key: jax.Array,
    state: RankKModeState,
    spatial_cache: SpatialKernelCache,
) -> RankKModeState:
    """Execute one partially collapsed, approximate Gibbs sweep for one mode ``k`` and one chain.

    Sweep order:
      1. Prognostic forest ``mu_k`` on ``base_mask`` (``I_0`` under
         ``baseline_cells='untreated'``), with the diagonal working noise
         variance ``sigma_k^2 + s_zeta^2 / 2``.
      2. Conjugate blocks ``(intercept, theta_k, f_k, s_f^2, ell_k, s_zeta^2,
         sigma_k^2, zeta_{k,I_0}, zeta_{k,I_1})`` on ``I_0`` (see
         :func:`_update_conjugate_blocks`), then the residual update on all units
         immediately before Step 3. With forests the treated residual holds
         ``zeta_{k,I_1}`` at its kriging mean ``A_g E[zeta_{k,I_0} | control
         residuals]``; the deviation ``delta_zeta1`` of the kriging draw from that
         mean enters only the scalar treatment-intercept draw of Step 3.
      3. Treatment forest ``tau_k`` on treated units ``I_1``, with the diagonal
         working noise variance ``sigma_k^2 + s_zeta^2``, followed by the
         treatment-intercept draw.

    The working variances stand in for the spatially correlated likelihood, so
    with forests the sweep is an approximate Gibbs step and the chain is not an
    exact posterior sampler. With ``state.use_forests=False`` (the
    linear-Gaussian verification limit used in the tests) both forests are
    constants and the kriging draw of ``zeta_{k,I_1}`` enters the residual.
    """
    k_mu, k_conj, k_nu = random.split(key, 3)
    R = state.resid
    # Working likelihood of the mu_k step (approximation): a diagonal noise variance
    # inflated by half the spatial variance replaces the correlated field.
    sigma2_tree = state.sigma2 + 0.5 * state.s_zeta2
    sqrt_w = jnp.sqrt(jnp.maximum(state.w_unit, 1e-8))
    base_mask = state.ctrl_mask if state.baseline_untreated_only else jnp.ones_like(state.ctrl_mask)
    mask_mu = base_mask

    # ------------------------------------------------------------------
    # 1. Prognostic forest mu_k
    # ------------------------------------------------------------------
    if state.use_forests:
        w_mu = jnp.where(mask_mu, sqrt_w, 0.0)
        attrs_mu, _ = compute_treatment_scale_attrs(w_mu, mask_mu, eps=1e-6)
        mu_raw_old = current_forest_fit(state.forest)
        R_for_mu = (
            R
            + (state.mu_fit - mu_raw_old)
            + state.theta_f * (state.df_vec - state.df_point)
        )
        R_mu_for_bartz = jnp.where(mask_mu, R_for_mu, 0.0)
        view_mu = State(
            _chain_anchor=state._chain_anchor,
            X=state.X,
            y=state.y,
            z=None,
            binary_indices=None,
            resid=_load(R_mu_for_bartz, state.resid_unit),
            resid_unit=state.resid_unit,
            resid_eff_scale=state.resid_eff_scale,
            resid_inexact_integral=state.resid_inexact_integral,
            error_cov_inv=Wishart(nu=None, rate=None, value=jnp.reciprocal(sigma2_tree)),
            error_scale=None,
            prec_scale=attrs_mu.prec_scale,
            inv_sdev_scale=attrs_mu.inv_sdev_scale,
            inv_sdev_unit=attrs_mu.inv_sdev_unit,
            n_non_missing=attrs_mu.n_non_missing,
            sum_diag_prec_scale=attrs_mu.sum_diag_prec_scale,
            forest=state.forest,
            config=state.config,
        )
        view_mu = refresh_prec_tree(view_mu)
        view_mu = bartz_step(k_mu, view_mu)
        mu_fit = current_forest_fit(view_mu.forest)
        forest_mu = view_mu.forest
        R = R - (mu_fit - state.mu_fit)
    else:
        forest_mu = state.forest
        mu_fit = state.mu_fit

    # ------------------------------------------------------------------
    # 2. Control-only conjugate blocks + conditional kriging to I_1
    # ------------------------------------------------------------------
    (
        forest_mu,
        mu_fit,
        intercept_new,
        theta_new,
        f_new,
        s_f2_new,
        zeta_new,
        s_zeta2_new,
        ell_idx_new,
        sigma2_new,
        R,
        delta_zeta1,
    ) = _update_conjugate_blocks(k_conj, state, R, mu_fit, forest_mu, spatial_cache)

    # ------------------------------------------------------------------
    # 3. Treatment forest tau_k on treated units I_1
    # ------------------------------------------------------------------
    if state.use_forests:
        k_nu_tree, k_tau_int = random.split(k_nu)
        if state.treated_extra_variance:
            k_tau_int, k_eta = random.split(k_tau_int)
        if state.learn_sigma_nu:
            k_tau_int, k_leaf = random.split(k_tau_int)
        # Working likelihood of the tau_k step (approximation): zeta_{k,I_1} sits at its
        # kriging mean in R and the spatial variance is added to a diagonal noise variance
        # (weighted by ``nu_working_zeta_weight``, 1.0 in the paper), plus the treated-arm
        # excess variance s_eta^2 when ``treated_extra_variance`` is on.
        sigma2_nu_tree = sigma2_new + state.nu_working_zeta_weight * s_zeta2_new
        if state.treated_extra_variance:
            sigma2_nu_tree = sigma2_nu_tree + state.s_eta2
        w_trt = jnp.where(state.trt_mask, sqrt_w, 0.0)
        attrs_nu, mask_nu = compute_treatment_scale_attrs(w_trt, state.trt_mask, eps=1e-6)
        nu_raw_old = current_forest_fit(state.forest_nu)
        R_for_nu = R + state.z_vec * (state.nu_fit - nu_raw_old)
        resid_nu_data = jnp.where(mask_nu, R_for_nu, 0.0)
        view_nu = State(
            _chain_anchor=state._chain_anchor,
            X=view_mu.X,
            y=view_mu.y,
            z=None,
            binary_indices=None,
            resid=_load(resid_nu_data, state.resid_unit_nu),
            resid_unit=state.resid_unit_nu,
            resid_eff_scale=state.resid_eff_scale_nu,
            resid_inexact_integral=state.resid_inexact_integral_nu,
            error_cov_inv=Wishart(nu=None, rate=None, value=jnp.reciprocal(sigma2_nu_tree)),
            error_scale=None,
            prec_scale=attrs_nu.prec_scale,
            inv_sdev_scale=attrs_nu.inv_sdev_scale,
            inv_sdev_unit=attrs_nu.inv_sdev_unit,
            n_non_missing=attrs_nu.n_non_missing,
            sum_diag_prec_scale=attrs_nu.sum_diag_prec_scale,
            forest=state.forest_nu,
            config=view_mu.config,
        )
        view_nu = refresh_prec_tree(view_nu)
        nu_before = state.nu_fit
        view_nu = bartz_step(k_nu_tree, view_nu)
        nu_raw_new = current_forest_fit(view_nu.forest)
        forest_nu = view_nu.forest
        with enable_x64():
            f64 = jnp.float64
            w_trt64 = jnp.where(state.trt_mask, jnp.asarray(state.w_unit, dtype=f64), 0.0)
            w_trt_sum = jnp.maximum(jnp.sum(w_trt64), 1e-12)
            if state.att_centering:
                # Re-center the treatment forest with the ATT weights a_i (sum_i a_i = 1 on
                # I_1), so that tau_intercept == ATT at every sweep and the prior /
                # truncation on the intercept acts exactly on the ATT (Issue 1a).
                c_trt64 = jnp.where(state.trt_mask, jnp.asarray(state.w_att, dtype=f64), 0.0)
                nu_mean_trt = jnp.sum(c_trt64 * jnp.asarray(nu_raw_new, dtype=f64)) / jnp.maximum(
                    jnp.sum(c_trt64), 1e-12
                )
            else:
                nu_mean_trt = jnp.sum(w_trt64 * jnp.asarray(nu_raw_new, dtype=f64)) / w_trt_sum
            _, _, is_leaf_nu = compute_active_leaf_stats(forest_nu.split_tree, forest_nu.leaf_tree)
            m_nu = jnp.float64(forest_nu.leaf_tree.shape[0])
            leaf_unit_nu = jnp.asarray(forest_nu.leaf_unit, dtype=f64)
            leaf_shift_nu = jnp.where(
                is_leaf_nu,
                (-nu_mean_trt / (m_nu * jnp.maximum(leaf_unit_nu, 1e-12))).astype(jnp.float32),
                0.0,
            )
            forest_nu = eqx.tree_at(
                lambda f: f.leaf_tree, forest_nu, (forest_nu.leaf_tree + leaf_shift_nu).astype(jnp.float32)
            )
            nu_fit = current_forest_fit(forest_nu)
            R = R - state.z_vec * (nu_fit - nu_before)
            sig2_64 = jnp.maximum(jnp.asarray(sigma2_new, dtype=f64), 1e-10)
            if state.treated_extra_variance:
                # Treated-arm residual variance sigma_1^2 = sigma_0^2 + s_eta^2 (Issue 3c).
                sig2_trt_64 = sig2_64 + jnp.asarray(state.s_eta2, dtype=f64)
            else:
                sig2_trt_64 = sig2_64
            tau_clean_prev = (
                jnp.sum(w_trt64 * jnp.asarray(R, dtype=f64)) / w_trt_sum
                if False
                else jnp.asarray(state.tau_intercept, dtype=f64)
            )
            r_tau = jnp.asarray(R, dtype=f64) + jnp.asarray(state.z_vec, dtype=f64) * tau_clean_prev
            tau_prior_prec = 1.0 / jnp.maximum(jnp.square(jnp.asarray(state.tau_prior_sd, dtype=f64)), 1e-12)
            tau_prior_m64 = jnp.asarray(state.tau_prior_mean, dtype=f64)
            prec_tau = w_trt_sum / sig2_trt_64 + tau_prior_prec
            w1_vec = jnp.asarray(state.w_unit[state.trt_idx], dtype=f64)
            dz1_64 = jnp.asarray(delta_zeta1, dtype=f64)
            lik_sum = jnp.sum(w_trt64 * r_tau) - jnp.sum(w1_vec * dz1_64)
            rhs_tau = lik_sum / sig2_trt_64 + tau_prior_prec * tau_prior_m64
            mean_tau = rhs_tau / prec_tau
            sd_tau = 1.0 / jnp.sqrt(prec_tau)
            if state.tau_prior_family == 0:
                unconstrained_tau = mean_tau + random.normal(k_tau_int, (), dtype=jnp.float32).astype(f64) * sd_tau
                trunc_tau = _sample_truncated_normal_64(k_tau_int, mean_tau, sd_tau, state.tau_min, state.tau_max)
                has_tau_clip = jnp.isfinite(jnp.asarray(state.tau_min, dtype=f64)) | jnp.isfinite(
                    jnp.asarray(state.tau_max, dtype=f64)
                )
                tau_scalar = jnp.where(has_tau_clip, trunc_tau, unconstrained_tau)
            else:
                # Non-conjugate ATT prior (half-Student-t or Gamma spike): exact scalar slice
                # sampling of the intercept's full conditional (Gaussian likelihood x prior).
                prec_lik = w_trt_sum / sig2_trt_64
                mean_lik = lik_sum / w_trt_sum
                family = int(state.tau_prior_family)
                p_mean, p_sd, p1, p2 = state.tau_prior_mean, state.tau_prior_sd, state.tau_prior_p1, state.tau_prior_p2

                def _log_cond(x):
                    return -0.5 * prec_lik * jnp.square(x - mean_lik) + _log_prior_tau_64(
                        x, family, p_mean, p_sd, p1, p2
                    )

                lo_support = 0.0
                tau_lo = jnp.maximum(jnp.asarray(state.tau_min, dtype=f64), lo_support)
                tau_hi = jnp.asarray(state.tau_max, dtype=f64)
                x_start = jnp.clip(jnp.asarray(state.tau_intercept, dtype=f64), tau_lo + 1e-9, tau_hi)
                tau_scalar = _slice_sample_scalar_64(
                    k_tau_int, x_start, _log_cond, tau_lo, tau_hi, 3.0 / jnp.sqrt(prec_lik)
                )
            tau_intercept_new = tau_scalar.astype(jnp.float32)
            R = (r_tau - jnp.asarray(state.z_vec, dtype=f64) * tau_scalar).astype(jnp.float32)
            if state.treated_extra_variance:
                # Slice-sample the treated-arm excess variance s_eta^2 (half-Cauchy(0, sigma_0)
                # prior on s_eta) from the treated residuals, which hold zeta_{k,I_1} at its
                # kriging mean (delta_zeta1 removed, as in the intercept draw).
                r_trt = jnp.asarray(R, dtype=f64)[state.trt_idx] - dz1_64
                ssr_1 = jnp.sum(w1_vec * jnp.square(r_trt))
                n_1 = jnp.float64(state.N_trt_dof)
                sig0 = jnp.sqrt(sig2_64)

                def _log_cond_log_s2(u):  # u = log s_eta^2 (Jacobian included)
                    s2 = jnp.exp(u)
                    v = sig2_64 + s2
                    log_lik = -0.5 * n_1 * jnp.log(v) - 0.5 * ssr_1 / v
                    log_prior = -jnp.log1p(s2 / sig2_64) + 0.5 * u - jnp.log(sig0)  # half-Cauchy on s_eta, times ds/du
                    return log_lik + log_prior

                u0 = jnp.log(jnp.maximum(jnp.asarray(state.s_eta2, dtype=f64), 1e-8 * sig2_64))
                u_new = _slice_sample_scalar_64(k_eta, u0, _log_cond_log_s2, jnp.log(1e-8 * sig2_64), jnp.log(1e4 * sig2_64), 1.0)
                s_eta2_new = jnp.exp(u_new).astype(jnp.float32)
            else:
                s_eta2_new = state.s_eta2
            if state.learn_sigma_nu:
                # Half-Cauchy(0, A) prior on the per-leaf SD of the treatment forest via the
                # inverse-gamma parameter expansion (Issue 3b): sigma_leaf^2 | xi ~ IG(1/2, 1/xi),
                # xi ~ IG(1/2, 1/A^2).
                k_leaf_a, k_leaf_b = random.split(k_leaf)
                n_leaf, ss_leaf, _ = compute_active_leaf_stats(forest_nu.split_tree, forest_nu.leaf_tree)
                ss_leaf_data = jnp.asarray(ss_leaf, dtype=f64) * jnp.square(leaf_unit_nu)
                xi_64 = jnp.maximum(jnp.asarray(state.xi_nu, dtype=f64), 1e-30)
                sigma_leaf2 = _sample_inv_gamma_64(
                    k_leaf_a, 0.5 * (1.0 + jnp.asarray(n_leaf, dtype=f64)), 1.0 / xi_64 + 0.5 * ss_leaf_data
                )
                sigma_leaf2 = jnp.clip(sigma_leaf2, 1e-10, 1e6)
                a2 = jnp.maximum(jnp.asarray(state.sigma_nu_prior_scale2, dtype=f64), 1e-30)
                xi_new = _sample_inv_gamma_64(k_leaf_b, jnp.float64(1.0), 1.0 / a2 + 1.0 / sigma_leaf2)
                forest_nu = eqx.tree_at(
                    lambda f: f.leaf_prior_cov_inv,
                    forest_nu,
                    jnp.asarray(1.0 / sigma_leaf2, dtype=forest_nu.leaf_prior_cov_inv.dtype),
                )
                xi_nu_new = xi_new.astype(jnp.float32)
            else:
                xi_nu_new = state.xi_nu
    else:
        # Forests-off linear-Gaussian constant treatment effect draw on I_1
        with enable_x64():
            f64 = jnp.float64
            w_trt64 = jnp.where(state.trt_mask, jnp.asarray(state.w_unit, dtype=f64), 0.0)
            sig2_64 = jnp.maximum(jnp.asarray(sigma2_new, dtype=f64), 1e-10)
            r_tau = jnp.asarray(R, dtype=f64) + jnp.asarray(state.z_vec, dtype=f64) * jnp.asarray(state.nu_fit, dtype=f64)
            tau_prior_prec = 1.0 / jnp.maximum(jnp.square(jnp.asarray(state.tau_prior_sd, dtype=f64)), 1e-12)
            tau_prior_m64 = jnp.asarray(state.tau_prior_mean, dtype=f64)
            prec_tau = jnp.sum(w_trt64) / sig2_64 + tau_prior_prec
            rhs_tau = jnp.sum(w_trt64 * r_tau) / sig2_64 + tau_prior_prec * tau_prior_m64
            mean_tau = rhs_tau / prec_tau
            sd_tau = 1.0 / jnp.sqrt(prec_tau)
            unconstrained_tau = mean_tau + random.normal(k_nu, (), dtype=jnp.float32).astype(f64) * sd_tau
            trunc_tau = _sample_truncated_normal_64(k_nu, mean_tau, sd_tau, state.tau_min, state.tau_max)
            has_tau_clip = jnp.isfinite(jnp.asarray(state.tau_min, dtype=f64)) | jnp.isfinite(jnp.asarray(state.tau_max, dtype=f64))
            tau_scalar = jnp.where(has_tau_clip, trunc_tau, unconstrained_tau)
            nu_fit = jnp.full(state.N_units, tau_scalar.astype(jnp.float32), dtype=jnp.float32)
            tau_intercept_new = jnp.float32(0.0)
            R = (r_tau - jnp.asarray(state.z_vec, dtype=f64) * tau_scalar).astype(jnp.float32)
            forest_nu = state.forest_nu
            s_eta2_new = state.s_eta2
            xi_nu_new = state.xi_nu

    # If baseline_cells == "all", refresh sigma^2 from the full residual after tau_k update
    if not state.baseline_untreated_only:
        with enable_x64():
            f64 = jnp.float64
            ssr_all = jnp.sum(jnp.asarray(state.w_unit, dtype=f64) * jnp.square(jnp.asarray(R, dtype=f64)))
            sig2_all = _sample_inv_gamma_64(
                random.fold_in(key, 9901),
                jnp.asarray(state.sigma_prior_a, dtype=f64) + 0.5 * state.N_all_dof,
                jnp.asarray(state.sigma_prior_b, dtype=f64) + 0.5 * ssr_all,
            )
            sigma2_new = jnp.clip(sig2_all, 1e-8, 1e6).astype(jnp.float32)

    error_cov_inv = eqx.tree_at(
        lambda w_: w_.value,
        view_mu.error_cov_inv if state.use_forests else state.error_cov_inv,
        jnp.reciprocal(sigma2_new).astype(jnp.float32),
    )
    if state.use_forests:
        return eqx.tree_at(
            lambda s: (
                s.resid,
                s.forest,
                s.resid_eff_scale,
                s.resid_inexact_integral,
                s.error_cov_inv,
                s.prec_scale,
                s.inv_sdev_scale,
                s.inv_sdev_unit,
                s.forest_nu,
                s.resid_nu,
                s.resid_eff_scale_nu,
                s.resid_inexact_integral_nu,
                s.prec_scale_nu,
                s.inv_sdev_scale_nu,
                s.inv_sdev_unit_nu,
                s.mu_fit,
                s.nu_fit,
                s.intercept,
                s.tau_intercept,
                s.theta_f,
                s.f_loadings,
                s.s_f2,
                s.zeta,
                s.s_zeta2,
                s.ell_idx,
                s.sigma2,
                s.s_eta2,
                s.xi_nu,
            ),
            state,
            (
                R,
                forest_mu,
                view_mu.resid_eff_scale,
                view_mu.resid_inexact_integral,
                error_cov_inv,
                view_mu.prec_scale,
                view_mu.inv_sdev_scale,
                view_mu.inv_sdev_unit,
                forest_nu,
                view_nu.resid,
                view_nu.resid_eff_scale,
                view_nu.resid_inexact_integral,
                view_nu.prec_scale,
                view_nu.inv_sdev_scale,
                view_nu.inv_sdev_unit,
                mu_fit,
                nu_fit,
                intercept_new,
                tau_intercept_new,
                theta_new,
                f_new,
                s_f2_new,
                zeta_new,
                s_zeta2_new,
                ell_idx_new,
                sigma2_new,
                s_eta2_new,
                xi_nu_new,
            ),
            is_leaf=lambda x: x is None,
        )
    return eqx.tree_at(
        lambda s: (
            s.resid,
            s.forest,
            s.forest_nu,
            s.error_cov_inv,
            s.mu_fit,
            s.nu_fit,
            s.intercept,
            s.tau_intercept,
            s.theta_f,
            s.f_loadings,
            s.s_f2,
            s.zeta,
            s.s_zeta2,
            s.ell_idx,
            s.sigma2,
        ),
        state,
        (
            R,
            forest_mu,
            forest_nu,
            error_cov_inv,
            mu_fit,
            nu_fit,
            intercept_new,
            tau_intercept_new,
            theta_new,
            f_new,
            s_f2_new,
            zeta_new,
            s_zeta2_new,
            ell_idx_new,
            sigma2_new,
        ),
    )


def closed_form_linear_gaussian_posterior(
    *,
    dy_k: np.ndarray,
    df_k: np.ndarray,
    z_vec: np.ndarray,
    w_unit: np.ndarray,
    lambda_mat: np.ndarray,
    ctrl_idx: np.ndarray,
    trt_idx: np.ndarray,
    sigma2: float,
    s_f2: float,
    s_zeta2: float,
    C_corr: np.ndarray | None = None,
    intercept_prior_sd: float = 10.0,
    theta_prior_mean: float = 1.0,
    theta_prior_sd: float = 0.5,
    tau_prior_mean: float = 0.0,
    tau_prior_sd: float = 0.5,
) -> dict[str, Any]:
    """Exact closed-form Gaussian posterior for ``(intercept, theta_k, f_k, zeta_k, tau_k)``
    in the linear-Gaussian limit (``use_forests=False``, fixed variance components and
    untruncated Normal prior on ``theta_k``) under the control-only cut (``baseline_cells='untreated'``).

    Used by the Section 5.1 unit test to verify that ``rank_k_mode_single_step``
    with ``use_forests=False`` recovers the analytical mean and covariance.
    """
    dy = np.asarray(dy_k, dtype=np.float64)
    df = np.asarray(df_k, dtype=np.float64)
    w = np.asarray(w_unit, dtype=np.float64)
    lam = np.asarray(lambda_mat, dtype=np.float64)
    c_idx = np.asarray(ctrl_idx, dtype=np.int32)
    t_idx = np.asarray(trt_idx, dtype=np.int32)
    n0 = int(c_idx.size)
    n1 = int(t_idx.size)
    N = int(dy.size)
    q = int(lam.shape[1])

    # Control-group design for beta_0 = (intercept, theta_k, f_k)
    X0 = np.column_stack([np.ones(n0, dtype=np.float64), df[c_idx], lam[c_idx]])
    prior_mean_0 = np.zeros(2 + q, dtype=np.float64)
    prior_mean_0[1] = float(theta_prior_mean)
    prior_var_0 = np.array(
        [float(intercept_prior_sd) ** 2, float(theta_prior_sd) ** 2] + [float(s_f2)] * q,
        dtype=np.float64,
    )
    Sigma_beta0 = np.diag(prior_var_0)

    if C_corr is not None and s_zeta2 > 0.0:
        C00 = np.asarray(C_corr[np.ix_(c_idx, c_idx)], dtype=np.float64)
        C10 = np.asarray(C_corr[np.ix_(t_idx, c_idx)], dtype=np.float64)
        C11 = np.asarray(C_corr[np.ix_(t_idx, t_idx)], dtype=np.float64)
    else:
        C00 = np.zeros((n0, n0), dtype=np.float64)
        C10 = np.zeros((n1, n0), dtype=np.float64)
        C11 = np.zeros((n1, n1), dtype=np.float64)

    # Joint Gaussian distribution of (beta_0, zeta_0, zeta_1, dy_0) on I_0:
    # Marginal covariance of dy_0 = X0 beta_0 + zeta_0 + eps_0
    V_y0 = (
        X0 @ Sigma_beta0 @ X0.T
        + float(s_zeta2) * C00
        + np.diag(float(sigma2) / w[c_idx])
    )
    m_y0 = X0 @ prior_mean_0
    dy0_centered = dy[c_idx] - m_y0
    V_y0_inv_dy0 = np.linalg.solve(V_y0, dy0_centered)

    # Posterior of beta_0 = (intercept, theta_k, f_k) | dy_0
    Cov_beta0_y0 = Sigma_beta0 @ X0.T
    post_mean_beta0 = prior_mean_0 + Cov_beta0_y0 @ V_y0_inv_dy0
    post_cov_beta0 = Sigma_beta0 - Cov_beta0_y0 @ np.linalg.solve(V_y0, Cov_beta0_y0.T)

    # Posterior of counterfactual baseline on treated units:
    # eta_1 = X1 beta_0 + zeta_1, where Cov(eta_1, dy_0) = X1 Sigma_beta0 X0^T + s_zeta^2 C10
    X1 = np.column_stack([np.ones(n1, dtype=np.float64), df[t_idx], lam[t_idx]])
    Cov_eta1_y0 = X1 @ Sigma_beta0 @ X0.T + float(s_zeta2) * C10
    Var_eta1_prior = X1 @ Sigma_beta0 @ X1.T + float(s_zeta2) * C11
    post_mean_eta1 = X1 @ prior_mean_0 + Cov_eta1_y0 @ V_y0_inv_dy0
    post_cov_eta1 = Var_eta1_prior - Cov_eta1_y0 @ np.linalg.solve(V_y0, Cov_eta1_y0.T)

    # Treatment scalar tau_k | eta_1, dy_1:
    # Conditional on eta_1: prec_tau = sum_{i in I_1} w_i / sigma^2 + 1 / tau_prior_sd^2,
    # mean_tau(eta_1) = prec_tau^{-1} (sum_{i in I_1} (w_i / sigma^2) (dy_{1,i} - eta_{1,i}) + tau_prior_mean / tau_prior_sd^2).
    w1 = w[t_idx]
    tau_prior_prec = 1.0 / (float(tau_prior_sd) ** 2)
    prec_tau = float(np.sum(w1) / float(sigma2) + tau_prior_prec)
    alpha_1 = (w1 / float(sigma2)) / prec_tau  # (n1,)
    tau_mean = float(alpha_1 @ (dy[t_idx] - post_mean_eta1) + (tau_prior_prec / prec_tau) * float(tau_prior_mean))
    tau_var_within = 1.0 / prec_tau
    tau_var_between = float(alpha_1 @ post_cov_eta1 @ alpha_1)
    tau_var_total = tau_var_within + tau_var_between

    return {
        "beta0_mean": post_mean_beta0,
        "beta0_cov": post_cov_beta0,
        "eta1_mean": post_mean_eta1,
        "eta1_cov": post_cov_eta1,
        "tau_mean": tau_mean,
        "tau_sd": math.sqrt(tau_var_total),
        "tau_var_within": tau_var_within,
        "tau_var_between": tau_var_between,
        "N": N,
    }


def _bin_mode_predictors(
    x_mod: np.ndarray,
    ps_vec: np.ndarray,
    df_std: np.ndarray,
    has_forecast: bool,
    key: jax.Array,
    max_bins: int = 255,
) -> tuple[jax.Array, jax.Array, jax.Array]:
    """Bin cross-sectional predictors into a single ``(p, N)`` uint8 matrix for ``bartz``.

    Layout of rows:
      * ``0 .. P_x - 1``: ``x_mod`` (visible to both ``mu_k`` and ``tau_k``)
      * ``P_x``: ``df_std`` (when ``has_forecast``; visible to ``mu_k`` only)
      * ``P_x + 1`` (or ``P_x``): ``ps_vec`` (visible to ``mu_k`` only)
    """
    cols = [np.asarray(x_mod, dtype=np.float32)]
    if has_forecast:
        cols.append(np.asarray(df_std, dtype=np.float32)[:, None])
    cols.append(np.asarray(ps_vec, dtype=np.float32)[:, None])
    raw = np.column_stack(cols).T  # (p, N)

    binner = UniqueQuantileBinner(jnp.asarray(raw, dtype=jnp.float32), max_bins=max_bins, key=key)
    X_binned = binner.bin(jnp.asarray(raw, dtype=jnp.float32))
    max_split_all = np.asarray(binner.max_split, dtype=np.uint8)

    p_x = int(x_mod.shape[1])
    max_split_mu = max_split_all.copy()
    max_split_nu = np.zeros_like(max_split_all)
    max_split_nu[:p_x] = max_split_all[:p_x]
    return (
        X_binned,
        jnp.asarray(max_split_mu, dtype=np.uint8),
        jnp.asarray(max_split_nu, dtype=np.uint8),
    )


def _scan_sweeps(
    key: jax.Array,
    st0: RankKModeState,
    spatial_cache: SpatialKernelCache,
    num_burnin: int,
    num_draws: int,
    thin: int,
) -> tuple[RankKModeState, dict[str, jax.Array]]:
    """Run ``num_burnin + num_draws * thin`` sweeps from ``st0`` in one ``lax.scan``; keep every ``thin``-th post-burn-in sweep."""

    def scan_body(carry: tuple[jax.Array, RankKModeState], _: None) -> tuple[tuple[jax.Array, RankKModeState], dict[str, jax.Array]]:
        k_c, st = carry
        k_c, k_s = random.split(k_c)
        st = rank_k_mode_single_step(k_s, st, spatial_cache)
        out = {
            "tau": st.tau_intercept + st.nu_fit,
            "mu": st.mu_fit + st.intercept,
            "theta_f": st.theta_f,
            "sigma": jnp.sqrt(st.sigma2),
            "s_zeta": jnp.sqrt(st.s_zeta2),
            "s_f": jnp.sqrt(st.s_f2),
            "ell_idx": st.ell_idx,
            "tau_intercept": st.tau_intercept,
            "s_eta2": st.s_eta2,
            "sigma_nu_leaf2": jnp.reciprocal(jnp.maximum(st.forest_nu.leaf_prior_cov_inv, 1e-30)),
        }
        return (k_c, st), out

    thin = max(int(thin), 1)
    total_steps = int(num_burnin) + int(num_draws) * thin
    (_, final_state), full_trace = lax.scan(scan_body, (key, st0), None, length=total_steps)
    b = int(num_burnin)
    trace = {k: v[b:][::thin] for k, v in full_trace.items()}
    return final_state, trace


def _run_chain_loop(
    key: jax.Array,
    df_vec: jax.Array,
    init_state: RankKModeState,
    spatial_cache: SpatialKernelCache,
    num_burnin: int,
    num_draws: int,
    thin: int = 1,
) -> tuple[RankKModeState, dict[str, jax.Array]]:
    """Run ``num_burnin + num_draws * thin`` sweeps in one ``lax.scan``; keep every ``thin``-th post-burn-in sweep."""
    R_init = init_state.resid + init_state.theta_f * (init_state.df_vec - df_vec)
    st0 = eqx.tree_at(lambda s: (s.df_vec, s.resid), init_state, (df_vec, R_init))
    return _scan_sweeps(key, st0, spatial_cache, num_burnin, num_draws, thin)


def _run_chain_loop_shifted(
    key: jax.Array,
    df_vec: jax.Array,
    dy_shift: jax.Array,
    init_state: RankKModeState,
    spatial_cache: SpatialKernelCache,
    num_burnin: int,
    num_draws: int,
    thin: int = 1,
) -> tuple[RankKModeState, dict[str, jax.Array]]:
    """:func:`_run_chain_loop` for a track whose mode outcome is ``dy_k + dy_shift`` (scaled units).

    The sampler only ever reads the outcome through the residual, so a per-track
    outcome shift is a per-track shift of the initial residual. Used to propagate
    higher-mode posterior draws into the Mode-1 target of masked panels, exactly as
    forecast draws are propagated through ``df_vec``.
    """
    R_init = init_state.resid + dy_shift + init_state.theta_f * (init_state.df_vec - df_vec)
    st0 = eqx.tree_at(lambda s: (s.df_vec, s.resid), init_state, (df_vec, R_init))
    return _scan_sweeps(key, st0, spatial_cache, num_burnin, num_draws, thin)


_run_chains_vmap = enable_x64()(
    jax.jit(
        enable_x64()(jax.vmap(_run_chain_loop, in_axes=(0, 0, None, None, None, None, None))),
        static_argnums=(4, 5, 6),
    )
)

_run_chains_vmap_shifted = enable_x64()(
    jax.jit(
        enable_x64()(jax.vmap(_run_chain_loop_shifted, in_axes=(0, 0, 0, None, None, None, None, None))),
        static_argnums=(5, 6, 7),
    )
)


def sample_rank_k_mode_jax(
    *,
    dy_k: np.ndarray,
    df_k_point: np.ndarray,
    df_k_draws: np.ndarray | None,
    x_mod: np.ndarray,
    ps_vec: np.ndarray,
    z_vec: np.ndarray,
    w_unit: np.ndarray,
    lambda_mat: np.ndarray,
    ctrl_idx: np.ndarray,
    trt_idx: np.ndarray,
    spatial_cache: SpatialKernelCache,
    num_chains: int = 4,
    num_burnin: int = 200,
    num_draws_per_chain: int = 150,
    warm_burnin: int = 50,
    num_trees_mu: int = 75,
    num_trees_nu: int = 35,
    has_forecast: bool = True,
    use_forests: bool = True,
    use_factor_loadings: bool = True,
    use_spatial_gp: bool = True,
    baseline_untreated_only: bool = True,
    theta_f_clip: tuple[float, float] = (0.0, 1.5),
    theta_prior_mean: float = 1.0,
    theta_prior_sd: float = 0.5,
    tau_prior_mean: float = 0.0,
    tau_prior_sd: float | None = None,
    tau_clip: tuple[float, float] = (-float("inf"), float("inf")),
    random_seed: int = 2026,
    att_unit_weights: np.ndarray | None = None,
    tau_prior_family: str = "normal",
    tau_prior_params: tuple[float, float] | None = None,
    sigma_nu_total: float = 0.25,
    alpha_nu: float = 0.25,
    beta_nu: float = 3.0,
    nu_working_zeta_weight: float = 1.0,
    learn_sigma_nu: bool = False,
    sigma_nu_prior_scale: float = 0.5,
    treated_extra_variance: bool = False,
    thin: int = 1,
    dy_k_draws: np.ndarray | None = None,
) -> dict[str, Any]:
    """Sample one mode ``k`` across either a single point forecast or ``M`` posterior forecast draws.

    When ``df_k_draws`` is ``None`` (or ``M == 1``), runs ``num_chains`` parallel
    chains on ``df_k_point`` via ``jax.vmap``. When ``df_k_draws`` has shape
    ``(M, N)`` with ``M > 1``, runs ``C * M`` parallel chain-by-forecast tracks
    via the same depth-1 ``jax.vmap`` with unbatched ``init_state``.

    Burn-in: every chain (or chain-by-forecast track) starts from the same
    ``init_state`` and runs its own burn-in, ``num_burnin`` sweeps in the
    single-forecast path and ``max(num_burnin, warm_burnin)`` sweeps in the
    multi-draw path; ``warm_burnin`` is therefore not a warm start and is ignored
    in the single-forecast path. :func:`longbet.fit_rank_k_longbet` passes
    ``warm_burnin=0`` unless the deprecated ``RankKLongBetConfig.warm_burnin`` is
    set, and passes its own tree counts; the keyword defaults here are not the
    package defaults.

    The trailing keyword block (``att_unit_weights`` ... ``thin``) holds the
    remediation extensions; at their defaults the sampler is the paper's,
    bit-for-bit. ``att_unit_weights`` (ATT weights ``a_i = L_{G_i} / |O_trt|``)
    are only used when a prior or truncation acts on the mode-1 intercept, in
    which case the treatment forest is re-centered with them so that
    ``tau_intercept`` is the ATT. ``thin`` keeps every ``thin``-th post-burn-in
    sweep (``num_draws_per_chain`` retained draws per track either way).

    ``dy_k_draws`` (``(M, N)``, optional) gives every track its own mode outcome
    in place of ``dy_k`` -- the Mode-1 targets of a masked panel under ``M``
    higher-mode posterior draws. Tracks are paired with the forecast draws when
    both are given (equal ``M``); the outcome enters as a per-track shift of the
    initial residual, so a single track or ``None`` reproduces the plain path.
    Units with ``w_unit == 0`` take no part in any likelihood (their forest
    predictions are still returned) and are excluded from the variance updates'
    sample counts and from the outcome scale.
    """
    base_key = random.key(int(random_seed))
    k_bin, k_main = random.split(base_key)

    # Standardize point forecast for tree splitting
    if has_forecast:
        df_std_point = (df_k_point - np.mean(df_k_point)) / max(float(np.std(df_k_point)), 1e-6)
    else:
        df_std_point = np.zeros_like(df_k_point)

    X_binned_pt, max_split_mu, max_split_nu = _bin_mode_predictors(
        x_mod=x_mod,
        ps_vec=ps_vec,
        df_std=df_std_point,
        has_forecast=has_forecast,
        key=k_bin,
    )

    # Standardize dy_k by the control-group post-forecast residual SD so leaf priors
    # operate on exact post-CUPED residual units (strictly pre/control-only). Units
    # with weight 0 (unsupported modes of masked panels) carry no information.
    c_idx_np = np.asarray(ctrl_idx, dtype=np.int32)
    w_pos_c = np.asarray(w_unit[c_idx_np], dtype=np.float64) > 0.0
    if has_forecast and c_idx_np.size > 0:
        w_c = np.asarray(w_unit[c_idx_np], dtype=np.float64)
        df_c = np.asarray(df_k_point[c_idx_np], dtype=np.float64)
        dy_c = np.asarray(dy_k[c_idx_np], dtype=np.float64)
        df_c_cent = df_c - np.average(df_c, weights=w_c)
        dy_c_cent = dy_c - np.average(dy_c, weights=w_c)
        denom_c = float(np.sum(w_c * np.square(df_c_cent))) + 1e-4
        num_c = float(np.sum(w_c * df_c_cent * dy_c_cent))
        th_c = float(np.clip(num_c / denom_c, float(theta_f_clip[0]), float(theta_f_clip[1])))
        resid_ctrl = dy_c - th_c * df_c
    else:
        resid_ctrl = np.asarray(dy_k[c_idx_np], dtype=np.float64)
    ctrl_sd = float(np.std(resid_ctrl[w_pos_c]))
    y_scale = max(ctrl_sd, 1e-4)
    dy_scaled = dy_k / y_scale
    df_pt_scaled = df_k_point / y_scale

    tau_prior_mean_scaled = float(tau_prior_mean) / y_scale
    tau_prior_sd_scaled = (float(tau_prior_sd) / y_scale) if tau_prior_sd is not None else None
    tau_clip_scaled = (float(tau_clip[0]) / y_scale, float(tau_clip[1]) / y_scale)

    init_state = init_rank_k_mode(
        X_binned=X_binned_pt,
        max_split_mu=max_split_mu,
        max_split_nu=max_split_nu,
        dy_k=jnp.asarray(dy_scaled, dtype=jnp.float32),
        df_k=jnp.asarray(df_pt_scaled, dtype=jnp.float32),
        z_vec=jnp.asarray(z_vec, dtype=jnp.float32),
        w_unit=jnp.asarray(w_unit, dtype=jnp.float32),
        lambda_mat=jnp.asarray(lambda_mat, dtype=jnp.float32),
        ctrl_idx=ctrl_idx,
        trt_idx=trt_idx,
        num_trees_mu=num_trees_mu,
        num_trees_nu=num_trees_nu,
        has_forecast=has_forecast,
        use_forests=use_forests,
        use_factor_loadings=use_factor_loadings,
        use_spatial_gp=use_spatial_gp,
        baseline_untreated_only=baseline_untreated_only,
        theta_f_clip=theta_f_clip,
        theta_prior_mean=theta_prior_mean,
        theta_prior_sd=theta_prior_sd,
        tau_prior_mean=tau_prior_mean_scaled,
        tau_prior_sd=tau_prior_sd_scaled,
        tau_clip=tau_clip_scaled,
        att_unit_weights=att_unit_weights,
        tau_prior_family=tau_prior_family,
        tau_prior_params=tau_prior_params,
        sigma_nu_total=float(sigma_nu_total),
        alpha_nu=float(alpha_nu),
        beta_nu=float(beta_nu),
        nu_working_zeta_weight=float(nu_working_zeta_weight),
        learn_sigma_nu=bool(learn_sigma_nu),
        sigma_nu_prior_scale=float(sigma_nu_prior_scale),
        treated_extra_variance=bool(treated_extra_variance),
    )
    thin_eff = max(int(thin), 1)

    C = int(num_chains)
    m_df = int(df_k_draws.shape[0]) if (df_k_draws is not None and has_forecast) else 1
    m_dy = int(dy_k_draws.shape[0]) if dy_k_draws is not None else 1
    if m_df > 1 and m_dy > 1 and m_df != m_dy:
        raise ValueError(f"df_k_draws ({m_df}) and dy_k_draws ({m_dy}) must hold the same number of draws")
    if m_df <= 1 and m_dy <= 1:
        chain_keys = random.split(k_main, C)
        if df_k_draws is not None and has_forecast and df_k_draws.shape[0] == 1:
            df_single = jnp.asarray(np.asarray(df_k_draws[0], dtype=np.float64) / y_scale, dtype=jnp.float32)
        else:
            df_single = init_state.df_vec
        df_batch = jnp.repeat(df_single[None, :], C, axis=0)
        if dy_k_draws is not None and dy_k_draws.shape[0] == 1:
            dy_shift = jnp.asarray((np.asarray(dy_k_draws[0], dtype=np.float64) - np.asarray(dy_k, dtype=np.float64)) / y_scale, dtype=jnp.float32)
            _, trace = _run_chains_vmap_shifted(
                chain_keys,
                df_batch,
                jnp.repeat(dy_shift[None, :], C, axis=0),
                init_state,
                spatial_cache,
                int(num_burnin),
                int(num_draws_per_chain),
                thin_eff,
            )
        else:
            _, trace = _run_chains_vmap(
                chain_keys,
                df_batch,
                init_state,
                spatial_cache,
                int(num_burnin),
                int(num_draws_per_chain),
                thin_eff,
            )
        # trace["tau"] has shape (C, S, N) -> rescale by y_scale
        tau_chains = np.asarray(trace["tau"], dtype=np.float64) * y_scale  # (C, S, N)
        theta_chains = np.asarray(trace["theta_f"], dtype=np.float64)      # (C, S)
        sigma_chains = np.asarray(trace["sigma"], dtype=np.float64) * y_scale
        s_zeta_chains = np.asarray(trace["s_zeta"], dtype=np.float64) * y_scale
        ell_idx_chains = np.asarray(trace["ell_idx"], dtype=np.int32)
        ell_grid_np = np.asarray(spatial_cache.ell_grid, dtype=np.float64)
        ell_chains = ell_grid_np[np.clip(ell_idx_chains, 0, len(ell_grid_np) - 1)]
        tau_int_chains = np.asarray(trace["tau_intercept"], dtype=np.float64) * y_scale
        s_eta2_chains = np.asarray(trace["s_eta2"], dtype=np.float64) * (y_scale**2)
        sigma_nu_chains = np.sqrt(np.asarray(trace["sigma_nu_leaf2"], dtype=np.float64)) * y_scale * math.sqrt(
            max(int(num_trees_nu), 1)
        )

        C, S, N = tau_chains.shape
        tau_flat = tau_chains.reshape(C * S, N).T  # (N, M_total)
        return {
            "tau_draws": tau_flat,
            "tau_chains": tau_chains,
            "comp_tau_draws": tau_chains.reshape(1, C * S, N),
            "theta_f_draws": theta_chains.reshape(-1),
            "sigma_draws": sigma_chains.reshape(-1),
            "s_zeta_draws": s_zeta_chains.reshape(-1),
            "ell_draws": ell_chains.reshape(-1),
            "tau_intercept_draws": tau_int_chains.reshape(-1),
            "s_eta2_draws": s_eta2_chains.reshape(-1),
            "sigma_nu_total_draws": sigma_nu_chains.reshape(-1),
            "y_scale": float(y_scale),
        }

    # Multi-component propagation over M forecast draws (and / or M Mode-1 targets) via flat depth-1 jax.vmap
    M = max(m_df, m_dy)
    if m_df > 1:
        df_draws_scaled = jnp.asarray(np.asarray(df_k_draws, dtype=np.float64) / y_scale, dtype=jnp.float32)
    else:
        df_draws_scaled = jnp.repeat(init_state.df_vec[None, :], M, axis=0)
    comp_keys_flat = random.split(k_main, C * M)
    df_batch = jnp.tile(df_draws_scaled, (C, 1))
    total_burn = max(int(num_burnin), int(warm_burnin))

    if m_dy > 1:
        dy_shift_scaled = jnp.asarray(
            (np.asarray(dy_k_draws, dtype=np.float64) - np.asarray(dy_k, dtype=np.float64)[None, :]) / y_scale,
            dtype=jnp.float32,
        )
        _, trace_flat = _run_chains_vmap_shifted(
            comp_keys_flat,
            df_batch,
            jnp.tile(dy_shift_scaled, (C, 1)),
            init_state,
            spatial_cache,
            int(total_burn),
            int(num_draws_per_chain),
            thin_eff,
        )
    else:
        _, trace_flat = _run_chains_vmap(
            comp_keys_flat,
            df_batch,
            init_state,
            spatial_cache,
            int(total_burn),
            int(num_draws_per_chain),
            thin_eff,
        )
    S = int(num_draws_per_chain)
    N = int(dy_k.shape[0])
    tau_cm = np.asarray(trace_flat["tau"], dtype=np.float64).reshape(C, M, S, N) * y_scale
    # Component-major array (M, C * S, N) for exact within/between variance decomposition
    comp_tau_arr = np.transpose(tau_cm, (1, 0, 2, 3)).reshape(M, C * S, N)
    # Chain-major array (C, M * S, N) ordered so split-half R-hat compares
    # the first half of every component against the second half of every component
    half = S // 2
    if half >= 1:
        tau_chains = np.concatenate(
            [
                tau_cm[:, :, :half, :].reshape(C, M * half, N),
                tau_cm[:, :, half:, :].reshape(C, M * (S - half), N),
            ],
            axis=1,
        )
    else:
        tau_chains = tau_cm.reshape(C, M * S, N)
    tau_flat = tau_chains.reshape(C * M * S, N).T  # (N, C * M * S), chain-major

    theta_cm = np.transpose(np.asarray(trace_flat["theta_f"], dtype=np.float64).reshape(C, M, S), (1, 0, 2)).reshape(-1)
    sigma_cm = (np.transpose(np.asarray(trace_flat["sigma"], dtype=np.float64).reshape(C, M, S), (1, 0, 2)) * y_scale).reshape(-1)
    szeta_cm = (np.transpose(np.asarray(trace_flat["s_zeta"], dtype=np.float64).reshape(C, M, S), (1, 0, 2)) * y_scale).reshape(-1)
    e_idx = np.transpose(np.asarray(trace_flat["ell_idx"], dtype=np.int32).reshape(C, M, S), (1, 0, 2)).reshape(-1)
    ell_grid_np = np.asarray(spatial_cache.ell_grid, dtype=np.float64)
    ell_cm = ell_grid_np[np.clip(e_idx, 0, len(ell_grid_np) - 1)]

    def _chain_major_scalar(name: str, scale: float) -> np.ndarray:
        """Scalar trace in the draw order of ``tau_draws`` (chain-major with the split-half reordering)."""
        arr = np.asarray(trace_flat[name], dtype=np.float64).reshape(C, M, S) * scale
        if half >= 1:
            arr = np.concatenate([arr[:, :, :half].reshape(C, M * half), arr[:, :, half:].reshape(C, M * (S - half))], axis=1)
        else:
            arr = arr.reshape(C, M * S)
        return arr.reshape(-1)

    return {
        "tau_draws": tau_flat,
        "tau_chains": tau_chains,
        "comp_tau_draws": comp_tau_arr,
        "theta_f_draws": theta_cm,
        "sigma_draws": sigma_cm,
        "s_zeta_draws": szeta_cm,
        "ell_draws": ell_cm,
        # Aligned with the columns of ``tau_draws`` (unlike the component-major traces above)
        "tau_intercept_draws": _chain_major_scalar("tau_intercept", y_scale),
        "s_eta2_draws": _chain_major_scalar("s_eta2", y_scale**2),
        "sigma_nu_total_draws": np.sqrt(_chain_major_scalar("sigma_nu_leaf2", 1.0)) * y_scale * math.sqrt(max(int(num_trees_nu), 1)),
        "y_scale": float(y_scale),
    }
