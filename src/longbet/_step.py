# Copyright 2026 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""One full Gibbs sweep for LongBet.

Two invariants govern this module.

**Units.** ``bartz`` stores its residual scaled: ``resid_unit * resid == y -
sum(trees)``, and ``resid_unit`` differs between the two forests because their
leaf priors and tree counts differ.  LongBet therefore carries the full-model
residual ``R`` in **data units** at all times and converts only at the boundary
of a ``bartz`` call, via :func:`_load` and :func:`_read`. Forest fits are summed
from cached leaf memberships, without traversing trees again. Differencing
residuals after division by a small treatment multiplier loses leaf updates
to cancellation; those errors would then enter the GP and coding conditionals.

**Chains.** :func:`longbet_single_step` always operates on a single chain.
Multiple chains are run by :func:`longbet_step`, which splits the state into its
per-chain and shared parts and vmaps over the former only, so the design matrix
is never replicated (see ``longbet._state``).
"""

from __future__ import annotations

from functools import partial
from typing import Any

import equinox as eqx
import jax
import jax.numpy as jnp
from jax import lax, random
import jax.scipy.linalg  # noqa: F401  (used by the linear-moderator step)
from jaxtyping import Array, Float32, Key

from bartz.mcmcstep import State, Wishart
from bartz.mcmcstep import step as bartz_step
from bartz.mcmcstep._step import step_z

import longbet._ar1 as ar1_mod
from longbet._change_internal_move import change_internal_step
from longbet._change_move import change_step
from longbet._regrow_move import regrow_step
from longbet._gp import beta_prior_quadform, sample_beta_gp
from longbet._ordinal import sample_cutpoints_marginalized, sample_ordinal_latents
from longbet._forest_cache import refresh_prec_tree, current_forest_fit
from longbet._x64 import enable_x64
from longbet._inter_ensemble_move import (
    inter_ensemble_beta_gamma_step,
    inter_ensemble_transfer_step,
)
from longbet._ridge import ridge_scale_step, scale_group_step
from longbet._scales import compute_treatment_scale_attrs
from longbet._trend_basis import (
    horseshoe_col_scales,
    sample_horseshoe,
    trend_block_step,
)
from longbet._state import LongBetState, chain_filter_spec

#: Weights at or below this magnitude are treated as negligible in the ordinary
#: treatment update. The stable collapsed fallback retains nonzero weights.
WEIGHT_EPS = 1e-6


def _load(
    r: Float32[Array, ' n'], resid_unit: Float32[Array, '']
) -> Float32[Array, ' n']:
    """Convert a residual from data units to bartz's state units."""
    return (r / resid_unit).astype(jnp.float32)


def _read(
    resid: Float32[Array, ' n'], resid_unit: Float32[Array, '']
) -> Float32[Array, ' n']:
    """Convert a residual from bartz's state units to data units."""
    return resid * resid_unit


def _sample_inv_gamma(
    key: Key[Array, ''],
    shape: Float32[Array, ''],
    rate: Float32[Array, ''],
) -> Float32[Array, '']:
    """Draw from Inv-Gamma(shape, rate) via a Gamma draw."""
    g = jnp.maximum(random.gamma(key, shape, dtype=jnp.float32), 1e-30)
    return (rate / g).astype(jnp.float32)


def _tempered_single_step(key: Key[Array, ''], state: LongBetState) -> LongBetState:
    """One sweep against ``p(theta) L(theta)^beta`` with ``beta = state.temperature``.

    Every mean, latent and coefficient update runs with the per-cell conditional
    precision ``beta / sigma2``; a continuous outcome's ``sigma2`` is then drawn
    from its tempered conditional ``IG(a + beta m / 2, b + beta SSR / 2)``. See
    ``longbet._tempering``.
    """
    if state.baseline_untreated_only:
        raise NotImplementedError(
            "baseline_cells='untreated' cannot be combined with parallel tempering"
        )
    beta = jnp.asarray(state.temperature, dtype=jnp.float32)
    prec = jnp.where(state.obs_mask, beta / state.sigma2, 0.0).astype(jnp.float32)
    k_step, k_sig = random.split(key)
    new = longbet_single_step(k_step, state, prec)
    if state.num_categories < 2:
        m_obs = jnp.sum(state.obs_mask.astype(jnp.float32))
        resid_innov = (ar1_mod.transform(jnp.where(new.obs_mask, new.resid, 0.0),
                                         new.ar1_phi, new.ar1_first)
                       if state.ar1 else new.resid)
        ssr = jnp.sum(jnp.square(jnp.where(new.obs_mask, resid_innov, 0.0)))
        sigma2_new = _sample_inv_gamma(
            k_sig,
            jnp.float32(state.sigma_prior_a) + 0.5 * beta * m_obs,
            jnp.float32(state.sigma_prior_b) + 0.5 * beta * ssr,
        )
        error_cov_inv = eqx.tree_at(
            lambda w: w.value, new.error_cov_inv, jnp.reciprocal(sigma2_new).astype(jnp.float32)
        )
        new = eqx.tree_at(lambda s: (s.sigma2, s.error_cov_inv), new, (sigma2_new, error_cov_inv))
    return new


#: Proposal standard deviation of ``log c`` in the beta-nu scale ridge move.
RIDGE_PROPOSAL_SIGMA = 0.2


def _stabilize_collapsed_nu_stats(
    kappa, target, second_moment, denominator, u, w, obs_mask, unit_idx,
    num_units, prior_precision, cell_precision=None,
):
    """Repair cancellation-prone collapsed treatment statistics.

    The ordinary float32 path is retained where its denominators and precision
    factor are well conditioned. In the formerly floored regime, use weighted
    centered sums in float64. For P=sum(p), c=prior_precision and weighted means
    wbar, ubar, the exact marginal precision and score are

        Q = sum(p*(w-wbar)**2) + P*wbar**2*c/(P+c),
        h = sum(p*(w-wbar)*(u-ubar)) + P*wbar*ubar*c/(P+c).

    Here p=1 and c=sigma2/sigma_gamma2 for the scalar path; otherwise p is
    conditional observation precision and c=1/sigma_gamma2. Missing cells and
    zero-weight units contribute zero. No positive lower bound is imposed on Q.
    The returned mask identifies units needing this numerically stable path.
    """
    # Never-treated or fully missing units already contribute exactly zero;
    # they must not force every ordinary sweep onto the float64 path. Detect
    # signal before squaring so nonzero tiny weights still receive the repair.
    has_signal = jnp.zeros(num_units, dtype=jnp.bool_).at[unit_idx].max(
        obs_mask & (w != 0.0)
    )
    repair = has_signal & (
        (kappa <= 1e-4) | (second_moment <= 1e-12) | (denominator <= 1e-6)
        | (kappa * second_moment <= 1e-12)
        | ~jnp.isfinite(kappa) | ~jnp.isfinite(target)
    )

    @enable_x64()
    def stable(_):
        p = (jnp.ones_like(w, dtype=jnp.float64) if cell_precision is None
             else jnp.asarray(cell_precision, dtype=jnp.float64))
        p = jnp.where(obs_mask, p, 0.0)
        weights = jnp.where(obs_mask, jnp.asarray(w, jnp.float64), 0.0)
        response = jnp.where(obs_mask, jnp.asarray(u, jnp.float64), 0.0)

        def unit_sum(values):
            return jnp.zeros(num_units, dtype=jnp.float64).at[unit_idx].add(values)

        P = unit_sum(p)
        P_safe = jnp.where(P > 0.0, P, 1.0)
        wbar = unit_sum(p * weights) / P_safe
        ubar = unit_sum(p * response) / P_safe
        dw = weights - wbar[unit_idx]
        du = response - ubar[unit_idx]
        c = jnp.asarray(prior_precision, jnp.float64)
        mean_precision = P * (c / (P + c))
        precision = unit_sum(p * jnp.square(dw)) + mean_precision * jnp.square(wbar)
        score = unit_sum(p * dw * du) + mean_precision * wbar * ubar
        second = unit_sum(p * jnp.square(weights))
        stable_kappa = precision / jnp.where(second > 0.0, second, 1.0)
        stable_target = score / jnp.where(precision > 0.0, precision, 1.0)
        return (
            jnp.where(repair, stable_kappa.astype(kappa.dtype), kappa),
            jnp.where(repair, stable_target.astype(target.dtype), target),
        )

    kappa, target = lax.cond(jnp.any(repair), stable, lambda _: (kappa, target), None)
    return kappa, target, repair


def _draw_trajectory(
    key: Key[Array, ''],
    R: Float32[Array, ' n'],
    beta: Float32[Array, ' S_max_plus_1'],
    nu_fit: Float32[Array, ' n'],
    forest_nu: Any,
    b_z: Float32[Array, ' n'],
    state: LongBetState,
    sigma2: Float32[Array, ''],
    conditional_precision: Float32[Array, ' n'] | None,
    K_chol: Any | None = None,
    lin_fit: Float32[Array, ' n'] | None = None,
    delta: Any | None = None,
) -> tuple[Any, Any, Any, Any, Any, Any]:
    """Draw the exposure trajectory ``beta_S`` from its GP conditional.

    Conditions on the current treatment forest, unit intercepts and trend; the
    multiplier ``d_it = b_{Z_it} (nu_it + lin_it)`` is zero on untreated cells,
    so only treated cells enter the sufficient statistics. The sign of the
    trajectory is then mapped onto the canonical half-space
    ``sum_{s>=1} beta_s >= 0`` (the ``beta -> -beta, nu -> -nu`` quotient), the
    forest leaves (and, when present, the linear moderator ``delta`` and its
    fit ``lin_fit``) flipped to match, and the residual re-formed from the
    rounded draw.

    ``K_chol`` defaults to ``state.K_chol``; callers pass the effective kernel
    factor when ``learn_beta_scale`` rescales it within the sweep. ``lin_fit``
    and ``delta`` are the encouragement-model linear moderator (``None`` when
    ``linear_trt`` is off) and are returned unchanged in that case.

    Returns ``(beta_new, nu_fit, forest_nu, R, delta, lin_fit)``.
    """
    obs_mask = state.obs_mask
    if K_chol is None:
        K_chol = state.K_chol
    use_lin = lin_fit is not None
    d_vec = b_z * (nu_fit + lin_fit if use_lin else nu_fit)  # treatment factor with beta removed
    r_partial = R + d_vec * beta[state.exposure_idx]

    d_obs = jnp.where(obs_mask, d_vec, 0.0)
    r_obs = jnp.where(obs_mask, r_partial, 0.0)

    A = (
        jnp.zeros(state.S_max + 1, dtype=jnp.float32)
        .at[state.exposure_idx]
        .add(jnp.square(d_obs) / sigma2 if conditional_precision is None else
             jnp.square(d_obs) * conditional_precision)
    )
    C = (
        jnp.zeros(state.S_max + 1, dtype=jnp.float32)
        .at[state.exposure_idx]
        .add(d_obs * r_obs / sigma2 if conditional_precision is None else
             d_obs * r_obs * conditional_precision)
    )

    if state.sample_beta:
        beta_new = sample_beta_gp(key, A, C, K_chol)
        if not state.anchor_moderator:
            sign_beta = jnp.where(jnp.sum(beta_new[1:]) < 0.0, -1.0, 1.0)
            beta_new = sign_beta * beta_new
            nu_fit = sign_beta * nu_fit
            forest_nu = eqx.tree_at(
                lambda f: f.leaf_tree, forest_nu, sign_beta * forest_nu.leaf_tree
            )
            if use_lin:
                delta = sign_beta * delta
                lin_fit = sign_beta * lin_fit
                d_vec = b_z * (nu_fit + lin_fit)
            else:
                d_vec = b_z * nu_fit
    else:
        beta_new = beta

    beta_expanded = beta_new[state.exposure_idx]
    R = jnp.where(obs_mask, r_partial - d_vec * beta_expanded, 0.0)
    return beta_new, nu_fit, forest_nu, R, delta, lin_fit


# ---------------------------------------------------------------------------
# Distributed-lag exposure kernel (``LongBetState.use_lag_kernel``)
# ---------------------------------------------------------------------------
# The treatment term of cell (i, t) is ``w_it nu(x_i)`` with
# ``w_it = sum_l D_{i,t-l} kappa_l = (lag_design @ kappa)_it`` in place of
# ``b_z beta_{S_it}``. Everything that consumes ``w`` -- both forests, the
# linear moderator, the collapsed unit-intercept statistics -- is untouched;
# only the draw of the exposure coefficients changes, from a per-exposure
# diagonal to a dense (L + 1)-dimensional Gaussian, because the lags of one
# cell co-occur. The ridge and inter-ensemble moves scatter by
# ``exposure_idx`` and are skipped: with the joint ``(c, gamma, kappa)`` draw
# of the trend block there is no trend/level/exposure ridge left for them to
# cross, and under the anchored moderator the beta-nu scale ridge is gone.


def _lag_kernel_factor(state: LongBetState) -> Any:
    """The kernel factor at the current scale, ``L_kappa * s``."""
    L = state.lag_K_chol
    if state.lag_learn_scale:
        L = (L * jnp.sqrt(state.lag_scale2)).astype(jnp.float32)
    return L


def _lag_sign(kappa: Any) -> Any:
    """Canonical-orientation sign for the ``(kappa, nu) -> (-kappa, -nu)`` quotient:
    the cumulative response to a unit pulse, ``sum_l kappa_l``, is non-negative."""
    return jnp.where(jnp.sum(kappa) < 0.0, -1.0, 1.0)


def _draw_lag_kernel(
    key: Key[Array, ''],
    R: Float32[Array, ' n'],
    kappa: Any,
    nu_fit: Float32[Array, ' n'],
    forest_nu: Any,
    state: LongBetState,
    sigma2: Float32[Array, ''],
    conditional_precision: Float32[Array, ' n'] | None,
    L_kappa: Any,
    lin_fit: Float32[Array, ' n'] | None = None,
    delta: Any | None = None,
    ar1: Any = None,
) -> tuple[Any, Any, Any, Any, Any, Any]:
    """Draw the lag kernel ``kappa`` from its dense Gaussian conditional.

    The analogue of :func:`_draw_trajectory`: with ``d_it = nu_it + lin_it``
    the treatment term is ``d_it (lag_design @ kappa)_it``, linear in
    ``kappa`` with design rows ``d_it lag_design[it]``, so the conditional is
    ``N(P^-1 C, P^-1)``, ``P = (L L')^-1 + A``, with ``A = sum_it p_it d_it^2
    D_it D_it'`` dense and accumulated in float64; the draw is the whitened
    solve of ``longbet._ar1.sample_beta_gp_full``. ``ar1 = (phi, first, lam)``
    runs it on the Prais-Winsten transformed rows and residual. Without the
    anchored moderator the sign is then mapped onto ``sum(kappa) >= 0`` and
    the forest leaves (and linear part) flipped to match.

    Returns ``(kappa_new, nu_fit, forest_nu, R, delta, lin_fit)``.
    """
    obs_mask = state.obs_mask
    use_lin = lin_fit is not None
    d_vec = (nu_fit + lin_fit) if use_lin else nu_fit
    r_partial = R + d_vec * (state.lag_design @ kappa)
    if state.sample_beta:
        with enable_x64():
            f64 = jnp.float64
            rows = jnp.where(obs_mask[:, None],
                             jnp.asarray(d_vec, f64)[:, None] * jnp.asarray(state.lag_design, f64),
                             0.0)
            r = jnp.where(obs_mask, jnp.asarray(r_partial, f64), 0.0)
            if ar1 is not None:
                phi, first, lam = ar1
                rows = ar1_mod.transform(rows, jnp.asarray(phi, f64), first)
                r = ar1_mod.transform(r, jnp.asarray(phi, f64), first)
                p = jnp.where(obs_mask, jnp.asarray(lam, f64), 0.0)
            elif conditional_precision is None:
                p = jnp.where(obs_mask, jnp.reciprocal(jnp.asarray(sigma2, f64)), 0.0)
            else:
                p = jnp.where(obs_mask, jnp.asarray(conditional_precision, f64), 0.0)
            rows_p = rows * p[:, None]
            A_mat = rows.T @ rows_p
            C_vec = rows_p.T @ r
        kappa_new = ar1_mod.sample_beta_gp_full(key, A_mat, C_vec, L_kappa)
        if not state.anchor_moderator:
            sign = _lag_sign(kappa_new)
            kappa_new = sign * kappa_new
            nu_fit = sign * nu_fit
            forest_nu = eqx.tree_at(lambda f: f.leaf_tree, forest_nu, sign * forest_nu.leaf_tree)
            if use_lin:
                delta = sign * delta
                lin_fit = sign * lin_fit
            d_vec = (nu_fit + lin_fit) if use_lin else nu_fit
    else:
        kappa_new = kappa
    R = jnp.where(obs_mask, r_partial - d_vec * (state.lag_design @ kappa_new), 0.0)
    return kappa_new, nu_fit, forest_nu, R, delta, lin_fit


@enable_x64()
def _inv_gamma64(key, shape_a, rate_b):
    """Inverse-gamma draw with the division in float64 (as the horseshoe does)."""
    rate = jnp.asarray(rate_b, dtype=jnp.float64)
    g = random.gamma(key, jnp.asarray(shape_a, dtype=jnp.float32), dtype=jnp.float32)
    return rate / jnp.maximum(jnp.asarray(g, dtype=jnp.float64), 1e-30)


#: Bounds on sampled half-Cauchy variances, for float32 safety.
SCALE_VAR_MIN, SCALE_VAR_MAX = 1e-10, 1e10


@enable_x64()
def _half_cauchy_scale_step(key, sum_sq, dim, scale2, aux, global_scale):
    """Gibbs step for ``s^2`` with ``s ~ C+(0, A)`` given ``dim`` Gaussian
    coordinates with prior variance ``s^2`` times a fixed matrix and quadratic
    form ``sum_sq`` (Makalic and Schmidt's inverse-gamma mixture):
    ``s^2 | xi ~ IG((dim + 1)/2, 1/xi + sum_sq/2)``,
    ``xi | s^2 ~ IG(1, 1/A^2 + 1/s^2)``."""
    k1, k2 = random.split(key)
    A2 = float(global_scale) ** 2
    s2 = _inv_gamma64(k1, 0.5 * (dim + 1.0),
                      jnp.reciprocal(jnp.asarray(aux, jnp.float64))
                      + 0.5 * jnp.asarray(sum_sq, jnp.float64))
    s2 = jnp.clip(s2, SCALE_VAR_MIN, SCALE_VAR_MAX)
    xi = _inv_gamma64(k2, 1.0, 1.0 / A2 + jnp.reciprocal(s2))
    xi = jnp.clip(xi, SCALE_VAR_MIN, SCALE_VAR_MAX)
    return s2.astype(jnp.float32), xi.astype(jnp.float32)


@enable_x64()
def _linear_moderator_step(key, R, w, X, delta, lin_fit, tau2, tau_aux, global_scale,
                           obs_mask, unit_idx, N_units, gamma, sigma2, sigma_gamma2,
                           conditional_precision, collapse, ar1=None, *, local2, local_aux):
    """Exact draw of the moderator's linear coefficients and their horseshoe scales.

    Given everything else the treatment term is ``w_it (x_i' delta + trees)``,
    so the partial residual ``r = R + w * lin_fit`` is linear in ``delta`` with
    regressors ``w_it x_i``. With ``collapse`` the unit intercepts are
    integrated out (the same partially collapsed scheme as the moderator
    forest, and the sweep redraws ``gamma`` right after): per unit the
    residual covariance is ``diag(1/p) + sigma_gamma^2 1 1'``, whose inverse is
    ``diag(p) - p p' / (1/sigma_gamma^2 + sum p)``. Sums in float64.

    ``ar1 = (phi, first, lam)`` runs the same draw on the Prais-Winsten
    transformed rows and residual (``longbet._ar1``) with innovation precision
    ``lam``; ``gamma_i`` then enters transformed cell ``(i, t)`` with weight
    ``g_t``, the transform of a unit constant.
    """
    k_d, k_s = random.split(key)
    f64 = jnp.float64
    p = X.shape[1]
    wv = jnp.where(obs_mask, jnp.asarray(w, f64), 0.0)
    D = wv[:, None] * jnp.asarray(X, f64)
    r = jnp.asarray(R, f64) + wv * jnp.asarray(lin_fit, f64)
    if collapse:
        r = r + jnp.asarray(gamma, f64)[unit_idx]
    r = jnp.where(obs_mask, r, 0.0)
    gw = jnp.ones_like(r)
    if ar1 is not None:
        phi, first, lam = ar1
        D = ar1_mod.transform(D, jnp.asarray(phi, f64), first)
        r = ar1_mod.transform(r, jnp.asarray(phi, f64), first)
        gw = ar1_mod.transform(gw, jnp.asarray(phi, f64), first)
        pc = jnp.where(obs_mask, jnp.asarray(lam, f64), 0.0)
    elif conditional_precision is None:
        pc = jnp.where(obs_mask, jnp.reciprocal(jnp.asarray(sigma2, f64)), 0.0)
    else:
        pc = jnp.where(obs_mask, jnp.asarray(conditional_precision, f64), 0.0)
    Dp = D * pc[:, None]
    Q = D.T @ Dp
    b = Dp.T @ r
    if collapse:
        S = jnp.zeros((N_units, p), f64).at[unit_idx].add(Dp * gw[:, None])
        P = jnp.zeros(N_units, f64).at[unit_idx].add(pc * gw * gw)
        rs = jnp.zeros(N_units, f64).at[unit_idx].add(pc * gw * r)
        denom = jnp.reciprocal(jnp.asarray(sigma_gamma2, f64)) + P
        Q = Q - S.T @ (S / denom[:, None])
        b = b - S.T @ (rs / denom)
    prior_var = jnp.asarray(tau2, f64) * jnp.asarray(local2, f64)
    Q = Q + jnp.diag(jnp.reciprocal(prior_var))
    Q = 0.5 * (Q + Q.T)
    L = jnp.linalg.cholesky(Q)
    eta = random.normal(k_d, (p,), jnp.float32).astype(f64)
    delta_new = (jax.scipy.linalg.cho_solve((L, True), b)
                 + jax.scipy.linalg.solve_triangular(L.T, eta, lower=False)).astype(jnp.float32)
    lin_new = (jnp.asarray(X, jnp.float32) @ delta_new).astype(jnp.float32)
    R_new = jnp.where(obs_mask, R + w * (lin_fit - lin_new), 0.0).astype(jnp.float32)
    # Horseshoe: delta_j ~ N(0, tau^2 lambda_j^2), lambda_j ~ C+(0, 1), tau ~ C+(0, A)
    # (the trend block's Makalic-Schmidt sampler).
    local2_new, local_aux_new, tau2_new, aux_new = sample_horseshoe(
        k_s, delta_new, local2, local_aux, tau2, tau_aux, global_scale)
    return delta_new, lin_new, R_new, tau2_new, aux_new, local2_new, local_aux_new



@enable_x64()
def _interweave_gamma_scale(key, gamma, sigma_gamma2, P_unit, h_unit, prior_a, prior_b):
    """Ancillary update of ``s = sigma_gamma`` with ``eta = gamma / s`` fixed.

    ``P_unit`` is each unit's likelihood precision for its intercept and
    ``h_unit`` the precision-weighted sum of its residual with the intercept
    added back, so with ``eta`` fixed the log-likelihood in ``s`` is
    ``-A s^2 / 2 + B s`` (``A = sum eta^2 P``, ``B = sum eta h``). With the
    IG(a, b) prior on ``s^2`` and ``u = log s`` the exact target is

    ``l(u) = -A e^{2u} / 2 + B e^u - 2 a u - b e^{-2u}``

    (the ``s^N`` Jacobian of ``gamma = s eta`` cancels the normal's
    normalizer). Independence Metropolis from a Laplace fit (Newton mode), as
    for the AR(1) coefficient. Returns ``(gamma, sigma_gamma2, accepted)``.
    """
    f64 = jnp.float64
    s0 = jnp.sqrt(jnp.asarray(sigma_gamma2, f64))
    eta = jnp.asarray(gamma, f64) / s0
    A = jnp.sum(jnp.square(eta) * jnp.asarray(P_unit, f64))
    B = jnp.sum(eta * jnp.asarray(h_unit, f64))
    a, b = float(prior_a), float(prior_b)

    def ell(u):
        return -0.5 * A * jnp.exp(2.0 * u) + B * jnp.exp(u) - 2.0 * a * u - b * jnp.exp(-2.0 * u)

    def d12(u):
        g1 = -A * jnp.exp(2.0 * u) + B * jnp.exp(u) - 2.0 * a + 2.0 * b * jnp.exp(-2.0 * u)
        h1 = -2.0 * A * jnp.exp(2.0 * u) + B * jnp.exp(u) - 4.0 * b * jnp.exp(-2.0 * u)
        return g1, h1

    u0 = jnp.log(s0)
    m = u0
    for _ in range(8):
        g1, h1 = d12(m)
        step = jnp.where(h1 < -1e-12, -g1 / jnp.where(h1 < -1e-12, h1, -1.0), 0.5 * jnp.sign(g1))
        m = m + jnp.clip(step, -1.0, 1.0)
    _, h1 = d12(m)
    sd = jnp.where(h1 < -1e-8, 1.25 / jnp.sqrt(jnp.where(h1 < -1e-8, -h1, 1.0)), 0.5)
    k_p, k_a = random.split(key)
    prop = m + sd * random.normal(k_p, (), jnp.float32).astype(f64)

    def log_q(x):
        return -0.5 * jnp.square((x - m) / sd)

    log_ratio = ell(prop) + log_q(u0) - ell(u0) - log_q(prop)
    ok = jnp.isfinite(log_ratio) & (jnp.log(random.uniform(k_a, (), jnp.float32).astype(f64) + 1e-300) < log_ratio)
    s_new = jnp.exp(jnp.where(ok, prop, u0))
    return (s_new * eta).astype(jnp.float32), jnp.square(s_new).astype(jnp.float32), ok


def _ar1_single_step(key, state, conditional_precision=None):
    """One Gibbs sweep with stationary AR(1) innovations (``ar1_errors``).

    The same blocks, in the same order, as :func:`longbet_single_step` for the
    configuration AR(1) supports (continuous outcome, random intercept with
    both forests collapsed over it, balanced panel), with every likelihood
    term taken on the Prais-Winsten transformed residual (``longbet._ar1``),
    then ``sigma^2`` from the innovations and ``phi`` by its independence
    Metropolis step. ``conditional_precision`` may only be the tempered
    ``temperature / sigma^2`` of :func:`_tempered_single_step`; the sweep uses
    the scalar ``lam = temperature / sigma^2`` directly.
    """
    keys = random.split(key, 10)
    R = state.resid
    obs_mask = state.obs_mask
    first = state.ar1_first
    unit_idx, N_units, T = state.unit_idx, state.N_units, state.T_periods
    sigma2 = state.sigma2
    temp = (jnp.float32(1.0) if state.temperature is None
            else jnp.asarray(state.temperature, jnp.float32))
    lam = (temp / sigma2).astype(jnp.float32)
    phi = state.ar1_phi
    alpha = state.alpha
    alpha_safe = jnp.where(jnp.abs(alpha) < 1e-6, 1.0, alpha)

    def tf(v):
        return ar1_mod.transform(v, phi, first)

    def usum(v):
        return ar1_mod.unit_sum(v, unit_idx, N_units)

    g = tf(jnp.ones_like(R))                               # transformed unit constant
    G2 = usum(g * g)                                       # (N,), equal on a balanced panel
    G2_bar = jnp.mean(G2)

    K_chol = state.K_chol
    beta_scale2_cur = state.beta_scale2
    if state.learn_beta_scale and not state.anchor_moderator:
        K_chol = (state.K_chol * jnp.sqrt(beta_scale2_cur)).astype(jnp.float32)
    use_lin = state.use_linear_trt
    delta_new, lin_fit = state.delta, state.lin_fit
    tau_delta2_new, tau_delta_aux_new = state.tau_delta2, state.tau_delta_aux
    # Distributed-lag kernel: its factor at the current scale, and the
    # trajectory the block draws (kappa in place of beta).
    lag = state.use_lag_kernel
    lag_scale2_cur = state.lag_scale2
    L_kappa = _lag_kernel_factor(state) if lag else None
    traj_cur = state.kappa if lag else state.beta
    extra_coef_new = state.extra_coef

    def refresh_gamma(k, R_now, gamma_now, sg2):
        e = jnp.where(obs_mask, R_now + gamma_now[unit_idx], 0.0)
        V = jnp.reciprocal(lam * G2 + jnp.reciprocal(sg2))
        gam = V * lam * usum(g * tf(e)) + random.normal(k, (N_units,), jnp.float32) * jnp.sqrt(V)
        return gam.astype(jnp.float32), jnp.where(obs_mask, e - gam[unit_idx], 0.0)

    # 1. Prognostic forest, collapsed over gamma: per unit, g' T u / |g|^2 is an
    # observation of alpha mu(x_i) with variance 1/(lam |g|^2) + sigma_gamma^2.
    u_mu = jnp.where(obs_mask, R + state.gamma[unit_idx] + alpha * state.mu_fit, 0.0)
    u_mu_bar = usum(g * tf(u_mu)) / jnp.maximum(G2, 1e-12)
    R_mu_for_bartz = jnp.where(obs_mask, u_mu_bar[unit_idx] - alpha * state.mu_fit, 0.0)
    cell_prec_mu = (lam * G2_bar / T) / (1.0 + lam * G2_bar * state.sigma_gamma2)
    view_mu = State(
        _chain_anchor=state._chain_anchor, X=state.X, y=state.y, z=None, binary_indices=None,
        resid=_load(R_mu_for_bartz / alpha_safe, state.resid_unit),
        resid_unit=state.resid_unit, resid_eff_scale=state.resid_eff_scale,
        resid_inexact_integral=state.resid_inexact_integral,
        error_cov_inv=Wishart(nu=None, rate=None, value=jnp.square(alpha_safe) * cell_prec_mu),
        error_scale=None, prec_scale=state.prec_scale, inv_sdev_scale=state.inv_sdev_scale,
        inv_sdev_unit=state.inv_sdev_unit, n_non_missing=state.n_non_missing,
        sum_diag_prec_scale=state.sum_diag_prec_scale, forest=state.forest, config=state.config,
    )
    view_mu = bartz_step(keys[1], view_mu)
    view_mu = change_step(random.fold_in(key, 7101), view_mu)
    view_mu = change_internal_step(random.fold_in(key, 7105), view_mu)
    view_mu = regrow_step(random.fold_in(key, 7103), view_mu, 1)
    mu_fit = current_forest_fit(view_mu.forest)
    R = jnp.where(obs_mask, R - alpha * (mu_fit - state.mu_fit), 0.0)

    # 3. Treatment forest, collapsed over gamma, on the transformed weights.
    b_z = jnp.where(state.z_vec == 1.0, state.b1, state.b0)
    w = (state.lag_design @ state.kappa) if lag else b_z * state.beta[state.exposure_idx]
    nu_fit_curr = state.nu_fit
    if state.anchor_moderator:
        mean_nu_curr = jnp.mean(nu_fit_curr)
        nu_tilde_curr = nu_fit_curr - mean_nu_curr
        m_nu = state.forest_nu.leaf_tree.shape[0]
        forest_nu_in = eqx.tree_at(
            lambda f: f.leaf_tree, state.forest_nu,
            state.forest_nu.leaf_tree - mean_nu_curr / (m_nu * state.forest_nu.leaf_unit))
    else:
        nu_tilde_curr = nu_fit_curr
        forest_nu_in = state.forest_nu
    with enable_x64():
        f64 = jnp.float64
        w64 = jnp.where(obs_mask, jnp.asarray(w, f64), 0.0)
        g64 = ar1_mod.transform(jnp.ones(R.shape, f64), jnp.asarray(phi, f64), first)
        wt = ar1_mod.transform(w64, jnp.asarray(phi, f64), first)
        u_nu = jnp.where(obs_mask, jnp.asarray(R, f64) + jnp.asarray(state.gamma, f64)[unit_idx]
                         + w64 * jnp.asarray(nu_tilde_curr, f64), 0.0)
        ut = ar1_mod.transform(u_nu, jnp.asarray(phi, f64), first)
        W1 = ar1_mod.unit_sum(g64 * wt, unit_idx, N_units)
        W2 = ar1_mod.unit_sum(wt * wt, unit_idx, N_units)
        d_i = ar1_mod.unit_sum(g64 * g64, unit_idx, N_units) + jnp.reciprocal(
            jnp.asarray(lam, f64) * jnp.asarray(state.sigma_gamma2, f64))
        gamma_hat0 = ar1_mod.unit_sum(g64 * ut, unit_idx, N_units) / d_i
        score = ar1_mod.unit_sum(wt * ut, unit_idx, N_units) - W1 * gamma_hat0
        prec_i = W2 - jnp.square(W1) / d_i
        has = prec_i > 1e-12 * jnp.maximum(W2, 1e-300)
        kappa = jnp.where(has, prec_i / jnp.where(W2 > 0, W2, 1.0), 0.0)
        target = jnp.where(has, score / jnp.where(has, prec_i, 1.0), 0.0)
        kappa = kappa.astype(jnp.float32)
        target = target.astype(jnp.float32)
        wt32 = wt.astype(jnp.float32)
    w_eff_nu = wt32 * jnp.sqrt(kappa[unit_idx] * temp)
    attrs_nu, mask_nu = compute_treatment_scale_attrs(w_eff_nu, obs_mask, eps=WEIGHT_EPS)
    resid_nu_data = jnp.where(mask_nu, target[unit_idx] - nu_tilde_curr, 0.0)
    view_nu = State(
        _chain_anchor=state._chain_anchor, X=view_mu.X, y=view_mu.y, z=None, binary_indices=None,
        resid=_load(resid_nu_data, state.resid_unit_nu), resid_unit=state.resid_unit_nu,
        resid_eff_scale=state.resid_eff_scale_nu,
        resid_inexact_integral=state.resid_inexact_integral_nu,
        error_cov_inv=Wishart(nu=None, rate=None, value=jnp.reciprocal(sigma2)),
        error_scale=None, prec_scale=attrs_nu.prec_scale, inv_sdev_scale=attrs_nu.inv_sdev_scale,
        inv_sdev_unit=attrs_nu.inv_sdev_unit, n_non_missing=attrs_nu.n_non_missing,
        sum_diag_prec_scale=attrs_nu.sum_diag_prec_scale, forest=forest_nu_in,
        config=view_mu.config,
    )
    view_nu = refresh_prec_tree(view_nu)
    view_nu = bartz_step(keys[4], view_nu)
    view_nu = change_step(random.fold_in(key, 7102), view_nu)
    view_nu = change_internal_step(random.fold_in(key, 7106), view_nu)
    view_nu = regrow_step(random.fold_in(key, 7104), view_nu, 1)
    if state.anchor_moderator:
        nu_tilde_raw = current_forest_fit(view_nu.forest)
        shift_new = 1.0 - jnp.mean(nu_tilde_raw)
        m_nu = view_nu.forest.leaf_tree.shape[0]
        forest_nu = eqx.tree_at(
            lambda f: f.leaf_tree, view_nu.forest,
            view_nu.forest.leaf_tree + shift_new / (m_nu * view_nu.forest.leaf_unit))
        nu_fit = nu_tilde_raw + shift_new
    else:
        nu_fit = current_forest_fit(view_nu.forest)
        forest_nu = view_nu.forest
    R = jnp.where(obs_mask, R - w * (nu_fit - nu_fit_curr), 0.0)

    # 3b. Linear part of the moderator (collapsed over gamma, transformed).
    lin_local2_new, lin_local_aux_new = state.lin_local2, state.lin_local_aux
    if use_lin:
        (delta_new, lin_fit, R, tau_delta2_new, tau_delta_aux_new,
         lin_local2_new, lin_local_aux_new) = _linear_moderator_step(
            random.fold_in(key, 7201), R, w, state.lin_design, state.delta, state.lin_fit,
            state.tau_delta2, state.tau_delta_aux, state.linear_trt_global_scale,
            obs_mask, unit_idx, N_units, state.gamma, sigma2, state.sigma_gamma2,
            None, True, ar1=(phi, first, lam),
            local2=state.lin_local2, local_aux=state.lin_local_aux,
        )

    ridge_accepts_new = state.ridge_accepts
    if state.anchor_moderator:
        gamma_new = state.gamma
        beta_new = traj_cur
        sigma_gamma2_new = state.sigma_gamma2
        if lag and not state.use_trend_block:
            # Dense kappa draw (no sign quotient under the anchor), with the
            # unit intercepts refreshed on either side as for beta below.
            gamma_mid, R = refresh_gamma(random.fold_in(key, 7109), R, state.gamma, state.sigma_gamma2)
            beta_new, nu_fit, forest_nu, R, delta_new, lin_fit = _draw_lag_kernel(
                keys[3], R, state.kappa, nu_fit, forest_nu, state, sigma2, None, L_kappa,
                lin_fit=lin_fit if use_lin else None, delta=delta_new if use_lin else None,
                ar1=(phi, first, lam))
            gamma_new, R = refresh_gamma(keys[7], R, gamma_mid, state.sigma_gamma2)
            sigma_gamma2_new = _sample_inv_gamma(
                keys[8], jnp.float32(state.gamma_prior_a + N_units / 2.0),
                state.gamma_prior_b + 0.5 * jnp.sum(jnp.square(gamma_new)))
        elif not state.use_trend_block:
            gamma_mid, R = refresh_gamma(random.fold_in(key, 7109), R, state.gamma, state.sigma_gamma2)
            d_vec = b_z * (nu_fit + lin_fit if use_lin else nu_fit)
            r_partial = jnp.where(obs_mask, R + d_vec * state.beta[state.exposure_idx], 0.0)
            if state.sample_beta:
                rows = jnp.where(obs_mask[:, None],
                                 d_vec[:, None] * jax.nn.one_hot(state.exposure_idx, state.S_max + 1,
                                                                 dtype=jnp.float32), 0.0)
                with enable_x64():
                    rows_t = ar1_mod.transform(jnp.asarray(rows, jnp.float64), jnp.asarray(phi, jnp.float64), first)
                    r_t = ar1_mod.transform(jnp.asarray(r_partial, jnp.float64), jnp.asarray(phi, jnp.float64), first)
                    lam64 = jnp.asarray(lam, jnp.float64)
                    A_mat = lam64 * (rows_t.T @ rows_t)
                    C_vec = lam64 * (rows_t.T @ r_t)
                beta_new = ar1_mod.sample_beta_gp_full(keys[3], A_mat, C_vec, K_chol)
            R = jnp.where(obs_mask, r_partial - d_vec * beta_new[state.exposure_idx], 0.0)
            gamma_new, R = refresh_gamma(keys[7], R, gamma_mid, state.sigma_gamma2)
            sigma_gamma2_new = _sample_inv_gamma(
                keys[8], jnp.float32(state.gamma_prior_a + N_units / 2.0),
                state.gamma_prior_b + 0.5 * jnp.sum(jnp.square(gamma_new)))
    elif lag:
        # Unanchored lag kernel: gamma, the dense kappa draw with its sign
        # quotient, gamma again, and the scale-group move. When use_trend_block
        # is active, the joint (c, gamma, kappa) block removes the
        # trend/level/exposure ridge so ridge_scale_step and the inter-ensemble
        # moves (which scatter by exposure_idx) are skipped.
        gamma_mid, R = refresh_gamma(random.fold_in(key, 7109), R, state.gamma, state.sigma_gamma2)
        beta_new, nu_fit, forest_nu, R, delta_new, lin_fit = _draw_lag_kernel(
            keys[3], R, state.kappa, nu_fit, forest_nu, state, sigma2, None, L_kappa,
            lin_fit=lin_fit if use_lin else None, delta=delta_new if use_lin else None,
            ar1=(phi, first, lam))
        gamma_new, R = refresh_gamma(keys[7], R, gamma_mid, state.sigma_gamma2)
        sigma_gamma2_new = _sample_inv_gamma(
            keys[8], jnp.float32(state.gamma_prior_a + N_units / 2.0),
            state.gamma_prior_b + 0.5 * jnp.sum(jnp.square(gamma_new)))
        if state.sample_beta and not state.use_trend_block:
            beta_new, forest_nu, nu_fit, delta_r, lin_r, ridge_ok = ridge_scale_step(
                random.fold_in(key, 888), beta_new, forest_nu, nu_fit, L_kappa,
                state.leaf_prior_cov_inv_nu, proposal_sigma=RIDGE_PROPOSAL_SIGMA,
                delta=delta_new if use_lin else None, lin_fit=lin_fit if use_lin else None,
                delta_prior_prec=(jnp.reciprocal(tau_delta2_new * lin_local2_new)
                                  if use_lin else None),
                return_extras=True)
            if use_lin:
                delta_new, lin_fit = delta_r, lin_r
            if ridge_accepts_new is not None:
                ridge_accepts_new = ridge_accepts_new + ridge_ok.astype(jnp.float32)
        if state.sample_beta and (use_lin or state.lag_learn_scale):
            (beta_new, forest_nu, nu_fit, ls2_g, delta_g, lin_g, tau2_g,
             group_ok) = scale_group_step(
                random.fold_in(key, 7501), beta_new, forest_nu, nu_fit, state.lag_K_chol_sq,
                state.leaf_prior_cov_inv_nu,
                beta_scale2=lag_scale2_cur if state.lag_learn_scale else None,
                beta_scale_aux=state.lag_scale_aux if state.lag_learn_scale else None,
                delta=delta_new if use_lin else None, lin_fit=lin_fit if use_lin else None,
                tau_delta2=tau_delta2_new if use_lin else None,
                tau_delta_aux=tau_delta_aux_new if use_lin else None)
            if state.lag_learn_scale:
                lag_scale2_cur = ls2_g
                L_kappa = (state.lag_K_chol * jnp.sqrt(lag_scale2_cur)).astype(state.lag_K_chol.dtype)
            if use_lin:
                delta_new, lin_fit, tau_delta2_new = delta_g, lin_g, tau2_g
    else:
        # gamma given the collapsed updates, before anything conditions on it.
        gamma_mid, R = refresh_gamma(random.fold_in(key, 7109), R, state.gamma, state.sigma_gamma2)

        # 4. Trajectory beta: the transformed treatment rows touch two exposure
        # coordinates (S_it and S_{i,t-1}), so the data precision is banded.
        d_vec = b_z * (nu_fit + lin_fit if use_lin else nu_fit)
        r_partial = jnp.where(obs_mask, R + d_vec * state.beta[state.exposure_idx], 0.0)
        if state.sample_beta:
            rows = jnp.where(obs_mask[:, None],
                             d_vec[:, None] * jax.nn.one_hot(state.exposure_idx, state.S_max + 1,
                                                             dtype=jnp.float32), 0.0)
            with enable_x64():
                rows_t = ar1_mod.transform(jnp.asarray(rows, jnp.float64), jnp.asarray(phi, jnp.float64), first)
                r_t = ar1_mod.transform(jnp.asarray(r_partial, jnp.float64), jnp.asarray(phi, jnp.float64), first)
                lam64 = jnp.asarray(lam, jnp.float64)
                A_mat = lam64 * (rows_t.T @ rows_t)
                C_vec = lam64 * (rows_t.T @ r_t)
            beta_new = ar1_mod.sample_beta_gp_full(keys[3], A_mat, C_vec, K_chol)
            sign_beta_4 = jnp.where(jnp.sum(beta_new[1:]) < 0.0, -1.0, 1.0)
            beta_new = sign_beta_4 * beta_new
            nu_fit = sign_beta_4 * nu_fit
            forest_nu = eqx.tree_at(lambda f: f.leaf_tree, forest_nu, sign_beta_4 * forest_nu.leaf_tree)
            if use_lin:
                delta_new = sign_beta_4 * delta_new
                lin_fit = sign_beta_4 * lin_fit
            d_vec = b_z * (nu_fit + lin_fit if use_lin else nu_fit)
        else:
            beta_new = state.beta
        R = jnp.where(obs_mask, r_partial - d_vec * beta_new[state.exposure_idx], 0.0)

        # 6. gamma and sigma_gamma^2.
        gamma_new, R = refresh_gamma(keys[7], R, gamma_mid, state.sigma_gamma2)
        sigma_gamma2_new = _sample_inv_gamma(
            keys[8], jnp.float32(state.gamma_prior_a + N_units / 2.0),
            state.gamma_prior_b + 0.5 * jnp.sum(jnp.square(gamma_new)))

        # 7. Ridge move (the likelihood is invariant along it).
        if state.sample_beta:
            beta_new, forest_nu, nu_fit, delta_r, lin_r, ridge_ok = ridge_scale_step(
                random.fold_in(key, 888), beta_new, forest_nu, nu_fit, K_chol,
                state.leaf_prior_cov_inv_nu, proposal_sigma=RIDGE_PROPOSAL_SIGMA,
                delta=delta_new if use_lin else None, lin_fit=lin_fit if use_lin else None,
                delta_prior_prec=(jnp.reciprocal(tau_delta2_new * lin_local2_new)
                                  if use_lin else None),
                return_extras=True)
            if use_lin:
                delta_new, lin_fit = delta_r, lin_r
            if ridge_accepts_new is not None:
                ridge_accepts_new = ridge_accepts_new + ridge_ok.astype(jnp.float32)

        # 7b. Ridge move carrying the learned scales along (encourage-v2 only).
        if state.sample_beta and (use_lin or state.learn_beta_scale):
            (beta_new, forest_nu, nu_fit, bs2_g, delta_g, lin_g, tau2_g,
             group_ok) = scale_group_step(
                random.fold_in(key, 7501), beta_new, forest_nu, nu_fit, state.K_chol,
                state.leaf_prior_cov_inv_nu,
                beta_scale2=beta_scale2_cur if state.learn_beta_scale else None,
                beta_scale_aux=state.beta_scale_aux if state.learn_beta_scale else None,
                delta=delta_new if use_lin else None, lin_fit=lin_fit if use_lin else None,
                tau_delta2=tau_delta2_new if use_lin else None,
                tau_delta_aux=tau_delta_aux_new if use_lin else None)
            if state.learn_beta_scale:
                beta_scale2_cur = bs2_g
                K_chol = (state.K_chol * jnp.sqrt(beta_scale2_cur)).astype(jnp.float32)
            if use_lin:
                delta_new, lin_fit, tau_delta2_new = delta_g, lin_g, tau2_g

        # 8. Inter-ensemble moves on the transformed likelihood.
        if state.use_inter_ensemble_move:
            view_mu_forest, mu_fit, gamma_new, beta_new, R = inter_ensemble_beta_gamma_step(
                random.fold_in(key, 887), view_mu.forest, mu_fit, gamma_new, beta_new, R, alpha,
                nu_fit + lin_fit if use_lin else nu_fit, b_z, obs_mask, state.z_vec, state.time_idx,
                state.T_periods, unit_idx, N_units, state.exposure_idx, state.S_max, sigma2,
                view_mu.forest.leaf_prior_cov_inv, sigma_gamma2_new, K_chol, state.random_intercept,
                state.sample_beta, conditional_precision=None, temperature=temp,
                X_unified=state.X, y_vec=state.y, ar1=(phi, first))
            view_mu = eqx.tree_at(lambda v: v.forest, view_mu, view_mu_forest)
            w_transfer = b_z * beta_new[state.exposure_idx]
            view_mu_forest, mu_fit, forest_nu, nu_fit, R = inter_ensemble_transfer_step(
                random.fold_in(key, 889), view_mu.forest, mu_fit, forest_nu, nu_fit, R, alpha,
                w_transfer, obs_mask, state.z_vec, state.time_idx, state.T_periods, sigma2,
                view_mu.forest.leaf_prior_cov_inv, state.leaf_prior_cov_inv_nu,
                conditional_precision=None, temperature=temp, unit_idx=unit_idx, N_units=N_units,
                X_unified=state.X, y_vec=state.y, ar1=(phi, first))
            view_mu = eqx.tree_at(lambda v: v.forest, view_mu, view_mu_forest)

    # 8b. Trend block, jointly with gamma and beta (kappa under the lag
    # kernel), and with the fixed extra columns appended to the design.
    trend_coef_new = state.trend_coef
    trend_local2_new, trend_local_aux_new = state.trend_local2, state.trend_local_aux
    trend_global2_new, trend_global_aux_new = state.trend_global2, state.trend_global_aux
    if state.use_trend_block:
        if state.use_trend_horseshoe:
            col_scale = horseshoe_col_scales(
                jnp.asarray(state.trend_prior_scale_vec, dtype=jnp.float32),
                state.trend_num_period, state.trend_local2, state.trend_global2)
        else:
            col_scale = jnp.asarray(state.trend_prior_scale_vec, dtype=jnp.float32)
        q_tot = state.trend_design.shape[1]
        if state.num_extra_cols:
            H_blk = jnp.concatenate([state.trend_design, state.extra_design], axis=1)
            col_scale = jnp.concatenate([col_scale, jnp.asarray(state.extra_prior_scale_vec, jnp.float32)])
            coef_blk = jnp.concatenate([state.trend_coef, state.extra_coef])
        else:
            H_blk, coef_blk = state.trend_design, state.trend_coef
        d_blk = (nu_fit + lin_fit if use_lin else nu_fit) if lag else b_z * (nu_fit + lin_fit if use_lin else nu_fit)
        coef_blk, gamma_new, beta_new, R = ar1_mod.trend_block_step_ar1(
            random.fold_in(key, 9311), coef_blk, gamma_new, beta_new, R,
            H_blk, d_blk, obs_mask,
            unit_idx, N_units, state.exposure_idx, lam, sigma_gamma2_new,
            L_kappa if lag else K_chol, col_scale,
            phi, first, sample_beta=state.sample_beta, stats=state.ar1_trend_stats,
            lag_design=state.lag_design if lag else None,
            D_lag_K_chol=state.D_lag_K_chol if lag else None,
            lag_scale=jnp.sqrt(lag_scale2_cur) if (lag and state.lag_learn_scale) else None)
        if state.num_extra_cols:
            trend_coef_new, extra_coef_new = coef_blk[:q_tot], coef_blk[q_tot:]
        else:
            trend_coef_new = coef_blk
        if state.use_trend_horseshoe:
            (trend_local2_new, trend_local_aux_new, trend_global2_new,
             trend_global_aux_new) = sample_horseshoe(
                random.fold_in(key, 9313), trend_coef_new[state.trend_num_period:],
                state.trend_local2, state.trend_local_aux, state.trend_global2,
                state.trend_global_aux, state.trend_hs_global_scale)
        sigma_gamma2_new = _sample_inv_gamma(
            random.fold_in(key, 9312), jnp.float32(state.gamma_prior_a + N_units / 2.0),
            state.gamma_prior_b + 0.5 * jnp.sum(jnp.square(gamma_new)))

    if state.sample_beta and not state.anchor_moderator:
        sign_beta_end = _lag_sign(beta_new) if lag else jnp.where(jnp.sum(beta_new[1:]) < 0.0, -1.0, 1.0)
        beta_new = sign_beta_end * beta_new
        nu_fit = sign_beta_end * nu_fit
        forest_nu = eqx.tree_at(lambda f: f.leaf_tree, forest_nu, sign_beta_end * forest_nu.leaf_tree)
        if use_lin:
            delta_new = sign_beta_end * delta_new
            lin_fit = sign_beta_end * lin_fit

    beta_scale2_new, beta_scale_aux_new = state.beta_scale2, state.beta_scale_aux
    if state.learn_beta_scale and not state.anchor_moderator and not lag:
        beta_scale2_new, beta_scale_aux_new = _half_cauchy_scale_step(
            random.fold_in(key, 7301), beta_prior_quadform(beta_new, state.K_chol),
            float(state.S_max + 1), state.beta_scale2, state.beta_scale_aux,
            state.beta_scale_prior)
    # 8c'. Lag-kernel scale s^2 given kappa (prior-only; runs anchored or not).
    lag_scale2_new, lag_scale_aux_new = state.lag_scale2, state.lag_scale_aux
    if lag and state.lag_learn_scale and state.sample_beta:
        lag_scale2_new, lag_scale_aux_new = _half_cauchy_scale_step(
            random.fold_in(key, 7302), beta_prior_quadform(beta_new, state.lag_K_chol_sq),
            float(beta_new.shape[0]), state.lag_scale2, state.lag_scale_aux,
            state.lag_scale_prior)

    # 8d. Interweaving update of sigma_gamma on the transformed likelihood.
    interweave_accepts_new = state.interweave_accepts
    if state.interweave_gamma:
        e_g = jnp.where(obs_mask, R + gamma_new[unit_idx], 0.0)
        gamma_iw, sigma_gamma2_new, iw_ok = _interweave_gamma_scale(
            random.fold_in(key, 7601), gamma_new, sigma_gamma2_new, lam * G2,
            lam * usum(g * tf(e_g)), state.gamma_prior_a, state.gamma_prior_b)
        R = jnp.where(obs_mask, e_g - gamma_iw[unit_idx], 0.0)
        gamma_new = gamma_iw
        interweave_accepts_new = interweave_accepts_new + iw_ok.astype(jnp.float32)

    # 9. sigma^2 from the innovations (the tempered caller draws its own).
    if conditional_precision is None:
        m_obs = jnp.sum(obs_mask.astype(jnp.float32))
        ssr = jnp.sum(jnp.square(jnp.where(obs_mask, tf(R), 0.0)))
        sigma2_new = _sample_inv_gamma(keys[9], state.sigma_prior_a + m_obs / 2.0,
                                       state.sigma_prior_b + 0.5 * ssr)
        lam_phi = jnp.reciprocal(sigma2_new)
    else:
        sigma2_new = sigma2
        lam_phi = lam

    # 10. phi given the residual.
    phi_new, phi_ok = ar1_mod.sample_phi(random.fold_in(key, 7401), phi, R, first, obs_mask,
                                          lam_phi, temp, state.ar1_phi_max)

    error_cov_inv = eqx.tree_at(lambda w_: w_.value, state.error_cov_inv,
                                jnp.reciprocal(sigma2_new).astype(jnp.float32))
    extra_names = ["ar1_phi", "ar1_accepts"]
    extra_values = [phi_new, state.ar1_accepts + phi_ok.astype(jnp.float32)]
    if state.interweave_gamma:
        extra_names.append("interweave_accepts"); extra_values.append(interweave_accepts_new)
    if ridge_accepts_new is not None:
        extra_names.append("ridge_accepts"); extra_values.append(ridge_accepts_new)
    if use_lin:
        extra_names += ["delta", "lin_fit", "tau_delta2", "tau_delta_aux", "lin_local2",
                        "lin_local_aux"]
        extra_values += [delta_new, lin_fit, tau_delta2_new, tau_delta_aux_new, lin_local2_new,
                         lin_local_aux_new]
    if state.learn_beta_scale:
        extra_names += ["beta_scale2", "beta_scale_aux"]
        extra_values += [beta_scale2_new, beta_scale_aux_new]
    if lag:
        # ``beta_new`` carried the kernel through the sweep; the placeholder
        # ``beta`` is left as it was.
        extra_names.append("kappa"); extra_values.append(beta_new)
        beta_new = state.beta
        if state.lag_learn_scale:
            extra_names += ["lag_scale2", "lag_scale_aux"]
            extra_values += [lag_scale2_new, lag_scale_aux_new]
    if state.num_extra_cols:
        extra_names.append("extra_coef"); extra_values.append(extra_coef_new)
    state = eqx.tree_at(lambda s: tuple(getattr(s, n) for n in extra_names), state,
                        tuple(extra_values))
    return eqx.tree_at(
        lambda s: (s.X, s.y, s.config, s.resid, s.forest, s.resid_unit, s.resid_eff_scale,
                   s.resid_inexact_integral, s.error_cov_inv, s.n_non_missing,
                   s.forest_nu, s.resid_nu, s.resid_unit_nu, s.resid_eff_scale_nu,
                   s.resid_inexact_integral_nu, s.prec_scale_nu, s.inv_sdev_scale_nu,
                   s.inv_sdev_unit_nu, s.leaf_prior_cov_inv_nu, s.beta, s.gamma, s.sigma2,
                   s.sigma_gamma2, s.mu_fit, s.nu_fit, s.trend_coef, s.trend_local2,
                   s.trend_local_aux, s.trend_global2, s.trend_global_aux),
        state,
        (view_nu.X, view_nu.y, view_nu.config, R, view_mu.forest, view_mu.resid_unit,
         view_mu.resid_eff_scale, view_mu.resid_inexact_integral, error_cov_inv,
         view_mu.n_non_missing, forest_nu, view_nu.resid, view_nu.resid_unit,
         view_nu.resid_eff_scale, view_nu.resid_inexact_integral, view_nu.prec_scale,
         view_nu.inv_sdev_scale, view_nu.inv_sdev_unit, forest_nu.leaf_prior_cov_inv,
         beta_new, gamma_new, sigma2_new, sigma_gamma2_new, mu_fit, nu_fit, trend_coef_new,
         trend_local2_new, trend_local_aux_new, trend_global2_new, trend_global_aux_new),
        is_leaf=lambda x: x is None,
    )


@enable_x64()
@jax.jit
@enable_x64(False)
def longbet_single_step(
    key: Key[Array, ''],
    state: LongBetState,
    conditional_precision: Float32[Array, ' n'] | None = None,
) -> LongBetState:
    """Execute one Gibbs sweep on a single chain.

    Jitted at the boundary. ``bartz.mcmcstep.step`` donates its argument's
    buffers, so an un-jitted sweep would delete arrays the caller still holds --
    including ``X``, which is shared between chains. Inside ``jit`` the donation
    is an XLA aliasing hint and stays contained.

    Sweep order (plan section 4.9)::

        z (binary only) -> mu -> nu -> beta -> gamma -> sigma_gamma^2
                        -> ridge move -> subspace transfers
                        -> (c, gamma, beta) joint draw -> sigma^2

    Any order is valid provided each draw conditions on current values; this one
    keeps ``R`` consistent with the least bookkeeping and puts the variance draws
    last so they see the fully updated residual.  In particular ``gamma`` is drawn
    *after* the treatment fit, so a treated unit's post-adoption periods do not
    drag its baseline upward.

    Under ``baseline_cells='untreated'`` the same sweep runs as the alternating
    sampler of a two-module (cut) posterior: every baseline conditional --
    ``mu``, the pre- and post-forest ``gamma`` refreshes, the trend block,
    ``sigma_gamma^2`` and ``sigma^2`` -- is accumulated over the untreated
    cells only, and ``nu`` and ``beta`` are drawn from the treated cells given
    the baseline. The inter-ensemble transfers are skipped. The residual is
    still maintained on every observed cell, because the treatment module
    reads it on the treated ones.

    When ``conditional_precision`` is supplied, ``state.resid`` is the
    conditional (pseudo-response) residual from the full SUR likelihood.
    Every mean/latent update uses that cell-specific precision. The structural
    innovation variance ``state.sigma2`` is then LEFT UNCHANGED: its caller
    must update it from the structural innovation, not this pseudo-residual.
    None selects the scalar likelihood with its existing random-key schedule.
    Both paths refresh leaf caches when their observation weights change.
    """
    if state.ar1:
        if bool(state.baseline_untreated_only):
            raise NotImplementedError(
                "baseline_cells='untreated' is not available with ar1_errors=True "
                "(the AR(1) sweep fits the baseline on all cells)"
            )
        return _ar1_single_step(key, state, conditional_precision)
    keys = random.split(key, 10)

    # encouragement-v2 extensions: the effective trajectory kernel and the
    # moderator's linear part (both static no-ops when their flags are off).
    K_chol = state.K_chol
    beta_scale2_cur = state.beta_scale2
    if state.learn_beta_scale and not state.anchor_moderator:
        K_chol = (state.K_chol * jnp.sqrt(beta_scale2_cur)).astype(jnp.float32)
    use_lin = state.use_linear_trt
    delta_new = state.delta
    lin_fit = state.lin_fit
    tau_delta2_new = state.tau_delta2
    tau_delta_aux_new = state.tau_delta_aux
    lin_local2_new, lin_local_aux_new = state.lin_local2, state.lin_local_aux
    # Distributed-lag kernel (static no-op when off): its factor at the
    # current scale and the trajectory variable the sweep carries.
    lag = state.use_lag_kernel
    lag_scale2_cur = state.lag_scale2
    L_kappa = _lag_kernel_factor(state) if lag else None
    extra_coef_new = state.extra_coef

    R = state.resid  # full-model residual, data units, shape (n,)
    obs_mask = state.obs_mask
    sigma2 = state.sigma2
    alpha = state.alpha
    is_binary = state.num_categories == 2
    is_ordered = state.num_categories > 2
    cutpoints = state.cutpoints

    # Baseline module support (see LongBetConfig.baseline_cells). Every
    # conditional of mu, the trend block, gamma_i and sigma^2 is accumulated
    # over ``base_mask`` / ``base_counts``; the residual R is maintained on
    # ``obs_mask`` throughout, since the treatment module still reads it on the
    # treated cells. Under 'all' these are the same arrays and the sweep below
    # is, operation for operation, the one it was before the option existed.
    cut = bool(state.baseline_untreated_only)
    if cut and conditional_precision is not None:
        raise NotImplementedError(
            "baseline_cells='untreated' is not available with a coupled "
            "multi-outcome (SUR) likelihood"
        )
    base_mask = state.base_mask if cut else obs_mask
    base_counts = state.unit_counts_base if cut else state.unit_counts

    # Keep the global innovation precision factor and express conditional
    # weights relative to it. These temporary, chain-dependent attributes must
    # not replace the shared baseline observation mask in the returned state.
    conditional_attrs = None
    if conditional_precision is not None:
        conditional_attrs, _ = compute_treatment_scale_attrs(
            jnp.sqrt(conditional_precision * sigma2), obs_mask, eps=0.0
        )

    def weighted_sum(value):
        if conditional_precision is None:
            return jnp.sum(value) / sigma2
        return jnp.sum(jnp.where(obs_mask, value * conditional_precision, 0.0))

    # ------------------------------------------------------------------
    # 0. Probit latent z (binary outcomes only)
    # ------------------------------------------------------------------
    if is_ordered:
        # Distinct ordinal-only fold-in tags leave random.split(key, 10) and
        # every legacy binary/continuous transition unchanged.
        sd = 1. if conditional_precision is None else lax.rsqrt(conditional_precision)
        mean = state.z - R
        cutpoints = sample_cutpoints_marginalized(
            random.fold_in(key, 8112), cutpoints, state.y, mean, sd, obs_mask, state.cutpoint_prior_scale
        )
        current_z = sample_ordinal_latents(
            random.fold_in(key, 8111), state.y, mean, sd, cutpoints, obs_mask, state.z
        )
        R = jnp.where(obs_mask, R + current_z - state.z, 0.)
    elif is_binary:
        # step_z reconstructs the mean as z - resid * resid_unit, so it recovers
        # alpha*mu + b*beta*nu + gamma correctly without knowing the model has
        # more terms in it -- provided the residual it sees is the full-model one.
        view_z = State(
            _chain_anchor=state._chain_anchor,
            X=state.X,
            y=state.y,
            z=state.z,
            binary_indices=state.binary_indices,
            resid=_load(R, state.resid_unit),
            resid_unit=state.resid_unit,
            resid_eff_scale=state.resid_eff_scale,
            resid_inexact_integral=state.resid_inexact_integral,
            error_cov_inv=Wishart(nu=None, rate=None, value=jnp.float32(1.0)),
            # Binary structural sigma2 is 1, but its conditional latent
            # variance given continuous outcomes is generally smaller than 1.
            error_scale=(state.error_scale if conditional_attrs is None else
                         lax.rsqrt(conditional_precision)),
            prec_scale=(state.prec_scale if conditional_attrs is None else
                        conditional_attrs.prec_scale),
            inv_sdev_scale=(state.inv_sdev_scale if conditional_attrs is None else
                            conditional_attrs.inv_sdev_scale),
            inv_sdev_unit=(state.inv_sdev_unit if conditional_attrs is None else
                           conditional_attrs.inv_sdev_unit),
            n_non_missing=state.n_non_missing,
            sum_diag_prec_scale=(state.sum_diag_prec_scale if conditional_attrs is None else
                                 conditional_attrs.sum_diag_prec_scale),
            forest=state.forest,
            config=state.config,
        )
        view_z = step_z(keys[0], view_z)
        current_z = view_z.z
        R = R + jnp.where(obs_mask, current_z - state.z, 0.0)
    else:
        current_z = state.z

    # ------------------------------------------------------------------
    # 1. Prognostic forest mu (collapsed over unit random intercepts gamma_i)
    # ------------------------------------------------------------------
    # Target is Y - b*beta*nu - trend = alpha*mu + gamma_i + eps_it.
    # Conditioning on fixed gamma_i pins mu(X_i) and nu(X_i) to the previous
    # sweep's unit intercepts with precision T/sigma^2 (~20x the marginal
    # precision 1/(sigma_gamma^2 + sigma^2/T)), freezing tree partitions.
    # Integrating out gamma_i ~ N(0, sigma_gamma^2) analytically yields exact
    # marginal sufficient statistics for every unit-level tree move.
    alpha_safe = jnp.where(jnp.abs(alpha) < 1e-6, 1.0, alpha)
    collapse_gamma_mu = state.random_intercept and state.collapse_gamma_mu
    # Under the cut gamma_i is a module-1 parameter: its conditional is the
    # untreated-cell posterior and does not involve nu, so the treatment forest
    # conditions on the current draw instead of integrating it out. The
    # mu <-> gamma coupling that motivates collapsing stays inside module 1.
    collapse_gamma_nu = state.random_intercept and state.collapse_gamma_nu and not cut
    conditional_attrs_mu = conditional_attrs

    if collapse_gamma_mu and (conditional_precision is None):
        R_no_gamma = jnp.where(base_mask, R + state.gamma[state.unit_idx], 0.0)
        u_mu = jnp.where(base_mask, R_no_gamma + alpha * state.mu_fit, 0.0)
        u_mu_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(u_mu)
        )
        n_i_safe = jnp.maximum(base_counts, 1.0)
        u_mu_bar = u_mu_sum / n_i_safe
        kappa_mu = jnp.reciprocal(1.0 + base_counts * (state.sigma_gamma2 / sigma2))
        if state.has_unbalanced_units or cut:
            conditional_attrs_mu, _ = compute_treatment_scale_attrs(
                jnp.sqrt(kappa_mu[state.unit_idx]),
                base_mask,
                eps=0.0,
            )
            R_mu_for_bartz = jnp.where(
                base_mask, u_mu_bar[state.unit_idx] - alpha * state.mu_fit, 0.0
            )
            error_cov_inv_mu = jnp.square(alpha_safe) / sigma2
        else:
            n_bar = jnp.mean(n_i_safe)
            kappa_mu_bar = jnp.reciprocal(1.0 + n_bar * (state.sigma_gamma2 / sigma2))
            # Observed units in this branch have identical counts, hence this
            # ratio is exactly one, including at very small marginal precision.
            ratio_mu = jnp.ones_like(R)
            R_mu_for_bartz = jnp.where(
                obs_mask, ratio_mu * (u_mu_bar[state.unit_idx] - alpha * state.mu_fit), 0.0
            )
            error_cov_inv_mu = jnp.square(alpha_safe) * kappa_mu_bar / sigma2
    elif collapse_gamma_mu and (conditional_precision is not None):
        cell_prec_mu = jnp.where(obs_mask, conditional_precision, 0.0)
        P_i_mu = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_mu)
        )
        R_no_gamma = jnp.where(obs_mask, R + state.gamma[state.unit_idx], 0.0)
        u_mu = jnp.where(obs_mask, R_no_gamma + alpha * state.mu_fit, 0.0)
        u_mu_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_mu * u_mu)
        )
        u_mu_bar = u_mu_sum / jnp.maximum(P_i_mu, 1e-12)
        kappa_mu = jnp.reciprocal(1.0 + P_i_mu * state.sigma_gamma2)
        conditional_attrs_mu, _ = compute_treatment_scale_attrs(
            jnp.sqrt(conditional_precision * kappa_mu[state.unit_idx] * sigma2),
            obs_mask,
            eps=0.0,
        )
        R_mu_for_bartz = jnp.where(
            obs_mask, u_mu_bar[state.unit_idx] - alpha * state.mu_fit, 0.0
        )
        error_cov_inv_mu = jnp.square(alpha_safe) / sigma2
    elif cut:
        # Conditional on gamma_i (mu splits on time-varying inputs); the
        # treated cells are masked out of the forest's precision instead of
        # relying on the init-time mask, which spans every observed cell.
        conditional_attrs_mu, _ = compute_treatment_scale_attrs(
            jnp.ones_like(R), base_mask, eps=0.0
        )
        R_mu_for_bartz = jnp.where(base_mask, R, 0.0)
        error_cov_inv_mu = jnp.square(alpha_safe) / sigma2
    else:
        R_mu_for_bartz = R
        error_cov_inv_mu = jnp.square(alpha_safe) / sigma2

    view_mu = State(
        _chain_anchor=state._chain_anchor,
        X=state.X,
        y=state.y,
        # The full-model latent response was sampled above. bartz.step would
        # sample it AGAIN after the tree update if z were non-None, corrupting
        # fit-by-residual-differencing and leaving state.z inconsistent with R.
        z=None,
        binary_indices=None,
        resid=_load(R_mu_for_bartz / alpha_safe, state.resid_unit),
        resid_unit=state.resid_unit,
        resid_eff_scale=state.resid_eff_scale,
        resid_inexact_integral=state.resid_inexact_integral,
        error_cov_inv=Wishart(
            nu=None, rate=None, value=error_cov_inv_mu
        ),
        error_scale=None,
        prec_scale=(state.prec_scale if conditional_attrs_mu is None else
                    conditional_attrs_mu.prec_scale),
        inv_sdev_scale=(state.inv_sdev_scale if conditional_attrs_mu is None else
                        conditional_attrs_mu.inv_sdev_scale),
        inv_sdev_unit=(state.inv_sdev_unit if conditional_attrs_mu is None else
                       conditional_attrs_mu.inv_sdev_unit),
        n_non_missing=state.n_non_missing,
        sum_diag_prec_scale=(state.sum_diag_prec_scale if conditional_attrs_mu is None else
                             conditional_attrs_mu.sum_diag_prec_scale),
        forest=state.forest,
        config=state.config,
    )

    if conditional_attrs_mu is not None:
        view_mu = refresh_prec_tree(view_mu)
    view_mu = bartz_step(keys[1], view_mu)
    view_mu = change_step(random.fold_in(key, 7101), view_mu)
    view_mu = change_internal_step(random.fold_in(key, 7105), view_mu)
    view_mu = regrow_step(random.fold_in(key, 7103), view_mu, 1)
    mu_fit = current_forest_fit(view_mu.forest)
    R = jnp.where(obs_mask, R-alpha*(mu_fit-state.mu_fit), 0.)

    # When mu collapses gamma_i (collapse_gamma_mu=True) but nu updates conditionally
    # on gamma_i (collapse_gamma_nu=False, e.g., with time-varying moderators x_trt_tv),
    # refresh gamma_i from p(gamma | mu_new, nu_old, ...) before Step 3 so that
    # (mu_new, gamma_pre_nu) is an exact joint draw from p(mu, gamma | nu_old, ...)
    # and Step 3 conditions on a valid gamma draw (Van Dyk & Park, 2008).
    if collapse_gamma_mu and (not collapse_gamma_nu) and (conditional_precision is None):
        e_mid = R + state.gamma[state.unit_idx]
        e_mid_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(jnp.where(base_mask, e_mid, 0.0))
        ) / sigma2
        V_mid = jnp.reciprocal(base_counts / sigma2 + jnp.reciprocal(state.sigma_gamma2))
        gamma_pre_nu = V_mid * e_mid_sum + random.normal(
            random.fold_in(key, 7108), (state.N_units,), dtype=jnp.float32
        ) * jnp.sqrt(V_mid)
        R = jnp.where(obs_mask, e_mid - gamma_pre_nu[state.unit_idx], 0.0)
        state = eqx.tree_at(lambda s: s.gamma, state, gamma_pre_nu)
    elif collapse_gamma_mu and (not collapse_gamma_nu) and (conditional_precision is not None):
        e_mid = R + state.gamma[state.unit_idx]
        cell_prec_mid = jnp.where(obs_mask, conditional_precision, 0.0)
        e_mid_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_mid * e_mid)
        )
        P_i_mid = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_mid)
        )
        V_mid = jnp.reciprocal(P_i_mid + jnp.reciprocal(state.sigma_gamma2))
        gamma_pre_nu = V_mid * e_mid_sum + random.normal(
            random.fold_in(key, 7108), (state.N_units,), dtype=jnp.float32
        ) * jnp.sqrt(V_mid)
        R = jnp.where(obs_mask, e_mid - gamma_pre_nu[state.unit_idx], 0.0)
        state = eqx.tree_at(lambda s: s.gamma, state, gamma_pre_nu)

    # ------------------------------------------------------------------
    # 2. Prognostic scale alpha
    # ------------------------------------------------------------------
    alpha_new = alpha  # fixed at one; kept as a state field for the archive layout

    # ------------------------------------------------------------------
    # 3. Treatment forest nu (collapsed over unit random intercepts gamma_i)
    # ------------------------------------------------------------------
    b_z = jnp.where(state.z_vec == 1.0, state.b1, state.b0)
    if lag:
        w = state.lag_design @ state.kappa
    else:
        beta_expanded = state.beta[state.exposure_idx]
        w = b_z * beta_expanded

    gamma_curr = state.gamma
    nu_fit_curr = state.nu_fit
    if state.anchor_moderator:
        mean_nu_curr = jnp.mean(nu_fit_curr)
        nu_tilde_curr = nu_fit_curr - mean_nu_curr
        m_nu = state.forest_nu.leaf_tree.shape[0]
        forest_nu_curr = eqx.tree_at(
            lambda f: f.leaf_tree,
            state.forest_nu,
            state.forest_nu.leaf_tree - mean_nu_curr / (m_nu * state.forest_nu.leaf_unit),
        )
    else:
        nu_tilde_curr = nu_fit_curr
        forest_nu_curr = state.forest_nu

    if collapse_gamma_nu and (conditional_precision is None):
        R_no_gamma = jnp.where(obs_mask, R + gamma_curr[state.unit_idx], 0.0)
        u_nu = jnp.where(obs_mask, R_no_gamma + w * nu_tilde_curr, 0.0)
        w_obs = jnp.where(obs_mask, w, 0.0)
        W1_i = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(w_obs)
        )
        W2_i = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(jnp.square(w_obs))
        )
        d_i = state.unit_counts + sigma2 / state.sigma_gamma2
        kappa_nu_i = jnp.minimum(
            1.0 - jnp.square(W1_i) / (jnp.maximum(W2_i, 1e-12) * jnp.maximum(d_i, 1e-6)),
            1.0,
        )
        u_nu_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(u_nu)
        )
        gamma_hat0_i = u_nu_sum / jnp.maximum(d_i, 1e-6)
        u_nu_tilde = jnp.where(obs_mask, u_nu - gamma_hat0_i[state.unit_idx], 0.0)
        score_nu_i = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(w_obs * u_nu_tilde)
        )
        prec_nu_i = kappa_nu_i * W2_i
        nu_target_i = score_nu_i / jnp.maximum(prec_nu_i, 1e-12)
        # Compute the prior precision in float64 too, so an extreme but finite
        # variance ratio does not underflow before entering the stable branch.
        with enable_x64():
            prior_precision = jnp.asarray(sigma2, jnp.float64) / jnp.asarray(state.sigma_gamma2, jnp.float64)
        kappa_nu_i, nu_target_i, repaired_nu = _stabilize_collapsed_nu_stats(
            kappa_nu_i, nu_target_i, W2_i, d_i, u_nu, w, obs_mask,
            state.unit_idx, state.N_units, prior_precision,
        )
        w_eff_nu = w * jnp.sqrt(kappa_nu_i[state.unit_idx])
        attrs_nu, mask_nu = compute_treatment_scale_attrs(
            w_eff_nu, obs_mask,
            eps=jnp.where(repaired_nu[state.unit_idx], 0.0, WEIGHT_EPS),
        )
        resid_nu_data = jnp.where(mask_nu, nu_target_i[state.unit_idx] - nu_tilde_curr, 0.0)
    elif collapse_gamma_nu and (conditional_precision is not None):
        cell_prec_nu = jnp.where(obs_mask, conditional_precision, 0.0)
        R_no_gamma = jnp.where(obs_mask, R + gamma_curr[state.unit_idx], 0.0)
        u_nu = jnp.where(obs_mask, R_no_gamma + w * nu_tilde_curr, 0.0)
        w_obs = jnp.where(obs_mask, w, 0.0)
        P_i_nu = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_nu)
        )
        W1_i = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_nu * w_obs)
        )
        W2_i = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_nu * jnp.square(w_obs))
        )
        d_i = P_i_nu + jnp.reciprocal(state.sigma_gamma2)
        kappa_nu_i = jnp.minimum(
            1.0 - jnp.square(W1_i) / (jnp.maximum(W2_i, 1e-12) * jnp.maximum(d_i, 1e-6)),
            1.0,
        )
        u_nu_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_nu * u_nu)
        )
        gamma_hat0_i = u_nu_sum / jnp.maximum(d_i, 1e-6)
        u_nu_tilde = jnp.where(obs_mask, u_nu - gamma_hat0_i[state.unit_idx], 0.0)
        score_nu_i = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_nu * w_obs * u_nu_tilde)
        )
        prec_nu_i = kappa_nu_i * W2_i
        nu_target_i = score_nu_i / jnp.maximum(prec_nu_i, 1e-12)
        with enable_x64():
            prior_precision = jnp.reciprocal(jnp.asarray(state.sigma_gamma2, jnp.float64))
        kappa_nu_i, nu_target_i, repaired_nu = _stabilize_collapsed_nu_stats(
            kappa_nu_i, nu_target_i, W2_i, d_i, u_nu, w, obs_mask,
            state.unit_idx, state.N_units, prior_precision, cell_prec_nu,
        )
        w_eff_nu = w * jnp.sqrt(kappa_nu_i[state.unit_idx] * cell_prec_nu * sigma2)
        attrs_nu, mask_nu = compute_treatment_scale_attrs(
            w_eff_nu, obs_mask,
            eps=jnp.where(repaired_nu[state.unit_idx], 0.0, WEIGHT_EPS),
        )
        resid_nu_data = jnp.where(mask_nu, nu_target_i[state.unit_idx] - nu_tilde_curr, 0.0)
    else:
        attrs_nu, mask_nu = compute_treatment_scale_attrs(w, obs_mask, eps=WEIGHT_EPS)
        if conditional_precision is not None:
            attrs_nu, _ = compute_treatment_scale_attrs(
                w * jnp.sqrt(conditional_precision * sigma2), mask_nu, eps=0.0
            )
        w_safe = jnp.where(mask_nu, w, 1.0)
        resid_nu_data = jnp.where(mask_nu, R / w_safe, 0.0)

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
        error_cov_inv=Wishart(nu=None, rate=None, value=jnp.reciprocal(sigma2)),
        error_scale=None,
        prec_scale=attrs_nu.prec_scale,
        inv_sdev_scale=attrs_nu.inv_sdev_scale,
        inv_sdev_unit=attrs_nu.inv_sdev_unit,
        n_non_missing=attrs_nu.n_non_missing,
        sum_diag_prec_scale=attrs_nu.sum_diag_prec_scale,
        forest=forest_nu_curr,
        config=view_mu.config,
    )

    view_nu = refresh_prec_tree(view_nu)
    view_nu = bartz_step(keys[4], view_nu)
    view_nu = change_step(random.fold_in(key, 7102), view_nu)
    view_nu = change_internal_step(random.fold_in(key, 7106), view_nu)
    view_nu = regrow_step(random.fold_in(key, 7104), view_nu, 1)
    if state.anchor_moderator:
        nu_tilde_raw = current_forest_fit(view_nu.forest)
        shift_new = 1.0 - jnp.mean(nu_tilde_raw)
        m_nu = view_nu.forest.leaf_tree.shape[0]
        forest_nu = eqx.tree_at(
            lambda f: f.leaf_tree,
            view_nu.forest,
            view_nu.forest.leaf_tree + shift_new / (m_nu * view_nu.forest.leaf_unit),
        )
        nu_fit = nu_tilde_raw + shift_new
    else:
        nu_fit = current_forest_fit(view_nu.forest)
        forest_nu = view_nu.forest
    private_delta = nu_fit - nu_fit_curr

    # Apply the physical fit delta even on zero/negligible-weight cells. Avoid
    # subtracting large R/w values to recover a small leaf change.
    R = jnp.where(obs_mask, R-w*private_delta, 0.0)

    # 3b. Linear part of the moderator (collapsed over gamma when nu is; the
    # gamma refresh below then restores a valid joint draw).
    if use_lin:
        (delta_new, lin_fit, R, tau_delta2_new, tau_delta_aux_new,
         lin_local2_new, lin_local_aux_new) = _linear_moderator_step(
            random.fold_in(key, 7201), R, w, state.lin_design, state.delta, state.lin_fit,
            state.tau_delta2, state.tau_delta_aux, state.linear_trt_global_scale,
            obs_mask, state.unit_idx, state.N_units, state.gamma, sigma2,
            state.sigma_gamma2, conditional_precision,
            bool(collapse_gamma_nu),
            local2=state.lin_local2, local_aux=state.lin_local_aux,
        )

    # Refresh gamma_i immediately after the collapsed forest updates so that Step 4
    # (sample_beta_gp) conditions on unit intercepts consistent with (mu_fit, nu_fit).
    if (collapse_gamma_mu or collapse_gamma_nu) and (conditional_precision is None):
        e_pre = R + state.gamma[state.unit_idx]
        e_pre_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(jnp.where(base_mask, e_pre, 0.0))
        ) / sigma2
        V_pre = jnp.reciprocal(base_counts / sigma2 + jnp.reciprocal(state.sigma_gamma2))
        gamma_mid = V_pre * e_pre_sum + random.normal(
            random.fold_in(key, 7109), (state.N_units,), dtype=jnp.float32
        ) * jnp.sqrt(V_pre)
        R = jnp.where(obs_mask, e_pre - gamma_mid[state.unit_idx], 0.0)
        state = eqx.tree_at(lambda s: s.gamma, state, gamma_mid)
    elif (collapse_gamma_mu or collapse_gamma_nu) and (conditional_precision is not None):
        e_pre = R + state.gamma[state.unit_idx]
        cell_prec_pre = jnp.where(obs_mask, conditional_precision, 0.0)
        e_pre_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_pre * e_pre)
        )
        unit_prec_pre = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(cell_prec_pre)
        )
        V_pre = jnp.reciprocal(unit_prec_pre + jnp.reciprocal(state.sigma_gamma2))
        gamma_mid = V_pre * e_pre_sum + random.normal(
            random.fold_in(key, 7109), (state.N_units,), dtype=jnp.float32
        ) * jnp.sqrt(V_pre)
        R = jnp.where(obs_mask, e_pre - gamma_mid[state.unit_idx], 0.0)
        state = eqx.tree_at(lambda s: s.gamma, state, gamma_mid)

    # ------------------------------------------------------------------
    # 4. GP trajectory beta_S (dense lag kernel kappa under exposure_lags)
    # ------------------------------------------------------------------
    delta_prev, lin_prev = delta_new, lin_fit
    if lag:
        beta_new, nu_fit, forest_nu, R, delta_new, lin_fit = _draw_lag_kernel(
            keys[3], R, state.kappa, nu_fit, forest_nu, state, sigma2,
            conditional_precision, L_kappa,
            lin_fit=lin_fit if use_lin else None, delta=delta_new if use_lin else None,
        )
    else:
        beta_new, nu_fit, forest_nu, R, delta_new, lin_fit = _draw_trajectory(
            keys[3], R, state.beta, nu_fit, forest_nu, b_z, state, sigma2,
            conditional_precision, K_chol=K_chol,
            lin_fit=lin_fit if use_lin else None, delta=delta_new if use_lin else None,
        )
    if not use_lin:
        delta_new, lin_fit = delta_prev, lin_prev

    # ------------------------------------------------------------------
    # 5. Coding scales: fixed at b0 = 0, b1 = 1 (treatment-only coding)
    # ------------------------------------------------------------------
    b0_new, b1_new = state.b0, state.b1

    # ------------------------------------------------------------------
    # 6. Unit intercepts gamma_i and sigma_gamma^2
    # ------------------------------------------------------------------
    gamma_new, sigma_gamma2_new = state.gamma, state.sigma_gamma2
    if state.random_intercept:
        e = R + state.gamma[state.unit_idx]
        e_sum = (
            jnp.zeros(state.N_units, dtype=jnp.float32)
            .at[state.unit_idx]
            .add(jnp.where(base_mask, e, 0.0))
        )
        if conditional_precision is not None:
            e_sum = jnp.zeros(state.N_units, dtype=jnp.float32).at[state.unit_idx].add(
                jnp.where(obs_mask, e * conditional_precision, 0.0)
            )
            unit_precision = jnp.zeros(state.N_units, dtype=jnp.float32).at[state.unit_idx].add(
                jnp.where(obs_mask, conditional_precision, 0.0)
            )
        else:
            e_sum = e_sum / sigma2
            unit_precision = base_counts / sigma2
        V = jnp.reciprocal(
            unit_precision + jnp.reciprocal(state.sigma_gamma2)
        )
        gamma_new = V * e_sum + random.normal(
            keys[7], (state.N_units,), dtype=jnp.float32
        ) * jnp.sqrt(V)
        R = jnp.where(obs_mask, e - gamma_new[state.unit_idx], 0.0)

        sigma_gamma2_new = _sample_inv_gamma(
            keys[8],
            jnp.float32(state.gamma_prior_a + state.N_units / 2.0),
            state.gamma_prior_b + 0.5 * jnp.sum(jnp.square(gamma_new)),
        )

    # ------------------------------------------------------------------
    # 7. Metropolis move along the beta-nu scale ridge
    # ------------------------------------------------------------------
    # Skipped under the lag kernel when use_trend_block=True: the joint
    # (c, gamma, kappa) draw below already updates kappa jointly with the
    # baseline, while scale_group_step below carries the learned scales along.
    ridge_accepts_new = state.ridge_accepts
    if state.sample_beta and not state.anchor_moderator and not (lag and state.use_trend_block):
        beta_new, forest_nu, nu_fit, delta_r, lin_r, ridge_ok = ridge_scale_step(
            random.fold_in(key, 888),
            beta_new,
            forest_nu,
            nu_fit,
            L_kappa if lag else K_chol,
            state.leaf_prior_cov_inv_nu,
            proposal_sigma=RIDGE_PROPOSAL_SIGMA,
            delta=delta_new if use_lin else None,
            lin_fit=lin_fit if use_lin else None,
            delta_prior_prec=(jnp.reciprocal(tau_delta2_new * lin_local2_new)
                              if use_lin else None),
            return_extras=True,
        )
        if use_lin:
            delta_new, lin_fit = delta_r, lin_r
        if ridge_accepts_new is not None:
            ridge_accepts_new = ridge_accepts_new + ridge_ok.astype(jnp.float32)

    # 7b. Ridge move carrying the learned scales along (encourage-v2 only).
    if lag:
        if (state.sample_beta and not state.anchor_moderator
                and (use_lin or state.lag_learn_scale)):
            (beta_new, forest_nu, nu_fit, ls2_g, delta_g, lin_g, tau2_g,
             group_ok) = scale_group_step(
                random.fold_in(key, 7501), beta_new, forest_nu, nu_fit, state.lag_K_chol_sq,
                state.leaf_prior_cov_inv_nu,
                beta_scale2=lag_scale2_cur if state.lag_learn_scale else None,
                beta_scale_aux=state.lag_scale_aux if state.lag_learn_scale else None,
                delta=delta_new if use_lin else None, lin_fit=lin_fit if use_lin else None,
                tau_delta2=tau_delta2_new if use_lin else None,
                tau_delta_aux=tau_delta_aux_new if use_lin else None)
            if state.lag_learn_scale:
                lag_scale2_cur = ls2_g
                L_kappa = (state.lag_K_chol * jnp.sqrt(lag_scale2_cur)).astype(state.lag_K_chol.dtype)
            if use_lin:
                delta_new, lin_fit, tau_delta2_new = delta_g, lin_g, tau2_g
    elif (state.sample_beta and not state.anchor_moderator
            and (use_lin or state.learn_beta_scale)):
        (beta_new, forest_nu, nu_fit, bs2_g, delta_g, lin_g, tau2_g,
         group_ok) = scale_group_step(
            random.fold_in(key, 7501), beta_new, forest_nu, nu_fit, state.K_chol,
            state.leaf_prior_cov_inv_nu,
            beta_scale2=beta_scale2_cur if state.learn_beta_scale else None,
            beta_scale_aux=state.beta_scale_aux if state.learn_beta_scale else None,
            delta=delta_new if use_lin else None, lin_fit=lin_fit if use_lin else None,
            tau_delta2=tau_delta2_new if use_lin else None,
            tau_delta_aux=tau_delta_aux_new if use_lin else None)
        if state.learn_beta_scale:
            beta_scale2_cur = bs2_g
            K_chol = (state.K_chol * jnp.sqrt(beta_scale2_cur)).astype(jnp.float32)
        if use_lin:
            delta_new, lin_fit, tau_delta2_new = delta_g, lin_g, tau2_g

    # ------------------------------------------------------------------
    # 8. Inter-ensemble residual-transfer Metropolis move (mu <-> gamma, beta, nu)
    # ------------------------------------------------------------------
    # Not run under the cut: both moves rebalance fitted mass between the
    # baseline and the treatment module using every cell at once.
    if state.use_inter_ensemble_move and not cut and not state.anchor_moderator and not lag:
        temp = (
            jnp.float32(1.0)
            if state.temperature is None
            else jnp.asarray(state.temperature, dtype=jnp.float32)
        )
        view_mu_forest, mu_fit, gamma_new, beta_new, R = inter_ensemble_beta_gamma_step(
            random.fold_in(key, 887),
            view_mu.forest,
            mu_fit,
            gamma_new,
            beta_new,
            R,
            alpha_new,
            nu_fit + lin_fit if use_lin else nu_fit,
            b_z,
            obs_mask,
            state.z_vec,
            state.time_idx,
            state.T_periods,
            state.unit_idx,
            state.N_units,
            state.exposure_idx,
            state.S_max,
            sigma2,
            view_mu.forest.leaf_prior_cov_inv,
            sigma_gamma2_new,
            K_chol,
            state.random_intercept,
            state.sample_beta,
            conditional_precision=conditional_precision,
            temperature=temp,
            X_unified=state.X,
            y_vec=state.y,
        )
        view_mu = eqx.tree_at(lambda v: v.forest, view_mu, view_mu_forest)

        w_transfer = b_z * beta_new[state.exposure_idx]
        view_mu_forest, mu_fit, forest_nu, nu_fit, R = inter_ensemble_transfer_step(
            random.fold_in(key, 889),
            view_mu.forest,
            mu_fit,
            forest_nu,
            nu_fit,
            R,
            alpha_new,
            w_transfer,
            obs_mask,
            state.z_vec,
            state.time_idx,
            state.T_periods,
            sigma2,
            view_mu.forest.leaf_prior_cov_inv,
            state.leaf_prior_cov_inv_nu,
            conditional_precision=conditional_precision,
            temperature=temp,
            unit_idx=state.unit_idx,
            N_units=state.N_units,
            X_unified=state.X,
            y_vec=state.y,
        )
        view_mu = eqx.tree_at(lambda v: v.forest, view_mu, view_mu_forest)

    # ------------------------------------------------------------------
    # 8b. Smooth trend block, drawn jointly with gamma and beta
    # ------------------------------------------------------------------
    # The forests are held fixed here, and conditional on them the model is
    # linear-Gaussian in (c, gamma, beta), so this draw is from their exact
    # joint conditional: acceptance one, no tuning. Drawing the three together
    # rather than in sequence is what removes the ridge between the baseline
    # trend, the unit levels and the exposure profile -- a coordinate-wise
    # sweep through blocks this dependent crawls along it. See
    # ``longbet._trend_basis``.
    trend_coef_new = state.trend_coef
    trend_local2_new = state.trend_local2
    trend_local_aux_new = state.trend_local_aux
    trend_global2_new = state.trend_global2
    trend_global_aux_new = state.trend_global_aux
    if state.use_trend_block:
        temp_tb = (
            jnp.float32(1.0)
            if state.temperature is None
            else jnp.asarray(state.temperature, dtype=jnp.float32)
        )
        # Under the horseshoe the interaction columns' prior scale is the
        # current tau * lambda_j; under the ridge it is the fixed vector. Either
        # way the block below is the same exact Gaussian draw -- the horseshoe
        # is conjugate precisely so that it does not have to change.
        if state.use_trend_horseshoe:
            col_scale = horseshoe_col_scales(
                jnp.asarray(state.trend_prior_scale_vec, dtype=jnp.float32),
                state.trend_num_period,
                state.trend_local2,
                state.trend_global2,
            )
        else:
            col_scale = jnp.asarray(
                state.trend_prior_scale_vec, dtype=jnp.float32
            )
        d_trt = b_z * (nu_fit + lin_fit if use_lin else nu_fit)
        if cut:
            # Module 1: (c, gamma) from the untreated cells alone. The
            # trajectory is not part of this draw -- it is a module-2
            # parameter and its cells are not in base_mask -- but the treated
            # cells' residuals still have to absorb the new baseline, hence
            # resid_mask=obs_mask. beta is then redrawn from the treated cells
            # given the new baseline, which is its exact module-2 conditional.
            trend_coef_new, gamma_new, _, R = trend_block_step(
                random.fold_in(key, 9311),
                state.trend_coef,
                gamma_new,
                beta_new,
                R,
                state.trend_design,
                d_trt,
                base_mask,
                state.unit_idx,
                state.N_units,
                state.exposure_idx,
                sigma2,
                sigma_gamma2_new,
                K_chol,
                col_scale,
                state.random_intercept,
                False,
                conditional_precision=None,
                temperature=temp_tb,
                trend_gram=state.trend_gram,
                trend_unit_sum=state.trend_unit_sum,
                resid_mask=obs_mask,
            )
            delta_prev, lin_prev = delta_new, lin_fit
            beta_new, nu_fit, forest_nu, R, delta_new, lin_fit = _draw_trajectory(
                random.fold_in(key, 9314), R, beta_new, nu_fit, forest_nu, b_z,
                state, sigma2, None, K_chol=K_chol,
                lin_fit=lin_fit if use_lin else None,
                delta=delta_new if use_lin else None,
            )
            if not use_lin:
                delta_new, lin_fit = delta_prev, lin_prev
        else:
            # The fixed extra columns ride along in the same conjugate draw:
            # appended to the design with their own prior scales, their
            # coefficients come out of the same solve as (c, gamma, beta).
            q_tot = state.trend_design.shape[1]
            if state.num_extra_cols:
                H_blk = jnp.concatenate([state.trend_design, state.extra_design], axis=1)
                col_scale = jnp.concatenate(
                    [col_scale, jnp.asarray(state.extra_prior_scale_vec, jnp.float32)])
                coef_blk = jnp.concatenate([state.trend_coef, state.extra_coef])
            else:
                H_blk, coef_blk = state.trend_design, state.trend_coef
            if lag:
                d_trt = nu_fit + lin_fit if use_lin else nu_fit
            coef_blk, gamma_new, beta_new, R = trend_block_step(
                random.fold_in(key, 9311),
                coef_blk,
                gamma_new,
                beta_new,
                R,
                H_blk,
                d_trt,
                obs_mask,
                state.unit_idx,
                state.N_units,
                state.exposure_idx,
                sigma2,
                sigma_gamma2_new,
                L_kappa if lag else K_chol,
                col_scale,
                state.random_intercept,
                state.sample_beta,
                conditional_precision=conditional_precision,
                temperature=temp_tb,
                trend_gram=state.trend_gram,
                trend_unit_sum=state.trend_unit_sum,
                lag_design=state.lag_design if lag else None,
                D_lag_K_chol=state.D_lag_K_chol if lag else None,
                lag_scale=jnp.sqrt(lag_scale2_cur) if (lag and state.lag_learn_scale) else None,
            )
            if state.num_extra_cols:
                trend_coef_new, extra_coef_new = coef_blk[:q_tot], coef_blk[q_tot:]
            else:
                trend_coef_new = coef_blk
        if state.use_trend_horseshoe:
            # Scales given the coefficients just drawn. Ordering matters only
            # in that both directions are valid Gibbs sweeps; drawing the
            # scales last means the state carries scales consistent with the
            # coefficients it also carries.
            (
                trend_local2_new,
                trend_local_aux_new,
                trend_global2_new,
                trend_global_aux_new,
            ) = sample_horseshoe(
                random.fold_in(key, 9313),
                trend_coef_new[state.trend_num_period:],
                state.trend_local2,
                state.trend_local_aux,
                state.trend_global2,
                state.trend_global_aux,
                state.trend_hs_global_scale,
            )
        if state.random_intercept:
            # gamma moved, so its variance is redrawn from the conditional that
            # matches the values now in the state.
            sigma_gamma2_new = _sample_inv_gamma(
                random.fold_in(key, 9312),
                jnp.float32(state.gamma_prior_a + state.N_units / 2.0),
                state.gamma_prior_b + 0.5 * jnp.sum(jnp.square(gamma_new)),
            )

    if state.sample_beta and not state.anchor_moderator:
        sign_beta_end = (_lag_sign(beta_new) if lag
                         else jnp.where(jnp.sum(beta_new[1:]) < 0.0, -1.0, 1.0))
        beta_new = sign_beta_end * beta_new
        nu_fit = sign_beta_end * nu_fit
        forest_nu = eqx.tree_at(
            lambda f: f.leaf_tree, forest_nu, sign_beta_end * forest_nu.leaf_tree
        )
        if use_lin:
            delta_new = sign_beta_end * delta_new
            lin_fit = sign_beta_end * lin_fit

    # 8c. Trajectory scale s^2 given beta (prior-only, so never tempered).
    beta_scale2_new = state.beta_scale2
    beta_scale_aux_new = state.beta_scale_aux
    if state.learn_beta_scale and not state.anchor_moderator and not lag:
        beta_scale2_new, beta_scale_aux_new = _half_cauchy_scale_step(
            random.fold_in(key, 7301),
            beta_prior_quadform(beta_new, state.K_chol),
            float(state.S_max + 1),
            state.beta_scale2, state.beta_scale_aux, state.beta_scale_prior,
        )
    # 8c'. Lag-kernel scale s^2 given kappa (prior-only; anchored or not).
    lag_scale2_new, lag_scale_aux_new = state.lag_scale2, state.lag_scale_aux
    if lag and state.lag_learn_scale and state.sample_beta:
        lag_scale2_new, lag_scale_aux_new = _half_cauchy_scale_step(
            random.fold_in(key, 7302),
            beta_prior_quadform(beta_new, state.lag_K_chol_sq),
            float(beta_new.shape[0]),
            state.lag_scale2, state.lag_scale_aux, state.lag_scale_prior,
        )

    # ------------------------------------------------------------------
    # 8d. Interweaving update of sigma_gamma (encourage-v2, interweave_gamma)
    # ------------------------------------------------------------------
    interweave_accepts_new = state.interweave_accepts
    if state.interweave_gamma:
        e_g = jnp.where(obs_mask, R + gamma_new[state.unit_idx], 0.0)
        if conditional_precision is None:
            p_c = jnp.where(obs_mask, jnp.reciprocal(sigma2), 0.0)
        else:
            p_c = jnp.where(obs_mask, conditional_precision, 0.0)
        P_u = jnp.zeros(state.N_units, jnp.float32).at[state.unit_idx].add(p_c)
        h_u = jnp.zeros(state.N_units, jnp.float32).at[state.unit_idx].add(p_c * e_g)
        gamma_iw, sigma_gamma2_new, iw_ok = _interweave_gamma_scale(
            random.fold_in(key, 7601), gamma_new, sigma_gamma2_new, P_u, h_u,
            state.gamma_prior_a, state.gamma_prior_b)
        R = jnp.where(obs_mask, e_g - gamma_iw[state.unit_idx], 0.0)
        gamma_new = gamma_iw
        interweave_accepts_new = interweave_accepts_new + iw_ok.astype(jnp.float32)

    # ------------------------------------------------------------------
    # 9. Error variance sigma^2
    # ------------------------------------------------------------------
    if conditional_precision is not None:
        sigma2_new = sigma2
    elif is_binary or is_ordered:
        sigma2_new = jnp.float32(1.0)
    else:
        # Module-1 parameter under the cut: the noise level is learned from
        # the untreated cells, as every other baseline quantity is.
        m_obs = jnp.sum(base_mask.astype(jnp.float32))
        sigma2_new = _sample_inv_gamma(
            keys[9],
            state.sigma_prior_a + m_obs / 2.0,
            state.sigma_prior_b + 0.5 * jnp.sum(jnp.square(jnp.where(base_mask, R, 0.0))),
        )

    # ------------------------------------------------------------------
    # Assemble the updated state
    # ------------------------------------------------------------------
    # error_cov_inv records the *sampled* error precision so the saved trace
    # carries a meaningful value; the sweep itself always reads state.sigma2.
    error_cov_inv = eqx.tree_at(
        lambda w_: w_.value,
        state.error_cov_inv,
        jnp.reciprocal(sigma2_new).astype(jnp.float32),
    )

    extra_names, extra_values = ["ridge_accepts"], [ridge_accepts_new]
    if state.interweave_gamma:
        extra_names.append("interweave_accepts"); extra_values.append(interweave_accepts_new)
    if use_lin:
        extra_names += ["delta", "lin_fit", "tau_delta2", "tau_delta_aux", "lin_local2",
                        "lin_local_aux"]
        extra_values += [delta_new, lin_fit, tau_delta2_new, tau_delta_aux_new, lin_local2_new,
                         lin_local_aux_new]
    if state.learn_beta_scale:
        extra_names += ["beta_scale2", "beta_scale_aux"]
        extra_values += [beta_scale2_new, beta_scale_aux_new]
    if lag:
        # ``beta_new`` carried the kernel through the sweep; the placeholder
        # ``beta`` is left as it was.
        extra_names.append("kappa"); extra_values.append(beta_new)
        beta_new = state.beta
        if state.lag_learn_scale:
            extra_names += ["lag_scale2", "lag_scale_aux"]
            extra_values += [lag_scale2_new, lag_scale_aux_new]
    if state.num_extra_cols:
        extra_names.append("extra_coef"); extra_values.append(extra_coef_new)
    if state.ridge_accepts is None:
        extra_names, extra_values = extra_names[1:], extra_values[1:]
    if extra_names:
        state = eqx.tree_at(
            lambda s: tuple(getattr(s, n) for n in extra_names), state, tuple(extra_values))

    return eqx.tree_at(
        lambda s: (
            s.X,
            s.y,
            s.config,
            s.z,
            s.cutpoints,
            s.resid,
            s.forest,
            s.resid_unit,
            s.resid_eff_scale,
            s.resid_inexact_integral,
            s.error_cov_inv,
            s.prec_scale,
            s.inv_sdev_scale,
            s.inv_sdev_unit,
            s.n_non_missing,
            s.sum_diag_prec_scale,
            s.forest_nu,
            s.resid_nu,
            s.resid_unit_nu,
            s.resid_eff_scale_nu,
            s.resid_inexact_integral_nu,
            s.prec_scale_nu,
            s.inv_sdev_scale_nu,
            s.inv_sdev_unit_nu,
            s.leaf_prior_cov_inv_nu,
            s.alpha,
            s.beta,
            s.b0,
            s.b1,
            s.gamma,
            s.sigma2,
            s.sigma_gamma2,
            s.mu_fit,
            s.nu_fit,
            s.trend_coef,
            s.trend_local2,
            s.trend_local_aux,
            s.trend_global2,
            s.trend_global_aux,
        ),
        state,
        (
            view_nu.X,
            view_nu.y,
            view_nu.config,
            current_z,
            cutpoints,
            R,
            view_mu.forest,
            view_mu.resid_unit,
            view_mu.resid_eff_scale,
            view_mu.resid_inexact_integral,
            error_cov_inv,
            state.prec_scale if conditional_attrs_mu is not None else view_mu.prec_scale,
            state.inv_sdev_scale if conditional_attrs_mu is not None else view_mu.inv_sdev_scale,
            state.inv_sdev_unit if conditional_attrs_mu is not None else view_mu.inv_sdev_unit,
            view_mu.n_non_missing,
            state.sum_diag_prec_scale if conditional_attrs_mu is not None else view_mu.sum_diag_prec_scale,
            forest_nu,
            view_nu.resid,
            view_nu.resid_unit,
            view_nu.resid_eff_scale,
            view_nu.resid_inexact_integral,
            view_nu.prec_scale,
            view_nu.inv_sdev_scale,
            view_nu.inv_sdev_unit,
            forest_nu.leaf_prior_cov_inv,
            alpha_new,
            beta_new,
            b0_new,
            b1_new,
            gamma_new,
            sigma2_new,
            sigma_gamma2_new,
            mu_fit,
            nu_fit,
            trend_coef_new,
            trend_local2_new,
            trend_local_aux_new,
            trend_global2_new,
            trend_global_aux_new,
        ),
        is_leaf=lambda x: x is None,
    )


@jax.jit
@enable_x64()
def longbet_step(key: Key[Array, ''], state: LongBetState) -> LongBetState:
    """Execute one Gibbs sweep, across chains if the state carries a chain axis.

    Chains are mapped with ``jax.vmap`` over the per-chain partition of the
    state only; ``X``, ``y``, the panel indices and the GP kernel factor are
    closed over, so a ``k``-chain fit does not hold ``k`` copies of the design
    matrix.  ``bartz``'s own ``num_chains`` cannot be used here because its
    weighted-observation path does not accept a per-chain ``prec_scale``, and
    LongBet's treatment weights differ across chains by construction.
    """
    sweep = _tempered_single_step if state.is_tempered else longbet_single_step
    if not state.has_chain_axis:
        return sweep(key, state)

    spec = chain_filter_spec(state)
    per_chain, shared = eqx.partition(state, spec)
    keys = random.split(key, state.num_chains)

    def one(kp: tuple[Key[Array, ''], Any]) -> Any:
        k, p = kp
        out = sweep(k, eqx.combine(p, shared))
        # Return only the per-chain part, so map does not add a chain axis to
        # the shared arrays.
        return eqx.filter(out, spec)

    out_per = lax.map(one, (keys, per_chain))
    return eqx.combine(out_per, shared)
