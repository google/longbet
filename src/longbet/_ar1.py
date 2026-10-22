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

"""Stationary AR(1) innovations within a unit's series (``ar1_errors``).

The model is ``eps_it = phi eps_{i,t-1} + xi_it`` with ``xi ~ N(0, sigma^2)``
i.i.d. and a stationary first period, ``eps_i1 ~ N(0, sigma^2 / (1 - phi^2))``.
The Prais-Winsten transform of each unit's series,

.. code::

    (T v)_{i1} = sqrt(1 - phi^2) v_{i1},    (T v)_{it} = v_{it} - phi v_{i,t-1},

maps the errors to the i.i.d. innovations, so the log-likelihood is
``-(m/2) log sigma^2 + (N/2) log(1 - phi^2) - |T R|^2 / (2 sigma^2)`` with ``R``
the untransformed residual. Every Gibbs block that is linear in its parameters
keeps its exact Gaussian conditional in the transformed space: the sensitivity
rows are transformed with ``T`` and so is the residual.

For a unit-level component (the prognostic forest with no calendar splits,
the unit intercept) the transformed contribution of the unit's value ``a_i`` is
``a_i g`` with ``g = (sqrt(1 - phi^2), 1 - phi, ..., 1 - phi)`` on a balanced
panel; for the treatment term ``w_it nu(x_i)`` it is ``nu(x_i) (T w)_it``. So
both forests keep the partially collapsed update over ``gamma_i``, with
``g``-weighted unit sums in place of plain ones (see ``longbet._step``).

Scope: continuous outcome, balanced panel, unit-level partitions (the
production configuration). ``phi`` is drawn by an independence Metropolis step
whose proposal is a Laplace fit of its exact conditional under a uniform prior
on ``(-phi_max, phi_max)``.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import jax.scipy.linalg as jla
from jax import random
from jaxtyping import Array, Bool, Float32, Key

from longbet._x64 import enable_x64


def transform(v: Array, phi: Array, first: Bool[Array, ' n']) -> Array:
    """``T_phi`` along the cell axis (axis 0) of a unit-major balanced panel.

    ``first`` marks each unit's first period. ``v`` is ``(n,)`` or ``(n, k)``;
    masked cells must already be zero.
    """
    lag = jnp.roll(v, 1, axis=0)
    phi = jnp.asarray(phi, dtype=v.dtype)
    root = jnp.sqrt(jnp.maximum(1.0 - jnp.square(phi), 0.0)).astype(v.dtype)
    f = first if v.ndim == 1 else first[:, None]
    return jnp.where(f, root * v, v - phi * lag)


def unit_weights(phi: Array, first: Bool[Array, ' n']) -> Array:
    """``g``: the transformed coefficient of a unit-level constant, per cell."""
    return transform(jnp.ones(first.shape, jnp.float32), phi, first)


def unit_sum(values: Array, unit_idx: Array, num_units: int) -> Array:
    """Sum over a unit's cells (axis 0)."""
    shape = (num_units,) + values.shape[1:]
    return jnp.zeros(shape, values.dtype).at[unit_idx].add(values)


def sample_beta_gp_full(
    key: Key[Array, ''],
    A: Float32[Array, 'd d'],
    C: Float32[Array, ' d'],
    K_chol: Float32[Array, 'd d'],
) -> Float32[Array, ' d']:
    """Draw ``beta ~ N(P^-1 C, P^-1)``, ``P = K^-1 + A``, for a full data precision ``A``.

    The whitened form of ``longbet._gp.sample_beta_gp``: ``beta = L u`` with
    ``u ~ N(M^-1 L' C, M^-1)``, ``M = I + L' A L``, computed in float64.
    """
    with enable_x64():
        L = jnp.asarray(K_chol, jnp.float64)
        A64 = jnp.asarray(A, jnp.float64)
        d = L.shape[-1]
        M = jnp.eye(d, dtype=jnp.float64) + L.T @ A64 @ L
        M = 0.5 * (M + M.T)
        LM = jnp.linalg.cholesky(M)
        mean_u = jla.cho_solve((LM, True), L.T @ jnp.asarray(C, jnp.float64))
        eta = random.normal(key, (d,), jnp.float32).astype(jnp.float64)
        u = mean_u + jla.solve_triangular(LM.T, eta, lower=False)
        return (L @ u).astype(jnp.float32)


def residual_stats(R: Array, first: Bool[Array, ' n'], obs_mask: Bool[Array, ' n']):
    """Sufficient statistics of ``phi``: ``(A0, B, C, D, N)`` with ``A0 = sum eps_{i1}^2``,
    ``B = sum eps_t eps_{t-1}``, ``C = sum eps_{t-1}^2``, ``D = sum eps_t^2`` over
    ``t >= 2``, and ``N`` the number of units."""
    with enable_x64():
        e = jnp.where(obs_mask, jnp.asarray(R, jnp.float64), 0.0)
        lag = jnp.roll(e, 1)
        later = ~first
        A0 = jnp.sum(jnp.where(first, e * e, 0.0))
        B = jnp.sum(jnp.where(later, e * lag, 0.0))
        C = jnp.sum(jnp.where(later, lag * lag, 0.0))
        D = jnp.sum(jnp.where(later, e * e, 0.0))
        n_units = jnp.sum(first.astype(jnp.float64))
        return A0, B, C, D, n_units


def phi_log_conditional(phi, stats, lam, temperature):
    """Tempered log conditional of ``phi`` (uniform prior), up to a constant.

    ``temperature * (N/2) log(1 - phi^2) - lam/2 * SSR(phi)`` with ``lam`` the
    innovation precision times the temperature and
    ``SSR(phi) = (1 - phi^2) A0 + D - 2 phi B + phi^2 C``.
    """
    A0, B, C, D, N = stats
    ssr = (1.0 - phi * phi) * A0 + D - 2.0 * phi * B + phi * phi * C
    return 0.5 * temperature * N * jnp.log1p(-phi * phi) - 0.5 * lam * ssr


@enable_x64()
def sample_phi(
    key: Key[Array, ''],
    phi: Float32[Array, ''],
    R: Float32[Array, ' n'],
    first: Bool[Array, ' n'],
    obs_mask: Bool[Array, ' n'],
    lam: Float32[Array, ''],
    temperature: Float32[Array, ''],
    phi_max: float,
) -> tuple[Float32[Array, ''], Bool[Array, '']]:
    """Independence Metropolis draw of ``phi`` from its exact conditional.

    The proposal is ``N(m, (1.25 / sqrt(h))^2)`` with ``m`` the conditional's
    mode (six Newton steps from the least-squares estimate, clipped to the
    support) and ``h`` its negative curvature there. The log conditional is
    concave on the support, so the fit is close and the acceptance high; the
    ratio corrects exactly for the difference. Out-of-support proposals are
    rejected.
    """
    stats = residual_stats(R, first, obs_mask)
    A0, B, C, D, N = stats
    lam = jnp.asarray(lam, jnp.float64)
    temp = jnp.asarray(temperature, jnp.float64)
    lim = float(phi_max)

    def grad_hess(p):
        g = (-temp * N * p / (1.0 - p * p)
             - 0.5 * lam * (-2.0 * p * A0 - 2.0 * B + 2.0 * p * C))
        h = (-temp * N * (1.0 + p * p) / jnp.square(1.0 - p * p)
             - 0.5 * lam * (-2.0 * A0 + 2.0 * C))
        return g, h

    m = jnp.clip(B / jnp.maximum(C, 1e-12), -lim, lim)
    for _ in range(6):
        g, h = grad_hess(m)
        m = jnp.clip(m - g / jnp.minimum(h, -1e-12), -lim + 1e-6, lim - 1e-6)
    _, h = grad_hess(m)
    sd = 1.25 / jnp.sqrt(jnp.maximum(-h, 1e-12))
    k_prop, k_acc = random.split(key)
    prop = m + sd * random.normal(k_prop, (), jnp.float32).astype(jnp.float64)
    cur = jnp.asarray(phi, jnp.float64)

    def log_q(x):
        return -0.5 * jnp.square((x - m) / sd)

    in_support = jnp.abs(prop) < lim
    safe_prop = jnp.where(in_support, prop, 0.0)
    log_ratio = (phi_log_conditional(safe_prop, stats, lam, temp) + log_q(cur)
                 - phi_log_conditional(cur, stats, lam, temp) - log_q(safe_prop))
    u = random.uniform(k_acc, (), jnp.float32).astype(jnp.float64)
    accept = in_support & (jnp.log(u + 1e-300) < log_ratio)
    return jnp.where(accept, prop, cur).astype(jnp.float32), accept


@enable_x64()
def trend_ar1_statistics(H, obs_mask, first, unit_idx, N_units):
    """Fixed pieces of the transformed trend design's cross-products.

    With ``Ht = T_phi H`` on a balanced unit-major panel,

    ``Ht' Ht = G_later - phi (C + C') + phi^2 G_lag + (1 - phi^2) G_first``

    (``G_first``: first-period rows, ``G_later``: later rows, ``G_lag``: their
    lagged rows, ``C = sum_{t >= 2} H_t H_{t-1}'``), and a unit's ``g``-weighted
    sum of transformed rows is
    ``(1 - phi^2) U_first + (1 - phi)(U_later - phi U_lag)``. Computed once, in
    float64, so the per-sweep cost of the trend block no longer scales with
    ``n q^2``.
    """
    H = jnp.where(obs_mask[:, None], jnp.asarray(H, jnp.float64), 0.0)
    lag = jnp.roll(H, 1, axis=0)
    f = first[:, None]
    Hf = jnp.where(f, H, 0.0)
    Hl = jnp.where(f, 0.0, H)
    Lg = jnp.where(f, 0.0, lag)
    return {
        "G_first": Hf.T @ Hf, "G_later": Hl.T @ Hl, "G_lag": Lg.T @ Lg, "C": Hl.T @ Lg,
        "U_first": unit_sum(Hf, unit_idx, N_units), "U_later": unit_sum(Hl, unit_idx, N_units),
        "U_lag": unit_sum(Lg, unit_idx, N_units),
    }


@enable_x64()
def trend_block_step_ar1(
    key: Key[Array, ''],
    trend_coef: Float32[Array, ' q_tot'],
    gamma: Float32[Array, ' N_units'],
    beta: Float32[Array, ' S_max_plus_1'],
    R: Float32[Array, ' n'],
    H: Float32[Array, 'n q_tot'],
    d_vec: Float32[Array, ' n'],
    obs_mask: Bool[Array, ' n'],
    unit_idx,
    N_units: int,
    exposure_idx,
    lam: Float32[Array, ''],
    sigma_gamma2: Float32[Array, ''],
    K_chol: Float32[Array, 'S_max_plus_1 S_max_plus_1'],
    trend_col_scale: Float32[Array, ' q_tot'],
    phi: Float32[Array, ''],
    first: Bool[Array, ' n'],
    sample_beta: bool = True,
    stats: dict | None = None,
    lag_design: Float32[Array, 'n L_plus_1'] | None = None,
    D_lag_K_chol: Array | None = None,
    lag_scale: Array | float | None = None,
):
    """Exact joint draw of ``(c, gamma, beta)`` under AR(1) innovations.

    ``longbet._trend_basis.trend_block_step`` with the design rows, the
    trajectory rows ``d_it L[S_it]`` and the residual transformed by ``T_phi``,
    and with ``gamma_i`` entering transformed cell ``(i, t)`` with weight
    ``g_t``. ``lam`` is the (tempered) innovation precision. The returned
    residual is untransformed.

    With ``lag_design`` (the distributed-lag kernel) ``beta`` is the kernel
    ``kappa``, ``K_chol`` its factor ``L_kappa`` of shape ``(L + 1, r)``, and
    the trajectory rows are ``d_it (lag_design @ L_kappa)[it]`` (or
    ``(lag_scale * d_it) * D_lag_K_chol[it]`` when ``D_lag_K_chol = lag_design
    @ K_chol`` is precomputed) -- the dense analogue of the one-hot lookup,
    whitened the same way; ``exposure_idx`` is then unused.
    """
    f64 = jnp.float64
    H = jnp.where(obs_mask[:, None], jnp.asarray(H, f64), 0.0)
    R = jnp.asarray(R, f64)
    d_vec = jnp.asarray(d_vec, f64)
    trend_coef = jnp.asarray(trend_coef, f64)
    gamma = jnp.asarray(gamma, f64)
    beta = jnp.asarray(beta, f64)
    lam = jnp.asarray(lam, f64)
    phi64 = jnp.asarray(phi, f64)
    L = jnp.asarray(K_chol, f64)
    S1 = L.shape[-1]
    q_tot = H.shape[1]
    s_col = jnp.asarray(trend_col_scale, f64)
    lag64 = None if lag_design is None else jnp.asarray(lag_design, f64)

    def traj_fit(coef):
        return d_vec * (coef[exposure_idx] if lag64 is None else lag64 @ coef)

    R_tilde = R + H @ trend_coef + gamma[unit_idx]
    if sample_beta:
        R_tilde = R_tilde + traj_fit(beta)
    R_tilde = jnp.where(obs_mask, R_tilde, 0.0)

    Rt = transform(R_tilde, phi64, first)
    Hs_t = transform(H * s_col[None, :], phi64, first)
    g = transform(jnp.ones(first.shape, f64), phi64, first)

    if stats is not None:
        gram = (stats["G_later"] - phi64 * (stats["C"] + stats["C"].T)
                + phi64 * phi64 * stats["G_lag"] + (1.0 - phi64 * phi64) * stats["G_first"])
        P_cc = lam * (s_col[:, None] * gram * s_col[None, :]) + jnp.eye(q_tot, dtype=f64)
    else:
        P_cc = lam * (Hs_t.T @ Hs_t) + jnp.eye(q_tot, dtype=f64)
    h_c = lam * (Hs_t.T @ Rt)
    if sample_beta:
        if D_lag_K_chol is not None:
            scale64 = jnp.float64(1.0) if lag_scale is None else jnp.asarray(lag_scale, f64)
            Dt = jnp.where(
                obs_mask[:, None],
                (scale64 * d_vec)[:, None] * jnp.asarray(D_lag_K_chol, f64),
                0.0,
            )
        else:
            rows_w = L[exposure_idx] if lag64 is None else lag64 @ L
            Dt = jnp.where(obs_mask[:, None], d_vec[:, None] * rows_w, 0.0)
        Dt_t = transform(Dt, phi64, first)
        P_uu = lam * (Dt_t.T @ Dt_t) + jnp.eye(S1, dtype=f64)
        P_cu = lam * (Hs_t.T @ Dt_t)
        h_u = lam * (Dt_t.T @ Rt)
        P_ww = jnp.block([[P_cc, P_cu], [P_cu.T, P_uu]])
        h_w = jnp.concatenate([h_c, h_u])
        rows = None if stats is not None else jnp.concatenate([Hs_t, Dt_t], axis=1)
    else:
        P_ww, h_w, rows = P_cc, h_c, Hs_t

    p_gg = unit_sum(lam * g * g, unit_idx, N_units) + jnp.reciprocal(jnp.asarray(sigma_gamma2, f64))
    h_g = unit_sum(lam * g * Rt, unit_idx, N_units)
    if stats is not None:
        U = ((1.0 - phi64 * phi64) * stats["U_first"]
             + (1.0 - phi64) * (stats["U_later"] - phi64 * stats["U_lag"]))          # (N, q)
        P_cg = lam * (U * s_col[None, :]).T
        if sample_beta:
            P_ug = unit_sum(lam * g[:, None] * Dt_t, unit_idx, N_units).T
            P_wg = jnp.concatenate([P_cg, P_ug], axis=0)
        else:
            P_wg = P_cg
    else:
        P_wg = unit_sum(lam * g[:, None] * rows, unit_idx, N_units).T          # (k, N)
    scaled = P_wg / p_gg[None, :]
    P_eff = P_ww - scaled @ P_wg.T
    h_eff = h_w - scaled @ h_g

    P_eff = 0.5 * (P_eff + P_eff.T)
    L_eff = jnp.linalg.cholesky(P_eff)
    eta = random.normal(key, h_eff.shape, jnp.float32).astype(f64)
    w_new = jla.cho_solve((L_eff, True), h_eff) + jla.solve_triangular(L_eff.T, eta, lower=False)

    trend_new = (s_col * w_new[:q_tot]).astype(jnp.float32)
    beta_new = (L @ w_new[q_tot:]).astype(jnp.float32) if sample_beta else beta.astype(jnp.float32)
    v_gg = jnp.reciprocal(p_gg)
    mean_g = v_gg * (h_g - P_wg.T @ w_new)
    gamma_new = (mean_g + random.normal(random.fold_in(key, 1), (N_units,), jnp.float32)
                 * jnp.sqrt(v_gg)).astype(jnp.float32)

    fitted = H @ jnp.asarray(trend_new, f64) + jnp.asarray(gamma_new, f64)[unit_idx]
    if sample_beta:
        fitted = fitted + traj_fit(jnp.asarray(beta_new, f64))
    R_new = jnp.where(obs_mask, R_tilde - fitted, 0.0)
    return trend_new, gamma_new, beta_new, R_new.astype(jnp.float32)
