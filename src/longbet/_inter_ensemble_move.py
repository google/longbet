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

"""Exact two-stage conjugate subspace Gibbs inter-ensemble residual-transfer step.

On staggered panels where treated and control cohorts follow diverging baseline
trends, the prognostic ensemble ``mu(X_i)``, unit intercepts ``gamma_i``, and
treatment component ``beta_{S_it} * nu(X_i, t)`` can trade off post-adoption
level, calendar-time slope/curvature, and covariate-moderation profiles on
treated cells (Z_it = 1). Sequential Gibbs updates update one component holding
the others fixed, which can slow mixing along their shared affine subspace.

This module implements an exact two-stage conjugate Gaussian subspace Gibbs step:
  - Stage 1 (8-D subspace ``delta in R^8``): jointly shifts the prognostic leaf
    values ``mu_l``, unit intercepts ``gamma_i``, and treatment exposure curve
    ``beta_s`` along unit-average polynomials ``(1, bar_t_i, bar_q_i, bar_w_i)``
    and within-unit demeaned profiles ``(t_it - bar_t_i, q_it - bar_q_i,
    w_it - bar_w_i)``, weighted on prognostic leaves by squared treated fraction
    ``p_trt^2`` so pure control leaves remain untouched.
  - Stage 2 (6-D subspace ``delta_het in R^6``): jointly shifts the prognostic
    leaves ``mu_l`` and treatment leaves ``nu_l`` along the tensor product of
    time polynomials ``(1, t_norm, q_norm)`` and standardized prognostic/treatment
    covariate scores.

Because both the Gaussian likelihood (or probit/ordinal latent Gaussian
likelihood) and the leaf/intercept/GP priors are quadratic in the shifted
parameters, each stage's conditional posterior distribution is an exact
low-dimensional Gaussian distribution. Drawing from this conditional has 100%
acceptance probability, requires no step-size tuning, and leaves the target
joint posterior distribution strictly invariant.
"""

from __future__ import annotations

from typing import Any

import equinox as eqx
import jax
import jax.numpy as jnp
from jaxtyping import Array, Bool, Float32, Int32, Key

from bartz.mcmcstep._step import apply_moves_to_leaf_indices
from longbet._linalg import sample_from_precision
from longbet._ridge import compute_active_leaf_stats


def _leaf_stats_vec(
    ids: Int32[Array, ' n'],
    Phi: Float32[Array, 'n K'],
    mask: Bool[Array, ' n'],
    max_slots: int,
) -> tuple[Float32[Array, 'max_slots K'], Float32[Array, ' max_slots']]:
    """Compute per-slot mean of ``Phi`` (shape ``(n, K)``) and active count."""
    m_f32 = mask.astype(jnp.float32)
    Phi_masked = jnp.where(mask[:, None], Phi, 0.0).astype(jnp.float32)
    sum_Phi = jnp.zeros((max_slots, Phi.shape[1]), dtype=jnp.float32).at[ids].add(Phi_masked)
    cnt = jnp.zeros(max_slots, dtype=jnp.float32).at[ids].add(m_f32)
    mean_Phi = jnp.where(cnt[:, None] > 0, sum_Phi / jnp.maximum(cnt[:, None], 1.0), 0.0)
    return mean_Phi, cnt


def _compute_unit_prognostic_score(
    X_unified: Any,
    y_vec: Float32[Array, ' n'] | None,
    z_vec: Float32[Array, ' n'],
    obs_mask: Bool[Array, ' n'],
    unit_idx: Int32[Array, ' n'] | None,
    N_units: int | None,
) -> Float32[Array, ' n']:
    """Compute standardized unit-level pre-treatment prognostic score g_i projected onto X_i."""
    if X_unified is None or y_vec is None or unit_idx is None or N_units is None:
        return jnp.zeros_like(z_vec)
    ctrl_obs = obs_mask & (z_vec == 0.0)
    ctrl_f32 = ctrl_obs.astype(jnp.float32)
    sum_y0 = (
        jnp.zeros(N_units, dtype=jnp.float32)
        .at[unit_idx]
        .add(jnp.where(ctrl_obs, y_vec, 0.0))
    )
    cnt_y0 = jnp.zeros(N_units, dtype=jnp.float32).at[unit_idx].add(ctrl_f32)
    y0_mean = jnp.where(cnt_y0 > 0, sum_y0 / jnp.maximum(cnt_y0, 1.0), 0.0)

    p_take = min(int(X_unified.shape[0]), 12)
    X_cell = X_unified[:p_take, :].T.astype(jnp.float32)
    sum_X = (
        jnp.zeros((N_units, p_take), dtype=jnp.float32)
        .at[unit_idx]
        .add(X_cell)
    )
    cnt_X = (
        jnp.zeros(N_units, dtype=jnp.float32)
        .at[unit_idx]
        .add(jnp.ones(X_cell.shape[0], dtype=jnp.float32))
    )
    X_unit = sum_X / jnp.maximum(cnt_X[:, None], 1.0)
    X_std = (X_unit - jnp.mean(X_unit, axis=0, keepdims=True)) / (
        jnp.std(X_unit, axis=0, keepdims=True) + 1e-4
    )
    idx_i, idx_j = jnp.triu_indices(p_take)
    X_pairs = X_std[:, idx_i] * X_std[:, idx_j]
    X_des = jnp.concatenate(
        [jnp.ones((N_units, 1), dtype=jnp.float32), X_std, X_pairs], axis=1
    )
    d_des = X_des.shape[1]
    gram = X_des.T @ X_des + 10.0 * jnp.eye(d_des, dtype=jnp.float32)
    rhs = X_des.T @ y0_mean
    beta_ols = jax.scipy.linalg.solve(gram, rhs, assume_a='pos')
    g_hat = X_des @ beta_ols
    g_tilde = (g_hat - jnp.mean(g_hat)) / (jnp.std(g_hat) + 1e-4)
    return g_tilde[unit_idx]


def inter_ensemble_transfer_step(
    key: Key[Array, ''],
    forest_mu: Any,
    mu_fit: Float32[Array, ' n'],
    forest_nu: Any,
    nu_fit: Float32[Array, ' n'],
    R: Float32[Array, ' n'],
    alpha: Float32[Array, ''],
    w: Float32[Array, ' n'],
    obs_mask: Bool[Array, ' n'],
    z_vec: Float32[Array, ' n'],
    time_idx: Int32[Array, ' n'],
    T_periods: int,
    sigma2: Float32[Array, ''],
    leaf_prior_cov_inv_mu: Float32[Array, ''],
    leaf_prior_cov_inv_nu: Float32[Array, ''],
    conditional_precision: Float32[Array, ' n'] | None = None,
    temperature: Float32[Array, ''] = jnp.float32(1.0),
    min_abs_beta: float = 0.05,
    unit_idx: Int32[Array, ' n'] | None = None,
    N_units: int | None = None,
    X_unified: Any = None,
    y_vec: Float32[Array, ' n'] | None = None,
    ar1: Any = None,
) -> tuple[Any, Float32[Array, ' n'], Any, Float32[Array, ' n'], Float32[Array, ' n']]:
    """Execute one exact conjugate subspace Gibbs transfer step between mu and nu leaves."""
    denom_t = jnp.float32(max(1, T_periods - 1))
    t_norm = (time_idx.astype(jnp.float32) - 0.5 * denom_t) / denom_t
    t_quad = jnp.square(t_norm) - jnp.float32(1.0 / 12.0)
    w_norm = w / jnp.maximum(jnp.max(jnp.abs(w)), 1.0)
    g_cell = _compute_unit_prognostic_score(
        X_unified, y_vec, z_vec, obs_mask, unit_idx, N_units
    )

    Phi = jnp.stack(
        [
            jnp.ones_like(t_norm),
            t_norm,
            t_quad,
            w_norm,
            g_cell * t_norm,
            g_cell * w_norm,
        ],
        axis=-1,
    )  # shape (n, K=6)
    K = Phi.shape[1]

    ids_mu = apply_moves_to_leaf_indices(
        forest_mu.leaf_indices, forest_mu.to_prune, forest_mu.move_node
    )
    ids_nu = apply_moves_to_leaf_indices(
        forest_nu.leaf_indices, forest_nu.to_prune, forest_nu.move_node
    )

    _, _, is_leaf_mu = compute_active_leaf_stats(
        forest_mu.split_tree, forest_mu.leaf_tree
    )
    _, _, is_leaf_nu = compute_active_leaf_stats(
        forest_nu.split_tree, forest_nu.leaf_tree
    )

    max_slots_mu = forest_mu.leaf_tree.shape[-1]
    max_slots_nu = forest_nu.leaf_tree.shape[-1]
    m_mu = jnp.float32(forest_mu.leaf_tree.shape[0])
    m_nu = jnp.float32(forest_nu.leaf_tree.shape[0])

    # 1. Compute leaf basis directions A_mu (data units), shape (m_mu, max_slots_mu, K)
    trt_obs_mask = obs_mask & (z_vec == 1.0)
    Phi_mean_mu, cnt_all_mu = jax.vmap(
        _leaf_stats_vec, in_axes=(0, None, None, None)
    )(ids_mu, Phi, obs_mask, max_slots_mu)
    _, cnt_trt_mu = jax.vmap(
        _leaf_stats_vec, in_axes=(0, None, None, None)
    )(ids_mu, Phi[:, :1], trt_obs_mask, max_slots_mu)

    p_trt_mu = jnp.where(
        cnt_all_mu > 0, cnt_trt_mu / jnp.maximum(cnt_all_mu, 1.0), 0.0
    )
    active_mu = is_leaf_mu & (cnt_trt_mu > 0)
    weight_mu = jnp.where(active_mu, jnp.square(p_trt_mu) / m_mu, 0.0)
    A_mu = weight_mu[:, :, None] * Phi_mean_mu  # (m_mu, max_slots_mu, K)

    # Induced change J_mu in mu_fit per unit of delta: shape (n, K)
    J_mu_per_tree = jax.vmap(lambda a_tree, idx: a_tree[idx])(A_mu, ids_mu)
    J_mu = jnp.sum(J_mu_per_tree, axis=0)  # (n, K)

    # 2. Compute compensating leaf basis directions A_nu (data units), shape (m_nu, max_slots_nu, K)
    target_nu_cell = alpha * J_mu  # (n, K)
    target_mean_nu, cnt_trt_nu = jax.vmap(
        _leaf_stats_vec, in_axes=(0, None, None, None)
    )(ids_nu, target_nu_cell, trt_obs_mask, max_slots_nu)
    w_mean_nu, _ = jax.vmap(
        _leaf_stats_vec, in_axes=(0, None, None, None)
    )(ids_nu, w[:, None], trt_obs_mask, max_slots_nu)
    w_mean_nu = w_mean_nu[:, :, 0]

    active_nu = (
        is_leaf_nu
        & (cnt_trt_nu > 0)
        & (jnp.abs(w_mean_nu) > min_abs_beta)
    )
    safe_w_nu = jnp.where(active_nu, w_mean_nu, 1.0)
    A_nu = jnp.where(
        active_nu[:, :, None],
        -target_mean_nu / (m_nu * safe_w_nu[:, :, None]),
        0.0,
    )  # (m_nu, max_slots_nu, K)

    # Induced change J_nu in nu_fit per unit of delta: shape (n, K)
    J_nu_per_tree = jax.vmap(lambda a_tree, idx: a_tree[idx])(A_nu, ids_nu)
    J_nu = jnp.sum(J_nu_per_tree, axis=0)  # (n, K)

    # Total residual sensitivity H = alpha * J_mu + w * J_nu: shape (n, K)
    H = jnp.where(
        obs_mask[:, None],
        alpha * J_mu + w[:, None] * J_nu,
        0.0,
    )

    # 3. Assemble exact Gaussian conditional precision (Sigma_inv) and linear term (h)
    if conditional_precision is None:
        cell_prec = jnp.where(obs_mask, temperature / sigma2, 0.0)
    else:
        cell_prec = jnp.where(obs_mask, conditional_precision, 0.0)

    if ar1 is None:
        H_lik, R_lik = H, R
    else:
        # AR(1) innovations: the likelihood is Gaussian in the Prais-Winsten
        # transformed sensitivities and residual (longbet._ar1); the update
        # itself is applied to the untransformed parameters below.
        from longbet._ar1 import transform as _ar1_transform
        phi_ar1, first_ar1 = ar1
        H_lik = _ar1_transform(H, phi_ar1, first_ar1)
        R_lik = _ar1_transform(jnp.where(obs_mask, R, 0.0), phi_ar1, first_ar1)
    H_weighted = H_lik * cell_prec[:, None]
    prec_lik = H_lik.T @ H_weighted  # (K, K)
    h_lik = H_weighted.T @ R_lik  # (K,)

    # Prior contribution from mu leaves
    curr_mu_data = forest_mu.leaf_tree * forest_mu.leaf_unit  # (m_mu, max_slots_mu)
    A_mu_flat = jnp.where(is_leaf_mu[:, :, None], A_mu, 0.0).reshape(-1, K)
    curr_mu_flat = jnp.where(is_leaf_mu, curr_mu_data, 0.0).reshape(-1)
    prec_mu = leaf_prior_cov_inv_mu * (A_mu_flat.T @ A_mu_flat)
    h_mu = -leaf_prior_cov_inv_mu * (A_mu_flat.T @ curr_mu_flat)

    # Prior contribution from nu leaves
    curr_nu_data = forest_nu.leaf_tree * forest_nu.leaf_unit  # (m_nu, max_slots_nu)
    A_nu_flat = jnp.where(is_leaf_nu[:, :, None], A_nu, 0.0).reshape(-1, K)
    curr_nu_flat = jnp.where(is_leaf_nu, curr_nu_data, 0.0).reshape(-1)
    prec_nu = leaf_prior_cov_inv_nu * (A_nu_flat.T @ A_nu_flat)
    h_nu = -leaf_prior_cov_inv_nu * (A_nu_flat.T @ curr_nu_flat)

    prec_total = prec_lik + prec_mu + prec_nu
    h_total = h_lik + h_mu + h_nu

    # 4. Sample delta ~ N(Sigma_delta @ h_total, Sigma_delta)
    delta = sample_from_precision(key, prec_total, h_total)

    # 5. Apply exact subspace update to leaf trees, forest fits, and residual
    shift_mu_data = jnp.einsum('tsk,k->ts', A_mu, delta)
    new_leaf_tree_mu = forest_mu.leaf_tree + shift_mu_data / forest_mu.leaf_unit
    forest_mu_new = eqx.tree_at(lambda f: f.leaf_tree, forest_mu, new_leaf_tree_mu)
    mu_fit_new = mu_fit + J_mu @ delta

    shift_nu_data = jnp.einsum('tsk,k->ts', A_nu, delta)
    new_leaf_tree_nu = forest_nu.leaf_tree + shift_nu_data / forest_nu.leaf_unit
    forest_nu_new = eqx.tree_at(lambda f: f.leaf_tree, forest_nu, new_leaf_tree_nu)
    nu_fit_new = nu_fit + J_nu @ delta

    R_new = jnp.where(obs_mask, R - H @ delta, 0.0)

    return forest_mu_new, mu_fit_new, forest_nu_new, nu_fit_new, R_new


def inter_ensemble_beta_gamma_step(
    key: Key[Array, ''],
    forest_mu: Any,
    mu_fit: Float32[Array, ' n'],
    gamma: Float32[Array, ' N_units'],
    beta: Float32[Array, ' S_max_plus_1'],
    R: Float32[Array, ' n'],
    alpha: Float32[Array, ''],
    nu_fit: Float32[Array, ' n'],
    b_z: Float32[Array, ' n'],
    obs_mask: Bool[Array, ' n'],
    z_vec: Float32[Array, ' n'],
    time_idx: Int32[Array, ' n'],
    T_periods: int,
    unit_idx: Int32[Array, ' n'],
    N_units: int,
    exposure_idx: Int32[Array, ' n'],
    S_max: int,
    sigma2: Float32[Array, ''],
    leaf_prior_cov_inv_mu: Float32[Array, ''],
    sigma_gamma2: Float32[Array, ''],
    K_chol: Float32[Array, 'S_max_plus_1 S_max_plus_1'],
    random_intercept: bool,
    sample_beta: bool,
    conditional_precision: Float32[Array, ' n'] | None = None,
    temperature: Float32[Array, ''] = jnp.float32(1.0),
    X_unified: Any = None,
    y_vec: Float32[Array, ' n'] | None = None,
    ar1: Any = None,
) -> tuple[
    Any,
    Float32[Array, ' n'],
    Float32[Array, ' N_units'],
    Float32[Array, ' S_max_plus_1'],
    Float32[Array, ' n'],
]:
    """Execute exact conjugate subspace Gibbs step jointly transferring trend between mu, gamma, and beta."""
    denom_t = jnp.float32(max(1, T_periods - 1))
    t_norm = (time_idx.astype(jnp.float32) - 0.5 * denom_t) / denom_t
    s_norm = exposure_idx.astype(jnp.float32) / jnp.float32(max(1, S_max))
    g_cell = _compute_unit_prognostic_score(
        X_unified, y_vec, z_vec, obs_mask, unit_idx, N_units
    )

    Phi_base = jnp.stack(
        [
            jnp.ones_like(t_norm),
            t_norm,
            s_norm,
            g_cell * t_norm,
        ],
        axis=-1,
    )  # shape (n, 4)

    ids_mu = apply_moves_to_leaf_indices(
        forest_mu.leaf_indices, forest_mu.to_prune, forest_mu.move_node
    )
    _, _, is_leaf_mu = compute_active_leaf_stats(
        forest_mu.split_tree, forest_mu.leaf_tree
    )

    max_slots_mu = forest_mu.leaf_tree.shape[-1]
    m_mu = jnp.float32(forest_mu.leaf_tree.shape[0])

    # 1. Prognostic leaf basis directions A_mu (data units), shape (m_mu, max_slots_mu, K=8)
    trt_obs_mask = obs_mask & (z_vec == 1.0)
    Phi_mean_mu, cnt_all_mu = jax.vmap(
        _leaf_stats_vec, in_axes=(0, None, None, None)
    )(ids_mu, Phi_base, obs_mask, max_slots_mu)
    _, cnt_trt_mu = jax.vmap(
        _leaf_stats_vec, in_axes=(0, None, None, None)
    )(ids_mu, Phi_base[:, :1], trt_obs_mask, max_slots_mu)

    p_trt_mu = jnp.where(
        cnt_all_mu > 0, cnt_trt_mu / jnp.maximum(cnt_all_mu, 1.0), 0.0
    )
    active_mu = is_leaf_mu & (cnt_trt_mu > 0)
    w_prop = jnp.where(active_mu, p_trt_mu / m_mu, 0.0)
    w_unif = jnp.where(active_mu, 1.0 / m_mu, 0.0)
    A_mu = jnp.concatenate(
        [
            w_prop[:, :, None] * Phi_mean_mu,
            w_unif[:, :, None] * Phi_mean_mu,
        ],
        axis=-1,
    )  # (m_mu, max_slots_mu, K=8)
    K = A_mu.shape[-1]

    J_mu_per_tree = jax.vmap(lambda a_tree, idx: a_tree[idx])(A_mu, ids_mu)
    J_mu = jnp.sum(J_mu_per_tree, axis=0)  # (n, K)

    # 2. Unit random intercept basis G anchored strictly on untreated cells (Z_it = 0)
    if random_intercept:
        ctrl_obs_mask = obs_mask & (z_vec == 0.0)
        ctrl_mask_f32 = ctrl_obs_mask.astype(jnp.float32)
        alpha_J_mu_ctrl = jnp.where(ctrl_obs_mask[:, None], alpha * J_mu, 0.0)
        sum_ctrl_J = (
            jnp.zeros((N_units, K), dtype=jnp.float32)
            .at[unit_idx]
            .add(alpha_J_mu_ctrl)
        )
        cnt_ctrl = (
            jnp.zeros(N_units, dtype=jnp.float32)
            .at[unit_idx]
            .add(ctrl_mask_f32)
        )
        G = jnp.where(
            cnt_ctrl[:, None] > 0,
            -sum_ctrl_J / jnp.maximum(cnt_ctrl[:, None], 1.0),
            0.0,
        )
        J_gamma = G[unit_idx]  # (n, K)
    else:
        G = jnp.zeros((N_units, K), dtype=jnp.float32)
        J_gamma = jnp.zeros((J_mu.shape[0], K), dtype=jnp.float32)

    # 3. Exposure GP trajectory basis B, shape (S_max + 1, K)
    d_vec = b_z * nu_fit  # (n,)
    if sample_beta:
        J_base = alpha * J_mu + J_gamma  # (n, K)
        d_trt = jnp.where(trt_obs_mask, d_vec, 0.0)
        num_B = (
            jnp.zeros((S_max + 1, K), dtype=jnp.float32)
            .at[exposure_idx]
            .add(d_trt[:, None] * J_base)
        )
        den_B = (
            jnp.zeros(S_max + 1, dtype=jnp.float32)
            .at[exposure_idx]
            .add(jnp.square(d_trt))
        )
        valid_s = (jnp.arange(S_max + 1) >= 1) & (den_B > 1e-6)
        B = jnp.where(
            valid_s[:, None],
            -num_B / (den_B[:, None] + 1e-4),
            0.0,
        )
        J_beta = d_vec[:, None] * B[exposure_idx]  # (n, K)
    else:
        B = jnp.zeros((S_max + 1, K), dtype=jnp.float32)
        J_beta = jnp.zeros((J_mu.shape[0], K), dtype=jnp.float32)

    # Total residual sensitivity H = alpha * J_mu + J_gamma + J_beta: shape (n, K)
    H = jnp.where(
        obs_mask[:, None],
        alpha * J_mu + J_gamma + J_beta,
        0.0,
    )

    # 4. Assemble exact Gaussian conditional precision and linear term
    if conditional_precision is None:
        cell_prec = jnp.where(obs_mask, temperature / sigma2, 0.0)
    else:
        cell_prec = jnp.where(obs_mask, conditional_precision, 0.0)

    if ar1 is None:
        H_lik, R_lik = H, R
    else:
        # AR(1) innovations: the likelihood is Gaussian in the Prais-Winsten
        # transformed sensitivities and residual (longbet._ar1); the update
        # itself is applied to the untransformed parameters below.
        from longbet._ar1 import transform as _ar1_transform
        phi_ar1, first_ar1 = ar1
        H_lik = _ar1_transform(H, phi_ar1, first_ar1)
        R_lik = _ar1_transform(jnp.where(obs_mask, R, 0.0), phi_ar1, first_ar1)
    H_weighted = H_lik * cell_prec[:, None]
    prec_lik = H_lik.T @ H_weighted  # (K, K)
    h_lik = H_weighted.T @ R_lik  # (K,)

    # Prior contribution from mu leaves
    curr_mu_data = forest_mu.leaf_tree * forest_mu.leaf_unit  # (m_mu, max_slots_mu)
    A_mu_flat = jnp.where(is_leaf_mu[:, :, None], A_mu, 0.0).reshape(-1, K)
    curr_mu_flat = jnp.where(is_leaf_mu, curr_mu_data, 0.0).reshape(-1)
    prec_mu = leaf_prior_cov_inv_mu * (A_mu_flat.T @ A_mu_flat)
    h_mu = -leaf_prior_cov_inv_mu * (A_mu_flat.T @ curr_mu_flat)

    # Prior contribution from gamma
    if random_intercept:
        inv_sigma_gamma2 = jnp.reciprocal(sigma_gamma2)
        prec_gamma = inv_sigma_gamma2 * (G.T @ G)
        h_gamma = -inv_sigma_gamma2 * (G.T @ gamma)
    else:
        prec_gamma = jnp.zeros((K, K), dtype=jnp.float32)
        h_gamma = jnp.zeros((K,), dtype=jnp.float32)

    # Prior contribution from beta GP: || L_gp^{-1} (beta + B theta) ||^2
    if sample_beta:
        L_gp = jnp.asarray(K_chol, dtype=jnp.float32)
        v_beta = jax.scipy.linalg.solve_triangular(
            L_gp, beta[:, None], lower=True
        ).squeeze(-1)  # (S_max + 1,)
        V_B = jax.scipy.linalg.solve_triangular(L_gp, B, lower=True)  # (S_max + 1, K)
        prec_beta = V_B.T @ V_B  # (K, K)
        h_beta = -(V_B.T @ v_beta)  # (K,)
    else:
        prec_beta = jnp.zeros((K, K), dtype=jnp.float32)
        h_beta = jnp.zeros((K,), dtype=jnp.float32)

    prec_total = (
        prec_lik
        + prec_mu
        + prec_gamma
        + prec_beta
    )
    h_total = h_lik + h_mu + h_gamma + h_beta

    # 5. Sample theta ~ N(Sigma_theta @ h_total, Sigma_theta)
    #
    # The subspace is not guaranteed to have full rank. On a single-period
    # panel the normalized calendar coordinate is constant, so the {1, t}
    # directions coincide and two of the eight directions vanish. See
    # ``longbet._linalg``.
    theta = sample_from_precision(key, prec_total, h_total)

    # 6. Apply exact subspace updates
    shift_mu_data = jnp.einsum('tsk,k->ts', A_mu, theta)
    new_leaf_tree_mu = forest_mu.leaf_tree + shift_mu_data / forest_mu.leaf_unit
    forest_mu_new = eqx.tree_at(lambda f: f.leaf_tree, forest_mu, new_leaf_tree_mu)
    mu_fit_new = mu_fit + J_mu @ theta
    gamma_new = gamma + G @ theta
    beta_new = beta + B @ theta
    R_new = jnp.where(obs_mask, R - H @ theta, 0.0)

    return forest_mu_new, mu_fit_new, gamma_new, beta_new, R_new

