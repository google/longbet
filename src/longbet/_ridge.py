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

"""Global-scale Metropolis move along the beta-nu ridge.

The treatment term enters as a product ``beta_S * nu(X_i, S, t)`` and the
likelihood does not pin down the split: for any positive ``c``,
``(beta, nu) -> (c * beta, nu / c)`` reproduces the fit exactly.  A Metropolis
grow/prune sampler traverses that ridge slowly with local tree moves alone, so
this module gives it an explicit move along the manifold.

Because the likelihood is invariant, the acceptance ratio needs only the two
priors and the Jacobian.  The move costs ``O(S_max^2 + L)`` and can cross the
ridge in a single accepted step.
"""

from __future__ import annotations

from typing import Any

import equinox as eqx
import jax
import jax.numpy as jnp
from jaxtyping import Array, Bool, Float, Float32, Int32, Key, UInt

from bartz.grove._grove import is_actual_leaf

from longbet._gp import beta_prior_quadform


def compute_active_leaf_stats(
    split_tree: UInt[Array, 'num_trees half_tree_size'],
    leaf_tree: Float[Array, 'num_trees 2*half_tree_size'],
) -> tuple[Int32[Array, ''], Float[Array, ''], Bool[Array, 'num_trees 2*half_tree_size']]:
    """Count the forest's actual leaves and sum their squares.

    Returns the count, the sum of squared leaf values (in storage units), and
    the boolean mask of actual-leaf slots.  The mask is returned so the caller
    can rescale exactly the values the Jacobian accounts for.
    """
    st_flat = split_tree.reshape(-1, split_tree.shape[-1])
    is_leaf_flat = jax.vmap(lambda s: is_actual_leaf(s, add_bottom_level=True))(st_flat)
    is_leaf = is_leaf_flat.reshape((*split_tree.shape[:-1], is_leaf_flat.shape[-1]))

    num_active = jnp.sum(is_leaf, axis=(-2, -1))
    sum_sq = jnp.sum(jnp.square(leaf_tree) * is_leaf, axis=(-2, -1))
    return num_active, sum_sq, is_leaf


def ridge_scale_step(
    key: Key[Array, ''],
    beta: Float32[Array, ' S_max_plus_1'],
    forest_nu: Any,
    nu_fit: Float32[Array, ' n'],
    K_chol: Float32[Array, 'S_max_plus_1 S_max_plus_1'],
    leaf_prior_cov_inv_nu: Float32[Array, ''],
    proposal_sigma: float = 0.2,
    *,
    delta: Float32[Array, ' p'] | None = None,
    lin_fit: Float32[Array, ' n'] | None = None,
    delta_prior_prec: Float32[Array, ''] | None = None,
    return_extras: bool = False,
) -> Any:
    """One Metropolis step on ``log c`` for ``(beta, ell) -> (c beta, ell / c)``.

    With ``beta ~ N(0, K_tilde)`` over ``d = S_max + 1`` points and ``L`` active
    leaves ``ell ~ N(0, 1 / leaf_prior_cov_inv)``, the log acceptance ratio is

    .. code::

        -1/2 (c^2 - 1) beta^T K_tilde^-1 beta
        -1/2 leaf_prec (c^-2 - 1) sum(ell^2)
        + (d - L) log c

    the last term being the Jacobian of the map.  The proposal on ``log c`` is
    symmetric, so it drops out.

    Only the *actual leaves* are rescaled, matching the ``L`` counted in the
    Jacobian; unreachable leaf slots are left alone (`bartz` redraws leaf values
    when a node is grown, so their contents never enter the model).

    ``nu_fit`` is rescaled alongside, and the full-model residual is left
    untouched -- it is invariant by construction, which is the point of the move.

    Parameters
    ----------
    key
        PRNG key.
    beta
        Current exposure trajectory, shape ``(S_max + 1,)``.
    forest_nu
        Current treatment forest.
    nu_fit
        Current treatment forest fit, shape ``(n,)``.
    K_chol
        Lower Cholesky factor of the marginalized GP kernel.
    leaf_prior_cov_inv_nu
        Prior precision of the treatment forest's leaves, in data units.
    proposal_sigma
        Standard deviation of the ``log c`` proposal.
    delta, lin_fit, delta_prior_prec
        The moderator's linear part (``LongBetConfig.linear_trt``): its
        coefficients, its fitted values and the prior precision
        ``1 / tau_delta^2`` of each coefficient. The map then rescales
        ``delta -> delta / c`` too, which adds
        ``-1/2 delta_prec (c^-2 - 1) |delta|^2 - p log c`` to the log ratio.
    return_extras
        Also return ``(delta_new, lin_fit_new, accepted)``.

    Returns
    -------
    beta_new, forest_nu_new, nu_fit_new[, delta_new, lin_fit_new, accepted]
    """
    k_prop, k_acc = jax.random.split(key)

    log_c = jax.random.normal(k_prop, shape=(), dtype=jnp.float32) * proposal_sigma
    c = jnp.exp(log_c)

    # beta^T K_tilde^-1 beta by triangular solve; never forms the precision.
    beta_quad = beta_prior_quadform(beta, K_chol)
    log_p_beta_diff = -0.5 * (jnp.square(c) - 1.0) * beta_quad

    num_active, sum_sq, is_leaf = compute_active_leaf_stats(
        forest_nu.split_tree, forest_nu.leaf_tree
    )
    sum_sq_data = sum_sq * jnp.square(forest_nu.leaf_unit)
    log_p_leaves_diff = (
        -0.5 * leaf_prior_cov_inv_nu * (jnp.reciprocal(jnp.square(c)) - 1.0) * sum_sq_data
    )

    d = beta.shape[-1]
    log_jacobian = (d - num_active.astype(jnp.float32)) * log_c

    log_alpha = log_p_beta_diff + log_p_leaves_diff + log_jacobian
    if delta is not None:
        inv_c2_minus_1 = jnp.reciprocal(jnp.square(c)) - 1.0
        log_alpha = (log_alpha
                     - 0.5 * inv_c2_minus_1 * jnp.sum(delta_prior_prec * jnp.square(delta))
                     - jnp.float32(delta.shape[-1]) * log_c)
    log_u = jnp.log(
        jax.random.uniform(k_acc, shape=(), dtype=jnp.float32) + jnp.finfo(jnp.float32).tiny
    )
    accept = log_u < log_alpha
    c_eff = jnp.where(accept, c, 1.0)

    beta_new = beta * c_eff
    new_leaf_tree = jnp.where(is_leaf, forest_nu.leaf_tree / c_eff, forest_nu.leaf_tree)
    forest_nu_new = eqx.tree_at(lambda f: f.leaf_tree, forest_nu, new_leaf_tree)
    nu_fit_new = nu_fit / c_eff

    if not return_extras:
        return beta_new, forest_nu_new, nu_fit_new
    delta_new = None if delta is None else delta / c_eff
    lin_fit_new = None if lin_fit is None else lin_fit / c_eff
    return beta_new, forest_nu_new, nu_fit_new, delta_new, lin_fit_new, accept


def scale_group_step(
    key: Key[Array, ''],
    beta: Float32[Array, ' S_max_plus_1'],
    forest_nu: Any,
    nu_fit: Float32[Array, ' n'],
    K_chol_base: Float32[Array, 'S_max_plus_1 S_max_plus_1'],
    leaf_prior_cov_inv_nu: Float32[Array, ''],
    proposal_sigma: float = 0.3,
    *,
    beta_scale2: Float32[Array, ''] | None = None,
    beta_scale_aux: Float32[Array, ''] | None = None,
    delta: Float32[Array, ' p'] | None = None,
    lin_fit: Float32[Array, ' n'] | None = None,
    tau_delta2: Float32[Array, ''] | None = None,
    tau_delta_aux: Float32[Array, ''] | None = None,
) -> tuple[Any, ...]:
    """Metropolis move along the ridge *with* the learned scales (encourage-v2).

    With a learned trajectory scale ``s`` (``beta ~ N(0, s^2 K)``) and/or a
    learned scale ``tau`` of the moderator's linear coefficients, the ridge
    ``(beta, nu) -> (c beta, nu / c)`` gains a funnel: the scales can trade
    against each other along it at nearly constant density, and the plain
    ridge move, which holds them fixed, crawls. This move carries them along:

    .. code::

        beta -> c beta,  s -> c s,  tree leaves -> leaves / c,
        delta -> delta / c,  tau -> tau / c

    so ``beta / s`` and ``delta / tau`` are unchanged and the fit is invariant.
    With the half-Cauchy scales in their inverse-gamma mixture form
    (``s^2 | xi ~ IG(1/2, 1/xi)``), the log acceptance ratio, priors plus the
    Jacobian, is

    .. code::

        -1/2 leaf_prec (c^-2 - 1) sum(ell^2) - L log c
        - (c^-2 - 1) / (xi_s s^2) - log c                 (learned s)
        - 1/2 (c^2 - 1) beta' K^-1 beta + d log c          (fixed kernel)
        - (c^2 - 1) / (xi_tau tau^2) + log c              (linear part)

    Returns ``(beta, forest_nu, nu_fit, beta_scale2, delta, lin_fit, tau_delta2, accepted)``.
    """
    k_prop, k_acc = jax.random.split(key)
    log_c = jax.random.normal(k_prop, shape=(), dtype=jnp.float32) * proposal_sigma
    c = jnp.exp(log_c)
    inv_c2_m1 = jnp.reciprocal(jnp.square(c)) - 1.0
    c2_m1 = jnp.square(c) - 1.0

    num_active, sum_sq, is_leaf = compute_active_leaf_stats(forest_nu.split_tree, forest_nu.leaf_tree)
    sum_sq_data = sum_sq * jnp.square(forest_nu.leaf_unit)
    log_alpha = (-0.5 * leaf_prior_cov_inv_nu * inv_c2_m1 * sum_sq_data
                 - num_active.astype(jnp.float32) * log_c)
    if beta_scale2 is not None:
        log_alpha = log_alpha - inv_c2_m1 / (beta_scale_aux * beta_scale2) - log_c
    else:
        quad = beta_prior_quadform(beta, K_chol_base)
        log_alpha = log_alpha - 0.5 * c2_m1 * quad + jnp.float32(beta.shape[-1]) * log_c
    if delta is not None:
        log_alpha = log_alpha - c2_m1 / (tau_delta_aux * tau_delta2) + log_c
    log_u = jnp.log(jax.random.uniform(k_acc, shape=(), dtype=jnp.float32) + jnp.finfo(jnp.float32).tiny)
    accept = log_u < log_alpha
    c_eff = jnp.where(accept, c, 1.0)

    beta_new = beta * c_eff
    new_leaf_tree = jnp.where(is_leaf, forest_nu.leaf_tree / c_eff, forest_nu.leaf_tree)
    forest_nu_new = eqx.tree_at(lambda f: f.leaf_tree, forest_nu, new_leaf_tree)
    nu_fit_new = nu_fit / c_eff
    scale2_new = None if beta_scale2 is None else beta_scale2 * jnp.square(c_eff)
    delta_new = None if delta is None else delta / c_eff
    lin_new = None if lin_fit is None else lin_fit / c_eff
    tau2_new = None if tau_delta2 is None else tau_delta2 / jnp.square(c_eff)
    return beta_new, forest_nu_new, nu_fit_new, scale2_new, delta_new, lin_new, tau2_new, accept
