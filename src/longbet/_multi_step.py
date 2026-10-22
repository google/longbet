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

"""One coupled Gibbs sweep for multiple outcomes in LongBet."""

from __future__ import annotations

from typing import Any

import equinox as eqx
import jax

from longbet._x64 import enable_x64
import jax.numpy as jnp
from jax import random
from jaxtyping import Array, Float32, Int32, Key

from longbet._linalg import sample_from_precision
from longbet._multi_state import (
    MultiLongBetState,
    multi_chain_filter_spec,
    split_multi_chain_fields,
)
from longbet._state import LongBetState
from longbet._step import longbet_single_step, _sample_inv_gamma, _tempered_single_step
from longbet._sur import conditional_residual, innovation_residual


def sample_loadings(
    key: Key[Array, ''],
    *,
    predecessors: list[LongBetState],
    response: LongBetState,
    sigma2: Float32[Array, ''],
    prior_var: float,
) -> Float32[Array, ' m']:
    """Draw one row of Gamma from its conjugate Gaussian posterior.

    Parameters
    ----------
    key
        PRNG key.
    predecessors
        List of predecessor ``LongBetState`` objects (internal order 0..m-1).
    response
        Current outcome's ``LongBetState``.
    sigma2
        Innovation variance of the current equation.
    prior_var
        Marginal prior variance for each element of Gamma[m, :m].

    Returns
    -------
    Float32 array of shape ``(m,)`` representing Gamma[m, :m].
    """
    m = len(predecessors)
    if m == 0:
        return jnp.zeros((0,), dtype=jnp.float32)
    if prior_var <= 0.0:
        return jnp.zeros((m,), dtype=jnp.float32)

    # Stack predecessor raw residuals: shape (Q, m)
    X_raw = jnp.column_stack([p.resid for p in predecessors])
    r = response.resid
    mask = response.obs_mask

    # Mask unobserved cells
    Xo = jnp.where(mask[:, None], X_raw, 0.0)
    ro = jnp.where(mask, r, 0.0)

    # Posterior precision: Xo.T @ Xo / sigma2 + I_m / prior_var
    P = (Xo.T @ Xo) / sigma2 + jnp.eye(m, dtype=jnp.float32) / prior_var
    # Symmetrize against roundoff
    P = 0.5 * (P + P.T)

    h = (Xo.T @ ro) / sigma2

    # The prior term makes P positive definite in exact arithmetic, but two
    # predecessors with near-identical residuals and a small sigma2 put the
    # likelihood term many orders of magnitude above it, which is enough for a
    # float32 Cholesky to fail. See ``longbet._linalg``.
    return sample_from_precision(key, P, h)


@enable_x64()
@jax.jit
@enable_x64(False)
def multi_single_step(
    outcome_step_keys: Key[Array, 'M ...'],
    loading_keys: Key[Array, 'M ...'],
    multi_state: MultiLongBetState,
    sweep_index: Int32[Array, ''],
) -> MultiLongBetState:
    """Execute one coupled Gibbs transition on a single chain.

    Full-precision SUR update, including downstream likelihood feedback into
    earlier means and binary latent responses. Store raw residuals at state
    boundaries, use temporary conditional residuals/precisions for local mean
    updates, and draw structural innovation variances from their own equations.
    """
    states = list(multi_state.states)
    Gamma = multi_state.gamma_loadings
    M = multi_state.M
    coupled = multi_state.sur_active and any(multi_state.continuous_mask[1:])
    observed = jnp.stack([st.obs_mask for st in states])
    variances = jnp.stack([st.sigma2 for st in states])
    # Parallel tempering: the replica's inverse temperature scales
    # every likelihood precision; see longbet._tempering.
    tempered = states[0].is_tempered
    beta = jnp.asarray(states[0].temperature, dtype=jnp.float32) if tempered else None

    for m in range(M):
        st = states[m]
        mask = st.obs_mask
        if not coupled:
            # Preserve the scalar sampler exactly, including its key schedule.
            states[m] = (_tempered_single_step(outcome_step_keys[m], st) if tempered
                         else longbet_single_step(outcome_step_keys[m], st))
            continue

        raw = jnp.stack([child.resid for child in states])
        resid, precision = conditional_residual(raw, observed, Gamma, variances, m)
        if tempered:
            precision = jnp.where(mask, precision * beta, precision)
        offset = jnp.where(mask, st.resid - resid, 0.0)
        adjusted = eqx.tree_at(
            lambda s: s.resid, st, resid
        )
        updated = longbet_single_step(outcome_step_keys[m], adjusted, precision)
        states[m] = eqx.tree_at(
            lambda s: s.resid, updated, jnp.where(mask, updated.resid + offset, 0.0)
        )

    # Mean/latent blocks are now complete. Gamma and innovation variances each
    # condition on these current raw residuals. In particular, do not estimate
    # a structural variance from an early equation's conditional residual.
    if coupled:
        raw = jnp.stack([st.resid for st in states])
        for m in range(M):
            st = states[m]
            if not multi_state.continuous_mask[m]:
                continue  # identifying variance stays one; no incoming loadings
            if m > 0:
                row_draw = sample_loadings(
                    loading_keys[m], predecessors=states[:m], response=st,
                    sigma2=(st.sigma2 / beta if tempered else st.sigma2),
                    prior_var=multi_state.sur_prior_var,
                )
                # Retain the historical initialization convention only on the
                # first burn-in sweep; all subsequent sweeps update loadings.
                Gamma = Gamma.at[m, :m].set(
                    jnp.where(sweep_index > 0, row_draw, Gamma[m, :m])
                )
            innovation = jnp.where(st.obs_mask, innovation_residual(raw, Gamma, m), 0.0)
            scale = beta if tempered else 1.0
            sigma2 = _sample_inv_gamma(
                random.split(outcome_step_keys[m], 10)[9],
                st.sigma_prior_a + scale * jnp.sum(st.obs_mask) / 2.0,
                st.sigma_prior_b + 0.5 * scale * jnp.sum(innovation**2),
            )
            states[m] = eqx.tree_at(
                lambda s: (s.sigma2, s.error_cov_inv.value), st,
                (sigma2, jnp.reciprocal(sigma2).astype(jnp.float32)),
            )

    return eqx.tree_at(
        lambda s: (s.states, s.gamma_loadings),
        multi_state,
        (tuple(states), Gamma),
    )


@enable_x64()
def multi_step(
    key: Key[Array, ''],
    multi_state: MultiLongBetState,
    sweep_index: Int32[Array, ''] = jnp.int32(1),
) -> MultiLongBetState:
    """Execute one Gibbs sweep across one or multiple chains.

    Derives independent keys for outcome equations and loading draws according to
    Section 6.1. Vmaps over chains when C > 1.
    """
    k_outcome = random.fold_in(key, 0)
    k_loading = random.fold_in(key, 1)
    M = multi_state.M

    if multi_state.has_chain_axis:
        C = multi_state.num_chains
        outcome_chains = random.split(k_outcome, C)
        loading_chains = random.split(k_loading, C)

        # Per-chain key arrays shaped (C, M, 2)
        outcome_keys = jnp.stack(
            [
                jnp.stack([random.fold_in(outcome_chains[c], m) for m in range(M)])
                for c in range(C)
            ]
        )
        loading_keys = jnp.stack(
            [
                jnp.stack([random.fold_in(loading_chains[c], m) for m in range(M)])
                for c in range(C)
            ]
        )

        per_chain, shared = split_multi_chain_fields(multi_state)

        def _vmapped_step(per: Any, k_out: Any, k_load: Any) -> Any:
            single_st = eqx.combine(per, shared)
            updated = multi_single_step(k_out, k_load, single_st, sweep_index)
            per_updated, _ = split_multi_chain_fields(updated)
            return per_updated

        new_per = jax.vmap(_vmapped_step, in_axes=(0, 0, 0))(
            per_chain, outcome_keys, loading_keys
        )
        return eqx.combine(new_per, shared)

    # Single chain
    outcome_keys = jnp.stack([random.fold_in(k_outcome, m) for m in range(M)])
    loading_keys = jnp.stack([random.fold_in(k_loading, m) for m in range(M)])
    return multi_single_step(outcome_keys, loading_keys, multi_state, sweep_index)
