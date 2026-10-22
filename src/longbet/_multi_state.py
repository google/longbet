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

"""Multi-outcome state and initialization for coupled LongBet.

Child states are held in INTERNAL outcome order (binary first, preserving caller
order within each type, followed by continuous).
"""

from __future__ import annotations

import dataclasses
from typing import Any

import equinox as eqx
import jax
import jax.numpy as jnp
from jaxtyping import Array, Bool, Float32, Int32, Key, UInt

from bartz._jaxext import field
from bartz.mcmcstep._axes import CHAIN_AXIS

from longbet._config import LongBetConfig
from longbet._multi_input import NormalizedMultiInput
from longbet._state import (
    LongBetState,
    _overdisperse_chains,
    chain_filter_spec,
    init_longbet,
)


class MultiLongBetState(eqx.Module):
    """The full MCMC state for coupled multi-outcome LongBet.

    Attributes
    ----------
    states
        Tuple of child ``LongBetState`` objects in INTERNAL outcome order.
    gamma_loadings
        SUR regression loadings matrix, shaped ``(M, M)`` for a single chain or
        ``(C, M, M)`` for C parallel chains. Strictly lower-triangular with a
        zero diagonal.
    M
        Number of outcomes.
    sur_active
        Whether SUR coupling is enabled and active.
    sur_prior_var
        Prior variance on SUR regression loadings.
    continuous_mask
        Tuple of bools of length M indicating which equations are continuous.
    """

    states: tuple[LongBetState, ...]
    gamma_loadings: Float32[Array, '*chains M M'] = field(chains=CHAIN_AXIS)

    M: int = field(static=True)
    sur_active: bool = field(static=True)
    sur_prior_var: float = field(static=True)
    continuous_mask: tuple[bool, ...] = field(static=True)

    @property
    def has_chain_axis(self) -> bool:
        """Whether this state carries a leading chain axis."""
        return self.gamma_loadings.ndim > 2

    @property
    def num_chains(self) -> int:
        """Number of chains represented by this state."""
        return self.gamma_loadings.shape[0] if self.has_chain_axis else 1


def multi_chain_filter_spec(state: MultiLongBetState) -> Any:
    """Boolean pytree marking the leaves of MultiLongBetState that carry a chain axis."""
    return chain_filter_spec(state)


def split_multi_chain_fields(state: MultiLongBetState) -> tuple[Any, Any]:
    """Split into ``(per_chain, shared)`` pytrees."""
    return eqx.partition(state, multi_chain_filter_spec(state))


def broadcast_multi_to_chains(
    state: MultiLongBetState,
    num_chains: int,
    key: Any = None,
) -> MultiLongBetState:
    """Replicate the chain-axis fields of a single-chain MultiLongBetState.

    Shared quantities (X, indices, K_chol) remain un-replicated across chains.
    When ``key`` is given, each child's parameters are overdispersed according to
    its prior and fixed settings.
    """
    if num_chains <= 1:
        return state
    if state.has_chain_axis:
        raise ValueError("MultiLongBetState already carries a chain axis")

    per, shared = split_multi_chain_fields(state)
    per = jax.tree.map(lambda x: jnp.broadcast_to(x, (num_chains, *x.shape)), per)
    state = eqx.combine(per, shared)

    if key is None:
        return state

    return _overdisperse_multi_chains(state, key, num_chains)


def _overdisperse_multi_chains(
    state: MultiLongBetState, key: Any, num_chains: int
) -> MultiLongBetState:
    """Overdisperse each child's chain start respecting fixed options."""
    new_states = []
    for m in range(state.M):
        k_m = jax.random.fold_in(key, m)
        child = state.states[m]
        new_child = _overdisperse_chains(child, k_m, num_chains)
        new_states.append(new_child)

    return eqx.tree_at(
        lambda s: s.states,
        state,
        tuple(new_states),
    )


def init_multi_longbet(
    *,
    X_unified: UInt[Array, 'p n'],
    unit_idx: Int32[Array, ' n'],
    time_idx: Int32[Array, ' n'],
    exposure_idx: Int32[Array, ' n'],
    z_vec: Float32[Array, ' n'],
    max_split_mu: UInt[Array, ' p'],
    max_split_nu: UInt[Array, ' p'],
    norm_input: NormalizedMultiInput,
    config: LongBetConfig,
    trend_design: Float32[Array, 'n q_tot'] | None = None,
    trend_num_period: int = 0,
    num_chains: int | None = None,
    init_key: Any = None,
    mesh: Any = None,
    **kwargs: Any,
) -> MultiLongBetState:
    """Construct the initial ``MultiLongBetState``.

    Shares the unified design matrix X, panel index vectors, and split bounds
    across all M outcome equations.

    ``trend_design`` is the calendar-time block of
    :func:`longbet._trend_basis.build_trend_design`. It depends only on ``(x,
    t)``, which every equation shares, so the *same* array is handed to all M
    children -- exactly as the scalar fit hands one array to its single
    equation. Each child then masks it to its own observed cells inside
    :func:`longbet._state.init_longbet`, which is what makes an outcome with
    missing cells behave like the scalar fit on that same outcome.

    Passing it is not optional in practice. With the shipped defaults
    :attr:`longbet.LongBetConfig.split_calendar_mu` resolves to ``False``, so
    the prognostic forest cannot split on calendar time; if the block were not
    wired in as well, nothing in the model could represent a calendar-time
    movement at all.
    """
    M = norm_input.M
    continuous_mask = tuple(t == "continuous" for t in norm_input.internal_outcomes)
    if config.baseline_cells != "all":
        raise NotImplementedError(
            "baseline_cells='untreated' is not available for LongBetMulti: the "
            "coupled likelihood shares the treated cells across equations"
        )

    child_states = []
    for m in range(M):
        otype = norm_input.internal_outcomes[m]
        child_cfg = dataclasses.replace(config, outcome=otype,
            num_categories=norm_input.internal_num_categories[m])
        y_m = norm_input.y_prepared[m].ravel()
        obs_m = norm_input.obs_masks[m].ravel()
        offset_m = norm_input.offset_[m]

        st = init_longbet(
            X_unified=X_unified,
            y=y_m,
            unit_idx=unit_idx,
            time_idx=time_idx,
            exposure_idx=exposure_idx,
            z_vec=z_vec,
            obs_mask=obs_m,
            max_split_mu=max_split_mu,
            max_split_nu=max_split_nu,
            config=child_cfg,
            offset=offset_m,
            trend_design=trend_design,
            trend_num_period=trend_num_period,
            num_chains=None,
            chain_key=(jax.random.fold_in(init_key, m)
                       if init_key is not None and (num_chains is None or num_chains <= 1) else None),
            mesh=mesh,
            **kwargs,
        )
        child_states.append(st)

    gamma_loadings = jnp.zeros((M, M), dtype=jnp.float32)

    multi_state = MultiLongBetState(
        states=tuple(child_states),
        gamma_loadings=gamma_loadings,
        M=M,
        sur_active=config.sur_active,
        sur_prior_var=float(config.sur_prior_var),
        continuous_mask=continuous_mask,
    )

    if config.tempering_levels > 1 and (num_chains is None or num_chains % config.tempering_levels):
        raise ValueError("tempering needs num_chains * tempering_levels replicas")
    if num_chains is not None and num_chains > 1:
        multi_state = broadcast_multi_to_chains(
            multi_state, int(num_chains), key=init_key
        )

    if config.tempering_levels > 1:
        from longbet._tempering import ladder_temperatures, multi_set_temperature
        multi_state = multi_set_temperature(
            multi_state, ladder_temperatures(multi_state.num_chains, config.tempering_levels,
                                             config.tempering_beta_min))
    return multi_state
