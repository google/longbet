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

"""MCMC driver for coupled multi-outcome LongBet."""

from __future__ import annotations

from collections.abc import Callable
from typing import NamedTuple

import equinox as eqx
import jax
import jax.numpy as jnp
from jax import lax, random
from jaxtyping import Array, Bool, Float32, Int32, Key

from bartz.mcmcloop._loop import _empty_trace, _set
from bartz.mcmcloop._trace import BurninTrace, MainTrace

from longbet._loop import (
    LongBetBurninTrace,
    LongBetTrace,
    _make_views,
    _set_param,
)
from longbet._multi_state import MultiLongBetState
from longbet._multi_step import multi_step
from longbet._tempering import multi_select_replicas, multi_swap_step


class MultiLongBetTrace(eqx.Module):
    """Posterior draws for coupled multi-outcome LongBet.

    Attributes
    ----------
    traces
        Tuple of child ``LongBetTrace`` objects in INTERNAL outcome order.
    gamma_loadings
        SUR regression loadings trace, shape ``(*chains, n_save, M, M)``.
        Strictly lower-triangular with zero diagonal.
    """

    traces: tuple[LongBetTrace, ...]
    gamma_loadings: Float32[Array, '...']


def _multi_views(state: MultiLongBetState, m: int):
    """Prediction-ready views of outcome ``m``'s two forests."""
    return _make_views(state.states[m])


class MultiLongBetBurninTrace(eqx.Module):
    """Diagnostic burn-in traces for coupled multi-outcome LongBet."""

    traces: tuple[LongBetBurninTrace, ...]


class RunMultiLongBetResult(NamedTuple):
    """Return value of :func:`run_multi_longbet_mcmc`."""

    final_state: MultiLongBetState
    main_trace: MultiLongBetTrace
    burnin_trace: MultiLongBetBurninTrace | None
    #: Per replica, accepted exchanges with the next level (tempered runs only).
    swap_accepts: Any = None


class _MultiMCMCCarry(eqx.Module):
    """Carry for the coupled MCMC ``lax.while_loop``."""

    state: MultiLongBetState
    key: Key[Array, '']
    i_total: Int32[Array, '']
    mu_burnins: tuple[BurninTrace, ...]
    nu_burnins: tuple[BurninTrace, ...]
    mu_mains: tuple[MainTrace, ...]
    nu_mains: tuple[MainTrace, ...]
    beta_traces: tuple[Float32[Array, '...'], ...]
    gamma_traces: tuple[Float32[Array, '...'], ...]
    b0_traces: tuple[Float32[Array, '...'], ...]
    b1_traces: tuple[Float32[Array, '...'], ...]
    alpha_traces: tuple[Float32[Array, '...'], ...]
    sigma2_traces: tuple[Float32[Array, '...'], ...]
    sigma_gamma2_traces: tuple[Float32[Array, '...'], ...]
    trend_coef_traces: tuple[Float32[Array, '...'], ...]
    cutpoint_traces: tuple[Float32[Array, '...'] | None, ...]
    gamma_loadings_trace: Float32[Array, '...']
    swap_accepts: Float32[Array, '...'] | None


def run_multi_longbet_mcmc(
    key: Key[Array, ''],
    state: MultiLongBetState,
    n_burn: int,
    n_save: int,
    n_skip: int = 1,
    inner_loop_length: int | None = None,
    callback: Callable[[int, int, MultiLongBetState], None] | None = None,
) -> RunMultiLongBetResult:
    return _run_multi_longbet_mcmc(key, state, n_burn, n_save, n_skip,
                                  inner_loop_length, callback)


def _run_multi_longbet_mcmc(
    key, state, n_burn, n_save, n_skip=1, inner_loop_length=None, callback=None,
) -> RunMultiLongBetResult:
    """Run the coupled multi-outcome LongBet sampler.

    Parameters
    ----------
    key
        PRNG master key.
    state
        Initial MultiLongBetState, single-chain or multi-chain.
    n_burn
        Burn-in sweeps to discard.
    n_save
        Posterior draws to retain per chain.
    n_skip
        Thinning interval.
    inner_loop_length
        Sweeps per outer dispatch.
    callback
        Optional host callable ``(outer_step, n_outer, state)`` called between batches.

    Returns
    -------
    RunMultiLongBetResult
    """
    n_burn_i = int(n_burn)
    n_save_i = int(n_save)
    n_skip_i = int(n_skip)
    n_iters = n_burn_i + n_skip_i * n_save_i

    M = state.M
    levels = int(state.states[0].tempering_levels)
    tempered = state.has_chain_axis and levels > 1
    cold_idx = jnp.arange(0, state.num_chains, levels) if tempered else None
    cold = (lambda st: multi_select_replicas(st, cold_idx)) if tempered else (lambda st: st)
    cold_init = cold(state)
    chains = cold_init.gamma_loadings.shape[:-2]
    sample_axis = 1 if chains else 0

    coupled = state.sur_active and any(state.continuous_mask[1:])
    if not coupled and not tempered and callback is None:
        from longbet import _loop
        from longbet._state import split_chain_fields
        from longbet._step import longbet_single_step

        final_states = []
        main_traces = []
        burnin_traces = []
        orig_step = _loop.longbet_step
        try:
            for m in range(M):
                child = state.states[m]

                def _outcome_step(k, st, _m=m):
                    k_outcome = random.fold_in(k, 0)
                    if not st.has_chain_axis:
                        return longbet_single_step(random.fold_in(k_outcome, _m), st)
                    keys = random.split(k_outcome, st.num_chains)
                    per, shared = split_chain_fields(st)

                    def _one(ck, p):
                        updated = longbet_single_step(
                            random.fold_in(ck, _m), eqx.combine(p, shared)
                        )
                        return split_chain_fields(updated)[0]

                    return eqx.combine(jax.vmap(_one)(keys, per), shared)

                _loop.longbet_step = _outcome_step
                res_m = _loop.run_longbet_mcmc(
                    key=key,
                    state=child,
                    n_burn=n_burn_i,
                    n_save=n_save_i,
                    n_skip=n_skip_i,
                    inner_loop_length=inner_loop_length,
                )
                final_states.append(res_m.final_state)
                main_traces.append(res_m.main_trace)
                burnin_traces.append(res_m.burnin_trace)
        finally:
            _loop.longbet_step = orig_step

        final_multi = eqx.tree_at(lambda s: s.states, state, tuple(final_states))
        return RunMultiLongBetResult(
            final_state=final_multi,
            main_trace=MultiLongBetTrace(
                traces=tuple(main_traces),
                gamma_loadings=jnp.zeros((*chains, n_save_i, M, M), jnp.float32),
            ),
            burnin_trace=(
                MultiLongBetBurninTrace(traces=tuple(burnin_traces))
                if n_burn_i > 0
                else None
            ),
            swap_accepts=None,
        )

    # Preallocate traces for each equation
    mu_burnins = []
    nu_burnins = []
    mu_mains = []
    nu_mains = []
    beta_traces = []
    gamma_traces = []
    b0_traces = []
    b1_traces = []
    alpha_traces = []
    sigma2_traces = []
    sigma_gamma2_traces = []
    trend_coef_traces = []
    cutpoint_traces = []

    for m in range(M):
        child = cold_init.states[m]
        v_mu, v_nu = _multi_views(cold_init, m)
        mu_burnins.append(_empty_trace(n_burn_i, v_mu, BurninTrace))
        nu_burnins.append(_empty_trace(n_burn_i, v_nu, BurninTrace))
        mu_mains.append(_empty_trace(n_save_i, v_mu, MainTrace))
        nu_mains.append(_empty_trace(n_save_i, v_nu, MainTrace))
        beta_traces.append(
            jnp.zeros((*chains, n_save_i, child.S_max + 1), jnp.float32)
        )
        gamma_traces.append(
            jnp.zeros((*chains, n_save_i, child.N_units), jnp.float32)
        )
        b0_traces.append(jnp.zeros((*chains, n_save_i), jnp.float32))
        b1_traces.append(jnp.zeros((*chains, n_save_i), jnp.float32))
        alpha_traces.append(jnp.zeros((*chains, n_save_i), jnp.float32))
        sigma2_traces.append(jnp.zeros((*chains, n_save_i), jnp.float32))
        sigma_gamma2_traces.append(jnp.zeros((*chains, n_save_i), jnp.float32))
        # Calendar-time trend coefficients, exactly as ``longbet._loop`` stores
        # them for the scalar model. When the block is off ``q_tot`` is zero and
        # this is a zero-width array rather than ``None``, so predict(),
        # save/load and the diagnostics need no multi-specific special case.
        trend_coef_traces.append(
            jnp.zeros((*chains, n_save_i, child.trend_design.shape[1]), jnp.float32)
        )
        cutpoint_traces.append(jnp.zeros((*chains, n_save_i, child.num_categories - 2), jnp.float32)
                               if child.outcome_type_str == "ordinal" else None)

    gamma_loadings_trace = jnp.zeros((*chains, n_save_i, M, M), jnp.float32)

    carry = _MultiMCMCCarry(
        state=state,
        key=key,
        i_total=jnp.int32(0),
        mu_burnins=tuple(mu_burnins),
        nu_burnins=tuple(nu_burnins),
        mu_mains=tuple(mu_mains),
        nu_mains=tuple(nu_mains),
        beta_traces=tuple(beta_traces),
        gamma_traces=tuple(gamma_traces),
        b0_traces=tuple(b0_traces),
        b1_traces=tuple(b1_traces),
        alpha_traces=tuple(alpha_traces),
        sigma2_traces=tuple(sigma2_traces),
        sigma_gamma2_traces=tuple(sigma_gamma2_traces),
        trend_coef_traces=tuple(trend_coef_traces),
        cutpoint_traces=tuple(cutpoint_traces),
        gamma_loadings_trace=gamma_loadings_trace,
        swap_accepts=jnp.zeros((state.num_chains,), jnp.float32) if tempered else None,
    )

    inner_length = (
        n_iters if inner_loop_length is None else max(1, int(inner_loop_length))
    )
    n_outer = max(1, -(-n_iters // inner_length))
    noop_idx = jnp.iinfo(jnp.int32).max

    def run_inner_batch(c: _MultiMCMCCarry, i_target: Int32[Array, '']) -> _MultiMCMCCarry:
        def cond_fn(carry_in: _MultiMCMCCarry) -> Bool[Array, '']:
            return carry_in.i_total < i_target

        def body_fn(carry_in: _MultiMCMCCarry) -> _MultiMCMCCarry:
            key_step, key_next = random.split(carry_in.key)
            key_swap = random.fold_in(key_step, 4242)  # keeps the untempered stream unchanged
            i = carry_in.i_total
            new_state = multi_step(key_step, carry_in.state, sweep_index=i)
            swap_accepts = carry_in.swap_accepts
            if tempered:
                new_state, accepted = multi_swap_step(key_swap, new_state, i)
                swap_accepts = swap_accepts + accepted.astype(jnp.float32)
            cold_state = cold(new_state)

            is_burn = i < n_burn_i
            burnin_idx = jnp.where(is_burn, i, noop_idx)

            i_from_burn = i - n_burn_i
            is_save = (~is_burn) & (((i_from_burn + 1) % n_skip_i) == 0)
            main_idx = jnp.where(is_save, i_from_burn // n_skip_i, noop_idx)

            new_mu_burnins = []
            new_nu_burnins = []
            new_mu_mains = []
            new_nu_mains = []
            new_beta_traces = []
            new_gamma_traces = []
            new_b0_traces = []
            new_b1_traces = []
            new_alpha_traces = []
            new_sigma2_traces = []
            new_sigma_gamma2_traces = []
            new_trend_coef_traces = []
            new_cutpoint_traces = []

            for m_idx in range(M):
                child_st = cold_state.states[m_idx]
                cp_trace = carry_in.cutpoint_traces[m_idx]
                new_cutpoint_traces.append(_set_param(cp_trace, main_idx, child_st.cutpoints, sample_axis)
                                           if cp_trace is not None else None)
                v_mu, v_nu = _multi_views(cold_state, m_idx)
                new_mu_burnins.append(
                    _set(carry_in.mu_burnins[m_idx], burnin_idx, BurninTrace.from_state(v_mu))
                )
                new_nu_burnins.append(
                    _set(carry_in.nu_burnins[m_idx], burnin_idx, BurninTrace.from_state(v_nu))
                )
                new_mu_mains.append(
                    _set(carry_in.mu_mains[m_idx], main_idx, MainTrace.from_state(v_mu))
                )
                new_nu_mains.append(
                    _set(carry_in.nu_mains[m_idx], main_idx, MainTrace.from_state(v_nu))
                )
                new_beta_traces.append(
                    _set_param(carry_in.beta_traces[m_idx], main_idx, child_st.beta, sample_axis)
                )
                new_gamma_traces.append(
                    _set_param(carry_in.gamma_traces[m_idx], main_idx, child_st.gamma, sample_axis)
                )
                new_b0_traces.append(
                    _set_param(carry_in.b0_traces[m_idx], main_idx, child_st.b0, sample_axis)
                )
                new_b1_traces.append(
                    _set_param(carry_in.b1_traces[m_idx], main_idx, child_st.b1, sample_axis)
                )
                new_alpha_traces.append(
                    _set_param(carry_in.alpha_traces[m_idx], main_idx, child_st.alpha, sample_axis)
                )
                new_sigma2_traces.append(
                    _set_param(carry_in.sigma2_traces[m_idx], main_idx, child_st.sigma2, sample_axis)
                )
                new_sigma_gamma2_traces.append(
                    _set_param(
                        carry_in.sigma_gamma2_traces[m_idx],
                        main_idx,
                        child_st.sigma_gamma2,
                        sample_axis,
                    )
                )
                new_trend_coef_traces.append(
                    _set_param(
                        carry_in.trend_coef_traces[m_idx],
                        main_idx,
                        child_st.trend_coef,
                        sample_axis,
                    )
                )

            new_gamma_loadings_trace = _set_param(
                carry_in.gamma_loadings_trace,
                main_idx,
                cold_state.gamma_loadings,
                sample_axis,
            )

            return _MultiMCMCCarry(
                state=new_state,
                key=key_next,
                i_total=i + 1,
                mu_burnins=tuple(new_mu_burnins),
                nu_burnins=tuple(new_nu_burnins),
                mu_mains=tuple(new_mu_mains),
                nu_mains=tuple(new_nu_mains),
                beta_traces=tuple(new_beta_traces),
                gamma_traces=tuple(new_gamma_traces),
                b0_traces=tuple(new_b0_traces),
                b1_traces=tuple(new_b1_traces),
                alpha_traces=tuple(new_alpha_traces),
                sigma2_traces=tuple(new_sigma2_traces),
                sigma_gamma2_traces=tuple(new_sigma_gamma2_traces),
                trend_coef_traces=tuple(new_trend_coef_traces),
                cutpoint_traces=tuple(new_cutpoint_traces),
                gamma_loadings_trace=new_gamma_loadings_trace,
                swap_accepts=swap_accepts,
            )

        return lax.while_loop(cond_fn, body_fn, c)

    jitted_inner = jax.jit(run_inner_batch)

    for outer_step in range(n_outer):
        i_target = min((outer_step + 1) * inner_length, n_iters)
        carry = jitted_inner(carry, jnp.int32(i_target))
        if callback is not None:
            callback(outer_step, n_outer, carry.state)

    child_traces = []
    child_burnin_traces = []
    for m in range(M):
        tr = LongBetTrace(
            mu_trace=carry.mu_mains[m],
            nu_trace=carry.nu_mains[m],
            beta=carry.beta_traces[m],
            gamma=carry.gamma_traces[m],
            b0=carry.b0_traces[m],
            b1=carry.b1_traces[m],
            alpha=carry.alpha_traces[m],
            sigma2=carry.sigma2_traces[m],
            sigma_gamma2=carry.sigma_gamma2_traces[m],
            trend_coef=carry.trend_coef_traces[m],
            cutpoints=carry.cutpoint_traces[m],
        )
        child_traces.append(tr)
        b_tr = LongBetBurninTrace(
            mu_trace=carry.mu_burnins[m],
            nu_trace=carry.nu_burnins[m],
        )
        child_burnin_traces.append(b_tr)

    main_trace = MultiLongBetTrace(
        traces=tuple(child_traces),
        gamma_loadings=carry.gamma_loadings_trace,
    )
    burnin_trace = (
        MultiLongBetBurninTrace(traces=tuple(child_burnin_traces))
        if n_burn_i > 0
        else None
    )

    return RunMultiLongBetResult(
        final_state=carry.state,
        main_trace=main_trace,
        burnin_trace=burnin_trace,
        swap_accepts=carry.swap_accepts)
