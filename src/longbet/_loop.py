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

"""MCMC driver for LongBet."""

from __future__ import annotations

from collections.abc import Callable
from functools import partial
from typing import Any, NamedTuple

import equinox as eqx
import jax
import jax.numpy as jnp
from jax import random
from jaxtyping import Array, Float32, Int32, Key

from bartz.mcmcloop._loop import _empty_trace, add_dummy_axis, trace_sample_axes
from bartz.mcmcloop._trace import BurninTrace, MainTrace
from bartz.mcmcstep import State

from longbet._state import LongBetState
from longbet._step import _load, enable_x64, longbet_step
from longbet._tempering import select_replicas, swap_step

_DEFAULT_LONGBET_STEP = longbet_step


@partial(jax.jit, static_argnames=('custom_step_fn',))
@enable_x64()
def _custom_step_dispatch(
    custom_step_fn: Callable[[Key[Array, ''], LongBetState], LongBetState],
    key: Key[Array, ''],
    state: LongBetState,
) -> LongBetState:
    return custom_step_fn(key, state)


class LongBetTrace(eqx.Module):
    """Posterior draws saved during the main sampling phase.

    Parameter arrays are shaped ``(*chains, n_save, ...)``; the forest traces
    follow ``bartz``'s own layout, which is the same convention.
    """

    mu_trace: MainTrace
    nu_trace: MainTrace
    beta: Float32[Array, '...']
    gamma: Float32[Array, '...']
    b0: Float32[Array, '...']
    b1: Float32[Array, '...']
    alpha: Float32[Array, '...']
    sigma2: Float32[Array, '...']
    sigma_gamma2: Float32[Array, '...']
    trend_coef: Float32[Array, '...'] | None = None
    cutpoints: Float32[Array, '...'] | None = None
    #: encouragement-v2 extensions, ``None`` when off: the moderator's linear
    #: coefficients ``(*chains, n_save, p)`` and scale, and the trajectory scale.
    delta: Float32[Array, '...'] | None = None
    tau_delta2: Float32[Array, '...'] | None = None
    beta_scale2: Float32[Array, '...'] | None = None
    ar1_phi: Float32[Array, '...'] | None = None
    #: Distributed-lag kernel (``exposure_lags``): ``kappa`` ``(*chains,
    #: n_save, L + 1)`` and its learned scale; and the coefficients of the
    #: fixed extra columns ``(*chains, n_save, Q)``. ``None`` when absent.
    kappa: Float32[Array, '...'] | None = None
    lag_scale2: Float32[Array, '...'] | None = None
    extra_coef: Float32[Array, '...'] | None = None


class LongBetBurninTrace(eqx.Module):
    """Diagnostic traces saved during burn-in."""

    mu_trace: BurninTrace
    nu_trace: BurninTrace


class RunLongBetResult(NamedTuple):
    """Return value of :func:`run_longbet_mcmc`."""

    final_state: LongBetState
    main_trace: LongBetTrace
    burnin_trace: LongBetBurninTrace | None
    #: Per replica, accepted exchanges with the next level (tempered runs only).
    swap_accepts: Any = None


def _make_views(state: LongBetState) -> tuple[State, State]:
    """Build plain ``State`` views of the two forests, for the trace machinery.

    Only the forest, the error precision and the move counters are read by
    ``MainTrace.from_state`` / ``BurninTrace.from_state``, but the residuals are
    converted into each forest's own units anyway so the views are internally
    consistent if upstream ever starts reading them.
    """
    view_mu = State(
        _chain_anchor=state._chain_anchor,
        X=state.X,
        y=state.y,
        z=state.z,
        binary_indices=state.binary_indices,
        resid=_load(state.resid, state.resid_unit),
        resid_unit=state.resid_unit,
        resid_eff_scale=state.resid_eff_scale,
        resid_inexact_integral=state.resid_inexact_integral,
        error_cov_inv=state.error_cov_inv,
        error_scale=state.error_scale,
        prec_scale=state.prec_scale,
        inv_sdev_scale=state.inv_sdev_scale,
        inv_sdev_unit=state.inv_sdev_unit,
        n_non_missing=state.n_non_missing,
        sum_diag_prec_scale=state.sum_diag_prec_scale,
        forest=state.forest,
        config=state.config,
    )
    view_nu = State(
        _chain_anchor=state._chain_anchor,
        X=state.X,
        y=state.y,
        z=None,
        binary_indices=None,
        resid=state.resid_nu,
        resid_unit=state.resid_unit_nu,
        resid_eff_scale=state.resid_eff_scale_nu,
        resid_inexact_integral=state.resid_inexact_integral_nu,
        error_cov_inv=state.error_cov_inv,
        error_scale=None,
        prec_scale=state.prec_scale_nu,
        inv_sdev_scale=state.inv_sdev_scale_nu,
        inv_sdev_unit=state.inv_sdev_unit_nu,
        n_non_missing=state.n_non_missing,
        sum_diag_prec_scale=state.sum_diag_prec_scale,
        forest=state.forest_nu,
        config=state.config,
    )
    return view_mu, view_nu


def _set_param(
    trace: Float32[Array, '...'],
    index: Int32[Array, ''],
    val: Float32[Array, '...'],
    sample_axis: int,
) -> Float32[Array, '...']:
    """Write ``val`` at ``index`` along ``sample_axis``; out-of-range is dropped."""
    ndindex = (slice(None),) * sample_axis + (index, ...)
    return trace.at[ndindex].set(val, mode='drop')


@jax.jit
def _extract_burnin(state: LongBetState) -> tuple[BurninTrace, BurninTrace]:
    v_mu, v_nu = _make_views(state)
    return BurninTrace.from_state(v_mu), BurninTrace.from_state(v_nu)


@jax.jit
def _extract_main(state: LongBetState) -> tuple[MainTrace, MainTrace]:
    v_mu, v_nu = _make_views(state)
    return MainTrace.from_state(v_mu), MainTrace.from_state(v_nu)


@partial(jax.jit, static_argnames=('levels',))
def _tempered_post_step(
    key_step: Key[Array, ''],
    state: LongBetState,
    i: Int32[Array, ''],
    swap_accepts: Float32[Array, '...'],
    levels: int,
) -> tuple[LongBetState, Float32[Array, '...'], LongBetState]:
    key_swap = random.fold_in(key_step, 4242)  # keeps the untempered stream unchanged
    new_state, accepted = swap_step(key_swap, state, i)
    swap_accepts = swap_accepts + accepted.astype(jnp.float32)
    cold_idx = jnp.arange(0, new_state.num_chains, levels)
    cold_state = select_replicas(new_state, cold_idx)
    return new_state, swap_accepts, cold_state


def _stack_forest_traces(
    items: list[Any],
    length: int,
    view_init: State,
    trace_cls: type[BurninTrace] | type[MainTrace],
) -> Any:
    if length == 0 or not items:
        return _empty_trace(0, view_init, trace_cls)
    axes = trace_sample_axes(add_dummy_axis(items[0]))

    def _combine(ax: int | None, *xs: Any) -> Any:
        return xs[0] if ax is None else jnp.stack(xs, axis=ax)

    return jax.tree.map(_combine, axes, *items, is_leaf=lambda x: x is None)


def _stack_param(
    items: list[Float32[Array, '...']],
    sample_axis: int,
    empty_shape: tuple[int, ...],
) -> Float32[Array, '...']:
    if not items:
        return jnp.zeros(empty_shape, jnp.float32)
    return jnp.stack(items, axis=sample_axis)


def run_longbet_mcmc(
    key: Key[Array, ''],
    state: LongBetState,
    n_burn: int,
    n_save: int,
    n_skip: int = 1,
    inner_loop_length: int | None = None,
    callback: Callable[[int, int, LongBetState], None] | None = None,
) -> RunLongBetResult:
    """Run the LongBet sampler.

    The thinning convention matches ``bartz`` exactly: the total number of
    iterations is ``n_burn + n_skip * n_save``, and ``n_skip = 1`` means no
    thinning.

    Parameters
    ----------
    key
        PRNG key.
    state
        Initial state, possibly carrying a chain axis.
    n_burn
        Burn-in sweeps, discarded.
    n_save
        Posterior draws to retain.
    n_skip
        Thinning interval.
    inner_loop_length
        Sweeps per outer dispatch. Bounds compile time and lets ``callback``
        run between batches. The loop bound is passed as a traced value, so
        varying it does **not** trigger recompilation.
    callback
        Optional host callable ``(outer_step, n_outer, state)`` run between
        batches.
    """
    n_burn_i, n_save_i, n_skip_i = int(n_burn), int(n_save), int(n_skip)
    n_iters = n_burn_i + n_skip_i * n_save_i

    # Tempered runs carry every replica in the state but save only the
    # posterior (beta = 1) replicas, so memory does not scale with the ladder.
    levels = int(state.tempering_levels)
    tempered = state.has_chain_axis and levels > 1
    cold_idx = jnp.arange(0, state.num_chains, levels) if tempered else None
    cold = (lambda st: select_replicas(st, cold_idx)) if tempered else (lambda st: st)
    chains = cold(state).beta.shape[:-1]
    sample_axis = 1 if chains else 0

    view_mu_init, view_nu_init = _make_views(cold(state))
    swap_accepts = jnp.zeros((state.num_chains,), jnp.float32) if tempered else None
    has_cutpoints = state.outcome_type_str == "ordinal"

    mu_burnin_list: list[BurninTrace] = []
    nu_burnin_list: list[BurninTrace] = []
    mu_main_list: list[MainTrace] = []
    nu_main_list: list[MainTrace] = []
    beta_list: list[Float32[Array, '...']] = []
    gamma_list: list[Float32[Array, '...']] = []
    b0_list: list[Float32[Array, '...']] = []
    b1_list: list[Float32[Array, '...']] = []
    alpha_list: list[Float32[Array, '...']] = []
    sigma2_list: list[Float32[Array, '...']] = []
    sigma_gamma2_list: list[Float32[Array, '...']] = []
    trend_coef_list: list[Float32[Array, '...']] = []
    cutpoints_list: list[Float32[Array, '...']] = []
    delta_list: list[Float32[Array, '...']] = []
    tau_delta2_list: list[Float32[Array, '...']] = []
    beta_scale2_list: list[Float32[Array, '...']] = []
    ar1_phi_list: list[Float32[Array, '...']] = []
    kappa_list: list[Float32[Array, '...']] = []
    lag_scale2_list: list[Float32[Array, '...']] = []
    extra_coef_list: list[Float32[Array, '...']] = []

    step_fn = (
        longbet_step
        if longbet_step is _DEFAULT_LONGBET_STEP
        else partial(_custom_step_dispatch, longbet_step)
    )

    inner_length = (
        max(1, n_iters) if inner_loop_length is None else max(1, int(inner_loop_length))
    )
    n_outer = max(1, -(-n_iters // inner_length))

    cur_state = state
    cur_key = key
    for outer_step in range(n_outer):
        i_start = outer_step * inner_length
        i_target = min((outer_step + 1) * inner_length, n_iters)
        for i in range(i_start, i_target):
            key_step, cur_key = random.split(cur_key)
            cur_state = step_fn(key_step, cur_state)
            if tempered:
                assert swap_accepts is not None
                cur_state, swap_accepts, cold_state = _tempered_post_step(
                    key_step, cur_state, jnp.int32(i), swap_accepts, levels
                )
            else:
                cold_state = cur_state

            if i < n_burn_i:
                b_mu, b_nu = _extract_burnin(cold_state)
                mu_burnin_list.append(b_mu)
                nu_burnin_list.append(b_nu)
            elif ((i - n_burn_i + 1) % n_skip_i) == 0:
                m_mu, m_nu = _extract_main(cold_state)
                mu_main_list.append(m_mu)
                nu_main_list.append(m_nu)
                beta_list.append(cold_state.beta)
                gamma_list.append(cold_state.gamma)
                b0_list.append(cold_state.b0)
                b1_list.append(cold_state.b1)
                alpha_list.append(cold_state.alpha)
                sigma2_list.append(cold_state.sigma2)
                sigma_gamma2_list.append(cold_state.sigma_gamma2)
                trend_coef_list.append(cold_state.trend_coef)
                if has_cutpoints:
                    cutpoints_list.append(cold_state.cutpoints)
                if state.use_linear_trt:
                    delta_list.append(cold_state.delta)
                    tau_delta2_list.append(cold_state.tau_delta2)
                if state.learn_beta_scale:
                    beta_scale2_list.append(cold_state.beta_scale2)
                if state.ar1:
                    ar1_phi_list.append(cold_state.ar1_phi)
                if state.use_lag_kernel:
                    kappa_list.append(cold_state.kappa)
                    if state.lag_learn_scale:
                        lag_scale2_list.append(cold_state.lag_scale2)
                if state.num_extra_cols:
                    extra_coef_list.append(cold_state.extra_coef)
        if callback is not None:
            callback(outer_step, n_outer, cur_state)

    main_trace = LongBetTrace(
        mu_trace=_stack_forest_traces(mu_main_list, n_save_i, view_mu_init, MainTrace),
        nu_trace=_stack_forest_traces(nu_main_list, n_save_i, view_nu_init, MainTrace),
        beta=_stack_param(beta_list, sample_axis, (*chains, n_save_i, state.S_max + 1)),
        gamma=_stack_param(gamma_list, sample_axis, (*chains, n_save_i, state.N_units)),
        b0=_stack_param(b0_list, sample_axis, (*chains, n_save_i)),
        b1=_stack_param(b1_list, sample_axis, (*chains, n_save_i)),
        alpha=_stack_param(alpha_list, sample_axis, (*chains, n_save_i)),
        sigma2=_stack_param(sigma2_list, sample_axis, (*chains, n_save_i)),
        sigma_gamma2=_stack_param(sigma_gamma2_list, sample_axis, (*chains, n_save_i)),
        trend_coef=_stack_param(
            trend_coef_list, sample_axis, (*chains, n_save_i, state.trend_design.shape[1])
        ),
        cutpoints=(
            _stack_param(
                cutpoints_list, sample_axis, (*chains, n_save_i, state.num_categories - 2)
            )
            if has_cutpoints
            else None
        ),
        delta=(
            _stack_param(delta_list, sample_axis, (*chains, n_save_i, state.lin_design.shape[1]))
            if state.use_linear_trt
            else None
        ),
        tau_delta2=(
            _stack_param(tau_delta2_list, sample_axis, (*chains, n_save_i))
            if state.use_linear_trt
            else None
        ),
        beta_scale2=(
            _stack_param(beta_scale2_list, sample_axis, (*chains, n_save_i))
            if state.learn_beta_scale
            else None
        ),
        ar1_phi=(
            _stack_param(ar1_phi_list, sample_axis, (*chains, n_save_i))
            if state.ar1
            else None
        ),
        kappa=(
            _stack_param(kappa_list, sample_axis, (*chains, n_save_i, state.lag_design.shape[1]))
            if state.use_lag_kernel
            else None
        ),
        lag_scale2=(
            _stack_param(lag_scale2_list, sample_axis, (*chains, n_save_i))
            if state.use_lag_kernel and state.lag_learn_scale
            else None
        ),
        extra_coef=(
            _stack_param(extra_coef_list, sample_axis, (*chains, n_save_i, state.num_extra_cols))
            if state.num_extra_cols
            else None
        ),
    )
    burnin_trace = (
        LongBetBurninTrace(
            mu_trace=_stack_forest_traces(mu_burnin_list, n_burn_i, view_mu_init, BurninTrace),
            nu_trace=_stack_forest_traces(nu_burnin_list, n_burn_i, view_nu_init, BurninTrace),
        )
        if n_burn_i > 0
        else None
    )
    return RunLongBetResult(cur_state, main_trace, burnin_trace, swap_accepts)
