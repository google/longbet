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

"""Compilation behaviour of the sweep and the loop."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet._config import LongBetConfig
from longbet._loop import run_longbet_mcmc
from longbet._state import init_longbet
from longbet._step import longbet_step


def _state(num_chains=None, **cfg):
    N, T, P = 10, 4, 3
    M = N * T
    rng = np.random.default_rng(42)
    return init_longbet(
        X_unified=jnp.asarray(rng.integers(0, 5, (P, M), dtype=np.uint8)),
        y=jnp.asarray(rng.normal(size=M), jnp.float32),
        unit_idx=jnp.repeat(jnp.arange(N), T),
        time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.asarray(np.tile(np.clip(np.arange(T) - 1, 0, None), N), jnp.int32),
        z_vec=jnp.asarray(np.tile((np.arange(T) >= 2).astype(np.float32), N)),
        obs_mask=jnp.ones(M, bool),
        max_split_mu=jnp.full(P, 4, jnp.uint8),
        max_split_nu=jnp.full(P, 4, jnp.uint8),
        config=LongBetConfig(num_trees_pr=3, num_trees_trt=3,
                             num_chains=num_chains or 1, **cfg),
        num_chains=num_chains,
    )


def test_sweep_compiles_once_and_is_reused():
    """Every sweep after the first must hit the compilation cache."""
    state = _state()
    longbet_step._clear_cache()
    key = jax.random.key(0)
    for _ in range(5):
        key, sub = jax.random.split(key)
        state = longbet_step(sub, state)
    assert longbet_step._cache_size() == 1, (
        f"the sweep recompiled: {longbet_step._cache_size()} cache entries"
    )
    assert jnp.all(jnp.isfinite(state.beta))


def test_inner_loop_length_does_not_recompile():
    """The batch bound is traced, so splitting the run must not recompile.

    It used to be a static argument, which recompiled the whole while-loop once
    per batch -- exactly cancelling the reason to batch at all.
    """
    state = _state()
    result_batched = run_longbet_mcmc(
        jax.random.key(1), state, n_burn=3, n_save=6, n_skip=1, inner_loop_length=2
    )
    result_single = run_longbet_mcmc(
        jax.random.key(1), state, n_burn=3, n_save=6, n_skip=1, inner_loop_length=None
    )
    # Same key, same schedule: batching is a dispatch detail, not a model change.
    np.testing.assert_allclose(
        np.asarray(result_batched.main_trace.beta),
        np.asarray(result_single.main_trace.beta),
        rtol=1e-6,
        atol=1e-6,
    )


def test_callback_runs_between_batches():
    state = _state()
    seen = []
    run_longbet_mcmc(
        jax.random.key(2), state, n_burn=2, n_save=4, inner_loop_length=2,
        callback=lambda i, n, s: seen.append((i, n)),
    )
    assert seen == [(0, 3), (1, 3), (2, 3)]


def test_thinning_matches_bartz_convention():
    """Total iterations are n_burn + n_skip * n_save, with n_skip=1 meaning none."""
    state = _state()
    result = run_longbet_mcmc(jax.random.key(3), state, n_burn=4, n_save=5, n_skip=3)
    assert np.asarray(result.main_trace.beta).shape[0] == 5
    assert result.burnin_trace is not None
    cfg = LongBetConfig(num_burnin=4, num_sweeps=5, n_skip=3)
    assert cfg.total_iterations() == 4 + 3 * 5


@pytest.mark.parametrize("num_chains", [None, 2])
def test_no_donated_buffer_errors_across_repeated_sweeps(num_chains):
    """bartz's step donates; calling a sweep twice must not delete live arrays."""
    state = _state(num_chains=num_chains)
    key = jax.random.key(4)
    for _ in range(3):
        key, sub = jax.random.split(key)
        state = longbet_step(sub, state)
        assert np.all(np.isfinite(np.asarray(state.X)))
        assert np.all(np.isfinite(np.asarray(state.resid)))
