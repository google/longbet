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

"""Unit conversions, fit-by-differencing, and the full-model residual invariant.

``bartz`` stores its residual scaled (``resid_unit * resid == y - sum(trees)``)
and ``resid_unit`` differs between the two forests. LongBet carries the
full-model residual in data units and converts only at the boundary of a bartz
call. These tests pin both halves of that contract.
"""

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from bartz.grove import evaluate_forest
from bartz.mcmcstep import Wishart, init, make_p_nonterminal, step

from longbet._config import LongBetConfig
from longbet._state import init_longbet
from longbet._step import _load, _read, longbet_step


def _bartz_state(n=30, num_trees=6, seed=1):
    return init(
        X=jnp.zeros((2, n), dtype=jnp.uint8),
        y=jax.random.normal(jax.random.key(seed), (n,)),
        outcome_type='continuous',
        offset=0.0,
        max_split=jnp.array([1, 1], dtype=jnp.uint8),
        num_trees=num_trees,
        p_nonterminal=make_p_nonterminal(4, 0.95, 2.0),
        leaf_prior_cov_inv=1.0,
        error_cov_inv=Wishart(nu=jnp.array(3.0), rate=jnp.array(3.0), value=jnp.array(1.0)),
    )


def test_load_read_round_trip():
    state = _bartz_state()
    r = jax.random.normal(jax.random.key(2), (state.y.shape[0],))
    assert jnp.allclose(_read(_load(r, state.resid_unit), state.resid_unit), r, atol=1e-6)


def test_resid_unit_is_not_one():
    """If it were, a missing conversion would pass unnoticed."""
    assert float(_bartz_state().resid_unit) != 1.0


def test_fit_by_differencing_equals_evaluate_forest():
    """Recovering a forest's fit by differencing the residual must be exact.

    LongBet never re-evaluates the trees to get a fit; it differences the
    residual across the bartz call, which is only valid because bartz applies
    its leaf delta to every observation, masked ones included.
    """
    state = _bartz_state(seed=10)
    before_fit = evaluate_forest(state.X, state.forest, sum_batch_axis=-1)
    before_resid = _read(state.resid, state.resid_unit)

    state_next = step(jax.random.key(42), state)
    after_fit = evaluate_forest(state_next.X, state_next.forest, sum_batch_axis=-1)
    after_resid = _read(state_next.resid, state_next.resid_unit)

    np.testing.assert_allclose(
        np.asarray(after_fit - before_fit),
        np.asarray(before_resid - after_resid),
        atol=1e-5,
    )


@pytest.mark.parametrize("binary", [False, True])
def test_full_model_residual_invariant(binary):
    """After every sweep, ``resid == y - fitted`` on observed cells.

    This is the single strongest end-to-end check on the sweep's bookkeeping: it
    fails on a wrong ``resid_unit`` conversion, a missed residual update after
    any of the seven blocks, or a mis-scaled treatment weight. Missing cells are
    included in the panel so the masked path is exercised too.
    """
    N, T, P = 12, 5, 3
    M = N * T
    rng = np.random.default_rng(0)
    obs = np.ones(M, dtype=bool)
    obs[rng.choice(M, 8, replace=False)] = False
    exposure = np.tile(np.clip(np.arange(T) - 1, 0, None), N).astype(np.int32)

    state = init_longbet(
        X_unified=jnp.asarray(rng.integers(0, 5, (P, M), dtype=np.uint8)),
        y=jnp.asarray(rng.binomial(1, 0.3, size=M) if binary else rng.normal(size=M), jnp.float32),
        unit_idx=jnp.repeat(jnp.arange(N), T),
        time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.asarray(exposure),
        z_vec=jnp.asarray(np.tile((np.arange(T) >= 2).astype(np.float32), N)),
        obs_mask=jnp.asarray(obs),
        max_split_mu=jnp.full(P, 4, jnp.uint8),
        max_split_nu=jnp.full(P, 4, jnp.uint8),
        config=LongBetConfig(num_trees_pr=4, num_trees_trt=4,
                             outcome="binary" if binary else "continuous"),
        offset=-0.7 if binary else 0.0,
    )

    key = jax.random.key(0)
    for sweep in range(6):
        key, sub = jax.random.split(key)
        state = longbet_step(sub, state)

        b_z = np.where(np.asarray(state.z_vec) == 1.0, float(state.b1), float(state.b0))
        fitted = (
            float(state.alpha) * np.asarray(state.mu_fit)
            + b_z
            * np.asarray(state.beta)[np.asarray(state.exposure_idx)]
            * np.asarray(state.nu_fit)
            + np.asarray(state.gamma)[np.asarray(state.unit_idx)]
        )
        # This independent forest evaluation catches latent-response changes
        # incorrectly being counted as tree-fit deltas.
        np.testing.assert_allclose(
            state.mu_fit,
            evaluate_forest(state.X, state.forest, sum_batch_axis=-1) + (-0.7 if binary else 0.0), atol=1e-4,
        )
        np.testing.assert_allclose(
            state.nu_fit,
            evaluate_forest(state.X, state.forest_nu, sum_batch_axis=-1), atol=1e-4,
        )
        implied = np.asarray(state.z if binary else state.y) - fitted
        err = np.abs(np.where(obs, implied - np.asarray(state.resid), 0.0)).max()
        assert err < 1e-4, f"residual bookkeeping drifted at sweep {sweep}: {err:.2e}"

    # Masked cells carry a zero residual so nothing non-finite propagates.
    assert np.allclose(np.asarray(state.resid)[~obs], 0.0)
