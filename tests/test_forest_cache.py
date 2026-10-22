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

"""Scalar sampler regressions for changing treatment weights and fixed settings."""

import dataclasses

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet import LongBetConfig
from longbet._forest_cache import refresh_prec_tree
from longbet._state import init_longbet
from longbet._step import longbet_single_step


def _constant_state(y, observed, variance=1.0, **overrides):
    """One constant mu leaf, with nu switched off by a fixed zero trajectory."""
    n = len(y)
    cfg = LongBetConfig(
        num_trees_pr=1, num_trees_trt=1,
        max_depth_pr=1, max_depth_trt=1,
        min_points_per_leaf_pr=1, min_points_per_leaf_trt=1,
        sample_beta=False, random_intercept=False,
        use_inter_ensemble_move=False, use_trend_block=False,
        standardize=False, num_chains=1,
    )
    cfg = dataclasses.replace(cfg, **overrides)
    state = init_longbet(
        X_unified=jnp.zeros((1, n), jnp.uint8),
        y=jnp.array(y, jnp.float32), unit_idx=jnp.zeros(n, jnp.int32),
        time_idx=jnp.arange(n), exposure_idx=jnp.zeros(n, jnp.int32),
        z_vec=jnp.ones(n), obs_mask=jnp.array(observed),
        max_split_mu=jnp.zeros(1, jnp.uint8),
        max_split_nu=jnp.zeros(1, jnp.uint8), config=cfg,
    )
    return eqx.tree_at(
        lambda s: (s.beta, s.sigma2), state,
        (jnp.zeros_like(state.beta), jnp.float32(variance)),
    )


def test_scalar_treatment_draw_refreshes_changing_precision():
    """Changing weights must preserve the analytical conditional variance.

    Conditional on the updated mu, a constant nu leaf has a known Gaussian
    law. Standardize its draws by that law while changing the fixed beta each
    iteration, including zero weights. A stale precision cache fails this test.
    """
    y = jnp.array([-1.0, 0.2, 1.4, 2.0, 0.4, -0.3])
    variance = 1.7
    initial = _constant_state(y, np.ones(len(y), bool), variance)
    amplitudes = jnp.array([0.0, 0.2, 0.75, 3.0])
    # Read the leaf prior off the state rather than assuming it: the treatment
    # forest carries the BCF total prior SD of 0.5, so one tree gives 1/0.5**2.
    leaf_prec = float(initial.forest_nu.leaf_prior_cov_inv)

    @jax.jit
    def run(state):
        def step(s, inputs):
            key, i = inputs
            w = amplitudes[i % len(amplitudes)]
            # This fixture has fixed b0=b1=1 and one exposure bin.
            s = eqx.tree_at(
                lambda s: (s.beta, s.sigma2, s.resid), s,
                (jnp.full_like(s.beta, w), jnp.float32(variance),
                 y - s.mu_fit - w * s.nu_fit),
            )
            updated = longbet_single_step(key, s)
            v = 1 / (leaf_prec + len(y) * w**2 / variance)
            mean = v * w * jnp.sum(y - updated.mu_fit) / variance
            score = (updated.nu_fit[0] - mean) / jnp.sqrt(v)
            return updated, score

        return jax.lax.scan(
            step, state,
            (jax.random.split(jax.random.key(882), 8000), jnp.arange(8000)),
        )

    _, scores = run(initial)
    for group in range(len(amplitudes)):
        z = np.asarray(scores)[group::len(amplitudes)]
        np.testing.assert_allclose(z.mean(), 0, atol=0.07)
        np.testing.assert_allclose(z.var(), 1, rtol=0.09)


def test_precision_cache_respects_pending_prunes_and_untouched_leaves():
    state = _constant_state(np.zeros(4), np.ones(4, bool))
    indices = jnp.array([[4, 5, 3, 3], [2, 3, 2, 3]], jnp.uint8)
    state = eqx.tree_at(
        lambda s: (s.forest.leaf_tree, s.forest.leaf_indices,
                   s.forest.to_prune, s.forest.move_node, s.prec_scale), state,
        (jnp.zeros((2, 8)), indices, jnp.array([True, False]),
         jnp.array([2, 1]), jnp.array([1.0, 2.0, 0.0, 4.0])),
        is_leaf=lambda x: x is None,
    )
    result = refresh_prec_tree(state)
    expected = np.zeros((2, 8))
    expected[0, [2, 3]] = [3, 4]
    expected[1, [2, 3]] = [1, 6]
    np.testing.assert_array_equal(result.forest.prec_tree, expected)
    np.testing.assert_array_equal(result.forest.leaf_indices, indices)
    np.testing.assert_array_equal(result.forest.to_prune, state.forest.to_prune)
    np.testing.assert_array_equal(result.forest.count_tree, state.forest.count_tree)


def test_precision_refresh_never_builds_dense_membership_arrays():
    """An all-leaf refresh must not allocate observations-by-leaf-slot arrays."""
    n, trees, slots = 2048, 3, 1024
    state = _constant_state(np.zeros(n), np.ones(n, bool))
    indices = jnp.tile(
        1 + jnp.arange(n, dtype=jnp.uint16) % (slots - 1), (trees, 1),
    )
    state = eqx.tree_at(
        lambda s: (s.forest.leaf_tree, s.forest.leaf_indices,
                   s.forest.to_prune, s.forest.move_node), state,
        (jnp.zeros((trees, slots)), indices,
         jnp.zeros(trees, bool), jnp.ones(trees, jnp.int32)),
    )
    traced = jax.make_jaxpr(lambda s: refresh_prec_tree(s).forest.prec_tree)(state)

    def shapes(value):
        if hasattr(value, "eqns"):
            for eqn in value.eqns:
                for var in eqn.outvars:
                    yield getattr(var.aval, "shape", ())
                yield from shapes(eqn.params)
        elif hasattr(value, "jaxpr"):
            yield from shapes(value.jaxpr)
        elif isinstance(value, dict):
            for child in value.values():
                yield from shapes(child)
        elif isinstance(value, (tuple, list)):
            for child in value:
                yield from shapes(child)

    bound = 8 * (trees * n + trees * slots * 16)
    for shape in shapes(traced):
        assert np.prod(shape) <= bound, f"Dense cache-refresh intermediate: {shape}"
    result = refresh_prec_tree(state)
    expected = np.bincount(np.asarray(indices[0]), minlength=slots)
    np.testing.assert_allclose(result.forest.prec_tree, np.tile(expected, (trees, 1)))


def test_fixed_beta_is_not_rescaled_by_ridge_move():
    state = _constant_state(np.zeros(4), np.ones(4, bool))
    state = eqx.tree_at(lambda s: s.beta, state, jnp.ones_like(state.beta))
    for i in range(20):
        state = longbet_single_step(jax.random.key(i), state)
        np.testing.assert_array_equal(state.beta, 1)


@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("num_chains", [1, 3])
def test_chain_initialization_preserves_fixed_parameters(binary, num_chains):
    cfg = LongBetConfig(
        num_trees_pr=2, num_trees_trt=2, num_chains=num_chains,
        sample_beta=False, random_intercept=False,
        outcome="binary" if binary else "continuous",
    )
    state = init_longbet(
        X_unified=jnp.zeros((1, 4), jnp.uint8),
        y=jnp.array([0.0, 1.0, 0.0, 1.0]),
        unit_idx=jnp.array([0, 0, 1, 1]), time_idx=jnp.array([0, 1, 0, 1]),
        exposure_idx=jnp.array([0, 1, 0, 0]), z_vec=jnp.array([0.0, 1.0, 0.0, 0.0]),
        obs_mask=jnp.ones(4, bool), max_split_mu=jnp.zeros(1, jnp.uint8),
        max_split_nu=jnp.zeros(1, jnp.uint8), config=cfg,
        num_chains=num_chains, chain_key=jax.random.key(43),
    )
    np.testing.assert_array_equal(state.beta, 1)
    np.testing.assert_array_equal(state.b0, 0)
    np.testing.assert_array_equal(state.b1, 1)
    np.testing.assert_array_equal(state.gamma, 0)
    if binary:
        np.testing.assert_array_equal(state.sigma2, 1)
