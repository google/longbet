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

"""Tests for multi-outcome single step, residual protocol, and binary invariance."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet import LongBetConfig
from longbet._multi_input import normalize_multi_inputs
from longbet._multi_state import init_multi_longbet
from longbet._multi_step import multi_step


def test_residual_protocol_and_binary_sigma2():
    """Verify that stored residuals are raw, adjustments are temporary,

    and binary sigma2 is held at 1.0 throughout the step.
    """
    N, T = 6, 4
    rng = np.random.default_rng(20)

    # 1 binary outcome, 1 continuous outcome
    y_bin = rng.choice([0.0, 1.0], size=(N, T)).astype(np.float32)
    y_cont = rng.standard_normal((N, T)).astype(np.float32)

    X = rng.standard_normal((N, 3)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[:, 2:] = 1.0  # absorbing

    cfg = LongBetConfig(
        num_trees_pr=5,
        num_trees_trt=3,  # unequal tree counts
        sur=True,
        sur_prior_var=1.0,
        random_intercept=True,
    )

    norm = normalize_multi_inputs(
        {"bin": y_bin, "cont": y_cont},
        outcome={"bin": "binary", "cont": "continuous"},
        outcome_names=None,
        config=cfg,
    )

    # Order: bin is internal 0, cont is internal 1
    assert norm.internal_outcomes == ("binary", "continuous")

    P = 3
    X_unified = jnp.zeros((P, N * T), dtype=jnp.uint8)
    unit_idx = jnp.repeat(jnp.arange(N), T)
    time_idx = jnp.tile(jnp.arange(T), N)
    exposure_idx = jnp.zeros(N * T, dtype=jnp.int32)
    z_vec = jnp.asarray(z.ravel(order="C"), dtype=jnp.float32)
    max_split = jnp.full(P, 10, dtype=jnp.uint8)

    key = jax.random.key(42)
    k_init, k_step = jax.random.split(key)

    state = init_multi_longbet(
        X_unified=X_unified,
        unit_idx=unit_idx,
        time_idx=time_idx,
        exposure_idx=exposure_idx,
        z_vec=z_vec,
        max_split_mu=max_split,
        max_split_nu=max_split,
        norm_input=norm,
        config=cfg,
        num_chains=1,
        init_key=k_init,
    )

    # Invariant: Binary sigma2 must be exactly 1.0 before sweep 0
    assert float(state.states[0].sigma2) == 1.0

    # Execute one step at sweep_index = 0
    k_step0, k_step1 = jax.random.split(k_step)
    step0_out = multi_step(k_step0, state, sweep_index=jnp.int32(0))

    # At sweep_index 0, Gamma remains 0
    np.testing.assert_allclose(step0_out.gamma_loadings, 0.0)
    # Binary sigma2 remains 1.0
    assert float(step0_out.states[0].sigma2) == 1.0
    # Binary row in Gamma stays 0
    np.testing.assert_allclose(step0_out.gamma_loadings[0, :], 0.0)

    # Execute second step at sweep_index = 1
    step1_out = multi_step(k_step1, step0_out, sweep_index=jnp.int32(1))

    # Binary row 0 stays 0
    np.testing.assert_allclose(step1_out.gamma_loadings[0, :], 0.0)
    # Diagonal stays 0
    assert float(step1_out.gamma_loadings[1, 1]) == 0.0
    # Upper triangle stays 0
    assert float(step1_out.gamma_loadings[0, 1]) == 0.0
    # Binary sigma2 remains 1.0
    assert float(step1_out.states[0].sigma2) == 1.0

    # Stored residuals in step1_out must be finite raw residuals
    assert np.all(np.isfinite(step1_out.states[0].resid))
    assert np.all(np.isfinite(step1_out.states[1].resid))

    # Check raw residuals against independently evaluated forests through
    # several NONZERO-loading sweeps, not merely that residuals are finite.
    from bartz.grove import evaluate_forest
    state = step1_out
    for sweep in range(2, 6):
        sweep_key = jax.random.fold_in(k_step, sweep)
        state = multi_step(sweep_key, state, jnp.int32(sweep))
        for child in state.states:
            mu = evaluate_forest(child.X, child.forest, sum_batch_axis=-1) + child.forest.offset
            nu = evaluate_forest(child.X, child.forest_nu, sum_batch_axis=-1) + child.forest_nu.offset
            fitted = (child.alpha * mu + jnp.where(child.z_vec == 1, child.b1, child.b0)
                      * child.beta[child.exposure_idx] * nu + child.gamma[child.unit_idx])
            response = child.z if child.outcome_type_str == "binary" else child.y
            expected = jnp.where(child.obs_mask, response - fitted, 0)
            np.testing.assert_allclose(child.resid, expected, atol=2e-4)

        # Structural variance uses e_m = r_m - Gamma[m] @ r, never the
        # conditional pseudo-residual passed temporarily to its forest step.
        from longbet._step import _sample_inv_gamma
        child = state.states[1]
        innovation = child.resid - state.gamma_loadings[1] @ jnp.stack(
            [s.resid for s in state.states])
        outcome_key = jax.random.fold_in(jax.random.fold_in(sweep_key, 0), 1)
        expected_sigma2 = _sample_inv_gamma(
            jax.random.split(outcome_key, 10)[9],
            child.sigma_prior_a + jnp.sum(child.obs_mask)/2,
            child.sigma_prior_b + .5*jnp.sum(jnp.where(child.obs_mask, innovation, 0)**2),
        )
        np.testing.assert_allclose(child.sigma2, expected_sigma2, rtol=1e-6)
