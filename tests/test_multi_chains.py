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

"""Tests for multi-chain vectorization, key schedule, and SUR toggling invariants."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, LongBetMulti
from longbet._multi_input import normalize_multi_inputs
from longbet._multi_state import broadcast_multi_to_chains, init_multi_longbet
from longbet._multi_step import multi_single_step, multi_step


def test_vmap_equals_single_chain_step():
    """Verify vmap across chains gives identical transitions to sequential single chain calls."""
    N, T, P = 6, 4, 2
    Q = N * T
    rng = np.random.default_rng(101)

    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    cfg = LongBetConfig(
        num_trees_pr=3,
        num_trees_trt=3,
        num_chains=2,
        sur=True,
        sur_prior_var=1.0,
    )
    norm = normalize_multi_inputs(
        {"y1": y1, "y2": y2},
        outcome="continuous",
        outcome_names=None,
        config=cfg,
    )

    X_unified = jnp.asarray(rng.integers(0, 4, (P, Q), dtype=np.uint8))
    unit_idx = jnp.repeat(jnp.arange(N), T)
    time_idx = jnp.tile(jnp.arange(T), N)
    exposure_idx = jnp.zeros(Q, dtype=jnp.int32)
    z_vec = jnp.zeros(Q, dtype=jnp.float32)
    max_split = jnp.full(P, 5, dtype=jnp.uint8)

    key = jax.random.key(123)
    k_init, k_step = jax.random.split(key)

    # Initial state with 2 chains
    state2 = init_multi_longbet(
        X_unified=X_unified,
        unit_idx=unit_idx,
        time_idx=time_idx,
        exposure_idx=exposure_idx,
        z_vec=z_vec,
        max_split_mu=max_split,
        max_split_nu=max_split,
        norm_input=norm,
        config=cfg,
        num_chains=2,
        init_key=k_init,
    )

    assert state2.has_chain_axis
    assert state2.gamma_loadings.shape == (2, 2, 2)

    # Step via vmap (multi_step)
    stepped_vmap = multi_step(k_step, state2, sweep_index=jnp.int32(1))

    # Reconstruct single chains and step individually
    # Outcome and loading keys as derived in multi_step
    k_outcome = jax.random.fold_in(k_step, 0)
    k_loading = jax.random.fold_in(k_step, 1)
    outcome_chains = jax.random.split(k_outcome, 2)
    loading_chains = jax.random.split(k_loading, 2)

    import equinox as eqx
    from longbet._multi_state import split_multi_chain_fields

    per_chain, shared = split_multi_chain_fields(state2)
    per_vmap, _ = split_multi_chain_fields(stepped_vmap)

    for c in range(2):
        per_c = jax.tree.map(lambda x: x[c], per_chain)
        single_c = eqx.combine(per_c, shared)

        out_keys_c = jnp.stack([jax.random.fold_in(outcome_chains[c], m) for m in range(2)])
        load_keys_c = jnp.stack([jax.random.fold_in(loading_chains[c], m) for m in range(2)])

        stepped_c = multi_single_step(out_keys_c, load_keys_c, single_c, jnp.int32(1))
        per_stepped_c, _ = split_multi_chain_fields(stepped_c)
        per_vmap_c = jax.tree.map(lambda x: x[c], per_vmap)

        for leaf_vmap, leaf_c in zip(jax.tree.leaves(per_vmap_c), jax.tree.leaves(per_stepped_c)):
            np.testing.assert_allclose(leaf_vmap, leaf_c, rtol=1e-5, atol=1e-5)


def test_sur_feedback_changes_binary_draws_but_not_identifying_variance():
    """Full SUR changes the binary fit while preserving marginal latent scale."""
    N, T = 6, 4
    rng = np.random.default_rng(202)
    y_bin = rng.choice([0.0, 1.0], size=(N, T)).astype(np.float32)
    y_cont = rng.standard_normal((N, T)).astype(np.float32)
    X = rng.standard_normal((N, 3)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[0, :] = 1.0

    y_data = {"bin": y_bin, "cont": y_cont}
    outcome_types = {"bin": "binary", "cont": "continuous"}

    # Fit with SUR on
    cfg_on = LongBetConfig(
        sigma_prior_a=2, sigma_prior_b=1,
        num_sweeps=4,
        num_burnin=2,
        num_chains=1,
        num_trees_pr=3,
        num_trees_trt=3,
        random_seed=42,
        sur=True,
        sur_prior_var=1.0,
    )
    m_on = LongBetMulti(cfg_on).fit(y=y_data, x=X, z=z, outcome=outcome_types)

    # Fit with SUR off (same seed)
    cfg_off = LongBetConfig(
        num_sweeps=4,
        num_burnin=2,
        num_chains=1,
        num_trees_pr=3,
        num_trees_trt=3,
        random_seed=42,
        sur=False,
        sur_prior_var=1.0,
    )
    m_off = LongBetMulti(cfg_off).fit(y=y_data, x=X, z=z, outcome=outcome_types)

    # Full SUR includes downstream information in the binary mean/latent
    # updates. Unlike the old recursive approximation, it must NOT reproduce
    # the independent binary trace when nonzero coupling is learned.
    assert not np.allclose(m_on.trace.traces[0].beta, m_off.trace.traces[0].beta)
    # Binary sigma2 must be 1.0 in both
    np.testing.assert_allclose(m_on.trace.traces[0].sigma2, 1.0)
    np.testing.assert_allclose(m_off.trace.traces[0].sigma2, 1.0)
