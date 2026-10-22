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

"""Tests for exact saved iteration schedules, chain-major flattening, and thinning."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet import LongBetConfig, LongBetMulti
from longbet._multi_io import load_multi_npz, save_multi_npz
from longbet._multi_loop import run_multi_longbet_mcmc


def test_chain_major_flattening_and_thinning():
    """Verify that draws are saved chain-major (d = c * K + k) and thinning schedule is exact."""
    N, T = 6, 4
    rng = np.random.default_rng(77)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    X = rng.standard_normal((N, 3)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[0, :] = 1.0

    # 2 chains, 6 sweeps, burnin=2, skip=2
    # Number of saved sweeps per chain = 6
    cfg = LongBetConfig(
        sigma_prior_a=2, sigma_prior_b=1,
        num_sweeps=6,
        num_burnin=2,
        n_skip=2,
        num_chains=2,
        num_trees_pr=3,
        num_trees_trt=3,
        random_seed=42,
        sur=True,
    )

    model = LongBetMulti(cfg).fit(y={"a": y1, "b": y2}, x=X, z=z)

    assert model.trace is not None
    # gamma_loadings shape is (C, K, M, M) = (2, 6, 2, 2)
    assert model.trace.gamma_loadings.shape == (2, 6, 2, 2)

    # Gamma_draws shape must be (M*M, D) where D = C * K = 12
    gd = model.Gamma_draws
    assert gd.shape == (4, 12)

    # Verify that Gamma_draws columns correspond to chain-major draws
    gl = np.asarray(model.trace.gamma_loadings)
    for c in range(2):
        for k in range(6):
            d = c * 6 + k
            expected_mat = gl[c, k].reshape(4)
            np.testing.assert_allclose(gd[:, d], expected_mat)


def test_zero_burnin_and_single_chain():
    """Zero burn-in and C=1 must correctly record sweep 0 with zero Gamma."""
    N, T = 6, 4
    rng = np.random.default_rng(88)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    X = rng.standard_normal((N, 3)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)

    cfg = LongBetConfig(
        sigma_prior_a=2, sigma_prior_b=1,
        num_sweeps=3,
        num_burnin=0,
        n_skip=1,
        num_chains=1,
        num_trees_pr=3,
        num_trees_trt=3,
        random_seed=42,
    )

    model = LongBetMulti(cfg).fit(y={"a": y1, "b": y2}, x=X, z=z)
    assert model.trace.gamma_loadings.shape == (3, 2, 2)
    # The first draw (sweep 0) must have Gamma == 0
    np.testing.assert_allclose(model.trace.gamma_loadings[0], 0.0)
