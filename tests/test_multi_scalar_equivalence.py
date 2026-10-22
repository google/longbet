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

"""Uncoupled multi fits must reduce to the scalar sampler, not just look plausible.

Match initialization and sweep keys explicitly: the public APIs deliberately
use different streams for scalar and multi fits, even with the same base seed.
Compare full traces/predictions, not estimates against truth at a different
exposure for never-treated units.
"""

import dataclasses

import equinox as eqx
import jax
import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, LongBetMulti, effect_draws
from longbet import _loop, _model, _multi_model
from longbet._state import split_chain_fields
from longbet._step import longbet_single_step


@pytest.mark.parametrize("chains", [1, 2])
def test_uncoupled_fit_matches_scalar_with_matched_keys(monkeypatch, chains):
    rng = np.random.default_rng(871)
    n, t = 24, 6
    x = np.asfortranarray(rng.normal(size=(n, 2)))
    z = np.zeros((n, t))
    z[:16, 2:] = 1
    y = {
        "large": np.asfortranarray(40 + 7 * rng.normal(size=(n, t))),
        "binary": (rng.normal(size=(n, t)) > 0.8).astype(float),
        "small": 0.1 * rng.normal(size=(n, t)),
    }
    # Non-nested masks are valid when coupling is off; neither fill values nor
    # R/Fortran memory layout may change the scalar likelihood.
    y["large"][0, 0] = np.nan
    y["small"][1, 1] = np.nan
    types = {"large": "continuous", "binary": "binary", "small": "continuous"}
    cfg = LongBetConfig(num_burnin=2, num_sweeps=3, n_skip=2,
                        num_chains=chains, num_trees_pr=3, num_trees_trt=2,
                        max_depth_pr=3, max_depth_trt=3,
                        min_points_per_leaf_pr=2, min_points_per_leaf_trt=2,
                        sur=False, random_seed=63, device="cpu")
    starts = []
    original_multi_init = _multi_model.init_multi_longbet
    def capture_multi_init(**kwargs):
        state = original_multi_init(**kwargs)
        starts.append(state)
        return state
    with monkeypatch.context() as patch:
        patch.setattr(_multi_model, "init_multi_longbet", capture_multi_init)
        multi = LongBetMulti(cfg).fit(y, x, z, outcome=types)
    ze = np.broadcast_to(z[0], z.shape).copy()
    mp = multi.predict(x, ze)
    original_init = _model.init_longbet

    for user_index, name in enumerate(y):
        m = multi.inverse_order[user_index]

        def scalar_init(*, chain_key, **kwargs):
            state = original_init(chain_key=jax.random.fold_in(chain_key, m), **kwargs)
            for (path, actual), expected in zip(jax.tree_util.tree_flatten_with_path(starts[0].states[m])[0], jax.tree.leaves(state)):
                np.testing.assert_allclose(actual, expected, atol=1e-6, rtol=1e-6,
                                           err_msg=f"init {name}: {path}")
            return state

        def scalar_step(key, state):
            k_outcome = jax.random.fold_in(key, 0)
            if not state.has_chain_axis:
                return longbet_single_step(jax.random.fold_in(k_outcome, m), state)
            keys = jax.random.split(k_outcome, chains)
            per, shared = split_chain_fields(state)

            def one(k, p):
                updated = longbet_single_step(jax.random.fold_in(k, m), eqx.combine(p, shared))
                return split_chain_fields(updated)[0]

            return eqx.combine(jax.vmap(one)(keys, per), shared)

        with monkeypatch.context() as patch:
            patch.setattr(_model, "init_longbet", scalar_init)
            patch.setattr(_loop, "longbet_step", scalar_step)
            scalar = LongBet(dataclasses.replace(cfg, outcome=types[name])).fit(y[name], x, z)
        for (path, actual), expected in zip(jax.tree_util.tree_flatten_with_path(multi[name].trace)[0],
                                    jax.tree.leaves(scalar.trace)):
            np.testing.assert_allclose(actual, expected, atol=3e-5, rtol=3e-5,
                                       err_msg=f"trace {name}: {path}")
        sp = scalar.predict(x, ze)
        np.testing.assert_allclose(effect_draws(mp, name), effect_draws(sp),
                                   atol=3e-5, rtol=3e-5)
