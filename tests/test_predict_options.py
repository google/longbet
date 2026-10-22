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

"""Prediction options added for the chapter's forecast sweep.

None of these had tests when they were written; each covers a way the feature
can be wrong without looking wrong.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from _panels import make_staggered_panel
from longbet import LongBet, LongBetConfig

FAST = dict(num_sweeps=10, num_burnin=6, num_trees_pr=5, num_trees_trt=5,
            num_chains=1, random_seed=3)


def _fitted(**cfg):
    data = make_staggered_panel(N=50, T=8)
    model = LongBet(LongBetConfig(**{**FAST, **cfg}))
    model.fit(y=data["y"], x=data["x"], z=data["z"], t=data["t"])
    return model, data


def _extended(data, extra=5):
    """A panel that runs past the fitted exposure horizon, forcing a projection."""
    z, t = data["z"], data["t"]
    n, n_t = z.shape
    z2 = np.concatenate([z, np.repeat(z[:, [-1]], extra, axis=1)], axis=1)
    t2 = np.arange(1, n_t + extra + 1, dtype=np.float32)
    return z2, t2


class TestProjectionKernelOverride:
    def test_override_changes_only_the_projection(self):
        """The in-sample block must be untouched; only the tail may move.

        This is the whole claim the option rests on. If the override leaked into
        the fitted draws it would be a silent refit, and a sweep over it would
        not be the sensitivity analysis it is described as.
        """
        model, data = _fitted()
        z2, t2 = _extended(data)
        base = model.predict(x=data["x"], z=z2, t=t2, summary_only=True)
        wide = model.predict(x=data["x"], z=z2, t=t2, summary_only=True,
                             lambda_knl=20.0, sig_knl=3.0)

        fitted_block = model.S_max_ + 1
        np.testing.assert_allclose(
            base.beta_values[:, :fitted_block],
            wide.beta_values[:, :fitted_block],
            rtol=0, atol=0,
            err_msg="overriding the projection kernel altered the fitted beta draws",
        )
        assert base.beta_values.shape[1] > fitted_block, "no projection happened"
        assert not np.allclose(
            base.beta_values[:, fitted_block:], wide.beta_values[:, fitted_block:]
        ), "the override did not reach the projection"

    def test_default_matches_the_fitted_kernel(self):
        model, data = _fitted(sig_knl=1.0, lambda_knl=2.0)
        z2, t2 = _extended(data)
        implicit = model.predict(x=data["x"], z=z2, t=t2, summary_only=True)
        explicit = model.predict(x=data["x"], z=z2, t=t2, summary_only=True,
                                 sig_knl=1.0, lambda_knl=2.0)
        np.testing.assert_array_equal(implicit.beta_values, explicit.beta_values)

    def test_a_short_lengthscale_reverts_to_the_estimated_level(self):
        """The limiting common level is a conditional mean, not mean(beta).

        Test known trajectories against an independent Gaussian conditional.
        Ten noisy forecasts from an unconverged fit cannot establish a universal
        ordering of realized long/short-lengthscale sample means.
        """
        import jax
        from longbet._gp import forecast_beta_gp

        s = np.arange(6, dtype=float)
        observed = 1 + .3 * s
        draws = 6000
        means = []
        for length in (1., 25.):
            K = np.exp(-.5 * ((s[:, None] - s) / length)**2) + 1 + 1e-6*np.eye(len(s))
            cross = np.exp(-.5 * ((30 - s) / length)**2) + 1
            expected_mean = cross @ np.linalg.solve(K, observed)
            expected_var = 2 + 2e-6 - cross @ np.linalg.solve(K, cross)
            result = np.asarray(forecast_beta_gp(
                jax.random.key(173), np.tile(observed, (draws, 1)), s,
                np.array([30.]), sig_knl=1., lambda_knl=length))[:, -1]
            np.testing.assert_allclose(result.mean(), expected_mean,
                                       atol=4*np.sqrt(expected_var/draws))
            np.testing.assert_allclose(result.var(), expected_var, rtol=.06)
            if length == 1:
                common_level = np.ones(len(s)) @ np.linalg.solve(K, observed)
                np.testing.assert_allclose(expected_mean, common_level, atol=1e-10)
            means.append(expected_mean)
        assert abs(means[1] - means[0]) > .5


class TestForestEvaluationCache:
    def test_off_by_default(self):
        """summary_only promises bounded memory; the cache must not break it."""
        model, data = _fitted()
        model.predict(x=data["x"], z=data["z"], t=data["t"], summary_only=True)
        assert not hasattr(model, "_eval_cache") or model._eval_cache[0] is None

    def test_cached_and_uncached_agree(self):
        model, data = _fitted()
        kw = dict(x=data["x"], z=data["z"], t=data["t"], summary_only=True)
        plain = model.predict(**kw)
        model.predict(**kw, cache_forest_evaluations=True)      # populate
        cached = model.predict(**kw, cache_forest_evaluations=True)  # hit
        np.testing.assert_allclose(plain.tau_summary.mean, cached.tau_summary.mean,
                                   rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(plain.att_full, cached.att_full,
                                   rtol=1e-6, atol=1e-6)

    def test_a_new_design_evicts_rather_than_accumulates(self):
        """One entry only: an unbounded cache of (draws x cells) arrays would be
        a memory leak on exactly the panels this package targets."""
        model, data = _fitted()
        z_alt = data["z"].copy()
        z_alt[0, :] = 0
        for z in (data["z"], z_alt, data["z"], z_alt):
            model.predict(x=data["x"], z=z, t=data["t"], summary_only=True,
                          cache_forest_evaluations=True)
        key, blocks = model._eval_cache
        assert key is not None and blocks is not None
        # A single design's worth of blocks, not four.
        n_blocks = len(blocks)
        expected = -(-data["z"].size // (model._eval_cache[0][5] or data["z"].size))
        assert n_blocks == expected

    def test_a_changed_design_is_not_served_from_cache(self):
        """The dangerous failure: returning one design's forests for another."""
        model, data = _fitted()
        z_alt = data["z"].copy()
        z_alt[: len(z_alt) // 2] = 0          # materially different treatment
        first = model.predict(x=data["x"], z=data["z"], t=data["t"],
                              summary_only=True, cache_forest_evaluations=True)
        second = model.predict(x=data["x"], z=z_alt, t=data["t"],
                               summary_only=True, cache_forest_evaluations=True)
        again = model.predict(x=data["x"], z=data["z"], t=data["t"],
                              summary_only=True, cache_forest_evaluations=True)
        assert not np.allclose(first.tau_summary.mean, second.tau_summary.mean)
        np.testing.assert_allclose(first.tau_summary.mean, again.tau_summary.mean,
                                   rtol=1e-6, atol=1e-6)


def test_multichain_save_load_round_trip():
    """Regression: the rebuilt trace's move counters need the chain axis.

    The bug this covers was found and fixed without a test; a single-chain
    round trip does not exercise it.
    """
    model, data = _fitted(num_chains=3)
    before = model.predict(x=data["x"], z=data["z"], t=data["t"])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "m.npz"
        model.save(path)
        reloaded = LongBet.load(path)
        after = reloaded.predict(x=data["x"], z=data["z"], t=data["t"])
    assert before.tauhats.shape[-1] == 3 * FAST["num_sweeps"]
    np.testing.assert_allclose(before.tauhats, after.tauhats, atol=1e-5)
    np.testing.assert_allclose(before.att_full, after.att_full, atol=1e-5)
