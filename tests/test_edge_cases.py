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

"""Degenerate panels that should degrade gracefully rather than crash."""

import warnings

import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, derive_exposure, get_att

FAST = dict(num_sweeps=6, num_burnin=3, num_trees_pr=4, num_trees_trt=4,
            num_chains=1, random_seed=1)


def _covariates(N=30, P=2, seed=0):
    return np.random.default_rng(seed).normal(size=(N, P)).astype(np.float32)


def test_panel_with_no_treated_units():
    """No treatment anywhere: the ATT is undefined, not an exception."""
    N, T = 30, 5
    x = _covariates(N)
    z = np.zeros((N, T), dtype=np.float32)
    y = np.random.default_rng(1).normal(size=(N, T)).astype(np.float32)
    t = np.arange(1, T + 1, dtype=np.float32)

    model = LongBet(LongBetConfig(**FAST)).fit(y=y, x=x, z=z, t=t)
    pred = model.predict(x=x, z=z, t=t)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert np.all(np.isnan(get_att(pred)["att"]))
        assert np.isnan(pred.stability(warn=False).summary["ess_median"])
    # The effect itself is still estimated, and identically zero is plausible.
    assert np.all(np.isfinite(pred.tauhats))


def test_single_period_panel():
    N = 30
    x = _covariates(N)
    rng = np.random.default_rng(2)
    z = rng.binomial(1, 0.5, (N, 1)).astype(np.float32)
    y = rng.normal(size=(N, 1)).astype(np.float32)
    t = np.array([1.0], dtype=np.float32)

    model = LongBet(LongBetConfig(**FAST)).fit(y=y, x=x, z=z, t=t)
    pred = model.predict(x=x, z=z, t=t)
    assert pred.tauhats.shape == (N, 1, FAST["num_sweeps"])
    assert np.all(np.isfinite(pred.tauhats))


def test_unit_with_every_cell_missing():
    """Its intercept falls back to the prior instead of dividing by zero."""
    N, T = 30, 5
    x = _covariates(N)
    rng = np.random.default_rng(3)
    z = np.zeros((N, T), dtype=np.float32)
    z[:15, 2:] = 1.0
    y = rng.normal(size=(N, T)).astype(np.float32)
    y[3, :] = np.nan
    t = np.arange(1, T + 1, dtype=np.float32)

    model = LongBet(LongBetConfig(**FAST)).fit(y=y, x=x, z=z, t=t)
    assert np.all(np.isfinite(np.asarray(model.trace.gamma)))
    pred = model.predict(x=x, z=z, t=t)
    assert np.all(np.isfinite(pred.tauhats))


def test_all_missing_outcomes_is_refused():
    N, T = 10, 4
    x = _covariates(N)
    y = np.full((N, T), np.nan, dtype=np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    with pytest.raises(ValueError, match="no observed cells"):
        LongBet(LongBetConfig(**FAST)).fit(y=y, x=x, z=z)


def test_constant_outcome_is_refused():
    N, T = 10, 4
    x = _covariates(N)
    y = np.ones((N, T), dtype=np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    with pytest.raises(ValueError, match="constant"):
        LongBet(LongBetConfig(**FAST)).fit(y=y, x=x, z=z)


def test_non_integer_time_spacing_warns():
    """The exposure index is an integer grid, so uneven time must be flagged."""
    N = 30
    x = _covariates(N)
    z = np.zeros((N, 4), dtype=np.float32)
    z[:15, 2:] = 1.0
    y = np.random.default_rng(4).normal(size=(N, 4)).astype(np.float32)
    with pytest.warns(UserWarning, match="integer-spaced"):
        LongBet(LongBetConfig(**FAST)).fit(
            y=y, x=x, z=z, t=np.array([0.0, 0.5, 1.7, 2.3], dtype=np.float32)
        )


def test_decreasing_time_is_refused():
    N = 10
    x = _covariates(N)
    z = np.zeros((N, 3), dtype=np.float32)
    y = np.random.default_rng(5).normal(size=(N, 3)).astype(np.float32)
    with pytest.raises(ValueError, match="strictly increasing"):
        LongBet(LongBetConfig(**FAST)).fit(
            y=y, x=x, z=z, t=np.array([1.0, 3.0, 2.0], dtype=np.float32)
        )


def test_mismatched_shapes_are_refused():
    x = _covariates(10)
    z = np.zeros((10, 4), dtype=np.float32)
    y = np.zeros((10, 5), dtype=np.float32)
    with pytest.raises(ValueError, match="expected"):
        LongBet(LongBetConfig(**FAST)).fit(y=y, x=x, z=z)
    with pytest.raises(ValueError, match="x must be"):
        LongBet(LongBetConfig(**FAST)).fit(
            y=np.zeros((10, 4), dtype=np.float32), x=_covariates(7), z=z
        )


class TestDeriveExposure:
    """The exposure index counts the first treated period as 1."""

    def test_simultaneous_adoption(self):
        z = np.array([[0, 0, 1, 1], [0, 0, 0, 0]], dtype=np.float32)
        s = derive_exposure(z, np.arange(1, 5, dtype=np.float32))
        np.testing.assert_array_equal(s, [[0, 0, 1, 2], [0, 0, 0, 0]])

    def test_treated_from_the_first_period(self):
        z = np.array([[1, 1, 1]], dtype=np.float32)
        s = derive_exposure(z, np.arange(1, 4, dtype=np.float32))
        np.testing.assert_array_equal(s, [[1, 2, 3]])

    def test_staggered_adoption(self):
        z = np.array([[0, 1, 1, 1], [0, 0, 0, 1]], dtype=np.float32)
        s = derive_exposure(z, np.arange(1, 5, dtype=np.float32))
        np.testing.assert_array_equal(s, [[0, 1, 2, 3], [0, 0, 0, 1]])

    def test_exposure_counts_elapsed_time_not_periods(self):
        """With a gap in the calendar, exposure follows the clock."""
        z = np.array([[0, 1, 1]], dtype=np.float32)
        s = derive_exposure(z, np.array([1.0, 2.0, 5.0], dtype=np.float32))
        np.testing.assert_array_equal(s, [[0, 1, 4]])
