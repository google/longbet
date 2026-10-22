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

"""Numerical diagnostics for reported means and individual target checks."""

import arviz as az
import numpy as np
import pytest

from longbet import att_stability, compute_ess


def test_mean_mcse_uses_mean_ess_for_skewed_autocorrelated_draws():
    rng = np.random.default_rng(823)
    log_draws = np.zeros((4, 1000))
    for draw in range(1, log_draws.shape[1]):
        log_draws[:, draw] = 0.95 * log_draws[:, draw - 1] + rng.normal(size=4)
    draws = np.exp(log_draws)
    result = att_stability(draws[..., None], warn=False)

    expected_mcse = float(az.mcse(draws, method="mean"))
    np.testing.assert_allclose(result.by_exposure["mcse"], [expected_mcse])
    np.testing.assert_allclose(
        result.by_exposure["ess_mean"], [float(az.ess(draws, method="mean"))]
    )
    bulk_approximation = np.std(draws) / np.sqrt(float(az.ess(draws, method="bulk")))
    # Rank normalization changes the estimand of ESS. This example would fail
    # substantially if MCSE of a raw mean were approximated using bulk ESS.
    assert abs(bulk_approximation / expected_mcse - 1) > 0.25


def test_one_poorly_mixed_exposure_fails_even_when_median_ess_passes():
    draws = np.random.default_rng(72).normal(size=(4, 1000, 3))
    draws[:, :, 2] += np.array([-4, -2, 2, 4])[:, None]
    result = att_stability(draws, warn=False)
    assert result.summary["ess_median"] > 400
    assert result.summary["ess_min"] < 400
    assert result.summary["ess_ok"] is False
    assert result.summary["rhat_ok"] is False


def test_tail_ess_must_also_pass_at_each_exposure():
    draws = np.random.default_rng(1).normal(size=(4, 1000))
    bulk = float(az.ess(draws, method="bulk"))
    tail = float(az.ess(draws, method="tail", prob=0.05))
    assert tail < bulk
    threshold = (bulk + tail) / 2
    result = att_stability(draws[..., None], min_ess=threshold, warn=False)
    assert result.summary["ess_min"] > threshold
    assert result.summary["ess_tail_min"] < threshold
    assert result.summary["ess_ok"] is False


@pytest.mark.parametrize("undefined", [np.nan, 0.0])
def test_undefined_exposure_cannot_be_hidden_by_other_good_diagnostics(undefined):
    draws = np.random.default_rng(56).normal(size=(4, 1000, 3))
    draws[:, :, 2] = undefined
    result = att_stability(draws, warn=False)
    assert result.summary["ess_median"] > 400
    assert result.summary["rhat_max"] < 1.01
    assert result.summary["ess_ok"] is False
    assert result.summary["rhat_ok"] is False
    assert np.isnan(result.by_exposure["mcse"][2])


def test_small_effect_variation_is_not_misclassified_as_constant():
    draws = 1.0 + 1e-7 * np.random.default_rng(9).normal(size=(4, 1000))
    result = att_stability(draws[..., None], warn=False)
    assert np.isfinite(result.by_exposure["ess_bulk"][0])
    assert np.isfinite(result.by_exposure["mcse"][0])
    np.testing.assert_allclose(
        result.by_exposure["mcse"], [float(az.mcse(draws, method="mean"))]
    )


@pytest.mark.parametrize("shape", [(0, 10, 2), (4, 0, 2), (4, 10, 0)])
def test_empty_draw_axes_fail_with_an_explanatory_error(shape):
    with pytest.raises(ValueError, match="nonempty"):
        att_stability(np.empty(shape), warn=False)


def test_unknown_ess_method_is_not_silently_replaced_with_bulk():
    with pytest.raises(ValueError, match="Unknown ESS method"):
        compute_ess(np.arange(10.0), method="typo")
