# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""v0.10 item B3: robust (Huber IRLS) cut-space projection (``projection="huber"``)."""

from __future__ import annotations

import math
import pathlib
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import _rank_k_v010_mocks as mocks  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet._rank_k_multi_arm as rkma  # noqa: E402  pylint: disable=g-import-not-at-top

_ARMS4 = ("c", "b1", "b2", "b3")
_PAIRS4 = mocks.all_pairs(4)
_MU = np.array([0.0, 0.35, 0.55, 0.65])
_V4 = np.array([0.0020, 0.0022, 0.0024, 0.0025])


def _fits(seed: int, pair_bias: dict | None = None):
    rng = np.random.default_rng(seed)
    return mocks.simulate_shared_arm_fits(
        rng, _MU, _V4, _PAIRS4, _ARMS4, n_draws=300, n_units=4, t_post=2, pair_bias=pair_bias
    )


def test_huber_equals_ols_under_exact_shared_arm_model() -> None:
    fits, _ = _fits(1)
    kw = dict(arms=_ARMS4, prior_scale=1.0, n_proposals=20_000, seed=3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ols = longbet.reconcile_multi_arm_draws(fits, **kw)
        hub = longbet.reconcile_multi_arm_draws(fits, projection="huber", **kw)
    np.testing.assert_array_equal(hub.mu_hat, ols.mu_hat)
    np.testing.assert_array_equal(hub.tau_mu_draws, ols.tau_mu_draws)
    assert hub.flagged_pairs == ()
    assert all(w == 1.0 for w in hub.projection_weights.values())
    assert ols.projection == "ols" and ols.projection_weights is None


@pytest.mark.parametrize("k_sd", [2.0, 4.0])
def test_huber_flags_and_removes_one_corrupted_fit(k_sd: float) -> None:
    corrupted = (3, 1)
    sd = math.sqrt(_V4[3] + _V4[1])
    fits, _ = _fits(2, pair_bias={corrupted: k_sd * sd})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ols = longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, prior_scale=1.0, n_proposals=20_000)
        hub = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_scale=1.0, projection="huber", n_proposals=20_000
        )
    assert hub.flagged_pairs[0] == ("b3", "b1")
    assert hub.projection_weights[("b3", "b1")] < 0.1
    # The robust projection reproduces the projection of the five clean fits.
    keep = [p for p in _PAIRS4 if p != corrupted]
    b_keep = np.zeros((5, 3))
    y_keep = np.zeros(5)
    for i, (j, l) in enumerate(keep):
        if j:
            b_keep[i, j - 1] += 1
        if l:
            b_keep[i, l - 1] -= 1
        y_keep[i] = float(np.mean(fits[(_ARMS4[j], _ARMS4[l])].att_draws))
    mu_clean = np.linalg.lstsq(b_keep, y_keep, rcond=None)[0]
    np.testing.assert_allclose(hub.mu_hat, mu_clean, atol=1e-6)
    assert np.max(np.abs(ols.mu_hat - mu_clean)) > 0.25 * k_sd * sd / 2.0


def test_huber_irls_helper_bounded_weights() -> None:
    b = np.array([[1, 0], [0, 1], [1, -1], [-1, 1]], dtype=float)
    y = np.array([1.0, 2.0, -1.0, 5.0])
    mu, w, scale = rkma._huber_cut_space_projection(b, y, np.ones(4))
    assert np.all((w > 0) & (w <= 1))
    assert scale > 0
    assert mu.shape == (2,)


def test_huber_threads_through_multi_outcome() -> None:
    fits_x, _ = _fits(4, pair_bias={(2, 0): 0.3})
    fits_y, _ = _fits(5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = longbet.reconcile_multi_arm_multi_outcome_draws(
            {"x": fits_x, "y": fits_y},
            arms=_ARMS4,
            prior_scales=1.0,
            projection="huber",
            composite_outcomes={"t": ("x", "y")},
            n_proposals=20_000,
        )
    assert res.by_outcome["x"].flagged_pairs[0] == ("b2", "c")
    assert res.by_outcome["y"].flagged_pairs == ()
    assert res.max_abs_composite_discrepancy["t"] < 1e-12
