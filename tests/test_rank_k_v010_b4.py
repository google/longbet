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

"""v0.10 item B4: arm variances from per-arm outcome variances (``arm_variance_source="direct"``)."""

from __future__ import annotations

import pathlib
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import _rank_k_v010_mocks as mocks  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet  # noqa: E402  pylint: disable=g-import-not-at-top

_ARMS4 = ("c", "b1", "b2", "b3")
_STAR = [(1, 0), (2, 0), (3, 0)]
_MU = np.array([0.0, 0.35, 0.55, 0.65])
_V_HET = np.array([0.0012, 0.0020, 0.0035, 0.0060])
_N = np.array([200.0, 200.0, 200.0, 200.0])


def _star_fits(seed: int):
    rng = np.random.default_rng(seed)
    return mocks.simulate_shared_arm_fits(rng, _MU, _V_HET, _STAR, _ARMS4, n_draws=300, n_units=4, t_post=2)[0]


def test_direct_recovers_heteroscedastic_arm_variances_on_a_star() -> None:
    fits = _star_fits(1)
    sigma2 = {a: float(v * n) for a, v, n in zip(_ARMS4, _V_HET, _N)}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rec = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_scale=1.0, arm_variance_source="direct",
            arm_outcome_variances=sigma2, arm_sample_sizes=_N, n_proposals=20_000,
        )
    assert not any(w.category is longbet.RankKDesignWarning for w in caught)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        back = longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, prior_scale=1.0, n_proposals=20_000)
    v_dir = np.asarray(rec.arm_variances, dtype=np.float64)
    v_back = np.asarray(back.arm_variances, dtype=np.float64)
    err_dir = np.linalg.norm(v_dir - _V_HET) / np.linalg.norm(_V_HET)
    err_back = np.linalg.norm(v_back - _V_HET) / np.linalg.norm(_V_HET)
    assert err_dir < 1e-8
    assert err_back > 0.2
    assert rec.arm_variance_source == "direct"
    assert rec.bipartite_policy is None
    np.testing.assert_allclose(rec.sigma_lik, np.diag(v_dir[1:]) + v_dir[0], rtol=1e-12)


def test_direct_scale_matches_pairwise_variances_and_validates_inputs() -> None:
    fits = _star_fits(2)
    # Outcome variances known only up to a common factor: c absorbs it.
    sigma2 = [float(3.0 * v * n) for v, n in zip(_V_HET, _N)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_scale=1.0, arm_variance_source="direct",
            arm_outcome_variances=sigma2, arm_sample_sizes=_N, n_proposals=10_000,
        )
    np.testing.assert_allclose(np.asarray(rec.arm_variances), _V_HET, rtol=1e-8)
    assert rec.nnls_residual < 1e-10
    with pytest.raises(ValueError, match="requires arm_outcome_variances"):
        longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, arm_variance_source="direct")
    with pytest.raises(ValueError, match="arm_variance_source"):
        longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, arm_variance_source="oracle")
    with pytest.raises(ValueError, match="finite positive"):
        longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, arm_variance_source="direct", arm_outcome_variances=[1.0, -1.0, 1.0, 1.0]
        )


def test_direct_threads_through_multi_outcome() -> None:
    fits = {"x": _star_fits(3), "y": _star_fits(4)}
    sigma2 = {a: float(v * n) for a, v, n in zip(_ARMS4, _V_HET, _N)}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = longbet.reconcile_multi_arm_multi_outcome_draws(
            fits, arms=_ARMS4, prior_scales=1.0, arm_variance_source="direct",
            arm_outcome_variances={"x": sigma2, "y": sigma2}, sample_sizes=_N,
            composite_outcomes={"t": ("x", "y")}, n_proposals=10_000,
        )
    for m in ("x", "y"):
        np.testing.assert_allclose(np.asarray(res.by_outcome[m].arm_variances), _V_HET, rtol=1e-8)
        assert res.by_outcome[m].arm_variance_source == "direct"
    assert res.max_abs_composite_discrepancy["t"] < 1e-12
