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

"""v0.10 item B1: spike-and-slab bundle prior (``prior_mode="bundle_spike_halfnormal"``)."""

from __future__ import annotations

import math
import pathlib
import sys
import warnings

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import _rank_k_v010_mocks as mocks  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet  # noqa: E402  pylint: disable=g-import-not-at-top

_ARMS4 = ("c", "b1", "b2", "b3")
_PAIRS4 = mocks.all_pairs(4)
_V4 = np.array([0.0020, 0.0022, 0.0024, 0.0025])


def _fits_from_mu(rng: np.random.Generator, mu_full: np.ndarray, v: np.ndarray, n_draws: int = 300):
    return mocks.simulate_shared_arm_fits(
        rng, mu_full, v, _PAIRS4, _ARMS4, n_draws=n_draws, n_units=4, t_post=2
    )[0]


def test_spike_slab_point_mass_and_closed_form_weight() -> None:
    rng = np.random.default_rng(1)
    fits = _fits_from_mu(rng, np.zeros(4), _V4)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_mode="bundle_spike_halfnormal", prior_scale=0.3,
            n_proposals=40_000, seed=2,
        )
    w = rec.spike_posterior_prob
    assert 0.0 < w < 1.0
    at_zero = np.all(rec.tau_mu_draws == 0.0, axis=1)
    # Draw-wise point mass at zero with the closed-form posterior weight (binomial tolerance).
    s = rec.tau_mu_draws.shape[0]
    assert abs(float(np.mean(at_zero)) - w) <= 4.0 * math.sqrt(w * (1 - w) / s) + 1e-12
    # Zero Mode-1 draws give exactly zero trajectories on every contrast.
    traj = rec.contrast_trajectory_draws[("b3", "c")]
    assert np.all(traj[at_zero] == 0.0)
    assert rec.prior_mode == "bundle_spike_halfnormal"


def test_spike_weight_small_on_active_effect() -> None:
    rng = np.random.default_rng(3)
    fits = _fits_from_mu(rng, np.array([0.0, 0.35, 0.55, 0.65]), _V4)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_mode="bundle_spike_halfnormal", prior_scale=1.0,
            n_proposals=40_000, seed=4,
        )
    assert rec.spike_posterior_prob < 1e-6
    assert np.all(rec.tau_mu_draws[:, -1] > 0.0)


def test_spike_requires_zero_in_support() -> None:
    rng = np.random.default_rng(5)
    fits = _fits_from_mu(rng, np.zeros(4), _V4)
    with pytest.raises(ValueError, match="spike"):
        longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_mode="bundle_spike_halfnormal", prior_scale=1.0, tau_clip=(0.01, np.inf)
        )
    with pytest.raises(ValueError, match="spike_prob"):
        longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_mode="bundle_spike_halfnormal", prior_scale=1.0, spike_prob=1.0
        )


def test_spike_multi_outcome_coherent_total_composite_vanishes_with_spike() -> None:
    rng = np.random.default_rng(6)
    fits = {m: _fits_from_mu(rng, np.zeros(4), _V4) for m in ("dense", "sparse")}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = longbet.reconcile_multi_arm_multi_outcome_draws(
            fits, arms=_ARMS4, prior_mode="bundle_spike_halfnormal", prior_scales=0.3,
            composite_outcomes={"total": ("dense", "sparse")}, n_proposals=40_000, seed=7,
        )
    assert res.max_abs_composite_discrepancy["total"] < 1e-12
    total = res.by_outcome["total"].tau_mu_draws
    assert np.mean(np.all(total == 0.0, axis=1)) > 0.05


def _randomized_rank(post: np.ndarray, truth: float, rng: np.random.Generator) -> float:
    below = int(np.sum(post < truth))
    ties = int(np.sum(post == truth))
    return (below + int(rng.integers(0, ties + 1))) / float(post.size + 1)


def _sbc_spike_ranks(n_trials: int, n_draws: int, seed: int) -> np.ndarray:
    """SBC: mu* ~ 0.5 delta_0 + 0.5 bundle half-normal slab, y = B mu* + E eps, reconcile."""
    rng = np.random.default_rng(seed)
    scale = 0.3
    ranks = np.zeros(n_trials)
    for t in range(n_trials):
        if rng.random() < 0.5:
            mu_star = np.zeros(3)
        else:
            tot = abs(rng.standard_normal()) * scale
            mu_star = np.cumsum(rng.dirichlet(np.ones(3)) * tot)
        fits = _fits_from_mu(rng, np.concatenate([[0.0], mu_star]), _V4, n_draws=n_draws)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rec = longbet.reconcile_multi_arm_draws(
                fits, arms=_ARMS4, prior_mode="bundle_spike_halfnormal", prior_scale=scale,
                n_proposals=20_000, seed=int(rng.integers(1, 10**8)),
            )
        ranks[t] = _randomized_rank(rec.tau_mu_draws[:, -1], float(mu_star[-1]), rng)
    return ranks


def _chi2_p(ranks: np.ndarray, n_bins: int) -> float:
    counts, _ = np.histogram(ranks, bins=np.linspace(0.0, 1.0, n_bins + 1))
    return float(stats.chisquare(counts).pvalue)


def test_sbc_spike_slab_fast() -> None:
    assert _chi2_p(_sbc_spike_ranks(80, 200, 808), 5) > 0.01


@pytest.mark.slow
def test_sbc_spike_slab_slow() -> None:
    assert _chi2_p(_sbc_spike_ranks(400, 300, 909), 10) >= 0.05
