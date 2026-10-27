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

"""v0.10 item B2: empirical-Bayes prior scale (``prior_scale="empirical_bayes"``)."""

from __future__ import annotations

import math
import pathlib
import sys
import warnings

import numpy as np
import pytest
from scipy import integrate, stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import _rank_k_v010_mocks as mocks  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet._rank_k_multi_arm as rkma  # noqa: E402  pylint: disable=g-import-not-at-top

_ARMS4 = ("c", "b1", "b2", "b3")
_PAIRS4 = mocks.all_pairs(4)
_V4 = np.array([0.0020, 0.0022, 0.0024, 0.0025])


def _fits(rng: np.random.Generator, mu_full: np.ndarray, n_draws: int = 300):
    return mocks.simulate_shared_arm_fits(
        rng, mu_full, _V4, _PAIRS4, _ARMS4, n_draws=n_draws, n_units=4, t_post=2
    )[0]


@pytest.mark.parametrize(
    ("family", "kwargs"),
    [
        ("bundle_halfnormal", {}),
        ("bundle_halfstudentt", {}),
        ("normal", {"low": -math.inf, "high": math.inf, "mean": 0.1}),
        ("normal", {"low": 0.0, "high": math.inf, "mean": 0.0}),
        ("bundle_halfnormal", {"spike_prob": 0.5}),
    ],
)
def test_log_marginal_matches_numerical_integration(family: str, kwargs: dict) -> None:
    sigma, t_hat, s_lik = 0.4, 0.3, 0.1

    def _prior_pdf(tau: float) -> float:
        if family == "bundle_halfnormal":
            return 2.0 * stats.norm.pdf(tau, scale=sigma) if tau > 0 else 0.0
        if family == "bundle_halfstudentt":
            sc = 0.5 * sigma
            return 2.0 * stats.t.pdf(tau / sc, df=2.5) / sc if tau > 0 else 0.0
        lo, hi, m = kwargs["low"], kwargs["high"], kwargs["mean"]
        mass = stats.norm.cdf((hi - m) / sigma) - stats.norm.cdf((lo - m) / sigma)
        return stats.norm.pdf(tau, loc=m, scale=sigma) / mass if lo <= tau <= hi else 0.0

    ref, _ = integrate.quad(
        lambda tau: _prior_pdf(tau) * stats.norm.pdf(t_hat, loc=tau, scale=s_lik), -5.0, 5.0, points=[0.0, t_hat],
        limit=400,
    )
    if "spike_prob" in kwargs:
        ref = 0.5 * stats.norm.pdf(t_hat, scale=s_lik) + 0.5 * ref
    got = rkma._log_marginal_1d(sigma, t_hat, s_lik, family, **kwargs)
    assert got == pytest.approx(math.log(ref), abs=2e-3)


def test_empirical_bayes_scale_floor_and_large_effect() -> None:
    # Effect far from zero: the half-normal Type-II ML scale is about sqrt(t^2 - s^2).
    sig = rkma._empirical_bayes_scale(0.65, 0.07, "bundle_halfnormal", floor=0.07)
    assert 0.5 < sig < 1.2
    # Effect at or below zero: the scale sits on the floor.
    assert rkma._empirical_bayes_scale(-0.05, 0.07, "bundle_halfnormal", floor=0.07) == pytest.approx(0.07)


def test_empirical_bayes_removes_in_sample_shrinkage() -> None:
    rng = np.random.default_rng(1)
    mu = np.array([0.0, 0.35, 0.55, 0.65])
    fits = _fits(rng, mu)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec_is = longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, prior_scale="in_sample", n_proposals=40_000)
        rec_eb = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_scale="empirical_bayes", n_proposals=40_000
        )
    top_is = float(np.mean(rec_is.tau_mu_draws[:, -1]))
    top_eb = float(np.mean(rec_eb.tau_mu_draws[:, -1]))
    assert rec_eb.prior_scale > 5.0 * rec_is.prior_scale
    assert abs(top_eb - rec_eb.mu_hat[-1]) < abs(top_is - rec_eb.mu_hat[-1])
    assert abs(top_eb - 0.65) < 0.1


def test_empirical_bayes_independent_family_and_multi_outcome() -> None:
    rng = np.random.default_rng(2)
    mu = np.array([0.0, 0.35, 0.55, 0.65])
    fits = _fits(rng, mu)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_mode="independent_halfnormal", prior_scale="empirical_bayes",
            n_proposals=20_000,
        )
        res = longbet.reconcile_multi_arm_multi_outcome_draws(
            {"x": fits, "y": _fits(rng, 0.5 * mu)},
            arms=_ARMS4,
            prior_scales="empirical_bayes",
            composite_outcomes={"t": ("x", "y")},
            n_proposals=20_000,
        )
    assert isinstance(rec.prior_scale, tuple) and len(rec.prior_scale) == 3
    assert res.composite_mode == "coherent_total"
    assert res.by_outcome["t"].prior_scale > 0.3
    assert res.max_abs_composite_discrepancy["t"] < 1e-12


def _sbc_eb_ranks(n_trials: int, n_draws: int, seed: int) -> np.ndarray:
    """SBC: sigma ~ log-uniform hyperprior, mu* ~ bundle half-normal(sigma), y = B mu* + E eps."""
    rng = np.random.default_rng(seed)
    s_top = math.sqrt(_V4[3] + _V4[0])
    ranks = np.zeros(n_trials)
    for t in range(n_trials):
        sigma = s_top * math.exp(rng.uniform(math.log(1.0), math.log(20.0)))
        tot = abs(rng.standard_normal()) * sigma
        mu_star = np.cumsum(rng.dirichlet(np.ones(3)) * tot)
        fits = _fits(rng, np.concatenate([[0.0], mu_star]), n_draws=n_draws)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rec = longbet.reconcile_multi_arm_draws(
                fits, arms=_ARMS4, prior_scale="empirical_bayes", n_proposals=20_000,
                seed=int(rng.integers(1, 10**8)),
            )
        ranks[t] = float(np.mean(rec.tau_mu_draws[:, -1] < mu_star[-1]))
    return ranks


def _chi2_p(ranks: np.ndarray, n_bins: int) -> float:
    counts, _ = np.histogram(ranks, bins=np.linspace(0.0, 1.0, n_bins + 1))
    return float(stats.chisquare(counts).pvalue)


def test_sbc_empirical_bayes_fast() -> None:
    assert _chi2_p(_sbc_eb_ranks(80, 200, 1201), 5) > 0.01


@pytest.mark.slow
def test_sbc_empirical_bayes_slow() -> None:
    assert _chi2_p(_sbc_eb_ranks(400, 300, 1301), 10) >= 0.05
