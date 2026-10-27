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

"""v0.10 item A5: block auto-continuation of the hit-and-run slice sampler for ``d >= 4``."""

from __future__ import annotations

import pathlib
import sys
import warnings

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import _rank_k_v010_mocks as mocks  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet._rank_k_multi_arm as rkma  # noqa: E402  pylint: disable=g-import-not-at-top


def _chain_fits(n_arms: int, seed: int, n_draws: int = 200) -> tuple[tuple[str, ...], dict]:
    arms = tuple(f"a{i}" for i in range(n_arms))
    mu = np.concatenate([[0.0], np.linspace(0.3, 0.7, n_arms - 1)])
    v = np.full(n_arms, 0.0036)
    rng = np.random.default_rng(seed)
    fits, _ = mocks.simulate_shared_arm_fits(
        rng, mu, v, mocks.all_pairs(n_arms), arms, n_draws=n_draws, n_units=4, t_post=2
    )
    return arms, fits


def test_first_block_is_bit_identical_when_targets_are_met_immediately() -> None:
    arms, fits = _chain_fits(5, 1)
    kwargs = dict(arms=arms, prior_scale=1.0, seed=3, num_burnin=60)
    off = longbet.reconcile_multi_arm_draws(fits, auto_continue=False, **kwargs)
    on = longbet.reconcile_multi_arm_draws(
        fits, auto_continue=True, target_rhat=10.0, target_ess=0.0, **kwargs
    )
    np.testing.assert_array_equal(on.tau_mu_draws, off.tau_mu_draws)
    assert on.sampler_diagnostics.n_iterations == off.sampler_diagnostics.n_iterations == 60 + 50 * 2
    assert on.sampler_diagnostics.converged is True
    assert off.sampler_diagnostics.converged is None


def test_cap_reached_warns_and_records_iterations() -> None:
    arms, fits = _chain_fits(5, 2)
    with pytest.warns(longbet.RankKConfigWarning, match="did not reach split R-hat"):
        rec = longbet.reconcile_multi_arm_draws(
            fits,
            arms=arms,
            prior_scale=1.0,
            seed=4,
            num_burnin=40,
            auto_continue=True,
            target_rhat=1.0,  # unattainable
            max_iter_multiplier=3.0,
        )
    diag = rec.sampler_diagnostics
    default_sweeps = 40 + 50 * 2
    assert diag.converged is False
    assert default_sweeps < diag.n_iterations <= 3 * default_sweeps
    assert rec.tau_mu_draws.shape == (200, 4)
    assert diag["n_iterations"] == diag.n_iterations


def test_iteration_cost_is_bounded_at_a5() -> None:
    """Median sweep count with auto-continuation stays within 3x of the default at A = 5."""
    ratios = []
    for seed in range(4):
        arms, fits = _chain_fits(5, 10 + seed, n_draws=400)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rec = longbet.reconcile_multi_arm_draws(fits, arms=arms, prior_scale=1.0, seed=seed)
        ratios.append(rec.sampler_diagnostics.n_iterations / (200.0 + 100 * 2))
    assert float(np.median(ratios)) <= 3.0


def _sbc_slice_ranks(
    n_arms: int,
    n_trials: int,
    n_post: int,
    seed: int,
    num_warmup: int = 100,
    max_iter_multiplier: float = 8.0,
) -> np.ndarray:
    """SBC on shared-arm summaries: mu* ~ bundle prior, y = B mu* + E eps, cut-space projection,
    auto-continued slice posterior; returns the normalized rank of the top contrast."""
    rng = np.random.default_rng(seed)
    dim = n_arms - 1
    pairs = mocks.all_pairs(n_arms)
    b_cut = np.zeros((len(pairs), dim))
    e_inc = np.zeros((len(pairs), n_arms))
    for p, (j, l) in enumerate(pairs):
        e_inc[p, j], e_inc[p, l] = 1.0, -1.0
    b_cut = e_inc[:, 1:]
    v_arm = np.linspace(0.02, 0.04, n_arms)
    sigma_lik = rkma._sigma_lik_from_arm_variances(v_arm)
    alpha = np.ones(dim)
    kw = dict(dim=dim, topology="chain", prior_mode="bundle_halfnormal", scale_val=0.6,
              alpha_vec=alpha, monotone=True, tau_clip=(0.0, np.inf),
              poset_edges_idx=[(a, a - 1) for a in range(1, n_arms)], non_top_arms=list(range(1, n_arms - 1)))

    def _p_draw(n: int) -> np.ndarray:
        return rkma._draw_prior_mu(rng, n=n, **kw)

    def _p_logpdf(pts: np.ndarray) -> np.ndarray:
        return rkma._prior_logpdf_mu(pts, **kw)

    ranks = np.zeros(n_trials)
    for t in range(n_trials):
        mu_star = _p_draw(1)[0]
        y = b_cut @ mu_star + e_inc @ (np.sqrt(v_arm) * rng.standard_normal(n_arms))
        mu_hat = np.linalg.solve(b_cut.T @ b_cut, b_cut.T @ y)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            post, _, _ = rkma._hit_and_run_slice_sampler(
                rng, prior_draw=_p_draw, prior_logpdf=_p_logpdf, mu_hat=mu_hat, sigma_lik=sigma_lik,
                n_draws=n_post, num_chains=4, num_warmup=num_warmup, thin=1, auto_continue=True,
                max_iter_multiplier=max_iter_multiplier,
            )
        ranks[t] = np.mean(post[:, -1] < mu_star[-1])
    return ranks


def _chi2_p(ranks: np.ndarray, n_bins: int) -> float:
    counts, _ = np.histogram(ranks, bins=np.linspace(0.0, 1.0, n_bins + 1))
    return float(stats.chisquare(counts).pvalue)


def test_sbc_auto_continued_slice_fast() -> None:
    ranks = _sbc_slice_ranks(5, n_trials=40, n_post=100, seed=505, num_warmup=60, max_iter_multiplier=3.0)
    assert _chi2_p(ranks, 4) > 0.01


@pytest.mark.slow
@pytest.mark.parametrize("n_arms", [5, 6])
def test_sbc_auto_continued_slice_slow(n_arms: int) -> None:
    ranks = _sbc_slice_ranks(n_arms, n_trials=300, n_post=200, seed=600 + n_arms)
    assert _chi2_p(ranks, 10) >= 0.05
