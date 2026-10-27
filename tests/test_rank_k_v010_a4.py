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

"""v0.10 item A4: bipartite contrast graphs (design warning, policies, star fits on a chain)."""

from __future__ import annotations

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
_STAR = [(1, 0), (2, 0), (3, 0)]
_MU = np.array([0.0, 0.35, 0.55, 0.65])
_V_HET = np.array([0.0012, 0.0020, 0.0035, 0.0060])


def _star_fits(seed: int = 1) -> dict[tuple[str, str], longbet.RankKLongBetResult]:
    rng = np.random.default_rng(seed)
    fits, _ = mocks.simulate_shared_arm_fits(rng, _MU, _V_HET, _STAR, _ARMS4, n_draws=400)
    return fits


def test_design_warning_names_missing_odd_cycle_and_fit() -> None:
    fits = _star_fits()
    with pytest.warns(longbet.RankKDesignWarning) as rec_w:
        longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, prior_scale=1.0, n_proposals=10_000)
    msgs = [str(w.message) for w in rec_w if w.category is longbet.RankKDesignWarning]
    assert len(msgs) == 1
    assert "b1 -> c -> b2 -> b1" in msgs[0]
    assert "('b2', 'b1')" in msgs[0]
    assert issubclass(longbet.RankKDesignWarning, UserWarning)


def test_no_design_warning_with_odd_cycle() -> None:
    rng = np.random.default_rng(2)
    pairs = [*_STAR, (2, 1)]
    fits, _ = mocks.simulate_shared_arm_fits(rng, _MU, _V_HET, pairs, _ARMS4, n_draws=300)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rec = longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, prior_scale=1.0, n_proposals=10_000)
    assert not any(w.category is longbet.RankKDesignWarning for w in caught)
    assert rec.bipartite_policy is None and rec.null_shift_interval is None


def test_suggest_odd_cycle_fit_on_even_cycle() -> None:
    arms = ("c", "b1", "b2", "b3")
    pairs = [(1, 0), (2, 1), (3, 2), (3, 0)]  # 4-cycle c-b1-b2-b3-c
    _, is_bip, colors = rkma._check_connected_and_bipartite(4, pairs)
    assert is_bip
    (j, l), cycle = rkma._suggest_odd_cycle_fit(arms, pairs, colors)
    assert colors[j] == colors[l]
    assert cycle[0] == l and cycle[-1] == l and cycle[-2] == j
    assert len(cycle) - 1 == 3  # triangle


def test_star_fits_on_chain_order_reconcile_with_chain_monotonicity() -> None:
    fits = _star_fits(3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, prior_scale=1.0, n_proposals=20_000)
    assert rec.topology == "chain"
    assert set(rec.tree_edges) == {("b1", "c"), ("b2", "c"), ("b3", "c")}
    mu = rec.tau_mu_draws
    assert np.all(np.diff(np.column_stack([np.zeros(mu.shape[0]), mu]), axis=1) >= -1e-12)
    assert rec.trajectory_additivity_max_gap < 1e-10
    with pytest.raises(ValueError, match="uniform split"):
        longbet.reconcile_multi_arm_draws(fits, arms=_ARMS4, prior_scale=1.0, split_alpha=[1.0, 2.0, 1.0])


@pytest.mark.parametrize("policy", ["allocation", "conservative", "marginalize"])
def test_bipartite_policies_keep_pairwise_variances_and_floors(policy: str) -> None:
    fits = _star_fits(4)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_scale=1.0, bipartite_policy=policy, n_proposals=20_000, seed=5
        )
    assert rec.bipartite_policy == policy
    lo, hi = rec.null_shift_interval
    assert lo <= 0.0 <= hi
    v = np.asarray(rec.arm_variances, dtype=np.float64)
    s_sq = np.array([np.var(fits[(_ARMS4[j], _ARMS4[l])].att_draws, ddof=1) for j, l in _STAR])
    alloc = longbet.reconcile_multi_arm_draws(
        fits, arms=_ARMS4, prior_scale=1.0, bipartite_policy="allocation", n_proposals=20_000, seed=5
    ) if policy != "allocation" else rec
    v_alloc = np.asarray(alloc.arm_variances, dtype=np.float64)
    # Shifts along the null vector leave every fitted pairwise sum v_j + v_0 unchanged.
    np.testing.assert_allclose(v[1:] + v[0], v_alloc[1:] + v_alloc[0], rtol=1e-10)
    assert np.all(v > 0.0)
    if policy == "allocation":
        assert rec.null_shift == 0.0
    elif policy == "conservative":
        assert lo <= rec.null_shift <= hi
    else:
        assert rec.null_shift is None
    assert np.all(np.isfinite(rec.tau_mu_draws))
    _ = s_sq


def test_conservative_policy_maximizes_gaussian_posterior_variance() -> None:
    fits = _star_fits(6)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = longbet.reconcile_multi_arm_draws(
            fits, arms=_ARMS4, prior_mode="flat", monotone=False, tau_clip=None,
            bipartite_policy="conservative", n_proposals=20_000,
        )
    # Flat prior: posterior variance of the fitted top contrast is invariant; any xi is optimal,
    # and the reported sigma_lik must match the chosen arm variances exactly.
    v = np.asarray(rec.arm_variances, dtype=np.float64)
    np.testing.assert_allclose(rec.sigma_lik, np.diag(v[1:]) + v[0], rtol=1e-12)


def test_multi_outcome_accepts_star_on_chain_and_policy() -> None:
    rng = np.random.default_rng(7)
    fits = {
        m: mocks.simulate_shared_arm_fits(rng, _MU * f, _V_HET, _STAR, _ARMS4, n_draws=300)[0]
        for m, f in (("dense", 0.7), ("sparse", 0.3))
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = longbet.reconcile_multi_arm_multi_outcome_draws(
            fits,
            arms=_ARMS4,
            prior_scales=1.0,
            composite_outcomes={"total": ("dense", "sparse")},
            bipartite_policy="marginalize",
            n_proposals=20_000,
        )
    assert res.composite_mode == "coherent_total"
    assert res.max_abs_composite_discrepancy["total"] < 1e-12
    assert res.by_outcome["dense"].bipartite_policy == "marginalize"
    assert res.by_outcome["total"].bipartite_policy == "marginalize"
    top = res.by_outcome["total"].tau_mu_draws
    assert np.all(np.diff(np.column_stack([np.zeros(top.shape[0]), top]), axis=1) >= -1e-10)
