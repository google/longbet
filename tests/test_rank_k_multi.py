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

"""Tests for v0.11 multi-outcome Bayesian Rank-K reconciliation (`_rank_k_multi.py`)."""

from __future__ import annotations

import numpy as np
import pytest

from longbet import (
    CrossOutcomeCorrelationResult,
    MultiArmMultiOutcomeReconciliationResult,
    check_rank_k_convergence,
    estimate_cross_outcome_correlation,
    ge_revenue_draws,
    multi_arm_decision_table,
    paid_winsor_cap,
    reconcile_multi_arm_multi_outcome_draws,
    reconcile_nested_3arm_multi_draws,
    roi_break_even_table,
)
from tests.test_rank_k_multi_arm import _make_stage1_fit


def test_utilities_winsor_convergence_ge_revenue_and_iw_correlation() -> None:
    rng = np.random.default_rng(11)
    y_pre = rng.exponential(scale=2.0, size=(50, 6))
    cap = paid_winsor_cap(y_pre, q=0.95)
    assert cap > 0.0
    assert paid_winsor_cap(np.zeros((10, 4))) == 0.0
    with pytest.raises(ValueError, match="q must lie"):
        paid_winsor_cap(y_pre, q=1.0)

    fit = _make_stage1_fit(rng, np.linspace(0.4, 0.6, 5))
    chk = check_rank_k_convergence(fit)
    assert chk.converged and len(chk.rhat_modes) == 2

    d1 = rng.normal(1.0, 0.1, size=(100, 4))
    d2 = rng.normal(0.5, 0.1, size=(100, 4))
    rev = ge_revenue_draws({"paid": d1, "organic": d2}, weights={"paid": 0.6, "organic": 1.0}, discount_rate=0.10)
    assert rev.shape == (100,) and np.all(np.isfinite(rev))

    n_u = 120
    cov_true = np.array([[1.0, 0.55], [0.55, 1.0]])
    raw_u = rng.multivariate_normal([0.0, 0.0], cov_true, size=n_u)
    corr_res = estimate_cross_outcome_correlation(
        {"paid": raw_u[:, 0], "organic": raw_u[:, 1]},
        arm_labels=np.repeat(["C", "B", "A"], 40),
    )
    assert isinstance(corr_res, CrossOutcomeCorrelationResult)
    assert corr_res.outcomes == ("paid", "organic")
    assert corr_res.correlation.shape == (2, 2)
    assert 0.25 < corr_res.correlation[0, 1] < 0.85
    assert set(corr_res.per_arm_correlation) == {"C", "B", "A"}


def test_multi_outcome_reconciliation_with_composite_and_iw_correlation() -> None:
    rng = np.random.default_rng(2026)
    h_post = 5
    paid_b = np.full(h_post, 0.45)
    paid_a = np.full(h_post, 0.95)
    org_b = np.full(h_post, 0.10)
    org_a = np.full(h_post, -0.15)  # Substitution on organic

    results_by_outcome = {
        "paid": {
            ("B", "C"): _make_stage1_fit(rng, paid_b),
            ("A", "B"): _make_stage1_fit(rng, paid_a - paid_b),
            ("A", "C"): _make_stage1_fit(rng, paid_a),
        },
        "organic": {
            ("B", "C"): _make_stage1_fit(rng, org_b),
            ("A", "B"): _make_stage1_fit(rng, org_a - org_b),
            ("A", "C"): _make_stage1_fit(rng, org_a),
        },
    }
    composite_results = {
        "total": {
            ("B", "C"): _make_stage1_fit(rng, paid_b + org_b),
            ("A", "B"): _make_stage1_fit(rng, (paid_a + org_a) - (paid_b + org_b)),
            ("A", "C"): _make_stage1_fit(rng, paid_a + org_a),
        }
    }
    resid = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.4], [0.4, 1.0]], size=100)

    mres = reconcile_multi_arm_multi_outcome_draws(
        results_by_outcome,
        composite_outcomes={"total": ("paid", "organic")},
        composite_results=composite_results,
        cross_outcome_residuals=resid,
        arms=("C", "B", "A"),
        order=("C", "B", "A"),
        monotone={"paid": True, "organic": False},
        seed=88,
    )
    assert isinstance(mres, MultiArmMultiOutcomeReconciliationResult)
    assert mres.outcomes == ("paid", "organic", "total")
    assert mres.structural_outcomes == ("paid", "organic")
    assert mres.additivity_violation_rate["total"] == 0.0
    assert mres.max_abs_composite_discrepancy["total"] < 1e-13
    assert mres.sampler_diagnostics.converged
    assert mres.sampler_diagnostics.r_hat <= 1.01
    assert mres.sampler_diagnostics.ess >= 400.0
    assert mres.posterior_correlation_draws.shape == (1600, 2, 2)
    assert mres.posterior_correlation_mean.shape == (2, 2)
    assert mres.ppc_composite_pvalue["total"] is not None

    # Paid is monotone-constrained; organic is unconstrained and recovers negative A vs B
    assert np.min(mres.contrast_trajectory_draws["paid"][("A", "B")]) >= -1e-12
    assert float(np.mean(mres.contrast_mode1_draws["organic"][("A", "B")])) < 0.0

    # Convenience 3-arm multi-outcome wrapper delegates to identical engine
    res_bc = {m: results_by_outcome[m][("B", "C")] for m in ("paid", "organic")}
    res_ab = {m: results_by_outcome[m][("A", "B")] for m in ("paid", "organic")}
    res_ac = {m: results_by_outcome[m][("A", "C")] for m in ("paid", "organic")}
    comp_3arm = {"total": {"BC": composite_results["total"][("B", "C")], "AB": composite_results["total"][("A", "B")], "AC": composite_results["total"][("A", "C")]}}
    wres = reconcile_nested_3arm_multi_draws(
        res_bc,
        res_ab,
        res_ac,
        composite_outcomes={"total": ("paid", "organic")},
        composite_results=comp_3arm,
        cross_outcome_residuals=resid,
        monotone={"paid": True, "organic": False},
        seed=88,
    )
    np.testing.assert_allclose(
        wres.contrast_trajectory_draws["total"][("A", "C")],
        mres.contrast_trajectory_draws["total"][("A", "C")],
        atol=1e-13,
    )

    # Decision table & ROI break-even table on multi-outcome result
    dt = multi_arm_decision_table(
        mres,
        costs={"C": 0.0, "B": 0.5, "A": 1.0},
        alpha={"paid": 0.7, "organic": 1.0},
    )
    assert list(dt["arm"]) == ["C", "B", "A"]
    roi_df = roi_break_even_table(mres, cost_grid=(0.5, 1.0), alpha_grid=(0.5, 1.0))
    assert not roi_df.empty and not roi_df.isna().any().any()
