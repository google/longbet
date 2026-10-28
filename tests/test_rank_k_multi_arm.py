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

"""Tests for v0.11 single-outcome multi-arm Bayesian Rank-K reconciliation (`_rank_k_multi_arm.py`)."""

from __future__ import annotations

import types
import warnings

import numpy as np
import pytest

from longbet import (
    MultiArmReconciliationResult,
    Nested3ArmReconciliationResult,
    RankKDesignWarning,
    multi_arm_decision_table,
    orthonormal_exposure_basis,
    reconcile_multi_arm_draws,
    reconcile_nested_3arm_draws,
    roi_break_even_table,
)


def _make_stage1_fit(
    rng: np.random.Generator,
    true_traj: np.ndarray,
    *,
    n_draws: int = 240,
    n_trt: int = 50,
    n_ctrl: int = 50,
    t_pre: int = 4,
    sigma_arm_trt: float = 0.05,
    sigma_arm_ctrl: float = 0.05,
    sigma_forest: float = 0.04,
    include_arm_metadata: bool = True,
) -> types.SimpleNamespace:
    """Create a synthetic Stage-1 RankKLongBetResult mock with realistic arm_metadata."""
    h_post = int(true_traj.shape[0])
    t_total = t_pre + h_post
    n_units = n_trt + n_ctrl
    cohorts = np.full(n_units, t_total, dtype=np.int32)
    cohorts[:n_trt] = t_pre

    cov_trt = (sigma_arm_trt**2) * np.eye(h_post, dtype=np.float64)
    cov_ctrl = (sigma_arm_ctrl**2) * np.eye(h_post, dtype=np.float64)
    cov_total = cov_trt + cov_ctrl + (sigma_forest**2) * np.eye(h_post, dtype=np.float64)

    mean_shift = rng.multivariate_normal(np.zeros(h_post), cov_total)
    center_traj = true_traj + mean_shift
    draws_post = rng.multivariate_normal(center_traj, cov_total, size=n_draws)

    catt_draws = np.zeros((n_draws, n_units, t_total), dtype=np.float64)
    catt_draws[:, :n_trt, t_pre:] = draws_post[:, None, :]
    att_draws = np.mean(draws_post, axis=1)

    basis = orthonormal_exposure_basis(h_post, min(2, h_post))
    mode0 = np.tile(att_draws[None, :], (n_units, 1))
    mode0[n_trt:, :] = 0.0
    if h_post >= 2:
        slope_draws = (draws_post @ basis[:, 1]) / float(h_post)
        mode1 = np.tile(slope_draws[None, :], (n_units, 1))
        mode1[n_trt:, :] = 0.0
        mode_tau_draws = (mode0, mode1)
    else:
        mode_tau_draws = (mode0,)

    s_sq_val = float((sigma_arm_trt**2) * n_trt + (sigma_arm_ctrl**2) * n_ctrl) / 2.0
    arm_meta = (
        {
            "n_treated": n_trt,
            "n_control": n_ctrl,
            "post_residual_cov_treated": cov_trt,
            "post_residual_cov_control": cov_ctrl,
            "pre_residual_var_treated": float(sigma_arm_trt**2 * n_trt),
            "pre_residual_var_control": float(sigma_arm_ctrl**2 * n_ctrl),
            "trt_trajectory_draws": draws_post,
        }
        if include_arm_metadata
        else None
    )

    return types.SimpleNamespace(
        att_draws=att_draws,
        att_hat=float(np.mean(att_draws)),
        att_se=float(np.std(att_draws, ddof=1)),
        catt_draws=catt_draws,
        cohorts=cohorts,
        bases={t_pre: basis},
        mode_tau_draws=mode_tau_draws,
        s_sq=np.full(n_draws, max(s_sq_val, 1e-4), dtype=np.float64),
        rhat_att=1.001,
        ess_att=float(n_draws),
        arm_metadata=arm_meta,
    )


def test_chain_3arm_zero_kwarg_and_nested_wrapper_equivalence() -> None:
    rng = np.random.default_rng(101)
    h_post = 6
    traj_b = np.linspace(0.35, 0.45, h_post)
    traj_a = np.linspace(0.80, 0.95, h_post)
    results = {
        ("B", "C"): _make_stage1_fit(rng, traj_b),
        ("A", "B"): _make_stage1_fit(rng, traj_a - traj_b),
        ("A", "C"): _make_stage1_fit(rng, traj_a),
    }

    rec = reconcile_multi_arm_draws(results, arms=("C", "B", "A"), order=("C", "B", "A"), seed=2026)
    assert isinstance(rec, MultiArmReconciliationResult)
    assert rec.arms == ("C", "B", "A")
    assert rec.topology == "chain"
    assert rec.max_cycle_discrepancy < 1e-13
    assert rec.sampler_diagnostics.converged
    assert rec.sampler_diagnostics.r_hat <= 1.01
    assert rec.sampler_diagnostics.ess >= 400.0
    assert rec.ppc_cycle_pvalue is not None and 0.0 <= rec.ppc_cycle_pvalue <= 1.0
    assert rec.ppc_basis_pvalue is not None and 0.0 <= rec.ppc_basis_pvalue <= 1.0

    # Every cover-edge trajectory is non-negative under default monotone=True, period_positivity=True
    assert np.min(rec.contrast_trajectory_draws[("B", "C")]) >= -1e-12
    assert np.min(rec.contrast_trajectory_draws[("A", "B")]) >= -1e-12
    np.testing.assert_allclose(
        rec.contrast_trajectory_draws[("A", "C")],
        rec.contrast_trajectory_draws[("B", "C")] + rec.contrast_trajectory_draws[("A", "B")],
        atol=1e-13,
    )

    # Convenience wrapper delegates to the exact same engine
    nested = reconcile_nested_3arm_draws(
        results[("B", "C")], results[("A", "B")], results[("A", "C")], seed=2026
    )
    assert isinstance(nested, Nested3ArmReconciliationResult)
    np.testing.assert_allclose(nested.trt_catt_bc, rec.trt_catt_bc, atol=1e-13)
    np.testing.assert_allclose(nested.trt_catt_ab, rec.trt_catt_ab, atol=1e-13)
    np.testing.assert_allclose(nested.trt_catt_ac, rec.trt_catt_ac, atol=1e-13)
    bc, ab, ac = nested
    np.testing.assert_allclose(ac, bc + ab, atol=1e-13)


def test_diamond_poset_4arm_reconciliation_and_decision_table() -> None:
    rng = np.random.default_rng(202)
    h_post = 5
    arms = ("00", "10", "01", "11")
    poset = (("11", "10"), ("11", "01"), ("10", "00"), ("01", "00"))
    true_arm = {
        "00": np.zeros(h_post),
        "10": np.full(h_post, 0.40),
        "01": np.full(h_post, 0.35),
        "11": np.full(h_post, 0.85),
    }
    pairs = [("10", "00"), ("01", "00"), ("11", "10"), ("11", "01"), ("11", "00")]
    results = {
        (hi, lo): _make_stage1_fit(rng, true_arm[hi] - true_arm[lo])
        for hi, lo in pairs
    }

    rec = reconcile_multi_arm_draws(results, arms=arms, poset=poset, seed=77)
    assert rec.topology == "poset"
    assert rec.max_cycle_discrepancy < 1e-13
    assert rec.sampler_diagnostics.converged
    for edge in poset:
        assert np.min(rec.contrast_trajectory_draws[edge]) >= -1e-12

    # Decision table & ROI break-even table
    dt = multi_arm_decision_table(rec, costs={"00": 0.0, "10": 0.5, "01": 0.5, "11": 1.0})
    assert list(dt["arm"]) == list(arms)
    assert np.isclose(dt["pr_best_arm"].sum(), 1.0, atol=1e-6)
    roi_df = roi_break_even_table(rec, cost_grid=(0.5, 1.0), alpha_grid=(0.5, 1.0))
    assert not roi_df.empty and not roi_df.isna().any().any()


def test_bipartite_star_without_arm_metadata_emits_design_warning() -> None:
    rng = np.random.default_rng(303)
    h_post = 5
    arms = ("C", "T1", "T2", "T3")
    poset = (("T1", "C"), ("T2", "C"), ("T3", "C"))
    results = {
        ("T1", "C"): _make_stage1_fit(rng, np.full(h_post, 0.3), include_arm_metadata=False),
        ("T2", "C"): _make_stage1_fit(rng, np.full(h_post, 0.5), include_arm_metadata=False),
        ("T3", "C"): _make_stage1_fit(rng, np.full(h_post, 0.7), include_arm_metadata=False),
    }
    with pytest.warns(RankKDesignWarning, match="bipartite"):
        rec = reconcile_multi_arm_draws(results, arms=arms, poset=poset, seed=42)
    assert rec.sampler_diagnostics.converged
    assert rec.max_cycle_discrepancy < 1e-13


def test_unconstrained_and_horizon_truncation() -> None:
    rng = np.random.default_rng(404)
    h_post = 6
    results = {
        ("B", "C"): _make_stage1_fit(rng, np.full(h_post, -0.25)),
        ("A", "B"): _make_stage1_fit(rng, np.full(h_post, 0.60)),
        ("A", "C"): _make_stage1_fit(rng, np.full(h_post, 0.35)),
    }
    rec = reconcile_multi_arm_draws(
        results, arms=("C", "B", "A"), order=("C", "B", "A"), monotone=False, horizon=4, seed=55
    )
    assert rec.horizon == 4
    assert rec.contrast_trajectory_draws[("B", "C")].shape[1] == 4
    assert float(np.mean(rec.contrast_mode1_draws[("B", "C")])) < 0.0
    assert rec.max_cycle_discrepancy < 1e-13
