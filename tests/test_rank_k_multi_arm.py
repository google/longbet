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


def test_bipartite_star_with_arm_metadata_is_quiet() -> None:
    rng = np.random.default_rng(304)
    h_post = 5
    arms = ("C", "T1", "T2", "T3")
    poset = (("T1", "C"), ("T2", "C"), ("T3", "C"))
    results = {
        ("T1", "C"): _make_stage1_fit(rng, np.full(h_post, 0.3), include_arm_metadata=True),
        ("T2", "C"): _make_stage1_fit(rng, np.full(h_post, 0.5), include_arm_metadata=True),
        ("T3", "C"): _make_stage1_fit(rng, np.full(h_post, 0.7), include_arm_metadata=True),
    }
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rec = reconcile_multi_arm_draws(results, arms=arms, poset=poset, seed=42)
    design_warnings = [w for w in caught if issubclass(w.category, RankKDesignWarning)]
    assert not design_warnings
    assert rec.sampler_diagnostics.converged


def test_mixed_rank_stage1_partial_mode_observation_preserves_dynamic_slope() -> None:
    rng = np.random.default_rng(505)
    h_post = 6
    traj_b = np.linspace(0.25, 0.75, h_post)
    traj_a = np.linspace(0.60, 1.30, h_post)

    # Fit ("B", "C") at rank_k=1 (constant draws in Mode 0 only), while ("A", "B") and ("A", "C") are rank_k=2
    fit_bc = _make_stage1_fit(rng, traj_b, sigma_arm_trt=0.03, sigma_arm_ctrl=0.03, sigma_forest=0.02)
    bc_mean_draws = np.mean(fit_bc.arm_metadata["trt_trajectory_draws"], axis=1, keepdims=True)
    bc_rank1_draws = np.repeat(bc_mean_draws, h_post, axis=1)
    fit_bc.arm_metadata["trt_trajectory_draws"] = bc_rank1_draws
    fit_bc.catt_draws[:, :50, 4:] = bc_rank1_draws[:, None, :]
    fit_bc.mode_tau_draws = (fit_bc.mode_tau_draws[0],)
    fit_bc.rank_k = 1

    fit_ab = _make_stage1_fit(rng, traj_a - traj_b, sigma_arm_trt=0.03, sigma_arm_ctrl=0.03, sigma_forest=0.02)
    fit_ac = _make_stage1_fit(rng, traj_a, sigma_arm_trt=0.03, sigma_arm_ctrl=0.03, sigma_forest=0.02)
    fit_ab.rank_k = 2
    fit_ac.rank_k = 2

    results = {("B", "C"): fit_bc, ("A", "B"): fit_ab, ("A", "C"): fit_ac}
    rec = reconcile_multi_arm_draws(results, arms=("C", "B", "A"), order=("C", "B", "A"), n_modes=2, seed=91)
    assert rec.sampler_diagnostics.converged
    assert rec.max_cycle_discrepancy < 1e-13
    assert rec.ppc_cycle_pvalue is not None and 0.0 <= rec.ppc_cycle_pvalue <= 1.0

    # Dynamic slope on ("B", "C") is recovered via ("A", "C") - ("A", "B") instead of clamped to zero!
    bc_post_mean = np.mean(rec.contrast_trajectory_draws[("B", "C")], axis=0)
    assert float(bc_post_mean[-1] - bc_post_mean[0]) > 0.20


def test_basis_ppc_returns_none_when_draws_lie_in_col_psi_and_s0_is_leak_free() -> None:
    rng = np.random.default_rng(606)
    h_post = 6
    psi = orthonormal_exposure_basis(h_post, 2) / np.sqrt(float(h_post))

    def _in_subspace_fit(level: float, post_noise_mult: float = 1.0) -> types.SimpleNamespace:
        fit = _make_stage1_fit(
            rng,
            np.linspace(level, level + 0.2, h_post),
            sigma_arm_trt=0.04 * post_noise_mult,
            sigma_arm_ctrl=0.04 * post_noise_mult,
            sigma_forest=0.03 * post_noise_mult,
        )
        # Fix pre-period variances to an exact constant independent of post_noise_mult
        fit.arm_metadata["pre_residual_var_treated"] = 0.25
        fit.arm_metadata["pre_residual_var_control"] = 0.25
        # Project draws strictly onto col(Psi)
        raw_d = fit.arm_metadata["trt_trajectory_draws"]
        in_col = (raw_d @ psi) @ psi.T
        fit.arm_metadata["trt_trajectory_draws"] = in_col
        fit.catt_draws[:, :50, 4:] = in_col[:, None, :]
        return fit

    res_low_noise = {
        ("B", "C"): _in_subspace_fit(0.3, 1.0),
        ("A", "B"): _in_subspace_fit(0.4, 1.0),
        ("A", "C"): _in_subspace_fit(0.7, 1.0),
    }
    res_high_post_noise = {
        ("B", "C"): _in_subspace_fit(0.3, 50.0),
        ("A", "B"): _in_subspace_fit(0.4, 50.0),
        ("A", "C"): _in_subspace_fit(0.7, 50.0),
    }
    rec_lo = reconcile_multi_arm_draws(res_low_noise, arms=("C", "B", "A"), order=("C", "B", "A"), seed=12)
    rec_hi = reconcile_multi_arm_draws(res_high_post_noise, arms=("C", "B", "A"), order=("C", "B", "A"), seed=12)
    assert rec_lo.ppc_basis_pvalue is None
    assert rec_lo.s0_scale == rec_hi.s0_scale == pytest.approx(0.5)


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
