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

"""Unit tests for ``longbet._rank_k_multi`` and multi-outcome extensions to ``roi_break_even_table``."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import longbet
import longbet._rank_k as rk
import longbet._rank_k_multi as rkm
import longbet._rank_k_plots as rkp


def _make_mock_rank_k_result(
    *,
    att_mean: float,
    att_sd: float,
    n_draws: int = 800,
    n_units: int = 40,
    t_pre: int = 4,
    t_post: int = 12,
    rank_k: int = 2,
    mode2_mean: float = 0.0,
    mode2_sd: float = 0.02,
    seed: int = 0,
    rhat_att: float = 1.002,
    ess_att: float = 600.0,
) -> rk.RankKLongBetResult:
    """Build a lightweight synthetic :class:`RankKLongBetResult` for fast unit testing."""
    rng = np.random.default_rng(seed)
    n_periods = t_pre + t_post
    cohorts = np.full(n_units, n_periods, dtype=np.int32)
    cohorts[: n_units // 2] = t_pre
    ever = cohorts < n_periods
    unit_weights = np.where(ever, 1.0 / np.sum(ever), 0.0)

    att_draws = rng.normal(loc=att_mean, scale=att_sd, size=n_draws)
    mode1_draws = np.zeros((n_units, n_draws), dtype=np.float64)
    mode1_draws[ever] = att_draws[None, :] + rng.normal(0.0, 0.01, size=(int(np.sum(ever)), n_draws))
    # Force exact treated mean match
    mode1_draws[ever] += (att_draws - mode1_draws[ever].mean(axis=0))[None, :]

    mode_tau_draws = [mode1_draws]
    basis = rk.orthonormal_exposure_basis(t_post, rank_k)
    if rank_k >= 2:
        m2_trace = rng.normal(loc=mode2_mean, scale=mode2_sd, size=n_draws)
        mode2_draws = np.zeros((n_units, n_draws), dtype=np.float64)
        mode2_draws[ever] = m2_trace[None, :]
        mode_tau_draws.append(mode2_draws)

    catt_draws = np.zeros((n_draws, n_units, n_periods), dtype=np.float64)
    for k_idx, m_draws in enumerate(mode_tau_draws):
        psi_k = basis[:, k_idx]
        catt_draws[:, ever, t_pre:] += m_draws[ever].T[:, :, None] * psi_k[None, None, :]

    catt_hat = catt_draws.mean(axis=0)
    catt_ci_low = np.quantile(catt_draws, 0.025, axis=0)
    catt_ci_high = np.quantile(catt_draws, 0.975, axis=0)
    return rk.RankKLongBetResult(
        catt_draws=catt_draws,
        catt_hat=catt_hat,
        catt_ci_low=catt_ci_low,
        catt_ci_high=catt_ci_high,
        att_draws=att_draws,
        att_hat=float(np.mean(att_draws)),
        att_se=float(np.std(att_draws, ddof=1)),
        att_ci_low=float(np.quantile(att_draws, 0.025)),
        att_ci_high=float(np.quantile(att_draws, 0.975)),
        mode_tau_draws=mode_tau_draws,
        theta_f_modes=[1.0] * rank_k,
        sigma_modes=[att_sd] * rank_k,
        bases={t_pre: basis},
        runtime_sec=0.01,
        rhat_att=rhat_att,
        ess_att=ess_att,
        cohorts=cohorts,
        unit_weights=unit_weights,
    )


def test_public_exports_in_rank_k_and_top_level() -> None:
    expected = {
        "CrossOutcomeCorrelationResult",
        "ImanConoverResult",
        "Nested3ArmMultiReconciliationResult",
        "RankKConvergenceCheck",
        "check_rank_k_convergence",
        "estimate_cross_outcome_correlation",
        "ge_revenue_draws",
        "iman_conover_couple",
        "paid_winsor_cap",
        "reconcile_nested_3arm_multi_draws",
    }
    assert expected <= set(rkm.__all__)
    assert expected <= set(longbet.__all__)
    for name in expected:
        assert getattr(longbet, name) is getattr(rkm, name)
        assert getattr(rk, name) is getattr(rkm, name)


def test_estimate_cross_outcome_correlation_recovers_structure_and_floors_eigenvalues() -> None:
    rng = np.random.default_rng(42)
    n_units = 600
    x_covs = rng.normal(size=(n_units, 3))
    y_pre_1 = rng.normal(size=n_units)
    y_pre_2 = rng.normal(size=n_units)
    y_pre_3 = rng.normal(size=n_units)

    true_corr = np.array(
        [
            [1.0, 0.25, 0.65],
            [0.25, 1.0, 0.15],
            [0.65, 0.15, 1.0],
        ],
        dtype=np.float64,
    )
    chol = np.linalg.cholesky(true_corr)
    shocks = rng.normal(size=(n_units, 3)) @ chol.T
    arms = np.array(["C"] * 200 + ["B"] * 200 + ["A"] * 200)
    arm_shift = np.where(arms == "A", 0.8, np.where(arms == "B", 0.4, 0.0))

    y_post = {
        "outcome_1": 0.9 * y_pre_1 + 0.3 * x_covs[:, 0] + arm_shift + shocks[:, 0],
        "outcome_2": 0.7 * y_pre_2 - 0.2 * x_covs[:, 1] + 0.5 * arm_shift + shocks[:, 1],
        "outcome_3": 1.1 * y_pre_3 + 0.4 * x_covs[:, 2] + 1.2 * arm_shift + shocks[:, 2],
    }
    y_pre = {"outcome_1": y_pre_1, "outcome_2": y_pre_2, "outcome_3": y_pre_3}

    res = longbet.estimate_cross_outcome_correlation(
        y_post, y_pre, x_covs=x_covs, arm_labels=arms, min_eigenvalue=1e-4
    )
    assert res.outcomes == ("outcome_1", "outcome_2", "outcome_3")
    assert set(res.per_arm_correlation.keys()) == {"C", "B", "A"}
    assert res.correlation.shape == (3, 3)
    np.testing.assert_allclose(np.diag(res.correlation), 1.0, atol=1e-12)
    np.testing.assert_allclose(res.correlation, true_corr, atol=0.08)
    assert res.min_eigenvalue_after_floor >= 1e-4 - 1e-7
    assert all(0.0 <= lam <= 1.0 for lam in res.shrinkage_intensities.values())

    # Near-collinear small-sample stress test: shrinkage + eigenvalue floor keeps matrix SPD
    small_x = rng.normal(size=(12, 1))
    base = rng.normal(size=12)
    collinear_post = np.column_stack([base, base + 1e-6 * rng.normal(size=12), -base + 1e-6 * rng.normal(size=12)])
    res_small = longbet.estimate_cross_outcome_correlation(
        collinear_post, x_covs=small_x, outcomes=("m1", "m2", "m3"), min_eigenvalue=1e-3
    )
    eigvals = np.linalg.eigvalsh(res_small.correlation)
    assert float(np.min(eigvals)) >= 1e-3 - 1e-6
    np.testing.assert_allclose(np.diag(res_small.correlation), 1.0, atol=1e-12)

    with pytest.raises(ValueError, match="min_eigenvalue"):
        longbet.estimate_cross_outcome_correlation(collinear_post, min_eigenvalue=0.0)


def test_iman_conover_couple_preserves_exact_marginals_and_matches_spearman() -> None:
    rng = np.random.default_rng(123)
    n_draws = 2000
    draws = {
        "outcome_1": np.column_stack([rng.gamma(2.0, 0.5, size=n_draws), rng.gamma(1.5, 0.3, size=n_draws)]),
        "outcome_2": np.column_stack([rng.exponential(0.2, size=n_draws), rng.exponential(0.1, size=n_draws)]),
        "outcome_3": np.column_stack([rng.lognormal(0.0, 0.4, size=n_draws), rng.lognormal(-0.2, 0.3, size=n_draws)]),
    }
    target = np.array(
        [
            [1.0, 0.30, 0.65],
            [0.30, 1.0, 0.20],
            [0.65, 0.20, 1.0],
        ],
        dtype=np.float64,
    )
    coupled = longbet.iman_conover_couple(draws, target, key="sum", seed=77)
    assert coupled.outcomes == ("outcome_1", "outcome_2", "outcome_3")
    # Outcome 0 is anchored to identity permutation
    np.testing.assert_array_equal(coupled.permutations["outcome_1"], np.arange(n_draws))
    for name in draws:
        orig = draws[name]
        reordered = coupled.coupled_draws[name]
        assert reordered.shape == orig.shape
        # Exact marginal preservation: sorted columns are identical bit-for-bit
        np.testing.assert_array_equal(np.sort(reordered[:, 0]), np.sort(orig[:, 0]))
        np.testing.assert_array_equal(np.sort(reordered[:, 1]), np.sort(orig[:, 1]))
        np.testing.assert_array_equal(
            np.sort(reordered.sum(axis=1)), np.sort(orig.sum(axis=1))
        )
    np.testing.assert_allclose(coupled.realized_spearman, target, atol=0.025)

    # Test col0 key and callable key
    coupled_col0 = longbet.iman_conover_couple(draws, target, key="col0", seed=77)
    np.testing.assert_allclose(coupled_col0.realized_spearman, target, atol=0.025)
    coupled_fn = longbet.iman_conover_couple(draws, target, key=lambda a: a[:, 0], seed=77)
    np.testing.assert_allclose(coupled_fn.realized_spearman, coupled_col0.realized_spearman, atol=1e-12)


def test_reconcile_nested_3arm_multi_draws_bottom_up_exact_additivity_and_marginals() -> None:
    n_draws = 1000
    fits_bc = {
        "outcome_1": _make_mock_rank_k_result(att_mean=0.22, att_sd=0.06, n_draws=n_draws, rank_k=2, seed=11),
        "outcome_2": _make_mock_rank_k_result(att_mean=0.08, att_sd=0.04, n_draws=n_draws, rank_k=1, seed=12),
    }
    fits_ab = {
        "outcome_1": _make_mock_rank_k_result(att_mean=0.14, att_sd=0.05, n_draws=n_draws, rank_k=2, seed=21),
        "outcome_2": _make_mock_rank_k_result(att_mean=0.05, att_sd=0.03, n_draws=n_draws, rank_k=1, seed=22),
    }
    fits_ac = {
        "outcome_1": _make_mock_rank_k_result(att_mean=0.38, att_sd=0.06, n_draws=n_draws, rank_k=2, seed=31),
        "outcome_2": _make_mock_rank_k_result(att_mean=0.12, att_sd=0.04, n_draws=n_draws, rank_k=1, seed=32),
    }
    prior_scales = {"outcome_1": 0.30, "outcome_2": 0.15}
    target_corr = np.array([[1.0, 0.35], [0.35, 1.0]], dtype=np.float64)

    multi = longbet.reconcile_nested_3arm_multi_draws(
        fits_bc,
        fits_ab,
        fits_ac,
        prior_scales=prior_scales,
        cross_outcome_corr=target_corr,
        composite_outcomes={"total": ("outcome_1", "outcome_2")},
        composite_mode="bottom_up",
        n_proposals=60_000,
        seed=2026,
    )
    assert multi.outcomes == ("outcome_1", "outcome_2", "total")
    assert multi.structural_outcomes == ("outcome_1", "outcome_2")
    assert multi.composite_outcomes == {"total": ("outcome_1", "outcome_2")}
    assert multi.max_abs_composite_discrepancy["total"] < 1e-12
    assert multi.additivity_violation_rate["total"] == 0.0
    np.testing.assert_allclose(multi.realized_spearman, target_corr, atol=0.04)

    # Exact draw-wise additivity across components and across contrasts
    np.testing.assert_allclose(
        multi.post_att_bc["total"],
        multi.post_att_bc["outcome_1"] + multi.post_att_bc["outcome_2"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        multi.post_att_ab["total"],
        multi.post_att_ab["outcome_1"] + multi.post_att_ab["outcome_2"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        multi.post_att_ac["total"],
        multi.post_att_ac["outcome_1"] + multi.post_att_ac["outcome_2"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        multi.trt_catt_ac["total"],
        multi.trt_catt_ac["outcome_1"] + multi.trt_catt_ac["outcome_2"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        multi.cum_att_ac["total"],
        multi.cum_att_ac["outcome_1"] + multi.cum_att_ac["outcome_2"],
        atol=1e-12,
    )
    for name in multi.outcomes:
        np.testing.assert_allclose(
            multi.trt_catt_ac[name],
            multi.trt_catt_bc[name] + multi.trt_catt_ab[name],
            atol=1e-12,
        )

    # Exact marginal preservation vs single-outcome reconcile_nested_3arm_draws
    for idx, m in enumerate(("outcome_1", "outcome_2")):
        single = longbet.reconcile_nested_3arm_draws(
            fits_bc[m],
            fits_ab[m],
            fits_ac[m],
            prior_scale=prior_scales[m],
            n_proposals=60_000,
            seed=2026 + 1009 * idx,
        )
        np.testing.assert_allclose(
            np.sort(multi.post_att_ac[m]),
            np.sort(single.post_att_ac),
            atol=1e-14,
        )
        np.testing.assert_allclose(
            np.sort(multi.cum_att_ac[m]),
            np.sort(single.trt_catt_ac.sum(axis=1)),
            atol=1e-13,
        )


def test_reconcile_nested_3arm_multi_draws_coherent_total_preserves_null_floor_and_additivity() -> None:
    """Under a global null, coherent_total prevents sqrt(2) sum-of-half-normals inflation while keeping AVR == 0."""
    n_draws = 1200
    # Global null: all pairwise point estimates near 0 with SD = 0.10
    fits_bc = {
        "outcome_1": _make_mock_rank_k_result(att_mean=0.0, att_sd=0.10, n_draws=n_draws, rank_k=1, seed=101),
        "outcome_2": _make_mock_rank_k_result(att_mean=0.0, att_sd=0.10, n_draws=n_draws, rank_k=1, seed=102),
    }
    fits_ab = {
        "outcome_1": _make_mock_rank_k_result(att_mean=0.0, att_sd=0.10, n_draws=n_draws, rank_k=1, seed=201),
        "outcome_2": _make_mock_rank_k_result(att_mean=0.0, att_sd=0.10, n_draws=n_draws, rank_k=1, seed=202),
    }
    fits_ac = {
        "outcome_1": _make_mock_rank_k_result(att_mean=0.0, att_sd=0.10, n_draws=n_draws, rank_k=1, seed=301),
        "outcome_2": _make_mock_rank_k_result(att_mean=0.0, att_sd=0.10, n_draws=n_draws, rank_k=1, seed=302),
    }
    # Direct composite fits under the same global null (variance = 0.10^2 + 0.10^2 -> SD = sqrt(2)*0.10)
    sd_comp = float(np.sqrt(2.0) * 0.10)
    comp_fits = {
        "BC": _make_mock_rank_k_result(att_mean=0.0, att_sd=sd_comp, n_draws=n_draws, rank_k=1, seed=401),
        "AB": _make_mock_rank_k_result(att_mean=0.0, att_sd=sd_comp, n_draws=n_draws, rank_k=1, seed=402),
        "AC": _make_mock_rank_k_result(att_mean=0.0, att_sd=sd_comp, n_draws=n_draws, rank_k=1, seed=403),
    }
    direct_comp = longbet.reconcile_nested_3arm_draws(
        comp_fits["BC"],
        comp_fits["AB"],
        comp_fits["AC"],
        prior_scale=sd_comp,
        n_proposals=80_000,
        seed=999,
    )
    direct_null_mean = float(np.mean(direct_comp.post_att_ac))

    bottom_up = longbet.reconcile_nested_3arm_multi_draws(
        fits_bc,
        fits_ab,
        fits_ac,
        prior_scales={"outcome_1": 0.10, "outcome_2": 0.10},
        composite_outcomes={"total": ("outcome_1", "outcome_2")},
        composite_mode="bottom_up",
        n_proposals=80_000,
        seed=500,
    )
    coherent = longbet.reconcile_nested_3arm_multi_draws(
        fits_bc,
        fits_ab,
        fits_ac,
        prior_scales={"outcome_1": 0.10, "outcome_2": 0.10},
        composite_outcomes={"total": ("outcome_1", "outcome_2")},
        composite_mode="coherent_total",
        n_proposals=80_000,
        seed=500,
    )

    # Bottom-up inflates the null floor by ~sqrt(2) = 1.41x
    bu_ratio = float(np.mean(bottom_up.post_att_ac["total"])) / direct_null_mean
    coh_ratio = float(np.mean(coherent.post_att_ac["total"])) / direct_null_mean
    assert bu_ratio > 1.30
    assert coh_ratio <= 1.05

    # And coherent_total STILL satisfies exact draw-wise additivity (AVR == 0, gap < 1e-12)
    assert coherent.max_abs_composite_discrepancy["total"] < 1e-12
    assert coherent.additivity_violation_rate["total"] == 0.0
    np.testing.assert_allclose(
        coherent.post_att_ac["total"],
        coherent.post_att_ac["outcome_1"] + coherent.post_att_ac["outcome_2"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        coherent.trt_catt_ac["total"],
        coherent.trt_catt_ac["outcome_1"] + coherent.trt_catt_ac["outcome_2"],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        coherent.cum_att_ac["total"],
        coherent.cum_att_ac["outcome_1"] + coherent.cum_att_ac["outcome_2"],
        atol=1e-12,
    )

    # Also verify explicit composite_results passed into coherent_total
    coherent_explicit = longbet.reconcile_nested_3arm_multi_draws(
        fits_bc,
        fits_ab,
        fits_ac,
        prior_scales={"outcome_1": 0.10, "outcome_2": 0.10},
        composite_outcomes={"total": ("outcome_1", "outcome_2")},
        composite_mode="coherent_total",
        composite_results={"total": comp_fits},
        composite_prior_scales={"total": sd_comp},
        n_proposals=80_000,
        seed=500,
    )
    assert coherent_explicit.max_abs_composite_discrepancy["total"] < 1e-12
    assert coherent_explicit.additivity_violation_rate["total"] == 0.0
    assert float(np.mean(coherent_explicit.post_att_ac["total"])) / direct_null_mean <= 1.05

    # Overlapping composite component sets under coherent_total must raise ValueError
    with pytest.raises(ValueError, match="disjoint composite component sets"):
        longbet.reconcile_nested_3arm_multi_draws(
            fits_bc,
            fits_ab,
            fits_ac,
            prior_scales={"outcome_1": 0.10, "outcome_2": 0.10},
            composite_outcomes={"total_1": ("outcome_1", "outcome_2"), "total_2": ("outcome_1", "outcome_2")},
            composite_mode="coherent_total",
            n_proposals=10_000,
            seed=500,
        )


def test_ge_revenue_draws_and_vector_alpha_roi_break_even_table() -> None:
    rng = np.random.default_rng(314)
    n_draws = 1000
    cum_1 = rng.uniform(1.0, 3.0, size=n_draws)
    cum_2 = rng.uniform(0.5, 1.5, size=n_draws)
    cum_total = cum_1 + cum_2

    # ge_revenue_draws on 1D cumulative draws
    ge_half = longbet.ge_revenue_draws(
        {"outcome_1": cum_1, "outcome_2": cum_2},
        weights={"outcome_1": 0.5, "outcome_2": 1.0},
    )
    np.testing.assert_allclose(ge_half, 0.5 * cum_1 + 1.0 * cum_2, atol=1e-12)

    # ge_revenue_draws on 2D trajectories with horizon and discount_rate
    traj_1 = np.tile(cum_1[:, None] / 12.0, (1, 12))
    traj_2 = np.tile(cum_2[:, None] / 12.0, (1, 12))
    ge_2d_undisc = longbet.ge_revenue_draws(
        {"outcome_1": traj_1, "outcome_2": traj_2},
        weights={"outcome_1": 0.5, "outcome_2": 1.0},
        horizon=12,
        discount_rate=0.0,
    )
    np.testing.assert_allclose(ge_2d_undisc, ge_half, atol=1e-12)
    ge_2d_disc = longbet.ge_revenue_draws(
        {"outcome_1": traj_1, "outcome_2": traj_2},
        weights={"outcome_1": 0.5, "outcome_2": 1.0},
        horizon=6,
        discount_rate=0.10,
    )
    assert ge_2d_disc.shape == (n_draws,)
    assert np.all(ge_2d_disc < ge_2d_undisc)

    # When alpha_1 == alpha_2 == alpha, multi-outcome roi_break_even_table matches single-outcome on total
    alphas = (0.25, 0.50, 0.75, 1.00)
    costs = (1.5, 2.5)
    single_df = rkp.roi_break_even_table(
        {"AC": cum_total},
        cost_grid=costs,
        alpha_grid=alphas,
        horizon_periods=12,
    )
    uniform_weights = [{"outcome_1": a, "outcome_2": a} for a in alphas]
    multi_uniform_df = rkp.roi_break_even_table(
        cost_grid=costs,
        alpha_grid=alphas,
        horizon_periods=12,
        component_cum_draws={"outcome_1": cum_1, "outcome_2": cum_2},
        component_weights=uniform_weights,
    )
    np.testing.assert_allclose(
        multi_uniform_df["pr_roi_positive"].to_numpy(),
        single_df["pr_roi_positive"].to_numpy(),
        atol=1e-12,
    )
    # When non-focal weights vary across rows (outcome_2 weight = a), break_even_alpha_50 reflects each row's weight
    c25_rows = multi_uniform_df[multi_uniform_df["cost"] == 2.5]
    assert len(set(np.round(c25_rows["break_even_alpha_50"].to_numpy(), 6))) > 1

    # With alpha_2 = 1.0 fixed and alpha_1 varying over alphas:
    vec_df = rkp.roi_break_even_table(
        cost_grid=costs,
        alpha_grid=alphas,
        horizon_periods=12,
        component_cum_draws={"outcome_1": cum_1, "outcome_2": cum_2},
    )
    assert list(vec_df.columns) == [
        "tier",
        "contrast",
        "horizon_periods",
        "cost",
        "alpha",
        "pr_roi_positive",
        "break_even_alpha_50",
    ]
    for row in vec_df.itertuples(index=False):
        expected_prob = float(np.mean(row.alpha * cum_1 + 1.0 * cum_2 > row.cost))
        assert row.pr_roi_positive == pytest.approx(expected_prob, abs=1e-12)
        if row.cost == 1.5:
            med_at_be = float(np.median(row.break_even_alpha_50 * cum_1 + cum_2))
            assert med_at_be == pytest.approx(row.cost, abs=1e-6)


def test_check_rank_k_convergence_and_paid_winsor_cap() -> None:
    rng = np.random.default_rng(99)
    good_res = _make_mock_rank_k_result(
        att_mean=0.25, att_sd=0.05, n_draws=800, rank_k=2, seed=55, rhat_att=1.001, ess_att=500.0
    )
    chk_ok = longbet.check_rank_k_convergence(good_res, num_chains=4)
    assert chk_ok.converged is True
    assert chk_ok.failed_modes == ()
    assert len(chk_ok.rhat_modes) == 2 and len(chk_ok.ess_modes) == 2

    # Corrupt mode 1 (index 1, i.e. second mode) with divergent chains (large shift between halves + noise)
    bad_mode2 = rng.normal(0.0, 0.1, size=good_res.mode_tau_draws[1].shape)
    bad_mode2[:, :400] -= 5.0
    bad_mode2[:, 400:] += 5.0
    good_res.mode_tau_draws[1] = bad_mode2
    chk_bad = longbet.check_rank_k_convergence(good_res, num_chains=4, rhat_max=1.05, ess_min=100.0)
    assert chk_bad.converged is False
    assert 1 in chk_bad.failed_modes

    # paid_winsor_cap on zero-inflated pre-period data
    y_pre = np.zeros((200, 12), dtype=np.float64)
    pos_vals = rng.lognormal(mean=1.5, sigma=0.8, size=300)
    y_pre.flat[:300] = pos_vals
    cap = longbet.paid_winsor_cap(y_pre, q=0.99)
    assert cap == pytest.approx(float(np.quantile(pos_vals, 0.99)), rel=1e-12)
    assert longbet.paid_winsor_cap(np.zeros((20, 5))) == 0.0
    with pytest.raises(ValueError, match="q must lie"):
        longbet.paid_winsor_cap(y_pre, q=0.0)
