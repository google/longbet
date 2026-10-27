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

"""Unit tests for ``longbet._rank_k_multi_arm`` and multi-arm extensions in ``_rank_k_multi`` / ``_rank_k_plots``."""

from __future__ import annotations

import dataclasses
import itertools
import math
import warnings

import numpy as np
import pytest
from scipy import stats

import longbet
import longbet._rank_k as rk
import longbet._rank_k_multi as rkm
import longbet._rank_k_multi_arm as rkma
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
    exact_moments: bool = False,
) -> rk.RankKLongBetResult:
    """Build a synthetic :class:`RankKLongBetResult` for multi-arm testing."""
    rng = np.random.default_rng(seed)
    n_periods = t_pre + t_post
    cohorts = np.full(n_units, n_periods, dtype=np.int32)
    cohorts[: n_units // 2] = t_pre
    ever = cohorts < n_periods
    unit_weights = np.where(ever, 1.0 / np.sum(ever), 0.0)

    raw_z = rng.normal(loc=0.0, scale=1.0, size=n_draws)
    if exact_moments:
        raw_z = (raw_z - float(np.mean(raw_z))) / float(np.std(raw_z, ddof=1))
    att_draws = float(att_mean) + float(att_sd) * raw_z

    mode1_draws = np.zeros((n_units, n_draws), dtype=np.float64)
    mode1_draws[ever] = att_draws[None, :] + rng.normal(0.0, 0.01, size=(int(np.sum(ever)), n_draws))
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
        rhat_att=1.001,
        ess_att=float(n_draws * 0.8),
        cohorts=cohorts,
        unit_weights=unit_weights,
    )


def test_public_exports_multi_arm() -> None:
    expected_ma = {
        "ContrastReconciliationResult",
        "MultiArmReconciliationResult",
        "SamplerDiagnostics",
        "reconcile_multi_arm_draws",
    }
    assert expected_ma <= set(rkma.__all__)
    assert expected_ma <= set(longbet.__all__)
    for name in expected_ma:
        assert getattr(longbet, name) is getattr(rkma, name)
        assert getattr(rk, name) is getattr(rkma, name)

    assert "MultiArmMultiOutcomeReconciliationResult" in longbet.__all__
    assert "reconcile_multi_arm_multi_outcome_draws" in longbet.__all__
    assert "multi_arm_decision_table" in longbet.__all__


# ---- Row 1: A = 3 delegation (bit-for-bit equality to reconcile_nested_3arm_draws) ---------------


@pytest.mark.parametrize(
    "prior_mode",
    ["bundle_halfnormal", "bundle_halfstudentt", "bundle_gamma_spike", "independent_halfnormal"],
)
def test_row1_a3_delegation_bit_for_bit_equality(prior_mode: str) -> None:
    res_bc = _make_mock_rank_k_result(att_mean=0.45, att_sd=0.14, n_draws=600, seed=101)
    res_ab = _make_mock_rank_k_result(att_mean=0.30, att_sd=0.16, n_draws=600, seed=102)
    res_ac = _make_mock_rank_k_result(att_mean=0.78, att_sd=0.15, n_draws=600, seed=103)

    legacy = rk.reconcile_nested_3arm_draws(
        res_bc,
        res_ab,
        res_ac,
        prior_scale=0.50,
        prior_mode=prior_mode,
        n_proposals=80_000,
        seed=2026,
        return_prior_draws=True,
    )
    multi = rkma.reconcile_multi_arm_draws(
        {("B", "C"): res_bc, ("A", "B"): res_ab, ("A", "C"): res_ac},
        arms=["C", "B", "A"],
        topology="chain",
        prior_scale=0.50,
        prior_mode=prior_mode,
        n_proposals=80_000,
        seed=2026,
        return_prior_draws=True,
    )

    assert np.array_equal(multi.post_att_bc, legacy.post_att_bc)
    assert np.array_equal(multi.post_att_ab, legacy.post_att_ab)
    assert np.array_equal(multi.post_att_ac, legacy.post_att_ac)
    assert np.array_equal(multi.trt_catt_bc, legacy.trt_catt_bc)
    assert np.array_equal(multi.trt_catt_ab, legacy.trt_catt_ab)
    assert np.array_equal(multi.trt_catt_ac, legacy.trt_catt_ac)
    assert np.array_equal(multi.prior_att_bc, legacy.prior_att_bc)
    assert np.array_equal(multi.prior_att_ab, legacy.prior_att_ab)
    assert np.array_equal(multi.prior_att_ac, legacy.prior_att_ac)
    assert np.array_equal(multi.prior_catt_bc, legacy.prior_catt_bc)
    assert np.array_equal(multi.prior_catt_ab, legacy.prior_catt_ab)
    assert np.array_equal(multi.prior_catt_ac, legacy.prior_catt_ac)
    assert np.array_equal(multi.contrast_mode1_draws[("B", "C")], legacy.post_att_bc)
    assert np.array_equal(multi.contrast_mode1_draws[("A", "B")], legacy.post_att_ab)
    assert np.array_equal(multi.contrast_mode1_draws[("A", "C")], legacy.post_att_ac)
    assert multi.sampler_diagnostics.ess == legacy.sir_ess
    assert multi.sampler_diagnostics["method"] == "sir"
    assert not multi.is_bipartite


# ---- Row 2: A = 3 general path (_force_general=True matches legacy within MC tolerance) ----------


def test_row2_a3_general_path_matches_legacy_within_mc_tolerance() -> None:
    res_bc = _make_mock_rank_k_result(att_mean=0.50, att_sd=0.12, n_draws=4000, seed=201)
    res_ab = _make_mock_rank_k_result(att_mean=0.35, att_sd=0.14, n_draws=4000, seed=202)
    res_ac = _make_mock_rank_k_result(att_mean=0.82, att_sd=0.13, n_draws=4000, seed=203)

    legacy = rk.reconcile_nested_3arm_draws(
        res_bc,
        res_ab,
        res_ac,
        prior_scale=0.60,
        prior_mode="bundle_halfnormal",
        n_proposals=200_000,
        seed=42,
    )
    general = rkma.reconcile_multi_arm_draws(
        {("B", "C"): res_bc, ("A", "B"): res_ab, ("A", "C"): res_ac},
        arms=["C", "B", "A"],
        topology="chain",
        prior_scale=0.60,
        prior_mode="bundle_halfnormal",
        n_proposals=200_000,
        seed=42,
        _force_general=True,
    )

    # Cut-space GLS moments and arm variances match identically
    np.testing.assert_allclose(
        general.mu_hat,
        np.array([legacy.lik_mean[0], legacy.lik_mean[0] + legacy.lik_mean[1]]),
        atol=1e-12,
    )
    for pair, leg_draws in [
        (("B", "C"), legacy.post_att_bc),
        (("A", "B"), legacy.post_att_ab),
        (("A", "C"), legacy.post_att_ac),
    ]:
        gen_draws = general.contrast_mode1_draws[pair]
        n_eff = min(float(leg_draws.size), float(legacy.sir_ess), float(general.sir_ess))
        mc_se = math.sqrt(2.0) * float(np.std(leg_draws, ddof=1)) / math.sqrt(n_eff)
        tol = max(3.0 * mc_se, 0.01)
        assert abs(float(np.mean(gen_draws)) - float(np.mean(leg_draws))) <= tol
        for q in (0.30, 0.50, 0.70):
            assert abs(float(np.quantile(gen_draws, q)) - float(np.quantile(leg_draws, q))) <= 1.5 * tol


# ---- Row 3: Flat-prior closed form (A in {3, 4, 5}) ----------------------------------------------


@pytest.mark.parametrize("n_arms", [3, 4, 5])
def test_row3_flat_prior_closed_form_matches_mu_hat_and_sigma_lik(n_arms: int) -> None:
    arms = [f"arm_{i}" for i in range(n_arms)]
    true_v = np.linspace(0.02, 0.04, n_arms)
    true_mu = np.linspace(0.3, 1.2, n_arms - 1)
    mu_full = np.concatenate([[0.0], true_mu])

    # Build complete contrast graph with exact moments
    results = {}
    seed_idx = 1000 * n_arms
    for j in range(1, n_arms):
        for l_idx in range(j):
            pair = (arms[j], arms[l_idx])
            y_mean = float(mu_full[j] - mu_full[l_idx])
            y_sd = math.sqrt(float(true_v[j] + true_v[l_idx]))
            results[pair] = _make_mock_rank_k_result(
                att_mean=y_mean,
                att_sd=y_sd,
                n_draws=2000,
                seed=seed_idx,
                exact_moments=True,
            )
            seed_idx += 1

    rec = rkma.reconcile_multi_arm_draws(
        results,
        arms=arms,
        topology="chain",
        prior_mode="flat",
        monotone=False,
        tau_clip=None,
        seed=777 + n_arms,
    )

    np.testing.assert_allclose(rec.mu_hat, true_mu, atol=1e-10)
    sample_mean = np.mean(rec.tau_mu_draws, axis=0)
    sample_cov = np.cov(rec.tau_mu_draws, rowvar=False, ddof=1)
    n_draws = rec.tau_mu_draws.shape[0]
    se_mean = np.sqrt(np.diag(rec.sigma_lik) / float(n_draws))
    assert np.all(np.abs(sample_mean - rec.mu_hat) <= 3.0 * se_mean)
    diag_se = np.diag(rec.sigma_lik) * math.sqrt(2.0 / float(n_draws - 1))
    assert np.all(np.abs(np.diag(sample_cov) - np.diag(rec.sigma_lik)) <= 3.5 * diag_se)


# ---- Row 4: SBC rank uniformity across A in {3, 4, 5} --------------------------------------------


@pytest.mark.parametrize(
    ("n_arms", "topology", "prior_mode", "sampler"),
    [
        (3, "chain", "bundle_halfnormal", "sir"),
        (4, "chain", "bundle_halfstudentt", "sir"),
        (4, "poset", "bundle_halfnormal", "sir"),
        (5, "chain", "bundle_gamma_spike", "slice"),
    ],
)
def test_row4_sbc_rank_uniformity(
    n_arms: int, topology: str, prior_mode: str, sampler: str
) -> None:
    rng = np.random.default_rng(4000 + 100 * n_arms + len(prior_mode))
    dim = n_arms - 1
    arms = [f"arm_{i}" for i in range(n_arms)]
    if topology == "poset" and n_arms == 4:
        # 2x2 factorial diamond: 0 -> 1, 0 -> 2, 1 -> 3, 2 -> 3 (spanning tree edges: (1,0), (2,0), (3,1))
        poset_edges_idx = [(1, 0), (2, 0), (3, 1), (3, 2)]
        non_top_arms = [1, 2]
        lin_exts = rkma._enumerate_linear_extensions(non_top_arms, poset_edges_idx)
    else:
        poset_edges_idx = [(a, a - 1) for a in range(1, n_arms)]
        non_top_arms = list(range(1, n_arms - 1))
        lin_exts = None

    scale_val = 0.60
    alpha_vec = np.ones(dim, dtype=np.float64)
    tau_clip = (0.0, float("inf"))
    sigma_lik = 0.04 * np.eye(dim, dtype=np.float64) + 0.015 * np.ones((dim, dim), dtype=np.float64)
    chol_lik = np.linalg.cholesky(sigma_lik)

    def _p_draw(n: int) -> np.ndarray:
        return rkma._draw_prior_mu(
            rng,
            n=n,
            dim=dim,
            topology=topology,
            prior_mode=prior_mode,
            scale_val=scale_val,
            alpha_vec=alpha_vec,
            monotone=True,
            tau_clip=tau_clip,
            poset_edges_idx=poset_edges_idx,
            non_top_arms=non_top_arms,
            linear_extensions=lin_exts,
        )

    def _p_logpdf(pts: np.ndarray) -> np.ndarray:
        return rkma._prior_logpdf_mu(
            pts,
            dim=dim,
            topology=topology,
            prior_mode=prior_mode,
            scale_val=scale_val,
            alpha_vec=alpha_vec,
            monotone=True,
            tau_clip=tau_clip,
            poset_edges_idx=poset_edges_idx,
            non_top_arms=non_top_arms,
            linear_extensions=lin_exts,
        )

    n_sbc = 48
    n_bins = 4
    n_post = 120
    ranks = np.zeros((n_sbc, dim), dtype=np.float64)
    for rep in range(n_sbc):
        mu_star = _p_draw(1)[0]
        mu_hat = mu_star + chol_lik @ rng.standard_normal(dim)
        if sampler == "sir":
            post_mu, _, _ = rkma._sir_posterior_nd(
                rng,
                prior_draw=_p_draw,
                prior_logpdf=_p_logpdf,
                mu_hat=mu_hat,
                sigma_lik=sigma_lik,
                n_draws=n_post,
                n_prop=15_000,
            )
        else:
            post_mu, _, _ = rkma._hit_and_run_slice_sampler(
                rng,
                prior_draw=_p_draw,
                prior_logpdf=_p_logpdf,
                mu_hat=mu_hat,
                sigma_lik=sigma_lik,
                n_draws=n_post,
                num_chains=2,
                num_warmup=80,
                thin=1,
            )
        ranks[rep] = np.mean(post_mu < mu_star[None, :], axis=0)

    for j in range(dim):
        counts, _ = np.histogram(ranks[:, j], bins=np.linspace(0.0, 1.0, n_bins + 1))
        stat = stats.chisquare(counts)
        assert float(stat.pvalue) > 0.01, f"SBC rank uniformity failed for coord {j}: p={stat.pvalue:.4f}"
    _ = arms


# ---- Row 5: Variance identification (3-cycle, 4-cycle, tree, all-pairs A=5) ----------------------


def test_row5_variance_identification_across_topologies() -> None:
    # 1. Complete 3-cycle (triangle): exact recovery when no floor binds
    arms3 = ["C", "B", "A"]
    v_true3 = np.array([0.02, 0.03, 0.04], dtype=np.float64)
    pairs3 = {
        ("B", "C"): _make_mock_rank_k_result(
            att_mean=0.3, att_sd=math.sqrt(v_true3[1] + v_true3[0]), seed=501, exact_moments=True
        ),
        ("A", "B"): _make_mock_rank_k_result(
            att_mean=0.2, att_sd=math.sqrt(v_true3[2] + v_true3[1]), seed=502, exact_moments=True
        ),
        ("A", "C"): _make_mock_rank_k_result(
            att_mean=0.5, att_sd=math.sqrt(v_true3[2] + v_true3[0]), seed=503, exact_moments=True
        ),
    }
    rec3 = rkma.reconcile_multi_arm_draws(pairs3, arms=arms3, _force_general=True, seed=510)
    assert not rec3.is_bipartite
    np.testing.assert_allclose(np.asarray(rec3.arm_variances), v_true3, rtol=1e-6)
    assert rec3.arm_variances["C"] == pytest.approx(v_true3[0], rel=1e-6)
    assert rec3.arm_variances["B"] == pytest.approx(v_true3[1], rel=1e-6)
    assert rec3.arm_variances["A"] == pytest.approx(v_true3[2], rel=1e-6)
    assert rec3.nnls_residual < 1e-10

    # 2. 4-cycle (bipartite: 0-1-2-3-0): emits UserWarning and uses ridge path
    arms4 = ["a0", "a1", "a2", "a3"]
    pairs4_cycle = {
        ("a1", "a0"): _make_mock_rank_k_result(att_mean=0.2, att_sd=0.20, seed=511, exact_moments=True),
        ("a2", "a1"): _make_mock_rank_k_result(att_mean=0.2, att_sd=0.22, seed=512, exact_moments=True),
        ("a3", "a2"): _make_mock_rank_k_result(att_mean=0.2, att_sd=0.21, seed=513, exact_moments=True),
        ("a3", "a0"): _make_mock_rank_k_result(att_mean=0.6, att_sd=0.19, seed=514, exact_moments=True),
    }
    with pytest.warns(UserWarning, match="bipartite"):
        rec4_cycle = rkma.reconcile_multi_arm_draws(
            pairs4_cycle,
            arms=arms4,
            sample_sizes={"a0": 200, "a1": 150, "a2": 150, "a3": 200},
            seed=520,
        )
    assert rec4_cycle.is_bipartite
    assert np.all(np.asarray(rec4_cycle.arm_variances) > 0.0)

    # 3. Tree (star/chain of A=4 with 3 edges, bipartite): emits UserWarning
    pairs4_tree = {
        k: v for k, v in pairs4_cycle.items() if k != ("a3", "a0")
    }
    with pytest.warns(UserWarning, match="bipartite"):
        rec4_tree = rkma.reconcile_multi_arm_draws(pairs4_tree, arms=arms4, seed=530)
    assert rec4_tree.is_bipartite
    assert len(rec4_tree.redundant_pairs) == 0

    # 4. A=5 all-pairs (10 contrasts, overdetermined non-bipartite): NNLS residual reported
    arms5 = [f"arm_{i}" for i in range(5)]
    v_true5 = np.array([0.015, 0.020, 0.025, 0.030, 0.035], dtype=np.float64)
    pairs5 = {}
    for j in range(1, 5):
        for l_idx in range(j):
            s_jl = math.sqrt(float(v_true5[j] + v_true5[l_idx]))
            pairs5[(arms5[j], arms5[l_idx])] = _make_mock_rank_k_result(
                att_mean=0.15 * (j - l_idx),
                att_sd=s_jl,
                n_draws=400,
                seed=540 + 10 * j + l_idx,
                exact_moments=True,
            )
    rec5 = rkma.reconcile_multi_arm_draws(
        pairs5, arms=arms5, num_warmup=150, seed=599
    )
    assert not rec5.is_bipartite
    assert np.isfinite(rec5.nnls_residual) and rec5.nnls_residual < 1e-6
    np.testing.assert_allclose(np.asarray(rec5.arm_variances), v_true5, rtol=0.05)


# ---- Row 6: Cycle coherence across all triples, draws, and periods -------------------------------


def test_row6_cycle_coherence_complete_a4_and_poset_a5() -> None:
    # A=4 complete graph (6 pairs)
    arms4 = ["c", "d1", "d2", "d3"]
    pairs4 = {}
    for j in range(1, 4):
        for l_idx in range(j):
            pairs4[(arms4[j], arms4[l_idx])] = _make_mock_rank_k_result(
                att_mean=0.25 * (j - l_idx) + 0.02 * j,
                att_sd=0.15,
                n_draws=500,
                seed=600 + 10 * j + l_idx,
            )
    rec4 = rkma.reconcile_multi_arm_draws(pairs4, arms=arms4, topology="chain", seed=610)
    assert rec4.trajectory_additivity_max_gap < 1e-12
    for a_i, a_j, a_k in itertools.permutations(arms4, 3):
        t_ik = rec4.get_contrast_result(a_i, a_k).trt_catt
        t_ij = rec4.get_contrast_result(a_i, a_j).trt_catt
        t_jk = rec4.get_contrast_result(a_j, a_k).trt_catt
        assert float(np.max(np.abs(t_ik - (t_ij + t_jk)))) < 1e-12
        m_ik = rec4.get_contrast_result(a_i, a_k).post_att
        m_ij = rec4.get_contrast_result(a_i, a_j).post_att
        m_jk = rec4.get_contrast_result(a_j, a_k).post_att
        assert float(np.max(np.abs(m_ik - (m_ij + m_jk)))) < 1e-12

    # A=5 poset (control + 2x2 factorial with top arm)
    arms5 = ["ctrl", "f1", "f2", "f1_hi", "f12_top"]
    edges5 = [("f1", "ctrl"), ("f2", "ctrl"), ("f1_hi", "f1"), ("f12_top", "f1_hi"), ("f12_top", "f2")]
    pairs5 = {
        edge: _make_mock_rank_k_result(att_mean=0.30, att_sd=0.14, n_draws=400, seed=650 + idx)
        for idx, edge in enumerate(edges5)
    }
    pairs5[("f12_top", "ctrl")] = _make_mock_rank_k_result(
        att_mean=0.90, att_sd=0.16, n_draws=400, seed=660
    )
    rec5 = rkma.reconcile_multi_arm_draws(
        pairs5,
        arms=arms5,
        topology="poset",
        edges=edges5,
        num_warmup=150,
        seed=670,
    )
    assert rec5.trajectory_additivity_max_gap < 1e-12
    assert set(rec5.own_vs_implied_max_gap.keys()) == set(rec5.redundant_pairs)
    for a_i, a_j, a_k in itertools.permutations(arms5, 3):
        t_ik = rec5.get_contrast_result(a_i, a_k).trt_catt
        t_ij = rec5.get_contrast_result(a_i, a_j).trt_catt
        t_jk = rec5.get_contrast_result(a_j, a_k).trt_catt
        assert float(np.max(np.abs(t_ik - (t_ij + t_jk)))) < 1e-12


# ---- Row 7: Order-polytope prior (chain spacings vs Dirichlet(1); 2x2 independent Uniform) -------


def test_row7_order_polytope_prior_chain_and_2x2_factorial() -> None:
    rng = np.random.default_rng(707)
    n_samples = 8000

    # 1. Chain poset on A=4 arms (0 < 1 < 2 < 3): w_1 <= w_2 in [0, 1]^2 -> uniform order statistics
    poset_chain_edges = [(1, 0), (2, 1), (3, 2)]
    non_top_chain = [1, 2]
    lin_chain = rkma._enumerate_linear_extensions(non_top_chain, poset_chain_edges)
    assert lin_chain.shape == (1, 2)
    mu_poset_chain = rkma._draw_prior_mu(
        rng,
        n=n_samples,
        dim=3,
        topology="poset",
        prior_mode="bundle_halfnormal",
        scale_val=1.0,
        alpha_vec=np.ones(3),
        monotone=True,
        tau_clip=(0.0, float("inf")),
        poset_edges_idx=poset_chain_edges,
        non_top_arms=non_top_chain,
        linear_extensions=lin_chain,
    )
    # Normalized increments (theta_1, theta_2, theta_3) / tau_top must follow Dirichlet(1, 1, 1),
    # whose marginal is Beta(1, 2)
    tau_top = mu_poset_chain[:, 2]
    u1 = mu_poset_chain[:, 0] / tau_top
    u2 = (mu_poset_chain[:, 1] - mu_poset_chain[:, 0]) / tau_top
    u3 = (mu_poset_chain[:, 2] - mu_poset_chain[:, 1]) / tau_top
    for u_vec in (u1, u2, u3):
        ks_res = stats.kstest(u_vec, stats.beta(1.0, 2.0).cdf)
        assert float(ks_res.pvalue) > 0.05

    # 2. 2x2 factorial diamond on A=4 arms (0 < 1 < 3, 0 < 2 < 3; 1 and 2 incomparable):
    # w_1 and w_2 are independent Uniform(0, 1)
    poset_2x2_edges = [(1, 0), (2, 0), (3, 1), (3, 2)]
    non_top_2x2 = [1, 2]
    lin_2x2 = rkma._enumerate_linear_extensions(non_top_2x2, poset_2x2_edges)
    assert lin_2x2.shape == (2, 2)
    mu_2x2 = rkma._draw_prior_mu(
        rng,
        n=n_samples,
        dim=3,
        topology="poset",
        prior_mode="bundle_halfnormal",
        scale_val=1.0,
        alpha_vec=np.ones(3),
        monotone=True,
        tau_clip=(0.0, float("inf")),
        poset_edges_idx=poset_2x2_edges,
        non_top_arms=non_top_2x2,
        linear_extensions=lin_2x2,
    )
    w1 = mu_2x2[:, 0] / mu_2x2[:, 2]
    w2 = mu_2x2[:, 1] / mu_2x2[:, 2]
    assert float(stats.kstest(w1, "uniform").pvalue) > 0.05
    assert float(stats.kstest(w2, "uniform").pvalue) > 0.05
    assert abs(float(np.corrcoef(w1, w2)[0, 1])) < 0.04


# ---- Row 8: Multi-outcome A=4, M=3 (bottom_up & coherent_total + decision tables) ----------------


def test_row8_multi_outcome_a4_m3_bottom_up_and_coherent_total() -> None:
    arms = ["ctrl", "tier_1", "tier_2", "tier_3"]
    outcomes = ("out_a", "out_b", "out_c")
    results_by_outcome: dict[str, dict[tuple[str, str], rk.RankKLongBetResult]] = {
        m: {} for m in outcomes
    }
    base_means = {"out_a": 0.08, "out_b": 0.06, "out_c": 0.20}
    base_sds = {"out_a": 0.12, "out_b": 0.11, "out_c": 0.10}

    seed_counter = 800
    for m in outcomes:
        for j in range(1, 4):
            for l_idx in range(j):
                pair = (arms[j], arms[l_idx])
                results_by_outcome[m][pair] = _make_mock_rank_k_result(
                    att_mean=base_means[m] * (j - l_idx),
                    att_sd=base_sds[m] * math.sqrt(1.0 + 0.15 * (j - l_idx)),
                    n_draws=600,
                    seed=seed_counter,
                )
                seed_counter += 1

    target_corr = np.array(
        [
            [1.0, 0.55, 0.30],
            [0.55, 1.0, 0.40],
            [0.30, 0.40, 1.0],
        ],
        dtype=np.float64,
    )

    # 1. bottom_up mode
    bu = rkm.reconcile_multi_arm_multi_outcome_draws(
        results_by_outcome,
        arms=arms,
        topology="chain",
        prior_scales=0.50,
        cross_outcome_corr=target_corr,
        composite_outcomes={"sub_ab": ("out_a", "out_b")},
        composite_mode="bottom_up",
        n_proposals=60_000,
        seed=888,
        return_prior_draws=True,
    )
    assert bu.outcomes == ("out_a", "out_b", "out_c", "sub_ab")
    assert bu.additivity_violation_rate["sub_ab"] == 0.0
    assert bu.max_abs_composite_discrepancy["sub_ab"] < 1e-12
    np.testing.assert_allclose(bu.realized_spearman, target_corr, atol=0.05)

    # Verify exact marginal preservation under bottom_up compared to standalone per-outcome reconciliation
    standalone_a = rkma.reconcile_multi_arm_draws(
        results_by_outcome["out_a"],
        arms=arms,
        topology="chain",
        prior_scale=0.50,
        n_proposals=60_000,
        seed=888,
    )
    for pair in standalone_a.fitted_pairs:
        np.testing.assert_allclose(
            np.sort(bu["out_a"].contrast_mode1_draws[pair]),
            np.sort(standalone_a.contrast_mode1_draws[pair]),
            atol=1e-12,
        )

    # 2. coherent_total mode: AVR=0 across all 6 contrasts AND null-floor ratio <= 1.05
    ct = rkm.reconcile_multi_arm_multi_outcome_draws(
        results_by_outcome,
        arms=arms,
        topology="chain",
        prior_scales=0.50,
        cross_outcome_corr=target_corr,
        composite_outcomes={"sub_ab": ("out_a", "out_b")},
        composite_mode="coherent_total",
        n_proposals=60_000,
        seed=888,
    )
    assert ct.additivity_violation_rate["sub_ab"] == 0.0
    assert ct.max_abs_composite_discrepancy["sub_ab"] < 1e-12

    # Check all 6 unordered contrasts have exact draw-wise composite additivity
    for j in range(1, 4):
        for l_idx in range(j):
            pair = (arms[j], arms[l_idx])
            comp_t = ct["sub_ab"].contrast_trajectory_draws[pair]
            sum_t = (
                ct["out_a"].contrast_trajectory_draws[pair]
                + ct["out_b"].contrast_trajectory_draws[pair]
            )
            assert float(np.max(np.abs(comp_t - sum_t))) < 1e-12

    # Single-truncation composite benchmark on key_contrast = ("tier_3", "ctrl")
    single_trunc_mu, _, _, _, _, _ = rkm._sample_synthetic_composite_multi_arm(
        ("out_a", "out_b"),
        ("out_a", "out_b", "out_c"),
        tuple(arms),
        "chain",
        ct.tree_edges,
        {m: bu[m] for m in ("out_a", "out_b", "out_c")},
        {a: target_corr for a in arms},
        target_corr,
        comp_prior_scale=None,
        prior_mode="bundle_halfnormal",
        split_alpha=None,
        monotone=True,
        tau_clip=(0.0, float("inf")),
        n_draws=600,
        n_proposals=60_000,
        num_chains=4,
        num_warmup=200,
        thin=2,
        seed=888 + 20011,
    )
    key_pair = ("tier_3", "ctrl")
    ct_key_mean = float(np.mean(ct["sub_ab"].contrast_mode1_draws[key_pair]))
    st_key_mean = float(np.mean(single_trunc_mu[:, -1]))
    assert ct_key_mean / max(st_key_mean, 1e-12) <= 1.05

    # 3. Also test decision tables and GE revenue draws on multi-arm results
    dec_df = rkp.multi_arm_decision_table(
        ct,
        step_costs={"tier_1": 0.5, "tier_2": 0.5, "tier_3": 0.5},
        alpha={"sub_ab": 1.0, "out_c": 0.8},
    )
    assert list(dec_df["arm"]) == arms
    assert float(dec_df["pr_best_arm"].sum()) == pytest.approx(1.0, abs=1e-12)
    assert float(dec_df["pr_best_treated_arm"].sum()) == pytest.approx(1.0, abs=1e-12)
    assert np.all(np.isfinite(dec_df.loc[~dec_df["is_control"], "pr_incremental_step_worthwhile"]))

    roi_df = rkp.roi_break_even_table(
        ct["sub_ab"],
        contrast=[("tier_1", "ctrl"), ("tier_3", "ctrl")],
        cost_grid=(1.0, 2.0),
        alpha_grid=(0.5, 1.0),
    )
    assert set(roi_df["contrast"]) == {("tier_1", "ctrl"), ("tier_3", "ctrl")}

    ge_draws = rkm.ge_revenue_draws(
        ct,
        weights={"sub_ab": 1.0, "out_c": 0.5},
        contrast=("tier_3", "ctrl"),
    )
    assert ge_draws.shape == (600,) and np.all(np.isfinite(ge_draws))


# ---- Row 9: Input validation (disconnected graph, missing tree edge, draw mismatch, non-unique top)


def test_row9_input_validation_raises_clear_value_errors() -> None:
    fit_600 = _make_mock_rank_k_result(att_mean=0.3, att_sd=0.1, n_draws=600, seed=901)
    fit_400 = _make_mock_rank_k_result(att_mean=0.4, att_sd=0.1, n_draws=400, seed=902)

    # 1. Disconnected contrast graph
    with pytest.raises(ValueError, match="connected"):
        rkma.reconcile_multi_arm_draws(
            {("a1", "a0"): fit_600, ("a3", "a2"): fit_600},
            arms=["a0", "a1", "a2", "a3"],
        )

    # 2. Missing spanning-tree edge fit
    with pytest.raises(ValueError, match="spanning-tree edge"):
        rkma.reconcile_multi_arm_draws(
            {("a1", "a0"): fit_600, ("a2", "a0"): fit_600, ("a3", "a0"): fit_600},
            arms=["a0", "a1", "a2", "a3"],
            topology="chain",
        )

    # 3. Inconsistent draw counts across fitted pairs
    with pytest.raises(ValueError, match="draw count"):
        rkma.reconcile_multi_arm_draws(
            {("a1", "a0"): fit_600, ("a2", "a1"): fit_400},
            arms=["a0", "a1", "a2"],
            topology="chain",
        )

    # 4. Non-unique top arm in poset with bundle prior (two maximal arms a1 and a2)
    with pytest.raises(ValueError, match="unique top arm"):
        rkma.reconcile_multi_arm_draws(
            {("a1", "a0"): fit_600, ("a2", "a0"): fit_600},
            arms=["a0", "a1", "a2"],
            topology="poset",
            edges=[("a1", "a0"), ("a2", "a0")],
            prior_mode="bundle_halfnormal",
        )

    # Note: non-unique top poset IS supported when prior_mode="independent"
    with pytest.warns(UserWarning, match="bipartite"):
        rec_star = rkma.reconcile_multi_arm_draws(
            {("a1", "a0"): fit_600, ("a2", "a0"): fit_600},
            arms=["a0", "a1", "a2"],
            topology="poset",
            edges=[("a1", "a0"), ("a2", "a0")],
            prior_mode="independent",
            seed=905,
        )
    assert rec_star.tau_mu_draws.shape == (600, 2)


def test_reverse_oriented_tree_edge_in_reconcile_multi_arm_draws() -> None:
    arms = ["C", "B1", "B2", "A"]
    fit_b1_c = _make_mock_rank_k_result(att_mean=0.20, att_sd=0.10, n_draws=300, seed=1001)
    fit_b2_b1 = _make_mock_rank_k_result(att_mean=0.25, att_sd=0.11, n_draws=300, seed=1002)
    fit_a_b2 = _make_mock_rank_k_result(att_mean=0.30, att_sd=0.12, n_draws=300, seed=1003)
    fit_a_c = _make_mock_rank_k_result(att_mean=0.75, att_sd=0.13, n_draws=300, seed=1004)

    # Create reverse-oriented fit ("B1", "B2") with exact negated draws
    fit_b1_b2 = dataclasses.replace(
        fit_b2_b1,
        att_draws=-fit_b2_b1.att_draws,
        att_hat=-fit_b2_b1.att_hat,
        att_ci_low=-fit_b2_b1.att_ci_high,
        att_ci_high=-fit_b2_b1.att_ci_low,
        catt_draws=-fit_b2_b1.catt_draws,
        catt_hat=-fit_b2_b1.catt_hat,
        mode_tau_draws=tuple(-m for m in fit_b2_b1.mode_tau_draws),
    )

    rec_fwd = rkma.reconcile_multi_arm_draws(
        {
            ("B1", "C"): fit_b1_c,
            ("B2", "B1"): fit_b2_b1,
            ("A", "B2"): fit_a_b2,
            ("A", "C"): fit_a_c,
        },
        arms=arms,
        topology="chain",
        prior_scale=0.50,
        return_prior_draws=True,
        seed=1010,
    )
    rec_rev = rkma.reconcile_multi_arm_draws(
        {
            ("B1", "C"): fit_b1_c,
            ("B1", "B2"): fit_b1_b2,
            ("A", "B2"): fit_a_b2,
            ("A", "C"): fit_a_c,
        },
        arms=arms,
        topology="chain",
        prior_scale=0.50,
        return_prior_draws=True,
        seed=1010,
    )

    np.testing.assert_allclose(rec_rev.mu_hat, rec_fwd.mu_hat, atol=1e-12)
    np.testing.assert_allclose(rec_rev.tau_increment_draws, rec_fwd.tau_increment_draws, atol=1e-12)
    np.testing.assert_allclose(
        rec_rev.get_contrast_result("B2", "B1").trajectory_draws,
        rec_fwd.get_contrast_result("B2", "B1").trajectory_draws,
        atol=1e-12,
    )
    assert rec_rev.prior_contrast_trajectory_draws is not None
    assert rec_fwd.prior_contrast_trajectory_draws is not None
    np.testing.assert_allclose(
        rec_rev.prior_contrast_trajectory_draws[("B2", "B1")],
        rec_fwd.prior_contrast_trajectory_draws[("B2", "B1")],
        atol=1e-12,
    )


def test_bottom_up_composite_with_flat_and_independent_priors_and_reversed_key_contrast() -> None:
    arms = ["ctrl", "tier_1", "tier_2"]
    outcomes = ("out_a", "out_b")
    results_by_outcome: dict[str, dict[tuple[str, str], rk.RankKLongBetResult]] = {
        m: {} for m in outcomes
    }
    seed_idx = 1100
    for m in outcomes:
        for j in range(1, 3):
            for l_idx in range(j):
                results_by_outcome[m][(arms[j], arms[l_idx])] = _make_mock_rank_k_result(
                    att_mean=0.15 * (j - l_idx),
                    att_sd=0.10,
                    n_draws=200,
                    seed=seed_idx,
                )
                seed_idx += 1

    # 1. prior_mode="flat" (prior_scale is None) + reverse-oriented key_contrast=("ctrl", "tier_2")
    bu_flat = rkm.reconcile_multi_arm_multi_outcome_draws(
        results_by_outcome,
        arms=arms,
        topology="chain",
        prior_mode="flat",
        monotone=False,
        tau_clip=None,
        key_contrast=("ctrl", "tier_2"),
        composite_outcomes={"tot": ("out_a", "out_b")},
        composite_mode="bottom_up",
        seed=1111,
    )
    assert bu_flat["tot"].prior_scale is None
    assert bu_flat.additivity_violation_rate["tot"] == 0.0
    assert bu_flat.max_abs_composite_discrepancy["tot"] < 1e-12

    # 2. prior_mode="independent" (prior_scale is a tuple of floats) + bottom_up
    bu_indep = rkm.reconcile_multi_arm_multi_outcome_draws(
        results_by_outcome,
        arms=arms,
        topology="chain",
        prior_mode="independent",
        prior_scales={"out_a": (0.3, 0.4), "out_b": (0.2, 0.3)},
        composite_outcomes={"tot": ("out_a", "out_b")},
        composite_mode="bottom_up",
        seed=1112,
    )
    assert isinstance(bu_indep["tot"].prior_scale, tuple)
    assert len(bu_indep["tot"].prior_scale) == 2
    assert bu_indep.additivity_violation_rate["tot"] == 0.0
    assert bu_indep.max_abs_composite_discrepancy["tot"] < 1e-12


def test_multi_arm_decision_table_alias_kwargs_and_columns() -> None:
    arms = ["C", "B1", "B2", "A"]
    pairs = {}
    for j in range(1, 4):
        for l_idx in range(j):
            pairs[(arms[j], arms[l_idx])] = _make_mock_rank_k_result(
                att_mean=0.25 * (j - l_idx),
                att_sd=0.10,
                n_draws=300,
                seed=1200 + 10 * j + l_idx,
            )
    rec = rkma.reconcile_multi_arm_draws(pairs, arms=arms, prior_scale=0.50, seed=1250)
    dec_df = rkp.multi_arm_decision_table(
        rec,
        arm_costs={"C": 0.0, "B1": 1.0, "B2": 2.0, "A": 3.5},
        incremental_costs={("B1", "C"): 1.0, ("B2", "B1"): 1.0, ("A", "B2"): 1.5},
        outcome_weights=1.0,
    )
    for col in (
        "pr_best_arm",
        "pr_best_treated_arm",
        "pr_positive_net_value",
        "pr_incremental_step_worthwhile",
        "expected_net_value",
        "p_best_gross",
        "p_best_net",
        "p_beats_control",
        "p_positive_step",
        "p_step_roi_positive",
    ):
        assert col in dec_df.columns
    np.testing.assert_allclose(dec_df["p_best_net"], dec_df["pr_best_arm"], atol=1e-12)
    np.testing.assert_allclose(
        dec_df.loc[~dec_df["is_control"], "p_step_roi_positive"],
        dec_df.loc[~dec_df["is_control"], "pr_incremental_step_worthwhile"],
        atol=1e-12,
    )


def test_independent_prior_logpdf_exact_truncnorm_normalization_and_monotone() -> None:
    pts = np.array([[0.2, 0.5], [0.6, 0.3], [-0.1, 0.4]], dtype=np.float64)
    sd = 0.8
    tau_clip = (-0.5, 1.5)
    logp = rkma._prior_logpdf_mu(
        pts,
        dim=2,
        topology="chain",
        prior_mode="independent",
        scale_val=sd,
        alpha_vec=np.ones(2),
        monotone=True,
        tau_clip=tau_clip,
        poset_edges_idx=[(1, 0), (2, 1)],
        non_top_arms=[1],
    )
    # Row 0 (0.2 <= 0.5 and >= 0) is valid and matches exact scipy.stats.truncnorm.logpdf sum
    a_std, b_std = tau_clip[0] / sd, tau_clip[1] / sd
    expected_row0 = float(
        np.sum(stats.truncnorm.logpdf(pts[0], a_std, b_std, loc=0.0, scale=sd))
    )
    assert logp[0] == pytest.approx(expected_row0, abs=1e-12)
    # Row 1 violates chain monotonicity (0.6 > 0.3) and Row 2 violates mu_1 >= 0 (-0.1 < 0)
    assert logp[1] == -np.inf
    assert logp[2] == -np.inf

