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

"""Stage 1 API hygiene, mathematical property, and sampler exactness tests for v0.11 (`methodology_v3_final_plan.md`)."""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest
from scipy import integrate, stats

import longbet
import longbet._rank_k as rk_mod
import longbet._rank_k_multi as rkm_mod
import longbet._rank_k_multi_arm as rkma_mod
from longbet import (
    RankKLongBetConfig,
    orthonormal_exposure_basis,
    reconcile_multi_arm_draws,
    reconcile_multi_arm_multi_outcome_draws,
    reconcile_nested_3arm_draws,
    reconcile_nested_3arm_multi_draws,
)
from longbet._rank_k_multi_arm import (
    _OutcomeConePriorSpec,
    _asis_interweave_outcome_scales,
    _resolve_arms_and_poset,
    _sample_invgamma,
    _sample_truncated_mvn,
    _update_scales_outcome,
)
from tests.test_rank_k_multi_arm import _make_stage1_fit


PRUNED_SYMBOLS = (
    "iman_conover_couple",
    "ImanConoverResult",
    "_empirical_bayes_scale",
    "_log_marginal_1d",
    "_huber_cut_space_projection",
    "_spike_slab_mix",
    "_enumerate_linear_extensions",
    "_sir_posterior_nd",
    "_find_interior_starting_points",
    "_slice_posterior_nd",
    "_project_cohort_damped_trajectory",
    "_rescale_multi_arm_outcome",
    "_sample_synthetic_composite_3arm",
    "_sample_synthetic_composite_multi_arm",
    "_ledoit_wolf_correlation_shrinkage",
    "_recover_arm_variances",
    "_bipartite_null_interval",
)

PRUNED_KWARGS = (
    "prior_mode",
    "spike_prob",
    "projection",
    "huber_delta",
    "sampler",
    "n_proposals",
    "n_draws",
    "n_burnin",
    "num_burnin",
    "n_chains",
    "tau_clip",
    "catt_floor",
    "composite_mode",
    "key_contrast",
    "arm_variance_source",
    "bipartite_policy",
    "bipartite_ridge",
    "bipartite_ridge_fraction",
    "min_group_var_share",
    "sample_sizes",
    "return_prior_draws",
    "share_mean_logit",
    "share_sd_logit",
    "prior_sd_bc",
    "prior_sd_ab",
    "prior_sd_ac",
)


@pytest.mark.parametrize("sym", PRUNED_SYMBOLS)
def test_pruned_symbols_do_not_exist_anywhere(sym: str) -> None:
    for mod in (longbet, rk_mod, rkma_mod, rkm_mod):
        assert not hasattr(mod, sym), f"{sym} must be deleted from {mod.__name__}"


@pytest.mark.parametrize("kw", PRUNED_KWARGS)
def test_pruned_kwargs_raise_type_error_on_all_reconciliation_entrypoints(kw: str) -> None:
    for fn in (
        reconcile_multi_arm_draws,
        reconcile_multi_arm_multi_outcome_draws,
        reconcile_nested_3arm_draws,
        reconcile_nested_3arm_multi_draws,
    ):
        sig = inspect.signature(fn)
        assert kw not in sig.parameters, f"Pruned kwarg {kw!r} still in {fn.__name__} signature"


def test_public_reconciliation_signatures_and_thin_3arm_wrappers() -> None:
    sig_ma = tuple(inspect.signature(reconcile_multi_arm_draws).parameters)
    assert sig_ma == (
        "results",
        "arms",
        "order",
        "poset",
        "monotone",
        "prior_scale",
        "horizon",
        "n_modes",
        "seed",
        "_model_options",
    )
    sig_mamo = tuple(inspect.signature(reconcile_multi_arm_multi_outcome_draws).parameters)
    assert sig_mamo == (
        "results_by_outcome",
        "composite_outcomes",
        "composite_results",
        "cross_outcome_residuals",
        "arms",
        "order",
        "poset",
        "monotone",
        "prior_scale",
        "horizon",
        "n_modes",
        "seed",
        "_model_options",
    )

    for wrapper in (reconcile_nested_3arm_draws, reconcile_nested_3arm_multi_draws):
        lines = [
            ln
            for ln in inspect.getsource(wrapper).splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        assert len(lines) <= 45, f"{wrapper.__name__} must remain a thin convenience wrapper"


@pytest.mark.parametrize("bad_scale", ["empirical_bayes", "in_sample", "auto", 0.0, -0.5, True, False, float("nan")])
def test_string_or_nonpositive_prior_scale_raises_value_error(bad_scale) -> None:
    rng = np.random.default_rng(1)
    fit = _make_stage1_fit(rng, np.ones(4) * 0.5)
    results = {("B", "C"): fit, ("A", "B"): fit, ("A", "C"): fit}
    with pytest.raises(ValueError, match="prior_scale"):
        reconcile_multi_arm_draws(results, arms=("C", "B", "A"), order=("C", "B", "A"), prior_scale=bad_scale)
    with pytest.raises(ValueError, match="prior_scale"):
        reconcile_multi_arm_multi_outcome_draws(
            {"m1": results}, arms=("C", "B", "A"), order=("C", "B", "A"), prior_scale={"m1": bad_scale}
        )


def test_stage1_treated_extra_variance_default_policy() -> None:
    assert RankKLongBetConfig().treated_extra_variance is True
    assert RankKLongBetConfig(baseline_cells="untreated").treated_extra_variance is True
    assert RankKLongBetConfig(baseline_cells="all").treated_extra_variance is False
    assert RankKLongBetConfig(baseline_cells="untreated", treated_extra_variance=False).treated_extra_variance is False
    with pytest.raises(ValueError, match="baseline_cells='untreated'"):
        RankKLongBetConfig(baseline_cells="all", treated_extra_variance=True)


def test_pakman_paninski_tn_hmc_matches_exact_rejection_sampling_in_3d() -> None:
    rng = np.random.default_rng(202605)
    mean = np.array([0.25, 0.55, -0.10], dtype=np.float64)
    cov = np.array(
        [
            [0.30, 0.12, -0.05],
            [0.12, 0.40, 0.08],
            [-0.05, 0.08, 0.25],
        ],
        dtype=np.float64,
    )
    prec = np.linalg.inv(cov)
    a_cone = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 1.0, 0.0],
            [1.0, 0.0, 0.8],
            [1.0, 0.0, -0.8],
        ],
        dtype=np.float64,
    )
    interior_ray = np.array([1.0, 2.0, 0.0], dtype=np.float64)

    # Exact rejection sample
    raw = rng.multivariate_normal(mean, cov, size=250_000)
    keep = np.all((raw @ a_cone.T) >= 0.0, axis=1)
    exact_draws = raw[keep][:8000]
    assert exact_draws.shape[0] == 8000

    # Pakman-Paninski TN-HMC chain
    th = interior_ray.copy()
    hmc_draws = np.empty((4000, 3), dtype=np.float64)
    for s in range(4400):
        th = _sample_truncated_mvn(mean, prec, a_cone, th, interior_ray, rng, n_hmc_steps=2)
        if s >= 400:
            hmc_draws[s - 400] = th

    assert np.all((hmc_draws @ a_cone.T) >= -1e-13)
    np.testing.assert_allclose(hmc_draws.mean(axis=0), exact_draws.mean(axis=0), atol=0.025)
    np.testing.assert_allclose(hmc_draws.std(axis=0), exact_draws.std(axis=0), atol=0.025)
    for j in range(3):
        ks_stat, ks_p = stats.ks_2samp(hmc_draws[::4, j], exact_draws[::4, j])
        assert ks_p > 0.01, f"Coord {j} KS p={ks_p:.4f} (stat={ks_stat:.4f})"


def test_2x2_diamond_symmetry_under_10_01_label_swap() -> None:
    h_post = 5
    psi = orthonormal_exposure_basis(h_post, 2) / math.sqrt(float(h_post))
    c_cover = np.array(
        [
            [1.0, 0.0, 0.0],   # 10 - 00
            [0.0, 1.0, 0.0],   # 01 - 00
            [-1.0, 0.0, 1.0],  # 11 - 10
            [0.0, -1.0, 1.0],  # 11 - 01
        ],
        dtype=np.float64,
    )
    spec = _OutcomeConePriorSpec(
        c_cover=c_cover,
        psi=psi,
        horizon=h_post,
        monotone=True,
        period_positivity=True,
        scale_prior="local_halft",
    )
    assert spec.is_2x2_diamond
    # Prior covariance at unit scales gives identical marginal variance to 10 and 01, and 1.0 to top arm 11
    cov0 = np.linalg.inv(spec.lambda0_unit)
    assert math.isclose(cov0[0, 0], cov0[1, 1], rel_tol=1e-12)
    assert math.isclose(cov0[2, 2], 1.0, rel_tol=1e-12)
    # Normalizer is invariant under swapping edge scales (e_10_00 <-> e_01_00 and e_11_10 <-> e_11_01)
    lam_a = np.array([0.7, 1.4, 0.9, 1.8])
    lam_b = np.array([1.4, 0.7, 1.8, 0.9])
    assert math.isclose(spec.log_normalizer(0.8, lam_a), spec.log_normalizer(0.8, lam_b), rel_tol=1e-12)


def test_huang_wand_scale_mixture_matches_exact_halft3_density_and_quadrature() -> None:
    s0 = 0.65
    # Verify marginal density of sigma ~ Half-t_3(s0) via numerical quadrature of p(sigma | a) p(a)
    # sigma^2 | a ~ InvGamma(3/2, 3/a), a ~ InvGamma(1/2, 1/s0^2)
    ht3 = stats.halfnorm  # compare against exact Student-t df=3 formula
    for sig in (0.2, 0.65, 1.3):
        sig2 = sig * sig
        exact_pdf = 2.0 * stats.t.pdf(sig / s0, df=3.0) / s0

        def integrand(a: float) -> float:
            # p(sigma | a) = 2 * sigma * InvGamma(sigma^2; 3/2, 3/a)
            p_sig2_given_a = stats.invgamma.pdf(sig2, a=1.5, scale=3.0 / a)
            p_a = stats.invgamma.pdf(a, a=0.5, scale=1.0 / (s0 * s0))
            return 2.0 * sig * p_sig2_given_a * p_a

        quad_pdf, _ = integrate.quad(integrand, 0.0, np.inf, epsabs=1e-10, epsrel=1e-10)
        assert math.isclose(quad_pdf, exact_pdf, rel_tol=1e-7)

    # Empirical check of _sample_invgamma Huang-Wand hierarchy
    rng = np.random.default_rng(99)
    n_mc = 40_000
    a_draws = np.array([_sample_invgamma(0.5, 1.0 / (s0 * s0), rng) for _ in range(n_mc)])
    sig_draws = np.array([math.sqrt(_sample_invgamma(1.5, 3.0 / a, rng)) for a in a_draws])
    ref_draws = s0 * np.abs(rng.standard_t(df=3.0, size=n_mc))
    ks_stat, ks_p = stats.ks_2samp(sig_draws, ref_draws)
    assert ks_p > 0.05, f"Huang-Wand Half-t3 KS p={ks_p:.4f} (stat={ks_stat:.4f})"


def test_asis_interweaving_preserves_exact_2d_cone_posterior_stationarity() -> None:
    """Verify centered + ASIS interweaved sampler matches exact rejection sampling on a 2-D cone."""
    rng = np.random.default_rng(202604)
    h_post = 4
    psi = orthonormal_exposure_basis(h_post, 2) / math.sqrt(float(h_post))
    c_cover = np.array([[1.0]], dtype=np.float64)  # A = 2 (1 non-control arm, K' = 2)
    spec = _OutcomeConePriorSpec(
        c_cover=c_cover,
        psi=psi,
        horizon=h_post,
        monotone=True,
        period_positivity=True,
        scale_prior="global_halft",
    )
    s0 = 0.50
    y_obs = np.array([0.35, 0.08], dtype=np.float64)
    sd_obs = 0.25
    wx = np.eye(2, dtype=np.float64) / sd_obs
    wy = y_obs / sd_obs
    prec_lik = wx.T @ wx
    h_lik = wx.T @ wy
    interior_ray = np.array([1.0, 0.0], dtype=np.float64)

    # 1. Exact prior-proposal importance reference draws from joint posterior p(theta, sigma, kappa | y)
    n_prop = 200_000
    sig_p = s0 * np.abs(rng.standard_t(df=1.0, size=n_prop))
    kap_p = s0 * np.abs(rng.standard_t(df=3.0, size=n_prop))
    r_ratio = np.maximum(kap_p / np.maximum(sig_p, 1e-30), 1e-30)
    phi_lo = -np.arctan(1.0 / np.maximum(r_ratio * spec.m_plus, 1e-30))
    phi_hi = np.arctan(1.0 / np.maximum(r_ratio * spec.m_minus, 1e-30))
    rad = np.sqrt(-2.0 * np.log(np.clip(rng.uniform(size=n_prop), 1e-15, 1.0)))
    ang = rng.uniform(phi_lo, phi_hi)
    scale_fac = math.sqrt(float(h_post) / spec.c_scale)
    th0_p = rad * np.cos(ang) * sig_p * scale_fac
    th1_p = rad * np.sin(ang) * kap_p * scale_fac

    log_lik = -0.5 * (((th0_p - y_obs[0]) / sd_obs) ** 2 + ((th1_p - y_obs[1]) / sd_obs) ** 2)
    log_lik -= np.max(log_lik)
    weights = np.exp(log_lik)
    weights /= np.sum(weights)
    ref_idx = rng.choice(n_prop, size=4000, replace=True, p=weights)
    ref_th0 = th0_p[ref_idx]
    ref_sig = sig_p[ref_idx]
    ref_kap = kap_p[ref_idx]

    # 2. Run centered + ASIS interweaved chain
    theta = np.array([0.4, 0.0], dtype=np.float64)
    sig_sq = s0 * s0
    kap_sq = s0 * s0
    a_sig = s0 * s0
    a_kap = s0 * s0
    lam = np.ones(1, dtype=np.float64)
    b_loc = np.ones(1, dtype=np.float64)

    mcmc_th0 = np.empty(2500, dtype=np.float64)
    mcmc_sig = np.empty(2500, dtype=np.float64)
    mcmc_kap = np.empty(2500, dtype=np.float64)
    for step in range(2800):
        l0, _ = spec.precision_and_logdet(lam)
        prec_prior = np.diag([float(l0[0, 0]) / (float(h_post) * sig_sq), float(l0[0, 0]) / (float(h_post) * kap_sq)])
        prec_post = prec_lik + prec_prior
        mean_post = np.linalg.solve(prec_post, h_lik)
        theta = _sample_truncated_mvn(mean_post, prec_post, spec.a_cone, theta, interior_ray, rng, n_hmc_steps=1)
        w_res = wy - wx @ theta
        sig_sq, kap_sq, a_sig, a_kap, lam, b_loc = _update_scales_outcome(
            theta, spec, sig_sq, kap_sq, a_sig, a_kap, lam, b_loc, s0, None, True, rng
        )
        sig_sq, kap_sq, a_sig, a_kap, w_res = _asis_interweave_outcome_scales(
            m_idx=0,
            theta=theta,
            w_res=w_res,
            wx=wx,
            spec=spec,
            sigma_sq_m=sig_sq,
            kappa_sq_m=kap_sq,
            a_sigma_m=a_sig,
            a_kappa_m=a_kap,
            local_scales_m=lam,
            s0_scale_m=s0,
            fixed_sigma_m=None,
            dynamic_shrinkage=True,
            rng=rng,
        )
        if step >= 300:
            mcmc_th0[step - 300] = theta[0]
            mcmc_sig[step - 300] = math.sqrt(sig_sq)
            mcmc_kap[step - 300] = math.sqrt(kap_sq)

    assert stats.ks_2samp(mcmc_th0[::2], ref_th0[::2]).pvalue > 0.01
    assert stats.ks_2samp(mcmc_sig[::2], ref_sig[::2]).pvalue > 0.01
    assert stats.ks_2samp(mcmc_kap[::2], ref_kap[::2]).pvalue > 0.01


def test_polar_tree_and_diamond_prior_sampler_match_rejection() -> None:
    """Verify exact polar tree prior sampler and diamond TN-HMC fallback match naive rejection at r_m=1.5."""
    rng = np.random.default_rng(202605)
    h = 4
    psi = orthonormal_exposure_basis(h, 2) / math.sqrt(float(h))
    sig_m, kap_m = 0.4, 0.6  # r_m = 1.5

    # 1. Polar tree sampler on 5-arm chain
    topo_tree = _resolve_arms_and_poset(
        (("A1", "A0"), ("A2", "A1"), ("A3", "A2"), ("A4", "A3")),
        ("A0", "A1", "A2", "A3", "A4"),
        ("A0", "A1", "A2", "A3", "A4"),
        None,
    )
    spec_tree = _OutcomeConePriorSpec(
        c_cover=topo_tree["c_cover"],
        psi=psi,
        horizon=h,
        monotone=True,
        period_positivity=True,
        scale_prior="local_halft",
    )
    lam_m = np.array([0.8, 1.3, 0.9, 1.5], dtype=np.float64)
    prec_l0, _ = spec_tree.precision_and_logdet(lam_m)
    chol_l0 = np.linalg.cholesky(np.linalg.inv(prec_l0 / float(h)))

    rej_draws: list[np.ndarray] = []
    while len(rej_draws) < 2500:
        z0 = rng.standard_normal((4096, 4)) @ chol_l0.T * sig_m
        z1 = rng.standard_normal((4096, 4)) @ chol_l0.T * kap_m
        cand = np.hstack([z0, z1])
        ok = np.all((cand @ spec_tree.a_cone.T) >= 0.0, axis=1)
        rej_draws.extend(cand[ok])
    rej_tree = np.asarray(rej_draws[:2500])

    assert spec_tree.c_cover_inv is not None
    edge_sd0 = sig_m * lam_m * math.sqrt(float(h) / spec_tree.c_scale)
    edge_sdk = kap_m * lam_m * math.sqrt(float(h) / spec_tree.c_scale)
    r_m = kap_m / sig_m
    phi_lo = -math.atan(1.0 / (r_m * spec_tree.m_plus))
    phi_hi = math.atan(1.0 / (r_m * spec_tree.m_minus))
    polar_tree = np.empty((2500, 8), dtype=np.float64)
    for i in range(2500):
        rho = np.sqrt(-2.0 * np.log(rng.uniform(size=4)))
        phi = rng.uniform(phi_lo, phi_hi, size=4)
        u_modes = np.vstack([rho * np.cos(phi) * edge_sd0, rho * np.sin(phi) * edge_sdk])
        polar_tree[i] = (u_modes @ spec_tree.c_cover_inv.T).reshape(-1)

    for j in range(8):
        assert stats.ks_2samp(rej_tree[:, j], polar_tree[:, j]).pvalue > 0.01

    # 2. Diamond TN-HMC fallback vs naive rejection
    topo_diam = _resolve_arms_and_poset(
        (("10", "00"), ("01", "00"), ("11", "10"), ("11", "01")),
        ("00", "10", "01", "11"),
        None,
        (("11", "10"), ("11", "01"), ("10", "00"), ("01", "00")),
    )
    spec_diam = _OutcomeConePriorSpec(
        c_cover=topo_diam["c_cover"],
        psi=psi,
        horizon=h,
        monotone=True,
        period_positivity=True,
        scale_prior="global_halft",
    )
    prec_d0, _ = spec_diam.precision_and_logdet(np.ones(4))
    chol_d0 = np.linalg.cholesky(np.linalg.inv(prec_d0 / float(h)))
    rej_d: list[np.ndarray] = []
    while len(rej_d) < 2000:
        z0 = rng.standard_normal((4096, 3)) @ chol_d0.T * sig_m
        z1 = rng.standard_normal((4096, 3)) @ chol_d0.T * kap_m
        cand = np.hstack([z0, z1])
        ok = np.all((cand @ spec_diam.a_cone.T) >= 0.0, axis=1)
        rej_d.extend(cand[ok])
    rej_diam = np.asarray(rej_d[:2000])

    prec_prior = np.zeros((6, 6), dtype=np.float64)
    prec_prior[:3, :3] = prec_d0 / (float(h) * sig_m**2)
    prec_prior[3:, 3:] = prec_d0 / (float(h) * kap_m**2)
    ray_d = np.array([sig_m, sig_m, 2.0 * sig_m, 0.0, 0.0, 0.0], dtype=np.float64)
    th_cur = ray_d.copy()
    hmc_diam = np.empty((1000, 6), dtype=np.float64)
    for s in range(4200):
        th_cur = _sample_truncated_mvn(
            np.zeros(6, dtype=np.float64), prec_prior, spec_diam.a_cone, th_cur, ray_d, rng, n_hmc_steps=3
        )
        if s >= 200 and (s - 200) % 4 == 0:
            hmc_diam[(s - 200) // 4] = th_cur
    for j in range(6):
        assert stats.ks_2samp(rej_diam[:, j], hmc_diam[:, j]).pvalue > 0.01


