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
    _slice_sample_1d_with_logp,
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
        use_local_scales=False,
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
            theta_m=theta,
            spec=spec,
            sigma_sq=sig_sq,
            kappa_sq=kap_sq,
            a_sigma=a_sig,
            a_kappa=a_kap,
            local_scales=lam,
            b_local=b_loc,
            log_s0_eff=math.log(s0),
            log_sk_eff=math.log(s0),
            nu_sig=1.0,
            nu_kap=3.0,
            nu_loc=3.0,
            fixed_sigma=None,
            dynamic_shrinkage=True,
            rng=rng,
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
            log_s0_eff_m=math.log(s0),
            log_sk_eff_m=math.log(s0),
            nu_sig=1.0,
            nu_kap=3.0,
            nu_loc=3.0,
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



def test_pakman_paninski_bounce_cap_rejection_keeps_exact_truncated_normal() -> None:
    """A tiny fixed bounce budget forces frequent rejections; reject-and-stay must stay exact."""
    rng = np.random.default_rng(20261029)
    # Narrow wedge {z : f_j . z + g_j >= 0} around the +z_1 axis with apex at (-0.5, 0):
    # trajectories bounce between its two walls many times, so max_events=3 rejects often.
    half_angle = 0.15
    f_walls = np.array(
        [
            [math.sin(half_angle), math.cos(half_angle)],
            [math.sin(half_angle), -math.cos(half_angle)],
        ],
        dtype=np.float64,
    )
    g_walls = np.full(2, 0.5 * math.sin(half_angle), dtype=np.float64)

    raw = rng.standard_normal((2_000_000, 2))
    exact = raw[np.all(raw @ f_walls.T + g_walls >= 0.0, axis=1)][:20_000]
    assert exact.shape[0] == 20_000

    z_cur = np.array([1.0, 0.0], dtype=np.float64)
    n_keep = 30_000
    chain = np.empty((n_keep, 2), dtype=np.float64)
    n_stay = 0
    for s in range(n_keep + 500):
        z_new = rkma_mod._pakman_paninski_hmc_whitened(z_cur, f_walls, g_walls, rng, n_steps=1, max_events=3)
        n_stay += int(np.array_equal(z_new, z_cur))
        z_cur = z_new
        if s >= 500:
            chain[s - 500] = z_cur

    reject_frac = n_stay / float(n_keep + 500)
    assert 0.05 < reject_frac < 0.95, f"bounce cap should bind for a meaningful test; got {reject_frac:.3f}"
    assert np.all(chain @ f_walls.T + g_walls >= -1e-12)
    np.testing.assert_allclose(chain.mean(axis=0), exact.mean(axis=0), atol=0.03)
    np.testing.assert_allclose(chain.std(axis=0), exact.std(axis=0), atol=0.03)
    for j in range(2):
        assert stats.ks_2samp(chain[::10, j], exact[::7, j]).pvalue > 0.01


def test_slice_sampler_logp_reuse_is_bitwise_identical_to_reevaluation() -> None:
    def log_prob(x: float) -> float:
        return -0.5 * (x - 0.3) ** 2 / 0.04 + 2.0 * math.log(max(x, 1e-300)) if x > 0 else -math.inf

    for seed in range(25):
        x0 = 0.05 + 0.04 * seed
        x_plain = rkma_mod._slice_sample_1d(log_prob, x0, np.random.default_rng(seed), width=0.35, max_steps=3, bounds=(1e-5, 2.0))
        x_reuse, lp_reuse = rkma_mod._slice_sample_1d_with_logp(
            log_prob,
            x0,
            np.random.default_rng(seed),
            width=0.35,
            max_steps=3,
            bounds=(1e-5, 2.0),
            logp_init=log_prob(x0),
        )
        assert x_plain == x_reuse
        assert lp_reuse == log_prob(x_reuse)


def test_lapack_collapsed_gaussian_loglik_matches_scipy_logpdf() -> None:
    """dpotrf/dtrtrs collapsed log-likelihood (used for cached slice log densities) is exact."""
    rng = np.random.default_rng(7)
    n_q, k_dim, n_xi = 9, 2, 8
    g_xi = rng.standard_normal((n_q * k_dim, n_xi))
    a_xi = rng.standard_normal((n_xi, n_xi))
    sigma_xi = a_xi @ a_xi.T / n_xi
    f_blocks = np.stack([np.diag(rng.uniform(0.05, 0.4, k_dim)) for _ in range(n_q)])
    omega = rng.gamma(2.0, 0.5, n_q)
    omega_cov = rkma_mod._build_omega_cov(g_xi, sigma_xi, f_blocks, omega, k_dim)
    resid = rng.standard_normal(n_q * k_dim)

    chol = rkma_mod._fast_cholesky(omega_cov, floor=1e-13)
    v_sol = rkma_mod._solve_lower(chol, resid)
    loglik = -float(np.sum(np.log(np.diag(chol)))) - 0.5 * float(v_sol @ v_sol)
    ref = float(stats.multivariate_normal(mean=np.zeros(n_q * k_dim), cov=omega_cov).logpdf(resid))
    ref += 0.5 * (n_q * k_dim) * math.log(2.0 * math.pi)
    assert abs(loglik - ref) < 1e-10

    x_back = rkma_mod._solve_lower(chol, rkma_mod._solve_lower(chol, resid), transpose=True)
    np.testing.assert_allclose(omega_cov @ x_back, resid, atol=1e-10)
    np.testing.assert_allclose(rkma_mod._inv_lower(chol) @ chol, np.eye(n_q * k_dim), atol=1e-10)


def test_iw2_marginal_z_slice_and_imh_matches_bartlett_correlation() -> None:
    """Verify exact 1-D IW_2(nu, Psi) marginal z=atanh(r) log-density and IMH + z-slice stationarity."""
    # 1. Quadrature check of _log_iw2_marginal_z across small and large df (no overflow at nu=90)
    for df_iw, r_psi in ((5.0, 0.75), (25.0, -0.60), (90.0, 0.82)):
        iw2_t2, iw2_w_vt = rkma_mod._precompute_iw2_quadrature(df_iw)
        z_ref = 0.0
        lp_ref = rkma_mod._log_iw2_marginal_z(z_ref, df_iw, r_psi, iw2_t2, iw2_w_vt)

        def _exact_unnorm_logp_z(z_val: float) -> float:
            r_val = math.tanh(z_val)
            a_val = r_psi * r_val
            log_1ma = math.log1p(-a_val)
            integ, _ = integrate.quad(
                lambda psi: math.exp(
                    (df_iw - 1.0) * math.log(max(math.sin(psi), 1e-300))
                    - df_iw * (math.log1p(-a_val * math.sin(psi)) - log_1ma)
                ),
                0.0,
                math.pi,
                epsabs=1e-13,
                epsrel=1e-13,
            )
            # p(r) propto (1 - r^2)^{(df - 3)/2} * integ; |dr/dz| = 1 - r^2 => p(z) propto sech^{df - 1}(z) * integ
            return (
                -(df_iw - 1.0) * np.logaddexp(z_val, -z_val)
                + (df_iw - 1.0) * math.log(2.0)
                - df_iw * log_1ma
                + math.log(integ)
            )

        q_ref = _exact_unnorm_logp_z(z_ref)
        for z_test in (-1.35, -0.4, 0.55, 1.50):
            diff_fast = rkma_mod._log_iw2_marginal_z(z_test, df_iw, r_psi, iw2_t2, iw2_w_vt) - lp_ref
            diff_quad = _exact_unnorm_logp_z(z_test) - q_ref
            assert abs(diff_fast - diff_quad) < 1e-10, (df_iw, r_psi, z_test, diff_fast, diff_quad)

    # 2. Stationarity of IMH + z-slice kernel under both flat likelihood and tilted likelihood
    rng = np.random.default_rng(20261030)
    for df_iw, r_psi, tilt in ((6.0, 0.75, 0.0), (25.0, -0.60, 1.6)):
        psi_mat = np.array([[4.0, 4.0 * r_psi * 1.5], [4.0 * r_psi * 1.5, 9.0]], dtype=np.float64)
        chol_inv = np.linalg.cholesky(np.linalg.inv(psi_mat))
        iw2_t2, iw2_w_vt = rkma_mod._precompute_iw2_quadrature(df_iw)

        def _loglik_r(r_val: float) -> float:
            return -tilt * ((r_val - 0.25) ** 2)

        # Reference draws via Bartlett IW_2 + importance resampling
        n_prop = 60_000
        bart_r = np.array(
            [rkma_mod._sample_iw_correlation(df_iw, psi_mat, rng, chol_inv=chol_inv)[0, 1] for _ in range(n_prop)]
        )
        bart_r = np.clip(bart_r, -math.tanh(6.5), math.tanh(6.5))
        if tilt == 0.0:
            ref_r = bart_r[:5000]
        else:
            lw = -tilt * ((bart_r - 0.25) ** 2)
            lw -= np.max(lw)
            w = np.exp(lw)
            w /= np.sum(w)
            ref_r = bart_r[rng.choice(n_prop, size=5000, replace=True, p=w)]

        r_cur = 0.0
        ll_cur = _loglik_r(r_cur)
        chain_r = np.empty(3500, dtype=np.float64)
        for step in range(3800):
            # Exact Step 6a: IMH from IW prior
            r_prop = float(rkma_mod._sample_iw_correlation(df_iw, psi_mat, rng, chol_inv=chol_inv)[0, 1])
            ll_prop = _loglik_r(r_prop)
            if math.log(max(float(rng.uniform()), 1e-300)) < ll_prop - ll_cur:
                r_cur = r_prop
                ll_cur = ll_prop
            # Exact Step 6b: 1-D slice step on z = atanh(r)
            z_cur = math.atanh(max(-1.0 + 1e-12, min(1.0 - 1e-12, r_cur)))

            def _logp_z(z_val: float) -> float:
                rv = math.tanh(z_val)
                return rkma_mod._log_iw2_marginal_z(z_val, df_iw, r_psi, iw2_t2, iw2_w_vt) + _loglik_r(rv)

            z_new, _ = rkma_mod._slice_sample_1d_with_logp(
                _logp_z,
                z_cur,
                rng,
                width=0.5,
                max_steps=6,
                bounds=(-6.5, 6.5),
                logp_init=rkma_mod._log_iw2_marginal_z(z_cur, df_iw, r_psi, iw2_t2, iw2_w_vt) + ll_cur,
            )
            r_cur = math.tanh(z_new)
            ll_cur = _loglik_r(r_cur)
            if step >= 300:
                chain_r[step - 300] = r_cur

        ks_res = stats.ks_2samp(chain_r[::2], ref_r[::2])
        assert ks_res.pvalue > 0.01, f"IW2 z-slice KS failed for df={df_iw}, r={r_psi}, tilt={tilt}: p={ks_res.pvalue:.4f}"


def test_hierarchical_global_local_scale_prior_stationarity() -> None:
    """Verify per-outcome (tau_m, sigma_m = tau_m * s_{0,m}, kappa_m = gamma_m * sigma_m, theta_m) stationarity."""
    rng = np.random.default_rng(20261031)
    h_post = 4
    psi = orthonormal_exposure_basis(h_post, 2) / math.sqrt(float(h_post))
    c_cover = np.eye(2, dtype=np.float64)  # A = 3 (2 non-control increments, K' = 2, M = 2)
    specs = [
        _OutcomeConePriorSpec(
            c_cover=c_cover,
            psi=psi,
            horizon=h_post,
            monotone=True,
            period_positivity=True,
            scale_prior="global_halft",
            use_local_scales=False,
        )
        for _ in range(2)
    ]
    s0_scales = [0.50, 0.80]
    log_s0 = [math.log(s) for s in s0_scales]
    d_m = 4
    d_tot = 8
    y_obs = np.array([0.55, 0.70, 0.05, -0.04, 0.04, 0.03, 0.01, -0.01], dtype=np.float64)
    sd_obs = 0.35
    wx = np.eye(d_tot, dtype=np.float64) / sd_obs
    wy = y_obs / sd_obs
    prec_lik = wx.T @ wx
    h_lik = wx.T @ wy

    a_cone = np.zeros((2 * specs[0].a_cone.shape[0], d_tot), dtype=np.float64)
    nw = specs[0].a_cone.shape[0]
    a_cone[:nw, :d_m] = specs[0].a_cone
    a_cone[nw:, d_m:] = specs[1].a_cone
    interior_ray = np.array([0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0], dtype=np.float64)

    tau0 = rkma_mod._SCALE_PRIOR_TAU0
    nu_gamma = rkma_mod._SCALE_PRIOR_NU_GAMMA
    gamma0 = rkma_mod._SCALE_PRIOR_GAMMA0
    wt_b = rkma_mod._SCALE_PRIOR_W_TAU_BOUND
    ug_b = rkma_mod._SCALE_PRIOR_U_GAMMA_BOUND
    nu_loc = rkma_mod._SCALE_PRIOR_NU_LOCAL

    # Verify normalized two-group per-outcome log-scale density integrates to 1
    w_grid = np.linspace(-25.0, 25.0, 20001)
    p_grid = np.array([math.exp(rkma_mod._log_outcome_tau_prior(float(w))) for w in w_grid])
    assert abs(float(np.trapezoid(p_grid, w_grid)) - 1.0) < 1e-5

    def _draw_box_halft(scale_val: float, nu: float, n_s: int, lo: float, hi: float) -> np.ndarray:
        out = np.empty(n_s, dtype=np.float64)
        rem = np.arange(n_s)
        while rem.size > 0:
            prop = scale_val * np.abs(rng.standard_t(df=nu, size=rem.size))
            w_p = np.log(np.maximum(prop, 1e-300))
            ok = (w_p > lo) & (w_p < hi)
            out[rem[ok]] = prop[ok]
            rem = rem[~ok]
        return out

    def _draw_box_outcome_tau(n_s: int) -> np.ndarray:
        out = np.empty(n_s, dtype=np.float64)
        rem = np.arange(n_s)
        while rem.size > 0:
            prop = rkma_mod._SCALE_PRIOR_TAU0 * np.sqrt(
                np.maximum(
                    rng.f(rkma_mod._SCALE_PRIOR_ALPHA_TAU, rkma_mod._SCALE_PRIOR_NU_TAU, size=rem.size),
                    1e-300,
                )
            )
            w_p = np.log(np.maximum(prop, 1e-300))
            ok = (w_p > -wt_b) & (w_p < wt_b)
            out[rem[ok]] = prop[ok]
            rem = rem[~ok]
        return out

    # 1. Exact importance-sampled reference draws from joint posterior with independent tau_m
    n_prop = 220_000
    tau_p = np.vstack([_draw_box_outcome_tau(n_prop) for _ in range(2)])
    sig_p = np.vstack([tau_p[m] * s0_scales[m] for m in range(2)])
    gam_p = np.vstack([_draw_box_halft(gamma0, nu_gamma, n_prop, -ug_b, ug_b) for _ in range(2)])
    kap_p = sig_p * gam_p
    th_p = np.zeros((n_prop, d_tot), dtype=np.float64)
    scale_fac = math.sqrt(float(h_post) / specs[0].c_scale)
    for m in range(2):
        r_m = gam_p[m]
        phi_lo = -np.arctan(1.0 / np.maximum(r_m * specs[m].m_plus, 1e-30))
        phi_hi = np.arctan(1.0 / np.maximum(r_m * specs[m].m_minus, 1e-30))
        for e in range(2):
            rad = np.sqrt(-2.0 * np.log(np.clip(rng.uniform(size=n_prop), 1e-15, 1.0)))
            ang = rng.uniform(phi_lo, phi_hi)
            th_p[:, m * d_m + e] = rad * np.cos(ang) * sig_p[m] * scale_fac
            th_p[:, m * d_m + 2 + e] = rad * np.sin(ang) * kap_p[m] * scale_fac

    ll = -0.5 * np.sum(((th_p - y_obs[None, :]) / sd_obs) ** 2, axis=1)
    ll -= np.max(ll)
    w_is = np.exp(ll)
    w_is /= np.sum(w_is)
    ess_is = 1.0 / float(np.sum(w_is * w_is))
    assert ess_is > 1500.0, f"Importance reference ESS too small: {ess_is:.1f}"
    idx_ref = rng.choice(n_prop, size=3500, replace=True, p=w_is)

    # 2. Run per-outcome MCMC transitions + ASIS interweaving
    theta = interior_ray.copy()
    sigma_sq = np.array([(tau0 * s) ** 2 for s in s0_scales], dtype=np.float64)
    kappa_sq = np.array([(tau0 * gamma0 * s) ** 2 for s in s0_scales], dtype=np.float64)
    a_sig = sigma_sq.copy()
    a_kap = kappa_sq.copy()
    lam = [np.ones(2, dtype=np.float64), np.ones(2, dtype=np.float64)]
    b_loc = [np.ones(2, dtype=np.float64), np.ones(2, dtype=np.float64)]

    tr_tau0 = np.empty(3000, dtype=np.float64)
    tr_tau1 = np.empty(3000, dtype=np.float64)
    tr_sig0 = np.empty(3000, dtype=np.float64)
    tr_sig1 = np.empty(3000, dtype=np.float64)
    tr_kap0 = np.empty(3000, dtype=np.float64)
    tr_th0 = np.empty(3000, dtype=np.float64)

    for step in range(3400):
        prec_prior = np.zeros((d_tot, d_tot), dtype=np.float64)
        for m in range(2):
            l0_m, _ = specs[m].precision_and_logdet(lam[m])
            bm = m * d_m
            prec_prior[bm : bm + 2, bm : bm + 2] = l0_m / (float(h_post) * sigma_sq[m])
            prec_prior[bm + 2 : bm + 4, bm + 2 : bm + 4] = l0_m / (float(h_post) * kappa_sq[m])
        prec_post = prec_lik + prec_prior
        mean_post = np.linalg.solve(prec_post, h_lik)
        theta = _sample_truncated_mvn(mean_post, prec_post, a_cone, theta, interior_ray, rng, n_hmc_steps=1)

        w_res = wy - wx @ theta
        for m in range(2):
            sigma_sq[m], kappa_sq[m], a_sig[m], a_kap[m], lam[m], b_loc[m] = _update_scales_outcome(
                theta_m=theta[m * d_m : (m + 1) * d_m],
                spec=specs[m],
                sigma_sq=float(sigma_sq[m]),
                kappa_sq=float(kappa_sq[m]),
                a_sigma=float(a_sig[m]),
                a_kappa=float(a_kap[m]),
                local_scales=lam[m],
                b_local=b_loc[m],
                log_s0_eff=log_s0[m],
                log_sk_eff=log_s0[m],
                nu_sig=3.0,
                nu_kap=nu_gamma,
                nu_loc=nu_loc,
                fixed_sigma=None,
                dynamic_shrinkage=True,
                rng=rng,
                coupled_p1=True,
            )
            sigma_sq[m], kappa_sq[m], a_sig[m], a_kap[m], w_res = _asis_interweave_outcome_scales(
                m_idx=m,
                theta=theta,
                w_res=w_res,
                wx=wx,
                spec=specs[m],
                sigma_sq_m=float(sigma_sq[m]),
                kappa_sq_m=float(kappa_sq[m]),
                a_sigma_m=float(a_sig[m]),
                a_kappa_m=float(a_kap[m]),
                local_scales_m=lam[m],
                log_s0_eff_m=log_s0[m],
                log_sk_eff_m=log_s0[m],
                nu_sig=3.0,
                nu_kap=nu_gamma,
                nu_loc=nu_loc,
                fixed_sigma_m=None,
                dynamic_shrinkage=True,
                rng=rng,
                coupled_p1=True,
            )
        if step >= 400:
            tr_tau0[step - 400] = math.sqrt(sigma_sq[0]) / s0_scales[0]
            tr_tau1[step - 400] = math.sqrt(sigma_sq[1]) / s0_scales[1]
            tr_sig0[step - 400] = math.sqrt(sigma_sq[0])
            tr_sig1[step - 400] = math.sqrt(sigma_sq[1])
            tr_kap0[step - 400] = math.sqrt(kappa_sq[0])
            tr_th0[step - 400] = theta[0]

    assert stats.ks_2samp(tr_tau0[::2], tau_p[0, idx_ref][::2]).pvalue > 0.01
    assert stats.ks_2samp(tr_tau1[::2], tau_p[1, idx_ref][::2]).pvalue > 0.01
    assert stats.ks_2samp(tr_sig0[::2], sig_p[0, idx_ref][::2]).pvalue > 0.01
    assert stats.ks_2samp(tr_sig1[::2], sig_p[1, idx_ref][::2]).pvalue > 0.01
    assert stats.ks_2samp(tr_kap0[::2], kap_p[0, idx_ref][::2]).pvalue > 0.01
    assert stats.ks_2samp(tr_th0[::2], th_p[idx_ref, 0][::2]).pvalue > 0.01


def test_per_outcome_independent_scale_zero_contamination() -> None:
    """Verify an active outcome does not inflate a null outcome's scale when R = I and no composite links them."""
    arms = ("C", "B", "A")
    rng = np.random.default_rng(20261100)
    fits_active = {
        ("B", "C"): _make_stage1_fit(rng, np.ones(4) * 1.20),
        ("A", "B"): _make_stage1_fit(rng, np.ones(4) * 1.20),
        ("A", "C"): _make_stage1_fit(rng, np.ones(4) * 2.40),
    }
    fits_null = {
        ("B", "C"): _make_stage1_fit(rng, np.zeros(4)),
        ("A", "B"): _make_stage1_fit(rng, np.zeros(4)),
        ("A", "C"): _make_stage1_fit(rng, np.zeros(4)),
    }

    solo_null = reconcile_multi_arm_draws(
        fits_null, arms=arms, order=arms, monotone=True, seed=20261101
    )
    joint = reconcile_multi_arm_multi_outcome_draws(
        {"active": fits_active, "null": fits_null},
        cross_outcome_residuals=None,
        arms=arms,
        order=arms,
        monotone=True,
        seed=20261102,
    )

    tau_act = joint.by_outcome["active"].tau_draws
    tau_nul_joint = joint.by_outcome["null"].tau_draws
    tau_nul_solo = solo_null.tau_draws
    assert tau_act is not None and tau_nul_joint is not None and tau_nul_solo is not None
    assert float(np.median(tau_act)) > 4.0 * float(np.median(tau_nul_joint))
    assert abs(float(np.median(tau_nul_joint)) - float(np.median(tau_nul_solo))) < 0.04
    mean_top_joint = float(np.mean(joint.by_outcome["null"].contrast_mode1_draws[("A", "C")]))
    mean_top_solo = float(np.mean(solo_null.contrast_mode1_draws[("A", "C")]))
    assert abs(mean_top_joint - mean_top_solo) < 0.015


def test_slice_sample_1d_boundary_detailed_balance() -> None:
    """Verify bounded Neal (2003) slice sampler preserves exact stationarity near both endpoints."""
    rng = np.random.default_rng(20261105)
    bounds = (1e-5, 1.0 - 1e-5)

    for a_par, b_par, width in ((2.0, 10.0, 0.25), (5.0, 1.5, 0.35)):
        def log_beta(x: float, a: float = a_par, b: float = b_par) -> float:
            return (a - 1.0) * math.log(x) + (b - 1.0) * math.log1p(-x)

        x_cur = a_par / (a_par + b_par)
        lp_cur: float | None = None
        draws = np.empty(6000, dtype=np.float64)
        for step in range(25_000):
            x_cur, lp_cur = _slice_sample_1d_with_logp(
                log_beta, x_cur, rng, width=width, max_steps=3, bounds=bounds, logp_init=lp_cur
            )
            if step >= 1000 and (step - 1000) % 4 == 0:
                draws[(step - 1000) // 4] = x_cur

        ks_res = stats.kstest(draws, "beta", args=(a_par, b_par))
        assert ks_res.pvalue > 0.05, (
            f"Bounded slice sampler failed KS stationarity on Beta({a_par}, {b_par}): p={ks_res.pvalue:.4g}"
        )
        med = float(stats.beta.ppf(0.5, a_par, b_par))
        assert abs(float(np.mean(draws < med)) - 0.5) < 0.025


def test_sample_truncated_mvn_ultrasubmicro_scale_invariance() -> None:
    """Verify whitened cone-slack handling in `_sample_truncated_mvn` stays exact at `sigma = 1e-10`."""
    rng = np.random.default_rng(20261106)
    scale = 1e-10
    prec = np.eye(2, dtype=np.float64) / (scale * scale)
    mean = np.zeros(2, dtype=np.float64)
    a_cone = np.array([[1.0, -0.5], [1.0, 0.5]], dtype=np.float64)
    ray = np.array([1.0, 0.0], dtype=np.float64)

    th = np.array([1e-16, 0.0], dtype=np.float64)  # well below old 1e-11 unwhitened threshold
    draws = np.empty((2000, 2), dtype=np.float64)
    for step in range(4400):
        th = _sample_truncated_mvn(mean, prec, a_cone, th, ray, rng, n_hmc_steps=1)
        if step >= 400 and (step - 400) % 2 == 0:
            draws[(step - 400) // 2] = th / scale

    # Exact rejection reference at unit scale
    rej_list: list[np.ndarray] = []
    while len(rej_list) < 4000:
        cand = rng.standard_normal((4096, 2))
        ok = np.all((cand @ a_cone.T) >= 0.0, axis=1)
        rej_list.extend(cand[ok])
    ref = np.asarray(rej_list[:4000])
    assert stats.ks_2samp(draws[:, 0], ref[:, 0]).pvalue > 0.05
    assert stats.ks_2samp(draws[:, 1], ref[:, 1]).pvalue > 0.05


