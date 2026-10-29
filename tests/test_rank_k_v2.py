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

"""Section 5.1 Sampler Correctness Tests and Stage E Equivalence Gate for Rank-K BCF v2."""

from __future__ import annotations

import dataclasses

import equinox as eqx
import jax
import jax.numpy as jnp
from jax import random
from longbet._x64 import enable_x64
import numpy as np
import pytest
from scipy.stats import ks_2samp

from longbet import (
    RankKConfigWarning,
    RankKLongBetConfig,
    build_causal_features,
    fit_rank_k_longbet,
    rank_k_placebo_audit,
)
from longbet._rank_k_jax import (
    _bin_mode_predictors,
    _sample_inv_gamma_64,
    _sample_truncated_normal_64,
    build_spatial_kernel_cache,
    closed_form_linear_gaussian_posterior,
    init_rank_k_mode,
    matern32_corr,
    great_circle_distance_km,
    rank_k_mode_single_step,
)


def _synthetic_staggered_panel(
    seed: int = 101,
    N: int = 80,
    T: int = 10,
    T_hist: int = 24,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, 3))
    coords = np.column_stack(
        [35.0 + 4.0 * rng.uniform(size=N), -95.0 + 8.0 * rng.uniform(size=N)]
    )
    cohorts = np.full(N, T, dtype=np.int32)
    cohorts[: N // 4] = 4
    cohorts[N // 4 : N // 2] = 6
    z = (np.arange(T)[None, :] >= cohorts[:, None]).astype(np.int32)
    t_full = np.arange(T_hist + T, dtype=np.float64)
    season = (
        1.4
        * (1.0 + 0.35 * x[:, 0:1])
        * np.sin(2.0 * np.pi * t_full[None, :] / 12.0)
    )
    trend = 0.04 * (0.5 + 0.3 * x[:, 1:2]) * t_full[None, :]
    y_full = season + trend + 0.35 * rng.normal(size=(N, T_hist + T))
    yh = y_full[:, :T_hist].copy()
    y = y_full[:, T_hist:].copy()
    tau_unit = 1.1 + 0.45 * x[:, 0] - 0.25 * x[:, 1]
    for i in range(N):
        g = cohorts[i]
        if g < T:
            s = np.arange(1, T - g + 1, dtype=np.float64)
            y[i, g:] += tau_unit[i] * (1.0 - 0.4 * np.exp(-s / 2.5))
    return y, z, x, yh, coords, cohorts


def test_closed_form_linear_gaussian_limit() -> None:
    """Section 5.1 Test 2: Closed-form linear-Gaussian limit with use_forests=False.

    The closed-form posterior oracle covers only this forest-free limit (``use_forests``
    is an internal switch of the sampler, not a config option); the forest steps, which
    use working likelihoods, are not checked against it.
    """
    rng = np.random.default_rng(1234)
    N = 40
    q = 2
    n0 = 24
    ctrl_idx = np.arange(n0, dtype=np.int32)
    trt_idx = np.arange(n0, N, dtype=np.int32)
    z_vec = np.zeros(N, dtype=np.float64)
    z_vec[trt_idx] = 1.0
    w_unit = 0.8 + 0.4 * rng.uniform(size=N)
    w_unit /= np.mean(w_unit[ctrl_idx])

    df_k = rng.normal(size=N)
    lam_mat = rng.normal(size=(N, q))
    coords = np.column_stack(
        [36.0 + 3.0 * rng.uniform(size=N), -92.0 + 3.0 * rng.uniform(size=N)]
    )
    ell_fixed = 400.0
    dist_km = great_circle_distance_km(coords, coords)
    C_corr = matern32_corr(dist_km, ell_fixed)

    sigma2_true = 0.36
    s_f2_true = 0.25
    r_fixed = 0.5
    s_zeta2_true = r_fixed * sigma2_true
    intercept_prior_sd = 10.0
    theta_prior_mean = 1.0
    theta_prior_sd = 0.5
    tau_prior_sd = 0.5

    dy_k = (
        0.4
        + 0.85 * df_k
        + lam_mat @ np.array([0.3, -0.25])
        + 1.35 * z_vec
        + rng.normal(scale=np.sqrt(sigma2_true / w_unit), size=N)
    )

    exact = closed_form_linear_gaussian_posterior(
        dy_k=dy_k,
        df_k=df_k,
        z_vec=z_vec,
        w_unit=w_unit,
        lambda_mat=lam_mat,
        ctrl_idx=ctrl_idx,
        trt_idx=trt_idx,
        sigma2=sigma2_true,
        s_f2=s_f2_true,
        s_zeta2=s_zeta2_true,
        C_corr=C_corr,
        intercept_prior_sd=intercept_prior_sd,
        theta_prior_mean=theta_prior_mean,
        theta_prior_sd=theta_prior_sd,
        tau_prior_sd=tau_prior_sd,
    )

    # Build a spatial cache with a single lengthscale ell_fixed and single ratio r_fixed,
    # and concentrated inverse-gamma priors on sigma^2 and s_f^2 so variance components stay fixed.
    spatial_cache = build_spatial_kernel_cache(
        ctrl_idx=ctrl_idx,
        trt_idx=trt_idx,
        w_unit=w_unit,
        coords=coords,
        ell_grid_km=(ell_fixed,),
        ratio_grid=(r_fixed,),
        use_spatial_gp=True,
    )
    X_binned, max_split_mu, max_split_nu = _bin_mode_predictors(
        x_mod= rng.normal(size=(N, 2)),
        ps_vec=np.full(N, 0.4),
        df_std=df_k,
        has_forecast=True,
        key=random.key(77),
    )
    big_a = 1e7
    st = init_rank_k_mode(
        X_binned=X_binned,
        max_split_mu=max_split_mu,
        max_split_nu=max_split_nu,
        dy_k=jnp.asarray(dy_k, dtype=jnp.float32),
        df_k=jnp.asarray(df_k, dtype=jnp.float32),
        z_vec=jnp.asarray(z_vec, dtype=jnp.float32),
        w_unit=jnp.asarray(w_unit, dtype=jnp.float32),
        lambda_mat=jnp.asarray(lam_mat, dtype=jnp.float32),
        ctrl_idx=ctrl_idx,
        trt_idx=trt_idx,
        num_trees_mu=10,
        num_trees_nu=10,
        sigma_nu_total=tau_prior_sd,
        has_forecast=True,
        use_forests=False,
        use_factor_loadings=True,
        use_spatial_gp=True,
        baseline_untreated_only=True,
        theta_f_clip=(-100.0, 100.0),
        theta_prior_mean=theta_prior_mean,
        theta_prior_sd=theta_prior_sd,
        sigma_prior_a=big_a,
        sigma_prior_b=big_a * sigma2_true,
        sf_prior_a=big_a,
        sf_prior_b=big_a * s_f2_true,
        intercept_prior_sd=intercept_prior_sd,
    )
    st = eqx.tree_at(
        lambda s: (s.s_f2, s.sigma2),
        st,
        (
            jnp.float32(s_f2_true),
            jnp.float32(sigma2_true),
        ),
    )

    @enable_x64()
    def run_scan(st_init, keys_in):
        def body(carry, k_step):
            st_c = rank_k_mode_single_step(k_step, carry, spatial_cache)
            return st_c, (st_c.intercept, st_c.theta_f, st_c.f_loadings, st_c.nu_fit[0])

        return jax.lax.scan(body, st_init, keys_in)

    keys = random.split(random.key(2026), 2500)
    _, (int_tr, th_tr, f_tr, tau_tr) = run_scan(st, keys)
    burn = 500
    int_draws = np.asarray(int_tr[burn:], dtype=np.float64)
    th_draws = np.asarray(th_tr[burn:], dtype=np.float64)
    f_draws = np.asarray(f_tr[burn:], dtype=np.float64)
    tau_draws = np.asarray(tau_tr[burn:], dtype=np.float64)

    beta0_draws = np.column_stack([int_draws, th_draws, f_draws])
    np.testing.assert_allclose(beta0_draws.mean(axis=0), exact["beta0_mean"], atol=0.04)
    np.testing.assert_allclose(np.cov(beta0_draws, rowvar=False), exact["beta0_cov"], atol=0.03)

    tau_mcse = float(tau_draws.std(ddof=1) / np.sqrt(len(tau_draws) / 4.0))
    assert abs(float(tau_draws.mean()) - exact["tau_mean"]) < max(3.5 * tau_mcse, 0.025)
    assert abs(float(tau_draws.std(ddof=1)) - exact["tau_sd"]) / exact["tau_sd"] < 0.06


def test_geweke_joint_distribution_conjugate_and_spatial_blocks() -> None:
    """Section 5.1 Test 1: Geweke (2004) joint-distribution test for rank_k_mode_single_step."""
    rng = np.random.default_rng(99)
    N = 10
    n0 = 5
    q = 2
    ctrl_idx = np.arange(n0, dtype=np.int32)
    trt_idx = np.arange(n0, N, dtype=np.int32)
    z_vec = np.zeros(N, dtype=np.float64)
    z_vec[trt_idx] = 1.0
    w_unit = np.ones(N, dtype=np.float64)
    df_k = rng.normal(size=N)
    lam_mat = rng.normal(size=(N, q))
    coords = np.column_stack(
        [35.0 + 3.0 * rng.uniform(size=N), -95.0 + 3.0 * rng.uniform(size=N)]
    )

    # Note: for exact joint stationarity on all N units in Geweke's data-augmentation loop,
    # baseline_untreated_only=False ensures sigma^2 and beta_0 use the full likelihood p(y | Theta),
    # while use_spatial_gp=False tests the exact conjugate (intercept, theta_k, f_k, s_f^2, tau_k, sigma^2) transition.
    spatial_cache = build_spatial_kernel_cache(
        ctrl_idx=ctrl_idx,
        trt_idx=trt_idx,
        w_unit=w_unit,
        coords=coords,
        use_spatial_gp=False,
    )
    X_binned, max_split_mu, max_split_nu = _bin_mode_predictors(
        x_mod=rng.normal(size=(N, 2)),
        ps_vec=np.full(N, 0.45),
        df_std=df_k,
        has_forecast=True,
        key=random.key(11),
    )

    intercept_prior_sd = 0.5
    theta_prior_mean = 0.8
    theta_prior_sd = 0.4
    sigma_prior_a = 3.0
    sigma_prior_b = 1.5
    sf_prior_a = 3.0
    sf_prior_b = 0.8
    tau_prior_sd = 0.5

    base_state = init_rank_k_mode(
        X_binned=X_binned,
        max_split_mu=max_split_mu,
        max_split_nu=max_split_nu,
        dy_k=jnp.zeros(N, dtype=jnp.float32),
        df_k=jnp.asarray(df_k, dtype=jnp.float32),
        z_vec=jnp.asarray(z_vec, dtype=jnp.float32),
        w_unit=jnp.asarray(w_unit, dtype=jnp.float32),
        lambda_mat=jnp.asarray(lam_mat, dtype=jnp.float32),
        ctrl_idx=ctrl_idx,
        trt_idx=trt_idx,
        num_trees_mu=10,
        num_trees_nu=10,
        sigma_nu_total=tau_prior_sd,
        has_forecast=True,
        use_forests=False,
        use_factor_loadings=True,
        use_spatial_gp=False,
        baseline_untreated_only=False,
        theta_f_clip=(0.0, 1.5),
        theta_prior_mean=theta_prior_mean,
        theta_prior_sd=theta_prior_sd,
        sigma_prior_a=sigma_prior_a,
        sigma_prior_b=sigma_prior_b,
        sf_prior_a=sf_prior_a,
        sf_prior_b=sf_prior_b,
        intercept_prior_sd=intercept_prior_sd,
    )

    @enable_x64()
    def sample_prior_and_data(k_rng: jax.Array):
        k1, k2, k3, k4, k5, k6, k7 = random.split(k_rng, 7)
        int_draw = intercept_prior_sd * random.normal(k1, (), dtype=jnp.float32).astype(jnp.float64)
        th_draw = _sample_truncated_normal_64(k2, theta_prior_mean, theta_prior_sd, 0.0, 1.5)
        sf2_draw = _sample_inv_gamma_64(k3, sf_prior_a, sf_prior_b)
        f_draw = jnp.sqrt(sf2_draw) * random.normal(k4, (q,), dtype=jnp.float32).astype(jnp.float64)
        sig2_draw = _sample_inv_gamma_64(k5, sigma_prior_a, sigma_prior_b)
        tau_draw = tau_prior_sd * random.normal(k6, (), dtype=jnp.float32).astype(jnp.float64)
        mean_y = (
            int_draw
            + th_draw * jnp.asarray(df_k, dtype=jnp.float64)
            + jnp.asarray(lam_mat, dtype=jnp.float64) @ f_draw
            + tau_draw * jnp.asarray(z_vec, dtype=jnp.float64)
        )
        y_sim = mean_y + jnp.sqrt(sig2_draw) * random.normal(k7, (N,), dtype=jnp.float32).astype(jnp.float64)
        return int_draw, th_draw, sf2_draw, f_draw, sig2_draw, tau_draw, y_sim

    B = 350
    prior_keys = random.split(random.key(501), B)
    mc_int, mc_th, mc_sf2, _, mc_sig2, mc_tau, _ = jax.vmap(sample_prior_and_data)(prior_keys)

    # Successive-conditional Geweke chain
    int0, th0, sf20, f0, sig20, tau0, y0 = sample_prior_and_data(random.key(701))
    with enable_x64():
        R0 = (
            y0
            - int0
            - th0 * jnp.asarray(df_k, dtype=jnp.float64)
            - jnp.asarray(lam_mat, dtype=jnp.float64) @ f0
            - tau0 * jnp.asarray(z_vec, dtype=jnp.float64)
        ).astype(jnp.float32)
    st_curr = eqx.tree_at(
        lambda s: (s.y, s.resid, s.intercept, s.theta_f, s.s_f2, s.f_loadings, s.sigma2, s.nu_fit),
        base_state,
        (
            y0.astype(jnp.float32),
            R0,
            int0.astype(jnp.float32),
            th0.astype(jnp.float32),
            sf20.astype(jnp.float32),
            f0.astype(jnp.float32),
            sig20.astype(jnp.float32),
            jnp.full(N, tau0.astype(jnp.float32), dtype=jnp.float32),
        ),
    )

    @enable_x64()
    def run_geweke_scan(st_init, keys_in):
        df64 = jnp.asarray(df_k, dtype=jnp.float64)
        lam64 = jnp.asarray(lam_mat, dtype=jnp.float64)
        z64 = jnp.asarray(z_vec, dtype=jnp.float64)

        def outer_step(st_c, k_outer):
            sub_keys = random.split(k_outer, 6)

            def inner_step(i, st_i):
                k_y, k_g = random.split(sub_keys[i])
                mean_curr = (
                    jnp.asarray(st_i.intercept, dtype=jnp.float64)
                    + jnp.asarray(st_i.theta_f, dtype=jnp.float64) * df64
                    + lam64 @ jnp.asarray(st_i.f_loadings, dtype=jnp.float64)
                    + jnp.asarray(st_i.nu_fit, dtype=jnp.float64) * z64
                )
                eps_new = jnp.sqrt(jnp.asarray(st_i.sigma2, dtype=jnp.float64)) * random.normal(
                    k_y, (N,), dtype=jnp.float32
                ).astype(jnp.float64)
                y_new = (mean_curr + eps_new).astype(jnp.float32)
                R_new = eps_new.astype(jnp.float32)
                st_i = eqx.tree_at(lambda s: (s.y, s.resid), st_i, (y_new, R_new))
                return rank_k_mode_single_step(k_g, st_i, spatial_cache)

            st_next = jax.lax.fori_loop(0, 6, inner_step, st_c)
            return st_next, (
                st_next.intercept,
                st_next.theta_f,
                st_next.s_f2,
                st_next.sigma2,
                st_next.nu_fit[0],
            )

        return jax.lax.scan(outer_step, st_init, keys_in)

    scan_keys = random.split(random.key(909), B)
    _, (sc_int, sc_th, sc_sf2, sc_sig2, sc_tau) = run_geweke_scan(st_curr, scan_keys)

    for a, b in [
        (np.asarray(mc_int), np.asarray(sc_int)),
        (np.asarray(mc_th), np.asarray(sc_th)),
        (np.asarray(mc_sf2), np.asarray(sc_sf2)),
        (np.asarray(mc_sig2), np.asarray(sc_sig2)),
        (np.asarray(mc_tau), np.asarray(sc_tau)),
    ]:
        stat = ks_2samp(a, b)
        assert stat.pvalue > 0.002


def test_m1_reduction_and_variance_decomposition_identity() -> None:
    """Section 5.1 Tests 3 & 4: M=1 reduction and exact variance decomposition identity."""
    y, z, x, yh, coords, cohorts = _synthetic_staggered_panel(seed=202, N=60, T=9, T_hist=20)
    feats = build_causal_features(y, z, x=x, y_pre_history=yh, seed=11)

    # M = 1: a single forecast draw propagates no forecast uncertainty, so the JAX
    # sampler conditions on the point forecast and retains num_mcmc draws in total.
    # (engine="jax" with uncertainty_propagation=False is rejected by the config, and
    # uncertainty_propagation=False selects the stochtree plug-in variant.)
    cfg_m1 = RankKLongBetConfig(
        rank_k=2,
        uncertainty_propagation=True,
        n_forecast_draws=1,
        max_forecast_draws=1,
        num_chains=2,
        num_burnin=25,
        num_mcmc=16,
        num_trees_mu=25,
        num_trees_tau_mode1=15,
        num_trees_tau_higher=10,
        random_seed=77,
    )
    with pytest.warns(RankKConfigWarning, match="contradicts n_forecast_draws=1"):
        res_m1 = fit_rank_k_longbet(y, z, x, features=feats, coords=coords, cohorts=cohorts, cfg=cfg_m1)

    assert res_m1.n_components == 1
    assert res_m1.att_var_between == 0.0
    assert res_m1.att_draws.shape == (16,)  # num_chains x round(num_mcmc / num_chains)
    assert abs((res_m1.att_var_within + res_m1.att_var_between) - res_m1.att_se ** 2) < 1e-10

    # M = 1 reduction: n_forecast_draws=0 propagates no forecast draws either and gives the same posterior.
    cfg_m0 = dataclasses.replace(cfg_m1, n_forecast_draws=0, max_forecast_draws=None)
    with pytest.warns(RankKConfigWarning, match="contradicts n_forecast_draws=0"):
        res_m0 = fit_rank_k_longbet(y, z, x, features=feats, coords=coords, cohorts=cohorts, cfg=cfg_m0)
    assert res_m0.n_components == 1 and res_m0.att_var_between == 0.0
    np.testing.assert_allclose(res_m0.att_draws, res_m1.att_draws, atol=1e-12)

    # M = 6 with uncertainty propagation enabled
    cfg_m6 = RankKLongBetConfig(
        rank_k=2,
        uncertainty_propagation=True,
        n_forecast_draws=6,
        max_forecast_draws=6,
        num_chains=2,
        num_burnin=25,
        draws_per_component=8,
        num_trees_mu=25,
        num_trees_tau_mode1=15,
        num_trees_tau_higher=10,
        random_seed=77,
    )
    res_m6 = fit_rank_k_longbet(y, z, x, features=feats, coords=coords, cohorts=cohorts, cfg=cfg_m6)
    assert res_m6.n_components == 6
    assert res_m6.att_var_between > 0.0
    assert abs((res_m6.att_var_within + res_m6.att_var_between) - res_m6.att_se ** 2) < 1e-10


def test_proposition1_post_anchor_invariance_v2() -> None:
    """Section 5.1 Test 5: Proposition 1 post-anchor perturbation invariance of Module 1 & Rank-K v2."""
    y, z, x, yh, coords, cohorts = _synthetic_staggered_panel(seed=303, N=60, T=9, T_hist=20)
    feats_1 = build_causal_features(
        y, z, x=x, y_pre_history=yh, seed=19, n_forecast_draws=6, draw_method="posterior"
    )
    y_pert = y.copy()
    ever = cohorts < y.shape[1]
    y_pert[ever, feats_1.anchor_t :] += 250.0
    feats_2 = build_causal_features(
        y_pert, z, x=x, y_pre_history=yh, seed=19, n_forecast_draws=6, draw_method="posterior"
    )

    np.testing.assert_array_equal(feats_1.point_forecast, feats_2.point_forecast)
    np.testing.assert_array_equal(feats_1.point_forecast_draws, feats_2.point_forecast_draws)
    np.testing.assert_array_equal(feats_1.pre_anchor_loadings, feats_2.pre_anchor_loadings)

    cfg = RankKLongBetConfig(
        rank_k=2,
        uncertainty_propagation=True,
        n_forecast_draws=6,
        max_forecast_draws=6,
        num_chains=2,
        num_burnin=25,
        draws_per_component=8,
        num_trees_mu=25,
        num_trees_tau_mode1=15,
        num_trees_tau_higher=10,
        random_seed=88,
    )
    res_1 = fit_rank_k_longbet(y, z, x, features=feats_1, coords=coords, cohorts=cohorts, cfg=cfg)
    res_2 = fit_rank_k_longbet(y, z, x, features=feats_2, coords=coords, cohorts=cohorts, cfg=cfg)
    np.testing.assert_array_equal(res_1.att_draws, res_2.att_draws)
    np.testing.assert_array_equal(res_1.catt_draws, res_2.catt_draws)

    # Under baseline_cells="untreated", baseline parameters (theta_f, sigma, spatial GP)
    # are also bit-for-bit invariant to post-anchor perturbations of y[ever, anchor_t:]
    res_pert = fit_rank_k_longbet(y_pert, z, x, features=feats_2, coords=coords, cohorts=cohorts, cfg=cfg)
    for k in range(2):
        np.testing.assert_array_equal(res_1.theta_f_draws[k], res_pert.theta_f_draws[k])
        np.testing.assert_array_equal(res_1.gp_hyper_draws[k]["s_zeta"], res_pert.gp_hyper_draws[k]["s_zeta"])
        np.testing.assert_array_equal(res_1.gp_hyper_draws[k]["ell"], res_pert.gp_hyper_draws[k]["ell"])
    assert res_1.sigma_modes == res_pert.sigma_modes


def test_proposition2_k_invariance_v2() -> None:
    """Section 5.1 Test 6: Proposition 2 exact K-invariance of ATT draws and unit-averaged CATEs in v2."""
    y, z, x, yh, coords, cohorts = _synthetic_staggered_panel(seed=404, N=60, T=9, T_hist=20)
    feats = build_causal_features(y, z, x=x, y_pre_history=yh, seed=23)

    results = {}
    for K in (1, 2, 3):
        cfg = RankKLongBetConfig(
            rank_k=K,
            uncertainty_propagation=True,
            n_forecast_draws=6,
            max_forecast_draws=6,
            num_chains=2,
            num_burnin=25,
            draws_per_component=8,
            num_trees_mu=25,
            num_trees_tau_mode1=15,
            num_trees_tau_higher=10,
            random_seed=91,
        )
        results[K] = fit_rank_k_longbet(
            y, z, x, features=feats, coords=coords, cohorts=cohorts, cfg=cfg
        )

    np.testing.assert_allclose(results[1].att_draws, results[2].att_draws, atol=1e-12)
    np.testing.assert_allclose(results[1].att_draws, results[3].att_draws, atol=1e-12)

    ever = cohorts < y.shape[1]
    for i in np.where(ever)[0]:
        trt_i = z[i] == 1
        u1 = float(results[1].catt_hat[i, trt_i].mean())
        u2 = float(results[2].catt_hat[i, trt_i].mean())
        u3 = float(results[3].catt_hat[i, trt_i].mean())
        assert abs(u1 - u2) < 1e-12
        assert abs(u1 - u3) < 1e-12


def test_stage_e_jax_vs_stochtree_equivalence_and_placebo_audit() -> None:
    """Stage E Equivalence Gate: JAX sampler without forecast draws vs the stochtree plug-in variant."""
    pytest.importorskip("stochtree")
    y, z, x, yh, coords, cohorts = _synthetic_staggered_panel(seed=505, N=100, T=10, T_hist=24)
    feats = build_causal_features(y, z, x=x, y_pre_history=yh, seed=31)

    # A single forecast draw makes the JAX sampler condition on the point forecast (the
    # code path of the former engine="jax", uncertainty_propagation=False); the config
    # warns that this contradicts uncertainty_propagation=True.
    cfg_jax = RankKLongBetConfig(
        rank_k=2,
        engine="jax",
        uncertainty_propagation=True,
        n_forecast_draws=1,
        max_forecast_draws=1,
        use_spatial_gp=False,
        use_factor_loadings=False,
        baseline_cells="untreated",
        num_chains=4,
        num_burnin=150,
        num_mcmc=400,
        num_trees_mu=75,
        num_trees_tau_mode1=35,
        num_trees_tau_higher=20,
        random_seed=105,
    )
    cfg_st = RankKLongBetConfig(
        rank_k=2,
        engine="stochtree",
        num_gfr=10,
        num_burnin=150,
        num_mcmc=400,
        num_trees_mu=75,
        num_trees_tau_mode1=35,
        num_trees_tau_higher=20,
        random_seed=105,
    )

    with pytest.warns(RankKConfigWarning, match="contradicts n_forecast_draws=1"):
        res_jax = fit_rank_k_longbet(y, z, x, features=feats, cohorts=cohorts, cfg=cfg_jax)
    res_st = fit_rank_k_longbet(y, z, x, features=feats, cohorts=cohorts, cfg=cfg_st)

    # ATT point estimates agree closely
    assert abs(res_jax.att_hat - res_st.att_hat) < 0.15
    # Unit-averaged CATE correlation across treated units is high
    ever = cohorts < y.shape[1]
    u_jax = np.array([res_jax.catt_hat[i, z[i] == 1].mean() for i in np.where(ever)[0]])
    u_st = np.array([res_st.catt_hat[i, z[i] == 1].mean() for i in np.where(ever)[0]])
    corr = float(np.corrcoef(u_jax, u_st)[0, 1])
    assert corr > 0.85

    # Placebo audit covers zero on clean pre-adoption history (reuses N=60 panel so XLA cache hits)
    y60, z60, x60, yh60, coords60, _ = _synthetic_staggered_panel(seed=202, N=60, T=9, T_hist=20)
    pl = rank_k_placebo_audit(
        y60,
        z60,
        x60,
        y_pre_history=yh60,
        coords=coords60,
        cfg=RankKLongBetConfig(
            rank_k=2,
            uncertainty_propagation=True,
            n_forecast_draws=6,
            max_forecast_draws=6,
            num_chains=2,
            num_burnin=25,
            draws_per_component=8,
            num_trees_mu=25,
            num_trees_tau_mode1=15,
            num_trees_tau_higher=10,
            random_seed=105,
        ),
    )
    assert pl["covers_zero"] is True
    assert pl["flag_parallel_trends_violation"] is False


def test_informative_tau_priors_and_catt_floor() -> None:
    """Verify native informative Mode-1 priors (tau_prior_mean, tau_prior_sd, tau_clip) and catt_floor."""
    with pytest.raises(ValueError, match="tau_prior_sd must be > 0"):
        RankKLongBetConfig(tau_prior_sd=-0.5)
    with pytest.raises(ValueError, match="tau_clip must be a"):
        RankKLongBetConfig(tau_clip=(2.0, -1.0))

    rng = np.random.default_rng(777)
    N, T = 60, 7
    x = rng.normal(size=(N, 2))
    z = np.zeros((N, T), dtype=np.int32)
    z[:30, 1:] = 1
    # Noisy panel with a slight negative sample difference in means
    y = rng.normal(scale=2.0, size=(N, T))
    y[:30, 1:] -= 0.4

    cfg_flat = RankKLongBetConfig(
        rank_k=2,
        use_forecast=False,
        use_spatial_gp=False,
        use_factor_loadings=False,
        num_chains=2,
        num_burnin=30,
        num_mcmc=120,
        num_trees_mu=20,
        num_trees_tau_mode1=15,
        num_trees_tau_higher=10,
        random_seed=42,
    )
    cfg_nonneg = dataclasses.replace(
        cfg_flat,
        tau_prior_mean=0.0,
        tau_prior_sd=0.5,
        tau_clip=(0.0, float("inf")),
    )
    cfg_nonneg_k1 = dataclasses.replace(cfg_nonneg, rank_k=1)
    cfg_floor = dataclasses.replace(
        cfg_nonneg,
        catt_floor=0.0,
    )
    cfg_floor_k1 = dataclasses.replace(cfg_floor, rank_k=1)

    res_flat = fit_rank_k_longbet(y, z, x, cfg=cfg_flat)
    res_nonneg = fit_rank_k_longbet(y, z, x, cfg=cfg_nonneg)
    res_nonneg_k1 = fit_rank_k_longbet(y, z, x, cfg=cfg_nonneg_k1)
    res_floor = fit_rank_k_longbet(y, z, x, cfg=cfg_floor)
    res_floor_k1 = fit_rank_k_longbet(y, z, x, cfg=cfg_floor_k1)

    # 1. Non-negative Mode-1 prior guarantees all ATT draws >= 0 and reduces posterior SD
    assert np.all(res_nonneg.att_draws >= 0.0)
    assert res_nonneg.att_ci_low >= 0.0
    assert res_nonneg.att_se < res_flat.att_se

    # 2. Exact K-invariance holds under informative Mode-1 priors AND under catt_floor=0.0
    np.testing.assert_allclose(res_nonneg.att_draws, res_nonneg_k1.att_draws, atol=1e-12)
    np.testing.assert_allclose(res_floor.att_draws, res_floor_k1.att_draws, atol=1e-12)
    np.testing.assert_allclose(res_floor.mode_tau_draws[0], res_floor_k1.mode_tau_draws[0], atol=1e-12)
    np.testing.assert_allclose(res_floor.att_draws, res_nonneg.att_draws, atol=1e-12)
    np.testing.assert_allclose(res_floor.att_var_total, res_floor.att_se ** 2, rtol=1e-10, atol=1e-12)

    # 3. catt_floor=0.0 guarantees every post-adoption cell effect draw is >= 0
    assert np.all(res_floor.catt_draws[:, z == 1] >= 0.0)
    assert np.all(res_floor.catt_ci_low[z == 1] >= 0.0)


