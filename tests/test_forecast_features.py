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

"""Verification tests for the conjugate baseline_forecast channel and build_causal_features."""

from __future__ import annotations

from unittest.mock import patch

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from scipy import linalg

from longbet import (
    CausalFeatures,
    LongBet,
    LongBetConfig,
    build_causal_features,
)
from longbet._trend_basis import (
    apply_baseline_forecast,
    build_period_basis,
    build_time_basis,
    build_trend_design,
    build_unit_features,
    evaluate_trend,
    prepare_baseline_forecast,
    trend_block_step,
    trend_prior_scales,
    trend_sufficient_statistics,
)
from longbet._x64 import enable_x64


def _make_staggered_panel(seed: int = 42, N: int = 40, T: int = 10, T_hist: int = 24):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, 3)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.int32)
    # Staggered adoption starting at t = 4 (so anchor_t = 4)
    z[:10, 4:] = 1
    z[10:20, 6:] = 1
    # Units 20..39 are never-treated controls
    phase = rng.uniform(0.0, 2.0 * np.pi, size=N)
    slope = 0.25 * x[:, 0] + 0.15 * np.sin(x[:, 1])
    t_all = np.arange(-T_hist, T, dtype=np.float64)
    y_all = (
        1.5 * x[:, 0:1]
        + slope[:, None] * (t_all[None, :] / T)
        + 0.8 * np.cos(2.0 * np.pi * t_all[None, :] / 6.0 + phase[:, None])
        + 0.2 * rng.normal(size=(N, T_hist + T))
    )
    y_pre_history = y_all[:, :T_hist].astype(np.float32)
    y = y_all[:, T_hist:].copy().astype(np.float32)
    # Add true treatment effect on treated cells
    y[z == 1] += 2.0
    return y, z, x, y_pre_history


def test_prepare_baseline_forecast_orthogonality():
    """Prepared forecast columns must be orthogonal to unit intercepts and period basis."""
    rng = np.random.default_rng(101)
    N, T = 30, 12
    # Include arbitrary unit means, period shocks, and interactions
    unit_eff = rng.normal(scale=3.0, size=(N, 1, 1))
    period_eff = rng.normal(scale=2.0, size=(1, T, 1))
    raw_f = unit_eff + period_eff + rng.normal(size=(N, T, 2))

    f_cols, spec = prepare_baseline_forecast(raw_f, N, T)
    assert f_cols.shape == (N * T, 2)
    assert spec.num_cols == 2

    f_grid = f_cols.reshape(N, T, 2)
    # 1. Zero sum within every unit (orthogonal to gamma_i and unit-constant mu(X_i))
    np.testing.assert_allclose(f_grid.sum(axis=1), 0.0, atol=1e-5)
    # 2. Zero sum within every period (orthogonal to constant & period effects)
    np.testing.assert_allclose(f_grid.sum(axis=0), 0.0, atol=1e-5)
    # 3. Exact orthogonality to the complete (T - 1) DCT-II period basis
    psi = build_period_basis(T)  # (T, T - 1)
    for k in range(2):
        cross = (f_grid[:, :, k] @ psi).sum(axis=0)
        np.testing.assert_allclose(cross, 0.0, atol=1e-4)
    # 4. Unit RMS per column
    rms = np.sqrt(np.mean(f_cols.astype(np.float64) ** 2, axis=0))
    np.testing.assert_allclose(rms, 1.0, atol=1e-6)

    # 5. Round-trip through apply_baseline_forecast reproduces f_cols
    f_applied = apply_baseline_forecast(raw_f, spec, N, T)
    np.testing.assert_allclose(f_applied, f_cols, atol=1e-6)


def test_causal_features_post_treatment_invariance():
    """Perturbing treated outcomes y[z == 1] by +1000 must change nothing (zero leakage)."""
    y, z, x, y_hist = _make_staggered_panel(seed=77)
    feats_clean = build_causal_features(
        y, z, x=x, y_pre_history=y_hist, seed=11
    )
    assert isinstance(feats_clean, CausalFeatures)
    assert feats_clean.anchor_t == 4
    assert feats_clean.baseline_forecast.shape == y.shape
    assert feats_clean.unit_features.shape == (y.shape[0], 2)
    assert feats_clean.ps.shape == (y.shape[0],)

    # Corrupt all treated cells AND all post-anchor cells of ever-treated units
    y_corrupted = y.copy()
    y_corrupted[z == 1] += 1000.0
    ever_treated = np.any(z == 1, axis=1)
    y_corrupted[ever_treated, feats_clean.anchor_t :] -= 500.0

    feats_corrupted = build_causal_features(
        y_corrupted, z, x=x, y_pre_history=y_hist, seed=11
    )
    np.testing.assert_array_equal(
        feats_clean.baseline_forecast, feats_corrupted.baseline_forecast
    )
    np.testing.assert_array_equal(
        feats_clean.unit_features, feats_corrupted.unit_features
    )
    np.testing.assert_array_equal(feats_clean.ps, feats_corrupted.ps)

    # Verify that hazardous options (anchor_t, use_control_adapter, raw_zero_shot_forecast)
    # are not exposed on the public API
    import inspect

    sig = inspect.signature(build_causal_features)
    assert "anchor_t" not in sig.parameters
    assert "use_control_adapter" not in sig.parameters
    assert "forecast_draw_calibration" not in sig.parameters
    assert not hasattr(feats_clean, "raw_zero_shot_forecast")


def test_custom_forecaster_callable_protocol():
    """A callable fn(context, horizon) -> (N, horizon) is the only external-forecaster protocol."""
    y, z, x, y_hist = _make_staggered_panel(seed=19)
    seen: dict[str, int] = {}

    def my_forecaster(context: np.ndarray, horizon: int) -> np.ndarray:
        seen["len"] = int(context.shape[1])
        t_h = np.arange(horizon, dtype=np.float64)[None, :]
        return np.mean(context[:, -4:], axis=1, keepdims=True) + 0.8 * np.cos(2.0 * np.pi * t_h / 6.0) * np.std(
            context, axis=1, keepdims=True
        )

    feats = build_causal_features(y, z, x=x, y_pre_history=y_hist, forecaster=my_forecaster, seed=3)
    # Context length must be T_hist (24) + anchor_t (4) = 28: strictly pre-anchor
    assert seen["len"] == 24 + 4
    assert feats.baseline_forecast is not None
    assert np.all(np.isfinite(feats.baseline_forecast))
    # TimesFM-style model objects are not accepted any more: only callables and the two built-in names
    class NotCallable:
        def predict_batch(self, contexts, horizon, return_quantiles=False):  # pragma: no cover
            raise AssertionError("must not be called")

    with pytest.raises(TypeError):
        build_causal_features(y, z, x=x, y_pre_history=y_hist, forecaster=NotCallable(), seed=3)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        build_causal_features(y, z, x=x, y_pre_history=y_hist, forecaster=lambda c, h: c[:, :1], seed=3)


def test_collapse_gamma_mu_preserved_with_baseline_forecast_and_ps():
    """Conjugate baseline_forecast + unit-level ps preserves collapse_gamma_mu=True."""
    y, z, x, y_hist = _make_staggered_panel(seed=2026, N=32, T=8, T_hist=16)
    feats = build_causal_features(y, z, x=x, y_pre_history=y_hist, seed=5)
    x_aug = np.hstack([x, feats.unit_features])

    cfg = LongBetConfig(num_chains=2, num_burnin=4, num_sweeps=4, random_seed=0)
    model_conj = LongBet(cfg).fit(
        y=y,
        x=x_aug,
        z=z,
        ps=feats.ps,
        baseline_forecast=feats.baseline_forecast,
    )
    assert model_conj.state is not None
    assert model_conj.state.collapse_gamma_mu is True
    assert model_conj.state.collapse_gamma_nu is True
    # trend_num_period should be (T - 1) + 1 = 8
    assert model_conj.state.trend_num_period == (y.shape[1] - 1) + 1

    # Contrast with naive x_tv injection (Trap 1): disables collapse_gamma_mu
    model_xtv = LongBet(cfg).fit(
        y=y,
        x=x,
        z=z,
        x_tv=feats.baseline_forecast[:, :, None],
    )
    assert model_xtv.state is not None
    assert model_xtv.state.collapse_gamma_mu is False


def test_joint_conditional_matches_dense_oracle_with_forecast_cols():
    """Verify trend_block_step matches dense float64 oracle when forecast_cols are active."""
    rng = np.random.default_rng(909)
    N, T, S = 8, 6, 4
    n = N * T
    ui = np.repeat(np.arange(N), T).astype(np.int32)
    ti = np.tile(np.arange(T), N).astype(np.int32)

    x_unit = rng.normal(size=(N, 3)).astype(np.float32)
    features, _ = build_unit_features(x_unit, num_random_features=4)
    phi = build_time_basis(T, degree=2)
    period = build_period_basis(T)
    raw_forecast = rng.normal(size=(N, T, 2))
    f_cols, _ = prepare_baseline_forecast(raw_forecast, N, T)

    H = np.asarray(
        build_trend_design(
            features, phi, ui, ti, period=period, forecast_cols=f_cols
        ),
        dtype=np.float32,
    )
    q = H.shape[1]
    n_period = period.shape[1]
    n_forecast = f_cols.shape[1]
    n_inter = q - n_period - n_forecast
    sc = trend_prior_scales(
        n_period,
        n_inter,
        prior_scale=2.0,
        period_scale=2.0,
        num_forecast=n_forecast,
        forecast_scale=1.2,
    )

    exposure = np.tile([0, 0, 0, 1, 2, 3], N).astype(np.int32)
    exposure[ui >= 5] = 0
    d = ((exposure > 0) * rng.uniform(0.5, 1.3, n)).astype(np.float32)
    mask = np.ones(n, dtype=bool)
    mask[[2, 17]] = False

    A = rng.normal(size=(S, S))
    L = np.linalg.cholesky(A @ A.T + np.eye(S)).astype(np.float32)
    y = rng.normal(size=n).astype(np.float32)
    c = rng.normal(size=q).astype(np.float32)
    b = rng.normal(size=S).astype(np.float32)
    g = rng.normal(size=N).astype(np.float32)
    R = (y - H @ c - d * b[exposure] - g[ui]).astype(np.float32)
    sigma2, gs2, temp = map(np.float32, (0.15, 0.25, 1.0))
    precision = np.where(mask, float(temp) / float(sigma2), 0.0)

    restored = R.astype(float) + H.astype(float) @ c + g[ui] + d.astype(float) * b[exposure]
    blocks = [
        H.astype(float) * sc,
        d.astype(float)[:, None] * L[exposure].astype(float),
        np.eye(N)[ui] * np.sqrt(float(gs2)),
    ]
    X = np.column_stack(blocks)
    P = np.eye(X.shape[1]) + X.T @ (precision[:, None] * X)
    chol = linalg.cholesky(P, lower=True)
    covariance = linalg.cho_solve((chol, True), np.eye(len(P)))
    mean = linalg.cho_solve((chol, True), X.T @ (precision * restored))

    t_gram, t_unit_sum = trend_sufficient_statistics(
        jnp.asarray(np.where(mask[:, None], H, 0.0)),
        jnp.asarray(mask),
        jnp.asarray(ui),
        N,
    )
    args = dict(
        trend_coef=c,
        gamma=g,
        beta=b,
        R=R,
        H=np.where(mask[:, None], H, 0.0).astype(np.float32),
        d_vec=d,
        obs_mask=mask,
        unit_idx=ui,
        N_units=N,
        exposure_idx=exposure,
        sigma2=sigma2,
        sigma_gamma2=gs2,
        K_chol=L,
        trend_col_scale=sc,
        random_intercept=True,
        sample_beta=True,
        temperature=temp,
        conditional_precision=None,
        trend_gram=t_gram,
        trend_unit_sum=t_unit_sum,
    )
    k = q + S

    @enable_x64()
    def draw(noise):
        with patch("jax.random.normal", side_effect=[noise[:k], noise[k:]]):
            cn, gn, bn, _ = trend_block_step(jax.random.key(0), **args)
        return jnp.concatenate(
            [
                cn.astype(jnp.float64) / sc,
                jax.scipy.linalg.solve_triangular(
                    jnp.asarray(L, dtype=jnp.float64),
                    bn.astype(jnp.float64),
                    lower=True,
                ),
                gn.astype(jnp.float64) / np.sqrt(float(gs2)),
            ]
        )

    with enable_x64():
        noise = jnp.zeros(len(P), dtype=jnp.float64)
        actual_mean = np.asarray(draw(noise))
        factor = np.asarray(jax.jacfwd(draw)(noise), dtype=float)
    actual_covariance = factor @ factor.T
    relative_covariance = chol.T @ actual_covariance @ chol
    np.testing.assert_allclose(relative_covariance, np.eye(len(P)), atol=3e-4, rtol=0)
    np.testing.assert_allclose(
        (actual_mean - mean) / np.sqrt(np.diag(covariance)), 0.0, atol=3e-4, rtol=0
    )


def test_save_load_and_predict_roundtrip_with_baseline_forecast(tmp_path):
    """Save/load archive must preserve baseline_forecast spec and enforce predict contract."""
    y, z, x, y_hist = _make_staggered_panel(seed=314, N=24, T=8, T_hist=12)
    feats = build_causal_features(y, z, x=x, y_pre_history=y_hist, seed=7)

    cfg = LongBetConfig(num_chains=2, num_burnin=4, num_sweeps=4, random_seed=9)
    model = LongBet(cfg).fit(
        y=y,
        x=x,
        z=z,
        ps=feats.ps,
        baseline_forecast=feats.baseline_forecast,
    )
    pred1 = model.predict(
        x=x,
        z=z,
        ps=feats.ps,
        baseline_forecast=feats.baseline_forecast,
    )

    # Omitting baseline_forecast at predict time must raise ValueError
    with pytest.raises(ValueError, match="model was fitted with baseline_forecast"):
        model.predict(x=x, z=z, ps=feats.ps)

    archive = tmp_path / "model_with_forecast.npz"
    model.save(archive)
    loaded = LongBet.load(archive)
    assert loaded.trend_forecast_spec_ is not None

    with pytest.raises(ValueError, match="model was fitted with baseline_forecast"):
        loaded.predict(x=x, z=z, ps=feats.ps)

    pred2 = loaded.predict(
        x=x,
        z=z,
        ps=feats.ps,
        baseline_forecast=feats.baseline_forecast,
    )
    np.testing.assert_allclose(pred1.att_full, pred2.att_full, atol=1e-6)
    np.testing.assert_allclose(pred1.mu0_summary.mean, pred2.mu0_summary.mean, atol=1e-6)


def test_apply_baseline_forecast_batch_invariance_and_horizon_extension():
    """Subset predictions (N_new >= 1) must be batch-invariant and support T_new > T_fit."""
    rng = np.random.default_rng(404)
    N, T = 24, 8
    raw_f = rng.normal(size=(N, T, 2)).astype(np.float32)
    f_full, spec = prepare_baseline_forecast(raw_f, N, T)
    f_full_grid = f_full.reshape(N, T, 2)

    # Subset of 5 units and single unit must match full batch row-for-row
    subset_idx = np.array([0, 3, 7, 12, 19])
    f_sub = apply_baseline_forecast(raw_f[subset_idx], spec, len(subset_idx), T).reshape(
        len(subset_idx), T, 2
    )
    np.testing.assert_allclose(f_sub, f_full_grid[subset_idx], atol=1e-6)

    f_single = apply_baseline_forecast(raw_f[3:4], spec, 1, T).reshape(1, T, 2)
    np.testing.assert_allclose(f_single, f_full_grid[3:4], atol=1e-6)

    # Horizon extension T_new > T_fit must not raise IndexError and must preserve fit prefix
    T_new = T + 4
    raw_ext = np.concatenate(
        [raw_f, rng.normal(size=(N, 4, 2)).astype(np.float32)], axis=1
    )
    f_ext = apply_baseline_forecast(raw_ext, spec, N, T_new).reshape(N, T_new, 2)
    assert f_ext.shape == (N, T_new, 2)
    assert np.all(np.isfinite(f_ext))
    np.testing.assert_allclose(f_ext[:, :T, :], f_full_grid, atol=1e-6)


def test_save_load_with_trend_degree_zero_and_no_period_effects(tmp_path):
    """Archive round-trip must preserve 2D (T, 0) shapes when trend_degree=0."""
    y, z, x, y_hist = _make_staggered_panel(seed=515, N=28, T=8, T_hist=12)
    feats = build_causal_features(y, z, x=x, y_pre_history=y_hist, seed=13)
    cfg = LongBetConfig(
        num_chains=2,
        num_burnin=3,
        num_sweeps=3,
        trend_degree=0,
        trend_period_effects=False,
        random_seed=11,
    )
    model = LongBet(cfg).fit(
        y=y, x=x, z=z, ps=feats.ps, baseline_forecast=feats.baseline_forecast
    )
    archive = tmp_path / "model_deg0.npz"
    model.save(archive)
    loaded = LongBet.load(archive)
    assert loaded.trend_time_basis_.shape == (8, 0)
    assert loaded.trend_period_basis_.shape == (8, 0)
    pred = loaded.predict(
        x=x[:5],
        z=z[:5],
        ps=feats.ps[:5],
        baseline_forecast=feats.baseline_forecast[:5],
    )
    assert pred.mu0_summary.mean.shape == (5, 8)
    assert pred.muhats0.shape[:2] == (5, 8)


def test_cross_fit_propensity_pure_random_assignment_no_jackknife_inversion():
    """Under pure RCT assignment, cross-fitted propensity must not invert or overfit."""
    rng = np.random.default_rng(808)
    N, T, T_hist = 80, 10, 16
    x = rng.normal(size=(N, 4)).astype(np.float32)
    y_all = rng.normal(size=(N, T_hist + T)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.int32)
    # Pure coin-flip RCT independent of x and y_all
    treated = rng.uniform(size=N) < 0.4
    z[treated, 5:] = 1

    feats = build_causal_features(
        y_all[:, T_hist:], z, x=x, y_pre_history=y_all[:, :T_hist], seed=21
    )
    ever_treated = treated.astype(np.float64)
    # Must not have negative jackknife correlation or spurious dispersion
    assert float(np.std(feats.ps)) < 0.05
    if float(np.std(feats.ps)) > 1e-6:
        corr = float(np.corrcoef(ever_treated, feats.ps)[0, 1])
        assert corr >= -0.05


def test_beyond_poly_r2_and_recommended_forecast_gating():
    """beyond_poly_r2 must gate off pure polynomial panels and gate on seasonal panels."""
    rng = np.random.default_rng(909)
    N, T, T_hist = 40, 10, 24
    T_all = T_hist + T
    t_grid = np.arange(T_all, dtype=np.float64)[None, :]
    x = rng.normal(size=(N, 4)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.int32)
    z[:20, 4:] = 1

    # S1: Pure quadratic polynomial in time + noise (no non-polynomial harmonic dynamics)
    slope = 0.25 * x[:, 0:1]
    curve = 0.01 * x[:, 1:2]
    y_poly = (
        1.0 * x[:, 0:1]
        + slope * (t_grid - T_hist)
        + curve * ((t_grid - T_hist) ** 2)
        + rng.normal(0.0, 0.5, size=(N, T_all))
    ).astype(np.float32)
    feats_poly = build_causal_features(
        y_poly[:, T_hist:], z, x=x, y_pre_history=y_poly[:, :T_hist], seed=31
    )
    assert 0.0 <= feats_poly.beyond_poly_r2 < 0.05
    assert not feats_poly.gate_open
    assert feats_poly.baseline_forecast is None
    assert feats_poly.recommended_forecast is None
    assert feats_poly.recommended_ps is None
    assert feats_poly.augment_x(x) is x

    # S2: Strong multi-harmonic seasonality + local trend break (rich non-polynomial signal)
    y_seas, z_s, x_s, y_hist_s = _make_staggered_panel(seed=910, N=40, T=10, T_hist=24)
    feats_seas = build_causal_features(
        y_seas, z_s, x=x_s, y_pre_history=y_hist_s, seed=32
    )
    assert feats_seas.beyond_poly_r2 >= 0.20
    assert feats_seas.baseline_forecast is not None
    assert feats_seas.recommended_forecast is not None
    np.testing.assert_array_equal(feats_seas.recommended_forecast, feats_seas.baseline_forecast)

    # S3: Short context window (L_ctx = T_hist + anchor_t = 4 + 4 = 8 < 12) must gate off
    y_short, z_sh, x_sh, y_hist_sh = _make_staggered_panel(seed=911, N=40, T=10, T_hist=4)
    feats_short = build_causal_features(
        y_short, z_sh, x=x_sh, y_pre_history=y_hist_sh, seed=33
    )
    assert feats_short.baseline_forecast is None
    assert feats_short.recommended_forecast is None


def test_default_backbone_is_patch_ridge_and_reduces_to_spectral_on_small_panels():
    """forecaster=None is patch_ridge_forecaster; it equals the spectral backbone when the refinement is off."""
    from longbet import patch_ridge_forecaster
    from longbet._forecast import _spectral_damped_zero_shot

    # Large panel with a long context: the cross-series refinement is active and changes the forecast
    y, z, x, y_hist = _make_staggered_panel(seed=612, N=48, T=8, T_hist=24)
    f_def = build_causal_features(y, z, x=x, y_pre_history=y_hist, seed=41)
    f_pr = build_causal_features(y, z, x=x, y_pre_history=y_hist, forecaster=patch_ridge_forecaster, seed=41)
    f_sp = build_causal_features(y, z, x=x, y_pre_history=y_hist, forecaster="spectral", seed=41)
    np.testing.assert_array_equal(f_def.point_forecast, f_pr.point_forecast)
    np.testing.assert_array_equal(f_def.ps, f_pr.ps)
    assert f_def.point_forecast.shape == y.shape and np.all(np.isfinite(f_def.point_forecast))
    # Pre-anchor periods are the spectral fit under both backbones; post-anchor may differ
    np.testing.assert_allclose(f_def.point_forecast[:, : f_def.anchor_t], f_sp.point_forecast[:, : f_sp.anchor_t], atol=1e-5)

    # N < 12: the refinement is skipped and the backbone output is exactly the spectral forecast
    ctx = np.hstack([y_hist[:10], y[:10, :4]]).astype(np.float64)
    direct = patch_ridge_forecaster(ctx, 4)
    spectral, _ = _spectral_damped_zero_shot(ctx, 4, 0)
    np.testing.assert_array_equal(direct, spectral.astype(np.float32))
    assert patch_ridge_forecaster(ctx, 0).shape == (10, 0)

    # Forecast draws (bootstrap of the backbone + adapter) work with the default backbone
    f_draw = build_causal_features(y, z, x=x, y_pre_history=y_hist, n_forecast_draws=4, seed=41)
    assert f_draw._unfiltered_forecast_draws is not None
    assert f_draw._unfiltered_forecast_draws.shape == (4, y.shape[0], y.shape[1])
    np.testing.assert_allclose(f_draw._unfiltered_forecast_draws.mean(axis=0), f_draw.point_forecast, atol=1e-5)
    assert f_draw.forecast_draws is not None and float(f_draw.forecast_draws["post_anchor_sd_mean"]) > 1e-3
    assert "load_timesfm_forecaster" not in dir(__import__("longbet"))

    # Persistent innovations (as in monthly county employment) activate the refinement: the output then
    # differs from the spectral backbone on the first post-anchor steps and the adjustment decays with h
    rng = np.random.default_rng(7)
    N_ar, L_ar, rho = 60, 30, 0.9
    e = np.zeros((N_ar, L_ar))
    for t in range(1, L_ar):
        e[:, t] = rho * e[:, t - 1] + rng.normal(size=N_ar) * np.sqrt(1.0 - rho**2)
    tt = np.arange(L_ar)
    ctx_ar = 0.5 * rng.normal(size=(N_ar, 1)) + 0.8 * np.cos(2 * np.pi * tt / 6)[None, :] * (1 + 0.3 * rng.normal(size=(N_ar, 1))) + 0.6 * e
    pr = patch_ridge_forecaster(ctx_ar, 8).astype(np.float64)
    sp, _ = _spectral_damped_zero_shot(ctx_ar, 8, 0)
    gap = np.mean(np.abs(pr - sp), axis=0)
    assert gap[0] > 1e-4 and gap[-1] < gap[0]
