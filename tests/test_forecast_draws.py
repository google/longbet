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

"""Forecast draws (``n_forecast_draws``) and gate-aware recommended features."""

from __future__ import annotations

import numpy as np
import pytest

from longbet import LateContextSpec
from longbet._forecast import _moving_block_indices, build_causal_features

from test_latecontext import _panel  # noqa: E402  (same test directory)


@pytest.fixture(scope="module")
def panel():
    return _panel(1, 80, 10, 26, True)


def test_default_outputs_are_unchanged_when_draws_are_requested(panel):
    y, z, x, yh = panel
    f0 = build_causal_features(y, z, x=x, y_pre_history=yh, seed=1)
    f1 = build_causal_features(y, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=8)
    assert f0.baseline_forecast_draws is None and f0.forecast_draws is None
    np.testing.assert_array_equal(f0.point_forecast, f1.point_forecast)
    np.testing.assert_array_equal(f0.unit_features, f1.unit_features)
    np.testing.assert_array_equal(f0.ps, f1.ps)
    assert f0.beyond_poly_r2 == f1.beyond_poly_r2


def test_draws_are_centred_finite_and_dispersed(panel):
    y, z, x, yh = panel
    f = build_causal_features(y, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=16)
    d = f.baseline_forecast_draws
    assert f.gate_open and d is not None and d.shape == (16,) + y.shape and d.dtype == np.float32
    assert np.all(np.isfinite(d))
    # exactly centred on the point forecast (float32 rounding only)
    assert float(np.max(np.abs(d.mean(axis=0) - f.baseline_forecast))) < 1e-4
    sd_post = d.std(axis=0)[:, f.anchor_t:]
    assert 0.02 < float(sd_post.mean()) < 1.0
    diag = f.forecast_draws
    assert diag["n_draws"] == 16 and diag["gate_open"] is True and diag["block_length"] >= 2
    assert 0.0 <= diag["adapter_variance_share"] <= 1.0
    assert diag["post_anchor_sd_mean"] == pytest.approx(float(np.std(d - d.mean(axis=0), axis=0, ddof=1)[:, f.anchor_t:].mean()), rel=1e-3)


def test_draws_are_deterministic_in_seed(panel):
    y, z, x, yh = panel
    a = build_causal_features(y, z, x=x, y_pre_history=yh, seed=5, n_forecast_draws=6)
    b = build_causal_features(y, z, x=x, y_pre_history=yh, seed=5, n_forecast_draws=6)
    c = build_causal_features(y, z, x=x, y_pre_history=yh, seed=6, n_forecast_draws=6)
    np.testing.assert_array_equal(a.baseline_forecast_draws, b.baseline_forecast_draws)
    assert not np.array_equal(a.baseline_forecast_draws, c.baseline_forecast_draws)


@pytest.mark.parametrize("spec", [None, "none",
                                  LateContextSpec(break_kink=7.0, transient_centre=-2.5, transient_width=3.0)])
def test_draws_never_see_treated_post_anchor_outcomes(panel, spec):
    y, z, x, yh = panel
    f = build_causal_features(y, z, x=x, y_pre_history=yh, seed=2, n_forecast_draws=5, late_context=spec)
    y2 = y.copy()
    y2[z.any(axis=1), f.anchor_t:] += 50.0
    g = build_causal_features(y2, z, x=x, y_pre_history=yh, seed=2, n_forecast_draws=5, late_context=spec)
    np.testing.assert_array_equal(f.baseline_forecast_draws, g.baseline_forecast_draws)
    np.testing.assert_array_equal(f.point_forecast, g.point_forecast)


def test_gated_off_panel_has_no_draws_and_no_recommended_history_features():
    y, z, x, yh = _panel(3, 40, 10, 4, True)  # context of 8 periods: gate closed
    f = build_causal_features(y, z, x=x, y_pre_history=yh, seed=3, n_forecast_draws=4)
    assert not f.gate_open and f.baseline_forecast is None
    assert f.baseline_forecast_draws is None
    assert f.forecast_draws is not None and f.forecast_draws["gate_open"] is False
    assert f.recommended_unit_features is None and f.recommended_ps is None
    assert f.augment_x(x) is x
    assert f.augment_x(None) is None
    assert f.unit_features.shape == (40, 2) and f.ps.shape == (40,)


def test_gated_on_panel_recommends_its_history_features(panel):
    y, z, x, yh = panel
    f = build_causal_features(y, z, x=x, y_pre_history=yh, seed=1)
    assert f.gate_open
    np.testing.assert_array_equal(f.recommended_unit_features, f.unit_features)
    np.testing.assert_array_equal(f.recommended_ps, f.ps)
    np.testing.assert_array_equal(f.augment_x(x), np.hstack([x, f.unit_features]))
    np.testing.assert_array_equal(f.augment_x(None), f.unit_features)


def test_moving_block_indices_cover_length_and_respect_blocks():
    rng = np.random.default_rng(0)
    idx = _moving_block_indices(30, 4, rng)
    assert idx.shape == (30,) and idx.min() >= 0 and idx.max() < 30
    # inside each block indices are consecutive
    for s in range(0, 28, 4):
        blk = idx[s:s + 4]
        assert np.all(np.diff(blk) == 1)
    assert _moving_block_indices(5, 10, rng).shape == (5,)


def test_block_length_override_is_honoured(panel):
    y, z, x, yh = panel
    f = build_causal_features(y, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=3,
                                      forecast_draw_block_length=7)
    assert f.forecast_draws["block_length"] == 7


def test_control_calibration_inflates_and_is_leak_free(panel):
    y, z, x, yh = panel
    cal = build_causal_features(y, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=16)
    c = cal.forecast_draws["calibration"]
    assert c["applied"] is True and 1.0 <= c["kappa"] <= 3.0 and c["n_controls"] >= 6
    # the calibration uses only never-treated post-anchor outcomes: perturbing treated ones changes nothing,
    # perturbing a control's post-anchor outcomes changes kappa
    y_t = y.copy()
    y_t[z.any(axis=1), cal.anchor_t:] += 50.0
    again = build_causal_features(y_t, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=16)
    assert again.forecast_draws["calibration"]["kappa"] == c["kappa"]
    # (a per-unit constant shift is removed by the double-demeaning, so use a unit x time interaction)
    y_c = y.copy()
    ctl = np.flatnonzero(~z.any(axis=1))[:6]
    t_post = np.arange(y.shape[1] - cal.anchor_t, dtype=float)
    y_c[ctl, cal.anchor_t:] += 3.0 * np.outer(np.where(np.arange(6) % 2 == 0, 1.0, -1.0), t_post)
    moved = build_causal_features(y_c, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=16)
    m = moved.forecast_draws["calibration"]
    assert m["control_error_var_dd"] > c["control_error_var_dd"]
    assert m["kappa"] > c["kappa"]
    with pytest.raises(ValueError):
        build_causal_features(y, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=1)
    with pytest.raises(ValueError):
        build_causal_features(y, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=-1)
    with pytest.raises(ValueError):
        build_causal_features(y, z, x=x, y_pre_history=yh, seed=1, n_forecast_draws=4, draw_method="invalid")


def test_posterior_draws_centred_leak_free_and_match_closed_form_covariance(panel):
    y, z, x, yh = panel
    f = build_causal_features(
        y, z, x=x, y_pre_history=yh, seed=11, n_forecast_draws=1024, draw_method="posterior"
    )
    d = f.point_forecast_draws
    assert d is not None and d.shape == (1024,) + y.shape and d.dtype == np.float32
    assert float(np.max(np.abs(d.mean(axis=0) - f.point_forecast))) < 1e-4

    diag = f.forecast_draws
    assert diag["draw_method"] == "posterior"
    assert diag["calibration"]["applied"] is False and diag["calibration"]["kappa"] == 1.0

    # Empirical post-anchor variance and covariance match the exact closed-form NIG formulas
    ana_var = diag["analytical_post_anchor_var_mean"]
    emp_var = diag["empirical_post_anchor_var_mean"]
    assert ana_var > 0.0
    assert emp_var == pytest.approx(ana_var, rel=0.08)

    post_mean_draws = d[:, :, f.anchor_t:].astype(np.float64).mean(axis=2)  # (M, N)
    emp_cov = np.cov(post_mean_draws, rowvar=False, ddof=1)
    ana_cov = diag["analytical_cov_post_mean"]
    rel_frob = np.linalg.norm(emp_cov - ana_cov) / max(np.linalg.norm(ana_cov), 1e-12)
    assert rel_frob < 0.12

    # Proposition 1: perturbing treated post-anchor outcomes leaves posterior draws and loadings bit-identical
    y2 = y.copy()
    y2[z.any(axis=1), f.anchor_t:] += 100.0
    g = build_causal_features(
        y2, z, x=x, y_pre_history=yh, seed=11, n_forecast_draws=16, draw_method="posterior"
    )
    f16 = build_causal_features(
        y, z, x=x, y_pre_history=yh, seed=11, n_forecast_draws=16, draw_method="posterior"
    )
    np.testing.assert_array_equal(f16.point_forecast_draws, g.point_forecast_draws)
    np.testing.assert_array_equal(f16.pre_anchor_loadings, g.pre_anchor_loadings)


def test_pre_anchor_loadings_shape_and_standardization(panel):
    y, z, x, yh = panel
    f = build_causal_features(y, z, x=x, y_pre_history=yh, seed=1, n_factor_loadings=3)
    lam = f.pre_anchor_loadings
    assert lam is not None and lam.shape == (y.shape[0], 3) and lam.dtype == np.float32
    np.testing.assert_allclose(lam.mean(axis=0), np.zeros(3), atol=1e-5)
    np.testing.assert_allclose(lam.std(axis=0), np.ones(3), atol=1e-5)


