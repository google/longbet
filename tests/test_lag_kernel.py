# Copyright 2026 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""API and edge-case tests of the distributed-lag exposure kernel.

``LongBet.fit(exposure_lags=...)`` replaces the absorbing clock by
``w_it = sum_l D_{i,t-l} kappa_l``; these tests cover the public surface
(``kappa_draws_``, ``extra_coef_draws_``, ``lag_scale_draws_``,
``moderator_draws``, ``predict``, ``save``/``load``), the degenerate panels
(zero-dose cells, no dose at all, ``L = 0``, missing cells, float64 inputs)
and a short recovery check of a known kernel and extra-column coefficient.
"""

import warnings

import numpy as np
import pytest

from longbet import LongBet, LongBetConfig
from longbet._model import lag_dose_matrix

SMALL = dict(num_sweeps=30, num_burnin=30, num_chains=2, num_trees_pr=4, num_trees_trt=6,
             trend_num_features=4, random_seed=10011)


def _panel(seed=10001, N=30, T=5, L=2, missing=0.0, kappa=(0.6, 0.3, 0.1), rho=0.0,
           dose_rate=0.6, dtype=np.float64):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, 2))
    dose = rng.gamma(1.5, 1.0, size=(N, T)) * (rng.random((N, T)) < dose_rate)
    lags = lag_dose_matrix(dose, L)
    kap = np.asarray(kappa, dtype=np.float64)[: L + 1]
    w = lags @ kap
    nu = 1.0 + 0.3 * np.tanh(x[:, 1])
    cf = rng.normal(size=(N, T, 1))
    y = (0.5 * x[:, :1] + 0.4 * rng.normal(size=(N, 1)) + w * nu[:, None] + rho * cf[:, :, 0]
         + 0.3 * rng.normal(size=(N, T)))
    if missing:
        y = y.copy()
        y[rng.random((N, T)) < missing] = np.nan
    return dict(y=y.astype(dtype), x=x.astype(dtype), lags=lags.astype(dtype), dose=dose,
                cf=cf.astype(dtype), nu=nu, kappa=kap, w=w)


def _fit(p, extra=False, **cfg):
    kw = dict(SMALL, **cfg)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return LongBet(LongBetConfig(**kw)).fit(
            p["y"], p["x"], exposure_lags=p["lags"], extra_cols=p["cf"] if extra else None)


def test_lag_dose_matrix_pads_pre_panel_lags_with_zeros():
    d = np.arange(1, 7, dtype=float).reshape(2, 3)          # [[1,2,3],[4,5,6]]
    out = lag_dose_matrix(d, 2)
    assert out.shape == (2, 3, 3) and out.dtype == np.float32
    np.testing.assert_array_equal(out[0, :, 0], [1, 2, 3])
    np.testing.assert_array_equal(out[0, :, 1], [0, 1, 2])
    np.testing.assert_array_equal(out[0, :, 2], [0, 0, 1])
    # More lags than periods: the extra columns are all zero.
    assert np.all(lag_dose_matrix(d, 5)[:, :, 3:] == 0)
    with pytest.raises(ValueError):
        lag_dose_matrix(d, -1)


def test_fit_predict_accessors_and_roundtrip(tmp_path):
    p = _panel()
    m = _fit(p, extra=True, anchor_moderator=True, ar1_errors=True, linear_trt=True,
             split_calendar_mu=False)
    B = SMALL["num_sweeps"] * SMALL["num_chains"]
    assert m.lag_num_ == 3 and m.extra_num_ == 1
    kappa = m.kappa_draws_
    assert kappa.shape == (B, 3) and np.all(np.isfinite(kappa))
    assert m.extra_coef_draws_.shape == (B, 1)
    assert m.lag_scale_draws_.shape == (B,) and np.all(m.lag_scale_draws_ > 0)
    # Outcome units: the trace is standardized, the accessor rescales by sdy.
    np.testing.assert_allclose(kappa, np.asarray(m.trace.kappa).reshape(B, 3) * m.sdy, rtol=1e-6)
    nu = m.moderator_draws(p["x"])
    assert nu.shape == (B, p["x"].shape[0])
    # Anchored: the moderator averages to one over the fitted units, per draw.
    np.testing.assert_allclose(nu.mean(axis=1), 1.0, atol=1e-4)
    # New units of a different size go through the fitted design.
    assert m.moderator_draws(p["x"][:7]).shape == (B, 7)

    pred = m.predict(p["x"], exposure_lags=p["lags"], extra_cols=p["cf"])
    N, T = p["y"].shape
    assert pred.tauhats.shape == (N, T, B)
    # tau = (lags @ kappa) * nu, in outcome units.
    w = np.einsum("ntl,bl->ntb", p["lags"], kappa)
    np.testing.assert_allclose(pred.tauhats, w * nu.T[:, None, :], rtol=1e-4, atol=1e-4)
    # The extra columns enter the factual outcome but not the contrast.
    ex = np.einsum("ntq,bq->ntb", p["cf"], m.extra_coef_draws_)
    np.testing.assert_allclose(pred.yhats - pred.muhats0 - pred.tauhats, ex, rtol=1e-4, atol=1e-4)
    assert pred.att_full.shape == (1, B) and np.all(np.isfinite(pred.att_full))
    with pytest.raises(ValueError, match="exposure_lags"):
        m.predict(p["x"])
    with pytest.raises(ValueError, match="extra_cols"):
        m.predict(p["x"], exposure_lags=p["lags"])

    path = tmp_path / "lag.npz"
    m.save(path)
    m2 = LongBet.load(path)
    assert m2.lag_num_ == 3 and m2.extra_num_ == 1
    np.testing.assert_array_equal(m2.kappa_draws_, kappa)
    np.testing.assert_array_equal(m2.extra_coef_draws_, m.extra_coef_draws_)
    np.testing.assert_array_equal(m2.lag_scale_draws_, m.lag_scale_draws_)
    pred2 = m2.predict(p["x"], exposure_lags=p["lags"], extra_cols=p["cf"])
    np.testing.assert_array_equal(pred2.tauhats, pred.tauhats)
    np.testing.assert_allclose(m2.moderator_draws(p["x"]), nu, rtol=1e-6)


def test_no_extra_cols_gives_none_and_z_is_ignored():
    p = _panel()
    z = np.ones_like(p["y"])
    m = _fit(p)
    assert m.extra_coef_draws_ is None and m.extra_num_ == 0
    # A z passed alongside exposure_lags is ignored by the sampler.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m_z = LongBet(LongBetConfig(**SMALL)).fit(p["y"], p["x"], z, exposure_lags=p["lags"])
    np.testing.assert_array_equal(m_z.kappa_draws_, m.kappa_draws_)
    with pytest.raises(ValueError, match="z is required"):
        LongBet(LongBetConfig(**SMALL)).fit(p["y"], p["x"])


@pytest.mark.parametrize("case", ["zero_dose_cells", "all_zero", "L0", "unbalanced", "float32"])
def test_degenerate_panels(case):
    if case == "zero_dose_cells":
        p = _panel(dose_rate=0.25)                      # most cells carry no dose
    elif case == "all_zero":
        p = _panel(dose_rate=0.0)                       # no dose anywhere
    elif case == "L0":
        p = _panel(L=0)
    elif case == "unbalanced":
        p = _panel(missing=0.15)
    else:
        p = _panel(dtype=np.float32)
    m = _fit(p, anchor_moderator=True)
    kappa = m.kappa_draws_
    assert np.all(np.isfinite(kappa)) and kappa.shape[1] == p["lags"].shape[2]
    pred = m.predict(p["x"], exposure_lags=p["lags"], summary_only=True)
    assert np.all(np.isfinite(np.asarray(pred.tau_summary.mean)))
    if case == "all_zero":
        # Nothing to learn from: the kernel is a prior draw and the effect is zero.
        assert np.all(np.asarray(pred.tau_summary.mean) == 0.0)
        assert np.all(np.isnan(pred.att_full))
    if case == "zero_dose_cells":
        # Zero-dose cells have exactly zero effect for every draw.
        zero = p["lags"].sum(axis=2) == 0
        assert np.all(np.asarray(pred.tau_summary.mean)[zero] == 0.0)
    if case == "unbalanced":
        assert np.all(np.isfinite(np.asarray(pred.y_summary.mean)))


def test_recovers_kernel_and_extra_coefficient():
    """Known kernel and control-function coefficient, moderate N: posterior
    means within a few posterior standard deviations and well inside the
    prior's range (a loose recovery check; calibration lives in the Geweke /
    SBC tests)."""
    p = _panel(seed=10002, N=120, T=6, L=2, kappa=(0.6, 0.3, 0.1), rho=0.8)
    m = _fit(p, extra=True, anchor_moderator=True, ar1_errors=True, split_calendar_mu=False,
             num_sweeps=150, num_burnin=150)
    kappa, rho = m.kappa_draws_, m.extra_coef_draws_[:, 0]
    for l in range(3):
        err = abs(kappa[:, l].mean() - p["kappa"][l])
        assert err < max(4 * kappa[:, l].std(), 0.15), (l, kappa[:, l].mean(), kappa[:, l].std())
    assert abs(rho.mean() - 0.8) < max(4 * rho.std(), 0.15), (rho.mean(), rho.std())


def test_lag_kernel_rejects_unsupported_configurations():
    p = _panel()
    with pytest.raises(ValueError, match="continuous"):
        LongBet(LongBetConfig(**SMALL, outcome="binary")).fit(
            (p["y"] > 0).astype(float), p["x"], exposure_lags=p["lags"])
    with pytest.raises(ValueError, match="use_trend_block"):
        LongBet(LongBetConfig(**SMALL, use_trend_block=False, split_calendar_mu=True)).fit(
            p["y"], p["x"], exposure_lags=p["lags"], extra_cols=p["cf"])
    with pytest.raises(ValueError, match="exposure_lags"):
        LongBet(LongBetConfig(**SMALL)).fit(p["y"], p["x"], np.ones_like(p["y"]), extra_cols=p["cf"])
    # Signed lag inputs (propensity-centered T0 or pre-offer demeaned doses) are valid.
    signed = p["lags"] - p["lags"].mean(axis=(0, 1), keepdims=True)
    m_signed = _fit({**p, "lags": signed}, anchor_moderator=True)
    assert np.all(np.isfinite(m_signed.kappa_draws_))
    bad_nan = p["lags"].copy(); bad_nan[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        LongBet(LongBetConfig(**SMALL)).fit(p["y"], p["x"], exposure_lags=bad_nan)
    dead = np.zeros_like(p["cf"])
    with pytest.raises(ValueError, match="identically zero"):
        LongBet(LongBetConfig(**SMALL)).fit(p["y"], p["x"], exposure_lags=p["lags"], extra_cols=dead)
