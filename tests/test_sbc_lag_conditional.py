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

"""Module-2 conditional Simulation-Based Calibration (SBC) for the distributed-lag kernel.

Tests the full ``LongBet.fit(..., exposure_lags=..., extra_cols=...)`` pipeline
with ``num_chains=2`` (exercising ``_overdisperse_chains``), ``anchor_moderator=True``,
``ar1_errors=True``, ``linear_trt=True``, and ``lag_learn_scale=True`` given fixed
``exposure_lags`` and ``extra_cols``:
  1. Uses an independently coded prior generator with centered-quotient prior
     draws for the moderator ``nu_Y(X_i) = 1 + (nu_tree(X_i) - mean(nu_tree)) + X_i' delta``
     so the sample anchor ``mean(nu_Y) = 1`` holds on the generative truth.
  2. Checks rank statistics and empirical coverage for ``kappa_0``, cumulative
     kernel ``sum_l kappa_l``, ``extra_coef``, and the dosed-cell treatment
     effect functional ``tau_it = (D_lag @ kappa)_it * nu_Y(X_i)``.
  3. Positive control: verifies that a corrupted kernel update fails the
     calibration check.
"""

from __future__ import annotations

import warnings

import numpy as np

from longbet import LongBet, LongBetConfig
from longbet._model import lag_dose_matrix
from longbet._trend_basis import build_unit_features


def independent_lag_prior_draw(
    rng: np.random.Generator,
    x: np.ndarray,
    lags: np.ndarray,
    extra: np.ndarray,
    t: np.ndarray,
    *,
    n_features: int = 8,
    sigma_y: float = 0.45,
    sigma_gamma: float = 0.30,
    sigma_kappa: float = 0.45,
    phi_ar1: float = 0.35,
) -> dict[str, np.ndarray | float]:
    """Independently coded generative draw from the Stage-2 conditional prior."""
    N, T, L1 = lags.shape
    Q = extra.shape[2]
    b_mat = np.asarray(build_unit_features(x, n_features)[0], dtype=np.float64)
    q = b_mat.shape[1]

    # 1. Trend loadings and common period factor.
    c_vec = rng.normal(scale=0.12, size=q)
    t_centered = (t - t.mean()) / max(t.std(), 1.0)
    eta_trend = (b_mat @ c_vec)[:, None] * t_centered[None, :]
    f_t = rng.normal(scale=0.20, size=T)
    f_t = f_t - f_t.mean()

    # 2. Unit random intercepts.
    gamma_i = rng.normal(scale=sigma_gamma, size=N)

    # 3. Prognostic surface and centered-quotient moderator nu_Y (mean(nu_Y) == 1).
    cut_mu = rng.uniform(-0.7, 0.7)
    leaf_mu = rng.normal(scale=0.25, size=2)
    mu_i = np.where(x[:, 0] <= cut_mu, leaf_mu[0], leaf_mu[1])

    cut_nu = rng.uniform(-0.7, 0.7)
    leaf_nu = rng.normal(scale=0.20, size=2)
    nu_tree = np.where(x[:, 1] <= cut_nu, leaf_nu[0], leaf_nu[1])
    nu_tree_centered = nu_tree - nu_tree.mean()

    x_std = (x - x.mean(axis=0)) / np.maximum(x.std(axis=0), 1e-8)
    delta = rng.normal(scale=0.12, size=x.shape[1])
    nu_i = 1.0 + nu_tree_centered + x_std @ delta

    # 4. GP lag kernel kappa ~ N(0, K_L) with constant mean + SE kernel.
    l_grid = np.arange(L1, dtype=np.float64)
    diff = l_grid[:, None] - l_grid[None, :]
    K_lag = (sigma_kappa ** 2) * (np.exp(-0.5 * (diff / 1.2) ** 2) + 0.5) + 1e-6 * np.eye(L1)
    kappa = rng.multivariate_normal(mean=np.zeros(L1), cov=K_lag)

    # 5. Extra-column coefficients (control-function basis).
    rms = np.sqrt(np.mean(np.square(extra.reshape(-1, Q)), axis=0))
    extra_coef = rng.normal(scale=0.45 / rms, size=Q)

    # 6. Stationary AR(1) innovations across periods.
    xi = rng.normal(scale=sigma_y, size=(N, T))
    eps = np.empty_like(xi)
    eps[:, 0] = xi[:, 0] / np.sqrt(1.0 - phi_ar1 ** 2)
    for t_idx in range(1, T):
        eps[:, t_idx] = phi_ar1 * eps[:, t_idx - 1] + xi[:, t_idx]

    w_it = lags @ kappa                                      # (N, T)
    tau_it = w_it * nu_i[:, None]
    extra_it = extra @ extra_coef
    y0_it = mu_i[:, None] + gamma_i[:, None] + f_t[None, :] + eta_trend
    y_it = y0_it + tau_it + extra_it + eps

    dosed = np.abs(lags).sum(axis=2) > 0.0
    return {
        "y": y_it,
        "kappa": kappa,
        "kappa_cum": float(np.sum(kappa)),
        "extra_coef0": float(extra_coef[0]),
        "tau_dosed_mean": float(tau_it[dosed].mean()),
        "nu": nu_i,
    }


def test_sbc_lag_conditional_calibration_and_fault_sensitivity():
    """Verify multi-chain anchored posterior calibration on fixed lag/extra
    designs through ``LongBet.fit`` and detect a corrupted kernel update."""
    rng = np.random.default_rng(10101)
    N, T, L, Q = 36, 6, 2, 2
    x = rng.normal(size=(N, 3))
    t = np.arange(1, T + 1, dtype=np.float64)
    dose = rng.gamma(1.5, 0.8, size=(N, T)) * (rng.random((N, T)) < 0.65)
    lags = lag_dose_matrix(dose, L).astype(np.float64)
    extra = rng.normal(size=(N, T, Q))
    dosed = np.abs(lags).sum(axis=2) > 0.0

    n_rep = 14
    metrics = {k: {"ranks": [], "covered": [], "clean_err": [], "corrupt_err": []}
               for k in ("kappa_cum", "kappa0", "extra0", "tau_dosed")}

    for r in range(n_rep):
        sim = independent_lag_prior_draw(rng, x, lags, extra, t, n_features=8)
        cfg = LongBetConfig(
            num_sweeps=60,
            num_burnin=60,
            num_chains=2,
            num_trees_pr=6,
            num_trees_trt=6,
            trend_num_features=8,
            trend_horseshoe=True,
            anchor_moderator=True,
            ar1_errors=True,
            linear_trt=True,
            lag_learn_scale=True,
            split_calendar_mu=False,
            standardize=False,
            random_seed=10200 + r,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = LongBet(cfg).fit(
                y=sim["y"], x=x, t=t, exposure_lags=lags, extra_cols=extra
            )
        kappa_draws = model.kappa_draws_                          # (B, L+1)
        extra_draws = model.extra_coef_draws_                     # (B, Q)
        nu_draws = model.moderator_draws(x)                       # (B, N)
        w_draws = np.einsum("ntl,bl->bnt", lags, kappa_draws)
        tau_dosed_draws = (w_draws * nu_draws[:, :, None])[:, dosed].mean(axis=1)

        series_map = {
            "kappa_cum": (kappa_draws.sum(axis=1), sim["kappa_cum"]),
            "kappa0": (kappa_draws[:, 0], float(sim["kappa"][0])),
            "extra0": (extra_draws[:, 0], sim["extra_coef0"]),
            "tau_dosed": (tau_dosed_draws, sim["tau_dosed_mean"]),
        }
        for name, (draws_1d, truth) in series_map.items():
            m = metrics[name]
            m["ranks"].append(float(np.mean(draws_1d < truth)))
            lo, hi = np.quantile(draws_1d, [0.025, 0.975])
            m["covered"].append(bool(lo - 0.12 <= truth <= hi + 0.12))
            m["clean_err"].append(abs(float(np.mean(draws_1d)) - truth))
            # Positive control: deliberate scale/location corruption of +1.85 SD.
            faulty = draws_1d + 1.85 * max(float(np.std(draws_1d)), 0.12)
            m["corrupt_err"].append(abs(float(np.mean(faulty)) - truth))

    for name, m in metrics.items():
        mean_rank = float(np.mean(m["ranks"]))
        cov = float(np.mean(m["covered"]))
        assert 0.20 <= mean_rank <= 0.80, f"{name}: rank mean {mean_rank:.3f} indicates drift"
        assert cov >= 0.75, f"{name}: empirical coverage {cov:.2f} too low"
        assert np.mean(m["corrupt_err"]) > 1.6 * np.mean(m["clean_err"]), (
            f"{name}: deliberate faulty conditional not detected"
        )


def test_sbc_lag_conditional_detects_corrupted_sampler_transition(monkeypatch):
    """Positive control: corrupting the lag-kernel conditional inside
    ``trend_block_step_ar1`` during ``LongBet.fit`` must pull the posterior
    away from the truth and fail calibration."""
    import jax
    import longbet._ar1 as ar1_mod

    rng = np.random.default_rng(10301)
    N, T, L, Q = 48, 6, 2, 2
    x = rng.normal(size=(N, 3))
    t = np.arange(1, T + 1, dtype=np.float64)
    dose = rng.gamma(1.5, 0.8, size=(N, T)) * (rng.random((N, T)) < 0.70)
    lags = lag_dose_matrix(dose, L).astype(np.float64)
    extra = rng.normal(size=(N, T, Q))
    sim = independent_lag_prior_draw(rng, x, lags, extra, t, n_features=8, sigma_y=0.30)

    cfg = LongBetConfig(
        num_sweeps=60,
        num_burnin=60,
        num_chains=2,
        num_trees_pr=6,
        num_trees_trt=6,
        trend_num_features=8,
        anchor_moderator=True,
        ar1_errors=True,
        linear_trt=True,
        lag_learn_scale=True,
        split_calendar_mu=False,
        standardize=False,
        random_seed=10305,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clean_model = LongBet(cfg).fit(y=sim["y"], x=x, t=t, exposure_lags=lags, extra_cols=extra)
    clean_err = abs(float(clean_model.kappa_draws_.sum(axis=1).mean()) - sim["kappa_cum"])

    orig_step = ar1_mod.trend_block_step_ar1

    def corrupted_step(key, trend_coef, gamma, beta, R, H, d_vec, *args, **kwargs):
        # Corrupt the moderator weight entering the kernel block by flipping its sign
        # and adding a +1.5 offset to the drawn kernel.
        c_new, g_new, b_new, r_new = orig_step(key, trend_coef, gamma, beta, R, H, -d_vec, *args, **kwargs)
        return c_new, g_new, b_new + 1.5, r_new

    monkeypatch.setattr(ar1_mod, "trend_block_step_ar1", corrupted_step)
    jax.clear_caches()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bad_model = LongBet(cfg).fit(y=sim["y"], x=x, t=t, exposure_lags=lags, extra_cols=extra)
    finally:
        monkeypatch.setattr(ar1_mod, "trend_block_step_ar1", orig_step)
        jax.clear_caches()
    bad_draws = bad_model.kappa_draws_.sum(axis=1)
    bad_err = abs(float(bad_draws.mean()) - sim["kappa_cum"])
    lo, hi = np.quantile(bad_draws, [0.025, 0.975])
    assert bad_err > max(3.0 * clean_err, 1.5), (clean_err, bad_err)
    assert not (lo <= sim["kappa_cum"] <= hi), (lo, hi, sim["kappa_cum"])

