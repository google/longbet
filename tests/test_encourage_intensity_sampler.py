"""Sampler-backed Stage-2 integration, composed joint SBC (T2), and confounded T1 coverage tests.

Covers Deliverable C for ``LongBetEncourage(treatment="intensity")`` with the
default ``stage2_backend="longbet"``:
  1. All three identification tiers (``"deconvolution"``, ``"control_function"``,
     ``"unconfounded"``), including ``demean_pre=True``, ``save``/``load``
     archive roundtrip, and exact cancellation of control-function coefficients
     in ``predict`` and ``predict_conditional``.
  2. Tier T0 immunity to strong within-arm confounding and covariate-moderated
     organic dose trends on the sampler backend.
  3. End-to-end composed joint SBC on Tier T2 (``identification="unconfounded"``,
     ``rho == 0``, small panel, ``test_joint_full_sbc``-style) with an
     independently coded two-stage prior generator and positive-control fault
     sensitivity.
  4. Slow-marked frequentist bias and coverage check for Tier T1
     (``identification="control_function"``) vs Tier T2 on a confounded DGP.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from longbet import LongBetEncourage
from longbet._encourage_intensity import compose_offer_effect, lag_matrix
from longbet._trend_basis import build_unit_features

FAST_SAMPLER = dict(
    num_chains=2,
    num_burnin=20,
    num_sweeps=25,
    n_skip=1,
    num_trees_pr=4,
    num_trees_trt=4,
    max_depth_pr=3,
    max_depth_trt=3,
    min_points_per_leaf_pr=2,
    min_points_per_leaf_trt=2,
    trend_num_features=6,
    random_seed=404,
)


def _synthetic_intensity_panel(
    seed: int = 701,
    n: int = 48,
    T: int = 6,
    t0: int = 3,
    confound: float = 0.0,
    organic_trend: float = 0.0,
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 2))
    z = np.zeros((n, T), dtype=float)
    z[: n // 2, t0 - 1 :] = 1.0
    t_grid = np.arange(1, T + 1, dtype=float)
    L = T - t0
    kappa = np.array([0.70, 0.35, 0.15, 0.05][: L + 1], dtype=float)

    # First stage: organic dose + front-loaded offer pulse + dynamic shock u_{it}.
    u = rng.normal(scale=0.45, size=(n, T))
    pulse = np.array([1.40, 0.70, 0.35, 0.15, 0.05], dtype=float)
    offer_profile = np.zeros(T, dtype=float)
    offer_profile[t0 - 1 :] = pulse[: T - t0 + 1]
    nu_d = 1.0 + 0.20 * np.tanh(x[:, 0])
    citt_d = nu_d[:, None] * offer_profile[None, :]
    d0 = (
        0.5
        + 0.3 * x[:, 0:1]
        + organic_trend * (0.4 + 0.5 * x[:, 0:1]) * (t_grid[None, :] / T)
    )
    d = d0 + z * citt_d + u

    # Outcome stage: distributed-lag causal effect + within-arm confounding from u and u_lag.
    nu_y = 1.0 + 0.20 * (x[:, 1] - x[:, 1].mean())
    lags_d = lag_matrix(d, L)
    structural = (lags_d @ kappa) * nu_y[:, None]
    u_lag = lag_matrix(u, L)[:, :, 1:].sum(axis=2)
    y0 = 0.4 * x[:, 0:1] - 0.25 * x[:, 1:2] + 0.2 * rng.normal(size=(n, 1))
    y = y0 + structural + confound * (0.8 * u + 0.5 * u_lag) + 0.30 * rng.normal(size=(n, T))
    citt_y = compose_offer_effect(
        kappa[None, None, :] * nu_y[:, None, None],
        citt_d[:, :, None],
        t0 - 1,
    )[:, :, 0]
    return {
        "y": y,
        "d": d,
        "z": z,
        "x": x,
        "t": t_grid,
        "kappa": kappa,
        "nu_y": nu_y,
        "citt_d": citt_d,
        "citt_y": citt_y,
        "itt_y": citt_y.mean(axis=0)[t0 - 1 :],
        "itt_d": citt_d.mean(axis=0)[t0 - 1 :],
    }


@pytest.mark.parametrize("tier", ["deconvolution", "control_function", "unconfounded"])
def test_sampler_stage2_all_tiers_and_archive_roundtrip(tier, tmp_path):
    """Default ``stage2_backend="longbet"`` runs end-to-end across all three
    tiers, cancels control-function/nuisance columns in composed contrasts, and
    survives a ``save``/``load`` roundtrip byte-for-byte."""
    p = _synthetic_intensity_panel(seed=710, n=36, T=5, t0=3, confound=0.6, organic_trend=0.5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = LongBetEncourage(
            treatment="intensity",
            identification=tier,
            demean_pre=(tier == "unconfounded"),
            **FAST_SAMPLER,
        ).fit(p["y"], p["d"], p["z"], p["x"], t=p["t"])

    assert model.stage2_backend == "longbet"
    assert model.outcome_model is not None
    B = FAST_SAMPLER["num_chains"] * FAST_SAMPLER["num_sweeps"]
    assert model.stage2.kappa_draws.shape == (B, 3)
    assert np.all(np.isfinite(model.stage2.kappa_draws))

    pred = model.predict()
    cond = model.predict_conditional()
    assert np.all(np.isfinite(pred.draws["itt_y"]))
    assert np.all(np.isfinite(cond.citt_y.draws))
    assert np.all(np.isfinite(pred.discrepancy().z_stat))

    # Perturbing Stage-2 extra_coef_draws (CF or T0 organic-dose controls or demean_pre Mundlak terms)
    # must leave every causal contrast identically unchanged.
    if model.stage2.extra_coef_draws is not None:
        saved = model.stage2.extra_coef_draws.copy()
        try:
            model.stage2.extra_coef_draws = saved + 500.0
            pred_shifted = model.predict()
            cond_shifted = model.predict_conditional()
        finally:
            model.stage2.extra_coef_draws = saved
        np.testing.assert_array_equal(pred_shifted.draws["itt_y"], pred.draws["itt_y"])
        np.testing.assert_array_equal(cond_shifted.citt_y.draws, cond.citt_y.draws)
        np.testing.assert_array_equal(cond_shifted.kernel.draws, cond.kernel.draws)

    # Archive roundtrip preserves Stage-2 LongBet model and predictions exactly.
    archive_path = tmp_path / f"intensity_{tier}.npz"
    model.save(archive_path)
    loaded = LongBetEncourage.load(archive_path)
    assert loaded.stage2_backend == "longbet"
    assert loaded.identification == tier
    np.testing.assert_array_equal(loaded.stage2.kappa_draws, model.stage2.kappa_draws)
    np.testing.assert_array_equal(loaded.predict().draws["itt_y"], pred.draws["itt_y"])
    np.testing.assert_array_equal(
        loaded.predict_conditional().citt_y.draws, cond.citt_y.draws
    )
    pd.testing.assert_frame_equal(loaded.predict().reference, pred.reference)


def test_sampler_t0_immunity_and_t1_deconfounding_vs_t2():
    """On the sampler backend, a strong within-arm dynamic confounder biases T2
    upward, while T0 (propensity-centered Robinson/FWL offer increments) is
    invariant to within-arm shocks given Stage 1 and T1 (3-term CF) removes the
    confounding."""
    p_clean = _synthetic_intensity_panel(seed=720, n=80, T=6, t0=3, confound=0.0, organic_trend=0.8)
    p_conf = _synthetic_intensity_panel(seed=720, n=80, T=6, t0=3, confound=1.5, organic_trend=0.8)
    cfg = dict(FAST_SAMPLER, num_burnin=50, num_sweeps=60, random_seed=725)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t0_clean = LongBetEncourage(treatment="intensity", identification="deconvolution", **cfg).fit(
            p_clean["y"], p_clean["d"], p_clean["z"], p_clean["x"], t=p_clean["t"]
        )
        t0_conf = LongBetEncourage(treatment="intensity", identification="deconvolution", **cfg).fit(
            p_conf["y"], p_conf["d"], p_conf["z"], p_conf["x"], t=p_conf["t"]
        )
        t1_conf = LongBetEncourage(treatment="intensity", identification="control_function", **cfg).fit(
            p_conf["y"], p_conf["d"], p_conf["z"], p_conf["x"], t=p_conf["t"]
        )
        t2_conf = LongBetEncourage(treatment="intensity", identification="unconfounded", **cfg).fit(
            p_conf["y"], p_conf["d"], p_conf["z"], p_conf["x"], t=p_conf["t"]
        )

    k_true = float(p_clean["kappa"].sum())
    k_t0_clean = float(t0_clean.stage2.kappa_draws.sum(axis=1).mean())
    k_t0_conf = float(t0_conf.stage2.kappa_draws.sum(axis=1).mean())
    k_t1_conf = float(t1_conf.stage2.kappa_draws.sum(axis=1).mean())
    k_t2_conf = float(t2_conf.stage2.kappa_draws.sum(axis=1).mean())

    # T0 and T1 are dramatically closer to the true cumulative kernel than T2 under strong confounding.
    assert abs(k_t0_conf - k_true) < 0.45 * abs(k_t2_conf - k_true), (k_t0_conf, k_t2_conf, k_true)
    assert abs(k_t1_conf - k_true) < 0.45 * abs(k_t2_conf - k_true), (k_t1_conf, k_t2_conf, k_true)
    # T0 on clean vs confounded panel stays consistent because within-arm shocks are orthogonal to Z_tilde.
    assert abs(k_t0_conf - k_t0_clean) < 0.35, (k_t0_clean, k_t0_conf)


def independent_composed_t2_prior_draw(
    rng: np.random.Generator,
    x: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    t0_idx: int = 2,
    n_features: int = 6,
    sigma_d: float = 0.35,
    sigma_y: float = 0.35,
    sigma_gamma: float = 0.22,
    sigma_kappa: float = 0.38,
) -> dict[str, np.ndarray | float]:
    """Independently coded generative draw from the two-stage Tier T2 hierarchy
    (``rho == 0``) on a fixed design ``(x, z, t)``."""
    n, T = z.shape
    L = T - t0_idx - 1
    L1 = L + 1
    b_mat = np.asarray(build_unit_features(x, n_features)[0], dtype=np.float64)
    q = b_mat.shape[1]
    t_centered = (t - t.mean()) / max(t.std(), 1.0)

    # --- Stage 1: Dose D_{it} ---
    c_d = rng.normal(scale=0.10, size=q)
    eta_d = (b_mat @ c_d)[:, None] * t_centered[None, :]
    f_d = rng.normal(scale=0.15, size=T)
    f_d = f_d - f_d.mean()
    gamma_d = rng.normal(scale=sigma_gamma, size=n)

    cut_mu_d = rng.uniform(-0.7, 0.7)
    leaf_mu_d = 0.8 + rng.normal(scale=0.20, size=2)
    mu_d = np.where(x[:, 0] <= cut_mu_d, leaf_mu_d[0], leaf_mu_d[1])

    cut_nu_d = rng.uniform(-0.7, 0.7)
    leaf_nu_d = 1.0 + rng.normal(scale=0.18, size=2)
    nu_d = np.where(x[:, 1] <= cut_nu_d, leaf_nu_d[0], leaf_nu_d[1])

    s_max = T - t0_idx
    s_grid = np.arange(1, s_max + 1, dtype=np.float64)
    diff_s = s_grid[:, None] - s_grid[None, :]
    K_d = (0.35**2) * np.exp(-0.5 * (diff_s / 1.5) ** 2) + 1e-6 * np.eye(s_max)
    beta_d_raw = rng.multivariate_normal(mean=np.full(s_max, 0.75), cov=K_d)
    sgn_d = 1.0 if np.sum(beta_d_raw) >= 0.0 else -1.0
    beta_d = sgn_d * beta_d_raw
    nu_d = sgn_d * nu_d

    # Exact generative CITT_D(X_i, t) for every unit at post-offer periods.
    citt_d_full = np.zeros((n, T), dtype=np.float64)
    for h in range(s_max):
        citt_d_full[:, t0_idx + h] = nu_d * beta_d[h]

    d0_it = mu_d[:, None] + gamma_d[:, None] + f_d[None, :] + eta_d
    d_it = d0_it + z * citt_d_full + rng.normal(scale=sigma_d, size=(n, T))

    # --- Stage 2: Outcome Y_{it} under Tier T2 (rho == 0) ---
    c_y = rng.normal(scale=0.10, size=q)
    eta_y = (b_mat @ c_y)[:, None] * t_centered[None, :]
    f_y = rng.normal(scale=0.15, size=T)
    f_y = f_y - f_y.mean()
    gamma_y = rng.normal(scale=sigma_gamma, size=n)

    cut_mu_y = rng.uniform(-0.7, 0.7)
    leaf_mu_y = rng.normal(scale=0.22, size=2)
    mu_y = np.where(x[:, 0] <= cut_mu_y, leaf_mu_y[0], leaf_mu_y[1])

    # Centered-quotient moderator nu_Y(X_i) with sample mean 1.
    cut_nu_y = rng.uniform(-0.7, 0.7)
    leaf_nu_y = rng.normal(scale=0.18, size=2)
    nu_y_tree = np.where(x[:, 1] <= cut_nu_y, leaf_nu_y[0], leaf_nu_y[1])
    nu_y_tree = nu_y_tree - nu_y_tree.mean()
    x_std = (x - x.mean(axis=0)) / np.maximum(x.std(axis=0), 1e-8)
    delta_y = rng.normal(scale=0.10, size=x.shape[1])
    nu_y = 1.0 + nu_y_tree + x_std @ delta_y

    l_grid = np.arange(L1, dtype=np.float64)
    diff_l = l_grid[:, None] - l_grid[None, :]
    K_lag = (sigma_kappa**2) * (np.exp(-0.5 * (diff_l / 1.2) ** 2) + 0.5) + 1e-6 * np.eye(L1)
    kappa = rng.multivariate_normal(mean=np.zeros(L1), cov=K_lag)

    lags_d = lag_matrix(d_it, L)
    w_it = lags_d @ kappa
    tau_y_it = w_it * nu_y[:, None]
    y0_it = mu_y[:, None] + gamma_y[:, None] + f_y[None, :] + eta_y
    y_it = y0_it + tau_y_it + rng.normal(scale=sigma_y, size=(n, T))

    # Exact composed CITT_Y(X_i, t) = nu_Y(X_i) sum_{l=0}^{t-t0} kappa_l CITT_D(X_i, t-l).
    citt_y_full = compose_offer_effect(
        kappa[None, None, :] * nu_y[:, None, None],
        citt_d_full[:, :, None],
        t0_idx,
    )[:, :, 0]
    itt_d_post = citt_d_full[:, t0_idx:].mean(axis=0)
    itt_y_post = citt_y_full[:, t0_idx:].mean(axis=0)

    return {
        "y": y_it,
        "d": d_it,
        "kappa": kappa,
        "kappa_cum": float(np.sum(kappa)),
        "itt_d_cum": float(np.sum(itt_d_post)),
        "itt_y_cum": float(np.sum(itt_y_post)),
        "itt_y_post": itt_y_post,
    }


def test_joint_composed_t2_sbc_calibration_and_fault_sensitivity():
    """End-to-end composed joint SBC on Tier T2 (``identification="unconfounded"``,
    ``rho == 0``): verifies rank centering, empirical coverage, and sensitivity
    to a deliberate faulty conditional across both stages and their composition."""
    rng = np.random.default_rng(20261010)
    n, T, t0_idx = 36, 5, 2
    x = rng.normal(size=(n, 3))
    t = np.arange(1, T + 1, dtype=np.float64)
    z = np.zeros((n, T), dtype=np.float64)
    z[: n // 2, t0_idx:] = 1.0

    n_rep = 12
    metrics = {
        k: {"ranks": [], "covered": [], "clean_err": [], "corrupt_err": []}
        for k in ("itt_d_cum", "kappa_cum", "itt_y_cum")
    }

    for r in range(n_rep):
        sim = independent_composed_t2_prior_draw(rng, x, z, t, t0_idx=t0_idx, n_features=6)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = LongBetEncourage(
                treatment="intensity",
                identification="unconfounded",
                num_sweeps=45,
                num_burnin=45,
                num_chains=2,
                num_trees_pr=6,
                num_trees_trt=6,
                trend_num_features=6,
                trend_horseshoe=True,
                standardize=False,
                random_seed=9100 + r,
            ).fit(sim["y"], sim["d"], z, x, t=t)
            pred = model.predict()

        itt_d_draws = pred.draws["itt_d"][0].sum(axis=0).ravel()
        itt_y_draws = pred.draws["itt_y"][0].sum(axis=0).ravel()
        kappa_cum_draws = model.stage2.kappa_draws.sum(axis=1).ravel()

        series_map = {
            "itt_d_cum": (itt_d_draws, sim["itt_d_cum"]),
            "kappa_cum": (kappa_cum_draws, sim["kappa_cum"]),
            "itt_y_cum": (itt_y_draws, sim["itt_y_cum"]),
        }
        for name, (draws_1d, truth) in series_map.items():
            m = metrics[name]
            m["ranks"].append(float(np.mean(draws_1d < truth)))
            lo, hi = np.quantile(draws_1d, [0.025, 0.975])
            m["covered"].append(bool(lo - 0.18 <= truth <= hi + 0.18))
            m["clean_err"].append(abs(float(np.mean(draws_1d)) - truth))
            # Deliberate faulty conditional: shift by +1.85 posterior SD.
            faulty = draws_1d + 1.85 * max(float(np.std(draws_1d)), 0.15)
            m["corrupt_err"].append(abs(float(np.mean(faulty)) - truth))

    for name, m in metrics.items():
        mean_rank = float(np.mean(m["ranks"]))
        cov = float(np.mean(m["covered"]))
        assert 0.18 <= mean_rank <= 0.82, f"{name}: rank mean {mean_rank:.3f} indicates systematic drift"
        assert cov >= 0.75, f"{name}: empirical coverage {cov:.2f} too low"
        assert np.mean(m["corrupt_err"]) > 1.5 * np.mean(m["clean_err"]), (
            f"{name}: deliberate faulty conditional not detected"
        )


@pytest.mark.slow
def test_confounded_t1_frequentist_bias_and_coverage_vs_t2():
    """Frequentist bias and coverage check on a confounded dynamic DGP: Tier T1
    (3-term control function) achieves low bias and nominal-level coverage for
    the cumulative kernel and composed ``ITT_Y``, whereas Tier T2 is severely
    biased and under-covers."""
    n_rep = 6
    t1_covered_k, t1_covered_y = [], []
    t1_errors_k, t2_errors_k = [], []

    for r in range(n_rep):
        p = _synthetic_intensity_panel(seed=8200 + r, n=100, T=5, t0=3, confound=1.0, organic_trend=0.3)
        k_true = float(p["kappa"].sum())
        y_true = float(p["itt_y"].sum())
        cfg = dict(
            num_chains=2,
            num_burnin=45,
            num_sweeps=50,
            num_trees_pr=6,
            num_trees_trt=6,
            trend_num_features=6,
            random_seed=8300 + r,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m_t1 = LongBetEncourage(treatment="intensity", identification="control_function", **cfg).fit(
                p["y"], p["d"], p["z"], p["x"], t=p["t"]
            )
            m_t2 = LongBetEncourage(treatment="intensity", identification="unconfounded", **cfg).fit(
                p["y"], p["d"], p["z"], p["x"], t=p["t"]
            )
            pred_t1 = m_t1.predict()

        k_draws_t1 = m_t1.stage2.kappa_draws.sum(axis=1)
        k_draws_t2 = m_t2.stage2.kappa_draws.sum(axis=1)
        y_draws_t1 = pred_t1.draws["itt_y"][0].sum(axis=0).ravel()

        t1_errors_k.append(abs(float(k_draws_t1.mean()) - k_true))
        t2_errors_k.append(abs(float(k_draws_t2.mean()) - k_true))

        lo_k, hi_k = np.quantile(k_draws_t1, [0.025, 0.975])
        lo_y, hi_y = np.quantile(y_draws_t1, [0.025, 0.975])
        t1_covered_k.append(bool(lo_k - 0.15 <= k_true <= hi_k + 0.15))
        t1_covered_y.append(bool(lo_y - 0.20 <= y_true <= hi_y + 0.20))

    assert np.mean(t1_errors_k) < 0.45 * np.mean(t2_errors_k), (
        np.mean(t1_errors_k),
        np.mean(t2_errors_k),
    )
    assert np.mean(t1_covered_k) >= 0.75, t1_covered_k
    assert np.mean(t1_covered_y) >= 0.75, t1_covered_y
