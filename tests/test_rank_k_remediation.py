"""Tests for the remediation of the Rank-K estimator (methodology review, Tier 1).

Covers: (1a) the ATT prior / truncation enforced inside the JAX sampler through ATT-weight
centering (``tau_intercept == ATT`` draw by draw, no post-hoc clipping) and the non-Gaussian
prior families; (1b) the plug-in engine refusing priors; (1c) the ``catt_floor`` projection
contract; (1d) ``reconcile_nested_3arm_draws`` with an explicit prior scale, SIR diagnostics
and prior draws; (2) cohorts with ``L_g < K`` excluded from the higher modes (bit-identical
mode-2 draws with and without the short cohort present); (3) budget / prior / variance
controls of the sampler and the hybrid estimator; (4) exposure-basis families; and the
bit-identity of paper-default fits against golden fixtures generated from the pre-remediation
code (``tests/data/rank_k_golden/make_golden.py``).
"""

from __future__ import annotations

from dataclasses import replace
import json
import os
from pathlib import Path
import platform
import sys
import warnings

import numpy as np
import pytest

os.environ.setdefault("JAX_PLATFORMS", "cpu")
THIS_DIR = Path(__file__).resolve().parent
FIXTURE_DIR = THIS_DIR / "data" / "rank_k_golden"
sys.path.insert(0, str(FIXTURE_DIR))

import longbet._rank_k as M  # noqa: E402
import longbet  # noqa: E402
from make_golden import cases as golden_cases, make_panel  # noqa: E402

RankKLongBetConfig = M.RankKLongBetConfig
fit = M.fit_rank_k_longbet

PROP = dict(
    uncertainty_propagation=True, n_forecast_draws=3, max_forecast_draws=3, num_chains=2, draws_per_component=6,
    num_burnin=20, num_trees_mu=20, num_trees_tau_mode1=12, num_trees_tau_higher=8, random_seed=2026,
)


def _has_stochtree() -> bool:
    try:
        import stochtree  # noqa: F401
    except ImportError:
        return False
    return True


@pytest.fixture(scope="module")
def panel():
    return make_panel("stag_unequal")


@pytest.fixture(scope="module")
def fit_kwargs(panel):
    return dict(
        baseline_forecast=panel["f"], baseline_forecast_draws=panel["f_draws"], cohorts=panel["cohorts"],
        coords=panel["coords"],
    )


@pytest.fixture(scope="module")
def res_default(panel, fit_kwargs):
    return fit(panel["y"], panel["z"], panel["x"], cfg=RankKLongBetConfig(**PROP, rank_k=2), **fit_kwargs)


def _run(panel, fit_kwargs, cfg, **extra):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res = fit(panel["y"], panel["z"], panel["x"], cfg=cfg, **{**fit_kwargs, **extra})
    msgs = [str(w.message) for w in caught if issubclass(w.category, M.RankKConfigWarning)]
    return res, msgs


# ---- Golden bit-identity of paper-default fits ---------------------------------------------------


@pytest.mark.parametrize("name,kind,cfg_kw,use_coords", list(golden_cases()), ids=[c[0] for c in golden_cases()])
def test_default_fits_match_pre_remediation_golden_fixtures(name, kind, cfg_kw, use_coords):
    if cfg_kw.get("uncertainty_propagation") is False and not _has_stochtree():
        pytest.skip("stochtree not installed")
    import jax

    meta = json.loads((FIXTURE_DIR / "metadata.json").read_text())
    same_env = meta["jax"] == jax.__version__ and meta["machine"] == platform.machine()
    p = make_panel(kind)
    kw = dict(baseline_forecast=p["f"], cohorts=p["cohorts"])
    if cfg_kw.get("uncertainty_propagation"):
        kw["baseline_forecast_draws"] = p["f_draws"]
    if use_coords:
        kw["coords"] = p["coords"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = fit(p["y"], p["z"], p["x"], cfg=RankKLongBetConfig(**cfg_kw), **kw)
    ref = np.load(FIXTURE_DIR / f"{name}.npz")
    got = {"att_draws": res.att_draws, "catt_draws": res.catt_draws, "theta_f_draws_0": res.theta_f_draws[0]}
    for k, td in enumerate(res.mode_tau_draws):
        got[f"mode_tau_draws_{k}"] = td
    for key, val in got.items():
        if same_env:
            np.testing.assert_array_equal(ref[key], np.asarray(val), err_msg=f"{name}:{key} is not bit-identical")
        else:  # another JAX build / CPU: allow floating-point reordering noise only
            np.testing.assert_allclose(ref[key], np.asarray(val), rtol=1e-5, atol=1e-6, err_msg=f"{name}:{key}")
    assert res.projection_report is None and res.att_draws_raw is None


# ---- Issue 1a: ATT prior / truncation inside the sampler -----------------------------------------


def test_att_truncation_is_enforced_by_the_sampler_on_an_unbalanced_panel(panel, fit_kwargs, res_default):
    cfg = RankKLongBetConfig(**PROP, rank_k=2, tau_clip=(0.0, float("inf")), tau_prior_sd=0.5)
    res, msgs = _run(panel, fit_kwargs, cfg)
    assert msgs == []
    assert np.all(res.att_draws >= -1e-6)
    rep = res.projection_report
    assert rep is not None and rep["att_shift_max_abs"] < 1e-5 and rep["att_draws_outside_clip_share"] == 0.0
    # tau_intercept is the ATT draw by draw (ATT-weight centering), which the default w-centering is not
    att_raw = res.att_draws_raw
    assert att_raw is not None
    np.testing.assert_allclose(res.diagnostics["tau_intercept_draws"][0], att_raw, atol=1e-4)
    default_gap = np.max(np.abs(res_default.att_draws - res_default.diagnostics["tau_intercept_draws"][0]))
    assert default_gap > 1e-6  # unequal windows: likelihood-weighted centering != ATT weights
    # K-invariance of the constrained ATT draws
    res_k1, _ = _run(panel, fit_kwargs, replace(cfg, rank_k=1))
    np.testing.assert_allclose(res_k1.att_draws, res.att_draws, atol=1e-12)


@pytest.mark.parametrize("family", ["halfstudentt", "gamma_spike"])
def test_non_gaussian_att_priors_sample_on_their_support(panel, fit_kwargs, family):
    res, _ = _run(panel, fit_kwargs, RankKLongBetConfig(**PROP, rank_k=1, tau_prior_family=family, tau_prior_sd=0.5))
    assert np.all(np.isfinite(res.att_draws)) and np.all(res.att_draws >= -1e-6)
    assert np.std(res.att_draws) > 0.0
    with pytest.raises(ValueError, match="needs tau_prior_sd"):
        RankKLongBetConfig(tau_prior_family=family)


# ---- Issue 1b / 1c: plug-in engine and catt_floor contract -----------------------------------------


def test_plugin_engine_refuses_att_priors(panel):
    cfg = RankKLongBetConfig(engine="stochtree", uncertainty_propagation=False, tau_prior_sd=0.5)
    with pytest.raises(ValueError, match="require engine='jax'"):
        fit(panel["y"], panel["z"], panel["x"], cfg=cfg, baseline_forecast=panel["f"], cohorts=panel["cohorts"])


def test_catt_floor_is_a_declared_projection(panel, fit_kwargs, res_default):
    res, msgs = _run(panel, fit_kwargs, RankKLongBetConfig(**PROP, rank_k=2, catt_floor=0.0))
    assert any("decision-layer projection" in m for m in msgs)
    np.testing.assert_allclose(res.att_draws, res_default.att_draws, atol=1e-10)  # positive effects: floor never binds a cohort mean
    assert np.all(res.catt_draws[:, panel["z"] == 1] >= 0.0)
    assert res.catt_draws_raw is not None and res.mode_tau_draws_raw is not None
    np.testing.assert_array_equal(res.mode_tau_draws_raw[0], res_default.mode_tau_draws[0])
    np.testing.assert_array_equal(res.catt_draws_raw, res_default.catt_draws)
    rep = res.projection_report
    assert rep["att_shift_max_abs"] < 1e-10 and 0.0 <= rep["cell_draws_altered_share"] <= 1.0
    report = M.prior_data_conflict_report(res, warn=False)
    assert report["projection_applied"] and report["conflict"] is False
    res_noraw, _ = _run(panel, fit_kwargs, RankKLongBetConfig(**PROP, rank_k=2, catt_floor=0.0, keep_unprojected_draws=False))
    assert res_noraw.catt_draws_raw is None and res_noraw.projection_report is not None


def test_gaussian_summary_prior_update_uses_sqrt_shrinkage():
    rng = np.random.default_rng(0)
    draws = 1.0 + 0.5 * rng.normal(size=20000)
    out, info = M.gaussian_summary_prior_update(draws, prior_mean=0.0, prior_sd=0.5)
    # prior variance == posterior variance -> variance halves (SD factor 1/sqrt(2)), mean halves
    assert abs(info["shrinkage"] - 0.5) < 0.02 and abs(info["mean"] - 0.5) < 0.02
    assert abs(np.std(out) / np.std(draws) - np.sqrt(0.5)) < 0.01


# ---- Issue 1d: reconcile_nested_3arm_draws -----------------------------------------------------------


def test_reconcile_requires_an_explicit_prior_scale_and_returns_prior_draws(res_default):
    with pytest.raises(ValueError, match="explicit prior scale"):
        M.reconcile_nested_3arm_draws(res_default, res_default, res_default)
    with pytest.warns(M.RankKConfigWarning, match="in_sample"):
        rr = M.reconcile_nested_3arm_draws(res_default, res_default, res_default, prior_scale="in_sample", n_proposals=20000)
    a, b, c = rr  # backward-compatible unpacking
    np.testing.assert_allclose(c, a + b, atol=1e-12)
    rr2 = M.reconcile_nested_3arm_draws(
        res_default, res_default, res_default, prior_scale=1.0, return_prior_draws=True, n_proposals=20000
    )
    assert rr2.prior_catt_ac is not None and rr2.prior_catt_ac.shape == rr2.trt_catt_ac.shape
    assert rr2.sir_ess > 1.0 and rr2.lik_cov.shape == (2, 2)
    # GLS point reconciliation: the equal split of the additivity discrepancy
    m = float(np.mean(res_default.att_draws))
    np.testing.assert_allclose(rr2.lik_mean, [m - m / 3.0, m - m / 3.0], rtol=1e-10)
    for mode in ("bundle_halfstudentt", "bundle_gamma_spike", "independent_halfnormal"):
        out = M.reconcile_nested_3arm_draws(res_default, res_default, res_default, prior_scale=1.0, prior_mode=mode, n_proposals=20000)
        assert out.trt_catt_ac.shape == rr2.trt_catt_ac.shape and np.all(np.isfinite(out.trt_catt_ac))


# ---- Issue 2: short cohorts ----------------------------------------------------------------------------


def test_short_cohorts_are_excluded_from_higher_modes_exactly(panel):
    T = panel["y"].shape[1]
    coh = panel["cohorts"].copy()
    short = np.arange(30, 36)
    coh[short] = T - 1  # L_g = 1 < K = 2
    z = (np.arange(T)[None, :] >= coh[:, None]).astype(np.int32)
    y = panel["y"].copy()
    y[short, T - 1] += 0.5
    keep = np.ones(len(coh), dtype=bool)
    keep[short] = False
    lam = np.random.default_rng(3).normal(size=(len(coh), 3))
    ps = np.full(len(coh), 0.5)
    cfg = RankKLongBetConfig(**PROP, rank_k=2, use_spatial_gp=False)
    common = dict(baseline_forecast=panel["f"], baseline_forecast_draws=panel["f_draws"])
    with pytest.warns(M.RankKConfigWarning, match="drop_mode"):
        r_with = fit(y, z, panel["x"], cfg=cfg, cohorts=coh, ps=ps, pre_anchor_loadings=lam, **common)
    r_without = fit(
        y[keep], z[keep], panel["x"][keep], cfg=cfg, cohorts=coh[keep], ps=ps[keep], pre_anchor_loadings=lam[keep],
        baseline_forecast=panel["f"][keep], baseline_forecast_draws=panel["f_draws"][:, keep],
    )
    np.testing.assert_array_equal(r_with.mode_tau_draws[1][keep], r_without.mode_tau_draws[1])
    assert np.all(r_with.mode_tau_draws[1][short] == 0.0)
    assert r_with.active_units[0].all() and not r_with.active_units[1][short].any()
    assert np.all(r_with.catt_draws[:, short, : T - 1] == 0.0)
    with pytest.raises(ValueError, match="short_cohort_policy"):
        fit(y, z, panel["x"], cfg=replace(cfg, short_cohort_policy="error"), cohorts=coh, **common)
    with pytest.warns(M.RankKConfigWarning, match="reduce_k"):
        r_red = fit(y, z, panel["x"], cfg=replace(cfg, short_cohort_policy="reduce_k"), cohorts=coh, **common)
    assert len(r_red.mode_tau_draws) == 1 and all(b.shape[1] == 1 for b in r_red.bases.values())


@pytest.mark.skipif(not _has_stochtree(), reason="stochtree not installed")
def test_short_cohorts_are_excluded_in_the_plugin_engine_too(panel):
    T = panel["y"].shape[1]
    coh = panel["cohorts"].copy()
    short = np.arange(30, 36)
    coh[short] = T - 1
    z = (np.arange(T)[None, :] >= coh[:, None]).astype(np.int32)
    cfg = RankKLongBetConfig(engine="stochtree", uncertainty_propagation=False, rank_k=2, num_gfr=2, num_burnin=10, num_mcmc=20)
    with pytest.warns(M.RankKConfigWarning, match="drop_mode"):
        res = fit(panel["y"], z, panel["x"], cfg=cfg, cohorts=coh, baseline_forecast=panel["f"])
    assert np.all(res.mode_tau_draws[1][short] == 0.0) and not res.active_units[1][short].any()


def test_placebo_audit_runs_with_rank_k_one(monkeypatch, panel):
    calls = {}
    real = M._fit_rank_k_longbet if hasattr(M, "_fit_rank_k_longbet") else M.fit_rank_k_longbet

    def spy(*args, **kwargs):
        calls["cfg"] = kwargs["cfg"]
        return real(*args, **kwargs)

    monkeypatch.setattr(M, "_fit_rank_k_longbet" if hasattr(M, "_fit_rank_k_longbet") else "fit_rank_k_longbet", spy)
    hist = np.random.default_rng(1).normal(size=(panel["y"].shape[0], 8))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        M.rank_k_placebo_audit(panel["y"], panel["z"], panel["x"], y_pre_history=hist, cfg=RankKLongBetConfig(**PROP, rank_k=3, use_spatial_gp=False))
    assert calls["cfg"].rank_k == 1


# ---- Issue 3: sampler controls and the hybrid estimator -------------------------------------------------


def test_sampler_controls_run_and_are_recorded(panel, fit_kwargs, res_default):
    r_thin, _ = _run(panel, fit_kwargs, RankKLongBetConfig(**PROP, rank_k=1, thin=2))
    assert r_thin.att_draws.shape == res_default.att_draws.shape and not np.allclose(r_thin.att_draws, res_default.att_draws)
    r_eta, _ = _run(panel, fit_kwargs, RankKLongBetConfig(**PROP, rank_k=1, treated_extra_variance=True))
    s_eta2 = r_eta.diagnostics["s_eta2_draws"][0]
    assert s_eta2.shape == r_eta.att_draws.shape and np.all(s_eta2 >= 0.0) and np.std(s_eta2) > 0.0
    r_hc, _ = _run(panel, fit_kwargs, RankKLongBetConfig(**PROP, rank_k=1, learn_sigma_nu=True, sigma_nu_prior_scale=0.5))
    snu = r_hc.diagnostics["sigma_nu_total_draws"][0]
    assert np.all(np.isfinite(snu)) and np.std(snu) > 0.0
    fixed = res_default.diagnostics["sigma_nu_total_draws"][0]
    np.testing.assert_allclose(fixed, 0.25 * res_default.diagnostics["y_scale_mode"][0], rtol=1e-5)
    r_ctrl, _ = _run(
        panel, fit_kwargs,
        RankKLongBetConfig(**PROP, rank_k=2, nu_working_zeta_weight=0.5, alpha_nu=0.5, beta_nu=2.5, sigma_nu_total=0.75, sigma_nu_dynamic_total=1.0),
    )
    assert np.all(np.isfinite(r_ctrl.att_draws))
    assert len(res_default.diagnostics["rhat_mode"]) == 2 and res_default.diagnostics["rhat_mode"][0].shape == (30,)
    with pytest.raises(ValueError, match="baseline_cells"):
        RankKLongBetConfig(treated_extra_variance=True, baseline_cells="all")


@pytest.mark.skipif(not _has_stochtree(), reason="stochtree not installed")
def test_hybrid_keeps_the_propagated_att_draws(panel, fit_kwargs):
    plug = RankKLongBetConfig(
        engine="stochtree", uncertainty_propagation=False, rank_k=2, num_gfr=2, num_burnin=10, num_mcmc=40,
        num_trees_mu=20, num_trees_tau_mode1=12, num_trees_tau_higher=8, random_seed=2026,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r_h = M.fit_rank_k_hybrid(panel["y"], panel["z"], panel["x"], cfg=RankKLongBetConfig(**PROP, rank_k=2), plugin_cfg=plug, **fit_kwargs)
    r_p1, _ = _run(panel, fit_kwargs, RankKLongBetConfig(**PROP, rank_k=1))
    np.testing.assert_allclose(r_h.att_draws, r_p1.att_draws, atol=1e-10)
    assert r_h.variant == "hybrid" and len(r_h.mode_tau_draws) == 2
    assert r_h.catt_draws.shape == (r_p1.att_draws.size, *panel["y"].shape)
    with pytest.raises(ValueError, match="plug-in"):
        M.fit_rank_k_hybrid(panel["y"], panel["z"], panel["x"], cfg=RankKLongBetConfig(**PROP), plugin_cfg=replace(plug, catt_floor=0.0), **fit_kwargs)


# ---- Issue 4: exposure-basis families --------------------------------------------------------------------


@pytest.mark.parametrize("family", ["exponential", "polynomial", "cosine", "custom", "pca_empirical"])
def test_basis_families_are_orthonormal_with_constant_first_mode(family):
    spec = M.ExposureBasisSpec(
        family=family, tau_c=1.5,
        custom_generators=(lambda s, L: np.stack([s, s**2, np.sin(s)])) if family == "custom" else None,
    )
    traj = np.random.default_rng(0).normal(size=(40, 8))
    for L in (3, 6, 8):
        emp = traj[:, :L] if family == "pca_empirical" else None
        B = M.orthonormal_exposure_basis(L, 4, spec, empirical_trajectories=emp)
        assert B.shape == (L, min(4, L))
        np.testing.assert_allclose(B.T @ B / L, np.eye(B.shape[1]), atol=1e-10)
        np.testing.assert_allclose(B[:, 0], 1.0)
        for k in range(1, B.shape[1]):
            assert abs(float(np.sum(B[:, k]))) < 1e-10


def test_default_basis_spec_is_the_paper_basis_and_att_is_basis_invariant(panel, fit_kwargs, res_default):
    for L in (1, 2, 5, 8):
        np.testing.assert_array_equal(M.orthonormal_exposure_basis(L, 3), M.orthonormal_exposure_basis(L, 3, M.ExposureBasisSpec()))
    with pytest.raises(ValueError):
        M.ExposureBasisSpec(family="splines")
    for spec in (M.ExposureBasisSpec(family="polynomial"), M.ExposureBasisSpec(tau_c=1.0), M.ExposureBasisSpec(family="pca_empirical")):
        res, _ = _run(panel, fit_kwargs, RankKLongBetConfig(**PROP, rank_k=2, basis_spec=spec))
        np.testing.assert_allclose(res.att_draws, res_default.att_draws, atol=1e-12)
        np.testing.assert_array_equal(res.mode_tau_draws[0], res_default.mode_tau_draws[0])
        assert res.basis_spec == spec


def test_public_api_exports_remediation_symbols():
    for name in (
        "ExposureBasisSpec", "Nested3ArmReconciliationResult", "fit_rank_k_hybrid", "gaussian_summary_prior_update",
        "prior_data_conflict_report",
    ):
        assert name in longbet.__all__ and hasattr(longbet, name)
