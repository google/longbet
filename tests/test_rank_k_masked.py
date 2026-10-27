"""Tests of the masked event-time Rank-K path (0.7.0).

(a) bit-identity of complete single-wave (and staggered) panels against the reference
    generated from the unmodified 0.6.0 sources (``tests/data/make_masked_reference.py``),
    with ``mask=None`` and an all-``True`` ``mask``, both engines;
(b) validation errors of masked panels;
(c) unit tests of the event-time projection (known coordinates recovered under partial
    support, correction exactly 0 under full support, weights equal to the paper's);
(d) recovery of the population 12-period effect on a geometric-decay DGP with administrative
    censoring (fast smoke here, 100 replications under ``-m slow``);
(e) masked 3-arm reconciliation additivity;
(f) size and power of the transport audit.
"""

from __future__ import annotations

import os
from pathlib import Path
import platform
import sys
import warnings

import numpy as np
import pytest

os.environ.setdefault("JAX_PLATFORMS", "cpu")
THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR / "data"))

import longbet  # noqa: E402
import longbet._rank_k as M  # noqa: E402
import longbet._rank_k_event_time as ET  # noqa: E402
from longbet import (  # noqa: E402
    RankKLongBetConfig,
    cumulative_effect,
    event_time_trajectory,
    fit_rank_k_longbet,
    rank_k_transport_audit,
    reconcile_nested_3arm_draws,
)
from make_masked_reference import cases as reference_cases, make_panel, result_arrays  # noqa: E402

REFERENCE = THIS_DIR / "data" / "masked_reference_6f714d5.npz"

# Small propagated budget for the masked fits (JAX engine).
FAST = dict(
    uncertainty_propagation=True, use_forecast=False, n_forecast_draws=3, num_chains=2, draws_per_component=12,
    num_burnin=25, num_trees_mu=20, num_trees_tau_mode1=12, num_trees_tau_higher=8, use_spatial_gp=False,
    use_factor_loadings=False, random_seed=2026,
)


def _has_stochtree() -> bool:
    try:
        import stochtree  # noqa: F401
    except ImportError:
        return False
    return True


def _quiet_fit(y, z, x, cfg, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fit_rank_k_longbet(y, z, x, cfg=cfg, **kw)


def censored_dgp(seed, *, N=160, t0=5, H=15, rho=0.85, beta0=1.0, beta1=0.5, sigma=1.0, short_amp=1.0, p_short=0.4):
    """Event-time panel with geometric decay, unit / period effects and administrative censoring.

    ``L_i`` is 5..7 for a share ``p_short`` of the units and 8..15 otherwise; the short
    units' effect amplitude is multiplied by ``short_amp`` (1: no cohort-by-event-time
    interaction). Returns the panel and the population 12-period truth.
    """
    rng = np.random.default_rng(seed)
    T = t0 + H
    x = rng.normal(size=(N, 3))
    treated = np.zeros(N, dtype=bool)
    treated[: N // 2] = True
    rng.shuffle(treated)
    z = np.zeros((N, T), dtype=np.int32)
    z[treated, t0:] = 1
    alpha = rng.normal(scale=0.7, size=(N, 1))
    gamma = np.cumsum(rng.normal(scale=0.3, size=T))[None, :]
    s = np.arange(1, H + 1)
    L = np.where(rng.random(N) < p_short, rng.integers(5, 8, size=N), rng.integers(8, H + 1, size=N))
    # two units per arm observe the whole horizon, so every post period has both arms observed
    L[np.nonzero(treated)[0][:2]] = H
    L[np.nonzero(~treated)[0][:2]] = H
    beta = beta0 + beta1 * x[:, 0]
    amp = np.where(L <= 7, short_amp, 1.0)
    tau = (beta * amp)[:, None] * rho ** (s - 1)[None, :]
    y = alpha + gamma + 0.5 * x[:, :1] + sigma * rng.normal(size=(N, T))
    y[:, t0:] += z[:, t0:] * tau
    for i in range(N):
        y[i, t0 + L[i]:] = np.nan
    truth_12 = float(np.mean(tau[treated, :12].sum(axis=1)))
    return dict(y=y, z=z, x=x, L=L, treated=treated, tau=tau, truth_12=truth_12, t0=t0, H=H)


# ---- (a) Bit-identity on complete panels -------------------------------------------------------


@pytest.mark.parametrize("name,kind,cfg_kw,kw_builder", list(reference_cases()), ids=[c[0] for c in reference_cases()])
@pytest.mark.parametrize("mask_mode", ["none", "all_true"])
def test_complete_panels_are_bit_identical_to_0_6_0(name, kind, cfg_kw, kw_builder, mask_mode):
    if cfg_kw.get("uncertainty_propagation") is False and not _has_stochtree():
        pytest.skip("stochtree not installed")
    if not REFERENCE.exists():
        pytest.skip(f"reference {REFERENCE.name} not present (generate it from the 0.6.0 sources)")
    ref = np.load(REFERENCE)
    import jax

    same_env = str(ref["_version"]) == "0.6.0" and jax.__version__ == "0.11.2" and platform.machine() == "x86_64"
    panel = make_panel(kind)
    kw = kw_builder(panel)
    if mask_mode == "all_true":
        kw["mask"] = np.ones(panel["y"].shape, dtype=bool)
    res = _quiet_fit(panel["y"], panel["z"], panel["x"], RankKLongBetConfig(**cfg_kw), **kw)
    for key, val in result_arrays(res).items():
        r = ref[f"{name}/{key}"]
        if same_env:
            np.testing.assert_array_equal(r, np.asarray(val), err_msg=f"{name}:{key} ({mask_mode}) is not bit-identical")
        else:  # another JAX build / CPU: floating-point reordering noise only
            np.testing.assert_allclose(r, np.asarray(val), rtol=1e-5, atol=1e-6, err_msg=f"{name}:{key}")
    # The event-time result fields are populated on single-wave panels and consistent.
    if kind == "single":
        assert res.horizon == panel["y"].shape[1] - int(panel["cohorts"].min())
        assert res.catt_observed.sum() == int(np.sum(panel["z"] == 1))
        assert np.all(res.support_mask)
        np.testing.assert_allclose(res.unit_weights.sum(), 1.0)
    else:
        assert res.horizon is None and res.catt_observed is None


def test_unit_weighting_equals_cell_weighting_on_a_balanced_panel():
    panel = make_panel("single")
    cfg = RankKLongBetConfig(**dict(FAST, rank_k=2))
    res_cell = _quiet_fit(panel["y"], panel["z"], panel["x"], cfg)
    res_unit = _quiet_fit(panel["y"], panel["z"], panel["x"], RankKLongBetConfig(**dict(FAST, rank_k=2, att_weighting="unit")))
    np.testing.assert_array_equal(res_cell.att_draws, res_unit.att_draws)
    np.testing.assert_array_equal(res_cell.unit_weights, res_unit.unit_weights)
    assert res_unit.att_weighting == "unit" and res_cell.att_weighting == "cell"
    # the aggregates reduce to the ATT on a complete balanced panel
    np.testing.assert_allclose(cumulative_effect(res_unit, res_unit.horizon), res_unit.horizon * res_unit.att_draws, atol=1e-10)
    traj = event_time_trajectory(res_unit)
    assert traj.shape == (res_unit.att_draws.shape[0], res_unit.horizon)
    np.testing.assert_allclose(traj.mean(axis=1), res_unit.att_draws, atol=1e-10)


# ---- (b) Validation ----------------------------------------------------------------------------


def _small_masked_panel(seed=3, N=40, t0=4, H=8):
    d = censored_dgp(seed, N=N, t0=t0, H=H)
    return d


def test_validation_errors_of_masked_panels():
    d = _small_masked_panel()
    y, z, x = d["y"], d["z"], d["x"]
    cfg = RankKLongBetConfig(**dict(FAST, rank_k=2))
    validate = M._validate_rank_k_inputs
    # inf is never allowed
    y_inf = y.copy()
    y_inf[0, 0] = np.inf
    with pytest.raises(ValueError, match="infinite"):
        validate(y_inf, z, x, cfg=cfg)
    # mask must agree with the NaN pattern
    with pytest.raises(ValueError, match="mask and NaN pattern"):
        validate(y, z, x, cfg=cfg, mask=np.ones_like(y, dtype=bool))
    with pytest.raises(ValueError, match="mask must be a boolean"):
        validate(y, z, x, cfg=cfg, mask=np.ones((3, 3), dtype=bool))
    # a gap inside the observed post window
    i = int(np.nonzero(d["treated"] & (d["L"] >= 6))[0][0])
    y_gap = y.copy()
    y_gap[i, d["t0"] + 2] = np.nan
    with pytest.raises(ValueError, match="gap"):
        validate(y_gap, z, x, cfg=cfg)
    # no observed pre-period cell
    y_pre = y.copy()
    y_pre[i, : d["t0"]] = np.nan
    with pytest.raises(ValueError, match="no observed pre-period cell"):
        validate(y_pre, z, x, cfg=cfg)
    # fewer post periods than min_support_modes[0]
    with pytest.raises(ValueError, match="min_support_modes\\[0\\]"):
        validate(y, z, x, cfg=RankKLongBetConfig(**dict(FAST, rank_k=2, min_support_modes=(6, 6))))
    # a period with observed treated but no observed never-treated outcome
    y_ctrl = y.copy()
    y_ctrl[~d["treated"], d["t0"] + 5:] = np.nan
    with pytest.raises(ValueError, match="no observed never-treated outcome"):
        validate(y_ctrl, z, x, cfg=cfg)
    # masked panels need a single adoption period
    z_stag = z.copy()
    late = np.nonzero(d["treated"])[0][:5]
    z_stag[late, d["t0"]] = 0
    with pytest.raises(ValueError, match="single adoption period"):
        validate(y, z_stag, x, cfg=cfg)
    # horizon range and engine support
    with pytest.raises(ValueError, match="horizon must be"):
        validate(y, z, x, cfg=cfg, horizon=d["H"] + 1)
    with pytest.raises(ValueError, match="horizon must be"):
        validate(y, z, x, cfg=cfg, horizon=0)
    if _has_stochtree():
        with pytest.raises(NotImplementedError, match="JAX engine only"):
            validate(y, z, x, cfg=RankKLongBetConfig(uncertainty_propagation=False))
    with pytest.raises(ValueError, match="common_basis=False"):
        validate(y, z, x, cfg=RankKLongBetConfig(**dict(FAST, rank_k=2, common_basis=False)))
    # a valid masked panel returns its structure
    info = validate(y, z, x, cfg=cfg, warn=False)
    assert info["masked"] and info["common_basis"] and info["horizon"] == int(d["L"][d["treated"]].max())
    # the placebo audit and the hybrid estimator refuse masked panels
    with pytest.raises(ValueError, match="masked panels"):
        longbet.rank_k_placebo_audit(y, z, x, cfg=cfg)
    with pytest.raises(NotImplementedError, match="masked panels"):
        longbet.fit_rank_k_hybrid(y, z, x, cfg=cfg)


def test_config_validation_of_the_new_fields():
    with pytest.raises(ValueError, match="att_weighting"):
        RankKLongBetConfig(att_weighting="units")
    with pytest.raises(ValueError, match="projection"):
        RankKLongBetConfig(projection="joint")
    with pytest.raises(ValueError, match="mode_coupling"):
        RankKLongBetConfig(mode_coupling="draws")
    with pytest.raises(ValueError, match="min_support_modes"):
        RankKLongBetConfig(min_support_modes=(1, 8, 6))
    with pytest.raises(ValueError, match="min_support_modes"):
        RankKLongBetConfig(min_support_modes=())
    with pytest.raises(ValueError, match="min_support_modes"):
        RankKLongBetConfig(min_support_modes=(0, 6))
    cfg = RankKLongBetConfig(min_support_modes=[1.0, 6, 8.0], common_basis=1)
    assert cfg.min_support_modes == (1, 6, 8) and cfg.common_basis is True


def test_event_time_options_warn_on_the_per_cohort_path():
    panel = make_panel("staggered")
    cfg = RankKLongBetConfig(**dict(FAST, rank_k=2, min_support_modes=(1, 2)))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        M._validate_rank_k_inputs(panel["y"], panel["z"], panel["x"], cfg=cfg)
    msgs = [str(w.message) for w in caught if issubclass(w.category, longbet.RankKConfigWarning)]
    assert any("min_support_modes" in m and "per-cohort" in m for m in msgs)


# ---- (c) Projection unit tests -----------------------------------------------------------------


def _layout(obs, cohorts, t0, K, min_support=(1, 4, 6)):
    return ET.build_event_time_layout(
        obs, cohorts, t0, horizon=None, rank_k=K, min_support=min_support,
        basis_factory=lambda h, k: M.orthonormal_exposure_basis(h, k),
    )


def test_projection_matches_the_paper_on_complete_single_wave_panels():
    rng = np.random.default_rng(1)
    N, T, t0, K = 50, 12, 4, 3
    cohorts = np.full(N, T, dtype=np.int32)
    cohorts[:28] = t0
    ctrl = cohorts == T
    dy_dm = rng.normal(size=(N, T))
    stack = dy_dm[None] + 0.1 * rng.normal(size=(3, N, T))
    bases = {t0: M.orthonormal_exposure_basis(T - t0, K)}
    ww = {t0: 1.0}
    layout = _layout(np.ones((N, T), bool), cohorts, t0, K)
    assert not layout.masked and not layout.has_partial_windows()
    for k in range(K):
        np.testing.assert_array_equal(
            M._project_panel_to_mode(dy_dm, cohorts, ctrl, [t0], ww, bases, k), ET.project_event_time_mode(dy_dm, layout, k)
        )
        np.testing.assert_array_equal(
            M._project_panel_to_mode(stack, cohorts, ctrl, [t0], ww, bases, k), ET.project_event_time_mode(stack, layout, k)
        )
        np.testing.assert_array_equal(
            M._compute_mode_precision_weights(cohorts, ctrl, [t0], ww, bases, t0, T, k), ET.event_time_precision_weights(layout, k)
        )
    corr = ET.mode1_dynamic_correction(layout, [rng.normal(size=N), rng.normal(size=N)])
    assert np.all(corr == 0.0)


def test_two_pass_and_restricted_wls_recover_known_coordinates_under_partial_support():
    rng = np.random.default_rng(2)
    N, T, t0, K = 60, 14, 4, 3
    H = T - t0
    cohorts = np.full(N, T, dtype=np.int32)
    cohorts[:36] = t0
    ever = cohorts < T
    Psi = M.orthonormal_exposure_basis(H, K)
    L = rng.integers(4, H + 1, size=N)
    L[:4] = H
    obs = np.ones((N, T), dtype=bool)
    for i in range(N):
        obs[i, t0 + L[i]:] = False
    layout = _layout(obs, cohorts, t0, K)
    assert layout.horizon == H and layout.has_partial_windows()
    sup = layout.support_matrix()
    assert np.array_equal(sup[:, 1], L >= 4) and np.array_equal(sup[:, 2], L >= 6)

    def panel_from(theta):
        r = np.zeros((N, T))
        r[:, t0:] = theta @ Psi.T
        return np.where(obs, r, np.nan)

    # Higher modes: a level plus one dynamic mode is recovered exactly on every supported unit,
    # whatever the window (the window constant is partialled out on partial windows).
    for k in (1, 2):
        theta = np.zeros((N, K))
        theta[:, 0] = rng.normal(size=N)
        theta[:, k] = rng.normal(size=N)
        proj = ET.project_event_time_mode(panel_from(theta), layout, k)
        np.testing.assert_allclose(proj[sup[:, k]], theta[sup[:, k], k], atol=1e-10)
        assert np.all(proj[~sup[:, k]] == 0.0)
        w = ET.event_time_precision_weights(layout, k)
        assert np.all(w[~sup[:, k]] == 0.0) and np.all(w[sup[:, k]] > 0.0)
    # Without the centering the level would leak into the mode-2 projection on partial windows.
    theta = np.zeros((N, K))
    theta[:, 0] = 1.0
    r_nan = panel_from(theta)
    plain = np.array([np.nansum(r_nan[i, t0:] * np.where(obs[i, t0:], Psi[:, 1], 0.0)) for i in range(N)])
    assert np.max(np.abs(plain[ever & (L < H)])) > 0.1
    assert np.max(np.abs(ET.project_event_time_mode(r_nan, layout, 1))) < 1e-10

    # Mode 1: the uncorrected window mean is biased on partial windows ...
    theta = rng.normal(size=(N, K))
    r_nan = panel_from(theta)
    base = ET.project_event_time_mode(r_nan, layout, 0)
    partial_trt = ever & (L < H)
    assert np.max(np.abs(base[partial_trt] - theta[partial_trt, 0])) > 0.1
    # ... and exact after the correction with the true higher-mode coordinates (treated units)
    corr = ET.mode1_dynamic_correction(layout, [theta[:, 1], theta[:, 2]])
    np.testing.assert_allclose((base - corr)[ever], theta[ever, 0], atol=1e-10)
    assert np.all(corr[~ever] == 0.0) and np.all(corr[L == H] == 0.0)
    # restricted WLS: exact on every supported coordinate of the treated units (joint fit)
    th, vf, fell_back = ET.restricted_wls_projection(r_nan, layout, [theta[:, 1], theta[:, 2]])
    assert not np.any(fell_back)
    for k in range(K):
        sel = sup[:, k] & ever
        np.testing.assert_allclose(th[k][sel], theta[sel, k], atol=1e-9)
        assert np.all(np.isinf(vf[k][~sup[:, k]]))
    full = (L == H)
    np.testing.assert_allclose(vf[0][full], 1.0 / H)
    w0 = ET.restricted_wls_precision_weights(layout, vf, 0)
    assert np.all(w0 > 0.0)
    # weights: a shorter window gets a smaller mode-1 weight
    w1 = ET.event_time_precision_weights(layout, 0)
    trt = np.nonzero(ever)[0]
    order = np.argsort(L[trt])
    assert w1[trt][order[0]] < w1[trt][order[-1]]


def test_observation_structure_checks():
    N, T, t0 = 10, 9, 3
    cohorts = np.full(N, T, dtype=np.int32)
    cohorts[:5] = t0
    obs = np.ones((N, T), dtype=bool)
    info = ET.check_observation_structure(obs, cohorts, t0, horizon=None, min_support=(1, 6, 8))
    assert info["horizon"] == T - t0 and np.all(info["post_len"] == T - t0) and np.all(info["pre_count"] == t0)
    obs[0, 7:] = False
    obs[1, 5:] = False
    info = ET.check_observation_structure(obs, cohorts, t0, horizon=None, min_support=(1, 6, 8))
    assert info["post_len"][0] == 4 and info["post_len"][1] == 2
    info = ET.check_observation_structure(obs, cohorts, t0, horizon=3, min_support=(1, 6, 8))
    assert info["horizon"] == 3 and info["post_len"][0] == 3 and info["post_len"][1] == 2
    lengths, is_prefix = ET.observed_prefix_lengths(np.array([[True, True, False], [True, False, True]]))
    assert list(lengths) == [2, 2] and list(is_prefix) == [True, False]


# ---- (d) Recovery on a censored geometric-decay DGP ---------------------------------------------


def _recovery_stats(n_rep, *, K=2, coupling="propagated", seed0=1000, N=160):
    errs, covered, rel = [], [], []
    cfg = RankKLongBetConfig(**dict(FAST, rank_k=K, att_weighting="unit", mode_coupling=coupling, n_forecast_draws=4,
                                     draws_per_component=25, num_burnin=50, num_trees_mu=30, num_trees_tau_mode1=20,
                                     num_trees_tau_higher=10))
    for r in range(n_rep):
        d = censored_dgp(seed0 + r, N=N)
        res = _quiet_fit(d["y"], d["z"], np.column_stack([d["x"], d["L"]]), cfg)
        cum = cumulative_effect(res, 12)
        lo, hi = np.percentile(cum, [10, 90])
        errs.append(float(np.mean(cum)) - d["truth_12"])
        covered.append(bool(lo <= d["truth_12"] <= hi))
        rel.append(d["truth_12"])
    return float(np.mean(errs)), float(np.mean(covered)), float(np.mean(rel)), res


def test_masked_fit_recovers_the_population_12_period_effect_smoke():
    bias, coverage, truth, res = _recovery_stats(4)
    # loose smoke thresholds (4 replications); the calibrated check is the slow test below
    assert abs(bias) < 0.25 * truth
    assert coverage >= 0.5
    # structural checks of the result
    assert res.horizon == 15 and res.catt_draws.shape[2] == 20
    assert res.catt_observed.sum() < int(np.sum(res.catt_draws[0] != 0.0))  # predicted cells exist
    assert res.support_mask.shape == (res.catt_draws.shape[1], 2)
    assert res.variant == "propagated" and res.n_components == 4
    assert res.diagnostics["mode_coupling"] == "propagated"
    assert res.active_units[1].sum() == res.support_mask[:, 1].sum() < res.active_units[0].sum()
    assert np.all(np.isfinite(res.att_draws)) and res.att_se > 0


def test_mode_coupling_removes_the_short_window_bias():
    d = censored_dgp(7, N=160)
    x = np.column_stack([d["x"], d["L"]])
    cfg_none = RankKLongBetConfig(**dict(FAST, rank_k=2, att_weighting="unit", mode_coupling="none"))
    cfg_plug = RankKLongBetConfig(**dict(FAST, rank_k=2, att_weighting="unit", mode_coupling="plug_in"))
    res_none = _quiet_fit(d["y"], d["z"], x, cfg_none)
    res_plug = _quiet_fit(d["y"], d["z"], x, cfg_plug)
    err_none = abs(float(np.mean(cumulative_effect(res_none, 12))) - d["truth_12"])
    err_plug = abs(float(np.mean(cumulative_effect(res_plug, 12))) - d["truth_12"])
    assert err_plug < err_none
    assert res_none.diagnostics["mode_coupling"] == "none" and res_plug.diagnostics["mode_coupling"] == "plug_in"
    assert res_plug.variant == "jax" and res_plug.n_components == 1


def test_restricted_wls_runs_and_keeps_the_identity_sum_tau_equals_h_tau1():
    d = censored_dgp(11, N=160)
    cfg = RankKLongBetConfig(**dict(FAST, rank_k=3, att_weighting="unit", projection="restricted_wls", mode_coupling="plug_in"))
    res = _quiet_fit(d["y"], d["z"], d["x"], cfg)
    assert res.diagnostics["projection"] == "restricted_wls"
    trt = np.nonzero(res.treated_units)[0]
    t0 = int(min(res.bases))
    cum_h = res.catt_draws[:, trt, t0 : t0 + res.horizon].sum(axis=2)  # (S, N_trt)
    np.testing.assert_allclose(cum_h, res.horizon * res.mode_tau_draws[0][trt].T, atol=1e-8)


def test_forecast_path_on_a_masked_panel_with_complete_pre_period():
    d = censored_dgp(5, N=160)
    rng = np.random.default_rng(0)
    f = np.nan_to_num(d["y"], nan=0.0) * 0.0 + rng.normal(scale=0.3, size=d["y"].shape).cumsum(axis=1)
    f_draws = f[None] + 0.1 * rng.normal(size=(3, *f.shape))
    cfg = RankKLongBetConfig(**dict(FAST, use_forecast=True, rank_k=2, att_weighting="unit"))
    res = _quiet_fit(d["y"], d["z"], d["x"], cfg, baseline_forecast=f, baseline_forecast_draws=f_draws)
    assert res.n_components == 3 and res.variant == "propagated"
    assert np.all(np.isfinite(res.att_draws)) and np.all(np.isfinite(res.catt_draws))
    # heterogeneous pre coverage with a bare point forecast: the internal draw generator needs
    # complete outcomes, so the fit conditions on the point forecast and warns
    y2 = d["y"].copy()
    y2[:30, :2] = np.nan
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res2 = fit_rank_k_longbet(y2, d["z"], d["x"], cfg=cfg, baseline_forecast=f)
    msgs = [str(w.message) for w in caught if issubclass(w.category, longbet.RankKConfigWarning)]
    assert any("internal forecast-draw generator" in m for m in msgs)
    assert res2.diagnostics["n_forecast_components"] == 1 and np.all(np.isfinite(res2.att_draws))
    assert res2.diagnostics["n_mode_components"] == 3 and res2.n_components == 3  # higher-mode tracks only
    # ... while user-supplied draws are propagated as usual (paired with the higher-mode draws)
    res2b = _quiet_fit(y2, d["z"], d["x"], cfg, baseline_forecast=f, baseline_forecast_draws=f_draws)
    assert res2b.diagnostics["n_forecast_components"] == 3 and res2b.n_components == 3
    assert np.all(np.isfinite(res2b.att_draws))
    # explicit horizon shorter than the panel
    res3 = _quiet_fit(d["y"], d["z"], d["x"], RankKLongBetConfig(**dict(FAST, rank_k=2, att_weighting="unit")), horizon=6)
    assert res3.horizon == 6 and res3.bases[d["t0"]].shape == (6, 2)
    assert np.all(res3.catt_draws[:, :, d["t0"] + 6:] == 0.0)
    np.testing.assert_allclose(cumulative_effect(res3, 6), 6 * res3.att_draws, atol=1e-8)


@pytest.mark.slow
@pytest.mark.parametrize("coupling", ["propagated", "plug_in"])
def test_masked_fit_recovers_the_population_12_period_effect(coupling):
    bias, coverage, truth, _ = _recovery_stats(100, coupling=coupling)
    assert abs(bias) < 0.10 * truth, (bias, truth)
    assert 0.72 <= coverage <= 0.88, coverage


# ---- (e) Masked 3-arm reconciliation --------------------------------------------------------------


def test_masked_three_arm_reconciliation_is_additive():
    rng = np.random.default_rng(21)
    N, t0, H = 150, 4, 10
    T = t0 + H
    arm = np.repeat(np.arange(3), N // 3)  # 0 = C, 1 = B, 2 = A; equal arms so every pair has 100 units
    rng.shuffle(arm)
    x = rng.normal(size=(N, 2))
    L = np.where(rng.random(N) < 0.4, rng.integers(5, 8, size=N), rng.integers(8, H + 1, size=N))
    s = np.arange(1, H + 1)
    y = rng.normal(scale=0.5, size=(N, 1)) + rng.normal(size=(N, T))
    eff_b = 0.8 * 0.9 ** (s - 1)
    eff_a = 1.3 * 0.9 ** (s - 1)
    y[arm == 1, t0:] += eff_b
    y[arm == 2, t0:] += eff_a
    for i in range(N):
        y[i, t0 + L[i]:] = np.nan
    cfg = RankKLongBetConfig(**dict(FAST, rank_k=2, att_weighting="unit", mode_coupling="plug_in"))

    def pair(trt_arm, ctrl_arm):
        keep = (arm == trt_arm) | (arm == ctrl_arm)
        z = np.zeros((int(keep.sum()), T), dtype=np.int32)
        z[arm[keep] == trt_arm, t0:] = 1
        return _quiet_fit(y[keep], z, np.column_stack([x[keep], L[keep]]), cfg)

    res_bc, res_ab, res_ac = pair(1, 0), pair(2, 1), pair(2, 0)
    for r in (res_bc, res_ab, res_ac):
        assert r.att_weighting == "unit" and r.horizon == H and r.catt_observed is not None
    out = reconcile_nested_3arm_draws(res_bc, res_ab, res_ac, prior_scale=1.0, att_w_bc=res_bc.unit_weights, att_w_ab=res_ab.unit_weights)
    bc, ab, ac = out
    assert bc.shape == ab.shape == ac.shape == (res_bc.att_draws.shape[0], H)
    np.testing.assert_array_equal(ac, bc + ab)
    assert np.all(np.isfinite(ac))
    # the stored unit weights are picked up without passing them explicitly
    out2 = reconcile_nested_3arm_draws(res_bc, res_ab, res_ac, prior_scale=1.0)
    np.testing.assert_array_equal(out2.trt_catt_ac, ac)


# ---- (f) Transport audit --------------------------------------------------------------------------


def _audit(seed, short_amp, *, N=160):
    d = censored_dgp(seed, N=N, short_amp=short_amp)
    cfg = RankKLongBetConfig(**dict(FAST, rank_k=2, att_weighting="unit", mode_coupling="plug_in"))
    res = _quiet_fit(d["y"], d["z"], np.column_stack([d["x"], d["L"]]), cfg)
    strata = np.where(d["L"] >= 12, "long", np.where(d["L"] <= 7, "short", None))
    return rank_k_transport_audit(res, strata, overlap_horizon=5, cumulative_horizon=12), res, strata


def test_transport_audit_reports_and_detects_an_injected_interaction():
    null, res, strata = _audit(31, 1.0)
    assert null.labels == ("long", "short") and null.overlap_horizon == 5
    assert set(null.quantiles) == {"p5", "p30", "p50", "p70", "p95"}
    assert null.quantiles["p5"] <= null.quantiles["p50"] <= null.quantiles["p95"]
    assert null.diff_draws.shape == (res.att_draws.shape[0], 5)
    assert null.unobserved_cell_share["long"] == 0.0 and 0.0 < null.unobserved_cell_share["short"] < 1.0
    assert 0.0 < null.unobserved_effect_share["short"] < 1.0
    assert "Transport audit" in null.summary()
    # explicit pair and threshold
    same = rank_k_transport_audit(res, strata, overlap_horizon=5, pair=("short", "long"), threshold=100.0)
    np.testing.assert_allclose(same.mean_diff_draws, -null.mean_diff_draws)
    assert not same.flag
    with pytest.raises(ValueError, match="exactly two labels"):
        rank_k_transport_audit(res, np.where(strata == "long", "long", None))
    with pytest.raises(ValueError, match="overlap_horizon"):
        rank_k_transport_audit(res, strata, overlap_horizon=99)
    # a strong cohort-by-event-time interaction (short-window units with triple the effect) is flagged
    alt, _, _ = _audit(31, 3.0)
    assert alt.flag and alt.quantiles["p5"] > 0.0


@pytest.mark.slow
def test_transport_audit_size_and_power():
    flags_null = [_audit(100 + r, 1.0)[0].flag for r in range(40)]
    flags_alt = [_audit(200 + r, 2.5)[0].flag for r in range(40)]
    assert np.mean(flags_null) <= 0.25, np.mean(flags_null)   # nominal 10% for a 90% interval
    assert np.mean(flags_alt) >= 0.80, np.mean(flags_alt)


def test_masked_pre_period_causal_features_neutralizes_unobserved_cells():
    """Unobserved pre-anchor cells in build_causal_features do not inject harmonic-imputation artifacts into pre_anchor_loadings."""
    from longbet import build_causal_features

    d = censored_dgp(77, N=60)
    y_masked_pre = d["y"].copy()
    # Mask first 8 pre-period columns for half the units (late-entry cohorts)
    late = np.arange(60) >= 30
    y_masked_pre[late, :8] = np.nan
    feats = build_causal_features(y_masked_pre, d["z"], x=d["x"], n_factor_loadings=2, seed=77)
    assert feats.pre_anchor_loadings is not None
    assert feats.pre_anchor_loadings.shape == (60, 2)
    assert np.all(np.isfinite(feats.pre_anchor_loadings))
    # Verify complete-panel bit-identity is preserved
    feats_complete = build_causal_features(d["y"], d["z"], x=d["x"], n_factor_loadings=2, seed=77)
    assert np.all(np.isfinite(feats_complete.pre_anchor_loadings))

