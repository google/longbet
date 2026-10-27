"""Tests for ``longbet._rank_k_plots``: prior draws, prior-vs-posterior tables, sensitivity tiers, decay and figures.

The data layer is checked for shapes, additivity, the horizon-``L`` identity, quantile tables
without NaNs, ROI arithmetic, the standardized families against scipy and the prior-draw
projection; every plot function is smoke-tested under the ``Agg`` backend; and the lazy
matplotlib import is checked statically (``import longbet`` already pulls matplotlib in through
arviz, so ``sys.modules`` cannot tell).
"""

from __future__ import annotations

import ast
import os
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import stats

os.environ.setdefault("JAX_PLATFORMS", "cpu")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR / "data" / "rank_k_golden"))

import longbet  # noqa: E402
import longbet._rank_k as M  # noqa: E402
import longbet._rank_k_plots as P  # noqa: E402
from make_golden import make_panel  # noqa: E402

RankKConfigWarning = M.RankKConfigWarning
RankKLongBetConfig = M.RankKLongBetConfig
fit_rank_k_longbet = M.fit_rank_k_longbet

CFG = RankKLongBetConfig(
    uncertainty_propagation=True, n_forecast_draws=2, max_forecast_draws=2, num_chains=2, draws_per_component=20,
    num_burnin=20, num_trees_mu=20, num_trees_tau_mode1=12, num_trees_tau_higher=8, random_seed=2026,
    use_spatial_gp=False,
)
N_PROPOSALS = 20000


def _fit(kind: str):
    p = make_panel(kind)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RankKConfigWarning)
        return fit_rank_k_longbet(
            p["y"], p["z"], p["x"], cfg=CFG, baseline_forecast=p["f"], baseline_forecast_draws=p["f_draws"],
            cohorts=p["cohorts"],
        )


@pytest.fixture(scope="module")
def res_stag():
    return _fit("stag_unequal")


@pytest.fixture(scope="module")
def res_bal():
    return _fit("balanced")


@pytest.fixture(scope="module")
def tiers():
    return P.ytlr_prior_tiers(0.3)


@pytest.fixture(scope="module")
def sens(res_stag, res_bal, tiers):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RankKConfigWarning)
        return P.evaluate_prior_sensitivity(
            res_stag, res_bal, res_stag, tiers=tiers, rho=0.85, extrap_horizon=12, n_proposals=N_PROPOSALS, seed=2026
        )


@pytest.fixture(scope="module")
def synthetic_arms():
    rng = np.random.default_rng(7)
    bc = np.abs(rng.normal(0.5, 0.1, size=(80, 5)))
    ab = np.abs(rng.normal(0.3, 0.1, size=(80, 5)))
    return {"BC": bc, "AB": ab, "AC": bc + ab}


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ---- Static contract -----------------------------------------------------------------------------


def test_matplotlib_is_not_imported_at_module_level() -> None:
    tree = ast.parse(Path(P.__file__).read_text())
    top_level = []
    for node in tree.body:
        if isinstance(node, ast.Import):
            top_level.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            top_level.append(node.module or "")
    assert not any(name.split(".")[0] == "matplotlib" for name in top_level), top_level


def test_public_names_are_exported_from_the_package() -> None:
    expected = {
        "PriorTier", "RankKPriorDraws", "PriorSensitivityResult", "DecayExtrapolationResult", "ytlr_prior_tiers",
        "sample_rank_k_prior_draws", "prior_vs_posterior_summary", "assert_nested_additivity", "evaluate_prior_sensitivity",
        "extrapolate_rank_k_decay", "plot_rank_k_prior_vs_posterior", "plot_rank_k_prior_calibration",
        "plot_rank_k_prior_sensitivity", "plot_rank_k_decay_extrapolation",
    }
    assert expected <= set(P.__all__)
    for name in P.__all__:
        assert hasattr(P, name), name
    if P.__name__.startswith("longbet."):  # the package re-exports every non-constant name
        for name in P.__all__:
            if name.isupper():
                continue
            assert name in longbet.__all__, name
            assert getattr(longbet, name) is getattr(P, name)


# ---- Prior families and tiers --------------------------------------------------------------------


def test_ytlr_prior_tiers_match_the_specification(tiers) -> None:
    assert [(t.name, t.prior_mode, t.prior_scale) for t in tiers] == [
        ("Out-of-sample Half-Normal", "bundle_halfnormal", 0.3),
        ("In-sample SE Half-Normal", "bundle_halfnormal", "in_sample"),
        ("Half-Student-t", "bundle_halfstudentt", 0.3),
        ("Gamma spike", "bundle_gamma_spike", 0.3),
    ]
    assert tiers[1].is_in_sample and not tiers[0].is_in_sample
    with pytest.raises(ValueError):
        P.ytlr_prior_tiers(0.0)
    with pytest.raises(ValueError):
        P.PriorTier("x", "independent_halfnormal", 1.0)
    with pytest.raises(ValueError):
        P.PriorTier("x", "bundle_halfnormal", "posterior")


def test_standardized_families_match_scipy(tiers) -> None:
    df = P.standardized_prior_family_table(tiers)
    assert list(df["tier"]) == [t.name for t in tiers]
    assert not df.drop(columns=["scale_value"]).isna().any().any()
    hn = stats.halfnorm(scale=1.0)
    gm = stats.gamma(a=0.40, scale=1.60)
    row_hn = df.iloc[0]
    row_in = df.iloc[1]
    row_t = df.iloc[2]
    row_g = df.iloc[3]
    for q, col in ((0.3, "P30"), (0.5, "P50"), (0.7, "P70")):
        np.testing.assert_allclose(row_hn[col], hn.ppf(q), rtol=1e-12)
        np.testing.assert_allclose(row_in[col], hn.ppf(q), rtol=1e-12)
        np.testing.assert_allclose(row_t[col], 0.5 * stats.t.ppf(0.5 * (1 + q), df=2.5), rtol=1e-12)
        np.testing.assert_allclose(row_g[col], gm.ppf(q), rtol=1e-12)
    np.testing.assert_allclose(row_hn["Mean"], hn.mean(), rtol=1e-12)
    np.testing.assert_allclose(row_hn["SD"], hn.std(), rtol=1e-12)
    np.testing.assert_allclose(row_hn["Var"], hn.var(), rtol=1e-12)
    np.testing.assert_allclose(row_g["Mean"], gm.mean(), rtol=1e-12)
    np.testing.assert_allclose(row_g["SD"], gm.std(), rtol=1e-12)
    np.testing.assert_allclose(row_hn["pr_below_0.25_scale"], hn.cdf(0.25), rtol=1e-12)
    np.testing.assert_allclose(row_g["pr_below_0.5_scale"], gm.cdf(0.5), rtol=1e-12)
    np.testing.assert_allclose(row_t["pr_below_0.5_scale"], 2 * stats.t.cdf(0.5 / 0.5, df=2.5) - 1, rtol=1e-12)
    # Half-t moments: closed forms against a long Monte Carlo sample (mean only; the SD has a heavy tail).
    mc = np.abs(np.random.default_rng(0).standard_t(2.5, size=400_000)) * 0.5
    np.testing.assert_allclose(row_t["Mean"], mc.mean(), rtol=0.02)
    assert row_t["Var"] == pytest.approx(0.25 * 2.5 / 0.5 - row_t["Mean"] ** 2)
    assert row_in["scale_unit"] == "in-sample posterior SD" and np.isnan(row_in["scale_value"])
    assert row_hn["scale_unit"] == "sigma_hist" and row_hn["scale_value"] == 0.3
    with pytest.raises(ValueError, match="unique"):
        P.standardized_prior_family_table([tiers[0], tiers[0]])


# ---- Prior draws of a two-arm fit ----------------------------------------------------------------


@pytest.mark.parametrize("family", P.PRIOR_FAMILIES)
def test_prior_draw_projection_matches_posterior_trajectory_shape(res_stag, family) -> None:
    pri = P.sample_rank_k_prior_draws(res_stag, prior_family=family, prior_scale=1.0, prior_mean=0.4, seed=11)
    n_post = res_stag.att_draws.shape[0]
    t_post = res_stag.catt_draws.shape[2] - min(res_stag.bases)
    assert pri.att_draws.shape == (n_post,)
    assert pri.traj_draws.shape == (n_post, t_post)
    assert pri.posterior_traj_draws.shape == (n_post, t_post)
    assert pri.traj_draws.shape == pri.posterior_traj_draws.shape
    assert np.all(np.isfinite(pri.traj_draws)) and np.all(np.isfinite(pri.att_draws))
    assert np.all(pri.att_draws >= 0.0) and np.all(pri.traj_draws >= 0.0)  # tau_clip and catt_floor
    assert pri.family == family and pri.scale == 1.0 and pri.tau_clip == (0.0, float("inf"))
    np.testing.assert_array_equal(pri.posterior_att_draws, np.asarray(res_stag.att_draws))


def test_prior_draws_are_deterministic_and_follow_the_family(res_stag) -> None:
    a = P.sample_rank_k_prior_draws(res_stag, prior_scale=0.7, n_draws=4000, seed=3)
    b = P.sample_rank_k_prior_draws(res_stag, prior_scale=0.7, n_draws=4000, seed=3)
    np.testing.assert_array_equal(a.att_draws, b.att_draws)
    assert np.median(a.att_draws) == pytest.approx(stats.halfnorm(scale=0.7).ppf(0.5), rel=0.05)
    normal = P.sample_rank_k_prior_draws(
        res_stag, prior_family="normal", prior_mean=0.2, prior_scale=0.5, tau_clip=(-0.1, 0.9), n_draws=4000, seed=3
    )
    assert normal.att_draws.min() >= -0.1 and normal.att_draws.max() <= 0.9
    trunc = stats.truncnorm((-0.1 - 0.2) / 0.5, (0.9 - 0.2) / 0.5, loc=0.2, scale=0.5)
    assert normal.att_draws.mean() == pytest.approx(trunc.mean(), abs=0.03)
    assert np.all(np.isfinite(normal.traj_draws))


def test_prior_draw_input_validation(res_stag) -> None:
    with pytest.raises(ValueError, match="prior_family"):
        P.sample_rank_k_prior_draws(res_stag, prior_family="laplace", prior_scale=1.0)
    with pytest.raises(ValueError, match="positive"):
        P.sample_rank_k_prior_draws(res_stag, prior_scale=-1.0)
    with pytest.raises(ValueError, match="tau_clip"):
        P.sample_rank_k_prior_draws(res_stag, prior_scale=1.0, tau_clip=(-1.0, 1.0))
    with pytest.raises(ValueError, match="no mass"):
        P.sample_rank_k_prior_draws(res_stag, prior_family="normal", prior_mean=0.0, prior_scale=0.01, tau_clip=(5.0, 6.0))
    with pytest.raises(ValueError, match="n_draws"):
        P.sample_rank_k_prior_draws(res_stag, prior_scale=1.0, n_draws=1)


# ---- Prior-vs-posterior tables -------------------------------------------------------------------


def test_prior_vs_posterior_summary_columns_and_arithmetic(res_stag) -> None:
    pri = P.sample_rank_k_prior_draws(res_stag, prior_scale=1.0)
    df = P.prior_vs_posterior_summary(pri.traj_draws, pri.posterior_traj_draws, labels="ATT", ref_total=50.0, unit="$")
    assert len(df) == 1 and df.loc[0, "contrast"] == "ATT" and df.loc[0, "unit"] == "$"
    assert df.loc[0, "horizon_periods"] == pri.traj_draws.shape[1]
    expected = [
        "contrast", "horizon_periods", "unit", "prior_P30", "prior_P50", "prior_P70", "prior_Mean", "prior_SD",
        "post_P30", "post_P50", "post_P70", "post_Mean", "post_SD", "ref_total", "prior_lift_pct", "post_lift_pct",
        "width_change_pct", "interval_overlap",
    ]
    assert list(df.columns) == expected
    assert not df.drop(columns=["contrast", "unit", "interval_overlap"]).isna().any().any()
    post_total = pri.posterior_traj_draws.sum(axis=1)
    assert df.loc[0, "post_P50"] == pytest.approx(np.quantile(post_total, 0.5))
    assert df.loc[0, "post_Mean"] == pytest.approx(post_total.mean())
    assert df.loc[0, "post_SD"] == pytest.approx(post_total.std(ddof=1))
    assert df.loc[0, "post_lift_pct"] == pytest.approx(100.0 * np.quantile(post_total, 0.5) / 50.0)
    prior_width = df.loc[0, "prior_P70"] - df.loc[0, "prior_P30"]
    post_width = df.loc[0, "post_P70"] - df.loc[0, "post_P30"]
    assert df.loc[0, "width_change_pct"] == pytest.approx(100.0 * (post_width - prior_width) / prior_width)
    assert df.loc[0, "interval_overlap"] in ("PASS", "FAIL")


def test_prior_vs_posterior_summary_dict_inputs_and_overlap_assertion(synthetic_arms) -> None:
    prior = {k: v + 5.0 for k, v in synthetic_arms.items()}  # shifted prior: no overlap
    df = P.prior_vs_posterior_summary(prior, synthetic_arms, labels={"AC": "Combined"}, ref_total={"AC": 10.0})
    assert list(df["contrast"]) == ["BC", "AB", "Combined"]
    assert list(df["interval_overlap"]) == ["FAIL"] * 3
    assert np.isnan(df.loc[0, "ref_total"]) and df.loc[2, "ref_total"] == 10.0
    with pytest.raises(AssertionError, match="do not overlap"):
        P.prior_vs_posterior_summary(prior, synthetic_arms, assert_overlap=True)
    same = P.prior_vs_posterior_summary(synthetic_arms, synthetic_arms, assert_overlap=True)
    assert list(same["interval_overlap"]) == ["PASS"] * 3 and np.allclose(same["width_change_pct"], 0.0)
    with pytest.raises(ValueError, match="share their contrast keys"):
        P.prior_vs_posterior_summary({"BC": synthetic_arms["BC"]}, synthetic_arms)
    with pytest.raises(ValueError, match="both be mappings"):
        P.prior_vs_posterior_summary(synthetic_arms["BC"], synthetic_arms)
    with pytest.raises(ValueError, match="NaN"):
        bad = synthetic_arms["BC"].copy()
        bad[0, 0] = np.nan
        P.prior_vs_posterior_summary(bad, synthetic_arms["BC"])
    with pytest.raises(ValueError, match="odd length"):
        P.prior_vs_posterior_summary(synthetic_arms, synthetic_arms, quantiles=(0.3, 0.7))


def test_assert_nested_additivity(synthetic_arms) -> None:
    bc, ab, ac = synthetic_arms["BC"], synthetic_arms["AB"], synthetic_arms["AC"]
    P.assert_nested_additivity(bc, ab, ac)
    P.assert_nested_additivity(bc, ab, ac + 1e-9, atol=1e-8)
    with pytest.raises(AssertionError, match="draw by draw"):
        P.assert_nested_additivity(bc, ab, ac + np.where(np.arange(80)[:, None] == 0, 1e-3, -1e-3 / 79))
    with pytest.raises(AssertionError, match="mean"):
        P.assert_nested_additivity(bc, ab, ac + 0.1)
    # Means agree but draws are permuted: passes only without the elementwise check.
    perm = ac[np.random.default_rng(0).permutation(80)]
    P.assert_nested_additivity(bc, ab, perm, elementwise=False)
    with pytest.raises(AssertionError):
        P.assert_nested_additivity(bc, ab, perm)
    with pytest.raises(ValueError, match="same horizon"):
        P.assert_nested_additivity(bc, ab, ac[:, :3])


def test_nested_3arm_trajectories_requires_prior_draws(sens) -> None:
    rec = sens.draws_by_tier["Out-of-sample Half-Normal"]
    prior, post = P.nested_3arm_trajectories(rec)
    assert list(prior) == list(post) == ["BC", "AB", "AC"]
    assert prior["AC"].shape == post["AC"].shape == rec.trt_catt_ac.shape
    P.assert_nested_additivity(prior["BC"], prior["AB"], prior["AC"], atol=1e-10)
    P.assert_nested_additivity(post["BC"], post["AB"], post["AC"], atol=1e-10)
    df = P.prior_vs_posterior_summary(prior, post, ref_total=100.0)
    assert len(df) == 3 and not df["post_P50"].isna().any()
    assert df.loc[2, "post_Mean"] == pytest.approx(df.loc[0, "post_Mean"] + df.loc[1, "post_Mean"], abs=1e-8)
    stripped = M.Nested3ArmReconciliationResult(
        trt_catt_bc=rec.trt_catt_bc, trt_catt_ab=rec.trt_catt_ab, trt_catt_ac=rec.trt_catt_ac,
        post_att_bc=rec.post_att_bc, post_att_ab=rec.post_att_ab, post_att_ac=rec.post_att_ac,
    )
    with pytest.raises(ValueError, match="return_prior_draws"):
        P.nested_3arm_trajectories(stripped)


# ---- Prior sensitivity ---------------------------------------------------------------------------


def test_evaluate_prior_sensitivity_tables(sens, tiers) -> None:
    names = [t.name for t in tiers]
    assert list(sens.draws_by_tier) == names
    assert sens.n_observed == 8 and sens.extrap_horizon == 12 and sens.rho == 0.85
    comp = sens.posterior_comparison_df
    assert len(comp) == len(names) * 3 * 2
    assert set(comp["horizon"]) == {"observed", "extrapolated"} and set(comp["contrast"]) == {"BC", "AB", "AC"}
    assert not comp[["P30", "P50", "P70", "Mean", "SD", "prior_scale", "sir_ess"]].isna().any().any()
    assert (comp["P30"] <= comp["P50"]).all() and (comp["P50"] <= comp["P70"]).all()
    for name in names:
        rec = sens.draws_by_tier[name]
        assert rec.prior_catt_ac is not None
        obs = sens.observed_cumulative_by_tier[name]
        ext = sens.extrapolated_cumulative_by_tier[name]
        np.testing.assert_allclose(obs["AC"], obs["BC"] + obs["AB"], atol=1e-10)
        np.testing.assert_allclose(ext["AC"], ext["BC"] + ext["AB"], atol=1e-10)
        np.testing.assert_allclose(obs["AC"], rec.trt_catt_ac.sum(axis=1), atol=1e-12)
        tail = rec.trt_catt_ac[:, -1] * np.sum(0.85 ** np.arange(1, 5))
        np.testing.assert_allclose(ext["AC"], obs["AC"] + tail, atol=1e-10)
        row = comp.loc[(comp["tier"] == name) & (comp["contrast"] == "AC") & (comp["horizon"] == "extrapolated")].iloc[0]
        assert row["P50"] == pytest.approx(np.quantile(ext["AC"], 0.5)) and row["horizon_periods"] == 12
    in_sample = sens.standardized_families_df.set_index("tier").loc["In-sample SE Half-Normal"]
    assert in_sample["scale_value"] == pytest.approx(sens.draws_by_tier["In-sample SE Half-Normal"].prior_scale)
    assert sens.standardized_families_df.set_index("tier").loc["Gamma spike", "scale_value"] == 0.3


def test_roi_break_even_arithmetic(sens, tiers) -> None:
    roi = sens.roi_break_even_df
    assert len(roi) == len(tiers) * 3 * len(sens.cost_grid) * len(sens.alpha_grid)
    assert list(roi.columns) == [
        "tier", "contrast", "horizon_periods", "cost", "alpha", "pr_roi_positive", "break_even_alpha_50",
    ]
    assert roi["pr_roi_positive"].between(0.0, 1.0).all() and not roi.isna().any().any()
    for row in roi.itertuples(index=False):
        draws = sens.extrapolated_cumulative_by_tier[row.tier][row.contrast]
        assert row.pr_roi_positive == pytest.approx(np.mean(row.alpha * draws > row.cost))
        assert row.break_even_alpha_50 == pytest.approx(row.cost / np.median(draws))
    table = P.roi_break_even_table(
        {"AC": np.full(50, 4.0)}, cost_grid=(2.0,), alpha_grid=(0.25, 0.5, 0.75), horizon_periods=12
    )
    assert list(table["pr_roi_positive"]) == [0.0, 0.0, 1.0] and set(table["break_even_alpha_50"]) == {0.5}
    with pytest.raises(ValueError, match="alpha_grid"):
        P.roi_break_even_table({"AC": np.ones(10)}, cost_grid=(1.0,), alpha_grid=(1.5,), horizon_periods=1)


def test_evaluate_prior_sensitivity_silences_only_the_in_sample_warning(res_stag, res_bal) -> None:
    only_in_sample = (P.PriorTier("In-sample", "bundle_halfnormal", "in_sample"),)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        P.evaluate_prior_sensitivity(res_stag, res_bal, res_stag, tiers=only_in_sample, n_proposals=N_PROPOSALS)
    assert not any("in_sample" in str(w.message) for w in caught)
    # The filter must not leak: a direct call still warns.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        M.reconcile_nested_3arm_draws(res_stag, res_bal, res_stag, prior_scale="in_sample", n_proposals=5000)
    assert sum("in_sample" in str(w.message) for w in caught) == 1
    with pytest.raises(ValueError, match="extrap_horizon"):
        P.evaluate_prior_sensitivity(res_stag, res_bal, res_stag, tiers=only_in_sample, extrap_horizon=3, n_proposals=5000)
    with pytest.raises(ValueError, match="reconcile_kwargs"):
        P.evaluate_prior_sensitivity(res_stag, res_bal, res_stag, tiers=only_in_sample, reconcile_kwargs={"seed": 1})


# ---- Decay extrapolation -------------------------------------------------------------------------


def test_extrapolate_rank_k_decay_shapes_identities_and_certificate(synthetic_arms) -> None:
    d = P.extrapolate_rank_k_decay(synthetic_arms, horizon=24, rho=0.85, rho_sd=0.07, num_mc_per_draw=5, seed=1)
    assert d.arms == ("BC", "AB", "AC") and d.n_observed == 5 and d.horizon == 24
    for arm in d.arms:
        assert d.monthly_draws[arm].shape == (80, 24) and d.cumulative_draws[arm].shape == (80, 24)
        assert d.random_rho_cumulative_draws[arm].shape == (400, 24)
        np.testing.assert_array_equal(d.monthly_draws[arm][:, :5], synthetic_arms[arm])
        np.testing.assert_allclose(d.monthly_draws[arm][:, 5:], synthetic_arms[arm][:, -1:] * 0.85 ** np.arange(1, 20))
        np.testing.assert_allclose(d.cumulative_draws[arm][:, 4], synthetic_arms[arm].sum(axis=1), atol=1e-10)
        assert np.all(np.diff(d.cumulative_draws[arm], axis=1) >= 0.0)
        obs_total = synthetic_arms[arm].sum(axis=1)
        for rho, per_arm in d.rho_grid_cumulative.items():
            tail = synthetic_arms[arm][:, -1] * np.sum(rho ** np.arange(1, 20))
            np.testing.assert_allclose(per_arm[arm][:, -1], obs_total + tail)
    # Additivity survives the fixed-rho, grid and hierarchical extrapolations (shared decay draws).
    P.assert_nested_additivity(d.cumulative_draws["BC"], d.cumulative_draws["AB"], d.cumulative_draws["AC"], atol=1e-10)
    rand = d.random_rho_cumulative_draws
    P.assert_nested_additivity(rand["BC"], rand["AB"], rand["AC"], atol=1e-10)
    a, b = d.beta_params
    assert a / (a + b) == pytest.approx(0.85) and np.sqrt(a * b / ((a + b) ** 2 * (a + b + 1))) == pytest.approx(0.07)
    cert = d.data_health_certificate_df
    assert list(cert.columns) == ["arm", "month", "n_draws", "n_nan", "n_inf", "Min", "P25", "P50", "P95", "Max", "Mean"]
    assert list(cert.loc[cert["arm"] == "AC", "month"]) == [1, 3, 5, 12, 24]
    assert (cert["n_nan"] == 0).all() and (cert["n_inf"] == 0).all() and (cert["n_draws"] == 80).all()
    assert (cert["Min"] <= cert["P25"]).all() and (cert["P25"] <= cert["P50"]).all() and (cert["P95"] <= cert["Max"]).all()
    mile = d.milestone_table_df
    assert list(mile.columns) == [
        "arm", "month", "segment", "rate_P30", "rate_P50", "rate_P70", "rate_Mean", "cum_P30", "cum_P50", "cum_P70", "cum_Mean",
    ]
    assert list(mile.loc[mile["arm"] == "BC", "segment"]) == ["observed"] * 3 + ["extrapolated"] * 2
    row = mile.loc[(mile["arm"] == "AC") & (mile["month"] == 12)].iloc[0]
    assert row["cum_P50"] == pytest.approx(np.quantile(d.cumulative_draws["AC"][:, 11], 0.5))
    assert not mile.drop(columns=["arm", "segment"]).isna().any().any()
    cum = d.cumulative_summary_df
    assert len(cum) == 3 * 24 and not cum.isna().any().any()
    band = cum.loc[cum["arm"] == "AC"].sort_values("month")
    np.testing.assert_allclose(band["rho_band_P50"].to_numpy(), np.quantile(d.random_rho_cumulative_draws["AC"], 0.5, axis=0))
    assert list(d.monthly_summary_df.columns) == ["arm", "month", "segment", "P30", "P50", "P70", "Mean"]


def test_extrapolate_rank_k_decay_single_arm_and_options(synthetic_arms) -> None:
    arr = synthetic_arms["AC"]
    d = P.extrapolate_rank_k_decay(arr, horizon=5, rho_grid=(0.5,), milestones=(1, 3, 60))
    assert d.arms == ("effect",) and d.rho_sd is None and d.random_rho_cumulative_draws is None and d.beta_params is None
    np.testing.assert_array_equal(d.monthly_draws["effect"], arr)
    assert list(d.milestone_table_df["month"]) == [1, 3, 5]
    assert d.cumulative_summary_df["rho_band_P50"].isna().all()
    same = P.extrapolate_rank_k_decay(arr, horizon=12, rho_sd=0.07, num_mc_per_draw=3, seed=5)
    again = P.extrapolate_rank_k_decay(arr, horizon=12, rho_sd=0.07, num_mc_per_draw=3, seed=5)
    np.testing.assert_array_equal(
        same.random_rho_cumulative_draws["effect"], again.random_rho_cumulative_draws["effect"]
    )
    bad_options = (
        dict(horizon=3), dict(rho=1.5), dict(rho_sd=0.5), dict(num_mc_per_draw=0), dict(milestones=(0,)),
        dict(rho_grid=(2.0,)),
    )
    for bad in bad_options:
        with pytest.raises(ValueError):
            P.extrapolate_rank_k_decay(arr, **bad)
    with pytest.raises(ValueError, match="same"):
        P.extrapolate_rank_k_decay({"a": arr, "b": arr[:, :3]})
    with pytest.raises(ValueError, match="NaN"):
        P.extrapolate_rank_k_decay(np.full((10, 3), np.nan))


# ---- Figures (smoke tests under Agg) -------------------------------------------------------------


def _is_axes(obj) -> bool:
    return isinstance(obj, matplotlib.axes.Axes)


def test_plot_prior_vs_posterior_variants(res_stag, synthetic_arms) -> None:
    pri = P.sample_rank_k_prior_draws(res_stag, prior_scale=1.0)
    both = P.plot_rank_k_prior_vs_posterior(
        pri.traj_draws, pri.posterior_traj_draws, labels="ATT", unit="$", ref_totals=50.0, title="t"
    )
    assert isinstance(both, np.ndarray) and both.shape == (2,) and all(_is_axes(a) for a in both)
    shifted = {k: v + 0.1 for k, v in synthetic_arms.items()}
    forest = P.plot_rank_k_prior_vs_posterior(shifted, synthetic_arms, kind="forest")
    assert _is_axes(forest) and len(forest.get_yticklabels()) == 3
    _, ax = plt.subplots()
    density = P.plot_rank_k_prior_vs_posterior(
        synthetic_arms["AC"], synthetic_arms["AC"], kind="density", ax=ax, title="d"
    )
    assert density is ax and ax.get_title() == "d"
    with pytest.raises(ValueError, match="kind"):
        P.plot_rank_k_prior_vs_posterior(synthetic_arms, synthetic_arms, kind="violin")
    with pytest.raises(ValueError, match="expected 2 Axes"):
        P.plot_rank_k_prior_vs_posterior(synthetic_arms, synthetic_arms, ax=ax)


def test_plot_prior_calibration_variants() -> None:
    axes = P.plot_rank_k_prior_calibration(0.3, hist_mean=3.0, hist_se=0.5, horizon_periods=8, unit="$", title="c")
    assert isinstance(axes, np.ndarray) and axes.shape == (2,) and all(_is_axes(a) for a in axes)
    labels = [t.get_text() for t in axes[0].get_legend().get_texts()]
    assert any("Historical benchmark" in s for s in labels) and any("Implied single-component" in s for s in labels)
    single = P.plot_rank_k_prior_calibration(0.3, prior_mode="bundle_gamma_spike", show_rayleigh_comparison=False)
    assert _is_axes(single)
    for mode in P.BUNDLE_PRIOR_MODES:
        assert _is_axes(P.plot_rank_k_prior_calibration(1.0, prior_mode=mode, show_rayleigh_comparison=False, hist_mean=0.5))
    with pytest.raises(ValueError, match="hist_mean"):
        P.plot_rank_k_prior_calibration(0.3, hist_se=0.1)
    with pytest.raises(ValueError, match="prior_mode"):
        P.plot_rank_k_prior_calibration(0.3, prior_mode="independent_halfnormal")
    # Closed-form density of the sum of two iid half-normals integrates to one.
    grid = np.linspace(0.0, 40.0, 200_001)
    assert np.trapezoid(P._two_halfnormal_sum_pdf(grid, 2.0), grid) == pytest.approx(1.0, abs=1e-6)


def test_plot_prior_sensitivity(sens) -> None:
    axes = P.plot_rank_k_prior_sensitivity(sens, metric_label="revenue", unit="$", cost_highlight=1.10, title="s")
    assert isinstance(axes, np.ndarray) and axes.shape == (3,) and all(_is_axes(a) for a in axes)
    assert len(axes[1].get_yticklabels()) == len(sens.tiers)
    assert len(axes[2].get_lines()) >= len(sens.tiers)
    _, given = plt.subplots(1, 3)
    out = P.plot_rank_k_prior_sensitivity(sens, contrast="AB", axes=given)
    assert all(o is g for o, g in zip(out, given))
    with pytest.raises(ValueError, match="contrast"):
        P.plot_rank_k_prior_sensitivity(sens, contrast="CA")
    with pytest.raises(ValueError, match="expected 3 Axes"):
        P.plot_rank_k_prior_sensitivity(sens, axes=given[:2])


def test_plot_decay_extrapolation(sens) -> None:
    prior, post = P.nested_3arm_trajectories(sens.draws_by_tier["Half-Student-t"])
    d_post = P.extrapolate_rank_k_decay(post, horizon=36, num_mc_per_draw=4)
    d_prior = P.extrapolate_rank_k_decay(prior, horizon=36, num_mc_per_draw=4)
    axes = P.plot_rank_k_decay_extrapolation(d_post, prior_decay_result=d_prior, unit="$", title="decay")
    assert isinstance(axes, np.ndarray) and axes.shape == (2,) and all(_is_axes(a) for a in axes)
    no_band = P.extrapolate_rank_k_decay(post["AC"], horizon=12, rho_sd=None)
    _, given = plt.subplots(1, 2)
    out = P.plot_rank_k_decay_extrapolation(no_band, axes=given, milestones=(3, 6))
    assert all(o is g for o, g in zip(out, given))
    with pytest.raises(ValueError, match="share the arms"):
        P.plot_rank_k_decay_extrapolation(d_post, prior_decay_result=no_band)
