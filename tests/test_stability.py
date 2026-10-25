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

"""Tests for ``longbet._stability``: wear-out ANCOVA, retention, slope correlation, pooling, the rule and the report.

Exact identities (HC1 ANCOVA with and without covariates, the ``Var_R`` reductions and a hand ``Var_R`` with unequal
variances, DerSimonian-Laird against a hand computation, Morris shrinkage at ``tau^2 = 0``, NLS recovery and boundary
pile-up), the adversarial cap-censored outcome, degenerate windows under both policies, every rule clause, malformed
inputs and determinism all run on synthetic draws without MCMC; a single integration smoke test fits the golden
Rank-K fixture with a tiny configuration.
"""

from __future__ import annotations

import ast
import dataclasses
import os
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import stats

os.environ.setdefault("JAX_PLATFORMS", "cpu")

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR / "data" / "rank_k_golden"))

import longbet  # noqa: E402
import longbet._aggregates as A  # noqa: E402
import longbet._rank_k as M  # noqa: E402
import longbet._stability as S  # noqa: E402
from make_golden import make_panel  # noqa: E402

GRID_STEP = 0.0005
WINDOW = (4, 10)
Z_MDD = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)
ANCHOR, OBSERVED, HORIZON, PULSE_TAU = 7, 11, 48, 2.5  # pseudo-horizon L0, observed L, horizon H, pulse decay


def _paths(rho, *, amplitude=1.0, n_periods=10, first=4, pulse=2.0) -> np.ndarray:
    """``(S, L)`` trajectories: a pulse level before ``first``, then ``amplitude * rho^(m - first)``."""
    rho = np.atleast_1d(np.asarray(rho, dtype=np.float64))
    amp = np.broadcast_to(np.asarray(amplitude, dtype=np.float64), rho.shape)
    draws = np.full((rho.size, n_periods), float(pulse))
    steps = np.arange(n_periods - first + 1, dtype=np.float64)
    draws[:, first - 1 :] = amp[:, None] * rho[:, None] ** steps[None, :]
    return draws


def _exact_sample(mean: float, sd: float, n: int = 400, seed: int = 0) -> np.ndarray:
    """Draws whose sample mean and SD (ddof 1) are exactly ``mean`` and ``sd``."""
    z = np.random.default_rng(seed).standard_normal(n)
    z = (z - z.mean()) / z.std(ddof=1)
    return mean + sd * z


def _quiet(func, *args, **kwargs):
    """Call ``func`` with the module's StabilityWarnings silenced."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", S.StabilityWarning)
        return func(*args, **kwargs)


def _wearout_panel(*, seed: int = 3, n: int = 400, effect: float = 1.0, decay: float = 0.0, horizon: int = 12):
    """Two-arm panel: effect ``effect`` in periods 1-6 and ``effect - decay`` after; 3 strata; 25% censored units."""
    rng = np.random.default_rng(seed)
    treated = rng.random(n) < 0.5
    strata = rng.integers(1, 4, size=n)
    baseline = rng.gamma(2.0, 1.0, size=n)
    periods = np.arange(1, horizon + 1)
    effect_path = np.where(periods <= 6, effect, effect - decay)
    y = (
        0.8 * baseline[:, None] + 0.3 * strata[:, None] + treated[:, None] * effect_path[None, :]
        + rng.normal(scale=1.0, size=(n, horizon))
    )
    lengths = np.where(rng.random(n) < 0.25, rng.integers(6, horizon, size=n), horizon)
    y[np.arange(horizon)[None, :] >= lengths[:, None]] = np.nan
    return {"y": y, "treated": treated, "baseline": baseline, "strata": strata, "lengths": lengths}


def _hand_hc1(outcome: np.ndarray, design: np.ndarray) -> tuple[float, float]:
    """Treatment coefficient and HC1 SE by the normal equations (independent of the library's lstsq path)."""
    xtx = design.T @ design
    beta = np.linalg.solve(xtx, design.T @ outcome)
    resid = outcome - design @ beta
    bread = np.linalg.inv(xtx)
    meat = design.T @ (design * resid[:, None] ** 2)
    n_obs, n_params = design.shape
    cov = bread @ meat @ bread * n_obs / (n_obs - n_params)
    return float(beta[1]), float(np.sqrt(cov[1, 1]))


# ---- Static contract -----------------------------------------------------------------------------


def test_matplotlib_is_not_imported_at_module_level() -> None:
    tree = ast.parse(Path(S.__file__).read_text())
    imported = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.append(node.module or "")
    assert not any(name.split(".")[0] == "matplotlib" for name in imported), imported


def test_public_names_are_documented_and_no_package_module_imports_stability() -> None:
    for name in S.__all__:
        assert hasattr(S, name), name
        obj = getattr(S, name)
        if not name.isupper():
            assert name in longbet.__all__ and getattr(longbet, name) is obj, name
        if callable(obj):
            doc = (obj.__doc__ or "").strip()
            assert doc and not doc.startswith(f"{name}("), name  # dataclasses auto-generate "Name(field, ...)"
    package_dir = Path(S.__file__).parent
    for path in package_dir.glob("*.py"):
        if path.name in ("__init__.py", "_stability.py"):
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom):
                modules = [node.module or ""] + [f"{node.module}.{alias.name}" for alias in node.names]
            elif isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            else:
                continue
            assert not any(module == "longbet._stability" for module in modules), path.name


# ---- Wear-out ANCOVA -------------------------------------------------------------------------


def test_wearout_ancova_matches_a_hand_ols_with_hc1_on_the_balanced_subpanel() -> None:
    panel = _wearout_panel()
    result = S.wearout_ancova(panel["y"], panel["treated"], panel["baseline"], strata=panel["strata"])
    balanced = panel["lengths"] >= 12
    assert result.n_units == int(balanced.sum()) and result.min_post_periods == 12
    assert result.n_treated == int(panel["treated"][balanced].sum())
    y_bal = panel["y"][balanced]
    base = panel["baseline"][balanced]
    strata = panel["strata"][balanced]
    design = np.column_stack([
        np.ones(balanced.sum()), panel["treated"][balanced].astype(float), base - base.mean(),
        (strata == 2).astype(float), (strata == 3).astype(float),
    ])
    early, late = y_bal[:, :6].mean(axis=1), y_bal[:, 6:12].mean(axis=1)
    expected = {
        "level_early": early,
        "level_late": late,
        "change_late_minus_early": late - early,
        "change_late_split": y_bal[:, 9:12].mean(axis=1) - y_bal[:, 6:9].mean(axis=1),
    }
    summary = result.summary_df.set_index("measure")
    assert list(summary.index) == list(expected)
    for measure, values in expected.items():
        est, se = _hand_hc1(values, design)
        row = summary.loc[measure]
        assert row["est"] == pytest.approx(est, abs=1e-10) and row["se_hc1"] == pytest.approx(se, rel=1e-10)
        assert row["t_stat"] == pytest.approx(est / se, rel=1e-10)
        assert row["p_value"] == pytest.approx(2.0 * stats.norm.sf(abs(est / se)), rel=1e-10)
        assert row["mdd"] == pytest.approx(Z_MDD * se, rel=1e-12)
        assert row["mdd_pct_of_early"] == pytest.approx(100.0 * Z_MDD * se / abs(summary.loc["level_early", "est"]))
    assert result.mdd_multiplier == pytest.approx(2.8016, abs=1e-4)
    change = summary.loc["change_late_minus_early"]
    assert result.change_mdd == change["mdd"] and result.early_estimate == summary.loc["level_early", "est"]
    assert result.powered_for_half_wearout == bool(change["mdd"] <= 0.5 * max(result.early_estimate, 0.0))
    assert summary.loc["level_early", "powered_for_half_wearout"] is None
    # Passing the window lengths explicitly selects the same balanced subpanel as the NaN prefix.
    explicit = S.wearout_ancova(
        panel["y"], panel["treated"].astype(int), panel["baseline"], strata=panel["strata"],
        post_window_length=panel["lengths"],
    )
    pd.testing.assert_frame_equal(explicit.summary_df, result.summary_df)
    cert = result.data_health_certificate_df
    assert list(cert.columns) == ["kind", "quantity", "n", "n_nan", "n_inf", "Min", "P25", "P50", "P95", "Max", "Mean"]
    window_row = cert.loc[cert["quantity"].str.startswith("y_event periods")].iloc[0]
    assert window_row["n_nan"] == 0 and window_row["n"] == result.n_units * 12
    assert set(cert["kind"]) == {"row_count", "distribution"}


def test_wearout_ancova_power_flags_and_covariates() -> None:
    strong = S.wearout_ancova(**_wearout_kwargs(_wearout_panel(n=4000, effect=1.0, decay=1.0, seed=5)))
    assert strong.powered_for_half_wearout
    assert strong.change_estimate == pytest.approx(-1.0, abs=0.15) and strong.mdd_pct_of_early < 50.0
    weak = S.wearout_ancova(**_wearout_kwargs(_wearout_panel(n=60, effect=0.2, seed=6)))
    assert not weak.powered_for_half_wearout and weak.mdd_pct_of_early > 50.0
    panel = _wearout_panel(seed=7)
    covariate = np.random.default_rng(1).normal(size=panel["y"].shape[0])
    with_x = S.wearout_ancova(**_wearout_kwargs(panel), x=covariate)
    without_split = S.wearout_ancova(**_wearout_kwargs(panel), late_split=None)
    assert list(without_split.summary_df["measure"]) == ["level_early", "level_late", "change_late_minus_early"]
    assert with_x.summary_df.loc[0, "est"] != pytest.approx(without_split.summary_df.loc[0, "est"], abs=1e-12)


def test_wearout_ancova_covariates_match_a_hand_hc1_and_ignore_nan_outside_the_subpanel() -> None:
    panel = _wearout_panel(seed=7)
    balanced = panel["lengths"] >= 12
    covariates = np.random.default_rng(1).normal(size=(panel["y"].shape[0], 2))
    covariates[~balanced] = np.nan  # unobserved outside the balanced subpanel: harmless
    baseline = np.where(balanced, panel["baseline"], np.nan)
    result = S.wearout_ancova(panel["y"], panel["treated"], baseline, strata=panel["strata"], x=covariates)
    y_bal, base, strata = panel["y"][balanced], panel["baseline"][balanced], panel["strata"][balanced]
    design = np.column_stack([
        np.ones(balanced.sum()), panel["treated"][balanced].astype(float), base - base.mean(),
        (strata == 2).astype(float), (strata == 3).astype(float), covariates[balanced],
    ])
    early, late = y_bal[:, :6].mean(axis=1), y_bal[:, 6:12].mean(axis=1)
    expected = {
        "level_early": early,
        "level_late": late,
        "change_late_minus_early": late - early,
        "change_late_split": y_bal[:, 9:12].mean(axis=1) - y_bal[:, 6:9].mean(axis=1),
    }
    summary = result.summary_df.set_index("measure")
    for measure, values in expected.items():
        est, se = _hand_hc1(values, design)
        assert summary.loc[measure, "est"] == pytest.approx(est, abs=1e-10)
        assert summary.loc[measure, "se_hc1"] == pytest.approx(se, rel=1e-10)
    cert = result.data_health_certificate_df.set_index("quantity")
    assert cert.loc["baseline (balanced)", "n_nan"] == 0 and cert.loc["baseline (balanced)", "n"] == balanced.sum()
    # Mixed-type object strata fall back cleanly to repr-sorted dummies, and differential attrition emits a
    # StabilityWarning.
    mixed_strata = np.array([1 if s == 1 else f"s{s}" for s in panel["strata"]], dtype=object)
    mixed = S.wearout_ancova(panel["y"], panel["treated"], panel["baseline"], strata=mixed_strata)
    pd.testing.assert_frame_equal(mixed.summary_df, S.wearout_ancova(**_wearout_kwargs(panel)).summary_df)
    attrit_lengths = np.where(panel["treated"] & (np.arange(panel["y"].shape[0]) < 150), 6, panel["lengths"])
    with pytest.warns(S.StabilityWarning, match="treated share in the balanced subpanel"):
        attrit = S.wearout_ancova(
            panel["y"], panel["treated"], panel["baseline"], strata=panel["strata"], post_window_length=attrit_lengths
        )
    assert len(attrit.warnings) == 1 and "administrative censoring" in attrit.warnings[0]


def _wearout_kwargs(panel: dict) -> dict:
    return {
        "y_event": panel["y"], "treated": panel["treated"], "baseline": panel["baseline"], "strata": panel["strata"],
    }


def test_wearout_ancova_rejects_malformed_inputs() -> None:
    panel = _wearout_panel()
    good = _wearout_kwargs(panel)
    n_units = panel["y"].shape[0]
    first_balanced = int(np.flatnonzero(panel["lengths"] >= 12)[0])
    nan_baseline = panel["baseline"].copy()
    nan_baseline[first_balanced] = np.nan
    nan_x = np.random.default_rng(2).normal(size=n_units)
    nan_x[first_balanced] = np.nan
    infinite = panel["y"].copy()
    infinite[0, 0] = np.inf
    lying_lengths = np.full(n_units, 12)
    cases = [
        dict(good, y_event=panel["y"][:, 0]),
        dict(good, y_event=infinite),
        dict(good, y_event=np.ones_like(panel["y"])),  # every measure is constant
        dict(good, treated=np.full(n_units, 2)),
        dict(good, baseline=nan_baseline),  # NaN inside the balanced subpanel
        dict(good, x=nan_x),
        dict(good, x=np.r_[np.inf, np.zeros(n_units - 1)]),
        dict(good, strata=np.where(np.arange(n_units) == 0, np.nan, 1.0)),
        dict(good, min_post_periods=11),
        dict(good, late=(7, 13)),
        dict(good, early=(0, 6)),
        dict(good, late_split=((7, 9),)),
        dict(good, alpha=1.5),
        dict(good, power=0.0),
        dict(good, treated=np.zeros(n_units, dtype=bool)),
        dict(good, x=np.column_stack([panel["baseline"], panel["baseline"]])),
        dict(good, post_window_length=lying_lengths),
    ]
    for case in cases:
        with pytest.raises(ValueError):
            S.wearout_ancova(**case)
    with pytest.raises(ValueError, match="baseline of the balanced subpanel"):
        S.wearout_ancova(**dict(good, baseline=nan_baseline))
    with pytest.raises(ValueError, match="covariates of the balanced subpanel"):
        S.wearout_ancova(**dict(good, x=nan_x))
    with pytest.raises(ValueError, match="constant on the balanced subpanel.*early / late windows"):
        S.wearout_ancova(**dict(good, y_event=np.ones_like(panel["y"])))
    split_constant = panel["y"].copy()
    split_constant[:, 9:12] = split_constant[:, 6:9]  # identical late sub-windows: a zero late-split change
    with pytest.raises(ValueError, match="late_split=None"):
        S.wearout_ancova(**dict(good, y_event=split_constant))
    without_split = S.wearout_ancova(**dict(good, y_event=split_constant), late_split=None)
    assert "change_late_split" not in set(without_split.summary_df["measure"])


# ---- Post-pulse retention --------------------------------------------------------------------


def test_retention_recovers_rho_from_noiseless_geometric_paths() -> None:
    on_grid = np.array([0.80, 0.85, 0.9, 0.95, 0.9995, 0.6125])
    off_grid = np.array([0.91234, 0.77777, 0.98765])
    amplitude = np.linspace(0.5, 3.0, on_grid.size + off_grid.size)
    result = S.estimate_post_pulse_retention(_paths(np.r_[on_grid, off_grid], amplitude=amplitude), window=WINDOW)
    np.testing.assert_allclose(result.rho_draws[: on_grid.size], on_grid, atol=1e-12)
    assert np.all(np.abs(result.rho_draws[on_grid.size :] - off_grid) <= GRID_STEP)
    np.testing.assert_allclose(result.lambda_draws, -np.log(result.rho_draws))
    np.testing.assert_allclose(result.rho_unconstrained_draws, result.rho_draws)
    assert result.share_at_upper == 0.0 and result.share_degenerate == 0.0 and not result.boundary_pileup_flag
    assert result.grid_step == pytest.approx(GRID_STEP) and result.unconstrained_upper == pytest.approx(1.5)
    assert result.support_source == "uniform" and result.trajectory_source == "draws"
    row = result.summary_df.iloc[0]
    assert row["rho_P10"] == pytest.approx(np.quantile(result.rho_draws, 0.10))
    assert row["n_identified"] == on_grid.size + off_grid.size
    assert list(result.support_df.columns) == ["event_time", "support", "in_window", "weight", "weight_share"]
    assert result.support_df["weight_share"].sum() == pytest.approx(1.0)
    # Periods before the window (the pulse) do not enter the fit.
    shifted = _paths(on_grid, pulse=50.0)
    np.testing.assert_allclose(S.estimate_post_pulse_retention(shifted, window=WINDOW).rho_draws, on_grid, atol=1e-12)


def test_flat_and_growing_windows_pile_up_on_the_rho_bounds() -> None:
    flat = S.estimate_post_pulse_retention(_paths(np.ones(40), amplitude=np.linspace(0.1, 2.0, 40)), window=WINDOW)
    assert np.all(flat.rho_draws == 1.0) and flat.share_at_upper == 1.0 and flat.boundary_pileup_flag
    np.testing.assert_allclose(flat.rho_unconstrained_draws, 1.0)
    assert flat.share_unconstrained_above_one == 0.0
    growing = S.estimate_post_pulse_retention(_paths(np.full(10, 1.1)), window=WINDOW)
    assert np.all(growing.rho_draws == 1.0) and growing.share_at_upper == 1.0
    assert np.all(np.abs(growing.rho_unconstrained_draws - 1.1) <= GRID_STEP)
    assert growing.share_unconstrained_above_one == 1.0
    fast = S.estimate_post_pulse_retention(_paths(np.full(10, 0.3)), window=WINDOW)
    assert np.all(fast.rho_draws == 0.5) and fast.share_at_lower == 1.0 and fast.share_at_upper == 0.0
    mixed = S.estimate_post_pulse_retention(_paths(np.r_[np.ones(30), np.full(70, 0.9)]), window=WINDOW)
    assert mixed.share_at_upper == pytest.approx(0.30) and mixed.boundary_pileup_flag
    strict = S.estimate_post_pulse_retention(
        _paths(np.r_[np.ones(30), np.full(70, 0.9)]), window=WINDOW, pileup_threshold=0.5
    )
    assert not strict.boundary_pileup_flag


def test_retention_support_weights_and_input_types() -> None:
    rho = np.array([0.8, 0.9, 0.95])
    draws = _paths(rho)
    corrupted = draws.copy()
    corrupted[:, -1] = 7.0
    full = np.r_[np.ones(9), 0.0]
    by_full = S.estimate_post_pulse_retention(corrupted, window=WINDOW, support=full)
    by_window = S.estimate_post_pulse_retention(corrupted, window=WINDOW, support=full[WINDOW[0] - 1 :])
    np.testing.assert_allclose(by_full.rho_draws, rho, atol=1e-12)
    np.testing.assert_array_equal(by_window.rho_draws, by_full.rho_draws)
    assert by_full.support_source == "user" and by_full.weights.tolist() == [1.0] * 6 + [0.0]
    nested = M.Nested3ArmReconciliationResult(
        trt_catt_bc=draws * 0.25, trt_catt_ab=draws * 0.75, trt_catt_ac=draws, post_att_bc=rho, post_att_ab=rho,
        post_att_ac=rho,
    )
    from_nested = S.estimate_post_pulse_retention(nested, window=WINDOW)
    from_mapping = S.estimate_post_pulse_retention({"AC": draws, "BC": draws}, window=WINDOW)
    np.testing.assert_allclose(from_nested.rho_draws, rho, atol=1e-12)
    np.testing.assert_array_equal(from_mapping.rho_draws, from_nested.rho_draws)
    assert from_nested.trajectory_source.endswith("trt_catt_ac")
    bc = S.estimate_post_pulse_retention(nested, window=WINDOW, contrast="BC")
    np.testing.assert_allclose(bc.rho_draws, rho, atol=1e-12)


def test_grid_search_does_not_depend_on_the_block_size(monkeypatch: pytest.MonkeyPatch) -> None:
    rng = np.random.default_rng(17)
    draws = _paths(rng.uniform(0.6, 1.05, size=152), amplitude=rng.uniform(0.2, 2.0, size=152))
    draws[:5, WINDOW[0] - 1 :] = 0.0  # degenerate draws must be flagged identically
    single = S.estimate_post_pulse_retention(draws, window=WINDOW)  # one block of 152 draws
    # Blocks of 7 draws on the 1001-point grid and 3 on the 2001-point grid: 152 = 21 * 7 + 5 = 50 * 3 + 2, so
    # both grids end with a partial block.
    monkeypatch.setattr(S, "_BLOCK_CELLS", 7 * 1001 + 5)
    blocked = S.estimate_post_pulse_retention(draws, window=WINDOW)
    np.testing.assert_array_equal(blocked.rho_draws, single.rho_draws)
    np.testing.assert_array_equal(blocked.rho_unconstrained_draws, single.rho_unconstrained_draws)
    np.testing.assert_array_equal(blocked.degenerate, single.degenerate)
    assert single.degenerate[:5].all() and not single.degenerate[5:].any()
    assert single.share_at_upper > 0.0 and single.share_unconstrained_above_one > 0.0


def test_degenerate_windows_under_both_policies() -> None:
    rho = np.tile([0.8, 0.9, 1.0, 0.95], 10)
    draws = _paths(rho)
    draws[:8, WINDOW[0] - 1 :] = 0.0  # identically zero window
    draws[8:10, WINDOW[0] - 1 :] = -0.5  # no positive projection: amplitude 0 for every rho
    nan_policy = S.estimate_post_pulse_retention(draws, window=WINDOW)
    assert nan_policy.share_degenerate == pytest.approx(10 / 40) and nan_policy.degenerate[:10].all()
    assert np.isnan(nan_policy.rho_draws[:10]).all() and np.isnan(nan_policy.lambda_draws[:10]).all()
    assert np.isnan(nan_policy.rho_unconstrained_draws[:10]).all()
    np.testing.assert_allclose(nan_policy.rho_draws[10:], rho[10:], atol=1e-12)
    assert nan_policy.summary_df.iloc[0]["rho_P10"] == pytest.approx(np.quantile(rho[10:], 0.10))
    assert nan_policy.share_at_upper == pytest.approx(np.mean(rho[10:] == 1.0) * 30 / 40)
    cert = nan_policy.data_health_certificate_df.set_index("quantity")
    assert cert.loc["rho", "n_nan"] == 10 and cert.loc["degenerate draws", "n"] == 10
    resampled = S.estimate_post_pulse_retention(draws, window=WINDOW, degenerate_policy="resample", seed=3)
    again = S.estimate_post_pulse_retention(draws, window=WINDOW, degenerate_policy="resample", seed=3)
    other = S.estimate_post_pulse_retention(draws, window=WINDOW, degenerate_policy="resample", seed=4)
    np.testing.assert_array_equal(resampled.rho_draws, again.rho_draws)
    assert set(resampled.rho_draws[:10]) <= set(rho[10:]) and np.isfinite(resampled.lambda_draws).all()
    np.testing.assert_array_equal(resampled.rho_draws[10:], nan_policy.rho_draws[10:])
    assert resampled.share_degenerate == nan_policy.share_degenerate
    assert not np.array_equal(other.rho_draws[:10], resampled.rho_draws[:10])
    zero = np.zeros((5, 10))
    all_nan = S.estimate_post_pulse_retention(zero, window=WINDOW)
    assert all_nan.share_degenerate == 1.0 and np.isnan(all_nan.rho_draws).all()
    assert np.isnan(all_nan.summary_df.iloc[0]["rho_P50"])
    with pytest.raises(ValueError, match="nothing to resample"):
        S.estimate_post_pulse_retention(zero, window=WINDOW, degenerate_policy="resample")


def test_retention_rejects_malformed_inputs() -> None:
    draws = _paths([0.9, 0.8, 0.95])
    with_nan = draws.copy()
    with_nan[0, 5] = np.nan
    nested = M.Nested3ArmReconciliationResult(
        trt_catt_bc=draws, trt_catt_ab=draws, trt_catt_ac=draws, post_att_bc=draws[:, 0], post_att_ab=draws[:, 0],
        post_att_ac=draws[:, 0],
    )
    cases = [
        dict(traj=with_nan),
        dict(traj=draws[0]),
        dict(traj={"BC": draws}),
        dict(traj=nested, contrast="CA"),
        dict(window=(4, 4)),
        dict(window=(4, 11)),
        dict(window=(0, 10)),
        dict(window=4),
        dict(rho_bounds=(0.0, 1.0)),
        dict(rho_bounds=(0.9, 0.8)),
        dict(unconstrained_upper=0.9),
        dict(grid_points=2),
        dict(degenerate_policy="drop"),
        dict(pileup_threshold=0.0),
        dict(pileup_threshold=1.5),
        dict(support=-np.ones(7)),
        dict(support=np.zeros(10)),
        dict(support=np.ones(5)),
        dict(seed=-1),
        dict(quantiles=(0.1, 0.9)),
    ]
    for case in cases:
        kwargs = {"traj": draws, "window": WINDOW, **case}
        with pytest.raises(ValueError):
            S.estimate_post_pulse_retention(kwargs.pop("traj"), **kwargs)


# ---- Unit slope correlation ------------------------------------------------------------------


def _slope_panel(seed: int = 9):
    """Noiseless linear windows (periods 3-9) with correlated unit slopes, an arm shift, a third arm and censoring."""
    rng = np.random.default_rng(seed)
    corr = np.array([[1.0, 0.6, 0.3], [0.6, 1.0, 0.5], [0.3, 0.5, 1.0]])
    arms = np.array(["A"] * 300 + ["C"] * 300 + ["B"] * 100, dtype=object)
    slopes = rng.multivariate_normal(np.zeros(3), corr, size=arms.size) + np.where(arms == "A", 2.0, 0.0)[:, None]
    periods = np.arange(3, 10, dtype=np.float64)
    observed = np.where(rng.random(arms.size) < 0.2, rng.integers(3, 7, size=arms.size), 7)
    windows = {}
    for j, name in enumerate(("y1", "y2", "y3")):
        level = rng.normal(size=arms.size)[:, None]
        window = level + slopes[:, [j]] * (periods[None, :] - 6.0)
        window[np.arange(7)[None, :] >= observed[:, None]] = np.nan
        windows[name] = window
    return windows, arms, slopes, observed, periods, corr


def test_unit_slope_correlation_pools_within_arm_correlations() -> None:
    windows, arms, slopes, observed, periods, corr = _slope_panel()
    result = S.unit_slope_correlation(windows, arms, periods=periods)
    eligible = np.isin(arms, ["A", "C"]) & (observed >= 5)
    np.testing.assert_array_equal(result.eligible, eligible)
    np.testing.assert_allclose(result.slopes, slopes[eligible], atol=1e-10)
    assert result.n_units == 700 and result.n_in_arms == 600 and result.n_eligible == int(eligible.sum())
    assert result.n_eligible_by_arm == {arm: int(np.sum(eligible & (arms == arm))) for arm in ("A", "C")}
    z_sum, weight_sum = np.zeros((3, 3)), 0.0
    for arm in ("A", "C"):
        arm_slopes = slopes[eligible & (arms == arm)]
        matrix = np.corrcoef(arm_slopes, rowvar=False)
        np.testing.assert_allclose(result.per_arm_matrix[arm], matrix, atol=1e-10)
        assert result.per_arm_r_hat[arm] == pytest.approx(matrix[np.triu_indices(3, 1)].mean(), abs=1e-10)
        z_sum += (arm_slopes.shape[0] - 3) * np.arctanh(np.clip(matrix, -1 + 1e-12, 1 - 1e-12))
        weight_sum += arm_slopes.shape[0] - 3
    expected = np.tanh(z_sum / weight_sum)
    np.fill_diagonal(expected, 1.0)
    np.testing.assert_allclose(result.correlation_matrix, expected, atol=1e-10)
    np.testing.assert_allclose(result.correlation_matrix, result.correlation_matrix.T)
    assert np.all(np.diag(result.correlation_matrix) == 1.0)
    assert result.r_hat == pytest.approx(expected[np.triu_indices(3, 1)].mean())
    # The arm shift of the slopes does not leak into the within-arm correlation.
    assert np.max(np.abs(result.correlation_matrix - corr)) < 0.12
    assert list(result.summary_df.columns) == ["outcome_a", "outcome_b", "r_pooled", "r_A", "r_C"]
    assert len(result.summary_df) == 3
    # On noiseless linear windows the delta slope is 3.5x the OLS slope, so every correlation agrees.
    delta = S.unit_slope_correlation(windows, arms, periods=periods, slope_method="delta", delta_early_periods=3)
    full = observed == 7
    np.testing.assert_allclose(delta.slopes[full[eligible]], 3.5 * slopes[eligible & full], atol=1e-10)
    spearman = S.unit_slope_correlation(windows, arms, periods=periods, correlation="spearman")
    rank_corr = np.corrcoef(stats.rankdata(slopes[eligible & (arms == "C")], axis=0), rowvar=False)
    np.testing.assert_allclose(spearman.per_arm_matrix["C"], rank_corr, atol=1e-10)
    cert = result.data_health_certificate_df.set_index("quantity")
    assert cert.loc["units with finite slopes (eligible)", "n"] == result.n_eligible
    assert cert.loc["slope (ols): y1", "Mean"] == pytest.approx(slopes[eligible, 0].mean())


def test_unit_slope_correlation_drops_nonfinite_delta_slopes() -> None:
    windows, arms, _, _, periods, _ = _slope_panel(seed=10)
    gappy = {name: window.copy() for name, window in windows.items()}
    for window in gappy.values():
        window[:5, :3] = np.nan  # five units observe the late half only
    result = S.unit_slope_correlation(gappy, arms, periods=periods, slope_method="delta", min_obs=4)
    assert result.n_dropped_nonfinite >= 1 and not result.eligible[:5].any()
    assert np.isfinite(result.slopes).all()


def test_unit_slope_correlation_rejects_malformed_inputs() -> None:
    windows, arms, _, _, periods, _ = _slope_panel()
    constant = {name: np.tile(np.arange(7.0), (arms.size, 1)) for name in windows}
    infinite = {**windows, "y1": np.where(np.arange(7) == 0, np.inf, windows["y1"])}
    cases = [
        dict(y_event_by_outcome={"y1": windows["y1"]}),
        dict(y_event_by_outcome={**windows, "y2": windows["y2"][:, :5]}),
        dict(y_event_by_outcome=infinite),
        dict(slope_method="median"),
        dict(correlation="kendall"),
        dict(arms_used=("A", "Z")),
        dict(arms_used=("A", "A")),
        dict(arms=arms[:-1]),
        dict(min_obs=1),
        dict(min_obs=8),
        dict(periods=periods[::-1]),
        dict(delta_early_periods=7),
        dict(y_event_by_outcome=constant),
        dict(arms=np.where(np.arange(arms.size) < 3, "C", "A").astype(object), arms_used=("A", "C")),
    ]
    for case in cases:
        kwargs = {"y_event_by_outcome": windows, "arms": arms, "periods": periods, **case}
        with pytest.raises(ValueError):
            S.unit_slope_correlation(kwargs.pop("y_event_by_outcome"), kwargs.pop("arms"), **kwargs)


# ---- Multi-outcome pooling -------------------------------------------------------------------


def _hand_dersimonian_laird(means, variances):
    """Textbook DerSimonian-Laird with plain Python floats."""
    w = [1.0 / v for v in variances]
    total = sum(w)
    mu_fe = sum(wi * mi for wi, mi in zip(w, means)) / total
    q = sum(wi * (mi - mu_fe) ** 2 for wi, mi in zip(w, means))
    df = len(means) - 1
    tau2 = max(0.0, (q - df) / (total - sum(wi * wi for wi in w) / total))
    w_re = [1.0 / (v + tau2) for v in variances]
    mu = sum(wi * mi for wi, mi in zip(w_re, means)) / sum(w_re)
    return {"q": q, "df": df, "tau2": tau2, "w": w_re, "mu": mu, "var": 1.0 / sum(w_re)}


def test_var_r_reduces_to_the_independent_variance_at_r_identity() -> None:
    draws = {f"o{j}": _exact_sample(0.02 + 0.003 * j, 0.01 + 0.002 * j, seed=j) for j in range(4)}
    focal = _exact_sample(0.03, 0.02, seed=9)
    identity = _quiet(S.pool_multi_outcome_retention, draws, focal=focal, correlation=np.eye(4))
    zero = _quiet(S.pool_multi_outcome_retention, draws, focal=focal, correlation=0.0)
    assert identity.var_mu_correlated == pytest.approx(identity.var_mu_independent, rel=1e-12)
    assert identity.var_mu_common == pytest.approx(identity.var_mu_independent, rel=1e-12) and identity.r_common == 0.0
    assert zero.var_mu_common == pytest.approx(zero.var_mu_independent, rel=1e-12)
    np.testing.assert_allclose(identity.pooled_rho_draws_adjusted["focal"], identity.pooled_rho_draws["focal"])
    assert identity.adjusted_p10 == pytest.approx(identity.literal_p10, rel=1e-12)


@pytest.mark.parametrize("r", [0.25, 0.5, 1.0])
def test_equicorrelation_inflates_var_mu_by_one_plus_m_minus_one_r(r: float) -> None:
    base = _exact_sample(0.04, 0.015, seed=2)
    draws = {f"o{j}": base.copy() for j in range(5)}  # equal variances and weights, Q = 0
    common = _quiet(S.pool_multi_outcome_retention, draws, correlation=r)
    matrix = np.full((5, 5), r)
    np.fill_diagonal(matrix, 1.0)
    full = _quiet(S.pool_multi_outcome_retention, draws, correlation=matrix)
    assert common.tau2 == 0.0 and common.q_statistic == pytest.approx(0.0, abs=1e-20)
    np.testing.assert_allclose(list(common.weight_shares.values()), 0.2)
    assert common.var_mu_common == pytest.approx(common.var_mu_independent * (1.0 + 4.0 * r), rel=1e-12)
    assert full.var_mu_correlated == pytest.approx(common.var_mu_common, rel=1e-12)
    assert full.adjustment == "matrix" and common.adjustment == "common"


def test_dersimonian_laird_matches_a_hand_computation_on_three_outcomes() -> None:
    means, sds = [0.10, 0.05, 0.20], [0.02, 0.03, 0.04]
    draws = {name: _exact_sample(m, s, seed=k) for k, (name, m, s) in enumerate(zip("abc", means, sds))}
    pooled = _quiet(S.pool_multi_outcome_retention, draws)
    hand = _hand_dersimonian_laird(means, [s * s for s in sds])
    assert pooled.q_statistic == pytest.approx(hand["q"], rel=1e-10) and pooled.q_dof == 2
    assert pooled.q_statistic == pytest.approx(9.01639, abs=1e-5)
    assert pooled.tau2 == pytest.approx(hand["tau2"], rel=1e-10) and pooled.tau2 == pytest.approx(0.0029517, abs=1e-7)
    assert pooled.mu_hat == pytest.approx(hand["mu"], rel=1e-10)
    assert pooled.var_mu_independent == pytest.approx(hand["var"], rel=1e-10)
    np.testing.assert_allclose(list(pooled.weight_shares.values()), np.array(hand["w"]) / sum(hand["w"]), rtol=1e-10)
    assert pooled.q_p_value_upper == pytest.approx(stats.chi2.sf(hand["q"], 2), rel=1e-10)
    assert pooled.q_p_value_lower == pytest.approx(stats.chi2.cdf(hand["q"], 2), rel=1e-10)
    assert pooled.i2 == pytest.approx((hand["q"] - 2) / hand["q"], rel=1e-10)
    assert pooled.rho_at_mu == pytest.approx(np.exp(-hand["mu"]))
    assert pooled.dominant_outcome == "a" and pooled.rule_outcomes == ("a", "b", "c") and pooled.focal_label is None
    variants = pooled.rule_variants_df.set_index("label")
    loo = _hand_dersimonian_laird(means[1:], [s * s for s in sds[1:]])
    assert variants.loc["Leave one out: without a", "mu_hat"] == pytest.approx(loo["mu"], rel=1e-10)
    assert variants.loc["Leave one out: without a", "tau2"] == pytest.approx(loo["tau2"], rel=1e-10)
    equal = pooled.rule_variants_df.set_index("variant").loc["equal_weights"]
    assert equal["mu_hat"] == pytest.approx(np.mean(means), rel=1e-12)
    assert equal["var_mu"] == pytest.approx(sum(s * s + hand["tau2"] for s in sds) / 9.0, rel=1e-10)
    q_moment = pooled.rule_variants_df.set_index("variant").loc["q_moment"]
    assert q_moment["correlation"] == 0.0 and pooled.q_moment_r == 0.0  # Q > df: no moment correlation
    assert list(pooled.rule_variants_df["variant"]) == ["literal", "q_moment"] + ["leave_one_out"] * 3 + [
        "equal_weights"
    ]


def test_morris_shrinkage_with_zero_tau2_matches_the_closed_form() -> None:
    sds = {"x": 0.02, "y": 0.03, "z": 0.04}
    pool = {name: _exact_sample(0.05, sd, seed=k) for k, (name, sd) in enumerate(sds.items())}
    focal = _exact_sample(0.08, 0.05, seed=7)
    with pytest.warns(S.StabilityWarning) as caught:
        pooled = S.pool_multi_outcome_retention(pool, focal=focal, focal_label="primary", correlation=0.5)
    messages = [str(item.message) for item in caught]
    assert any(message.startswith("tau^2-hat = 0, so B = 1") for message in messages)
    assert any("lower chi2 tail" in message for message in messages)  # Q = 0: identical means
    assert pooled.tau2 == 0.0 and pooled.shrinkage_b["primary"] == 1.0
    var_ind = 1.0 / sum(1.0 / sd**2 for sd in sds.values())
    assert pooled.var_mu_independent == pytest.approx(var_ind, rel=1e-12)
    kappa = np.sqrt(var_ind / 0.05**2)
    assert pooled.kappa["primary"] == pytest.approx(kappa, rel=1e-12)
    expected = np.maximum(pooled.mu_hat + kappa * (focal - 0.08), 0.0)
    np.testing.assert_allclose(pooled.pooled_lambda_draws["primary"], expected, rtol=0, atol=1e-15)
    assert pooled.clip_share["primary"] == pytest.approx(np.mean(pooled.mu_hat + kappa * (focal - 0.08) < 0.0))
    kappa_adjusted = np.sqrt(pooled.var_mu_common / 0.05**2)
    adjusted = np.exp(-np.maximum(pooled.mu_hat + kappa_adjusted * (focal - 0.08), 0.0))
    np.testing.assert_allclose(pooled.pooled_rho_draws_adjusted["primary"], adjusted)
    assert pooled.literal_p10 == pytest.approx(np.quantile(np.exp(-expected), 0.10))
    assert pooled.rule_outcomes == ("primary",) and "B = 1" in pooled.warnings[0]
    summary = pooled.summary_df.set_index("outcome")
    assert summary.loc["primary", "role"] == "focal" and np.isnan(summary.loc["primary", "weight_share"])
    no_mu_term = _quiet(S.pool_multi_outcome_retention, pool, focal=focal, with_mu_uncertainty=False)
    np.testing.assert_allclose(no_mu_term.pooled_rho_draws["focal"], np.exp(-pooled.mu_hat))
    assert no_mu_term.kappa["focal"] == 0.0 and no_mu_term.literal_p10 == pytest.approx(np.exp(-pooled.mu_hat))
    unclipped = _quiet(S.pool_multi_outcome_retention, pool, focal=focal, clip_lambda_at_zero=False)
    np.testing.assert_allclose(unclipped.pooled_lambda_draws["focal"], pooled.mu_hat + kappa * (focal - 0.08))


def test_breakeven_correlation_and_correlated_variants() -> None:
    pool = {f"o{j}": _exact_sample(0.02 + 0.004 * j, 0.012, seed=j) for j in range(4)}
    focal = _exact_sample(0.03, 0.03, seed=11)
    literal = _quiet(S.pool_multi_outcome_retention, pool, focal=focal)
    threshold = literal.literal_p10 - 0.004
    pooled = _quiet(S.pool_multi_outcome_retention, pool, focal=focal, p10_threshold=threshold)
    assert pooled.breakeven_r is not None and 0.0 < pooled.breakeven_r < 1.0
    at_breakeven = _quiet(S.pool_multi_outcome_retention, pool, focal=focal, correlation=pooled.breakeven_r)
    assert _common_row(at_breakeven)["rho_P10"] == pytest.approx(threshold, abs=1e-9)
    above = _quiet(S.pool_multi_outcome_retention, pool, focal=focal, correlation=min(1.0, pooled.breakeven_r + 0.1),
                   p10_threshold=threshold)
    assert not _common_row(above)["passed"]
    assert _quiet(S.pool_multi_outcome_retention, pool, focal=focal, p10_threshold=0.999).breakeven_r is None
    assert _quiet(S.pool_multi_outcome_retention, pool, focal=focal, p10_threshold=0.5).breakeven_r is None
    variants = list(above.rule_variants_df["variant"])
    assert variants == ["literal", "correlation_common", "q_moment"] + ["leave_one_out"] * 4 + [
        "equal_weights", "equal_weights_correlated"
    ]
    variant_table = above.rule_variants_df
    np.testing.assert_allclose(variant_table["se_mu"].to_numpy() ** 2, variant_table["var_mu"].to_numpy())


def _common_row(pooled: S.PooledRetentionResult) -> pd.Series:
    return pooled.rule_variants_df.set_index("variant").loc["correlation_common"]


def _equicorrelation(r: float, size: int) -> np.ndarray:
    matrix = np.full((size, size), r)
    np.fill_diagonal(matrix, 1.0)
    return matrix


def test_correlation_matrix_variants_match_a_hand_var_r_with_unequal_variances() -> None:
    means, sds = [0.01, 0.05, 0.03], [0.010, 0.015, 0.020]
    draws = {name: _exact_sample(m, s, seed=k) for k, (name, m, s) in enumerate(zip("abc", means, sds))}
    focal = _exact_sample(0.025, 0.02, seed=5)
    matrix = np.array([[1.0, 0.5, 0.2], [0.5, 1.0, 0.4], [0.2, 0.4, 1.0]])
    pooled = _quiet(S.pool_multi_outcome_retention, draws, focal=focal, correlation=matrix)
    hand = _hand_dersimonian_laird(means, [s * s for s in sds])
    assert hand["tau2"] > 0.0 and pooled.tau2 == pytest.approx(hand["tau2"], rel=1e-10)
    variances = np.array(sds) ** 2

    def var_r(weights: np.ndarray, correlation: np.ndarray) -> float:
        covariance = correlation * np.sqrt(np.outer(variances, variances)) + hand["tau2"] * np.eye(3)
        return float(weights @ covariance @ weights / weights.sum() ** 2)

    weights = np.array(hand["w"])
    variants = pooled.rule_variants_df.set_index("variant")
    assert pooled.var_mu_correlated == pytest.approx(var_r(weights, matrix), rel=1e-10)
    assert variants.loc["correlation_matrix", "var_mu"] == pytest.approx(var_r(weights, matrix), rel=1e-10)
    r_bar = (0.5 + 0.2 + 0.4) / 3.0
    assert pooled.r_common == pytest.approx(r_bar, rel=1e-12)
    assert pooled.var_mu_common == pytest.approx(var_r(weights, _equicorrelation(r_bar, 3)), rel=1e-10)
    equal = variants.loc["equal_weights_correlated"]
    assert equal["var_mu"] == pytest.approx(var_r(np.ones(3), matrix), rel=1e-10)
    assert equal["mu_hat"] == pytest.approx(np.mean(means), rel=1e-12)
    # The adjusted draws and adjusted_p10 use the matrix, not its common r.
    shrink = 0.02**2 / (0.02**2 + hand["tau2"])
    kappa = np.sqrt(((1.0 - shrink) * 0.02**2 + shrink**2 * var_r(weights, matrix)) / 0.02**2)
    lam = np.maximum((1.0 - shrink) * 0.025 + shrink * hand["mu"] + kappa * (focal - 0.025), 0.0)
    assert pooled.adjustment == "matrix"
    np.testing.assert_allclose(pooled.pooled_rho_draws_adjusted["focal"], np.exp(-lam), rtol=1e-10)
    assert pooled.adjusted_p10 == pytest.approx(np.quantile(np.exp(-lam), 0.10), rel=1e-10)
    assert pooled.adjusted_p10 == pytest.approx(variants.loc["correlation_matrix", "rho_P10"], rel=1e-12)


def test_under_dispersed_outcomes_fail_the_q_clause_and_set_a_positive_q_moment_r() -> None:
    means, sds = [0.020, 0.022, 0.024], [0.010, 0.015, 0.020]
    draws = {name: _exact_sample(m, s, seed=k) for k, (name, m, s) in enumerate(zip("abc", means, sds))}
    with pytest.warns(S.StabilityWarning) as caught:
        pooled = S.pool_multi_outcome_retention(draws, focal=_exact_sample(0.02, 0.01, seed=8), correlation=0.3)
    messages = [str(item.message) for item in caught]
    assert any("lower chi2 tail" in message for message in messages)
    assert any(message.startswith("tau^2-hat = 0") for message in messages)
    hand = _hand_dersimonian_laird(means, [s * s for s in sds])
    assert hand["tau2"] == 0.0 and pooled.q_p_value_lower == pytest.approx(stats.chi2.cdf(hand["q"], 2), rel=1e-10)
    assert pooled.q_p_value_lower < 0.05
    r_q = 1.0 - hand["q"] / 2.0
    assert pooled.q_moment_r == pytest.approx(r_q, rel=1e-10) and r_q > 0.9
    variances = np.array(sds) ** 2
    weights = 1.0 / variances
    covariance = _equicorrelation(r_q, 3) * np.sqrt(np.outer(variances, variances))
    q_row = pooled.rule_variants_df.set_index("variant").loc["q_moment"]
    assert q_row["correlation"] == pytest.approx(r_q, rel=1e-10)
    assert q_row["var_mu"] == pytest.approx(float(weights @ covariance @ weights / weights.sum() ** 2), rel=1e-10)
    evaluation = S.evaluate_retention_rule(pooled, S.RetentionRuleSpec(max_cap_hit_share=None))
    clause = evaluation.summary_df.set_index("clause").loc["q_lower_tail"]
    assert clause["enabled"] and not clause["passed"] and clause["value"] == pooled.q_p_value_lower
    assert evaluation.fallback_recommended and any(reason.startswith("q_lower_tail") for reason in evaluation.reasons)
    lenient = S.evaluate_retention_rule(pooled, S.RetentionRuleSpec(max_cap_hit_share=None, max_q_lower_tail_p=None))
    assert lenient.summary_df.set_index("clause").loc["q_lower_tail", "passed"] is None


def test_pooling_reorders_a_slope_correlation_to_the_pooling_outcomes() -> None:
    windows, arms, _, _, periods, _ = _slope_panel()
    slope = S.unit_slope_correlation(windows, arms, periods=periods)
    order = ("y3", "y1", "y2")
    draws = {name: _exact_sample(0.02 + 0.004 * k, 0.01 + 0.003 * k, seed=k) for k, name in enumerate(order)}
    pooled = _quiet(S.pool_multi_outcome_retention, draws, correlation=slope)
    index = [slope.outcomes.index(name) for name in order]
    reordered = slope.correlation_matrix[np.ix_(index, index)]
    assert not np.allclose(reordered, slope.correlation_matrix)  # the order matters
    np.testing.assert_array_equal(pooled.correlation_matrix, reordered)
    direct = _quiet(S.pool_multi_outcome_retention, draws, correlation=reordered)
    assert pooled.var_mu_correlated == direct.var_mu_correlated
    assert pooled.r_common == pytest.approx(slope.r_hat, rel=1e-12)
    not_psd = dataclasses.replace(slope, correlation_matrix=_equicorrelation(-0.9, 3))
    with pytest.raises(ValueError, match="positive semi-definite"):
        _quiet(S.pool_multi_outcome_retention, draws, correlation=not_psd)


def test_degenerate_shares_feed_the_degenerate_clause() -> None:
    rho = np.tile([0.90, 0.93, 0.96, 0.99], 10)
    paths = {f"o{j}": _paths(rho * (1.0 - 0.02 * j)) for j in range(3)}
    paths["o0"][:10, WINDOW[0] - 1 :] = 0.0  # 25% of the draws have an identically zero window
    results = {name: S.estimate_post_pulse_retention(path, window=WINDOW, label=name) for name, path in paths.items()}
    pooled = _quiet(S.pool_multi_outcome_retention, results)
    assert pooled.degenerate_shares == {"o0": 0.25, "o1": 0.0, "o2": 0.0}
    assert pooled.summary_df.set_index("outcome").loc["o0", "degenerate_share"] == 0.25
    resampled = S.estimate_post_pulse_retention(paths["o0"], window=WINDOW, degenerate_policy="resample", label="o0")
    from_resampled = _quiet(S.pool_multi_outcome_retention, {**results, "o0": resampled})
    assert np.isfinite(resampled.lambda_draws).all() and from_resampled.degenerate_shares["o0"] == 0.25
    spec = S.RetentionRuleSpec(require_corr_adjusted=False, max_cap_hit_share=None)
    clause = S.evaluate_retention_rule(pooled, spec).summary_df.set_index("clause").loc["max_degenerate_share"]
    assert clause["enabled"] and not clause["passed"] and clause["value"] == 0.25 and "'o0'" in clause["detail"]
    relaxed = S.evaluate_retention_rule(pooled, dataclasses.replace(spec, max_degenerate_share=0.30))
    assert relaxed.summary_df.set_index("clause").loc["max_degenerate_share", "passed"]
    disabled = S.evaluate_retention_rule(pooled, dataclasses.replace(spec, max_degenerate_share=None))
    assert disabled.summary_df.set_index("clause").loc["max_degenerate_share", "passed"] is None
    with_nan = np.r_[np.full(4, np.nan), results["o1"].lambda_draws[4:]]
    arrays = _quiet(S.pool_multi_outcome_retention, {"a": with_nan, "b": results["o2"].lambda_draws})
    assert arrays.degenerate_shares == {"a": 0.1, "b": 0.0}


def test_pooling_with_nan_draws_retention_inputs_and_no_focal() -> None:
    results = {}
    for j in range(3):
        rho = np.r_[np.linspace(0.85, 1.0, 39), 0.97 + 0.01 * j]
        results[f"o{j}"] = S.estimate_post_pulse_retention(_paths(rho), window=WINDOW, label=f"o{j}")
    with_nan = np.r_[np.nan, results["o0"].lambda_draws[1:]]
    pooled = _quiet(S.pool_multi_outcome_retention, {"o0": with_nan, "o1": results["o1"], "o2": results["o2"]})
    assert pooled.lambda_means["o0"] == pytest.approx(np.mean(with_nan[1:]))
    assert np.isnan(pooled.pooled_lambda_draws["o0"][0]) and np.isfinite(pooled.pooled_lambda_draws["o0"][1:]).all()
    assert pooled.cap_hit_shares == {"o1": results["o1"].share_at_upper, "o2": results["o2"].share_at_upper}
    assert pooled.degenerate_shares == {"o0": 1 / 40, "o1": 0.0, "o2": 0.0}
    per_outcome = [np.quantile(np.exp(-d[np.isfinite(d)]), 0.10) for d in pooled.pooled_lambda_draws.values()]
    assert pooled.literal_p10 == pytest.approx(min(per_outcome))
    summary = pooled.summary_df.set_index("outcome")
    assert summary.loc["o0", "n_finite"] == 39 and np.isnan(summary.loc["o0", "cap_hit_share"])
    from_results = _quiet(S.pool_multi_outcome_retention, {k: results[k] for k in ("o0", "o1")}, focal=results["o2"])
    assert from_results.focal_label == "o2" and set(from_results.cap_hit_shares) == {"o0", "o1", "o2"}


def test_pooling_rejects_malformed_inputs() -> None:
    pool = {"a": _exact_sample(0.02, 0.01, seed=1), "b": _exact_sample(0.03, 0.01, seed=2)}
    slope_like = S.unit_slope_correlation(*_slope_panel()[:2], periods=np.arange(3, 10))
    cases = [
        dict(lambda_draws={"a": pool["a"]}),
        dict(lambda_draws={"a": pool["a"], "b": np.full(400, 0.02)}),
        dict(lambda_draws={"a": pool["a"], "b": np.r_[np.inf, pool["b"][1:]]}),
        dict(lambda_draws={"a": pool["a"], "b": np.r_[0.1, np.full(399, np.nan)]}),
        dict(focal=pool["a"], focal_label="a"),
        dict(focal_label="primary"),
        dict(correlation=1.5),
        dict(correlation=-1.5),
        dict(correlation=np.eye(3)),
        dict(correlation=np.array([[1.0, 0.2], [0.3, 1.0]])),
        dict(correlation=np.array([[0.9, 0.2], [0.2, 1.0]])),
        dict(correlation=slope_like),
        dict(correlation="high"),
        dict(p10_threshold=0.0),
        dict(p10_threshold=1.2),
        dict(quantiles=(0.5, 0.1, 0.9)),
    ]
    for case in cases:
        kwargs = {"lambda_draws": pool, **case}
        with pytest.raises(ValueError):
            _quiet(S.pool_multi_outcome_retention, kwargs.pop("lambda_draws"), **kwargs)
    three = {name: _exact_sample(0.02, 0.01, seed=k) for k, name in enumerate("xyz")}
    with pytest.raises(ValueError, match="positive semi-definite"):
        _quiet(S.pool_multi_outcome_retention, three, correlation=_equicorrelation(-0.9, 3))
    with pytest.raises(ValueError, match="positive semi-definite"):
        _quiet(S.pool_multi_outcome_retention, three, correlation=-0.6)  # below -1 / (M - 1) = -0.5
    unequal = {name: _exact_sample(0.02, sd, seed=k) for k, (name, sd) in enumerate(zip("xyz", (0.01, 0.02, 0.03)))}
    assert _quiet(S.pool_multi_outcome_retention, unequal, correlation=-0.5).r_common == -0.5


# ---- Retention rule --------------------------------------------------------------------------


def _agreeing_lambda_draws() -> dict[str, np.ndarray]:
    """Five outcomes with equal variances, ``Q`` near its df, small lambda (P10 well above 0.90)."""
    offsets = (-0.007, -0.003, 0.0, 0.003, 0.007)
    return {f"o{j}": _exact_sample(0.02 + offset, 0.005, seed=j) for j, offset in enumerate(offsets)}


def test_rule_spec_validation_and_clause_switches() -> None:
    for bad in (
        dict(p10_threshold=0.0), dict(p10_threshold=1.5), dict(require_corr_adjusted="yes"),
        dict(require_leave_one_out=1), dict(max_single_outcome_weight=0.0), dict(max_cap_hit_share=1.2),
        dict(max_degenerate_share=1.5), dict(max_degenerate_share=-0.1), dict(max_q_lower_tail_p=1.0),
    ):
        with pytest.raises(ValueError):
            S.RetentionRuleSpec(**bad)
    with pytest.raises(TypeError):
        S.RetentionRuleSpec(0.9)  # keyword-only: a pre-registered spec cannot be built positionally
    pooled = _quiet(S.pool_multi_outcome_retention, _agreeing_lambda_draws(), focal=_exact_sample(0.02, 0.005, seed=9))
    default = S.evaluate_retention_rule(pooled)
    clauses = default.summary_df.set_index("clause")
    assert default.literal_passed and not default.robust_passed and default.fallback_recommended
    assert default.headline_regimes == S.BOUND_REGIMES and default.warnings == ()
    assert not clauses.loc["correlation_adjusted_p10", "passed"]  # no correlation supplied
    assert not clauses.loc["max_cap_hit_share", "passed"]  # arrays carry no cap-hit shares
    assert clauses.loc["leave_one_out_p10", "passed"] and clauses.loc["max_single_outcome_weight", "passed"]
    assert clauses.loc["max_degenerate_share", "passed"] and clauses.loc["max_degenerate_share", "value"] == 0.0
    assert clauses.loc["q_lower_tail", "passed"]
    assert len(default.reasons) == 2 and default.reasons[0].startswith("correlation_adjusted_p10")
    relaxed = S.RetentionRuleSpec(require_corr_adjusted=False, max_cap_hit_share=None)
    passed = S.evaluate_retention_rule(pooled, relaxed)
    assert passed.robust_passed and not passed.fallback_recommended and passed.reasons == ()
    assert passed.headline_regimes == S.PERSISTENCE_REGIMES
    disabled = passed.summary_df.set_index("clause")
    assert not disabled.loc["max_cap_hit_share", "enabled"] and disabled.loc["max_cap_hit_share", "passed"] is None
    correlated = _quiet(
        S.pool_multi_outcome_retention, _agreeing_lambda_draws(), focal=_exact_sample(0.02, 0.005, seed=9),
        correlation=0.5,
    )
    assert S.evaluate_retention_rule(correlated, S.RetentionRuleSpec(max_cap_hit_share=None)).robust_passed
    with pytest.warns(S.StabilityWarning, match="differs from the pooling's p10_threshold"):
        strict = S.evaluate_retention_rule(correlated, S.RetentionRuleSpec(p10_threshold=0.999, max_cap_hit_share=None))
    assert not strict.literal_passed and strict.reasons[0].startswith("literal_p10")
    assert len(strict.warnings) == 1 and "breakeven_r" in strict.warnings[0]
    with pytest.raises(ValueError):
        S.evaluate_retention_rule("pooled")
    with pytest.raises(ValueError):
        S.evaluate_retention_rule(pooled, spec={"p10_threshold": 0.9})


def _cap_censored_panel(seed: int = 11, n: int = 400) -> dict[str, np.ndarray]:
    """One cap-censored outcome (85% of draws flat at rho = 1) and five well-spread outcomes."""
    rng = np.random.default_rng(seed)
    grid = np.linspace(0.5, 1.0, 1001)

    def snap(rho: np.ndarray) -> np.ndarray:
        return grid[np.clip(np.round((rho - 0.5) / GRID_STEP).astype(int), 0, 1000)]

    def lognormal_rho(mean: float, sd: float) -> np.ndarray:
        sigma2 = np.log(1.0 + (sd / mean) ** 2)
        return snap(np.exp(-rng.lognormal(np.log(mean) - sigma2 / 2.0, np.sqrt(sigma2), size=n)))

    capped = np.ones(n)
    capped[: n - int(round(0.85 * n))] = snap(rng.uniform(0.96, 0.999, size=n - int(round(0.85 * n))))
    rhos = {"capped": capped}
    for j in range(4):
        rhos[f"outcome_{j}"] = lognormal_rho(0.025, 0.025)
    rhos["primary"] = lognormal_rho(0.025, 0.03)
    return {name: _paths(rho) for name, rho in rhos.items()}


def test_cap_censored_dominant_outcome_triggers_pileup_weight_flags_and_fallback() -> None:
    with pytest.warns(S.StabilityWarning) as caught:
        report = S.diagnose_trajectory_stability(
            _cap_censored_panel(), window=WINDOW, focal="primary", correlation=0.3
        )
    capped = report.retention["capped"]
    assert capped.share_at_upper == pytest.approx(0.85) and capped.boundary_pileup_flag
    assert not any(report.retention[name].boundary_pileup_flag for name in report.outcomes if name != "capped")
    assert report.pooled.dominant_outcome == "capped" and 0.55 < report.pooled.dominant_weight < 0.70
    clauses = report.rule.summary_df.set_index("clause")
    assert not clauses.loc["max_single_outcome_weight", "passed"]
    assert not clauses.loc["max_cap_hit_share", "passed"]
    assert clauses.loc["max_cap_hit_share", "value"] == pytest.approx(0.85)
    assert report.rule.literal_passed and report.rule.fallback_recommended and not report.rule.robust_passed
    assert report.headline_regimes == ("U", "L") and report.verdict is None
    reasons = " ".join(report.rule.reasons)
    assert "max_single_outcome_weight" in reasons and "max_cap_hit_share" in reasons
    messages = [str(item.message) for item in caught]
    assert any("boundary pile-up" in message and "'capped'" in message for message in messages)
    assert any(message.startswith("dominant outcome: 'capped'") for message in messages)
    assert any(message.startswith("retention rule fallback") for message in messages)
    assert set(report.warnings) == set(messages)


# ---- Report ----------------------------------------------------------------------------------


def test_diagnose_trajectory_stability_report_contract() -> None:
    panel = _cap_censored_panel(seed=12, n=120)
    wearout = S.wearout_ancova(**_wearout_kwargs(_wearout_panel(n=60, effect=0.2, seed=6)))
    audit = A.TransportAuditResult(
        labels=("short", "long"), overlap_horizon=5, n_units={"short": 3, "long": 4}, profile_draws={},
        diff_draws=np.zeros((2, 5)), mean_diff_draws=np.zeros(2), quantiles={}, prob_positive=0.5, flag=False,
        threshold=None, unobserved_cell_share={"all": 0.3}, unobserved_effect_share={"all": 0.4},
        cumulative_horizon=10,
    )
    supports = {name: np.r_[np.ones(9), 0.5] for name in panel}
    with pytest.warns(S.StabilityWarning):
        report = S.diagnose_trajectory_stability(
            panel, window=WINDOW, focal="primary", pooling_outcomes=["capped", "outcome_0", "outcome_1"],
            support=supports, wearout={"A-C": wearout}, transport_audit=audit,
            rule=S.RetentionRuleSpec(max_cap_hit_share=None),
        )
    assert report.verdict is None and report.outcomes == ("capped", "outcome_0", "outcome_1", "primary")
    assert report.focal == "primary" and report.window == WINDOW
    assert report.unobserved_effect_share == {"all": 0.4} and report.unobserved_cell_share == {"all": 0.3}
    assert report.wearout == {"A-C": wearout} and report.slope_correlation is None
    assert report.retention["capped"].support_source == "user"
    assert list(report.summary_df["outcome"]) == list(report.outcomes)
    assert {"share_at_upper", "boundary_pileup_flag", "pooled_rho_P10", "role"} <= set(report.summary_df.columns)
    assert report.data_health_certificate_df.columns[0] == "component"
    assert set(report.data_health_certificate_df["component"]) == {
        "retention[capped]", "retention[outcome_0]", "retention[outcome_1]", "retention[primary]", "pooling",
        "wearout[A-C]",
    }
    assert any(message.startswith("wear-out test 'A-C' is underpowered") for message in report.warnings)
    for bad in (
        dict(focal="unknown"), dict(pooling_outcomes=["capped"]), dict(pooling_outcomes=["capped", "primary"]),
        dict(pooling_outcomes=["capped", "nope"]), dict(wearout="A-C"), dict(transport_audit="audit"),
        dict(rule={"p10_threshold": 0.9}), dict(support={"nope": np.ones(10)}),
    ):
        kwargs = {"window": WINDOW, "focal": "primary", **bad}
        with pytest.raises(ValueError):
            _quiet(S.diagnose_trajectory_stability, panel, **kwargs)
    with pytest.raises(ValueError, match="support keys"):
        _quiet(S.diagnose_trajectory_stability, panel, window=WINDOW, support={"primary": np.ones(10), "Primary": 1})
    with pytest.raises(ValueError):
        S.diagnose_trajectory_stability({}, window=WINDOW)


def test_report_with_a_slope_correlation_two_pooling_outcomes_and_every_wearout_message() -> None:
    windows, arms, _, _, periods, _ = _slope_panel()
    slope = S.unit_slope_correlation({name: windows[name] for name in ("y1", "y2")}, arms, periods=periods)
    rng = np.random.default_rng(21)
    traj = {
        "y1": _paths(rng.uniform(0.90, 0.99, size=200)),
        "y2": _paths(rng.uniform(0.85, 0.99, size=200)),
        "primary": _paths(rng.uniform(0.88, 0.99, size=200)),
    }
    base = S.wearout_ancova(**_wearout_kwargs(_wearout_panel(n=60, effect=0.2, seed=6)))
    assert not base.powered_for_half_wearout and base.early_estimate > 0.0
    wearouts = {
        "underpowered": base,
        "powered": dataclasses.replace(base, powered_for_half_wearout=True),
        "negative": dataclasses.replace(base, early_estimate=-0.1),
        "exact": dataclasses.replace(base, change_mdd=0.0),
    }
    with pytest.warns(S.StabilityWarning) as caught:
        report = S.diagnose_trajectory_stability(
            traj, window=WINDOW, focal="primary", correlation=slope, wearout=wearouts
        )
    assert report.slope_correlation is slope
    assert "slope_correlation" in set(report.data_health_certificate_df["component"])
    np.testing.assert_array_equal(report.pooled.correlation_matrix, slope.correlation_matrix)
    variants = report.pooled.rule_variants_df
    loo = variants.loc[variants["variant"] == "leave_one_out"].set_index("left_out")
    for left_out, kept in (("y1", "y2"), ("y2", "y1")):  # one outcome left: its own mean and variance
        assert loo.loc[left_out, "tau2"] == 0.0
        assert loo.loc[left_out, "mu_hat"] == pytest.approx(report.pooled.lambda_means[kept], rel=1e-12)
        assert loo.loc[left_out, "var_mu"] == pytest.approx(report.pooled.lambda_variances[kept], rel=1e-12)
    weight = report.rule.summary_df.set_index("clause").loc["max_single_outcome_weight"]
    assert not weight["passed"] and "unreachable" in weight["detail"]
    messages = [str(item.message) for item in caught]
    assert any(message.startswith("the max_single_outcome_weight clause cannot pass") for message in messages)
    assert any(message.startswith("wear-out test 'underpowered' is underpowered") for message in messages)
    assert any(message.startswith("wear-out test 'negative' has a non-positive early-level") for message in messages)
    assert any(message.startswith("wear-out test 'exact' has a zero standard error") for message in messages)
    assert not any("'powered'" in message for message in messages)
    assert set(report.warnings) == set(messages)
    with pytest.raises(ValueError, match="withholds a verdict"):
        dataclasses.replace(report, verdict="persistent")


# ---- Pseudo-horizon separation ---------------------------------------------------------------


def _weibull_tail(observed: int, horizon: int, scale: float = 40.0, shape: float = 2.0) -> np.ndarray:
    """``(H,)`` survival weights: 1 on ``t <= L`` and ``S(t) / S(L)`` beyond, for a Weibull ``S``."""
    t = np.arange(1, horizon + 1, dtype=float)
    weights = np.exp(-((t / scale) ** shape)) / np.exp(-((observed / scale) ** shape))
    weights[:observed] = 1.0
    return weights


def test_pseudo_horizon_separation_matches_closed_forms() -> None:
    rng = np.random.default_rng(4)
    plateau = rng.uniform(0.1, 0.3, size=(200, OBSERVED))
    pulse = rng.uniform(0.0, 0.2, size=(200, OBSERVED))
    arms = {"BC": plateau, "AB": pulse, "AC": plateau + pulse}
    survival = _weibull_tail(OBSERVED, HORIZON)
    result = S.pseudo_horizon_separation(
        arms, pseudo_horizon=ANCHOR, observed_horizon=OBSERVED, horizon=HORIZON, pulse_decay_tau=PULSE_TAU,
        survival=survival, targets={"model-free": (0.9, 0.5), "posterior": 0.4},
    )
    steps = np.arange(1, OBSERVED - ANCHOR + 1, dtype=float)
    at_anchor = ANCHOR - 1
    pulse_path = pulse[:, [at_anchor]] * np.exp(-steps / PULSE_TAU)
    sums = {
        "U": (plateau[:, [at_anchor]] + pulse_path).sum(axis=1),
        "L": (plateau[:, [at_anchor]] * 0.85**steps + pulse_path).sum(axis=1),
    }
    np.testing.assert_allclose(result.gap_draws, plateau[:, at_anchor] * np.sum(1.0 - 0.85**steps))
    assert result.delta == pytest.approx(np.mean(sums["U"]) - np.mean(sums["L"]))
    window_sum = (plateau + pulse)[:, ANCHOR:OBSERVED].sum(axis=1)
    assert result.observed_window_sd == pytest.approx(np.std(window_sum, ddof=1))
    power = result.power_df.set_index("target")
    assert list(power.index) == ["posterior SD of the observed cumulative (same draws)", "model-free", "posterior"]
    z_alpha, z_power = stats.norm.ppf(0.95), stats.norm.ppf(0.80)
    for target, sd in zip(power.index, (np.std(window_sum, ddof=1), 0.5, 0.4)):
        combined = np.sqrt(np.var(sums["L"], ddof=1) + sd**2)
        assert power.loc[target, "combined_sd"] == pytest.approx(combined)
        assert power.loc[target, "power"] == pytest.approx(stats.norm.cdf(result.delta / combined - z_alpha))
        assert power.loc[target, "mdd"] == pytest.approx((z_alpha + z_power) * combined)
    z_u = (0.9 - np.mean(sums["U"])) / np.sqrt(np.var(sums["U"], ddof=1) + 0.25)
    assert power.loc["model-free", "z_target_minus_U"] == pytest.approx(z_u)
    assert np.isnan(power.loc["posterior", "z_target_minus_L"])
    assert result.min_detectable_gap == power["mdd"].iloc[0]
    assert result.share_gap_above_min == pytest.approx(np.mean(result.gap_draws > result.min_detectable_gap))
    tail = np.arange(1, HORIZON - OBSERVED + 1, dtype=float)
    expected_lever = np.sum(1.0 - 0.85**steps) / np.sum((1.0 - 0.85**tail) * survival[OBSERVED:])
    assert result.lever_ratio == pytest.approx(expected_lever)
    assert result.tail_lever_multiple == pytest.approx(1.0 / expected_lever)
    assert set(result.lever_df["segment"]) == {"pseudo", "tail"} and len(result.lever_df) == steps.size + tail.size
    whole = S.pseudo_horizon_separation(
        arms["AC"][:, :ANCHOR], pseudo_horizon=ANCHOR, observed_horizon=OBSERVED, horizon=HORIZON,
        targets={"external": 0.3}, min_detectable_gap=0.05,
    )
    np.testing.assert_allclose(whole.gap_draws, arms["AC"][:, at_anchor] * np.sum(1.0 - 0.85**steps))
    assert whole.observed_window_sd is None and whole.min_detectable_gap == 0.05
    assert whole.lever_ratio == pytest.approx(np.sum(1.0 - 0.85**steps) / np.sum(1.0 - 0.85**tail))
    nested = M.Nested3ArmReconciliationResult(
        trt_catt_bc=plateau, trt_catt_ab=pulse, trt_catt_ac=plateau + pulse, post_att_bc=plateau[:, 0],
        post_att_ab=pulse[:, 0], post_att_ac=plateau[:, 0],
    )
    from_nested = S.pseudo_horizon_separation(
        nested, pseudo_horizon=ANCHOR, observed_horizon=OBSERVED, horizon=HORIZON, pulse_decay_tau=PULSE_TAU,
        survival=survival,
    )
    np.testing.assert_allclose(from_nested.gap_draws, result.gap_draws)


def test_pseudo_horizon_separation_rejects_malformed_inputs() -> None:
    draws = np.random.default_rng(0).uniform(size=(50, OBSERVED))
    n_tail = HORIZON - OBSERVED
    cases = [
        dict(pseudo_horizon=OBSERVED),
        dict(observed_horizon=HORIZON),
        dict(pseudo_horizon=OBSERVED + 1, observed_horizon=OBSERVED + 2),  # anchor beyond the draws
        dict(rho_bounds=(1.0, 0.85)),
        dict(rho_bounds=(0.85, 1.2)),
        dict(pulse_decay_tau=PULSE_TAU),  # the plateau / pulse split needs arm trajectories
        dict(pulse_decay_tau=-1.0, traj_draws={"BC": draws, "AB": draws}),
        dict(pulse_decay_tau=PULSE_TAU, traj_draws={"BC": draws, "AB": draws, "AC": draws}),  # AC != BC + AB
        dict(survival=np.ones(HORIZON - 1)),
        dict(survival=np.r_[np.ones(OBSERVED - 1), 0.9, np.ones(n_tail)]),
        dict(survival=np.r_[np.ones(OBSERVED - 1), 1.0 - 1e-13, np.ones(n_tail)]),  # strict 1.0 on t <= L
        dict(survival=np.r_[np.ones(OBSERVED), np.zeros(n_tail)]),  # zero tail lever
        dict(alpha=0.0),
        dict(targets={"x": -1.0}),
        dict(targets={"x": (1.0, 2.0, 3.0)}),
        dict(traj_draws=draws[:, : OBSERVED - 2]),  # the draws stop before L and no target is given
        dict(traj_draws={"AB": draws}),
    ]
    for case in cases:
        kwargs = {
            "traj_draws": draws, "pseudo_horizon": ANCHOR, "observed_horizon": OBSERVED, "horizon": HORIZON, **case,
        }
        with pytest.raises(ValueError):
            S.pseudo_horizon_separation(kwargs.pop("traj_draws"), **kwargs)


# ---- Determinism -----------------------------------------------------------------------------------


def test_every_function_is_deterministic() -> None:
    panel = _cap_censored_panel(seed=13, n=80)

    def run():
        report = _quiet(S.diagnose_trajectory_stability, panel, window=WINDOW, focal="primary", correlation=0.2,
                        degenerate_policy="resample", seed=5)
        wear = S.wearout_ancova(**_wearout_kwargs(_wearout_panel(seed=8)))
        windows, arms, _, _, periods, _ = _slope_panel(seed=14)
        slope = S.unit_slope_correlation(windows, arms, periods=periods)
        separation = S.pseudo_horizon_separation(panel["capped"], pseudo_horizon=6, observed_horizon=10, horizon=40)
        return report, wear, slope, separation

    first, second = run(), run()
    pd.testing.assert_frame_equal(first[0].summary_df, second[0].summary_df)
    pd.testing.assert_frame_equal(first[0].pooled.rule_variants_df, second[0].pooled.rule_variants_df)
    for name in first[0].retention:
        np.testing.assert_array_equal(first[0].retention[name].rho_draws, second[0].retention[name].rho_draws)
    pd.testing.assert_frame_equal(first[1].summary_df, second[1].summary_df)
    np.testing.assert_array_equal(first[2].correlation_matrix, second[2].correlation_matrix)
    pd.testing.assert_frame_equal(first[3].power_df, second[3].power_df)
    assert first[0].warnings == second[0].warnings


# ---- Integration smoke (golden fixture, tiny configuration) ---------------------------------------


CFG = M.RankKLongBetConfig(
    uncertainty_propagation=True, n_forecast_draws=2, max_forecast_draws=2, num_chains=2, draws_per_component=20,
    num_burnin=20, num_trees_mu=20, num_trees_tau_mode1=12, num_trees_tau_higher=8, random_seed=2026,
    use_spatial_gp=False,
)


def test_integration_smoke_on_the_golden_fixture() -> None:
    panel = make_panel("balanced")
    y = panel["y"].copy()
    short = np.r_[0:8, 30:38]  # 8 treated and 8 control units observe 6 of the 8 post periods
    y[short, 4 + 6 :] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = M.fit_rank_k_longbet(
            y, panel["z"], panel["x"], cfg=CFG, baseline_forecast=panel["f"],
            baseline_forecast_draws=panel["f_draws"], cohorts=panel["cohorts"],
        )
    result = S.estimate_post_pulse_retention(fit, window=(3, 8), label="golden")
    assert result.support_source == "post_window_length" and result.trajectory_source.startswith("event_time")
    assert result.weights.tolist() == [60.0, 60.0, 60.0, 60.0, 44.0, 44.0]
    traj = A.event_time_trajectory(fit, weighting="unit")
    direct = S.estimate_post_pulse_retention(traj, window=(3, 8), support=result.weights)
    np.testing.assert_array_equal(direct.rho_draws, result.rho_draws)
    assert np.isfinite(result.rho_draws).all() and np.all((result.rho_draws >= 0.5) & (result.rho_draws <= 1.0))
    assert result.rho_draws.shape == (traj.shape[0],)
    cert = result.data_health_certificate_df
    assert (cert.loc[cert["kind"] == "distribution", "n_inf"] == 0).all()
