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

"""Tests of the trajectory layer upgrade and the required-persistence decision layer of ``_rank_k_plots``.

Synthetic draws (no MCMC) check the survival and discount helpers against their closed forms, the
plateau / pulse split (additivity, horizon-``L`` identity, observed window kept), the ``rho_sd``
deprecation path, the required persistence ``X`` (a hand-computed case, an independent per-draw
recomputation, the identity ``Pr(X <= k) = Pr(payback_U <= L + k)`` and the regime ordering), the
``alpha_share = 1`` reduction to :func:`roi_break_even_table`, the break-even retention against a
brute-force characterization, the placebo gate, the labelled extra regimes (a fixed retention of 1
reproduces Bound U bit for bit, a per-draw retention is paired by index, headline / diagnostic roles,
placebo retentions or the skip note), the negative-draw guard and input validation; one smoke test
runs the golden fixture end to end and the figure is smoke-tested under ``Agg``.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import optimize

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

BAND_LABEL = "assumed Beta(rho, sd) scenario band — not identified by the data"
CAVEAT = "pre-registered scenario, not a statistical bound; faster decay is not excluded"
CFG = M.RankKLongBetConfig(
    uncertainty_propagation=True, n_forecast_draws=2, max_forecast_draws=2, num_chains=2, draws_per_component=20,
    num_burnin=20, num_trees_mu=20, num_trees_tau_mode1=12, num_trees_tau_higher=8, random_seed=2026,
    use_spatial_gp=False,
)


def _nested_arms(seed: int, n_draws: int, n_obs: int, plateau: float = 0.2, pulse: float = 0.1) -> dict:
    """Non-negative BC / AB draws and AC = BC + AB."""
    rng = np.random.default_rng(seed)
    bc = np.abs(rng.normal(plateau, 0.1, size=(n_draws, n_obs)))
    ab = np.abs(rng.normal(pulse, 0.05, size=(n_draws, n_obs)))
    return {"BC": bc, "AB": ab, "AC": bc + ab}


def _value_weights(result: P.RequiredPersistenceResult, alpha: float) -> np.ndarray:
    """``m_t(alpha) delta_t`` recomputed from the result's inputs."""
    return (alpha * result.alpha_share + (1.0 - result.alpha_share)) * result.discount


def _reference_persistence(values: np.ndarray, cost: float, n_obs: int) -> np.ndarray:
    """Required persistence from its definition, one draw at a time (first crossing, linear interpolation)."""
    out = np.full(values.shape[0], np.inf)
    for draw, row in enumerate(values):
        cumulative = 0.0
        for month, value in enumerate(row, start=1):
            if cumulative + value > cost:
                out[draw] = max(0.0, month - 1 - n_obs + (cost - cumulative) / value)
                break
            cumulative += value
    return out


def _value_at_horizon(
    arms: dict, rho: np.ndarray, weights: np.ndarray, survival: np.ndarray, pulse_decay_tau: float
) -> np.ndarray:
    """``V_rho(H)`` per draw written out directly: observed value + pulse value + plateau value at ``rho``."""
    n_obs = arms["AC"].shape[1]
    steps = np.arange(1, weights.size - n_obs + 1)
    tail_weights = weights[n_obs:] * survival[n_obs:]
    observed = (arms["AC"] * weights[:n_obs]).sum(axis=1)
    pulse = arms["AB"][:, -1] * np.sum(np.exp(-steps / pulse_decay_tau) * tail_weights)
    plateau = arms["BC"][:, -1] * (np.power(rho[:, None], steps[None, :]) @ tail_weights)
    return observed + pulse + plateau


def _year_options() -> dict:
    """Options of the ``year_result`` fixture: monthly periods, five-year horizon, split, survival and NPV."""
    return dict(
        cost_grid=(1.0, 3.0, 6.0), alpha_grid=(1.0, 0.5), alpha_share=np.array([0.8] * 12 + [0.7] * 48),
        pulse_decay_tau=3.0, survival=P.weibull_survival_weights(60, 12), discount_rate=0.10,
    )


def _regime_rows(frame: pd.DataFrame, regime: str) -> pd.DataFrame:
    """The rows of one regime without the ``regime`` and ``role`` columns, re-indexed from 0."""
    return frame.loc[frame["regime"] == regime].drop(columns=["regime", "role"]).reset_index(drop=True)


@pytest.fixture(scope="module")
def short_arms():
    return _nested_arms(seed=7, n_draws=60, n_obs=6)


@pytest.fixture(scope="module")
def year_arms():
    return _nested_arms(seed=11, n_draws=300, n_obs=12)


@pytest.fixture(scope="module")
def year_result(year_arms):
    return P.compute_required_persistence(year_arms, **_year_options())


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ---- Survival and discounting --------------------------------------------------------------------


def test_survival_and_discount_weights_match_closed_forms() -> None:
    weights = P.weibull_survival_weights(60, 12)
    scale = 60.0 / math.log(2.0) ** (1.0 / 2.5)
    survival = np.exp(-((np.arange(1, 61) / scale) ** 2.5))
    assert scale == pytest.approx(69.47, abs=5e-3)
    assert survival[11] == pytest.approx(0.9877, abs=5e-5) and survival[59] == pytest.approx(0.5, abs=1e-12)
    np.testing.assert_array_equal(weights[:12], 1.0)
    np.testing.assert_allclose(weights[12:], survival[12:] / survival[11], rtol=1e-13)
    assert np.all(np.diff(weights) <= 0.0) and weights.min() > 0.0
    # shape 1 is the exponential life curve: S(t) / S(L) = 2^(-(t - L) / median).
    exponential = P.weibull_survival_weights(24, 5, median_periods=30.0, shape=1.0)
    np.testing.assert_allclose(exponential[5:], 2.0 ** (-np.arange(1, 20) / 30.0), rtol=1e-12)
    discount = P.discount_weights(60)
    assert discount[0] == 1.0 and discount[12] == pytest.approx(1.0 / (1.0 + 0.10), rel=1e-14)
    np.testing.assert_allclose(discount, (1.0 + 0.10) ** (-np.arange(60) / 12.0), rtol=1e-13)
    quarterly = P.discount_weights(8, annual_rate=0.05, periods_per_year=4)
    np.testing.assert_allclose(quarterly, 1.05 ** (-np.arange(8) / 4.0), rtol=1e-13)
    np.testing.assert_array_equal(P.discount_weights(5, annual_rate=0.0), np.ones(5))
    for bad in (
        dict(horizon=10, n_observed=12), dict(horizon=60, n_observed=0), dict(horizon=12.5, n_observed=12),
        dict(horizon=60, n_observed=12, median_periods=0.0), dict(horizon=60, n_observed=12, shape=-1.0),
        dict(horizon=60, n_observed=12, median_periods=1e-3),
    ):
        with pytest.raises(ValueError):
            P.weibull_survival_weights(**bad)
    for bad in (
        dict(horizon=0), dict(horizon=12, annual_rate=-1.0), dict(horizon=12, annual_rate=float("nan")),
        dict(horizon=12, periods_per_year=0),
    ):
        with pytest.raises(ValueError):
            P.discount_weights(**bad)


# ---- Trajectory layer ----------------------------------------------------------------------------


def test_split_tail_is_additive_and_keeps_the_observed_window(short_arms) -> None:
    survival = P.weibull_survival_weights(30, 6, median_periods=24.0)
    d = P.extrapolate_rank_k_decay(
        short_arms, horizon=30, rho=0.9, rho_sd=0.05, rho_grid=(0.8, 0.95), num_mc_per_draw=4, seed=3,
        pulse_decay_tau=3.0, survival=survival,
    )
    assert d.arm_roles == {"BC": "plateau", "AB": "pulse", "AC": "total"}
    assert d.pulse_decay_tau == 3.0 and np.array_equal(d.survival, survival)
    steps = np.arange(1, 25)
    sigma = survival[6:]
    for arm, observed in short_arms.items():
        np.testing.assert_array_equal(d.monthly_draws[arm][:, :6], observed)
        np.testing.assert_allclose(d.cumulative_draws[arm][:, 5], observed.sum(axis=1), atol=1e-10)
    np.testing.assert_allclose(d.monthly_draws["BC"][:, 6:], short_arms["BC"][:, -1:] * 0.9**steps * sigma, rtol=1e-13)
    pulse_tail = short_arms["AB"][:, -1:] * np.exp(-steps / 3.0) * sigma
    np.testing.assert_allclose(d.monthly_draws["AB"][:, 6:], pulse_tail, rtol=1e-13)
    split_sum = d.monthly_draws["BC"][:, 6:] + d.monthly_draws["AB"][:, 6:]
    np.testing.assert_array_equal(d.monthly_draws["AC"][:, 6:], split_sum)
    P.assert_nested_additivity(d.monthly_draws["BC"], d.monthly_draws["AB"], d.monthly_draws["AC"], atol=1e-12)
    cumulative = d.cumulative_draws
    P.assert_nested_additivity(cumulative["BC"], cumulative["AB"], cumulative["AC"], atol=1e-10)
    for per_arm in d.rho_grid_cumulative.values():  # the grid moves the plateau only
        P.assert_nested_additivity(per_arm["BC"], per_arm["AB"], per_arm["AC"], atol=1e-10)
        np.testing.assert_array_equal(per_arm["AB"], d.cumulative_draws["AB"])
    band = d.random_rho_cumulative_draws  # so does the Beta band
    P.assert_nested_additivity(band["BC"], band["AB"], band["AC"], atol=1e-10)
    np.testing.assert_array_equal(band["AB"], np.repeat(d.cumulative_draws["AB"], 4, axis=0))
    assert "exp(-h / 3)" in d.scenario_note and "not estimated from the data" in d.scenario_note
    # A reconciliation result is accepted and maps BC / AB / AC onto the plateau / pulse / total roles.
    rec = M.Nested3ArmReconciliationResult(
        trt_catt_bc=short_arms["BC"], trt_catt_ab=short_arms["AB"], trt_catt_ac=short_arms["AC"],
        post_att_bc=short_arms["BC"].mean(axis=1), post_att_ab=short_arms["AB"].mean(axis=1),
        post_att_ac=short_arms["AC"].mean(axis=1),
    )
    from_rec = P.extrapolate_rank_k_decay(
        rec, horizon=30, rho=0.9, rho_sd=None, pulse_decay_tau=3.0, survival=survival
    )
    for arm in short_arms:
        np.testing.assert_array_equal(from_rec.monthly_draws[arm], d.monthly_draws[arm])
    # Without the split every arm decays geometrically and the survival weights multiply its tail.
    unsplit = P.extrapolate_rank_k_decay(short_arms, horizon=30, rho=0.9, rho_sd=None, survival=survival)
    assert set(unsplit.arm_roles.values()) == {"geometric"} and unsplit.pulse_decay_tau is None
    for arm, observed in short_arms.items():
        geometric_tail = observed[:, -1:] * 0.9**steps * sigma
        np.testing.assert_allclose(unsplit.monthly_draws[arm][:, 6:], geometric_tail, rtol=1e-13)
    # survival=None keeps the 0.7.0 extrapolation bit for bit.
    legacy = P.extrapolate_rank_k_decay(short_arms, horizon=30, rho=0.9, rho_sd=None)
    assert legacy.survival is None
    for arm, observed in short_arms.items():
        np.testing.assert_array_equal(legacy.monthly_draws[arm][:, 6:], observed[:, -1:] * np.power(0.9, steps))


def test_implicit_rho_sd_is_none_without_warning_and_band_is_labelled_as_assumed(short_arms) -> None:
    effect = short_arms["AC"]
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        implicit = P.extrapolate_rank_k_decay(effect, horizon=12, num_mc_per_draw=2)
        explicit = P.extrapolate_rank_k_decay(effect, horizon=12, rho_sd=0.07, num_mc_per_draw=2)
        no_band = P.extrapolate_rank_k_decay(effect, horizon=12, rho_sd=None)
    assert implicit.rho_sd is None
    assert implicit.random_rho_cumulative_draws is None
    np.testing.assert_array_equal(
        implicit.cumulative_draws["effect"], no_band.cumulative_draws["effect"]
    )
    assert set(explicit.cumulative_summary_df["rho_band_label"]) == {BAND_LABEL}
    assert set(no_band.cumulative_summary_df["rho_band_label"]) == {""}
    assert set(implicit.cumulative_summary_df["rho_band_label"]) == {""}
    assert BAND_LABEL in explicit.scenario_note and "Beta" not in no_band.scenario_note
    assert "Beta" not in implicit.scenario_note
    axes = P.plot_rank_k_decay_extrapolation(explicit)
    legend = [text.get_text() for text in axes[1].get_legend().get_texts()]
    assert any(BAND_LABEL in text for text in legend)
    assert not any("rho ~ Beta band" in text for text in legend)
    # Under the plateau / pulse split, _decay_tail_label distinguishes the arms and the pulse arm skips the rho band
    # and rho-grid lines.
    split_decay = P.extrapolate_rank_k_decay(
        short_arms, horizon=18, rho=0.9, rho_sd=0.05, rho_grid=(0.8, 0.95), num_mc_per_draw=2, pulse_decay_tau=3.0
    )
    assert P._decay_tail_label("pulse", split_decay) == "exp(-h / 3) tail"
    assert P._decay_tail_label("total", split_decay) == "plateau rho = 0.9 + pulse tail"
    assert P._decay_tail_label("plateau", split_decay) == "rho = 0.9 tail"
    split_axes = P.plot_rank_k_decay_extrapolation(split_decay)
    split_legend = [text.get_text() for text in split_axes[1].get_legend().get_texts()]
    assert any("BC rho grid" in text for text in split_legend)
    assert any("AC rho grid" in text for text in split_legend)
    assert not any("AB rho grid" in text for text in split_legend)
    assert sum(BAND_LABEL in text for text in split_legend) == 1


# ---- Required persistence ------------------------------------------------------------------------


def test_required_persistence_matches_a_hand_computed_example() -> None:
    # Two draws, L = 2, H = 6, whole level persisting, no discounting or survival, alpha = 1.
    draws = np.array([[1.0, 1.0], [2.0, 2.0]])
    res = P.compute_required_persistence(
        draws, horizon=6, cost_grid=(3.5,), alpha_grid=(1.0,), rho_bounds=(0.5, 1.0), checkpoints=(1, 4)
    )
    # Draw 1: V_U = 1, 2, 3, 4, ... exceeds 3.5 in period 4 -> X = (4 - 1 - 2) + (3.5 - 3) / 1 = 1.5.
    # Draw 2: V = 2, 4 already exceeds 3.5 in period 2 -> X = 0.
    np.testing.assert_allclose(res.required_persistence_draws[(3.5, 1.0)], [1.5, 0.0])
    np.testing.assert_array_equal(res.payback_month_draws[("U", 3.5, 1.0)], [4, 2])
    np.testing.assert_array_equal(res.payback_month_draws[("L", 3.5, 1.0)], [7, 2])  # 2 + 0.5 + ... < 3.5
    np.testing.assert_array_equal(res.payback_month_draws[("rho0", 3.5, 1.0)], [7, 2])
    # Break-even: rho + rho^2 + rho^3 + rho^4 > 1.5 for draw 1, rounded up to the 5e-5 grid; 0 for draw 2.
    root = optimize.brentq(lambda r: r + r**2 + r**3 + r**4 - 1.5, 0.0, 1.0)
    rho_be = res.breakeven_rho_draws[(3.5, 1.0)]
    assert root <= rho_be[0] <= root + 5e-5 and rho_be[1] == 0.0
    row = res.required_persistence_df.iloc[0]
    assert row["pr_X_never"] == 0.0 and row["pr_paid_back_in_window"] == 0.5
    assert row["pr_X_le_1"] == 0.5 and row["pr_X_le_4"] == 1.0
    assert res.arm_roles == {"effect": "total"} and res.pulse_decay_tau is None
    assert "whole level" in row["persisting_level"] and "in full" in row["persisting_level"]
    # An exact tie (the value equals the cost at the end of period L) keeps the strict payback rule.
    tie = P.compute_required_persistence(draws, horizon=6, cost_grid=(2.0,), alpha_grid=(1.0,))
    x_tie = tie.required_persistence_draws[(2.0, 1.0)]
    assert 0.0 < x_tie[0] < 1e-300 and x_tie[1] == 0.0
    assert tie.payback_month_draws[("U", 2.0, 1.0)][0] == 3 and tie.cdf_identity_checked
    # So does a near tie whose fraction is below half an ulp: V(12) = 3.0999999999999996 < 3.1 pays back in
    # period 13, so X must stay above 0 (it used to round onto 0 and trip the CDF identity).
    row_near = [0.7] * 3 + [0.0] * 8 + [1.0]
    near = P.compute_required_persistence(
        np.array([row_near, row_near]), horizon=24, cost_grid=(3.1,), alpha_grid=(1.0,)
    )
    x_near = near.required_persistence_draws[(3.1, 1.0)]
    assert np.all(x_near > 0.0) and np.all(x_near < 1e-15) and near.cdf_identity_checked
    np.testing.assert_array_equal(near.payback_month_draws[("U", 3.1, 1.0)], [13, 13])
    assert near.required_persistence_df.loc[0, "pr_paid_back_in_window"] == 0.0
    # With rho_U < 1, X is persistence at that retention and the labels say so.
    partial = P.compute_required_persistence(
        draws, horizon=6, cost_grid=(3.5,), alpha_grid=(1.0,), rho_bounds=(0.5, 0.95)
    )
    assert "at rho = 0.95" in partial.required_persistence_df.loc[0, "persisting_level"]
    assert "persistence at rho = 0.95 (Bound U)" in partial.scenario_note
    assert "full persistence" not in partial.scenario_note


def test_required_persistence_cdf_identity_and_regime_ordering(year_result, year_arms) -> None:
    res = year_result
    assert res.cdf_identity_checked and res.warnings == ()
    assert res.n_observed == 12 and res.horizon == 60 and res.report_months == (12, 24, 36, 60)
    for (cost, alpha), persistence in res.required_persistence_draws.items():
        upper = res.payback_month_draws[("U", cost, alpha)]
        for k in range(49):
            assert np.mean(persistence <= k) == np.mean(upper <= 12 + k), (cost, alpha, k)
        table = res.payback_cdf_df
        cdf = table.loc[(table["cost"] == cost) & (table["alpha"] == alpha)]
        by_regime = {regime: cdf.loc[cdf["regime"] == regime, "pr_paid_back"].to_numpy() for regime in P._REGIMES}
        assert np.all(by_regime["U"] >= by_regime["L"]) and np.all(by_regime["L"] >= by_regime["rho0"])
        # Independent recomputation of X from its definition.
        values = res.monthly_draws["U"] * _value_weights(res, alpha)[None, :]
        expected = _reference_persistence(values, cost, 12)
        np.testing.assert_array_equal(np.isinf(expected), np.isinf(persistence))
        finite = np.isfinite(expected)
        np.testing.assert_allclose(persistence[finite], expected[finite], rtol=1e-12, atol=1e-12)
        row = res.required_persistence_df.loc[
            (res.required_persistence_df["cost"] == cost) & (res.required_persistence_df["alpha"] == alpha)
        ].iloc[0]
        for k in (12, 24, 48):
            assert row[f"pr_X_le_{k}"] == np.mean(persistence <= k)
        assert row["X_P50"] == np.quantile(persistence, 0.5, method="inverted_cdf")
    # The cumulative effect at L is the observed sum in every regime; the split holds draw by draw.
    for regime, monthly in res.monthly_draws.items():
        np.testing.assert_allclose(monthly[:, :12].sum(axis=1), year_arms["AC"].sum(axis=1), atol=1e-12)
        assert regime in res.regime_labels
    assert CAVEAT in res.regime_labels["U"] and CAVEAT in res.regime_labels["L"]
    certificate = res.data_health_certificate_df
    assert set(zip(certificate["cost"], certificate["alpha"])) == set(res.required_persistence_draws)
    assert (certificate["n_nan"] == 0).all() and (certificate["n_inf"] == 0).all()
    never = certificate.loc[certificate["quantity"] == "required_persistence_X", "never_code"]
    assert (never == 49.0).all()


def test_discount_array_matches_discount_rate(year_arms, year_result) -> None:
    options = _year_options()
    del options["discount_rate"]
    explicit = P.compute_required_persistence(year_arms, discount=P.discount_weights(60, annual_rate=0.10), **options)
    np.testing.assert_array_equal(explicit.discount, year_result.discount)
    for key, persistence in year_result.required_persistence_draws.items():
        np.testing.assert_array_equal(explicit.required_persistence_draws[key], persistence)
        np.testing.assert_array_equal(explicit.breakeven_rho_draws[key], year_result.breakeven_rho_draws[key])
    pd.testing.assert_frame_equal(explicit.decision_df, year_result.decision_df)
    assert "discounted" in explicit.scenario_note and "undiscounted" not in explicit.scenario_note
    # Discounting lowers every cumulative value after the first period, so Bound U pays back later.
    undiscounted = P.compute_required_persistence(year_arms, **options)
    for key, months in undiscounted.payback_month_draws.items():
        assert np.all(year_result.payback_month_draws[key] >= months)


def test_alpha_share_one_reproduces_roi_break_even_table(short_arms) -> None:
    effect = short_arms["AC"]
    costs, alphas = (1.0, 3.0), (0.25, 0.5, 1.0)
    res = P.compute_required_persistence(
        effect, horizon=24, cost_grid=costs, alpha_grid=alphas, rho_bounds=(0.85, 1.0)
    )
    decision = res.decision_df
    for regime, rho in (("U", 1.0), ("L", 0.85), ("rho0", 0.0)):
        cumulative = {"effect": P._extend_with_decay(effect, rho, 24).sum(axis=1)}
        table = P.roi_break_even_table(cumulative, cost_grid=costs, alpha_grid=alphas, horizon_periods=24)
        for row in table.itertuples(index=False):
            cell = decision.loc[
                (decision["regime"] == regime) & (decision["horizon"] == 24) & (decision["cost"] == row.cost)
                & (decision["alpha"] == row.alpha)
            ]
            assert cell["pr_roi_positive"].item() == pytest.approx(row.pr_roi_positive, abs=1e-12)
    observed = P.roi_break_even_table(
        {"effect": effect.sum(axis=1)}, cost_grid=costs, alpha_grid=alphas, horizon_periods=6
    )
    at_l = decision.loc[decision["horizon"] == 6]
    assert set(at_l["regime"]) == {"observed"}
    np.testing.assert_allclose(at_l["pr_roi_positive"].to_numpy(), observed["pr_roi_positive"].to_numpy(), atol=1e-12)


def test_breakeven_rho_matches_brute_force(short_arms) -> None:
    survival = P.weibull_survival_weights(24, 6, median_periods=30.0)
    res = P.compute_required_persistence(
        short_arms, horizon=24, cost_grid=(1.5, 3.0, 4.5), alpha_grid=(1.0, 0.5), alpha_share=0.8,
        pulse_decay_tau=2.0, survival=survival, discount_rate=0.10,
    )
    n_interior = 0
    for (cost, alpha), rho_be in res.breakeven_rho_draws.items():
        weights = (alpha * 0.8 + 0.2) * P.discount_weights(24, annual_rate=0.10)
        n_draws = len(rho_be)
        at_zero = _value_at_horizon(short_arms, np.zeros(n_draws), weights, survival, 2.0)
        at_one = _value_at_horizon(short_arms, np.ones(n_draws), weights, survival, 2.0)
        np.testing.assert_array_equal(rho_be == 0.0, at_zero > cost)
        np.testing.assert_array_equal(np.isinf(rho_be), at_one <= cost)
        interior = (rho_be > 0.0) & np.isfinite(rho_be)
        n_interior += int(interior.sum())
        # The smallest grid rho that pays back: it does, and one grid step below does not.
        on_grid = np.where(interior, rho_be, 0.0)
        assert np.all(_value_at_horizon(short_arms, on_grid, weights, survival, 2.0)[interior] > cost)
        below = np.where(interior, rho_be - 5e-5, 0.0)
        assert np.all(_value_at_horizon(short_arms, below, weights, survival, 2.0)[interior] <= cost)
    assert n_interior > 20


# ---- Placebo gate and guards ---------------------------------------------------------------------


def test_placebo_gate_flags_a_null_shifted_to_pr_055(year_arms) -> None:
    rng = np.random.default_rng(5)
    n_draws = year_arms["AC"].shape[0]
    paying = rng.permutation(n_draws) < int(0.55 * n_draws)  # 55% of the draws carry an effect, 45% none
    bc = np.where(paying[:, None], 0.2, 0.0) * np.ones((n_draws, 12))
    ab = np.where(paying[:, None], 0.05, 0.0) * np.ones((n_draws, 12))
    shifted = {"BC": bc, "AB": ab, "AC": bc + ab}
    null_bc = rng.normal(0.0, 0.005, size=(n_draws, 12))
    null_ab = rng.normal(0.0, 0.005, size=(n_draws, 12))
    null = {"BC": null_bc, "AB": null_ab, "AC": null_bc + null_ab}
    survival = P.weibull_survival_weights(60, 12)
    options = dict(cost_grid=(1.0, 3.0), alpha_grid=(1.0, 0.5), pulse_decay_tau=3.0, survival=survival)
    flagged = P.compute_required_persistence(year_arms, placebo_draws=shifted, **options)
    decision = flagged.decision_df
    cell = decision.loc[
        (decision["regime"] == "U") & (decision["horizon"] == 60) & (decision["cost"] == 1.0)
        & (decision["alpha"] == 1.0)
    ].iloc[0]
    assert cell["pr_roi_positive_placebo"] == pytest.approx(0.55) and not cell["decision_immune"]
    assert cell["immunity_margin"] == pytest.approx(-0.05) and cell["pr_roi_positive"] > 0.9
    summary = flagged.placebo_summary_df.set_index("scope")
    assert summary.loc["all", "n_cells"] == len(decision)
    assert summary.loc["all", "n_immune"] == int(decision["decision_immune"].sum())
    assert summary.loc["observed+headline", "n_cells"] == 4 * 7 and summary.loc["U", "n_failing"] >= 3
    assert summary.loc["reference", "n_cells"] == 4 * 3 and "diagnostic" not in summary.index
    assert list(flagged.placebo_summary_df["kind"]) == ["regime"] * 4 + ["role"] * 2 + ["all"]
    clean = P.compute_required_persistence(year_arms, placebo_draws=null, **options)
    assert clean.decision_df["decision_immune"].all()
    assert (clean.placebo_summary_df["n_failing"] == 0).all()
    assert P.compute_required_persistence(year_arms, **options).placebo_summary_df is None


# ---- Extra regimes and roles ---------------------------------------------------------------------


def test_fixed_extra_regimes_reproduce_bound_u_and_the_floor_bitwise(year_arms, year_result) -> None:
    res = P.compute_required_persistence(year_arms, extra_regimes={"full": 1.0, "none": 0.0}, **_year_options())
    for twin, builtin in (("full", "U"), ("none", "rho0")):
        np.testing.assert_array_equal(res.monthly_draws[twin], res.monthly_draws[builtin])
        for alpha in res.alpha_grid:
            np.testing.assert_array_equal(
                res.cumulative_value_draws[(twin, alpha)], res.cumulative_value_draws[(builtin, alpha)]
            )
            for cost in res.cost_grid:
                np.testing.assert_array_equal(
                    res.payback_month_draws[(twin, cost, alpha)], res.payback_month_draws[(builtin, cost, alpha)]
                )
        for frame in (
            res.payback_df, res.payback_cdf_df, res.monthly_summary_df, res.cumulative_value_df,
            res.milestone_table_df, res.decision_df,
        ):
            pd.testing.assert_frame_equal(_regime_rows(frame, twin), _regime_rows(frame, builtin))
    # The built-in regimes, X and the break-even retention do not move.
    for key, months in year_result.payback_month_draws.items():
        np.testing.assert_array_equal(res.payback_month_draws[key], months)
    for key, persistence in year_result.required_persistence_draws.items():
        np.testing.assert_array_equal(res.required_persistence_draws[key], persistence)
        np.testing.assert_array_equal(res.breakeven_rho_draws[key], year_result.breakeven_rho_draws[key])
    builtin_rows = res.payback_df["regime"].isin(P._REGIMES)
    pd.testing.assert_frame_equal(res.payback_df.loc[builtin_rows].reset_index(drop=True), year_result.payback_df)
    pd.testing.assert_frame_equal(res.required_persistence_df, year_result.required_persistence_df)
    breakeven = res.breakeven_rho_df
    pd.testing.assert_frame_equal(breakeven[list(year_result.breakeven_rho_df.columns)], year_result.breakeven_rho_df)
    np.testing.assert_allclose(breakeven["pr_rho_ge_breakeven_full"], 1.0 - breakeven["pr_rho_be_never"], atol=1e-12)
    np.testing.assert_array_equal(breakeven["pr_rho_ge_breakeven_none"], breakeven["pr_rho_be_zero"])
    assert isinstance(res.regime_rhos["full"], float) and res.regime_rhos["full"] == 1.0
    assert res.regime_roles["full"] == res.regime_roles["none"] == "diagnostic"
    assert res.regime_labels["full"] == "full: plateau rho = 1 (diagnostic regime)"
    assert "each extra regime (full, none; fixed or one per draw)" in res.scenario_note
    certificate = res.data_health_certificate_df
    regime_free = certificate["quantity"].isin(["required_persistence_X", "breakeven_rho"])
    assert certificate.loc[regime_free, ["regime", "role"]].isna().all().all()
    regime_rows = certificate.loc[~regime_free]
    assert (regime_rows["role"] == regime_rows["regime"].map(res.regime_roles)).all()
    # Without the plateau / pulse split the whole level persists, and the twin is still bit-identical.
    whole = P.compute_required_persistence(
        year_arms["AC"], cost_grid=(1.0, 3.0), alpha_grid=(1.0, 0.5), extra_regimes={"full": 1.0}
    )
    np.testing.assert_array_equal(whole.monthly_draws["full"], whole.monthly_draws["U"])
    for (regime, cost, alpha), months in whole.payback_month_draws.items():
        if regime == "full":
            np.testing.assert_array_equal(months, whole.payback_month_draws[("U", cost, alpha)])


def test_per_draw_extra_regime_is_paired_by_index(year_arms) -> None:
    rng = np.random.default_rng(13)
    n_draws = year_arms["AC"].shape[0]
    rho_draws = rng.uniform(0.8, 1.0, size=n_draws)
    survival = P.weibull_survival_weights(60, 12)
    options = dict(
        cost_grid=(1.0, 3.0), alpha_grid=(1.0, 0.5), pulse_decay_tau=3.0, survival=survival, discount_rate=0.10
    )
    res = P.compute_required_persistence(year_arms, extra_regimes={"P": rho_draws}, **options)
    steps = np.arange(1, 49)
    plateau = year_arms["BC"][:, -1:] * rho_draws[:, None] ** steps * survival[12:]
    pulse = year_arms["AB"][:, -1:] * np.exp(-steps / 3.0) * survival[12:]
    np.testing.assert_array_equal(res.monthly_draws["P"][:, :12], year_arms["AC"])
    np.testing.assert_allclose(res.monthly_draws["P"][:, 12:], plateau + pulse, rtol=1e-13)
    # Draw i under the per-draw retention is draw i under the fixed retention rho_draws[i].
    for draw in (0, 17, n_draws - 1):
        fixed = P.compute_required_persistence(year_arms, extra_regimes={"P": float(rho_draws[draw])}, **options)
        np.testing.assert_array_equal(res.monthly_draws["P"][draw], fixed.monthly_draws["P"][draw])
        for key, months in res.payback_month_draws.items():
            assert months[draw] == fixed.payback_month_draws[key][draw], key
    reversed_rho = P.compute_required_persistence(year_arms, extra_regimes={"P": rho_draws[::-1]}, **options)
    assert not np.array_equal(reversed_rho.monthly_draws["P"], res.monthly_draws["P"])
    assert res.regime_rhos["P"] is not rho_draws
    np.testing.assert_array_equal(res.regime_rhos["P"], rho_draws)
    assert (res.payback_df.loc[res.payback_df["regime"] == "P", "rho"] == np.median(rho_draws)).all()
    assert res.regime_labels["P"].startswith("P: plateau rho per draw (median ")
    for (cost, alpha), rho_be in res.breakeven_rho_draws.items():
        table = res.breakeven_rho_df
        row = table.loc[(table["cost"] == cost) & (table["alpha"] == alpha)].iloc[0]
        assert row["pr_rho_ge_breakeven_P"] == np.mean(rho_draws >= rho_be)
    # A retention between the floor and Bound U, draw by draw, pays back between them.
    assert res.cdf_identity_checked
    for (regime, cost, alpha), months in res.payback_month_draws.items():
        if regime == "P":
            assert np.all(months >= res.payback_month_draws[("U", cost, alpha)])
            assert np.all(months <= res.payback_month_draws[("rho0", cost, alpha)])


def test_role_columns_and_headline_regimes(year_arms) -> None:
    options = dict(cost_grid=(1.0, 3.0), alpha_grid=(1.0,), pulse_decay_tau=3.0)
    res = P.compute_required_persistence(
        year_arms, extra_regimes={"P": 0.95, "Q": 0.9}, headline_regimes=("U", "P"), **options
    )
    assert res.regime_roles == {
        "observed": "headline", "U": "headline", "L": "diagnostic", "rho0": "reference", "P": "headline",
        "Q": "diagnostic",
    }
    for frame in (
        res.payback_df, res.payback_cdf_df, res.monthly_summary_df, res.cumulative_value_df, res.milestone_table_df,
        res.decision_df,
    ):
        columns = list(frame.columns)
        assert columns[columns.index("regime") + 1] == "role"
        assert (frame["role"] == frame["regime"].map(res.regime_roles)).all()
        assert {"U", "L", "rho0", "P", "Q"} <= set(frame["regime"])
    assert set(res.decision_df.loc[res.decision_df["horizon"] == 12, "role"]) == {"headline"}
    assert res.regime_labels["P"].endswith("(headline regime)")
    assert res.regime_labels["Q"].endswith("(diagnostic regime)")
    default = P.compute_required_persistence(year_arms, extra_regimes={"P": 0.95}, **options)
    assert default.regime_roles["P"] == "diagnostic" and default.regime_roles["U"] == "headline"
    single = P.compute_required_persistence(year_arms, extra_regimes={"P": 0.95}, headline_regimes="P", **options)
    assert single.regime_roles["U"] == single.regime_roles["L"] == "diagnostic"
    assert single.regime_roles["P"] == "headline"
    bounds_only = P.compute_required_persistence(year_arms, headline_regimes=(), **options)
    assert bounds_only.regime_roles == {
        "observed": "headline", "U": "diagnostic", "L": "diagnostic", "rho0": "reference",
    }
    for bad in ("rho0", "observed", "R", ("U", "U"), None, 3):
        with pytest.raises(ValueError, match="headline_regimes"):
            P.compute_required_persistence(year_arms, extra_regimes={"P": 0.95}, headline_regimes=bad, **options)


def test_extra_regime_labels_and_retentions_are_validated(short_arms) -> None:
    effect = short_arms["AC"]
    grids = dict(horizon=24, cost_grid=(1.0,), alpha_grid=(1.0,))
    for label in ("U", "L", "rho0", "observed", "observed+headline", "diagnostic", "reference", "all"):
        with pytest.raises(ValueError, match="reserved"):
            P.compute_required_persistence(effect, extra_regimes={label: 0.9}, **grids)
    n_draws = effect.shape[0]
    bad_rhos = (
        1.2, -0.1, float("nan"), float("inf"), None, "high", np.r_[np.full(n_draws - 1, 0.9), 1.01],
        np.r_[np.full(n_draws - 1, 0.9), np.nan], np.full(n_draws + 1, 0.9), np.full((n_draws, 1), 0.9),
    )
    for rho in bad_rhos:
        with pytest.raises(ValueError, match="extra_regimes"):
            P.compute_required_persistence(effect, extra_regimes={"P": rho}, **grids)
    for bad in ({"": 0.9}, {1: 0.9}, [("P", 0.9)]):
        with pytest.raises(ValueError, match="extra_regimes"):
            P.compute_required_persistence(effect, extra_regimes=bad, **grids)
    edges = P.compute_required_persistence(effect, extra_regimes={"zero": 0, "one": np.ones(n_draws)}, **grids)
    assert edges.regime_rhos["zero"] == 0.0 and edges.regime_rhos["one"].shape == (n_draws,)
    # An empty mapping is the same as no extra regimes.
    empty = P.compute_required_persistence(effect, extra_regimes={}, **grids)
    plain = P.compute_required_persistence(effect, **grids)
    assert empty.regime_roles == plain.regime_roles and empty.regime_rhos == plain.regime_rhos
    for name in ("payback_df", "decision_df", "breakeven_rho_df", "data_health_certificate_df"):
        pd.testing.assert_frame_equal(getattr(empty, name), getattr(plain, name))


def test_placebo_extra_regimes_gate_the_extra_regimes_or_skip_them_with_a_note(year_arms) -> None:
    rng = np.random.default_rng(21)
    n_draws = year_arms["AC"].shape[0]
    null_bc = rng.normal(0.0, 0.005, size=(n_draws, 12))
    null_ab = rng.normal(0.0, 0.005, size=(n_draws, 12))
    null = {"BC": null_bc, "AB": null_ab, "AC": null_bc + null_ab}
    rho_draws = rng.uniform(0.85, 1.0, size=n_draws)
    options = dict(cost_grid=(1.0, 3.0), alpha_grid=(1.0, 0.5), pulse_decay_tau=3.0, extra_regimes={"P": rho_draws})
    with pytest.warns(UserWarning, match="placebo_extra_regimes was not given"):
        skipped = P.compute_required_persistence(year_arms, placebo_draws=null, **options)
    assert len(skipped.warnings) == 1 and "skip the placebo gate" in skipped.warnings[0]
    assert set(skipped.placebo_regime_rhos) == {"U", "L", "rho0"}
    decision = skipped.decision_df
    p_rows = decision.loc[decision["regime"] == "P"]
    assert len(p_rows) == 3 * 2 * 2 and decision["decision_immune"].dtype == "boolean"
    assert p_rows["pr_roi_positive_placebo"].isna().all() and p_rows["immunity_margin"].isna().all()
    assert p_rows["decision_immune"].isna().all() and decision["decision_immune"].notna().sum() == len(decision) - 12
    summary = skipped.placebo_summary_df.set_index("scope")
    assert summary.loc["P", "n_cells"] == summary.loc["diagnostic", "n_cells"] == 0
    assert summary.loc["P", "n_cells_ungated"] == summary.loc["diagnostic", "n_cells_ungated"] == 12
    assert np.isnan(summary.loc["P", "max_pr_roi_positive_placebo"]) and pd.isna(summary.loc["P", "max_at"])
    assert summary.loc["all", "n_cells"] == len(decision) - len(p_rows) and summary.loc["all", "n_cells_ungated"] == 12
    assert summary.loc["observed+headline", "n_cells"] == 4 * 7
    assert summary.loc["observed+headline", "n_cells_ungated"] == 0
    assert "placebo_value_at_H_P" not in set(skipped.data_health_certificate_df["quantity"])
    # An ungated headline extra regime shows up in the headline count as ungated cells, not as immune ones.
    with pytest.warns(UserWarning, match="placebo_extra_regimes was not given"):
        ungated_headline = P.compute_required_persistence(
            year_arms, placebo_draws=null, headline_regimes=("U", "L", "P"), **options
        )
    summary = ungated_headline.placebo_summary_df.set_index("scope")
    assert summary.loc["observed+headline", "n_cells"] == 4 * 7
    assert summary.loc["observed+headline", "n_cells_ungated"] == 12
    # With placebo retentions the extra regime is gated and counted with the diagnostic regimes.
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        gated = P.compute_required_persistence(
            year_arms, placebo_draws=null, placebo_extra_regimes={"P": rho_draws}, **options
        )
    assert gated.warnings == ()
    np.testing.assert_array_equal(gated.placebo_regime_rhos["P"], rho_draws)
    assert gated.placebo_regime_rhos["P"] is not rho_draws and gated.placebo_regime_rhos["U"] == 1.0
    table = gated.placebo_summary_df
    assert list(table["kind"]) == ["regime"] * 5 + ["role"] * 3 + ["all"]
    assert table["role"].iloc[:-1].tolist() == [
        "headline", "headline", "headline", "reference", "diagnostic", "headline", "diagnostic", "reference",
    ]
    assert pd.isna(table["role"].iloc[-1])
    assert (table["n_cells_ungated"] == 0).all()
    summary = table.set_index("scope")
    p_rows = gated.decision_df.loc[gated.decision_df["regime"] == "P"]
    assert summary.loc["diagnostic", "n_cells"] == summary.loc["P", "n_cells"] == len(p_rows) == 12
    assert summary.loc["diagnostic", "n_immune"] == int(p_rows["decision_immune"].sum()) == 12
    assert summary.loc["observed+headline", "n_cells"] == 4 * 7
    assert summary.loc["all", "n_cells"] == len(gated.decision_df)
    certificate = gated.data_health_certificate_df
    placebo_p = certificate.loc[certificate["quantity"] == "placebo_value_at_H_P"]
    assert len(placebo_p) == 4 and set(placebo_p["regime"]) == {"P"} and set(placebo_p["role"]) == {"diagnostic"}
    # The placebo retention is paired with the placebo draws: the gate is the placebo run as traj_draws.
    with pytest.warns(UserWarning, match="negative"):
        as_draws = P.compute_required_persistence(null, **options)
    expected = as_draws.decision_df.loc[as_draws.decision_df["regime"] == "P", "pr_roi_positive"].to_numpy()
    np.testing.assert_array_equal(p_rows["pr_roi_positive_placebo"].to_numpy(), expected)
    # So is a fixed placebo retention.
    fixed = dict(options, extra_regimes={"P": 0.95})
    scalar = P.compute_required_persistence(year_arms, placebo_draws=null, placebo_extra_regimes={"P": 0.95}, **fixed)
    assert scalar.placebo_regime_rhos["P"] == 0.95
    with pytest.warns(UserWarning, match="negative"):
        scalar_as_draws = P.compute_required_persistence(null, **fixed)
    np.testing.assert_array_equal(
        scalar.decision_df.loc[scalar.decision_df["regime"] == "P", "pr_roi_positive_placebo"].to_numpy(),
        scalar_as_draws.decision_df.loc[scalar_as_draws.decision_df["regime"] == "P", "pr_roi_positive"].to_numpy(),
    )
    headline = P.compute_required_persistence(
        year_arms, placebo_draws=null, placebo_extra_regimes={"P": rho_draws}, headline_regimes=("U", "L", "P"),
        **options,
    )
    summary = headline.placebo_summary_df.set_index("scope")
    assert summary.loc["observed+headline", "n_cells"] == 4 * 10 and "diagnostic" not in summary.index
    assert P.compute_required_persistence(year_arms, **options).placebo_regime_rhos is None
    for bad in (
        dict(placebo_extra_regimes={"P": rho_draws}),
        dict(placebo_draws=null, placebo_extra_regimes={"Q": rho_draws}),
        dict(placebo_draws=null, placebo_extra_regimes={}),
        dict(placebo_draws=null, placebo_extra_regimes={"P": rho_draws[:-1]}),
        dict(placebo_draws=null, placebo_extra_regimes={"P": 1.5}),
    ):
        with pytest.raises(ValueError, match="placebo_extra_regimes"):
            P.compute_required_persistence(year_arms, **options, **bad)


def test_negative_monthly_effects_disable_the_cdf_identity_with_a_warning(year_arms) -> None:
    arms = {key: value.copy() for key, value in year_arms.items()}
    arms["AB"][:, 4] -= 0.5
    arms["AC"] = arms["BC"] + arms["AB"]
    options = dict(cost_grid=(1.0, 3.0), alpha_grid=(1.0, 0.5), pulse_decay_tau=3.0)
    with pytest.warns(UserWarning, match="not checked"):
        res = P.compute_required_persistence(arms, **options)
    assert not res.cdf_identity_checked and len(res.warnings) == 1 and "negative" in res.warnings[0]
    assert not res.required_persistence_df.drop(columns=["persisting_level"]).isna().any().any()
    assert (res.data_health_certificate_df["n_nan"] == 0).all()
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        assert P.compute_required_persistence(year_arms, **options).cdf_identity_checked


def test_checkpoints_beyond_extrapolated_window_are_nan_with_warning() -> None:
    arms = _nested_arms(seed=19, n_draws=200, n_obs=15)  # L = 15: at H = 60, H - L = 45 < 48
    with pytest.warns(UserWarning, match="checkpoints \\[48\\] exceed the extrapolated window"):
        res_60 = P.compute_required_persistence(
            arms, horizon=60, cost_grid=(6.0, 12.0), alpha_grid=(1.0,), pulse_decay_tau=3.0
        )
    assert len(res_60.warnings) == 1 and "horizon - L = 45" in res_60.warnings[0]
    assert res_60.required_persistence_df["pr_X_le_48"].isna().all()
    assert res_60.required_persistence_df[["pr_X_le_12", "pr_X_le_24"]].notna().all().all()
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        res_63 = P.compute_required_persistence(
            arms, horizon=63, cost_grid=(6.0, 12.0), alpha_grid=(1.0,), pulse_decay_tau=3.0
        )
    assert res_63.required_persistence_df["pr_X_le_48"].notna().all()
    pd.testing.assert_series_equal(
        res_60.required_persistence_df["pr_X_le_12"], res_63.required_persistence_df["pr_X_le_12"]
    )
    pd.testing.assert_series_equal(
        res_60.required_persistence_df["pr_X_le_24"], res_63.required_persistence_df["pr_X_le_24"]
    )


def test_malformed_inputs_raise_value_error(short_arms) -> None:
    effect = short_arms["AC"]
    bad_decay = (
        (effect, dict(pulse_decay_tau=3.0)),
        ({"BC": short_arms["BC"], "AB": short_arms["AB"]}, dict(pulse_decay_tau=3.0)),
        (dict(short_arms, extra=effect), dict(pulse_decay_tau=3.0)),
        (dict(short_arms, AC=short_arms["AC"] + 0.01), dict(pulse_decay_tau=3.0)),
        (short_arms, dict(pulse_decay_tau=3.0, plateau_arm="AC")),
        (short_arms, dict(pulse_decay_tau=0.0)),
        (short_arms, dict(pulse_decay_tau=float("nan"))),
        (short_arms, dict(survival=np.ones(12))),
        (short_arms, dict(survival=np.full(24, 1.2))),
        (short_arms, dict(survival=np.r_[np.full(6, 0.9), np.full(18, 0.5)])),
        (short_arms, dict(survival=np.r_[np.ones(6), np.full(18, np.nan)])),
    )
    for draws, options in bad_decay:
        with pytest.raises(ValueError):
            P.extrapolate_rank_k_decay(draws, horizon=24, rho_sd=None, **options)
    grids = dict(cost_grid=(1.0,), alpha_grid=(1.0,))
    bad_value = (
        (effect, dict(cost_grid=(), alpha_grid=(1.0,))),
        (effect, dict(cost_grid=(0.0,), alpha_grid=(1.0,))),
        (effect, dict(cost_grid=(1.0, 1.0), alpha_grid=(1.0,))),
        (effect, dict(cost_grid=(1.0,), alpha_grid=(1.5,))),
        (effect, dict(grids, alpha_share=1.2)),
        (effect, dict(grids, alpha_share=np.ones(5))),
        (effect, dict(grids, alpha_share=float("nan"))),
        (effect, dict(grids, rho_bounds=(0.9, 0.8))),
        (effect, dict(grids, rho_bounds=(0.5, 1.2))),
        (effect, dict(grids, rho_bounds=(0.5, 0.9, 1.0))),
        (effect, dict(grids, discount=np.ones(24), discount_rate=0.1)),
        (effect, dict(grids, discount=np.ones(10))),
        (effect, dict(grids, discount=np.zeros(24))),
        (effect, dict(grids, horizon=4)),
        (effect, dict(grids, checkpoints=(0,))),
        (effect, dict(grids, checkpoints=(1.5,))),
        (effect, dict(grids, checkpoints=None)),
        (effect, dict(grids, quantiles=(0.3, 0.7))),
        (effect, dict(grids, placebo_draws=effect[:, :5])),
        (short_arms, dict(grids, pulse_decay_tau=3.0, placebo_draws=effect)),
        ({"BC": short_arms["BC"], "AB": short_arms["AB"]}, grids),
        (np.full((10, 6), np.nan), grids),
    )
    for draws, options in bad_value:
        with pytest.raises(ValueError):
            P.compute_required_persistence(draws, **{"horizon": 24, **options})


# ---- Figure (smoke test under Agg) ---------------------------------------------------------------


def test_plot_required_persistence_and_payback(year_result, short_arms) -> None:
    axes = P.plot_required_persistence_and_payback(
        year_result, cost_highlight=3.0, alpha_highlight=1.0, unit="$", title="required persistence"
    )
    assert isinstance(axes, np.ndarray) and axes.shape == (3,)
    assert all(isinstance(ax, matplotlib.axes.Axes) for ax in axes)
    legend = [text.get_text() for text in axes[0].get_legend().get_texts()]
    assert any(CAVEAT in text for text in legend)
    annotations = [child.get_text() for child in axes[1].texts]
    assert any(text.startswith("X P50") for text in annotations)
    panel_b = [text.get_text() for text in axes[1].get_legend().get_texts()]
    assert any("Pr(X <= x)" in text for text in panel_b) and any("U: plateau rho = 1" in text for text in panel_b)
    _, given = plt.subplots(1, 3)
    out = P.plot_required_persistence_and_payback(year_result, cost_highlight=1.0, alpha_highlight=0.5, axes=given)
    assert all(o is g for o, g in zip(out, given, strict=True))
    with pytest.raises(ValueError, match="cost_highlight"):
        P.plot_required_persistence_and_payback(year_result, cost_highlight=2.0, alpha_highlight=1.0)
    with pytest.raises(ValueError, match="alpha_highlight"):
        P.plot_required_persistence_and_payback(year_result, cost_highlight=3.0, alpha_highlight=0.7)
    with pytest.raises(ValueError, match="expected 3 Axes"):
        P.plot_required_persistence_and_payback(year_result, cost_highlight=3.0, alpha_highlight=1.0, axes=given[:2])
    # Extra regimes stay in the tables; the figure keeps to Bounds U and L and the floor.
    with_extra = P.compute_required_persistence(
        short_arms, horizon=24, cost_grid=(1.0,), alpha_grid=(1.0,), pulse_decay_tau=3.0, extra_regimes={"P": 0.9}
    )
    extra_axes = P.plot_required_persistence_and_payback(with_extra, cost_highlight=1.0, alpha_highlight=1.0)
    for ax in extra_axes:
        assert not any(text.get_text().startswith("P") for text in ax.get_legend().get_texts())
    # The whole-level envelope and a horizon equal to L also draw.
    envelope = P.compute_required_persistence(short_arms["AC"], horizon=6, cost_grid=(1.0,), alpha_grid=(1.0,))
    assert envelope.report_months == (6,) and set(envelope.decision_df["regime"]) == {"observed"}
    envelope_axes = P.plot_required_persistence_and_payback(envelope, cost_highlight=1.0, alpha_highlight=1.0)
    assert envelope_axes.shape == (3,)
    panel_b = [text.get_text() for text in envelope_axes[1].get_legend().get_texts()]
    assert any("whole-level rho" in text for text in panel_b) and not any("plateau" in text for text in panel_b)


# ---- Integration smoke test on the golden fixture ------------------------------------------------


@pytest.fixture(scope="module")
def fitted_stag():
    panel = make_panel("stag_unequal")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", M.RankKConfigWarning)
        return M.fit_rank_k_longbet(
            panel["y"], panel["z"], panel["x"], cfg=CFG, baseline_forecast=panel["f"],
            baseline_forecast_draws=panel["f_draws"], cohorts=panel["cohorts"],
        )


def test_integration_from_the_golden_fixture(fitted_stag) -> None:
    rec = M.reconcile_nested_3arm_draws(fitted_stag, fitted_stag, fitted_stag, prior_scale=0.3)
    n_obs = rec.trt_catt_ac.shape[1]
    survival = P.weibull_survival_weights(36, n_obs, median_periods=24.0)
    decay = P.extrapolate_rank_k_decay(rec, horizon=36, rho=0.9, rho_sd=None, pulse_decay_tau=3.0, survival=survival)
    P.assert_nested_additivity(
        decay.cumulative_draws["BC"], decay.cumulative_draws["AB"], decay.cumulative_draws["AC"], atol=1e-10
    )
    scale = float(np.median(rec.trt_catt_ac.sum(axis=1)))
    assert scale > 0.0
    costs = (0.5 * scale, 1.5 * scale)
    placebo = {"BC": 0.1 * rec.trt_catt_bc, "AB": 0.1 * rec.trt_catt_ab, "AC": 0.1 * rec.trt_catt_ac}
    res = longbet.compute_required_persistence(
        rec, horizon=36, cost_grid=costs, alpha_grid=(1.0, 0.5), pulse_decay_tau=3.0, survival=survival,
        discount_rate=0.10, placebo_draws=placebo,
    )
    certificate_values = res.data_health_certificate_df.drop(columns=["never_code", "regime", "role"])
    assert res.n_observed == n_obs and not certificate_values.isna().any().any()
    assert res.placebo_summary_df is not None and res.arm_roles == {"BC": "plateau", "AB": "pulse", "AC": "total"}
    trajectory = longbet.event_time_trajectory(fitted_stag, weighting="unit")
    envelope = longbet.compute_required_persistence(
        trajectory, horizon=36, cost_grid=costs, alpha_grid=(1.0,), survival=None, discount_rate=0.10
    )
    assert envelope.n_observed == trajectory.shape[1]
    axes = longbet.plot_required_persistence_and_payback(res, cost_highlight=costs[1], alpha_highlight=0.5, unit="$")
    assert axes.shape == (3,)
