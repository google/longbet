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

"""The estimator API: fit, predict, summaries, persistence."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, get_att, get_catt
from _panels import make_staggered_panel

# Deliberately tiny, and explicitly single-chain: these tests assert the
# single-chain trace layout. Multi-chain is covered separately below.
FAST = dict(num_sweeps=12, num_burnin=6, num_trees_pr=6, num_trees_trt=6,
            num_chains=1, random_seed=42)


def _fit(data, **cfg):
    model = LongBet(LongBetConfig(**{**FAST, **cfg}))
    return model.fit(y=data["y"], x=data["x"], z=data["z"], t=data["t"])


def test_fit_and_predict_shapes(staggered_panel):
    data = staggered_panel
    model = _fit(data)
    N, T = data["y"].shape

    assert model.trace.beta.shape == (FAST["num_sweeps"], model.S_max_ + 1)
    assert model.trace.gamma.shape == (FAST["num_sweeps"], N)

    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])
    for arr in (pred.tauhats, pred.muhats0, pred.yhats):
        assert arr.shape == (N, T, FAST["num_sweeps"])
        assert np.all(np.isfinite(arr))
    for summary in (pred.tau_summary, pred.mu0_summary, pred.y_summary):
        assert summary.mean.shape == (N, T)
        assert np.all(summary.lower <= summary.upper)


def test_summary_only_matches_full_draws_and_keeps_the_att(staggered_panel):
    """summary_only must be a memory strategy, not a different answer."""
    data = staggered_panel
    model = _fit(data)
    full = model.predict(x=data["x"], z=data["z"], t=data["t"], summary_only=False)
    lean = model.predict(x=data["x"], z=data["z"], t=data["t"], summary_only=True)

    assert lean.tauhats is None and lean.muhats0 is None and lean.yhats is None
    for a, b in ((full.tau_summary, lean.tau_summary), (full.mu0_summary, lean.mu0_summary)):
        for name in ("mean", "std", "lower", "upper"):
            np.testing.assert_allclose(getattr(a, name), getattr(b, name), rtol=1e-5, atol=1e-6)

    # The ATT is small, so it survives summary_only; get_att and stability work.
    np.testing.assert_allclose(full.att_full, lean.att_full, rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(get_att(full)["att"], get_att(lean)["att"], rtol=1e-5, atol=1e-6)
    assert model.stability(lean, warn=False).summary["total_draws"] == FAST["num_sweeps"]


def test_block_size_does_not_change_the_answer(staggered_panel):
    data = staggered_panel
    model = _fit(data)
    a = model.predict(x=data["x"], z=data["z"], t=data["t"], summary_only=True)
    b = model.predict(x=data["x"], z=data["z"], t=data["t"], summary_only=True, block_size=7)
    np.testing.assert_allclose(a.tau_summary.mean, b.tau_summary.mean, rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(a.tau_summary.lower, b.tau_summary.lower, rtol=1e-5, atol=1e-6)


def test_untreated_outcome_includes_the_unit_intercept(staggered_panel):
    """muhats0 is a counterfactual *outcome*, so it must carry gamma_i.

    Without it the per-unit level is systematically wrong in a model whose
    headline feature is unit random intercepts.
    """
    data = staggered_panel
    model = _fit(data, num_sweeps=25, num_burnin=15, num_trees_pr=10, num_trees_trt=10)
    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])

    never_treated = data["z"].sum(axis=1) == 0
    assert never_treated.sum() > 5

    # On never-treated units the untreated outcome is the observed outcome, so
    # dropping gamma_i must make the per-unit fit strictly worse. Comparing the
    # two directly avoids leaning on an arbitrary correlation threshold.
    gamma_hat = np.asarray(model.trace.gamma).mean(axis=0)
    with_gamma = pred.mu0_summary.mean[never_treated]
    without_gamma = with_gamma - gamma_hat[never_treated, None]
    observed = data["y"][never_treated]

    rmse = lambda a: float(np.sqrt(np.mean((a - observed) ** 2)))
    assert rmse(with_gamma) < rmse(without_gamma), (
        f"muhats0 does not include the unit intercept: "
        f"RMSE with gamma {rmse(with_gamma):.3f} vs without {rmse(without_gamma):.3f}"
    )
    assert np.abs(with_gamma - observed).mean() < 0.5


def test_factual_predictions_track_the_outcome(staggered_panel):
    data = staggered_panel
    model = _fit(data, num_sweeps=25, num_burnin=15, num_trees_pr=10, num_trees_trt=10)
    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])
    assert np.corrcoef(pred.y_summary.mean.ravel(), data["y"].ravel())[0, 1] > 0.9


def test_unbalanced_panel_is_marginalized_not_imputed(unbalanced_panel):
    data = unbalanced_panel
    assert np.isnan(data["y"]).any()
    model = _fit(data, num_sweeps=20, num_burnin=10)
    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])
    assert np.all(np.isfinite(pred.tauhats))
    assert np.all(np.isfinite(np.asarray(model.trace.sigma2)))


def test_att_and_catt_helpers_agree_with_methods(staggered_panel):
    data = staggered_panel
    model = _fit(data)
    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])
    np.testing.assert_allclose(get_att(pred)["att"], pred.att()["att"])
    np.testing.assert_allclose(get_catt(pred)["catt"], pred.catt()["catt"])
    assert get_att(pred)["exposure"].tolist() == list(range(1, pred.att_full.shape[0] + 1))


def test_save_load_round_trip(staggered_panel):
    data = staggered_panel
    model = _fit(data)
    before = model.predict(x=data["x"], z=data["z"], t=data["t"])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "model.npz"
        model.save(path)
        reloaded = LongBet.load(path)
        after = reloaded.predict(x=data["x"], z=data["z"], t=data["t"])
    np.testing.assert_allclose(before.tauhats, after.tauhats, atol=1e-5)
    np.testing.assert_allclose(before.muhats0, after.muhats0, atol=1e-5)
    np.testing.assert_allclose(before.att_full, after.att_full, atol=1e-5)


def test_multichain_flattens_chains_into_draws(staggered_panel):
    data = staggered_panel
    model = _fit(data, num_chains=3, num_sweeps=10)
    assert model.trace.beta.shape[0] == 3
    assert model.state.X.ndim == 2, "the design matrix must not be replicated per chain"
    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])
    assert pred.tauhats.shape[-1] == 3 * 10
    stab = model.stability(pred, warn=False)
    assert stab.summary["num_chains"] == 3
    assert np.isfinite(stab.summary["rhat_max"])


def test_stability_withholds_the_verdict(staggered_panel):
    data = staggered_panel
    model = _fit(data, num_chains=2, num_sweeps=10)
    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])
    stab = model.stability(pred, warn=False)
    assert stab.reliable is None
    assert stab.summary["verdict"] == "withheld"
    assert "coverage" in stab.summary["verdict_note"]
    assert isinstance(stab.summary["ess_ok"], bool)


def test_single_chain_warns_that_rhat_needs_two(staggered_panel):
    data = staggered_panel
    model = _fit(data)
    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])
    with pytest.warns(UserWarning, match="at least two chains"):
        model.stability(pred, min_ess=1.0, warn=True)


def test_forecast_beyond_the_fitted_horizon(staggered_panel):
    """Predicting past the fitted exposure horizon extends beta by the GP."""
    data = staggered_panel
    model = _fit(data)
    N, T = data["y"].shape

    T2 = T + 4
    z2 = np.zeros((N, T2), dtype=np.float32)
    z2[:, :T] = data["z"]
    z2[:, T:] = data["z"][:, [-1]]
    t2 = np.arange(1, T2 + 1, dtype=np.float32)

    pred = model.predict(x=data["x"], z=z2, t=t2)
    assert pred.beta_values.shape[1] > model.S_max_ + 1
    assert np.all(np.isfinite(pred.tauhats))

    # Uncertainty must widen past the fitted horizon: beyond it the trajectory
    # is the GP prior conditioned on the fitted block, not a reading of data.
    width = pred.tau_summary.upper - pred.tau_summary.lower
    inside = (pred.s >= 1) & (pred.s <= model.S_max_)
    beyond = pred.s > model.S_max_
    assert beyond.any(), "the test panel must actually run past the fitted horizon"
    assert width[beyond].mean() > width[inside].mean(), (
        f"projected interval width {width[beyond].mean():.3f} is not wider than "
        f"the in-sample {width[inside].mean():.3f}"
    )


def test_projection_is_reproducible_from_the_seed(staggered_panel):
    data = staggered_panel
    model = _fit(data)
    N, T = data["y"].shape
    T2 = T + 3
    z2 = np.zeros((N, T2), dtype=np.float32)
    z2[:, :T] = data["z"]
    z2[:, T:] = data["z"][:, [-1]]
    t2 = np.arange(1, T2 + 1, dtype=np.float32)
    a = model.predict(x=data["x"], z=z2, t=t2)
    b = model.predict(x=data["x"], z=z2, t=t2)
    np.testing.assert_array_equal(a.beta_values, b.beta_values)


def test_non_absorbing_treatment_is_refused():
    data = make_staggered_panel(N=10, T=5)
    z = data["z"].copy()
    z[0] = [0, 1, 1, 0, 0]
    with pytest.raises(ValueError, match="absorbing"):
        LongBet(LongBetConfig(**FAST)).fit(y=data["y"], x=data["x"], z=z, t=data["t"])


def test_always_treated_units_warn():
    data = make_staggered_panel(N=20, T=5)
    z = data["z"].copy()
    z[0] = 1.0
    with pytest.warns(UserWarning, match="treated in every period"):
        LongBet(LongBetConfig(**FAST)).fit(y=data["y"], x=data["x"], z=z, t=data["t"])


def test_binary_outcome_uses_the_probit_offset(binary_panel):
    """Regression: the probit intercept used to be computed and then dropped."""
    data = binary_panel
    model = LongBet(LongBetConfig(outcome="binary", num_sweeps=20, num_burnin=12,
                                  num_trees_pr=8, num_trees_trt=6, random_seed=3))
    model.fit(y=data["y"], x=data["x"], z=data["z"], t=data["t"])

    from scipy.stats import norm
    rate = float(data["y"].mean())
    assert abs(model.offset_ - norm.ppf(rate)) < 1e-5
    # The offset reaches bartz and comes back through evaluate_trace.
    assert np.allclose(np.asarray(model.trace.mu_trace.offset), model.offset_, atol=1e-4)
    # sigma^2 is held at 1 for a probit model.
    assert np.allclose(np.asarray(model.trace.sigma2), 1.0)

    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])
    implied = norm.cdf(pred.y_summary.mean)
    assert abs(implied.mean() - rate) < 0.15, (
        f"implied event rate {implied.mean():.3f} vs observed {rate:.3f}"
    )


def test_binary_outcome_rejects_non_binary_y(staggered_panel):
    data = staggered_panel
    with pytest.raises(ValueError, match="requires y in"):
        LongBet(LongBetConfig(outcome="binary", **FAST)).fit(
            y=data["y"], x=data["x"], z=data["z"], t=data["t"]
        )


def test_predict_before_fit_is_refused(staggered_panel):
    with pytest.raises(RuntimeError, match="fit"):
        LongBet(LongBetConfig(**FAST)).predict(
            x=staggered_panel["x"], z=staggered_panel["z"]
        )


def test_config_validates_its_arguments():
    with pytest.raises(ValueError, match="num_chains"):
        LongBetConfig(num_chains=0)
    with pytest.raises(ValueError, match="kernel_type"):
        LongBetConfig(kernel_type="rbf")
    with pytest.raises(ValueError, match="device"):
        LongBetConfig(device="tpu")
    with pytest.raises(ValueError, match="lambda_knl"):
        LongBetConfig(lambda_knl=0.0)


def test_stability_reports_the_support_behind_each_exposure_time(staggered_panel):
    """A bad R-hat at an exposure time backed by 20 units means something
    different from one backed by 2,000, so the support is reported alongside."""
    data = staggered_panel
    model = _fit(data, num_chains=2, num_sweeps=10)
    pred = model.predict(x=data["x"], z=data["z"], t=data["t"])

    counts = pred.att_counts
    assert counts is not None
    assert len(counts) == pred.att_full.shape[0]
    # A staggered rollout is reached by fewer units at later exposure times.
    assert counts[0] >= counts[-1]
    # Every exposure time with a finite ATT must have at least one treated cell.
    finite = np.isfinite(pred.att_full).any(axis=1)
    assert np.all(counts[finite] > 0)

    stab = model.stability(pred, warn=False)
    np.testing.assert_array_equal(stab.by_exposure["n_treated"], counts)
    assert stab.summary["min_treated_cells"] == float(counts.min())
    # The worst R-hat is located, so the reader can see whether it sits where
    # the data is thin.
    worst = stab.summary["rhat_max_at_exposure"]
    assert 1 <= worst <= len(counts)
    assert stab.summary["rhat_max_n_treated"] == float(counts[worst - 1])


def test_rhat_warning_names_where_the_problem_is(staggered_panel):
    """'R-hat is 1.09' is not actionable; 'at exposure time 17, backed by 40
    cells' tells the reader whether to sample longer or to stop asking."""
    rng = np.random.default_rng(0)
    # Three chains that disagree badly at the last, thinly supported exposure.
    att = rng.normal(size=(3, 200, 5))
    att[:, :, 4] += np.array([0.0, 3.0, -3.0])[:, None]

    from longbet import att_stability

    with pytest.warns(UserWarning, match=r"exposure time 5, which 25 treated"):
        result = att_stability(
            att, min_ess=1.0, warn=True, n_treated=np.array([900, 700, 400, 120, 25])
        )
    assert result.summary["rhat_max_at_exposure"] == 5
    assert result.summary["rhat_max_n_treated"] == 25.0
    assert result.summary["rhat_ok"] is False
