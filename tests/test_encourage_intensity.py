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

"""``LongBetEncourage(treatment="intensity")``: the two-stage distributed-lag model.

Fast tests use the closed-form ``"linear"`` Stage-2 backend; Stage 1 is the
real sampler on a tiny panel. The sampler-backed Stage 2 is exercised in
``test_encourage_intensity_sampler.py``.
"""

import inspect
import warnings

import numpy as np
import pandas as pd
import pytest

from longbet import (
    IntensityConditionalPrediction,
    IntensityEncouragementPrediction,
    LongBet,
    LongBetConfig,
    LongBetEncourage,
    lag_matrix,
)
from longbet._encourage_intensity import LinearLagStage2, control_function_columns, standardized_covariates

SMALL = dict(num_chains=2, num_burnin=4, num_sweeps=8, n_skip=1, num_trees_pr=2, num_trees_trt=2,
             max_depth_pr=3, max_depth_trt=3, min_points_per_leaf_pr=2, min_points_per_leaf_trt=2,
             random_seed=5)
TRUTH = np.array([-0.3, 0.5, 0.4, 0.2])
PULSE = np.array([1.0, 0.6, 0.3, 0.1])


def intensity_panel(seed=21, n=40, T=7, t0=4):
    """A small continuous-dose encouragement panel with a known lag kernel (``t0`` one-based)."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 2))
    z = np.zeros((n, T)); z[: n // 2, t0 - 1:] = 1
    u = rng.normal(size=(n, 1))
    offer = np.zeros((n, T)); offer[: n // 2, t0 - 1:] = PULSE[: T - t0 + 1]
    d = np.maximum(0.5 + 0.3 * rng.normal(size=(n, T)) + 0.5 * u + offer, 0)
    y = lag_matrix(d, 3) @ TRUTH + u + 0.5 * x[:, :1] + 0.3 * rng.normal(size=(n, T))
    strata = np.where(np.arange(n) % 2 == 0, "north", "south")
    return dict(y=y, d=d, z=z, x=x, t=np.arange(1.0, T + 1), strata=strata)


def fit_tier(identification, **overrides):
    p = intensity_panel()
    model = LongBetEncourage(treatment="intensity", identification=identification, stage2_backend="linear",
                             **{**SMALL, **overrides})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        model.fit(p["y"], p["d"], p["z"], p["x"], t=p["t"], strata=p["strata"])
    return model, p


@pytest.fixture(scope="module")
def fitted_t1():
    return fit_tier("control_function")


def test_constructor_validation():
    model = LongBetEncourage(treatment="intensity")
    assert model.identification == "control_function" and model.stage2_backend == "longbet"
    assert LongBetEncourage().treatment == "absorbing" and LongBetEncourage().identification is None
    with pytest.raises(ValueError, match="treatment must be"):
        LongBetEncourage(treatment="dose")
    with pytest.raises(ValueError, match="identification applies"):
        LongBetEncourage(identification="unconfounded")
    with pytest.raises(ValueError, match="apply to treatment='intensity'"):
        LongBetEncourage(lag_horizon=2)
    with pytest.raises(ValueError, match="identification must be one of"):
        LongBetEncourage(treatment="intensity", identification="iv")
    with pytest.raises(ValueError, match="continuous outcome only"):
        LongBetEncourage(treatment="intensity", outcome="binary")
    with pytest.raises(ValueError, match="lag_horizon must be"):
        LongBetEncourage(treatment="intensity", lag_horizon=-1)
    with pytest.raises(ValueError, match="stage2_backend"):
        LongBetEncourage(treatment="intensity", stage2_backend="stan")
    sig = inspect.signature(LongBetEncourage.fit)
    assert sig.parameters["strata"].default is None
    p = intensity_panel()
    with pytest.raises(ValueError, match="strata is used by"):
        LongBetEncourage(**SMALL).fit(p["y"], (p["d"] > 0.8).astype(float), p["z"], p["x"], strata=p["strata"])


def test_two_stages_are_fitted_as_designed(fitted_t1):
    model, p = fitted_t1
    assert model.dose_model.config.outcome == "continuous"
    assert model.dose_model.config.random_intercept and model.dose_model.config.ar1_errors
    assert model.dose_model.config.anchor_moderator
    assert isinstance(model.stage2, LinearLagStage2) and model.outcome_model is None
    meta = model.metadata
    assert meta["engine"] == "distributed_lag" and meta["inference_version"] == "distributed_lag_v1"
    assert meta["treatment"] == "intensity" and meta["tier"] == "T1"
    assert meta["lag_horizon"] == 3 and meta["n_strata"] == 2
    assert meta["component_streams"] == "independent" and isinstance(meta["engine_sha"], str)
    assert meta["dose_nonabsorbing"] and meta["encouragement_index"] == 3
    # The 3-term control function interacted with [1, x]: 3 * (1 + P) columns, never two terms.
    assert len(meta["extra_columns"]) == 3 * (1 + p["x"].shape[1])
    assert {"v", "v_lag", "v_bar"} <= set(meta["extra_columns"])
    assert model.stage2.kappa_draws.shape == (SMALL["num_chains"] * SMALL["num_sweeps"], 4)
    assert model.stage2.extra_coef_draws.shape[1] == len(meta["extra_columns"])


def test_tiers_build_the_right_stage2_design(fitted_t1):
    from longbet._encourage_intensity import organic_dose_columns
    model, p = fitted_t1
    d_hat, citt_d_mean = model._stage1_summaries(model._data)
    assert np.all(citt_d_mean[:, :3] == 0)
    x_tilde = standardized_covariates(p["x"])
    model.identification = "unconfounded"
    lags, extra, names = model._stage2_design(model._data, d_hat, citt_d_mean)
    assert lags.shape == (40, 7, 4) and extra is None and names == []
    np.testing.assert_allclose(lags[:, :, 0], p["d"])
    model.identification = "deconvolution"
    lags, extra, names = model._stage2_design(model._data, d_hat, citt_d_mean)
    assert lags.shape == (40, 7, 4) and np.all(lags[:, :3, :] == 0)
    np.testing.assert_allclose(lags[:20, 3:, 0], 0.5 * citt_d_mean[:20, 3:])
    np.testing.assert_allclose(lags[20:, 3:, 0], -0.5 * citt_d_mean[20:, 3:])
    expected_d0, expected_d0_names = organic_dose_columns(model._d0_hat, 3, x_tilde)
    np.testing.assert_allclose(extra, expected_d0)
    assert names == expected_d0_names and {"d0", "d0_lag"} <= set(names)
    model.identification = "control_function"
    lags, extra, names = model._stage2_design(model._data, d_hat, citt_d_mean)
    expected, expected_names = control_function_columns(p["d"] - d_hat, 3, x_tilde, 3)
    np.testing.assert_allclose(extra, expected)
    assert names == expected_names


def test_control_function_terms_cancel_exactly_in_the_offer_effect(fitted_t1):
    """The offer contrast is composed from kappa, nu_Y and CITT_D only: the control-function
    coefficients do not enter it at all (E[v | X, Z, t] = 0 in both arms)."""
    model, _ = fitted_t1
    before = model.predict()
    cond_before = model.predict_conditional()
    saved = model.stage2.extra_coef_draws
    try:
        model.stage2.extra_coef_draws = saved + 1e6 * np.random.default_rng(0).normal(size=saved.shape)
        after = model.predict()
        cond_after = model.predict_conditional()
    finally:
        model.stage2.extra_coef_draws = saved
    for name in ("itt_y", "itt_d", "omega_cum"):
        np.testing.assert_array_equal(after.draws[name], before.draws[name])
    np.testing.assert_array_equal(cond_after.citt_y.draws, cond_before.citt_y.draws)
    np.testing.assert_array_equal(cond_after.kernel.draws, cond_before.kernel.draws)


def test_stages_draw_independent_random_streams(monkeypatch):
    import jax
    from longbet import _encourage_model as em

    keys = []
    original = em.LongBet.fit

    def recording_fit(self, *args, key=None, **kwargs):
        keys.append(np.asarray(jax.random.key_data(key)))
        return original(self, *args, key=key, **kwargs)

    monkeypatch.setattr(em.LongBet, "fit", recording_fit)
    model, _ = fit_tier("unconfounded")
    key_y, key_d = model._split_keys(None, SMALL["random_seed"])
    assert len(keys) == 1                                       # the linear backend never calls the sampler
    np.testing.assert_array_equal(keys[0], np.asarray(jax.random.key_data(key_d)))
    assert not np.array_equal(np.asarray(jax.random.key_data(key_y)), keys[0])
    np.testing.assert_array_equal(np.asarray(jax.random.key_data(key_d)),
                                  np.asarray(jax.random.key_data(jax.random.fold_in(key_y, 1))))


def test_population_prediction_structure_and_reference(fitted_t1):
    model, p = fitted_t1
    pred = model.predict(groups=np.where(p["x"][:, 1] > 0, "hi", "lo"))
    assert isinstance(pred, IntensityEncouragementPrediction)
    assert pred.group_labels == ("all", "hi", "lo")
    for name in ("itt_y", "itt_d", "wald", "omega_cum"):
        assert pred.draws[name].shape == (3, 4, SMALL["num_chains"], SMALL["num_sweeps"])
    assert np.all(pred.draws["itt_d"][:, 0] != 0)
    np.testing.assert_allclose(
        pred.draws["omega_cum"],
        np.cumsum(pred.draws["itt_y"], axis=1) / np.cumsum(pred.draws["itt_d"], axis=1),
    )
    eff = pred.effects()
    assert set(eff.quantity) == {"itt_y", "itt_d", "wald", "omega_cum"}
    assert eff.loc[eff.quantity == "omega_cum", "posterior_mean"].isna().all()
    assert eff.loc[eff.quantity == "itt_y", "posterior_mean"].notna().all()
    ref = pred.reference
    assert len(ref) == 12 and {"itt_y_se", "itt_d_se", "itt_y_d_cov", "omega_cum_set_type", "condition_number",
                               "z_dip", "fosd_pvalue", "itt_y_cum_se", "n_units"} <= set(ref.columns)
    np.testing.assert_allclose(ref.loc[ref.group == "all", "itt_d"],
                               p["d"][:20, 3:].mean(axis=0) - p["d"][20:, 3:].mean(axis=0), atol=1e-12)
    assert pred.references["all"].metadata["n_strata"] == 2
    kernel = pred.kernel()
    assert set(kernel.quantity) == {"kappa", "kappa_acr"} and set(kernel.lag) == {0, 1, 2, 3}
    row = kernel.query("group == 'all' and quantity == 'kappa' and lag == 0").iloc[0]
    assert row.posterior_mean == pytest.approx(model.stage2.kappa_draws[:, 0].mean())
    assert np.isfinite(row.reference_kappa) and np.isfinite(row.reference_kappa_cum_se)
    disc = pred.discrepancy()
    assert set(disc.label) == {"T1"} and len(disc) == 12 and np.isfinite(disc.z_stat).all()
    diag = pred.diagnostics()
    assert list(diag.group) == ["all", "hi", "lo"] and np.isfinite(diag.condition_number).all()
    cum = pred.cumulative_return()
    assert len(cum) == 12 and {"omega_cum_lower_1", "posterior_median"} <= set(cum.columns)
    assert pred.metadata["engine"] == "distributed_lag" and "omega_cum" in pred.metadata["ratio_interpretation"]
    assert pred.first_stage().posterior_probability_positive.between(0, 1).all()
    comp = pred.comparison()
    assert comp.difference.notna().all() and set(comp.quantity) == {"itt_y", "itt_d"}
    cum_disc = pred.cumulative_discrepancy()
    assert np.isfinite(cum_disc.z_stat).all()
    with pytest.raises(ValueError, match="block_size"):
        model.predict(block_size=0)
    same = model.predict(block_size=7)
    np.testing.assert_allclose(same.draws["itt_y"], model.predict().draws["itt_y"])


def test_missing_group_arm_is_retained(fitted_t1):
    model, p = fitted_t1
    # Units 0..19 are encouraged: 'both' keeps 5 encouraged and all 20 controls, the other group has no control.
    groups = np.array(["only_encouraged"] * 20 + ["both"] * 20)
    groups[:5] = "both"
    pred = model.predict(groups=groups)
    ref = pred.reference
    assert (ref.loc[ref.group == "only_encouraged", "wald_reason"] == "missing_group_arm").all()
    assert ref.loc[ref.group == "only_encouraged", "itt_y"].isna().all()
    assert pred.references["only_encouraged"] is None
    assert pred.effects().query("group == 'only_encouraged'").model_extrapolation.all()
    assert pred.diagnostics().set_index("group").loc["only_encouraged"].isna().all()


def test_conditional_prediction_in_and_out_of_sample(fitted_t1):
    model, p = fitted_t1
    cond = model.predict_conditional()
    assert isinstance(cond, IntensityConditionalPrediction)
    n_draws = SMALL["num_chains"] * SMALL["num_sweeps"]
    assert cond.citt_y.draws.shape == (40, 4, n_draws) and cond.citt_d.draws.shape == (40, 4, n_draws)
    assert cond.kernel.draws.shape == (40, 4, n_draws) and cond.kernel_cum.draws.shape == (40, 4, n_draws)
    # The linear backend has nu_Y = 1, so every unit carries the population kernel.
    np.testing.assert_allclose(cond.kernel.draws[0], model.stage2.kappa_draws.T)
    np.testing.assert_allclose(cond.kernel_cum.draws, np.cumsum(cond.kernel.draws, axis=1))
    # Composition identity per unit and draw: CITT_Y(t) = sum_l kappa_l CITT_D(t - l).
    for b in (0, n_draws - 1):
        for i in (0, 39):
            expected = [sum(model.stage2.kappa_draws[b, lag] * cond.citt_d.draws[i, h - lag, b]
                            for lag in range(h + 1)) for h in range(4)]
            np.testing.assert_allclose(cond.citt_y.draws[i, :, b], expected, rtol=1e-6, atol=1e-8)
    assert cond.cace is cond.omega_cum and cond.omega_cum.draws.shape == (40, 4, n_draws)
    assert cond.always_takers is None and cond.p0_draws is None
    with pytest.raises(ValueError, match="continuous dose"):
        cond.principal_strata()
    with pytest.raises(ValueError, match="continuous dose"):
        cond.strata_summary()
    summary = cond.kernel_summary()
    assert list(summary.columns[:2]) == ["lag", "kappa_mean"] and len(summary) == 4
    # Decision layer is inherited unchanged.
    assert cond.optimal_policy(cost=0.1).dtype == bool
    assert cond.knapsack_policy(cost=0.1, capacity=3).sum() <= 3
    assert np.isfinite(cond.policy_value(cost=0.1))
    # The population ITT equals the equal-weight average of the conditional effects of the study units.
    pop = model.predict()
    np.testing.assert_allclose(pop.draws["itt_y"][0].reshape(4, -1), cond.citt_y.draws.mean(axis=0), rtol=1e-6)
    # New units with a counterfactual offer schedule and calendar.
    new = model.predict_conditional(x=np.zeros((3, 2)), z=np.array([[0, 0, 1, 1, 1, 1, 1]] * 3), t=p["t"])
    assert new.citt_y.draws.shape == (3, 5, n_draws) and list(new.horizons) == [1, 2, 3, 4, 5]
    assert np.all(new.citt_d.draws[:, 0, :] != 0)
    blocked = model.predict_conditional(block_size=7)
    np.testing.assert_allclose(blocked.citt_y.draws, cond.citt_y.draws)
    with pytest.raises(ValueError, match="x must be a 2-D array"):
        model.predict_conditional(x=np.zeros((3, 5)))


def test_lag_horizon_and_pre_offer_demeaning():
    model, p = fit_tier("unconfounded", lag_horizon=1)
    assert model.metadata["lag_horizon"] == 1 and model.stage2.kappa_draws.shape[1] == 2
    assert model.predict_conditional().kernel.draws.shape[1] == 2
    model, p = fit_tier("unconfounded", demean_pre=True)
    d_hat, citt = model._stage1_summaries(model._data)
    lags, _, _ = model._stage2_design(model._data, d_hat, citt)
    np.testing.assert_allclose(lags[:, :, 0], p["d"] - p["d"][:, :3].mean(axis=1, keepdims=True))
    with pytest.raises(ValueError, match="lag_horizon cannot exceed"):
        LongBetEncourage(treatment="intensity", lag_horizon=4, stage2_backend="linear", **SMALL).fit(
            p["y"], p["d"], p["z"], p["x"])


def test_roundtrip_archive(fitted_t1, tmp_path):
    model, p = fitted_t1
    path = tmp_path / "intensity.npz"
    model.save(path)
    loaded = LongBetEncourage.load(path)
    assert loaded.treatment == "intensity" and loaded.identification == "control_function"
    assert loaded.stage2_backend == "linear" and loaded.metadata == model.metadata
    assert "strata" in loaded._data and loaded._data["strata"].dtype.kind == "U"
    before, after = model.predict(), loaded.predict()
    for name in ("itt_y", "itt_d", "omega_cum"):
        np.testing.assert_array_equal(after.draws[name], before.draws[name])
    pd.testing.assert_frame_equal(after.reference, before.reference)
    np.testing.assert_array_equal(loaded.predict_conditional().citt_y.draws, model.predict_conditional().citt_y.draws)
    # Tampering with an array is caught by the digests.
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    arrays["d"] = arrays["d"] * 1.01
    tampered = tmp_path / "tampered.npz"
    with tampered.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    with pytest.raises(ValueError, match="digest mismatch"):
        LongBetEncourage.load(tampered)


def test_bootstrap_comparison_refits_both_stages(fitted_t1):
    model, p = fitted_t1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = model.bootstrap_comparison(replicates=2, random_seed=1)
    assert result.replicates.shape == (2, 1, 4, 2, 2)
    assert set(result.table.quantity) == {"itt_y", "itt_d"}
    assert result.metadata["resampling_unit"] == "entire_unit_record"
    assert np.isfinite(result.replicates).all()
    assert result.cumulative_discrepancy().disagreement.isna().all()      # fewer than 100 replicates


def test_absorbing_path_is_unchanged_by_the_intensity_extension():
    """``treatment="absorbing"`` must call the sampler exactly as before: both equations' draws equal
    direct LongBet fits with the same configurations and keys."""
    import jax
    rng = np.random.default_rng(918)
    n, T = 12, 4
    x = rng.normal(size=(n, 2))
    z = np.zeros((n, T)); z[: n // 2, 1:] = 1
    d = np.zeros((n, T)); d[: 4, 1:] = 1; d[4:7, 2:] = 1; d[8, :] = 1
    y = rng.normal(size=(n, T)) + np.cumsum(d, axis=1) + x[:, :1]
    model = LongBetEncourage(**SMALL)
    assert model.treatment == "absorbing" and model.identification is None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        model.fit(y, d, z, x)
        outcome_cfg, adoption_cfg = model._equation_configs()
        key_y = jax.random.key(outcome_cfg.random_seed)
        key_d = jax.random.fold_in(key_y, 1)
        direct_outcome = LongBet(outcome_cfg).fit(y, x, d, key=key_y)
        direct_adoption = LongBet(adoption_cfg).fit(LongBetEncourage._first_adoption_events(d), x, z, key=key_d)
    np.testing.assert_array_equal(np.asarray(model.outcome_model.trace.beta), np.asarray(direct_outcome.trace.beta))
    np.testing.assert_array_equal(np.asarray(model.outcome_model.trace.gamma), np.asarray(direct_outcome.trace.gamma))
    np.testing.assert_array_equal(np.asarray(model.adoption_model.trace.beta), np.asarray(direct_adoption.trace.beta))
    assert model.metadata["engine"] == "adoption_clock" and "treatment" not in model.metadata
    pred = model.predict()
    assert type(pred).__name__ == "EncouragementPrediction" and set(pred.draws) == {"itt_y", "itt_d", "wald"}
    assert type(model.predict_conditional()).__name__ == "ConditionalEncouragementPrediction"
