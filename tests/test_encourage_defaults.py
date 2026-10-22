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

"""The encouragement model's fixed components and the model-dependent
defaults of LongBetConfig: resolution, refusals, removed options, archives
and the tempering exchange warning."""

import json
import warnings

import numpy as np
import pytest
from scipy.special import ndtr

from longbet import LongBet, LongBetConfig, LongBetEncourage
from longbet._config import MODEL_DEFAULTS
from longbet._encourage_model import DEFAULT_N_SKIP, MODEL_COMPONENTS, _flat_draws, _predict_quiet
from longbet._tempering import SWAP_RATE_FLOOR, swap_rates, warn_if_not_exchanging

from test_encourage_model import encouragement_panel

SMALL = dict(num_chains=2, num_burnin=4, num_sweeps=6, num_trees_pr=2, num_trees_trt=2,
             max_depth_pr=3, max_depth_trt=3, min_points_per_leaf_pr=2,
             min_points_per_leaf_trt=2, random_seed=5)
MODEL = {**{name: True for name in MODEL_COMPONENTS}, "ar1_errors": True, "n_skip": DEFAULT_N_SKIP}


def _resolved(config):
    return {name: getattr(config, name) for name in MODEL_DEFAULTS}


def test_unset_options_are_left_to_the_model():
    assert LongBetConfig().unresolved == tuple(MODEL_DEFAULTS)
    assert _resolved(LongBet().config) == MODEL_DEFAULTS
    assert _resolved(LongBet(LongBetConfig()).config) == MODEL_DEFAULTS
    assert LongBet().config.unresolved == ()


def test_the_encouragement_model_fills_in_its_components():
    for model in (LongBetEncourage(), LongBetEncourage(LongBetConfig()),
                  LongBetEncourage(config=LongBetConfig(num_chains=2))):
        assert _resolved(model.config) == MODEL
    assert LongBetEncourage().config == LongBetEncourage.model_config()
    binary = LongBetEncourage(outcome="binary").config
    assert _resolved(binary) == {**MODEL, "ar1_errors": False}
    assert LongBetEncourage.model_config("binary") == binary
    assert LongBetEncourage(n_skip=1).config.n_skip == 1


@pytest.mark.parametrize("override", [dict(linear_trt=False), dict(learn_beta_scale=False),
                                      dict(interweave_gamma=False)])
def test_a_component_cannot_be_switched_off(override):
    with pytest.raises(ValueError, match="always uses"):
        LongBetEncourage(**override)
    with pytest.raises(ValueError, match="always uses"):
        LongBetEncourage(LongBetConfig(**override))


def test_ar1_errors_follow_the_outcome_and_need_their_prerequisites():
    with pytest.raises(ValueError, match="not for a binary one"):
        LongBetEncourage(outcome="binary", ar1_errors=True)
    with pytest.raises(ValueError, match="not for a binary one"):
        LongBetEncourage(ar1_errors=False)
    with pytest.raises(ValueError, match="unit intercepts"):
        LongBetEncourage(random_intercept=False)
    for override in (dict(split_exposure_trt=True), dict(split_calendar_trt=True),
                     dict(use_trend_block=False)):
        with pytest.raises(ValueError, match="unit-level partitions"):
            LongBetEncourage(**override)
    assert not LongBetEncourage(outcome="binary", split_exposure_trt=True).config.ar1_errors


@pytest.mark.parametrize("name, kept", [("linear_trt_horseshoe", True),
                                        ("linear_trt_outcome_only", False),
                                        ("adoption_frailty", False)])
def test_removed_options(name, kept):
    with pytest.warns(DeprecationWarning, match=name):
        config = LongBetConfig(**{name: kept})
    assert not hasattr(config, name)
    with pytest.raises(ValueError, match="no longer supported"):
        LongBetConfig(**{name: not kept})
    with pytest.raises(ValueError, match="no longer supported"):
        LongBetEncourage(**{name: not kept})


def test_a_configuration_saved_before_the_extensions_is_the_published_model():
    old = {k: v for k, v in LongBetConfig(num_chains=2, n_skip=2).to_dict().items()
           if k not in MODEL_DEFAULTS or k == "n_skip"}
    config = LongBetConfig.from_dict({**old, "adoption_frailty": True})
    assert _resolved(config) == {**MODEL_DEFAULTS, "n_skip": 2}
    assert LongBetConfig.from_dict({**old, "linear_trt_horseshoe": True}) == config


def test_swap_rates_and_the_exchange_warning():
    # Two chains of a three-level ladder over 11 sweeps: level 0 leads on the
    # six even sweeps, level 1 on the five odd ones; the top level never leads.
    counts = np.array([3.0, 0.0, 0.0, 6.0, 5.0, 0.0])
    rates = swap_rates(counts, num_chains=2, levels=3, total_sweeps=11)
    np.testing.assert_allclose(rates, [[0.5, 0.0], [1.0, 1.0]])
    with pytest.warns(UserWarning, match="levels 1 and 2 of chain 0"):
        warn_if_not_exchanging(rates)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warn_if_not_exchanging(np.full((2, 2), SWAP_RATE_FLOOR))


def test_a_tempered_fit_reports_its_exchange_rates():
    p = encouragement_panel()
    cfg = dict(SMALL, tempering_levels=2, tempering_beta_min=0.5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = LongBet(**cfg).fit(p["y"], p["x"], p["d"], t=p["t"])
        plain = LongBet(**SMALL).fit(p["y"], p["x"], p["d"], t=p["t"])
    assert fit.swap_rates_.shape == (2, 1)
    assert np.all((fit.swap_rates_ >= 0) & (fit.swap_rates_ <= 1))
    assert plain.swap_rates_ is None


@pytest.fixture(scope="module")
def model_fit():
    p = encouragement_panel()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return LongBetEncourage(**SMALL).fit(p["y"], p["d"], p["z"], p["x"], t=p["t"])


def test_the_fit_runs_every_component(model_fit):
    out, adopt = model_fit.outcome_model, model_fit.adoption_model
    assert out.config.ar1_errors and not adopt.config.ar1_errors
    assert out.config.random_intercept and not adopt.config.random_intercept
    for model in (out, adopt):
        assert model.trace.delta is not None and model.trace.beta_scale2 is not None
    assert out.trace.ar1_phi is not None
    assert out.trace.beta.shape[1] == SMALL["num_sweeps"]


def test_the_fit_round_trips(model_fit, tmp_path):
    path = tmp_path / "fit.npz"
    model_fit.save(path)
    loaded = LongBetEncourage.load(path)
    assert loaded.config == model_fit.config
    for name in ("outcome_model", "adoption_model"):
        a, b = getattr(model_fit, name), getattr(loaded, name)
        for key in ("mean", "sd", "keep"):
            np.testing.assert_array_equal(b.linear_trt_spec_[key], a.linear_trt_spec_[key])
        for key in ("delta", "tau_delta2", "beta_scale2"):
            np.testing.assert_array_equal(np.asarray(getattr(b.trace, key)),
                                          np.asarray(getattr(a.trace, key)))
    x_new = encouragement_panel(seed=3)["x"][:4]
    before = model_fit.predict_conditional(x=x_new)
    after = loaded.predict_conditional(x=x_new)
    np.testing.assert_array_equal(after.citt_y.draws, before.citt_y.draws)
    np.testing.assert_array_equal(after.citt_d.draws, before.citt_d.draws)


def test_an_archive_of_the_earlier_model_is_refused(model_fit, tmp_path):
    model_fit.save(tmp_path / "fit.npz")
    with np.load(tmp_path / "fit.npz", allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    meta = json.loads(str(arrays["metadata"]))
    meta["inference_version"] = "adoption_clock_v1"
    arrays["metadata"] = np.asarray(json.dumps(meta))
    np.savez(tmp_path / "old.npz", **arrays)
    with pytest.raises(ValueError, match="earlier release"):
        LongBetEncourage.load(tmp_path / "old.npz")


def test_a_binary_outcome_intercept_is_integrated_exactly():
    """The closed form against a 120-node Gauss-Hermite rule on the fitted draws."""
    p = encouragement_panel()
    y = (p["y"] > np.median(p["y"])).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = LongBetEncourage(outcome="binary", **SMALL).fit(y, p["d"], p["z"], p["x"], t=p["t"])
    from longbet._encourage_model import adoption_distribution, offer_effect_on_outcome
    model = fit.outcome_model
    x, t = p["x"][:3], p["t"]
    n, periods = 3, len(t)
    p1, _ = adoption_distribution(fit.adoption_model, x, np.ones((n, periods)), t)
    p0, _ = adoption_distribution(fit.adoption_model, x, np.zeros((n, periods)), t)
    exact = offer_effect_on_outcome(model, x, t, p1, p0)
    nodes, weights = np.polynomial.hermite_e.hermegauss(120)
    weights = weights / weights.sum()
    sd = np.sqrt(_flat_draws(model, "sigma_gamma2"))[None, None, :]
    rule = np.zeros_like(exact)
    for a in range(periods):
        schedule = np.zeros((n, periods), dtype=np.float32)
        schedule[:, a:] = 1.0
        pred = _predict_quiet(model, x, schedule, t, None)
        tau, mu0 = np.asarray(pred.tauhats, np.float64), np.asarray(pred.muhats0, np.float64)
        for node, w in zip(nodes, weights):
            rule += w * (p1 - p0)[:, a, None, :] * (ndtr(mu0 + tau + node * sd) - ndtr(mu0 + node * sd))
    np.testing.assert_allclose(exact, rule, rtol=0, atol=1e-7)
