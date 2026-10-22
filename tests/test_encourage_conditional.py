"""Conditional encouragement prediction: schedules, horizons and the marginal target."""

import numpy as np
import pytest

from longbet import LongBetEncourage
from longbet._encourage_model import ConditionalEncouragementPrediction
from test_encourage_model import encouragement_panel


@pytest.fixture(scope="module")
def fit_and_panel():
    p = encouragement_panel(seed=31)
    fit = LongBetEncourage(num_chains=2, num_burnin=4, num_sweeps=6, n_skip=1, num_trees_pr=2,
                           num_trees_trt=2, random_seed=9)
    return fit.fit(p["y"], p["d"], p["z"], p["x"], t=p["t"]), p


def test_default_schedule_offers_from_the_launch_period(fit_and_panel):
    fit, p = fit_and_panel
    default = fit.predict_conditional()
    explicit = fit.predict_conditional(z=np.where(np.arange(4)[None, :] >= 1, 1.0, 0.0) * np.ones((12, 1)))
    np.testing.assert_allclose(default.citt_y.draws, explicit.citt_y.draws)
    np.testing.assert_array_equal(default.periods, p["t"][1:])
    np.testing.assert_array_equal(default.horizons, p["t"][1:] - p["t"][1] + 1)
    assert isinstance(default, ConditionalEncouragementPrediction)


def test_never_offering_gives_no_effect(fit_and_panel):
    fit, p = fit_and_panel
    none = fit.predict_conditional(z=np.zeros((12, 4)))
    np.testing.assert_allclose(none.citt_y.draws, 0, atol=1e-12)
    np.testing.assert_allclose(none.citt_d.draws, 0, atol=1e-12)


def test_marginal_target_is_the_same_for_study_and_new_units(fit_and_panel):
    """A study unit and a copy of its covariates get identical marginal effects."""
    fit, p = fit_and_panel
    study = fit.predict_conditional()
    copy = fit.predict_conditional(x=np.asarray(p["x"][:3]))
    # Forest sums are accumulated in float32 in an order that depends on the batch of
    # units evaluated, so agreement is to single precision, not bit-for-bit.
    np.testing.assert_allclose(study.citt_y.draws[:3], copy.citt_y.draws, rtol=1e-5, atol=1e-7)
    np.testing.assert_allclose(study.citt_d.draws[:3], copy.citt_d.draws, rtol=1e-5, atol=1e-7)


def test_strata_are_probabilities_that_sum_to_one(fit_and_panel):
    fit, _ = fit_and_panel
    cond = fit.predict_conditional()
    total = cond.always_takers.draws + cond.compliers.draws + cond.never_takers.draws
    np.testing.assert_allclose(total, 1.0, atol=1e-6)
    assert np.all(cond.always_takers.draws >= 0) and np.all(cond.never_takers.draws >= 0)
