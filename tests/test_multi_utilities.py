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

"""Draw-scale and alignment regressions, independent of expensive model fitting."""

import numpy as np
import pytest
from scipy.special import ndtr

from longbet import (
    LongBetMultiPrediction, LongBetPrediction, effect_draws,
    effect_draws_from_arrays, joint_prob, reduce_joint_masks,
)


def _prediction(tau, outcome="continuous", mu0=None):
    tau = np.asarray(tau, dtype=float).reshape(1, 1, -1)
    return LongBetPrediction(
        tauhats=tau, muhats0=np.zeros_like(tau) if mu0 is None else mu0,
        yhats=tau, att_full=tau[:, 0], beta_values=np.zeros((tau.shape[-1], 1)),
        z=np.ones((1, 1)), s=np.ones((1, 1)), tau_summary=None,
        mu0_summary=None, y_summary=None, outcome=outcome, num_chains=2,
    )


def _multi():
    return LongBetMultiPrediction(
        preds=[_prediction([1, 1, -1, -1]), _prediction([-1, 1, -1, 1]),
               _prediction([-1, -1, -1, 1], "binary")],
        outcome_names=["gmv", "hours", "complaint"],
        outcome=["continuous", "continuous", "binary"],
        num_chains=2, num_sweeps=2, summary_only=False, provenance="synthetic",
    )


def test_three_outcome_event_uses_names_and_paired_draws():
    pred = _multi()
    conditions = {"complaint": lambda a: a < 0, "hours": lambda a: a < 0,
                  "gmv": lambda a: a > 0}
    np.testing.assert_array_equal(joint_prob(pred, conditions), [[0.25]])
    assert 0.5 * 0.5 * 0.75 == 3 / 16  # product differs from intersection
    np.testing.assert_array_equal(joint_prob(pred, conditions, np.ones((1, 1), bool)), [0.25])


@pytest.mark.parametrize("semantics", [None, "recursive_sur_v1"])
def test_joint_rejects_legacy_prediction_semantics(semantics):
    pred = _multi()
    pred.sampler_semantics = semantics
    with pytest.raises(ValueError, match="Refit from the original data"):
        joint_prob(pred, [lambda a: a > 0] * 3)


def test_binary_probability_transform_uses_baseline_from_same_draw():
    tau = np.array([0.1, 0.5, -0.3, -0.7]).reshape(1, 1, 4)
    mu0 = np.array([-2, 0, 0.7, 2]).reshape(1, 1, 4)
    p = _prediction(tau, "binary", mu0)
    np.testing.assert_allclose(effect_draws(p), ndtr(mu0 + tau) - ndtr(mu0))
    np.testing.assert_array_equal(effect_draws(_prediction(tau)), tau)
    with pytest.raises(ValueError, match="conflicts"):
        effect_draws(p, "continuous")
    p.muhats0 = mu0[:, :, :1]  # broadcasting would silently pair a wrong baseline
    with pytest.raises(ValueError, match="same.*shape"):
        effect_draws(p)


@pytest.mark.parametrize("bad", [np.zeros((1, 4)), np.zeros((1, 1, 0)), np.full((1, 1, 4), np.nan)])
def test_array_helper_rejects_bad_draws(bad):
    with pytest.raises(ValueError):
        effect_draws_from_arrays(bad, outcome="continuous")


def test_joint_validates_effect_dimensions_before_callback():
    pred = _multi()
    pred["hours"].tauhats = np.ones((1, 1, 2))
    with pytest.raises(ValueError, match="Effect draws.*shape"):
        joint_prob(pred, [lambda a: a > 0, lambda a: np.ones((1, 1, 4), bool), lambda a: a < 0])


def test_joint_rejects_inconsistent_chain_metadata():
    pred = _multi()
    pred["hours"].num_chains = 1
    with pytest.raises(ValueError, match="Draw/chain metadata"):
        joint_prob(pred, [lambda a: a > 0] * 3)


def test_exported_mask_reducer_and_invalid_conditions():
    with pytest.raises(ValueError, match="nonempty"):
        reduce_joint_masks([np.ones((1, 1, 0), bool)])
    with pytest.raises(ValueError, match="boolean"):
        joint_prob(_multi(), [lambda a: a] * 3)
    with pytest.raises(ValueError, match="shape"):
        joint_prob(_multi(), [lambda a: (a > 0).ravel()] * 3)
    pred = _multi()
    pred.summary_only = True
    with pytest.raises(ValueError, match="full posterior"):
        joint_prob(pred, [lambda a: a > 0] * 3)
