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

"""Tests for effect_draws, joint_prob, and outcome_correlation algebra."""

import numpy as np
import pytest

from longbet import LongBetConfig, LongBetMulti, effect_draws, joint_prob, outcome_correlation
from longbet._multi_model import LongBetMultiPrediction, reduce_joint_masks


def test_joint_prob_exact_events_and_marginal_bound():
    """Verify hand-built boolean events [T,T,F,F] and [T,F,T,F] give joint 1/4,

    and joint probability never exceeds marginal probabilities.
    """
    # Create synthetic boolean event masks
    # Outcome 1: [True, True, False, False]
    # Outcome 2: [True, False, True, False]
    # Joint:     [True, False, False, False] -> mean = 1/4 = 0.25
    N, T, D = 2, 2, 4
    mask1 = np.array([[[True, True, False, False]] * T] * N, dtype=bool)
    mask2 = np.array([[[True, False, True, False]] * T] * N, dtype=bool)

    p_joint = reduce_joint_masks([mask1, mask2])
    assert p_joint.shape == (N, T)
    np.testing.assert_allclose(p_joint, 0.25)

    p_marg1 = np.mean(mask1, axis=-1)
    p_marg2 = np.mean(mask2, axis=-1)
    assert np.all(p_joint <= p_marg1 + 1e-12)
    assert np.all(p_joint <= p_marg2 + 1e-12)


def test_joint_prob_cells_selection():
    """Verify cells argument: None -> (N, T), boolean -> length-N average over selected periods."""
    N, T, D = 3, 4, 10
    rng = np.random.default_rng(123)
    mask1 = rng.choice([True, False], size=(N, T, D))
    mask2 = rng.choice([True, False], size=(N, T, D))

    p_full = reduce_joint_masks([mask1, mask2])
    assert p_full.shape == (N, T)

    # Unit 0: select period 0 and 1
    # Unit 1: select period 2
    # Unit 2: select none (should give NaN)
    cells = np.zeros((N, T), dtype=bool)
    cells[0, 0] = True
    cells[0, 1] = True
    cells[1, 2] = True

    p_cells = reduce_joint_masks([mask1, mask2], cells=cells)
    assert p_cells.shape == (N,)
    assert p_cells[0] == pytest.approx(float(np.mean(p_full[0, :2])))
    assert p_cells[1] == pytest.approx(float(p_full[1, 2]))
    assert np.isnan(p_cells[2])


def test_outcome_correlation_analytic_two_equation():
    """Verify outcome_correlation reproduces analytic formula:

    With Gamma[1, 0] = g and variances v0, v1:
    Sigma = [[v0, g * v0], [g * v0, g^2 * v0 + v1]]
    rho = g * sqrt(v0) / sqrt(g^2 * v0 + v1)
    """
    g_val = 0.6
    v0_val = 1.4
    v1_val = 0.8

    expected_rho = g_val * np.sqrt(v0_val) / np.sqrt(g_val**2 * v0_val + v1_val)

    # Mock LongBetMulti with single draw
    class MockTrace:
        def __init__(self, s2):
            self.sigma2 = np.array([s2], dtype=np.float64)

    class MockMultiTrace:
        def __init__(self):
            self.gamma_loadings = np.array([[[[0.0, 0.0], [g_val, 0.0]]]], dtype=np.float64)
            self.traces = (MockTrace(v0_val), MockTrace(v1_val))

    cfg = LongBetConfig(sur=True, sur_prior_var=1.0)
    model = LongBetMulti(cfg)
    model.outcome_names = ("y0", "y1")
    model.outcome = ("continuous", "continuous")
    model.order = (0, 1)
    model.inverse_order = (0, 1)
    model.sur_active = True
    model.trace = MockMultiTrace()

    R = outcome_correlation(model)
    assert R.shape == (2, 2)
    assert R[0, 0] == 1.0
    assert R[1, 1] == 1.0
    assert R[0, 1] == pytest.approx(expected_rho, rel=1e-6)
    assert R[1, 0] == pytest.approx(expected_rho, rel=1e-6)


def test_outcome_correlation_negative_and_zero():
    """Verify negative loadings give negative correlation and zero loadings give identity."""
    # Negative loading
    g_val = -0.5
    v0_val = 1.0
    v1_val = 1.0
    expected_rho = g_val * np.sqrt(v0_val) / np.sqrt(g_val**2 * v0_val + v1_val)

    class MockTrace:
        def __init__(self, s2):
            self.sigma2 = np.array([s2], dtype=np.float64)

    class MockMultiTrace:
        def __init__(self, g):
            self.gamma_loadings = np.array([[[[0.0, 0.0], [g, 0.0]]]], dtype=np.float64)
            self.traces = (MockTrace(v0_val), MockTrace(v1_val))

    cfg = LongBetConfig(sur=True, sur_prior_var=1.0)
    model = LongBetMulti(cfg)
    model.outcome_names = ("y0", "y1")
    model.outcome = ("continuous", "continuous")
    model.order = (0, 1)
    model.inverse_order = (0, 1)
    model.sur_active = True
    model.trace = MockMultiTrace(g_val)

    R_neg = outcome_correlation(model)
    assert R_neg[0, 1] == pytest.approx(expected_rho, rel=1e-6)
    assert R_neg[0, 1] < 0.0

    # Inactive SUR -> identity
    model.sur_active = False
    R_ident = outcome_correlation(model)
    np.testing.assert_allclose(R_ident, np.eye(2))
