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

"""Milestone 0: with T = 1, beta fixed at 1 and no unit intercepts, LongBet is BCF.

This validates the composition machinery -- residual bookkeeping, unit
conversion, the weighted treatment forest, the disabled double sigma^2 draw --
independently of any panel-specific code.

The reference is the ``bartz.bcf`` module from the open pull request
bartz-org/bartz#189, which is not on PyPI. Point ``LONGBET_BARTZ_BCF_PATH`` at a
checkout of that branch's ``src/bartz`` to enable this test; without it the test
skips rather than failing, and the reduction is still checked against the
analytic truth below, which needs no external code.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from longbet import LongBet, LongBetConfig

# The faithful reduction. With T = 1 there is no exposure dynamic to model, so
# the treatment forest is blinded to the exposure index; beta is pinned at 1;
# and the treated/control contrast is carried by the adaptive coding weights,
# exactly as in XBCF. The model is then y = mu(x) + b_Z nu(x), whose estimand
# is (b1 - b0) nu(x) -- BCF's.
BCF_REDUCTION = dict(
    num_sweeps=100, num_burnin=60, num_trees_pr=20, num_trees_trt=20,
    num_chains=1,             # a cross-section needs no convergence story
    random_intercept=False,   # no panel structure
    sample_beta=False,        # beta stays at its initial value of 1
    standardize=False,
    random_seed=42,
)


def _cross_section(n=400, p=4, seed=42):
    """A BCF-shaped problem: one period, heterogeneous effect, confounding."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, p)).astype(np.float32)
    pi = 1.0 / (1.0 + np.exp(-x[:, 0]))
    z = rng.binomial(1, pi).astype(np.float32)
    mu = 1.0 * x[:, 0] + 0.8 * x[:, 1]
    tau = 2.0 + 1.0 * x[:, 2]
    y = (mu + tau * z + rng.normal(0, 0.3, size=n)).astype(np.float32)
    return x, y, z, pi.astype(np.float32), tau


def _fit_longbet(x, y, z, pi):
    model = LongBet(LongBetConfig(**BCF_REDUCTION))
    model.fit(
        y=y[:, None], x=x, z=z[:, None],
        t=np.array([1.0], dtype=np.float32), ps=pi,
    )
    # The effect is reported for the schedule passed to predict (zero before
    # treatment), so ask for everyone treated to read tau(x) for every unit.
    pred = model.predict(x=x, z=np.ones((len(x), 1), dtype=np.float32),
                         t=np.array([1.0], dtype=np.float32), ps=pi)
    return model, pred


def test_reduces_to_a_cross_sectional_causal_forest():
    """With T = 1 the model must recover a BCF-style heterogeneous effect."""
    x, y, z, pi, tau_true = _cross_section()
    _, pred = _fit_longbet(x, y, z, pi)
    tau_hat = pred.tau_summary.mean[:, 0]

    corr = np.corrcoef(tau_hat, tau_true)[0, 1]
    assert corr > 0.85, f"correlation with the true effect is only {corr:.3f}"
    bias = abs(tau_hat.mean() - tau_true.mean())
    assert bias < 0.25, f"average effect is off by {bias:.3f}"


def test_beta_is_held_fixed_in_the_reduction():
    x, y, z, pi, _ = _cross_section(n=200)
    model, _ = _fit_longbet(x, y, z, pi)
    assert np.allclose(np.asarray(model.trace.beta), 1.0), (
        "sample_beta=False must leave beta at its initial value"
    )
    assert np.all(np.asarray(model.trace.gamma) == 0.0), (
        "random_intercept=False must leave gamma at zero"
    )
    # Treatment-only coding: y = mu(x) + tau(x) z, which is BCF's contrast.
    assert np.all(np.asarray(model.trace.b0) == 0) and np.all(np.asarray(model.trace.b1) == 1)
    # The treatment forest cannot see the exposure index in this reduction.
    assert not {b.name: b for b in model.design_.blocks}["s"].nu_visible


@pytest.mark.slow
def test_matches_bartz_bcf_from_pr_189():
    """Agreement with the PR #189 BCF sampler, in posterior summaries.

    Not draw for draw: these are two different samplers of the same model, so
    they agree in distribution and nothing stronger is meaningful.
    """
    bcf_path = os.environ.get("LONGBET_BARTZ_BCF_PATH")
    if not bcf_path or not Path(bcf_path).is_dir():
        pytest.skip(
            "set LONGBET_BARTZ_BCF_PATH to a checkout of bartz-org/bartz#189's "
            "src/bartz to run the BCF comparison"
        )

    import bartz

    if str(bcf_path) not in bartz.__path__:
        bartz.__path__.append(str(bcf_path))
    try:
        from bartz.bcf import bcf
    except ImportError:  # pragma: no cover - depends on an external checkout
        pytest.skip(f"no bartz.bcf module under {bcf_path}")

    x, y, z, pi, tau_true = _cross_section()

    reference = bcf(
        x_train=x, y_train=y, z_train=z, pihat_train=pi,
        num_trees_mu=BCF_REDUCTION["num_trees_pr"],
        num_trees_tau=BCF_REDUCTION["num_trees_trt"],
        ndpost=BCF_REDUCTION["num_sweeps"], nskip=BCF_REDUCTION["num_burnin"],
        standardize=False, seed=42,
    )
    tau_bcf = np.mean(reference.predict(x, pihat_test=pi)["tau"], axis=0)

    _, pred = _fit_longbet(x, y, z, pi)
    tau_lb = pred.tau_summary.mean[:, 0]

    corr = np.corrcoef(tau_bcf, tau_lb)[0, 1]
    assert corr > 0.90, f"correlation with bartz.bcf is only {corr:.3f}"
    # Both should sit within Monte Carlo error of the same average effect.
    assert abs(tau_bcf.mean() - tau_lb.mean()) < 0.2
    assert abs(tau_bcf.mean() - tau_true.mean()) < 0.3
