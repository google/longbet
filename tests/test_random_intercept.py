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

"""Unit random intercepts, and the ordering that makes them identifiable."""

import numpy as np

from longbet import LongBet, LongBetConfig
from _panels import make_staggered_panel


def test_recovers_unobserved_unit_levels():
    data = make_staggered_panel(seed=20260905, N=150, T=8, gamma_sd=1.0, noise=0.25)
    model = LongBet(LongBetConfig(num_sweeps=100, num_burnin=100, num_trees_pr=10,
                                  num_trees_trt=10, random_intercept=True,
                                  num_chains=1, random_seed=42))
    model.fit(y=data["y"], x=data["x"], z=data["z"], t=data["t"])

    gamma_hat = np.asarray(model.trace.gamma).mean(axis=0)
    corr = np.corrcoef(gamma_hat, data["gamma_true"])[0, 1]
    # gamma_i competes with mu(X_i), which is also unit-constant, so the two are
    # separated only by their priors and this correlation has a ceiling below 1.
    assert corr > 0.85, f"gamma correlation with the truth is only {corr:.3f}"

    resid_sd = float(np.mean(np.asarray(model.trace.sigma2) ** 0.5) * model.sdy)
    assert abs(resid_sd - data["true_sigma"]) < 0.12, f"residual sd {resid_sd:.3f}"


def test_gamma_tracks_the_truth_not_adoption_timing():
    """Regression on sweep order.

    gamma_i is drawn conditional on the current treatment fit, so a treated
    unit's post-adoption periods must not drag its baseline upward. If the
    intercept block were moved ahead of the treatment forest, gamma would soak
    up the effect and correlate with adoption timing instead of with the truth.
    """
    data = make_staggered_panel(seed=5, N=200, T=10, effect=3.0, gamma_sd=1.0, noise=0.3)
    model = LongBet(LongBetConfig(num_sweeps=40, num_burnin=25, num_trees_pr=10,
                                  num_trees_trt=10, random_intercept=True,
                                  num_chains=1, random_seed=7))
    model.fit(y=data["y"], x=data["x"], z=data["z"], t=data["t"])
    gamma_hat = np.asarray(model.trace.gamma).mean(axis=0)

    exposure_total = data["s"].max(axis=1).astype(float)
    corr_truth = abs(np.corrcoef(gamma_hat, data["gamma_true"])[0, 1])
    corr_timing = abs(np.corrcoef(gamma_hat, exposure_total)[0, 1])

    assert corr_truth > 0.8, f"gamma barely tracks the truth ({corr_truth:.3f})"
    assert corr_truth > 3 * corr_timing, (
        f"gamma correlates with adoption timing ({corr_timing:.3f}) nearly as much as "
        f"with the truth ({corr_truth:.3f}); the intercept block is absorbing the effect"
    )


def test_random_intercept_off_leaves_gamma_at_zero():
    data = make_staggered_panel(seed=1, N=40, T=6)
    model = LongBet(LongBetConfig(num_sweeps=10, num_burnin=5, num_trees_pr=5,
                                  num_trees_trt=5, random_intercept=False,
                                  num_chains=1, random_seed=42))
    model.fit(y=data["y"], x=data["x"], z=data["z"], t=data["t"])
    assert np.all(np.asarray(model.trace.gamma) == 0.0)
    assert np.all(np.asarray(model.trace.sigma_gamma2) == model.trace.sigma_gamma2[0])
