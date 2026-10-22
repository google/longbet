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

"""Frequentist calibration of the reported intervals.

The Geweke test in ``test_geweke.py`` checks the conditionals exactly; this one
checks what a user actually reads off a fit. It replicates a fixed design many
times and asks how often the nominal 95 percent ATT interval covers the truth.

Coverage is the quantity the reference implementation's stability diagnostic was
calibrated against, and it is the only evidence that would justify reinstating a
reliability verdict (see ``longbet._diagnostics``).

The sampler settings below are deliberately shorter than the package defaults
for continuous integration. This is a regression check, not a guaranteed floor
on coverage or a substitute for convergence and calibration at deployed
settings. Earlier numerical coverage claims predate the dynamic precision-cache
repair and must not be quoted as validation of the current engine. Exposure
comparisons within one panel are dependent, so uncertainty in a formal coverage
study must be assessed across independently simulated panels.
"""

import warnings

import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, get_att

N_REPLICATES = 20
NOMINAL = 0.95


def _replicate(seed, n=120, T=8, t0=4, sigma=0.3):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 3)).astype(np.float32)
    gamma = rng.normal(0.0, 0.6, size=n).astype(np.float32)

    adopt = rng.integers(t0, T + 1, size=n)
    adopt = np.where(rng.uniform(size=n) < 0.4, 10**6, adopt)
    z = np.zeros((n, T), dtype=np.float32)
    for i in range(n):
        if adopt[i] <= T:
            z[i, adopt[i] - 1:] = 1.0

    t = np.arange(1, T + 1, dtype=np.float32)
    s = np.where(z == 1, np.maximum(t[None, :] - (adopt[:, None] - 1), 0), 0).astype(int)
    tau = 1.2 * np.sqrt(s) * (1.0 + 0.3 * x[:, 1])[:, None]
    y = (
        (0.5 * x[:, 0])[:, None] + gamma[:, None] + 0.1 * t[None, :]
        + tau + rng.normal(0.0, sigma, size=(n, T))
    ).astype(np.float32)

    s_max = int(s.max())
    true_att = np.array([
        tau[(s == k) & (z == 1)].mean() if np.any((s == k) & (z == 1)) else np.nan
        for k in range(1, s_max + 1)
    ])
    return x, y, z, t, true_att


@pytest.mark.slow
def test_att_interval_coverage_is_at_least_nominal():
    hits = []
    for rep in range(N_REPLICATES):
        x, y, z, t, true_att = _replicate(1000 + rep)
        model = LongBet(LongBetConfig(
            num_sweeps=60, num_burnin=40, num_chains=2,
            num_trees_pr=15, num_trees_trt=15, random_seed=rep,
        ))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            model.fit(y=y, x=x, z=z, t=t)
            att = get_att(model.predict(x=x, z=z, t=t, summary_only=True))

        ok = np.isfinite(true_att) & np.isfinite(att["att"])
        hits.append(
            (true_att[ok] >= att["intervals"][0][ok])
            & (true_att[ok] <= att["intervals"][1][ok])
        )

    covered = np.concatenate(hits)
    rate = float(covered.mean())
    # Binomial standard error at the nominal rate for this many comparisons.
    se = np.sqrt(NOMINAL * (1 - NOMINAL) / covered.size)
    assert rate >= NOMINAL - 3 * se, (
        f"ATT interval coverage is {rate:.3f} over {covered.size} comparisons, "
        f"below nominal {NOMINAL} by more than 3 standard errors ({se:.3f}). "
        "Intervals are too narrow."
    )
    assert rate <= 1.0


@pytest.mark.slow
def test_posterior_contracts_as_the_panel_grows():
    """More units must sharpen the ATT, or the posterior is not using the data."""
    widths = []
    for n in (60, 480):
        x, y, z, t, _ = _replicate(77, n=n)
        model = LongBet(LongBetConfig(num_sweeps=60, num_burnin=40, num_trees_pr=15,
                                      num_trees_trt=15, random_seed=0))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            model.fit(y=y, x=x, z=z, t=t)
            att = get_att(model.predict(x=x, z=z, t=t, summary_only=True))
        widths.append(float(np.nanmean(att["intervals"][1] - att["intervals"][0])))

    assert widths[1] < widths[0], (
        f"interval width did not shrink with n: {widths[0]:.3f} -> {widths[1]:.3f}"
    )
