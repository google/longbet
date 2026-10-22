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

"""End-to-end recovery on the reference simulation design.

Mirrors ``tests/simulation/dgp.R`` from the reference implementation: staggered
adoption driven by a propensity, an exposure-dependent effect that rises then
decays, unit heterogeneity, and a known residual scale.
"""

import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, get_att


def generate_dgp_panel(n=200, t0=5, T=10, seed=42, true_sigma=0.25):
    """Panel with h(s) = s * exp(-0.3 s) modulated by a covariate."""
    rng = np.random.default_rng(seed)

    x = np.column_stack([
        rng.normal(size=n), rng.normal(size=n), rng.normal(size=n),
        rng.binomial(1, 0.5, size=n), rng.choice([1.0, 2.0, 3.0], size=n),
    ]).astype(np.float32)

    gamma = rng.normal(0.0, 0.8, size=n).astype(np.float32)
    pi = np.clip(0.2 * (1 / (1 + np.exp(-0.5 * x[:, 0]))) ** 2 + rng.uniform(0, 0.1, n), 0.05, 0.95)

    z = np.zeros((n, T), dtype=np.float32)
    for j in range(t0 - 1, T):
        z[:, j] = np.where((z[:, j - 1] == 1) | rng.binomial(1, pi).astype(bool), 1.0, 0.0)

    s = np.zeros((n, T), dtype=np.int32)
    for i in range(n):
        treated = np.flatnonzero(z[i] == 1)
        if treated.size:
            s[i, treated[0]:] = np.arange(1, T - treated[0] + 1)

    t = np.arange(1, T + 1, dtype=np.float32)
    nu_x = 1.5 + 0.5 * x[:, 1]
    tau = np.where(s > 0, s * np.exp(-0.3 * s), 0.0) * nu_x[:, None]
    y = (
        (0.8 * x[:, 0] - 0.5 * x[:, 2] + gamma)[:, None]
        + 0.1 * t[None, :]
        + tau
        + rng.normal(0.0, true_sigma, size=(n, T))
    ).astype(np.float32)

    s_max = int(s.max())
    true_att = np.array([
        tau[(s == k) & (z == 1)].mean() if np.any((s == k) & (z == 1)) else np.nan
        for k in range(1, s_max + 1)
    ])

    return {"x": x, "y": y, "z": z, "t": t, "s": s, "gamma_true": gamma,
            "true_sigma": true_sigma, "true_att": true_att, "tau_true": tau, "ps": pi}


@pytest.mark.slow
def test_dgp_recovery():
    data = generate_dgp_panel(seed=101)
    model = LongBet(LongBetConfig(num_sweeps=120, num_burnin=60, num_trees_pr=20,
                                  num_trees_trt=20, num_chains=2, random_intercept=True,
                                  random_seed=42))
    model.fit(y=data["y"], x=data["x"], z=data["z"], t=data["t"], ps=data["ps"])

    gamma_hat = np.asarray(model.trace.gamma).reshape(-1, data["x"].shape[0]).mean(axis=0)
    corr_gamma = np.corrcoef(gamma_hat, data["gamma_true"])[0, 1]
    # gamma_i and mu(X_i) are both unit-constant, so only the priors separate
    # them; correlation with the truth plateaus here well below 1 no matter how
    # long the chain runs. Measured ceiling on this design is about 0.90, and
    # about 0.93 even with no prognostic signal at all.
    assert corr_gamma >= 0.85, f"gamma correlation {corr_gamma:.3f}"

    resid_sd = float(np.mean(np.asarray(model.trace.sigma2) ** 0.5) * model.sdy)
    assert abs(resid_sd - data["true_sigma"]) < 0.05, f"residual sd {resid_sd:.3f}"

    pred = model.predict(x=data["x"], z=data["z"], t=data["t"], ps=data["ps"])

    # Per-cell effects, not just their average.
    treated = data["z"] == 1
    corr_catt = np.corrcoef(
        pred.tau_summary.mean[treated], data["tau_true"][treated]
    )[0, 1]
    assert corr_catt >= 0.90, f"CATT correlation {corr_catt:.3f}"

    att = get_att(pred)
    finite = np.isfinite(data["true_att"]) & np.isfinite(att["att"])
    corr_att = np.corrcoef(att["att"][finite], data["true_att"][finite])[0, 1]
    assert corr_att >= 0.90, f"ATT correlation {corr_att:.3f}"

    # The effect peaks and then decays; the estimate must reproduce the shape.
    assert np.argmax(att["att"][finite]) <= np.argmax(data["true_att"][finite]) + 1

    covered = (
        (data["true_att"][finite] >= att["intervals"][0][finite])
        & (data["true_att"][finite] <= att["intervals"][1][finite])
    )
    assert covered.mean() >= 0.75, f"ATT interval coverage {covered.mean():.2f} on one draw"
