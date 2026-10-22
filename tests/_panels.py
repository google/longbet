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

"""Synthetic panels shared by the test suite.

A plain module rather than ``conftest.py``: importing helpers from a conftest
only works when the repository root happens to be on ``sys.path``, which is true
under ``python -m pytest`` and false under a bare ``pytest``. That difference
made the suite pass locally and fail in CI. pytest puts this file's directory on
``sys.path`` during collection, so ``from _panels import ...`` works either way.
"""

from __future__ import annotations

import numpy as np


def make_staggered_panel(seed=42, N=60, T=8, P=4, effect=1.5, noise=0.25,
                         gamma_sd=0.5, missing_rate=0.0, binary=False):
    """A staggered-adoption panel with a known exposure-dependent effect.

    The true effect grows with exposure as ``effect * sqrt(S)`` and is modulated
    by a covariate, so the treatment forest has something to find and the ATT is
    not flat.
    """
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, P)).astype(np.float32)
    t = np.arange(1, T + 1, dtype=np.float32)

    adopt = rng.choice(np.arange(3, T + 1).tolist() + [10**6], size=N)
    z = np.zeros((N, T), dtype=np.float32)
    for i in range(N):
        if adopt[i] <= T:
            z[i, adopt[i] - 1:] = 1.0
    s = np.where(z == 1, np.maximum(t[None, :] - (adopt[:, None] - 1), 0), 0).astype(int)

    gamma = rng.normal(0.0, gamma_sd, size=N).astype(np.float32)
    prognostic = 0.5 * x[:, 0] - 0.3 * x[:, 1]
    modifier = 1.0 + 0.4 * x[:, 2]
    tau = effect * np.sqrt(s) * modifier[:, None]

    y = (
        prognostic[:, None]
        + 0.1 * t[None, :]
        + gamma[:, None]
        + tau
        + rng.normal(0.0, noise, size=(N, T))
    ).astype(np.float32)

    true_att = np.array(
        [tau[(s == k) & (z == 1)].mean() if np.any((s == k) & (z == 1)) else np.nan
         for k in range(1, s.max() + 1)]
    )

    if binary:
        y = (y > np.median(y)).astype(np.float32)
    if missing_rate:
        y = y.copy()
        y[rng.uniform(size=y.shape) < missing_rate] = np.nan

    return {
        "x": x, "y": y, "z": z, "t": t, "s": s,
        "gamma_true": gamma, "tau_true": tau, "true_att": true_att,
        "true_sigma": noise,
    }
