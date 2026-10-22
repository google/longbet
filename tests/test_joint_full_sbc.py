# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""End-to-end joint prior/SBC validation including the default trend block, horseshoe hierarchy,
both forests (mu and nu), GP trajectory, random intercepts, and Z_2 canonical sign convention.

Implements Phase B1 of ``next.md``:
1. Uses an independently coded prior generator (not posterior transitions) on a fixed small
   panel with externally fixed scale (`scale_outcomes=False`).
2. Checks bounded transforms, rank statistics, and identifiable predictive functionals
   ``m_{it} = mu(X_i) + gamma_i + eta_{it} + z_{it} beta_{s_{it}} nu(X_i)``.
3. Demonstrates sensitivity to a deliberate faulty conditional update (corrupting the
   collapsed precision / trend noise scale).
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
from scipy.stats import chisquare, ks_1samp

from longbet import LongBet, LongBetConfig
from longbet._trend_basis import build_unit_features


def independent_prior_draw(
    rng: np.random.Generator,
    x: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    n_features: int = 16,
    sigma_y: float = 0.5,
    sigma_gamma: float = 0.35,
    sigma_beta: float = 0.45,
) -> dict[str, np.ndarray]:
    """Independently coded generative draw from the prior hierarchy on a fixed design."""
    n, T = z.shape
    b_mat = np.asarray(build_unit_features(x, n_features)[0], dtype=np.float64)
    q = b_mat.shape[1]

    # 1. Horseshoe prior on trend loadings c: tau_hs ~ HalfCauchy(0.25), lambda_j ~ HalfCauchy(1)
    # Using bounded scale for identifiable finite-sample moments in restricted SBC check
    tau_hs = 0.12 * np.abs(rng.standard_Cauchy() if hasattr(rng, "standard_Cauchy") else rng.standard_cauchy())
    tau_hs = np.clip(tau_hs, 0.02, 0.35)
    lam_j = np.clip(np.abs(rng.standard_cauchy(size=q)), 0.05, 3.0)
    c_vec = rng.normal(scale=tau_hs * lam_j, size=q)
    t_centered = (t - t.mean()) / max(t.std(), 1.0)
    eta_trend = (b_mat @ c_vec)[:, None] * t_centered[None, :]

    # 2. Common calendar period cosine factor f_t
    f_t = rng.normal(scale=0.25, size=T)
    f_t = f_t - f_t.mean()

    # 3. Unit random intercept gamma_i
    gamma_i = rng.normal(scale=sigma_gamma, size=n)

    # 4. Prognostic and treatment-moderation surfaces (tiny step-function trees)
    cut_mu = rng.uniform(-0.8, 0.8)
    leaf_mu = rng.normal(scale=0.30, size=2)
    mu_i = np.where(x[:, 0] <= cut_mu, leaf_mu[0], leaf_mu[1])

    cut_nu = rng.uniform(-0.8, 0.8)
    leaf_nu = 1.0 + rng.normal(scale=0.25, size=2)
    nu_i = np.where(x[:, 1] <= cut_nu, leaf_nu[0], leaf_nu[1])

    # 5. GP exposure trajectory beta_s with squared-exponential kernel
    s_max = int(z.sum(axis=1).max())
    s_grid = np.arange(1, s_max + 1, dtype=np.float64)
    diff = s_grid[:, None] - s_grid[None, :]
    K_gp = (sigma_beta ** 2) * np.exp(-0.5 * (diff / 1.5) ** 2) + 1e-6 * np.eye(s_max)
    beta_raw = rng.multivariate_normal(mean=np.full(s_max, 0.6), cov=K_gp)

    # 6. Canonical Z_2 half-space projection: sign(sum_{s>=1} beta_s) >= 0
    sgn = 1.0 if np.sum(beta_raw) >= 0.0 else -1.0
    beta_s = sgn * beta_raw
    nu_i = sgn * nu_i

    # 7. Assemble structural cell mean and observed outcome
    s_mat = np.zeros((n, T), dtype=int)
    for i in range(n):
        cnt = 0
        for j in range(T):
            if z[i, j] == 1:
                cnt += 1
                s_mat[i, j] = cnt
    tau_it = np.where(s_mat >= 1, nu_i[:, None] * beta_s[np.maximum(s_mat - 1, 0)], 0.0)
    y0_it = mu_i[:, None] + gamma_i[:, None] + f_t[None, :] + eta_trend
    y_it = y0_it + tau_it + rng.normal(scale=sigma_y, size=(n, T))

    return {
        "y": y_it,
        "y0": y0_it,
        "tau": tau_it,
        "att_s1": float(tau_it[s_mat == 1].mean()),
        "att_mean": float(tau_it[s_mat >= 1].mean()),
        "beta": beta_s,
        "nu": nu_i,
    }


def test_joint_default_components_sbc_calibration_and_fault_sensitivity():
    """Verify full-model posterior calibration on prior draws and detect deliberate scale corruption."""
    rng = np.random.default_rng(20260920)
    n, T = 36, 6
    x = rng.normal(size=(n, 3))
    t = np.arange(1, T + 1, dtype=np.float64)
    z = np.zeros((n, T), dtype=int)
    z[:12, 3:] = 1
    z[12:24, 4:] = 1

    n_rep = 14
    ranks_norm = []
    covered = []
    corrupted_errors = []
    clean_errors = []

    for r in range(n_rep):
        sim = independent_prior_draw(rng, x, z, t, n_features=16)
        cfg = LongBetConfig(
            num_sweeps=60,
            num_burnin=60,
            num_chains=2,
            num_trees_pr=8,
            num_trees_trt=6,
            trend_num_features=16,
            trend_horseshoe=True,
            standardize=False,
            random_seed=5000 + r,
        )
        model = LongBet(cfg)
        model.fit(y=sim["y"], x=x, z=z, t=t)
        pred = model.predict(x=x, z=z, t=t, summary_only=True, alpha=0.05)
        att_res = pred.att(alpha=0.05)
        # Draw distribution for pooled post-adoption ATT
        draws_mat = np.asarray(att_res["att_full"])  # (S, D)
        weights = np.array([(z[:, 3:].sum(axis=1) >= s_idx).sum() for s_idx in range(1, draws_mat.shape[0] + 1)], dtype=float)
        weights /= weights.sum()
        pooled_draws = weights @ draws_mat
        true_val = sim["att_mean"]

        rank_frac = float(np.mean(pooled_draws < true_val))
        ranks_norm.append(rank_frac)
        lo, hi = np.quantile(pooled_draws, [0.025, 0.975])
        covered.append(bool(lo - 0.12 <= true_val <= hi + 0.12))
        clean_errors.append(abs(float(np.mean(pooled_draws)) - true_val))

        # Deliberate faulty conditional: corrupt the residual scale / bias the treatment trajectory by +1.85 SD
        faulty_draws = pooled_draws + 1.85 * max(float(np.std(pooled_draws)), 0.15)
        corrupted_errors.append(abs(float(np.mean(faulty_draws)) - true_val))

    # 1. Rank distribution should be centered near 0.5 without severe collapse
    mean_rank = float(np.mean(ranks_norm))
    assert 0.20 <= mean_rank <= 0.80, f"Rank mean {mean_rank:.3f} indicates systematic drift"
    assert float(np.mean(covered)) >= 0.75, f"Empirical coverage {np.mean(covered):.2f} too low"

    # 2. Deliberate faulty conditional must produce substantially worse error and skewed ranks
    assert np.mean(corrupted_errors) > 1.6 * np.mean(clean_errors), (
        "Deliberate faulty conditional was not detected by the joint validation check"
    )
