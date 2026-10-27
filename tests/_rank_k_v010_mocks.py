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

"""Synthetic pairwise Rank-K fits shared by the v0.10 methodology tests.

The helpers build :class:`longbet.RankKLongBetResult` objects whose Mode-1 posterior
is ``N(y_hat_p, s_p^2)`` for a contrast ``p = (j, l)``, where ``y_hat_p`` follows the
shared-arm sampling model ``y = B mu + E eps`` (one realized arm mean per arm, shared
by every contrast that involves the arm). They are also used to generate the stored
v0.9.0 parity fixture ``tests/data/rank_k_v090_parity.npz``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import itertools
import math

import numpy as np

import longbet._rank_k as rk


def make_mock_fit(
    *,
    att_mean: float,
    att_sd: float,
    n_draws: int = 600,
    n_units: int = 40,
    t_pre: int = 4,
    t_post: int = 12,
    rank_k: int = 2,
    mode2_mean: float = 0.0,
    mode2_sd: float = 0.02,
    seed: int = 0,
    exact_moments: bool = False,
    att_draws: np.ndarray | None = None,
) -> rk.RankKLongBetResult:
    """Builds a synthetic :class:`RankKLongBetResult` with a Gaussian Mode-1 posterior.

    Args:
      att_mean: Posterior mean of the Mode-1 ATT draws.
      att_sd: Posterior standard deviation of the Mode-1 ATT draws.
      n_draws: Number of posterior draws ``S``.
      n_units: Number of panel units (half of them treated).
      t_pre: Number of pre-adoption periods.
      t_post: Number of post-adoption periods.
      rank_k: Number of exposure modes (1 or 2).
      mode2_mean: Mean of the Mode-2 (shape) amplitude draws.
      mode2_sd: Standard deviation of the Mode-2 amplitude draws.
      seed: Random seed.
      exact_moments: Standardize the raw normal draws so the sample mean and SD
        are exactly ``att_mean`` and ``att_sd``.
      att_draws: Optional explicit ``(S,)`` Mode-1 draws overriding the Gaussian.

    Returns:
      A :class:`RankKLongBetResult` with ``catt_draws``, ``mode_tau_draws``,
      ``bases``, ``cohorts`` and ``unit_weights`` filled in.
    """
    rng = np.random.default_rng(seed)
    n_periods = t_pre + t_post
    cohorts = np.full(n_units, n_periods, dtype=np.int32)
    cohorts[: n_units // 2] = t_pre
    ever = cohorts < n_periods
    unit_weights = np.where(ever, 1.0 / np.sum(ever), 0.0)

    if att_draws is None:
        raw_z = rng.normal(loc=0.0, scale=1.0, size=n_draws)
        if exact_moments:
            raw_z = (raw_z - float(np.mean(raw_z))) / float(np.std(raw_z, ddof=1))
        att = float(att_mean) + float(att_sd) * raw_z
    else:
        att = np.asarray(att_draws, dtype=np.float64)
        n_draws = int(att.size)

    mode1_draws = np.zeros((n_units, n_draws), dtype=np.float64)
    mode1_draws[ever] = att[None, :] + rng.normal(0.0, 0.01, size=(int(np.sum(ever)), n_draws))
    mode1_draws[ever] += (att - mode1_draws[ever].mean(axis=0))[None, :]

    mode_tau_draws = [mode1_draws]
    basis = rk.orthonormal_exposure_basis(t_post, rank_k)
    if rank_k >= 2:
        m2_trace = rng.normal(loc=mode2_mean, scale=mode2_sd, size=n_draws)
        mode2_draws = np.zeros((n_units, n_draws), dtype=np.float64)
        mode2_draws[ever] = m2_trace[None, :]
        mode_tau_draws.append(mode2_draws)

    catt_draws = np.zeros((n_draws, n_units, n_periods), dtype=np.float64)
    for k_idx, m_draws in enumerate(mode_tau_draws):
        psi_k = basis[:, k_idx]
        catt_draws[:, ever, t_pre:] += m_draws[ever].T[:, :, None] * psi_k[None, None, :]

    return rk.RankKLongBetResult(
        catt_draws=catt_draws,
        catt_hat=catt_draws.mean(axis=0),
        catt_ci_low=np.quantile(catt_draws, 0.025, axis=0),
        catt_ci_high=np.quantile(catt_draws, 0.975, axis=0),
        att_draws=att,
        att_hat=float(np.mean(att)),
        att_se=float(np.std(att, ddof=1)),
        att_ci_low=float(np.quantile(att, 0.025)),
        att_ci_high=float(np.quantile(att, 0.975)),
        mode_tau_draws=mode_tau_draws,
        theta_f_modes=[1.0] * rank_k,
        sigma_modes=[float(att_sd)] * rank_k,
        bases={t_pre: basis},
        runtime_sec=0.01,
        rhat_att=1.001,
        ess_att=float(n_draws * 0.8),
        cohorts=cohorts,
        unit_weights=unit_weights,
    )


def all_pairs(n_arms: int) -> list[tuple[int, int]]:
    """Returns every ordered pair ``(j, l)`` with ``j > l`` on ``n_arms`` arms."""
    return [(j, l) for l, j in itertools.combinations(range(n_arms), 2)]


def simulate_shared_arm_fits(
    rng: np.random.Generator,
    mu_full: np.ndarray,
    arm_vars: np.ndarray,
    pairs: Sequence[tuple[int, int]],
    arms: Sequence[str],
    *,
    n_draws: int = 600,
    n_units: int = 8,
    t_post: int = 6,
    pair_bias: Mapping[tuple[int, int], float] | None = None,
    z_bar: np.ndarray | None = None,
) -> tuple[dict[tuple[str, str], rk.RankKLongBetResult], np.ndarray]:
    """Simulates pairwise fits from the shared-arm model ``y = B mu + E eps``.

    Args:
      rng: Random generator.
      mu_full: ``(A,)`` true arm effects with ``mu_full[0] == 0`` for the control.
      arm_vars: ``(A,)`` true arm sampling variances ``v_a``.
      pairs: Fitted contrasts ``(j, l)`` as arm indices.
      arms: Arm labels.
      n_draws: Number of posterior draws per fit.
      n_units: Number of units in each mock fit (kept small for speed).
      t_post: Number of post-adoption periods in each mock fit.
      pair_bias: Optional additive bias on the realized summary of some pairs.
      z_bar: Optional pre-drawn realized arm means (overrides the simulation).

    Returns:
      ``(fits, y_hat)`` where ``fits`` maps ``(arms[j], arms[l])`` to a mock fit and
      ``y_hat`` holds the realized pairwise summaries in ``pairs`` order.
    """
    mu_full = np.asarray(mu_full, dtype=np.float64)
    arm_vars = np.asarray(arm_vars, dtype=np.float64)
    if z_bar is None:
        z_bar = mu_full + np.sqrt(arm_vars) * rng.standard_normal(mu_full.size)
    fits: dict[tuple[str, str], rk.RankKLongBetResult] = {}
    y_hat = np.zeros(len(pairs), dtype=np.float64)
    for p_idx, (j, l) in enumerate(pairs):
        bias = 0.0 if pair_bias is None else float(pair_bias.get((j, l), 0.0))
        sd_p = math.sqrt(float(arm_vars[j] + arm_vars[l]))
        y_p = float(z_bar[j] - z_bar[l]) + bias
        y_hat[p_idx] = y_p
        fits[(arms[j], arms[l])] = make_mock_fit(
            att_mean=y_p,
            att_sd=sd_p,
            n_draws=n_draws,
            n_units=n_units,
            t_post=t_post,
            seed=int(rng.integers(0, 2**31 - 1)),
            exact_moments=True,
        )
    return fits, y_hat
