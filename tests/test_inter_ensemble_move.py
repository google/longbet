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

"""Unit tests for the inter-ensemble residual-transfer move, support diagnostics, and API cleanup."""

import warnings

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, get_att
from longbet._diagnostics import compute_ess
from longbet._inter_ensemble_move import inter_ensemble_transfer_step
from longbet._state import init_longbet
from longbet._step import longbet_single_step


def test_inter_ensemble_step_preserves_residual_invariant():
    """Verify inter_ensemble_transfer_step preserves R == z - alpha*mu_fit - w*nu_fit - gamma."""
    N, T, P = 12, 6, 2
    M = N * T
    rng = np.random.default_rng(42)
    z_mat = np.zeros((N, T), dtype=np.float32)
    z_mat[:6, 3:] = 1.0
    exposure = np.where(z_mat == 1, np.tile(np.arange(1, T + 1), (N, 1)) - 3, 0).astype(np.int32)
    y = rng.normal(size=M).astype(np.float32)

    config = LongBetConfig(
        num_trees_pr=4,
        num_trees_trt=4,
        use_inter_ensemble_move=True,
    )
    state = init_longbet(
        X_unified=jnp.asarray(rng.integers(0, 4, (P, M), dtype=np.uint8)),
        y=jnp.asarray(y),
        unit_idx=jnp.repeat(jnp.arange(N, dtype=jnp.int32), T),
        time_idx=jnp.tile(jnp.arange(T, dtype=jnp.int32), N),
        exposure_idx=jnp.asarray(exposure.ravel()),
        z_vec=jnp.asarray(z_mat.ravel()),
        obs_mask=jnp.ones(M, dtype=jnp.bool_),
        max_split_mu=jnp.full(P, 3, dtype=jnp.uint8),
        max_split_nu=jnp.full(P, 3, dtype=jnp.uint8),
        config=config,
    )

    # Run a few sweeps so trees grow splits and leaves become non-trivial
    key = jax.random.key(101)
    for _ in range(5):
        key, sub = jax.random.split(key)
        state = longbet_single_step(sub, state)

    b_z = jnp.where(state.z_vec == 1.0, state.b1, state.b0)
    w = b_z * state.beta[state.exposure_idx]

    response = state.y if state.z is None else state.z

    trend_term = (
        state.trend_design @ state.trend_coef if state.use_trend_block else 0.0
    )

    # Check residual invariant before step
    expected_resid_before = (
        response
        - state.alpha * state.mu_fit
        - w * state.nu_fit
        - state.gamma[state.unit_idx]
        - trend_term
    )
    np.testing.assert_allclose(state.resid, expected_resid_before, atol=2e-5)

    # Execute inter_ensemble_transfer_step directly
    key, sub = jax.random.split(key)
    f_mu_new, mu_new, f_nu_new, nu_new, R_new = inter_ensemble_transfer_step(
        sub,
        state.forest,
        state.mu_fit,
        state.forest_nu,
        state.nu_fit,
        state.resid,
        state.alpha,
        w,
        state.obs_mask,
        state.z_vec,
        state.time_idx,
        state.T_periods,
        state.sigma2,
        state.forest.leaf_prior_cov_inv,
        state.leaf_prior_cov_inv_nu,
    )

    expected_resid_after = (
        response
        - state.alpha * mu_new
        - w * nu_new
        - state.gamma[state.unit_idx]
        - trend_term
    )
    np.testing.assert_allclose(R_new, expected_resid_after, atol=2e-5)


def test_geweke_invariance_with_inter_ensemble_move():
    """Geweke joint-distribution invariance check with use_inter_ensemble_move=True."""
    N, T, P = 10, 5, 2
    M = N * T
    rng = np.random.default_rng(7)
    z_mat = np.zeros((N, T), dtype=np.float32)
    z_mat[:5, 2:] = 1.0
    exposure = np.where(z_mat == 1, np.tile(np.arange(1, T + 1), (N, 1)) - 2, 0).astype(np.int32)

    sigma_a, sigma_b = 5.0, 3.0
    gamma_a, gamma_b = 4.0, 2.0
    config = LongBetConfig(
        num_trees_pr=2,
        num_trees_trt=2,
        random_intercept=True,
        gamma_prior_a=gamma_a,
        gamma_prior_b=gamma_b,
        sigma_prior_a=sigma_a,
        sigma_prior_b=sigma_b,
        use_inter_ensemble_move=True,
    )
    state = init_longbet(
        X_unified=jnp.asarray(rng.integers(0, 3, (P, M), dtype=np.uint8)),
        y=jnp.zeros(M, dtype=jnp.float32),
        unit_idx=jnp.repeat(jnp.arange(N, dtype=jnp.int32), T),
        time_idx=jnp.tile(jnp.arange(T, dtype=jnp.int32), N),
        exposure_idx=jnp.asarray(exposure.ravel()),
        z_vec=jnp.asarray(z_mat.ravel()),
        obs_mask=jnp.ones(M, dtype=jnp.bool_),
        max_split_mu=jnp.full(P, 2, dtype=jnp.uint8),
        max_split_nu=jnp.full(P, 2, dtype=jnp.uint8),
        config=config,
    )

    body = longbet_single_step.__wrapped__
    sweep = jax.jit(lambda k, s: body(k, s))

    key = jax.random.key(20260917)
    sigma2_draws = []
    sigma_gamma2_draws = []
    for it in range(1200):
        key, k_data, k_step = jax.random.split(key, 3)
        b_z = jnp.where(state.z_vec == 1.0, state.b1, state.b0)
        mean = (
            state.alpha * state.mu_fit
            + b_z * state.beta[state.exposure_idx] * state.nu_fit
            + state.gamma[state.unit_idx]
        )
        eps = jax.random.normal(k_data, mean.shape, jnp.float32) * jnp.sqrt(state.sigma2)
        state = eqx.tree_at(lambda s: (s.y, s.resid), state, (mean + eps, eps))
        state = sweep(k_step, state)
        if it >= 200:
            sigma2_draws.append(float(state.sigma2))
            sigma_gamma2_draws.append(float(state.sigma_gamma2))

    s2 = np.array(sigma2_draws)
    sg2 = np.array(sigma_gamma2_draws)

    prior_mean_s2 = sigma_b / (sigma_a - 1.0)
    ess_s2 = max(float(compute_ess(s2[None, :])), 2.0)
    mcse_s2 = float(np.std(s2, ddof=1)) / np.sqrt(ess_s2)
    z_s2 = abs(float(np.mean(s2)) - prior_mean_s2) / max(mcse_s2, 1e-12)
    assert z_s2 < 4.0, f"Geweke sigma^2 mean z-score {z_s2:.2f} exceeds 4.0"

    prior_mean_sg2 = gamma_b / (gamma_a - 1.0)
    ess_sg2 = max(float(compute_ess(sg2[None, :])), 2.0)
    mcse_sg2 = float(np.std(sg2, ddof=1)) / np.sqrt(ess_sg2)
    z_sg2 = abs(float(np.mean(sg2)) - prior_mean_sg2) / max(mcse_sg2, 1e-12)
    assert z_sg2 < 4.0, f"Geweke sigma_gamma^2 mean z-score {z_sg2:.2f} exceeds 4.0"


def test_zero_control_period_warning_and_att_diagnostics():
    """Verify UserWarning when post-treatment periods have 0 untreated controls and check get_att() diagnostics."""
    N, T = 10, 5
    rng = np.random.default_rng(123)
    x = rng.normal(size=(N, 2)).astype(np.float32)
    # Cohort 1 adopts at t=2 (index 1), Cohort 2 adopts at t=4 (index 3).
    # In periods t=4, 5 (indices 3, 4), ALL units are treated -> N_{0,t} = 0!
    z = np.zeros((N, T), dtype=np.float32)
    z[:5, 1:] = 1.0
    z[5:, 3:] = 1.0
    y = rng.normal(size=(N, T)).astype(np.float32)

    model = LongBet(num_burnin=10, num_sweeps=10, num_chains=2)
    with pytest.warns(UserWarning, match="contain zero untreated control observations"):
        model.fit(y=y, x=x, z=z)

    pred = model.predict(x=x, z=z)
    att_res = get_att(pred)

    assert "mean_concurrent_controls" in att_res
    assert "pct_zero_control_cells" in att_res

    # Exposure s=4 occurs only at period index 4 (where z[:5, 4] == 1 and s[:5, 4] == 4).
    # At period index 4, N_{0,4} == 0. So for exposure s=4 (index 3):
    # mean_concurrent_controls[3] == 0.0 and pct_zero_control_cells[3] == 1.0.
    assert att_res["mean_concurrent_controls"][3] == pytest.approx(0.0)
    assert att_res["pct_zero_control_cells"][3] == pytest.approx(1.0)
    # Exposure s=1 occurs at period index 1 for cohort 1 (where N_{0,1} = 5)
    # and at period index 3 for cohort 2 (where N_{0,3} = 0).
    # Mean concurrent controls for s=1 is (5*5 + 5*0)/10 = 2.5, pct_zero is 0.5.
    assert att_res["mean_concurrent_controls"][0] == pytest.approx(2.5)
    assert att_res["pct_zero_control_cells"][0] == pytest.approx(0.5)


def test_split_time_ps_deprecation_and_beta_draws_by_exposure():
    """Verify split_calendar_mu, deprecated split_time_ps warning, and beta_draws_by_exposure property."""
    with pytest.warns(DeprecationWarning, match="split_time_ps"):
        cfg = LongBetConfig(split_time_ps=False)
    assert cfg.split_calendar_mu is False

    with pytest.warns(DeprecationWarning, match="split_time_ps"):
        val = cfg.split_time_ps
    assert val is False

    N, T = 8, 4
    rng = np.random.default_rng(99)
    x = rng.normal(size=(N, 2)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[:4, 2:] = 1.0
    y = rng.normal(size=(N, T)).astype(np.float32)

    with pytest.warns(DeprecationWarning, match="split_time_ps"):
        model = LongBet(num_burnin=5, num_sweeps=6, num_chains=2, split_time_ps=True)
    model.fit(y=y, x=x, z=z)

    # Total draws = num_chains * num_sweeps = 12, S_max = 2 -> S_max + 1 = 3
    assert model.beta_draws_.shape == (12, 3)
    assert model.beta_draws_by_exposure.shape == (3, 12)
    assert model.beta_draws_exposure_major.shape == (3, 12)
    np.testing.assert_allclose(model.beta_draws_by_exposure, model.beta_draws_.T)

    pred = model.predict(x=x, z=z)
    assert pred.beta_draws_by_exposure.shape == (3, 12)
    np.testing.assert_allclose(pred.beta_draws_by_exposure, pred.beta_values.T)
