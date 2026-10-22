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

"""Tests for SUR loading regression posterior sampling (sample_loadings)."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest
import scipy.linalg

from longbet._multi_step import sample_loadings


def test_loading_posterior_vs_analytic():
    """Verify sample_loadings matches analytic Gaussian posterior mean and covariance.

    Uses non-orthogonal predictors, non-unit sigma2, and informative prior.
    Non-orthogonal X distinguishes L vs L.T (Cholesky orientation).
    """
    Q = 200
    m = 2
    rng = np.random.default_rng(42)

    # Correlated predecessors to form non-orthogonal X
    cov_x = np.array([[1.0, 0.6], [0.6, 1.5]])
    X_pred = rng.multivariate_normal([0.0, 0.0], cov_x, size=Q).astype(np.float32)

    # True loadings
    true_gamma = np.array([0.7, -0.4], dtype=np.float32)
    sigma2_val = 1.8
    noise = rng.normal(0.0, np.sqrt(sigma2_val), size=Q).astype(np.float32)
    y_resp = X_pred @ true_gamma + noise

    # Observations mask (drop a few cells)
    obs_mask = np.ones(Q, dtype=bool)
    obs_mask[::10] = False  # 10% missing

    prior_var = 1.5

    # Analytic posterior:
    # Prior precision: Lambda_0 = (1 / prior_var) * I
    # Likelihood precision: X_obs^T X_obs / sigma2
    X_obs = X_pred[obs_mask]
    y_obs = y_resp[obs_mask]

    Lambda_0 = (1.0 / prior_var) * np.eye(m)
    XtX = X_obs.T @ X_obs
    Xty = X_obs.T @ y_obs

    Lambda_post = Lambda_0 + XtX / sigma2_val
    cov_post = np.linalg.inv(Lambda_post)
    mu_post = cov_post @ (Xty / sigma2_val)

    # Mock predecessor states
    class MockState:
        def __init__(self, resid, mask):
            self.resid = jnp.asarray(resid)
            self.obs_mask = jnp.asarray(mask)

    preds = [
        MockState(X_pred[:, 0], obs_mask),
        MockState(X_pred[:, 1], obs_mask),
    ]
    resp = MockState(y_resp, obs_mask)

    # Draw many samples with sample_loadings
    num_samples = 4000
    key = jax.random.key(123)
    keys = jax.random.split(key, num_samples)

    draws = []
    for k in keys:
        d = sample_loadings(
            k,
            predecessors=preds,
            response=resp,
            sigma2=jnp.asarray(sigma2_val, dtype=jnp.float32),
            prior_var=prior_var,
        )
        draws.append(np.asarray(d))

    draws_arr = np.array(draws)  # shape (num_samples, m)
    emp_mean = np.mean(draws_arr, axis=0)
    emp_cov = np.cov(draws_arr, rowvar=False)

    # Monte Carlo standard error for mean
    mcse = np.sqrt(np.diag(cov_post) / num_samples)
    np.testing.assert_allclose(emp_mean, mu_post, atol=4.0 * np.max(mcse))
    # Covariance within 15% relative tolerance
    np.testing.assert_allclose(emp_cov, cov_post, rtol=0.15)


def test_loading_zero_prior_or_empty_preds():
    key = jax.random.key(0)

    class MockState:
        def __init__(self, n):
            self.resid = jnp.ones(n, dtype=jnp.float32)
            self.obs_mask = jnp.ones(n, dtype=bool)

    resp = MockState(10)

    # Empty predecessors (m = 0) returns empty array
    d0 = sample_loadings(key, predecessors=[], response=resp, sigma2=jnp.array(1.0), prior_var=1.0)
    assert d0.shape == (0,)

    # Zero prior_var returns zeros
    p1 = MockState(10)
    dz = sample_loadings(key, predecessors=[p1], response=resp, sigma2=jnp.array(1.0), prior_var=0.0)
    assert dz.shape == (1,)
    np.testing.assert_allclose(dz, 0.0)


def test_loading_m3_collinear():
    """Verify M=3 and collinear predecessors are handled stably by lower Cholesky with prior."""
    Q = 100
    rng = np.random.default_rng(99)
    x0 = rng.standard_normal(Q).astype(np.float32)
    x1 = (2.0 * x0).astype(np.float32)  # exactly collinear with x0
    x2 = rng.standard_normal(Q).astype(np.float32)
    y = (0.5 * x0 + 0.1 * x2 + rng.standard_normal(Q) * 0.5).astype(np.float32)

    class MockState:
        def __init__(self, arr):
            self.resid = jnp.asarray(arr)
            self.obs_mask = jnp.ones(Q, dtype=bool)

    preds = [MockState(x0), MockState(x1), MockState(x2)]
    resp = MockState(y)

    key = jax.random.key(7)
    draw = sample_loadings(key, predecessors=preds, response=resp, sigma2=jnp.array(0.25), prior_var=1.0)
    assert draw.shape == (3,)
    assert np.all(np.isfinite(draw))
