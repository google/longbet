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

"""GP kernels, the conditional draw, and float32 hygiene.

The squared-exponential Gram matrix on the integer exposure grid is severely
ill-conditioned, and the plan flags carrying its inverse in float32 as the
project's main numerical risk. The sampler avoids it entirely by working through
the Cholesky factor; these tests pin that it is both correct and better.
"""

import jax
import jax.numpy as jnp
import numpy as np
import pytest
import scipy.linalg as sla

from longbet._gp import (
    KERNEL_TYPES,
    beta_prior_quadform,
    build_kernel_matrix,
    compute_kernel_precision,
    forecast_beta_gp,
    kernel_cholesky,
    sample_beta_gp,
)


@pytest.mark.parametrize("kernel", ["se", "matern", "matern32", "matern52", "ar1"])
def test_kernels_are_symmetric_positive_definite(kernel):
    K = build_kernel_matrix(np.arange(12), sig_knl=1.0, lambda_knl=2.0, kernel_type=kernel)
    assert K.dtype == np.float64
    assert np.allclose(K, K.T)
    assert np.all(np.linalg.eigvalsh(K) > 0)


def test_matern_alias():
    a = build_kernel_matrix(np.arange(8), 1.0, 2.0, "matern")
    b = build_kernel_matrix(np.arange(8), 1.0, 2.0, "matern32")
    assert np.array_equal(a, b)
    assert set(KERNEL_TYPES) == {"se", "matern", "matern32", "matern52", "ar1"}


def test_unknown_kernel_is_refused():
    with pytest.raises(ValueError, match="Unknown kernel type"):
        build_kernel_matrix(np.arange(4), 1.0, 1.0, "rbf")


def _prior_covariance_error(sampler, K, n_draws=4000):
    d = K.shape[0]
    A = jnp.zeros(d, jnp.float32)
    C = jnp.zeros(d, jnp.float32)
    draws = np.stack([np.asarray(sampler(jax.random.key(i), A, C)) for i in range(n_draws)])
    return np.abs(np.cov(draws.T) - K).max() / np.abs(K).max()


def test_prior_draw_matches_the_kernel_where_float32_inverses_fail():
    """With no data the draw must reproduce the prior exactly.

    At ``lambda_knl = 3`` over 52 exposure levels, carrying ``K^-1`` in float32
    -- as an explicit-precision sampler must -- destroys the prior. The whitened
    parameterisation used by :func:`sample_beta_gp` does not.
    """
    K = build_kernel_matrix(np.arange(52), 1.0, 3.0, "se", 1.0, True, 1e-6)
    L32 = jnp.asarray(kernel_cholesky(K), jnp.float32)

    err = _prior_covariance_error(lambda k, A, C: sample_beta_gp(k, A, C, L32), K)

    # Monte Carlo noise floor from an exact float64 sampler with the same count.
    ref = np.random.default_rng(0).normal(size=(4000, K.shape[0])) @ sla.cholesky(K, lower=True).T
    floor = np.abs(np.cov(ref.T) - K).max() / np.abs(K).max()

    assert err < max(2.0 * floor, 0.1), (
        f"prior covariance error {err:.3f} against a Monte Carlo floor of {floor:.3f}"
    )


def test_posterior_draw_matches_an_exact_float64_solve():
    """With data, the draw must match the analytic posterior mean and covariance."""
    d = 24
    rng = np.random.default_rng(1)
    K = build_kernel_matrix(np.arange(d), 1.0, 2.0, "se", 1.0, True, 1e-6)
    L32 = jnp.asarray(kernel_cholesky(K), jnp.float32)

    A = jnp.asarray(rng.uniform(0.0, 30.0, d), jnp.float32)
    C = jnp.asarray(rng.normal(size=d) * 5.0, jnp.float32)

    P = np.linalg.inv(K) + np.diag(np.asarray(A, np.float64))
    want_mean = np.linalg.solve(P, np.asarray(C, np.float64))
    want_cov = np.linalg.inv(P)

    draws = np.stack(
        [np.asarray(sample_beta_gp(jax.random.key(i), A, C, L32)) for i in range(4000)]
    )
    np.testing.assert_allclose(draws.mean(0), want_mean, atol=0.02)
    np.testing.assert_allclose(np.cov(draws.T), want_cov, atol=0.02)


def test_sample_beta_gp_is_float32_and_batched():
    d, chains = 6, 3
    K = build_kernel_matrix(np.arange(d), 1.0, 1.5)
    L32 = jnp.asarray(kernel_cholesky(K), jnp.float32)
    out = sample_beta_gp(
        jax.random.key(0),
        jnp.ones((chains, d), jnp.float32) * 4.0,
        jnp.zeros((chains, d), jnp.float32),
        L32,
    )
    assert out.shape == (chains, d)
    # bartz annotates its fields Float32 and runtime-checks them; a float64 leak
    # here would fail deep inside the sampler instead of at the boundary.
    assert out.dtype == jnp.float32


def test_beta_prior_quadform_matches_the_explicit_precision():
    d = 16
    K = build_kernel_matrix(np.arange(d), 1.0, 2.0)
    L32 = jnp.asarray(kernel_cholesky(K), jnp.float32)
    beta = np.random.default_rng(2).normal(size=d).astype(np.float32)
    want = float(beta @ compute_kernel_precision(K) @ beta)
    got = float(beta_prior_quadform(jnp.asarray(beta), L32))
    assert abs(got - want) / abs(want) < 1e-3


def test_forecast_reverts_to_the_estimated_common_level():
    """The marginalized constant mean must pull the projection to the fitted level.

    Not to zero: that is what ``gp_constant_mean`` is for, and it comes for free
    from putting ``sigma_m^2 11'`` inside the kernel over the joint index set.
    """
    beta_obs = jnp.full((10,), 3.0, dtype=jnp.float32)
    full = forecast_beta_gp(
        key=jax.random.key(99),
        beta_obs=beta_obs,
        s_obs=np.arange(10),
        s_fut=np.arange(10, 40),
        sig_knl=1.0,
        lambda_knl=2.0,
        sigma_m=2.0,
        gp_constant_mean=True,
    )
    assert full.shape == (40,) and full.dtype == jnp.float32
    assert jnp.allclose(full[:10], beta_obs)
    assert float(jnp.mean(full[30:])) > 1.5

    to_zero = forecast_beta_gp(
        key=jax.random.key(99),
        beta_obs=beta_obs,
        s_obs=np.arange(10),
        s_fut=np.arange(10, 40),
        sig_knl=1.0,
        lambda_knl=2.0,
        gp_constant_mean=False,
    )
    assert float(jnp.mean(to_zero[30:])) < float(jnp.mean(full[30:]))


def test_forecast_is_seeded():
    kwargs = dict(
        beta_obs=jnp.ones((6,), jnp.float32),
        s_obs=np.arange(6),
        s_fut=np.arange(6, 12),
        sig_knl=1.0,
        lambda_knl=2.0,
    )
    a = forecast_beta_gp(key=jax.random.key(7), **kwargs)
    b = forecast_beta_gp(key=jax.random.key(7), **kwargs)
    c = forecast_beta_gp(key=jax.random.key(8), **kwargs)
    assert jnp.array_equal(a, b)
    assert not jnp.array_equal(a, c)


def test_no_future_points_is_a_no_op():
    beta = jnp.ones((5,), jnp.float32)
    out = forecast_beta_gp(
        key=jax.random.key(0), beta_obs=beta, s_obs=np.arange(5), s_fut=np.array([]),
        sig_knl=1.0, lambda_knl=1.0,
    )
    assert jnp.array_equal(out, beta)
