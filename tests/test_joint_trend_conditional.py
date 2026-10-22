"""Independent dense Gaussian oracle for the joint calendar/exposure/unit draw.

The oracle assembles one full design and inverts its posterior precision; it
neither eliminates unit intercepts nor calls the sampler's sufficient statistics.
Controlling the Gaussian innovations identifies mean and covariance exactly,
without Monte Carlo error or a tolerance that hides a precision ridge.
"""
from unittest.mock import patch

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from scipy import linalg

from longbet._trend_basis import trend_block_step, trend_sufficient_statistics
from longbet._x64 import enable_x64


@pytest.mark.parametrize("sigma2,cached,weighted,units,trajectory", [
    (0.1, False, False, True, True),
    (0.1, True, False, True, True),
    (1e-6, False, False, True, True),
    (1e-6, True, False, True, True),
    (0.1, False, True, True, True),
    (0.1, False, False, False, True),
    (0.1, True, False, True, False),
    (0.1, False, True, False, False),
])
def test_joint_conditional_matches_dense_oracle(sigma2, cached, weighted, units, trajectory):
    rng = np.random.default_rng(812)
    N, T, q, S = 9, 6, 10, 4
    n = N * T
    ui = np.repeat(np.arange(N), T).astype(np.int32)
    # Deliberately mix unit-constant and varying columns: eliminating the
    # former is the cancellation case that the previous ridge distorted.
    H = rng.normal(size=(n, q)).astype(np.float32)
    H[:, :4] = np.repeat(rng.normal(size=(N, 4)), T, axis=0)
    sc = np.linspace(.3, 2., q).astype(np.float32)
    exposure = np.tile([0, 0, 0, 1, 2, 3], N).astype(np.int32)
    exposure[ui >= 6] = 0
    d = ((exposure > 0) * rng.uniform(.4, 1.4, n)).astype(np.float32)
    mask = np.ones(n, bool)
    mask[[3, 9, 27, 28]] = False
    A = rng.normal(size=(S, S))
    L = np.linalg.cholesky(A @ A.T + np.eye(S)).astype(np.float32)
    y = rng.normal(size=n).astype(np.float32)
    c = rng.normal(size=q).astype(np.float32)
    b = rng.normal(size=S).astype(np.float32)
    g = rng.normal(size=N).astype(np.float32)
    R = (y - H @ c - (d * b[exposure] if trajectory else 0)
         - (g[ui] if units else 0)).astype(np.float32)
    sigma2, gs2, temp = map(np.float32, (sigma2, .2, .7))
    precision = (rng.uniform(.3, 3., n).astype(np.float32) if weighted
                 else np.full(n, float(temp) / float(sigma2)))
    precision = np.where(mask, precision, 0.)
    # Match the actual float32 inputs and residual restoration, in float64.
    restored = R.astype(float) + H.astype(float) @ c
    if units:
        restored += g[ui]
    if trajectory:
        restored += d.astype(float) * b[exposure]
    blocks = [H.astype(float) * sc]
    if trajectory:
        blocks.append(d.astype(float)[:, None] * L[exposure].astype(float))
    if units:
        blocks.append(np.eye(N)[ui] * np.sqrt(float(gs2)))
    X = np.column_stack(blocks)
    P = np.eye(X.shape[1]) + X.T @ (precision[:, None] * X)
    chol = linalg.cholesky(P, lower=True)
    covariance = linalg.cho_solve((chol, True), np.eye(len(P)))
    mean = linalg.cho_solve((chol, True), X.T @ (precision * restored))
    args = dict(trend_coef=c, gamma=g, beta=b, R=R, H=H, d_vec=d,
                obs_mask=mask, unit_idx=ui, N_units=N, exposure_idx=exposure,
                sigma2=sigma2, sigma_gamma2=gs2, K_chol=L,
                trend_col_scale=sc, random_intercept=units, sample_beta=trajectory,
                temperature=temp,
                conditional_precision=precision.astype(np.float32) if weighted else None)
    if cached:
        args['trend_gram'], args['trend_unit_sum'] = trend_sufficient_statistics(
            jnp.asarray(H), jnp.asarray(mask), jnp.asarray(ui), N)
        assert args['trend_gram'].dtype == np.float64
    k = q + (S if trajectory else 0)

    @enable_x64()
    def draw(noise):
        with patch('jax.random.normal', side_effect=[noise[:k], noise[k:]]):
            cn, gn, bn, rn = trend_block_step(jax.random.key(0), **args)
        output = [cn.astype(jnp.float64) / sc]
        if trajectory:
            output.append(jax.scipy.linalg.solve_triangular(
                jnp.asarray(L, dtype=jnp.float64), bn.astype(jnp.float64), lower=True))
        if units:
            output.append(gn.astype(jnp.float64) / np.sqrt(float(gs2)))
        return jnp.concatenate(output)

    with enable_x64():
        noise = jnp.zeros(len(P), dtype=jnp.float64)
        actual_mean = np.asarray(draw(noise))
        factor = np.asarray(jax.jacfwd(draw)(noise), dtype=float)
    actual_covariance = factor @ factor.T
    # Check covariance in posterior-standardized directions, not only its
    # large diagonal entries (which can conceal a badly distorted direction).
    relative_covariance = chol.T @ actual_covariance @ chol
    np.testing.assert_allclose(relative_covariance, np.eye(len(P)), atol=3e-4, rtol=0)
    np.testing.assert_allclose((actual_mean - mean) / np.sqrt(np.diag(covariance)),
                               0., atol=3e-4, rtol=0)
    cn, gn, bn, rn = trend_block_step(jax.random.key(93), **args)
    expected_resid = restored - H.astype(float) @ np.asarray(cn)
    if units:
        expected_resid -= np.asarray(gn)[ui]
    if trajectory:
        expected_resid -= d.astype(float) * np.asarray(bn)[exposure]
    np.testing.assert_allclose(rn, np.where(mask, expected_resid, 0), atol=1e-6, rtol=1e-6)
    assert all(v.dtype == np.float32 for v in (cn, gn, bn, rn))
