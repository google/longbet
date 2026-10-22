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

"""Unit checks of the AR(1) building blocks (``longbet._ar1``)."""

import jax
import jax.numpy as jnp
import numpy as np
from scipy import stats as st

from longbet import _ar1 as ar1
from longbet._x64 import enable_x64

N, T = 7, 5


def _first():
    return jnp.asarray(np.tile(np.arange(T), N) == 0)


def _dense_T(phi):
    """Dense Prais-Winsten matrix for one unit."""
    M = np.eye(T)
    M[0, 0] = np.sqrt(1 - phi ** 2)
    for t in range(1, T):
        M[t, t - 1] = -phi
    return M


def test_transform_matches_the_dense_matrix_and_whitens_the_ar1_covariance():
    phi = 0.7
    rng = np.random.default_rng(0)
    v = rng.normal(size=(N * T, 3))
    got = np.asarray(ar1.transform(jnp.asarray(v), jnp.asarray(phi), _first()))
    want = np.concatenate([_dense_T(phi) @ v[i * T:(i + 1) * T] for i in range(N)])
    np.testing.assert_allclose(got, want, rtol=1e-6, atol=1e-6)
    # T Sigma T' = I for the stationary AR(1) covariance with unit innovations.
    lags = np.abs(np.subtract.outer(np.arange(T), np.arange(T)))
    Sigma = phi ** lags / (1 - phi ** 2)
    Tm = _dense_T(phi)
    np.testing.assert_allclose(Tm @ Sigma @ Tm.T, np.eye(T), atol=1e-10)


@enable_x64()
def test_trend_block_precomputed_statistics_give_the_same_draw():
    rng = np.random.default_rng(1)
    n, q = N * T, 4
    unit_idx = jnp.repeat(jnp.arange(N), T)
    exposure_idx = jnp.asarray(np.minimum(np.tile(np.arange(T), N), 3))
    obs = jnp.ones(n, bool)
    first = _first()
    H = jnp.asarray(rng.normal(size=(n, q)), jnp.float32)
    L = jnp.asarray(np.linalg.cholesky(np.eye(4) + 0.5), jnp.float32)
    args = dict(trend_coef=jnp.asarray(rng.normal(size=q), jnp.float32),
                gamma=jnp.asarray(rng.normal(size=N), jnp.float32),
                beta=jnp.asarray(rng.normal(size=4), jnp.float32),
                R=jnp.asarray(rng.normal(size=n), jnp.float32), H=H,
                d_vec=jnp.asarray(rng.normal(size=n), jnp.float32), obs_mask=obs,
                unit_idx=unit_idx, N_units=N, exposure_idx=exposure_idx,
                lam=jnp.float32(1.7), sigma_gamma2=jnp.float32(0.4), K_chol=L,
                trend_col_scale=jnp.asarray(rng.uniform(0.5, 2.0, size=q), jnp.float32),
                phi=jnp.float32(0.6), first=first)
    key = jax.random.key(3)
    plain = ar1.trend_block_step_ar1(key, **args)
    stats = ar1.trend_ar1_statistics(H, obs, first, unit_idx, N)
    fast = ar1.trend_block_step_ar1(key, **args, stats=stats)
    for a, b in zip(plain, fast):
        np.testing.assert_allclose(np.asarray(a), np.asarray(b), rtol=2e-5, atol=2e-5)


def test_phi_step_targets_its_exact_conditional():
    """Many independent phi steps from one fixed residual reproduce the
    conditional computed on a fine grid (independence Metropolis is exact)."""
    rng = np.random.default_rng(2)
    phi_true, n_units, T_ = 0.8, 40, 6
    e = np.empty((n_units, T_))
    e[:, 0] = rng.normal(size=n_units) / np.sqrt(1 - phi_true ** 2)
    for t in range(1, T_):
        e[:, t] = phi_true * e[:, t - 1] + rng.normal(size=n_units)
    R = jnp.asarray(e.reshape(-1), jnp.float32)
    first = jnp.asarray(np.tile(np.arange(T_), n_units) == 0)
    obs = jnp.ones(R.shape, bool)

    @jax.jit
    def chain(key):
        def body(phi, k):
            new, _ = ar1.sample_phi(k, phi, R, first, obs, jnp.float32(1.0), jnp.float32(1.0), 0.95)
            return new, new
        _, draws = jax.lax.scan(body, jnp.float32(0.0), jax.random.split(key, 3000))
        return draws

    draws = np.asarray(chain(jax.random.key(4)))[200:]
    grid = np.linspace(-0.949, 0.949, 20001)
    stats_ = [float(v) for v in ar1.residual_stats(R, first, obs)]
    logp = np.asarray(ar1.phi_log_conditional(jnp.asarray(grid), tuple(stats_), 1.0, 1.0))
    w = np.exp(logp - logp.max()); w /= w.sum()
    mean, sd = float((grid * w).sum()), float(np.sqrt(((grid - (grid * w).sum()) ** 2 * w).sum()))
    assert abs(draws.mean() - mean) < 0.15 * sd, (draws.mean(), mean, sd)
    assert abs(draws.std() - sd) < 0.15 * sd, (draws.std(), sd)
    cdf = np.cumsum(w)
    p = st.kstest(draws[::5], lambda x: np.interp(x, grid, cdf)).pvalue
    assert p > 1e-3, p
