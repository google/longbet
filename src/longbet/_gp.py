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

"""Gaussian process trajectory prior, sampling, and extrapolation for beta.

Numerical strategy
------------------
The squared-exponential Gram matrix on the integer exposure grid is severely
ill-conditioned: at ``lambda_knl = 3`` over ``S_max = 51`` its condition number
exceeds ``1e7``, so the *explicit inverse* cannot be carried in float32 without
destroying the prior.  We therefore never form ``K_tilde^{-1}``.

Instead the kernel is built and factorised in float64 with NumPy/SciPy once, at
initialisation, and only the Cholesky factor ``L`` (``L L^T = K_tilde``) crosses
into JAX.  The conditional draw uses the whitened parameterisation
``beta = L u``, whose posterior precision is ``I + L^T diag(A) L`` -- a matrix
with eigenvalues bounded below by 1, hence well conditioned in float32 whatever
``K_tilde`` looks like.  ``L``'s entries are ``O(sigma_knl)`` rather than
``O(1 / jitter)``, so casting it to float32 costs almost nothing.

Measured on the prior (``A = C = 0``), maximum relative error of the sampled
covariance against ``K_tilde``, at ``S_max = 51`` and ``lambda_knl = 3``:
carrying ``K_tilde^{-1}`` in float32 gives 0.49; this parameterisation gives
0.034, against a Monte Carlo noise floor of 0.05.
"""

from __future__ import annotations

from typing import Literal

import jax
import jax.numpy as jnp
import jax.scipy.linalg as jla
import numpy as np
import scipy.linalg as sla
from jaxtyping import Array, Float32, Key

KernelType = Literal["se", "matern32", "matern52", "ar1"]

#: Kernel families accepted by :func:`build_kernel_matrix`.  ``"matern"`` is
#: retained as an alias for ``"matern32"``.
KERNEL_TYPES: tuple[str, ...] = ("se", "matern", "matern32", "matern52", "ar1")


def build_kernel_matrix(
    s: np.ndarray | Array,
    sig_knl: float,
    lambda_knl: float,
    kernel_type: str = "se",
    sigma_m: float = 1.0,
    gp_constant_mean: bool = True,
    jitter: float = 1e-6,
) -> np.ndarray:
    """Build the marginalized GP covariance matrix over grid `s`, in float64.

    The constant mean ``m ~ N(0, sigma_m^2)`` is marginalized into the kernel as
    ``sigma_m^2 * 1 1^T`` rather than sampled, so the forecast reverts to the
    estimated common level with no ``m`` in the state (plan section 4.4).

    Parameters
    ----------
    s
        Grid of exposure values (0, 1, ..., S_max), shape (d,).
    sig_knl
        Kernel standard deviation; the marginal variance is ``sig_knl ** 2``.
    lambda_knl
        Lengthscale parameter.
    kernel_type
        'se' (squared exponential), 'matern32' (Matern-3/2, alias 'matern'),
        'matern52' (Matern-5/2), or 'ar1'.
    sigma_m
        Prior standard deviation for the constant mean m ~ N(0, sigma_m^2).
    gp_constant_mean
        If True, marginalizes out m by adding sigma_m^2 * 1 1^T.
    jitter
        Relative jitter added to the diagonal: ``jitter * sig_knl^2 * I``.

    Returns
    -------
    K_tilde : np.ndarray
        Symmetric marginalized kernel matrix in float64, with diagonal jitter.
    """
    s_arr = np.asarray(s, dtype=np.float64)
    dist = np.abs(s_arr[:, None] - s_arr[None, :])
    var = float(sig_knl) ** 2

    if lambda_knl <= 0:
        raise ValueError(f"lambda_knl must be positive, got {lambda_knl}")

    if kernel_type == "se":
        K = var * np.exp(-0.5 * (dist / lambda_knl) ** 2)
    elif kernel_type in ("matern", "matern32"):
        scaled = np.sqrt(3.0) * dist / lambda_knl
        K = var * (1.0 + scaled) * np.exp(-scaled)
    elif kernel_type == "matern52":
        scaled = np.sqrt(5.0) * dist / lambda_knl
        K = var * (1.0 + scaled + scaled**2 / 3.0) * np.exp(-scaled)
    elif kernel_type == "ar1":
        rho = np.exp(-1.0 / lambda_knl)
        K = var * (rho**dist)
    else:
        raise ValueError(
            f"Unknown kernel type: {kernel_type!r}; expected one of {KERNEL_TYPES}"
        )

    d = s_arr.shape[0]
    jitter_val = max(float(jitter) * var, 1e-8)
    K = K + jitter_val * np.eye(d, dtype=np.float64)

    if gp_constant_mean:
        K = K + float(sigma_m) ** 2 * np.ones((d, d), dtype=np.float64)

    return 0.5 * (K + K.T)


def kernel_cholesky(K_tilde: np.ndarray | Array) -> np.ndarray:
    """Lower Cholesky factor of ``K_tilde``, computed in float64.

    Returns
    -------
    L : np.ndarray
        Lower-triangular float64 array with ``L @ L.T == K_tilde``.
    """
    K_np = np.asarray(K_tilde, dtype=np.float64)
    K_np = 0.5 * (K_np + K_np.T)
    return sla.cholesky(K_np, lower=True)


def compute_kernel_precision(K_tilde: np.ndarray | Array) -> np.ndarray:
    """Invert ``K_tilde`` in float64 via its Cholesky factor.

    Retained for diagnostics and tests. The sampler deliberately does **not**
    use this: an explicit precision matrix cannot survive the cast to float32
    (see the module docstring).
    """
    K_np = np.asarray(K_tilde, dtype=np.float64)
    L = kernel_cholesky(K_np)
    K_inv = sla.cho_solve((L, True), np.eye(K_np.shape[0], dtype=np.float64))
    return 0.5 * (K_inv + K_inv.T)


def sample_beta_gp(
    key: Key[Array, ''],
    A: Float32[Array, '*chains d'],
    C: Float32[Array, '*chains d'],
    K_chol: Float32[Array, 'd d'],
) -> Float32[Array, '*chains d']:
    """Draw ``beta`` from its Gaussian conditional posterior.

    The target is ``N(P^-1 C, P^-1)`` with ``P = K_tilde^-1 + diag(A)``. Writing
    ``beta = L u`` with ``L L^T = K_tilde`` turns this into

    .. code::

        M      = I + L^T diag(A) L
        u      ~ N(M^-1 L^T C, M^-1)
        beta   = L u

    which is algebraically identical but never forms ``K_tilde^-1``. ``M`` has
    eigenvalues at least 1, so its Cholesky is stable in float32 regardless of
    how ill-conditioned ``K_tilde`` is.

    Exposure values with no observations have ``A_s = C_s = 0`` and are drawn
    from the GP conditional prior, which is exactly right and needs no special
    case.

    Parameters
    ----------
    key
        PRNG key.
    A
        Diagonal precision contributions from the data, ``sum d_it^2 / sigma^2``.
    C
        Mean contributions from the data, ``sum d_it r_it / sigma^2``.
    K_chol
        Lower Cholesky factor of the marginalized kernel, float32, shape (d, d).

    Returns
    -------
    beta : Float32[Array, '*chains d']
        Sampled exposure trajectory, float32.
    """
    L = jnp.asarray(K_chol, dtype=jnp.float32)
    d = L.shape[-1]
    Lt = jnp.swapaxes(L, -1, -2)

    # M = I + L^T diag(A) L. Scaling L's rows by A gives diag(A) @ L.
    M = jnp.eye(d, dtype=jnp.float32) + jnp.matmul(Lt, A[..., :, None] * L)
    M = 0.5 * (M + jnp.swapaxes(M, -1, -2))

    LM = jla.cholesky(M, lower=True)

    mean_u = jla.cho_solve((LM, True), jnp.matmul(Lt, C[..., None])).squeeze(-1)
    eta = jax.random.normal(key, shape=C.shape, dtype=jnp.float32)
    noise_u = jla.solve_triangular(LM, eta[..., None], lower=True, trans=1).squeeze(-1)

    u = mean_u + noise_u
    return jnp.matmul(L, u[..., None]).squeeze(-1).astype(jnp.float32)


def beta_prior_quadform(
    beta: Float32[Array, '*chains d'],
    K_chol: Float32[Array, 'd d'],
) -> Float32[Array, '*chains']:
    """Return ``beta^T K_tilde^{-1} beta``, computed by triangular solve.

    Used by the ridge move. Solving ``L v = beta`` and returning ``v^T v`` is
    stable where multiplying by an explicit ``K_tilde^{-1}`` is not.
    """
    L = jnp.asarray(K_chol, dtype=jnp.float32)
    v = jla.solve_triangular(L, beta[..., None], lower=True).squeeze(-1)
    return jnp.sum(jnp.square(v), axis=-1)


def forecast_beta_gp(
    key: Key[Array, ''],
    beta_obs: Float32[Array, '*batch d_obs'],
    s_obs: np.ndarray | Array,
    s_fut: np.ndarray | Array,
    sig_knl: float,
    lambda_knl: float,
    kernel_type: str = "se",
    sigma_m: float = 1.0,
    gp_constant_mean: bool = True,
    jitter: float = 1e-6,
    noise_scale: np.ndarray | Array | None = None,
) -> Float32[Array, '*batch d_total']:
    """Project ``beta`` onto future exposure values, conditional on ``beta_obs``.

    Computes, per draw,

    .. math::

        \\beta_* \\mid \\beta \\sim
        N(K_{*,\\cdot} K^{-1} \\beta,\\; K_{*,*} - K_{*,\\cdot} K^{-1} K_{\\cdot,*})

    entirely in float64 NumPy/SciPy, over the *joint* index set so that the
    marginalized constant mean applies across both blocks and the projection
    reverts to the estimated common level rather than to zero.

    Each draw receives its own noise vector from a single split of ``key``, so
    the projection is reproducible from the model's seed.

    ``noise_scale`` (shape ``*batch``) multiplies each draw's conditional
    noise: with a learned trajectory scale the kernel is ``s^2 K_tilde``, whose
    conditional mean does not depend on ``s`` and whose conditional covariance
    is ``s^2`` times that of ``K_tilde``.

    Returns the concatenated trajectory ``[beta_obs, beta_*]`` as float32.
    """
    s_obs_np = np.asarray(s_obs, dtype=np.float64)
    s_fut_np = np.asarray(s_fut, dtype=np.float64)
    if s_fut_np.shape[0] == 0:
        return jnp.asarray(beta_obs, dtype=jnp.float32)

    d_obs = s_obs_np.shape[0]
    d_fut = s_fut_np.shape[0]

    K_all = build_kernel_matrix(
        np.concatenate([s_obs_np, s_fut_np]),
        sig_knl=sig_knl,
        lambda_knl=lambda_knl,
        kernel_type=kernel_type,
        sigma_m=sigma_m,
        gp_constant_mean=gp_constant_mean,
        jitter=jitter,
    )

    K11 = K_all[:d_obs, :d_obs]
    K21 = K_all[d_obs:, :d_obs]
    K22 = K_all[d_obs:, d_obs:]

    L11 = sla.cholesky(K11, lower=True)
    W_T = sla.cho_solve((L11, True), K21.T)  # K11^{-1} K12, shape (d_obs, d_fut)

    beta_obs_np = np.asarray(beta_obs, dtype=np.float64)
    mu_fut = beta_obs_np @ W_T  # (*batch, d_fut)

    Sigma_fut = K22 - K21 @ W_T
    jitter_val = max(float(jitter) * float(sig_knl) ** 2, 1e-8)
    Sigma_fut = 0.5 * (Sigma_fut + Sigma_fut.T) + jitter_val * np.eye(d_fut)
    L_fut = sla.cholesky(Sigma_fut, lower=True)

    eta = np.asarray(
        jax.random.normal(key, shape=(*beta_obs_np.shape[:-1], d_fut), dtype=jnp.float32),
        dtype=np.float64,
    )
    noise = eta @ L_fut.T
    if noise_scale is not None:
        noise = noise * np.asarray(noise_scale, dtype=np.float64)[..., None]
    beta_fut = mu_fut + noise

    return jnp.asarray(np.concatenate([beta_obs_np, beta_fut], axis=-1), dtype=jnp.float32)
