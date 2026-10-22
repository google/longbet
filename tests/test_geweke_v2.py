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

"""Geweke joint-distribution tests of the encourage-v2 extensions.

The successive-conditional simulator of ``test_geweke.py`` (draw data given
the parameters, then one Gibbs sweep) has the joint prior as its stationary
distribution, so every parameter's simulated marginal must match its prior.
Here the sweep runs with the moderator's linear part (``linear_trt``, a
horseshoe), the learned trajectory scale (``learn_beta_scale``), AR(1) errors,
the intercept-scale interweaving, and their combinations, and the checks add
their prior marginals: ``delta_j / (tau lambda_j) ~ N(0, 1)``, for each
half-Cauchy scale ``(2/pi) arctan(s / A) ~ U(0, 1)``, and ``phi`` uniform. The trajectory is
whitened with the *current* scale, ``u = L^-1 beta / s``. A sensitivity test
corrupts the linear coefficients' conditional and checks the harness sees it.
"""

import inspect

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import pytest
import scipy.stats as st

from longbet._config import LongBetConfig
from longbet._diagnostics import compute_ess
from longbet._gp import build_kernel_matrix
from longbet._state import init_longbet
from longbet._step import longbet_single_step
from longbet._x64 import enable_x64

GAMMA_A, GAMMA_B = 4.0, 2.0
SIGMA_A, SIGMA_B = 5.0, 3.0
SIG_KNL, LAMBDA_KNL, SIGMA_M = 1.0, 1.5, 0.8
Z_TOL = 4.0
PHI_MAX = 0.95


def _build_state(flags):
    N, T, P = 8, 4, 2
    M = N * T
    rng = np.random.default_rng(0)
    exposure = np.tile(np.arange(T), N).astype(np.int32)
    config = LongBetConfig(
        num_trees_pr=2, num_trees_trt=2, random_intercept=True,
        gamma_prior_a=GAMMA_A, gamma_prior_b=GAMMA_B,
        sigma_prior_a=SIGMA_A, sigma_prior_b=SIGMA_B,
        sig_knl=SIG_KNL, lambda_knl=LAMBDA_KNL, sigma_m=SIGMA_M,
        gp_constant_mean=True, **flags,
    )
    X_unit = rng.integers(0, 3, (P, N), dtype=np.uint8)
    lin_units = rng.normal(size=(N, 2))
    lin_units = (lin_units - lin_units.mean(axis=0)) / lin_units.std(axis=0)
    state = init_longbet(
        X_unified=jnp.asarray(np.repeat(X_unit, T, axis=1)),
        y=jnp.zeros(M, jnp.float32),
        unit_idx=jnp.repeat(jnp.arange(N), T),
        time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.asarray(exposure),
        z_vec=jnp.asarray((exposure > 0).astype(np.float32)),
        obs_mask=jnp.ones(M, bool),
        max_split_mu=jnp.full(P, 2, jnp.uint8),
        max_split_nu=jnp.full(P, 2, jnp.uint8),
        config=config,
        lin_design=jnp.asarray(np.repeat(lin_units, T, axis=0), jnp.float32),
    )
    return state, config


def _kernel_factor(state, config):
    K = build_kernel_matrix(np.arange(state.S_max + 1), SIG_KNL, LAMBDA_KNL,
                            config.kernel_type, SIGMA_M, config.gp_constant_mean)
    return np.linalg.cholesky(K)


def _canonical_mean0(L):
    c = np.ones(L.shape[0]); c[0] = 0.0
    v = L.T @ c
    return float(np.sqrt(2.0 / np.pi) * v[0] / np.linalg.norm(v))


def _run(flags, seed, n_iter=4000, warmup=500, corrupt_delta=None):
    import longbet._step as step_mod

    state, config = _build_state(flags)
    L = _kernel_factor(state, config)
    L_inv = np.linalg.inv(L)
    original = step_mod._linear_moderator_step
    if corrupt_delta is not None:
        # A units slip in the prior precision of delta (tau^2 read as c * tau^2).
        def broken(key, R, w, X, delta, lin_fit, tau2, *rest, **kw):
            return original(key, R, w, X, delta, lin_fit, corrupt_delta * tau2, *rest, **kw)
        step_mod._linear_moderator_step = broken
    try:
        body = inspect.unwrap(longbet_single_step)

        @jax.jit
        @enable_x64(False)
        def sweep(k, s):
            return body(k, s)

        key = jax.random.key(seed)
        out = {k: [] for k in ("sigma2", "sigma_gamma2", "gamma_std", "beta_white",
                               "delta_std", "tau_delta_u", "beta_scale_u", "phi", "local_u")}
        N, T = state.N_units, state.T_periods
        for it in range(n_iter):
            key, sub = jax.random.split(key)
            k_data, k_step = jax.random.split(sub)
            b_z = jnp.where(state.z_vec == 1.0, state.b1, state.b0)
            nu = state.nu_fit + (state.lin_fit if state.use_linear_trt else 0.0)
            mean = (state.alpha * state.mu_fit + b_z * state.beta[state.exposure_idx] * nu
                    + state.gamma[state.unit_idx])
            eps = jax.random.normal(k_data, mean.shape, jnp.float32) * jnp.sqrt(state.sigma2)
            if state.ar1:
                # Stationary AR(1) errors from the current phi: eps = T_phi^{-1} xi.
                phi = float(state.ar1_phi)
                xi = np.asarray(eps, dtype=np.float64).reshape(N, T)
                e = np.empty_like(xi)
                e[:, 0] = xi[:, 0] / np.sqrt(1.0 - phi ** 2)
                for t in range(1, T):
                    e[:, t] = phi * e[:, t - 1] + xi[:, t]
                eps = jnp.asarray(e.reshape(-1), jnp.float32)
            state = eqx.tree_at(lambda s: (s.y, s.resid), state, (mean + eps, eps))
            state = sweep(k_step, state)
            if it < warmup:
                continue
            sg2 = float(state.sigma_gamma2)
            out["sigma2"].append(float(state.sigma2))
            out["sigma_gamma2"].append(sg2)
            out["gamma_std"].append(np.asarray(state.gamma) / np.sqrt(sg2))
            s = np.sqrt(float(state.beta_scale2)) if state.learn_beta_scale else 1.0
            out["beta_white"].append(L_inv @ np.asarray(state.beta, dtype=np.float64) / s)
            if state.learn_beta_scale:
                out["beta_scale_u"].append(2.0 / np.pi * np.arctan(s / state.beta_scale_prior))
            if state.ar1:
                out["phi"].append(float(state.ar1_phi))
            if state.use_linear_trt:
                tau = np.sqrt(float(state.tau_delta2))
                local = np.sqrt(np.asarray(state.lin_local2, dtype=np.float64))
                out["delta_std"].append(np.asarray(state.delta, dtype=np.float64) / (tau * local))
                out["local_u"].append(2.0 / np.pi * np.arctan(local))
                out["tau_delta_u"].append(2.0 / np.pi * np.arctan(tau / state.linear_trt_global_scale))
    finally:
        step_mod._linear_moderator_step = original
    # Multi-valued blocks stay (iterations, values): the Monte Carlo error of a
    # mean is taken from the per-iteration average, whose ESS reflects the
    # chain's autocorrelation (an interleaved concatenation would not).
    res = {k: np.asarray(v) for k, v in out.items() if v}
    res["beta_white"] = np.stack(out["beta_white"])
    return res, L


def _z(sample, target, label):
    sample = np.asarray(sample, dtype=np.float64)
    if sample.ndim > 1:
        sample = sample.mean(axis=tuple(range(1, sample.ndim)))
    ess = max(float(compute_ess(sample[None, :])), 2.0)
    mcse = float(np.std(sample, ddof=1)) / np.sqrt(ess)
    return abs(float(np.mean(sample)) - target) / max(mcse, 1e-12), ess, mcse, label


def _thinned(sample):
    """Draws for a distributional test, thinned to roughly independent ones:
    every ceil(n / ESS)-th iteration (at least every 37th), ESS taken from the
    per-iteration average of a multi-valued block."""
    a = np.asarray(sample, dtype=np.float64)
    series = a if a.ndim == 1 else a.mean(axis=tuple(range(1, a.ndim)))
    ess = max(float(compute_ess(series[None, :])), 2.0)
    step = max(37, int(np.ceil(len(series) / ess)))
    return a[::step].ravel()


def _checks(d, L):
    c = [
        _z(d["sigma2"], SIGMA_B / (SIGMA_A - 1.0), "sigma^2 mean"),
        _z(1.0 / d["sigma2"], SIGMA_A / SIGMA_B, "1/sigma^2 mean"),
        _z(d["sigma_gamma2"], GAMMA_B / (GAMMA_A - 1.0), "sigma_gamma^2 mean"),
        _z(1.0 / d["sigma_gamma2"], GAMMA_A / GAMMA_B, "1/sigma_gamma^2 mean"),
        _z(d["gamma_std"], 0.0, "standardized gamma mean"),
        _z(d["beta_white"][:, 0], _canonical_mean0(L), "canonical whitened beta[0] mean"),
        _z(np.square(d["beta_white"][:, 0]), 1.0, "whitened beta[0] second moment"),
        _z(np.square(d["beta_white"][:, -1]), 1.0, "whitened beta[-1] second moment"),
    ]
    if "delta_std" in d:
        c += [_z(d["delta_std"], 0.0, "delta/tau mean"),
              _z(np.square(d["delta_std"]), 1.0, "delta/tau second moment"),
              _z(d["tau_delta_u"], 0.5, "arctan(tau_delta/A) mean"),
              _z(np.square(d["tau_delta_u"]), 1.0 / 3.0, "arctan(tau_delta/A)^2 mean")]
    if "beta_scale_u" in d:
        c += [_z(d["beta_scale_u"], 0.5, "arctan(s/A) mean"),
              _z(np.square(d["beta_scale_u"]), 1.0 / 3.0, "arctan(s/A)^2 mean")]
    if "local_u" in d:
        c += [_z(d["local_u"], 0.5, "arctan(lambda) mean"),
              _z(np.square(d["local_u"]), 1.0 / 3.0, "arctan(lambda)^2 mean")]
    if "phi" in d:
        m2 = PHI_MAX ** 2 / 3.0
        c += [_z(d["phi"], 0.0, "phi mean"), _z(np.square(d["phi"]), m2, "phi^2 mean")]
    return c


@pytest.mark.slow
@pytest.mark.parametrize("flags", [
    {},
    {"linear_trt": True},
    {"learn_beta_scale": True},
    {"linear_trt": True, "learn_beta_scale": True},
    {"ar1_errors": True},
    {"ar1_errors": True, "linear_trt": True, "learn_beta_scale": True},
    {"interweave_gamma": True},
    {"interweave_gamma": True, "linear_trt": True, "ar1_errors": True, "learn_beta_scale": True},
], ids=["frozen", "linear", "beta_scale", "both", "ar1", "ar1_all", "interweave", "all"])
def test_geweke_v2_joint_distribution(flags):
    d, L = _run(flags, seed=20260927)
    failures = [f"{label}: z = {z:.2f} (ESS {ess:.0f}, MCSE {mcse:.4g})"
                for z, ess, mcse, label in _checks(d, L) if not np.isfinite(z) or z > Z_TOL]
    assert not failures, "Geweke means disagree with the prior:\n  " + "\n  ".join(failures)
    for name in ("tau_delta_u", "beta_scale_u", "local_u"):
        if name in d:
            p = st.kstest(_thinned(d[name]), "uniform").pvalue
            assert p > 1e-3, f"{name} is not uniform: KS p = {p:.2e}"
    if "delta_std" in d:
        p = st.kstest(_thinned(d["delta_std"]), "norm").pvalue
        assert p > 1e-3, f"delta/tau is not N(0,1): KS p = {p:.2e}"
    if "phi" in d:
        p = st.kstest(_thinned(d["phi"]), "uniform", args=(-PHI_MAX, 2 * PHI_MAX)).pvalue
        assert p > 1e-3, f"phi is not U(-{PHI_MAX}, {PHI_MAX}): KS p = {p:.2e}"


@pytest.mark.slow
def test_geweke_v2_detects_a_broken_linear_conditional():
    """A prior precision off by a factor of four in the linear coefficients'
    conditional must fail the battery the correct sampler passes. The
    horseshoe's scales absorb part of the error, so the scale checks carry the
    signal (the correct sampler's largest z is below Z_TOL)."""
    bad, L = _run({"linear_trt": True}, seed=5, n_iter=6000, warmup=1000, corrupt_delta=4.0)
    worst = max(z for z, *_ in _checks(bad, L))
    p_scale = min(st.kstest(_thinned(bad[name]), "uniform").pvalue for name in ("tau_delta_u", "local_u"))
    assert worst > 3 * Z_TOL and p_scale < 1e-6, (
        f"broken delta conditional not detected: largest z {worst:.1f}, scale KS p {p_scale:.1e}")
