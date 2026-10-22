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

"""Geweke joint-distribution ("getting it right") test of the Gibbs sweep.

This checks the full transition against known prior marginals. Analytical
conditional/posterior-moment tests provide independent checks: good treatment
effect means alone can conceal incorrect conditional variances.

Method.  The successive-conditional simulator alternates

    1. draw data given the current parameters, and
    2. run one Gibbs sweep,

whose stationary distribution is the model's joint prior.  The marginal of any
parameter block under that joint is simply its own prior, which for
``sigma^2``, ``sigma_gamma^2``, ``gamma`` and ``beta`` is available in closed
form.  So the simulated marginals are compared against those analytic priors,
not against loose interval bounds.

Two different criteria are used, deliberately:

* For the **correct** sampler, deviations are measured in Monte Carlo standard
  errors computed from the *effective* sample size, since successive-conditional
  draws are autocorrelated.
* For a **broken** sampler that criterion is the wrong one -- a badly mixing
  chain has a tiny effective sample size, which inflates its own standard error
  and hides the bias.  The sensitivity test below therefore compares the point
  estimate against the correct sampler's, which is the comparison a human would
  make.
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

# Proper, informative priors: the reference prior on sigma^2 is improper and has
# no distribution to compare against.
GAMMA_A, GAMMA_B = 4.0, 2.0
SIGMA_A, SIGMA_B = 5.0, 3.0
SIG_KNL, LAMBDA_KNL, SIGMA_M = 1.0, 1.5, 0.8

N_ITER = 3000
N_WARMUP = 500
#: Deviation from the analytic prior mean, in Monte Carlo standard errors.
#: Four covers several simultaneous statistics without being lax: a conditional
#: off by a constant lands well beyond it.
Z_TOL = 4.0


def _build_state():
    """A deliberately small panel, so the prior rather than the data dominates."""
    N, T, P = 8, 4, 2
    M = N * T
    rng = np.random.default_rng(0)
    exposure = np.tile(np.arange(T), N).astype(np.int32)
    config = LongBetConfig(
        num_trees_pr=2,
        num_trees_trt=2,
        random_intercept=True,
        gamma_prior_a=GAMMA_A,
        gamma_prior_b=GAMMA_B,
        sigma_prior_a=SIGMA_A,
        sigma_prior_b=SIGMA_B,
        sig_knl=SIG_KNL,
        lambda_knl=LAMBDA_KNL,
        sigma_m=SIGMA_M,
        gp_constant_mean=True,
    )
    X_unit = rng.integers(0, 3, (P, N), dtype=np.uint8)
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
    )
    return state, config


def _whitener(state, config):
    """``L^-1`` for ``beta``'s prior, so the target becomes ``sign(v' u) u`` with ``u ~ N(0, I)``."""
    K = build_kernel_matrix(
        np.arange(state.S_max + 1), SIG_KNL, LAMBDA_KNL, config.kernel_type,
        SIGMA_M, config.gp_constant_mean,
    )
    return np.linalg.inv(np.linalg.cholesky(K))


def _canonical_beta_white0_mean(state, config) -> float:
    """Exact prior mean of ``tilde_u[0] = sign(sum_{s>=1} beta_s) * u[0]`` for ``u ~ N(0, I)``."""
    K = build_kernel_matrix(
        np.arange(state.S_max + 1), SIG_KNL, LAMBDA_KNL, config.kernel_type,
        SIGMA_M, config.gp_constant_mean,
    )
    L = np.linalg.cholesky(K)
    c = np.ones(state.S_max + 1, dtype=np.float64)
    c[0] = 0.0
    v = L.T @ c
    return float(np.sqrt(2.0 / np.pi) * v[0] / np.linalg.norm(v))


def _run_successive_conditional(seed, n_iter=N_ITER, warmup=N_WARMUP, beta_scale=1.0):
    """Alternate forward data generation and one Gibbs sweep.

    ``beta_scale != 1`` deliberately corrupts the ``beta`` conditional, which the
    sensitivity test uses as a known-bad reference.
    """
    import longbet._step as step_mod

    state, config = _build_state()
    L_inv = _whitener(state, config)
    original = step_mod.sample_beta_gp

    if beta_scale != 1.0:
        step_mod.sample_beta_gp = lambda k, A, C, L: original(k, A, C, beta_scale * L)
    try:
        # A fresh closure: jax.jit caches traces by function identity, and the
        # module-level sweep has already been traced with the real conditional
        # bound in, so re-using it would silently ignore the patch.
        # Unwrap every layer: the x64 context wrapper sits outside jax.jit,
        # so a single __wrapped__ would still call the previously cached JIT.
        body = inspect.unwrap(longbet_single_step)

        @jax.jit
        @enable_x64(False)
        def sweep(k, s):
            return body(k, s)

        key = jax.random.key(seed)
        out = {"sigma2": [], "sigma_gamma2": [], "gamma_std": [], "beta_white": []}
        for it in range(n_iter):
            key, sub = jax.random.split(key)
            k_data, k_step = jax.random.split(sub)

            b_z = jnp.where(state.z_vec == 1.0, state.b1, state.b0)
            mean = (
                state.alpha * state.mu_fit
                + b_z * state.beta[state.exposure_idx] * state.nu_fit
                + state.gamma[state.unit_idx]
            )
            eps = jax.random.normal(k_data, mean.shape, jnp.float32) * jnp.sqrt(state.sigma2)
            state = eqx.tree_at(lambda s: (s.y, s.resid), state, (mean + eps, eps))
            state = sweep(k_step, state)

            if it < warmup:
                continue
            sg2 = float(state.sigma_gamma2)
            out["sigma2"].append(float(state.sigma2))
            out["sigma_gamma2"].append(sg2)
            # gamma_i | sigma_gamma^2 ~ N(0, sigma_gamma^2): standard normal.
            out["gamma_std"].append(np.asarray(state.gamma) / np.sqrt(sg2))
            out["beta_white"].append(L_inv @ np.asarray(state.beta, dtype=np.float64))
    finally:
        step_mod.sample_beta_gp = original

    return {
        "sigma2": np.array(out["sigma2"]),
        "sigma_gamma2": np.array(out["sigma_gamma2"]),
        "gamma_std": np.concatenate(out["gamma_std"]),
        "beta_white": np.stack(out["beta_white"]),
    }


def _z_against(sample, prior_mean, label):
    """Standardized deviation of a sample mean from a known prior mean."""
    sample = np.asarray(sample, dtype=np.float64)
    ess = max(float(compute_ess(sample[None, :])), 2.0)
    mcse = float(np.std(sample, ddof=1)) / np.sqrt(ess)
    return abs(float(np.mean(sample)) - prior_mean) / max(mcse, 1e-12), ess, mcse, label


@pytest.mark.slow
def test_geweke_joint_distribution():
    state_ref, config_ref = _build_state()
    canonical_mean0 = _canonical_beta_white0_mean(state_ref, config_ref)
    draws = _run_successive_conditional(seed=20260906)

    checks = [
        _z_against(draws["sigma2"], SIGMA_B / (SIGMA_A - 1.0), "sigma^2 mean"),
        _z_against(1.0 / draws["sigma2"], SIGMA_A / SIGMA_B, "1/sigma^2 mean"),
        _z_against(draws["sigma_gamma2"], GAMMA_B / (GAMMA_A - 1.0), "sigma_gamma^2 mean"),
        _z_against(1.0 / draws["sigma_gamma2"], GAMMA_A / GAMMA_B, "1/sigma_gamma^2 mean"),
        _z_against(draws["gamma_std"], 0.0, "standardized gamma mean"),
        _z_against(draws["beta_white"][:, 0], canonical_mean0, "canonical whitened beta[0] mean"),
        _z_against(np.square(draws["beta_white"][:, 0]), 1.0, "whitened beta[0] second moment"),
    ]
    failures = [
        f"{label}: z = {z:.2f} (ESS {ess:.0f}, MCSE {mcse:.4g})"
        for z, ess, mcse, label in checks
        if not np.isfinite(z) or z > Z_TOL
    ]
    assert not failures, "Geweke means disagree with the prior:\n  " + "\n  ".join(failures)

    # Distributional shape, not just the first two moments. Thin hard, since
    # successive-conditional draws are autocorrelated and KS assumes independence.
    p_gamma = st.kstest(draws["gamma_std"][::37], "norm").pvalue
    assert p_gamma > 1e-3, f"standardized gamma is not N(0,1): KS p = {p_gamma:.2e}"

    p_sigma = st.kstest(
        draws["sigma2"][::37], lambda q: st.invgamma.cdf(q, a=SIGMA_A, scale=SIGMA_B)
    ).pvalue
    assert p_sigma > 1e-3, (
        f"sigma^2 is not Inv-Gamma({SIGMA_A}, {SIGMA_B}): KS p = {p_sigma:.2e}"
    )


@pytest.mark.slow
def test_geweke_detects_a_broken_conditional():
    """The harness must reject a wrong sampler, or passing it proves nothing.

    The mutation is realistic and non-compounding: the ``beta`` conditional is
    handed a kernel factor 50 percent too large, as a units slip would. Nothing
    diverges and the draws remain plausibly scaled, but the stationary marginal
    of ``beta`` is no longer its prior. Treatment-effect recovery alone would
    not directly test this conditional.

    The comparison is against the correct sampler's own estimate rather than
    against a Monte Carlo standard error, because a broken sampler mixes badly
    and its effective sample size collapses, which would inflate its standard
    error and hide the very bias being looked for.
    """
    # Keep the original seed, 50% mutation and rejection thresholds. The
    # corrected precision cache changes the trajectory: 3,000 sweeps gave
    # second moments .822 (correct) and 1.450 (mutated), detecting the mutation
    # but leaving too much Monte Carlo error for the strict deviation ratio.
    # Use more sweeps to reduce that Monte Carlo error; keep the ratio unchanged.
    good = _run_successive_conditional(seed=11, n_iter=12000, warmup=2000)
    bad = _run_successive_conditional(seed=11, n_iter=12000, warmup=2000,
                                      beta_scale=1.5)

    m_good = float(np.mean(np.square(good["beta_white"][:, 0])))
    m_bad = float(np.mean(np.square(bad["beta_white"][:, 0])))

    # The target is 1.0 by construction of the whitening.
    err_good, err_bad = abs(m_good - 1.0), abs(m_bad - 1.0)
    assert err_good < 0.25, (
        f"the correct sampler should sit near the prior: E[w^2] = {m_good:.3f}"
    )
    assert err_bad > 0.4, (
        f"the broken beta conditional was not detected: E[w^2] = {m_bad:.3f}, "
        f"correct sampler gives {m_good:.3f}. The Geweke harness is not sensitive "
        "enough to be worth running."
    )
    assert err_bad > 4.0 * err_good, (
        f"broken/correct deviation ratio is only {err_bad / max(err_good, 1e-9):.1f}"
    )
