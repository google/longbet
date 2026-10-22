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

"""Exact nesting of the absorbing clock inside the distributed-lag kernel.

For absorbing adoption ``D_it = 1{t >= A_i}`` with ``L = S_max - 1`` the lag
doses of a cell at exposure ``s`` are ``(1, ..., 1, 0, ..., 0)`` (``s`` ones),
so ``w_it = sum_{l < s} kappa_l = beta_s`` with ``beta_{1:} = C kappa`` and
``C`` the lower-triangular matrix of ones. The clock's prior
``beta ~ N(0, K_tilde)`` on ``S_max + 1`` coordinates (``beta_0`` is a dummy
slot the likelihood never touches) induces the *marginal* prior
``beta_{1:} ~ N(0, K[1:, 1:])``; the kernel reproduces it with the non-square
factor ``L_kappa = C^{-1} L_beta[1:, :]`` of shape ``(S_max, S_max + 1)``:
``C L_kappa L_kappa' C' = L_beta[1:, :] L_beta[1:, :]' = K[1:, 1:]``, and the
whitened draw keeps all ``S_max + 1`` coordinates, so the normal streams of the
two samplers line up with no slicing. (Conditioning the dummy slot out,
``L_kappa = C^{-1} L_beta[1:, 1:]``, would instead give
``Cov(beta_{1:} | beta_0 = 0)``, a different prior.)

With that factor every whitened design row the joint block sees is
``d_it L_beta[S_it, :]`` in both samplers, so the two Markov chains are the
same chain up to floating point: the kernel stores the *differenced*
parameter in float32, so ``cumsum(kappa)`` and ``beta[1:]`` agree to float32
rounding rather than bit for bit, and tree moves may eventually take a
different branch on a rounding difference. The tests therefore check (i)
draw-for-draw agreement at float tolerance over the first sweeps on the
production outcome-equation path (anchored moderator, AR(1), trend block,
horseshoe linear moderator, interweaving), reporting the first diverging
sweep, and (ii) agreement of long-chain posterior summaries; plus (iii) the
Abel-summation identity behind the composition formula.
"""

import warnings

import jax
import jax.numpy as jnp
import numpy as np
import pytest

import longbet._model as model_mod
from longbet import LongBet, LongBetConfig
from longbet._gp import build_kernel_matrix, kernel_cholesky
from longbet._model import lag_dose_matrix

N, T = 40, 6
#: Cohort adoption periods (0-indexed; ``T`` = never), giving S_max = 4.
ADOPT = np.array([2] * 12 + [3] * 10 + [T] * 18)


def _panel(seed=10021):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, 2))
    z = (np.arange(T)[None, :] >= ADOPT[:, None]).astype(np.float32)
    s = np.cumsum(z, axis=1) * z
    tau = np.array([0.0, 0.4, 0.7, 0.9, 1.0])[s.astype(int)]
    nu = 1.0 + 0.3 * np.tanh(x[:, 1])
    y = (0.5 * x[:, :1] + 0.4 * rng.normal(size=(N, 1)) + tau * nu[:, None]
         + 0.3 * rng.normal(size=(N, T)))
    return x, z, s, y


def _outcome_config(**over):
    """The production outcome-equation configuration of ``LongBetEncourage``
    (continuous outcome): anchored moderator, AR(1), trend block, horseshoe
    linear moderator, interweaving; the clock's scale is not rescaled there."""
    base = dict(num_trees_pr=4, num_trees_trt=6, trend_num_features=4, random_seed=10023,
                anchor_moderator=True, ar1_errors=True, linear_trt=True, learn_beta_scale=True,
                interweave_gamma=True, split_calendar_mu=False, lag_learn_scale=False)
    base.update(over)
    return LongBetConfig(**base)


def _nesting_factor(cfg, s_max):
    """``L_kappa = C^{-1} L_beta[1:, :]`` from the clock's float32 factor."""
    K = build_kernel_matrix(np.arange(s_max + 1, dtype=np.float64), cfg.sig_knl, cfg.lambda_knl,
                            cfg.kernel_type, cfg.sigma_m, cfg.gp_constant_mean)
    L_beta = kernel_cholesky(K).astype(np.float32).astype(np.float64)   # as the sampler stores it
    C = np.tril(np.ones((s_max, s_max)))
    return np.linalg.solve(C, L_beta[1:, :]), L_beta


def _patched_init(monkeypatch, lag_K_chol):
    """Route the nesting factor into ``init_longbet``. The kernel's default
    start (``kappa = (1 / mean current dose, 0, ...)``) already equals the
    clock's ``beta = 1`` start on a binary absorbing dose."""
    original = model_mod.init_longbet
    monkeypatch.setattr(model_mod, "init_longbet",
                        lambda **kw: original(**kw, lag_K_chol=lag_K_chol))


def _fit_pair(monkeypatch, cfg_kw, n_burn, n_save, chains):
    x, z, s, y = _panel()
    s_max = int(s.max())
    cfg = _outcome_config(num_burnin=n_burn, num_sweeps=n_save, num_chains=chains, **cfg_kw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clock = LongBet(cfg).fit(y, x, z)
        L_kappa, L_beta = _nesting_factor(cfg, s_max)
        lags = lag_dose_matrix(z, s_max - 1)
        _patched_init(monkeypatch, L_kappa)
        kernel = LongBet(cfg).fit(y, x, exposure_lags=lags)
        monkeypatch.undo()
    # Sanity on the factor itself: the lag rows reproduce the clock's whitened rows.
    for s_val in range(1, s_max + 1):
        row = lags.reshape(-1, s_max)[np.flatnonzero(s.ravel() == s_val)[0]]
        np.testing.assert_allclose(row @ L_kappa, L_beta[s_val], rtol=0, atol=1e-12)
    return clock, kernel, (x, z, s, y)


@pytest.mark.slow
def test_lag_kernel_reproduces_absorbing_clock_draw_for_draw(monkeypatch):
    """Single chain, no burn-in, identical keys: every retained sweep of the
    kernel sampler matches the clock sampler at float32 tolerance; the first
    sweep where any quantity drifts beyond it is reported."""
    n_sweeps = 12
    clock, kernel, _ = _fit_pair(monkeypatch, {}, n_burn=0, n_save=n_sweeps, chains=1)
    tc, tk = clock.trace, kernel.trace
    beta = np.asarray(tc.beta)[:, 1:]
    cum = np.cumsum(np.asarray(tk.kappa), axis=1)
    pairs = {
        "cumsum(kappa) vs beta[1:]": (cum, beta),
        "gamma": (np.asarray(tk.gamma), np.asarray(tc.gamma)),
        "trend_coef": (np.asarray(tk.trend_coef), np.asarray(tc.trend_coef)),
        "sigma2": (np.asarray(tk.sigma2), np.asarray(tc.sigma2)),
        "sigma_gamma2": (np.asarray(tk.sigma_gamma2), np.asarray(tc.sigma_gamma2)),
        "ar1_phi": (np.asarray(tk.ar1_phi), np.asarray(tc.ar1_phi)),
        "delta": (np.asarray(tk.delta), np.asarray(tc.delta)),
        "nu leaves": (np.asarray(tk.nu_trace.leaf_tree), np.asarray(tc.nu_trace.leaf_tree)),
        "mu leaves": (np.asarray(tk.mu_trace.leaf_tree), np.asarray(tc.mu_trace.leaf_tree)),
        "nu splits": (np.asarray(tk.nu_trace.split_tree), np.asarray(tc.nu_trace.split_tree)),
    }
    first_bad = {}
    for name, (a, b) in pairs.items():
        for sweep in range(n_sweeps):
            if not np.allclose(a[sweep], b[sweep], rtol=2e-4, atol=2e-4):
                first_bad[name] = (sweep, float(np.max(np.abs(a[sweep] - b[sweep]))))
                break
    assert not first_bad, (
        "lag-kernel sampler diverged from the absorbing-clock sampler: "
        + ", ".join(f"{k}: sweep {s_} (max abs diff {d:.2e})" for k, (s_, d) in first_bad.items()))


@pytest.mark.slow
def test_lag_kernel_and_absorbing_clock_agree_in_posterior(monkeypatch):
    """Two overdispersed chains, long runs, different starts: posterior means
    of the exposure response and of the treated cells' effects agree within
    Monte Carlo error."""
    clock, kernel, (x, z, s, y) = _fit_pair(monkeypatch, {}, n_burn=300, n_save=300, chains=2)
    beta = clock.beta_draws_[:, 1:] * clock.sdy
    cum = np.cumsum(kernel.kappa_draws_, axis=1)
    mcse = np.sqrt(beta.var(0) / 50.0 + cum.var(0) / 50.0)       # conservative ESS of 50
    assert np.all(np.abs(beta.mean(0) - cum.mean(0)) < 4.0 * mcse + 0.02), (beta.mean(0), cum.mean(0))
    s_max = int(s.max())
    pc = clock.predict(x, z, summary_only=True)
    pk = kernel.predict(x, exposure_lags=lag_dose_matrix(z, s_max - 1), summary_only=True)
    treated = s >= 1
    diff = np.asarray(pc.tau_summary.mean)[treated] - np.asarray(pk.tau_summary.mean)[treated]
    sd = np.asarray(pc.tau_summary.std)[treated]
    assert np.mean(np.abs(diff) < 0.5 * sd + 0.03) > 0.95, np.abs(diff).max()


def test_abel_summation_composition_identity():
    """For binary absorbing adoption the lag-kernel composition
    ``sum_l kappa_l CITT_D(t - l)`` with ``kappa_l = tau(l + 1) - tau(l)``
    equals the clock composition ``sum_a [pi_a(1) - pi_a(0)] tau(t - a + 1)``
    (summation by parts), for arbitrary adoption probabilities and offer
    periods."""
    rng = np.random.default_rng(7)
    T_, t0 = 9, 3
    tau = np.concatenate([[0.0], np.cumsum(rng.normal(size=T_))])        # tau(0) = 0, tau(1..T)
    kappa = np.diff(tau)                                                  # kappa_l = tau(l+1) - tau(l)
    # Adoption-period probabilities under offer / no offer, zero before t0.
    pi = {zz: np.concatenate([np.zeros(t0), rng.dirichlet(np.ones(T_ - t0 + 1))[:-1]]) for zz in (0, 1)}
    citt_d = np.cumsum(pi[1] - pi[0])                                     # P(adopted by t | 1) - P(.. | 0)
    for t in range(t0, T_):
        lag_side = sum(kappa[l] * citt_d[t - l] for l in range(0, t - t0 + 1))
        clock_side = sum((pi[1][a] - pi[0][a]) * tau[t - a + 1] for a in range(t0, t + 1))
        np.testing.assert_allclose(lag_side, clock_side, rtol=1e-12, atol=1e-12)
