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

"""Module-2 conditional Geweke joint-distribution tests for the distributed-lag kernel.

Tests the Stage-2 sampler given fixed ``lag_design`` (``D_lag_mat``) and fixed
``extra_design`` (``extra_cols``) across:
  - Conjugate block ``(trend_coef, extra_coef, gamma, kappa)`` in ``trend_block_step``
    and ``trend_block_step_ar1``,
  - Standalone dense ``_draw_lag_kernel`` when ``use_trend_block=False``,
  - Learned kernel scale ``lag_scale2``,
  - Linear treatment moderator ``delta`` (horseshoe),
  - Stationary AR(1) innovations ``ar1_phi`` and ``interweave_gamma``,
  - Both ``anchor_moderator=False`` (canonical half-space ``sum(kappa) >= 0``)
    and ``anchor_moderator=True`` (centered-quotient moderator ``nu_tree = 1``,
    signed ``kappa``),
  - Positive control: a corrupted kernel conditional is decisively rejected.
"""

from __future__ import annotations

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
from longbet._trend_basis import build_period_basis, build_time_basis, build_trend_design, build_unit_features
from longbet._x64 import enable_x64

GAMMA_A, GAMMA_B = 4.0, 2.0
SIGMA_A, SIGMA_B = 5.0, 3.0
LAG_SIG_KNL, LAG_LAMBDA_KNL, LAG_SIGMA_M = 0.9, 1.4, 0.7
EXTRA_PRIOR_SCALE = 0.8
Z_TOL = 4.0
PHI_MAX = 0.95


def _build_lag_state(flags: dict, *, with_extra: bool = True, split_nu: bool = True):
    N, T, P, L = 8, 4, 2, 2
    M = N * T
    rng = np.random.default_rng(10021)
    flags = dict(flags)
    use_tb = flags.pop("use_trend_block", True)
    config = LongBetConfig(
        num_trees_pr=2,
        num_trees_trt=2,
        random_intercept=True,
        gamma_prior_a=GAMMA_A,
        gamma_prior_b=GAMMA_B,
        sigma_prior_a=SIGMA_A,
        sigma_prior_b=SIGMA_B,
        lag_sig_knl=LAG_SIG_KNL,
        lag_lambda_knl=LAG_LAMBDA_KNL,
        lag_sigma_m=LAG_SIGMA_M,
        lag_gp_constant_mean=True,
        extra_cols_prior_scale=EXTRA_PRIOR_SCALE,
        trend_degree=1,
        trend_num_features=2,
        trend_horseshoe=False,
        split_calendar_mu=False,
        **flags,
    )
    X_unit = rng.integers(0, 3, (P, N), dtype=np.uint8)
    lin_units = rng.normal(size=(N, 2))
    lin_units = (lin_units - lin_units.mean(axis=0)) / lin_units.std(axis=0)

    # Fixed lag-dose matrix (N, T, L+1) with pre-panel zeros; scaled so the
    # prior rather than the data dominates on this small Geweke panel.
    dose = rng.uniform(0.0, 0.40, size=(N, T)) * (rng.random((N, T)) < 0.6)
    lags = np.zeros((N, T, L + 1), dtype=np.float32)
    for l in range(L + 1):
        if l < T:
            lags[:, l:, l] = dose[:, : T - l]
    lag_design = jnp.asarray(lags.reshape(M, L + 1), dtype=jnp.float32)

    unit_idx = np.repeat(np.arange(N, dtype=np.int32), T)
    time_idx = np.tile(np.arange(T, dtype=np.int32), N)

    trend_design = None
    trend_num_period = 0
    if use_tb:
        phi_t = build_time_basis(T, config.trend_degree)
        period_t = build_period_basis(T)
        x_raw = rng.normal(size=(N, 2)).astype(np.float32)
        features, _ = build_unit_features(x_raw, config.trend_num_features)
        trend_num_period = int(period_t.shape[1])
        trend_design = build_trend_design(features, phi_t, unit_idx, time_idx, period=period_t)

    extra_design = None
    if with_extra and use_tb:
        extra_raw = rng.normal(size=(M, 2)).astype(np.float32)
        extra_design = jnp.asarray(extra_raw, dtype=jnp.float32)

    max_split_nu = jnp.full(P, 2 if split_nu else 0, jnp.uint8)
    state = init_longbet(
        X_unified=jnp.asarray(np.repeat(X_unit, T, axis=1)),
        y=jnp.zeros(M, jnp.float32),
        unit_idx=jnp.asarray(unit_idx),
        time_idx=jnp.asarray(time_idx),
        exposure_idx=jnp.zeros(M, dtype=jnp.int32),
        z_vec=jnp.ones(M, dtype=jnp.float32),
        obs_mask=jnp.ones(M, bool),
        max_split_mu=jnp.full(P, 2, jnp.uint8),
        max_split_nu=max_split_nu,
        config=config,
        trend_design=trend_design,
        trend_num_period=trend_num_period,
        lin_design=jnp.asarray(np.repeat(lin_units, T, axis=0), jnp.float32),
        lag_design=lag_design,
        extra_design=extra_design,
    )
    return state, config


def _lag_kernel_chol(state, config):
    L1 = state.kappa.shape[0]
    K = build_kernel_matrix(
        np.arange(L1, dtype=np.float64),
        LAG_SIG_KNL,
        LAG_LAMBDA_KNL,
        config.lag_kernel_type,
        LAG_SIGMA_M,
        config.lag_gp_constant_mean,
    )
    return np.linalg.cholesky(K)


def _canonical_lag_mean0(L):
    """Exact prior mean of ``u[0]`` under the half-space projection ``sum(L u) >= 0``."""
    c = np.ones(L.shape[0], dtype=np.float64)
    v = L.T @ c
    return float(np.sqrt(2.0 / np.pi) * v[0] / np.linalg.norm(v))


def _run_lag_geweke(
    flags: dict,
    seed: int = 10027,
    n_iter: int = 4000,
    warmup: int = 500,
    with_extra: bool = True,
    split_nu: bool = True,
    corrupt_kernel_scale: float | None = None,
):
    import longbet._step as step_mod

    state, config = _build_lag_state(flags, with_extra=with_extra, split_nu=split_nu)
    L = _lag_kernel_chol(state, config)
    L_inv = np.linalg.inv(L)
    orig_factor = step_mod._lag_kernel_factor
    if corrupt_kernel_scale is not None:
        def broken_factor(st):
            return (orig_factor(st) * jnp.float32(corrupt_kernel_scale)).astype(jnp.float32)
        step_mod._lag_kernel_factor = broken_factor

    try:
        body = inspect.unwrap(longbet_single_step)
        N, T = state.N_units, state.T_periods

        @enable_x64()
        @jax.jit
        @enable_x64(False)
        def gewek_step(st, subkey):
            k_data, k_step = jax.random.split(subkey)
            nu = st.nu_fit + (st.lin_fit if st.use_linear_trt else 0.0)
            w = st.lag_design @ st.kappa
            trend = (st.trend_design @ st.trend_coef) if st.trend_design.shape[1] else 0.0
            extra = (st.extra_design @ st.extra_coef) if st.num_extra_cols else 0.0
            mean = st.alpha * st.mu_fit + trend + extra + w * nu + st.gamma[st.unit_idx]
            xi = jax.random.normal(k_data, (N, T), jnp.float32) * jnp.sqrt(st.sigma2)
            if st.ar1:
                phi = st.ar1_phi
                e0 = xi[:, 0] / jnp.sqrt(1.0 - jnp.square(phi))
                cols = [e0]
                for t_idx in range(1, T):
                    cols.append(phi * cols[-1] + xi[:, t_idx])
                eps = jnp.stack(cols, axis=1).reshape(-1)
            else:
                eps = xi.reshape(-1)
            st = eqx.tree_at(lambda s: (s.y, s.resid), st, (mean + eps, eps))
            st = body(k_step, st)
            rec = {
                "sigma2": st.sigma2,
                "sigma_gamma2": st.sigma_gamma2,
                "gamma": st.gamma,
                "kappa": st.kappa,
            }
            if st.lag_learn_scale:
                rec["lag_scale2"] = st.lag_scale2
            if st.trend_design.shape[1]:
                rec["trend_coef"] = st.trend_coef
            if st.num_extra_cols:
                rec["extra_coef"] = st.extra_coef
            if st.ar1:
                rec["phi"] = st.ar1_phi
            if st.use_linear_trt:
                rec["delta"] = st.delta
                rec["tau_delta2"] = st.tau_delta2
                rec["lin_local2"] = st.lin_local2
            return st, rec

        with enable_x64():
            keys = jax.random.split(jax.random.key(seed), n_iter)
            _, trace = jax.lax.scan(gewek_step, state, keys)
    finally:
        step_mod._lag_kernel_factor = orig_factor

    sl = slice(warmup, None)
    sigma2 = np.asarray(trace["sigma2"][sl], dtype=np.float64)
    sigma_gamma2 = np.asarray(trace["sigma_gamma2"][sl], dtype=np.float64)
    gamma = np.asarray(trace["gamma"][sl], dtype=np.float64)
    kappa = np.asarray(trace["kappa"][sl], dtype=np.float64)

    out = {
        "sigma2": sigma2,
        "sigma_gamma2": sigma_gamma2,
        "gamma_std": gamma / np.sqrt(sigma_gamma2)[:, None],
    }
    if state.lag_learn_scale:
        s = np.sqrt(np.asarray(trace["lag_scale2"][sl], dtype=np.float64))
        out["lag_scale_u"] = 2.0 / np.pi * np.arctan(s / state.lag_scale_prior)
    else:
        s = np.ones(len(sigma2), dtype=np.float64)
    out["kappa_white"] = (kappa @ L_inv.T) / s[:, None]

    if "trend_coef" in trace:
        tc = np.asarray(trace["trend_coef"][sl], dtype=np.float64)
        sc = np.asarray(state.trend_prior_scale_vec, dtype=np.float64)
        out["trend_std"] = tc / sc[None, :]
    if "extra_coef" in trace:
        ec = np.asarray(trace["extra_coef"][sl], dtype=np.float64)
        esc = np.asarray(state.extra_prior_scale_vec, dtype=np.float64)
        out["extra_std"] = ec / esc[None, :]
    if "phi" in trace:
        out["phi"] = np.asarray(trace["phi"][sl], dtype=np.float64)
    if "delta" in trace:
        tau = np.sqrt(np.asarray(trace["tau_delta2"][sl], dtype=np.float64))
        loc = np.sqrt(np.asarray(trace["lin_local2"][sl], dtype=np.float64))
        delta = np.asarray(trace["delta"][sl], dtype=np.float64)
        out["delta_std"] = delta / (tau[:, None] * loc)
        out["tau_delta_u"] = 2.0 / np.pi * np.arctan(tau / state.linear_trt_global_scale)
        out["local_u"] = 2.0 / np.pi * np.arctan(loc)
    return out, L, bool(state.anchor_moderator)


def _z(sample, target, label):
    sample = np.asarray(sample, dtype=np.float64)
    if sample.ndim > 1:
        sample = sample.mean(axis=tuple(range(1, sample.ndim)))
    ess = max(float(compute_ess(sample[None, :])), 2.0)
    mcse = float(np.std(sample, ddof=1)) / np.sqrt(ess)
    return abs(float(np.mean(sample)) - target) / max(mcse, 1e-12), ess, mcse, label


def _thinned(sample):
    a = np.asarray(sample, dtype=np.float64)
    series = a if a.ndim == 1 else a.mean(axis=tuple(range(1, a.ndim)))
    ess = max(float(compute_ess(series[None, :])), 2.0)
    step = max(37, int(np.ceil(len(series) / ess)))
    return a[::step].ravel()


def _checks(d, L, anchored: bool):
    mean0 = 0.0 if anchored else _canonical_lag_mean0(L)
    c = [
        _z(d["sigma2"], SIGMA_B / (SIGMA_A - 1.0), "sigma^2 mean"),
        _z(1.0 / d["sigma2"], SIGMA_A / SIGMA_B, "1/sigma^2 mean"),
        _z(d["sigma_gamma2"], GAMMA_B / (GAMMA_A - 1.0), "sigma_gamma^2 mean"),
        _z(1.0 / d["sigma_gamma2"], GAMMA_A / GAMMA_B, "1/sigma_gamma^2 mean"),
        _z(d["gamma_std"], 0.0, "standardized gamma mean"),
        _z(d["kappa_white"][:, 0], mean0, "whitened kappa[0] mean"),
        _z(np.square(d["kappa_white"][:, 0]), 1.0, "whitened kappa[0] second moment"),
        _z(np.square(d["kappa_white"][:, -1]), 1.0, "whitened kappa[-1] second moment"),
    ]
    if anchored:
        c.append(_z(d["kappa_white"][:, -1], 0.0, "whitened kappa[-1] mean"))
    if "lag_scale_u" in d:
        c += [
            _z(d["lag_scale_u"], 0.5, "arctan(s_lag/A) mean"),
            _z(np.square(d["lag_scale_u"]), 1.0 / 3.0, "arctan(s_lag/A)^2 mean"),
        ]
    if "trend_std" in d:
        c += [
            _z(d["trend_std"], 0.0, "standardized trend_coef mean"),
            _z(np.square(d["trend_std"]), 1.0, "standardized trend_coef second moment"),
        ]
    if "extra_std" in d:
        c += [
            _z(d["extra_std"], 0.0, "standardized extra_coef mean"),
            _z(np.square(d["extra_std"]), 1.0, "standardized extra_coef second moment"),
        ]
    if "delta_std" in d:
        c += [
            _z(d["delta_std"], 0.0, "delta/tau mean"),
            _z(np.square(d["delta_std"]), 1.0, "delta/tau second moment"),
            _z(d["tau_delta_u"], 0.5, "arctan(tau_delta/A) mean"),
            _z(np.square(d["tau_delta_u"]), 1.0 / 3.0, "arctan(tau_delta/A)^2 mean"),
        ]
    if "phi" in d:
        m2 = PHI_MAX ** 2 / 3.0
        c += [_z(d["phi"], 0.0, "phi mean"), _z(np.square(d["phi"]), m2, "phi^2 mean")]
    return c


@pytest.mark.parametrize(
    "flags,with_extra,split_nu",
    [
        (dict(use_trend_block=True, lag_learn_scale=True, linear_trt=True, anchor_moderator=False), True, True),
        (dict(use_trend_block=False, lag_learn_scale=True, linear_trt=True, anchor_moderator=False), False, True),
        (
            dict(
                use_trend_block=True,
                ar1_errors=True,
                interweave_gamma=True,
                lag_learn_scale=True,
                linear_trt=True,
                anchor_moderator=False,
            ),
            True,
            True,
        ),
        (
            dict(
                use_trend_block=True,
                ar1_errors=True,
                interweave_gamma=True,
                lag_learn_scale=True,
                linear_trt=True,
                anchor_moderator=True,
            ),
            True,
            False,
        ),
    ],
    ids=["iid_block", "iid_standalone", "ar1_block", "anchored_ar1_block"],
)
def test_geweke_lag_joint_distribution(flags, with_extra, split_nu):
    d, L, anchored = _run_lag_geweke(
        flags, seed=10027, n_iter=4000, warmup=500, with_extra=with_extra, split_nu=split_nu
    )
    failures = [
        f"{label}: z = {z:.2f} (ESS {ess:.0f}, MCSE {mcse:.4g})"
        for z, ess, mcse, label in _checks(d, L, anchored)
        if not np.isfinite(z) or z > Z_TOL
    ]
    assert not failures, "Lag-kernel Geweke means disagree with the prior:\n  " + "\n  ".join(failures)
    for name in ("lag_scale_u", "tau_delta_u", "local_u"):
        if name in d:
            p = st.kstest(_thinned(d[name]), "uniform").pvalue
            assert p > 1e-3, f"{name} is not uniform: KS p = {p:.2e}"
    for name in ("extra_std", "trend_std", "delta_std"):
        if name in d:
            p = st.kstest(_thinned(d[name]), "norm").pvalue
            assert p > 1e-3, f"{name} is not N(0,1): KS p = {p:.2e}"
    if anchored:
        p = st.kstest(_thinned(d["kappa_white"]), "norm").pvalue
        assert p > 1e-3, f"anchored whitened kappa is not N(0,1): KS p = {p:.2e}"
    if "phi" in d:
        p = st.kstest(_thinned(d["phi"]), "uniform", args=(-PHI_MAX, 2 * PHI_MAX)).pvalue
        assert p > 1e-3, f"phi is not U(-{PHI_MAX}, {PHI_MAX}): KS p = {p:.2e}"


def test_geweke_lag_detects_corrupted_kernel_update():
    """Positive control: a corrupted kernel factor in the conditional update must
    fail the Geweke battery that the uncorrupted sampler passes."""
    bad, L, anchored = _run_lag_geweke(
        dict(use_trend_block=True, ar1_errors=True, lag_learn_scale=True, anchor_moderator=True),
        seed=10029,
        n_iter=3500,
        warmup=500,
        with_extra=True,
        split_nu=False,
        corrupt_kernel_scale=2.5,
    )
    worst = max(z for z, *_ in _checks(bad, L, anchored))
    p_scale = st.kstest(_thinned(bad["lag_scale_u"]), "uniform").pvalue
    assert worst > 3 * Z_TOL and p_scale < 1e-6, (
        f"corrupted lag-kernel update not detected: largest z {worst:.1f}, scale KS p {p_scale:.1e}"
    )
