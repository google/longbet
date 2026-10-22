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

"""Initialization invariants of the sampler state.

The residual carried by the state must equal ``y`` minus the fitted mean at
initialization and after every sweep, because every update is incremental
(``R_new = R_old + old_component - new_component``): an offset introduced at
the start is never repaired and shifts the data the chain analyses.
"""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet._config import LongBetConfig
from longbet._state import init_longbet
from longbet._step import longbet_step

N, T, P = 12, 5, 2


def _fitted(st):
    bz = jnp.where(st.z_vec == 1.0, st.b1[..., None], st.b0[..., None])
    nu = st.nu_fit + (st.lin_fit if st.use_linear_trt else 0.0)
    trend = st.trend_coef @ st.trend_design.T if st.trend_design.shape[1] else 0.0
    if getattr(st, "use_lag_kernel", False):
        w = st.kappa @ st.lag_design.T
    else:
        w = bz * jnp.take(st.beta, st.exposure_idx, axis=-1)
    extra_design = getattr(st, "extra_design", None)
    extra = (st.extra_coef @ extra_design.T) if extra_design is not None else 0.0
    gamma = jnp.take(st.gamma, st.unit_idx, axis=-1)
    return st.alpha[..., None] * st.mu_fit + trend + w * nu + extra + gamma


def _build(anchor, ar1, num_chains, lag=False, extra=False):
    rng = np.random.default_rng(0)
    exposure = np.tile(np.arange(T), N).astype(np.int32)
    X_unit = rng.integers(0, 3, (P, N), dtype=np.uint8)
    y = jnp.asarray(rng.normal(size=N * T), jnp.float32)
    kw = dict(num_trees_pr=2, num_trees_trt=2, anchor_moderator=anchor, ar1_errors=ar1,
              learn_beta_scale=True, linear_trt=True)
    kw.update(use_trend_block=True, split_calendar_mu=False) if ar1 else kw.update(use_trend_block=False)
    cfg = LongBetConfig(**kw)
    lin = rng.normal(size=(N, 2))
    lin = (lin - lin.mean(0)) / lin.std(0)
    opt = {}
    if lag:
        dose = rng.gamma(1.0, 1.0, size=(N, T)) * (rng.random((N, T)) < 0.6)
        lags = np.stack([np.concatenate([np.zeros((N, l)), dose[:, :T - l]], axis=1)
                         for l in range(3)], axis=-1)
        opt["lag_design"] = jnp.asarray(lags.reshape(N * T, 3), jnp.float32)
    if extra:
        opt["extra_design"] = jnp.asarray(rng.normal(size=(N * T, 2)), jnp.float32)
    st = init_longbet(
        X_unified=jnp.asarray(np.repeat(X_unit, T, axis=1)), y=y,
        unit_idx=jnp.repeat(jnp.arange(N), T), time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.asarray(np.zeros_like(exposure) if lag else exposure),
        z_vec=jnp.asarray(np.ones_like(exposure, dtype=np.float32) if lag
                          else (exposure > 0).astype(np.float32)),
        obs_mask=jnp.ones(N * T, bool), max_split_mu=jnp.full(P, 2, jnp.uint8),
        max_split_nu=jnp.full(P, 2, jnp.uint8), config=cfg,
        lin_design=jnp.asarray(np.repeat(lin, T, axis=0), jnp.float32),
        num_chains=num_chains, chain_key=jax.random.key(1), **opt,
    )
    return st, y


@pytest.mark.parametrize("lag", [False, True], ids=["clock", "lag_extra"])
@pytest.mark.parametrize("anchor", [False, True], ids=["unanchored", "anchored"])
@pytest.mark.parametrize("ar1", [False, True], ids=["iid", "ar1"])
def test_overdispersed_start_keeps_resid_consistent(lag, anchor, ar1):
    """Regression: with ``anchor_moderator=True`` the dispersed starting
    ``beta`` (or ``kappa``) entered the initial fit through ``nu = 1`` and must
    be taken out of the residual so every chain analyses ``y`` without a shift."""
    st, y = _build(anchor, ar1, num_chains=3, lag=lag, extra=lag)
    key = jax.random.key(2)
    for sweep in range(4):
        err = float(jnp.max(jnp.abs(st.resid - (y[None, :] - _fitted(st)))))
        assert err < 5e-6, f"sweep {sweep}: max |resid - (y - fitted)| = {err:.3e}"
        key, k = jax.random.split(key)
        st = longbet_step(k, st)


def test_anchored_lag_overdisperse_warms_conjugate_block():
    """Engine Finding B regression: with ``use_lag_kernel=True``,
    ``anchor_moderator=True`` and ``use_trend_block=True``, ``_overdisperse_chains``
    draws ``(trend_coef, extra_coef, gamma, kappa)`` jointly at ``nu = 1,
    delta = 0`` before Sweep 1 so Step 3b (``_linear_moderator_step``) does not
    regress an unadjusted residual on a near-zero or negative initial ``w``."""
    st, y = _build(anchor=True, ar1=True, num_chains=4, lag=True, extra=True)
    assert jnp.all(jnp.abs(st.extra_coef) > 0.0)
    err0 = float(jnp.max(jnp.abs(st.resid - (y[None, :] - _fitted(st)))))
    assert err0 < 5e-6
    st1 = longbet_step(jax.random.key(99), st)
    assert float(jnp.max(st1.tau_delta2)) < 5.0

