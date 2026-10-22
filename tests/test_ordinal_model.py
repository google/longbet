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

"""Ordinal state, sweep, chain alignment, and legacy binary compatibility."""
import dataclasses

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, LongBetMulti
from longbet._loop import run_longbet_mcmc
from longbet._ordinal import prepare_ordinal
from longbet._state import init_longbet, split_chain_fields, chained_field_names
from longbet._step import longbet_single_step, longbet_step


def ordinal_config(K=4, **kwargs):
    settings = dict(outcome="ordinal", num_categories=K, num_chains=1,
        num_trees_pr=2, num_trees_trt=2, max_depth_pr=3, max_depth_trt=3,
        min_points_per_leaf_pr=2, min_points_per_leaf_trt=2,
        num_burnin=3, num_sweeps=5, n_skip=2, device="cpu")
    settings.update(kwargs)
    return LongBetConfig(**settings)


def ordinal_panel(K=4):
    rng = np.random.default_rng(223)
    N, T = 12, 4
    x = rng.normal(size=(N, 2))
    y = (np.arange(N*T) % K).reshape(N, T).astype(float)
    y[0, 0] = y[2, 3] = np.nan
    z = np.zeros((N, T))
    z[:6, 2:] = 1
    return dict(y=y, x=x, z=z, t=np.arange(1, T+1))


def ordinal_state(K=4, chains=1, key=None, placeholder=0.):
    data = ordinal_panel(K)
    y = data["y"].ravel()
    prepared = prepare_ordinal(y, K)
    N, T = data["y"].shape
    return init_longbet(X_unified=jnp.array(np.tile(np.arange(T), (2, N)), jnp.uint8),
        y=jnp.array(np.where(prepared.obs_mask, y, placeholder), jnp.float32),
        obs_mask=jnp.array(prepared.obs_mask), unit_idx=jnp.repeat(jnp.arange(N), T),
        time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.array(np.tile([0, 0, 1, 2], N)),
        z_vec=jnp.array(data["z"].ravel(), jnp.float32),
        max_split_mu=jnp.full(2, 3, jnp.uint8), max_split_nu=jnp.full(2, 3, jnp.uint8),
        config=ordinal_config(K), offset=prepared.offset,
        num_chains=chains, chain_key=key)


def assert_invariants(state, check_latents=True):
    mean = (np.asarray(state.alpha)[..., None]*np.asarray(state.mu_fit)
            + np.where(np.asarray(state.z_vec) == 1,
                       np.asarray(state.b1)[..., None], np.asarray(state.b0)[..., None])
            * np.asarray(state.beta)[..., np.asarray(state.exposure_idx)] * np.asarray(state.nu_fit)
            + np.asarray(state.gamma)[..., np.asarray(state.unit_idx)])
    if getattr(state, "use_trend_block", False):
        mean = mean + np.asarray(state.trend_coef) @ np.asarray(state.trend_design).T
    mask, labels = np.asarray(state.obs_mask), np.asarray(state.y, int)
    # The residual is maintained incrementally in float32; treated cells take an
    # add/subtract cycle of the treatment term in every sweep, so the drift from
    # the freshly recomputed mean reaches a few 1e-5 (2.5e-5 with jaxlib 0.11.2 on
    # the 2-chain K=2 case). A bookkeeping error would be of the order of the term
    # itself, i.e. >= 1e-2, which 1e-4 still catches.
    np.testing.assert_allclose(state.resid, np.where(mask, np.asarray(state.z)-mean, 0), atol=1e-4)
    assert (np.asarray(state.sigma2) == 1).all()
    assert (np.asarray(state.error_cov_inv.value) == 1).all()
    assert (np.asarray(state.resid)[..., ~mask] == 0).all()
    if check_latents:
        # K=2 deliberately dispatches the legacy binary kernel, whose unused
        # missing latents are random. Only the general ordinal path zeros them.
        if state.num_categories > 2:
            assert (np.asarray(state.z)[..., ~mask] == 0).all()
        cp = np.asarray(state.cutpoints)
        assert np.isfinite(cp).all()
        full = np.concatenate([np.full((*cp.shape[:-1], 1), -np.inf),
            np.zeros((*cp.shape[:-1], 1)), cp, np.full((*cp.shape[:-1], 1), np.inf)], axis=-1)
        assert (np.diff(full, axis=-1) > 0).all()
        assert (np.asarray(state.z)[..., mask] > full[..., labels[mask]]).all()
        assert (np.asarray(state.z)[..., mask] < full[..., labels[mask]+1]).all()


@pytest.mark.parametrize("K", [3, 5])
@pytest.mark.parametrize("chains", [1, 3])
def test_initialization_and_complete_sweep_invariants(K, chains):
    state = ordinal_state(K, chains, jax.random.key(202))
    assert_invariants(state)
    assert state.binary_indices is None
    assert state.num_categories == K
    assert state.X.ndim == 2 and state.y.ndim == 1 and state.obs_mask.ndim == 1
    assert "cutpoints" in chained_field_names(state)
    if chains > 1:
        assert np.unique(np.asarray(state.cutpoints), axis=0).shape[0] == chains
    for key in jax.random.split(jax.random.key(4), 4):
        state = longbet_step(key, state)
        assert_invariants(state)


def test_missing_placeholders_cannot_change_sweeps():
    a = ordinal_state(5, placeholder=0.)
    b = ordinal_state(5, placeholder=np.nan)
    c = ordinal_state(5, placeholder=999.)
    outputs = [longbet_single_step(jax.random.key(3), s) for s in [a, b, c]]
    for other in outputs[1:]:
        for name in ("z", "cutpoints", "resid", "beta", "gamma", "mu_fit", "nu_fit"):
            np.testing.assert_array_equal(getattr(outputs[0], name), getattr(other, name))


def test_ordinal_vmap_matches_independent_chain_steps():
    many = ordinal_state(4, 3, jax.random.key(222))
    per, shared = split_chain_fields(many)
    key = jax.random.key(787)
    together = longbet_step(key, many)
    for c, ck in enumerate(jax.random.split(key, 3)):
        single = eqx.combine(jax.tree.map(lambda x: x[c], per), shared)
        out = longbet_single_step(ck, single)
        for name in ("z", "cutpoints", "beta", "gamma", "b0", "b1", "alpha",
                     "sigma2", "resid", "mu_fit", "nu_fit"):
            # Batched and scalar float32 reductions in the forest kernels do not
            # agree bit for bit, and the GP solve amplifies the difference, so
            # this matches the tolerance of the sibling sweep-equivalence test
            # in test_chain_axis. A structural mismatch is order one, not 1e-3.
            np.testing.assert_allclose(np.asarray(getattr(together, name))[c], getattr(out, name),
                                       atol=2e-4, rtol=1e-3, err_msg=name)


@pytest.mark.parametrize("chains", [1, 2])
@pytest.mark.parametrize("label", [None, 0., 1.])
def test_binary_and_two_category_ordinal_are_identical(chains, label):
    data = ordinal_panel(2)
    if label is not None:
        data["y"][np.isfinite(data["y"])] = label
    cfg = ordinal_config(2, num_chains=chains)
    ordinal = LongBet(cfg).fit(**data)
    binary = LongBet(dataclasses.replace(cfg, outcome="binary", num_categories=None)).fit(**data)
    assert ordinal.offset_ == binary.offset_
    assert ordinal.meany == binary.meany == 0 and ordinal.sdy == binary.sdy == 1
    assert ordinal.trace.cutpoints.shape == ((chains, 5, 0) if chains > 1 else (5, 0))
    assert binary.trace.cutpoints is None
    # Only static outcome metadata and the optional empty threshold trace differ.
    for name in ("beta", "gamma", "b0", "b1", "alpha", "sigma2", "sigma_gamma2",
                 "mu_trace", "nu_trace"):
        for a, b in zip(jax.tree.leaves(getattr(ordinal.trace, name)),
                        jax.tree.leaves(getattr(binary.trace, name))):
            np.testing.assert_array_equal(a, b, err_msg=name)
    np.testing.assert_array_equal(ordinal.state.z, binary.state.z)
    assert_invariants(ordinal.state)


def test_cutpoint_trace_matches_saved_sweeps_after_thinning():
    initial = ordinal_state(3, 2, jax.random.key(35))
    key = jax.random.key(457)
    result = run_longbet_mcmc(key, initial, n_burn=2, n_save=3, n_skip=2,
                             inner_loop_length=3)
    state = initial
    kept_cp, kept_beta = [], []
    for i in range(8):
        step_key, key = jax.random.split(key)
        state = longbet_step(step_key, state)
        if i >= 2 and (i-2+1) % 2 == 0:
            kept_cp.append(np.asarray(state.cutpoints))
            kept_beta.append(np.asarray(state.beta))
    np.testing.assert_allclose(result.main_trace.cutpoints, np.stack(kept_cp, axis=1), atol=1e-6)
    np.testing.assert_allclose(result.main_trace.beta, np.stack(kept_beta, axis=1), atol=1e-6)
    assert (np.asarray(result.main_trace.sigma2) == 1).all()
    assert (np.asarray(result.main_trace.mu_trace.error_cov_inv) == 1).all()


@pytest.mark.parametrize("y", [["0", "1"], [0., np.inf], [np.nan, np.nan], [0, 1.5]])
def test_public_input_rejection(y):
    with pytest.raises(ValueError):
        LongBet(ordinal_config(3)).fit(y=[y], x=[[1]], z=[[0, 0]])


def test_public_fit_with_only_the_top_category_observed():
    data = ordinal_panel(5)
    data["y"][np.isfinite(data["y"])] = 4.
    model = LongBet(ordinal_config(5, num_chains=2, num_sweeps=3)).fit(**data)
    assert model.meany == 0 and model.sdy == 1
    assert model.offset_ > 0
    assert_invariants(model.state)
    assert model.trace.cutpoints.shape == (2,3,3)
