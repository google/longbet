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

"""Chain-axis safety.

The plan asked for a test that runs one sweep with ``num_chains=1`` and one with
``num_chains=3`` from the same key-per-chain and asserts the first chain
matches.  That exact form is not reachable through ``bartz``'s public API:
``bartz.mcmcstep.step`` splits the key internally, so chain 0 of a multi-chain
run is not the single-chain run whatever key is supplied.

What that test was for -- catching reduction-axis, broadcast and scatter
mistakes -- is covered here directly instead:

* the chain partition is derived from bartz's field metadata, and is checked to
  put the right leaves on each side (a wrong answer here silently breaks
  ``evaluate_trace``, which reads ``leaf_unit`` as a scalar);
* a ``k``-chain sweep is compared **value for value** against ``k`` independent
  single-chain sweeps built from the same per-chain states and driven with the
  same per-chain keys, which is the strongest available form of the intended
  assertion;
* the shared data is checked not to acquire a chain axis, which is what made
  the previous implementation hold ``k`` copies of the design matrix.
"""

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from jax import random

from longbet._config import LongBetConfig
from longbet._state import (
    broadcast_to_chains,
    chain_filter_spec,
    chained_field_names,
    init_longbet,
    split_chain_fields,
)
from longbet._step import longbet_single_step, longbet_step

SHARED_FIELDS = ("X", "y", "obs_mask", "unit_idx", "time_idx", "exposure_idx",
                 "z_vec", "unit_counts", "K_chol", "prec_scale", "resid_unit")


def _state(num_chains=None, **cfg):
    N, T, P = 10, 5, 3
    M = N * T
    rng = np.random.default_rng(0)
    exposure = np.tile(np.clip(np.arange(T) - 1, 0, None), N).astype(np.int32)
    config = LongBetConfig(num_trees_pr=4, num_trees_trt=4,
                           num_chains=num_chains or 1, **cfg)
    X = rng.integers(0, 5, (P, M), dtype=np.uint8)
    y = rng.normal(size=M)
    if config.outcome == "binary":
        y = (y > 0).astype(np.float32)
    return init_longbet(
        X_unified=jnp.asarray(X),
        y=jnp.asarray(y, jnp.float32),
        unit_idx=jnp.repeat(jnp.arange(N), T),
        time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.asarray(exposure),
        z_vec=jnp.asarray(np.tile((np.arange(T) >= 2).astype(np.float32), N)),
        obs_mask=jnp.ones(M, bool),
        max_split_mu=jnp.full(P, 4, jnp.uint8),
        max_split_nu=jnp.full(P, 4, jnp.uint8),
        config=config,
        num_chains=num_chains,
    )


def test_partition_follows_bartz_metadata():
    state = _state(num_chains=3)
    chained = chained_field_names(state)

    # Sampled quantities are per chain.
    for name in ("resid", "beta", "gamma", "sigma2", "mu_fit", "nu_fit",
                 "forest", "forest_nu", "prec_scale_nu", "_chain_anchor"):
        assert name in chained, f"{name} should carry a chain axis"

    # Data and fixed quantities are not.
    for name in SHARED_FIELDS:
        assert name not in chained, f"{name} must stay shared across chains"

    # bartz reads leaf_unit as a scalar even in its own multi-chain layout.
    assert state.forest.leaf_unit.ndim == 0
    assert state.forest.leaf_tree.shape[0] == 3


def test_shared_data_is_not_replicated_per_chain():
    single = _state()
    many = _state(num_chains=4)
    for name in SHARED_FIELDS:
        a, b = getattr(single, name), getattr(many, name)
        if a is None:
            continue
        assert a.shape == b.shape, f"{name} gained a chain axis: {a.shape} -> {b.shape}"
        assert jnp.array_equal(a, b)

    _, shared = split_chain_fields(many)
    assert shared.X.shape == single.X.shape


def test_multichain_sweep_equals_independent_single_chain_sweeps():
    """A k-chain sweep must equal k single-chain sweeps, value for value."""
    num_chains = 3
    base = _state()
    spec = chain_filter_spec(base)

    # Give each chain a genuinely different starting point, so an axis mistake
    # cannot hide behind identical chains.
    rng = np.random.default_rng(3)
    singles = []
    for c in range(num_chains):
        s = eqx.tree_at(
            lambda st: (st.beta, st.gamma, st.sigma2, st.b0, st.b1),
            base,
            (
                jnp.asarray(rng.uniform(0.5, 1.5, base.beta.shape), jnp.float32),
                jnp.asarray(rng.normal(0, 0.3, base.gamma.shape), jnp.float32),
                jnp.float32(0.5 + 0.4 * c),
                jnp.float32(-0.5 + 0.1 * c),
                jnp.float32(0.5 + 0.1 * c),
            ),
        )
        singles.append(s)

    # Stack them into one chained state.
    pers = [eqx.filter(s, spec) for s in singles]
    stacked = jax.tree.map(lambda *xs: jnp.stack(xs), *pers)
    _, shared = split_chain_fields(broadcast_to_chains(base, num_chains))
    chained = eqx.combine(stacked, shared)

    key = random.key(20260906)
    out_chained = longbet_step(key, chained)
    keys = random.split(key, num_chains)  # exactly what longbet_step uses
    outs_single = [longbet_single_step(keys[c], singles[c]) for c in range(num_chains)]

    for c, single_out in enumerate(outs_single):
        for name in ("beta", "gamma", "alpha", "b0", "b1", "sigma2",
                     "sigma_gamma2", "mu_fit", "nu_fit", "resid"):
            got = np.asarray(getattr(out_chained, name))[c]
            want = np.asarray(getattr(single_out, name))
            np.testing.assert_allclose(
                got, want, rtol=5e-4, atol=1e-4,
                err_msg=f"chain {c} disagrees with the single-chain sweep on {name!r}",
            )
        np.testing.assert_allclose(
            np.asarray(out_chained.forest_nu.leaf_tree)[c],
            np.asarray(single_out.forest_nu.leaf_tree),
            rtol=5e-4, atol=1e-4,
        )


@pytest.mark.parametrize("num_chains", [None, 2])
def test_sweep_is_finite_and_keeps_shapes(num_chains):
    state = _state(num_chains=num_chains)
    step = jax.jit(longbet_step)
    key = random.key(0)
    for _ in range(3):
        key, sub = random.split(key)
        state = step(sub, state)
    assert jnp.all(jnp.isfinite(state.beta))
    assert jnp.all(jnp.isfinite(state.gamma))
    assert jnp.all(jnp.isfinite(state.resid))
    if num_chains:
        assert state.beta.shape[0] == num_chains
        assert state.X.ndim == 2


def test_chains_start_overdispersed():
    """R-hat only detects a stuck chain if the chains could start apart.

    Chains separated solely by their random streams can agree with each other
    while all sitting in the same unvisited region, so the starting points are
    drawn from the priors instead.
    """
    common = _state(num_chains=None)
    identical = _state(num_chains=4)                      # no key: replicated
    dispersed = init_longbet(
        X_unified=common.X,
        y=common.y,
        unit_idx=common.unit_idx,
        time_idx=common.time_idx,
        exposure_idx=common.exposure_idx,
        z_vec=common.z_vec,
        obs_mask=common.obs_mask,
        max_split_mu=jnp.full(common.X.shape[0], 4, jnp.uint8),
        max_split_nu=jnp.full(common.X.shape[0], 4, jnp.uint8),
        config=LongBetConfig(num_trees_pr=4, num_trees_trt=4, num_chains=4),
        num_chains=4,
        chain_key=random.key(0),
    )

    for name in ("beta", "gamma", "sigma2", "sigma_gamma2"):
        same = np.asarray(getattr(identical, name))
        apart = np.asarray(getattr(dispersed, name))
        assert np.allclose(same[0], same[1]), f"{name} should be replicated without a key"
        assert not np.allclose(apart[0], apart[1]), f"{name} did not disperse"
    assert np.all(np.asarray(dispersed.sigma2) > 0)
    assert np.all(np.asarray(dispersed.sigma_gamma2) > 0)


def test_dispersed_start_keeps_the_residual_invariant():
    """Dispersing gamma changes the fitted values, so the residual must follow.

    At initialisation both forests are empty, so the fit is exactly gamma_i.
    """
    state = _state(num_chains=None)
    dispersed = init_longbet(
        X_unified=state.X,
        y=state.y,
        unit_idx=state.unit_idx,
        time_idx=state.time_idx,
        exposure_idx=state.exposure_idx,
        z_vec=state.z_vec,
        obs_mask=state.obs_mask,
        max_split_mu=jnp.full(state.X.shape[0], 4, jnp.uint8),
        max_split_nu=jnp.full(state.X.shape[0], 4, jnp.uint8),
        config=LongBetConfig(num_trees_pr=4, num_trees_trt=4, num_chains=3),
        num_chains=3,
        chain_key=random.key(1),
    )
    fitted = np.asarray(dispersed.gamma)[:, np.asarray(dispersed.unit_idx)]
    implied = np.asarray(dispersed.y)[None, :] - fitted
    np.testing.assert_allclose(implied, np.asarray(dispersed.resid), atol=1e-5)

    # And it survives a sweep.
    stepped = longbet_step(random.key(2), dispersed)
    assert jnp.all(jnp.isfinite(stepped.beta))
    assert jnp.all(jnp.isfinite(stepped.resid))


def test_gp_starts_cover_the_declared_prior():
    """Positive clipping can hide chain disagreement while still varying starts.

    The declared prior is ``beta ~ N(0, K)`` read on the identified quotient:
    the sweep maps every trajectory draw onto the canonical half-space
    ``sum_{s>=1} beta_s >= 0``, so the whitened starts are ``sign(v'u) u`` for
    ``u ~ N(0, I)`` and ``v = K_chol' c``, ``c = (0, 1, ..., 1)``. That law has
    closed-form first and second moments and a closed-form sign-crossing rate,
    all of which a start confined to the positive orthant fails.
    """
    base = _state(max_depth_pr=2, max_depth_trt=2)
    dispersed = broadcast_to_chains(base, 1024, key=random.key(781))
    beta = np.asarray(dispersed.beta, dtype=np.float64)
    K_chol = np.asarray(base.K_chol, dtype=np.float64)
    white = np.linalg.solve(K_chol, beta.T).T

    # Every start already satisfies the invariant the sweep maintains.
    assert np.all(beta[:, 1:].sum(axis=1) >= 0)

    c = np.ones(beta.shape[1])
    c[0] = 0.
    v = K_chol.T @ c
    rho = v[0] / np.linalg.norm(v)

    # E[sign(v'u) u] = sqrt(2/pi) v/||v||, and sign**2 == 1 leaves the second
    # moment at the identity, so the quotient costs no dispersion.
    np.testing.assert_allclose(white.mean(0),
                               np.sqrt(2/np.pi) * v / np.linalg.norm(v), atol=.12)
    np.testing.assert_allclose(white.T @ white / len(white),
                               np.eye(white.shape[1]), atol=.16)

    # beta_0 is outside the constrained sum, so it still crosses zero at the
    # rate P(sign(v'u) u_0 < 0) = 1/2 - arcsin(rho)/pi.
    np.testing.assert_allclose(np.mean(beta[:, 0] < 0),
                               .5 - np.arcsin(rho)/np.pi, atol=.05)


@pytest.mark.parametrize("outcome", ["continuous", "binary"])
def test_fixed_coding_preserves_both_arms(outcome):
    """The coding scales are fixed at b0 = 0, b1 = 1 and never move."""
    state = _state(max_depth_pr=2, max_depth_trt=2, outcome=outcome)
    for key in random.split(random.key(691), 3):
        assert float(state.b0) == 0
        assert float(state.b1) == 1
        state = longbet_single_step(key, state)
    assert float(state.b0) == 0 and float(state.b1) == 1
