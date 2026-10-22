"""CHANGE move: tree validity, cache and residual consistency, and target invariance."""

import dataclasses

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from bartz.grove import check_trace, evaluate_forest
from bartz.mcmcstep import State, Wishart
from bartz.mcmcstep._step import apply_moves_to_leaf_indices

from longbet import LongBetConfig
from longbet._change_move import change_step
from longbet._forest_cache import current_forest_fit
from longbet._state import init_longbet
from longbet._step import _load, longbet_single_step


def _single_chain_state(seed=0, N=30, T=6, P=3, **cfg):
    rng = np.random.default_rng(seed)
    M = N * T
    exposure = np.tile(np.concatenate([np.zeros(2), np.arange(1, T - 1)]), N).astype(np.int32)
    z_vec = (exposure > 0) & np.repeat(np.arange(N) < N // 2, T)
    exposure = np.where(z_vec, exposure, 0).astype(np.int32)
    options = dict(num_trees_pr=4, num_trees_trt=4, min_points_per_leaf_pr=3,
                   min_points_per_leaf_trt=3,
                   sigma_prior_a=3.0, sigma_prior_b=1.0)
    options.update(cfg)
    config = LongBetConfig(**options)
    X = rng.integers(0, 21, (P, M)).astype(np.uint8)
    y = (2.0 * (X[0] > 10) + 0.5 * z_vec * exposure + rng.normal(0, 0.5, M)).astype(np.float32)
    state = init_longbet(
        X_unified=jnp.asarray(X), y=jnp.asarray(y),
        unit_idx=jnp.repeat(jnp.arange(N), T), time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.asarray(exposure), z_vec=jnp.asarray(z_vec.astype(np.float32)),
        obs_mask=jnp.ones(M, bool), max_split_mu=jnp.full(P, 20, jnp.uint8),
        max_split_nu=jnp.full(P, 20, jnp.uint8), config=config,
    )
    return state, config


def _check_forest(forest, X, offset_expected=None):
    """Leaf memberships agree with the rules, caches with memberships, trees are valid."""
    ids = np.asarray(apply_moves_to_leaf_indices(forest.leaf_indices, forest.to_prune, forest.move_node))
    from_rules = np.asarray(evaluate_forest(X, forest, sum_batch_axis=0)) + float(forest.offset)
    from_members = np.asarray(current_forest_fit(forest))
    np.testing.assert_allclose(from_members, from_rules, rtol=1e-5, atol=1e-5)
    codes = np.asarray(check_trace(forest, forest.max_split))
    assert np.all(codes == 0), codes
    if forest.count_tree is not None:
        tree_size = forest.leaf_tree.shape[-1]
        for j in range(ids.shape[0]):
            counts = np.bincount(ids[j], minlength=tree_size)
            leaves = np.unique(ids[j])
            np.testing.assert_array_equal(np.asarray(forest.count_tree)[j, leaves], counts[leaves])


def test_sweeps_keep_forests_consistent():
    state, _ = _single_chain_state()
    key = jax.random.key(3)
    for _ in range(25):
        key, sub = jax.random.split(key)
        state = longbet_single_step(sub, state)
        _check_forest(state.forest, state.X)
        _check_forest(state.forest_nu, state.X)


def _mu_view(state):
    """The prognostic-forest view exactly as ``longbet_single_step`` builds it."""
    return State(
        _chain_anchor=state._chain_anchor, X=state.X, y=state.y, z=None, binary_indices=None,
        resid=_load(state.resid / state.alpha, state.resid_unit), resid_unit=state.resid_unit,
        resid_eff_scale=state.resid_eff_scale, resid_inexact_integral=state.resid_inexact_integral,
        error_cov_inv=Wishart(nu=None, rate=None, value=jnp.square(state.alpha) / state.sigma2),
        error_scale=None, prec_scale=state.prec_scale, inv_sdev_scale=state.inv_sdev_scale,
        inv_sdev_unit=state.inv_sdev_unit, n_non_missing=state.n_non_missing,
        sum_diag_prec_scale=state.sum_diag_prec_scale, forest=state.forest, config=state.config,
    )


def test_change_step_bookkeeping_and_acceptance():
    state, _ = _single_chain_state(seed=5)
    key = jax.random.key(11)
    for _ in range(8):  # grow some trees first
        key, sub = jax.random.split(key)
        state = longbet_single_step(sub, state)
    view = _mu_view(state)
    step = jax.jit(change_step)
    rules_changed = 0
    for _ in range(60):
        key, sub = jax.random.split(key)
        before_fit = np.asarray(current_forest_fit(view.forest))
        before_resid = np.asarray(view.resid, dtype=np.float64) * float(view.resid_unit)
        rules_before = (np.asarray(view.forest.var_tree).copy(), np.asarray(view.forest.split_tree).copy())
        new_view = step(sub, view)
        after_fit = np.asarray(current_forest_fit(new_view.forest))
        after_resid = np.asarray(new_view.resid, dtype=np.float64) * float(new_view.resid_unit)
        # resid + fit is invariant: the residual absorbs exactly the change in the forest
        np.testing.assert_allclose(after_resid + after_fit, before_resid + before_fit, rtol=1e-4, atol=1e-4)
        assert not np.any(np.asarray(new_view.forest.to_prune))
        _check_forest(new_view.forest, view.X)
        rules_changed += int(np.any(np.asarray(new_view.forest.var_tree) != rules_before[0])
                             or np.any(np.asarray(new_view.forest.split_tree) != rules_before[1]))
        view = new_view
    assert rules_changed > 0, "no CHANGE proposal was ever accepted; the test is vacuous"
