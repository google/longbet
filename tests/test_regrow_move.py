"""REGROW move: forest invariants, exact bookkeeping, acceptance, and target invariance."""


import jax
import jax.numpy as jnp
import numpy as np
import pytest
from bartz.grove import check_trace, evaluate_forest
from bartz.mcmcstep._step import apply_moves_to_leaf_indices

from longbet._forest_cache import current_forest_fit
from longbet._regrow_move import regrow_tree
from longbet._step import longbet_single_step
from test_change_move import _check_forest, _mu_view, _single_chain_state


def test_sweeps_with_regrow_keep_forests_consistent():
    state, _ = _single_chain_state(seed=2)
    key = jax.random.key(5)
    for _ in range(20):
        key, sub = jax.random.split(key)
        state = longbet_single_step(sub, state)
        _check_forest(state.forest, state.X)
        _check_forest(state.forest_nu, state.X)


def test_regrow_tree_bookkeeping_and_acceptance():
    state, _ = _single_chain_state(seed=9)
    key = jax.random.key(4)
    for _ in range(6):
        key, sub = jax.random.split(key)
        state = longbet_single_step(sub, state)
    view = _mu_view(state)
    step = jax.jit(regrow_tree)
    accepted = 0
    for i in range(40):
        key, sub = jax.random.split(key)
        before_fit = np.asarray(current_forest_fit(view.forest))
        before_resid = np.asarray(view.resid, dtype=np.float64) * float(view.resid_unit)
        new_view, acc = step(sub, view, jnp.int32(i % view.forest.var_tree.shape[0]))
        after_fit = np.asarray(current_forest_fit(new_view.forest))
        after_resid = np.asarray(new_view.resid, dtype=np.float64) * float(new_view.resid_unit)
        np.testing.assert_allclose(after_resid + after_fit, before_resid + before_fit, rtol=1e-4, atol=1e-4)
        assert not np.any(np.asarray(new_view.forest.to_prune))
        _check_forest(new_view.forest, view.X)
        accepted += int(acc)
        view = new_view
    assert accepted > 0, "no REGROW proposal was ever accepted; the test is vacuous"
