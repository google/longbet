"""Internal CHANGE move: tree validity, cache and residual consistency, prior evaluation, acceptance."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet._change_internal_move import change_internal_sweep, tree_log_prior
from longbet._forest_cache import current_forest_fit
from longbet._step import longbet_single_step
try:
    from tests.test_change_move import _check_forest, _mu_view, _single_chain_state
except ImportError:
    from test_change_move import _check_forest, _mu_view, _single_chain_state


def test_sweeps_keep_forests_consistent():
    state, _ = _single_chain_state(seed=3)
    key = jax.random.key(5)
    for i in range(12):
        state = longbet_single_step(jax.random.fold_in(key, i), state)
        _check_forest(state.forest, state.X)
        _check_forest(state.forest_nu, state.X)


def test_tree_log_prior_matches_direct_count():
    """A root split with two leaves: log p_nt[1] - log p - log(range) + 2 log(1 - p_nt[2])."""
    state, _ = _single_chain_state(seed=1)
    forest = state.forest
    half = forest.var_tree.shape[1]
    var_tree = jnp.zeros(half, forest.var_tree.dtype).at[1].set(0)
    split_tree = jnp.zeros(half, forest.split_tree.dtype).at[1].set(7)
    args = (forest.max_split, forest.blocked_vars, forest.p_nonterminal, half)
    lp = float(tree_log_prior(var_tree, split_tree, *args))
    p = int(forest.max_split.shape[0])
    p_nt = np.asarray(forest.p_nonterminal)
    expected = np.log(p_nt[1]) - np.log(p) - np.log(float(forest.max_split[0])) + 2 * np.log1p(-p_nt[2])
    assert np.isclose(lp, expected, rtol=1e-5), (lp, expected)
    # A descendant rule on the same variable outside its restricted range has zero prior.
    bad_var = var_tree.at[2].set(0)
    bad_split = split_tree.at[2].set(9)          # the left child of "x0 < 7" cannot split x0 at 9
    assert float(tree_log_prior(bad_var, bad_split, *args)) == -np.inf
    good_split = split_tree.at[2].set(3)
    lp2 = float(tree_log_prior(bad_var, good_split, *args))
    expected2 = (np.log(p_nt[1]) - np.log(p) - np.log(float(forest.max_split[0]))
                 + np.log(p_nt[2]) - np.log(p) - np.log(6.0)      # cutpoints 1..6 are admissible under "< 7"
                 + np.log1p(-p_nt[3]) + 2 * np.log1p(-p_nt[4]))
    assert np.isclose(lp2, expected2, rtol=1e-5), (lp2, expected2)


def _interaction_state(seed=7, N=80, T=6, P=3, **cfg):
    """A panel whose prognostic surface is a two-variable interaction, so trees need internal rules."""
    from longbet import LongBetConfig
    from longbet._state import init_longbet
    rng = np.random.default_rng(seed)
    M = N * T
    exposure = np.tile(np.concatenate([np.zeros(2), np.arange(1, T - 1)]), N).astype(np.int32)
    z_vec = (exposure > 0) & np.repeat(np.arange(N) < N // 2, T)
    exposure = np.where(z_vec, exposure, 0).astype(np.int32)
    options = dict(num_trees_pr=2, num_trees_trt=2, min_points_per_leaf_pr=3, min_points_per_leaf_trt=3,
                   sigma_prior_a=3.0, sigma_prior_b=1.0,
                   use_inter_ensemble_move=False, use_trend_block=False)
    options.update(cfg)
    config = LongBetConfig(**options)
    X = rng.integers(0, 21, (P, M)).astype(np.uint8)
    # A moderate signal-to-noise ratio: strong enough to grow interaction trees, weak enough
    # that shifting an internal cutpoint by one bin is not a hopeless proposal.
    y = (1.0 * (X[0] > 10) * (X[1] > 10) - 0.7 * (X[0] <= 10) * (X[1] > 14) + rng.normal(0, 1.0, M)).astype(np.float32)
    state = init_longbet(
        X_unified=jnp.asarray(X), y=jnp.asarray(y),
        unit_idx=jnp.repeat(jnp.arange(N), T), time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.asarray(exposure), z_vec=jnp.asarray(z_vec.astype(np.float32)),
        obs_mask=jnp.ones(M, bool), max_split_mu=jnp.full(P, 20, jnp.uint8),
        max_split_nu=jnp.full(P, 20, jnp.uint8), config=config,
    )
    return state, config


def test_internal_change_bookkeeping_and_acceptance():
    state, _ = _interaction_state()
    key = jax.random.key(11)
    for _ in range(40):  # grow trees with internal structure first
        key, sub = jax.random.split(key)
        state = longbet_single_step(sub, state)
    view = _mu_view(state)
    split0 = np.asarray(view.forest.split_tree); half = split0.shape[1]
    deep = [(j, k) for j in range(split0.shape[0]) for k in range(1, half // 2)
            if split0[j, k] > 0 and (split0[j, 2 * k] > 0 or split0[j, 2 * k + 1] > 0)]
    assert deep, "the interaction panel did not produce a tree with nonterminal children"
    step = jax.jit(change_internal_sweep)
    internal_changed = 0
    for _ in range(150):
        key, sub = jax.random.split(key)
        before_fit = np.asarray(current_forest_fit(view.forest))
        before_resid = np.asarray(view.resid, dtype=np.float64) * float(view.resid_unit)
        var_b = np.asarray(view.forest.var_tree).copy(); split_b = np.asarray(view.forest.split_tree).copy()
        new_view, accepted = step(sub, view)
        after_fit = np.asarray(current_forest_fit(new_view.forest))
        after_resid = np.asarray(new_view.resid, dtype=np.float64) * float(new_view.resid_unit)
        np.testing.assert_allclose(after_resid + after_fit, before_resid + before_fit, rtol=1e-4, atol=1e-4)
        assert not np.any(np.asarray(new_view.forest.to_prune))
        _check_forest(new_view.forest, view.X)
        var_a = np.asarray(new_view.forest.var_tree); split_a = np.asarray(new_view.forest.split_tree)
        # The shape never changes: a node is nonterminal before iff it is after.
        np.testing.assert_array_equal(split_a > 0, split_b > 0)
        changed = (var_a != var_b) | (split_a != split_b)
        assert np.all(changed.any(axis=1) <= np.asarray(accepted)), "a rule changed in a tree whose proposal was rejected"
        # Count changes at nodes with nonterminal children: the case the leaf-parent CHANGE cannot reach.
        for j in np.flatnonzero(changed.any(axis=1)):
            for node in np.flatnonzero(changed[j]):
                left = 2 * node
                if left < half and (split_a[j, left] > 0 or split_a[j, left + 1] > 0):
                    internal_changed += 1
        view = new_view
    assert internal_changed > 0, "no rule at a node with nonterminal children ever changed; the test is vacuous"


def test_evaluate_change_is_exact_for_noop_and_routes_like_bartz():
    """A no-op proposal changes nothing; re-routed memberships agree with bartz's traversal; an
    inadmissible descendant rule gets zero prior."""
    from bartz import grove
    from bartz.mcmcstep._step import apply_moves_to_leaf_indices
    from longbet._change_internal_move import _forest_consts, evaluate_change

    state, _ = _interaction_state(seed=3)
    key = jax.random.key(2)
    for _ in range(40):
        key, sub = jax.random.split(key)
        state = longbet_single_step(sub, state)
    view = _mu_view(state)
    f = view.forest
    c = _forest_consts(view)
    idx_all = apply_moves_to_leaf_indices(f.leaf_indices, f.to_prune, f.move_node)
    split = np.asarray(f.split_tree); var = np.asarray(f.var_tree); half = split.shape[1]
    internal = [(j, k) for j in range(split.shape[0]) for k in range(1, half) if split[j, k] > 0]
    assert internal
    checked_inadmissible = False
    for j, k in internal:
        cases = [("noop", var[j, k], split[j, k]), ("shift", var[j, k], max(1, split[j, k] - 1)),
                 ("newvar", (var[j, k] + 1) % var.shape[0] if False else (int(var[j, k]) + 1) % int(f.max_split.shape[0]), 10)]
        for label, vn, sn in cases:
            ev = evaluate_change(c, f.leaf_tree[j], idx_all[j], view.resid, f.var_tree[j], f.split_tree[j],
                                 jnp.int32(k), jnp.int32(vn), jnp.int32(sn))
            ref = grove.traverse_forest(view.X, ev["var_tree_new"][None], ev["split_tree_new"][None])[0]
            np.testing.assert_array_equal(np.asarray(ev["new_idx"]), np.asarray(ref))
            if label == "noop":
                np.testing.assert_array_equal(np.asarray(ev["new_idx"]), np.asarray(idx_all[j]))
                assert float(ev["lp_new"]) == float(ev["lp_old"]) and np.isfinite(float(ev["lp_old"]))
                assert float(ev["ll_new"]) == float(ev["ll_old"])
                assert bool(ev["min_leaf_ok"])
            else:
                # The prior of the new tree is finite exactly when every descendant rule is still
                # admissible under the new ancestor restriction.
                v_new = np.asarray(ev["var_tree_new"]); s_new = np.asarray(ev["split_tree_new"])
                admissible = True
                for node in range(1, half):
                    if s_new[node] == 0:
                        continue
                    lo, hi = 1, int(f.max_split[v_new[node]]) + 1
                    a = node // 2; child = node
                    while a >= 1:
                        if s_new[a] > 0 and v_new[a] == v_new[node]:
                            if child == 2 * a:      # left child: x < split
                                hi = min(hi, int(s_new[a]))
                            else:
                                lo = max(lo, int(s_new[a]))
                        child, a = a, a // 2
                    if not (lo <= int(s_new[node]) < hi):
                        admissible = False
                assert np.isfinite(float(ev["lp_new"])) == admissible, (label, j, k)
                checked_inadmissible |= not admissible
    assert checked_inadmissible, "no inadmissible descendant case was exercised"


# ---------------------------------------------------------------------------
# Geweke-style check of the tree kernel: with the internal move in the sweep,
# the successive-conditional simulator must leave the tree prior invariant.
# ---------------------------------------------------------------------------

def _prior_panel_state(max_split=8):
    """A tiny panel with a wider cutpoint grid than test_geweke's, so cutpoint marginals are informative."""
    import equinox as eqx  # noqa: F401
    from longbet import LongBetConfig
    from longbet._state import init_longbet
    N, T, P = 8, 4, 2
    M = N * T
    rng = np.random.default_rng(0)
    exposure = np.tile(np.arange(T), N).astype(np.int32)
    config = LongBetConfig(num_trees_pr=2, num_trees_trt=2, random_intercept=True,
                           gamma_prior_a=4.0, gamma_prior_b=2.0, sigma_prior_a=5.0, sigma_prior_b=3.0,
                           min_points_per_leaf_pr=1, min_points_per_leaf_trt=1)
    state = init_longbet(
        X_unified=jnp.asarray(rng.integers(0, max_split + 1, (P, M), dtype=np.uint8)),
        y=jnp.zeros(M, jnp.float32), unit_idx=jnp.repeat(jnp.arange(N), T), time_idx=jnp.tile(jnp.arange(T), N),
        exposure_idx=jnp.asarray(exposure), z_vec=jnp.asarray((exposure > 0).astype(np.float32)),
        obs_mask=jnp.ones(M, bool), max_split_mu=jnp.full(P, max_split, jnp.uint8),
        max_split_nu=jnp.full(P, max_split, jnp.uint8), config=config,
    )
    return state


def _successive_conditional_forest(seed, n_iter, warmup, max_split=8, with_move=True):
    """Forest statistics along the successive-conditional chain, with or without the internal move.

    Count thresholds restrict the tree prior to feasible trees, so closed-form
    node probabilities are not the target; the established sweep without the
    move samples the same restricted prior and is the reference.
    """
    import equinox as eqx
    import longbet._step as step_mod

    original = step_mod.change_internal_step
    if not with_move:
        step_mod.change_internal_step = lambda k, s: s
    try:
        state = _prior_panel_state(max_split)
        body = longbet_single_step.__wrapped__      # fresh trace, so the patch is bound in
        sweep = jax.jit(lambda k, s: body(k, s))
        key = jax.random.key(seed)
        root_split, root_var, n_internal, child_split = [], [], [], []
        for it in range(n_iter):
            key, sub = jax.random.split(key)
            k_data, k_step = jax.random.split(sub)
            b_z = jnp.where(state.z_vec == 1.0, state.b1, state.b0)
            mean = (state.alpha * state.mu_fit + b_z * state.beta[state.exposure_idx] * state.nu_fit
                    + state.gamma[state.unit_idx])
            eps = jax.random.normal(k_data, mean.shape, jnp.float32) * jnp.sqrt(state.sigma2)
            state = eqx.tree_at(lambda s: (s.y, s.resid), state, (mean + eps, eps))
            state = sweep(k_step, state)
            if it < warmup:
                continue
            for forest in (state.forest, state.forest_nu):
                sp = np.asarray(forest.split_tree); va = np.asarray(forest.var_tree)
                root_split.append(sp[:, 1]); root_var.append(va[:, 1]); n_internal.append((sp > 0).sum(axis=1))
                child_split.append(sp[:, 2:4] > 0)
    finally:
        step_mod.change_internal_step = original
    return (np.stack(root_split), np.stack(root_var), np.stack(n_internal), np.stack(child_split))


def _forest_summaries(stats, max_split):
    """Per-tree summaries whose means the with/without runs must share: nonterminal root, root cutpoint and
    variable given a split, nonterminal children, and the number of nonterminal nodes."""
    root_split, root_var, n_internal, child_split = stats
    out = {}
    for j in range(root_split.shape[1]):
        is_split = root_split[:, j] > 0
        out[f"tree {j}: root nonterminal"] = is_split.astype(float)
        out[f"tree {j}: root cutpoint"] = np.where(is_split, root_split[:, j], (max_split + 1) / 2).astype(float)
        out[f"tree {j}: root variable"] = np.where(is_split, root_var[:, j], 0.5).astype(float)
        for c in range(2):
            out[f"tree {j}: child {c} nonterminal"] = child_split[:, j, c].astype(float)
        out[f"tree {j}: nonterminal nodes"] = n_internal[:, j].astype(float)
    return out


@pytest.mark.slow
def test_tree_prior_is_invariant_under_the_sweep():
    """With the internal move in the sweep, the successive-conditional simulator leaves the tree
    prior invariant: every forest summary agrees with the established sweep without the move."""
    from longbet._diagnostics import compute_ess
    max_split = 8
    with_move = _forest_summaries(_successive_conditional_forest(20260914, 8000, 500, max_split, True), max_split)
    without = _forest_summaries(_successive_conditional_forest(20260915, 8000, 500, max_split, False), max_split)

    def mean_mcse(sample):
        sample = np.asarray(sample, dtype=np.float64)
        ess = max(float(compute_ess(sample[None, :])), 2.0)
        return sample.mean(), float(np.std(sample, ddof=1)) / np.sqrt(ess), ess

    failures, report = [], []
    for label in with_move:
        m1, se1, ess1 = mean_mcse(with_move[label]); m2, se2, ess2 = mean_mcse(without[label])
        z = abs(m1 - m2) / max(np.hypot(se1, se2), 1e-12)
        report.append(f"{label}: with {m1:.3f} without {m2:.3f} z {z:.2f} (ESS {ess1:.0f}/{ess2:.0f})")
        if not np.isfinite(z) or z > 4.0:
            failures.append(report[-1])
    print("\n".join(report))
    assert not failures, "tree marginals differ with and without the internal move:\n  " + "\n  ".join(failures)
