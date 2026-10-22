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

"""Opt-in REGROW move: a data-driven whole-tree independence proposal.

Local GROW/PRUNE/CHANGE moves cannot travel between tree topologies that
represent nearly the same function when the likelihood is sharply peaked: every
intermediate tree fits far worse. XBART's grow-from-root recursion builds a
whole tree at once, choosing at every node between stopping and each admissible
split with probability proportional to the tree prior times the integrated
leaf likelihood of that choice. Used on its own that recursion is not a valid
MCMC kernel; used as a Metropolis-Hastings independence proposal it is, because
its density is a product of the categorical probabilities it sampled and can be
evaluated for the current tree under the same residual.

For one tree with the other trees fixed, the acceptance ratio is

    [p(T') L(T') / q(T')] / [p(T) L(T) / q(T)],

with ``p`` the bartz tree prior (nonterminal probabilities by depth and the
ancestor-restricted uniform rule prior), ``L`` the leaf-integrated likelihood,
and ``q`` the proposal density. Count thresholds restrict the proposal's
candidates and are applied identically when evaluating ``q(T)``. Accepted trees
get exact conditional leaf draws and an exact residual update. Univariate
outcomes with optional per-datapoint precision scales, complete or
precision-masked panels.
"""

from __future__ import annotations

from dataclasses import replace

import jax
import jax.numpy as jnp
from jax import lax, random
from jax.scipy.special import logsumexp
from jaxtyping import Array, Key

from bartz import grove
from bartz.mcmcstep import State
from bartz.mcmcstep._moves import split_range
from bartz.mcmcstep._state import scaled_error_cov_inv
from bartz.mcmcstep._step import apply_moves_to_leaf_indices


def _node_candidates(node, member, X, w, t, var_arr, split_arr, max_split, blocked,
                     p_nt, e, lam, min_leaf, num_bins, half):
    """Log-scores of stopping and of every admissible split at one node.

    Returns ``(score_stop, score_split[p, B], log_prior_leaf, log_prior_split[p, B],
    num_allowed_vars, growable)``; masked candidates carry ``-inf``.
    """
    p = X.shape[0]
    mw = member * w
    vals = jnp.stack([member, mw, mw * t], axis=-1)
    hist = jax.vmap(lambda xv: jnp.zeros((num_bins, 3), jnp.float32).at[xv].add(vals))(X)
    counts, sw, st = hist[..., 0], hist[..., 1], hist[..., 2]
    # bins < c go left: left sums for cut c are the cumulative sums up to c - 1
    lc = jnp.concatenate([jnp.zeros((p, 1), jnp.float32), jnp.cumsum(counts, axis=1)], axis=1)[:, :num_bins]
    lw = jnp.concatenate([jnp.zeros((p, 1), jnp.float32), jnp.cumsum(sw, axis=1)], axis=1)[:, :num_bins]
    lt = jnp.concatenate([jnp.zeros((p, 1), jnp.float32), jnp.cumsum(st, axis=1)], axis=1)[:, :num_bins]
    tc, tw, ttot = counts.sum(axis=1, keepdims=True), sw.sum(axis=1, keepdims=True), st.sum(axis=1, keepdims=True)
    rc, rw, rt = tc - lc, tw - lw, ttot - lt
    W_total = jnp.sum(mw)
    T_total = jnp.sum(mw * t)

    def logml(P, r):
        prec_post = lam + e * P
        return 0.5 * jnp.log(lam / prec_post) + 0.5 * e * e * r * r / prec_post

    ranges = jax.vmap(lambda v: split_range(var_arr, split_arr, max_split, node, v))(jnp.arange(p))
    l_v, r_v = ranges
    var_ok = (r_v > l_v)
    if blocked is not None:
        var_ok &= ~jnp.isin(jnp.arange(p), blocked)
    num_allowed = jnp.sum(var_ok)
    cuts = jnp.arange(num_bins)[None, :]
    admissible = var_ok[:, None] & (cuts >= l_v[:, None]) & (cuts < r_v[:, None])
    growable = (num_allowed > 0) & (node < half)
    pnt = jnp.where(growable, p_nt[node], 0.0)
    log_rule_prior = -jnp.log(jnp.maximum(num_allowed, 1).astype(jnp.float32)) - jnp.log(
        jnp.maximum(r_v - l_v, 1).astype(jnp.float32))[:, None]
    # The acceptance ratio flattens this array using (variable, cut) indices.
    # Match the complete score grid; broadcasting only in scores would leave
    # a (variables, 1) prior whose flattened indices select the wrong rule.
    log_prior_split = jnp.broadcast_to(
        jnp.log(pnt) + log_rule_prior, admissible.shape
    )
    log_prior_leaf = jnp.log1p(-pnt)
    feasible = admissible
    if min_leaf is not None:
        feasible &= (lc >= min_leaf) & (rc >= min_leaf)
    score_split = jnp.where(feasible, log_prior_split + logml(lw, lt) + logml(rw, rt), -jnp.inf)
    score_stop = log_prior_leaf + logml(W_total, T_total)
    return score_stop, score_split, log_prior_leaf, log_prior_split, growable


#: Maximum nodes visited by one grow-from-root traversal (a stack-ordered
#: depth-first walk). Trees that would need more are never proposed and never
#: left through this kernel, which keeps the kernel exact.
REGROW_MAX_VISITS = 48


def regrow_tree(key, state: State, tree: int):
    """Propose a regrown tree ``tree`` and accept or reject it. Returns (state, accepted)."""
    forest = state.forest
    num_trees, half = forest.var_tree.shape
    tree_size = 2 * half
    X = state.X
    n = X.shape[1]
    # Static bin count: cut indices live in [1, max_split] and X is stored as an
    # unsigned integer, so the type's range bounds every admissible cut.
    num_bins = int(jnp.iinfo(X.dtype).max) + 1

    leaf_indices = apply_moves_to_leaf_indices(forest.leaf_indices, forest.to_prune, forest.move_node)
    idx = leaf_indices[tree]
    var_old = forest.var_tree[tree]
    split_old = forest.split_tree[tree]
    leaf_old = forest.leaf_tree[tree]

    e = scaled_error_cov_inv(state).astype(jnp.float32)
    lam = forest.leaf_prior_cov_inv.astype(jnp.float32)
    leaf_unit = forest.leaf_unit.astype(jnp.float32)
    resid_unit = state.resid_unit.astype(jnp.float32)
    w = jnp.ones(n, jnp.float32) if state.prec_scale is None else state.prec_scale.astype(jnp.float32)
    resid32 = state.resid.astype(jnp.float32)
    prev_leaf = leaf_unit * leaf_old.astype(jnp.float32)
    t = resid32 * resid_unit + prev_leaf[idx]          # target for this tree, data units
    p_nt = forest.p_nonterminal
    depths = grove.tree_depths(tree_size).astype(jnp.int32)
    idx_depth = depths[idx]
    min_leaf = forest.min_points_per_leaf
    min_decision = forest.min_points_per_decision_node
    blocked = forest.blocked_vars
    capacity = REGROW_MAX_VISITS + 2

    def candidates(node, member, var_arr, split_arr):
        return _node_candidates(node, member, X, w, t, var_arr, split_arr, forest.max_split,
                                blocked, p_nt, e, lam, min_leaf, num_bins, half)

    def pop(stack, ptr):
        active = ptr > 0
        node = jnp.where(active, stack[jnp.maximum(ptr - 1, 0)], 0).astype(jnp.int32)
        return node, jnp.where(active, ptr - 1, ptr), active

    def push_children(stack, ptr, node, is_split):
        stack = stack.at[ptr].set(jnp.where(is_split, 2 * node + 1, stack[ptr]))
        stack = stack.at[ptr + 1].set(jnp.where(is_split, 2 * node, stack[ptr + 1]))
        return stack, jnp.where(is_split, ptr + 2, ptr)

    # ---- forward pass: grow a new tree from the root, recording log q and log prior
    def gen_body(carry, _):
        var_arr, split_arr, node_of, stack, ptr, log_q, log_p, key = carry
        key, k_choice = random.split(key)
        node, ptr_after_pop, active = pop(stack, ptr)

        def do_gen(_):
            member = (node_of == node).astype(jnp.float32)
            s_stop, s_split, lp_leaf, lp_split, _ = candidates(node, member, var_arr, split_arr)
            logits = jnp.concatenate([jnp.array([s_stop]), s_split.ravel()])
            choice = random.categorical(k_choice, logits)
            is_split = choice > 0
            flat = jnp.maximum(choice - 1, 0)
            v = flat // num_bins
            c = flat % num_bins
            lq = log_q + (logits[choice] - logsumexp(logits))
            lp = log_p + jnp.where(is_split, lp_split.ravel()[flat], lp_leaf)
            va = var_arr.at[node].set(jnp.where(is_split, v, var_arr[node]).astype(var_arr.dtype), mode="drop")
            sa = split_arr.at[node].set(jnp.where(is_split, c, split_arr[node]).astype(split_arr.dtype), mode="drop")
            go_right = (X[v].astype(jnp.int32) >= c).astype(jnp.int32)
            no = jnp.where(is_split & (node_of == node), 2 * node + go_right, node_of)
            st_new, pt_new = push_children(stack, ptr_after_pop, node, is_split)
            return va, sa, no, st_new, pt_new, lq, lp

        va, sa, no, st_new, pt_new, lq, lp = lax.cond(
            active, do_gen, lambda _: (var_arr, split_arr, node_of, stack, ptr_after_pop, log_q, log_p), None
        )
        return (va, sa, no, st_new, pt_new, lq, lp, key), None

    stack0 = jnp.zeros(capacity, jnp.int32).at[0].set(1)
    init = (jnp.zeros_like(var_old), jnp.zeros_like(split_old), jnp.ones(n, jnp.int32),
            stack0, jnp.int32(1), jnp.float32(0.0), jnp.float32(0.0), key)
    (var_new, split_new, node_of, _, ptr_new, log_q_new, log_p_new, key), _ = lax.scan(
        gen_body, init, None, length=REGROW_MAX_VISITS)
    overflow_new = ptr_new > 0
    node_of = node_of.astype(idx.dtype)

    # ---- reverse pass: density and prior of the current tree under the same proposal
    def rev_body(carry, _):
        stack, ptr, log_q, log_p = carry
        node, ptr_after_pop, active = pop(stack, ptr)

        def do_rev(_):
            d = depths[node]
            member = ((idx_depth >= d)
                      & (jnp.right_shift(idx.astype(jnp.int32), jnp.maximum(idx_depth - d, 0)) == node)).astype(jnp.float32)
            s_stop, s_split, lp_leaf, lp_split, _ = candidates(node, member, var_old, split_old)
            logits = jnp.concatenate([jnp.array([s_stop]), s_split.ravel()])
            is_split = split_old.at[node].get(mode="fill", fill_value=0) > 0
            flat = var_old.at[node].get(mode="fill", fill_value=0).astype(jnp.int32) * num_bins + split_old.at[node].get(
                mode="fill", fill_value=0).astype(jnp.int32)
            chosen = jnp.where(is_split, logits[1 + flat], logits[0])
            lq = log_q + (chosen - logsumexp(logits))
            lp = log_p + jnp.where(is_split, lp_split.ravel()[flat], lp_leaf)
            st_new, pt_new = push_children(stack, ptr_after_pop, node, is_split)
            return st_new, pt_new, lq, lp

        return lax.cond(active, do_rev, lambda _: (stack, ptr_after_pop, log_q, log_p), None), None

    (_, ptr_old, log_q_old, log_p_old), _ = lax.scan(
        rev_body, (stack0, jnp.int32(1), jnp.float32(0.0), jnp.float32(0.0)), None, length=REGROW_MAX_VISITS)
    overflow_old = ptr_old > 0

    # ---- leaf-integrated likelihoods of both trees
    def leaf_stats(index_vec):
        P = jnp.zeros(tree_size, jnp.float32).at[index_vec].add(w)
        r = jnp.zeros(tree_size, jnp.float32).at[index_vec].add(w * t)
        return P, r

    def logml(P, r):
        prec_post = lam + e * P
        return 0.5 * jnp.log(lam / prec_post) + 0.5 * e * e * r * r / prec_post

    P_old, r_old = leaf_stats(idx)
    P_new, r_new = leaf_stats(node_of)
    nodes = jnp.arange(tree_size)
    own_split_new = jnp.concatenate([split_new, jnp.zeros(half, split_new.dtype)])
    own_split_old = jnp.concatenate([split_old, jnp.zeros(half, split_old.dtype)])
    exists_new = (nodes == 1) | (own_split_new[nodes // 2] > 0)
    exists_old = (nodes == 1) | (own_split_old[nodes // 2] > 0)
    is_leaf_new = exists_new & (own_split_new == 0)
    is_leaf_old = exists_old & (own_split_old == 0)
    ll_new = jnp.sum(jnp.where(is_leaf_new, logml(P_new, r_new), 0.0))
    ll_old = jnp.sum(jnp.where(is_leaf_old, logml(P_old, r_old), 0.0))

    log_ratio = (log_p_new + ll_new - log_q_new) - (log_p_old + ll_old - log_q_old)
    key, k_u, k_leaf = random.split(key, 3)
    accept = (jnp.log(random.uniform(k_u)) <= log_ratio) & ~overflow_new & ~overflow_old

    # ---- apply: exact leaf draws, caches, affluence, residual
    var_post = jnp.reciprocal(lam + e * P_new)
    leaf_data = e * r_new * var_post + random.normal(k_leaf, (tree_size,), jnp.float32) * jnp.sqrt(var_post)
    leaf_store = jnp.where(is_leaf_new, leaf_data / leaf_unit, 0.0)
    quantization = state.config.leaf_quantization
    if quantization is not None:
        nmant = jnp.finfo(state.resid.dtype).nmant
        quantum = (state.resid_eff_scale / leaf_unit) * 2.0 ** (quantization - nmant)
        leaf_store = jnp.round(leaf_store / quantum) * quantum
    finfo = jnp.finfo(leaf_old.dtype)
    leaf_store = lax.reduce_precision(leaf_store, finfo.nexp, finfo.nmant).astype(leaf_old.dtype)

    leaf_out = jnp.where(accept, leaf_store, leaf_old)
    idx_out = jnp.where(accept, node_of, idx)
    var_out = jnp.where(accept, var_new, var_old)
    split_out = jnp.where(accept, split_new, split_old)
    new_fit = leaf_unit * leaf_out.astype(jnp.float32)
    delta = (prev_leaf[idx] - new_fit[idx_out]) / resid_unit
    resid_out = (resid32 + delta).astype(state.resid.dtype)

    counts_new = jnp.zeros(tree_size, jnp.float32).at[node_of].add(1.0)
    p = X.shape[0]

    def node_growable(node):
        lr = jax.vmap(lambda v: split_range(var_new, split_new, forest.max_split, node, v))(jnp.arange(p))
        return (node < half) & (jnp.sum(lr[1] > lr[0]) > 0)

    growable_new = jax.vmap(node_growable)(jnp.arange(half))
    aff_new = is_leaf_new[:half] & growable_new
    if min_decision is not None:
        aff_new &= counts_new[:half] >= min_decision
    affl_out = jnp.where(accept, aff_new, forest.affluence_tree[tree])

    forest_out = replace(
        forest,
        var_tree=forest.var_tree.at[tree].set(var_out),
        split_tree=forest.split_tree.at[tree].set(split_out),
        leaf_tree=forest.leaf_tree.at[tree].set(leaf_out),
        leaf_indices=leaf_indices.at[tree].set(idx_out),
        affluence_tree=forest.affluence_tree.at[tree].set(affl_out),
        to_prune=jnp.zeros_like(forest.to_prune),
    )
    if forest.count_tree is not None:
        forest_out = replace(forest_out, count_tree=forest.count_tree.at[tree].set(
            jnp.where(accept, counts_new.astype(forest.count_tree.dtype), forest.count_tree[tree])))
    if forest.prec_tree is not None:
        forest_out = replace(forest_out, prec_tree=forest.prec_tree.at[tree].set(
            jnp.where(accept, P_new.astype(forest.prec_tree.dtype), forest.prec_tree[tree])))
    return replace(state, resid=resid_out, forest=forest_out), accept


def regrow_step(key: Key[Array, ''], state: State, num_regrow: int) -> State:
    """Regrow ``num_regrow`` distinct random trees of a single-chain state, one at a time."""
    num_trees = state.forest.var_tree.shape[0]
    key, k_perm = random.split(key)
    order = random.permutation(k_perm, num_trees)[:num_regrow]
    keys = random.split(key, num_regrow)

    def body(state, xs):
        k, tree = xs
        state, acc = regrow_tree(k, state, tree)
        return state, acc

    state, _ = lax.scan(body, state, (keys, order))
    return state

