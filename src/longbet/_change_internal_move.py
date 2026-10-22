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

"""CHANGE move at any nonterminal node, keeping the tree shape.

The CHANGE move of ``_change_move`` re-draws a rule only at a parent of two
leaves. On threshold-shaped surfaces the decisive rule is often an internal
node with a subtree beneath it: an ``x1`` rule at the root with ``x2`` rules
below. Neither GROW/PRUNE nor the leaf-parent CHANGE can move that cutpoint
without first pruning the whole subtree, and the REGROW independence proposal
rarely reproduces an interaction exactly, so its acceptance for such trees is
low. Chains then freeze at whichever internal cutpoints they found first, each
chain at its own, and units that only those cutpoints separate never move.

This move picks a nonterminal node uniformly, proposes a new rule there, keeps
every other rule and the whole tree shape, re-routes the points of that node's
subtree, and accepts with the exact Metropolis-Hastings ratio of tree prior
times leaf-integrated likelihood. The rule proposal is the one used by the
leaf-parent CHANGE: with probability ``P_PERTURB`` a symmetric local shift of
the cutpoint, otherwise a fresh rule from the ancestor-restricted rule prior at
that node, whose density is divided out explicitly. The tree prior is evaluated
in full for the old and the new tree, which accounts for the descendants whose
admissible ranges (and hence rule-prior densities) change with the new
ancestor rule; a descendant rule that becomes inadmissible makes the proposal
invalid. The count threshold is applied to every leaf of the new tree, which
vetoes symmetrically because the current tree satisfies it. Accepted trees get
exact conditional leaf draws and an exact residual update, as in REGROW.
Univariate outcomes with optional per-datapoint precision scales.
"""

from __future__ import annotations

from dataclasses import replace
from typing import NamedTuple

import jax
import jax.numpy as jnp
from jax import lax, random
from jaxtyping import Array, Key

from bartz import grove
from bartz.mcmcstep import State
from bartz.mcmcstep._moves import choose_split, choose_variable, randint_masked, split_range
from bartz.mcmcstep._state import scaled_error_cov_inv
from bartz.mcmcstep._step import apply_moves_to_leaf_indices

from longbet._change_move import P_PERTURB, PERTURB_FRACTION


def tree_log_prior(var_tree, split_tree, max_split, blocked, p_nt, half):
    """Log prior of one tree under the bartz tree prior; ``-inf`` if a rule is inadmissible.

    Sums, over the nodes that exist, ``log p_nt`` plus the ancestor-restricted
    uniform rule prior at nonterminal nodes and ``log(1 - p_nt * growable)`` at
    leaves. Bottom-level nodes cannot grow and contribute zero.
    """
    p = max_split.shape[0]
    all_vars = jnp.arange(p)

    def node_terms(k):
        exists = (k == 1) | ((k > 0) & (split_tree[k // 2] > 0))

        def eval_node(_):
            l_v, r_v = jax.vmap(lambda v: split_range(var_tree, split_tree, max_split, k, v))(all_vars)
            var_ok = r_v > l_v
            if blocked is not None:
                var_ok &= ~jnp.isin(all_vars, blocked)
            num_avail = jnp.sum(var_ok)
            growable = num_avail > 0
            is_split = split_tree[k] > 0
            v = var_tree[k].astype(jnp.int32)
            c = split_tree[k].astype(jnp.int32)
            l_k, r_k = l_v[v], r_v[v]
            valid_rule = var_ok[v] & (c >= l_k) & (c < r_k)
            pnt = jnp.where(growable, p_nt[k], 0.0)
            lp_split = (jnp.log(jnp.maximum(pnt, 1e-30))
                        - jnp.log(jnp.maximum(num_avail, 1).astype(jnp.float32))
                        - jnp.log(jnp.maximum(r_k - l_k, 1).astype(jnp.float32)))
            lp_split = jnp.where(valid_rule & (pnt > 0), lp_split, -jnp.inf)
            lp_leaf = jnp.log1p(-pnt)
            return jnp.where(is_split, lp_split, lp_leaf)

        return lax.cond(exists, eval_node, lambda _: jnp.float32(0.0), None)

    return jnp.sum(jax.vmap(node_terms)(jnp.arange(half)))


def _propose(key, var_tree, split_tree, max_split, blocked_vars, log_s):
    """Propose a new rule at a uniformly chosen nonterminal node of one tree."""
    k_node, k_var, k_split, k_u, k_leaf, k_type, k_delta = random.split(key, 7)
    is_internal = split_tree > 0
    allowed = jnp.count_nonzero(is_internal) > 0
    node = jnp.where(allowed, randint_masked(k_node, is_internal), 1).astype(jnp.int32)
    var_old = var_tree[node].astype(jnp.int32)
    split_old = split_tree[node].astype(jnp.int32)
    l_old, r_old = split_range(var_tree, split_tree, max_split, node, var_old)

    var_draw, _num_avail = choose_variable(
        k_var, var_tree, split_tree, max_split, node, blocked_vars, log_s
    )
    split_draw, l_draw, r_draw = choose_split(k_split, var_draw, var_tree, split_tree, max_split, node)
    allowed_draw = l_draw < r_draw

    width = jnp.maximum(1, jnp.ceil(PERTURB_FRACTION * max_split[var_old].astype(jnp.float32))).astype(jnp.int32)
    delta = random.randint(k_delta, (), -width, width)
    delta = jnp.where(delta >= 0, delta + 1, delta)
    split_pert = split_old + delta
    allowed_pert = (l_old <= split_pert) & (split_pert < r_old)

    use_perturb = random.uniform(k_type) < P_PERTURB
    var_new = jnp.where(use_perturb, var_old, var_draw)
    split_new = jnp.where(use_perturb, split_pert, split_draw)
    allowed &= jnp.where(use_perturb, allowed_pert, allowed_draw)
    # Fresh rules are drawn from the rule prior at the node: q(T -> T') is the
    # prior density of the new rule and q(T' -> T) that of the old one. The
    # variable factor is the same in both directions (same ancestors).
    log_q_corr = jnp.where(
        use_perturb, 0.0,
        jnp.log(jnp.maximum(r_draw - l_draw, 1).astype(jnp.float32))
        - jnp.log(jnp.maximum(r_old - l_old, 1).astype(jnp.float32)),
    )
    u = random.uniform(k_u)
    logu = jnp.log(u)
    return node, var_new, split_new, allowed, log_q_corr, logu, k_leaf


class _TreeConsts(NamedTuple):
    """Per-forest constants shared by every tree of one internal-CHANGE sweep."""
    X: Array
    e: Array
    lam: Array
    leaf_unit: Array
    resid_unit: Array
    w: Array
    p_nt: Array
    max_split: Array
    blocked: Array | None
    half: int
    depths: Array
    min_leaf: int | None


def _forest_consts(state: State) -> _TreeConsts:
    forest = state.forest
    _num_trees, half = forest.var_tree.shape
    n = state.X.shape[1]
    return _TreeConsts(
        X=state.X,
        e=scaled_error_cov_inv(state).astype(jnp.float32),
        lam=forest.leaf_prior_cov_inv.astype(jnp.float32),
        leaf_unit=forest.leaf_unit.astype(jnp.float32),
        resid_unit=state.resid_unit.astype(jnp.float32),
        w=jnp.ones(n, jnp.float32) if state.prec_scale is None else state.prec_scale.astype(jnp.float32),
        p_nt=forest.p_nonterminal,
        max_split=forest.max_split,
        blocked=forest.blocked_vars,
        half=int(half),
        depths=grove.tree_depths(2 * half).astype(jnp.int32),
        min_leaf=forest.min_points_per_leaf,
    )


def route_subtree(c: _TreeConsts, node, var_arr, split_arr, idx):
    """Leaf of every point under ``var_arr``/``split_arr``; only points in the subtree of ``node`` move."""
    p, n = c.X.shape
    max_depth = int(c.half).bit_length()
    idx32 = idx.astype(jnp.int32)
    idx_depth = c.depths[idx32]
    d = c.depths[node]
    in_sub = (idx_depth >= d) & (jnp.right_shift(idx32, jnp.maximum(idx_depth - d, 0)) == node)
    point_ids = jnp.arange(n)

    def descend(_, k):
        sp = split_arr.at[k].get(mode="fill", fill_value=0).astype(jnp.int32)

        def step_down(k_in):
            v = var_arr.at[k_in].get(mode="fill", fill_value=0).astype(jnp.int32)
            xv = c.X[jnp.clip(v, 0, p - 1), point_ids].astype(jnp.int32)
            go_right = (xv >= sp).astype(jnp.int32)
            return jnp.where(sp > 0, 2 * k_in + go_right, k_in)

        return lax.cond(jnp.any(sp > 0), step_down, lambda k_in: k_in, k)

    k = lax.fori_loop(0, max_depth, descend, jnp.where(in_sub, node, idx32))
    return jnp.where(in_sub, k, idx32).astype(idx.dtype)


def evaluate_change(c: _TreeConsts, leaf_tree, idx, resid, var_tree, split_tree, node, var_new, split_new):
    """Everything the acceptance of one proposed rule change needs, for one tree.

    Returns a dict with the old and new tree log priors, leaf-integrated log
    likelihoods, the new leaf memberships and their per-node sufficient
    statistics, the leaf mask (shared by both trees) and the count-threshold flag.
    """
    tree_size = 2 * c.half
    var_tree_new = var_tree.at[node].set(var_new.astype(var_tree.dtype))
    split_tree_new = split_tree.at[node].set(split_new.astype(split_tree.dtype))
    new_idx = route_subtree(c, node, var_tree_new, split_tree_new, idx)
    lp_old = tree_log_prior(var_tree, split_tree, c.max_split, c.blocked, c.p_nt, c.half)
    lp_new = tree_log_prior(var_tree_new, split_tree_new, c.max_split, c.blocked, c.p_nt, c.half)

    resid32 = resid.astype(jnp.float32)
    prev_leaf = c.leaf_unit * leaf_tree.astype(jnp.float32)
    t = resid32 * c.resid_unit + prev_leaf[idx]
    wt = c.w * t

    def logml(P, r):
        prec_post = c.lam + c.e * P
        return 0.5 * jnp.log(c.lam / prec_post) + 0.5 * c.e * c.e * r * r / prec_post

    nodes = jnp.arange(tree_size)
    own_split = jnp.concatenate([split_tree, jnp.zeros(c.half, split_tree.dtype)])
    exists = (nodes == 1) | (own_split[nodes // 2] > 0)
    is_leaf = exists & (own_split == 0)          # the shape is unchanged by the move
    old_PR = jnp.zeros((tree_size, 2), jnp.float32).at[idx].add(jnp.stack([c.w, wt], axis=-1))
    P_old, r_old = old_PR[:, 0], old_PR[:, 1]
    new_PRC = jnp.zeros((tree_size, 3), jnp.float32).at[new_idx].add(
        jnp.stack([c.w, wt, jnp.ones_like(c.w)], axis=-1)
    )
    P_new, r_new, counts_new = new_PRC[:, 0], new_PRC[:, 1], new_PRC[:, 2]
    ll_old = jnp.sum(jnp.where(is_leaf, logml(P_old, r_old), 0.0))
    ll_new = jnp.sum(jnp.where(is_leaf, logml(P_new, r_new), 0.0))
    min_leaf_ok = jnp.bool_(True)
    if c.min_leaf is not None:
        min_leaf_ok = jnp.all(jnp.where(is_leaf, counts_new >= c.min_leaf, True))
    return dict(
        var_tree_new=var_tree_new, split_tree_new=split_tree_new, new_idx=new_idx,
        lp_old=lp_old, lp_new=lp_new, ll_old=ll_old, ll_new=ll_new,
        P_new=P_new, r_new=r_new, counts_new=counts_new, is_leaf=is_leaf, min_leaf_ok=min_leaf_ok,
        prev_leaf=prev_leaf, resid32=resid32,
    )


def change_internal_step(key: Key[Array, ''], state: State) -> State:
    """Apply one internal-CHANGE sweep to every tree of a single-chain ``bartz`` state."""
    return change_internal_sweep(key, state)[0]


def change_internal_sweep(key: Key[Array, ''], state: State):
    """As :func:`change_internal_step`, also returning per-tree acceptance flags."""
    forest = state.forest
    num_trees, half = forest.var_tree.shape
    tree_size = 2 * half
    p = state.X.shape[0]
    c = _forest_consts(state)

    leaf_indices = apply_moves_to_leaf_indices(
        forest.leaf_indices, forest.to_prune, forest.move_node
    )
    keys = random.split(key, num_trees)
    proposals = jax.vmap(_propose, in_axes=(0, 0, 0, None, None, None))(
        keys, forest.var_tree, forest.split_tree, forest.max_split, forest.blocked_vars, forest.log_s
    )

    has_count = forest.count_tree is not None
    has_prec = forest.prec_tree is not None
    min_decision = forest.min_points_per_decision_node
    leaf_dtype = forest.leaf_tree.dtype
    resid_dtype = state.resid.dtype
    quantization = state.config.leaf_quantization
    resid_eff_scale = state.resid_eff_scale

    def body(resid, xs):
        (leaf_tree, idx, var_tree, split_tree, affl, count_tree, prec_tree,
         node, var_new, split_new, allowed, log_q_corr, logu, k_leaf) = xs

        def do_allowed(_):
            ev = evaluate_change(c, leaf_tree, idx, resid, var_tree, split_tree, node, var_new, split_new)
            ok = ev["min_leaf_ok"] & jnp.isfinite(ev["lp_new"]) & jnp.isfinite(ev["lp_old"])
            log_ratio = jnp.where(
                ok, (ev["lp_new"] + ev["ll_new"]) - (ev["lp_old"] + ev["ll_old"]) + log_q_corr, -jnp.inf
            )
            accept = ok & (logu <= log_ratio)

            def do_accept(_):
                is_leaf, P_new, r_new, new_idx = ev["is_leaf"], ev["P_new"], ev["r_new"], ev["new_idx"]
                var_post = jnp.reciprocal(c.lam + c.e * P_new)
                leaf_data = c.e * r_new * var_post + random.normal(k_leaf, (tree_size,), jnp.float32) * jnp.sqrt(var_post)
                leaf_store = jnp.where(is_leaf, leaf_data / c.leaf_unit, 0.0)
                if quantization is not None:
                    nmant = jnp.finfo(resid_dtype).nmant
                    quantum = (resid_eff_scale / c.leaf_unit) * 2.0 ** (quantization - nmant)
                    leaf_store = jnp.round(leaf_store / quantum) * quantum
                finfo = jnp.finfo(leaf_dtype)
                leaf_out = lax.reduce_precision(leaf_store, finfo.nexp, finfo.nmant).astype(leaf_dtype)
                idx_out = new_idx
                var_out = ev["var_tree_new"]
                split_out = ev["split_tree_new"]
                new_fit = c.leaf_unit * leaf_out.astype(jnp.float32)
                delta = (ev["prev_leaf"][idx] - new_fit[idx_out]) / c.resid_unit
                resid_out = (ev["resid32"] + delta).astype(resid_dtype)

                def node_growable(k):
                    lr = jax.vmap(lambda v: split_range(var_out, split_out, forest.max_split, k, v))(jnp.arange(p))
                    ok_v = lr[1] > lr[0]
                    if forest.blocked_vars is not None:
                        ok_v &= ~jnp.isin(jnp.arange(p), forest.blocked_vars)
                    return (k < half) & (jnp.sum(ok_v) > 0)

                growable_out = jax.vmap(node_growable)(jnp.arange(half))
                counts_out = ev["counts_new"]
                aff_out = is_leaf[:half] & growable_out
                if min_decision is not None:
                    aff_out &= counts_out[:half] >= min_decision
                count_out = counts_out.astype(count_tree.dtype) if has_count else None
                prec_out = P_new.astype(prec_tree.dtype) if has_prec else None
                return resid_out, (leaf_out, idx_out, var_out, split_out, aff_out, count_out, prec_out, jnp.bool_(True))

            return lax.cond(
                accept,
                do_accept,
                lambda _: (resid, (leaf_tree, idx, var_tree, split_tree, affl, count_tree, prec_tree, jnp.bool_(False))),
                None,
            )

        return lax.cond(
            allowed,
            do_allowed,
            lambda _: (resid, (leaf_tree, idx, var_tree, split_tree, affl, count_tree, prec_tree, jnp.bool_(False))),
            None,
        )

    xs = (forest.leaf_tree, leaf_indices, forest.var_tree, forest.split_tree, forest.affluence_tree,
          forest.count_tree, forest.prec_tree) + tuple(proposals)
    resid, (leaf_tree, leaf_indices, var_tree, split_tree, affluence_tree, count_tree, prec_tree, accepted) = (
        lax.scan(body, state.resid, xs)
    )
    forest = replace(
        forest,
        leaf_tree=leaf_tree,
        leaf_indices=leaf_indices,
        var_tree=var_tree,
        split_tree=split_tree,
        affluence_tree=affluence_tree,
        count_tree=count_tree,
        prec_tree=prec_tree,
        to_prune=jnp.zeros_like(forest.to_prune),
    )
    return replace(state, resid=resid, forest=forest), accepted
