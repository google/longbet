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

"""Opt-in CHANGE move for tree decision rules.

``bartz`` proposes only GROW and PRUNE moves. Moving a decision rule from one
cutpoint to a neighbouring one therefore requires pruning the node and growing
it again, passing through a tree whose likelihood can be far lower. On strongly
identified threshold surfaces the chains freeze at whichever cutpoints they
found first, each chain at a different one. The CHANGE move proposes a fresh
rule at a random parent of two leaves, keeps the tree shape, and accepts with
the exact marginal-likelihood ratio of the two new leaves against the two old
ones. It is applied to every tree of a forest, one tree at a time, with the
residual updated exactly between trees, after the ordinary GROW/PRUNE step.

Validity. The new rule is drawn from the same ancestor-restricted rule prior
``bartz`` uses for GROW (uniform, or ``log_s``-weighted, over admissible
variables and uniform over the admissible cutpoints), so the proposal density
cancels the rule prior. The node stays nonterminal, so its terminality factor
cancels. The children's terminal-probability factors change only through their
rule admissibility, which is included. A count threshold vetoes the move
symmetrically: the reverse move starts from the current state, which satisfies
the threshold. The leaves of an accepted move are drawn from their exact
Gaussian conditionals, and the residual is updated with the stored (rounded)
leaf values, as in ``bartz``. Univariate outcomes only, which is all LongBet
uses; per-datapoint precision scales are supported.
"""

from __future__ import annotations

from dataclasses import replace

import jax
import jax.numpy as jnp
from jax import lax, random
from jaxtyping import Array, Key

from bartz import grove
from bartz.mcmcstep import State
from bartz.mcmcstep._moves import choose_split, choose_variable, randint_masked, split_range
from bartz.mcmcstep._state import scaled_error_cov_inv
from bartz.mcmcstep._step import apply_moves_to_leaf_indices


#: Probability of a local cutpoint perturbation instead of a fresh rule draw.
P_PERTURB = 0.5
#: Half-width of the perturbation window as a fraction of the variable's cutpoint count.
PERTURB_FRACTION = 0.05


def _propose(key, var_tree, split_tree, max_split, blocked_vars, log_s):
    """Propose a new rule at a random parent of two leaves of one tree.

    With probability ``P_PERTURB`` the proposal keeps the variable and shifts
    the cutpoint by a nonzero offset drawn uniformly from ``[-k, k]``, a
    symmetric kernel whose out-of-range proposals are rejected; otherwise the
    rule is re-drawn from the ancestor-restricted rule prior, whose density
    cancels against the prior. Both branches leave the tree prior's rule
    factor out of the acceptance ratio.
    """
    k_node, k_var, k_split, k_u, k_leaf, k_type, k_delta = random.split(key, 7)
    is_nog = grove.is_leaves_parent(split_tree)
    allowed = jnp.count_nonzero(is_nog) > 0
    node = jnp.where(allowed, randint_masked(k_node, is_nog), 1).astype(jnp.int32)
    var_old = var_tree[node].astype(jnp.int32)
    split_old = split_tree[node].astype(jnp.int32)
    l_old, r_old = split_range(var_tree, split_tree, max_split, node, var_old)

    # Fresh rule from the prior.
    var_draw, num_avail = choose_variable(
        k_var, var_tree, split_tree, max_split, node, blocked_vars, log_s
    )
    split_draw, l_draw, r_draw = choose_split(k_split, var_draw, var_tree, split_tree, max_split, node)
    allowed_draw = l_draw < r_draw

    # Local perturbation of the current cutpoint.
    width = jnp.maximum(1, jnp.ceil(PERTURB_FRACTION * max_split[var_old].astype(jnp.float32))).astype(jnp.int32)
    delta = random.randint(k_delta, (), -width, width)
    delta = jnp.where(delta >= 0, delta + 1, delta)
    split_pert = split_old + delta
    allowed_pert = (l_old <= split_pert) & (split_pert < r_old)

    use_perturb = random.uniform(k_type) < P_PERTURB
    var_new = jnp.where(use_perturb, var_old, var_draw)
    split_new = jnp.where(use_perturb, split_pert, split_draw)
    l_new = jnp.where(use_perturb, l_old, l_draw)
    r_new = jnp.where(use_perturb, r_old, r_draw)
    allowed &= jnp.where(use_perturb, allowed_pert, allowed_draw)

    other_vars = num_avail > 1
    growable_new = jnp.stack([other_vars | (l_new < split_new), other_vars | (split_new + 1 < r_new)])
    growable_old = jnp.stack([other_vars | (l_old < split_old), other_vars | (split_old + 1 < r_old)])
    u, exp1mlogu = random.uniform(k_u, (2,))
    logu = jnp.log1p(-exp1mlogu)
    normals = random.normal(k_leaf, (2,), jnp.float32)
    return node, var_new, split_new, allowed, growable_new, growable_old, logu, normals, use_perturb


def change_step(key: Key[Array, ''], state: State) -> State:
    """Apply one CHANGE sweep to every tree of a single-chain ``bartz`` state.

    The state is the forest view LongBet builds for one equation: residual in
    ``resid_unit`` units, optional ``prec_scale`` weights, and a fixed error
    precision. Pending prunes are resolved first, so the returned state carries
    none; ``current_forest_fit`` and the next GROW/PRUNE step remain valid.
    """
    return change_sweep(key, state)[0]


def change_sweep(key: Key[Array, ''], state: State):
    """As :func:`change_step`, also returning per-tree acceptance and move-type flags."""
    forest = state.forest
    num_trees, half = forest.var_tree.shape
    tree_size = 2 * half
    X = state.X
    p = X.shape[0]

    leaf_indices = apply_moves_to_leaf_indices(
        forest.leaf_indices, forest.to_prune, forest.move_node
    )

    keys = random.split(key, num_trees)
    log_s_axis = None if forest.log_s is None else None
    proposals = jax.vmap(_propose, in_axes=(0, 0, 0, None, None, log_s_axis))(
        keys, forest.var_tree, forest.split_tree, forest.max_split, forest.blocked_vars, forest.log_s
    )

    e = scaled_error_cov_inv(state).astype(jnp.float32)
    lam = forest.leaf_prior_cov_inv.astype(jnp.float32)
    leaf_unit = forest.leaf_unit.astype(jnp.float32)
    resid_unit = state.resid_unit.astype(jnp.float32)
    weights = None if state.prec_scale is None else state.prec_scale.astype(jnp.float32)
    has_count = forest.count_tree is not None
    has_prec = forest.prec_tree is not None
    min_leaf = forest.min_points_per_leaf
    min_decision = forest.min_points_per_decision_node
    leaf_dtype = forest.leaf_tree.dtype
    resid_dtype = state.resid.dtype
    quantization = state.config.leaf_quantization
    resid_eff_scale = state.resid_eff_scale

    def log_marginal(P, r):
        prec_post = lam + e * P
        return 0.5 * jnp.log(lam / prec_post) + 0.5 * e * e * r * r / prec_post

    def body(resid, xs):
        (leaf_tree, idx, var_tree, split_tree, affl, count_tree, prec_tree,
         node, var_new, split_new, allowed, g_new, g_old, logu, normals, _use_perturb) = xs

        def do_allowed(_):
            left = 2 * node
            children = jnp.stack([left, left + 1])
            in_sub = (idx >> 1).astype(jnp.int32) == node
            x_new = X[jnp.clip(var_new, 0, p - 1)]
            go_right = (x_new.astype(jnp.int32) >= split_new).astype(jnp.int32)
            new_idx = jnp.where(in_sub, (left + go_right).astype(idx.dtype), idx)

            resid32 = resid.astype(jnp.float32)
            prev_leaf = leaf_unit * leaf_tree.astype(jnp.float32)
            target = resid32 * resid_unit + prev_leaf[idx]
            w = jnp.ones_like(target) if weights is None else weights
            tw = jnp.where(in_sub, target * w, 0.0)
            w_sub = jnp.where(in_sub, w, 0.0)

            def leaf_sums(index_vec, values):
                return jnp.stack([jnp.sum(jnp.where(index_vec == c, values, 0.0)) for c in (left, left + 1)])

            r_old = leaf_sums(idx, tw)
            r_new = leaf_sums(new_idx, tw)
            P_old = leaf_sums(idx, w_sub)
            P_new = leaf_sums(new_idx, w_sub)
            ones_sub = jnp.where(in_sub, 1.0, 0.0)
            c_new = leaf_sums(new_idx, ones_sub)

            log_lik = jnp.sum(log_marginal(P_new, r_new) - log_marginal(P_old, r_old))
            pnt = forest.p_nonterminal[children]
            log_prior = jnp.sum(jnp.log1p(-pnt * g_new) - jnp.log1p(-pnt * g_old))
            ok = jnp.bool_(True)
            if min_leaf is not None:
                ok = jnp.min(c_new) >= min_leaf
            accept = ok & (logu <= log_lik + log_prior)

            # Exact conditional draw of the two new leaves (data units).
            var_post = jnp.reciprocal(lam + e * P_new)
            leaf_data = e * r_new * var_post + normals * jnp.sqrt(var_post)
            leaf_store = leaf_data / leaf_unit
            if quantization is not None:
                nmant = jnp.finfo(resid_dtype).nmant
                quantum = (resid_eff_scale / leaf_unit) * 2.0 ** (quantization - nmant)
                leaf_store = jnp.round(leaf_store / quantum) * quantum
            finfo = jnp.finfo(leaf_dtype)
            leaf_store = lax.reduce_precision(leaf_store, finfo.nexp, finfo.nmant).astype(leaf_dtype)
            new_leaf_tree = leaf_tree.at[children].set(jnp.where(accept, leaf_store, leaf_tree[children]))
            idx_out = jnp.where(accept, new_idx, idx)

            # Residual bookkeeping with the stored leaf values, as bartz does.
            new_leaf = leaf_unit * new_leaf_tree.astype(jnp.float32)
            delta = jnp.where(in_sub, (prev_leaf[idx] - new_leaf[idx_out]) / resid_unit, 0.0)
            resid_out = jnp.where(accept, (resid32 + delta).astype(resid_dtype), resid)

            var_out = var_tree.at[node].set(jnp.where(accept, var_new, var_tree[node]).astype(var_tree.dtype))
            split_out = split_tree.at[node].set(jnp.where(accept, split_new, split_tree[node]).astype(split_tree.dtype))
            aff_new = (children < half) & g_new
            if min_decision is not None:
                aff_new &= c_new >= min_decision
            affl_out = affl.at[children].set(jnp.where(accept, aff_new, affl.at[children].get(mode="fill", fill_value=False)), mode="drop")
            count_out = None
            if has_count:
                count_out = count_tree.at[children].set(jnp.where(accept, c_new.astype(count_tree.dtype), count_tree[children]))
            prec_out = None
            if has_prec:
                prec_out = prec_tree.at[children].set(jnp.where(accept, P_new.astype(prec_tree.dtype), prec_tree[children]))
            return resid_out, (new_leaf_tree, idx_out, var_out, split_out, affl_out, count_out, prec_out, accept)

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
    return replace(state, resid=resid, forest=forest), accepted, proposals[-1]
