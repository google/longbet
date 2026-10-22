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

"""Refresh bartz's leaf precision cache when observation weights change."""

import equinox as eqx
import jax
import jax.numpy as jnp
from bartz.mcmcstep import BatchedReduction
from bartz.mcmcstep._step import apply_moves_to_leaf_indices


def current_forest_fit(forest):
    """Evaluate single-chain stored leaves without residual differencing.

    Pending prunes must be applied to membership, just as for precision sums.
    This avoids cancellation in R/w when a coding or GP multiplier is tiny.
    Returns values for every design cell, including currently zero-weight cells.
    """
    ids = apply_moves_to_leaf_indices(forest.leaf_indices, forest.to_prune, forest.move_node)
    values = jax.vmap(lambda leaf, index: leaf[index])(forest.leaf_tree, ids)
    return forest.offset + forest.leaf_unit*jnp.sum(values.astype(jnp.float32), axis=0)


def refresh_prec_tree(state):
    """Return a forest view whose cached leaf sums match its current weights.

    bartz assumes fixed observation weights and only recomputes precision at
    the nodes involved in the proposed grow/prune. LongBet's treatment weights
    (and full-SUR prognostic weights) change between sweeps, so every occupied
    leaf must be refreshed. Group using the pending-prune-adjusted indices,
    without applying that prune a second time or altering the forest topology.
    Do NOT reuse the configured precision reduction: its default is a dense
    one-hot routine optimized for the two touched children. Applied to every
    leaf it creates an observations-by-leaf-slots array per tree and chain.
    Use an indexed reduction with at most 16 batches instead, retaining the
    configured data-axis collective without that quadratic-size temporary.
    """
    if state.prec_scale is None:
        return state
    forest = state.forest
    indices = apply_moves_to_leaf_indices(
        forest.leaf_indices, forest.to_prune, forest.move_node
    )
    reduction = BatchedReduction(
        num_batches=min(16, max(1, state.prec_scale.shape[-1] // 32))
    )
    prec_tree = jax.vmap(
        lambda ids: reduction._reduce(
            state.prec_scale, ids, size=forest.leaf_tree.shape[-1],
            dtype=jnp.float32, data_sharded=state.config.data_sharded,
        )
    )(indices)
    return eqx.tree_at(
        lambda s: s.forest.prec_tree, state, prec_tree,
        is_leaf=lambda x: x is None,
    )
