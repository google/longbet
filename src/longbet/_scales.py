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

"""Per-sweep observation scale reconstruction for the weighted treatment forest.

`bartz` derives ``prec_scale``, ``inv_sdev_scale``, ``inv_sdev_unit``,
``n_non_missing`` and ``sum_diag_prec_scale`` from ``error_scale`` inside
``init``, and ``init`` runs once.  The treatment forest's weights
``w = b_Z * beta_S`` change every sweep, so these have to be rebuilt each time.

This module mirrors ``bartz.mcmcstep._state.compute_scale_related_attrs``
rather than importing it, which keeps the private surface small and lets an
equivalence test (``tests/test_scales.py``) catch upstream drift on every bump.

Note that ``w ** 2`` is *not* used directly: if ``beta`` drifts small it
underflows in float32.  Following `bartz`, the magnitude is factored into a
power-of-two ``inv_sdev_unit`` first.
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp
from jaxtyping import Array, Bool, Float32, Int32


class ScaleRelatedAttrs(NamedTuple):
    """Scale attributes required by ``bartz.mcmcstep.State``.

    Field names, shapes and dtypes match
    ``bartz.mcmcstep._state.ScaleRelatedAttrs`` for the single-outcome case.
    """

    inv_sdev_scale: Float32[Array, '... n']
    prec_scale: Float32[Array, '... n']
    inv_sdev_unit: Float32[Array, '...']
    n_non_missing: Int32[Array, '...']
    sum_diag_prec_scale: Float32[Array, '...']


def round_to_pow2(x: Float32[Array, '...']) -> Float32[Array, '...']:
    """Round to the nearest power of two."""
    return 2.0 ** jnp.round(jnp.log2(jnp.maximum(x, 1e-30)))


def compute_treatment_scale_attrs(
    w: Float32[Array, '... n'],
    obs_mask: Bool[Array, ' n'],
    eps: float = 1e-6,
    prec_scale_dtype: jnp.dtype = jnp.float32,
) -> tuple[ScaleRelatedAttrs, Bool[Array, '... n']]:
    """Compute per-sweep scale attributes and mask for the treatment forest.

    Equivalent to ``bartz``'s ``init(error_scale=1/|w|, missing=~mask_nu)`` but
    computed from ``w`` directly and safe against underflow.

    Parameters
    ----------
    w
        Per-observation treatment multiplier ``w = b_Z * beta_S``, shape
        ``(..., n)``.
    obs_mask
        Boolean mask of observed cells (False where ``y`` is missing), shape ``(n,)``.
    eps
        Threshold below which a treatment weight is treated as zero, so the cell
        carries no information about ``nu`` and is masked out of its forest.
    prec_scale_dtype
        Dtype for ``prec_scale`` and ``inv_sdev_scale``.

    Returns
    -------
    attrs : ScaleRelatedAttrs
        Scale attributes to install on the treatment view.
    mask_nu : Bool[Array, '... n']
        Active mask for the treatment forest: ``obs_mask & (|w| > eps)``.
    """
    mask_nu = obs_mask & (jnp.abs(w) > eps)

    inv_sdev_scale = jnp.where(mask_nu, jnp.abs(w), 0.0)

    n_non_missing = jnp.sum(inv_sdev_scale != 0, axis=-1)
    sum_diag_prec_scale = jnp.einsum(
        '...n,...n->...', inv_sdev_scale, inv_sdev_scale
    )

    inv_sdev_unit = round_to_pow2(
        jnp.sqrt(sum_diag_prec_scale / jnp.maximum(n_non_missing, 1))
    )
    # Match bartz: the stored unit is the clamped one, never zero.
    inv_sdev_unit = jnp.where(inv_sdev_unit, inv_sdev_unit, 1.0)

    inv_sdev_scale = inv_sdev_scale / inv_sdev_unit[..., None]
    prec_scale = jnp.square(inv_sdev_scale)

    attrs = ScaleRelatedAttrs(
        inv_sdev_scale=inv_sdev_scale.astype(prec_scale_dtype),
        prec_scale=prec_scale.astype(prec_scale_dtype),
        inv_sdev_unit=inv_sdev_unit.astype(jnp.float32),
        n_non_missing=n_non_missing.astype(jnp.int32),
        sum_diag_prec_scale=sum_diag_prec_scale.astype(jnp.float32),
    )
    return attrs, mask_nu
