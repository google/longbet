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

"""Per-sweep scale reconstruction against bartz's own implementation.

``longbet._scales`` reimplements ``compute_scale_related_attrs`` rather than
importing it, to keep the private surface small. That is only safe with an
equivalence test, which also catches upstream drift on every version bump.
"""

import jax.numpy as jnp
import numpy as np
import pytest
from bartz.mcmcstep._state import compute_scale_related_attrs as bartz_scales

from longbet._scales import ScaleRelatedAttrs, compute_treatment_scale_attrs, round_to_pow2


def test_round_to_pow2():
    x = jnp.array([0.1, 0.7, 1.0, 1.9, 3.8, 8.1], dtype=jnp.float32)
    expected = jnp.array([0.125, 0.5, 1.0, 2.0, 4.0, 8.0], dtype=jnp.float32)
    assert jnp.allclose(round_to_pow2(x), expected, rtol=1e-5)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_matches_bartz_including_dtypes(seed):
    n = 120
    rng = np.random.default_rng(seed)
    w = rng.normal(size=n).astype(np.float32)
    w[:10] = 0.0                      # cells the treatment forest must not see
    w[10:14] *= 1e-9                  # weights below the masking threshold
    obs_mask = np.ones(n, dtype=bool)
    obs_mask[20:25] = False           # genuinely missing outcomes

    attrs, mask_nu = compute_treatment_scale_attrs(jnp.asarray(w), jnp.asarray(obs_mask))

    expected_mask = obs_mask & (np.abs(w) > 1e-6)
    assert np.array_equal(np.asarray(mask_nu), expected_mask)

    reference = bartz_scales(
        jnp.asarray(1.0 / np.where(expected_mask, np.abs(w), 1.0), jnp.float32),
        jnp.asarray(~expected_mask),
        jnp.float32,
        n,
    )
    for name in ScaleRelatedAttrs._fields:
        got, want = getattr(attrs, name), getattr(reference, name)
        assert jnp.allclose(got, want), f"{name} differs from bartz"
        assert got.dtype == want.dtype, (
            f"{name} dtype {got.dtype} != bartz's {want.dtype}"
        )


def test_small_weights_do_not_underflow():
    """prec_scale is w^2; computing it directly would underflow in float32."""
    n = 64
    w = np.full(n, 3e-20, dtype=np.float32)
    attrs, mask = compute_treatment_scale_attrs(jnp.asarray(w), jnp.ones(n, bool))
    assert not np.any(np.asarray(mask)), "weights this small must be masked out"
    assert np.all(np.isfinite(np.asarray(attrs.prec_scale)))

    w = np.full(n, 1e-3, dtype=np.float32)
    attrs, mask = compute_treatment_scale_attrs(jnp.asarray(w), jnp.ones(n, bool))
    assert np.all(np.asarray(mask))
    prec = np.asarray(attrs.prec_scale)
    assert np.all(prec > 0), "float32 underflow in prec_scale"
    # The unit is snapped to a power of two, so the normalised precision lands
    # within a factor of two of 1 however small the raw weights were. Squaring
    # 1e-3 directly would give 1e-6, and 1e-20 would flush to zero.
    assert np.all((prec > 0.5) & (prec < 2.0)), prec[:3]


def test_chain_axis_is_reduced_correctly():
    chains, n = 3, 50
    rng = np.random.default_rng(123)
    w = rng.normal(size=(chains, n)).astype(np.float32)
    obs_mask = np.ones(n, dtype=bool)

    attrs, _ = compute_treatment_scale_attrs(jnp.asarray(w), jnp.asarray(obs_mask))
    assert attrs.prec_scale.shape == (chains, n)
    assert attrs.inv_sdev_unit.shape == (chains,)

    for c in range(chains):
        single, _ = compute_treatment_scale_attrs(jnp.asarray(w[c]), jnp.asarray(obs_mask))
        for name in ScaleRelatedAttrs._fields:
            assert jnp.allclose(getattr(attrs, name)[c], getattr(single, name)), (
                f"{name} reduced over the wrong axis"
            )
