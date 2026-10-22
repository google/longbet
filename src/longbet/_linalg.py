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

"""Gaussian draws from a precision matrix that may be rank deficient.

Several of LongBet's moves are exact conjugate Gibbs steps in a small subspace:
assemble a precision matrix and a linear term, then draw from the corresponding
Gaussian. The subspaces are built from the data, so nothing guarantees they are
linearly independent. Two ways they collapse in practice:

* A single-period panel. The inter-ensemble subspace spans ``{1, t, s, g t}``,
  but with ``T = 1`` the normalized calendar coordinate is a constant, so the
  second direction is a multiple of the first and the subspace loses rank.
* Catastrophic cancellation. The calendar-time block marginalizes the unit
  intercepts through a Schur complement. The complement of a positive definite
  matrix is positive definite in exact arithmetic; in float32, with cell
  precisions of order ``1 / sigma^2``, the subtraction can return a matrix with
  small negative eigenvalues.

The fix is *not* a fixed jitter on the diagonal. A jitter of ``1e-5`` is
meaningless next to a precision of order ``1e6``: in float32, ``1e6 + 1e-5``
is ``1e6``. That is exactly how these moves used to produce NaN -- Cholesky of
a numerically singular matrix returns NaN, which then propagates through the
residual into every parameter in the model within a couple of sweeps.

Instead we diagonalize and invert only the directions the data actually
identify, with a tolerance *relative* to the largest eigenvalue. Directions
below the tolerance are dropped: the draw is projected onto the identified
subspace rather than being handed a near-infinite variance.

Dropping them is not an approximation that has to be argued away. In each of
these moves a null direction of the precision is a direction along which the
proposed update changes nothing -- the leaf shifts, the fits and the residual
all move by exactly zero -- so sampling it would only inject a huge coefficient
to be multiplied by a zero sensitivity, which is how float32 turns a harmless
degeneracy into garbage. The move also stays reversible, because the subspace
is determined by the tree structures, the design and the data, none of which
this move alters; the reverse move sees the same subspace.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
from jaxtyping import Array, Float32, Key

__all__ = ["sample_from_precision", "mean_from_precision"]

# Relative to the largest eigenvalue. float32 carries about 7 decimal digits,
# so 1e-6 keeps everything that survives the factorization and discards
# everything that is indistinguishable from rounding noise.
_DEFAULT_RTOL = 1e-6


def _pseudo_inverse_eigs(
    precision: Float32[Array, 'k k'], rtol: float
) -> tuple[Float32[Array, ' k'], Float32[Array, 'k k']]:
    """Return ``(inv_eigenvalues, eigenvectors)`` with null directions zeroed."""
    sym = 0.5 * (precision + precision.T)
    evals, evecs = jnp.linalg.eigh(sym)
    # eigh returns ascending eigenvalues, so the last is the largest. Clamping
    # at zero first means an all-negative matrix (which should not happen, but
    # would be a disaster if it did) yields a zero floor and therefore a
    # no-op move rather than a NaN.
    floor = rtol * jnp.maximum(evals[-1], 0.0)
    keep = evals > jnp.maximum(floor, jnp.finfo(jnp.float32).tiny)
    inv = jnp.where(keep, jnp.reciprocal(jnp.where(keep, evals, 1.0)), 0.0)
    return inv, evecs


def mean_from_precision(
    precision: Float32[Array, 'k k'],
    linear_term: Float32[Array, ' k'],
    *,
    rtol: float = _DEFAULT_RTOL,
) -> Float32[Array, ' k']:
    """Solve ``precision @ w = linear_term`` on the identified subspace."""
    if precision.shape[0] == 0:
        return jnp.zeros_like(linear_term)
    inv, evecs = _pseudo_inverse_eigs(precision, rtol)
    return evecs @ (inv * (evecs.T @ linear_term))


def sample_from_precision(
    key: Key[Array, ''],
    precision: Float32[Array, 'k k'],
    linear_term: Float32[Array, ' k'],
    *,
    rtol: float = _DEFAULT_RTOL,
) -> Float32[Array, ' k']:
    """Draw ``w ~ N(P^-1 h, P^-1)``, projected onto the identified subspace.

    Uses the canonical symmetric square root ``V diag(sqrt(inv)) V^T eta`` for
    the noise term so that the draw depends only on the matrix ``P`` and is
    strictly invariant to arbitrary eigenvector sign flips (``v_j -> -v_j``)
    and orthogonal rotations inside repeated-eigenvalue subspaces.
    """
    if precision.shape[0] == 0:
        return jnp.zeros_like(linear_term)
    inv, evecs = _pseudo_inverse_eigs(precision, rtol)
    eta = jax.random.normal(key, shape=linear_term.shape, dtype=jnp.float32)
    rotated_mean = inv * (evecs.T @ linear_term)
    rotated_noise = jnp.sqrt(inv) * (evecs.T @ eta)
    return evecs @ (rotated_mean + rotated_noise)


