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

"""Full-precision conditional likelihood for the identified triangular SUR model.

With raw residuals r=y*-f, B=I-Gamma, and independent innovations e=B r with
variances v, the likelihood precision is B.T diag(1/v) B. Every observed
innovation containing r_m contributes to outcome m's conditional, not only
its own equation. Binary and ordinal rows of Gamma are zero (marginal probit scale one).

Observation masks must be closed under predecessors: if a continuous equation
is observed, all its preceding outcomes are observed. Unobserved downstream
innovations then integrate to one. Discrete-only masks need not be nested,
because those innovation rows have no predecessors. Input validation enforces
this policy; the formulas below must NOT be used for arbitrary missingness.
"""
from __future__ import annotations

import jax.numpy as jnp

SAMPLER_SEMANTICS = "full_precision_sur_v1"


def conditional_residual(raw, observed, loadings, variances, outcome):
    """Return (conditional residual, cell precision) for a single equation.

    ``raw`` and ``observed`` are (M, n), loadings (M, M), variances (M,).
    Return benign residual=0, precision=1 for unobserved response cells, so
    downstream weighted operations never evaluate zero/zero or infinities.
    All inputs are in the outcomes' internal standardized/latent units.
    """
    B = jnp.eye(raw.shape[0], dtype=raw.dtype) - loadings
    innovations = B @ raw
    weights = observed.astype(raw.dtype) / variances[:, None]
    column = B[:, outcome, None]
    precision = jnp.sum(weights * column**2, axis=0)
    score = jnp.sum(weights * column * innovations, axis=0)
    precision = jnp.where(observed[outcome], precision, 1.0)
    residual = jnp.where(observed[outcome], score / precision, 0.0)
    return residual, precision


def innovation_residual(raw, loadings, outcome):
    """Structural innovation, NOT the full conditional pseudo-residual."""
    return raw[outcome] - loadings[outcome] @ raw
