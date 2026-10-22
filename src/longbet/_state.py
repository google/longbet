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

"""``LongBetState`` and initialization for the composed LongBet model.

The state subclasses ``bartz.mcmcstep.State`` so that ``X``, ``y`` and
``config`` (mesh, reduction strategies, sharding) are shared with the
prognostic forest rather than duplicated, and so the axis annotations `bartz`
reads are inherited.

Chain handling
--------------
`bartz` supports ``init(num_chains=k)``, but its weighted-observation path does
**not**: ``_compute_count_or_prec_tree`` assumes ``prec_scale`` carries no chain
axis (its second axis is the multi-outcome axis, not chains). LongBet's
treatment forest has per-observation weights ``w = b_Z * beta_S`` that differ
across chains by construction, so ``init(num_chains=k)`` cannot be used for it.

Chains are therefore run with ``jax.vmap`` over single-chain states. To avoid
replicating the design matrix -- which dominates device memory on a large panel
-- only the leaves `bartz` itself marks with a chain axis are replicated;
``X``, ``y``, the panel index vectors, the observation mask and the kernel
factor stay shared and are closed over by the vmapped function.  Which leaves
those are is read from `bartz`'s own field metadata, not hardcoded, so it cannot
drift from upstream (see :func:`chain_filter_spec`).
"""

from __future__ import annotations

import dataclasses
import math
from typing import Any

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
from jaxtyping import Array, Bool, Float32, Int32, UInt

from bartz._jaxext import field
from bartz.mcmcstep import Forest, State, Wishart, init, make_p_nonterminal
from bartz.mcmcstep._axes import CHAIN_AXIS

from longbet._config import LongBetConfig
from longbet._gp import build_kernel_matrix, kernel_cholesky
from longbet._ordinal import full_cutpoints, prepare_ordinal, sample_ordinal_latents
from longbet._trend_basis import (
    init_horseshoe,
    trend_block_step,
    trend_prior_scales,
    trend_sufficient_statistics,
)
from longbet._x64 import enable_x64

#: Which leaves carry a chain axis is **not** hardcoded here.  ``bartz``
#: annotates it on its own dataclass fields (``field(chains=CHAIN_AXIS)``), and
#: :func:`chain_filter_spec` reads that metadata, so LongBet's notion of a chain
#: axis cannot drift from upstream's.  ``Forest.leaf_tree`` is per chain while
#: ``Forest.leaf_unit`` is not, for instance, and getting that distinction wrong
#: silently breaks ``evaluate_trace``.  LongBet's own extra fields declare it the
#: same way in :class:`LongBetState` below.


class LongBetState(State):
    """The full MCMC state for LongBet.

    Inherits ``X``, ``y``, ``z``, ``binary_indices``, ``resid``, ``resid_unit``,
    ``resid_eff_scale``, ``resid_inexact_integral``, ``error_cov_inv``,
    ``error_scale``, ``prec_scale``, ``inv_sdev_scale``, ``inv_sdev_unit``,
    ``n_non_missing``, ``sum_diag_prec_scale``, ``forest`` (mu) and ``config``
    from ``bartz.mcmcstep.State``.

    ``resid`` holds the **full-model residual in data units**, not the scaled
    single-forest residual `bartz` stores there. Conversion happens only at the
    boundary of a `bartz` call (see ``longbet._step``).
    """

    # --- treatment forest nu and its observation scales ---------------------
    forest_nu: Forest = field()
    resid_nu: Float32[Array, '*chains n'] = field(chains=CHAIN_AXIS, data=-1)
    prec_scale_nu: Float32[Array, '*chains n'] | None = field(
        chains=CHAIN_AXIS, data=-1
    )
    inv_sdev_scale_nu: Float32[Array, '*chains n'] | None = field(
        chains=CHAIN_AXIS, data=-1
    )
    inv_sdev_unit_nu: Float32[Array, '*chains'] = field(chains=CHAIN_AXIS)
    resid_unit_nu: Float32[Array, ''] = field()
    leaf_prior_cov_inv_nu: Float32[Array, ''] = field()

    # --- panel indices and masks (M = N * T cells, unit-major) --------------
    unit_idx: Int32[Array, ' n'] = field(data=-1)
    time_idx: Int32[Array, ' n'] = field(data=-1)
    exposure_idx: Int32[Array, ' n'] = field(data=-1)
    z_vec: Float32[Array, ' n'] = field(data=-1)
    obs_mask: Bool[Array, ' n'] = field(data=-1)
    #: Number of observed cells per unit; precomputed since it never changes.
    unit_counts: Float32[Array, ' N'] = field()

    # --- panel dimensions (static) ------------------------------------------
    N_units: int = field(static=True)
    T_periods: int = field(static=True)
    S_max: int = field(static=True)

    # --- model parameters ----------------------------------------------------
    beta: Float32[Array, '*chains S_max_plus_1'] = field(chains=CHAIN_AXIS)
    gamma: Float32[Array, '*chains N'] = field(chains=CHAIN_AXIS)
    alpha: Float32[Array, '*chains'] = field(chains=CHAIN_AXIS)
    b0: Float32[Array, '*chains'] = field(chains=CHAIN_AXIS)
    b1: Float32[Array, '*chains'] = field(chains=CHAIN_AXIS)
    sigma2: Float32[Array, '*chains'] = field(chains=CHAIN_AXIS)
    sigma_gamma2: Float32[Array, '*chains'] = field(chains=CHAIN_AXIS)
    cutpoints: Float32[Array, '*chains K_free'] = field(chains=CHAIN_AXIS)

    # --- cached forest fits, in their own units ------------------------------
    mu_fit: Float32[Array, '*chains n'] = field(chains=CHAIN_AXIS, data=-1)
    nu_fit: Float32[Array, '*chains n'] = field(chains=CHAIN_AXIS, data=-1)

    # --- smooth calendar-time trend block ------------------------------------
    #: Row-wise Kronecker product ``phi(t) (x) b(X_i)``, shape ``(n, q_tot)``.
    #: Built once from fixed bases and **shared** across chains -- that it does
    #: not vary by chain is what makes the block cure the mismatch rather than
    #: relocate it (see ``longbet._trend_basis``). ``q_tot = 0`` disables it.
    trend_design: Float32[Array, 'n q_tot'] = field(data=-1)
    #: Trend coefficients, per chain.
    trend_coef: Float32[Array, '*chains q_tot'] = field(chains=CHAIN_AXIS)

    # --- GP kernel Cholesky factor over 0..S_max (shared, float32) -----------
    K_chol: Float32[Array, 'S_max_plus_1 S_max_plus_1'] = field()

    # --- fields with defaults ------------------------------------------------
    resid_eff_scale_nu: Any = field(chains=CHAIN_AXIS, default=None)
    resid_inexact_integral_nu: Any = field(chains=CHAIN_AXIS, default=None)
    sample_beta: bool = field(static=True, default=True)
    random_intercept: bool = field(static=True, default=True)
    gamma_prior_a: float = field(static=True, default=1.0)
    gamma_prior_b: float = field(static=True, default=0.1)
    sigma_prior_a: float = field(static=True, default=2.0)
    sigma_prior_b: float = field(static=True, default=1.0)
    outcome_type_str: str = field(static=True, default="continuous")
    num_categories: int = field(static=True, default=0)
    cutpoint_prior_scale: float = field(static=True, default=5.0)
    #: Inverse temperature of each replica (1 = the posterior). Per chain.
    temperature: Any = field(chains=CHAIN_AXIS, default=None)
    tempering_levels: int = field(static=True, default=1)
    tempering_beta_min: float = field(static=True, default=0.05)
    use_inter_ensemble_move: bool = field(static=True, default=True)
    use_trend_block: bool = field(static=True, default=False)
    #: Per-column prior standard deviation of ``trend_coef``, shape
    #: ``(q_tot,)``. Shared across chains: it is a fixed feature of the prior,
    #: not a sampled value. Under the horseshoe the interaction entries are the
    #: *global* scale ``A``, and the sampled ``tau lambda_j`` replaces them.
    trend_prior_scale_vec: Any = field(default=None)
    #: How many leading columns of ``trend_design`` are the free common-time
    #: factor. The horseshoe applies only past this point.
    trend_num_period: int = field(static=True, default=0)
    use_trend_horseshoe: bool = field(static=True, default=False)
    #: ``A``, the half-Cauchy scale of the horseshoe's global ``tau``. Static: a
    #: fixed hyperparameter, and :func:`longbet._trend_basis.sample_horseshoe`
    #: needs it as a Python float, not a traced array.
    trend_hs_global_scale: float = field(static=True, default=1.0)
    #: ``lambda_j^2`` and its auxiliary ``nu_j``, per chain, one entry per
    #: interaction column. Per chain because they are sampled.
    trend_local2: Any = field(chains=CHAIN_AXIS, default=None)
    trend_local_aux: Any = field(chains=CHAIN_AXIS, default=None)
    #: ``tau^2`` and its auxiliary ``xi``, per chain, scalars.
    trend_global2: Any = field(chains=CHAIN_AXIS, default=None)
    trend_global_aux: Any = field(chains=CHAIN_AXIS, default=None)
    trend_gram: Any = field(default=None)
    trend_unit_sum: Any = field(default=None)
    collapse_gamma_mu: bool = field(static=True, default=True)
    collapse_gamma_nu: bool = field(static=True, default=True)
    has_unbalanced_units: bool = field(static=True, default=False)
    #: Cells that inform the *baseline module* (``mu``, the trend block,
    #: ``gamma_i`` and ``sigma^2``): every observed cell under
    #: ``baseline_cells='all'``, the observed untreated cells under
    #: ``'untreated'``. The residual ``resid`` is maintained on ``obs_mask``
    #: either way. Shared across chains.
    base_mask: Any = field(default=None)
    #: Number of baseline-informing cells per unit, shape ``(N,)``.
    unit_counts_base: Any = field(default=None)
    #: Whether the baseline module is cut from the treated cells
    #: (``baseline_cells='untreated'``). Static: it selects code paths.
    baseline_untreated_only: bool = field(static=True, default=False)

    # --- encouragement-v2 extensions (None / False when off) -----------------
    #: Linear part of the moderator (``LongBetConfig.linear_trt``):
    #: ``nu(x) = x' delta + trees``. ``lin_design`` holds the standardized
    #: treatment covariates per cell, shape ``(n, p)``, shared across chains;
    #: ``lin_fit = lin_design @ delta`` is cached per chain. ``nu_fit`` stays
    #: the *tree* part only, so every forest update is unchanged.
    use_linear_trt: bool = field(static=True, default=False)
    lin_design: Any = field(data=0, default=None)
    delta: Any = field(chains=CHAIN_AXIS, default=None)
    lin_fit: Any = field(chains=CHAIN_AXIS, data=-1, default=None)
    #: ``tau_delta^2`` and its auxiliary in the half-Cauchy mixture
    #: ``tau^2 | xi ~ IG(1/2, 1/xi)``, ``xi ~ IG(1/2, 1/A^2)``.
    tau_delta2: Any = field(chains=CHAIN_AXIS, default=None)
    tau_delta_aux: Any = field(chains=CHAIN_AXIS, default=None)
    linear_trt_global_scale: float = field(static=True, default=1.0)
    #: Horseshoe local scales ``lambda_j^2`` and their auxiliaries.
    lin_local2: Any = field(chains=CHAIN_AXIS, default=None)
    lin_local_aux: Any = field(chains=CHAIN_AXIS, default=None)
    #: Learned trajectory scale (``LongBetConfig.learn_beta_scale``):
    #: ``beta ~ N(0, s^2 K_tilde)``; ``beta_scale2 = s^2`` and its auxiliary in
    #: the half-Cauchy mixture for ``s ~ C+(0, beta_scale_prior)``.
    learn_beta_scale: bool = field(static=True, default=False)
    beta_scale2: Any = field(chains=CHAIN_AXIS, default=None)
    beta_scale_aux: Any = field(chains=CHAIN_AXIS, default=None)
    beta_scale_prior: float = field(static=True, default=1.0)
    #: Accepted beta-nu ridge moves, per chain (instrumentation).
    ridge_accepts: Any = field(chains=CHAIN_AXIS, default=None)
    #: Stationary AR(1) innovations (``LongBetConfig.ar1_errors``, see
    #: ``longbet._ar1``): ``phi`` per chain, each unit's first-period cells
    #: (shared), and the accepted ``phi`` proposals per chain.
    ar1: bool = field(static=True, default=False)
    ar1_phi: Any = field(chains=CHAIN_AXIS, default=None)
    ar1_first: Any = field(data=-1, default=None)
    ar1_phi_max: float = field(static=True, default=0.95)
    ar1_accepts: Any = field(chains=CHAIN_AXIS, default=None)
    #: Fixed pieces of the AR(1)-transformed trend cross-products (shared).
    ar1_trend_stats: Any = field(default=None)
    #: Interweaving update of sigma_gamma (``interweave_gamma``) and its
    #: accepted proposals per chain.
    interweave_gamma: bool = field(static=True, default=False)
    interweave_accepts: Any = field(chains=CHAIN_AXIS, default=None)
    anchor_moderator: bool = field(static=True, default=False)

    # --- distributed-lag exposure kernel (``LongBet.fit(exposure_lags=...)``) --
    #: Whether the cell weight is ``w = lag_design @ kappa`` in place of
    #: ``b_z beta[exposure_idx]``. Static: it selects code paths. ``beta``,
    #: ``exposure_idx`` and ``K_chol`` are then inert placeholders.
    use_lag_kernel: bool = field(static=True, default=False)
    #: Lag-dose rows ``(D_{i,t}, ..., D_{i,t-L})`` per cell, shape ``(n, L+1)``,
    #: float32, shared across chains, masked rows zero.
    lag_design: Any = field(data=0, default=None)
    #: The lag kernel ``kappa``, per chain, shape ``(L+1,)``.
    kappa: Any = field(chains=CHAIN_AXIS, default=None)
    #: Factor of the kernel prior: ``kappa = s L_kappa u``, ``u ~ N(0, I_r)``,
    #: shape ``(L+1, r)``, shared. The lower Cholesky factor of ``K_L`` unless
    #: a caller supplies its own (possibly non-square) factor.
    lag_K_chol: Any = field(default=None)
    #: Precomputed whitened lag-dose design ``D_lag_mat @ K_chol`` (at unit
    #: scale ``s = 1``), shape ``(n, r)``, float64, shared across chains.
    D_lag_K_chol: Any = field(data=0, default=None)
    #: Lower Cholesky factor of ``L_kappa L_kappa'`` (the prior covariance at
    #: ``s = 1``), for the prior quadratic form of the scale step.
    lag_K_chol_sq: Any = field(default=None)
    #: Learned kernel scale ``s^2`` and its half-Cauchy auxiliary, per chain.
    lag_learn_scale: bool = field(static=True, default=False)
    lag_scale2: Any = field(chains=CHAIN_AXIS, default=None)
    lag_scale_aux: Any = field(chains=CHAIN_AXIS, default=None)
    lag_scale_prior: float = field(static=True, default=1.0)
    # --- fixed extra columns of the conjugate block --------------------------
    #: Fixed design ``(n, Q)`` whose Gaussian coefficients are drawn together
    #: with ``(trend_coef, gamma, kappa)``; shared, masked rows zero. ``None``
    #: when the fit has none.
    extra_design: Any = field(data=0, default=None)
    extra_coef: Any = field(chains=CHAIN_AXIS, default=None)
    #: Prior standard deviation per extra coefficient, shape ``(Q,)``, shared.
    extra_prior_scale_vec: Any = field(default=None)

    @property
    def num_extra_cols(self) -> int:
        """Number of fixed extra columns in the conjugate block."""
        return 0 if self.extra_design is None else int(self.extra_design.shape[-1])

    @property
    def is_tempered(self) -> bool:
        """Whether the sweep must run against a tempered likelihood."""
        return self.tempering_levels > 1

    @property
    def has_chain_axis(self) -> bool:
        """Whether this state carries a leading chain axis."""
        return self.beta.ndim > 1

    @property
    def num_chains(self) -> int:
        """Number of chains represented by this state."""
        return self.beta.shape[0] if self.has_chain_axis else 1


def _leaf_is_chained(root: Any, path: tuple[Any, ...]) -> bool:
    """Whether the leaf at ``path`` carries a chain axis.

    Walks the path through the actual objects and reads each dataclass field's
    ``chains`` metadata; the innermost annotation wins.  A field with no
    annotation is shared, which is `bartz`'s convention (``X``, ``resid_unit``
    and ``Forest.leaf_unit`` are shared; ``resid``, ``Forest.leaf_tree`` and
    ``Wishart.value`` are not).
    """
    obj = root
    chained = False
    for entry in path:
        name = getattr(entry, "name", None)
        if name is None:  # a sequence index: keep the enclosing decision
            idx = getattr(entry, "idx", None)
            if idx is not None and isinstance(obj, (tuple, list)) and idx < len(obj):
                obj = obj[idx]
            continue
        try:
            fields_by_name = {f.name: f for f in dataclasses.fields(obj)}
        except TypeError:  # not a dataclass; nothing more to learn
            return chained
        fld = fields_by_name.get(name)
        if fld is None:
            return chained
        chained = fld.metadata.get("chains", None) == CHAIN_AXIS
        obj = getattr(obj, name)
    return chained


def chain_filter_spec(state: LongBetState) -> Any:
    """Boolean pytree marking the leaves that carry a chain axis.

    Used with :func:`equinox.partition` to split a chained state into the part
    ``jax.vmap`` should map over and the part it should close over.
    """
    leaves_with_path, treedef = jax.tree_util.tree_flatten_with_path(state)
    flags = [_leaf_is_chained(state, path) for path, _ in leaves_with_path]
    return jax.tree_util.tree_unflatten(treedef, flags)


def chained_field_names(state: LongBetState) -> set[str]:
    """Top-level field names that contain at least one chain-axis leaf."""
    names: set[str] = set()
    for path, _ in jax.tree_util.tree_flatten_with_path(state)[0]:
        if _leaf_is_chained(state, path):
            top = next((getattr(e, "name") for e in path if hasattr(e, "name")), None)
            if top is not None:
                names.add(top)
    return names


def split_chain_fields(state: LongBetState) -> tuple[Any, Any]:
    """Split into ``(per_chain, shared)`` pytrees."""
    return eqx.partition(state, chain_filter_spec(state))


def broadcast_to_chains(
    state: LongBetState,
    num_chains: int,
    key: Any = None,
) -> LongBetState:
    """Replicate the chain-axis fields of a single-chain state.

    ``X``, ``y``, the index vectors and the kernel factor are left untouched, so
    a ``k``-chain fit costs ``k`` copies of the parameters and residuals rather
    than ``k`` copies of the design matrix.

    When ``key`` is given the chains are additionally **overdispersed**: each
    one's scalar and vector parameters are drawn from their priors rather than
    started at a common point. This matters for what R-hat means. Chains started
    at identical values and separated only by their random streams can agree
    with each other while all of them sit in the same unvisited region, so R-hat
    understates non-convergence; the statistic is designed for starts that are
    dispersed relative to the target.

    Only quantities that leave the initial fit unchanged are dispersed freely --
    at initialisation both forests are empty, so ``beta``, ``b0``, ``b1`` and the
    variances do not enter the fitted values at all. ``gamma`` does, so the
    residual is corrected for it and the model's invariant
    ``resid == y - fitted`` continues to hold on the first sweep.
    """
    if num_chains <= 1:
        return state
    if state.has_chain_axis:
        raise ValueError("state already carries a chain axis")

    per, shared = split_chain_fields(state)
    per = jax.tree.map(lambda x: jnp.broadcast_to(x, (num_chains, *x.shape)), per)
    state = eqx.combine(per, shared)
    if key is None:
        return state
    return _overdisperse_chains(state, key, num_chains)


def _overdisperse_chains(
    state: LongBetState, key: Any, num_chains: int
) -> LongBetState:
    """Draw each chain's starting parameters from their priors."""
    k_beta, k_gamma, k_sig, k_siggam, k_b = jax.random.split(key, 5)

    # beta ~ N(0, K_tilde), via the same Cholesky factor the sampler uses, then
    # mapped onto the canonical half-space sum_{s>=1} beta_s >= 0 by the same
    # idempotent Z_2 quotient the sweep applies after every trajectory draw. The
    # starts therefore cover the declared prior on the identified quotient
    # instead of a narrow positive band, which is what R-hat needs; confining
    # them to the positive orthant would let chains agree for a reason that has
    # nothing to do with convergence.
    if state.sample_beta:
        eta = jax.random.normal(
            k_beta, (num_chains, state.K_chol.shape[-1]), jnp.float32
        )
        raw_beta = jnp.einsum('ij,cj->ci', state.K_chol, eta)
        orientation = jnp.where(jnp.sum(raw_beta[:, 1:], axis=-1) < 0.0, -1.0, 1.0)
        beta = orientation[:, None] * raw_beta
    else:
        beta = state.beta

    # The variances start from their priors when those are proper, and from a
    # lognormal spread around 1 under the improper reference prior, which has no
    # scale of its own to draw from.
    def _inv_gamma_or_spread(k, a, b):
        if a > 0 and b > 0:
            return (b / jnp.maximum(jax.random.gamma(k, a, (num_chains,)), 1e-30)).astype(
                jnp.float32
            )
        return jnp.exp(jax.random.normal(k, (num_chains,), jnp.float32)).astype(jnp.float32)

    if state.outcome_type_str in ("binary", "ordinal"):
        sigma2 = jnp.ones((num_chains,), dtype=jnp.float32)
    else:
        sigma2 = _inv_gamma_or_spread(k_sig, state.sigma_prior_a, state.sigma_prior_b)

    if state.random_intercept:
        sigma_gamma2 = _inv_gamma_or_spread(
            k_siggam, state.gamma_prior_a, state.gamma_prior_b
        )
        # gamma_i | sigma_gamma^2 ~ N(0, sigma_gamma^2); damped by 0.25 at Sweep 0
        # so initial random intercepts do not corrupt the first forest sweep before
        # gamma is updated at Step 6.
        gamma = 0.25 * jax.random.normal(
            k_gamma, (num_chains, state.N_units), jnp.float32
        ) * jnp.sqrt(sigma_gamma2)[:, None]
        resid = state.resid - jnp.where(
            state.obs_mask, gamma[:, state.unit_idx], 0.0
        )
    else:
        sigma_gamma2 = state.sigma_gamma2
        gamma = state.gamma
        resid = state.resid

    # The coding scales are fixed (b0 = 0, b1 = 1): nothing to disperse.
    del k_b
    b0 = state.b0
    b1 = state.b1

    # The smooth trend block is dispersed as well (damped by 0.25 at Sweep 0 so
    # the initial un-updated period wave does not overwhelm the standardized
    # outcome before trend_block_step runs at the end of Sweep 1).
    if state.use_trend_block and state.trend_design.shape[1] > 0:
        sigma_c = jnp.asarray(state.trend_prior_scale_vec, dtype=jnp.float32)
        trend_coef = 0.25 * sigma_c[None, :] * jax.random.normal(
            jax.random.fold_in(key, 9301),
            (num_chains, state.trend_design.shape[1]),
            jnp.float32,
        )
        resid = resid - jnp.where(
            state.obs_mask, trend_coef @ state.trend_design.T, 0.0
        )
    else:
        trend_coef = state.trend_coef

    # With the anchored moderator the treatment forest starts at ``nu = 1``,
    # so ``beta`` *does* enter the initial fit (``w = b_z beta_S``) and the
    # dispersed trajectory has to be taken out of the residual, exactly as
    # ``gamma`` is above. Without the anchor ``nu`` starts at zero and nothing
    # changes; the branch is skipped so that path stays bit-identical.
    if state.sample_beta and state.anchor_moderator and not state.use_lag_kernel:
        b_z_init = jnp.where(state.z_vec == 1.0, b1[:, None], b0[:, None])
        delta_w = b_z_init * (beta - state.beta)[:, state.exposure_idx]
        resid = resid - jnp.where(state.obs_mask, delta_w * state.nu_fit, 0.0)

    state = eqx.tree_at(
        lambda s: (s.beta, s.gamma, s.sigma2, s.sigma_gamma2, s.b0, s.b1,
                   s.trend_coef, s.resid),
        state,
        (beta, gamma, sigma2, sigma_gamma2, b0, b1, trend_coef, resid),
    )
    # The lag kernel starts from its prior at unit scale, mapped onto the
    # canonical half-space ``sum(kappa) >= 0`` when the moderator is not
    # anchored (the sweep's own quotient). With the anchored moderator the
    # kernel enters the initial fit through ``w = lag_design @ kappa`` with
    # ``nu = 1``: damp the perturbation around the data-scaled start so no
    # chain begins with a near-zero or negative exposure weight before the
    # first moderator draw, and when the conjugate block is active draw
    # ``(trend_coef, extra_coef, gamma, kappa)`` jointly at ``nu = 1, delta = 0``
    # before Sweep 1's moderator steps run.
    if state.use_lag_kernel:
        eta_k = jax.random.normal(
            jax.random.fold_in(key, 9401), (num_chains, state.lag_K_chol.shape[-1]),
            jnp.float32)
        raw_k = jnp.einsum('ij,cj->ci', jnp.asarray(state.lag_K_chol, jnp.float32), eta_k)
        if state.anchor_moderator:
            kappa = (state.kappa + 0.25 * raw_k).astype(jnp.float32)
        else:
            orient_k = jnp.where(jnp.sum(raw_k, axis=-1) < 0.0, -1.0, 1.0)
            kappa = (orient_k[:, None] * raw_k).astype(jnp.float32)
        delta_w = (kappa - state.kappa) @ state.lag_design.T              # (chains, n)
        resid = state.resid - jnp.where(state.obs_mask, delta_w * state.nu_fit, 0.0)
        state = eqx.tree_at(lambda s: (s.kappa, s.resid), state, (kappa, resid))
        if (state.anchor_moderator and state.use_trend_block
                and state.outcome_type_str == "continuous"
                and not state.baseline_untreated_only):
            q_tot = state.trend_design.shape[1]
            col_scale = jnp.asarray(state.trend_prior_scale_vec, jnp.float32)
            if state.num_extra_cols:
                H_blk = jnp.concatenate([state.trend_design, state.extra_design], axis=1)
                col_scale = jnp.concatenate(
                    [col_scale, jnp.asarray(state.extra_prior_scale_vec, jnp.float32)])
                coef_blk = jnp.concatenate([state.trend_coef, state.extra_coef], axis=1)
            else:
                H_blk, coef_blk = state.trend_design, state.trend_coef
            tb_keys = jax.random.split(jax.random.fold_in(key, 9402), num_chains)
            with enable_x64():
                coef_blk, gamma_w, kappa_w, resid_w = jax.vmap(
                    lambda k_c, c_c, g_c, b_c, r_c, d_c, s2_c, sg2_c: trend_block_step(
                        k_c, c_c, g_c, b_c, r_c, H_blk, d_c, state.obs_mask,
                        state.unit_idx, state.N_units, state.exposure_idx, s2_c, sg2_c,
                        state.lag_K_chol, col_scale, state.random_intercept,
                        state.sample_beta, trend_gram=state.trend_gram,
                        trend_unit_sum=state.trend_unit_sum, lag_design=state.lag_design,
                        D_lag_K_chol=state.D_lag_K_chol,
                    )
                )(tb_keys, coef_blk, state.gamma, state.kappa, state.resid,
                  state.nu_fit, state.sigma2, state.sigma_gamma2)
            if state.num_extra_cols:
                tc_w, ec_w = coef_blk[:, :q_tot], coef_blk[:, q_tot:]
                state = eqx.tree_at(
                    lambda s: (s.trend_coef, s.extra_coef, s.gamma, s.kappa, s.resid),
                    state, (tc_w, ec_w, gamma_w, kappa_w, resid_w))
            else:
                state = eqx.tree_at(
                    lambda s: (s.trend_coef, s.gamma, s.kappa, s.resid),
                    state, (coef_blk, gamma_w, kappa_w, resid_w))
    if state.num_categories > 2:
        # Ordinal-only streams; preserve the five legacy initialization keys.
        gap_key = jax.random.fold_in(key, 8101)
        latent_keys = jax.random.split(jax.random.fold_in(key, 8102), num_chains)
        gaps = jnp.diff(jnp.concatenate((jnp.zeros((num_chains, 1)),
                                         state.cutpoints), axis=1), axis=1)
        gaps *= jnp.exp(.25 * jax.random.normal(gap_key, gaps.shape, jnp.float32))
        cutpoints = jnp.cumsum(gaps, axis=1)
        mean = (state.alpha[:, None] * state.mu_fit
                + trend_coef @ state.trend_design.T
                + jnp.where(state.z_vec == 1, b1[:, None], b0[:, None])
                * beta[:, state.exposure_idx] * state.nu_fit
                + gamma[:, state.unit_idx])
        z = jax.vmap(lambda k, m, cp, old: sample_ordinal_latents(
            k, state.y, m, 1., cp, state.obs_mask, old))(
                latent_keys, mean, cutpoints, state.z)
        state = eqx.tree_at(lambda s: (s.cutpoints, s.z, s.resid), state,
                            (cutpoints, z, jnp.where(state.obs_mask, z - mean, 0.)))
    return state


@enable_x64(False)
def init_longbet(
    *,
    X_unified: UInt[Array, 'p n'],
    y: Float32[Array, ' n'],
    unit_idx: Int32[Array, ' n'],
    time_idx: Int32[Array, ' n'],
    exposure_idx: Int32[Array, ' n'],
    z_vec: Float32[Array, ' n'],
    obs_mask: Bool[Array, ' n'],
    max_split_mu: UInt[Array, ' p'],
    max_split_nu: UInt[Array, ' p'],
    config: LongBetConfig,
    offset: float = 0.0,
    trend_design: Float32[Array, 'n q_tot'] | None = None,
    trend_num_period: int = 0,
    trend_num_forecast: int = 0,
    num_chains: int | None = None,
    chain_key: Any = None,
    mesh: Any = None,
    lin_design: Float32[Array, 'n p'] | None = None,
    lag_design: Float32[Array, 'n L_plus_1'] | None = None,
    extra_design: Float32[Array, 'n Q'] | None = None,
    lag_K_chol: Any = None,
    **kwargs: Any,
) -> LongBetState:
    """Build the initial ``LongBetState``.

    Parameters
    ----------
    X_unified
        Unified binned predictor matrix, shape ``(p, n)``, shared by both
        forests. Which columns each forest may split on is controlled by its own
        ``max_split`` vector, not by a second copy of the matrix.
    y
        Response vector of length ``n`` (already standardized for a continuous
        outcome; raw 0/1 for a binary one). Masked cells must be finite.
    unit_idx, time_idx, exposure_idx, z_vec, obs_mask
        Panel indices and indicators, length ``n``.
    max_split_mu, max_split_nu
        Per-forest split bounds, length ``p``. A zero entry blocks that column.
    config
        Sampler configuration.
    offset
        Additive offset handed to `bartz`. For a binary outcome this is the
        probit intercept ``Phi^-1(rate)``; for a standardized continuous outcome
        it is 0.
    trend_design
        Conjugate calendar-time design ``[ psi(t) | phi(t) (x) b(X_i) ]``, shape
        ``(n, q_tot)``, from :func:`longbet._trend_basis.build_trend_design`.
        ``None`` (or zero columns) disables the block. It is shared across
        chains and never resampled.
    trend_num_period
        How many leading columns of ``trend_design`` are the free common-time
        block. The prior shrinks those columns loosely and the remaining
        interaction columns in proportion to their number, so getting this wrong
        misprices the prior rather than breaking the sampler.
    num_chains
        Number of chains to replicate the parameter block into.
    chain_key
        PRNG key used to overdisperse the chains' starting points. Without it
        every chain starts identically, which weakens R-hat.
    mesh
        Optional JAX device mesh.
    lin_design
        Standardized treatment covariates per cell, shape ``(n, p)``, for the
        moderator's linear part; used only when ``config.linear_trt``.
    lag_design
        Lag-dose rows per cell, shape ``(n, L + 1)``: column ``l`` holds the
        dose ``l`` periods before the cell's period (zero before the panel
        starts). Switches the exposure representation to the distributed-lag
        kernel ``w = lag_design @ kappa``; ``exposure_idx`` and ``z_vec`` are
        then inert for the sampler (pass zeros and ones).
    extra_design
        Fixed columns ``(n, Q)`` whose Gaussian coefficients are drawn in the
        conjugate block together with the trend, the unit intercepts and the
        kernel. The block runs with these columns alone if the calendar-time
        design is empty.
    lag_K_chol
        Optional factor ``L_kappa`` of the kernel prior at unit scale, shape
        ``(L + 1, r)`` (``kappa = s L_kappa u``, ``u ~ N(0, I_r)``), replacing
        the one built from ``config``. Used by the exact-nesting tests.
    """
    y = jnp.asarray(y, dtype=jnp.float32)
    offset = float(offset)
    unit_idx = jnp.asarray(unit_idx, dtype=jnp.int32)
    time_idx = jnp.asarray(time_idx, dtype=jnp.int32)
    exposure_idx = jnp.asarray(exposure_idx, dtype=jnp.int32)
    z_vec = jnp.asarray(z_vec, dtype=jnp.float32)
    obs_mask = jnp.asarray(obs_mask, dtype=jnp.bool_)

    n = y.shape[0]
    N_units = int(jnp.max(unit_idx)) + 1
    T_periods = int(jnp.max(time_idx)) + 1
    S_max = int(jnp.max(exposure_idx))

    p_nonterminal_mu = make_p_nonterminal(
        config.max_depth_pr, config.alpha_split_pr, config.beta_split_pr
    )
    p_nonterminal_nu = make_p_nonterminal(
        config.max_depth_trt, config.alpha_split_trt, config.beta_split_trt
    )

    filter_splitless_mu = int(jnp.sum(max_split_mu == 0))
    filter_splitless_nu = int(jnp.sum(max_split_nu == 0))

    # Standard BART/BCF leaf-variance priors on a standardized response:
    # prognostic forest has total prior SD 1.5, moderating forest has total prior SD 0.5.
    sigma_mu = 3.0 / (2.0 * math.sqrt(config.num_trees_pr))
    sigma_nu = 0.5 / math.sqrt(config.num_trees_trt)
    leaf_prior_cov_inv_mu = jnp.array(1.0 / sigma_mu**2, dtype=jnp.float32)
    leaf_prior_cov_inv_nu = jnp.array(1.0 / sigma_nu**2, dtype=jnp.float32)

    # Treatment-only coding, fixed: the treatment term beta_S * nu enters
    # treated cells only, the untreated surface is mu + gamma, and beta_0 is
    # prior-only. (Adaptive coding scales were removed: they add an untreated
    # term b0 * beta_0 * nu(x, t) that competes with mu on every untreated cell
    # and enters the reported effect, a slow, effect-relevant ridge.)
    b0_init, b1_init = 0.0, 1.0

    beta_init = jnp.ones(S_max + 1, dtype=jnp.float32)
    b_z_init = jnp.where(z_vec == 1.0, b1_init, b0_init)

    # --- distributed-lag exposure kernel: validate, factor the prior, and
    # choose a start. The kernel starts flat at the level that gives the
    # average dosed cell a unit weight -- the analogue of ``beta = 1`` on the
    # absorbing clock -- so the first sweep's moderator step sees the dosed
    # cells; a panel with no dose at all starts at ``kappa = 1``.
    use_lag_kernel = lag_design is not None
    lag_fields: dict[str, Any] = {}
    if use_lag_kernel:
        if config.outcome != "continuous":
            raise ValueError("exposure_lags (the distributed-lag kernel) supports a continuous "
                             f"outcome only, got outcome={config.outcome!r}")
        if config.baseline_cells != "all":
            raise ValueError("exposure_lags cannot be combined with baseline_cells='untreated'")
        if config.tempering_levels > 1:
            raise ValueError("exposure_lags cannot be combined with parallel tempering")
        lag_design = jnp.asarray(lag_design, dtype=jnp.float32)
        if lag_design.ndim != 2 or lag_design.shape[0] != n or lag_design.shape[1] < 1:
            raise ValueError(f"lag_design must be (n, L + 1) with n={n}, got {tuple(lag_design.shape)}")
        lag_np = np.asarray(lag_design, dtype=np.float64)
        obs_np = np.asarray(obs_mask)
        if not np.all(np.isfinite(lag_np[obs_np])):
            raise ValueError("exposure_lags must be finite on every observed cell")
        lag_design = jnp.where(obs_mask[:, None], lag_design, 0.0)
        L1 = int(lag_design.shape[1])
        if lag_K_chol is None:
            K_lag = build_kernel_matrix(
                s=np.arange(L1, dtype=np.float64),
                sig_knl=config.lag_sig_knl,
                lambda_knl=config.lag_lambda_knl,
                kernel_type=config.lag_kernel_type,
                sigma_m=config.lag_sigma_m,
                gp_constant_mean=config.lag_gp_constant_mean,
            )
            L_lag = kernel_cholesky(K_lag)
            L_lag_sq = L_lag
        else:
            L_lag = np.asarray(lag_K_chol, dtype=np.float64)
            if L_lag.ndim != 2 or L_lag.shape[0] != L1 or L_lag.shape[1] < L1:
                raise ValueError(f"lag_K_chol must be (L + 1, r) with r >= L + 1 = {L1}, "
                                 f"got {tuple(L_lag.shape)}")
            L_lag_sq = L_lag if (L_lag.shape[1] == L1 and np.allclose(np.triu(L_lag, 1), 0.0)
                                 and np.all(np.diag(L_lag) > 0)) else kernel_cholesky(L_lag @ L_lag.T)
        # Start: kappa_0 = 1 / (mean |current dose| over the currently dosed
        # cells), kappa_{l >= 1} = 0, so the average currently dosed cell has
        # unit weight -- for a binary absorbing dose this is exactly the
        # clock's ``beta = 1`` start (w = 1 on treated cells). Signed lag
        # inputs (propensity-centered T0 or pre-offer demeaned doses) are
        # supported. A panel with no current dose anywhere falls back to the
        # same level over the lag sums, and one with no dose at all to 1.
        cur = np.where(obs_np, np.abs(lag_np[:, 0]), 0.0)
        if np.any(cur > 0.0):
            kappa_init_np = np.zeros(L1)
            kappa_init_np[0] = 1.0 / float(cur[cur > 0.0].mean())
        else:
            row_dose = np.where(obs_np, np.abs(lag_np).sum(axis=1), 0.0)
            dosed = row_dose > 0.0
            kappa_init_np = np.full(L1, 1.0 / float(row_dose[dosed].mean()) if dosed.any() else 1.0)
        kappa_init = jnp.asarray(kappa_init_np, dtype=jnp.float32)
        lag_scale_prior = float(config.lag_scale_prior)
        if lag_K_chol is None or L_lag.dtype != np.float64:
            lag_K_chol_arr = jnp.asarray(L_lag, dtype=jnp.float32)
        else:
            with enable_x64():
                lag_K_chol_arr = jnp.asarray(L_lag, dtype=jnp.float64)
        with enable_x64():
            d_lag_k_chol_arr = jnp.asarray(
                np.asarray(lag_design, dtype=np.float64)
                @ np.asarray(lag_K_chol_arr, dtype=np.float64),
                dtype=jnp.float64,
            )
        lag_fields = dict(
            use_lag_kernel=True, lag_design=lag_design, kappa=kappa_init,
            # The built-in factor is float32 like ``K_chol``; a supplied one is
            # kept at its own precision (float64 keeps an exact-nesting factor
            # exact inside the float64 block).
            lag_K_chol=lag_K_chol_arr,
            D_lag_K_chol=d_lag_k_chol_arr,
            lag_K_chol_sq=jnp.asarray(L_lag_sq, dtype=jnp.float32),
            lag_learn_scale=bool(config.lag_learn_scale),
            lag_scale2=jnp.array(1.0, dtype=jnp.float32),
            lag_scale_aux=jnp.array(1.0 / lag_scale_prior ** 2, dtype=jnp.float32),
            lag_scale_prior=lag_scale_prior,
        )
        w_init = lag_design @ kappa_init
    else:
        w_init = b_z_init * beta_init[exposure_idx]

    mask_nu_init = obs_mask & (jnp.abs(w_init) > 1e-6)
    error_scale_nu = jnp.reciprocal(jnp.where(mask_nu_init, jnp.abs(w_init), 1.0))

    is_binary = config.outcome == "binary" or (
        config.outcome == "ordinal" and config.num_categories == 2)
    is_ordered = config.outcome == "ordinal" and config.num_categories > 2
    cutpoints = jnp.empty((0,), jnp.float32)
    working_y = y
    if config.outcome == "ordinal":
        prepared = prepare_ordinal(np.where(np.asarray(obs_mask), np.asarray(y), np.nan),
                                   config.num_categories)
        y = jnp.asarray(prepared.labels)
        cutpoints = jnp.asarray(prepared.cutpoints)
        if is_ordered:
            full = full_cutpoints(cutpoints)
            labels = y.astype(jnp.int32)
            lo, hi = full[labels], full[labels + 1]
            # A finite interior working response, never category labels passed
            # through bartz's binary branch. No hidden initialization RNG.
            lo_f = jnp.where(jnp.isfinite(lo), lo, hi - 2.)
            hi_f = jnp.where(jnp.isfinite(hi), hi, lo + 2.)
            working_y = jnp.where(obs_mask, lo_f / 2 + hi_f / 2, 0.)
        else:
            working_y = y

    def _err_cov() -> Wishart:
        # A fresh object per call: bartz's init donates and deletes its buffers.
        return Wishart(nu=jnp.array(3.0), rate=jnp.array(3.0), value=jnp.array(1.0))

    # X and y are copied because bartz's init donates and deletes its buffers.
    state_mu = init(
        X=jnp.copy(X_unified),
        y=jnp.copy(working_y),
        outcome_type="binary" if is_binary else "continuous",
        offset=offset,
        max_split=jnp.copy(max_split_mu),
        num_trees=config.num_trees_pr,
        p_nonterminal=p_nonterminal_mu,
        leaf_prior_cov_inv=leaf_prior_cov_inv_mu,
        error_cov_inv=None if is_binary else _err_cov(),
        missing=~obs_mask,
        min_points_per_leaf=config.min_points_per_leaf_pr,
        filter_splitless_vars=filter_splitless_mu,
        num_chains=None,
        mesh=mesh,
        **kwargs,
    )

    state_nu = init(
        X=jnp.copy(X_unified),
        y=jnp.copy(working_y),
        outcome_type="continuous",
        offset=0.0,
        max_split=jnp.copy(max_split_nu),
        num_trees=config.num_trees_trt,
        p_nonterminal=p_nonterminal_nu,
        leaf_prior_cov_inv=leaf_prior_cov_inv_nu,
        error_cov_inv=_err_cov(),
        error_scale=error_scale_nu,
        missing=~mask_nu_init,
        min_points_per_leaf=config.min_points_per_leaf_trt,
        filter_splitless_vars=filter_splitless_nu,
        num_chains=None,
        mesh=mesh,
        **kwargs,
    )

    # Disable bartz's own error-variance draw on *both* views: sigma^2 is drawn
    # once per sweep from the full-model residual. Leaving it on either view
    # would draw it twice from two different conditionals.
    def _freeze_variance(st: State) -> State:
        frozen = eqx.tree_at(
            lambda w: (w.nu, w.rate),
            st.error_cov_inv,
            (None, None),
            is_leaf=lambda x: x is None,
        )
        return eqx.tree_at(lambda s: s.error_cov_inv, st, frozen)

    state_mu = _freeze_variance(state_mu)
    state_nu = _freeze_variance(state_nu)

    # GP kernel: build and factor in float64, hand only the Cholesky to JAX.
    K_tilde = build_kernel_matrix(
        s=np.arange(S_max + 1, dtype=np.float64),
        sig_knl=config.sig_knl,
        lambda_knl=config.lambda_knl,
        kernel_type=config.kernel_type,
        sigma_m=config.sigma_m,
        gp_constant_mean=config.gp_constant_mean,
    )
    K_chol = jnp.asarray(kernel_cholesky(K_tilde), dtype=jnp.float32)

    init_sigma_gamma2 = float(config.gamma_prior_b / max(config.gamma_prior_a, 1e-4))
    unit_counts = (
        jnp.zeros((N_units,), dtype=jnp.float32)
        .at[unit_idx]
        .add(obs_mask.astype(jnp.float32))
    )

    # Baseline module support. Under baseline_cells='untreated' the prognostic
    # forest, the trend block, the unit intercepts and sigma^2 see only the
    # observed untreated cells (z = 0); the treated cells are, to them, missing.
    # Under 'all' the two masks are the same object, so the default sweep
    # traces to exactly the same program as before the option existed.
    baseline_untreated_only = config.baseline_cells == "untreated"
    if baseline_untreated_only:
        base_mask = obs_mask & (z_vec == 0.0)
        unit_counts_base = (
            jnp.zeros((N_units,), dtype=jnp.float32)
            .at[unit_idx]
            .add(base_mask.astype(jnp.float32))
        )
    else:
        base_mask = obs_mask
        unit_counts_base = unit_counts

    # Smooth calendar-time block. Zero columns is the disabled case and is kept
    # as a real (n, 0) array rather than None so that every pytree operation --
    # partition, vmap, trace -- sees the same structure either way.
    if trend_design is None:
        trend_design = jnp.zeros((n, 0), dtype=jnp.float32)
    trend_design = jnp.asarray(trend_design, dtype=jnp.float32)
    if trend_design.shape[0] != n:
        raise ValueError(
            f"trend_design has {trend_design.shape[0]} rows, expected {n}"
        )
    # A masked cell contributes nothing to any conditional; zeroing its row here
    # means the block's sufficient statistics never have to re-apply the mask.
    trend_design = jnp.where(obs_mask[:, None], trend_design, 0.0)
    q_tot = int(trend_design.shape[1])

    # --- fixed extra columns of the conjugate block --------------------------
    # Validated and whitened here; the block's cross-products below are taken
    # over the concatenated design ``[trend | extra]`` once, so the per-sweep
    # cost of the extra columns is their share of a (q_tot + Q)-dimensional
    # solve and nothing else.
    extra_fields: dict[str, Any] = {}
    num_extra = 0
    if extra_design is not None and int(np.shape(extra_design)[-1]) > 0:
        extra_design = jnp.asarray(extra_design, dtype=jnp.float32)
        if extra_design.ndim != 2 or extra_design.shape[0] != n:
            raise ValueError(f"extra_design must be (n, Q) with n={n}, got {tuple(extra_design.shape)}")
        if config.outcome != "continuous":
            raise ValueError("extra_cols supports a continuous outcome only, "
                             f"got outcome={config.outcome!r}")
        if config.baseline_cells != "all":
            raise ValueError("extra_cols cannot be combined with baseline_cells='untreated'")
        if config.tempering_levels > 1:
            raise ValueError("extra_cols cannot be combined with parallel tempering")
        extra_np = np.asarray(extra_design, dtype=np.float64)
        obs_np_e = np.asarray(obs_mask)
        if not np.all(np.isfinite(extra_np[obs_np_e])):
            raise ValueError("extra_cols must be finite on every observed cell")
        rms = np.sqrt(np.mean(np.square(extra_np[obs_np_e]), axis=0))
        if np.any(rms <= 1e-12):
            dead = np.flatnonzero(rms <= 1e-12).tolist()
            raise ValueError(f"extra_cols column(s) {dead} are identically zero on the observed "
                             "cells; drop them (their coefficients would be prior-only)")
        extra_design = jnp.where(obs_mask[:, None], extra_design, 0.0)
        num_extra = int(extra_design.shape[1])
        extra_fields = dict(
            extra_design=extra_design,
            extra_coef=jnp.zeros((num_extra,), dtype=jnp.float32),
            extra_prior_scale_vec=jnp.asarray(
                float(config.extra_cols_prior_scale) / rms, dtype=jnp.float32),
        )
    use_trend_block = (bool(config.use_trend_block) and q_tot > 0) or num_extra > 0
    # The design is laid out as [ period contrasts | forecast | time x feature ],
    # so the prior splits at n_fixed = n_period + n_forecast. The blocks are
    # shrunk differently: see trend_prior_scales.
    n_period = min(int(trend_num_period), q_tot)
    n_forecast = min(max(int(trend_num_forecast), 0), max(q_tot - n_period, 0))
    n_fixed = n_period + n_forecast
    n_inter = q_tot - n_fixed
    if q_tot:
        trend_prior_scale_vec = jnp.asarray(
            trend_prior_scales(
                n_period,
                n_inter,
                config.trend_prior_scale,
                config.trend_period_scale,
                num_forecast=n_forecast,
                forecast_scale=config.trend_forecast_scale,
            ),
            dtype=jnp.float32,
        )
    else:
        trend_prior_scale_vec = jnp.zeros((0,), dtype=jnp.float32)
    use_trend_horseshoe = (
        bool(config.trend_horseshoe) and use_trend_block and n_inter > 0
    )
    # Started at the ridge, so the horseshoe has to earn any departure from it
    # out of the data rather than out of its initialization.
    hs_global_scale = float(config.trend_prior_scale) / max(n_inter, 1) ** 0.5
    hs_init = init_horseshoe(n_inter, hs_global_scale)
    trend_local2 = jnp.asarray(hs_init["local2"], dtype=jnp.float32)
    trend_local_aux = jnp.asarray(hs_init["local_aux"], dtype=jnp.float32)
    trend_global2 = jnp.asarray(hs_init["global2"], dtype=jnp.float32)
    trend_global_aux = jnp.asarray(hs_init["global_aux"], dtype=jnp.float32)

    # The initial full-model residual is whatever bartz's own init left on the
    # prognostic view, read back into data units. For a standardized continuous
    # outcome that is y - offset; for a binary outcome the latent z starts at
    # the offset and the residual starts at zero. Deriving it rather than
    # writing y - offset keeps the binary path correct.
    resid_init = state_mu.resid * state_mu.resid_unit
    if config.anchor_moderator:
        m_nu_init = state_nu.forest.leaf_tree.shape[0]
        forest_nu_init = eqx.tree_at(
            lambda f: f.leaf_tree,
            state_nu.forest,
            state_nu.forest.leaf_tree + 1.0 / (m_nu_init * state_nu.forest.leaf_unit),
        )
        nu_fit_init = jnp.ones(n, dtype=jnp.float32)
        resid_init = resid_init - w_init
    else:
        forest_nu_init = state_nu.forest
        nu_fit_init = jnp.zeros(n, dtype=jnp.float32)
    resid_init = jnp.where(obs_mask, resid_init, 0.0).astype(jnp.float32)

    # The block's cross-products are accumulated over the cells that inform the
    # baseline module; under the cut that excludes the treated cells. With
    # extra columns they are taken over the concatenated design the block
    # actually solves for.
    block_design = (jnp.concatenate([trend_design, extra_fields["extra_design"]], axis=1)
                    if num_extra else trend_design)
    trend_gram, trend_unit_sum = (
        trend_sufficient_statistics(block_design, base_mask, unit_idx, N_units)
        if use_trend_block else (None, None)
    )

    def _predictors_unit_constant(
        X_arr: np.ndarray, unit_arr: np.ndarray, max_split_arr: np.ndarray, n_u: int
    ) -> bool:
        active = np.flatnonzero(np.asarray(max_split_arr) > 0)
        if active.size == 0:
            return True
        for col in active:
            vals = np.asarray(X_arr[col])
            u_min = np.full(n_u, np.iinfo(vals.dtype).max, dtype=vals.dtype)
            u_max = np.zeros(n_u, dtype=vals.dtype)
            np.minimum.at(u_min, unit_arr, vals)
            np.maximum.at(u_max, unit_arr, vals)
            seen = u_min <= u_max
            if not np.all(u_min[seen] == u_max[seen]):
                return False
        return True

    X_np_init = np.asarray(X_unified)
    unit_idx_np_init = np.asarray(unit_idx)
    unit_counts_np_init = np.asarray(unit_counts)
    has_unbalanced_units = bool(
        unit_counts_np_init.size > 0
        and not np.all(unit_counts_np_init == unit_counts_np_init[0])
    )
    collapse_gamma_mu = bool(
        config.random_intercept
        and not config.split_calendar_mu
        and _predictors_unit_constant(X_np_init, unit_idx_np_init, np.asarray(max_split_mu), N_units)
    )
    collapse_gamma_nu = bool(
        config.random_intercept
        and not config.split_calendar_trt
        and not config.split_exposure_trt
        and _predictors_unit_constant(X_np_init, unit_idx_np_init, np.asarray(max_split_nu), N_units)
    )

    use_linear_trt = bool(config.linear_trt) and lin_design is not None and np.shape(lin_design)[1] > 0
    if use_linear_trt:
        lin_design = jnp.where(obs_mask[:, None], jnp.asarray(lin_design, dtype=jnp.float32), 0.0)
        p_lin = int(lin_design.shape[1])
        lin_scale = float(config.linear_trt_scale) / math.sqrt(p_lin)
        lin_fields = dict(
            use_linear_trt=True, lin_design=lin_design,
            delta=jnp.zeros(p_lin, dtype=jnp.float32),
            lin_fit=jnp.zeros(n, dtype=jnp.float32),
            tau_delta2=jnp.array(lin_scale ** 2, dtype=jnp.float32),
            tau_delta_aux=jnp.array(1.0 / lin_scale ** 2, dtype=jnp.float32),
            linear_trt_global_scale=lin_scale,
            lin_local2=jnp.ones(p_lin, dtype=jnp.float32),
            lin_local_aux=jnp.ones(p_lin, dtype=jnp.float32),
        )
    else:
        lin_fields = {}
    ar1_fields = {}
    if config.ar1_errors:
        tidx = np.asarray(time_idx)
        if not (config.random_intercept and collapse_gamma_mu and collapse_gamma_nu):
            raise ValueError(
                "ar1_errors needs random_intercept=True and unit-level partitions (covariates "
                "constant within unit, no calendar or exposure splits in either forest)")
        if has_unbalanced_units or not bool(np.all(np.asarray(obs_mask))):
            raise ValueError("ar1_errors needs a balanced panel with every cell observed")
        if T_periods < 2 or not np.array_equal(tidx, np.tile(np.arange(T_periods), N_units)):
            raise ValueError("ar1_errors needs a unit-major panel with at least two periods")
        ar1_fields = dict(
            ar1=True, ar1_phi=jnp.array(0.0, dtype=jnp.float32),
            ar1_first=jnp.asarray(tidx == 0), ar1_phi_max=float(config.ar1_phi_max),
            ar1_accepts=jnp.array(0.0, dtype=jnp.float32),
        )
        if use_trend_block:
            from longbet._ar1 import trend_ar1_statistics
            ar1_fields["ar1_trend_stats"] = trend_ar1_statistics(
                block_design, obs_mask, jnp.asarray(tidx == 0), unit_idx, N_units)
    if config.learn_beta_scale:
        scale_fields = dict(
            learn_beta_scale=True,
            beta_scale2=jnp.array(1.0, dtype=jnp.float32),
            beta_scale_aux=jnp.array(1.0 / float(config.beta_scale_prior) ** 2, dtype=jnp.float32),
            beta_scale_prior=float(config.beta_scale_prior),
        )
    else:
        scale_fields = {}

    state = LongBetState(
        _chain_anchor=state_mu._chain_anchor,
        X=state_mu.X,
        y=y if is_ordered else state_mu.y,
        z=working_y if is_ordered else state_mu.z,
        binary_indices=state_mu.binary_indices,
        resid=resid_init,
        resid_unit=state_mu.resid_unit,
        resid_eff_scale=state_mu.resid_eff_scale,
        resid_inexact_integral=state_mu.resid_inexact_integral,
        error_cov_inv=state_mu.error_cov_inv,
        error_scale=state_mu.error_scale,
        prec_scale=state_mu.prec_scale,
        inv_sdev_scale=state_mu.inv_sdev_scale,
        inv_sdev_unit=state_mu.inv_sdev_unit,
        n_non_missing=state_mu.n_non_missing,
        sum_diag_prec_scale=state_mu.sum_diag_prec_scale,
        forest=state_mu.forest,
        config=state_mu.config,
        forest_nu=forest_nu_init,
        resid_nu=state_nu.resid,
        prec_scale_nu=state_nu.prec_scale,
        inv_sdev_scale_nu=state_nu.inv_sdev_scale,
        inv_sdev_unit_nu=state_nu.inv_sdev_unit,
        resid_unit_nu=state_nu.resid_unit,
        leaf_prior_cov_inv_nu=state_nu.forest.leaf_prior_cov_inv,
        unit_idx=unit_idx,
        time_idx=time_idx,
        exposure_idx=exposure_idx,
        z_vec=z_vec,
        obs_mask=obs_mask,
        unit_counts=unit_counts,
        N_units=N_units,
        T_periods=T_periods,
        S_max=S_max,
        beta=beta_init,
        gamma=jnp.zeros(N_units, dtype=jnp.float32),
        alpha=jnp.array(1.0, dtype=jnp.float32),
        b0=jnp.array(b0_init, dtype=jnp.float32),
        b1=jnp.array(b1_init, dtype=jnp.float32),
        sigma2=jnp.array(1.0, dtype=jnp.float32),
        sigma_gamma2=jnp.array(init_sigma_gamma2, dtype=jnp.float32),
        cutpoints=cutpoints,
        # Include the forest offset: prediction includes it, and the alpha
        # conditional must scale the same prognostic mean as the likelihood.
        mu_fit=jnp.full(n, offset, dtype=jnp.float32),
        nu_fit=nu_fit_init,
        trend_design=trend_design,
        trend_coef=jnp.zeros(q_tot, dtype=jnp.float32),
        K_chol=K_chol,
        resid_eff_scale_nu=state_nu.resid_eff_scale,
        resid_inexact_integral_nu=state_nu.resid_inexact_integral,
        sample_beta=config.sample_beta,
        random_intercept=config.random_intercept,
        temperature=jnp.array(1.0, dtype=jnp.float32),
        tempering_levels=int(config.tempering_levels),
        tempering_beta_min=float(config.tempering_beta_min),
        # The transfer moves shift fitted mass between mu and (gamma, beta, nu)
        # using every cell; under the cut that would re-couple the modules.
        use_inter_ensemble_move=bool(config.use_inter_ensemble_move) and not baseline_untreated_only,
        use_trend_block=use_trend_block,
        trend_prior_scale_vec=trend_prior_scale_vec,
        trend_num_period=n_fixed,
        use_trend_horseshoe=use_trend_horseshoe,
        trend_hs_global_scale=hs_global_scale,
        trend_local2=trend_local2,
        trend_local_aux=trend_local_aux,
        trend_global2=trend_global2,
        trend_global_aux=trend_global_aux,
        trend_gram=trend_gram,
        trend_unit_sum=trend_unit_sum,
        collapse_gamma_mu=collapse_gamma_mu,
        collapse_gamma_nu=collapse_gamma_nu,
        has_unbalanced_units=has_unbalanced_units,
        base_mask=base_mask,
        unit_counts_base=unit_counts_base,
        baseline_untreated_only=baseline_untreated_only,
        gamma_prior_a=config.gamma_prior_a,
        gamma_prior_b=config.gamma_prior_b,
        sigma_prior_a=config.sigma_prior_a,
        sigma_prior_b=config.sigma_prior_b,
        outcome_type_str=config.outcome,
        num_categories=(config.num_categories if config.outcome == "ordinal"
                        else 2 if is_binary else 0),
        cutpoint_prior_scale=config.cutpoint_prior_scale,
        ridge_accepts=jnp.array(0.0, dtype=jnp.float32),
        anchor_moderator=bool(config.anchor_moderator),
        **(dict(interweave_gamma=True, interweave_accepts=jnp.array(0.0, dtype=jnp.float32))
           if (config.interweave_gamma and config.random_intercept) else {}),
        **lin_fields,
        **scale_fields,
        **ar1_fields,
        **lag_fields,
        **extra_fields,
    )

    if is_ordered:
        latent = working_y
        if chain_key is not None and (num_chains is None or num_chains <= 1):
            latent = sample_ordinal_latents(jax.random.fold_in(chain_key, 8103),
                y, state.mu_fit, 1., cutpoints, obs_mask, working_y)
        state = eqx.tree_at(lambda s: (s.z, s.resid), state,
                            (latent, jnp.where(obs_mask, latent - state.mu_fit, 0.)))

    if num_chains is not None and num_chains > 1:
        state = broadcast_to_chains(state, int(num_chains), key=chain_key)
    if config.tempering_levels > 1 and state.has_chain_axis:
        # A single-chain state (e.g. a child of the coupled sampler) gets its
        # ladder once the caller has broadcast it to replicas.
        from longbet._tempering import ladder_temperatures
        if state.num_chains % config.tempering_levels:
            raise ValueError("tempering needs num_chains * tempering_levels replicas")
        state = eqx.tree_at(lambda s: s.temperature, state,
                            ladder_temperatures(state.num_chains, config.tempering_levels,
                                                config.tempering_beta_min))
    return state
