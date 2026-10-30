# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Forecast-Augmented Rank-K Exposure-Basis Bayesian Causal Forest for longitudinal experiments.

The estimator of the companion paper for randomized (or stratified) staggered
rollouts with heterogeneous dynamic effects. Post-adoption outcomes are
projected, within unit, onto an orthonormal exposure basis whose first mode is
the constant, so mode 1 is exactly the unit's time-averaged effect and every
higher mode is a dynamic contrast that integrates to zero.

One estimator, two uncertainty variants, selected by
`RankKLongBetConfig.uncertainty_propagation` (the engine follows from it):

* "plug-in" (the default; `uncertainty_propagation=False`, engine
  `"stochtree"`): one ``stochtree`` BCF per mode on a plug-in OLS-ridge
  forecast slope fitted on never-treated units, without forecast draws, factor
  loadings or spatial field. The default for unit-level effects, targeting and
  the ATT under unit-level randomization or covariate-stratified rollouts: at
  the validated budgets it recovered unit-level effects better than the
  propagated variant, and both variants met every coverage rule of the paper's
  second frozen confirmation.
* "propagated" (`uncertainty_propagation=True`, engine `"jax"`): each mode is
  fitted by a JAX/``bartz`` sampler under the control-only cut
  (`baseline_cells="untreated"`) that updates the prognostic forest `mu_k`, the
  treatment forest `tau_k`, the CUPED forecast slope `theta_k`, the pre-anchor
  factor loading shocks `f_k` and the Matérn-3/2 spatial GP / cluster random
  field `zeta_k` (kriged from control to treated units before the `tau_k`
  update), and pools Module 1 posterior forecast draws through mode 1 into a
  modular (cut) posterior. The conservative option when assignment or shocks
  may be spatially or regionally clustered (geographically cluster-randomized
  designs, regional rollouts), at a cost in unit-level accuracy. The sampler is
  an approximate partially collapsed Gibbs sampler, not an exact posterior
  sampler: the two forest steps use diagonal working likelihoods (see
  :mod:`longbet._rank_k_jax`).

The defaults of :class:`RankKLongBetConfig` are, for each variant, the budget
validated for it in the companion paper. Unit-level (CATE) intervals are
approximate and below nominal in the paper's benchmarks; treat them as
descriptive. Invalid inputs raise ``ValueError`` before any computation, and
:func:`fit_rank_k_longbet` emits a :class:`RankKConfigWarning` only for explicit
inputs that it would silently ignore or that contradict each other, for the
deprecated ``warm_burnin`` and for cohorts too short for ``rank_k`` modes.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field as dc_field, fields as dc_fields, replace as dc_replace
import math
import os
import sys
import time
from typing import Any, NoReturn
import warnings

import numpy as np

try:
    from longbet._forecast import (
        _cross_fit_control_adapter,
        _forecast_draws_posterior,
        build_causal_features,
        extract_pre_anchor_loadings,
    )
    from longbet._latecontext import resolve_late_context
except ImportError:
    from forecast_v2 import (  # type: ignore[no-redef]
        _cross_fit_control_adapter,
        _forecast_draws_posterior,
        build_causal_features,
        extract_pre_anchor_loadings,
    )
    from longbet._latecontext import resolve_late_context

from longbet._rank_k_event_time import (
    EventTimeLayout,
    build_event_time_layout,
    check_observation_structure,
    event_time_precision_weights,
    mode1_dynamic_correction,
    project_event_time_mode,
    restricted_wls_precision_weights,
    restricted_wls_projection,
)


def _get_rhat_ess():
    try:
        from longbet._diagnostics import compute_ess, compute_rhat
        return compute_rhat, compute_ess
    except ImportError:
        import arviz as az

        def _rhat(x: np.ndarray) -> float:
            arr = np.asarray(x, dtype=np.float64)
            if arr.ndim == 1 or arr.shape[0] < 2 or arr.shape[1] < 4:
                return float("nan")
            if float(np.std(arr)) < 1e-15:
                return 1.0
            return float(az.rhat(az.convert_to_dataset(arr[None, ...])).to_array().values.item())

        def _ess(x: np.ndarray, method: str = "bulk") -> float:
            arr = np.asarray(x, dtype=np.float64)
            if arr.ndim == 1:
                arr = arr[None, :]
            if arr.shape[1] < 4 or float(np.std(arr)) < 1e-15:
                return float(arr.size)
            return float(az.ess(az.convert_to_dataset(arr[None, ...]), method=method).to_array().values.item())

        return _rhat, _ess


def _get_rank_k_jax():
    try:
        from longbet._rank_k_jax import build_spatial_kernel_cache, sample_rank_k_mode_jax
        return build_spatial_kernel_cache, sample_rank_k_mode_jax
    except ImportError:
        from rank_k_jax import build_spatial_kernel_cache, sample_rank_k_mode_jax  # type: ignore[no-redef]
        return build_spatial_kernel_cache, sample_rank_k_mode_jax


class RankKConfigWarning(UserWarning):
    """An explicit Rank-K input or option is ignored by the fit or contradicts another one.

    Emitted by :func:`fit_rank_k_longbet` (and :func:`rank_k_placebo_audit`) when
    an argument or a non-default :class:`RankKLongBetConfig` field would be
    silently ignored (for example ``coords`` or ``n_forecast_draws`` with the
    plug-in variant, ``cluster_ids`` next to ``coords``, or ``num_mcmc``
    while forecast draws are propagated), when two settings contradict each
    other (``uncertainty_propagation=True`` with ``n_forecast_draws=1``), when
    the deprecated ``warm_burnin`` is set, and when cohorts are too short for
    ``rank_k`` modes. Choosing a supported variant or option (for example
    ``uncertainty_propagation=False``) never warns; contradictory variant
    settings (``engine="jax"`` with ``uncertainty_propagation=False`` or
    ``engine="stochtree"`` with ``uncertainty_propagation=True``) raise
    ``ValueError`` instead. Each message names the affected inputs and how to
    resolve the conflict; silence deliberate cases with
    ``warnings.filterwarnings("ignore", category=longbet.RankKConfigWarning)``.
    """


_BASIS_FAMILIES = ("exponential", "polynomial", "cosine", "pca_empirical", "custom")


@dataclass(frozen=True)
class ExposureBasisSpec:
    """Generators of the orthonormal exposure basis (Issue 4 of the methodology review).

    Every family keeps ``psi_1(s) == 1`` and is orthonormalized by modified
    Gram--Schmidt under ``<u, v>_L = (1 / L) sum_{s=1}^L u(s) v(s)``, so the ATT
    and the unit-averaged effects are invariant to the family and to ``K``
    (Proposition 2 of the paper) and every higher mode sums to zero.

    Attributes:
      family: ``"exponential"`` (the paper's generators ``{1, 1 - exp(-s / tau_c),
        s / L}`` plus cosine harmonics for ``k >= 4``; the default ``tau_c = 3.0``
        reproduces the paper bit-for-bit), ``"polynomial"`` (discrete Legendre
        polynomials on the window), ``"cosine"`` (DCT-II harmonics),
        ``"pca_empirical"`` (leading singular vectors of the never-treated units'
        double-demeaned post-window trajectories; no treated outcome enters), or
        ``"custom"`` (``custom_generators(s, L)`` returning ``(K - 1, L)`` rows).
      tau_c: Time constant (periods) of the exponential saturation generator.
      custom_generators: Callable for the ``"custom"`` family.
    """

    family: str = "exponential"
    tau_c: float = 3.0
    custom_generators: Callable[[np.ndarray, int], np.ndarray] | None = None

    def __post_init__(self) -> None:
        if self.family not in _BASIS_FAMILIES:
            raise ValueError(f"ExposureBasisSpec.family must be one of {_BASIS_FAMILIES}; got {self.family!r}")
        if not (isinstance(self.tau_c, (int, float)) and math.isfinite(float(self.tau_c)) and float(self.tau_c) > 0):
            raise ValueError(f"ExposureBasisSpec.tau_c must be a positive finite number; got {self.tau_c!r}")
        if self.family == "custom" and self.custom_generators is None:
            raise ValueError("ExposureBasisSpec(family='custom') needs custom_generators")

    @property
    def is_default(self) -> bool:
        return self.family == "exponential" and float(self.tau_c) == 3.0


DEFAULT_BASIS_SPEC = ExposureBasisSpec()
_SHORT_COHORT_POLICIES = ("drop_mode", "reduce_k", "error")
_TAU_PRIOR_FAMILIES = ("normal", "halfstudentt", "gamma_spike")
_ATT_WEIGHTINGS = ("cell", "unit")
_PROJECTIONS = ("two_pass", "restricted_wls")
_MODE_COUPLINGS = ("propagated", "plug_in", "none")
# Config fields that only the common-basis (masked / event-time) path reads.
_EVENT_TIME_ONLY_OPTIONS = ("projection", "min_support_modes", "mode_coupling", "support_covariate")


# Validated budget of each uncertainty variant: RankKLongBetConfig.__post_init__
# fills these in for the variant-dependent fields that are left at None.
_PLUG_IN_DEFAULTS: dict[str, int] = {
    "num_chains": 1,
    "num_burnin": 200,
    "num_trees_mu": 75,
    "num_trees_tau_mode1": 35,
    "num_trees_tau_higher": 20,
}
_PROPAGATED_DEFAULTS: dict[str, int] = {
    "num_chains": 2,
    "num_burnin": 100,
    "num_trees_mu": 50,
    "num_trees_tau_mode1": 30,
    "num_trees_tau_higher": 15,
}


@dataclass(frozen=True)
class RankKLongBetConfig:
    """Configuration for Forecast-Augmented Rank-K Exposure-Basis BCF (v2).

    One estimator, two uncertainty variants, selected by
    ``uncertainty_propagation``; ``engine`` follows from it.

    * "propagated" (the default: ``uncertainty_propagation=True``,
      ``engine="jax"``) fits the baseline components (prognostic forest,
      forecast slope, factor shocks, spatial field) on never-treated units only
      (``baseline_cells="untreated"``), samples them jointly with the treatment
      forests, and pools ``n_forecast_draws`` posterior forecast draws. It is
      the headline Bayesian estimator validated in the companion paper: it
      keeps the baseline free of treatment-effect leakage and propagates
      forecast, slope, factor-loading, and spatial/cluster uncertainty into
      every equal-tailed posterior percentile interval. Its sampler is an
      approximate partially collapsed Gibbs sampler (working likelihoods in the
      forest steps).
    * "plug-in" (``uncertainty_propagation=False``, ``engine="stochtree"``) fits
      one ``stochtree`` BCF per mode across all units on a plug-in OLS-ridge
      forecast slope fitted on never-treated units, without forecast draws,
      factor loadings or spatial field. Use it as a fast plug-in approximation
      or when unit-level CATE point accuracy on small covariate-stratified
      panels takes priority over joint uncertainty propagation.

    Evidence (companion paper, at the validated budgets below). Both variants
    met every ATT coverage rule of the second frozen confirmation (pooled 95%
    quantile-interval coverage 94.3% propagated, 93.0% plug-in across 768
    replications in Sets E--G), whereas the plug-in variant missed the 90%
    threshold by one interval on Set A of the first confirmation (230/256 =
    89.8%). The propagated default eliminates treatment-effect leakage into the
    baseline (null and effect ATT errors match to SD 0.005 vs 0.068 for the
    plug-in variant) and passes four of six pre-committed criteria on a six-DGP
    synthetic study with known spatial correlation, factor shocks, and slope
    uncertainty, at a +4% (+0.010) unit-PEHE trade-off on N=250 panels because
    its prognostic forest is learned on the N_0=125 never-treated controls only.

    Variant-dependent defaults. ``engine``, ``num_chains``, ``num_burnin`` and
    the three tree counts default to ``None`` and are resolved in
    ``__post_init__`` to the budget validated for the selected variant, so a
    constructed config always holds concrete values:

    * propagated (default): ``engine="jax"`` with the benchmark budget of the
      paper's development benchmark and second frozen hold-out confirmation:
      ``num_chains=2``, ``num_burnin=100``, ``draws_per_component=25``,
      ``n_forecast_draws=8`` without adaptive expansion, 50 / 30 / 15 trees,
      ``use_spatial_gp=True``, ``use_factor_loadings=True`` with
      ``n_factor_loadings=3``, and ``baseline_cells="untreated"``. For small
      panels use the paper's application budget,
      ``RankKLongBetConfig(uncertainty_propagation=True, num_chains=4,
      num_burnin=300, draws_per_component=100)``. (The paper's runs also set
      ``warm_burnin`` to 60 and 150, below ``num_burnin``, where it has no
      effect.)
    * plug-in: ``engine="stochtree"``, one ``stochtree`` chain per mode
      (``num_chains=1``) with ``num_gfr=10`` grow-from-root sweeps,
      ``num_burnin=200`` burn-in sweeps and ``num_mcmc=600`` retained draws,
      and 75 / 35 / 20 trees for ``mu_k`` / ``tau_1`` / ``tau_k`` (``k >= 2``),
      the slope clipped to ``theta_f_clip=(0, 1.5)``: the configuration behind
      every plug-in result of the paper.

    ``max_forecast_draws`` defaults to ``None``, which means
    ``n_forecast_draws`` (no adaptive expansion); it stays ``None``, so
    ``dataclasses.replace(cfg, n_forecast_draws=...)`` never turns the
    expansion on. Every other field has a single default; an option that only
    one variant reads is ignored by the other, and a non-default value of it
    then triggers a :class:`RankKConfigWarning`. ``engine`` is kept for
    backward compatibility: ``engine="jax"`` alone selects the propagated
    variant and ``engine="stochtree"`` alone the plug-in variant, whereas
    ``engine="jax"`` with ``uncertainty_propagation=False`` and
    ``engine="stochtree"`` with ``uncertainty_propagation=True`` raise
    ``ValueError``. Because the resolved engine and budget are stored,
    ``dataclasses.replace`` cannot switch variants: replacing
    ``uncertainty_propagation`` alone contradicts the stored engine and raises
    ``ValueError``, and replacing it together with ``engine=None`` keeps the
    other variant's budget (the fit then warns). Build a new config to switch
    variants.

    Posterior draws. The plug-in variant retains ``num_mcmc`` draws per mode.
    In the propagated variant every chain is run once per forecast draw: the
    ``num_chains x n_forecast_draws`` tracks start from the same initial state,
    each runs its own ``num_burnin`` burn-in sweeps and retains
    ``draws_per_component`` draws, so the posterior has ``num_chains x
    n_forecast_draws x draws_per_component`` draws (``2 x 8 x 25 = 400`` at the
    defaults, ``4 x 8 x 100 = 3200`` at the application budget). When it
    propagates no forecast draws (no forecast, ``use_forecast=False`` or
    ``n_forecast_draws=1``) it retains ``num_mcmc`` draws in total, split evenly
    across chains. Forecast draws are propagated through mode 1 only; the
    higher (dynamic) modes condition on the point forecast. The ATT and the
    unit-averaged effects depend on mode 1 alone, and because mode 1 is seeded
    with ``random_seed`` whatever ``rank_k`` is, they are invariant to ``K`` by
    construction (Proposition 2) in both variants.

    Intervals (both variants). ``att_ci_low`` / ``att_ci_high`` are the
    2.5% / 97.5% posterior percentiles of the pooled ATT draws, the interval
    whose coverage the paper evaluated; ``catt_ci_low`` / ``catt_ci_high`` are
    the 2.5% / 97.5% posterior percentiles of each cell effect. Unit-level
    (CATE) intervals are approximate and below nominal in the paper's
    benchmarks (95% unit-level intervals covered 88-93% with the plug-in and
    70-92% with the propagated variant); treat them as descriptive.

    Options and warnings. In the propagated variant ``use_spatial_gp=False``,
    ``use_factor_loadings=False`` and ``baseline_cells="all"`` reproduce
    ablations of the paper (each removes one component, and ``"all"`` lets the
    treatment effect leak into the baseline); use them for sensitivity
    analyses. ``max_forecast_draws > n_forecast_draws`` (the adaptive
    expansion) and ``rank_k > 3`` were not validated. Choosing any of these
    never warns: :func:`fit_rank_k_longbet` emits a :class:`RankKConfigWarning`
    only when an explicit input or a non-default option would be silently
    ignored or contradicts another one (including a budget equal to the other
    variant's default budget, as left by ``dataclasses.replace``), when the
    deprecated ``warm_burnin`` is set, and when cohorts are too short for
    ``rank_k`` modes. Invalid values, including counts that are not whole
    numbers, flags that are not booleans and non-finite reals, raise
    ``ValueError`` at construction.

    Attributes:
      rank_k: Number ``K >= 1`` of exposure modes (default 2: level + concave
        ramp). Only ``K`` in {1, 2, 3} was validated. A cohort with ``L_g``
        post-adoption periods gets ``min(K, L_g)`` modes.
      use_forecast: Use the counterfactual forecast (``features.point_forecast``
        or ``baseline_forecast``) as the mode-wise control variate and as a
        prognostic-forest covariate (default True). Without a forecast the model
        has no control variate and no forecast draws.
      engine: ``None`` (default: derived from ``uncertainty_propagation``),
        ``"jax"`` (the propagated default) or ``"stochtree"`` (the plug-in
        variant); resolved to ``"jax"`` or ``"stochtree"``. The plug-in variant
        uses ``rank_k``, ``use_forecast``, ``num_gfr``, ``num_burnin``,
        ``num_mcmc``, the tree counts, ``theta_f_clip``, ``tau_prior_mean``,
        ``tau_prior_sd``, ``tau_clip``, ``catt_floor``, ``ridge_eps``,
        ``num_threads`` and ``random_seed``, and ignores the other options (a
        non-default value triggers a :class:`RankKConfigWarning`) as well as the
        ``coords``, ``cluster_ids``, ``pre_anchor_loadings`` and
        ``baseline_forecast_draws`` arguments.
      uncertainty_propagation: Selects the variant: ``True`` (the default) the
        propagated variant, which pools ``n_forecast_draws`` posterior forecast
        draws into a modular (cut) posterior; ``False`` the fast plug-in
        variant. ``None`` (default) means ``False`` if ``engine="stochtree"``
        and ``True`` otherwise; resolved to a bool.
      n_forecast_draws: Number ``M >= 0`` of posterior forecast draws of the
        propagated variant (default 8, as validated). ``M = 0`` or ``1``
        propagates no forecast uncertainty (with ``uncertainty_propagation=True``
        it triggers a :class:`RankKConfigWarning`).
      max_forecast_draws: Upper limit ``>= n_forecast_draws`` of the adaptive
        expansion (default ``None``: ``n_forecast_draws``, no expansion, as
        validated). Larger values enable the expansion, which the paper did not
        validate: mode 1 is refitted with ``max_forecast_draws`` draws whenever
        the between-forecast share of the mode-1 ATT variance exceeds
        ``between_var_threshold``.
      between_var_threshold: Between-forecast variance share in ``[0, 1]`` above
        which the adaptive expansion is triggered (default 0.20; 1 never
        triggers it); unused unless ``max_forecast_draws > n_forecast_draws``.
      num_chains: Independent chains per forecast draw of the propagated
        variant, ``>= 1`` (variant-dependent default: 2 propagated; 1 plug-in,
        whose ``stochtree`` sampler runs one chain per mode).
      draws_per_component: Retained draws per chain and forecast draw of the
        propagated variant, ``>= 2`` (default 25; application budget 100).
      num_burnin: Burn-in sweeps of the ``stochtree`` sampler and of every
        chain-by-forecast-draw track of the JAX sampler, ``>= 0``
        (variant-dependent default: 200 plug-in, 100 propagated).
      warm_burnin: Deprecated (default ``None``). It is not a warm start: when
        forecast draws are propagated, every mode-1 track runs ``max(num_burnin,
        warm_burnin)`` burn-in sweeps, so it has no effect unless it exceeds
        ``num_burnin``, and it is ignored everywhere else. Setting it triggers a
        :class:`RankKConfigWarning`; set ``num_burnin`` instead.
      num_gfr: Grow-from-root sweeps of the plug-in variant, ``>= 0`` (default
        10).
      num_mcmc: Retained posterior draws, ``>= 2`` (default 600): per mode in the
        plug-in variant; in total (``round(num_mcmc / num_chains)`` per chain) in
        the propagated variant when it propagates no forecast draws, and ignored
        when it does.
      num_trees_mu: Trees of each mode's prognostic forest ``mu_k``, ``>= 1``
        (variant-dependent default: 75 plug-in, 50 propagated).
      num_trees_tau_mode1: Trees of the mode-1 treatment forest, ``>= 1``
        (variant-dependent default: 35 plug-in, 30 propagated).
      num_trees_tau_higher: Trees of each higher-mode treatment forest, ``>= 1``
        (variant-dependent default: 20 plug-in, 15 propagated).
      theta_f_clip: Truncation interval ``(low, high)``, ``low <= high``, of the
        forecast control-variate slope ``theta_k`` (default ``(0.0, 1.5)``);
        ``low == high`` fixes the slope. The plug-in variant clips its OLS-ridge
        slope to this interval.
      theta_prior_mean: Mean of the Normal slope prior before truncation (default
        1.0: the forecast taken at face value; propagated variant).
      theta_prior_sd: Standard deviation ``> 0`` of the slope prior (default 0.5;
        propagated variant).
      tau_prior_mean: Mean of the Normal prior on the Mode-1 time-averaged
        treatment effect intercept (in per-period outcome units) before
        truncation (default 0.0).
      tau_prior_sd: Standard deviation ``> 0`` of the Normal prior on the Mode-1
        treatment effect intercept (in per-period outcome units), or ``None``
        (default: diffuse prior of 10.0 on standardized control-residual units).
      tau_clip: Truncation interval ``(low, high)``, ``low <= high`` (endpoints
        may be ``+/-inf``), for the Mode-1 time-averaged treatment effect
        intercept in per-period outcome units (default ``(-inf, +inf)``; pass
        ``(0.0, float('inf'))`` for a domain-informed non-negative ATT prior).
      catt_floor: Optional lower bound on post-adoption cell treatment effects
        ``tau_{i,t}`` for ``z_{i,t} == 1`` in per-period outcome units (default
        ``None``; pass ``0.0`` to enforce per-period non-negativity
        ``tau_{i,t} >= 0``). Enforced via cohort-mean-preserving Mode-1 unit
        damping and zero-sum higher-mode amplitude damping
        ``tau_{i,t} = tau_{i,1} + alpha_i * sum_{k>=2} psi_k(t) tau_{i,k}``,
        preserving exact ``K``-invariance of ``mode_tau_draws[0]`` and
        ``att_draws`` (Proposition 2) without Jensen's inequality inflation.
      use_factor_loadings: Include the pre-anchor factor term ``lambda_i' f_k``
        (default True; propagated variant).
      n_factor_loadings: Number ``q >= 0`` of loadings extracted from the
        in-panel pre-period residuals when neither ``features`` nor
        ``pre_anchor_loadings`` supplies them (default 3; propagated variant).
      use_spatial_gp: Include the spatial random field ``zeta_k`` (default True;
        propagated variant): a Matérn-3/2 GP on great-circle distances when
        ``coords`` are supplied, or a cluster random effect when only
        ``cluster_ids`` are, estimated on never-treated units and kriged to
        treated units. Without ``coords`` or ``cluster_ids`` no spatial field is
        fitted.
      lengthscale_grid_km: Non-empty grid of positive Matérn ranges in km, with a
        discrete uniform prior (default ``(100, 200, 400, 800, 1600)``); used with
        ``coords`` only.
      baseline_cells: ``"untreated"`` (default; the control-only cut: every
        baseline component is fitted on never-treated units) or ``"all"`` (lets the
        treatment effect leak into the baseline; sensitivity analyses only;
        propagated variant).
      ridge_eps: Ridge term ``>= 0`` of the plug-in slope (default 1e-4; plug-in
        variant).
      num_threads: Threads of the ``stochtree`` sampler, ``>= 1`` or ``-1`` for
        all available threads (default 1; plug-in variant).
      random_seed: Base seed ``>= 0`` (default 2026). Mode ``k`` uses
        ``random_seed + 101 * (k - 1)``, so mode 1 is seeded identically for
        every ``rank_k``.
      tau_prior_family: Family of the prior on the mode-1 ATT when ``tau_prior_sd``
        is set: ``"normal"`` (default; conjugate, truncated to ``tau_clip``),
        ``"halfstudentt"`` (``T^+_{2.5}(0, (0.5 tau_prior_sd)^2)``) or
        ``"gamma_spike"`` (``Gamma(0.40, 1.60 tau_prior_sd)``). With any prior or a
        finite ``tau_clip`` the mode-1 treatment forest is re-centered with the ATT
        weights, so the sampled intercept is the ATT and the prior / truncation acts
        on it exactly inside the sampler (JAX engine only; the plug-in engine raises).
      short_cohort_policy: What to do with cohorts whose post-adoption length ``L_g``
        is below ``rank_k``: ``"drop_mode"`` (default; they take part in modes
        ``k <= L_g`` only and are excluded from the fit of the higher modes),
        ``"reduce_k"`` (lower ``K`` to ``min_g L_g``) or ``"error"``.
      thin: Keep every ``thin``-th post-burn-in sweep of the JAX sampler, ``>= 1``
        (default 1: consecutive sweeps, the paper's budget).
      sigma_nu_total: Total leaf prior scale of the treatment forests in
        standardized units (default 0.25, the paper's value; JAX engine).
      sigma_nu_dynamic_total: Optional separate leaf scale for the dynamic modes
        ``k >= 2`` (default ``None``: ``sigma_nu_total``; JAX engine).
      alpha_nu: Depth-prior base of the treatment forests, in ``(0, 1)`` (default 0.25;
        JAX engine).
      beta_nu: Depth-prior power of the treatment forests, ``>= 0`` (default 3.0; JAX
        engine).
      nu_working_zeta_weight: Weight ``c`` of the spatial variance in the working noise
        variance ``sigma_k^2 + c s_zeta^2`` of the treatment-forest step, in ``[0, 1]``
        (default 1.0, the paper's value; JAX engine).
      learn_sigma_nu: Learn the treatment-forest leaf scale with a half-Cauchy
        hyperprior of scale ``sigma_nu_prior_scale`` instead of fixing it (default
        False; JAX engine).
      sigma_nu_prior_scale: Scale ``> 0`` of that half-Cauchy hyperprior (default 0.5).
      treated_extra_variance: Add a treated-arm excess residual variance ``s_eta^2``
        (idiosyncratic unit effects), learned from the treated residuals; its draws are
        in ``RankKLongBetResult.diagnostics["s_eta2_draws"]`` (default ``None``: ``True``
        when ``baseline_cells="untreated"`` and ``False`` otherwise; JAX engine,
        ``baseline_cells="untreated"`` only).
      basis_spec: :class:`ExposureBasisSpec` of the exposure-basis generators
        (default ``None``: the paper's ``{1, 1 - exp(-s / 3), s / L}`` plus cosines).
      keep_unprojected_draws: Keep the unprojected posterior draws (``*_raw``) whenever
        ``catt_floor`` or ``tau_clip`` is active (default True).
      att_weighting: Weights of the treated units in ``att_draws``: ``"cell"``
        (default; ``L_i / |O_trt|``, the paper's convention, ``L_i`` being the
        observed post-adoption cells under a mask) or ``"unit"`` (``1 / N_trt``,
        the population estimand of masked event-time panels).
      projection: Per-unit projection on the common event-time basis when windows
        differ (masked panels): ``"two_pass"`` (default; single-mode projections of
        the higher modes on their supported units, Mode 1 as the window mean after
        the fitted dynamics are removed) or ``"restricted_wls"`` (joint
        ``(Psi_O' Psi_O)^{-1} Psi_O' r_O`` of the supported modes, condition-number
        guarded; a sensitivity option).
      min_support_modes: Minimum observed post-adoption periods a unit needs to
        project mode ``k`` (index ``k``, 0-based; the last entry applies to higher
        modes), non-decreasing whole numbers (default ``(1, 6, 8)``). Units below
        support get precision weight 0 for that mode and a cross-sectional
        posterior prediction; units below ``min_support_modes[0]`` are an error.
      common_basis: ``None`` (default: one common basis on ``s = 1..H`` for masked
        and single-wave panels, the paper's per-cohort bases otherwise), ``True`` or
        ``False`` (not available for masked panels).
      mode_coupling: How the fitted higher modes enter the Mode-1 targets of
        partially observed units: ``"propagated"`` (default; ``n_forecast_draws``
        Mode-1 tracks, one higher-mode posterior draw each, paired with the
        forecast draws when present), ``"plug_in"`` (posterior means) or ``"none"``
        (no correction; an ablation). Irrelevant on fully observed panels.
      support_covariate: Append the observed window length ``L_i`` to the forest
        covariates when windows differ (default True).
    """

    rank_k: int = 2
    use_forecast: bool = True
    engine: str | None = None
    uncertainty_propagation: bool | None = None
    n_forecast_draws: int = 8
    max_forecast_draws: int | None = None
    between_var_threshold: float = 0.20
    num_chains: int | None = None
    draws_per_component: int = 25
    num_burnin: int | None = None
    warm_burnin: int | None = None
    num_gfr: int = 10
    num_mcmc: int = 600
    num_trees_mu: int | None = None
    num_trees_tau_mode1: int | None = None
    num_trees_tau_higher: int | None = None
    theta_f_clip: tuple[float, float] = (0.0, 1.5)
    theta_prior_mean: float = 1.0
    theta_prior_sd: float = 0.5
    tau_prior_mean: float = 0.0
    tau_prior_sd: float | None = None
    tau_clip: tuple[float, float] = (-float("inf"), float("inf"))
    catt_floor: float | None = None
    use_factor_loadings: bool = True
    n_factor_loadings: int = 3
    use_spatial_gp: bool = True
    lengthscale_grid_km: tuple[float, ...] = (100.0, 200.0, 400.0, 800.0, 1600.0)
    baseline_cells: str = "untreated"
    ridge_eps: float = 1e-4
    num_threads: int = 1
    random_seed: int = 2026
    # Remediation options (defaults reproduce the paper's estimator bit-for-bit; see the class docstring)
    tau_prior_family: str = "normal"
    short_cohort_policy: str = "drop_mode"
    thin: int = 1
    sigma_nu_total: float = 0.25
    sigma_nu_dynamic_total: float | None = None
    alpha_nu: float = 0.25
    beta_nu: float = 3.0
    nu_working_zeta_weight: float = 1.0
    learn_sigma_nu: bool = False
    sigma_nu_prior_scale: float = 0.5
    treated_extra_variance: bool | None = None
    basis_spec: ExposureBasisSpec | None = None
    keep_unprojected_draws: bool = True
    # Masked event-time panels (0.7.0): population estimand, per-unit windows, common basis
    att_weighting: str = "cell"
    projection: str = "two_pass"
    min_support_modes: tuple[int, ...] = (1, 6, 8)
    common_basis: bool | None = None
    mode_coupling: str = "propagated"
    support_covariate: bool = True

    def __post_init__(self) -> None:
        """Resolve the variant and its defaults; reject invalid values (raises ``ValueError``)."""

        def _fail(message: str) -> NoReturn:
            raise ValueError(f"RankKLongBetConfig: {message}")

        def _set(name: str, value: Any) -> None:
            object.__setattr__(self, name, value)  # the dataclass is frozen

        def _flag(name: str, value: Any, allowed: str = "True or False") -> bool:
            if isinstance(value, (bool, np.bool_)) or (isinstance(value, (int, np.integer)) and value in (0, 1)):
                return bool(value)
            _fail(f"{name} must be {allowed}; got {value!r}")

        def _count(name: str, minimum: int, *, optional: bool = False) -> None:
            # Whole numbers only: int() would silently truncate 2.7 to 2 in the samplers.
            value = getattr(self, name)
            if value is None and optional:
                return
            whole = isinstance(value, (int, np.integer)) or (
                isinstance(value, (float, np.floating)) and bool(np.isfinite(value)) and float(value).is_integer()
            )
            if isinstance(value, (bool, np.bool_)) or not whole:
                _fail(f"{name} must be an integer >= {minimum}; got {value!r}")
            if value < minimum:
                _fail(f"{name} must be >= {minimum}; got {value!r}")
            _set(name, int(value))

        def _real(name: str) -> float:
            value = getattr(self, name)
            is_number = isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(
                value, (bool, np.bool_)
            )
            if not is_number or not np.isfinite(value):
                _fail(f"{name} must be a finite number; got {value!r}")
            _set(name, float(value))
            return float(value)

        def _reals(name: str) -> tuple[float, ...] | None:
            # A non-string sequence of numbers as a tuple of floats; None if it is not one.
            value = getattr(self, name)
            if isinstance(value, (str, bytes)):
                return None
            try:
                items = list(value)
            except TypeError:
                return None
            if not all(
                isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))
                for v in items
            ):
                return None
            return tuple(float(v) for v in items)

        # ---- Variant: uncertainty_propagation selects it, engine follows ----------------
        if self.engine not in (None, "jax", "stochtree"):
            _fail(f"engine must be None, 'jax' or 'stochtree'; got {self.engine!r}")
        propagation = self.uncertainty_propagation
        if propagation is not None:
            propagation = _flag("uncertainty_propagation", propagation, allowed="None, True or False")
        switch = (
            " To switch the variant of an existing config, build a new RankKLongBetConfig instead of "
            "calling dataclasses.replace, which keeps the stored engine and budget."
        )
        if self.engine == "jax" and propagation is False:
            _fail(
                "engine='jax' with uncertainty_propagation=False is not a supported variant: the JAX "
                "sampler without forecast draws is only a synthetic-study analogue of the plug-in "
                "variant, and its intervals under-covered (60-84%) under spatially correlated shocks. "
                "Use RankKLongBetConfig(uncertainty_propagation=False) for the plug-in variant (engine "
                "'stochtree') or RankKLongBetConfig(uncertainty_propagation=True) for the propagated "
                "variant." + switch
            )
        if self.engine == "stochtree" and propagation is True:
            _fail(
                "engine='stochtree' with uncertainty_propagation=True is not implemented: the "
                "stochtree (plug-in) variant does not propagate forecast draws. Use "
                "RankKLongBetConfig(uncertainty_propagation=True) for the propagated variant (engine "
                "'jax') or RankKLongBetConfig(uncertainty_propagation=False) for the plug-in variant."
                + switch
            )
        if propagation is None:
            propagation = self.engine != "stochtree"
        _set("uncertainty_propagation", propagation)
        _set("engine", "jax" if propagation else "stochtree")
        for name, value in (_PROPAGATED_DEFAULTS if propagation else _PLUG_IN_DEFAULTS).items():
            if getattr(self, name) is None:
                _set(name, value)

        # ---- Types and values: counts are whole numbers, flags booleans, reals finite ----
        _count("rank_k", 1)
        if self.baseline_cells not in ("untreated", "all"):
            _fail(f"baseline_cells must be 'untreated' or 'all'; got {self.baseline_cells!r}")
        for name in ("use_forecast", "use_factor_loadings", "use_spatial_gp"):
            _set(name, _flag(name, getattr(self, name)))
        # 0 and 1 both mean that no forecast draws are propagated (as in 0.4.0).
        _count("n_forecast_draws", 0)
        # max_forecast_draws=None stays None (it means n_forecast_draws), so that
        # dataclasses.replace(cfg, n_forecast_draws=...) never turns the expansion on.
        _count("max_forecast_draws", 0, optional=True)
        if self.max_forecast_draws is not None and self.max_forecast_draws < self.n_forecast_draws:
            _fail(
                f"max_forecast_draws ({self.max_forecast_draws!r}) must be >= n_forecast_draws "
                f"({self.n_forecast_draws!r}); leave it at None (or set it equal) to disable the "
                "adaptive expansion"
            )
        _count("num_chains", 1)
        _count("num_burnin", 0)
        _count("warm_burnin", 0, optional=True)
        _count("draws_per_component", 2)
        _count("num_mcmc", 2)
        for name in ("num_trees_mu", "num_trees_tau_mode1", "num_trees_tau_higher"):
            _count(name, 1)
        _count("num_gfr", 0)
        _count("n_factor_loadings", 0)
        _count("num_threads", -1)
        if self.num_threads == 0:
            _fail("num_threads must be >= 1, or -1 for all available threads; got 0")
        _count("random_seed", 0)
        clip = _reals("theta_f_clip")
        if clip is None or len(clip) != 2 or not clip[0] <= clip[1]:
            _fail(f"theta_f_clip must be a (low, high) pair with low <= high; got {self.theta_f_clip!r}")
        _set("theta_f_clip", clip)
        _real("theta_prior_mean")
        if not _real("theta_prior_sd") > 0:
            _fail(f"theta_prior_sd must be > 0; got {self.theta_prior_sd!r}")
        _real("tau_prior_mean")
        if self.tau_prior_sd is not None:
            if not _real("tau_prior_sd") > 0:
                _fail(f"tau_prior_sd must be > 0 or None; got {self.tau_prior_sd!r}")
        t_clip = _reals("tau_clip")
        if (
            t_clip is None
            or len(t_clip) != 2
            or np.isnan(t_clip[0])
            or np.isnan(t_clip[1])
            or not t_clip[0] <= t_clip[1]
            or (np.isinf(t_clip[0]) and t_clip[0] == t_clip[1])
        ):
            _fail(f"tau_clip must be a (low, high) pair with low <= high; got {self.tau_clip!r}")
        _set("tau_clip", t_clip)
        if self.catt_floor is not None:
            _real("catt_floor")
        if not _real("ridge_eps") >= 0:
            _fail(f"ridge_eps must be >= 0; got {self.ridge_eps!r}")
        grid = _reals("lengthscale_grid_km")
        if not grid or not all(np.isfinite(v) and v > 0.0 for v in grid):
            _fail(
                "lengthscale_grid_km must be a non-empty sequence of positive, finite ranges "
                f"in km; got {self.lengthscale_grid_km!r}"
            )
        _set("lengthscale_grid_km", grid)
        # share > threshold triggers the expansion: 1 never does, 0 on any between-draw variance.
        if not 0.0 <= _real("between_var_threshold") <= 1.0:
            _fail(f"between_var_threshold must lie in [0, 1]; got {self.between_var_threshold!r}")

        # ---- Remediation options -------------------------------------------------------------
        if self.tau_prior_family not in _TAU_PRIOR_FAMILIES:
            _fail(f"tau_prior_family must be one of {_TAU_PRIOR_FAMILIES}; got {self.tau_prior_family!r}")
        if self.tau_prior_family != "normal" and self.tau_prior_sd is None:
            _fail(f"tau_prior_family={self.tau_prior_family!r} needs tau_prior_sd (its scale)")
        if self.short_cohort_policy not in _SHORT_COHORT_POLICIES:
            _fail(f"short_cohort_policy must be one of {_SHORT_COHORT_POLICIES}; got {self.short_cohort_policy!r}")
        _count("thin", 1)
        if not _real("sigma_nu_total") > 0:
            _fail(f"sigma_nu_total must be > 0; got {self.sigma_nu_total!r}")
        if self.sigma_nu_dynamic_total is not None and not _real("sigma_nu_dynamic_total") > 0:
            _fail(f"sigma_nu_dynamic_total must be > 0 or None; got {self.sigma_nu_dynamic_total!r}")
        if not 0.0 < _real("alpha_nu") < 1.0:
            _fail(f"alpha_nu must lie in (0, 1); got {self.alpha_nu!r}")
        if not _real("beta_nu") >= 0.0:
            _fail(f"beta_nu must be >= 0; got {self.beta_nu!r}")
        if not 0.0 <= _real("nu_working_zeta_weight") <= 1.0:
            _fail(f"nu_working_zeta_weight must lie in [0, 1]; got {self.nu_working_zeta_weight!r}")
        if not _real("sigma_nu_prior_scale") > 0:
            _fail(f"sigma_nu_prior_scale must be > 0; got {self.sigma_nu_prior_scale!r}")
        for name in ("learn_sigma_nu", "keep_unprojected_draws"):
            _set(name, _flag(name, getattr(self, name)))
        if self.treated_extra_variance is None:
            _set("treated_extra_variance", self.baseline_cells == "untreated")
        else:
            _set(
                "treated_extra_variance",
                _flag("treated_extra_variance", self.treated_extra_variance, allowed="None, True or False"),
            )
        if self.treated_extra_variance and self.baseline_cells != "untreated":
            _fail("treated_extra_variance=True requires baseline_cells='untreated'")
        if self.basis_spec is not None and not isinstance(self.basis_spec, ExposureBasisSpec):
            _fail("basis_spec must be an ExposureBasisSpec or None")

        # ---- Masked event-time options ------------------------------------------------------
        if self.att_weighting not in _ATT_WEIGHTINGS:
            _fail(f"att_weighting must be one of {_ATT_WEIGHTINGS}; got {self.att_weighting!r}")
        if self.projection not in _PROJECTIONS:
            _fail(f"projection must be one of {_PROJECTIONS}; got {self.projection!r}")
        if self.mode_coupling not in _MODE_COUPLINGS:
            _fail(f"mode_coupling must be one of {_MODE_COUPLINGS}; got {self.mode_coupling!r}")
        support = _reals("min_support_modes")
        if (
            support is None
            or len(support) == 0
            or any(not float(v).is_integer() or v < 1 for v in support)
            or any(b < a for a, b in zip(support, support[1:]))
        ):
            _fail(
                "min_support_modes must be a non-empty, non-decreasing sequence of whole numbers >= 1 (the minimum "
                f"observed post-adoption periods needed to project mode k); got {self.min_support_modes!r}"
            )
        _set("min_support_modes", tuple(int(v) for v in support))
        if self.common_basis is not None:
            _set("common_basis", _flag("common_basis", self.common_basis, allowed="None, True or False"))
        _set("support_covariate", _flag("support_covariate", self.support_covariate))

    @property
    def has_att_constraint(self) -> bool:
        """Whether a prior or truncation acts on the ATT (then the sampler uses ATT-weight centering)."""
        return bool(
            self.tau_prior_sd is not None
            or math.isfinite(float(self.tau_clip[0]))
            or math.isfinite(float(self.tau_clip[1]))
            or self.tau_prior_family != "normal"
        )

    @property
    def resolved_basis_spec(self) -> ExposureBasisSpec:
        return DEFAULT_BASIS_SPEC if self.basis_spec is None else self.basis_spec


@dataclass
class RankKLongBetResult:
    """Posterior summaries, draws, and variance-decomposition diagnostics from Rank-K BCF.

    ``S`` is the number of retained posterior draws: ``num_mcmc`` in the plug-in
    variant; in the propagated variant ``num_chains x n_components x
    draws_per_component`` when forecast draws are propagated (``n_components``
    is ``n_forecast_draws``, or ``max_forecast_draws`` when the adaptive
    expansion triggers), otherwise ``num_chains x round(num_mcmc /
    num_chains)``. In both variants
    ``att_ci_low`` / ``att_ci_high`` are the 2.5% / 97.5% posterior percentiles
    of ``att_draws`` (the interval whose coverage the paper evaluated), and
    ``catt_ci_low`` / ``catt_ci_high`` are 2.5% / 97.5% posterior percentiles of
    each cell effect (zero on untreated cells). Unit-level intervals computed
    from ``mode_tau_draws`` are approximate and below nominal in the paper's
    benchmarks and the cell-level intervals are approximate as well, so treat
    both as descriptive (see :func:`fit_rank_k_longbet`). Row ``i`` of
    ``mode_tau_draws[0]`` holds the draws of unit ``i``'s time-averaged effect
    over its post-adoption periods; it is meaningful for treated units only.
    ``rhat_att`` needs at least two chains, so it is NaN for the plug-in
    variant (one ``stochtree`` chain per mode); ``ess_att`` is defined for both.
    """

    catt_draws: np.ndarray          # (S, N, T) posterior draws of cell treatment effects
    catt_hat: np.ndarray            # (N, T) posterior mean cell treatment effects
    catt_ci_low: np.ndarray         # (N, T) 2.5% posterior quantile
    catt_ci_high: np.ndarray        # (N, T) 97.5% posterior quantile
    att_draws: np.ndarray           # (S,) posterior draws of overall ATT
    att_hat: float
    att_se: float
    att_ci_low: float
    att_ci_high: float
    mode_tau_draws: list[np.ndarray]  # length K list of (N, S) posterior draws per mode
    theta_f_modes: list[float]        # length K list of posterior mean CUPED forecast slopes
    sigma_modes: list[float]          # length K list of posterior mean mode error SDs
    bases: dict[int, np.ndarray]      # wave g -> (L_g, K_eff) orthonormal exposure basis
    runtime_sec: float
    theta_f_draws: list[np.ndarray] = dc_field(default_factory=list)
    gp_hyper_draws: list[dict[str, np.ndarray]] = dc_field(default_factory=list)
    att_var_within: float = 0.0
    att_var_between: float = 0.0
    att_var_slope: float = 0.0
    n_components: int = 1
    rhat_att: float = 1.0
    ess_att: float = 0.0
    # Remediation additions. ``*_raw`` hold the unprojected posterior draws whenever a
    # post-processing projection (``catt_floor``, or the no-op-up-to-rounding
    # ``tau_clip`` guard) was applied and ``keep_unprojected_draws`` is set;
    # ``projection_report`` quantifies the projection; ``diagnostics`` holds per-mode
    # convergence summaries and the extension traces (``s_eta2_draws``,
    # ``sigma_nu_total_draws``, ``tau_intercept_draws``); ``active_units[k]`` is the
    # boolean mask of units that took part in the mode-``k`` fit.
    att_draws_raw: np.ndarray | None = None
    mode_tau_draws_raw: list[np.ndarray] | None = None
    catt_draws_raw: np.ndarray | None = None
    projection_report: dict[str, Any] | None = None
    diagnostics: dict[str, Any] = dc_field(default_factory=dict)
    variant: str = "propagated"
    basis_spec: Any = None
    active_units: list[np.ndarray] = dc_field(default_factory=list)
    # Masked event-time additions (0.7.0). ``catt_observed[i, t]`` marks the treated
    # post cells inside the horizon that were observed (the other treated cells of
    # ``catt_draws`` are posterior-predictive); ``horizon`` is the event-time horizon
    # ``H`` of the common basis (``None`` for per-cohort bases); ``unit_weights`` are
    # the ATT weights behind ``att_draws`` (``att_weighting`` "cell" or "unit");
    # ``support_mask[i, k]`` tells whether unit ``i``'s observed window supported the
    # projection of mode ``k``; ``post_window_length`` is ``L_i``; ``cohorts`` the
    # first-adoption period of every unit (``T`` for never-treated units).
    catt_observed: np.ndarray | None = None
    horizon: int | None = None
    att_weighting: str = "cell"
    unit_weights: np.ndarray | None = None
    support_mask: np.ndarray | None = None
    post_window_length: np.ndarray | None = None
    cohorts: np.ndarray | None = None
    arm_metadata: dict[str, Any] | None = None
    rank_k: int = 0

    def __post_init__(self) -> None:
        if self.rank_k <= 0:
            self.rank_k = max(len(self.mode_tau_draws), 1)

    @property
    def treated_units(self) -> np.ndarray:
        """``(N,)`` ever-treated indicator (units with a positive ATT weight)."""
        if self.unit_weights is not None:
            return np.asarray(self.unit_weights) > 0.0
        return np.any(np.abs(np.asarray(self.catt_draws)) > 0.0, axis=(0, 2))

    @property
    def att_var_total(self) -> float:
        return float(self.att_var_within + self.att_var_between)

    @property
    def att_between_share(self) -> float:
        tot = self.att_var_total
        return float(self.att_var_between / tot) if tot > 1e-12 else 0.0

    @property
    def s_zeta_modes(self) -> list[float]:
        return [float(np.mean(g["s_zeta"])) for g in self.gp_hyper_draws]

    @property
    def ell_modes(self) -> list[float]:
        return [float(np.mean(g["ell"])) for g in self.gp_hyper_draws]


def _basis_generators(
    L: int,
    k_eff: int,
    spec: ExposureBasisSpec,
    empirical_trajectories: np.ndarray | None,
) -> list[np.ndarray]:
    """Raw (pre-orthonormalization) generators on ``s = 1..L``; the first is the constant."""
    s = np.arange(1, L + 1, dtype=np.float64)
    ones = np.ones(L, dtype=np.float64)
    if spec.family == "exponential":
        gens = [ones, 1.0 - np.exp(-s / float(spec.tau_c)), s / float(max(L, 1))]
        for j in range(4, k_eff + 1):
            gens.append(np.cos((j - 2) * np.pi * (s - 0.5) / float(L)))
        return gens[:k_eff]
    if spec.family == "polynomial":
        u = (2.0 * s - float(L + 1)) / float(max(L - 1, 1))
        return [ones] + [np.polynomial.legendre.Legendre.basis(j)(u) for j in range(1, k_eff)]
    if spec.family == "cosine":
        return [ones] + [np.cos(j * np.pi * (s - 0.5) / float(L)) for j in range(1, k_eff)]
    if spec.family == "pca_empirical":
        if empirical_trajectories is None:
            raise ValueError("ExposureBasisSpec(family='pca_empirical') needs empirical_trajectories (n, L)")
        traj = np.asarray(empirical_trajectories, dtype=np.float64)
        if traj.ndim != 2 or traj.shape[1] != L:
            raise ValueError(f"empirical_trajectories must have shape (n, L={L}); got {traj.shape}")
        shape_only = traj - traj.mean(axis=1, keepdims=True)
        shape_only = shape_only - shape_only.mean(axis=0, keepdims=True)
        _, _, vt = np.linalg.svd(shape_only, full_matrices=False)
        gens = [ones]
        for j in range(k_eff - 1):
            if j < vt.shape[0]:
                v = vt[j]
                # Deterministic sign: positive end-point so the mode reads as "late minus early"
                gens.append(v if v[-1] >= v[0] else -v)
            else:
                gens.append(np.cos((j + 1) * np.pi * (s - 0.5) / float(L)))
        return gens
    # custom
    custom = np.asarray(spec.custom_generators(s, L), dtype=np.float64)  # type: ignore[misc]
    if custom.ndim == 1:
        custom = custom[None, :]
    if custom.ndim != 2 or custom.shape[1] != L:
        raise ValueError(f"custom_generators(s, L) must return an array of shape (K - 1, L={L}); got {custom.shape}")
    if custom.shape[0] < k_eff - 1:
        raise ValueError(f"custom_generators returned {custom.shape[0]} generator(s); {k_eff - 1} are needed")
    return [ones] + [custom[j] for j in range(k_eff - 1)]


def orthonormal_exposure_basis(
    L: int,
    K: int,
    spec: ExposureBasisSpec | None = None,
    empirical_trajectories: np.ndarray | None = None,
) -> np.ndarray:
    """Construct an ``(L, min(K, L))`` orthonormal exposure basis on ``s = 1, ..., L``.

    Basis functions are orthonormal with respect to the empirical exposure inner product
    ``<u, v>_L = (1 / L) sum_{s=1}^L u(s) v(s)``. With the default generators
    (``spec=None`` or ``ExposureBasisSpec()``):
      * ``k = 1``: ``psi_1(s) = 1`` (time-averaged level mode)
      * ``k = 2``: orthonormalized exponential saturation ``1 - exp(-s / 3)`` (concave ramp mode)
      * ``k = 3``: orthonormalized linear slope ``s / L`` (monotone tilt / wear-out contrast)
      * ``k >= 4``: orthonormalized cosine harmonics ``cos((k - 2) * pi * (s - 0.5) / L)``
    Other :class:`ExposureBasisSpec` families replace the generators ``k >= 2``; the
    constant first mode and the Gram--Schmidt orthonormalization are the same for all.
    """
    if L <= 0:
        raise ValueError(f"Horizon L must be positive; got {L}")
    if K <= 0:
        raise ValueError(f"Rank K must be positive; got {K}")
    spec = DEFAULT_BASIS_SPEC if spec is None else spec
    k_eff = min(int(K), int(L))
    raw_generators = _basis_generators(int(L), k_eff, spec, empirical_trajectories)

    cols: list[np.ndarray] = []
    for k in range(k_eff):
        v = raw_generators[k].copy()
        for q in cols:
            proj = float(np.mean(v * q))
            v = v - proj * q
        norm = float(np.sqrt(np.mean(v ** 2)))
        if norm > 1e-8:
            cols.append(v / norm)
        else:
            cols.append(np.zeros(L, dtype=np.float64))
    return np.column_stack(cols)


def _infer_cohorts(z: np.ndarray) -> np.ndarray:
    """Infer first-adoption period ``G_i in {0, ..., T}`` from binary absorbing matrix ``z``."""
    N, T = z.shape
    cohorts = np.full(N, T, dtype=np.int32)
    for i in range(N):
        ones = np.nonzero(z[i] == 1)[0]
        if ones.size > 0:
            cohorts[i] = int(ones[0])
    return cohorts


def _nanmean_rows(a: np.ndarray) -> np.ndarray:
    """Row means over the observed (non-NaN) cells, shape ``(n, 1)``; NaN for all-NaN rows (silently)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(np.asarray(a, dtype=np.float64), axis=1, keepdims=True)


def _nanmean_cols(a: np.ndarray) -> np.ndarray:
    """Column means over the observed (non-NaN) cells, shape ``(1, m)``; NaN for all-NaN columns (silently)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(np.asarray(a, dtype=np.float64), axis=0, keepdims=True)


def _att_unit_weights(
    cohorts_arr: np.ndarray,
    waves: list[int],
    T: int,
    total_trt_cells: float,
    weighting: str,
    layout: EventTimeLayout | None,
) -> np.ndarray:
    """ATT weights ``a_i`` on the treated units (0 on never-treated units).

    ``"cell"``: ``L_i / |O_trt|`` -- the paper's convention (``L_{G_i} / |O_trt|``),
    ``L_i`` being the observed post-adoption cells inside the horizon under a mask.
    ``"unit"``: ``1 / N_trt``, the population estimand (every treated unit counts once).
    """
    N = int(cohorts_arr.shape[0])
    att_w = np.zeros(N, dtype=np.float64)
    ever = cohorts_arr < T
    if weighting == "unit":
        att_w[ever] = 1.0 / max(float(np.sum(ever)), 1.0)
        return att_w
    for g in waves:
        rows = cohorts_arr == g
        if layout is None:
            att_w[rows] = float(T - g) / max(total_trt_cells, 1.0)
        else:
            att_w[rows] = layout.window_len[g][rows].astype(np.float64) / max(total_trt_cells, 1.0)
    return att_w


def _project_panel_to_mode(
    panel_dm: np.ndarray,
    cohorts_arr: np.ndarray,
    ctrl: np.ndarray,
    active_waves: list[int],
    wave_weights: dict[int, float],
    bases: dict[int, np.ndarray],
    k: int,
) -> np.ndarray:
    """Project a double-demeaned panel ``(N, T)`` or draw stack ``(M, N, T)`` onto mode ``k``."""
    w_sum = max(sum(wave_weights[g] for g in active_waves), 1e-12)
    if panel_dm.ndim == 2:
        N, T = panel_dm.shape
        out = np.zeros(N, dtype=np.float64)
        for g in active_waves:
            psi_k = bases[g][:, k]
            mask_g = cohorts_arr == g
            L_g = float(T - g)
            out[mask_g] = (panel_dm[mask_g, g:] @ psi_k) / L_g
            w_norm = wave_weights[g] / w_sum
            out[ctrl] += w_norm * ((panel_dm[ctrl, g:] @ psi_k) / L_g)
        return out

    M, N, T = panel_dm.shape
    out_m = np.zeros((M, N), dtype=np.float64)
    for g in active_waves:
        psi_k = bases[g][:, k]
        mask_g = cohorts_arr == g
        L_g = float(T - g)
        out_m[:, mask_g] = np.einsum("mit,t->mi", panel_dm[:, mask_g, g:], psi_k) / L_g
        w_norm = wave_weights[g] / w_sum
        out_m[:, ctrl] += w_norm * (np.einsum("mit,t->mi", panel_dm[:, ctrl, g:], psi_k) / L_g)
    return out_m


def _compute_mode_precision_weights(
    cohorts_arr: np.ndarray,
    ctrl: np.ndarray,
    active_waves: list[int],
    wave_weights: dict[int, float],
    bases: dict[int, np.ndarray],
    t0: int,
    T: int,
    k: int,
) -> np.ndarray:
    """Compute relative precision weights ``w_{ik}`` inversely proportional to cell-projection variance."""
    N = cohorts_arr.shape[0]
    w_sum = max(sum(wave_weights[g] for g in active_waves), 1e-12)
    # Time weight vector on t = 0..T-1 for control units
    c_ctrl = np.zeros(T, dtype=np.float64)
    for g in active_waves:
        psi_k = bases[g][:, k]
        L_g = float(T - g)
        w_norm = wave_weights[g] / w_sum
        c_ctrl[g:] += w_norm * (psi_k / L_g)

    pre_var = (1.0 / float(max(t0, 1))) if k == 0 else 0.0
    var_ctrl = float(np.sum(c_ctrl ** 2)) + pre_var
    var_unit = np.full(N, max(var_ctrl, 1e-8), dtype=np.float64)
    for g in active_waves:
        mask_g = cohorts_arr == g
        L_g = float(T - g)
        # Since (1/L_g) sum_s psi_k(s)^2 == 1, sum_t (psi_k(t)/L_g)^2 == 1 / L_g
        var_g = (1.0 / L_g) + pre_var
        var_unit[mask_g] = max(var_g, 1e-8)

    inv_var = 1.0 / var_unit
    ctrl_mean = float(np.mean(inv_var[ctrl]))
    return inv_var / max(ctrl_mean, 1e-8)


def _fallback_posterior_forecast_draws(
    *,
    y_arr: np.ndarray,
    f_arr: np.ndarray,
    x_mod: np.ndarray,
    ctrl_idx: np.ndarray,
    t0: int,
    n_draws: int,
    seed: int,
) -> np.ndarray:
    """Construct leak-free NIG posterior forecast draws when raw `baseline_forecast` was passed without `features`."""
    N, T = y_arr.shape
    ctx = y_arr[:, :t0].copy()
    f_ctx = f_arr[:, :t0].copy()
    H_embed = np.column_stack([x_mod, ctx - f_ctx])
    y_ctrl = y_arr[ctrl_idx].copy()
    rng_ad = np.random.default_rng(int(seed) + 7907)
    _, adapter_fit = _cross_fit_control_adapter(
        y_ctrl=y_ctrl,
        ctx_mean_ctrl=np.mean(ctx[ctrl_idx], axis=1, keepdims=True),
        H_all=H_embed,
        ctrl_idx=ctrl_idx,
        zero_shot_all=f_arr,
        n_folds=5,
        rng=rng_ad,
        return_fit=True,
    )
    y_safe = y_arr.copy()
    ever_mask = np.ones(N, dtype=bool)
    ever_mask[ctrl_idx] = False
    y_safe[ever_mask, t0:] = np.nan
    draws, _ = _forecast_draws_posterior(
        n_draws=int(n_draws),
        context=ctx,
        unfiltered_forecast=f_arr,
        adapter_fit=adapter_fit,
        T=T,
        anchor_t=t0,
        lc_spec=resolve_late_context("fixed"),
        forecaster=None,
        y_safe=y_safe,
        ctrl_idx=ctrl_idx,
        rng=np.random.default_rng(int(seed) + 9973),
    )
    return draws


def _resolve_forecast_draws(
    *,
    n_draws: int,
    features: Any | None,
    baseline_forecast_draws: np.ndarray | None,
    y_arr: np.ndarray,
    f_arr: np.ndarray,
    x_mod: np.ndarray,
    ctrl_idx: np.ndarray,
    t0: int,
    seed: int,
    allow_fallback: bool = True,
) -> np.ndarray | None:
    """Obtain ``n_draws`` leak-free posterior forecast panels of shape ``(n_draws, N, T)``.

    ``allow_fallback=False`` (masked panels) skips the internal Normal-inverse-gamma
    generator, which needs complete control outcomes, and warns that no forecast
    draws are propagated.
    """
    if n_draws <= 1:
        return None
    if baseline_forecast_draws is not None:
        d_arr = np.asarray(baseline_forecast_draws, dtype=np.float64)
        if d_arr.shape[0] >= n_draws:
            return d_arr[:n_draws]
    if features is not None:
        pf_draws = getattr(features, "point_forecast_draws", None)
        f_diag = getattr(features, "forecast_draws", None) or {}
        if (
            pf_draws is not None
            and pf_draws.shape[0] >= n_draws
            and f_diag.get("draw_method") == "posterior"
        ):
            return np.asarray(pf_draws[:n_draws], dtype=np.float64)
        if hasattr(features, "draw_posterior_forecasts"):
            out = features.draw_posterior_forecasts(n_draws)
            if out is not None:
                return np.asarray(out[0], dtype=np.float64)
    if ctrl_idx.size >= 6 and t0 >= 1:
        if not allow_fallback:
            _warn_rank_k(
                "Masked panel: the internal forecast-draw generator needs complete outcomes, so no forecast draws "
                "are propagated (the sampler conditions on the point forecast). Pass baseline_forecast_draws, or "
                "features built on the complete common pre window, to propagate forecast uncertainty."
            )
            return None
        return _fallback_posterior_forecast_draws(
            y_arr=y_arr,
            f_arr=f_arr,
            x_mod=x_mod,
            ctrl_idx=ctrl_idx,
            t0=t0,
            n_draws=n_draws,
            seed=seed,
        )
    return None


def _project_mode1_tau_draws(
    tau_0: np.ndarray,
    *,
    cohorts_arr: np.ndarray,
    ever: np.ndarray,
    waves: list[int],
    att_unit_weights: np.ndarray,
    tau_clip: tuple[float, float],
    catt_floor: float | None,
) -> np.ndarray:
    """Enforce ``tau_clip`` and ``catt_floor`` on Mode-1 unit draws while preserving cohort means and K-invariance."""
    out = tau_0.copy()
    if np.isfinite(tau_clip[0]) or np.isfinite(tau_clip[1]):
        att_raw = np.einsum("ns,n->s", out, att_unit_weights)
        att_clipped = np.clip(att_raw, float(tau_clip[0]), float(tau_clip[1]))
        out[ever] += (att_clipped - att_raw)[None, :]
    if catt_floor is not None:
        floor_val = float(catt_floor)
        for g in waves:
            mask_g = cohorts_arr == g
            tau_g = out[mask_g]
            mean_g = tau_g.mean(axis=0, keepdims=True)
            dev_g = tau_g - mean_g
            mean_g_floored = np.maximum(mean_g, floor_val)
            max_neg_dev = np.max(np.maximum(-dev_g, 0.0), axis=0, keepdims=True)
            beta_g = np.where(
                max_neg_dev > 1e-12,
                np.minimum(1.0, np.maximum(mean_g_floored - floor_val, 0.0) / np.maximum(max_neg_dev, 1e-12)),
                1.0,
            )
            out[mask_g] = np.maximum(mean_g_floored + beta_g * dev_g, floor_val)
    return out


def _assemble_catt_draws(
    mode_tau_draws: list[np.ndarray],
    *,
    cohorts_arr: np.ndarray,
    waves: list[int],
    bases: dict[int, np.ndarray],
    n_units: int,
    n_periods: int,
    catt_floor: float | None,
) -> np.ndarray:
    """Assemble ``(S, N, T)`` cell treatment effects with zero-sum higher-mode amplitude damping under ``catt_floor``."""
    n_total_draws = mode_tau_draws[0].shape[1]
    catt_level = np.zeros((n_total_draws, n_units, n_periods), dtype=np.float64)
    catt_dyn = np.zeros((n_total_draws, n_units, n_periods), dtype=np.float64)
    for k, tau_k in enumerate(mode_tau_draws):
        active_waves = [g for g in waves if k < bases[g].shape[1]]
        for g in active_waves:
            psi_k = bases[g][:, k]
            mask_g = cohorts_arr == g
            contrib = tau_k[mask_g].T[:, :, None] * psi_k[None, None, :]
            if k == 0:
                catt_level[:, mask_g, g : g + psi_k.shape[0]] += contrib
            else:
                catt_dyn[:, mask_g, g : g + psi_k.shape[0]] += contrib
    if catt_floor is None:
        return catt_level + catt_dyn
    floor_val = float(catt_floor)
    catt_draws = catt_level.copy()
    for g in waves:
        mask_g = cohorts_arr == g
        end = g + bases[g].shape[0]
        level_g = catt_level[:, mask_g, g:end]
        dyn_g = catt_dyn[:, mask_g, g:end]
        max_neg_dyn = np.max(np.maximum(-dyn_g, 0.0), axis=2, keepdims=True)
        headroom = np.maximum(level_g - floor_val, 0.0)
        alpha_g = np.where(
            max_neg_dyn > 1e-12,
            np.minimum(1.0, headroom / np.maximum(max_neg_dyn, 1e-12)),
            1.0,
        )
        catt_draws[:, mask_g, g:end] = np.maximum(level_g + alpha_g * dyn_g, floor_val)
    return catt_draws


# ----------------------------------------------------------------------------------------------
# Remediation helpers: short-cohort handling, projection accounting, prior diagnostics
# ----------------------------------------------------------------------------------------------


def _mode_active_mask(cohorts_arr: np.ndarray, ctrl: np.ndarray, active_waves: list[int]) -> np.ndarray:
    """Units taking part in a mode's fit: every control plus the treated cohorts that have the mode."""
    active = np.asarray(ctrl, dtype=bool).copy()
    for g in active_waves:
        active |= cohorts_arr == g
    return active


def _resolve_short_cohorts(
    cfg: RankKLongBetConfig, waves: list[int], T: int, warn: bool = True
) -> int:
    """Apply ``cfg.short_cohort_policy`` to the cohorts with ``L_g < rank_k``; return the effective ``K``.

    Raises under ``"error"``; warns (when ``warn``) under ``"reduce_k"`` and ``"drop_mode"``.
    """
    K = int(cfg.rank_k)
    short = [g for g in waves if T - g < K]
    if not short:
        return K
    lengths = [T - g for g in short]
    if cfg.short_cohort_policy == "error":
        raise ValueError(
            f"rank_k={K} exceeds the post-adoption length L_g of {len(short)} treated cohort(s) "
            f"(first-adoption period(s) {short} with L_g = {lengths}); set short_cohort_policy to "
            "'drop_mode' (fit each higher mode on the cohorts that have it) or 'reduce_k'"
        )
    if cfg.short_cohort_policy == "reduce_k":
        k_eff = max(1, min(K, min(T - g for g in waves)))
        if warn:
            _warn_rank_k(
                f"short_cohort_policy='reduce_k': rank_k lowered from {K} to {k_eff} because cohort(s) {short} "
                f"have only L_g = {lengths} post-adoption period(s). The ATT and unit-averaged effects are "
                "unchanged (Proposition 2)."
            )
        return k_eff
    if warn:
        _warn_rank_k(
            f"rank_k={K} exceeds the post-adoption length L_g of {len(short)} treated cohort(s) (first-adoption "
            f"period(s) {short} with L_g = {lengths}). short_cohort_policy='drop_mode': these units take part in "
            "modes k <= L_g only and are excluded from the fit of the higher modes (never padded with zero "
            "outcomes). Their unit-averaged effects and the ATT are unaffected."
        )
    return K


def _scatter_rows(values: np.ndarray, idx: np.ndarray, n_units: int, axis: int = 0) -> np.ndarray:
    """Scatter an array whose ``axis`` indexes ``idx`` back onto ``n_units`` units (zeros elsewhere)."""
    shape = list(values.shape)
    shape[axis] = n_units
    out = np.zeros(shape, dtype=values.dtype)
    index = [slice(None)] * values.ndim
    index[axis] = idx
    out[tuple(index)] = values
    return out


def _projection_report(
    *,
    mode1_raw: np.ndarray,
    mode1_proj: np.ndarray,
    catt_raw: np.ndarray,
    catt_proj: np.ndarray,
    att_unit_weights: np.ndarray,
    z_arr: np.ndarray,
    tau_clip: tuple[float, float],
    catt_floor: float | None,
) -> dict[str, Any]:
    """Quantify a post-processing projection of the posterior draws (``catt_floor`` / ``tau_clip`` guard)."""
    att_raw = np.einsum("ns,n->s", mode1_raw, att_unit_weights)
    att_proj = np.einsum("ns,n->s", mode1_proj, att_unit_weights)
    trt_cells = np.asarray(z_arr) == 1
    d_mode1 = np.abs(mode1_proj - mode1_raw)[att_unit_weights > 0]
    d_catt = np.abs(catt_proj - catt_raw)[:, trt_cells]
    report: dict[str, Any] = {
        "tau_clip": (float(tau_clip[0]), float(tau_clip[1])),
        "catt_floor": None if catt_floor is None else float(catt_floor),
        "att_shift_max_abs": float(np.max(np.abs(att_proj - att_raw))) if att_raw.size else 0.0,
        "att_mean_shift": float(np.mean(att_proj - att_raw)) if att_raw.size else 0.0,
        "att_sd_ratio": float(np.std(att_proj, ddof=1) / max(np.std(att_raw, ddof=1), 1e-12)) if att_raw.size > 1 else 1.0,
        "att_draws_outside_clip_share": float(
            np.mean((att_raw < float(tau_clip[0]) - 1e-9) | (att_raw > float(tau_clip[1]) + 1e-9))
        ),
        "mode1_draws_altered_share": float(np.mean(d_mode1 > 1e-12)) if d_mode1.size else 0.0,
        "mode1_shift_mean_abs": float(np.mean(d_mode1)) if d_mode1.size else 0.0,
        "mode1_shift_max_abs": float(np.max(d_mode1)) if d_mode1.size else 0.0,
        "cell_draws_altered_share": float(np.mean(d_catt > 1e-12)) if d_catt.size else 0.0,
        "cell_shift_mean_abs": float(np.mean(d_catt)) if d_catt.size else 0.0,
    }
    if catt_floor is not None:
        report["unprojected_mode1_mass_below_floor"] = float(np.mean(mode1_raw[att_unit_weights > 0] < float(catt_floor)))
        report["unprojected_cell_mass_below_floor"] = float(np.mean(catt_raw[:, trt_cells] < float(catt_floor)))
        report["cells_at_floor_share"] = float(np.mean(np.abs(catt_proj[:, trt_cells] - float(catt_floor)) <= 1e-12))
    return report


def prior_data_conflict_report(
    res: "RankKLongBetResult",
    *,
    unconstrained: "RankKLongBetResult | None" = None,
    warn: bool = True,
) -> dict[str, Any]:
    """Summarize how much posterior mass a prior / truncation / floor removed (prior-data conflict check).

    ``catt_floor`` is a post-sampling projection, so its removed mass is read from the
    fit's own unprojected draws. A ``tau_clip`` truncation or a ``tau_prior_*`` prior
    is enforced *inside* the sampler, so measuring what it removed needs the matching
    ``unconstrained`` fit (same data and budget, no prior / clip): the report then
    gives the share of its ATT draws outside ``tau_clip``, the 30-70% interval
    overlap of the two ATT posteriors and the ratio of their posterior SDs. Warns
    (``RankKConfigWarning``) when more than half of the unconstrained ATT mass lies
    outside ``tau_clip``, when the two central intervals do not overlap, or when more
    than half of the unprojected unit-level mass lies below ``catt_floor``.
    """
    rep = dict(res.projection_report or {})
    out: dict[str, Any] = {
        "projection_applied": bool(rep),
        "mode1_mass_below_floor": rep.get("unprojected_mode1_mass_below_floor", 0.0),
        "cell_mass_below_floor": rep.get("unprojected_cell_mass_below_floor", 0.0),
        "att_mass_outside_clip": None,
        "central_interval_overlap": None,
        "posterior_sd_ratio": None,
    }
    lo, hi = (float(rep.get("tau_clip", (-np.inf, np.inf))[0]), float(rep.get("tau_clip", (-np.inf, np.inf))[1]))
    out["att_mean_shift_in_sd"] = (
        float(abs(rep.get("att_mean_shift", 0.0)) / max(float(np.std(np.asarray(res.att_draws), ddof=1)), 1e-12))
        if rep else 0.0
    )
    conflict = out["mode1_mass_below_floor"] > 0.5 or out["att_mean_shift_in_sd"] > 0.5
    if unconstrained is not None:
        u = np.asarray(unconstrained.att_draws, dtype=np.float64)
        c = np.asarray(res.att_draws, dtype=np.float64)
        out["att_mass_outside_clip"] = float(np.mean((u < lo) | (u > hi)))
        u30, u70 = np.percentile(u, [30, 70])
        c30, c70 = np.percentile(c, [30, 70])
        out["central_interval_overlap"] = bool(min(u70, c70) >= max(u30, c30))
        out["posterior_sd_ratio"] = float(np.std(c, ddof=1) / max(np.std(u, ddof=1), 1e-12))
        conflict = conflict or out["att_mass_outside_clip"] > 0.5 or not out["central_interval_overlap"]
    out["conflict"] = bool(conflict)
    if conflict and warn:
        _warn_rank_k(
            "Prior-data conflict: the constraint removed most of the unconstrained posterior mass "
            f"(ATT mass outside tau_clip: {out['att_mass_outside_clip']}; central intervals overlap: "
            f"{out['central_interval_overlap']}; unit-level mass below catt_floor: "
            f"{out['mode1_mass_below_floor']:.0%}). Report the unconstrained fit alongside."
        )
    return out


def gaussian_summary_prior_update(
    draws: np.ndarray, prior_mean: float, prior_sd: float
) -> tuple[np.ndarray, dict[str, float]]:
    """Post-processing normal-approximation prior update of a scalar posterior sample (sensitivity use only).

    Treats the draws as a Gaussian likelihood summary ``N(m, v)`` and combines it with
    ``N(prior_mean, prior_sd^2)`` by the exact conjugate affine map: posterior mean
    ``(m / v + mu_0 / s_0^2) / (1 / v + 1 / s_0^2)`` and SD factor ``sqrt(v' / v)`` (not
    the precision-weight ``lambda`` of the removed plug-in code). It double counts the
    estimator's own weak prior and is outside every coverage guarantee: use it only to
    tabulate prior sensitivity, never as the reported posterior.
    """
    d = np.asarray(draws, dtype=np.float64)
    m, v = float(np.mean(d)), float(np.var(d, ddof=1))
    prec = 1.0 / max(v, 1e-18) + 1.0 / max(float(prior_sd) ** 2, 1e-18)
    v_new = 1.0 / prec
    m_new = v_new * (m / max(v, 1e-18) + float(prior_mean) / max(float(prior_sd) ** 2, 1e-18))
    out = m_new + np.sqrt(v_new / max(v, 1e-18)) * (d - m)
    return out, {"mean": m_new, "sd": float(np.sqrt(v_new)), "shrinkage": float(v_new / max(v, 1e-18))}


def _mode1_diagnostics(chains_proj: np.ndarray, compute_rhat: Any, compute_ess: Any) -> dict[str, Any]:
    """Per-unit split R-hat and bulk / tail ESS of the mode-1 draws (chains, draws, units)."""
    C, S, N = chains_proj.shape
    rhat = np.full(N, np.nan)
    ess_bulk = np.full(N, np.nan)
    ess_tail = np.full(N, np.nan)
    if C >= 1 and S >= 4:
        for i in range(N):
            x = chains_proj[:, :, i]
            if float(np.std(x)) < 1e-15:
                continue
            try:
                rhat[i] = float(compute_rhat(x))
                ess_bulk[i] = float(compute_ess(x, method="bulk"))
                ess_tail[i] = float(compute_ess(x, method="tail"))
            except (ArithmeticError, ValueError, np.linalg.LinAlgError):  # pragma: no cover - diagnostic only
                pass
    return {"rhat_unit": rhat, "ess_bulk_unit": ess_bulk, "ess_tail_unit": ess_tail}



def _fit_rank_k_stochtree(
    *,
    y_arr: np.ndarray,
    z_arr: np.ndarray,
    x_mod: np.ndarray,
    ps_vec: np.ndarray,
    dy_dm: np.ndarray,
    df_dm: np.ndarray,
    has_forecast: bool,
    cohorts_arr: np.ndarray,
    ever: np.ndarray,
    ctrl: np.ndarray,
    waves: list[int],
    wave_weights: dict[int, float],
    bases: dict[int, np.ndarray],
    cfg: RankKLongBetConfig,
    t_start: float,
    att_w: np.ndarray | None = None,
    result_fields: dict[str, Any] | None = None,
) -> RankKLongBetResult:
    """Legacy plug-in ``stochtree.BCFModel`` implementation (``engine='stochtree'``).

    ``att_w`` are the ATT weights (default: the paper's ``L_{G_i} / |O_trt|``) and
    ``result_fields`` extra :class:`RankKLongBetResult` fields of the caller.
    """
    try:
        from stochtree import BCFModel
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "The plug-in variant of fit_rank_k_longbet (engine 'stochtree', "
            "uncertainty_propagation=False) requires 'stochtree', a core dependency of "
            "longbet-jax >= 0.5.0; reinstall with  pip install longbet-jax  (or use the "
            "default propagated variant with RankKLongBetConfig(uncertainty_propagation=True))"
        ) from exc

    N, T = y_arr.shape
    K = int(cfg.rank_k)
    if cfg.has_att_constraint:
        raise ValueError(
            "tau_prior_sd, tau_clip and tau_prior_family require engine='jax' (the propagated variant), where "
            "the prior and the truncation act on the ATT inside the Gibbs sampler. The plug-in stochtree BCF "
            "has no intercept parameter to carry them, and the former post-hoc shrinkage of its draws was not "
            "a valid posterior update (it double-counted the prior and contracted the SD by lambda instead of "
            "sqrt(lambda)). Use gaussian_summary_prior_update() for an explicit, labeled sensitivity table."
        )
    n_draws = int(cfg.num_mcmc)
    mode_tau_draws: list[np.ndarray] = []
    theta_f_modes: list[float] = []
    theta_f_draws: list[np.ndarray] = []
    sigma_modes: list[float] = []
    gp_hyper_draws: list[dict[str, np.ndarray]] = []
    active_units: list[np.ndarray] = []
    resid_means: list[np.ndarray] = []
    resid_vars: list[np.ndarray] = []
    w_units_st: list[np.ndarray] = []
    bayes_var_ctrl: list[float] = []
    t0_st = min(waves) if waves else 1
    trt_keep_vars = list(range(x_mod.shape[1]))
    total_trt_cells_st = float(np.sum(z_arr == 1))
    if att_w is None:
        att_w = np.zeros(N, dtype=np.float64)
        for g in waves:
            att_w[cohorts_arr == g] = float(T - g) / max(total_trt_cells_st, 1.0)
    else:
        att_w = np.asarray(att_w, dtype=np.float64)

    for k in range(K):
        active_waves = [g for g in waves if k < bases[g].shape[1]]
        if not active_waves:
            break
        dy_k = _project_panel_to_mode(dy_dm, cohorts_arr, ctrl, active_waves, wave_weights, bases, k)
        df_k = _project_panel_to_mode(df_dm, cohorts_arr, ctrl, active_waves, wave_weights, bases, k)
        w_unit_k = _compute_mode_precision_weights(
            cohorts_arr, ctrl, active_waves, wave_weights, bases, t0_st, T, k
        )
        # Units that have mode k: every control and the treated cohorts with L_g >= k + 1
        # (short cohorts are excluded from the fit instead of being padded with zero outcomes).
        active = _mode_active_mask(cohorts_arr, ctrl, active_waves)
        active_units.append(active)
        w_units_st.append(np.where(active, w_unit_k, 0.0))
        idx = np.nonzero(active)[0]
        full = bool(np.all(active))
        dy_k_a = dy_k if full else dy_k[idx]
        df_k_a = df_k if full else df_k[idx]
        ctrl_a = ctrl if full else ctrl[idx]
        ever_a = ever if full else ever[idx]
        x_mod_a = x_mod if full else x_mod[idx]
        ps_a = ps_vec if full else ps_vec[idx]

        if has_forecast:
            df_kc = df_k_a[ctrl_a] - df_k_a[ctrl_a].mean()
            dy_kc = dy_k_a[ctrl_a] - dy_k_a[ctrl_a].mean()
            var_f = float(np.sum(df_kc ** 2))
            if var_f > 1e-8:
                raw_theta = float(np.sum(df_kc * dy_kc) / (var_f + float(cfg.ridge_eps)))
                theta_fk = float(np.clip(raw_theta, cfg.theta_f_clip[0], cfg.theta_f_clip[1]))
            else:
                theta_fk = 0.0
            y_k = dy_k_a - theta_fk * df_k_a
            df_std = (df_k_a - df_k_a.mean()) / max(float(df_k_a.std()), 1e-6)
            x_k = np.column_stack([x_mod_a, df_std])
        else:
            theta_fk = 0.0
            y_k = dy_k_a
            x_k = x_mod_a

        theta_f_modes.append(theta_fk)
        theta_f_draws.append(np.full(n_draws, theta_fk, dtype=np.float64))
        n_trees_tau = int(cfg.num_trees_tau_mode1 if k == 0 else cfg.num_trees_tau_higher)
        model = BCFModel()
        model.sample(
            X_train=x_k,
            Z_train=ever_a.astype(np.float64),
            y_train=y_k,
            propensity_train=ps_a,
            num_gfr=int(cfg.num_gfr),
            num_burnin=int(cfg.num_burnin),
            num_mcmc=n_draws,
            general_params={
                "random_seed": int(cfg.random_seed + 101 * k),
                "keep_gfr": False,
                "num_threads": int(cfg.num_threads),
            },
            prognostic_forest_params={"num_trees": int(cfg.num_trees_mu)},
            treatment_effect_forest_params={
                "num_trees": n_trees_tau,
                "keep_vars": trt_keep_vars,
            },
        )
        pred_k = model.predict(
            X=x_k, Z=ever_a.astype(np.float64), propensity=ps_a, type="posterior", terms=["tau", "y_hat"]
        )
        if isinstance(pred_k, dict):
            tau_k = pred_k.get("tau_hat", pred_k.get("tau"))
            yhat_k = pred_k.get("y_hat")
        else:
            tau_k = pred_k
            yhat_k = None
        tau_k = np.asarray(tau_k, dtype=np.float64)
        if tau_k.ndim == 3:
            tau_k = tau_k[:, :, 0]
        if yhat_k is not None:
            yhat_arr = np.asarray(yhat_k, dtype=np.float64)
            if yhat_arr.ndim == 3:
                yhat_arr = yhat_arr[:, :, 0]
            resid_draws_k = y_k[:, None] - yhat_arr
            r_mean_a = resid_draws_k.mean(axis=1)
            r_var_a = resid_draws_k.var(axis=1, ddof=1) if resid_draws_k.shape[1] > 1 else np.zeros_like(y_k)
        else:
            r_mean_a = y_k - ever_a.astype(np.float64) * tau_k.mean(axis=1)
            r_var_a = np.zeros_like(y_k)
        if not full:
            tau_k = _scatter_rows(tau_k, idx, N, axis=0)
            r_mean_a = _scatter_rows(r_mean_a, idx, N, axis=0)
            r_var_a = _scatter_rows(r_var_a, idx, N, axis=0)
        resid_means.append(r_mean_a)
        resid_vars.append(r_var_a)
        mode_tau_draws.append(tau_k)
        gvar_mean = float(np.mean(model.global_var_samples))
        bayes_var_ctrl.append(gvar_mean)
        sigma_modes.append(float(np.sqrt(gvar_mean)))
        gp_hyper_draws.append(
            {
                "s_zeta": np.zeros(n_draws, dtype=np.float64),
                "ell": np.zeros(n_draws, dtype=np.float64),
            }
        )

    compute_rhat, compute_ess = _get_rhat_ess()
    t_post_st = int((result_fields or {}).get("horizon") or (T - t0_st))
    arm_meta_inputs = {
        "y_arr": y_arr,
        "f_arr": np.zeros_like(y_arr),
        "has_forecast": False,
        "dy_dm": dy_dm,
        "df_dm": df_dm,
        "x_mod": x_mod,
        "obs": np.isfinite(y_arr),
        "ctrl": ctrl,
        "t0": t0_st,
        "t_post": t_post_st,
        "resid_means": resid_means,
        "resid_vars": resid_vars,
        "w_units": w_units_st,
        "bayes_var_ctrl": bayes_var_ctrl,
        "bayes_var_trt": bayes_var_ctrl,
        "treated_extra_variance": bool(cfg.treated_extra_variance),
    }
    return _finalize_result(
        mode_tau_draws_raw=mode_tau_draws,
        cohorts_arr=cohorts_arr,
        ever=ever,
        waves=waves,
        bases=bases,
        att_unit_weights=att_w,
        z_arr=z_arr,
        cfg=cfg,
        theta_f_modes=theta_f_modes,
        sigma_modes=sigma_modes,
        theta_f_draws=theta_f_draws,
        gp_hyper_draws=gp_hyper_draws,
        att_var_within=None,
        att_var_between=0.0,
        att_var_slope=0.0,
        n_components=1,
        rhat_att=None,
        ess_att=None,
        diagnostics={"engine": "stochtree"},
        variant="plug-in",
        active_units=active_units,
        t_start=t_start,
        compute_rhat=compute_rhat,
        compute_ess=compute_ess,
        result_fields=result_fields,
        arm_meta_inputs=arm_meta_inputs,
    )


def _extract_event_time_trt_trajectory(
    catt_draws: np.ndarray,
    ever: np.ndarray,
    cohorts_arr: np.ndarray,
    waves: list[int],
    t0: int,
    t_post: int,
) -> np.ndarray:
    """Return ``(S, t_post)`` treated-unit-average trajectory draws over event time ``1..t_post``."""
    n_draws, _, n_periods = catt_draws.shape
    t_post_eff = max(int(t_post), 1)
    if len(waves) <= 1:
        end_t = min(t0 + t_post_eff, n_periods)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            traj = np.nanmean(catt_draws[:, ever, t0:end_t], axis=1)
        if traj.shape[1] < t_post_eff:
            pad = np.repeat(traj[:, -1:], t_post_eff - traj.shape[1], axis=1)
            traj = np.concatenate([traj, pad], axis=1)
        return np.asarray(traj, dtype=np.float64)

    traj = np.zeros((n_draws, t_post_eff), dtype=np.float64)
    trt_idx = np.nonzero(ever)[0]
    trt_cohorts = cohorts_arr[trt_idx]
    for s in range(t_post_eff):
        cal_t = trt_cohorts + s
        valid = cal_t < n_periods
        if np.any(valid):
            vals = catt_draws[:, trt_idx[valid], cal_t[valid]]
            traj[:, s] = np.nanmean(vals, axis=1)
        elif s > 0:
            traj[:, s] = traj[:, s - 1]
    return traj


def _psd_clip_matrix(mat: np.ndarray, floor: float = 0.0) -> np.ndarray:
    """Symmetrize ``mat`` and clip eigenvalues below ``floor``."""
    sym = 0.5 * (np.asarray(mat, dtype=np.float64) + np.asarray(mat, dtype=np.float64).T)
    if sym.size == 0:
        return sym
    eigvals, eigvecs = np.linalg.eigh(sym)
    eigvals = np.maximum(eigvals, float(floor))
    out = (eigvecs * eigvals[None, :]) @ eigvecs.T
    return 0.5 * (out + out.T)


def _pairwise_mean_cov(u_arm: np.ndarray) -> np.ndarray:
    """Estimate covariance of the column-wise sample mean across observed cells of ``u_arm``."""
    _, h_len = u_arm.shape
    if u_arm.shape[0] <= 1 or h_len == 0:
        return np.zeros((h_len, h_len), dtype=np.float64)
    finite = np.isfinite(u_arm)
    counts = finite.sum(axis=0).astype(np.float64)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        col_means = np.nanmean(u_arm, axis=0, keepdims=True)
    col_means = np.where(np.isfinite(col_means), col_means, 0.0)
    centered = np.where(finite, u_arm - col_means, 0.0)
    pair_n = finite.astype(np.float64).T @ finite.astype(np.float64)
    denom_cov = np.maximum(pair_n - 1.0, 1.0)
    raw_cov = (centered.T @ centered) / denom_cov
    raw_cov = np.where(pair_n >= 2.0, raw_cov, 0.0)
    n_eff = np.maximum(np.maximum(counts[:, None], counts[None, :]), 1.0)
    return _psd_clip_matrix(raw_cov / n_eff, floor=0.0)


def _build_stage1_arm_metadata(
    *,
    catt_draws: np.ndarray,
    catt_hat: np.ndarray,
    ever: np.ndarray,
    cohorts_arr: np.ndarray,
    waves: list[int],
    bases: dict[int, np.ndarray],
    theta_f_modes: list[float],
    arm_meta_inputs: dict[str, Any],
) -> dict[str, Any]:
    """Assemble lightweight Stage-1 arm metadata for Bayesian multi-arm reconciliation (§3)."""
    y_arr = np.asarray(arm_meta_inputs["y_arr"], dtype=np.float64)
    f_arr = np.asarray(arm_meta_inputs["f_arr"], dtype=np.float64)
    has_forecast = bool(arm_meta_inputs["has_forecast"])
    dy_dm = np.asarray(arm_meta_inputs["dy_dm"], dtype=np.float64)
    df_dm = np.asarray(arm_meta_inputs["df_dm"], dtype=np.float64)
    obs = np.asarray(arm_meta_inputs["obs"], dtype=bool)
    ctrl = np.asarray(arm_meta_inputs["ctrl"], dtype=bool)
    t0 = int(arm_meta_inputs["t0"])
    t_post = max(int(arm_meta_inputs["t_post"]), 1)
    resid_means: list[np.ndarray] = list(arm_meta_inputs["resid_means"])
    resid_vars: list[np.ndarray] = list(arm_meta_inputs["resid_vars"])
    w_units: list[np.ndarray] = list(arm_meta_inputs["w_units"])
    bayes_var_ctrl: list[float] = list(arm_meta_inputs["bayes_var_ctrl"])
    bayes_var_trt: list[float] = list(arm_meta_inputs["bayes_var_trt"])

    x_mod_raw = arm_meta_inputs.get("x_mod")
    x_mod = np.asarray(x_mod_raw, dtype=np.float64) if x_mod_raw is not None else None
    if x_mod is not None and x_mod.ndim == 1:
        x_mod = x_mod[:, None]
    tev_enabled = bool(arm_meta_inputs.get("treated_extra_variance", True))

    n_units, n_periods = y_arr.shape
    theta_0 = float(theta_f_modes[0]) if len(theta_f_modes) > 0 else 0.0

    # Pre-period unit residual variance per arm
    pre_u = y_arr[:, :t0] - (theta_0 * f_arr[:, :t0] if has_forecast else 0.0)
    pre_masked = np.where(obs[:, :t0], pre_u, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        unit_pre_mean = np.nanmean(pre_masked, axis=1)
        unit_pre_mean = np.where(np.isfinite(unit_pre_mean), unit_pre_mean, 0.0)
        ctrl_pre_col = np.nanmean(pre_masked[ctrl], axis=0, keepdims=True)
        ctrl_pre_col = np.where(np.isfinite(ctrl_pre_col), ctrl_pre_col, 0.0)
        pre_dm = pre_masked - ctrl_pre_col
        if t0 >= 2:
            pre_within = pre_dm - np.nanmean(pre_dm, axis=1, keepdims=True)
            df_scale = float(t0) / float(t0 - 1)
        else:
            pre_within = pre_dm - np.nanmean(pre_dm, axis=0, keepdims=True)
            df_scale = 1.0

    def _arm_pre_var(mask_a: np.ndarray) -> float:
        vals = pre_within[mask_a]
        fin = vals[np.isfinite(vals)]
        return float(np.var(fin, ddof=1) * df_scale) if fin.size > 1 else 0.0

    pre_var_trt = _arm_pre_var(ever)
    pre_var_ctrl = _arm_pre_var(ctrl)

    # Reference exposure basis on event time 1..t_post
    k_eff = len(resid_means)
    if t0 in bases and bases[t0].shape[0] >= t_post and k_eff >= 1:
        b_post = np.asarray(bases[t0][:t_post, : min(k_eff, bases[t0].shape[1])], dtype=np.float64)
    elif k_eff >= 1:
        b_post = orthonormal_exposure_basis(t_post, k_eff)[:, :k_eff]
    else:
        b_post = np.zeros((t_post, 0), dtype=np.float64)
    k_use = int(b_post.shape[1])

    # Post-adoption outcomes adjusted for forecast and posterior mean treatment effect
    y_adj = y_arr - (theta_0 * f_arr if has_forecast else 0.0) - np.where(ever[:, None], catt_hat, 0.0)
    y_adj_obs = np.where(obs, y_adj, np.nan)

    def _event_time_post_rows(mask_a: np.ndarray, is_trt: bool) -> tuple[np.ndarray, np.ndarray]:
        idx_a = np.nonzero(mask_a)[0]
        if not is_trt or len(waves) <= 1:
            end_t = min(t0 + t_post, n_periods)
            sub = y_adj_obs[idx_a, t0:end_t]
            if sub.shape[1] < t_post:
                pad = np.full((idx_a.size, t_post - sub.shape[1]), np.nan, dtype=np.float64)
                sub = np.concatenate([sub, pad], axis=1)
            return idx_a, sub
        out_u = np.full((idx_a.size, t_post), np.nan, dtype=np.float64)
        for r_i, u_i in enumerate(idx_a):
            g_i = int(cohorts_arr[u_i])
            len_i = min(t_post, max(n_periods - g_i, 0))
            if len_i > 0:
                out_u[r_i, :len_i] = y_adj_obs[u_i, g_i : g_i + len_i]
        return idx_a, out_u

    def _arm_unit_innovations(mask_a: np.ndarray, is_trt: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """Project out within-arm prognostic predictors ``[1, x_mod, y_pre_mean]`` from post-period cells."""
        idx_a, post_a = _event_time_post_rows(mask_a, is_trt=is_trt)
        n_a = int(idx_a.size)
        finite_a = np.isfinite(post_a)
        counts_a = np.sum(finite_a, axis=0).astype(np.float64)
        if n_a <= 1 or t_post == 0:
            return np.zeros((n_a, t_post), dtype=np.float64), finite_a, counts_a, 1

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            col_means = np.nanmean(post_a, axis=0, keepdims=True)
        col_means = np.where(np.isfinite(col_means), col_means, 0.0)
        post_fill = np.where(finite_a, post_a, col_means)

        des_cols = [np.ones((n_a, 1), dtype=np.float64)]
        if n_a >= 4:
            des_cols.append(unit_pre_mean[idx_a, None])
        if x_mod is not None and x_mod.shape[1] > 0 and n_a >= x_mod.shape[1] + 6:
            des_cols.append(x_mod[idx_a])
        des_a = np.hstack(des_cols)
        coef_a, _, rank_d, _ = np.linalg.lstsq(des_a, post_fill, rcond=None)
        err_a = np.where(finite_a, post_a - des_a @ coef_a, 0.0)
        return err_a, finite_a, counts_a, max(int(rank_d), 1)

    err_trt, fin_trt, cnt_trt, rank_trt = _arm_unit_innovations(ever, is_trt=True)
    err_ctrl, fin_ctrl, cnt_ctrl, rank_ctrl = _arm_unit_innovations(ctrl, is_trt=False)

    pair_n_trt = fin_trt.astype(np.float64).T @ fin_trt.astype(np.float64)
    pair_n_ctrl = fin_ctrl.astype(np.float64).T @ fin_ctrl.astype(np.float64)
    denom_trt = np.sqrt(np.maximum(cnt_trt[:, None], 1.0) * np.maximum(cnt_trt[None, :], 1.0))
    denom_ctrl = np.sqrt(np.maximum(cnt_ctrl[:, None], 1.0) * np.maximum(cnt_ctrl[None, :], 1.0))

    if tev_enabled:
        unit_cov_ctrl = np.where(
            pair_n_ctrl >= 2.0,
            (err_ctrl.T @ err_ctrl) / np.maximum(pair_n_ctrl - float(rank_ctrl), 1.0),
            0.0,
        )
        unit_cov_trt = np.where(
            pair_n_trt >= 2.0,
            (err_trt.T @ err_trt) / np.maximum(pair_n_trt - float(rank_trt), 1.0),
            unit_cov_ctrl,
        )
        cov_ctrl = _psd_clip_matrix(unit_cov_ctrl / denom_ctrl, floor=0.0)
        cov_trt = _psd_clip_matrix(unit_cov_trt / denom_trt, floor=0.0)
    else:
        pair_n_pool = pair_n_trt + pair_n_ctrl
        df_pool = np.maximum(pair_n_pool - float(rank_trt + rank_ctrl), 1.0)
        unit_cov_pool = np.where(
            pair_n_pool >= 2.0,
            (err_trt.T @ err_trt + err_ctrl.T @ err_ctrl) / df_pool,
            0.0,
        )
        cov_ctrl = _psd_clip_matrix(unit_cov_pool / denom_ctrl, floor=0.0)
        cov_trt = _psd_clip_matrix(unit_cov_pool / denom_trt, floor=0.0)

    trt_traj_draws = _extract_event_time_trt_trajectory(catt_draws, ever, cohorts_arr, waves, t0, t_post)
    return {
        "n_treated": int(np.sum(ever)),
        "n_control": int(np.sum(ctrl)),
        "t0": t0,
        "t_post": t_post,
        "post_residual_cov_treated": cov_trt,
        "post_residual_cov_control": cov_ctrl,
        "pre_residual_var_treated": pre_var_trt,
        "pre_residual_var_control": pre_var_ctrl,
        "trt_trajectory_draws": trt_traj_draws,
    }


def _finalize_result(
    *,
    mode_tau_draws_raw: list[np.ndarray],
    cohorts_arr: np.ndarray,
    ever: np.ndarray,
    waves: list[int],
    bases: dict[int, np.ndarray],
    att_unit_weights: np.ndarray,
    z_arr: np.ndarray,
    cfg: RankKLongBetConfig,
    theta_f_modes: list[float],
    sigma_modes: list[float],
    theta_f_draws: list[np.ndarray],
    gp_hyper_draws: list[dict[str, np.ndarray]],
    att_var_within: float | None,
    att_var_between: float,
    att_var_slope: float,
    n_components: int,
    rhat_att: float | None,
    ess_att: float | None,
    diagnostics: dict[str, Any],
    variant: str,
    active_units: list[np.ndarray],
    t_start: float,
    compute_rhat: Any,
    compute_ess: Any,
    result_fields: dict[str, Any] | None = None,
    treated_cells: np.ndarray | None = None,
    arm_meta_inputs: dict[str, Any] | None = None,
    arm_metadata_override: dict[str, Any] | None = None,
) -> RankKLongBetResult:
    """Apply the (declared) post-processing projections, assemble cell draws and build the result.

    ``mode_tau_draws_raw`` are the sampler's draws. Mode 1 goes through
    :func:`_project_mode1_tau_draws` (``tau_clip`` guard, a no-op up to float32
    rounding since the sampler enforces the truncation itself, and the
    ``catt_floor`` projection); the cell draws are assembled with the zero-sum
    higher-mode damping of ``catt_floor``. Whenever a projection is active the
    unprojected draws are kept (``cfg.keep_unprojected_draws``) and quantified in
    ``projection_report``. ``result_fields`` are extra :class:`RankKLongBetResult`
    fields (the masked event-time additions) and ``treated_cells`` the ``(N, T)``
    cells the projection report averages over (default: ``z == 1``).
    """
    N, T = np.asarray(z_arr).shape
    mode1_raw = mode_tau_draws_raw[0]
    projecting = cfg.catt_floor is not None or math.isfinite(float(cfg.tau_clip[0])) or math.isfinite(float(cfg.tau_clip[1]))
    mode1_proj = _project_mode1_tau_draws(
        mode1_raw,
        cohorts_arr=cohorts_arr,
        ever=ever,
        waves=waves,
        att_unit_weights=att_unit_weights,
        tau_clip=cfg.tau_clip,
        catt_floor=cfg.catt_floor,
    )
    mode_tau_draws = [mode1_proj, *mode_tau_draws_raw[1:]]
    catt_draws = _assemble_catt_draws(
        mode_tau_draws, cohorts_arr=cohorts_arr, waves=waves, bases=bases, n_units=N, n_periods=T,
        catt_floor=cfg.catt_floor,
    )
    catt_hat = catt_draws.mean(axis=0)
    catt_ci_low, catt_ci_high = np.percentile(catt_draws, [2.5, 97.5], axis=0)
    # By Proposition 2, overall ATT draws depend only on mode 1 (since sum_s psi_k(s) == 0 for k >= 2)
    att_draws = np.einsum("ns,n->s", mode_tau_draws[0], att_unit_weights)
    att_hat = float(att_draws.mean())
    att_se = float(att_draws.std(ddof=1))
    q_low, q_high = np.percentile(att_draws, [2.5, 97.5])

    att_draws_raw = mode_tau_draws_raw_cat = catt_draws_raw = None
    report = None
    if projecting:
        catt_draws_raw = _assemble_catt_draws(
            mode_tau_draws_raw, cohorts_arr=cohorts_arr, waves=waves, bases=bases, n_units=N, n_periods=T,
            catt_floor=None,
        )
        report = _projection_report(
            mode1_raw=mode1_raw,
            mode1_proj=mode1_proj,
            catt_raw=catt_draws_raw,
            catt_proj=catt_draws,
            att_unit_weights=att_unit_weights,
            z_arr=z_arr if treated_cells is None else treated_cells.astype(np.int32),
            tau_clip=cfg.tau_clip,
            catt_floor=cfg.catt_floor,
        )
        if cfg.catt_floor is not None:
            _warn_rank_k(f"catt_floor={cfg.catt_floor!r} is a decision-layer projection applied after sampling, outside "
                "the coverage guarantees of the estimator: it moved "
                f"{report['mode1_draws_altered_share']:.1%} of the treated unit-level mode-1 draws (mean |shift| "
                f"{report['mode1_shift_mean_abs']:.3g}) and {report['cell_draws_altered_share']:.1%} of the treated cell "
                f"draws, and shifted the ATT draws by up to {report['att_shift_max_abs']:.3g} (a cohort whose mean "
                "effect is below the floor is raised to it). The unprojected draws are in att_draws_raw / "
                "mode_tau_draws_raw / catt_draws_raw (keep_unprojected_draws=True); report both.")
        if cfg.keep_unprojected_draws:
            att_draws_raw = np.einsum("ns,n->s", mode1_raw, att_unit_weights)
            mode_tau_draws_raw_cat = list(mode_tau_draws_raw)
        else:
            catt_draws_raw = None
    runtime_sec = float(time.perf_counter() - t_start)
    diag = dict(diagnostics)
    diag.setdefault("variant", variant)
    diag.setdefault("basis_spec", cfg.resolved_basis_spec)
    extra = {"unit_weights": np.asarray(att_unit_weights, dtype=np.float64), "cohorts": np.asarray(cohorts_arr)}
    extra.update(result_fields or {})
    if arm_metadata_override is not None:
        arm_meta_out = dict(arm_metadata_override)
        t0_ov = int(arm_meta_out.get("t0", min(waves) if waves else 1))
        t_post_ov = int(arm_meta_out.get("t_post", max(T - t0_ov, 1)))
        arm_meta_out["trt_trajectory_draws"] = _extract_event_time_trt_trajectory(
            catt_draws, ever, cohorts_arr, waves, t0_ov, t_post_ov
        )
        extra["arm_metadata"] = arm_meta_out
    elif arm_meta_inputs is not None:
        extra["arm_metadata"] = _build_stage1_arm_metadata(
            catt_draws=catt_draws,
            catt_hat=catt_hat,
            ever=ever,
            cohorts_arr=cohorts_arr,
            waves=waves,
            bases=bases,
            theta_f_modes=theta_f_modes,
            arm_meta_inputs=arm_meta_inputs,
        )
    return RankKLongBetResult(
        catt_draws=catt_draws,
        catt_hat=catt_hat,
        catt_ci_low=catt_ci_low,
        catt_ci_high=catt_ci_high,
        att_draws=att_draws,
        att_hat=att_hat,
        att_se=att_se,
        att_ci_low=float(q_low),
        att_ci_high=float(q_high),
        mode_tau_draws=mode_tau_draws,
        theta_f_modes=theta_f_modes,
        sigma_modes=sigma_modes,
        bases=bases,
        runtime_sec=runtime_sec,
        theta_f_draws=theta_f_draws,
        gp_hyper_draws=gp_hyper_draws,
        att_var_within=float(att_se**2) if att_var_within is None else float(att_var_within),
        att_var_between=float(att_var_between),
        att_var_slope=float(att_var_slope),
        n_components=int(n_components),
        rhat_att=float(compute_rhat(att_draws)) if rhat_att is None else float(rhat_att),
        ess_att=float(compute_ess(att_draws, method="bulk")) if ess_att is None else float(ess_att),
        att_draws_raw=att_draws_raw,
        mode_tau_draws_raw=mode_tau_draws_raw_cat,
        catt_draws_raw=catt_draws_raw,
        projection_report=report,
        diagnostics=diag,
        variant=variant,
        basis_spec=cfg.resolved_basis_spec,
        active_units=active_units,
        rank_k=int(len(mode_tau_draws)),
        **extra,
    )


_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))


def _warn_rank_k(
    message: str,
    category: type[Warning] = RankKConfigWarning,
    *,
    stacklevel: int = 3,
) -> None:
    """Issue a warning attributed to the first caller outside the ``longbet`` package.

    Args:
      message: Warning text.
      category: Warning class (default :class:`RankKConfigWarning`).
      stacklevel: Fallback caller stack depth on Python < 3.12 where
        ``skip_file_prefixes`` is unavailable.
    """
    if sys.version_info >= (3, 12):
        warnings.warn(message, category, skip_file_prefixes=(_PACKAGE_DIR,))
    else:  # pragma: no cover
        warnings.warn(message, category, stacklevel=stacklevel)


def _as_float_array(value: Any, name: str) -> np.ndarray:
    """Convert ``value`` to a float64 array, raising ``ValueError`` naming ``name`` on failure."""
    try:
        return np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric array; got {type(value).__name__}") from exc


def _require_finite(arr: np.ndarray, name: str) -> None:
    """Raise ``ValueError`` if ``arr`` contains NaN or infinite values."""
    n_bad = int(arr.size - np.count_nonzero(np.isfinite(arr)))
    if n_bad:
        raise ValueError(
            f"{name} contains {n_bad} non-finite value(s) (NaN or +/-inf); impute or drop them before fitting"
        )


def _require_unit_rows(arr: np.ndarray, name: str, n_units: int) -> np.ndarray:
    """Return ``arr`` as an ``(N, Q)`` matrix, raising ``ValueError`` unless it has one row per unit."""
    mat = arr[:, None] if arr.ndim == 1 else arr
    if mat.ndim != 2 or mat.shape[0] != n_units:
        raise ValueError(f"{name} must have one row per unit (N={n_units}); got shape {arr.shape}")
    return mat


# Options read by only one engine: a non-default value under the other engine is ignored.
_JAX_ONLY_OPTIONS = (
    "n_forecast_draws",
    "max_forecast_draws",
    "between_var_threshold",
    "num_chains",
    "draws_per_component",
    "theta_prior_mean",
    "theta_prior_sd",
    "use_factor_loadings",
    "n_factor_loadings",
    "use_spatial_gp",
    "lengthscale_grid_km",
    "baseline_cells",
    "thin",
    "sigma_nu_total",
    "sigma_nu_dynamic_total",
    "alpha_nu",
    "beta_nu",
    "nu_working_zeta_weight",
    "learn_sigma_nu",
    "sigma_nu_prior_scale",
    "treated_extra_variance",
)
_STOCHTREE_ONLY_OPTIONS = ("num_gfr", "ridge_eps", "num_threads")


def _max_forecast_draws(cfg: RankKLongBetConfig) -> int:
    """Return the effective expansion limit: ``cfg.max_forecast_draws``, ``None`` meaning ``n_forecast_draws``."""
    return int(cfg.n_forecast_draws if cfg.max_forecast_draws is None else cfg.max_forecast_draws)


def _non_default_options(cfg: RankKLongBetConfig, names: tuple[str, ...]) -> list[str]:
    """Describe the options in ``names`` that ``cfg`` sets away from their defaults.

    A frozen dataclass does not record which fields were passed explicitly, so a
    value that differs from the default stands in for an explicit setting. The
    variant-dependent fields are compared with the default of ``cfg``'s variant,
    the values that :class:`RankKLongBetConfig` resolves them to, and the
    effective ``max_forecast_draws`` (``None`` meaning ``n_forecast_draws``)
    with ``n_forecast_draws``.

    Args:
      cfg: Configuration to inspect.
      names: Field names of :class:`RankKLongBetConfig`.

    Returns:
      ``"cfg.<name>=<value>"`` for every changed option, in the order of ``names``.
    """
    defaults = {f.name: f.default for f in dc_fields(RankKLongBetConfig)}
    defaults.update(_PROPAGATED_DEFAULTS if cfg.uncertainty_propagation else _PLUG_IN_DEFAULTS)
    defaults["max_forecast_draws"] = cfg.n_forecast_draws
    defaults["treated_extra_variance"] = cfg.baseline_cells == "untreated"
    changed = []
    for name in names:
        value = _max_forecast_draws(cfg) if name == "max_forecast_draws" else getattr(cfg, name)
        default = defaults[name]
        if isinstance(default, tuple):
            try:
                differs = tuple(float(v) for v in value) != tuple(float(v) for v in default)
            except (TypeError, ValueError):
                differs = True
        else:
            differs = value != default
        if differs:
            changed.append(f"cfg.{name}={value!r}")
    return changed


def _ignored(items: list[str]) -> str:
    """Return ``"a, b are ignored"`` (``"a is ignored"`` for one item) for warning texts."""
    return f"{', '.join(items)} {'is' if len(items) == 1 else 'are'} ignored"


def _validate_rank_k_inputs(
    y: Any,
    z: Any,
    x: Any,
    *,
    cfg: RankKLongBetConfig,
    features: Any | None = None,
    unit_features: Any | None = None,
    ps: Any | None = None,
    baseline_forecast: Any | None = None,
    baseline_forecast_draws: Any | None = None,
    pre_anchor_loadings: Any | None = None,
    coords: Any | None = None,
    cluster_ids: Any | None = None,
    cohorts: Any | None = None,
    t_0_panel: Any | None = None,
    mask: Any | None = None,
    horizon: Any | None = None,
    warn: bool = True,
    warn_short_cohorts: bool = True,
) -> dict[str, Any]:
    """Validate the inputs of :func:`fit_rank_k_longbet`; warn about ignored or conflicting ones.

    Cheap (no model fitting), so it runs before any heavy computation and can be
    unit-tested directly. The arguments are those passed to
    :func:`fit_rank_k_longbet`, before ``unit_features``, ``ps``,
    ``baseline_forecast`` and ``pre_anchor_loadings`` are filled in from
    ``features``; a value taken from ``features`` is validated under its
    ``features.<name>`` name. All errors are checked first. The warnings follow,
    one per issue, and only for explicit inputs that the fit would silently
    ignore or that contradict each other (a config field set away from its
    default, or from its variant's default, counts as explicit), for the
    deprecated ``warm_burnin`` and for cohorts too short for ``rank_k`` modes.
    Choosing a supported variant or option, such as
    ``uncertainty_propagation=True``, never warns.

    Masked panels: ``y`` may hold ``NaN`` where a cell is unobserved (never
    ``inf``), or ``mask`` (``(N, T)``, ``True`` = observed) may be given; the two
    must agree exactly. Unobserved cells are allowed only as trailing
    (administrative) censoring of every unit's post-adoption window and as missing
    pre-period cells (at least one pre-period cell per unit must be observed); see
    :func:`longbet._rank_k_event_time.check_observation_structure`. Masked panels
    need a single adoption period (align the panel in event time) and the JAX
    engine (``NotImplementedError`` for the plug-in ``stochtree`` engine).

    Args:
      y: Outcome panel.
      z: Treatment panel.
      x: Unit covariates.
      cfg: Configuration of the fit.
      features: Optional :class:`~longbet.CausalFeatures` (or compatible object).
      unit_features: Explicit extra unit covariates.
      ps: Explicit selection score.
      baseline_forecast: Explicit counterfactual forecast.
      baseline_forecast_draws: Explicit posterior forecast draws.
      pre_anchor_loadings: Explicit pre-anchor factor loadings.
      coords: Unit (latitude, longitude) in degrees.
      cluster_ids: Unit cluster labels.
      cohorts: First-adoption periods.
      t_0_panel: Pre-period length override.
      mask: Explicit observation mask.
      horizon: Event-time horizon override.
      warn: Emit the warnings (errors are always raised).
      warn_short_cohorts: Warn about cohorts with fewer post-adoption periods than
        ``cfg.rank_k``; disabled for the placebo panel, whose capped
        pseudo-cohorts do not affect the audited ATT.

    Returns:
      ``{"obs": (N, T) observed-cell indicator, "masked": bool, "common_basis": bool,
      "horizon": int | None}``: the observation structure the fit uses.

    Raises:
      ValueError: If an input is malformed or would leak treated outcomes into
        the baseline (see :func:`fit_rank_k_longbet`).
      NotImplementedError: Masked panel with the plug-in (``stochtree``) engine.
    """
    # ---- Panel ------------------------------------------------------------------
    y_arr = _as_float_array(y, "y")
    if y_arr.ndim != 2 or y_arr.shape[1] < 2:
        raise ValueError(f"y must be a 2-D (N, T) outcome panel with T >= 2; got shape {y_arr.shape}")
    N, T = y_arr.shape
    n_inf = int(np.count_nonzero(np.isinf(y_arr)))
    if n_inf:
        raise ValueError(f"y contains {n_inf} infinite value(s); only NaN marks an unobserved cell")
    obs = ~np.isnan(y_arr)
    if mask is not None:
        m_arr = np.asarray(mask)
        if m_arr.shape != (N, T):
            raise ValueError(f"mask must be a boolean (N, T) = {(N, T)} array (True = observed); got shape {m_arr.shape}")
        if m_arr.dtype != bool and not np.all((m_arr == 0) | (m_arr == 1)):
            raise ValueError("mask must hold booleans (True = observed)")
        m_arr = m_arr.astype(bool)
        if not np.array_equal(m_arr, obs):
            n_bad = int(np.count_nonzero(m_arr != obs))
            raise ValueError(
                f"mask and NaN pattern of y disagree on {n_bad} cell(s): a cell must be NaN in y if and only if "
                "mask is False (set y to NaN where you want a cell unobserved)"
            )
    masked = bool(not np.all(obs))

    z_arr = _as_float_array(z, "z")
    if z_arr.shape != (N, T):
        raise ValueError(f"z shape {z_arr.shape} must match y shape {(N, T)}")
    if not np.all((z_arr == 0.0) | (z_arr == 1.0)):
        raise ValueError("z must contain only 0/1 treatment indicators")
    switch_back = np.any(np.diff(z_arr, axis=1) < 0.0, axis=1)
    if np.any(switch_back):
        raise ValueError(
            "z must be absorbing (non-decreasing in t: once treated, always treated); "
            f"{int(switch_back.sum())} unit(s) return to z = 0, e.g. unit {int(np.argmax(switch_back))}"
        )

    x_arr = _require_unit_rows(_as_float_array(x, "x"), "x", N)
    _require_finite(x_arr, "x")

    # ---- Adoption timing ------------------------------------------------------------
    first_adopt = np.where(np.any(z_arr == 1.0, axis=1), np.argmax(z_arr, axis=1), T)
    if cohorts is not None:
        c_arr = _as_float_array(cohorts, "cohorts")
        if c_arr.shape != (N,):
            raise ValueError(
                f"cohorts must be a 1-D array with one first-adoption period per unit (N={N}); "
                f"got shape {c_arr.shape}"
            )
        if (
            not np.all(np.isfinite(c_arr))
            or np.any(c_arr != np.round(c_arr))
            or np.any(c_arr < 0)
            or np.any(c_arr > T)
        ):
            raise ValueError(
                f"cohorts must be integer first-adoption periods in [0, T] = [0, {T}], "
                f"with cohorts[i] = T = {T} for never-treated units"
            )
        mismatch = c_arr.astype(np.int64) != first_adopt
        if np.any(mismatch):
            i = int(np.argmax(mismatch))
            seen = "z is never 1" if first_adopt[i] == T else f"z first equals 1 at t = {int(first_adopt[i])}"
            raise ValueError(
                f"cohorts is inconsistent with z for {int(mismatch.sum())} unit(s): z[i, t] must equal 1 "
                f"if and only if t >= cohorts[i] (cohorts[i] = T = {T} for never-treated units), but "
                f"cohorts[{i}] = {int(c_arr[i])} while {seen}"
            )
    ever = first_adopt < T
    if not np.any(ever) or np.all(ever):
        raise ValueError(
            "Panel must contain both ever-treated and never-treated units: the never-treated units "
            "identify every baseline component."
        )
    first_wave = int(first_adopt[ever].min())
    if first_wave < 1:
        raise ValueError(
            "Every ever-treated unit needs at least one pre-adoption period, but some units are "
            "treated from period 0; drop them or prepend their pre-adoption periods to y"
        )
    t0 = first_wave
    if t_0_panel is not None:
        t0_arr = np.asarray(t_0_panel)
        t0_val = float(t0_arr.reshape(-1)[0]) if t0_arr.size == 1 and t0_arr.dtype.kind in "iuf" else np.nan
        if not np.isfinite(t0_val) or t0_val != round(t0_val) or not 1 <= t0_val <= first_wave:
            raise ValueError(
                f"t_0_panel must be a single integer in [1, {first_wave}] (at most the first adoption "
                f"period); got {t_0_panel!r}. Periods from the first adoption on contain treated "
                "outcomes, which would leak into the pre-period means."
            )
        t0 = int(t0_val)

    # ---- Observation structure (masked panels and the common event-time basis) ------------
    waves_all = sorted({int(g) for g in first_adopt[ever]})
    single_wave = len(waves_all) == 1
    if cfg.common_basis is None:
        common_basis = bool(masked or single_wave)
    else:
        common_basis = bool(cfg.common_basis)
        if masked and not common_basis:
            raise ValueError(
                "common_basis=False is not available for masked panels: unobserved cells are handled on the "
                "common event-time basis (leave common_basis=None or set it to True)"
            )
    if masked:
        if cfg.engine == "stochtree":
            raise NotImplementedError(
                "Masked panels (NaN in y or a mask) are implemented for the JAX engine only "
                "(RankKLongBetConfig(uncertainty_propagation=True)); the plug-in stochtree engine does not "
                "support unobserved cells yet"
            )
        if not single_wave:
            raise ValueError(
                f"Masked panels need a single adoption period, but the treated units adopt in periods {waves_all}: "
                "align the panel in event time (every unit's first post-adoption period in the same column, "
                "unobserved cells set to NaN) before fitting"
            )
        if cfg.resolved_basis_spec.family == "pca_empirical":
            raise ValueError("ExposureBasisSpec(family='pca_empirical') is not available for masked panels")
    if horizon is not None:
        h_arr = np.asarray(horizon)
        h_val = float(h_arr.reshape(-1)[0]) if h_arr.size == 1 and h_arr.dtype.kind in "iuf" else np.nan
        if not np.isfinite(h_val) or h_val != round(h_val) or not 1 <= h_val <= T - first_wave:
            raise ValueError(
                f"horizon must be a single integer in [1, {T - first_wave}] (at most the post-adoption length of "
                f"the first wave); got {horizon!r}. To predict beyond the last panel period, append all-NaN columns to y."
            )
        if not common_basis:
            raise ValueError("horizon applies to the common event-time basis only; set common_basis=True")
    horizon_eff: int | None = None
    if common_basis:
        info = check_observation_structure(
            obs, first_adopt, t0, horizon=None if horizon is None else int(float(np.asarray(horizon).reshape(-1)[0])),
            min_support=cfg.min_support_modes,
        )
        horizon_eff = int(info["horizon"])
    structure = {"obs": obs, "masked": masked, "common_basis": common_basis, "horizon": horizon_eff}

    # ---- Features built on the same panel -------------------------------------------
    if features is not None:
        try:
            anchor = int(getattr(features, "anchor_t"))
        except (AttributeError, TypeError, ValueError):
            anchor = None
        if anchor is not None and anchor > first_wave:
            raise ValueError(
                f"features.anchor_t = {anchor} is after the first adoption period {first_wave}, so these "
                "features used post-adoption outcomes of treated units; build them from the same "
                "(y, z) with build_causal_features(y=y, z=z, ...)"
            )
        pf = getattr(features, "point_forecast", None)
        if pf is not None and np.shape(pf) not in ((N, T), (N, T, 1)):
            raise ValueError(
                f"features.point_forecast has shape {np.shape(pf)}, but y has shape {(N, T)}: build the "
                "features from the same (y, z) panel that is passed to the fit"
            )
        pf_draws = getattr(features, "point_forecast_draws", None)
        if pf_draws is not None and (np.ndim(pf_draws) != 3 or tuple(np.shape(pf_draws)[1:]) != (N, T)):
            raise ValueError(
                f"features.point_forecast_draws has shape {np.shape(pf_draws)}; expected (M, {N}, {T}): "
                "build the features from the same (y, z) panel that is passed to the fit"
            )

    # ---- Optional unit-level inputs (explicit arguments take precedence over features) ------
    def _resolve(explicit: Any | None, name: str, attr: str) -> tuple[Any | None, str]:
        """Return the value the fit will use and the name to report it under."""
        if explicit is not None or features is None:
            return explicit, name
        return getattr(features, attr, None), f"features.{attr}"

    uf_val, uf_name = _resolve(unit_features, "unit_features", "unit_features")
    ps_val, ps_name = _resolve(ps, "ps", "ps")
    bf_val, bf_name = _resolve(baseline_forecast, "baseline_forecast", "point_forecast")
    lam_val, lam_name = _resolve(pre_anchor_loadings, "pre_anchor_loadings", "pre_anchor_loadings")
    if uf_val is not None:
        uf = _require_unit_rows(_as_float_array(uf_val, uf_name), uf_name, N)
        _require_finite(uf, uf_name)
    if ps_val is not None:
        ps_arr = _as_float_array(ps_val, ps_name).reshape(-1)
        if ps_arr.shape[0] != N:
            raise ValueError(
                f"{ps_name} must hold one selection score per unit (N={N}); got shape {np.shape(ps_val)}"
            )
        _require_finite(ps_arr, ps_name)
        if np.any((ps_arr < 0.0) | (ps_arr > 1.0)):
            raise ValueError(f"{ps_name} must hold selection (propensity) scores in [0, 1]")
    if bf_val is not None:
        bf = _as_float_array(bf_val, bf_name)
        if bf.shape not in ((N, T), (N, T, 1)):
            raise ValueError(f"{bf_name} must have shape (N, T) = {(N, T)} or (N, T, 1); got {bf.shape}")
        _require_finite(bf, bf_name)
    n_draws_given = 0
    if baseline_forecast_draws is not None:
        bfd = _as_float_array(baseline_forecast_draws, "baseline_forecast_draws")
        if bfd.ndim != 3 or bfd.shape[0] < 1 or bfd.shape[1:] != (N, T):
            raise ValueError(
                f"baseline_forecast_draws must have shape (M, N, T) = (M, {N}, {T}); got {bfd.shape}"
            )
        _require_finite(bfd, "baseline_forecast_draws")
        n_draws_given = int(bfd.shape[0])
    n_loadings = None
    if lam_val is not None:
        lam = _require_unit_rows(_as_float_array(lam_val, lam_name), lam_name, N)
        _require_finite(lam, lam_name)
        n_loadings = int(lam.shape[1])

    # ---- Spatial inputs -------------------------------------------------------------------
    if coords is not None:
        c_mat = _as_float_array(coords, "coords")
        if c_mat.shape != (N, 2):
            raise ValueError(
                f"coords must have shape (N, 2) = ({N}, 2) with columns (latitude, longitude) in "
                f"degrees; got {c_mat.shape}"
            )
        _require_finite(c_mat, "coords")
        lat, lon = c_mat[:, 0], c_mat[:, 1]
        if np.any(np.abs(lat) > 90.0):
            raise ValueError(
                "coords[:, 0] must be latitudes in [-90, 90] degrees (columns are (latitude, "
                f"longitude)); got range [{lat.min():g}, {lat.max():g}]"
            )
        if np.any(lon < -180.0) or np.any(lon > 360.0):
            raise ValueError(
                f"coords[:, 1] must be longitudes in [-180, 360] degrees; got range [{lon.min():g}, {lon.max():g}]"
            )
    if cluster_ids is not None:
        cid = np.asarray(cluster_ids)
        if cid.shape != (N,):
            raise ValueError(
                f"cluster_ids must be a 1-D array with one cluster label per unit (N={N}); got shape {cid.shape}"
            )
        if cid.dtype.kind in "fc" and not np.all(np.isfinite(cid)):
            raise ValueError("cluster_ids contains NaN or infinite labels; give every unit a finite cluster label")

    if not warn:
        return structure

    is_jax = cfg.engine == "jax"
    has_forecast = bool(cfg.use_forecast and bf_val is not None)
    propagates = bool(is_jax and has_forecast and cfg.uncertainty_propagation and cfg.n_forecast_draws >= 2)
    no_forecast = (
        "use_forecast=False"
        if not cfg.use_forecast
        else "neither features nor baseline_forecast supplies a forecast"
    )

    # ---- Deprecated option -------------------------------------------------------------------
    if cfg.warm_burnin is not None:
        if propagates and cfg.warm_burnin > cfg.num_burnin:
            effect = f"raises the mode-1 burn-in from num_burnin={cfg.num_burnin} to {cfg.warm_burnin} sweeps"
        else:
            effect = "has no effect"
        _warn_rank_k(
            f"cfg.warm_burnin={cfg.warm_burnin!r} is deprecated and will be removed; in this fit it {effect}. "
            "It is not a warm start: when forecast draws are propagated, every chain-by-forecast-draw track "
            "of mode 1 runs its own burn-in of max(num_burnin, warm_burnin) sweeps, and warm_burnin is "
            "ignored everywhere else. Set num_burnin instead."
        )

    # ---- Budget of the other variant (as left by dataclasses.replace) ------------------------
    variant, other, other_budget = (
        ("propagated", "plug-in", _PLUG_IN_DEFAULTS)
        if cfg.uncertainty_propagation
        else ("plug-in", "propagated", _PROPAGATED_DEFAULTS)
    )
    if all(getattr(cfg, name) == value for name, value in other_budget.items()):
        budget = ", ".join(f"{name}={value}" for name, value in other_budget.items())
        _warn_rank_k(
            f"cfg selects the {variant} variant (uncertainty_propagation={cfg.uncertainty_propagation}) but "
            f"holds the default budget of the {other} variant ({budget}), as dataclasses.replace(cfg, "
            f"engine=None, uncertainty_propagation={cfg.uncertainty_propagation}) leaves it: replace keeps the "
            "stored budget. Build a new "
            f"RankKLongBetConfig(uncertainty_propagation={cfg.uncertainty_propagation}, ...) instead for the "
            f"budget validated for the {variant} variant."
        )

    # ---- Inputs and options the chosen variant does not read --------------------------------
    if is_jax:
        ignored = _non_default_options(cfg, _STOCHTREE_ONLY_OPTIONS)
        if ignored:
            _warn_rank_k(
                f"With the propagated variant (uncertainty_propagation=True, engine 'jax') {_ignored(ignored)}: "
                "these options apply only to the plug-in variant (engine 'stochtree')."
            )
    else:
        explicit = (
            ("coords", coords),
            ("cluster_ids", cluster_ids),
            ("pre_anchor_loadings", pre_anchor_loadings),
            ("baseline_forecast_draws", baseline_forecast_draws),
        )
        ignored = [name for name, value in explicit if value is not None]
        ignored += _non_default_options(cfg, _JAX_ONLY_OPTIONS)
        if ignored:
            _warn_rank_k(
                f"With the plug-in variant (uncertainty_propagation=False, engine 'stochtree') "
                f"{_ignored(ignored)}: it fits one stochtree BCF per mode on a plug-in OLS-ridge forecast "
                "slope and has no forecast draws, factor loadings, spatial random field or JAX chain "
                "settings. Drop these inputs, or set uncertainty_propagation=True for the propagated "
                "variant."
            )

    # ---- Forecast and forecast-slope options ---------------------------------------------------
    if not has_forecast:
        ignored = ["baseline_forecast"] if baseline_forecast is not None else []
        slope_options = (
            ("theta_f_clip", "theta_prior_mean", "theta_prior_sd") if is_jax else ("theta_f_clip", "ridge_eps")
        )
        ignored += _non_default_options(cfg, slope_options)
        if ignored:
            _warn_rank_k(f"No forecast is used ({no_forecast}), so {_ignored(ignored)}.")

    if is_jax:
        # ---- Forecast draws ----------------------------------------------------------------------
        if has_forecast and cfg.uncertainty_propagation and cfg.n_forecast_draws < 2:
            _warn_rank_k(
                f"uncertainty_propagation=True contradicts n_forecast_draws={cfg.n_forecast_draws}: with fewer "
                "than 2 forecast draws no forecast uncertainty is propagated, so the JAX sampler conditions on the "
                "point forecast. That is not a validated variant (without forecast draws its intervals "
                "under-covered under spatially correlated shocks); it is kept for software tests. Set "
                "n_forecast_draws >= 2 (validated: 8), or use the plug-in variant (uncertainty_propagation=False)."
            )
        m_max = _max_forecast_draws(cfg)
        if propagates:
            if baseline_forecast_draws is not None and n_draws_given < m_max:
                _warn_rank_k(
                    f"baseline_forecast_draws holds M={n_draws_given} draws, fewer than "
                    f"max_forecast_draws={m_max}, so it is ignored and the posterior forecast "
                    "draws are generated internally. Pass at least max_forecast_draws draws, or omit "
                    "baseline_forecast_draws."
                )
            if _non_default_options(cfg, ("num_mcmc",)):
                n_total = cfg.num_chains * cfg.n_forecast_draws * cfg.draws_per_component
                _warn_rank_k(
                    f"cfg.num_mcmc={cfg.num_mcmc!r} is ignored because forecast draws are propagated: the "
                    "posterior has num_chains x n_forecast_draws x draws_per_component = "
                    f"{cfg.num_chains} x {cfg.n_forecast_draws} x {cfg.draws_per_component} = {n_total} "
                    "draws. Change draws_per_component or num_chains instead."
                )
            if m_max == cfg.n_forecast_draws and _non_default_options(cfg, ("between_var_threshold",)):
                _warn_rank_k(
                    f"cfg.between_var_threshold={cfg.between_var_threshold!r} is ignored: it only triggers "
                    "the adaptive forecast-draw expansion, which is off because max_forecast_draws == "
                    f"n_forecast_draws = {cfg.n_forecast_draws}."
                )
        else:
            # engine 'jax' implies uncertainty_propagation=True (enforced by the config).
            reason = no_forecast if not has_forecast else f"n_forecast_draws={cfg.n_forecast_draws}"
            ignored = ["baseline_forecast_draws"] if baseline_forecast_draws is not None else []
            # On a masked panel with mode_coupling="propagated" the higher-mode draws run
            # through n_forecast_draws Mode-1 tracks of draws_per_component draws, so those
            # two options are used even without forecast draws.
            modes_propagated = bool(masked and cfg.mode_coupling == "propagated" and cfg.rank_k >= 2 and cfg.n_forecast_draws >= 2)
            draw_counts = {"n_forecast_draws": cfg.n_forecast_draws, "max_forecast_draws": m_max}
            if modes_propagated:
                draw_counts.pop("n_forecast_draws")
            ignored += _non_default_options(cfg, tuple(name for name, m in draw_counts.items() if m > 1))
            ignored += _non_default_options(
                cfg, ("between_var_threshold",) if modes_propagated else ("between_var_threshold", "draws_per_component")
            )
            if ignored:
                _warn_rank_k(
                    f"No forecast draws are propagated ({reason}), so {_ignored(ignored)}; num_mcmc sets "
                    "the number of posterior draws."
                )

        # ---- Factor loadings ---------------------------------------------------------------------
        if not cfg.use_factor_loadings:
            ignored = ["pre_anchor_loadings"] if pre_anchor_loadings is not None else []
            ignored += _non_default_options(cfg, ("n_factor_loadings",))
            if ignored:
                _warn_rank_k(f"use_factor_loadings=False, so {_ignored(ignored)}.")
        elif _non_default_options(cfg, ("n_factor_loadings",)):
            # Supplied loadings are used as they are; a matching column count is consistent.
            if n_loadings is not None and n_loadings != cfg.n_factor_loadings:
                hint = (
                    " Set n_factor_loadings in build_causal_features instead."
                    if lam_name.startswith("features.")
                    else ""
                )
                _warn_rank_k(
                    f"cfg.n_factor_loadings={cfg.n_factor_loadings!r} is ignored because the factor loadings "
                    f"come from {lam_name} ({n_loadings} column(s)); it applies only when the loadings are "
                    f"extracted from the in-panel pre-period.{hint}"
                )
            elif n_loadings is None and t0 < 2:
                _warn_rank_k(
                    f"cfg.n_factor_loadings={cfg.n_factor_loadings!r} is ignored: a pre-period of {t0} "
                    "period is too short to extract factor loadings, so the factor term is dropped."
                )

        # ---- Spatial random field ----------------------------------------------------------------
        has_coords = coords is not None
        has_clusters = cluster_ids is not None
        if not cfg.use_spatial_gp:
            spatial_inputs = (("coords", coords), ("cluster_ids", cluster_ids))
            ignored = [name for name, value in spatial_inputs if value is not None]
            ignored += _non_default_options(cfg, ("lengthscale_grid_km",))
            if ignored:
                _warn_rank_k(
                    f"use_spatial_gp=False, so {_ignored(ignored)} and no spatial random field is fitted."
                )
        else:
            if has_coords and has_clusters:
                _warn_rank_k(
                    "Both coords and cluster_ids were supplied: coords take precedence (Matérn-3/2 field on "
                    "great-circle distances) and cluster_ids is ignored. Pass only one of them."
                )
            if not has_coords and _non_default_options(cfg, ("lengthscale_grid_km",)):
                spatial = (
                    "the cluster random effect on cluster_ids has no range"
                    if has_clusters
                    else "without coords or cluster_ids no spatial random field is fitted"
                )
                _warn_rank_k(
                    f"cfg.lengthscale_grid_km={tuple(cfg.lengthscale_grid_km)!r} is ignored: it is the range "
                    f"grid of the Matérn-3/2 field on coords, and {spatial}."
                )

    # ---- Cohorts too short for rank_k modes ----------------------------------------------------
    if warn_short_cohorts:
        waves = sorted({int(g) for g in first_adopt[ever]})
        _resolve_short_cohorts(cfg, waves, T, warn=True)

    # ---- Options of the common event-time basis on the per-cohort path --------------------------
    if not common_basis:
        ignored = _non_default_options(cfg, _EVENT_TIME_ONLY_OPTIONS)
        if ignored:
            _warn_rank_k(
                f"{_ignored(ignored)}: with per-cohort exposure bases (staggered adoption without a mask and "
                "common_basis=None or False) every unit is fully observed on its own window, so the two-pass "
                "projection, its support thresholds and the mode coupling do not apply. Set common_basis=True "
                "for the common event-time basis."
            )
    return structure


def fit_rank_k_longbet(
    y: np.ndarray,
    z: np.ndarray,
    x: np.ndarray,
    *,
    features: Any | None = None,
    unit_features: np.ndarray | None = None,
    ps: np.ndarray | None = None,
    baseline_forecast: np.ndarray | None = None,
    baseline_forecast_draws: np.ndarray | None = None,
    pre_anchor_loadings: np.ndarray | None = None,
    coords: np.ndarray | None = None,
    cluster_ids: np.ndarray | None = None,
    cohorts: np.ndarray | None = None,
    t_0_panel: int | None = None,
    mask: np.ndarray | None = None,
    horizon: int | None = None,
    cfg: RankKLongBetConfig | None = None,
) -> RankKLongBetResult:
    """Fit Forecast-Augmented Rank-K Exposure-Basis BCF on a panel ``(y, z, x)``.

    Every treated unit's post-adoption outcomes are projected onto an
    orthonormal exposure basis whose first mode is the constant, and one
    cross-sectional BCF is fitted per mode, with the counterfactual forecast as
    a control variate whose slope is learned from never-treated units. Mode 1 is
    each treated unit's time-averaged effect, so the ATT depends on mode 1 only.

    One estimator, two uncertainty variants, selected by
    ``cfg.uncertainty_propagation`` (see :class:`RankKLongBetConfig`):
    "propagated" (the default; engine ``"jax"``, ``uncertainty_propagation=True``)
    learns every baseline component (prognostic forest, slope, factor shocks,
    spatial field, noise variance) from never-treated units, krigs the spatial
    field to the treated units and pools posterior forecast draws through mode 1.
    "Plug-in" (engine ``"stochtree"``, ``uncertainty_propagation=False``) fits
    one ``stochtree`` BCF per mode across all units on a plug-in OLS-ridge slope
    and is available as a fast plug-in approximation. The JAX sampler is an
    approximate partially collapsed Gibbs sampler (working likelihoods in the
    forest steps), not an exact posterior sampler.

    Recommended calls (the configurations validated in the companion paper)::

        feats = build_causal_features(y=y, z=z, x=x, y_pre_history=y_hist, seed=1)
        res = fit_rank_k_longbet(y, z, x, features=feats, coords=latlon)  # propagated default
        # or with cluster labels: fit_rank_k_longbet(y, z, x, features=feats, cluster_ids=region)
        # fast plug-in variant:
        cfg_plug = RankKLongBetConfig(uncertainty_propagation=False)
        res_plug = fit_rank_k_longbet(y, z, x, features=feats, cfg=cfg_plug)

    Leak-free requirements: ``features`` must be built from the same ``(y, z)``
    (its forecasts, trajectory summaries, selection score and factor loadings
    then use only outcomes before the first adoption period and never-treated
    outcomes); arrays passed directly (``baseline_forecast``,
    ``baseline_forecast_draws``, ``unit_features``, ``ps``,
    ``pre_anchor_loadings``) must likewise not depend on treated outcomes
    ``y[z == 1]``; ``x`` must hold pre-treatment covariates; and the pre-period
    ``y[:, :t_0_panel]`` may not reach the first adoption period.

    Interval coverage: in both variants the 95% ATT interval is the 2.5% /
    97.5% posterior percentiles of the pooled ATT draws. Both variants met
    every coverage rule of the paper's second frozen confirmation (pooled
    coverage 93.0% plug-in and 94.3% propagated, with propagated intervals 23%
    longer); the propagated variant also reached pooled coverage of 93.8% on
    the temperature panel and 90.6% on county employment in the development
    benchmark. Unit-level (CATE) intervals, such as those computed from
    ``mode_tau_draws``, are approximate and below nominal in the paper's
    benchmarks (95% unit-level intervals covered 88-93% with the plug-in and
    70-92% with the propagated variant), and the cell-level intervals
    ``catt_ci_low`` / ``catt_ci_high`` (2.5% / 97.5% posterior percentiles) are
    approximate as well; treat both as descriptive.

    Masked event-time panels (0.7.0). ``y`` may hold ``NaN`` where a unit is not
    observed (or an explicit ``mask``): each unit's observed post-adoption cells
    must be a contiguous prefix ``s = 1, ..., L_i`` of the event-time horizon
    ``H`` (administrative censoring) and at least one pre-period cell must be
    observed. All treated units must then share one adoption period (align the
    panel in event time). Every unit is projected on the *common* orthonormal
    basis on ``s = 1..H`` (``cfg.common_basis``, the default for masked and
    single-wave panels) with the two-pass projection of
    :mod:`longbet._rank_k_event_time`: higher modes on the units whose window
    supports them (``cfg.min_support_modes``), mode 1 on every unit as the
    observed-window mean of the residual after the fitted treatment dynamics are
    removed (``cfg.mode_coupling``: higher-mode posterior draws propagated through
    Mode-1 tracks, or the posterior mean), with precision weights from the
    observed window. ``catt_draws`` then holds the posterior-predictive cell
    effect ``tau_1(x_i) + sum_k tau_k(x_i) psi_k(s)`` on *every* post cell
    ``s <= H`` of every treated unit (``catt_observed`` marks the observed ones),
    ``att_draws`` the ``cfg.att_weighting`` ("cell": observed-cell weights, the
    paper's convention; "unit": one weight per treated unit, the population
    estimand) average of the mode-1 levels, and :func:`longbet.cumulative_effect`
    / :func:`longbet.event_time_trajectory` / :func:`longbet.rank_k_transport_audit`
    the population aggregates and the transport diagnostic. The observed window
    length ``L_i`` is appended to the forest covariates (``cfg.support_covariate``).
    Under a mask the ATT depends on the fitted dynamics (``rank_k``, basis family)
    through the correction, so report a ``rank_k`` sensitivity; on fully observed
    panels the correction is exactly zero and the fit reproduces the per-cohort
    estimator bit for bit. The plug-in ``stochtree`` engine raises
    ``NotImplementedError`` for masked panels.

    Args:
      y: Outcome panel, shape ``(N, T)``; ``NaN`` marks an unobserved cell (see
        above), ``inf`` is never allowed.
      z: Binary absorbing treatment panel, shape ``(N, T)``: ``z[i, t] = 1`` from
        unit ``i``'s first adoption period on. The panel needs both ever-treated
        and never-treated units, and every treated unit at least one
        pre-adoption period.
      x: Pre-treatment unit covariates, shape ``(N, P)`` or ``(N,)``, finite;
        used by both forests.
      features: Output of :func:`build_causal_features` for the same ``(y, z)``
        (recommended). Supplies ``unit_features``, ``ps``, ``baseline_forecast``
        (its ``point_forecast``) and, for the propagated variant,
        ``pre_anchor_loadings`` and the posterior forecast draws (drawn on
        demand, so the builder's own ``n_forecast_draws`` can stay 0). Explicit
        arguments take precedence.
      unit_features: Extra unit covariates, shape ``(N, Q)``, appended to ``x``.
      ps: Selection (propensity) score in ``[0, 1]``, shape ``(N,)``, a covariate
        of the prognostic forest (clipped to ``[0.02, 0.98]``). Default: the
        ever-treated share.
      baseline_forecast: Counterfactual forecast, shape ``(N, T)`` or
        ``(N, T, 1)``; the control variate and a prognostic covariate of every
        mode.
      baseline_forecast_draws: Posterior forecast draws of the propagated
        variant, shape ``(M, N, T)``; the first ``max_forecast_draws`` are used
        when ``M >= max_forecast_draws``, otherwise the draws are ignored (with a
        warning) and generated internally: from ``features`` when given, else
        from a control-group Normal-inverse-gamma model around
        ``baseline_forecast``.
      pre_anchor_loadings: Pre-anchor factor loadings ``lambda_i`` of the
        propagated variant, shape ``(N, q)``. Default: from ``features``, else
        extracted from the pre-period residuals ``y[:, :t0]`` minus the forecast.
      coords: Shape ``(N, 2)``: (latitude, longitude) in degrees, for the
        Matérn-3/2 spatial field on great-circle distances (propagated variant).
        Without ``coords`` or ``cluster_ids`` no spatial field is fitted; with
        the propagated variant, pass one of them whenever units can share
        regional shocks.
      cluster_ids: Length-``N`` cluster labels for a cluster random effect
        (propagated variant); used only when ``coords`` is ``None``.
      cohorts: First-adoption period of every unit, length ``N`` (``T`` for
        never-treated units); must agree with ``z``. Default: inferred from ``z``.
      t_0_panel: Number ``t0`` of leading periods forming the pre-period, an
        integer in ``[1, first adoption period]``. Default: the first adoption
        period.
      mask: Optional boolean ``(N, T)`` observation indicator (``True`` =
        observed); must agree with the ``NaN`` pattern of ``y`` exactly. ``NaN``
        in ``y`` alone is sufficient.
      horizon: Event-time horizon ``H`` of the common basis, an integer in
        ``[1, T - first adoption period]``. Default: the longest observed
        post-adoption window of a treated unit. Observed cells beyond ``H``
        leave the likelihood; to predict beyond the panel, append all-``NaN``
        columns to ``y``.
      cfg: :class:`RankKLongBetConfig`; default ``RankKLongBetConfig()``, the
        propagated variant at its validated benchmark budget.

    Returns:
      A :class:`RankKLongBetResult`: the ATT summaries ``att_hat``, ``att_se``
      and the interval ``att_ci_low`` / ``att_ci_high`` (2.5% / 97.5% posterior
      percentiles of ``att_draws``); cell effects ``catt_draws`` / ``catt_hat`` /
      ``catt_ci_low`` / ``catt_ci_high`` of shape ``(S, N, T)`` / ``(N, T)``;
      ``mode_tau_draws`` (``K`` arrays of shape ``(N, S)``, row ``i`` of
      ``mode_tau_draws[0]`` being treated unit ``i``'s time-averaged effect);
      slope, spatial and variance-decomposition diagnostics; the convergence
      diagnostics ``rhat_att`` (NaN for the plug-in variant, which runs one
      chain) and ``ess_att``; and, for the common event-time basis, the
      observation fields ``catt_observed``, ``horizon``, ``unit_weights``,
      ``att_weighting``, ``support_mask``, ``post_window_length`` and ``cohorts``.

    Raises:
      ValueError: If ``y`` holds ``inf`` or ``x`` is non-finite or ``x`` does not
        have ``N`` rows; ``z`` is not 0/1 or not absorbing; the panel lacks
        treated or never-treated units or has units treated in period 0;
        ``cohorts`` disagrees with ``z``; ``t_0_panel`` is outside ``[1, first
        adoption period]``; ``features`` was built on a panel of another shape or
        with an anchor after the first adoption period; ``ps``,
        ``unit_features``, ``baseline_forecast``, ``baseline_forecast_draws`` or
        ``pre_anchor_loadings`` is misshapen or non-finite (or ``ps`` is outside
        ``[0, 1]``); ``coords`` is not ``(N, 2)``, non-finite or outside
        ``[-90, 90] x [-180, 360]``; ``cluster_ids`` does not have length ``N``;
        or, for masked panels, ``mask`` disagrees with the ``NaN`` pattern, a
        unit has a gap inside its post-adoption window or no observed pre-period
        cell, a unit observes fewer than ``min_support_modes[0]`` post-adoption
        periods, a period has observed treated but no observed never-treated
        outcomes, the treated units adopt in several periods, or ``horizon`` is
        out of range.
      NotImplementedError: Masked panel with the plug-in ``stochtree`` engine.

    Warns:
      RankKConfigWarning: Only for explicit inputs that the fit would silently
        ignore or that contradict each other (a config field set away from its
        variant's default counts as explicit), never for choosing a supported
        variant or option. Triggers: ``coords``, ``cluster_ids``,
        ``pre_anchor_loadings``, ``baseline_forecast_draws`` or non-default
        propagated-only options (``n_forecast_draws``, ``num_chains``,
        ``draws_per_component``, the slope-prior, factor-loading and spatial
        options, ``baseline_cells``, ...) with the plug-in variant; non-default
        ``num_gfr``, ``ridge_eps`` or ``num_threads`` with the propagated
        variant; both ``coords`` and ``cluster_ids`` (``coords`` take
        precedence); ``coords`` / ``cluster_ids`` with ``use_spatial_gp=False``;
        a non-default ``lengthscale_grid_km`` without ``coords``;
        ``pre_anchor_loadings`` or a non-default ``n_factor_loadings`` that is
        not used (supplied loadings with another number of columns, or a
        pre-period too short to extract them); ``baseline_forecast`` with
        ``use_forecast=False`` or non-default slope options without a forecast;
        ``uncertainty_propagation=True`` with ``n_forecast_draws`` 0 or 1; in
        the propagated variant without forecast draws,
        ``baseline_forecast_draws``, ``n_forecast_draws > 1``,
        ``max_forecast_draws > 1``, ``draws_per_component`` or
        ``between_var_threshold``; and, when forecast draws are propagated, a
        non-default ``num_mcmc``, a non-default ``between_var_threshold``
        without adaptive expansion, or fewer than ``max_forecast_draws``
        ``baseline_forecast_draws``. Also when ``num_chains``, ``num_burnin``
        and the three tree counts all equal the other variant's default budget
        (as ``dataclasses.replace`` leaves them when the variant is switched),
        for the deprecated ``cfg.warm_burnin`` and for cohorts with fewer
        post-adoption periods than ``rank_k``.
    """
    return _fit_rank_k_longbet(
        y,
        z,
        x,
        features=features,
        unit_features=unit_features,
        ps=ps,
        baseline_forecast=baseline_forecast,
        baseline_forecast_draws=baseline_forecast_draws,
        pre_anchor_loadings=pre_anchor_loadings,
        coords=coords,
        cluster_ids=cluster_ids,
        cohorts=cohorts,
        t_0_panel=t_0_panel,
        mask=mask,
        horizon=horizon,
        cfg=cfg,
    )


def _fit_rank_k_longbet(
    y: np.ndarray,
    z: np.ndarray,
    x: np.ndarray,
    *,
    features: Any | None = None,
    unit_features: np.ndarray | None = None,
    ps: np.ndarray | None = None,
    baseline_forecast: np.ndarray | None = None,
    baseline_forecast_draws: np.ndarray | None = None,
    pre_anchor_loadings: np.ndarray | None = None,
    coords: np.ndarray | None = None,
    cluster_ids: np.ndarray | None = None,
    cohorts: np.ndarray | None = None,
    t_0_panel: int | None = None,
    mask: np.ndarray | None = None,
    horizon: int | None = None,
    cfg: RankKLongBetConfig | None = None,
    warn_short_cohorts: bool = True,
) -> RankKLongBetResult:
    """Implementation of :func:`fit_rank_k_longbet` (``warn_short_cohorts`` as in the validator)."""
    t_start = time.perf_counter()
    if cfg is None:
        cfg = RankKLongBetConfig()
    if cfg.engine not in ("jax", "stochtree"):
        raise ValueError(f"cfg.engine must be 'jax' or 'stochtree'; got {cfg.engine!r}")
    if cfg.baseline_cells not in ("untreated", "all"):
        raise ValueError(f"cfg.baseline_cells must be 'untreated' or 'all'; got {cfg.baseline_cells!r}")

    # Validate the explicit arguments before they are filled in from ``features``, so
    # the warnings can tell explicit inputs from the ones that ``features`` supplies.
    structure = _validate_rank_k_inputs(
        y,
        z,
        x,
        cfg=cfg,
        features=features,
        unit_features=unit_features,
        ps=ps,
        baseline_forecast=baseline_forecast,
        baseline_forecast_draws=baseline_forecast_draws,
        pre_anchor_loadings=pre_anchor_loadings,
        coords=coords,
        cluster_ids=cluster_ids,
        cohorts=cohorts,
        t_0_panel=t_0_panel,
        mask=mask,
        horizon=horizon,
        warn_short_cohorts=warn_short_cohorts,
    )
    if t_0_panel is not None:
        t_0_panel = int(np.asarray(t_0_panel, dtype=np.float64).reshape(-1)[0])

    if features is not None:
        if unit_features is None:
            unit_features = getattr(features, "unit_features", None)
        if ps is None:
            ps = getattr(features, "ps", None)
        if baseline_forecast is None:
            baseline_forecast = getattr(features, "point_forecast", None)
        if pre_anchor_loadings is None:
            pre_anchor_loadings = getattr(features, "pre_anchor_loadings", None)

    y_arr = np.asarray(y, dtype=np.float64)
    z_arr = np.asarray(z, dtype=np.int32)
    x_base = np.asarray(x, dtype=np.float64)
    if x_base.ndim == 1:
        x_base = x_base[:, None]
    N, T = y_arr.shape
    if z_arr.shape != (N, T):
        raise ValueError(f"z shape {z_arr.shape} must match y shape {(N, T)}")

    cohorts_arr = _infer_cohorts(z_arr) if cohorts is None else np.asarray(cohorts, dtype=np.int32)
    ever = cohorts_arr < T
    ctrl = ~ever
    if not np.any(ever) or not np.any(ctrl):
        raise ValueError("Panel must contain both ever-treated and never-treated units.")
    ctrl_idx = np.nonzero(ctrl)[0].astype(np.int32)
    trt_idx = np.nonzero(ever)[0].astype(np.int32)

    t0 = int(np.min(cohorts_arr[ever])) if t_0_panel is None else int(t_0_panel)
    if t0 <= 0:
        raise ValueError(f"First adoption period t0 must be >= 1; got {t0}")

    obs = np.asarray(structure["obs"], dtype=bool)
    masked = bool(structure["masked"])
    # Complete panels keep the paper's arithmetic verbatim; masked panels average the
    # observed cells only (unit pre-period anchor and control period means).
    if masked:
        y_pre = _nanmean_rows(y_arr[:, :t0])
        dy = y_arr - y_pre
        dy_dm = dy - _nanmean_cols(dy[ctrl])
    else:
        y_pre = y_arr[:, :t0].mean(axis=1, keepdims=True)
        dy = y_arr - y_pre
        dy_dm = dy - dy[ctrl].mean(axis=0, keepdims=True)

    has_forecast = bool(cfg.use_forecast and baseline_forecast is not None)
    if has_forecast:
        f_arr = np.asarray(baseline_forecast, dtype=np.float64)
        if f_arr.ndim == 3:
            f_arr = f_arr[:, :, 0]
        if masked:
            f_pre = _nanmean_rows(np.where(obs[:, :t0], f_arr[:, :t0], np.nan))
            df = f_arr - f_pre
            df_dm = df - _nanmean_cols(np.where(obs[ctrl], df[ctrl], np.nan))
        else:
            f_pre = f_arr[:, :t0].mean(axis=1, keepdims=True)
            df = f_arr - f_pre
            df_dm = df - df[ctrl].mean(axis=0, keepdims=True)
    else:
        f_arr = np.zeros_like(y_arr)
        df_dm = np.zeros_like(dy_dm)

    if unit_features is not None:
        u_feat = np.asarray(unit_features, dtype=np.float64)
        if u_feat.ndim == 1:
            u_feat = u_feat[:, None]
        x_mod = np.column_stack([x_base, u_feat])
    else:
        x_mod = x_base

    if ps is not None:
        ps_vec = np.clip(np.asarray(ps, dtype=np.float64).reshape(-1), 0.02, 0.98)
    else:
        ps_vec = np.full(N, float(np.mean(ever)), dtype=np.float64)

    waves = sorted(int(g) for g in np.unique(cohorts_arr[ever]))
    K = _resolve_short_cohorts(cfg, waves, T, warn=False)
    spec = cfg.resolved_basis_spec
    layout: EventTimeLayout | None = None
    if structure["common_basis"]:
        # Common event-time basis on s = 1..H; per-unit observed windows (masked panels,
        # single-wave panels, or staggered panels with common_basis=True).
        def _basis_factory(h: int, k_req: int) -> np.ndarray:
            if spec.is_default:
                return orthonormal_exposure_basis(h, k_req)
            return orthonormal_exposure_basis(
                h, k_req, spec,
                empirical_trajectories=dy_dm[ctrl, t0 : t0 + h] if spec.family == "pca_empirical" else None,
            )

        layout = build_event_time_layout(
            obs, cohorts_arr, t0, horizon=structure["horizon"], rank_k=K,
            min_support=cfg.min_support_modes, basis_factory=_basis_factory,
        )
        bases = layout.bases
        wave_weights = layout.wave_weights
        total_trt_cells = float(sum(int(np.sum(layout.window_len[g][cohorts_arr == g])) for g in waves))
        if layout.has_partial_windows() and cfg.support_covariate:
            x_mod = np.column_stack([x_mod, layout.post_len.astype(np.float64)])
    else:
        total_trt_cells = float(np.sum(z_arr == 1))
        wave_weights = {
            g: float(np.sum(cohorts_arr == g) * (T - g)) / max(total_trt_cells, 1.0)
            for g in waves
        }
        if spec.is_default:
            bases = {g: orthonormal_exposure_basis(T - g, K) for g in waves}
        else:
            bases = {
                g: orthonormal_exposure_basis(
                    T - g, K, spec,
                    empirical_trajectories=dy_dm[ctrl, g:] if spec.family == "pca_empirical" else None,
                )
                for g in waves
            }
    cfg = dc_replace(cfg, rank_k=K) if K != int(cfg.rank_k) else cfg

    # Unit weights of the ATT: "cell" (the paper's L_{G_i} / |O_trt|, observed cells under a
    # mask) or "unit" (1 / N_trt, the population estimand).
    att_unit_weights = _att_unit_weights(cohorts_arr, waves, T, total_trt_cells, cfg.att_weighting, layout)
    event_time_fields: dict[str, Any] | None = None
    treated_cells_modeled: np.ndarray | None = None
    if layout is not None:
        event_time_fields = {
            "catt_observed": layout.treated_cells_observed(),
            "horizon": int(layout.horizon),
            "att_weighting": cfg.att_weighting,
            "support_mask": layout.support_matrix(),
            "post_window_length": layout.post_len.copy(),
        }
        treated_cells_modeled = layout.treated_cells_modeled()
    else:
        event_time_fields = {"att_weighting": cfg.att_weighting}

    if cfg.engine == "stochtree":
        if layout is not None and layout.has_partial_windows():
            raise NotImplementedError(
                "The common event-time basis with heterogeneous observed windows requires engine='jax' "
                "(RankKLongBetConfig(uncertainty_propagation=True)); the plug-in stochtree engine fits "
                "fully observed per-cohort windows only"
            )
        return _fit_rank_k_stochtree(
            y_arr=y_arr,
            z_arr=z_arr,
            x_mod=x_mod,
            ps_vec=ps_vec,
            dy_dm=dy_dm,
            df_dm=df_dm,
            has_forecast=has_forecast,
            cohorts_arr=cohorts_arr,
            ever=ever,
            ctrl=ctrl,
            waves=waves,
            wave_weights=wave_weights,
            bases=bases,
            cfg=cfg,
            t_start=t_start,
            att_w=att_unit_weights,
            result_fields=event_time_fields,
        )

    build_spatial_kernel_cache, sample_rank_k_mode_jax = _get_rank_k_jax()
    compute_rhat, compute_ess = _get_rhat_ess()

    # Pre-anchor factor loadings lambda_i (strictly pre-anchor: invariant to y[ever, t0:])
    if cfg.use_factor_loadings:
        if pre_anchor_loadings is not None:
            lam_mat = np.asarray(pre_anchor_loadings, dtype=np.float64)
            if lam_mat.ndim == 1:
                lam_mat = lam_mat[:, None]
        elif t0 >= 2:
            pre_resid = y_arr[:, :t0] - (f_arr[:, :t0] if has_forecast else 0.0)
            if masked:
                # Unobserved pre-period residuals enter the SVD at the unit's own observed
                # pre-period mean, i.e. with a neutral (zero) deviation pattern.
                pre_resid = np.where(obs[:, :t0], pre_resid, _nanmean_rows(pre_resid))
            lam_mat = np.asarray(
                extract_pre_anchor_loadings(pre_resid, q=int(cfg.n_factor_loadings)),
                dtype=np.float64,
            )
        else:
            lam_mat = np.zeros((N, 0), dtype=np.float64)
    else:
        lam_mat = np.zeros((N, 0), dtype=np.float64)

    # Determine initial component count M and chain draw count S_comp
    num_chains = max(int(cfg.num_chains), 1)
    # Deprecated: the sampler burns in max(num_burnin, warm_burnin) sweeps in the
    # mode-1 forecast-draw path, so ``None`` (0) leaves the burn-in at ``num_burnin``.
    warm_burnin_eff = 0 if cfg.warm_burnin is None else int(cfg.warm_burnin)
    if cfg.uncertainty_propagation and has_forecast and int(cfg.n_forecast_draws) >= 2:
        m_init = int(cfg.n_forecast_draws)
        m_max = max(_max_forecast_draws(cfg), m_init)
        draws_per_chain = max(int(cfg.draws_per_component), 2)
    else:
        m_init = 1
        m_max = 1
        draws_per_chain = max(int(round(float(cfg.num_mcmc) / float(num_chains))), 2)

    f_draws_all = _resolve_forecast_draws(
        n_draws=m_max,
        features=features,
        baseline_forecast_draws=baseline_forecast_draws,
        y_arr=y_arr,
        f_arr=f_arr,
        x_mod=x_mod,
        ctrl_idx=ctrl_idx,
        t0=t0,
        seed=int(cfg.random_seed),
        allow_fallback=not masked,
    )
    if f_draws_all is None:
        m_init = 1
        m_max = 1
        draws_per_chain = max(int(round(float(cfg.num_mcmc) / float(num_chains))), 2)
        df_dm_draws_all = None
    elif masked:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            f_pre_draws = np.nanmean(np.where(obs[None, :, :t0], f_draws_all[:, :, :t0], np.nan), axis=2, keepdims=True)
            df_draws = f_draws_all - f_pre_draws
            ctrl_means = np.nanmean(np.where(obs[None, ctrl, :], df_draws[:, ctrl, :], np.nan), axis=1, keepdims=True)
        df_dm_draws_all = df_draws - ctrl_means
    else:
        f_pre_draws = f_draws_all[:, :, :t0].mean(axis=2, keepdims=True)
        df_draws = f_draws_all - f_pre_draws
        df_dm_draws_all = df_draws - df_draws[:, ctrl, :].mean(axis=1, keepdims=True)

    def _decompose_mode1_att(comp_tau: np.ndarray) -> tuple[float, float, float]:
        """Exact algebraic within/between variance decomposition of ATT draws."""
        # comp_tau has shape (M, K_c, N) -> att_comp has shape (M, K_c)
        att_comp = np.einsum("msn,n->ms", comp_tau, att_unit_weights)
        M_c, K_c = att_comp.shape
        total_draws = M_c * K_c
        if total_draws <= 1:
            return 0.0, 0.0, 0.0
        grand_mean = float(np.mean(att_comp))
        comp_means = np.mean(att_comp, axis=1, keepdims=True)  # (M, 1)
        ss_within = float(np.sum((att_comp - comp_means) ** 2))
        ss_between = float(K_c * np.sum((comp_means[:, 0] - grand_mean) ** 2)) if M_c > 1 else 0.0
        denom = float(total_draws - 1)
        v_within = ss_within / denom
        v_between = ss_between / denom
        v_total = v_within + v_between
        return v_within, v_between, v_total

    mode_tau_draws: list[np.ndarray] = []
    theta_f_modes: list[float] = []
    theta_f_draws: list[np.ndarray] = []
    sigma_modes: list[float] = []
    gp_hyper_draws: list[dict[str, np.ndarray]] = []
    active_units: list[np.ndarray] = []
    diagnostics: dict[str, Any] = {
        "engine": "jax",
        "s_eta2_draws": [],
        "sigma_nu_total_draws": [],
        "tau_intercept_draws": [],
        "rhat_mode": [],
        "ess_bulk_mode": [],
        "ess_tail_mode": [],
        "y_scale_mode": [],
    }

    m_eff = m_init
    att_var_within = 0.0
    att_var_between = 0.0
    att_var_slope = 0.0
    rhat_att = 1.0
    ess_att = 0.0

    def _sampler_kwargs(k: int) -> dict[str, Any]:
        """Issue-3 controls of the JAX sampler for mode ``k`` (defaults reproduce the paper)."""
        sigma_nu_k = (
            float(cfg.sigma_nu_total)
            if (k == 0 or cfg.sigma_nu_dynamic_total is None)
            else float(cfg.sigma_nu_dynamic_total)
        )
        return dict(
            att_unit_weights=att_unit_weights if k == 0 else None,
            tau_prior_family=(cfg.tau_prior_family if k == 0 else "normal"),
            sigma_nu_total=sigma_nu_k,
            alpha_nu=float(cfg.alpha_nu),
            beta_nu=float(cfg.beta_nu),
            nu_working_zeta_weight=float(cfg.nu_working_zeta_weight),
            learn_sigma_nu=bool(cfg.learn_sigma_nu),
            sigma_nu_prior_scale=float(cfg.sigma_nu_prior_scale),
            treated_extra_variance=bool(cfg.treated_extra_variance),
            thin=int(cfg.thin),
        )

    def _run_mode(
        k: int,
        *,
        dy_k: np.ndarray,
        df_k: np.ndarray,
        df_k_draws_m: np.ndarray | None,
        w_unit: np.ndarray,
        n_trees_tau: int,
        draws_per_chain_k: int,
        idx: np.ndarray | None,
        dy_k_draws_m: np.ndarray | None = None,
    ) -> dict[str, Any]:
        """Sample mode ``k``; with ``idx`` the problem is restricted to the active units and scattered back."""
        if idx is None:
            sub = dict(
                dy_k=dy_k, df_k=df_k, df_k_draws=df_k_draws_m, x_mod=x_mod, ps_vec=ps_vec,
                z_vec=ever.astype(np.float64), w_unit=w_unit, lambda_mat=lam_mat, ctrl_idx=ctrl_idx,
                trt_idx=trt_idx, coords_s=coords, cluster_s=cluster_ids, dy_k_draws=dy_k_draws_m,
            )
        else:
            pos = np.full(N, -1, dtype=np.int64)
            pos[idx] = np.arange(idx.size)
            sub = dict(
                dy_k=dy_k[idx], df_k=df_k[idx],
                df_k_draws=None if df_k_draws_m is None else df_k_draws_m[:, idx],
                x_mod=x_mod[idx], ps_vec=ps_vec[idx], z_vec=ever[idx].astype(np.float64), w_unit=w_unit[idx],
                lambda_mat=lam_mat[idx], ctrl_idx=pos[ctrl_idx].astype(np.int32),
                trt_idx=pos[trt_idx[np.isin(trt_idx, idx)]].astype(np.int32),
                coords_s=None if coords is None else np.asarray(coords)[idx],
                cluster_s=None if cluster_ids is None else np.asarray(cluster_ids)[idx],
                dy_k_draws=None if dy_k_draws_m is None else dy_k_draws_m[:, idx],
            )
        spatial_cache = build_spatial_kernel_cache(
            coords=sub["coords_s"],
            cluster_ids=sub["cluster_s"],
            ctrl_idx=sub["ctrl_idx"],
            trt_idx=sub["trt_idx"],
            w_unit=sub["w_unit"],
            ell_grid_km=cfg.lengthscale_grid_km,
            use_spatial_gp=bool(cfg.use_spatial_gp),
        )
        extra_draws = {} if sub["dy_k_draws"] is None else {"dy_k_draws": sub["dy_k_draws"]}
        res = sample_rank_k_mode_jax(
            dy_k=sub["dy_k"],
            df_k_point=sub["df_k"],
            df_k_draws=sub["df_k_draws"],
            x_mod=sub["x_mod"],
            ps_vec=sub["ps_vec"],
            z_vec=sub["z_vec"],
            w_unit=sub["w_unit"],
            lambda_mat=sub["lambda_mat"],
            ctrl_idx=sub["ctrl_idx"],
            trt_idx=sub["trt_idx"],
            spatial_cache=spatial_cache,
            num_chains=num_chains,
            num_burnin=int(cfg.num_burnin),
            num_draws_per_chain=draws_per_chain_k,
            warm_burnin=warm_burnin_eff,
            num_trees_mu=int(cfg.num_trees_mu),
            num_trees_nu=n_trees_tau,
            has_forecast=has_forecast,
            use_forests=True,
            use_factor_loadings=bool(cfg.use_factor_loadings),
            use_spatial_gp=bool(cfg.use_spatial_gp),
            baseline_untreated_only=(cfg.baseline_cells == "untreated"),
            theta_f_clip=cfg.theta_f_clip,
            theta_prior_mean=float(cfg.theta_prior_mean),
            theta_prior_sd=float(cfg.theta_prior_sd),
            tau_prior_mean=(float(cfg.tau_prior_mean) if k == 0 else 0.0),
            tau_prior_sd=(cfg.tau_prior_sd if k == 0 else None),
            tau_clip=(cfg.tau_clip if k == 0 else (-float("inf"), float("inf"))),
            random_seed=int(cfg.random_seed + 101 * k),
            **_sampler_kwargs(k),
            **extra_draws,
        )
        res = dict(res)
        if "resid_mean" not in res:
            res["resid_mean"] = np.zeros_like(sub["dy_k"], dtype=np.float64)
        if "resid_var" not in res:
            res["resid_var"] = np.zeros_like(sub["dy_k"], dtype=np.float64)
        if idx is not None:
            res["tau_draws"] = _scatter_rows(res["tau_draws"], idx, N, axis=0)
            res["tau_chains"] = _scatter_rows(res["tau_chains"], idx, N, axis=2)
            res["comp_tau_draws"] = _scatter_rows(res["comp_tau_draws"], idx, N, axis=2)
            res["resid_mean"] = _scatter_rows(res["resid_mean"], idx, N, axis=0)
            res["resid_var"] = _scatter_rows(res["resid_var"], idx, N, axis=0)
        return res

    if layout is not None:
        # ---- Common event-time basis: higher modes first, then the coupled Mode 1 -----------
        K_eff = layout.k_eff
        partial = layout.has_partial_windows()
        coupling = cfg.mode_coupling if (partial and K_eff >= 2) else "none"
        propagate_modes = coupling == "propagated" and int(cfg.n_forecast_draws) >= 2
        if propagate_modes and m_init <= 1:
            # Higher-mode draws are propagated through Mode-1 tracks exactly like forecast
            # draws: n_forecast_draws components of draws_per_component draws per chain.
            m_init = m_max = int(cfg.n_forecast_draws)
            draws_per_chain = max(int(cfg.draws_per_component), 2)
        m_eff = m_init
        wls_df: list[np.ndarray] | None = None
        if cfg.projection == "restricted_wls":
            wls_df = restricted_wls_projection(df_dm, layout, None)[0]
        fitted: dict[int, dict[str, Any]] = {}
        fitted_m: dict[int, int] = {}
        w_units: dict[int, np.ndarray] = {}
        # Without a live correction the modes run in the paper's order (mode 1 first, so the
        # adaptive forecast-draw expansion is decided exactly as before); with partial windows
        # the higher modes come first because Mode 1 needs their posteriors.
        order = list(range(K_eff)) if coupling == "none" else list(range(K_eff - 1, 0, -1)) + [0]

        def _tau_tilde_mean() -> list[np.ndarray]:
            return [fitted[j]["tau_draws"].mean(axis=1) for j in range(1, K_eff)]

        def _stratum_order_and_means(m: int) -> tuple[list[np.ndarray], list[list[np.ndarray]]]:
            """Partition higher-mode posterior draws into ``m`` equal strata sorted by Mode-1 correction."""
            tt_all = [fitted[j]["tau_draws"].T for j in range(1, K_eff)]  # each (S_total, N)
            corr_all = mode1_dynamic_correction(layout, tt_all)  # (S_total, N)
            score = np.einsum("sn,n->s", corr_all, att_unit_weights)
            strata = [np.sort(np.asarray(idx, dtype=np.int64)) for idx in np.array_split(np.argsort(score), m)]
            comps: list[list[np.ndarray]] = []
            for idx in strata:
                comps.append([fitted[j]["tau_draws"][:, idx].mean(axis=1) for j in range(1, K_eff)])
            return strata, comps

        def _tau_tilde_draws(m: int) -> list[list[np.ndarray]]:
            """``m`` quantile-stratum means of every higher mode, per component (exact grand mean)."""
            return _stratum_order_and_means(m)[1]

        def _mode1_targets(m: int | None) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
            """Mode-1 outcome, forecast projection and (optionally) ``m`` propagated targets."""
            tt_mean = _tau_tilde_mean() if coupling != "none" else None
            if cfg.projection == "two_pass":
                base = project_event_time_mode(dy_dm, layout, 0)
                df0 = project_event_time_mode(df_dm, layout, 0)
                point = base if tt_mean is None else base - mode1_dynamic_correction(layout, tt_mean)
                draws = None
                if m is not None:
                    comps = _tau_tilde_draws(m)
                    draws = np.stack([base - mode1_dynamic_correction(layout, comp) for comp in comps], axis=0)
                return point, df0, draws
            theta, _, _ = restricted_wls_projection(dy_dm, layout, tt_mean)
            draws = None
            if m is not None:
                comps = _tau_tilde_draws(m)
                draws = np.stack([restricted_wls_projection(dy_dm, layout, comp)[0][0] for comp in comps], axis=0)
            return theta[0], wls_df[0], draws

        for k in order:
            n_trees_tau = int(cfg.num_trees_tau_mode1 if k == 0 else cfg.num_trees_tau_higher)
            dy_k_draws_m = None
            df_k_draws_m = None
            if k >= 1:
                if cfg.projection == "two_pass":
                    dy_k = project_event_time_mode(dy_dm, layout, k)
                    df_k = project_event_time_mode(df_dm, layout, k)
                    w_unit = event_time_precision_weights(layout, k)
                else:
                    tt = [fitted[j]["tau_draws"].mean(axis=1) if j in fitted else np.zeros(N) for j in range(1, K_eff)]
                    theta, var_factor, _ = restricted_wls_projection(dy_dm, layout, tt)
                    dy_k = theta[k]
                    df_k = wls_df[k]
                    w_unit = restricted_wls_precision_weights(layout, var_factor, k)
                draws_per_chain_k = m_eff * draws_per_chain
            else:
                dy_k, df_k, dy_k_draws_m = _mode1_targets(m_eff if propagate_modes else None)
                if cfg.projection == "two_pass":
                    w_unit = event_time_precision_weights(layout, 0)
                else:
                    _, var_factor, _ = restricted_wls_projection(dy_dm, layout, None)
                    w_unit = restricted_wls_precision_weights(layout, var_factor, 0)
                if df_dm_draws_all is not None and m_eff > 1:
                    df_k_draws_m = project_event_time_mode(df_dm_draws_all[:m_eff], layout, 0)
                multi = df_k_draws_m is not None or dy_k_draws_m is not None
                draws_per_chain_k = draws_per_chain if multi else (m_eff * draws_per_chain)
            if not np.all(np.isfinite(dy_k[w_unit > 0])) or not np.all(np.isfinite(df_k)):
                raise RuntimeError(f"non-finite mode-{k + 1} projection on a supported unit; this is a bug")
            dy_k = np.where(w_unit > 0, dy_k, 0.0)
            w_units[k] = w_unit
            res_k = _run_mode(
                k, dy_k=dy_k, df_k=df_k, df_k_draws_m=df_k_draws_m, w_unit=w_unit, n_trees_tau=n_trees_tau,
                draws_per_chain_k=draws_per_chain_k, idx=None, dy_k_draws_m=dy_k_draws_m,
            )

            if k == 0:
                comp_raw = res_k["comp_tau_draws"]
                comp_proj = _project_mode1_tau_draws(
                    comp_raw.reshape(-1, N).T, cohorts_arr=cohorts_arr, ever=ever, waves=waves,
                    att_unit_weights=att_unit_weights, tau_clip=cfg.tau_clip, catt_floor=cfg.catt_floor,
                ).T.reshape(comp_raw.shape)
                v_w, v_b, v_tot = _decompose_mode1_att(comp_proj)
                if (
                    df_dm_draws_all is not None
                    and m_max > m_eff
                    and v_tot > 1e-12
                    and (v_b / v_tot) > float(cfg.between_var_threshold)
                ):
                    m_eff = m_max
                    df_k_draws_m = project_event_time_mode(df_dm_draws_all[:m_eff], layout, 0)
                    if propagate_modes:
                        dy_k, _, dy_k_draws_m = _mode1_targets(m_eff)
                        dy_k = np.where(w_unit > 0, dy_k, 0.0)
                    res_k = _run_mode(
                        k, dy_k=dy_k, df_k=df_k, df_k_draws_m=df_k_draws_m, w_unit=w_unit, n_trees_tau=n_trees_tau,
                        draws_per_chain_k=draws_per_chain, idx=None, dy_k_draws_m=dy_k_draws_m,
                    )
                    comp_raw = res_k["comp_tau_draws"]
                    comp_proj = _project_mode1_tau_draws(
                        comp_raw.reshape(-1, N).T, cohorts_arr=cohorts_arr, ever=ever, waves=waves,
                        att_unit_weights=att_unit_weights, tau_clip=cfg.tau_clip, catt_floor=cfg.catt_floor,
                    ).T.reshape(comp_raw.shape)
                    v_w, v_b, v_tot = _decompose_mode1_att(comp_proj)
                att_var_within = v_w
                att_var_between = v_b
                if has_forecast:
                    f_gap_1 = float(np.sum(att_unit_weights * df_k) - np.mean(df_k[ctrl]))
                    th_var = float(np.var(res_k["theta_f_draws"], ddof=1)) if res_k["theta_f_draws"].size > 1 else 0.0
                    att_var_slope = float((f_gap_1 ** 2) * th_var)
            fitted[k] = res_k
            fitted_m[k] = int(m_eff)

        # Higher modes fitted before Mode 1 expanded the component count would hold fewer
        # draws than Mode 1: refit them at the expanded budget (adaptive expansion only).
        for k in range(1, K_eff):
            if fitted_m.get(k) == int(m_eff):
                continue
            if cfg.projection == "two_pass":
                dy_k = project_event_time_mode(dy_dm, layout, k)
                df_k = project_event_time_mode(df_dm, layout, k)
            else:
                tt = [fitted[j]["tau_draws"].mean(axis=1) for j in range(1, K_eff)]
                dy_k = restricted_wls_projection(dy_dm, layout, tt)[0][k]
                df_k = wls_df[k]
            w_unit = w_units[k]
            dy_k = np.where(w_unit > 0, dy_k, 0.0)
            fitted[k] = _run_mode(
                k, dy_k=dy_k, df_k=df_k, df_k_draws_m=None, w_unit=w_unit,
                n_trees_tau=int(cfg.num_trees_tau_higher), draws_per_chain_k=m_eff * draws_per_chain, idx=None,
            )
            fitted_m[k] = int(m_eff)

        # Pair each Mode-1 posterior draw with its within-stratum higher-mode draw and propagate
        # the mean-zero within-stratum dynamic residual when mode_coupling == "propagated".
        half_s = draws_per_chain // 2
        comp_grid = np.broadcast_to(
            np.arange(m_eff, dtype=np.int64)[None, :, None], (num_chains, m_eff, draws_per_chain)
        )
        if m_eff > 1 and half_s >= 1:
            comp_chains = np.concatenate(
                [
                    comp_grid[:, :, :half_s].reshape(num_chains, m_eff * half_s),
                    comp_grid[:, :, half_s:].reshape(num_chains, m_eff * (draws_per_chain - half_s)),
                ],
                axis=1,
            )
        else:
            comp_chains = comp_grid.reshape(num_chains, m_eff * draws_per_chain)
        comp_of_col = comp_chains.reshape(-1)

        chains_raw = fitted[0]["tau_chains"]
        chains_proj = _project_mode1_tau_draws(
            chains_raw.reshape(-1, N).T, cohorts_arr=cohorts_arr, ever=ever, waves=waves,
            att_unit_weights=att_unit_weights, tau_clip=cfg.tau_clip, catt_floor=cfg.catt_floor,
        ).T.reshape(chains_raw.shape)
        att_chains = np.einsum("csn,n->cs", chains_proj, att_unit_weights)
        rhat_att = float(compute_rhat(att_chains))
        if propagate_modes and m_eff > 1:
            att_chains_ess = att_chains.copy()
            grand_att = float(np.mean(att_chains))
            for c_idx in range(m_eff):
                sel_c = comp_chains[0] == c_idx
                att_chains_ess[:, sel_c] -= float(np.mean(att_chains[:, sel_c])) - grand_att
            ess_att = float(compute_ess(att_chains_ess, method="bulk"))
        else:
            ess_att = float(compute_ess(att_chains, method="bulk"))

        if propagate_modes and m_eff > 1:
            strata, comps = _stratum_order_and_means(m_eff)
            inv_perm = np.empty_like(comp_of_col)
            stratum_of_col = np.empty_like(comp_of_col)
            for c_idx in range(m_eff):
                inv_perm[strata[c_idx]] = np.nonzero(comp_of_col == c_idx)[0]
                stratum_of_col[strata[c_idx]] = c_idx
            res0 = dict(fitted[0])
            res0["tau_draws"] = res0["tau_draws"][:, inv_perm]
            if cfg.projection == "two_pass":
                resid_tt = [
                    fitted[j]["tau_draws"].T
                    - np.stack([comps[c_idx][j - 1] for c_idx in range(m_eff)], axis=0)[stratum_of_col]
                    for j in range(1, K_eff)
                ]
                corr_resid = mode1_dynamic_correction(layout, resid_tt).T  # (N, S_total)
                res0["tau_draws"] = res0["tau_draws"] - corr_resid
                comp_raw_new = np.empty_like(res0["comp_tau_draws"])
                for c_idx in range(m_eff):
                    comp_raw_new[c_idx] = res0["tau_draws"][:, strata[c_idx]].T
                res0["comp_tau_draws"] = comp_raw_new
            fitted[0] = res0

        for k in range(K_eff):
            res_k = fitted[k]
            trt_chains = res_k["tau_chains"][:, :, trt_idx]
            if trt_chains.shape[0] >= 2 and trt_chains.shape[1] >= 4:
                diag_k = _mode1_diagnostics(trt_chains, compute_rhat, compute_ess)
                diagnostics["rhat_mode"].append(diag_k["rhat_unit"])
                diagnostics["ess_bulk_mode"].append(diag_k["ess_bulk_unit"])
                diagnostics["ess_tail_mode"].append(diag_k["ess_tail_unit"])
            else:
                diagnostics["rhat_mode"].append(np.full(trt_idx.size, np.nan))
                diagnostics["ess_bulk_mode"].append(np.full(trt_idx.size, np.nan))
                diagnostics["ess_tail_mode"].append(np.full(trt_idx.size, np.nan))
            diagnostics["s_eta2_draws"].append(np.asarray(res_k.get("s_eta2_draws", np.zeros(0))))
            diagnostics["sigma_nu_total_draws"].append(np.asarray(res_k.get("sigma_nu_total_draws", np.zeros(0))))
            diagnostics["tau_intercept_draws"].append(np.asarray(res_k.get("tau_intercept_draws", np.zeros(0))))
            diagnostics["y_scale_mode"].append(float(res_k.get("y_scale", float("nan"))))
            mode_tau_draws.append(res_k["tau_draws"])  # (N, S_total), unprojected
            theta_f_draws.append(res_k["theta_f_draws"])
            theta_f_modes.append(float(np.mean(res_k["theta_f_draws"])))
            sigma_modes.append(float(np.mean(res_k["sigma_draws"])))
            gp_hyper_draws.append({"s_zeta": res_k["s_zeta_draws"], "ell": res_k["ell_draws"]})
        active_units = [w_units[k] > 0 for k in range(K_eff)]
        diagnostics["mode_coupling"] = coupling
        diagnostics["projection"] = cfg.projection
        diagnostics["precision_weights"] = [w_units[k] for k in range(K_eff)]
        diagnostics["n_mode_components"] = int(m_eff if propagate_modes else 1)
        diagnostics["n_forecast_components"] = int(m_eff if df_dm_draws_all is not None else 1)
        b_ctrl_et: list[float] = []
        b_trt_et: list[float] = []
        for k in range(K_eff):
            sig2_k = float(np.mean(np.asarray(fitted[k]["sigma_draws"], dtype=np.float64) ** 2))
            eta2_arr = np.asarray(fitted[k].get("s_eta2_draws", np.zeros(0)), dtype=np.float64)
            eta2_k = float(np.mean(eta2_arr)) if eta2_arr.size > 0 else 0.0
            b_ctrl_et.append(sig2_k)
            b_trt_et.append(sig2_k + eta2_k)
        arm_meta_inputs_et = {
            "y_arr": y_arr,
            "f_arr": f_arr,
            "has_forecast": has_forecast,
            "dy_dm": dy_dm,
            "df_dm": df_dm,
            "x_mod": x_mod,
            "obs": obs,
            "ctrl": ctrl,
            "t0": t0,
            "t_post": int(layout.horizon),
            "resid_means": [fitted[k]["resid_mean"] for k in range(K_eff)],
            "resid_vars": [fitted[k]["resid_var"] for k in range(K_eff)],
            "w_units": [w_units[k] for k in range(K_eff)],
            "bayes_var_ctrl": b_ctrl_et,
            "bayes_var_trt": b_trt_et,
            "treated_extra_variance": bool(cfg.treated_extra_variance),
        }
        return _finalize_result(
            mode_tau_draws_raw=mode_tau_draws,
            cohorts_arr=cohorts_arr,
            ever=ever,
            waves=waves,
            bases=bases,
            att_unit_weights=att_unit_weights,
            z_arr=z_arr,
            cfg=cfg,
            theta_f_modes=theta_f_modes,
            sigma_modes=sigma_modes,
            theta_f_draws=theta_f_draws,
            gp_hyper_draws=gp_hyper_draws,
            att_var_within=att_var_within,
            att_var_between=att_var_between,
            att_var_slope=att_var_slope,
            n_components=int(m_eff),
            rhat_att=rhat_att,
            ess_att=ess_att,
            diagnostics=diagnostics,
            variant="propagated" if (cfg.uncertainty_propagation and m_eff > 1) else "jax",
            active_units=active_units,
            t_start=t_start,
            compute_rhat=compute_rhat,
            compute_ess=compute_ess,
            result_fields=event_time_fields,
            treated_cells=treated_cells_modeled,
            arm_meta_inputs=arm_meta_inputs_et,
        )

    resid_means_pc: list[np.ndarray] = []
    resid_vars_pc: list[np.ndarray] = []
    w_units_pc: list[np.ndarray] = []
    b_ctrl_pc: list[float] = []
    b_trt_pc: list[float] = []
    for k in range(K):
        active_waves = [g for g in waves if k < bases[g].shape[1]]
        if not active_waves:
            break
        dy_k = _project_panel_to_mode(dy_dm, cohorts_arr, ctrl, active_waves, wave_weights, bases, k)
        df_k = _project_panel_to_mode(df_dm, cohorts_arr, ctrl, active_waves, wave_weights, bases, k)
        w_unit = _compute_mode_precision_weights(
            cohorts_arr, ctrl, active_waves, wave_weights, bases, t0, T, k
        )
        # Units that have mode k (every control + treated cohorts with L_g >= k + 1). Short
        # cohorts are excluded from the mode-k fit (Issue 2); mode 1 always has every unit.
        active = _mode_active_mask(cohorts_arr, ctrl, active_waves)
        active_units.append(active)
        idx = None if bool(np.all(active)) else np.nonzero(active)[0]
        n_trees_tau = int(cfg.num_trees_tau_mode1 if k == 0 else cfg.num_trees_tau_higher)

        df_k_draws_m = (
            _project_panel_to_mode(
                df_dm_draws_all[:m_eff], cohorts_arr, ctrl, active_waves, wave_weights, bases, k
            )
            if (k == 0 and df_dm_draws_all is not None and m_eff > 1)
            else None
        )
        draws_per_chain_k = draws_per_chain if df_k_draws_m is not None else (m_eff * draws_per_chain)

        res_k = _run_mode(
            k, dy_k=dy_k, df_k=df_k, df_k_draws_m=df_k_draws_m, w_unit=w_unit, n_trees_tau=n_trees_tau,
            draws_per_chain_k=draws_per_chain_k, idx=idx,
        )

        # Adaptive component expansion determined on Mode 1 (k == 0) only,
        # preserving exact K-invariance of ATT and unit-averaged effects (Proposition 2)
        if k == 0:
            comp_raw = res_k["comp_tau_draws"]
            comp_proj = _project_mode1_tau_draws(
                comp_raw.reshape(-1, N).T,
                cohorts_arr=cohorts_arr,
                ever=ever,
                waves=waves,
                att_unit_weights=att_unit_weights,
                tau_clip=cfg.tau_clip,
                catt_floor=cfg.catt_floor,
            ).T.reshape(comp_raw.shape)
            v_w, v_b, v_tot = _decompose_mode1_att(comp_proj)
            if (
                df_dm_draws_all is not None
                and m_max > m_eff
                and v_tot > 1e-12
                and (v_b / v_tot) > float(cfg.between_var_threshold)
            ):
                m_eff = m_max
                df_k_draws_m = _project_panel_to_mode(
                    df_dm_draws_all[:m_eff], cohorts_arr, ctrl, active_waves, wave_weights, bases, k
                )
                res_k = _run_mode(
                    k, dy_k=dy_k, df_k=df_k, df_k_draws_m=df_k_draws_m, w_unit=w_unit, n_trees_tau=n_trees_tau,
                    draws_per_chain_k=draws_per_chain, idx=idx,
                )
                comp_raw = res_k["comp_tau_draws"]
                comp_proj = _project_mode1_tau_draws(
                    comp_raw.reshape(-1, N).T,
                    cohorts_arr=cohorts_arr,
                    ever=ever,
                    waves=waves,
                    att_unit_weights=att_unit_weights,
                    tau_clip=cfg.tau_clip,
                    catt_floor=cfg.catt_floor,
                ).T.reshape(comp_raw.shape)
                v_w, v_b, v_tot = _decompose_mode1_att(comp_proj)

            att_var_within = v_w
            att_var_between = v_b
            if has_forecast:
                f_gap_1 = float(np.sum(att_unit_weights * df_k) - np.mean(df_k[ctrl]))
                th_var = float(np.var(res_k["theta_f_draws"], ddof=1)) if res_k["theta_f_draws"].size > 1 else 0.0
                att_var_slope = float((f_gap_1 ** 2) * th_var)
            chains_raw = res_k["tau_chains"]
            chains_proj = _project_mode1_tau_draws(
                chains_raw.reshape(-1, N).T,
                cohorts_arr=cohorts_arr,
                ever=ever,
                waves=waves,
                att_unit_weights=att_unit_weights,
                tau_clip=cfg.tau_clip,
                catt_floor=cfg.catt_floor,
            ).T.reshape(chains_raw.shape)
            att_chains = np.einsum("csn,n->cs", chains_proj, att_unit_weights)
            rhat_att = float(compute_rhat(att_chains))
            ess_att = float(compute_ess(att_chains, method="bulk"))

        # Per-mode convergence summaries of the treated units' draws (chains, draws, units)
        trt_chains = res_k["tau_chains"][:, :, trt_idx]
        if trt_chains.shape[0] >= 2 and trt_chains.shape[1] >= 4:
            diag_k = _mode1_diagnostics(trt_chains, compute_rhat, compute_ess)
            diagnostics["rhat_mode"].append(diag_k["rhat_unit"])
            diagnostics["ess_bulk_mode"].append(diag_k["ess_bulk_unit"])
            diagnostics["ess_tail_mode"].append(diag_k["ess_tail_unit"])
        else:
            diagnostics["rhat_mode"].append(np.full(trt_idx.size, np.nan))
            diagnostics["ess_bulk_mode"].append(np.full(trt_idx.size, np.nan))
            diagnostics["ess_tail_mode"].append(np.full(trt_idx.size, np.nan))
        diagnostics["s_eta2_draws"].append(np.asarray(res_k.get("s_eta2_draws", np.zeros(0))))
        diagnostics["sigma_nu_total_draws"].append(np.asarray(res_k.get("sigma_nu_total_draws", np.zeros(0))))
        diagnostics["tau_intercept_draws"].append(np.asarray(res_k.get("tau_intercept_draws", np.zeros(0))))
        diagnostics["y_scale_mode"].append(float(res_k.get("y_scale", float("nan"))))

        mode_tau_draws.append(res_k["tau_draws"])  # (N, S_total), unprojected
        theta_f_draws.append(res_k["theta_f_draws"])
        theta_f_modes.append(float(np.mean(res_k["theta_f_draws"])))
        sigma_modes.append(float(np.mean(res_k["sigma_draws"])))
        resid_means_pc.append(res_k["resid_mean"])
        resid_vars_pc.append(res_k["resid_var"])
        w_units_pc.append(np.where(active, w_unit, 0.0))
        sig2_pc_k = float(np.mean(np.asarray(res_k["sigma_draws"], dtype=np.float64) ** 2))
        eta2_pc_arr = np.asarray(res_k.get("s_eta2_draws", np.zeros(0)), dtype=np.float64)
        eta2_pc_k = float(np.mean(eta2_pc_arr)) if eta2_pc_arr.size > 0 else 0.0
        b_ctrl_pc.append(sig2_pc_k)
        b_trt_pc.append(sig2_pc_k + eta2_pc_k)
        gp_hyper_draws.append(
            {
                "s_zeta": res_k["s_zeta_draws"],
                "ell": res_k["ell_draws"],
            }
        )

    arm_meta_inputs_pc = {
        "y_arr": y_arr,
        "f_arr": f_arr,
        "has_forecast": has_forecast,
        "dy_dm": dy_dm,
        "df_dm": df_dm,
        "x_mod": x_mod,
        "obs": obs,
        "ctrl": ctrl,
        "t0": t0,
        "t_post": T - t0,
        "resid_means": resid_means_pc,
        "resid_vars": resid_vars_pc,
        "w_units": w_units_pc,
        "bayes_var_ctrl": b_ctrl_pc,
        "bayes_var_trt": b_trt_pc,
        "treated_extra_variance": bool(cfg.treated_extra_variance),
    }
    return _finalize_result(
        mode_tau_draws_raw=mode_tau_draws,
        cohorts_arr=cohorts_arr,
        ever=ever,
        waves=waves,
        bases=bases,
        att_unit_weights=att_unit_weights,
        z_arr=z_arr,
        cfg=cfg,
        theta_f_modes=theta_f_modes,
        sigma_modes=sigma_modes,
        theta_f_draws=theta_f_draws,
        gp_hyper_draws=gp_hyper_draws,
        att_var_within=att_var_within,
        att_var_between=att_var_between,
        att_var_slope=att_var_slope,
        n_components=int(m_eff),
        rhat_att=rhat_att,
        ess_att=ess_att,
        diagnostics=diagnostics,
        variant="propagated" if (cfg.uncertainty_propagation and m_eff > 1) else "jax",
        active_units=active_units,
        t_start=t_start,
        compute_rhat=compute_rhat,
        compute_ess=compute_ess,
        arm_meta_inputs=arm_meta_inputs_pc,
    )



@dataclass
class PlaceboAuditResult:
    """Result of :func:`rank_k_placebo_audit` (supports both attribute and dict-style indexing).

    Attributes:
      placebo_att_hat: Posterior mean of the placebo ATT (true value zero).
      placebo_att_se: Posterior SD of the placebo ATT.
      placebo_ci_low: Lower end of the 95% placebo ATT interval (2.5% posterior percentile).
      placebo_ci_high: Upper end of the 95% placebo ATT interval (97.5% posterior percentile).
      placebo_z_score: ``placebo_att_hat / placebo_att_se`` (alias ``placebo_z_stat``).
      covers_zero: Whether the 95% interval covers 0 (aliases ``covers_zero_95`` and
        ``passes_audit``).
      flag_parallel_trends_violation: ``not covers_zero``.
      result: The full :class:`RankKLongBetResult` of the placebo fit.
      placebo_anchor_t: Pseudo first-adoption period within the placebo panel.
      placebo_horizon: Number of placebo post-adoption periods of the first wave,
        ``max(pseudo_anchor_offset, 2)``.
    """

    placebo_att_hat: float
    placebo_att_se: float
    placebo_ci_low: float
    placebo_ci_high: float
    placebo_z_score: float
    covers_zero: bool
    flag_parallel_trends_violation: bool
    result: RankKLongBetResult
    placebo_anchor_t: int = 0
    placebo_horizon: int = 0

    @property
    def placebo_z_stat(self) -> float:
        return self.placebo_z_score

    @property
    def covers_zero_95(self) -> bool:
        return self.covers_zero

    @property
    def passes_audit(self) -> bool:
        return self.covers_zero

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def get(self, key: str, default: Any = None) -> Any:
        return getattr(self, key, default)

    def to_dict(self) -> dict[str, Any]:
        return {
            "placebo_anchor_t": self.placebo_anchor_t,
            "placebo_horizon": self.placebo_horizon,
            "placebo_att_hat": self.placebo_att_hat,
            "placebo_att_se": self.placebo_att_se,
            "placebo_ci_low": self.placebo_ci_low,
            "placebo_ci_high": self.placebo_ci_high,
            "placebo_z_score": self.placebo_z_score,
            "placebo_z_stat": self.placebo_z_score,
            "covers_zero": self.covers_zero,
            "covers_zero_95": self.covers_zero,
            "passes_audit": self.covers_zero,
            "flag_parallel_trends_violation": self.flag_parallel_trends_violation,
            "result": self.result,
        }


def rank_k_placebo_audit(
    y: np.ndarray,
    z: np.ndarray,
    x: np.ndarray,
    *,
    y_pre_history: np.ndarray | None = None,
    features: Any | None = None,
    coords: np.ndarray | None = None,
    cluster_ids: np.ndarray | None = None,
    cohorts: np.ndarray | None = None,
    pseudo_anchor_offset: int = 4,
    cfg: RankKLongBetConfig | None = None,
) -> PlaceboAuditResult:
    """Pre-adoption placebo audit: refit Rank-K on a pseudo-rollout whose true effect is zero.

    Let ``t0`` be the first adoption period, ``Y_pre = [y_pre_history | y[:, :t0]]``
    the ``L_pre`` pre-adoption periods of every unit, and
    ``h = max(pseudo_anchor_offset, 2)``. Then:

    * the placebo panel is the last ``T_pl = min(L_pre, h + max(t0, 2))``
      periods of ``Y_pre`` (the last ``pseudo_anchor_offset + t0`` periods
      whenever ``t0 >= 2``, ``pseudo_anchor_offset >= 2`` and the history is long
      enough), and the earlier periods of ``Y_pre`` become its history;
    * the pseudo first adoption is ``t0_pl = T_pl - h``: the final ``h``
      pre-adoption periods form the placebo post-period of the first wave;
    * every ever-treated unit keeps its adoption lag relative to the first wave,
      ``G_i - t0`` (capped at the last placebo period); never-treated units stay
      never-treated.

    All placebo outcomes precede the real rollout, so the true effect is zero by
    construction. Features are rebuilt on the placebo panel with
    :func:`build_causal_features` (default forecaster, ``seed=cfg.random_seed``,
    ``n_factor_loadings=cfg.n_factor_loadings``) and :func:`fit_rank_k_longbet`
    is run with ``cfg``, ``coords`` and ``cluster_ids``. The audit passes
    (``passes_audit``) iff the 95% placebo ATT interval covers 0; otherwise
    ``flag_parallel_trends_violation`` is set, pointing at a parallel-trends or
    trajectory-selection problem. The placebo ATT depends on mode 1 only, so it
    is invariant to ``cfg.rank_k`` by construction (``rank_k=1`` gives the same
    placebo ATT at lower cost).

    Args:
      y: Outcome panel ``(N, T)`` passed to :func:`fit_rank_k_longbet`.
      z: Binary absorbing treatment panel ``(N, T)``.
      x: Pre-treatment unit covariates ``(N, P)``.
      y_pre_history: Optional ``(N, T_hist)`` outcomes observed before period 0
        of ``y``; needed when ``L_pre = T_hist + t0 < h + 2``.
      features: Ignored; accepted only for API symmetry with
        :func:`fit_rank_k_longbet`. The features are always rebuilt on the
        placebo panel with the default forecaster, so a custom forecaster or
        other options used to build ``features`` are not exercised by the audit.
        Passing a non-``None`` value emits a :class:`RankKConfigWarning`.
      coords: Optional ``(N, 2)`` (latitude, longitude) in degrees.
      cluster_ids: Optional length-``N`` cluster labels.
      cohorts: Optional first-adoption periods, as in :func:`fit_rank_k_longbet`.
      pseudo_anchor_offset: Whole number of final pre-adoption periods held out as
        the placebo post-period (default 4; values below 2 are raised to 2).
      cfg: Configuration of the placebo fit (default ``RankKLongBetConfig()``,
        the propagated variant); use the configuration of the audited fit.

    Returns:
      A :class:`PlaceboAuditResult` with ``placebo_anchor_t = t0_pl`` and
      ``placebo_horizon = h``.

    Raises:
      ValueError: For invalid inputs (as in :func:`fit_rank_k_longbet`, plus a
        misshapen or non-finite ``y_pre_history`` and a ``pseudo_anchor_offset``
        that is not a whole number) or when ``L_pre < h + 2``.
    """
    if cfg is None:
        cfg = RankKLongBetConfig()
    offset = pseudo_anchor_offset
    if isinstance(offset, (bool, np.bool_)) or not (
        isinstance(offset, (int, np.integer))
        or (isinstance(offset, (float, np.floating)) and float(offset).is_integer())
    ):
        # int() would silently truncate 2.5 to 2.
        raise ValueError(f"pseudo_anchor_offset must be a whole number of periods; got {offset!r}")
    # Fail fast, before the feature rebuild; warnings come from the placebo fit.
    structure = _validate_rank_k_inputs(
        y, z, x, cfg=cfg, coords=coords, cluster_ids=cluster_ids, cohorts=cohorts, warn=False
    )
    if structure["masked"]:
        raise ValueError(
            "rank_k_placebo_audit does not support masked panels (NaN in y): the placebo features are rebuilt "
            "with build_causal_features on the pre-adoption window, which needs complete outcomes. Run the "
            "audit on the complete pre-period rows (e.g. the common observed pre window)."
        )
    if features is not None:
        # The placebo rebuilds features strictly on the pre-adoption window.
        _warn_rank_k(
            "rank_k_placebo_audit ignores `features`: it always rebuilds the features on the placebo "
            "window with build_causal_features and the default forecaster, so the options used to "
            "build `features` (for example a custom forecaster) are not exercised by the audit. "
            "Omit `features` to silence this warning."
        )
    y_arr = np.asarray(y, dtype=np.float64)
    z_arr = np.asarray(z, dtype=np.int32)
    N, T = y_arr.shape
    cohorts_arr = _infer_cohorts(z_arr) if cohorts is None else np.asarray(cohorts, dtype=np.int32)
    ever = cohorts_arr < T
    if not np.any(ever) or np.all(ever):
        raise ValueError("Panel must contain both ever-treated and never-treated units.")
    t0 = int(np.min(cohorts_arr[ever]))
    pre_y = y_arr[:, :t0]
    if y_pre_history is not None:
        yh = _as_float_array(y_pre_history, "y_pre_history")
        if yh.ndim != 2 or yh.shape[0] != N:
            raise ValueError(f"y_pre_history must be a 2-D (N, T_hist) panel with N={N} rows; got shape {yh.shape}")
        _require_finite(yh, "y_pre_history")
        full_pre = np.concatenate([yh, pre_y], axis=1)
    else:
        full_pre = pre_y

    L_pre = full_pre.shape[1]
    h_post = max(int(pseudo_anchor_offset), 2)
    if L_pre < h_post + 2:
        raise ValueError(
            f"Pre-adoption history length L_pre={L_pre} (T_hist + first adoption period) is too short "
            f"for pseudo_anchor_offset={h_post}: at least {h_post + 2} periods are needed; pass a longer "
            "y_pre_history or a smaller pseudo_anchor_offset."
        )
    t_panel_pl = min(L_pre, h_post + max(t0, 2))
    t0_pl = t_panel_pl - h_post
    hist_len_pl = L_pre - t_panel_pl

    y_pl = full_pre[:, hist_len_pl:]
    yh_pl = full_pre[:, :hist_len_pl] if hist_len_pl > 0 else None

    cohorts_pl = np.full(N, t_panel_pl, dtype=np.int32)
    for i in np.nonzero(ever)[0]:
        offset_i = int(cohorts_arr[i]) - t0
        cohorts_pl[i] = min(t0_pl + offset_i, t_panel_pl - 1)
    z_pl = (np.arange(t_panel_pl)[None, :] >= cohorts_pl[:, None]).astype(np.int32)

    feats_pl = build_causal_features(
        y_pl,
        z_pl,
        x=x,
        y_pre_history=yh_pl,
        seed=int(cfg.random_seed),
        n_factor_loadings=int(cfg.n_factor_loadings),
    )
    # The placebo ATT depends on mode 1 only (Proposition 2), so the audit runs with
    # rank_k=1: identical placebo ATT draws at lower cost, and the capped pseudo-cohorts
    # (a single placebo period) can no longer fall short of the number of modes.
    cfg_pl = dc_replace(cfg, rank_k=1) if int(cfg.rank_k) != 1 else cfg
    res_pl = _fit_rank_k_longbet(
        y_pl,
        z_pl,
        x,
        features=feats_pl,
        coords=coords,
        cluster_ids=cluster_ids,
        cohorts=cohorts_pl,
        t_0_panel=t0_pl,
        cfg=cfg_pl,
        warn_short_cohorts=False,
    )
    covers_zero = bool(res_pl.att_ci_low <= 0.0 <= res_pl.att_ci_high)
    z_score = float(res_pl.att_hat / max(res_pl.att_se, 1e-12))
    return PlaceboAuditResult(
        placebo_att_hat=float(res_pl.att_hat),
        placebo_att_se=float(res_pl.att_se),
        placebo_ci_low=float(res_pl.att_ci_low),
        placebo_ci_high=float(res_pl.att_ci_high),
        placebo_z_score=z_score,
        covers_zero=covers_zero,
        flag_parallel_trends_violation=bool(not covers_zero),
        result=res_pl,
        placebo_anchor_t=int(t0_pl),
        placebo_horizon=int(h_post),
    )


def _infer_panel_masks_from_result(
    res: RankKLongBetResult,
    cohorts: np.ndarray | None,
    ever: np.ndarray | None,
    att_w: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resolve ``(cohorts, ever, att_w)`` from explicit arguments, the result's own fields, or ``res``.

    Results of 0.7.0 and later carry ``cohorts`` and ``unit_weights`` (the ATT
    weights behind ``att_draws``, cell- or unit-weighted); they are used unless the
    caller passes explicit values. Older results fall back to the single-wave
    inference from ``catt_draws``.
    """
    n_units = int(res.catt_draws.shape[1])
    n_periods = int(res.catt_draws.shape[2])
    waves = sorted(int(g) for g in res.bases.keys())
    t0 = min(waves) if waves else 1
    stored_cohorts = getattr(res, "cohorts", None)
    stored_w = getattr(res, "unit_weights", None)
    if cohorts is None and stored_cohorts is not None:
        cohorts = np.asarray(stored_cohorts, dtype=np.int32)
    if att_w is None and stored_w is not None:
        att_w = np.asarray(stored_w, dtype=np.float64)
    if ever is None:
        if cohorts is not None:
            ever_arr = np.asarray(cohorts, dtype=np.int32) < n_periods
        elif att_w is not None:
            ever_arr = np.asarray(att_w, dtype=np.float64) > 0.0
        else:
            ever_arr = np.any(np.abs(res.catt_draws) > 0.0, axis=(0, 2))
    else:
        ever_arr = np.asarray(ever, dtype=bool)
    if cohorts is None:
        cohorts_arr = np.where(ever_arr, t0, n_periods).astype(np.int32)
    else:
        cohorts_arr = np.asarray(cohorts, dtype=np.int32)
    if att_w is None:
        att_w_arr = np.zeros(n_units, dtype=np.float64)
        total_trt_cells = float(np.sum(np.maximum(n_periods - cohorts_arr[ever_arr], 0)))
        for g in waves:
            mask_g = cohorts_arr == g
            att_w_arr[mask_g] = float(n_periods - g) / max(total_trt_cells, 1.0)
    else:
        att_w_arr = np.asarray(att_w, dtype=np.float64)
    return cohorts_arr, ever_arr, att_w_arr


@dataclass
class Nested3ArmReconciliationResult:
    """Output of :func:`reconcile_nested_3arm_draws`.

    Unpacks as the 3-tuple ``(trt_catt_bc, trt_catt_ab, trt_catt_ac)`` of
    ``(S, H)`` posterior trajectory draws and exposes the window-mean Mode-1
    posterior draws ``post_att_*`` and the underlying
    :class:`MultiArmReconciliationResult` (``multi_arm_result``).
    """

    trt_catt_bc: np.ndarray
    trt_catt_ab: np.ndarray
    trt_catt_ac: np.ndarray
    post_att_bc: np.ndarray
    post_att_ab: np.ndarray
    post_att_ac: np.ndarray
    prior_scale: float | None = None
    multi_arm_result: Any = None

    def __iter__(self):
        yield self.trt_catt_bc
        yield self.trt_catt_ab
        yield self.trt_catt_ac

    def __getitem__(self, i: int) -> np.ndarray:
        return (self.trt_catt_bc, self.trt_catt_ab, self.trt_catt_ac)[i]

    def __len__(self) -> int:
        return 3

    @property
    def sampler_diagnostics(self) -> Any:
        """Delegate :attr:`sampler_diagnostics` to :attr:`multi_arm_result`."""
        return getattr(self.multi_arm_result, "sampler_diagnostics", None)

    @property
    def diagnostics(self) -> Any:
        """Delegate :attr:`diagnostics` to :attr:`multi_arm_result`."""
        return getattr(self.multi_arm_result, "diagnostics", None)

    @property
    def max_cycle_discrepancy(self) -> float:
        """Delegate :attr:`max_cycle_discrepancy` to :attr:`multi_arm_result`."""
        return float(getattr(self.multi_arm_result, "max_cycle_discrepancy", 0.0))

    @property
    def ppc_cycle_pvalue(self) -> float | None:
        """Delegate :attr:`ppc_cycle_pvalue` to :attr:`multi_arm_result`."""
        return getattr(self.multi_arm_result, "ppc_cycle_pvalue", None)

    @property
    def ppc_basis_pvalue(self) -> float | None:
        """Delegate :attr:`ppc_basis_pvalue` to :attr:`multi_arm_result`."""
        return getattr(self.multi_arm_result, "ppc_basis_pvalue", None)


def reconcile_nested_3arm_draws(
    res_bc: RankKLongBetResult,
    res_ab: RankKLongBetResult,
    res_ac: RankKLongBetResult,
    *,
    monotone: bool = True,
    prior_scale: float | None = None,
    horizon: int | None = None,
    n_modes: int = 2,
    seed: int = 2026,
    _model_options: Any = None,
) -> Nested3ArmReconciliationResult:
    """Convenience wrapper reconciling nested 3-arm fits ``(C <= B <= A)`` via :func:`reconcile_multi_arm_draws`."""
    rec = reconcile_multi_arm_draws(
        {("B", "C"): res_bc, ("A", "B"): res_ab, ("A", "C"): res_ac},
        arms=("C", "B", "A"),
        order=("C", "B", "A"),
        monotone=monotone,
        prior_scale=prior_scale,
        horizon=horizon,
        n_modes=n_modes,
        seed=seed,
        _model_options=_model_options,
    )
    return Nested3ArmReconciliationResult(
        trt_catt_bc=rec.contrast_trajectory_draws[("B", "C")],
        trt_catt_ab=rec.contrast_trajectory_draws[("A", "B")],
        trt_catt_ac=rec.contrast_trajectory_draws[("A", "C")],
        post_att_bc=rec.contrast_mode1_draws[("B", "C")],
        post_att_ab=rec.contrast_mode1_draws[("A", "B")],
        post_att_ac=rec.contrast_mode1_draws[("A", "C")],
        prior_scale=rec.prior_scale,
        multi_arm_result=rec,
    )


reconcile_nested_3arm_monthly_draws = reconcile_nested_3arm_draws


def fit_rank_k_hybrid(
    y: np.ndarray,
    z: np.ndarray,
    x: np.ndarray,
    *,
    features: Any | None = None,
    unit_features: np.ndarray | None = None,
    ps: np.ndarray | None = None,
    baseline_forecast: np.ndarray | None = None,
    baseline_forecast_draws: np.ndarray | None = None,
    pre_anchor_loadings: np.ndarray | None = None,
    coords: np.ndarray | None = None,
    cluster_ids: np.ndarray | None = None,
    cohorts: np.ndarray | None = None,
    t_0_panel: int | None = None,
    cfg: RankKLongBetConfig | None = None,
    plugin_cfg: RankKLongBetConfig | None = None,
) -> RankKLongBetResult:
    """Hybrid Rank-K BCF: propagated mode-1 level (ATT) with the plug-in moderation surface and dynamics.

    Runs the propagated JAX sampler for mode 1 (``cfg`` with ``rank_k=1``: the ATT
    draws, their nuisance-propagated intervals and any ATT prior / truncation) and
    the plug-in ``stochtree`` fit (``plugin_cfg``, default: the paper's plug-in
    budget with ``cfg.rank_k`` modes) for the unit-level heterogeneity. The hybrid
    mode-1 draws are ``ATT^(s) + [tau^plug_{i,1} - sum_j a_j tau^plug_{j,1}]`` with
    the plug-in draws randomly paired to the propagated draws (the two fits are
    independent), so the hybrid ATT draws equal the propagated ones exactly and
    the unit-level deviations and higher modes come from the plug-in forests.
    ``catt_floor`` / ``tau_clip`` of ``cfg`` are applied to the combined draws by
    the usual projection accounting.
    """
    if cfg is None:
        cfg = RankKLongBetConfig()
    if cfg.engine != "jax":
        raise ValueError("fit_rank_k_hybrid needs a propagated (engine='jax') cfg for the mode-1 level")
    cfg_p = dc_replace(cfg, rank_k=1, catt_floor=None)
    if np.any(np.isnan(np.asarray(y, dtype=np.float64))):
        raise NotImplementedError(
            "fit_rank_k_hybrid does not support masked panels (NaN in y): its heterogeneity surface comes from "
            "the plug-in stochtree engine, which fits fully observed windows only. Use fit_rank_k_longbet with "
            "the JAX engine for masked panels."
        )
    if plugin_cfg is None:
        plugin_cfg = RankKLongBetConfig(
            rank_k=int(cfg.rank_k),
            use_forecast=cfg.use_forecast,
            uncertainty_propagation=False,
            num_gfr=10,
            num_mcmc=600,
            theta_f_clip=cfg.theta_f_clip,
            random_seed=int(cfg.random_seed),
            short_cohort_policy=cfg.short_cohort_policy,
            basis_spec=cfg.basis_spec,
        )
    if plugin_cfg.engine != "stochtree" or plugin_cfg.has_att_constraint or plugin_cfg.catt_floor is not None:
        raise ValueError(
            "plugin_cfg must be a plug-in (engine='stochtree') config without tau prior / clip / catt_floor; "
            "those act on the propagated mode-1 level (cfg)"
        )
    common = dict(
        features=features, unit_features=unit_features, ps=ps, baseline_forecast=baseline_forecast,
        pre_anchor_loadings=pre_anchor_loadings, cohorts=cohorts, t_0_panel=t_0_panel,
    )
    t_start = time.perf_counter()
    res_p = fit_rank_k_longbet(
        y, z, x, baseline_forecast_draws=baseline_forecast_draws, coords=coords, cluster_ids=cluster_ids,
        cfg=cfg_p, **common,
    )
    res_q = fit_rank_k_longbet(y, z, x, cfg=plugin_cfg, **common)

    z_arr = np.asarray(z, dtype=np.int32)
    N, T = z_arr.shape
    cohorts_arr = _infer_cohorts(z_arr) if cohorts is None else np.asarray(cohorts, dtype=np.int32)
    ever = cohorts_arr < T
    waves = sorted(int(g) for g in np.unique(cohorts_arr[ever]))
    total_trt_cells = float(np.sum(z_arr == 1))
    att_w = np.zeros(N, dtype=np.float64)
    for g in waves:
        att_w[cohorts_arr == g] = float(T - g) / max(total_trt_cells, 1.0)
    if res_p.unit_weights is not None and cfg.att_weighting != "cell":
        att_w = np.asarray(res_p.unit_weights, dtype=np.float64)

    att_p = np.asarray(res_p.att_draws, dtype=np.float64)
    s_p = int(att_p.shape[0])
    q1 = np.asarray(res_q.mode_tau_draws[0], dtype=np.float64)
    s_q = int(q1.shape[1])
    rng = np.random.default_rng(int(cfg.random_seed) + 7)
    if s_q >= s_p:
        sel = rng.permutation(s_q)[:s_p]
    else:
        _warn_rank_k(f"fit_rank_k_hybrid: the plug-in fit has {s_q} draws, fewer than the {s_p} propagated draws; its "
            "draws are recycled with replacement. Raise plugin_cfg.num_mcmc.")
        sel = rng.choice(s_q, size=s_p, replace=True)
    q1_sel = q1[:, sel]
    q1_centered = q1_sel - np.einsum("ns,n->s", q1_sel, att_w)[None, :]
    mode1_hybrid = q1_centered + att_p[None, :]
    mode_tau_raw = [mode1_hybrid, *[np.asarray(td, dtype=np.float64)[:, sel] for td in res_q.mode_tau_draws[1:]]]
    compute_rhat, compute_ess = _get_rhat_ess()
    diagnostics = {
        "engine": "hybrid",
        "propagated": res_p.diagnostics,
        "plugin_draw_index": sel,
        "propagated_result": res_p,
        "plugin_result": res_q,
    }
    return _finalize_result(
        mode_tau_draws_raw=mode_tau_raw,
        cohorts_arr=cohorts_arr,
        ever=ever,
        waves=waves,
        bases=res_q.bases,
        att_unit_weights=att_w,
        z_arr=z_arr,
        cfg=cfg,
        theta_f_modes=[res_p.theta_f_modes[0], *res_q.theta_f_modes[1:]],
        sigma_modes=[res_p.sigma_modes[0], *res_q.sigma_modes[1:]],
        theta_f_draws=[res_p.theta_f_draws[0], *[np.asarray(d)[sel] for d in res_q.theta_f_draws[1:]]],
        gp_hyper_draws=[
            res_p.gp_hyper_draws[0],
            *[{k_: np.asarray(v_)[sel] for k_, v_ in g.items()} for g in res_q.gp_hyper_draws[1:]],
        ],
        att_var_within=res_p.att_var_within,
        att_var_between=res_p.att_var_between,
        att_var_slope=res_p.att_var_slope,
        n_components=res_p.n_components,
        rhat_att=res_p.rhat_att,
        ess_att=res_p.ess_att,
        diagnostics=diagnostics,
        variant="hybrid",
        active_units=list(res_q.active_units),
        t_start=t_start,
        compute_rhat=compute_rhat,
        compute_ess=compute_ess,
        arm_metadata_override=res_p.arm_metadata or res_q.arm_metadata,
    )


from longbet._rank_k_multi import (  # noqa: E402
    CrossOutcomeCorrelationResult,
    MultiArmMultiOutcomeReconciliationResult,
    Nested3ArmMultiReconciliationResult,
    RankKConvergenceCheck,
    check_rank_k_convergence,
    estimate_cross_outcome_correlation,
    ge_revenue_draws,
    paid_winsor_cap,
    reconcile_multi_arm_multi_outcome_draws,
    reconcile_nested_3arm_multi_draws,
)
from longbet._rank_k_multi_arm import (  # noqa: E402
    ContrastReconciliationResult,
    MultiArmReconciliationResult,
    RankKDesignWarning,
    SamplerDiagnostics,
    reconcile_multi_arm_draws,
)
