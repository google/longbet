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

"""Configuration for the LongBet sampler."""

from __future__ import annotations

import dataclasses
import math
import warnings
from dataclasses import dataclass
from numbers import Integral, Real
from typing import Any, Literal


@dataclass(frozen=True)
class LongBetConfig:
    """Configuration of the LongBet panel model and its sampler.

    The model
    ---------
    For unit ``i`` in period ``t`` with exposure ``S`` (periods since the
    absorbing treatment started, 0 before it), the mean is

    .. code::

        alpha(mu(x, t)) + beta_S * nu(x, t) * 1{S >= 1} + gamma_i

    with ``mu`` a prognostic tree ensemble, ``nu`` a treatment tree ensemble
    for covariate (and optionally calendar-time) heterogeneity, ``beta_S`` a
    Gaussian-process trajectory over the exposure clock that gives the effect
    its time profile, and ``gamma_i`` a unit random intercept. The treatment
    term enters treated cells only, so the untreated surface is ``mu + gamma``
    and the effect at exposure ``S`` is ``beta_S * nu``. Continuous outcomes
    have Gaussian innovations; binary and ordinal outcomes use a probit latent
    response.

    The sampler
    -----------
    Every sweep updates both forests with GROW, PRUNE, CHANGE and one
    data-driven REGROW proposal per forest, then the exposure trajectory, the
    unit intercepts, the variances, a Metropolis move along the exact scale
    ridge between ``beta`` and the treatment leaves, and an exact conjugate
    subspace Gibbs transfer between ``mu`` and ``nu``. When the calendar-time
    block is on, the trend coefficients, the unit intercepts and the exposure
    trajectory are then drawn together from their exact joint Gaussian
    conditional. Parallel tempering (``tempering_levels``) is available for
    posteriors whose forest modes the local moves connect too slowly. Chains
    start overdispersed from the prior so that R-hat measures convergence
    rather than seed agreement.

    Defaults
    --------
    Four chains, 2,000 burn-in sweeps and 1,000 retained sweeps per chain.
    Options left unset (``None``) take the fitted model's own default:
    ``n_skip`` and the encouragement-model extensions (``linear_trt``,
    ``learn_beta_scale``, ``ar1_errors``, ``interweave_gamma``). ``LongBet``
    and ``LongBetMulti`` resolve them to the published model
    (:data:`MODEL_DEFAULTS`: no extension, no thinning); in
    ``LongBetEncourage`` the extensions are part of the model and ``n_skip``
    thins by 2.
    Four thousand retained draws are what the package's own convergence gate
    (bulk and tail ESS of at least 400, rank R-hat at most 1.01 on every
    reported quantity) needs to be attainable; check ``att_stability`` or the
    encouragement model's ``stability()`` on your own panel rather than
    assuming it. ``jax`` dispatches asynchronously, so time a fit with
    ``jax.block_until_ready``.

    Parameters
    ----------
    num_sweeps
        Retained draws per chain after burn-in.
    num_burnin
        Burn-in sweeps discarded per chain.
    n_skip
        Thinning interval; ``1`` keeps every post-burn-in sweep (the default of
        ``LongBet``; ``LongBetEncourage`` thins by 2). Total sweeps
        are ``num_burnin + n_skip * num_sweeps``.
    num_chains
        Parallel chains; at least two are needed for R-hat.
    inner_loop_length
        Sweeps per outer dispatch. Bounds compile time; ``None`` runs the whole
        chain in one dispatch.
    tempering_levels
        Parallel-tempering replicas per chain. ``1`` (default) runs the plain
        sampler. With ``K > 1`` every chain is a ladder of ``K`` replicas whose
        likelihood is raised to powers equally spaced in ``sqrt(beta)`` from 1
        down to ``tempering_beta_min``; replicas run the same sweep, adjacent
        levels exchange states after every sweep, and only the posterior
        replicas are retained. Exact, and lets chains cross between forest
        modes that single-tree moves cannot connect, at ``K`` times the cost.
        A ladder of about ``2 * sd(log-likelihood) * (1 - sqrt(beta_min))``
        levels keeps the exchange rate near a third. See ``longbet._tempering``.
    tempering_beta_min
        Inverse temperature of the hottest replica, in ``(0, 1]``.
    use_inter_ensemble_move
        Whether to run the exact conjugate subspace Gibbs transfer between
        ``mu`` and ``nu`` each sweep. The two ensembles are only jointly
        identified through the treatment indicator, so a component of the fit
        can move between them at nearly constant likelihood; this step samples
        along that subspace in closed form instead of leaving the chain to
        random-walk down it.
    use_trend_block
        Whether the prognostic surface carries an explicit conjugate
        calendar-time component,

        ``sum_k psi_k(t) f_k + sum_m phi_m(t) b(X_i)' c_m``,

        drawn jointly with the unit intercepts and the exposure trajectory from
        their exact Gaussian conditional. The first term is a free common-time
        factor, the second the heterogeneous trend; they are orthogonal by
        construction. See ``longbet._trend_basis``.

        The block exists because a sum of trees is piecewise constant in
        calendar time while the trajectory ``beta_S`` has one free coordinate
        per period, and on a staggered panel the identity
        ``t = (E_i - 1) + S_it`` ties the two together. The mismatch leaves a
        chain-specific remainder that no leaf shift can cancel. Note that the
        block only *removes* the mismatch if the forest is not also fitting
        calendar time -- otherwise the two compete for the same signal and a new
        ridge appears in place of the old one. Use it together with
        ``split_calendar_mu = False``.
    trend_degree
        Number of orthonormal polynomial basis functions in the heterogeneous
        trend. ``1`` is a linear trend, ``2`` (the default) adds curvature. This
        governs only how a unit's trend may *deviate* from the common one; the
        common time factor is unrestricted whenever
        ``trend_period_effects`` is on.
    trend_num_features
        Number of nonlinear dictionary columns ``F`` (default ``64``) appended
        to the standardized linear covariates in the unit feature map ``b(X_i)``
        (1D piecewise-linear hinges, 2D bilinear products, and directed 2D hinge
        interactions, truncated or padded with isotropic random Fourier features
        from fixed seed ``20260918`` and blended with ``+0.05`` RFF). They let
        the trend block carry a trend whose covariate loading is nonlinear.
        Setting ``trend_num_features = 0`` restricts ``b(X_i)`` to the linear
        standardized covariates: when true differential trends are linear in
        ``X_i`` and selection is concentrated in the tails of ``X_i``, ``0``
        cuts finite-sample regularization bias in the ATT by more than half and
        mixes ``2x--3x`` faster, at the cost of bias when differential trends
        load nonlinearly on ``X_i``.
    trend_prior_scale
        Prior standard deviation of a unit's trend *deviation* on the
        standardized response. The per-coefficient prior is rescaled by the
        basis size, so this stays interpretable as the basis grows.
    trend_period_effects
        Whether to include the free common-time factor ``sum_k psi_k(t) f_k``,
        a complete centred basis over calendar periods. This is what replaces
        the forest's calendar-time step functions when ``split_calendar_mu`` is
        off, so turning it off there leaves the model with no way to fit a
        common shock.
    trend_period_scale
        Prior standard deviation of the common time factor on the standardized
        response. Deliberately looser than ``trend_prior_scale`` and *not*
        rescaled by the basis size: every period has data, so the factor is well
        identified and should not be shrunk hard.
    trend_horseshoe
        Whether the *interaction* half of the trend block gets a horseshoe
        prior instead of a ridge. The block asks order 10^2 columns which
        covariates carry a diverging baseline trend, and typically only a
        handful do; one ridge scale must either shrink those or let the rest
        absorb noise. The horseshoe's global scale is ``trend_prior_scale``
        divided by the square root of the basis size -- the ridge value -- so
        the two priors agree on typical magnitude and differ only in tail
        weight. The common-time factor is never horseshoe'd. Sampled by the
        inverse-Gamma scheme of Makalic and Schmidt (2016), which keeps the
        model conditionally Gaussian in the coefficients and so leaves the
        exact joint ``(c, gamma, beta)`` draw intact.
    baseline_cells
        Which observed cells inform the *baseline module* -- the prognostic
        forest ``mu``, the calendar-time block, the unit intercepts
        ``gamma_i`` and the innovation variance ``sigma^2``.

        ``'all'`` (the default) is the single joint posterior. Every observed
        cell enters the baseline conditionals (treated cells through the
        partial residual ``y - beta_S nu``), and ``trend_block_step`` draws the
        trend coefficients ``c``, unit intercepts ``gamma``, and exposure
        trajectory ``beta`` from their exact joint Gaussian conditional given
        the trees. When adoption is strongly selective in the tails of ``X_i``
        and effects are non-null, horseshoe shrinkage on the ``d * (p_0 + F)``
        trend-interaction columns and BART shrinkage on ``mu`` leave a modest
        finite-sample regularization bias (reduced when ``trend_num_features =
        0``).

        ``'untreated'`` fits the baseline module **only from untreated cells**
        (``z = 0``: never-treated units and the pre-adoption periods of
        adopters) and conditions the treatment module ``(beta, nu)`` on the
        resulting baseline draw, targeting the modular (cut) posterior
        ``p(baseline | Y_0) p(beta, nu | baseline, Y_1)`` (Plummer, 2015;
        cf. Borusyak, Jaravel and Spiess, 2024; Liu, Wang and Xu, 2024). The
        sweep alternates baseline draws on ``z = 0`` cells with treatment draws
        on ``z = 1`` cells and disables inter-ensemble transfer moves. Note
        that because treated cohorts contribute only their pre-adoption periods
        to ``eta`` and ``gamma_i``, their post-adoption baseline is an
        out-of-sample polynomial extrapolation over ``t >= G_i``; on short
        panels with covariate-driven selection or unit-specific seasonality
        without a forecast channel, shrinkage on the pre-adoption slope
        estimate increases rather than reduces finite-sample ATT bias, so
        ``'untreated'`` is provided for modular sensitivity analysis rather than
        as the default. Requires a continuous outcome and no tempering, and is
        not available for ``LongBetMulti``.
    num_trees_pr, num_trees_trt
        Trees in the prognostic and treatment forests. The treatment forest
        needs more trees than a prognostic one of the same size: with 20 trees
        each tree carries a large share of a threshold surface and the chains
        freeze at their own cutpoints; with 60 the same panel mixes.
    min_points_per_leaf_pr, min_points_per_leaf_trt
        Minimum observed cells in a leaf.
    max_depth_pr, max_depth_trt
        Maximum tree depth. Trees are stored as heaps of ``2 ** depth`` nodes,
        so this bounds the archive size; posterior trees under the depth prior
        rarely exceed depth 5.
    alpha_split_pr, beta_split_pr, alpha_split_trt, beta_split_trt
        Tree-depth prior: a node at depth ``d`` splits with probability
        ``alpha / (1 + d) ** beta``. The treatment prior is shallower than the
        prognostic one, as in Bayesian causal forests.
    num_cutpoints
        Maximum quantile bins per continuous covariate.
    sig_knl, lambda_knl, kernel_type
        Marginal standard deviation, lengthscale and family (``'se'``,
        ``'matern32'`` (alias ``'matern'``), ``'matern52'``, ``'ar1'``) of the
        exposure-trajectory kernel.
    sigma_m, gp_constant_mean
        Prior standard deviation of the trajectory's constant mean, and whether
        that marginalized mean is included, so projections beyond the fitted
        horizon revert to an estimated common level rather than to zero.
    split_calendar_mu
        Whether the prognostic forest may split on calendar time. ``None``
        (the default) derives it, and deriving it is the point: this is not an
        independent choice but a consequence of whether anything else in the
        model can carry calendar time. With the calendar-time block on, the
        period basis already spans every function of ``t``, so tree splits on
        ``t`` are redundant and actively harmful -- the two overlap, the design
        goes near-collinear, and the chains disagree about which one holds the
        trend. With the block off, the trees are the only thing that can fit a
        common shock and must keep the splits. Setting this explicitly is
        supported for ablation; the combination that leaves nothing able to fit
        calendar time raises.
    split_calendar_trt
        Whether the treatment forest may split on calendar time.
    split_exposure_trt
        Whether the treatment forest may also split on the exposure clock.
        Off by default: the trajectory ``beta_S`` carries the whole exposure
        profile, which keeps the product ``beta_S * nu`` identified. Turn it
        on when different kinds of units need genuinely different shapes over
        exposure, and read the diagnostics: the extra freedom is an
        exposure-shape ridge between the trajectory and the forest.
    random_intercept
        Whether to fit unit random intercepts ``gamma_i``.
    gamma_prior_a, gamma_prior_b
        Inverse-gamma prior on the unit-intercept variance.
    sigma_prior_a, sigma_prior_b
        Inverse-gamma prior on the innovation variance of a continuous outcome,
        on the standardized scale. The default IG(2, 1) is proper with prior
        mean 1; an improper reference prior can give an improper posterior in
        coupled models and is not offered.
    outcome
        ``'continuous'``, ``'binary'`` or ``'ordinal'``.
    num_categories
        Number of ordered categories for an ordinal outcome.
    cutpoint_prior_scale
        Ordered-normal prior scale of the ordinal cutpoints, latent units.
    standardize
        Whether to centre and scale a continuous response internally; priors
        are on the standardized scale and outputs are rescaled.
    random_seed
        Base PRNG seed.
    device
        ``'auto'``, ``'cpu'`` or ``'gpu'``.
    sur
        Whether ``LongBetMulti`` couples the outcomes' within-period innovations
        through a triangular seemingly-unrelated-regression likelihood. Inert
        for a scalar fit.
    sur_prior_var
        Prior variance of the SUR loadings.
    """

    # --- sampler -----------------------------------------------------------
    num_sweeps: int = 1000
    num_burnin: int = 2000
    n_skip: int | None = None
    num_chains: int = 4
    inner_loop_length: int | None = None
    tempering_levels: int = 1
    tempering_beta_min: float = 0.05
    use_inter_ensemble_move: bool = True
    use_trend_block: bool = True
    trend_degree: int = 2
    trend_num_features: int = 64
    trend_prior_scale: float = 2.0
    trend_period_effects: bool = True
    trend_period_scale: float = 2.0
    trend_forecast_scale: float = 1.0
    trend_horseshoe: bool = True
    baseline_cells: Literal["all", "untreated"] = "all"

    # --- forests -----------------------------------------------------------
    num_trees_pr: int = 20
    num_trees_trt: int = 60
    min_points_per_leaf_pr: int = 10
    min_points_per_leaf_trt: int = 10
    max_depth_pr: int = 8
    max_depth_trt: int = 8
    alpha_split_pr: float = 0.95
    beta_split_pr: float = 2.0
    alpha_split_trt: float = 0.25
    beta_split_trt: float = 3.0
    num_cutpoints: int = 100

    # --- exposure trajectory -----------------------------------------------
    sig_knl: float = 1.0
    lambda_knl: float = 1.0
    kernel_type: Literal["se", "matern", "matern32", "matern52", "ar1"] = "se"
    sigma_m: float = 1.0
    gp_constant_mean: bool = True

    # --- structure ---------------------------------------------------------
    split_calendar_mu: bool | None = None
    split_calendar_trt: bool = False
    split_exposure_trt: bool = False

    # --- priors ------------------------------------------------------------
    random_intercept: bool = True
    gamma_prior_a: float = 1.0
    gamma_prior_b: float = 0.1
    sigma_prior_a: float = 2.0
    sigma_prior_b: float = 1.0

    # --- outcome -----------------------------------------------------------
    outcome: Literal["continuous", "binary", "ordinal"] = "continuous"
    num_categories: int | None = None
    cutpoint_prior_scale: float = 5.0
    #: Internal: hold the exposure trajectory at one (used by the BCF
    #: reduction test where there is no exposure clock).
    sample_beta: bool = True

    # --- misc --------------------------------------------------------------
    standardize: bool = True
    random_seed: int = 0
    device: Literal["auto", "cpu", "gpu"] = "auto"

    # --- multiple outcomes -------------------------------------------------
    sur: bool = True
    sur_prior_var: float = 1.0

    # --- encouragement-v2 extensions ----------------------------------------
    # Off here, so that a plain LongBet fit is the published model (Wang,
    # Martinez and Hahn, 2024); LongBetEncourage turns them on by default
    # (longbet._encourage_model.MODEL_COMPONENTS).
    #: Add a linear term to the treatment moderator,
    #: ``nu(x) = x' delta + sum of trees``, on standardized time-invariant
    #: treatment covariates, with a horseshoe prior
    #: ``delta_j ~ N(0, tau^2 lambda_j^2)``, ``lambda_j ~ C+(0, 1)``,
    #: ``tau ~ C+(0, linear_trt_scale / sqrt(p))``: the shallow tree prior makes
    #: heterogeneity a priori tiny relative to the effect's level, and a single
    #: global scale would shrink one strong moderator with the null ones.
    linear_trt: bool | None = None
    linear_trt_scale: float = 0.5
    #: Learn the scale of the exposure trajectory: ``beta ~ N(0, s^2 K_tilde)``
    #: with ``s ~ C+(0, beta_scale_prior)``; ``s = 1`` is the fixed kernel.
    learn_beta_scale: bool | None = None
    beta_scale_prior: float = 1.0
    #: Stationary AR(1) innovations within each unit's series (continuous
    #: outcome, balanced panel, unit-level partitions): ``eps_it = phi
    #: eps_{i,t-1} + xi_it``, ``phi ~ U(-ar1_phi_max, ar1_phi_max)``.
    ar1_errors: bool | None = None
    ar1_phi_max: float = 0.95
    #: Ancillarity-sufficiency interweaving (Yu and Meng, 2011) for the unit
    #: intercept variance: after its conjugate draw given gamma, redraw the
    #: scale sigma_gamma with the standardized intercepts gamma / sigma_gamma
    #: held fixed, from its exact conditional. Same target; it mixes where the
    #: intercepts are weakly informed.
    interweave_gamma: bool | None = None
    #: Anchor the treatment moderator ensemble at unit empirical mean,
    #: ``nu(x) = 1 + (nu_trees(x) - mean(nu_trees)) + x' delta``, so ``beta_S``
    #: directly represents the population-average exposure trajectory and avoids
    #: bilinear errors-in-variables Gibbs attenuation under moderate SNR.
    anchor_moderator: bool = False

    # --- distributed-lag exposure kernel (``LongBet.fit(exposure_lags=...)``) --
    # Inert unless ``exposure_lags`` is passed to ``fit``. The cell weight is
    # then ``w_it = sum_l D_{i,t-l} kappa_l`` in place of ``b_z beta_{S_it}``,
    # with ``kappa ~ N(0, s^2 K_L)`` a Gaussian process over the lag index
    # ``l = 0..L`` built from the same kernel families as ``beta``'s.
    #: Prior standard deviation of ``kappa_l`` (standardized response per unit
    #: dose) when the scale is not learned; the typical magnitude of ``s``'s
    #: half-Cauchy prior when it is.
    lag_sig_knl: float = 1.0
    #: Lengthscale of the lag kernel, in lag units.
    lag_lambda_knl: float = 1.0
    lag_kernel_type: Literal["se", "matern", "matern32", "matern52", "ar1"] = "se"
    #: Prior standard deviation of the kernel's constant mean, and whether that
    #: marginalized mean is included (as ``sigma_m`` / ``gp_constant_mean``).
    lag_sigma_m: float = 1.0
    lag_gp_constant_mean: bool = True
    #: Learn the kernel's scale: ``kappa ~ N(0, s^2 K_L)``, ``s ~ C+(0,
    #: lag_scale_prior)``. A plain Gibbs step given ``kappa``, so it runs with
    #: or without ``anchor_moderator``.
    lag_learn_scale: bool = True
    lag_scale_prior: float = 1.0
    #: Prior standard deviation of each fixed extra column's *contribution*
    #: (``LongBet.fit(extra_cols=...)``): coefficient ``j`` has prior
    #: ``N(0, (extra_cols_prior_scale / sd_j)^2)`` with ``sd_j`` the column's
    #: standard deviation over the observed cells, and is reported on the raw
    #: column's scale.
    extra_cols_prior_scale: float = 1.0

    def __post_init__(self) -> None:
        # --- resolve split_calendar_mu before anything validates booleans ----
        # Whether the prognostic forest should be allowed to cut calendar time
        # is not a free choice: it depends on whether something else in the
        # model can carry calendar time. Getting it wrong is silently
        # catastrophic in both directions, so the default is derived rather
        # than fixed.
        #
        #   block on  -> the DCT period basis spans every function of calendar
        #                time, so tree splits on t are redundant *and* harmful:
        #                the two span overlapping subspaces, the design goes
        #                near-collinear, and chains disagree about the split
        #                (R-hat 1.10, and NaN draws on one replicate).
        #   block off -> the trees are the only thing left that can fit a
        #                common time shock, so they must keep the splits
        #                (without them, RMSE degrades ~30x at zero coverage).
        can_carry_time = bool(self.use_trend_block) and bool(self.trend_period_effects)
        if self.split_calendar_mu is None:
            object.__setattr__(self, "split_calendar_mu", not can_carry_time)
        elif self.split_calendar_mu is False and not can_carry_time:
            raise ValueError(
                "split_calendar_mu=False removes calendar time from the "
                "prognostic forest, but use_trend_block="
                f"{self.use_trend_block} and trend_period_effects="
                f"{self.trend_period_effects} leave nothing able to fit a "
                "common time shock in its place. Either leave "
                "split_calendar_mu=None (the default, which resolves this for "
                "you) or enable the calendar-time block."
            )

        for name in ("num_sweeps", "num_burnin", "n_skip", "num_chains", "num_trees_pr",
                     "num_trees_trt", "min_points_per_leaf_pr", "min_points_per_leaf_trt",
                     "max_depth_pr", "max_depth_trt", "num_cutpoints", "random_seed"):
            value = getattr(self, name)
            if value is None and name in MODEL_DEFAULTS:
                continue
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise ValueError(f"{name} must be an integer, got {value!r}")
            object.__setattr__(self, name, int(value))
        if self.inner_loop_length is not None and (
                isinstance(self.inner_loop_length, bool)
                or not isinstance(self.inner_loop_length, Integral) or self.inner_loop_length < 1):
            raise ValueError("inner_loop_length must be None or a positive integer")
        if self.num_sweeps < 1:
            raise ValueError(f"num_sweeps must be >= 1, got {self.num_sweeps}")
        if self.num_burnin < 0:
            raise ValueError(f"num_burnin must be >= 0, got {self.num_burnin}")
        if self.n_skip is not None and self.n_skip < 1:
            raise ValueError(f"n_skip must be >= 1, got {self.n_skip}")
        if self.num_chains < 1:
            raise ValueError(f"num_chains must be >= 1, got {self.num_chains}")
        if self.tempering_levels < 1:
            raise ValueError("tempering_levels must be an integer >= 1")
        if not (isinstance(self.tempering_beta_min, Real) and not isinstance(self.tempering_beta_min, bool)
                and 0 < self.tempering_beta_min <= 1):
            raise ValueError("tempering_beta_min must be in (0, 1]")
        for name in ("trend_degree", "trend_num_features"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer, got {value!r}")
            object.__setattr__(self, name, int(value))
        for name in ("trend_prior_scale", "trend_period_scale", "trend_forecast_scale"):
            value = getattr(self, name)
            if (isinstance(value, bool) or not isinstance(value, Real)
                    or not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name} must be a finite positive scalar")
            object.__setattr__(self, name, float(value))
        for name in ("num_trees_pr", "num_trees_trt", "min_points_per_leaf_pr",
                     "min_points_per_leaf_trt", "max_depth_pr", "max_depth_trt", "num_cutpoints"):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be >= 1, got {getattr(self, name)}")
        if self.num_cutpoints > 255:
            raise ValueError(
                f"num_cutpoints must be at most 255, got {self.num_cutpoints}: the design "
                "stores split indices as unsigned 8-bit integers and the REGROW proposal "
                "scores every cutpoint of every variable at every node"
            )
        for name in ("alpha_split_pr", "alpha_split_trt"):
            value = getattr(self, name)
            if not (isinstance(value, Real) and 0 < value < 1):
                raise ValueError(f"{name} must lie strictly between 0 and 1, got {value!r}")
        for name in ("beta_split_pr", "beta_split_trt"):
            value = getattr(self, name)
            if not (isinstance(value, Real) and math.isfinite(value) and value >= 0):
                raise ValueError(f"{name} must be a finite nonnegative scalar, got {value!r}")
        for name in ("split_calendar_mu", "split_calendar_trt", "split_exposure_trt",
                     "random_intercept", "gp_constant_mean", "sample_beta",
                     "standardize", "sur", "use_inter_ensemble_move",
                     "use_trend_block", "trend_period_effects",
                     "trend_horseshoe", "anchor_moderator"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a boolean")
        for name in ("gamma_prior_a", "gamma_prior_b", "sigma_prior_a", "sigma_prior_b"):
            value = getattr(self, name)
            if (isinstance(value, bool) or not isinstance(value, Real)
                    or not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name} must be a finite positive scalar, got {value!r}")
        if self.outcome not in ("continuous", "binary", "ordinal"):
            raise ValueError(
                f"outcome must be 'continuous', 'binary', or 'ordinal', got {self.outcome!r}"
            )
        if self.outcome == "ordinal":
            if (isinstance(self.num_categories, bool)
                    or not isinstance(self.num_categories, Integral)
                    or self.num_categories < 2):
                raise ValueError("num_categories must be an integer >=2 for ordinal outcomes")
            object.__setattr__(self, "num_categories", int(self.num_categories))
        elif self.num_categories is not None:
            raise ValueError("num_categories must be None for nonordinal configurations")
        if self.baseline_cells not in ("all", "untreated"):
            raise ValueError(
                f"baseline_cells must be 'all' or 'untreated', got {self.baseline_cells!r}"
            )
        if self.baseline_cells == "untreated":
            if self.outcome != "continuous":
                raise ValueError(
                    "baseline_cells='untreated' is implemented for continuous outcomes "
                    f"only, got outcome={self.outcome!r}"
                )
            if self.tempering_levels > 1:
                raise ValueError(
                    "baseline_cells='untreated' cannot be combined with parallel "
                    "tempering (tempering_levels > 1)"
                )
        if (isinstance(self.cutpoint_prior_scale, bool)
                or not isinstance(self.cutpoint_prior_scale, Real)
                or not math.isfinite(self.cutpoint_prior_scale)
                or self.cutpoint_prior_scale <= 0):
            raise ValueError("cutpoint_prior_scale must be finite and positive")
        object.__setattr__(self, "cutpoint_prior_scale", float(self.cutpoint_prior_scale))
        if self.device not in ("auto", "cpu", "gpu"):
            raise ValueError(
                f"device must be 'auto', 'cpu' or 'gpu', got {self.device!r}"
            )
        if self.kernel_type not in ("se", "matern", "matern32", "matern52", "ar1"):
            raise ValueError(f"unknown kernel_type {self.kernel_type!r}")
        for name in ("lambda_knl", "sig_knl", "sigma_m"):
            value = getattr(self, name)
            if not (isinstance(value, Real) and math.isfinite(value) and value > 0):
                raise ValueError(f"{name} must be positive, got {value!r}")
        if (isinstance(self.sur_prior_var, bool) or not isinstance(self.sur_prior_var, Real)
                or not math.isfinite(self.sur_prior_var) or self.sur_prior_var < 0.0):
            raise ValueError(
                f"sur_prior_var must be a finite nonnegative scalar, got {self.sur_prior_var!r}"
            )
        for name in ("linear_trt", "learn_beta_scale", "ar1_errors", "interweave_gamma"):
            if getattr(self, name) is not None and not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a boolean or None (the model's default)")
        for name in ("linear_trt_scale", "beta_scale_prior"):
            value = getattr(self, name)
            if (isinstance(value, bool) or not isinstance(value, Real)
                    or not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name} must be a finite positive scalar, got {value!r}")
            object.__setattr__(self, name, float(value))
        if (isinstance(self.ar1_phi_max, bool) or not isinstance(self.ar1_phi_max, Real)
                or not 0 < self.ar1_phi_max < 1):
            raise ValueError("ar1_phi_max must lie strictly between 0 and 1")
        object.__setattr__(self, "ar1_phi_max", float(self.ar1_phi_max))
        if self.learn_beta_scale and not self.sample_beta:
            raise ValueError("learn_beta_scale needs sample_beta=True")
        if self.anchor_moderator and self.outcome != "continuous":
            raise ValueError("anchor_moderator=True supports a continuous outcome only")
        if self.ar1_errors and self.outcome != "continuous":
            raise ValueError("ar1_errors supports a continuous outcome only")
        if self.ar1_errors and not (self.random_intercept and not self.split_calendar_mu
                                    and not self.split_calendar_trt and not self.split_exposure_trt):
            raise ValueError(
                "ar1_errors needs random_intercept=True and unit-level partitions: no calendar "
                "or exposure splits in either forest (split_calendar_mu, split_calendar_trt, "
                "split_exposure_trt all False)")
        if self.ar1_errors and self.baseline_cells == "untreated":
            raise ValueError(
                "baseline_cells='untreated' cannot be combined with ar1_errors=True "
                "(the AR(1) sweep fits the baseline on all cells)")
        # --- distributed-lag kernel -------------------------------------------
        if self.lag_kernel_type not in ("se", "matern", "matern32", "matern52", "ar1"):
            raise ValueError(f"unknown lag_kernel_type {self.lag_kernel_type!r}")
        for name in ("lag_sig_knl", "lag_lambda_knl", "lag_sigma_m", "lag_scale_prior",
                     "extra_cols_prior_scale"):
            value = getattr(self, name)
            if (isinstance(value, bool) or not isinstance(value, Real)
                    or not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name} must be a finite positive scalar, got {value!r}")
            object.__setattr__(self, name, float(value))
        for name in ("lag_gp_constant_mean", "lag_learn_scale"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a boolean")

    @property
    def split_time_ps(self) -> bool:
        """Deprecated alias for :attr:`split_calendar_mu`."""
        warnings.warn(
            "'split_time_ps' is deprecated and will be removed in a future release; "
            "use 'split_calendar_mu' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.split_calendar_mu

    @property
    def sur_active(self) -> bool:
        """Whether SUR coupling is active (sur is True and sur_prior_var > 0)."""
        return bool(self.sur and self.sur_prior_var > 0.0)

    def total_iterations(self) -> int:
        """Total Gibbs sweeps to run (burn-in plus saved draws with thinning);
        an unset ``n_skip`` counts as the plain model's 1."""
        return self.num_burnin + self.num_sweeps * (1 if self.n_skip is None else self.n_skip)

    @property
    def unresolved(self) -> tuple[str, ...]:
        """Options still left to the fitted model's default (``None``)."""
        return tuple(name for name in MODEL_DEFAULTS if getattr(self, name) is None)

    def to_dict(self) -> dict[str, Any]:
        """Convert configuration to a plain dictionary."""
        d = dataclasses.asdict(self)
        d["split_time_ps"] = self.split_calendar_mu
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> LongBetConfig:
        """Construct a configuration from a dictionary, ignoring unknown keys.

        A key the dictionary lacks takes the dataclass default, except the
        encouragement-model extensions: a configuration saved before they
        existed describes the published model, so they take
        :data:`MODEL_DEFAULTS` rather than being left to a model's default."""
        d = dict(d)
        if "split_time_ps" in d and "split_calendar_mu" not in d:
            d["split_calendar_mu"] = d.pop("split_time_ps")
        for name in _EXTENSIONS:
            d.setdefault(name, MODEL_DEFAULTS[name])
        field_names = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in field_names})


#: What an unset option means for a plain ``LongBet`` or ``LongBetMulti`` fit:
#: the published model (Wang, Martinez and Hahn, 2024), without the
#: encouragement-model extensions and without thinning.
MODEL_DEFAULTS: dict[str, Any] = {
    "n_skip": 1, "linear_trt": False, "learn_beta_scale": False, "ar1_errors": False,
    "interweave_gamma": False,
}
_EXTENSIONS = ("linear_trt", "learn_beta_scale", "ar1_errors", "interweave_gamma")


def resolve_model_defaults(config: LongBetConfig,
                           defaults: dict[str, Any] | None = None) -> LongBetConfig:
    """``config`` with every unset option set from ``defaults`` (by default
    :data:`MODEL_DEFAULTS`); options set explicitly are kept."""
    defaults = MODEL_DEFAULTS if defaults is None else defaults
    unset = {name: defaults[name] for name in config.unresolved}
    return dataclasses.replace(config, **unset) if unset else config


_orig_longbet_config_init = LongBetConfig.__init__


#: Options removed after the encourage-v2 study, with the reason; passing the
#: removed behaviour raises, passing the surviving one is accepted and ignored.
_REMOVED = {
    "linear_trt_horseshoe": (True, "the linear moderator term always uses the horseshoe prior; a "
                             "single global scale shrank a strong moderator together with the null ones"),
    "linear_trt_outcome_only": (False, "the linear moderator term applies to both equations of "
                                "LongBetEncourage; it improved the adoption hazard, not only the outcome"),
    "adoption_frailty": (False, "LongBetEncourage's adoption hazard is population-averaged; with one "
                         "adoption event per unit a frailty's variance is weakly identified, mixes "
                         "poorly and flattens the offer's effect on adoption"),
}


def _longbet_config_init(self: LongBetConfig, *args: Any, **kwargs: Any) -> None:
    for name, (kept, reason) in _REMOVED.items():
        if name in kwargs:
            value = kwargs.pop(name)
            if bool(value) != kept:
                raise ValueError(f"{name}={value!r} is no longer supported: {reason}.")
            warnings.warn(f"'{name}' is no longer needed ({reason}); ignoring it.",
                          DeprecationWarning, stacklevel=2)
    if "split_time_ps" in kwargs:
        warnings.warn(
            "'split_time_ps' is deprecated and will be removed in a future release; "
            "use 'split_calendar_mu' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        val = kwargs.pop("split_time_ps")
        if "split_calendar_mu" not in kwargs:
            kwargs["split_calendar_mu"] = val
    _orig_longbet_config_init(self, *args, **kwargs)


LongBetConfig.__init__ = _longbet_config_init  # type: ignore[method-assign]
