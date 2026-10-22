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

#' Build leak-free counterfactual baseline forecasts, trajectory summaries, and a selection score
#'
#' Computes, strictly from pre-treatment data (every unit's context ends at the
#' common anchor `anchor_t = min G_i`, the first period in which any unit is
#' treated), a counterfactual baseline forecast from the built-in guarded-core
#' spectral forecaster with cross-series patch-ridge refinement (the backbone of
#' the companion paper), refined by a K-fold out-of-fold ridge adapter fitted on
#' never-treated units; unit-constant trajectory summaries; a cross-fitted
#' selection score; and pre-anchor residual factor loadings. The gated
#' `baseline_forecast` is `NULL` when the pre-treatment context has fewer than
#' 12 periods or `beyond_poly_r2 < 0.05`; the ungated `point_forecast` is always
#' returned and is what [fit_rank_k_longbet()] uses. [fit_rank_k_longbet()]
#' draws its posterior forecast panels on demand from the returned object, so
#' forecast draws are optional here (`n_forecast_draws = 0L` by default).
#'
#' @param y Outcome panel, `[N x T]`.
#' @param z Binary absorbing treatment indicator, `[N x T]`.
#' @param x Optional baseline covariates, `[N x P]`.
#' @param y_pre_history Optional strictly pre-treatment history panel, `[N x T_hist]`.
#' @param forecaster `NULL` (default; the built-in patch-ridge forecaster),
#'   `"spectral"` (the spectral backbone alone), or a Python callable
#'   `fn(context, horizon)` returning an `[N x horizon]` array.
#' @param n_folds Number of cross-fitting folds.
#' @param seed Random seed for deterministic fold assignment.
#' @param late_context Late-context structure of the built-in forecaster:
#'   `"fixed"` (default) or `"none"`.
#' @param n_forecast_draws Number `M` of forecast draws returned in
#'   `point_forecast_draws` and `baseline_forecast_draws`: `0L` (default; none)
#'   or at least `2L`. Not needed for [fit_rank_k_longbet()]: only its
#'   propagated variant uses forecast draws, and it reuses these only when
#'   `draw_method = "posterior"` and `M` is at least the number of draws it
#'   uses (`max_forecast_draws`, by default `n_forecast_draws`).
#' @param draw_method `"bootstrap"` (default; moving-block bootstrap whose
#'   cell-level channel is inflated by a control-error factor `kappa` in
#'   `[1, 3]`) or `"posterior"` (conjugate Normal-inverse-gamma posterior draws
#'   without `kappa` inflation; the law [fit_rank_k_longbet()] uses). Only used
#'   when `n_forecast_draws > 0`.
#' @param forecast_draw_block_length Optional block length of the moving-block
#'   bootstrap (default `NULL`: `round(L_ctx^(1/3)) + 1`).
#' @param n_factor_loadings Number of pre-anchor residual factor loadings
#'   returned in `pre_anchor_loadings` (default `3L`).
#' @return A list with `baseline_forecast`, `point_forecast`, `point_forecast_draws`,
#'   `unit_features`, `ps`, `pre_anchor_loadings`, `anchor_t`, `beyond_poly_r2`,
#'   `gate_open`, the gate-aware accessors `recommended_forecast`,
#'   `recommended_unit_features`, `recommended_ps`, `augment_x`,
#'   `baseline_forecast_draws`, `forecast_draws`, and `.py_features` (the Python
#'   object, which does not survive `saveRDS()`).
#' @export
build_causal_features <- function(y, z, x = NULL,
                                  y_pre_history = NULL,
                                  forecaster = NULL,
                                  n_folds = 5L,
                                  seed = 0L,
                                  late_context = "fixed",
                                  n_forecast_draws = 0L,
                                  draw_method = "bootstrap",
                                  forecast_draw_block_length = NULL,
                                  n_factor_loadings = 3L) {
  y <- as.matrix(y)
  z <- as.matrix(z)
  if (!is.null(x)) x <- as.matrix(x)
  if (!is.null(y_pre_history)) y_pre_history <- as.matrix(y_pre_history)

  lb <- longbet_py()
  np <- reticulate::import("numpy", convert = TRUE)
  res <- lb$build_causal_features(
    y = .as_np_matrix(y, "y"),
    z = .as_np_matrix(z, "z"),
    x = if (is.null(x)) NULL else .as_np_matrix(x, "x"),
    y_pre_history = if (is.null(y_pre_history)) NULL else .as_np_matrix(y_pre_history, "y_pre_history"),
    forecaster = forecaster,
    n_folds = as.integer(n_folds),
    seed = as.integer(seed),
    late_context = late_context,
    n_forecast_draws = as.integer(n_forecast_draws),
    forecast_draw_block_length = if (is.null(forecast_draw_block_length)) NULL else as.integer(forecast_draw_block_length),
    draw_method = as.character(draw_method),
    n_factor_loadings = as.integer(n_factor_loadings)
  )

  bf <- if (is.null(res$baseline_forecast)) NULL else as.matrix(np$asarray(res$baseline_forecast))
  pf <- as.matrix(np$asarray(res$point_forecast))
  uf <- as.matrix(np$asarray(res$unit_features))
  ps <- as.numeric(np$asarray(res$ps))
  r2 <- as.numeric(res$beyond_poly_r2)
  gate_open <- !is.null(bf)
  draws <- if (is.null(res$baseline_forecast_draws)) NULL else np$asarray(res$baseline_forecast_draws)
  pdraws <- if (is.null(res$point_forecast_draws)) NULL else np$asarray(res$point_forecast_draws)
  loadings <- if (is.null(res$pre_anchor_loadings)) NULL else as.matrix(np$asarray(res$pre_anchor_loadings))

  list(
    baseline_forecast = bf,
    point_forecast = pf,
    point_forecast_draws = pdraws,
    unit_features = uf,
    ps = ps,
    pre_anchor_loadings = loadings,
    anchor_t = as.integer(res$anchor_t),
    beyond_poly_r2 = r2,
    gate_open = gate_open,
    recommended_forecast = bf,
    recommended_unit_features = if (gate_open) uf else NULL,
    recommended_ps = if (gate_open) ps else NULL,
    augment_x = function(x_mat) {
      x_mat <- as.matrix(x_mat)
      if (gate_open) cbind(x_mat, uf) else x_mat
    },
    baseline_forecast_draws = draws,
    forecast_draws = res$forecast_draws,
    .py_features = res
  )
}

#' Deprecated: use [build_causal_features()]
#'
#' Kept for backward compatibility; the pretrained-TimesFM loader was removed in longbet-jax 0.3.0.
#' @inheritParams build_causal_features
#' @param ... Further arguments passed to [build_causal_features()].
#' @return See [build_causal_features()].
#' @export
build_causal_timesfm_features <- function(y, z, x = NULL, y_pre_history = NULL, forecaster = NULL,
                                          n_folds = 5L, seed = 0L, late_context = "fixed",
                                          n_forecast_draws = 0L, ...) {
  .Deprecated("build_causal_features", package = "longbet")
  build_causal_features(y = y, z = z, x = x, y_pre_history = y_pre_history, forecaster = forecaster,
                        n_folds = n_folds, seed = seed, late_context = late_context,
                        n_forecast_draws = n_forecast_draws, ...)
}

# Builds the Python RankKLongBetConfig. The defaults are those of
# RankKLongBetConfig() in Python and must stay identical to them: the Python
# fit treats a field that differs from its default as an explicit choice and
# warns when such a choice is ignored. As in Python, the variant-dependent
# arguments (`engine`, `uncertainty_propagation`, `max_forecast_draws`,
# `num_chains`, `num_burnin` and the three tree counts) default to NULL, which
# is passed on as Python None: RankKLongBetConfig then derives the engine from
# `uncertainty_propagation` (the propagated variant unless it is FALSE or
# engine = "stochtree") and fills in the budget validated for that variant.
# Every Rank-K front door in this file forwards its configuration arguments here.
.rank_k_config <- function(rank_k = 2L,
                           engine = NULL,
                           use_forecast = TRUE,
                           uncertainty_propagation = NULL,
                           n_forecast_draws = 8L,
                           max_forecast_draws = NULL,
                           use_spatial_gp = TRUE,
                           use_factor_loadings = TRUE,
                           n_factor_loadings = 3L,
                           num_chains = NULL,
                           num_burnin = NULL,
                           draws_per_component = 25L,
                           num_mcmc = 600L,
                           num_trees_mu = NULL,
                           num_trees_tau_mode1 = NULL,
                           num_trees_tau_higher = NULL,
                           random_seed = 2026L,
                           num_threads = 1L) {
  longbet_py()$RankKLongBetConfig(
    rank_k = .rank_k_count(rank_k, "rank_k"),
    engine = .rank_k_or_null(engine, as.character),
    use_forecast = .rank_k_flag(use_forecast, "use_forecast"),
    uncertainty_propagation = .rank_k_flag(uncertainty_propagation, "uncertainty_propagation",
                                           allow_null = TRUE),
    n_forecast_draws = .rank_k_count(n_forecast_draws, "n_forecast_draws"),
    max_forecast_draws = .rank_k_count(max_forecast_draws, "max_forecast_draws", allow_null = TRUE),
    use_spatial_gp = .rank_k_flag(use_spatial_gp, "use_spatial_gp"),
    use_factor_loadings = .rank_k_flag(use_factor_loadings, "use_factor_loadings"),
    n_factor_loadings = .rank_k_count(n_factor_loadings, "n_factor_loadings"),
    num_chains = .rank_k_count(num_chains, "num_chains", allow_null = TRUE),
    num_burnin = .rank_k_count(num_burnin, "num_burnin", allow_null = TRUE),
    draws_per_component = .rank_k_count(draws_per_component, "draws_per_component"),
    num_mcmc = .rank_k_count(num_mcmc, "num_mcmc"),
    num_trees_mu = .rank_k_count(num_trees_mu, "num_trees_mu", allow_null = TRUE),
    num_trees_tau_mode1 = .rank_k_count(num_trees_tau_mode1, "num_trees_tau_mode1", allow_null = TRUE),
    num_trees_tau_higher = .rank_k_count(num_trees_tau_higher, "num_trees_tau_higher", allow_null = TRUE),
    random_seed = .rank_k_count(random_seed, "random_seed"),
    num_threads = .rank_k_count(num_threads, "num_threads")
  )
}

# NULL (passed on as Python None: the selected variant's default) or the converted value.
.rank_k_or_null <- function(value, convert) {
  if (is.null(value)) NULL else convert(value)
}

# A single whole number as an integer (or NULL when allowed). Checked here
# because as.integer() would silently truncate 2.7 to 2 and pass NA on; the
# ranges (for example rank_k >= 1) are checked by RankKLongBetConfig.
.rank_k_count <- function(value, name, allow_null = FALSE) {
  if (allow_null && is.null(value)) return(NULL)
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value != round(value) ||
      abs(value) > .Machine$integer.max) {
    allowed <- if (allow_null) "NULL or a single whole number" else "a single whole number"
    stop(sprintf("`%s` must be %s.", name, allowed), call. = FALSE)
  }
  as.integer(value)
}

# A single TRUE or FALSE (or NULL when allowed). Checked here because isTRUE()
# would silently turn NA or a vector into FALSE.
.rank_k_flag <- function(value, name, allow_null = FALSE) {
  if (allow_null && is.null(value)) return(NULL)
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    allowed <- if (allow_null) "NULL, TRUE or FALSE" else "TRUE or FALSE"
    stop(sprintf("`%s` must be %s.", name, allowed), call. = FALSE)
  }
  value
}

# Converts the optional cohorts / t_0_panel / cluster_ids arguments shared by
# fit_rank_k_longbet() and rank_k_placebo_audit(). Numbers are passed as doubles
# so that the Python validator rejects non-integer values instead of R
# truncating them.
.rank_k_cohorts <- function(cohorts) {
  if (is.null(cohorts)) return(NULL)
  if (anyNA(cohorts)) {
    stop("`cohorts` contains NA; use ncol(y) for never-treated units.", call. = FALSE)
  }
  .as_np_vector(as.numeric(cohorts))
}

.rank_k_t0 <- function(t_0_panel) {
  if (is.null(t_0_panel)) return(NULL)
  if (length(t_0_panel) != 1L || is.na(t_0_panel)) {
    stop("`t_0_panel` must be a single integer: the number of leading pre-period columns of `y`.",
         call. = FALSE)
  }
  as.numeric(t_0_panel)
}

.rank_k_clusters <- function(cluster_ids) {
  if (is.null(cluster_ids)) return(NULL)
  if (anyNA(cluster_ids)) {
    stop("`cluster_ids` contains NA; give every unit a cluster label.", call. = FALSE)
  }
  .as_np_vector(as.character(cluster_ids))
}

#' Forecast-Augmented Rank-K Exposure-Basis BCF for longitudinal experiments
#'
#' One estimator with two uncertainty variants, selected by
#' `uncertainty_propagation` (the engine follows from it). Post-adoption
#' outcomes are projected within unit onto an orthonormal exposure basis whose
#' first mode is the constant (so mode 1 is each treated unit's time-averaged
#' effect, and the ATT depends on mode 1 only), and one cross-sectional BCF is
#' fitted per mode, with the counterfactual forecast as a control variate whose
#' slope is learned from never-treated units.
#'
#' * "propagated" (the default: `uncertainty_propagation = TRUE`, engine
#'   `"jax"`): the prognostic and treatment forests, the forecast
#'   control-variate slope `theta_k`, the pre-anchor factor shocks and the
#'   Matern-3/2 spatial field (or cluster random effect, kriged from
#'   never-treated to treated units) are learned from never-treated units and
#'   sampled together, and `n_forecast_draws` posterior forecast draws are
#'   pooled through mode 1 into a modular (cut) posterior. Pass `coords` or
#'   `cluster_ids` whenever units can share regional shocks. The sampler is an
#'   approximate partially collapsed Gibbs sampler (working likelihoods in the
#'   forest steps), not an exact posterior sampler.
#' * "plug-in" (`uncertainty_propagation = FALSE`, engine `"stochtree"`): one
#'   `stochtree` BCF per mode on a plug-in OLS-ridge forecast slope, without
#'   forecast draws, factor loadings or spatial field. Use it as a fast plug-in
#'   approximation for unit-level screening and targeting under unit-level
#'   randomization without spatial clustering. `stochtree` is installed with
#'   the Python package `longbet-jax` (a core dependency since 0.5.0).
#'
#' In the companion paper's RCT study the propagated default passed the
#' pre-registered ATT coverage rule on every development panel (93.8% on
#' temperature, 90.6% on county employment) and every hold-out set of the
#' second frozen confirmation (pooled 95% coverage 94.4%, SD ratios 0.93-1.02;
#' rules H1'-H7 passed), whereas the fast plug-in variant missed the pooled
#' coverage rule on one Stage-1 hold-out set (89.8% against 90%) and
#' under-covered in synthetic panels with spatially clustered shocks,
#' loading-correlated assignment or uncertain forecast slopes, while achieving
#' slightly sharper unit-level point recovery (unit PEHE +0.010 lower on
#' covariate-stratified cells; unit sign accuracy 83.4% vs 71.7% on temperature
#' and 70.0% vs 61.0% on county employment in the sign-flip regime).
#'
#' The defaults mirror `RankKLongBetConfig()` in Python. The variant-dependent
#' arguments (`engine`, `uncertainty_propagation`, `max_forecast_draws`,
#' `num_chains`, `num_burnin` and the three tree counts) default to `NULL` and
#' are resolved in Python to the budget validated for the selected variant:
#'
#' * propagated (default): engine `"jax"`, `num_chains = 2`, `num_burnin = 100`,
#'   `n_forecast_draws = 8` without adaptive expansion and
#'   `draws_per_component = 25` (2 x 8 x 25 = 400 posterior draws), 50 / 30 /
#'   15 trees, spatial field and 3 factor loadings. For small panels use the
#'   paper's application budget, `num_chains = 4L, num_burnin = 300L,
#'   draws_per_component = 100L`;
#' * plug-in: engine `"stochtree"`, one chain per mode (`num_chains = 1`), 10
#'   grow-from-root sweeps, `num_burnin = 200`, `num_mcmc = 600` retained draws
#'   per mode, and 75 / 35 / 20 trees for `mu_k` / `tau_1` / `tau_k`.
#'
#' In both variants the 95% ATT interval (`att_ci_low`, `att_ci_high`) is the
#' 2.5% / 97.5% posterior percentiles of the ATT draws. Unit-level (CATE)
#' intervals are approximate and below nominal in the paper's benchmarks (95%
#' unit-level intervals covered 70-92% with the propagated and 88-93% with the
#' plug-in variant), and the cell-level intervals `catt_ci_low` /
#' `catt_ci_high` (2.5% / 97.5% posterior percentiles) are approximate as well;
#' treat both as descriptive.
#'
#' Invalid inputs raise an error before any computation, as do the
#' contradictory combinations `engine = "jax"` with
#' `uncertainty_propagation = FALSE` and `engine = "stochtree"` with
#' `uncertainty_propagation = TRUE`; counts must be single whole numbers and
#' flags single `TRUE` / `FALSE` values. Explicit inputs that the fit would
#' ignore (for example `coords` with the plug-in variant), or that contradict
#' each other, trigger a Python `RankKConfigWarning` (issued through Python's
#' `warnings` module, not as an R condition); choosing a supported variant or
#' option never does. The arguments added in 0.5.0 (`pre_anchor_loadings`,
#' `cohorts`, `max_forecast_draws` and `draws_per_component`) come after
#' `num_threads`, so 0.4.0 positional calls keep their meaning; pass them by
#' name.
#'
#' @param y Outcome panel, `[N x T]`.
#' @param z Binary absorbing treatment indicator, `[N x T]`.
#' @param x Pre-treatment baseline covariates, `[N x P]`.
#' @param features The list returned by [build_causal_features()] for the same
#'   `y` and `z` (recommended). Its Python object supplies `unit_features`, `ps`,
#'   the forecast (`point_forecast`), `pre_anchor_loadings` and the posterior
#'   forecast draws, drawn on demand. A list without the Python object (for
#'   example one restored with `readRDS()`) supplies its stored arrays instead.
#' @param unit_features,ps,baseline_forecast,baseline_forecast_draws,pre_anchor_loadings
#'   Optional explicit inputs that take precedence over `features`: extra unit
#'   covariates `[N x Q]`, a selection score in `[0, 1]` (length `N`), the
#'   counterfactual forecast `[N x T]`, posterior forecast draws `[M x N x T]`,
#'   and pre-anchor factor loadings `[N x q]`. None of them may depend on
#'   treated outcomes. The forecast draws and the loadings are read by the
#'   propagated variant only (the plug-in variant ignores them with a warning);
#'   the draws are used only when `M` is at least `max_forecast_draws` (by
#'   default `n_forecast_draws`).
#' @param coords Optional `[N x 2]` matrix of `(latitude, longitude)` in degrees
#'   for the Matern-3/2 spatial field on great-circle distances (propagated
#'   variant; the plug-in variant ignores it with a warning). Without `coords`
#'   or `cluster_ids` the propagated variant fits no spatial field; pass one of
#'   them with it whenever units can share regional shocks.
#' @param cluster_ids Optional length-`N` vector of cluster labels for a cluster
#'   random effect (propagated variant); used only when `coords` is `NULL`.
#' @param cohorts Optional length-`N` integer vector: the number of periods before
#'   each unit's first treated period (the zero-based index of its first treated
#'   column), `ncol(y)` for never-treated units. Must agree with `z`; inferred
#'   from `z` by default.
#' @param t_0_panel Optional single integer: the number of leading columns of `y`
#'   that form the pre-period, between 1 and the first adoption period (default:
#'   the first adoption period).
#' @param rank_k Number of exposure modes (default `2L`: level + concave ramp;
#'   only 1 to 3 were validated).
#' @param engine `NULL` (default: derived from `uncertainty_propagation`),
#'   `"jax"` (the propagated variant) or `"stochtree"` (the plug-in variant).
#'   Kept for backward compatibility: `engine = "stochtree"` alone selects the
#'   plug-in variant, whereas `engine = "jax"` with
#'   `uncertainty_propagation = FALSE` and `engine = "stochtree"` with
#'   `uncertainty_propagation = TRUE` are errors.
#' @param use_forecast Use the forecast as the control variate and as a
#'   prognostic covariate (default `TRUE`).
#' @param uncertainty_propagation Selects the variant: `TRUE` (the default) the
#'   propagated variant, which pools `n_forecast_draws` posterior forecast draws
#'   into a modular (cut) posterior, `FALSE` the fast plug-in variant. `NULL`
#'   (default) means `FALSE` if `engine = "stochtree"` and `TRUE` otherwise.
#' @param n_forecast_draws Number `M` of posterior forecast draws of the
#'   propagated variant (default `8L`, as validated). `0L` or `1L` propagate
#'   no forecast uncertainty (with `uncertainty_propagation = TRUE` this
#'   triggers a warning).
#' @param max_forecast_draws Upper limit of the adaptive expansion (default
#'   `NULL`: `n_forecast_draws`, no expansion, as validated). Larger values
#'   enable the expansion, which was not validated.
#' @param use_spatial_gp Include the spatial field when `coords` or `cluster_ids`
#'   are supplied (default `TRUE`; propagated variant).
#' @param use_factor_loadings Include the pre-anchor factor term (default `TRUE`;
#'   propagated variant).
#' @param n_factor_loadings Number `q` of factor loadings extracted when neither
#'   `features` nor `pre_anchor_loadings` supplies them (default `3L`;
#'   propagated variant).
#' @param num_chains Chains per forecast draw of the propagated variant
#'   (default `NULL`: `2L` propagated; `1L` plug-in, whose `stochtree` sampler
#'   runs one chain per mode).
#' @param num_burnin Burn-in sweeps (default `NULL`: `200L` plug-in, `100L`
#'   propagated).
#' @param draws_per_component Retained draws per chain and forecast draw of the
#'   propagated variant (default `25L`): its posterior has
#'   `num_chains * n_forecast_draws * draws_per_component` draws (400 by
#'   default).
#' @param num_mcmc Retained draws (default `600L`): per mode in the plug-in
#'   variant; in total, split evenly across chains, in the propagated variant
#'   when it propagates no forecast draws, and ignored when it does.
#' @param num_trees_mu,num_trees_tau_mode1,num_trees_tau_higher Trees of each
#'   prognostic forest, of the mode-1 treatment forest and of each higher-mode
#'   treatment forest (default `NULL`: 75 / 35 / 20 plug-in, 50 / 30 / 15
#'   propagated).
#' @param random_seed Base random seed (default `2026L`); mode 1 is seeded
#'   identically for every `rank_k`.
#' @param num_threads Threads of the `stochtree` sampler (default `1L`; plug-in
#'   variant).
#' @return A list with the ATT summaries `att_hat`, `att_se`, `att_ci_low`,
#'   `att_ci_high` and `att_draws`; the cell effects `catt_hat`, `catt_ci_low`,
#'   `catt_ci_high` (`[N x T]`) and `catt_draws` (`[S x N x T]`);
#'   `mode_tau_draws`, a list of `[N x S]` matrices, one per mode (row `i` of the
#'   first is treated unit `i`'s time-averaged effect); the forecast slopes
#'   `theta_f_modes` and `theta_f_draws`; `sigma_modes`; the within- and
#'   between-forecast-draw ATT variances `att_within_var`, `att_between_var` and
#'   `att_between_share`; the spatial-field scale and range per mode,
#'   `s_zeta_modes` and `ell_modes` (meaningful only when a spatial field is
#'   fitted); `n_components`; the convergence diagnostics `rhat_att` (`NaN`
#'   for the plug-in variant, which runs one chain) and `ess_att`; and
#'   `runtime_sec`.
#' @export
fit_rank_k_longbet <- function(y, z, x,
                               features = NULL,
                               unit_features = NULL,
                               ps = NULL,
                               baseline_forecast = NULL,
                               baseline_forecast_draws = NULL,
                               coords = NULL,
                               cluster_ids = NULL,
                               t_0_panel = NULL,
                               rank_k = 2L,
                               engine = NULL,
                               use_forecast = TRUE,
                               uncertainty_propagation = NULL,
                               n_forecast_draws = 8L,
                               use_spatial_gp = TRUE,
                               use_factor_loadings = TRUE,
                               n_factor_loadings = 3L,
                               num_chains = NULL,
                               num_burnin = NULL,
                               num_mcmc = 600L,
                               num_trees_mu = NULL,
                               num_trees_tau_mode1 = NULL,
                               num_trees_tau_higher = NULL,
                               random_seed = 2026L,
                               num_threads = 1L,
                               # New in 0.5.0, appended so that 0.4.0 positional calls keep their meaning.
                               pre_anchor_loadings = NULL,
                               cohorts = NULL,
                               max_forecast_draws = NULL,
                               draws_per_component = 25L) {
  y <- as.matrix(y)
  z <- as.matrix(z)
  x <- as.matrix(x)
  cfg <- .rank_k_config(
    rank_k = rank_k,
    engine = engine,
    use_forecast = use_forecast,
    uncertainty_propagation = uncertainty_propagation,
    n_forecast_draws = n_forecast_draws,
    max_forecast_draws = max_forecast_draws,
    use_spatial_gp = use_spatial_gp,
    use_factor_loadings = use_factor_loadings,
    n_factor_loadings = n_factor_loadings,
    num_chains = num_chains,
    num_burnin = num_burnin,
    draws_per_component = draws_per_component,
    num_mcmc = num_mcmc,
    num_trees_mu = num_trees_mu,
    num_trees_tau_mode1 = num_trees_tau_mode1,
    num_trees_tau_higher = num_trees_tau_higher,
    random_seed = random_seed,
    num_threads = num_threads
  )

  py_feats <- if (!is.null(features)) features$.py_features else NULL
  if (!is.null(py_feats) && reticulate::py_is_null_xptr(py_feats)) py_feats <- NULL
  if (!is.null(features) && is.null(py_feats)) {
    # A features list without a live Python object (e.g. restored with readRDS()):
    # fill in only the stored arrays that this configuration uses, so that they
    # are not reported as ignored explicit inputs. The resolved configuration
    # gives the variant (the engine derived from `uncertainty_propagation`) and
    # the number of forecast draws (`max_forecast_draws = NULL` means
    # `n_forecast_draws`; fewer than 2 propagate none).
    propagated <- isTRUE(cfg$uncertainty_propagation)
    forecast_on <- isTRUE(cfg$use_forecast)
    m_max <- if (is.null(cfg$max_forecast_draws)) cfg$n_forecast_draws else cfg$max_forecast_draws
    if (is.null(unit_features)) unit_features <- features$unit_features
    if (is.null(ps)) ps <- features$ps
    if (is.null(baseline_forecast) && forecast_on) baseline_forecast <- features$point_forecast
    if (is.null(pre_anchor_loadings) && propagated && isTRUE(cfg$use_factor_loadings)) {
      pre_anchor_loadings <- features$pre_anchor_loadings
    }
    stored_draws <- features$point_forecast_draws
    if (is.null(baseline_forecast_draws) && !is.null(stored_draws) && propagated &&
        forecast_on && !is.null(baseline_forecast) && cfg$n_forecast_draws >= 2L &&
        identical(features$forecast_draws$draw_method, "posterior") &&
        dim(stored_draws)[1] >= m_max) {
      baseline_forecast_draws <- stored_draws
    }
  }

  lb <- longbet_py()
  np <- reticulate::import("numpy", convert = TRUE)
  res <- lb$fit_rank_k_longbet(
    y = .as_np_matrix(y, "y"),
    z = .as_np_matrix(z, "z"),
    x = .as_np_matrix(x, "x"),
    features = py_feats,
    unit_features = if (is.null(unit_features)) NULL else .as_np_matrix(unit_features, "unit_features"),
    ps = if (is.null(ps)) NULL else .as_np_vector(as.numeric(ps)),
    baseline_forecast = if (is.null(baseline_forecast)) NULL else .as_np_matrix(baseline_forecast, "baseline_forecast"),
    baseline_forecast_draws = if (is.null(baseline_forecast_draws)) NULL else np$asarray(baseline_forecast_draws),
    pre_anchor_loadings = if (is.null(pre_anchor_loadings)) NULL else .as_np_matrix(pre_anchor_loadings, "pre_anchor_loadings"),
    coords = if (is.null(coords)) NULL else .as_np_matrix(as.matrix(coords), "coords"),
    cluster_ids = .rank_k_clusters(cluster_ids),
    cohorts = .rank_k_cohorts(cohorts),
    t_0_panel = .rank_k_t0(t_0_panel),
    cfg = cfg
  )

  list(
    att_hat = as.numeric(res$att_hat),
    att_se = as.numeric(res$att_se),
    att_ci_low = as.numeric(res$att_ci_low),
    att_ci_high = as.numeric(res$att_ci_high),
    att_draws = as.numeric(np$asarray(res$att_draws)),
    catt_hat = as.matrix(np$asarray(res$catt_hat)),
    catt_ci_low = as.matrix(np$asarray(res$catt_ci_low)),
    catt_ci_high = as.matrix(np$asarray(res$catt_ci_high)),
    catt_draws = np$asarray(res$catt_draws),
    mode_tau_draws = lapply(res$mode_tau_draws, as.matrix),
    theta_f_modes = as.numeric(unlist(res$theta_f_modes)),
    theta_f_draws = lapply(res$theta_f_draws, as.numeric),
    sigma_modes = as.numeric(unlist(res$sigma_modes)),
    att_within_var = as.numeric(res$att_var_within),
    att_between_var = as.numeric(res$att_var_between),
    att_between_share = as.numeric(res$att_between_share),
    s_zeta_modes = as.numeric(unlist(res$s_zeta_modes)),
    ell_modes = as.numeric(unlist(res$ell_modes)),
    n_components = as.integer(res$n_components),
    rhat_att = as.numeric(res$rhat_att),
    ess_att = as.numeric(res$ess_att),
    runtime_sec = as.numeric(res$runtime_sec)
  )
}

#' Pre-adoption placebo audit for Rank-K Exposure-Basis BCF
#'
#' Holds out the final `pseudo_anchor_offset` pre-adoption periods (minimum 2) as
#' a placebo post-period, gives each ever-treated unit a pseudo-cohort that keeps
#' its adoption lag relative to the first wave, rebuilds the features on the
#' placebo panel with [build_causal_features()], its default forecaster and the
#' configured `n_factor_loadings`, and runs [fit_rank_k_longbet()] on these
#' strictly pre-adoption outcomes, where the true effect is zero. The audit
#' passes when the 95% placebo ATT interval covers zero. The placebo ATT depends
#' on mode 1 only, so `rank_k = 1L` gives the same placebo ATT at lower cost.
#'
#' @inheritParams fit_rank_k_longbet
#' @param y_pre_history Optional `[N x T_hist]` strictly pre-treatment history
#'   panel; needed when `T_hist` plus the first adoption period is shorter than
#'   `pseudo_anchor_offset + 2`.
#' @param features Ignored (accepted for symmetry with [fit_rank_k_longbet()]):
#'   the audit always rebuilds the features on the placebo panel. A non-`NULL`
#'   value triggers an R warning.
#' @param pseudo_anchor_offset A single whole number of pre-adoption periods held
#'   out for the placebo window (default `4L`; values below 2 are raised to 2).
#' @param ... Configuration arguments of [fit_rank_k_longbet()]
#'   (`uncertainty_propagation`, `rank_k`, `num_burnin`, ...), with the same
#'   defaults (the propagated variant); pass the configuration of the audited fit.
#' @return A list with `placebo_anchor_t`, `placebo_horizon`, `placebo_att_hat`,
#'   `placebo_att_se`, `placebo_ci_low`, `placebo_ci_high`, `placebo_z_stat`,
#'   `covers_zero_95`, `flag_parallel_trends_violation`, and `passes_audit`.
#' @export
rank_k_placebo_audit <- function(y, z, x,
                                 y_pre_history = NULL,
                                 features = NULL,
                                 coords = NULL,
                                 cluster_ids = NULL,
                                 cohorts = NULL,
                                 pseudo_anchor_offset = 4L,
                                 ...) {
  if ("placebo_horizon" %in% names(list(...))) {
    stop("`placebo_horizon` is not an argument; use `pseudo_anchor_offset`, the number of ",
         "pre-adoption periods held out.", call. = FALSE)
  }
  if (!is.null(features)) {
    warning("rank_k_placebo_audit() ignores `features`: it always rebuilds the features on the ",
            "placebo window with build_causal_features() and the default forecaster. Omit ",
            "`features` to silence this warning.", call. = FALSE)
  }
  y <- as.matrix(y)
  z <- as.matrix(z)
  x <- as.matrix(x)
  cfg <- .rank_k_config(...)
  lb <- longbet_py()
  res <- lb$rank_k_placebo_audit(
    y = .as_np_matrix(y, "y"),
    z = .as_np_matrix(z, "z"),
    x = .as_np_matrix(x, "x"),
    y_pre_history = if (is.null(y_pre_history)) NULL else .as_np_matrix(as.matrix(y_pre_history), "y_pre_history"),
    features = NULL,
    coords = if (is.null(coords)) NULL else .as_np_matrix(as.matrix(coords), "coords"),
    cluster_ids = .rank_k_clusters(cluster_ids),
    cohorts = .rank_k_cohorts(cohorts),
    pseudo_anchor_offset = .rank_k_count(pseudo_anchor_offset, "pseudo_anchor_offset"),
    cfg = cfg
  )
  list(
    placebo_anchor_t = as.integer(res$placebo_anchor_t),
    placebo_horizon = as.integer(res$placebo_horizon),
    placebo_att_hat = as.numeric(res$placebo_att_hat),
    placebo_att_se = as.numeric(res$placebo_att_se),
    placebo_ci_low = as.numeric(res$placebo_ci_low),
    placebo_ci_high = as.numeric(res$placebo_ci_high),
    placebo_z_stat = as.numeric(res$placebo_z_stat),
    covers_zero_95 = isTRUE(res$covers_zero_95),
    flag_parallel_trends_violation = isTRUE(res$flag_parallel_trends_violation),
    passes_audit = isTRUE(res$passes_audit)
  )
}
