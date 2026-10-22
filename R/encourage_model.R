# Plain R conversion is essential here: neither model nor prediction objects may
# retain reticulate handles across saveRDS or a new R process.
.encouragement_plain <- function(value) {
  if (inherits(value, "python.builtin.object") && reticulate::py_has_attr(value, "_asdict")) {
    value <- value$`_asdict`()
  }
  value <- .encouragement_as_r(value)
  if (is.data.frame(value)) {
    attr(value, "pandas.index") <- NULL
    rownames(value) <- NULL
    for (name in names(value)) value[[name]] <- .encouragement_plain(value[[name]])
  } else if (is.list(value)) {
    value <- lapply(value, .encouragement_plain)
  } else if (inherits(value, "python.builtin.object")) {
    stop("The engine returned an unsupported live Python value.", call. = FALSE)
  }
  value
}

.encouragement_load <- function(object) {
  if (!inherits(object, "longbet_encourage") || !is.raw(object$raw_archive)) {
    stop("object must be a fitted longbet_encourage with a raw archive.", call. = FALSE)
  }
  path <- tempfile(fileext = ".npz")
  on.exit(unlink(path), add = TRUE)
  writeBin(object$raw_archive, path)
  longbet_py()
  reticulate::import("longbet", convert = FALSE)$LongBetEncourage$load(path)
}

.encouragement_config_args <- function(config, lb) {
  if (!is.list(config) || (length(config) &&
      (is.null(names(config)) || anyNA(names(config)) ||
       any(!nzchar(names(config))) || anyDuplicated(names(config))))) {
    stop("config must be a list of uniquely named sampler options.", call. = FALSE)
  }
  defaults <- .encouragement_plain(reticulate::import("dataclasses")$asdict(lb$LongBetConfig()))
  if (any(!names(config) %in% names(defaults)) || "outcome" %in% names(config)) {
    stop("config names must be LongBetConfig options; set outcome separately.", call. = FALSE)
  }
  # R's ordinary numeric literals are doubles. Respect Python integer fields
  # without silently rounding noninteger values or maintaining a second config.
  # n_skip is unset (NULL) in LongBetConfig: the model supplies its default.
  for (name in names(config)) {
    if (is.integer(defaults[[name]]) || name %in% c("inner_loop_length", "n_skip")) {
      value <- config[[name]]
      if (is.null(value) && identical(name, "inner_loop_length")) next
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
          value != trunc(value) || abs(value) > .Machine$integer.max) {
        stop(sprintf("config$%s must be a finite integer.", name), call. = FALSE)
      }
      config[[name]] <- as.integer(value)
    }
  }
  config
}

#' The encouragement model: adoption clock or continuous treatment intensity
#'
#' With `treatment = "absorbing"` (default) this fits two LongBet equations and
#' composes them: the outcome on the observed adoption clock (periods since
#' adoption, both arms contributing adoption events) and adoption as a
#' population-averaged discrete-time hazard on the offer clock. The offer's
#' effects on the outcome and on adoption follow by composition. This uses the
#' exclusion restriction and assumes adoption timing is not confounded with the
#' outcome innovations given covariates and the outcome's unit intercept; the
#' design-based reference returned by `predict()` needs neither and is the
#' check on the population offer effects.
#'
#' With `treatment = "intensity"` the encouraged behaviour `d` is a continuous,
#' non-absorbing dose (hours, training-months). Stage 1 is a continuous LongBet
#' for the dose on the offer clock (unit intercept, AR(1) errors); Stage 2 is the
#' distributed-lag kernel for the outcome, `nu_Y(x) sum_l kappa_l D_{i,t-l}` with
#' a Gaussian-process prior over lags. `identification` chooses what Stage 2
#' learns the kernel from: `"deconvolution"` (T0: the Stage-1 fitted
#' offer-induced dose only; immune to within-arm dose endogeneity and
#' measurement error), `"control_function"` (T1, default: the observed dose plus
#' the three-term plug-in control function of the Stage-1 residuals) or
#' `"unconfounded"` (T2: the observed dose alone). Offer effects are composed
#' draw by draw, `CITT_Y(x,t) = nu_Y(x) sum_l kappa_l CITT_D(x,t-l)`.
#'
#' Both moderators carry a horseshoe linear term beside the trees, the
#' trajectory's scale is learned, and a continuous outcome has AR(1) errors;
#' these are part of the model. Draws are thinned by 2 unless `config` sets
#' `n_skip`.
#'
#' @param y Complete outcome matrix, `[N x T]`.
#' @param d,z Actual adoption and assigned encouragement, complete `[N x T]`
#'   matrices. `z` is binary and absorbing with a single wave and a permanent
#'   holdout; `d` is binary absorbing for `treatment = "absorbing"` and a finite
#'   continuous dose for `treatment = "intensity"`.
#' @param x Baseline covariates, `[N x P]`.
#' @param t Increasing calendar vector with whole-unit gaps; defaults to `1:T`.
#' @param x_trt Optional baseline treatment-forest covariates.
#' @param outcome `"continuous"` or `"binary"` for `y`; adoption is always a
#'   binary hazard. The intensity model supports a continuous outcome only.
#' @param treatment `"absorbing"` (adoption clock) or `"intensity"`
#'   (distributed-lag kernel for a continuous dose).
#' @param identification Identification tier of the intensity model:
#'   `"deconvolution"`, `"control_function"` (default when NULL) or
#'   `"unconfounded"`. Must be NULL for `treatment = "absorbing"`.
#' @param lag_horizon Number of lags beyond zero in the kernel; NULL means
#'   `T - t0`, one lag per post-offer period (intensity only).
#' @param demean_pre Demean the dose by each unit's pre-offer mean before
#'   building lags (intensity only).
#' @param stage2_backend `"longbet"` (the sampler) or `"linear"` (a closed-form
#'   two-way fixed-effects regression without covariate moderation, for smoke
#'   tests and fallbacks; intensity only).
#' @param strata Optional length-`N` baseline randomization strata (sites) used
#'   by the design-based reference and the dip screen of the intensity model.
#' @param config Named list of Python `LongBetConfig` sampler options applied to
#'   both equations. Set `num_chains`, `num_burnin`, `num_sweeps` and
#'   `random_seed` here; convergence still requires checking. The lag kernel's
#'   prior options (`lag_sig_knl`, `lag_lambda_knl`, ...) are set here too.
#' @param verbose Whether to print a fitting message.
#' @return A `longbet_encourage` object containing a versioned raw archive,
#'   metadata and configuration, with no live Python handles. The archive
#'   includes original outcomes, adoption, assignment and covariates so prediction
#'   replays the same units and target. It can be saved using `saveRDS()`.
#' @details Both ITTs average over all original study units with equal weights
#'   and each unit's own fitted intercepts. Binary contrasts are probability
#'   differences. The Wald ratio is not automatically CACE; exclusion,
#'   monotonicity, relevance and appropriate adoption-history assumptions are
#'   additional requirements. For the intensity model the summary ratio is the
#'   cumulative return per unit of offer-induced cumulative dose (`omega_cum`).
#' @export
longbet_encourage <- function(y, d, z, x, t = NULL, x_trt = NULL,
                             outcome = "continuous", treatment = "absorbing",
                             identification = NULL, lag_horizon = NULL, demean_pre = FALSE,
                             stage2_backend = "longbet", strata = NULL,
                             config = list(), verbose = FALSE) {
  outcome <- match.arg(outcome, c("continuous", "binary"))
  treatment <- match.arg(treatment, c("absorbing", "intensity"))
  stage2_backend <- match.arg(stage2_backend, c("longbet", "linear"))
  if (!is.null(identification)) {
    identification <- match.arg(identification, c("deconvolution", "control_function", "unconfounded"))
  }
  if (!is.null(lag_horizon)) {
    if (!is.numeric(lag_horizon) || length(lag_horizon) != 1L || !is.finite(lag_horizon) ||
        lag_horizon < 0 || lag_horizon != trunc(lag_horizon)) {
      stop("lag_horizon must be a nonnegative integer or NULL.", call. = FALSE)
    }
    lag_horizon <- as.integer(lag_horizon)
  }
  if (!is.logical(demean_pre) || length(demean_pre) != 1L || is.na(demean_pre)) {
    stop("demean_pre must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.factor(strata)) strata <- as.character(strata)
  lb <- longbet_py()
  options <- .encouragement_config_args(config, lb)
  model_args <- list(outcome = outcome, treatment = treatment)
  if (treatment == "intensity") {
    model_args <- c(model_args, list(identification = identification, lag_horizon = lag_horizon,
                                     demean_pre = demean_pre, stage2_backend = stage2_backend))
  } else if (!is.null(identification) || !is.null(lag_horizon) || isTRUE(demean_pre)) {
    stop("identification, lag_horizon and demean_pre apply to treatment = 'intensity' only.", call. = FALSE)
  }
  backend <- do.call(lb$LongBetEncourage, c(model_args, options))
  if (isTRUE(verbose)) {
    message(if (treatment == "intensity")
      "Fitting the dose on the offer clock and the outcome on the dose lags."
    else "Fitting the outcome on the adoption clock and the adoption hazard on the offer.")
  }
  backend$fit(y = as.matrix(y), d = as.matrix(d), z = as.matrix(z), x = as.matrix(x),
             t = if (is.null(t)) NULL else as.numeric(t),
             x_trt = if (is.null(x_trt)) NULL else as.matrix(x_trt),
             strata = if (is.null(strata)) NULL else as.character(strata))
  path <- tempfile(fileext = ".npz")
  on.exit(unlink(path), add = TRUE)
  backend$save(path)
  structure(list(
    raw_archive = readBin(path, what = "raw", n = file.info(path)$size),
    metadata = .encouragement_plain(backend$metadata),
    config = .encouragement_plain(reticulate::import("dataclasses")$asdict(backend$config)),
    treatment = treatment,
    call = match.call()
  ), class = "longbet_encourage")
}

#' Predict common-population encouragement effects
#'
#' Every study unit is evaluated under both schedules (offered from the launch
#' period onward, never offered) with its own fitted intercepts, so the target
#' is the finite study population that the design-based reference describes.
#'
#' @param object A fitted `longbet_encourage` object, including one read by
#'   `readRDS()` in another session.
#' @param summary_only Kept for API symmetry; the standardized draws are small
#'   and always retained.
#' @param groups Optional prespecified baseline group labels, length `N`.
#' @param weights Optional fixed nonnegative unit weights, normalized within each
#'   group; they must not depend on post-encouragement quantities.
#' @param block_size Optional positive integer number of units composed per block.
#' @param alpha Tail probability for pointwise intervals.
#' @param min_ess,max_rhat Diagnostic thresholds; passing is not a coverage claim.
#' @param practical_threshold Optional positive adoption risk difference for a
#'   separate posterior practical-relevance probability, not an IV strength cutoff.
#' @param ... Reserved; unsupported arguments raise an error.
#' @return A plain-R `longbet_encourage.pred` list with `effects`, `reference`,
#'   `wald`, `wald_unchecked`, `first_stage`, `comparison`, `diagnostics`, `draws`,
#'   `weights`, group information and metadata.
#'   Draw arrays retain axes `(group, horizon, chain, retained_draw)`.
#'   Wald summaries use medians, never a default ratio mean. `wald` withholds
#'   unconverged intervals; `wald_unchecked` allows inspection of those summaries.
#'   The `reference` contains separate normal Fieller/Anderson--Rubin-style
#'   confidence sets with infinite endpoints and unavailable cases preserved.
#' @export
predict.longbet_encourage <- function(object, summary_only = TRUE, groups = NULL,
                                      weights = NULL, block_size = NULL,
                                      alpha = 0.05, min_ess = 400, max_rhat = 1.01,
                                      practical_threshold = NULL, ...) {
  dots <- list(...)
  if (!is.null(dots$new_x)) {
    return(do.call(predict_conditional, c(list(object = object, alpha = alpha), dots)))
  }
  .reject_unsupported(...)
  if (!is.null(block_size)) {
    if (!is.numeric(block_size) || length(block_size) != 1L ||
        !is.finite(block_size) || block_size <= 0 || block_size != trunc(block_size) ||
        block_size > .Machine$integer.max) {
      stop("block_size must be a positive integer.", call. = FALSE)
    }
    block_size <- as.integer(block_size)
  }
  if (is.factor(groups)) groups <- as.character(groups)
  engine <- .encouragement_load(object)
  pred <- engine$predict(summary_only = summary_only, groups = groups,
                         weights = if (is.null(weights)) NULL else as.numeric(weights),
                         block_size = block_size, alpha = alpha)
  standardized <- pred$standardized
  covariance <- pred$itt_covariance()
  result <- list(
    effects = .encouragement_plain(pred$effects(min_ess = min_ess, max_rhat = max_rhat)),
    reference = .encouragement_plain(pred$reference),
    wald = .encouragement_plain(pred$wald(min_ess = min_ess, max_rhat = max_rhat)),
    wald_unchecked = .encouragement_plain(pred$wald(
      require_convergence = FALSE, min_ess = min_ess, max_rhat = max_rhat)),
    first_stage = .encouragement_plain(pred$first_stage(
      practical_threshold = practical_threshold, min_ess = min_ess, max_rhat = max_rhat)),
    comparison = .encouragement_plain(pred$comparison()),
    diagnostics = .encouragement_plain(pred$stability(min_ess = min_ess, max_rhat = max_rhat)),
    draws = .encouragement_plain(pred$draws),
    itt_covariance = .encouragement_plain(covariance$to_numpy()),
    itt_covariance_index = .encouragement_plain(covariance$index$to_frame(index = FALSE)),
    weights = .encouragement_plain(standardized$weights),
    group_labels = unlist(.encouragement_plain(standardized$group_labels), use.names = FALSE),
    group_counts = .encouragement_plain(standardized$group_counts),
    periods = .encouragement_plain(standardized$periods),
    period_indices = .encouragement_plain(standardized$period_indices),
    horizons = .encouragement_plain(standardized$horizons),
    metadata = .encouragement_plain(pred$metadata),
    alpha = alpha
  )
  if (identical(result$metadata$engine, "distributed_lag")) {
    # Continuous intensity: the lag kernel beside the design-based Toeplitz kernel,
    # the z_{tier - T0} discrepancy, design diagnostics and the cumulative return.
    result$kernel <- .encouragement_plain(pred$kernel())
    result$discrepancy <- .encouragement_plain(pred$discrepancy())
    result$design_diagnostics <- .encouragement_plain(pred$diagnostics())
    result$cumulative_return <- .encouragement_plain(pred$cumulative_return(
      min_ess = min_ess, max_rhat = max_rhat))
  }
  structure(result, class = "longbet_encourage.pred")
}

#' Extract encouragement posterior and comparison tables
#'
#' These return tables computed by the Python engine during prediction; they
#' require no live Python process and survive `saveRDS()`.
#' @param pred A `longbet_encourage.pred` object.
#' @param require_convergence Withhold posterior ratio intervals failing the
#'   diagnostic thresholds selected at prediction. Raw draws remain available.
#' @return A data frame. `encouragement_wald()` returns posterior ratio summaries;
#'   frequentist confidence sets are instead `pred$reference`.
#'   `encouragement_first_stage()` separates precision and practical relevance.
#'   `encouragement_comparison()` gives descriptive differences and does not
#'   fabricate a disagreement flag without paired sampling covariance.
#' @export
encouragement_wald <- function(pred, require_convergence = TRUE) {
  if (!inherits(pred, "longbet_encourage.pred")) stop("pred must be an encouragement prediction.", call. = FALSE)
  if (!is.logical(require_convergence) || length(require_convergence) != 1L || is.na(require_convergence)) {
    stop("require_convergence must be TRUE or FALSE.", call. = FALSE)
  }
  if (require_convergence) pred$wald else pred$wald_unchecked
}

#' @rdname encouragement_wald
#' @export
encouragement_first_stage <- function(pred) {
  if (!inherits(pred, "longbet_encourage.pred")) stop("pred must be an encouragement prediction.", call. = FALSE)
  pred$first_stage
}

#' @rdname encouragement_wald
#' @export
encouragement_comparison <- function(pred) {
  if (!inherits(pred, "longbet_encourage.pred")) stop("pred must be an encouragement prediction.", call. = FALSE)
  pred$comparison
}

#' @export
print.longbet_encourage <- function(x, ...) {
  if (identical(x$metadata$engine, "distributed_lag")) {
    cat(sprintf("LongBet encouragement model (continuous intensity, distributed lags, tier %s = %s): %s outcome\n",
                x$metadata$tier, x$metadata$identification, x$metadata$outcome))
    cat(sprintf("  lag horizon %d; Stage-2 backend %s\n", as.integer(x$metadata$lag_horizon),
                x$metadata$stage2_backend))
  } else {
    cat(sprintf("LongBet encouragement model (adoption clock): %s outcome\n",
                x$metadata$outcome))
  }
  cat(sprintf("  %d study units, %d periods; posterior calibration not established\n",
              x$metadata$n_units, x$metadata$n_periods))
  invisible(x)
}

#' @export
print.longbet_encourage.pred <- function(x, ...) {
  cat(sprintf("Encouragement prediction: %d groups, %d horizons (%s)\n",
              length(x$group_labels), length(x$horizons), x$metadata$standardization))
  cat("  Posterior summaries and design-based reference confidence sets are separate.\n")
  invisible(x)
}

#' Paired bootstrap comparison of model and reference encouragement effects
#'
#' Refits both estimators on the same whole-unit bootstrap samples, retaining
#' their paired covariance. This is model/reference comparison, not a new IV
#' identification assumption or an exact randomization test.
#' @param object A fitted `longbet_encourage` object.
#' @param replicates Number of refitted bootstrap samples; at least 100 and no
#'   failed fits/diagnostics are needed before any disagreement flag is issued.
#' @param groups Optional fixed baseline group labels.
#' @param random_seed Seed for bootstrap sampling and refit seed generation.
#' @param alpha,min_ess,max_rhat Interval and diagnostic settings.
#' @return Plain-R list containing `table`, `replicates`, `failures`, and
#'   `metadata`. Replicate axes are `(replicate, group, horizon, quantity,
#'   estimator)`, with quantities Y,D and estimators model,reference. Failures
#'   remain missing and are reported rather than removed.
#' @export
encouragement_bootstrap_comparison <- function(object, replicates = 200L, groups = NULL,
                                              random_seed = 0L, alpha = 0.05,
                                              min_ess = 400, max_rhat = 1.01) {
  for (name in c("replicates", "random_seed")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != trunc(value) || value < 0 || value > .Machine$integer.max) {
      stop(sprintf("%s must be a nonnegative integer.", name), call. = FALSE)
    }
  }
  if (is.factor(groups)) groups <- as.character(groups)
  result <- .encouragement_load(object)$bootstrap_comparison(
    replicates = as.integer(replicates), groups = groups,
    random_seed = as.integer(random_seed), alpha = alpha,
    min_ess = min_ess, max_rhat = max_rhat)
  list(table = .encouragement_plain(result$table),
       replicates = .encouragement_plain(result$replicates),
       failures = .encouragement_plain(result$failures),
       metadata = .encouragement_plain(result$metadata))
}

.encouragement_mapping <- function(value, name) {
  if (is.null(value)) return(NULL)
  if (is.null(names(value)) || anyNA(names(value)) || any(!nzchar(names(value))) ||
      anyDuplicated(names(value))) {
    stop(sprintf("%s must have unique nonempty names.", name), call. = FALSE)
  }
  as.list(value)
}

#' Declare the randomization scheme for encouragement inference
#'
#' @param assignment `"complete"`, `"bernoulli"`, or `"staggered"`. The last
#'   means independent categorical cohort assignment, not fixed cohort quotas.
#' @param blocks,clusters Baseline labels, length N. Each cluster must lie within
#'   one block and share an encouragement path.
#' @param target `"unit"` gives equal participant weight; `"cluster"` gives
#'   equal weight to included clusters and equal within-cluster subgroup weights.
#' @param probabilities Known scalar or vector Bernoulli assignment probabilities.
#'   For staggered assignment supply a named list mapping zero-based start column
#'   indices and optionally `"never"` to probabilities; they sum to one per
#'   independent randomization unit. Vectors may have length N or cluster count.
#' @param block_weights Optional named fixed nonnegative block weights summing
#'   to one. Every positive-weight block must contain the subgroup.
#' @param cohort_weights Required for staggered designs: named positive weights
#'   summing to one over the finite cohorts of interest, using zero-based indices.
#' @param control `"not_yet"` uses later cohorts and permanent holdouts;
#'   `"never"` uses holdouts. Both require no anticipation.
#' @return A plain R `encouragement_design` specification. Array and support
#'   checks occur in the Python engine when effects are computed.
#' @export
encouragement_design <- function(assignment = "complete", blocks = NULL, clusters = NULL,
                                 target = "unit", probabilities = NULL, block_weights = NULL,
                                 cohort_weights = NULL, control = "not_yet") {
  assignment <- match.arg(assignment, c("complete", "bernoulli", "staggered"))
  target <- match.arg(target, c("unit", "cluster"))
  control <- match.arg(control, c("not_yet", "never"))
  if (is.factor(blocks)) blocks <- as.character(blocks)
  if (is.factor(clusters)) clusters <- as.character(clusters)
  if (assignment == "staggered") probabilities <- .encouragement_mapping(probabilities, "probabilities")
  structure(list(assignment = assignment, blocks = blocks, clusters = clusters,
                 target = target, probabilities = probabilities,
                 block_weights = .encouragement_mapping(block_weights, "block_weights"),
                 cohort_weights = .encouragement_mapping(cohort_weights, "cohort_weights"),
                 control = control), class = "encouragement_design")
}

#' Encouragement effects under an explicit randomization design
#'
#' @param y,d,z Complete outcome, actual-adoption and encouragement panels `[N x T]`.
#' @param t Increasing calendar vector with whole-unit gaps; defaults to `1:T`.
#' @param design An [encouragement_design()] specification; NULL selects
#'   individual-unit complete randomization with equal unit weights.
#' @param groups Optional prespecified baseline subgroup labels, length N.
#' @param alpha Tail probability for pointwise normal intervals and AR-style sets.
#' @return A plain-R list with `table`, numeric `covariance`, `covariance_index`,
#'   and `metadata`. Rows and columns of `covariance` both follow the rows of
#'   `covariance_index`, whose columns are `contrast_id` and `quantity` (Y or D ITT).
#' @details Uses fixed-population Horvitz--Thompson contrasts and joint
#'   randomization-unit covariance, preserving shared controls, cross-horizon
#'   dependence and small or unsupported groups. Blocked complete designs use
#'   Neyman covariance. Bernoulli/categorical designs use a conservative uncentered
#'   covariance bound that can be loose and depends on outcome location.
#'   Subgroup HT estimates do not divide by random subgroup arm sizes.
#'   Missing cohort support makes a specified aggregate unavailable without
#'   renormalizing weights. Intervals are pointwise normal approximations, not
#'   simultaneous bands, exact permutation tests or posterior intervals.
#' @export
design_encouragement_effects <- function(y, d, z, t = NULL, design = NULL,
                                        groups = NULL, alpha = 0.05) {
  if (is.null(design)) design <- encouragement_design()
  if (!inherits(design, "encouragement_design")) {
    stop("design must be an encouragement_design() specification.", call. = FALSE)
  }
  if (is.factor(groups)) groups <- as.character(groups)
  longbet_py()
  lb <- reticulate::import("longbet", convert = FALSE)
  declared <- do.call(lb$EncouragementDesign, unclass(design))
  result <- lb$design_encouragement_effects(y = as.matrix(y), d = as.matrix(d), z = as.matrix(z),
    t = if (is.null(t)) NULL else as.numeric(t), design = declared, groups = groups, alpha = alpha)
  list(table = .encouragement_plain(result$table),
       covariance = .encouragement_plain(result$covariance$to_numpy()),
       covariance_index = .encouragement_plain(result$covariance$index$to_frame(index = FALSE)),
       metadata = .encouragement_plain(result$metadata))
}

#' Predict conditional encouragement effects and business decisions
#'
#' Unit-level offer effects on the outcome (CITT_Y) and on adoption (CITT_D),
#' their ratio (CACE), principal strata and risk-aware decision metrics, for the
#' study units or a new cohort. Each unit is described as a unit with its
#' covariate profile drawn from the fitted population: both unit intercepts are
#' integrated out, so the same rule applies to study and new units.
#'
#' @param object A fitted `longbet_encourage` object.
#' @param new_x Optional baseline covariates matrix for new units, `[N_new x P]`.
#'   If NULL, evaluates the covariate profiles of the study units.
#' @param new_z Optional counterfactual encouragement schedule, `[N_new x T]`.
#'   If NULL, assumes encouragement is assigned from the encouragement start period onward.
#' @param t Optional calendar time vector; defaults to the fitted calendar.
#' @param cost Optional positive scalar per-unit encouragement cost.
#' @param hurdle Decision certainty hurdle in (0, 1), default 0.50 (risk neutral).
#' @param alpha Tail probability for credible intervals (default 0.05 for 95% intervals).
#' @param cace_stabilization Positive floor added to the adoption contrast in the
#'   draw-wise CACE ratio (default 0.02).
#' @param monotonic_first_stage If TRUE (default), clips the adoption contrast at
#'   zero (no defiers).
#' @param budget Optional total monetary budget constraint on encouragement.
#' @param capacity Optional integer cap on the maximum number of accounts targeted.
#' @param ranking_metric Prioritization metric for resource-constrained selection:
#'   `"expected_net_value"` (default), `"certainty_adjusted"`, or `"breakeven_probability"`.
#' @param ... Reserved for future expansion.
#' @return A `longbet_conditional_pred` list containing `citt_y`, `citt_d`,
#'   `cace`, `cumulative_lift`, `breakeven_prob`, `decision`, `policy_value`,
#'   `periods`, `horizons`, `cost`, `hurdle`, `budget`, `capacity`, and `alpha`.
#' @export
predict_conditional <- function(object, ...) {
  UseMethod("predict_conditional")
}

#' @export
predict_conditional.longbet_encourage <- function(object, new_x = NULL, new_z = NULL, t = NULL,
                                                  cost = NULL, hurdle = 0.50, alpha = 0.05,
                                                  cace_stabilization = 0.02,
                                                  monotonic_first_stage = TRUE,
                                                  budget = NULL, capacity = NULL,
                                                  ranking_metric = "expected_net_value", ...) {
  .reject_unsupported(...)
  engine <- .encouragement_load(object)
  x_py <- if (is.null(new_x)) NULL else reticulate::r_to_py(as.matrix(new_x))
  z_py <- if (is.null(new_z)) NULL else reticulate::r_to_py(as.matrix(new_z))
  t_py <- if (is.null(t)) NULL else reticulate::r_to_py(as.numeric(t))
  pred <- engine$predict_conditional(
    x = x_py, z = z_py, t = t_py, alpha = alpha, cace_stabilization = cace_stabilization,
    monotonic_first_stage = monotonic_first_stage
  )

  citt_y_draws <- .encouragement_plain(pred$citt_y$draws)
  citt_y_mean <- .encouragement_plain(pred$citt_y$mean)
  citt_y_median <- .encouragement_plain(pred$citt_y$median)
  citt_y_lower <- .encouragement_plain(pred$citt_y$lower)
  citt_y_upper <- .encouragement_plain(pred$citt_y$upper)
  citt_y_sd <- .encouragement_plain(pred$citt_y$sd)

  citt_d_draws <- .encouragement_plain(pred$citt_d$draws)
  citt_d_mean <- .encouragement_plain(pred$citt_d$mean)
  citt_d_median <- .encouragement_plain(pred$citt_d$median)
  citt_d_lower <- .encouragement_plain(pred$citt_d$lower)
  citt_d_upper <- .encouragement_plain(pred$citt_d$upper)
  citt_d_sd <- .encouragement_plain(pred$citt_d$sd)

  cace_draws <- .encouragement_plain(pred$cace$draws)
  cace_mean <- .encouragement_plain(pred$cace$mean)
  cace_median <- .encouragement_plain(pred$cace$median)
  cace_lower <- .encouragement_plain(pred$cace$lower)
  cace_upper <- .encouragement_plain(pred$cace$upper)
  cace_sd <- .encouragement_plain(pred$cace$sd)

  periods <- unlist(.encouragement_plain(pred$periods), use.names = FALSE)
  horizons <- unlist(.encouragement_plain(pred$horizons), use.names = FALSE)

  # Cumulative lift across post-treatment horizons: [N x D]
  cum_draws <- apply(citt_y_draws, c(1, 3), sum)
  cum_mean <- rowMeans(cum_draws)
  cum_median <- apply(cum_draws, 1, stats::median)
  cum_lower <- apply(cum_draws, 1, stats::quantile, probs = alpha / 2)
  cum_upper <- apply(cum_draws, 1, stats::quantile, probs = 1 - alpha / 2)

  breakeven_prob <- if (!is.null(cost)) rowMeans(cum_draws > cost) else NULL

  decision <- if (!is.null(cost)) {
    if (!is.null(budget) || !is.null(capacity)) {
      as.logical(.encouragement_plain(pred$knapsack_policy(
        cost = as.numeric(cost),
        budget = if (!is.null(budget)) as.numeric(budget) else NULL,
        capacity = if (!is.null(capacity)) as.integer(capacity) else NULL,
        hurdle = as.numeric(hurdle),
        ranking_metric = ranking_metric
      )))
    } else {
      (breakeven_prob >= hurdle)
    }
  } else NULL

  policy_val <- if (!is.null(decision) && !is.null(cost)) {
    mean(ifelse(decision, cum_mean - as.numeric(cost), 0.0))
  } else NULL

  strata_df <- tryCatch({
    as.data.frame(.encouragement_plain(pred$principal_strata()))
  }, error = function(e) NULL)

  intensity <- reticulate::py_has_attr(pred, "kernel")
  effect_summary <- function(effect) {
    list(mean = .encouragement_plain(effect$mean), median = .encouragement_plain(effect$median),
         lower = .encouragement_plain(effect$lower), upper = .encouragement_plain(effect$upper),
         sd = .encouragement_plain(effect$sd), draws = .encouragement_plain(effect$draws))
  }

  structure(c(list(
    citt_y = list(mean = citt_y_mean, median = citt_y_median, lower = citt_y_lower,
                  upper = citt_y_upper, sd = citt_y_sd, draws = citt_y_draws),
    citt_d = list(mean = citt_d_mean, median = citt_d_median, lower = citt_d_lower,
                  upper = citt_d_upper, sd = citt_d_sd, draws = citt_d_draws),
    cace = list(mean = cace_mean, median = cace_median, lower = cace_lower,
                upper = cace_upper, sd = cace_sd, draws = cace_draws),
    cumulative_lift = list(draws = cum_draws, mean = cum_mean, median = cum_median,
                           lower = cum_lower, upper = cum_upper),
    breakeven_prob = breakeven_prob,
    decision = decision,
    policy_value = policy_val,
    cost = cost,
    hurdle = hurdle,
    budget = budget,
    capacity = capacity,
    ranking_metric = ranking_metric,
    periods = periods,
    horizons = horizons,
    alpha = alpha,
    strata = strata_df
  ), if (intensity) list(
    # Continuous intensity: the conditional lag kernel nu_Y(x) kappa_l, its running
    # sum and the cumulative return per unit of cumulative dose (`cace` aliases it).
    kernel = effect_summary(pred$kernel),
    kernel_cum = effect_summary(pred$kernel_cum),
    omega_cum = effect_summary(pred$omega_cum),
    kernel_summary = .encouragement_plain(pred$kernel_summary())
  ) else list()), class = "longbet_conditional_pred")
}

#' Multi-objective resource-constrained encouragement policy
#'
#' Prioritizes encouragement targeting under hard budget or account capacity constraints.
#'
#' @param pred A `longbet_conditional_pred` object from [predict_conditional()].
#' @param cost Per-unit encouragement cost.
#' @param budget Optional total monetary expenditure cap.
#' @param capacity Optional integer cap on maximum treated accounts.
#' @param hurdle Decision certainty hurdle in (0, 1), default 0.50.
#' @param ranking_metric Prioritization metric: `"expected_net_value"` (default),
#'   `"certainty_adjusted"`, or `"breakeven_probability"`.
#' @return Logical vector of length N indicating targeting decision.
#' @export
knapsack_policy <- function(pred, cost, budget = NULL, capacity = NULL,
                            hurdle = 0.50, ranking_metric = "expected_net_value") {
  if (!inherits(pred, "longbet_conditional_pred")) {
    stop("pred must be a longbet_conditional_pred object.", call. = FALSE)
  }
  cost <- as.numeric(cost)
  if (!is.finite(cost) || cost < 0) stop("cost must be a nonnegative numeric scalar.", call. = FALSE)
  
  expected_net <- pred$cumulative_lift$mean - cost
  p_breakeven <- if (!is.null(pred$breakeven_prob)) pred$breakeven_prob else rowMeans(pred$cumulative_lift$draws > cost)
  sd_lift <- apply(pred$cumulative_lift$draws, 1, stats::sd)

  n_units <- length(expected_net)
  eligible_idx <- which(p_breakeven >= hurdle & expected_net > 0)

  max_units <- n_units
  if (!is.null(capacity)) max_units <- min(max_units, as.integer(capacity))
  if (!is.null(budget)) max_units <- min(max_units, as.integer(floor(budget / cost)))

  decision <- rep(FALSE, n_units)
  if (max_units <= 0 || length(eligible_idx) == 0) return(decision)

  if (length(eligible_idx) <= max_units) {
    decision[eligible_idx] <- TRUE
    return(decision)
  }

  scores <- switch(
    ranking_metric,
    "expected_net_value" = expected_net,
    "certainty_adjusted" = expected_net / (sd_lift + 1e-6),
    "breakeven_probability" = p_breakeven,
    stop(sprintf("Unknown ranking_metric '%s'.", ranking_metric), call. = FALSE)
  )

  eligible_scores <- scores[eligible_idx]
  order_idx <- order(-eligible_scores)
  selected <- eligible_idx[order_idx[seq_len(max_units)]]
  decision[selected] <- TRUE
  decision
}

#' Estimate principal strata probabilities and counts (Compliers, Never-Takers, Always-Takers)
#'
#' @param object A `longbet_conditional_pred` object from [predict_conditional()].
#' @param horizon Optional exposure horizon to keep, matching the `horizon` column
#'   (1 is the first post-encouragement period when periods are unit-spaced).
#'   If NULL (default), returns all horizons.
#' @return A data.frame with columns `period`, `horizon`, `stratum`, `prob_mean`, `prob_median`,
#'   `prob_lower`, `prob_upper`, `count_mean`, `count_median`, `count_lower`, `count_upper`, `n_total`.
#' @export
principal_strata <- function(object, horizon = NULL) {
  UseMethod("principal_strata")
}

#' @export
principal_strata.longbet_conditional_pred <- function(object, horizon = NULL) {
  df <- object$strata
  if (is.null(df)) {
    stop("Principal strata estimates are not available on this object.")
  }
  if (!is.null(horizon)) {
    df <- df[df$horizon == as.integer(horizon), , drop = FALSE]
  }
  df
}

#' Print human-readable summary of principal strata estimates
#'
#' @param object A `longbet_conditional_pred` object from [predict_conditional()].
#' @param horizon Optional exposure horizon, matching the `horizon` column.
#'   Defaults to the first available horizon.
#' @export
strata_summary <- function(object, horizon = NULL) {
  UseMethod("strata_summary")
}

#' @export
strata_summary.longbet_conditional_pred <- function(object, horizon = NULL) {
  df <- object$strata
  if (is.null(df) || nrow(df) == 0) {
    cat("Principal strata estimates are not available.\n")
    return(invisible(NULL))
  }
  available_h <- sort(unique(df$horizon))
  target_h <- if (is.null(horizon)) available_h[1] else as.integer(horizon)
  sub <- df[df$horizon == target_h, ]
  if (nrow(sub) == 0) {
    stop(sprintf("Horizon %s not found. Available horizons: %s", target_h, paste(available_h, collapse = ", ")))
  }
  p_val <- sub$period[1]
  n_units <- sub$n_total[1]
  level <- if (is.numeric(object$alpha) && length(object$alpha) == 1L) 100 * (1 - object$alpha) else 95
  cat(sprintf("Principal Strata Estimates at Horizon %s (Period %.2f, N = %d):\n", target_h, p_val, n_units))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    st <- gsub("_", "-", tools::toTitleCase(as.character(r$stratum)))
    cat(sprintf("  - %-13s: %5.1f%% (%g%% CI: [%.1f%%, %.1f%%])  --  ~%.1f units (%g%% CI: [%.1f, %.1f])\n",
                st, r$prob_mean * 100, level, r$prob_lower * 100, r$prob_upper * 100,
                r$count_mean, level, r$count_lower, r$count_upper))
  }
  c_row <- sub[sub$stratum == "complier", ]
  if (nrow(c_row) > 0) {
    if (c_row$prob_mean < 0.10) {
      cat("  [WARNING: Weak Instrument detected - complier share < 10%]\n")
    } else {
      cat("  [Robust Instrument - complier share exceeds 10% threshold]\n")
    }
  }
  invisible(sub)
}



