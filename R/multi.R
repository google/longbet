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

#' LongBetMulti: Coupled Bayesian ensemble trees for multiple outcomes
#'
#' Fits outcome-specific prognostic forests on a shared panel,
#' using the full precision of a triangular SUR (seemingly unrelated regressions)
#' likelihood (`full_precision_sur_v1`). Downstream observations inform earlier
#' means and discrete latent responses. Binary and ordinal equations come first,
#' preserving their caller order, with zero incoming loadings and unit marginal
#' latent variance. Discrete outcomes have no free residual correlation parameter.
#' Note: aligned draws do not establish calibrated
#' joint inference. Treatment trees are outcome-specific; only the residual
#' likelihood couples the outcomes.
#' Check convergence of the effects and decision events before interpreting
#' their probabilities. Coupled fits containing continuous outcomes require
#' explicitly positive `sigma_prior_a` and `sigma_prior_b`; the default
#' improper reference prior can produce an improper joint posterior and is
#' rejected.
#'
#' @param y Outcome panels: a list of `[N x T]` matrices, or a 3-D numeric array
#'   `[N x T x M]`.
#' @param x Covariates for the prognostic forest, `[N x P]`.
#' @param z Binary treatment indicator, `[N x T]`. Must be absorbing.
#' @param t Calendar time vector of length `T`. Defaults to `1:ncol(y[[1]])`.
#' @param x_trt Covariates for the treatment forest, `[N x P_trt]`. Defaults to
#'   `x`.
#' @param x_tv Time-varying prognostic covariates, `[N x T x P]`.
#' @param x_trt_tv Time-varying treatment covariates, `[N x T x P]`.
#' @param ps Propensity score, length `N` or `[N x T]`.
#' @param num_sweeps Posterior draws saved per chain, after burn-in.
#' @param num_burnin Burn-in iterations discarded.
#' @param n_skip Thinning interval.
#' @param num_chains Parallel MCMC chains.
#' @param inner_loop_length Sweeps per dispatch.
#' @param num_trees_pr,num_trees_trt Trees in each forest.
#' @param min_points_per_leaf_pr,min_points_per_leaf_trt Minimum leaf sizes.
#' @param max_depth_pr,max_depth_trt Maximum tree depths.
#' @param alpha_split_pr,beta_split_pr Tree-depth prior, prognostic forest.
#' @param alpha_split_trt,beta_split_trt Tree-depth prior, treatment forest.
#' @param num_cutpoints Maximum quantile bins per continuous covariate.
#' @param sig_knl Kernel standard deviation for the exposure trajectory.
#' @param lambda_knl Kernel lengthscale over the exposure index.
#' @param kernel_type One of `"se"`, `"matern32"`, `"matern52"`, `"ar1"`.
#' @param sigma_m Prior standard deviation of the marginalized constant GP mean.
#' @param gp_constant_mean Whether the projection reverts to an estimated common level.
#' @param trend_num_features Random Fourier features in the unit feature map of
#'   the calendar-time block. More capacity to represent a trend whose covariate
#'   loading is nonlinear, at more cost per sweep. The horseshoe prior on the
#'   block shrinks each column on its own evidence, so raising this is close to
#'   free statistically; the cost is arithmetic.
#' @param random_intercept Whether to fit unit random intercepts.
#' @param gamma_prior_a,gamma_prior_b Inverse-gamma prior on unit-intercept variance.
#' @param sigma_prior_a,sigma_prior_b Inverse-gamma prior on innovation variance,
#'   proper by construction; the default IG(2, 1) has prior mean variance one on
#'   the standardized scale. When `standardize=FALSE`, choose these in outcome units.
#' @param outcome Outcome type(s): `"continuous"`, `"binary"`, `"ordinal"`, or a vector of length `M`.
#' @param num_categories Integer broadcast to ordinal outcomes, or a list in
#'   user order (optionally fully named), with `NULL` for nonordinal outcomes.
#'   Every ordinal outcome requires a count >=2 and explicit numeric labels
#'   `0,...,K-1`; empty categories are allowed. Counts are never inferred.
#' @param cutpoint_prior_scale Shared finite positive ordered-normal cutpoint
#'   prior scale in latent probit units. Inactive for nonordinal outcomes or K=2.
#' @param outcome_names Optional character vector of length `M` naming the outcomes.
#' @param standardize Whether to centre and scale continuous outcomes internally.
#' @param random_seed Base PRNG seed.
#' @param device `"auto"`, `"cpu"`, or `"gpu"`.
#' @param sur Whether to enable full-precision triangular SUR coupling.
#' @param sur_prior_var Prior variance for loading regression parameters.
#' @param verbose Whether to print progress messages.
#' @param ... Reserved; unsupported arguments raise an error.
#' @return An object of class `longbet_multi`.
#' @seealso [predict.longbet_multi()], [effect_draws()], [joint_prob()], [outcome_correlation()]
#' @examples
#' \dontrun{
#' n <- 100; Tn <- 8
#' x <- matrix(rnorm(n * 3), n, 3)
#' z <- matrix(0, n, Tn); z[1:50, 5:Tn] <- 1
#' y1 <- 0.5 * x[, 1] + 1.2 * z + matrix(rnorm(n * Tn, 0, 0.3), n, Tn)
#' y2 <- matrix(as.numeric(0.3 * x[, 2] - 0.8 * z + rnorm(n * Tn) > 0), n, Tn)
#' fit <- longbet_multi(
#'   y = list(rev = y1, churn = y2),
#'   x = x, z = z, t = 1:Tn,
#'   outcome = c(rev = "continuous", churn = "binary"),
#'   sigma_prior_a = 2, sigma_prior_b = 1,
#'   num_chains = 2
#' )
#' pred <- predict(fit, x = x, z = z, t = 1:Tn, summary_only = TRUE)
#' }
#' @export
longbet_multi <- function(y, x, z, t = NULL,
                          x_trt = NULL, x_tv = NULL, x_trt_tv = NULL, ps = NULL,
                          num_sweeps = 1000, num_burnin = 2000, n_skip = 1,
                          num_chains = 4, inner_loop_length = NULL,
                          num_trees_pr = 20, num_trees_trt = 60,
                          min_points_per_leaf_pr = 10, min_points_per_leaf_trt = 10,
                          max_depth_pr = 8, max_depth_trt = 8,
                          alpha_split_pr = 0.95, beta_split_pr = 2.0,
                          alpha_split_trt = 0.25, beta_split_trt = 3.0,
                          num_cutpoints = 100,
                          sig_knl = 1.0, lambda_knl = 1.0, kernel_type = "se",
                          sigma_m = 1.0, gp_constant_mean = TRUE,
                          trend_num_features = 64,
                          random_intercept = TRUE,
                          gamma_prior_a = 1.0, gamma_prior_b = 0.1,
                          sigma_prior_a = 2.0, sigma_prior_b = 1.0,
                          outcome = "continuous", outcome_names = NULL,
                          standardize = TRUE, random_seed = 0,
                          device = c("auto", "cpu", "gpu"),
                          sur = TRUE, sur_prior_var = 1.0,
                          verbose = FALSE, num_categories = NULL,
                          cutpoint_prior_scale = 5.0, ...) {

  .reject_unsupported(...)
  device <- match.arg(device)
  .check_cutpoint_prior(cutpoint_prior_scale)

  # ------------------------------------------------------------------
  # 1. Container validation and slice extraction
  # ------------------------------------------------------------------
  y_list <- list()
  container_names <- NULL

  if (is.array(y)) {
    if (length(dim(y)) == 2L) {
      stop("Bare 2-D panels are invalid for longbet_multi; use longbet() for single-outcome models.",
           call. = FALSE)
    }
    if (length(dim(y)) != 3L) {
      stop(sprintf("y must be a list of matrices or a 3-D numeric array [N x T x M], got ndim=%d.",
                   length(dim(y))), call. = FALSE)
    }
    d <- dim(y)
    N <- d[1]
    T_dim <- d[2]
    M <- d[3]
    if (M < 2L) {
      stop(sprintf("longbet_multi requires at least 2 outcomes, got %d.", M), call. = FALSE)
    }
    dn <- dimnames(y)
    if (!is.null(dn) && length(dn) >= 3L && !is.null(dn[[3]])) {
      container_names <- as.character(dn[[3]])
    }
    for (m in seq_len(M)) {
      sl <- y[, , m, drop = FALSE]
      dim(sl) <- c(N, T_dim)
      if (!is.numeric(sl)) {
        stop(sprintf("Outcome at index %d must be numeric.", m), call. = FALSE)
      }
      y_list[[m]] <- sl
    }
  } else if (is.list(y)) {
    M <- length(y)
    if (M < 2L) {
      stop(sprintf("longbet_multi requires at least 2 outcomes, got %d.", M), call. = FALSE)
    }
    nms <- names(y)
    if (!is.null(nms)) {
      if (any(nms == "") || anyNA(nms)) {
        stop("Partially named outcome lists are not supported; provide all names or none.",
             call. = FALSE)
      }
      if (length(unique(nms)) != length(nms)) {
        stop("Duplicate outcome names in list.", call. = FALSE)
      }
      container_names <- nms
    }
    for (m in seq_len(M)) {
      mat <- y[[m]]
      if (!is.matrix(mat) || !is.numeric(mat)) {
        stop(sprintf("Each outcome in list must be a numeric matrix; outcome %d is not.", m),
             call. = FALSE)
      }
      if (m == 1L) {
        N <- nrow(mat)
        T_dim <- ncol(mat)
      } else if (!identical(dim(mat), c(N, T_dim))) {
        stop(sprintf("All outcome panels must have matching dimensions [N x T]=[%d x %d]; outcome %d is [%d x %d].",
                     N, T_dim, m, nrow(mat), ncol(mat)), call. = FALSE)
      }
      y_list[[m]] <- mat
    }
  } else {
    stop("y must be a list of matrices or a 3-D numeric array [N x T x M].", call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 2. Resolve outcome names
  # ------------------------------------------------------------------
  if (!is.null(outcome_names)) {
    if (length(outcome_names) != M) {
      stop(sprintf("outcome_names has length %d, but there are %d outcomes.",
                   length(outcome_names), M), call. = FALSE)
    }
    outcome_names <- as.character(outcome_names)
    if (anyNA(outcome_names) || any(!nzchar(trimws(outcome_names)))) {
      stop("Outcome names must be non-empty strings.", call. = FALSE)
    }
    if (length(unique(outcome_names)) != M) {
      stop("Outcome names must be unique.", call. = FALSE)
    }
    if (!is.null(container_names) && !identical(outcome_names, container_names)) {
      stop("Explicit outcome_names conflict with names provided in y container.", call. = FALSE)
    }
    final_names <- outcome_names
  } else if (!is.null(container_names)) {
    final_names <- container_names
  } else {
    final_names <- paste0("outcome_", seq_len(M))
  }
  if (anyNA(final_names) || any(!nzchar(trimws(final_names))) || anyDuplicated(final_names)) {
    stop("Outcome names must be unique non-empty strings.", call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 3. Resolve outcome types per element
  # ------------------------------------------------------------------
  resolved_outcomes <- character(M)
  if (is.character(outcome) && length(outcome) == 1L) {
    if (!outcome %in% c("continuous", "binary", "ordinal")) {
      stop(sprintf("outcome must be 'continuous', 'binary', or 'ordinal', got '%s'.", outcome), call. = FALSE)
    }
    resolved_outcomes <- rep(outcome, M)
  } else if (is.character(outcome) && length(outcome) == M) {
    if (!is.null(names(outcome))) {
      if (!setequal(names(outcome), final_names)) {
        stop("Named outcome vector keys must exactly match outcome names.", call. = FALSE)
      }
      resolved_outcomes <- outcome[final_names]
    } else {
      resolved_outcomes <- outcome
    }
    for (i in seq_len(M)) {
      if (!resolved_outcomes[i] %in% c("continuous", "binary", "ordinal")) {
        stop(sprintf("outcome for '%s' must be 'continuous', 'binary', or 'ordinal', got '%s'.",
                     final_names[i], resolved_outcomes[i]), call. = FALSE)
      }
    }
  } else {
    stop("outcome must be a single string ('continuous', 'binary', or 'ordinal') or a vector of length M.",
         call. = FALSE)
  }

  category_counts <- .multi_category_counts(num_categories, resolved_outcomes, final_names)
  for (i in which(resolved_outcomes == "ordinal")) .check_ordinal_labels(y_list[[i]], category_counts[[i]])

  # ------------------------------------------------------------------
  # 4. Validate covariates and panel dimensions
  # ------------------------------------------------------------------
  x <- as.matrix(x)
  z <- as.matrix(z)
  if (nrow(x) != N) {
    stop(sprintf("x has %d rows but y has %d.", nrow(x), N), call. = FALSE)
  }
  if (!identical(dim(z), c(N, T_dim))) {
    stop(sprintf("z is %d x %d but y is %d x %d; they must match.",
                 nrow(z), ncol(z), N, T_dim), call. = FALSE)
  }
  if (is.null(t)) t <- seq_len(T_dim)
  t <- as.numeric(t)

  # ------------------------------------------------------------------
  # 5. Build config and fit Python model
  # ------------------------------------------------------------------
  lb <- longbet_py()
  config <- lb$LongBetConfig(
    num_sweeps = as.integer(num_sweeps),
    num_burnin = as.integer(num_burnin),
    n_skip = as.integer(n_skip),
    num_chains = as.integer(num_chains),
    inner_loop_length = if (is.null(inner_loop_length)) NULL else as.integer(inner_loop_length),
    num_trees_pr = as.integer(num_trees_pr),
    num_trees_trt = as.integer(num_trees_trt),
    min_points_per_leaf_pr = as.integer(min_points_per_leaf_pr),
    min_points_per_leaf_trt = as.integer(min_points_per_leaf_trt),
    max_depth_pr = as.integer(max_depth_pr),
    max_depth_trt = as.integer(max_depth_trt),
    alpha_split_pr = as.numeric(alpha_split_pr),
    beta_split_pr = as.numeric(beta_split_pr),
    alpha_split_trt = as.numeric(alpha_split_trt),
    beta_split_trt = as.numeric(beta_split_trt),
    num_cutpoints = as.integer(num_cutpoints),
    sig_knl = as.numeric(sig_knl),
    lambda_knl = as.numeric(lambda_knl),
    kernel_type = kernel_type,
    sigma_m = as.numeric(sigma_m),
    gp_constant_mean = as.logical(gp_constant_mean),
    trend_num_features = as.integer(trend_num_features),
    random_intercept = as.logical(random_intercept),
    gamma_prior_a = as.numeric(gamma_prior_a),
    gamma_prior_b = as.numeric(gamma_prior_b),
    sigma_prior_a = as.numeric(sigma_prior_a),
    sigma_prior_b = as.numeric(sigma_prior_b),
    outcome = resolved_outcomes[1],
    num_categories = category_counts[[1]],
    cutpoint_prior_scale = as.numeric(cutpoint_prior_scale),
    standardize = as.logical(standardize),
    random_seed = as.integer(random_seed),
    device = device,
    sur = as.logical(sur),
    sur_prior_var = as.numeric(sur_prior_var)
  )

  if (isTRUE(verbose)) {
    message(sprintf("longbet_multi: fitting %d outcomes on %d x %d panel on device '%s' (%d chain%s).",
                    M, N, T_dim, device, num_chains,
                    if (num_chains == 1) "" else "s"))
  }

  py_model <- lb$LongBetMulti(config)

  py_y <- stats::setNames(
    lapply(y_list, function(mat) .as_np_matrix(mat, "y")),
    final_names
  )

  py_model$fit(
    y = py_y,
    x = .as_np_matrix(x, "x"),
    z = .as_np_matrix(z, "z"),
    t = .as_np_vector(t),
    x_trt = .as_np_matrix(x_trt, "x_trt"),
    x_tv = .as_np_array3(x_tv, "x_tv"),
    x_trt_tv = .as_np_array3(x_trt_tv, "x_trt_tv"),
    ps = if (is.null(ps)) NULL else reticulate::r_to_py(if (is.matrix(ps)) as.matrix(ps) else as.numeric(ps)),
    outcome = reticulate::r_to_py(as.list(resolved_outcomes)),
    outcome_names = reticulate::r_to_py(as.list(final_names)),
    num_categories = if (any(resolved_outcomes == "ordinal")) reticulate::r_to_py(category_counts) else NULL
  )

  # Serialize parent model
  tmp <- tempfile(fileext = ".npz")
  on.exit(unlink(tmp), add = TRUE)
  py_model$save(tmp)
  raw_model <- readBin(tmp, what = "raw", n = file.info(tmp)$size)

  # ------------------------------------------------------------------
  # 6. Reconstitute child longbet fits in user order
  # ------------------------------------------------------------------
  np <- reticulate::import("numpy", convert = TRUE)
  as_arr <- function(v) as.array(reticulate::py_to_r(np$asarray(v)))

  flatten_draws <- function(a) {
    if (length(dim(a)) == 3L) {
      matrix(aperm(a, c(3, 2, 1)), nrow = dim(a)[3])
    } else {
      t(a)
    }
  }

  flatten_scalars <- function(a) {
    if (length(dim(a)) == 2L) {
      as.numeric(t(a))
    } else {
      as.numeric(a)
    }
  }

  child_fits <- list()
  for (u in seq_len(M)) {
    py_child <- py_model$fits[[as.integer(u - 1L)]]
    tmp_child <- tempfile(fileext = ".npz")
    py_child$save(tmp_child)
    child_raw <- readBin(tmp_child, what = "raw", n = file.info(tmp_child)$size)
    unlink(tmp_child)

    child_trace <- py_child$trace
    child_beta <- as_arr(child_trace$beta)
    child_gamma <- as_arr(child_trace$gamma)
    child_sigma2 <- flatten_scalars(as_arr(child_trace$sigma2))

    child_handle <- new.env(parent = emptyenv())
    child_handle$py_model <- py_child
    child_handle$pid <- Sys.getpid()

    child_fit <- list(
      beta_values = flatten_draws(child_beta),
      gamma_draws = flatten_draws(child_gamma),
      sigma_gamma_draws = flatten_scalars(as_arr(child_trace$sigma_gamma2)),
      sigma0_draws = matrix(sqrt(child_sigma2), nrow = 1),
      b0_draws = flatten_scalars(as_arr(child_trace$b0)),
      b1_draws = flatten_scalars(as_arr(child_trace$b1)),
      alpha_draws = flatten_scalars(as_arr(child_trace$alpha)),
      sdy = as.numeric(py_child$sdy),
      meany = as.numeric(py_child$meany),
      random_intercept = as.logical(random_intercept),
      outcome = as.character(py_child$config$outcome),
      num_categories = category_counts[[u]],
      cutpoint_prior_scale = as.numeric(cutpoint_prior_scale),
      multi_origin = reticulate::py_to_r(py_child$multi_origin),
      model_params = list(
        burnin = as.integer(num_burnin),
        num_sweeps = as.integer(num_sweeps),
        num_chains = as.integer(num_chains),
        num_trees_pr = as.integer(num_trees_pr),
        num_trees_trt = as.integer(num_trees_trt),
        device = device
      ),
      raw_model = child_raw,
      handle = child_handle,
      t = t
    )
    class(child_fit) <- "longbet"
    child_fits[[final_names[u]]] <- child_fit
  }

  handle <- new.env(parent = emptyenv())
  handle$py_model <- py_model
  handle$pid <- Sys.getpid()

  Gamma_draws <- as.matrix(reticulate::py_to_r(py_model$Gamma_draws))

  fit <- list(
    fits = child_fits,
    outcome_names = final_names,
    outcome = resolved_outcomes,
    num_categories = stats::setNames(category_counts, final_names),
    cutpoint_prior_scale = as.numeric(cutpoint_prior_scale),
    M = M,
    order = as.integer(reticulate::py_to_r(py_model$order)) + 1L,
    inverse_order = as.integer(reticulate::py_to_r(py_model$inverse_order)) + 1L,
    sur = as.logical(sur),
    sur_prior_var = as.numeric(sur_prior_var),
    sur_active = as.logical(py_model$sur_active),
    sampler_semantics = as.character(py_model$sampler_semantics),
    rng_scheme_version = as.integer(py_model$rng_scheme_version),
    provenance = as.character(py_model$provenance),
    Gamma_draws = Gamma_draws,
    model_params = list(
      burnin = as.integer(num_burnin),
      num_sweeps = as.integer(num_sweeps),
      num_chains = as.integer(num_chains),
      num_trees_pr = as.integer(num_trees_pr),
      num_trees_trt = as.integer(num_trees_trt),
      device = device
    ),
    raw_model = raw_model,
    handle = handle,
    t = t
  )
  class(fit) <- "longbet_multi"
  fit
}

#' @export
`[.longbet_multi` <- function(x, i) {
  x$fits[[i]]
}

#' @export
print.longbet_multi <- function(x, ...) {
  p <- x$model_params
  cat("LongBetMulti fit (JAX engine)\n")
  cat(sprintf("  Outcomes: %d (%s)\n", x$M, paste(x$outcome_names, collapse = ", ")))
  cat(sprintf("  Draws:    %d saved, %d burn-in, %d chain%s\n",
              p$num_sweeps, p$burnin, p$num_chains,
              if (p$num_chains == 1) "" else "s"))
  cat(sprintf("  SUR:      %s (active: %s, prior var: %g)\n",
              if (isTRUE(x$sur)) "yes" else "no",
              if (isTRUE(x$sur_active)) "yes" else "no",
              x$sur_prior_var))
  cat(sprintf("  Semantics: %s\n", x$sampler_semantics))
  invisible(x)
}

#' @export
summary.longbet_multi <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

#' Predict from a fitted multi-outcome LongBet model
#'
#' @param object A fitted `longbet_multi` object.
#' @param x Prognostic covariates, `[N x P]`.
#' @param z Treatment panel, `[N x T]`.
#' @param t Calendar time vector. Defaults to fitted time.
#' @param x_trt,x_tv,x_trt_tv,ps Design covariates as in [longbet_multi()].
#' @param summary_only If `TRUE`, return per-cell posterior means and credible
#'   bounds without materializing full `[N x T x draws]` arrays.
#' @param alpha Tail probability for credible intervals.
#' @param random_seed Seed for Gaussian-process projection beyond exposure horizon.
#' @param sig_knl,lambda_knl Optional projection kernel overrides.
#' @param cache_forest_evaluations Whether to cache prognostic forest evaluations.
#' @param ... Reserved; unsupported arguments raise an error.
#' @return An object of class `longbet_multi.pred`.
#' @seealso [longbet_multi()], [effect_draws()], [joint_prob()]
#' @examples
#' \dontrun{
#' pred <- predict(fit, x = x, z = z, t = 1:Tn)
#' pred_rev <- pred["rev"]
#' }
#' @export
predict.longbet_multi <- function(object, x, z, t = NULL,
                                  x_trt = NULL, x_tv = NULL, x_trt_tv = NULL, ps = NULL,
                                  summary_only = FALSE, alpha = 0.05,
                                  random_seed = NULL,
                                  sig_knl = NULL, lambda_knl = NULL,
                                  cache_forest_evaluations = FALSE, ...) {
  .reject_unsupported(...)
  py_model <- rehydrate_jax_model(object)
  np <- reticulate::import("numpy", convert = TRUE)

  x <- as.matrix(x)
  z <- as.matrix(z)
  if (is.null(t)) t <- object$t
  if (is.null(t)) t <- seq_len(ncol(z))
  t <- as.numeric(t)

  key <- NULL
  if (!is.null(random_seed)) {
    jax <- reticulate::import("jax", convert = FALSE)
    key <- jax$random$key(as.integer(random_seed))
  }

  py_pred <- py_model$predict(
    x = .as_np_matrix(x, "x"),
    z = .as_np_matrix(z, "z"),
    t = .as_np_vector(t),
    x_trt = .as_np_matrix(x_trt, "x_trt"),
    x_tv = .as_np_array3(x_tv, "x_tv"),
    x_trt_tv = .as_np_array3(x_trt_tv, "x_trt_tv"),
    ps = if (is.null(ps)) NULL else reticulate::r_to_py(if (is.matrix(ps)) as.matrix(ps) else as.numeric(ps)),
    summary_only = as.logical(summary_only),
    alpha = as.numeric(alpha),
    key = key,
    sig_knl = if (is.null(sig_knl)) NULL else as.numeric(sig_knl),
    lambda_knl = if (is.null(lambda_knl)) NULL else as.numeric(lambda_knl),
    cache_forest_evaluations = as.logical(cache_forest_evaluations)
  )

  arr <- function(v) as.array(np$asarray(v))
  summ <- function(s) list(mean = arr(s$mean), std = arr(s$std),
                           lower = arr(s$lower), upper = arr(s$upper))

  M <- length(object$outcome_names)
  child_preds <- list()

  for (u in seq_len(M)) {
    name <- object$outcome_names[u]
    py_child_pred <- py_pred$preds[[as.integer(u - 1L)]]

    tau <- summ(py_child_pred$tau_summary)
    mu0 <- summ(py_child_pred$mu0_summary)
    yhat <- summ(py_child_pred$y_summary)

    res_u <- list(
      tauhats = if (isTRUE(summary_only)) NULL else arr(py_child_pred$tauhats),
      muhats0 = if (isTRUE(summary_only)) NULL else arr(py_child_pred$muhats0),
      preds = if (isTRUE(summary_only)) NULL else arr(py_child_pred$yhats),
      tauhats.mean = tau$mean, tauhats.sd = tau$std,
      tauhats.lower = tau$lower, tauhats.upper = tau$upper,
      muhats0.mean = mu0$mean, muhats0.sd = mu0$std, muhats0.lower = mu0$lower, muhats0.upper = mu0$upper,
      preds.mean = yhat$mean, preds.sd = yhat$std, preds.lower = yhat$lower, preds.upper = yhat$upper,
      att_full = arr(py_child_pred$att_full),
      beta_values = arr(py_child_pred$beta_values),
      s = arr(py_child_pred$s),
      z = z,
      t = t,
      outcome = py_child_pred$outcome,
      alpha = alpha,
      summary_only = isTRUE(summary_only),
      num_chains = as.integer(py_child_pred$num_chains),
      num_sweeps = as.integer(object$model_params$num_sweeps),
      att_counts = as.integer(py_child_pred$att_counts),
      py_pred = py_child_pred
    )
    res_u <- c(res_u, .ordinal_prediction_fields(py_child_pred))
    class(res_u) <- "longbet.pred"
    child_preds[[name]] <- res_u
  }

  res <- list(
    preds = child_preds,
    outcome_names = object$outcome_names,
    outcome = object$outcome,
    num_categories = object$num_categories,
    summary_only = isTRUE(summary_only),
    alpha = alpha,
    num_chains = as.integer(py_pred$num_chains),
    num_sweeps = as.integer(py_pred$num_sweeps),
    provenance = as.character(py_pred$provenance),
    sampler_semantics = as.character(py_pred$sampler_semantics),
    py_pred = py_pred
  )
  class(res) <- "longbet_multi.pred"
  res
}

#' @export
`[.longbet_multi.pred` <- function(x, i) {
  x$preds[[i]]
}

#' @export
print.longbet_multi.pred <- function(x, ...) {
  cat("LongBetMulti predictions\n")
  cat(sprintf("  Outcomes:     %d (%s)\n", length(x$outcome_names), paste(x$outcome_names, collapse = ", ")))
  cat(sprintf("  Summary only: %s\n", if (isTRUE(x$summary_only)) "yes" else "no"))
  cat(sprintf("  Chains:       %d\n", x$num_chains))
  if (any(x$outcome == "ordinal")) cat("  Ordinal children include latent predictions and category probability summaries.\n")
  invisible(x)
}

#' Extract treatment effect draws on the natural outcome scale
#'
#' Continuous effects remain on the supplied response scale, including log
#' scales. Binary effects are probability differences in [-1, 1]: 0.01 means
#' one percentage point. Both transformations use the Python engine and remain
#' usable after saving and restoring ordinary R prediction arrays.
#' Ordinal selections require an explicit category or score and raise an error;
#' use [att_probabilities()] or [att_expected_score()] for ordinal ATTs.
#'
#' @param pred A `longbet.pred` or `longbet_multi.pred` object.
#' @param outcome Outcome name or 1-based index (required for multi predictions).
#' @return Numeric 3-D array `[N x T x D]`.
#' @seealso [longbet_multi()], [joint_prob()]
#' @examples
#' \dontrun{
#' tau_rev <- effect_draws(pred, outcome = "rev")
#' tau_churn <- effect_draws(pred, outcome = "churn")
#' }
#' @export
effect_draws <- function(pred, outcome = NULL) {
  if (inherits(pred, "longbet_multi.pred")) {
    if (is.null(outcome)) {
      stop("For longbet_multi.pred, `outcome` name or 1-based index must be specified.", call. = FALSE)
    }
    if (is.numeric(outcome) && length(outcome) == 1L) {
      if (!is.finite(outcome) || outcome != trunc(outcome) ||
          outcome < 1L || outcome > length(pred$outcome_names)) {
        stop("Outcome index must be a finite integer in [1, M].", call. = FALSE)
      }
      outcome <- pred$outcome_names[as.integer(outcome)]
    }
    if (length(outcome) != 1L || is.na(outcome) || !outcome %in% pred$outcome_names) {
      stop(sprintf("Unknown outcome '%s'. Expected one of: %s", outcome, paste(pred$outcome_names, collapse = ", ")), call. = FALSE)
    }
    child <- pred$preds[[outcome]]
  } else if (inherits(pred, "longbet.pred")) {
    if (!is.null(outcome) && outcome != pred$outcome) {
      stop(sprintf("Explicit outcome '%s' conflicts with prediction outcome '%s'.", outcome, pred$outcome), call. = FALSE)
    }
    child <- pred
  } else {
    stop("`pred` must be a longbet.pred or longbet_multi.pred object.", call. = FALSE)
  }

  if (identical(child$outcome, "ordinal")) {
    stop("Ordinal effects require selecting a category or score with att_probabilities() or att_expected_score().", call. = FALSE)
  }
  if (isTRUE(child$summary_only) || is.null(child$tauhats)) {
    stop("effect_draws requires full posterior draws. Re-run predict() with summary_only = FALSE.", call. = FALSE)
  }

  # Reconstruct from ordinary arrays, not a potentially dangling py_pred handle.
  # The scale conversion and shape checks have one implementation in Python.
  lb <- longbet_py()
  np <- reticulate::import("numpy", convert = TRUE)
  as.array(np$asarray(lb$effect_draws_from_arrays(
    tauhats = reticulate::r_to_py(child$tauhats),
    muhats0 = if (is.null(child$muhats0)) NULL else reticulate::r_to_py(child$muhats0),
    outcome = child$outcome,
    summary_only = isTRUE(child$summary_only)
  )))
}

#' Empirical joint probability of multiple treatment effect events
#'
#' This scalar-effect interface supports continuous and binary outcomes.
#' Ordinal selections raise an error; use the category/score ATT methods.
#'
#' @param pred A `longbet_multi.pred` object with full draws.
#' @param conditions Named list or sequence of functions in user outcome order.
#'   Each function receives an `[N x T x D]` numeric array of effect draws and
#'   must return an `[N x T x D]` logical array.
#' @param cells Optional logical panel of shape `[N x T]`. When `NULL`, returns
#'   an `[N x T]` matrix of cellwise joint probabilities. When provided, returns
#'   a length-N vector where each unit's probability is averaged over its
#'   selected periods (or `NA` if no periods are selected for that unit).
#' @return Numeric matrix `[N x T]` or length-N vector.
#' @seealso [longbet_multi()], [effect_draws()]
#' @examples
#' \dontrun{
#' p_joint <- joint_prob(pred, list(
#'   rev = function(e) e > 0,
#'   churn = function(e) e < 0
#' ))
#' }
#' @export
joint_prob <- function(pred, conditions, cells = NULL) {
  if (!inherits(pred, "longbet_multi.pred")) {
    stop("`pred` must be a longbet_multi.pred object.", call. = FALSE)
  }
  if (!(length(pred$sampler_semantics) == 1L && pred$sampler_semantics %in%
        c("full_precision_sur_v1", "shared_private_sur_v1"))) {
    stop("Joint prediction uses legacy or unsupported sampler semantics. Refit from the original data and regenerate predictions.", call. = FALSE)
  }
  if (isTRUE(pred$summary_only)) {
    stop("joint_prob requires full posterior draws; predict with summary_only = FALSE.", call. = FALSE)
  }

  M <- length(pred$outcome_names)
  cond_list <- list()

  if (is.list(conditions)) {
    nms <- names(conditions)
    if (!is.null(nms)) {
      if (length(nms) != M || anyNA(nms) || anyDuplicated(nms) ||
          !setequal(nms, pred$outcome_names)) {
        stop("Conditions named list keys must exactly match fitted outcome names.", call. = FALSE)
      }
      cond_list <- conditions[pred$outcome_names]
    } else {
      if (length(conditions) != M) {
        stop(sprintf("Number of conditions (%d) must match number of outcomes (%d).",
                     length(conditions), M), call. = FALSE)
      }
      cond_list <- conditions
    }
  } else {
    stop("`conditions` must be a list of functions.", call. = FALSE)
  }

  for (i in seq_len(M)) {
    if (!is.function(cond_list[[i]])) {
      stop(sprintf("Condition for outcome '%s' is not a function.", pred$outcome_names[i]), call. = FALSE)
    }
  }

  masks <- list()
  ref_dim <- NULL

  for (i in seq_len(M)) {
    name <- pred$outcome_names[i]
    fn <- cond_list[[i]]
    eff <- effect_draws(pred, outcome = name)

    if (is.null(ref_dim)) {
      ref_dim <- dim(eff)
    }
    if (!identical(dim(eff), ref_dim)) {
      stop("All outcomes must have identical [N x T x D] effect dimensions.", call. = FALSE)
    }
    child <- pred$preds[[name]]
    if (!identical(child$num_chains, pred$num_chains) ||
        (!is.null(pred$num_sweeps) && ref_dim[3] != pred$num_chains * pred$num_sweeps)) {
      stop("Draw/chain metadata is inconsistent with the joint prediction.", call. = FALSE)
    }
    if (anyNA(eff) || any(!is.finite(eff))) {
      stop(sprintf("effect_draws for outcome '%s' contains non-finite values.", name), call. = FALSE)
    }

    mask <- fn(eff)
    if (!is.logical(mask) || anyNA(mask)) {
      stop(sprintf("Condition for outcome '%s' must return a logical array without NA.", name), call. = FALSE)
    }
    if (!identical(dim(mask), ref_dim)) {
      stop(sprintf("Condition for outcome '%s' returned array of shape [%s], expected [%s].",
                   name, paste(dim(mask), collapse = " x "), paste(ref_dim, collapse = " x ")), call. = FALSE)
    }
    masks[[i]] <- mask
  }

  lb <- longbet_py()
  np <- reticulate::import("numpy", convert = TRUE)

  if (!is.null(cells) && (!is.logical(cells) || anyNA(cells) ||
                         !identical(dim(cells), ref_dim[1:2]))) {
    stop("cells must be a logical [N x T] matrix without NA.", call. = FALSE)
  }
  py_masks <- lapply(masks, function(m) reticulate::r_to_py(m))
  py_cells <- if (is.null(cells)) NULL else reticulate::r_to_py(as.matrix(cells))

  py_res <- lb$reduce_joint_masks(py_masks, cells = py_cells)
  as.array(np$asarray(py_res))
}

#' Mean posterior innovation correlation across outcomes
#'
#' @param model A fitted `longbet_multi` object.
#' @return Numeric matrix `[M x M]` with outcome dimnames.
#' @seealso [longbet_multi()]
#' @examples
#' \dontrun{
#' R <- outcome_correlation(fit)
#' }
#' @export
outcome_correlation <- function(model) {
  if (!inherits(model, "longbet_multi")) {
    stop("`model` must be a fitted longbet_multi object.", call. = FALSE)
  }
  py_model <- rehydrate_jax_model(model)
  lb <- longbet_py()
  np <- reticulate::import("numpy", convert = TRUE)

  mat <- as.matrix(np$asarray(lb$outcome_correlation(py_model)))
  dimnames(mat) <- list(model$outcome_names, model$outcome_names)
  mat
}
