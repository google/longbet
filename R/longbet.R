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

#' LongBet: Bayesian ensemble trees for causal inference on longitudinal data
#'
#' Fits time-varying heterogeneous treatment effects on panel data with separate
#' prognostic and treatment forests, a conjugate calendar-time block, a shared
#' Gaussian-process trajectory over the exposure index, and unit-level random
#' intercepts. The sampler is the JAX engine in the `longbet-jax` Python
#' package; this function marshals data in and results out.
#'
#' @section Estimand:
#' The reported effect is the contrast between being `S` periods into treatment
#' and not being treated at all,
#' \deqn{\tau_t(X_i, S) = \beta_S \nu(X_i),}
#' the exposure trajectory times the treatment forest. The treatment term enters
#' treated cells only, so the untreated surface is the prognostic forest plus the
#' calendar-time block and the unit intercept.
#'
#' @section Sampling defaults:
#' The defaults use 2,000 burn-in iterations and 1,000 retained draws, thinned by
#' 1, across four chains: 3,000 iterations per chain and 4,000 retained draws
#' overall. Check [att_stability()] and the subgroup or probability-scale
#' effects used in decisions.
#'
#' @section Categorical covariates:
#' Every column is treated as ordered numeric; there is no `pcat` argument.
#' One-hot encode unordered categorical variables with more than two levels
#' before calling this function. This is a real modelling difference from the
#' reference C++ implementation, not a detail.
#'
#' @param y Outcome panel, `[N x T]`. `NA` marks an unobserved cell; unbalanced
#'   panels are marginalized, not imputed.
#' @param x Covariates for the prognostic forest, `[N x P]`.
#' @param z Binary treatment indicator, `[N x T]`. Must be absorbing.
#' @param t Calendar time vector of length `T`. Defaults to `1:ncol(y)`.
#' @param x_trt Covariates for the treatment forest, `[N x P_trt]`. Defaults to
#'   `x`, in which case both forests share one block of columns.
#' @param x_tv Time-varying prognostic covariates, `[N x T x P]`.
#' @param x_trt_tv Time-varying treatment covariates, `[N x T x P]`.
#' @param ps Propensity score, length `N` or `[N x T]`. Shown to the prognostic
#'   forest only. For staggered adoption the principled analogue is a timing or
#'   hazard score rather than one scalar per unit.
#' @param baseline_forecast Optional counterfactual baseline forecast panel,
#'   `[N x T]` or `[N x T x K]`. Double-centered and appended to the conjugate
#'   calendar-time block while preserving collapsed unit-intercept updates.
#' @param num_sweeps Posterior draws saved per chain, after burn-in.
#' @param num_burnin Initial burn-in iterations discarded. This does not
#'   guarantee convergence or sufficient effective sample size.
#' @param n_skip Thinning interval; 1 saves every post-burn-in iteration.
#'   At fixed retained draw count, larger values add sampling work while keeping
#'   prediction cost and trace storage approximately fixed.
#' @param num_chains Parallel MCMC chains. At least 2 are needed for R-hat.
#' @param inner_loop_length Sweeps per dispatch; bounds compile time.
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
#' @param gp_constant_mean Whether the projection reverts to an estimated common
#'   level rather than to zero.
#' @param trend_num_features Random Fourier features in the unit feature map of
#'   the calendar-time block. More capacity to represent a trend whose covariate
#'   loading is nonlinear, at more cost per sweep. The horseshoe prior on the
#'   block shrinks each column on its own evidence, so raising this is close to
#'   free statistically; the cost is arithmetic.
#' @param random_intercept Whether to fit unit random intercepts.
#' @param gamma_prior_a,gamma_prior_b Inverse-gamma prior on the unit-intercept
#'   variance.
#' @param sigma_prior_a,sigma_prior_b Inverse-gamma prior on the error variance
#'   on the standardized scale; the default IG(2, 1) is proper with prior mean 1.
#' @param outcome `"continuous"`, `"binary"` (probit), or `"ordinal"` (ordered probit).
#' @param num_categories Declared integer category count, at least 2, for ordinal
#'   outcomes; `NULL` otherwise. Labels must be numeric `0,...,K-1`; empty
#'   categories are allowed. Ordinal outcomes are never standardized.
#' @param cutpoint_prior_scale Finite positive ordered-normal cutpoint prior
#'   scale in latent probit units, default 5. Inactive for nonordinal or K=2.
#' @param standardize Whether to centre and scale a continuous outcome internally.
#' @param random_seed Base PRNG seed.
#' @param device `"auto"`, `"cpu"` or `"gpu"`.
#' @param sur Whether to enable full-precision triangular SUR coupling (inert for scalar `longbet`).
#' @param sur_prior_var Prior variance for loading regression parameters.
#' @param verbose Whether to print progress messages.
#' @param ... Reserved. Passing an unsupported argument is an error rather than
#'   being silently ignored.
#' @return An object of class `longbet`.
#' @seealso [predict.longbet()], [get_att()], [att_stability()]
#' @export
longbet <- function(y, x, z, t = NULL,
                    x_trt = NULL, x_tv = NULL, x_trt_tv = NULL, ps = NULL,
                    baseline_forecast = NULL,
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
                    outcome = c("continuous", "binary", "ordinal"),
                    standardize = TRUE, random_seed = 0,
                    device = c("auto", "cpu", "gpu"),
                    sur = TRUE, sur_prior_var = 1.0,
                    verbose = FALSE, num_categories = NULL,
                    cutpoint_prior_scale = 5.0, ...) {

  .reject_unsupported(...)
  outcome <- match.arg(outcome)
  device <- match.arg(device)

  .check_cutpoint_prior(cutpoint_prior_scale)
  if (outcome == "ordinal") num_categories <- .check_category_count(num_categories)
  else if (!is.null(num_categories)) stop("num_categories must be NULL for nonordinal outcomes.", call. = FALSE)

  y <- as.matrix(y)
  if (outcome == "ordinal") .check_ordinal_labels(y, num_categories)
  x <- as.matrix(x)
  z <- as.matrix(z)
  if (is.null(t)) t <- seq_len(ncol(y))
  t <- as.numeric(t)

  if (!identical(dim(y), dim(z))) {
    stop(sprintf("y is %d x %d but z is %d x %d; they must match.",
                 nrow(y), ncol(y), nrow(z), ncol(z)), call. = FALSE)
  }
  if (nrow(x) != nrow(y)) {
    stop(sprintf("x has %d rows but y has %d.", nrow(x), nrow(y)), call. = FALSE)
  }

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
    outcome = outcome,
    num_categories = num_categories,
    cutpoint_prior_scale = as.numeric(cutpoint_prior_scale),
    standardize = as.logical(standardize),
    random_seed = as.integer(random_seed),
    device = device,
    sur = as.logical(sur),
    sur_prior_var = as.numeric(sur_prior_var)
  )

  if (isTRUE(verbose)) {
    message(sprintf("longbet: fitting %d x %d panel on device '%s' (%d chain%s).",
                    nrow(y), ncol(y), device, num_chains,
                    if (num_chains == 1) "" else "s"))
  }

  py_model <- lb$LongBet(config)
  py_model$fit(
    y = .as_np_matrix(y, "y"),
    x = .as_np_matrix(x, "x"),
    z = .as_np_matrix(z, "z"),
    t = .as_np_vector(t),
    x_trt = .as_np_matrix(x_trt, "x_trt"),
    x_tv = .as_np_array3(x_tv, "x_tv"),
    x_trt_tv = .as_np_array3(x_trt_tv, "x_trt_tv"),
    ps = if (is.null(ps)) NULL else reticulate::r_to_py(if (is.matrix(ps)) as.matrix(ps) else as.numeric(ps)),
    baseline_forecast = .as_np_forecast(baseline_forecast, "baseline_forecast")
  )

  # Serialize the engine state to a raw vector: the R object must hold plain R
  # data, never a live Python handle, or saveRDS() produces something that
  # cannot be reopened.
  tmp <- tempfile(fileext = ".npz")
  on.exit(unlink(tmp), add = TRUE)
  py_model$save(tmp)
  raw_model <- readBin(tmp, what = "raw", n = file.info(tmp)$size)

  trace <- py_model$trace
  as_arr <- function(v) as.array(reticulate::py_to_r(reticulate::import("numpy")$asarray(v)))

  # Draws go last in R and first in JAX; permute once, here.
  flatten_draws <- function(a) {
    if (length(dim(a)) == 3L) {           # (chains, draws, k) -> (k, chains*draws)
      matrix(aperm(a, c(3, 2, 1)), nrow = dim(a)[3])
    } else {                              # (draws, k) -> (k, draws)
      t(a)
    }
  }

  # (chains, draws) -> (chains*draws,) chain-major. In R, as.numeric() traverses
  # column-major, so transpose first to avoid interleaving chains.
  flatten_scalars <- function(a) {
    if (length(dim(a)) == 2L) {
      as.numeric(t(a))
    } else {
      as.numeric(a)
    }
  }

  beta_draws <- as_arr(trace$beta)
  gamma_draws <- as_arr(trace$gamma)
  sigma2_draws <- flatten_scalars(as_arr(trace$sigma2))

  # The live handle is cached in an environment, but the model object also
  # carries `raw_model`, which is what predict() falls back to after a
  # saveRDS/readRDS round trip leaves the cached pointer dangling.
  handle <- new.env(parent = emptyenv())
  handle$py_model <- py_model
  handle$pid <- Sys.getpid()

  fit <- list(
    beta_values = flatten_draws(beta_draws),
    gamma_draws = flatten_draws(gamma_draws),
    sigma_gamma_draws = flatten_scalars(as_arr(trace$sigma_gamma2)),
    # Standardized units, as the reference package reports it: multiply by sdy
    # to recover the residual standard deviation on the data scale.
    sigma0_draws = matrix(sqrt(sigma2_draws), nrow = 1),
    b0_draws = flatten_scalars(as_arr(trace$b0)),
    b1_draws = flatten_scalars(as_arr(trace$b1)),
    alpha_draws = flatten_scalars(as_arr(trace$alpha)),
    sdy = as.numeric(py_model$sdy),
    meany = as.numeric(py_model$meany),
    random_intercept = as.logical(random_intercept),
    outcome = outcome,
    num_categories = num_categories,
    cutpoint_prior_scale = as.numeric(cutpoint_prior_scale),
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
  class(fit) <- "longbet"
  fit
}

#' @export
print.longbet <- function(x, ...) {
  p <- x$model_params
  cat("LongBet fit (JAX engine)\n")
  cat(sprintf("  Draws:   %d saved, %d burn-in, %d chain%s\n",
              p$num_sweeps, p$burnin, p$num_chains,
              if (p$num_chains == 1) "" else "s"))
  cat(sprintf("  Trees:   %d prognostic / %d treatment\n",
              p$num_trees_pr, p$num_trees_trt))
  cat(sprintf("  Outcome: %s\n", x$outcome))
  if (identical(x$outcome, "ordinal")) cat(sprintf("  Categories: %d; effects and residual SD are on the latent probit scale.\n", x$num_categories))
  cat(sprintf("  Unit random intercept: %s\n",
              if (isTRUE(x$random_intercept)) "yes" else "no"))
  cat(sprintf("  Residual SD (posterior mean): %.4f\n",
              mean(x$sigma0_draws) * x$sdy))
  if (p$num_chains < 2) {
    cat("  Note: R-hat needs at least two chains; set num_chains >= 2.\n")
  }
  invisible(x)
}

#' @export
summary.longbet <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

#' Inspect available JAX compute devices
#'
#' @return Character vector such as `c("cpu:0", "gpu:0")`.
#' @export
longbet_devices <- function() {
  as.character(reticulate::py_to_r(longbet_py()$available_devices()))
}
