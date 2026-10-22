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

#' Predict from a fitted LongBet model
#'
#' Every input the model was fitted with must be supplied again. Omitting one --
#' a propensity score, say -- is an error, not a silently narrower design: the
#' fitted trees index columns by position, so a missing block would make them
#' read the wrong covariate.
#'
#' @param object A fitted `longbet` object.
#' @param x Prognostic covariates, `[N x P]`.
#' @param z Treatment panel, `[N x T]`.
#' @param t Calendar time vector. Defaults to the fitted one.
#' @param x_trt,x_tv,x_trt_tv,ps,baseline_forecast As in [longbet()]; required if used at fit time.
#' @param summary_only If `TRUE`, return per-cell posterior means and credible
#'   bounds without materializing the `[N x T x draws]` arrays. The panel is
#'   processed in blocks of cells, so the quantiles are still exact, and the ATT
#'   is returned in full either way.
#' @param alpha Tail probability for credible intervals.
#' @param random_seed Seed for the Gaussian-process projection used when the
#'   prediction panel runs past the fitted exposure horizon.
#' @param sig_knl,lambda_knl Optionally override the kernel used to **project**
#'   `beta` past the fitted exposure horizon. They do not refit anything: the
#'   posterior draws of `beta` inside the observed window are whatever the model
#'   was fitted with, and only the conditional projection beyond it changes.
#'   That makes a sweep over them a sensitivity analysis of the extrapolation
#'   prior, holding the in-sample estimates fixed -- which is the honest way to
#'   show how much a forecast owes to the prior rather than to the data. It is
#'   *not* the same as refitting with a different kernel, and should not be
#'   described as such. Defaults to the fitted values.
#' @param cache_forest_evaluations Whether to cache prognostic forest evaluations.
#' @param ... Reserved; unsupported arguments are an error.
#' @return An object of class `longbet.pred`. Ordinal predictions also contain
#'   `num_categories`, `categories`, `cutpoints_samples` `[draws x (K-2)]`,
#'   `prob_y`, `prob_mu0`, `prob_tau` `[N x T x K x draws]` (NULL in summary
#'   mode), their `*_summary` lists `[N x T x K]`, and `att_prob_full`
#'   `[S x K x draws]`. Existing prediction fields remain on the latent scale.
#' @export
predict.longbet <- function(object, x, z, t = NULL,
                            x_trt = NULL, x_tv = NULL, x_trt_tv = NULL, ps = NULL,
                            baseline_forecast = NULL,
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
    baseline_forecast = .as_np_forecast(baseline_forecast, "baseline_forecast"),
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

  tau <- summ(py_pred$tau_summary)
  mu0 <- summ(py_pred$mu0_summary)
  yhat <- summ(py_pred$y_summary)

  res <- list(
    tauhats = if (isTRUE(summary_only)) NULL else arr(py_pred$tauhats),
    muhats0 = if (isTRUE(summary_only)) NULL else arr(py_pred$muhats0),
    preds = if (isTRUE(summary_only)) NULL else arr(py_pred$yhats),
    tauhats.mean = tau$mean, tauhats.sd = tau$std,
    tauhats.lower = tau$lower, tauhats.upper = tau$upper,
    muhats0.mean = mu0$mean, muhats0.sd = mu0$std, muhats0.lower = mu0$lower, muhats0.upper = mu0$upper,
    preds.mean = yhat$mean, preds.sd = yhat$std, preds.lower = yhat$lower, preds.upper = yhat$upper,
    att_full = arr(py_pred$att_full),
    beta_values = arr(py_pred$beta_values),
    s = arr(py_pred$s),
    z = z,
    t = t,
    outcome = py_pred$outcome,
    alpha = alpha,
    summary_only = isTRUE(summary_only),
    num_chains = as.integer(py_pred$num_chains),
    att_counts = as.integer(py_pred$att_counts),
    py_pred = py_pred
  )
  res <- c(res, .ordinal_prediction_fields(py_pred))
  class(res) <- "longbet.pred"
  res
}

#' @export
print.longbet.pred <- function(x, ...) {
  d <- dim(x$tauhats.mean)
  cat(sprintf("LongBet predictions for a %d x %d panel\n", d[1], d[2]))
  cat(sprintf("  Exposure times: %d\n", nrow(x$att_full)))
  cat(sprintf("  Draws:          %d\n", ncol(x$att_full)))
  cat(sprintf("  Full draws kept: %s\n", if (x$summary_only) "no (summary_only)" else "yes"))
  if (identical(x$outcome, "binary")) {
    cat("  Scale: latent probit. Apply pnorm() to preds/muhats0 for probabilities.\n")
  }
  if (identical(x$outcome, "ordinal")) {
    cat(sprintf("  Scale: latent probit; category probabilities and effects available for %d categories.\n", x$num_categories))
  }
  invisible(x)
}
