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

# Shared ordinal validation and Python probability interfaces.
.check_cutpoint_prior <- function(scale) {
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("cutpoint_prior_scale must be finite and positive.", call. = FALSE)
  }
}

.check_category_count <- function(K) {
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) ||
      K != floor(K) || K < 2 || K > .Machine$integer.max) {
    stop("num_categories must be an integer >=2.", call. = FALSE)
  }
  as.integer(K)
}

.check_ordinal_labels <- function(y, K) {
  if (is.factor(y) || !is.numeric(y)) {
    stop("Ordinal y must contain explicit numeric category codes 0,...,K-1, not factors.", call. = FALSE)
  }
  if (any(is.infinite(y))) stop("Ordinal y cannot contain infinity; only NA/NaN means missing.", call. = FALSE)
  obs <- y[!is.na(y)]
  if (!length(obs)) stop("Ordinal y contains no observed cells.", call. = FALSE)
  if (any(obs != floor(obs) | obs < 0 | obs >= K)) {
    stop("Ordinal y requires integer labels in [0,K-1].", call. = FALSE)
  }
}

.multi_category_counts <- function(counts, outcomes, outcome_names) {
  ordinal <- outcomes == "ordinal"
  if (!any(ordinal)) {
    if (!is.null(counts)) stop("num_categories requires at least one ordinal outcome.", call. = FALSE)
    return(rep(list(NULL), length(outcomes)))
  }
  if (is.numeric(counts) && length(counts) == 1L) {
    K <- .check_category_count(counts)
    return(lapply(ordinal, function(o) if (o) K else NULL))
  }
  if (!is.list(counts) || length(counts) != length(outcomes)) {
    stop("num_categories must be an integer or a list in outcome order, with NULL for nonordinal outcomes.", call. = FALSE)
  }
  if (!is.null(names(counts))) {
    if (anyNA(names(counts)) || anyDuplicated(names(counts)) ||
        !setequal(names(counts), outcome_names)) {
      stop("num_categories names must exactly match all outcome names.", call. = FALSE)
    }
    counts <- counts[outcome_names]
  }
  for (i in seq_along(outcomes)) {
    if (ordinal[i]) counts[i] <- list(.check_category_count(counts[[i]]))
    else if (!is.null(counts[[i]])) stop("Nonordinal num_categories entries must be NULL.", call. = FALSE)
  }
  unname(counts)
}

.ordinal_prediction_fields <- function(pred) {
  names <- c("num_categories", "categories", "cutpoints_samples", "prob_y", "prob_mu0", "prob_tau",
             "prob_y_summary", "prob_mu0_summary", "prob_tau_summary", "att_prob_full")
  fields <- stats::setNames(rep(list(NULL), length(names)), names)
  if (!identical(as.character(pred$outcome), "ordinal")) return(fields)
  np <- reticulate::import("numpy", convert = TRUE)
  arr <- function(v) as.array(np$asarray(v))
  fields$num_categories <- as.integer(pred$num_categories)
  fields$categories <- as.integer(np$asarray(pred$categories))
  fields$cutpoints_samples <- arr(pred$cutpoints_samples)
  fields$att_prob_full <- arr(pred$att_prob_full)
  for (name in c("prob_y", "prob_mu0", "prob_tau")) {
    if (!isTRUE(as.logical(pred$summary_only))) fields[name] <- list(arr(reticulate::py_get_attr(pred, name)))
    s <- reticulate::py_get_attr(pred, paste0(name, "_summary"))
    fields[[paste0(name, "_summary")]] <- list(mean=arr(s$mean), std=arr(s$std),
                                               lower=arr(s$lower), upper=arr(s$upper))
  }
  fields
}

# Rebuild from ordinary arrays after saveRDS/readRDS invalidates py_pred. This
# invokes the same Python methods; R never duplicates CDF/ATT arithmetic.
.ordinal_prediction_handle <- function(object) {
  if (!inherits(object, "longbet.pred") || !identical(object$outcome, "ordinal")) {
    stop("This method requires an ordinal longbet prediction; select a child of a multi prediction.", call. = FALSE)
  }
  live <- tryCatch(!reticulate::py_is_null_xptr(object$py_pred), error=function(e) FALSE)
  if (isTRUE(live)) return(object$py_pred)
  lb <- longbet_py()
  py <- function(a) if (is.null(a)) NULL else reticulate::r_to_py(a)
  summ <- function(mean, std, lower, upper) lb$PosteriorSummary(py(mean), py(std), py(lower), py(upper))
  ps <- function(s) summ(s$mean, s$std, s$lower, s$upper)
  lb$LongBetPrediction(
    tauhats=py(object$tauhats), muhats0=py(object$muhats0), yhats=py(object$preds),
    att_full=py(object$att_full), beta_values=py(object$beta_values), z=py(object$z), s=py(object$s),
    tau_summary=summ(object$tauhats.mean, object$tauhats.sd, object$tauhats.lower, object$tauhats.upper),
    mu0_summary=summ(object$muhats0.mean, object$muhats0.sd, object$muhats0.lower, object$muhats0.upper),
    y_summary=summ(object$preds.mean, object$preds.sd, object$preds.lower, object$preds.upper),
    outcome="ordinal", num_chains=as.integer(object$num_chains), att_counts=py(object$att_counts),
    summary_only=isTRUE(object$summary_only), num_categories=as.integer(object$num_categories),
    cutpoints_samples=py(object$cutpoints_samples), prob_y=py(object$prob_y),
    prob_mu0=py(object$prob_mu0), prob_tau=py(object$prob_tau),
    prob_y_summary=ps(object$prob_y_summary), prob_mu0_summary=ps(object$prob_mu0_summary),
    prob_tau_summary=ps(object$prob_tau_summary), att_prob_full=py(object$att_prob_full))
}

#' Ordinal category probabilities and probability effects
#'
#' Probabilities integrate unit-variance observation noise. Effects are paired
#' treated-minus-control probability differences within each posterior draw.
#' Rows use fitted unit intercepts by position when unit counts agree; otherwise
#' they use zero intercepts, rather than marginalizing over a new random unit.
#'
#' @param object An ordinal `longbet.pred` object, including a selected multi child.
#' @param arm `"factual"`, `"control"`, or `"effect"`.
#' @param summary Return summaries (`TRUE`) or retained full draws (`FALSE`).
#' @return With `summary=TRUE`, a list of `mean`, `std`, `lower`, `upper`, each
#'   `[N x T x K]`, using the interval level chosen at prediction time. Otherwise
#'   returns `[N x T x K x draws]`; raises if full draws were discarded.
#' @export
predict_probabilities <- function(object, arm="factual", summary=TRUE) {
  pred <- .ordinal_prediction_handle(object)
  res <- pred$predict_probabilities(arm=arm, summary=summary)
  np <- reticulate::import("numpy", convert=TRUE)
  if (!isTRUE(summary)) return(as.array(np$asarray(res)))
  list(mean=as.array(np$asarray(res$mean)), std=as.array(np$asarray(res$std)),
       lower=as.array(np$asarray(res$lower)), upper=as.array(np$asarray(res$upper)))
}

#' Category-specific average treatment effects on the treated
#'
#' @param object An ordinal `longbet.pred`, including a selected multi child.
#' @param alpha Tail probability for credible intervals.
#' @return A list with `att` `[S x K]`, `intervals` `[2 x S x K]`,
#'   `att_full` `[S x K x draws]`, `exposure` and `categories`.
#'   Empty exposure groups contain `NA`/`NaN`.
#' @export
att_probabilities <- function(object, alpha=0.05) {
  res <- .ordinal_prediction_handle(object)$att_probabilities(alpha=alpha)
  np <- reticulate::import("numpy", convert=TRUE)
  list(att=as.array(np$asarray(res[["att"]])), intervals=as.array(np$asarray(res[["intervals"]])),
       att_full=as.array(np$asarray(res[["att_full"]])), exposure=as.integer(np$asarray(res[["exposure"]])),
       categories=as.integer(np$asarray(res[["categories"]])))
}

#' Average treatment effect on an ordinal score
#'
#' Forms weighted category effects within each draw before summarizing. Default
#' rank scores are a reporting convention, not an assumption that ordinal levels
#' have equal spacing. Indicator scores can request an exceedance probability.
#'
#' @param object An ordinal `longbet.pred`, including a selected multi child.
#' @param weights Finite numeric vector of length K; defaults to `0:(K-1)`.
#'   Weights need not be increasing.
#' @param alpha Tail probability for credible intervals.
#' @return A list with `att` `[S]`, `intervals` `[2 x S]`, `att_full` `[S x draws]`,
#'   `exposure` and the chosen `weights`. Works in summary mode.
#' @export
att_expected_score <- function(object, weights=NULL, alpha=0.05) {
  res <- .ordinal_prediction_handle(object)$att_expected_score(
    weights=if (is.null(weights)) NULL else reticulate::r_to_py(weights), alpha=alpha)
  np <- reticulate::import("numpy", convert=TRUE)
  list(att=as.numeric(np$asarray(res[["att"]])), intervals=as.array(np$asarray(res[["intervals"]])),
       att_full=as.array(np$asarray(res[["att_full"]])), exposure=as.integer(np$asarray(res[["exposure"]])),
       weights=as.numeric(np$asarray(res[["weights"]])))
}
