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

# Estimands and diagnostics.
#
# These all delegate to the Python engine rather than reimplementing the
# arithmetic in R. That is deliberate: an ATT aligned one way in Python and
# another way in R is exactly the silent drift a two-front-door package is prone
# to, and it cannot happen if there is only one implementation.

#' Average treatment effect on the treated, by exposure time
#'
#' @param object Output of [predict.longbet()].
#' @param alpha Tail probability for credible intervals.
#' @param ... Reserved.
#' @return A list with `att` (posterior mean by exposure time), `intervals`
#'   (a `2 x S` matrix of credible bounds), `att_full` (`S x draws`) and
#'   `exposure` (the exposure times, starting at 1).
#' @export
get_att <- function(object, alpha = 0.05, ...) {
  .reject_unsupported(...)
  if (!inherits(object, "longbet.pred")) {
    stop("get_att() requires the output of predict() on a longbet model.", call. = FALSE)
  }
  np <- reticulate::import("numpy", convert = TRUE)
  res <- object$py_pred$att(alpha = as.numeric(alpha))
  list(
    att = as.numeric(np$asarray(res[["att"]])),
    intervals = as.matrix(np$asarray(res[["intervals"]])),
    att_full = as.matrix(np$asarray(res[["att_full"]])),
    exposure = as.integer(np$asarray(res[["exposure"]]))
  )
}

#' Conditional average treatment effect on the treated, per cell
#'
#' @param object Output of [predict.longbet()].
#' @param alpha Unused; the bounds were fixed when `predict()` ran.
#' @param ... Reserved.
#' @return A list with `catt`, `sd`, `lower` and `upper`, each `[N x T]`.
#' @export
get_catt <- function(object, alpha = 0.05, ...) {
  .reject_unsupported(...)
  if (!inherits(object, "longbet.pred")) {
    stop("get_catt() requires the output of predict() on a longbet model.", call. = FALSE)
  }
  list(
    catt = object$tauhats.mean,
    sd = object$tauhats.sd,
    lower = object$tauhats.lower,
    upper = object$tauhats.upper
  )
}

#' Posterior mean treatment effect matrix
#'
#' @param object Output of [predict.longbet()].
#' @return An `[N x T]` matrix.
#' @export
getTaus <- function(object) {
  if (!inherits(object, "longbet.pred")) {
    stop("getTaus() requires the output of predict() on a longbet model.", call. = FALSE)
  }
  object$tauhats.mean
}

#' Posterior mean untreated-outcome matrix
#'
#' Includes the unit random intercept, so it is a counterfactual outcome rather
#' than a prognostic surface.
#'
#' @param object Output of [predict.longbet()].
#' @return An `[N x T]` matrix.
#' @export
getMus <- function(object) {
  if (!inherits(object, "longbet.pred")) {
    stop("getMus() requires the output of predict() on a longbet model.", call. = FALSE)
  }
  object$muhats0.mean
}

#' Sampling diagnostics for the ATT
#'
#' Reports rank-normalized effective sample size and split R-hat for the ATT at
#' each exposure time. This interface diagnoses the reported ATT. Proper priors
#' still define a parameter posterior despite sign/scale symmetries. Inspect
#' parameter traces when investigating poor mixing, and diagnose every effect
#' used in a decision.
#'
#' @section No reliability verdict:
#' The reference implementation reported a `reliable` column calibrated against
#' *measured coverage of the XBART sweep sampler*. That calibration does not
#' transfer to a Metropolis sampler with real chains, so this function reports
#' the numbers and withholds the verdict: `is_reliable` is `NA`. `ess_ok` and
#' `rhat_ok` are threshold checks on the reported numbers, not evidence about
#' coverage. The split-half interval-width ratio is deliberately not reported --
#' it read 0.93 for a fit that under-covered and 0.96 for one that did not.
#'
#' @param object Output of [predict.longbet()].
#' @param min_ess Threshold for both bulk and tail ESS at every exposure.
#'   Undefined diagnostics fail the check.
#' @param max_rhat R-hat threshold for `rhat_ok`.
#' @param warn Whether to emit warnings when a threshold is not met.
#' @param ... Reserved.
#' @return A list with `summary` (one-row data frame), `by_exposure` (data
#'   frame) and `is_reliable` (always `NA`). Mean MCSE uses mean ESS, not bulk ESS.
#' @export
att_stability <- function(object, min_ess = 400, max_rhat = 1.01, warn = TRUE, ...) {
  .reject_unsupported(...)
  if (!inherits(object, "longbet.pred")) {
    stop("att_stability() requires the output of predict() on a longbet model.", call. = FALSE)
  }
  # Delegated to the prediction object if live, or reconstructed from att_full
  py_pred_ok <- tryCatch(!reticulate::py_is_null_xptr(object$py_pred) &&
                         !is.null(object$py_pred$stability),
                         error = function(e) FALSE)
  if (isTRUE(py_pred_ok)) {
    res <- object$py_pred$stability(
      min_ess = as.numeric(min_ess),
      max_rhat = as.numeric(max_rhat),
      warn = as.logical(warn)
    )
  } else {
    lb <- longbet_py()
    att <- as.matrix(object$att_full)
    # Default to a single chain when the layout is unknown. Guessing 4 would
    # slice one chain into four and manufacture between-chain variance, which
    # is a fabricated R-hat rather than a missing one.
    num_chains <- if (!is.null(object$num_chains)) as.integer(object$num_chains) else 1L
    if (num_chains > 1L && (ncol(att) %% num_chains == 0L)) {
      per <- ncol(att) %/% num_chains
      arr <- aperm(array(t(att), dim = c(per, num_chains, nrow(att))), c(2, 1, 3))
    } else {
      arr <- array(t(att), dim = c(1L, ncol(att), nrow(att)))
    }
    n_treated <- if (!is.null(object$att_counts)) {
      as.integer(object$att_counts)
    } else if (!is.null(object$z) && !is.null(object$s)) {
      vapply(seq_len(nrow(att)), function(si) sum(object$z == 1 & object$s == si), integer(1))
    } else {
      NULL
    }
    res <- lb$att_stability(
      att_draws = reticulate::r_to_py(arr),
      min_ess = as.numeric(min_ess),
      max_rhat = as.numeric(max_rhat),
      warn = as.logical(warn),
      n_treated = if (!is.null(n_treated)) reticulate::r_to_py(n_treated) else NULL
    )
  }
  summary_list <- reticulate::py_to_r(res$summary)
  note <- summary_list$verdict_note
  summary_list$verdict_note <- NULL
  summary_list$reliable <- NULL
  summary_list <- lapply(summary_list, function(v) if (is.null(v)) NA else v)

  by_exp <- reticulate::py_to_r(res$by_exposure)
  list(
    summary = as.data.frame(summary_list, stringsAsFactors = FALSE),
    by_exposure = as.data.frame(lapply(by_exp, as.numeric)),
    is_reliable = NA,
    verdict_note = note
  )
}
