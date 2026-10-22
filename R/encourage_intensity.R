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

# Design-based (Tier T0) reference for a continuous, non-absorbing treatment
# intensity. Plain R values only: nothing here retains a reticulate handle.

#' Validate a single-wave encouragement panel with a continuous dose
#'
#' `z` is a complete binary absorbing `[N x T]` offer matrix with one
#' encouragement wave and a permanent control arm; `d` is a complete finite
#' `[N x T]` dose matrix that may rise, vary fractionally and switch off.
#' Periods must be consecutive and unit spaced.
#'
#' @param z Offer matrix, `[N x T]`.
#' @param d Dose matrix, `[N x T]`.
#' @param t Calendar vector with unit gaps; defaults to `1:T`.
#' @return A plain-R list of design metadata (dimensions, arm counts, offer
#'   period and index, dose range, share of zero doses, whether the dose is
#'   non-absorbing). Validation cannot establish randomization, exclusion or
#'   the linearity of the dose response.
#' @export
validate_intensity <- function(z, d, t = NULL) {
  longbet_py()
  lb <- reticulate::import("longbet", convert = FALSE)
  panel <- lb$validate_intensity(z = as.matrix(z), d = as.matrix(d),
                                 t = if (is.null(t)) NULL else as.numeric(t))
  .encouragement_plain(panel$metadata())
}

#' Design-based reference for a continuous-intensity encouragement experiment
#'
#' Model-free Tier T0 analysis: stratum-weighted offer contrasts on the dose
#' (`itt_d`) and the outcome (`itt_y`) with the joint Neyman covariance across
#' periods; the unregularized triangular Toeplitz deconvolution of the lag kernel
#' `kappa` with delta-method standard errors and the condition number of
#' `T_{ITT_D}`; the cumulative kernel; the per-period Wald ratio and the
#' cumulative return per unit of cumulative dose (`omega_cum`) with exact
#' Fieller/Anderson--Rubin confidence sets; the arm-wise empirical CDF shift of
#' cumulative dose with a first-order stochastic dominance check; and the
#' dose-history-conditioned pre-treatment dip screen `z_dip`.
#'
#' Nothing here uses within-arm dose variation to learn the kernel, so the
#' reference is immune to within-arm dose endogeneity and to classical
#' measurement error in the dose. Under randomization, exclusion, no
#' anticipation, linear calendar-homogeneous lag effects and time-invariant
#' conditional complier weights, the population deconvolution identifies the
#' compliance-weighted kernel.
#'
#' @param y,d,z Complete outcome, dose and offer matrices, `[N x T]`.
#' @param t Calendar vector with unit gaps; defaults to `1:T`.
#' @param x Optional baseline covariates `[N x P]`, used by the dip screen only.
#' @param strata Optional length-`N` baseline randomization strata (sites).
#' @param lag_horizon Number of lags beyond zero; NULL means `T - t0`.
#' @param demean_pre Demean the dose by each unit's pre-offer mean first.
#' @param alpha Tail probability for intervals and confidence sets.
#' @param kappa_draws Optional matrix of posterior kernel draws `[B x (L+1)]`
#'   from a fitted model; when supplied, the `discrepancy` table reports
#'   `z_{model - T0}` of the cumulative kernel against the reference.
#' @return A plain-R list with `itt`, `covariance`, `covariance_index`,
#'   `kernel`, `kernel_covariance`, `condition_number`, `ecdf`, `fosd`, `dip`,
#'   `anderson_rubin_projection`, `metadata` and, when `kappa_draws` is given,
#'   `discrepancy`. Intervals are pointwise asymptotic approximations.
#' @export
intensity_reference <- function(y, d, z, t = NULL, x = NULL, strata = NULL, lag_horizon = NULL,
                                demean_pre = FALSE, alpha = 0.05, kappa_draws = NULL) {
  if (!is.null(lag_horizon)) {
    if (!is.numeric(lag_horizon) || length(lag_horizon) != 1L || !is.finite(lag_horizon) ||
        lag_horizon < 0 || lag_horizon != trunc(lag_horizon)) {
      stop("lag_horizon must be a nonnegative integer or NULL.", call. = FALSE)
    }
    lag_horizon <- as.integer(lag_horizon)
  }
  if (is.factor(strata)) strata <- as.character(strata)
  longbet_py()
  lb <- reticulate::import("longbet", convert = FALSE)
  ref <- lb$intensity_reference(
    y = as.matrix(y), d = as.matrix(d), z = as.matrix(z),
    t = if (is.null(t)) NULL else as.numeric(t),
    x = if (is.null(x)) NULL else as.matrix(x),
    strata = if (is.null(strata)) NULL else as.character(strata),
    lag_horizon = lag_horizon, demean_pre = isTRUE(demean_pre), alpha = alpha)
  result <- list(
    itt = .encouragement_plain(ref$itt),
    covariance = .encouragement_plain(ref$covariance$to_numpy()),
    covariance_index = .encouragement_plain(ref$covariance$index$to_frame(index = FALSE)),
    kernel = .encouragement_plain(ref$kernel),
    kernel_covariance = .encouragement_plain(ref$kernel_covariance),
    condition_number = .encouragement_plain(ref$condition_number),
    ecdf = .encouragement_plain(ref$ecdf),
    fosd = .encouragement_plain(ref$fosd),
    dip = .encouragement_plain(ref$dip),
    metadata = .encouragement_plain(ref$metadata)
  )
  result$anderson_rubin_projection <- tryCatch(
    .encouragement_plain(ref$anderson_rubin_projection()), error = function(e) NULL)
  if (!is.null(kappa_draws)) {
    result$discrepancy <- .encouragement_plain(ref$discrepancy(as.matrix(kappa_draws)))
  }
  structure(result, class = "intensity_reference")
}

#' @export
print.intensity_reference <- function(x, ...) {
  post <- x$itt[x$itt$post_encouragement, , drop = FALSE]
  cat(sprintf("Design-based reference (T0): %d units, %d post-offer periods, lag horizon %d\n",
              x$metadata$n_units, nrow(post), as.integer(x$metadata$lag_horizon)))
  cat(sprintf("  cond(T_d) = %.3g; z_dip = %.2f; FOSD p = %.3f\n",
              x$condition_number, x$dip$z_stat, x$fosd$pvalue))
  cat(sprintf("  cumulative return per unit dose at the last horizon: %.3g (%s set)\n",
              post$omega_cum[nrow(post)], post$omega_cum_set_type[nrow(post)]))
  invisible(x)
}
