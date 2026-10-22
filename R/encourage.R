# Convert only when reticulate has not already converted the return value.
.encouragement_as_r <- function(value) {
  if (inherits(value, "python.builtin.object")) reticulate::py_to_r(value) else value
}

#' Validate a single-wave encouragement experiment
#'
#' Supports individual-unit complete randomization, a single encouragement wave,
#' and a permanent control arm. Validation checks arrays, not whether assignment
#' was actually randomized or whether exclusion and monotonicity hold. Blocked,
#' unequal-probability, and cluster assignment require different estimators.
#'
#' @param z Binary encouragement indicator, `[N x T]`, absorbing.
#' @param d Binary actual adoption indicator, `[N x T]`, absorbing and complete.
#' @param t Strictly increasing calendar times with whole-unit gaps; defaults
#'   to `1:ncol(z)`. Uneven gaps are allowed.
#' @return A list of dimensions, arm counts, encouragement period and zero-based
#'   column index, `has_pre_periods`, and `perfect_compliance`.
#' @details Perfect compliance, adoption before encouragement, negative sample
#'   first stages, and experiments without baseline periods are accepted.
#' @export
validate_encouragement <- function(z, d, t = NULL) {
  lb <- longbet_py()
  .encouragement_as_r(lb$validate_encouragement(
    z = reticulate::r_to_py(as.matrix(z)),
    d = reticulate::r_to_py(as.matrix(d)),
    t = if (is.null(t)) NULL else reticulate::r_to_py(as.numeric(t))
  ))
}

#' Describe adoption in assigned encouragement arms
#'
#' Includes baseline periods and preserves signed first-stage estimates. All
#' statistics delegate to the Python engine.
#'
#' @param z,d,t As in [validate_encouragement()].
#' @param alpha Significance level; intervals have level `1-alpha`.
#' @return A data frame with one row per calendar period, arm counts and take-up
#'   rates, first-stage differences and uncertainty, and observed adoption lags.
#' @details The standard error uses separate arm sample variances divided by
#'   independent unit counts. Intervals are pointwise normal confidence intervals
#'   under individual-unit complete randomization, not exact or Bayesian intervals.
#'   If an arm has fewer than two units, uncertainty is unavailable.
#'
#'   Lag summaries describe encouraged units observed adopted by each period.
#'   `share_immediate` is the fraction adopting at encouragement;
#'   `median_observed_lag` is adoption time minus encouragement time;
#'   `share_pre_encouragement_adoption` is the fraction adopting before the nudge.
#'   Undefined summaries are missing. None identifies compliance types or verifies
#'   counterfactual timing assumptions. `horizon` uses the model exposure clock;
#'   `period_index` is a zero-based column index.
#' @export
encouragement_summary <- function(z, d, t = NULL, alpha = 0.05) {
  lb <- longbet_py()
  result <- lb$encouragement_summary(
    z = reticulate::r_to_py(as.matrix(z)),
    d = reticulate::r_to_py(as.matrix(d)),
    t = if (is.null(t)) NULL else reticulate::r_to_py(as.numeric(t)),
    alpha = alpha
  )
  result <- as.data.frame(.encouragement_as_r(result))
  # pandas' index attribute carries a live Python handle through reticulate.
  # These results are ordinary tables and must survive saveRDS independently.
  attr(result, "pandas.index") <- NULL
  rownames(result) <- NULL
  result
}

#' Reference encouragement effects and Wald confidence sets
#'
#' Computes assigned-arm differences in means on the outcome and adoption scales,
#' their joint Neyman covariance estimate, and analytic normal-test inversion.
#' Assumes individual-unit complete randomization with one wave and a permanent
#' control arm. Every statistic is computed by the Python engine.
#'
#' @param y Complete finite numeric outcome panel, `[N x T]`. Binary outcomes
#'   yield risk differences without a latent-scale transformation.
#' @param d,z,t As in [validate_encouragement()].
#' @param alpha Significance level; confidence sets have level `1-alpha`.
#' @return A data frame with one row per observed post-encouragement period:
#'   `itt_y`, `itt_d`, their SEs and bounds, `itt_y_d_cov`, `wald`, and a confidence
#'   set in `wald_lower_1`, `wald_upper_1`, `wald_lower_2`, `wald_upper_2`.
#'   `wald_set_type` describes its form. Infinite endpoints remain infinite and
#'   unused components are missing. `wald_reason` identifies unavailable inference.
#' @details The covariance accounts for both outcomes measured on the same units.
#'   The confidence set inverts the studentized assigned-arm contrast of `Y-w*D`
#'   using a normal critical value. This is asymptotic Fieller/Anderson--Rubin-style
#'   inference, not an exact permutation procedure or Bayesian credible interval.
#'   Adequate arm sizes and moment conditions are still required. Intervals are
#'   pointwise, not simultaneous across horizons. No first-stage cutoff is used.
#'
#'   Sets may be bounded, disjoint, half-line, all-real, singleton, empty, or
#'   unavailable. Arms with fewer than two units retain point estimates but have
#'   unavailable uncertainty. Missing data are rejected rather than dropped.
#'
#'   The Wald ratio is not automatically CACE: that interpretation also needs
#'   exclusion, monotonicity, relevance, and appropriate adoption-history
#'   assumptions. Negative first stages are retained. Randomization alone does
#'   not establish those additional assumptions.
#' @export
encouragement_effects <- function(y, d, z, t = NULL, alpha = 0.05) {
  lb <- longbet_py()
  result <- lb$encouragement_effects(
    y = reticulate::r_to_py(as.matrix(y)),
    d = reticulate::r_to_py(as.matrix(d)),
    z = reticulate::r_to_py(as.matrix(z)),
    t = if (is.null(t)) NULL else reticulate::r_to_py(as.numeric(t)),
    alpha = alpha
  )
  result <- as.data.frame(.encouragement_as_r(result))
  attr(result, "pandas.index") <- NULL
  rownames(result) <- NULL
  result
}

#' Plot adoption by assigned encouragement arm
#'
#' Uses the Python-produced summary table; requires the optional ggplot2 package.
#' @param z,d,t As in [validate_encouragement()].
#' @return A ggplot with assigned-arm adoption curves, their gap shaded, and the
#'   encouragement time marked. Baseline periods are included.
#' @export
plot_encouragement <- function(z, d, t = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_encouragement() needs ggplot2; use encouragement_summary() ",
         "to obtain the table without a plotting dependency.", call. = FALSE)
  }
  df <- encouragement_summary(z, d, t)
  start <- df$period[which(df$post_encouragement)[1]]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$period)) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = pmin(.data$takeup_control, .data$takeup_encouraged),
      ymax = pmax(.data$takeup_control, .data$takeup_encouraged)), alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = .data$takeup_encouraged,
                                  colour = "Assigned encouragement")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$takeup_control,
                                  colour = "Assigned control")) +
    ggplot2::geom_vline(xintercept = start, linetype = "dashed", colour = "grey") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Period", y = "Adoption rate", colour = NULL) +
    ggplot2::theme_minimal()
}
