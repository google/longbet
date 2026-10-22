\name{encouragement_summary}
\alias{encouragement_summary}
\title{Describe Adoption in Assigned Encouragement Arms}
\usage{encouragement_summary(z, d, t = NULL, alpha = 0.05)}
\arguments{
\item{z,d,t}{As in \code{\link{validate_encouragement}}.}
\item{alpha}{Significance level; intervals have level \code{1-alpha}.}
}
\description{Describes assigned-arm adoption rates and observed lags, including baseline periods. All statistics delegate to the Python engine.}
\details{Inference assumes individual-unit complete randomization. Standard errors use separate arm sample variances divided by unit counts. Intervals are pointwise normal confidence intervals, not exact or Bayesian intervals. Uncertainty is unavailable with fewer than two units in either arm. Negative first stages are retained.

Lag columns describe encouraged units observed adopted by each period. \code{share_immediate} is the fraction adopting at encouragement, \code{median_observed_lag} is adoption time minus encouragement time, and \code{share_pre_encouragement_adoption} is the fraction adopting earlier. Undefined fractions are missing. These summaries cannot identify compliers or verify assumptions about both potential adoption times.}
\value{A data frame with one row per period: calendar time, zero-based \code{period_index}, model exposure \code{horizon}, post-encouragement indicator, unit counts, adoption rates, first stage and uncertainty, and observed lag summaries.}
\seealso{\code{\link{plot_encouragement}}, \code{\link{encouragement_effects}}}
