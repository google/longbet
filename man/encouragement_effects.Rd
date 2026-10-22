\name{encouragement_effects}
\alias{encouragement_effects}
\title{Reference Encouragement Effects and Wald Confidence Sets}
\usage{encouragement_effects(y, d, z, t = NULL, alpha = 0.05)}
\arguments{
\item{y}{Complete finite numeric outcome panel, \code{[N x T]}; binary outcomes yield risk differences.}
\item{d,z,t}{As in \code{\link{validate_encouragement}}.}
\item{alpha}{Significance level; confidence sets have level \code{1-alpha}.}
}
\description{Computes assigned-arm outcome and adoption differences, their joint Neyman covariance estimate, and normal-test-inversion Wald confidence sets. All statistics delegate to Python.}
\details{Assumes individual-unit complete randomization with one encouragement wave and a permanent control arm. Both contrasts use the same units. Missing cells are rejected.

The confidence set analytically inverts the studentized contrast of \code{Y-w*D} using a normal critical value (Fieller/Anderson--Rubin-style inference). This is an asymptotic approximation requiring adequate arm sizes and moment conditions, not an exact permutation procedure or a Bayesian credible interval. Intervals are pointwise across horizons. No arbitrary first-stage cutoff or truncation of infinite endpoints is used.

Randomization identifies the ITTs. A CACE interpretation of their ratio additionally requires exclusion, monotonicity, relevance, and appropriate adoption-history assumptions; these are not verified by the function.}
\value{A data frame with one row per observed post-encouragement period, \code{itt_y}, \code{itt_d}, SEs and bounds, \code{itt_y_d_cov}, and \code{wald}. The confidence set consists of up to two closed components, \code{wald_lower_1}/\code{wald_upper_1} and \code{wald_lower_2}/\code{wald_upper_2}. Infinite endpoints remain infinite; unused components are missing. \code{wald_set_type} is bounded, disjoint, half_line, all_real, singleton, empty, or unavailable. With fewer than two units in an arm, point estimates remain but inference is unavailable and \code{wald_reason} is insufficient_arm_size.}
\references{Aronow, Chang, and Lopatto (2024). Randomization-based Confidence Sets for the Local Average Treatment Effect. \url{https://arxiv.org/abs/2404.18786}.}
