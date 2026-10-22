% Generated from R/rollout.R -- regenerate with roxygen2::roxygenise().
\name{rollout_summary}
\alias{rollout_summary}
\title{Summarise a staggered rollout}
\usage{rollout_summary(z, t = NULL, labels = NULL,
  never_treated_label = "Never treated")}
\arguments{
\item{z}{Binary treatment indicator, \code{[N x T]}. Must be absorbing.}
\item{t}{Calendar time vector, length \code{T}. Defaults to \code{1:ncol(z)}.}
\item{labels}{Optional named vector of display names keyed by adoption time,
e.g. \code{c("4" = "W1")}. Unlabelled cohorts fall back to their adoption time.}
\item{never_treated_label}{Name for the cohort that is never treated.}
}
\value{A data frame with columns \code{cohort}, \code{first_treated},
\code{n_units}, \code{period}, \code{period_index}, \code{status} and
\code{exposure}.}
\description{
Describes who was treated and when, as one row per adoption cohort and period. A
cohort is the set of units sharing a first-treated period, derived from \code{z}
rather than declared, so this works whether or not the design has named launch
waves. It needs no fitted model, and the \code{exposure} column is computed by
the same code the sampler uses.
}
\seealso{\code{\link{plot_rollout}}}
