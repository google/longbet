\name{validate_encouragement}
\alias{validate_encouragement}
\title{Validate a Single-Wave Encouragement Experiment}
\usage{validate_encouragement(z, d, t = NULL)}
\arguments{
\item{z}{Complete binary absorbing encouragement panel, \code{[N x T]}.}
\item{d}{Complete binary absorbing actual adoption panel, \code{[N x T]}.}
\item{t}{Strictly increasing calendar times with whole-unit gaps; defaults to \code{1:ncol(z)}.}
}
\description{Checks the supported single-wave design with a permanent control arm and returns design metadata.}
\details{Inference assumes complete randomization across individual units. Array validation does not verify randomization, exclusion, monotonicity, or counterfactual timing. Blocked, unequal-probability, and cluster assignment need different estimators. Perfect compliance, adoption before encouragement, negative observed first stages, and experiments without baseline periods are accepted.}
\value{A list with design name, dimensions, arm counts, encouragement time and zero-based column index, \code{has_pre_periods}, and \code{perfect_compliance}.}
\seealso{\code{\link{encouragement_summary}}, \code{\link{encouragement_effects}}}
