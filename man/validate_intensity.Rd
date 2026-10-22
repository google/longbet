\name{validate_intensity}
\alias{validate_intensity}
\title{Validate a Single-Wave Encouragement Panel with a Continuous Dose}
\usage{validate_intensity(z, d, t = NULL)}
\arguments{
\item{z}{Complete binary absorbing offer panel, \code{[N x T]}, with one
encouragement wave and a permanent control arm.}
\item{d}{Complete finite dose panel, \code{[N x T]}; it may rise, vary
fractionally and switch off.}
\item{t}{Calendar times with unit gaps; defaults to \code{1:ncol(z)}.}
}
\description{Checks the supported single-wave design for a continuous,
non-absorbing treatment intensity and returns design metadata.}
\details{Array validation does not verify randomization, exclusion, no
anticipation or the linearity of the dose response. Periods must be consecutive
and unit spaced because the distributed-lag kernel is indexed by lag.}
\value{A list with the design name, \code{treatment = "intensity"}, dimensions,
arm counts, encouragement time and zero-based column index,
\code{has_pre_periods}, the dose range, the share of zero doses and whether the
dose is non-absorbing.}
\seealso{\code{\link{intensity_reference}}, \code{\link{longbet_encourage}}}
