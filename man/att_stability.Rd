% Generated from R/att.R -- regenerate with roxygen2::roxygenise().
\name{att_stability}
\alias{att_stability}
\title{Sampling diagnostics for the ATT}
\usage{att_stability(object, min_ess = 400, max_rhat = 1.01, warn = TRUE, ...)}
\arguments{
\item{object}{Output of \code{\link{predict.longbet}}.}
\item{min_ess}{Threshold for both bulk and tail ESS at every exposure.
Undefined diagnostics fail the check.}
\item{max_rhat}{R-hat threshold for \code{rhat_ok}.}
\item{warn}{Whether to emit warnings when a threshold is not met.}
\item{...}{Reserved.}
}
\value{A list with \code{summary} (one-row data frame), \code{by_exposure}
(data frame) and \code{is_reliable} (always \code{NA}). Mean MCSE uses mean ESS,
not bulk ESS.}
\description{
Reports rank-normalized effective sample size and split R-hat for the ATT at
each exposure time. This interface diagnoses the reported ATT. Proper priors
still define a parameter posterior despite sign/scale symmetries. Inspect
parameter traces when investigating poor mixing, and diagnose every effect
used in a decision.
}
\section{No reliability verdict}{
The reference implementation reported a \code{reliable} column calibrated
against measured coverage of the XBART sweep sampler. That calibration does not
transfer to a Metropolis sampler with real chains, so this function reports the
numbers and withholds the verdict: \code{is_reliable} is \code{NA}. The
split-half interval-width ratio is deliberately not reported -- it read 0.93 for
a fit that under-covered and 0.96 for one that did not.
}
