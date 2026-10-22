% Generated from R/att.R -- regenerate with roxygen2::roxygenise().
\name{getMus}
\alias{getMus}
\title{Posterior mean untreated-outcome matrix}
\usage{getMus(object)}
\arguments{\item{object}{Output of \code{\link{predict.longbet}}.}}
\value{An \code{[N x T]} matrix.}
\description{
Posterior mean untreated outcome. Includes the unit random intercept, so it is
a counterfactual outcome rather than a prognostic surface.
}
