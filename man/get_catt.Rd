% Generated from R/att.R -- regenerate with roxygen2::roxygenise().
\name{get_catt}
\alias{get_catt}
\title{Conditional average treatment effect on the treated, per cell}
\usage{get_catt(object, alpha = 0.05, ...)}
\arguments{
\item{object}{Output of \code{\link{predict.longbet}}.}
\item{alpha}{Unused; the bounds were fixed when \code{predict()} ran.}
\item{...}{Reserved.}
}
\value{A list with \code{catt}, \code{sd}, \code{lower} and \code{upper}, each
\code{[N x T]}.}
\description{Per-cell posterior summaries of the treatment effect.}
