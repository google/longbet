% Generated from R/att.R -- regenerate with roxygen2::roxygenise().
\name{get_att}
\alias{get_att}
\title{Average treatment effect on the treated, by exposure time}
\usage{get_att(object, alpha = 0.05, ...)}
\arguments{
\item{object}{Output of \code{\link{predict.longbet}}.}
\item{alpha}{Tail probability for credible intervals.}
\item{...}{Reserved.}
}
\value{A list with \code{att} (posterior mean by exposure time),
\code{intervals} (a \code{2 x S} matrix), \code{att_full} (\code{S x draws})
and \code{exposure}.}
\description{
Computes the event-time ATT. Delegates to the Python engine rather than
reimplementing the alignment, so the R and Python front doors cannot diverge.
}
