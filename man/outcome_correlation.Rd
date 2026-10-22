% Generated from R/multi.R -- regenerate with roxygen2::roxygenise().
\name{outcome_correlation}
\alias{outcome_correlation}
\title{Mean posterior innovation correlation across outcomes}
\usage{
outcome_correlation(model)
}
\arguments{
\item{model}{A fitted \code{longbet_multi} object.}
}
\value{Numeric matrix \code{[M x M]} with outcome dimnames.}
\description{
Returns the mean of per-draw implied innovation correlation matrices
across outcomes in user order.
}
\seealso{\code{\link{longbet_multi}}}
\examples{
\dontrun{
R <- outcome_correlation(fit)
}
}

