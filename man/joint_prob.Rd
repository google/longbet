% Generated from R/multi.R -- regenerate with roxygen2::roxygenise().
\name{joint_prob}
\alias{joint_prob}
\title{Empirical joint probability of multiple treatment effect events}
\usage{
joint_prob(pred, conditions, cells = NULL)
}
\arguments{
\item{pred}{A \code{longbet_multi.pred} object with full draws.}
\item{conditions}{Named list or sequence of functions in user outcome order.
Each function receives an \code{[N x T x D]} numeric array of effect draws and
must return an \code{[N x T x D]} logical array.}
\item{cells}{Optional logical panel of shape \code{[N x T]}. When \code{NULL},
returns an \code{[N x T]} matrix of cellwise joint probabilities. When provided,
returns a length-N vector where each unit's probability is averaged over its
selected periods (or \code{NA} if no periods are selected for that unit).}
}
\value{Numeric matrix \code{[N x T]} or length-N vector.}
\description{
Computes the empirical joint probability of multiple treatment effect events
by evaluating user-supplied condition functions on each outcome's effect draws,
taking the conjunction across outcomes, and averaging over draws.
This scalar-effect interface supports continuous and binary outcomes.
Ordinal selections raise an error; use the category/score ATT methods.
}
\seealso{\code{\link{longbet_multi}}, \code{\link{effect_draws}}}
\examples{
\dontrun{
p_joint <- joint_prob(pred, list(
  rev = function(e) e > 0,
  churn = function(e) e < 0
))
}
}
