% Generated from R/multi.R -- regenerate with roxygen2::roxygenise().
\name{effect_draws}
\alias{effect_draws}
\title{Extract treatment effect draws on the natural outcome scale}
\usage{
effect_draws(pred, outcome = NULL)
}
\arguments{
\item{pred}{A \code{longbet.pred} or \code{longbet_multi.pred} object.}
\item{outcome}{Outcome name or 1-based index (required for multi predictions).}
}
\value{Numeric 3-D array \code{[N x T x D]}.}
\description{
Returns treatment effect draws on the natural outcome scale: \code{tauhats} for
continuous outcomes (including log scales), and \code{Phi(mu0 + tau) - Phi(mu0)}
for binary outcomes. The latter is a probability difference in [-1, 1], so
\code{0.01} means one percentage point. Saved R prediction arrays can be used
without a live Python prediction handle.
Ordinal selections require an explicit category or score and raise an error;
use \code{\link{att_probabilities}} or \code{\link{att_expected_score}} for ordinal ATTs.
}
\seealso{\code{\link{longbet_multi}}, \code{\link{joint_prob}}}
\examples{
\dontrun{
tau_rev <- effect_draws(pred, outcome = "rev")
tau_churn <- effect_draws(pred, outcome = "churn")
}
}
