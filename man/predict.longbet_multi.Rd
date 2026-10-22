% Generated from R/multi.R -- regenerate with roxygen2::roxygenise().
\name{predict.longbet_multi}
\alias{predict.longbet_multi}
\title{Predict from a fitted multi-outcome LongBet model}
\usage{
\method{predict}{longbet_multi}(object, x, z, t = NULL, x_trt = NULL,
  x_tv = NULL, x_trt_tv = NULL, ps = NULL, summary_only = FALSE,
  alpha = 0.05, random_seed = NULL, sig_knl = NULL, lambda_knl = NULL,
  cache_forest_evaluations = FALSE, ...)
}
\arguments{
\item{object}{A fitted \code{longbet_multi} object.}
\item{x}{Prognostic covariates, \code{[N x P]}.}
\item{z}{Treatment panel, \code{[N x T]}.}
\item{t}{Calendar time vector. Defaults to fitted time.}
\item{x_trt, x_tv, x_trt_tv, ps}{Design covariates as in \code{\link{longbet_multi}}.}
\item{summary_only}{If \code{TRUE}, return per-cell posterior means and credible
bounds without materializing full \code{[N x T x draws]} arrays.}
\item{alpha}{Tail probability for credible intervals.}
\item{random_seed}{Seed for Gaussian-process projection beyond exposure horizon.}
\item{sig_knl, lambda_knl}{Optional projection kernel overrides.}
\item{cache_forest_evaluations}{Whether to cache prognostic forest evaluations.}
\item{...}{Reserved; unsupported arguments raise an error.}
}
\value{An object of class \code{longbet_multi.pred}. Ordinal children expose
category probability draws, summaries, and category/score ATT methods as
documented in \code{\link{predict.longbet}}.}
\description{
Predict treatment effects and outcomes for all fitted outcomes from a
\code{longbet_multi} model.
}
\seealso{\code{\link{longbet_multi}}, \code{\link{effect_draws}}, \code{\link{joint_prob}}}
\examples{
\dontrun{
pred <- predict(fit, x = x, z = z, t = 1:Tn)
pred_rev <- pred["rev"]
}
}
