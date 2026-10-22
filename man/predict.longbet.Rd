% Generated from R/predict.longbet.R -- regenerate with roxygen2::roxygenise().
\name{predict.longbet}
\alias{predict.longbet}
\title{Predict from a fitted LongBet model}
\usage{
\method{predict}{longbet}(object, x, z, t = NULL, x_trt = NULL, x_tv = NULL,
  x_trt_tv = NULL, ps = NULL, baseline_forecast = NULL, summary_only = FALSE,
  alpha = 0.05, random_seed = NULL, sig_knl = NULL, lambda_knl = NULL,
  cache_forest_evaluations = FALSE, ...)
}
\arguments{
\item{object}{A fitted \code{longbet} object.}
\item{x}{Prognostic covariates, \code{[N x P]}.}
\item{z}{Treatment panel, \code{[N x T]}.}
\item{t}{Calendar time vector. Defaults to the fitted one.}
\item{x_trt, x_tv, x_trt_tv, ps, baseline_forecast}{As in \code{\link{longbet}}; required if used
at fit time.}
\item{summary_only}{If \code{TRUE}, return per-cell posterior means and
credible bounds without materializing the \code{[N x T x draws]} arrays. The
panel is processed in blocks of cells, so quantiles are still exact and the
ATT is returned in full either way.}
\item{alpha}{Tail probability for credible intervals.}
\item{random_seed}{Seed for the Gaussian-process projection used when the
prediction panel runs past the fitted exposure horizon.}
\item{sig_knl, lambda_knl}{Optionally override the kernel used to \strong{project}
\code{beta} past the fitted exposure horizon. They do not refit anything: the
posterior draws of \code{beta} inside the observed window are whatever the model
was fitted with, and only the conditional projection beyond it changes. A sweep
over them is therefore a sensitivity analysis of the extrapolation prior holding
the in-sample estimates fixed, not a refit with a different kernel. Defaults to
the fitted values.}
\item{cache_forest_evaluations}{Whether to cache prognostic forest evaluations.}
\item{...}{Reserved; unsupported arguments are an error.}
}
\value{An object of class \code{longbet.pred}. Ordinal predictions also contain
\code{num_categories}, \code{categories}, \code{cutpoints_samples}
\code{[draws x (K-2)]}, \code{prob_y}, \code{prob_mu0}, \code{prob_tau}
\code{[N x T x K x draws]} (NULL in summary mode), their \code{*_summary}
lists \code{[N x T x K]}, and \code{att_prob_full} \code{[S x K x draws]}.
Existing prediction fields remain on the latent scale.}
\description{
Every input the model was fitted with must be supplied again. Omitting one is an
error, not a silently narrower design: the fitted trees index columns by
position, so a missing block would make them read the wrong covariate.
}
