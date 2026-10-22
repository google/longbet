\name{encouragement_bootstrap_comparison}
\alias{encouragement_bootstrap_comparison}
\title{Paired bootstrap comparison of encouragement effects}
\usage{
encouragement_bootstrap_comparison(object, replicates = 200L, groups = NULL,
  random_seed = 0L, alpha = 0.05, min_ess = 400, max_rhat = 1.01)
}
\arguments{
\item{object}{A fitted longbet_encourage object.}
\item{replicates}{Number of bootstrap refits; at least 100 and no failed fits or
diagnostics are needed before any disagreement flag is issued.}
\item{groups}{Optional prespecified baseline subgroup labels.}
\item{random_seed}{Seed for whole-unit resampling and refit seed generation.}
\item{alpha, min_ess, max_rhat}{Interval and diagnostic settings.}
}
\description{Refits model and reference estimators on the same whole-unit
bootstrap samples to retain their paired sampling covariance.}
\value{Plain-R list with table, replicates, failures and metadata. Replicate axes
are replicate, group, horizon, quantity, estimator; quantities are Y,D and
estimators are model,reference. Failures remain missing and are reported.}
\details{This compares two estimators on the same target. It is not an exact
randomization test or a new identification assumption. Repeated model refitting
can be expensive.}
