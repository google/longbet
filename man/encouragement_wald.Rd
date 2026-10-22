\name{encouragement_wald}
\alias{encouragement_wald}
\alias{encouragement_first_stage}
\alias{encouragement_comparison}
\title{Extract encouragement posterior and comparison tables}
\usage{
encouragement_wald(pred, require_convergence = TRUE)
encouragement_first_stage(pred)
encouragement_comparison(pred)
}
\arguments{
\item{pred}{A longbet_encourage.pred object.}
\item{require_convergence}{Withhold ratio intervals failing the diagnostic
thresholds selected during prediction. Raw draws remain available.}
}
\description{Returns tables computed by the Python engine during prediction,
without a live Python handle.}
\value{A data frame of posterior Wald summaries, first-stage precision and
practical relevance, or descriptive model/reference differences, respectively.}
\details{Frequentist confidence sets are pred$reference, separate from posterior
Wald intervals. A model/reference disagreement flag requires paired sampling
covariance and is not fabricated by summing independent standard errors.}
