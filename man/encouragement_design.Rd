\name{encouragement_design}
\alias{encouragement_design}
\title{Declare an encouragement randomization design}
\usage{
encouragement_design(assignment = "complete", blocks = NULL, clusters = NULL,
  target = "unit", probabilities = NULL, block_weights = NULL,
  cohort_weights = NULL, control = "not_yet")
}
\arguments{
\item{assignment}{Complete, Bernoulli, or staggered assignment. Staggered means
independent categorical cohort assignment, not fixed cohort quotas.}
\item{blocks, clusters}{Baseline labels, length N; clusters lie within blocks
and share an encouragement path.}
\item{target}{Unit gives equal participant weights; cluster gives equal included
cluster weights and equal within-cluster subgroup weights.}
\item{probabilities}{Known scalar or vector Bernoulli probabilities. For staggered
assignment a named list maps zero-based cohort start indices and optionally
never to probabilities summing to one per randomization unit. Vectors may have
length N or number of clusters.}
\item{block_weights}{Optional named nonnegative weights summing to one. Every
positive-weight block must contain the subgroup.}
\item{cohort_weights}{For staggered designs, required named positive fixed weights
summing to one over the finite cohorts of interest, using zero-based indices.}
\item{control}{Not_yet uses later cohorts and permanent holdouts; never uses
holdouts. Both require no anticipation.}
}
\description{Specifies the assignment mechanism and fixed analysis target for
design-based encouragement inference.}
\value{A plain-R encouragement_design specification. The Python engine validates
arrays and assignment support when effects are computed.}
