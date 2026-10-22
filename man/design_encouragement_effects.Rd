\name{design_encouragement_effects}
\alias{design_encouragement_effects}
\title{Encouragement inference under an explicit randomization design}
\usage{
design_encouragement_effects(y, d, z, t = NULL, design = NULL,
  groups = NULL, alpha = 0.05)
}
\arguments{
\item{y, d, z}{Complete outcome, adoption and encouragement panels, N by T.}
\item{t}{Increasing calendar times with whole-unit gaps; defaults to 1:T.}
\item{design}{An encouragement_design specification; NULL selects individual-unit
complete randomization and equal unit weights.}
\item{groups}{Optional prespecified baseline subgroup labels, length N.}
\item{alpha}{Tail probability for pointwise normal intervals and AR-style sets.}
}
\description{Computes fixed-population Horvitz--Thompson ITTs and their joint
randomization-unit covariance in Python.}
\value{A plain-R list with table, numeric covariance matrix, covariance_index
data frame and metadata. Both covariance axes follow covariance_index rows,
with columns contrast_id and quantity (itt_y or itt_d).}
\details{Preserves shared controls, cross-horizon dependence and unsupported
groups. Blocked complete designs use Neyman covariance. Bernoulli and categorical
designs use an uncentered conservative covariance bound that can be loose and
depends on outcome location. Subgroup HT estimates do not divide by random
subgroup arm sizes. Missing cohort support makes fixed-weight aggregates
unavailable without renormalization. Intervals are pointwise normal approximations,
not exact permutation tests, simultaneous bands or posterior intervals.}
