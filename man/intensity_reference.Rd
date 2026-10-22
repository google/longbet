\name{intensity_reference}
\alias{intensity_reference}
\alias{print.intensity_reference}
\title{Design-based reference for a continuous-intensity encouragement experiment}
\usage{
intensity_reference(y, d, z, t = NULL, x = NULL, strata = NULL,
  lag_horizon = NULL, demean_pre = FALSE, alpha = 0.05, kappa_draws = NULL)
}
\arguments{
\item{y, d, z}{Complete outcome, dose and offer matrices, \code{[N x T]}. \code{z}
is binary and absorbing with one wave and a permanent holdout; \code{d} is a
finite continuous dose.}
\item{t}{Calendar vector with unit gaps; defaults to \code{1:T}.}
\item{x}{Optional baseline covariates \code{[N x P]}, used by the dip screen only.}
\item{strata}{Optional length-N baseline randomization strata (sites). Contrasts
are stratum-share-weighted differences in arm means.}
\item{lag_horizon}{Number of lags beyond zero in the kernel; NULL means \code{T - t0}.}
\item{demean_pre}{Demean the dose by each unit's strictly pre-offer mean first.}
\item{alpha}{Tail probability for intervals and confidence sets.}
\item{kappa_draws}{Optional matrix of posterior kernel draws \code{[B x (L+1)]}
from a fitted model; when supplied, a \code{discrepancy} table reports
\code{z_{model - T0}} of the cumulative kernel against the reference.}
}
\description{Model-free Tier T0 analysis of a randomized encouragement experiment
whose encouraged behaviour is a continuous, non-absorbing dose: stratum-weighted
offer contrasts on the dose (\code{itt_d}) and the outcome (\code{itt_y}) with the
joint Neyman covariance across periods; the unregularized triangular Toeplitz
deconvolution of the lag kernel \code{kappa} with delta-method standard errors and
the condition number of \code{T_{ITT_D}}; the cumulative kernel; the per-period
Wald ratio and the cumulative return per unit of cumulative dose
(\code{omega_cum}) with exact Fieller/Anderson--Rubin confidence sets; the
arm-wise empirical CDF shift of cumulative dose with a first-order stochastic
dominance check; and the dose-history-conditioned pre-treatment dip screen
\code{z_dip}.}
\details{Nothing here uses within-arm dose variation to learn the kernel, so the
reference is immune to within-arm dose endogeneity and to classical measurement
error in the dose. Under randomization, exclusion, no anticipation, linear
calendar-homogeneous lag effects and time-invariant conditional complier weights,
the population deconvolution identifies the compliance-weighted kernel. Intervals
are pointwise asymptotic approximations, not simultaneous bands or posterior
intervals; the Toeplitz inverse amplifies noise when compliance ramps up
gradually, which \code{condition_number} reports.}
\value{A plain-R list of class \code{intensity_reference} with \code{itt},
\code{covariance}, \code{covariance_index}, \code{kernel}, \code{kernel_covariance},
\code{condition_number}, \code{ecdf}, \code{fosd}, \code{dip},
\code{anderson_rubin_projection}, \code{metadata} and, when \code{kappa_draws} is
given, \code{discrepancy}.}
\seealso{\code{\link{validate_intensity}}, \code{\link{longbet_encourage}}}
