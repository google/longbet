\name{longbet_encourage}
\alias{longbet_encourage}
\title{The encouragement model: adoption clock or continuous treatment intensity}
\usage{
longbet_encourage(y, d, z, x, t = NULL, x_trt = NULL,
  outcome = "continuous", treatment = "absorbing", identification = NULL,
  lag_horizon = NULL, demean_pre = FALSE, stage2_backend = "longbet",
  strata = NULL, config = list(), verbose = FALSE)
}
\arguments{
\item{y}{Complete numeric outcome matrix, N by T.}
\item{d, z}{Actual adoption and assigned encouragement, complete N by T matrices.
\code{z} is binary and absorbing with a single wave and a permanent holdout;
\code{d} is binary absorbing for \code{treatment = "absorbing"} and a finite
continuous dose for \code{treatment = "intensity"}.}
\item{x}{Baseline covariate matrix, N by P.}
\item{t}{Increasing calendar vector with whole-unit gaps; defaults to 1:T.}
\item{x_trt}{Optional baseline treatment-forest covariates.}
\item{outcome}{\code{"continuous"} or \code{"binary"} for y; adoption is always a
binary hazard. The intensity model supports a continuous outcome only.}
\item{treatment}{\code{"absorbing"} (adoption clock) or \code{"intensity"}
(distributed-lag kernel for a continuous dose).}
\item{identification}{Identification tier of the intensity model:
\code{"deconvolution"} (T0), \code{"control_function"} (T1, the default when
NULL) or \code{"unconfounded"} (T2). Must be NULL for \code{treatment = "absorbing"}.}
\item{lag_horizon}{Number of lags beyond zero in the kernel; NULL means
\code{T - t0}, one lag per post-offer period (intensity only).}
\item{demean_pre}{Demean the dose by each unit's pre-offer mean before building
lags (intensity only).}
\item{stage2_backend}{\code{"longbet"} (the sampler) or \code{"linear"} (a
closed-form two-way fixed-effects regression without covariate moderation, for
smoke tests and fallbacks; intensity only).}
\item{strata}{Optional length-N baseline randomization strata (sites) used by
the design-based reference and the dip screen of the intensity model.}
\item{config}{Named list of Python LongBetConfig sampler options applied to both
equations, such as num_chains, num_sweeps, num_burnin and random_seed. The
model's components (linear moderator term, learned trajectory scale, AR(1)
errors) cannot be switched off. The lag kernel's prior options
(\code{lag_sig_knl}, \code{lag_lambda_knl}, ...) are set here too.}
\item{verbose}{Whether to print a fitting message.}
}
\description{With \code{treatment = "absorbing"} this fits two LongBet equations
and composes them: the outcome on the observed adoption clock (periods since
adoption, both arms contributing adoption events) and adoption as a
population-averaged discrete-time hazard on the offer clock. With
\code{treatment = "intensity"} the encouraged behaviour is a continuous,
non-absorbing dose: Stage 1 is a continuous LongBet for the dose on the offer
clock and Stage 2 the distributed-lag kernel for the outcome,
\code{nu_Y(x) sum_l kappa_l D_{i,t-l}} with a Gaussian-process prior over lags.
The offer's effects follow by composition in both cases.}
\details{Both moderators carry a horseshoe linear term beside the trees, the
trajectory's scale is learned, and a continuous outcome has AR(1) errors; draws
are thinned by 2 unless config sets n_skip. The adoption-clock model uses the
exclusion restriction and assumes adoption timing is not confounded with the
outcome innovations given covariates and the outcome's unit intercept.
For the intensity model, \code{identification} chooses what Stage 2 learns the
kernel from: the Stage-1 fitted offer-induced dose only (T0, immune to within-arm
dose endogeneity and measurement error), the observed dose plus the three-term
plug-in control function of the Stage-1 residuals (T1), or the observed dose
alone (T2). The design-based reference returned by \code{predict()} needs none of
these assumptions and is the check on the population offer effects. Both ITTs
average over all original study units with equal weights and each unit's own
fitted intercepts. The Wald ratio needs exclusion, monotonicity, relevance and
appropriate adoption-history assumptions for a CACE interpretation; the
intensity model's summary ratio is the cumulative return per unit of
offer-induced cumulative dose (\code{omega_cum}).}
\value{A longbet_encourage object containing a versioned raw archive, configuration
and metadata. The archive includes both fitted equations and the original
outcomes, adoption, assignment, covariates and strata. No live Python handles are
stored; the object survives saveRDS and readRDS in a fresh R session.}
\seealso{\code{\link{predict.longbet_encourage}}, \code{\link{predict_conditional}},
\code{\link{encouragement_effects}}, \code{\link{intensity_reference}}}
