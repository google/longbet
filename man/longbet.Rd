% Generated from R/longbet.R -- regenerate with roxygen2::roxygenise().
\name{longbet}
\alias{longbet}
\title{LongBet: Bayesian ensemble trees for causal inference on longitudinal data}
\usage{
longbet(y, x, z, t = NULL, x_trt = NULL, x_tv = NULL, x_trt_tv = NULL,
  ps = NULL, baseline_forecast = NULL, num_sweeps = 1000, num_burnin = 2000,
  n_skip = 1, num_chains = 4, inner_loop_length = NULL, num_trees_pr = 20,
  num_trees_trt = 60, min_points_per_leaf_pr = 10, min_points_per_leaf_trt = 10,
  max_depth_pr = 8, max_depth_trt = 8, alpha_split_pr = 0.95,
  beta_split_pr = 2, alpha_split_trt = 0.25, beta_split_trt = 3,
  num_cutpoints = 100, sig_knl = 1, lambda_knl = 1, kernel_type = "se",
  sigma_m = 1, gp_constant_mean = TRUE, trend_num_features = 64,
  random_intercept = TRUE, gamma_prior_a = 1,
  gamma_prior_b = 0.1, sigma_prior_a = 2, sigma_prior_b = 1,
  outcome = c("continuous", "binary", "ordinal"),
  standardize = TRUE, random_seed = 0,
  device = c("auto", "cpu", "gpu"), sur = TRUE, sur_prior_var = 1,
  verbose = FALSE, num_categories = NULL, cutpoint_prior_scale = 5, ...)
}
\arguments{
\item{y}{Outcome panel, \code{[N x T]}. \code{NA} marks an unobserved cell;
unbalanced panels are marginalized, not imputed.}
\item{x}{Covariates for the prognostic forest, \code{[N x P]}.}
\item{z}{Binary treatment indicator, \code{[N x T]}. Must be absorbing.}
\item{t}{Calendar time vector of length \code{T}. Defaults to \code{1:ncol(y)}.}
\item{x_trt}{Covariates for the treatment forest, \code{[N x P_trt]}. Defaults
to \code{x}, in which case both forests share one block of columns.}
\item{x_tv}{Time-varying prognostic covariates, \code{[N x T x P]}.}
\item{x_trt_tv}{Time-varying treatment covariates, \code{[N x T x P]}.}
\item{ps}{Propensity score, length \code{N} or \code{[N x T]}. Shown to the
prognostic forest only.}
\item{baseline_forecast}{Optional counterfactual baseline forecast panel,
\code{[N x T]} or \code{[N x T x K]}. Double-centered and appended to the
conjugate calendar-time block while preserving collapsed unit-intercept
updates.}
\item{num_sweeps}{Posterior draws saved per chain, after burn-in.}
\item{num_burnin}{Initial burn-in iterations discarded. This does not
guarantee convergence or sufficient effective sample size.}
\item{n_skip}{Thinning interval; 1 saves every post-burn-in iteration.
At fixed retained draw count, larger values add sampling work while keeping
prediction cost and trace storage approximately fixed.}
\item{num_chains}{Parallel MCMC chains. At least 2 are needed for R-hat.}
\item{inner_loop_length}{Sweeps per dispatch; bounds compile time.}
\item{num_trees_pr, num_trees_trt}{Trees in each forest.}
\item{min_points_per_leaf_pr, min_points_per_leaf_trt}{Minimum leaf sizes.}
\item{max_depth_pr, max_depth_trt}{Maximum tree depths.}
\item{alpha_split_pr, beta_split_pr}{Tree-depth prior, prognostic forest.}
\item{alpha_split_trt, beta_split_trt}{Tree-depth prior, treatment forest.}
\item{num_cutpoints}{Maximum quantile bins per continuous covariate.}
\item{sig_knl}{Kernel standard deviation for the exposure trajectory.}
\item{lambda_knl}{Kernel lengthscale over the exposure index.}
\item{kernel_type}{One of \code{"se"}, \code{"matern32"}, \code{"matern52"},
\code{"ar1"}.}
\item{sigma_m}{Prior standard deviation of the marginalized constant GP mean.}
\item{gp_constant_mean}{Whether the projection reverts to an estimated common
level rather than to zero.}
\item{trend_num_features}{Random Fourier features in the unit feature map of
the calendar-time block. More capacity to represent a trend whose covariate
loading is nonlinear, at more cost per sweep. The horseshoe prior on the
block shrinks each column on its own evidence, so raising this is close to
free statistically; the cost is arithmetic.}

\item{random_intercept}{Whether to fit unit random intercepts.}
\item{gamma_prior_a, gamma_prior_b}{Inverse-gamma prior on the unit-intercept
variance.}
\item{sigma_prior_a, sigma_prior_b}{Inverse-gamma prior on the error variance on
the standardized scale; the default IG(2, 1) is proper with prior mean 1.}
\item{outcome}{\code{"continuous"}, \code{"binary"} (probit), or \code{"ordinal"} (ordered probit).}
\item{num_categories}{Declared integer category count >=2 for ordinal outcomes;
\code{NULL} otherwise. Use explicit numeric labels \code{0,...,K-1}.
Empty categories are allowed; ordinal outcomes are never standardized.}
\item{cutpoint_prior_scale}{Finite positive ordered-normal cutpoint prior scale
in latent probit units, default 5. Inactive for nonordinal outcomes or K=2.}
\item{standardize}{Whether to centre and scale a continuous outcome internally.}
\item{random_seed}{Base PRNG seed.}
\item{device}{\code{"auto"}, \code{"cpu"} or \code{"gpu"}.}
\item{sur}{Whether to enable full-precision triangular SUR coupling (inert for scalar \code{longbet}).}
\item{sur_prior_var}{Prior variance for loading regression parameters.}
\item{verbose}{Whether to print progress messages.}
\item{...}{Reserved. Passing an unsupported argument is an error rather than
being silently ignored.}
}
\value{An object of class \code{longbet}.}
\description{
Fits time-varying heterogeneous treatment effects on panel data with separate
prognostic and treatment forests, a shared Gaussian-process trajectory over the
exposure index, and unit-level random intercepts. The sampler is the JAX engine
in the \code{longbet-jax} Python package.
Fits made before the 2026-09-07 dynamic leaf-precision cache repair must be
refitted from the original data. Loading old posterior draws cannot repair
their uncertainty; archive loaders reject missing sampler provenance.
}
\section{Estimand}{
The reported effect is the contrast between being \code{S} periods into
treatment and not being treated at all,
\deqn{\tau_t(X_i, S) = \beta_S \nu(X_i, t),}
the exposure trajectory times the treatment forest. The treatment term enters
treated cells only, so the untreated surface is the prognostic forest plus the
unit intercept.
}
\section{Sampling defaults}{
The defaults use 2,000 burn-in iterations and 250 retained draws, thinned by
2, across four chains: 2,500 iterations per chain and 1,000 retained draws
overall. These are starting settings, not a convergence guarantee. Check
\code{\link{att_stability}} and the subgroup or probability-scale effects used
in decisions. More burn-in can remove an initial transient, but cannot by itself
repair persistent poor mixing; the retained chains also need adequate
effective sample size.
}

\section{Categorical covariates}{
Every column is treated as ordered numeric; there is no \code{pcat} argument.
One-hot encode unordered categorical variables with more than two levels before
calling this function. This is a real modelling difference from the reference
C++ implementation.
}
\examples{
\dontrun{
n <- 100; Tn <- 8
x <- matrix(rnorm(n * 3), n, 3)
z <- matrix(0, n, Tn); z[1:50, 5:Tn] <- 1
y <- 0.5 * x[, 1] + 1.2 * z + matrix(rnorm(n * Tn, 0, 0.3), n, Tn)
fit <- longbet(y = y, x = x, z = z, t = 1:Tn, num_chains = 2)
pred <- predict(fit, x = x, z = z, t = 1:Tn, summary_only = TRUE)
get_att(pred)
}
}
\seealso{\code{\link{predict.longbet}}, \code{\link{get_att}},
\code{\link{att_stability}}}
