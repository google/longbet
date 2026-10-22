% Generated from R/multi.R -- regenerate with roxygen2::roxygenise().
\name{longbet_multi}
\alias{longbet_multi}
\title{LongBetMulti: Coupled Bayesian ensemble trees for multiple outcomes}
\usage{
longbet_multi(y, x, z, t = NULL, x_trt = NULL, x_tv = NULL, x_trt_tv = NULL,
  ps = NULL, num_sweeps = 1000, num_burnin = 2000, n_skip = 1, num_chains = 4,
  inner_loop_length = NULL, num_trees_pr = 20, num_trees_trt = 60,
  min_points_per_leaf_pr = 10, min_points_per_leaf_trt = 10,
  max_depth_pr = 8, max_depth_trt = 8, alpha_split_pr = 0.95,
  beta_split_pr = 2, alpha_split_trt = 0.25, beta_split_trt = 3,
  num_cutpoints = 100, sig_knl = 1, lambda_knl = 1, kernel_type = "se",
  sigma_m = 1, gp_constant_mean = TRUE, trend_num_features = 64,
  random_intercept = TRUE, gamma_prior_a = 1,
  gamma_prior_b = 0.1, sigma_prior_a = 2, sigma_prior_b = 1,
  outcome = "continuous", outcome_names = NULL,
  standardize = TRUE, random_seed = 0,
  device = c("auto", "cpu", "gpu"), sur = TRUE, sur_prior_var = 1,
  verbose = FALSE, num_categories = NULL, cutpoint_prior_scale = 5, ...)
}
\arguments{
\item{y}{Outcome panels: a list of \code{[N x T]} matrices, or a 3-D numeric array
\code{[N x T x M]}.}
\item{x}{Covariates for the prognostic forest, \code{[N x P]}.}
\item{z}{Binary treatment indicator, \code{[N x T]}. Must be absorbing.}
\item{t}{Calendar time vector of length \code{T}. Defaults to \code{1:ncol(y[[1]])}.}
\item{x_trt}{Covariates for the treatment forest, \code{[N x P_trt]}. Defaults
to \code{x}.}
\item{x_tv}{Time-varying prognostic covariates, \code{[N x T x P]}.}
\item{x_trt_tv}{Time-varying treatment covariates, \code{[N x T x P]}.}
\item{ps}{Propensity score, length \code{N} or \code{[N x T]}.}
\item{num_sweeps}{Posterior draws saved per chain, after burn-in.}
\item{num_burnin}{Burn-in iterations discarded.}
\item{n_skip}{Thinning interval.}
\item{num_chains}{Parallel MCMC chains.}
\item{inner_loop_length}{Sweeps per dispatch.}
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
\item{gp_constant_mean}{Whether the projection reverts to an estimated common level.}
\item{trend_num_features}{Random Fourier features in the unit feature map of
the calendar-time block. More capacity to represent a trend whose covariate
loading is nonlinear, at more cost per sweep. The horseshoe prior on the
block shrinks each column on its own evidence, so raising this is close to
free statistically; the cost is arithmetic.}

\item{random_intercept}{Whether to fit unit random intercepts.}
\item{gamma_prior_a, gamma_prior_b}{Inverse-gamma prior on unit-intercept variance.}
\item{sigma_prior_a, sigma_prior_b}{Inverse-gamma prior on innovation variance,
proper by construction; the default IG(2, 1) has prior mean variance one on the
standardized scale.}
\item{outcome}{Outcome type(s): \code{"continuous"}, \code{"binary"}, \code{"ordinal"}, or a vector of length \code{M}.}
\item{num_categories}{Integer broadcast to ordinal outcomes or a list in user
order (optionally fully named), with \code{NULL} for nonordinal entries.
Every ordinal outcome requires K>=2 and explicit numeric labels \code{0,...,K-1}.
Empty categories are allowed; counts are never inferred.}
\item{cutpoint_prior_scale}{Shared finite positive ordered-normal cutpoint
prior scale in latent probit units, default 5. Inactive for nonordinal outcomes or K=2.}
\item{outcome_names}{Optional character vector of length \code{M} naming the outcomes.}
\item{standardize}{Whether to centre and scale continuous outcomes internally.}
\item{random_seed}{Base PRNG seed.}
\item{device}{\code{"auto"}, \code{"cpu"}, or \code{"gpu"}.}
\item{sur}{Whether to enable full-precision triangular SUR coupling.}
\item{sur_prior_var}{Prior variance for loading regression parameters.}
\item{verbose}{Whether to print progress messages.}
\item{...}{Reserved; unsupported arguments raise an error.}
}
\value{An object of class \code{longbet_multi}.}
\description{
Fits outcome-specific prognostic and treatment forests on a shared panel,
coupling the outcomes through the full precision of a triangular SUR
likelihood (\code{full_precision_sur_v1}). Binary and ordinal equations come
first with zero incoming loadings and unit latent variance. Aligned draws do
not establish calibrated joint inference; check convergence of the effects and
decision events before interpreting their probabilities.
}
\seealso{\code{\link{predict.longbet_multi}}, \code{\link{effect_draws}},
\code{\link{joint_prob}}, \code{\link{outcome_correlation}}}
\examples{
\dontrun{
n <- 100; Tn <- 8
x <- matrix(rnorm(n * 3), n, 3)
z <- matrix(0, n, Tn); z[1:50, 5:Tn] <- 1
y1 <- 0.5 * x[, 1] + 1.2 * z + matrix(rnorm(n * Tn, 0, 0.3), n, Tn)
y2 <- matrix(as.numeric(0.3 * x[, 2] - 0.8 * z + rnorm(n * Tn) > 0), n, Tn)
fit <- longbet_multi(
  y = list(rev = y1, churn = y2),
  x = x, z = z, t = 1:Tn,
  outcome = c(rev = "continuous", churn = "binary"),
  sigma_prior_a = 2, sigma_prior_b = 1,
  num_chains = 2
)
pred <- predict(fit, x = x, z = z, t = 1:Tn, summary_only = TRUE)
}
}
