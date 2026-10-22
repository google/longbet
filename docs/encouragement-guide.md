# Randomized encouragement in LongBet

This guide covers the encouragement tools of the package: the design-based
reference analysis, the adoption-clock model, its unit-level predictions and
decisions, model/reference comparison, and the generalized-design reference
estimators. Matching Python and R calls are shown throughout.

## Data and interpretation

Use complete `(N, T)` panels for outcomes `y`, encouragement `z`, and actual
adoption `d`. Both indicators are binary and absorbing. Calendar time `t` has
strictly positive whole-unit gaps. Baseline covariates `x` have shape `(N, P)`.
The model supports a single randomized encouragement wave and a permanent
holdout with complete individual assignment; adoption can precede the offer
or occur later. Experiments without baseline periods and perfect compliance
are valid.

Randomization identifies the offer's intent-to-treat effects on the outcome
and on adoption. Interpreting their ratio as a complier effect also requires
exclusion, monotonicity, relevance, and appropriate assumptions about adoption
history.

## Single-wave reference analysis

```python
from longbet import encouragement_effects, encouragement_summary

adoption = encouragement_summary(z, d, t)
reference = encouragement_effects(y, d, z, t)
```

```r
adoption <- encouragement_summary(z, d, t)
reference <- encouragement_effects(y, d, z, t)
```

These are assigned-arm differences in means with pointwise normal inference
and Fieller/Anderson-Rubin-style confidence sets for the ratio; the sets can be
bounded, unbounded, disjoint, all real numbers, or empty, and all of that is
preserved. They need no model and remain the inference of record for the
population offer effects.

## The adoption-clock model

`LongBetEncourage` fits two LongBet equations and composes them.

1. **Outcome on the adoption clock.** `y` on the observed adoption `d`, with
   LongBet's exposure clock counting periods since adoption. The equation has
   the adoption jump, its exposure profile (`beta_S`), its covariate
   heterogeneity (`nu(x)`), calendar-time prognostic structure, a unit
   intercept and, for a continuous outcome, AR(1) errors. Both arms
   contribute adoption events.
2. **Adoption as a discrete-time hazard.** A binary LongBet on the
   first-adoption indicator over the at-risk set (cells after adoption are
   masked), with the offer as the treatment and periods since the offer as its
   clock. The hazard is population-averaged: it carries no unit frailty.

For a unit with covariates `x`, the hazard equation gives the distribution of
its adoption period with and without the offer; the outcome equation gives the
expected outcome response `s` periods after adopting; the offer effect on the
outcome is the difference of the expected exposure responses, and the offer
effect on adoption is the difference in adoption probabilities. The fitted
hazards of a population-averaged hazard are survivor averages, and their
survival product is the marginal adoption-time distribution, which is all the
composition needs.

In both equations the moderator is a horseshoe linear term plus trees,
`nu(x) = x'delta + sum of trees` on standardized baseline covariates, and the
exposure profile's scale is learned (`beta ~ GP(0, s^2 K)`, `s ~ C+(0, 1)`).
These components are part of the model (`MODEL_COMPONENTS`), chosen on
development simulations and confirmed on held-out replicates; a configuration
that switches one off is refused.

| Component | Why |
| --- | --- |
| horseshoe linear moderator term (`linear_trt`) | the tree prior alone makes realistic heterogeneity a priori negligible |
| population-averaged hazard (no unit frailty) | with one adoption event per unit a frailty's variance is weakly identified, mixes poorly and flattens the offer's effect on adoption |
| AR(1) errors for a continuous outcome (`ar1_errors`) | independent errors overstate the information in a serially correlated series; needs the unit intercepts and unit-level forest partitions |
| learned trajectory scale (`learn_beta_scale`) | the effect can shrink under a null and grow under a strong signal |
| interweaving of the intercept variance (`interweave_gamma`) | mixing; same posterior |

Draws are thinned by 2 unless you set `n_skip`, the thinning of every
validated fit. `LongBetEncourage.model_config()` shows the configuration a fit
will use. A binary outcome's unit intercept is integrated out exactly,
`E[Phi(a + sigma Z)] = Phi(a / sqrt(1 + sigma^2))`. A plain `LongBet` fit is
unaffected and remains the published model. Archives of the earlier release
of the encouragement model are refused.

Assumptions beyond randomization: exclusion (the offer moves outcomes only
through adoption) and, given covariates and the outcome's unit intercept, no
confounding between the timing of adoption and the outcome innovations. Compare the
population effects with the reference; a disagreement is the sign that these
assumptions fail.

```python
from longbet import LongBetEncourage

fit = LongBetEncourage(num_chains=4, num_burnin=2000, num_sweeps=1000).fit(y, d, z, x, t=t)
pred = fit.predict(groups=baseline_group)
print(pred.table)        # posterior ITTs with diagnostics and support
print(pred.wald())       # ratio withheld when diagnostics fail
print(pred.reference)    # design-based confidence sets
```

```r
fit <- longbet_encourage(y, d, z, x, t = t,
                         config = list(num_chains = 4, num_burnin = 2000, num_sweeps = 1000))
pred <- predict(fit, groups = baseline_group)
pred$effects; encouragement_wald(pred); pred$reference
```

Sampler options (`config` in R, keyword arguments in Python) apply to both
equations. Binary outcomes use `outcome="binary"`; their contrasts are
probability differences. The package's convergence gate (bulk and tail ESS of
at least 400, rank R-hat at most 1.01 on every reported quantity) needs the
default four chains of 1,000 retained draws to be attainable; check
`pred.stability()` on your panel. On some panels a few conditional profiles
wander slowly between runs by more than their within-run Monte Carlo errors
suggest (a few tenths of a posterior SD); when a decision hinges on the
ranking of a handful of profiles, refit with a second `random_seed` and
compare. `tempering_levels` is available in the
Python engine, not the R API, for panels whose forest modes the local moves
connect too slowly; the fit warns when adjacent replicas exchange in fewer
than 5% of attempts (`swap_rates_` on each equation), because a ladder that
does not exchange only multiplies the cost.

## Population target and the reference

`predict()` evaluates every study unit under both schedules (offered from the
launch period onward, never offered) with its own fitted intercepts, then
averages with equal weights (or `weights`) within each prespecified `groups`
label. That is the finite study population the reference describes, so
`pred.comparison()` is on matched targets. `pred.draws["itt_y"]`,
`pred.draws["itt_d"]` and `pred.draws["wald"]` have axes
`(group, horizon, chain, retained_draw)`.

To estimate the sampling covariance of model and reference estimates, run the
paired bootstrap (`fit.bootstrap_comparison`, R
`encouragement_bootstrap_comparison`): whole unit records are resampled within
arm and group, both estimators are recomputed on each resample, and flags are
withheld if any fit fails diagnostics or fewer than 100 replications are run.

## Unit-level effects and decisions

`predict_conditional()` describes a unit with a given covariate profile drawn
from the fitted population: both unit intercepts are integrated out, so the
same rule applies to study units and to a new cohort.

```python
cond = fit.predict_conditional(x=new_cohort)       # or fit.predict_conditional() for the study profiles
cond.citt_y.mean, cond.citt_d.mean, cond.cace.median   # (N, H) arrays; draws in .draws
cond.principal_strata()                            # compliers, always-takers, never-takers by horizon
cond.cumulative_lift(); cond.breakeven_probability(cost=20)
cond.knapsack_policy(cost=20, capacity=80, budget=1600, hurdle=0.5)
```

```r
cond <- predict_conditional(fit, new_x = new_cohort, cost = 20, capacity = 80, budget = 1600)
cond$citt_y$mean; cond$strata; cond$decision
```

Adoption contrasts are clipped at zero under `monotonic_first_stage` (no
defiers); the CACE is the draw-wise ratio of the two offer effects with a
small floor on the adoption contrast. The decision quantities target the
offer's effect on the outcome, which is what an allocation of offers changes.

## Blocks, clusters, known probabilities and staggered cohorts

`EncouragementDesign` with `design_encouragement_effects()` (R:
`encouragement_design()`, `design_encouragement_effects()`) provides
design-based reference estimators for blocked, clustered, Bernoulli and
staggered-cohort randomization. These are reference analyses only; the model
above is for the single-wave individually randomized design.

## What was retired and why

Earlier versions offered a reduced-form model (outcome on the offer clock,
adoption as a probit or linear working response, coupled by SUR), a smooth
stump engine, an orthogonalized two-stage engine and a fixed-basis hazard.
An adoption jump at a random period that the outcome equation does not model
is unit-specific, time-structured residual noise, and flexible forests fit it;
the exact posterior of the reduced form was less accurate out of sample than
early-stopped chains, and none of the alternatives improved on the design-based
reference. The adoption-clock model recovers the offer effects with less than
half the error on the package's simulated panel and its chains agree. The
record is in `benchmarks/tempering_plan.md` and `benchmarks/encouragement/`.
