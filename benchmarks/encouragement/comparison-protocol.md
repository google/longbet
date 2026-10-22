# IV comparison protocol, 2026-09-10

This protocol is written before examining comparative results for the repaired
public direct-smoothing implementation. Its purpose is a continuation decision,
not to certify uniform validity or to select favorable examples for a vignette.

## Candidate and comparators

The candidate is `LongBetDirectSmooth`, after fixing assignment-standardization
and sampler correctness. It uses four baseline and four effect stump trees,
fixed GP length scale 2 and nugget .05, four chains, 500 warmup and 1,000 retained
draws per chain. The primary candidate uses independent intercepts; the repaired
correlated version needs separate calibration and is not silently substituted.
Sampler controls may be increased if diagnostic failures require it, with the
original results retained. Do not tune priors against evaluation truths.

Comparators are (1) unadjusted assigned-arm means with joint Fieller inversion,
(2) arm-interacted linear ANCOVA for both reduced forms, and (3) arm-interacted
cubic-spline ANCOVA. Both adjusted methods include pre-assignment outcome means
and baseline covariates, standardize to the same empirical population, and use
joint HC2 residual covariance for Fieller inversion. They are horizon-specific
IV/reduced-form estimators with different nuisance adjustments. A comparison
to naive OLS on actual adoption is not evidence of improvement over IV.

## Targets and data

Generate both potential adoption and outcome panels before complete random
assignment. Use 160 accounts, two baseline weeks and four follow-up weeks.
Preserve latent adoption/outcome confounding, and distinguish empirical-population
ITTs from CACE. Current-adoption CACE is evaluated only where specified exclusion,
monotonicity, and a positive finite-population first stage hold. Zero-stage ratios
are undefined, not assigned zero. All methods receive identical observed data.

Prespecified scenarios: linear constant effects; nonlinear smooth effects;
nonlinear abrupt effects; weak instruments; zero instruments; and serially
correlated, heteroskedastic errors. The last three assess limits, not opportunities
to inflate ratio RMSE by comparing against nearly zero sample denominators.
The duration estimator has a separate conventional-2SLS equivalence check; changing
the treatment definition does not by itself establish superiority over 2SLS.

## Evaluation

Use independent pilot seeds only to verify execution and diagnose computation.
The evaluation stream starts at seed 2026091100, with at least 100 independent
datasets per scenario if feasible. Store every run, failure, seed, runtime,
source hash, posterior diagnostics, point estimate, interval/set, and truth.
Report pointwise coverage at each horizon, not only an average over horizons;
retain joint-experiment coverage as a separate descriptive measure. No padding
intervals, conditioning coverage on success, dropping difficult datasets,
collapsing disjoint Fieller sets, or treating grid endpoints as infinite tails.

Primary accuracy criterion: paired mean squared error for CACE in the three
strong-stage linear/smooth/abrupt scenarios. Also report both ITTs, interval
coverage and width, failures, diagnostics, and runtime. Undefined/infinite
estimates are recorded separately. Paired bootstrap intervals resample whole
datasets, retaining all horizons together. A credible accuracy advantage requires
an upper 95% paired-bootstrap RMSE ratio below 1 against both adjusted comparators
in a prespecified scenario, without an evident coverage deficit. A reduction
of at least 10% in RMSE is considered practically meaningful.

For a continuation recommendation, the candidate must also have usable sampler
diagnostics on at least 95% of evaluation fits, and the pointwise coverage evidence
must not demonstrate undercoverage exceeding five percentage points. At this
sample size Monte Carlo uncertainty can preclude a verdict; report that outcome
as unestablished rather than converting uncertainty into success. A narrow
interval does not compensate for undercoverage. Any recommendation to continue
must name the settings with an observed advantage and its computational cost;
it is not a claim that LongBet dominates IV generally.
