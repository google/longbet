# Cross-fitted longitudinal IV prototype: evaluation protocol

Written 2026-09-10 before comparative evaluation. This is a new, bounded
experiment prompted by the negative direct-posterior comparison, which remains
available in `comparison-report.md`. This study does not replace those results.

## Question and scope

Can LongBet's longitudinal predictor improve estimation of the same
finite-population IV effects relative to spline and tree predictors when all
three use the same cross-fitted, randomization-aware inference procedure?
The proposed contribution is prediction and precision, not new identification.
The prototype does not develop a structural exposure-duration model or a joint
hazard/outcome causal likelihood.

## Common inference and information

The design is individual-unit complete randomization, one assignment wave,
complete balanced panels, and no interference. Within each assignment arm,
uniformly allocate accounts to two folds with prespecified arm counts. Conditional
on the split, the fold assignments are independent complete randomizations.
Train on one fold and predict the other, reversing roles. All observations from
an account stay together. Account-level baseline covariates and all observed
pre-assignment outcomes/adoption indicators are shared inputs to every learner.
Held-out post-assignment responses must never train or tune their own predictor.

For both outcomes and adoption, average the arm prediction difference plus the
inverse-probability-weighted observed residual correction, using each fold's
conditional assignment probability. Estimate the full outcome/adoption/horizon
covariance by summing fold-weighted, within-arm residual covariance matrices.
This is a finite-population Neyman construction, not the covariance of iid panel
rows. The ITT estimates are conditionally unbiased under the stated design;
the covariance and confidence statements are asymptotic and require prediction
stability and suitable moment/CLT conditions. They are not automatically exact
for arbitrary learners or small samples.

Invert the studentized average moment `ITT_Y - beta * ITT_D = 0` using that joint
covariance. Preserve bounded, disconnected, unbounded, all-real, singleton,
empty and unavailable sets. Do not divide first and apply delta-method intervals,
truncate first stages, clip final ITTs, or count finite posterior quantiles as
identification. This is asymptotic Fieller/AR-style inference for an average
moment, not exact permutation inference under an imputable sharp null.

## Predictors and tuning

Compare a LongBet GP-leaf forest predictor, a pooled spline/ridge predictor,
and an ExtraTrees trajectory predictor. The LongBet predictor must actually
reuse LongBet's temporal forest model and support prediction for new accounts.
All predictors use the same inner account holdout and three prespecified
hyperparameter candidates, scored by equally weighted outcome/adoption mean
squared prediction error after training-only response standardization.
The winning configuration is refit on the entire outer training fold.
LongBet's temporal choices include no pooling; interaction support and all
sampler settings are recorded in the configuration manifest.

Before freezing configuration, separate pilot seeds may be used for numerical
correctness, runtime, and prediction stability under a larger sampling budget.
No comparison with evaluation truths may guide configuration. Exact candidate
grids, seeds and source hashes will be frozen in the launch manifest before
evaluating the first comparative dataset. Failures are retained, not silently
replaced with a different learner.

## Experiments

Retain all six original generators: linear, smooth, abrupt, weak, zero, and
serial/heteroskedastic, with 160 accounts, two baseline and four follow-up
periods. Use 100 fresh experiments per scenario as guardrails.

Add two primary scenarios, `long_smooth` and `long_abrupt`, with 240 accounts,
six baseline and 16 follow-up periods, measured covariate interactions, latent
readiness affecting adoption and untreated outcomes, and heterogeneous absorbing
adoption timing. Use 500 independent experiments per primary scenario. This
gives 1,600 independent experiments in total. Potential outcomes and both
potential adoption paths are fixed before complete 50/50 assignment. The long
scenarios include always-adopters, accelerated adopters, encouragement-only
adopters and never-adopters, preserving a meaningful positive first stage across
the evaluated horizons. Smooth and abrupt outcome patterns are both required.

Outcomes in these experiments depend on current adoption, so horizon-specific
exclusion and monotonicity give a valid current-adoption CACE target. This does
not impose current-adoption exclusion on the vignette's duration example.
At zero relevance, the CACE target is undefined, not zero. Report reduced-form
coverage and null-moment rejection rates there for fixed beta values -2, 0, 2 and 10. Weak-stage ratio RMSE is
descriptive only and cannot justify continuation.

## Continuation rule

The primary score is CACE RMSE over all horizons and all experiments in each of
the two long scenarios. A case to continue requires at least a 10% reduction
against **both** competing predictors in at least one primary scenario. Each
corresponding paired-bootstrap upper bound must be below one: use one-sided
98.75% bounds (Bonferroni across two primary scenarios and two comparators).
Bootstrap 20,000 whole experiments, retaining their horizons together. Also
report ordinary 95% paired intervals for interpretation. Incomplete pairs are
reported; missing fits cannot be silently removed to establish a win.

Accuracy alone is insufficient. Require at least 99% completed fits. For each
primary scenario, every reported ITT/CACE horizon must have at least 92%
empirical coverage of its nominal 95% set; the average pointwise coverage for
each quantity must be at least 94%. Report Wilson intervals and Monte Carlo
uncertainty, including all unsuccessful calculations in coverage denominators.
These are operational evaluation thresholds, not mathematical validity claims.
At 500 repetitions, coverage near 95% has about one percentage point Monte Carlo
standard error. A guardrail failure means clear evidence of coverage below 90%
(Wilson upper limit below 90%), or a strong-stage RMSE ratio whose 95% lower
bound exceeds 1.10 against the best competing predictor.

Report interval/set availability and widths, prediction errors, tuning choices,
sampling stability, and total wall time. A precision-only result is secondary
and does not substitute for the primary RMSE criterion. Any positive conclusion
must name its scope and cost. If the criteria are not met, recommend stopping
the LongBet-specific IV model effort while retaining useful, tested generic
inference infrastructure.

## Records and references

Preserve every seed, source hash, configuration, failure, estimate, confidence-set
component and runtime. Keep the original posterior comparison and this new
experiment distinct. Freeze evaluation code before inspecting aggregate results.

- [Lu et al., conditional cross-fitting for randomized experiments](https://arxiv.org/abs/2508.15664).
- [Chernozhukov et al., double/debiased machine learning](https://academic.oup.com/ectj/article/21/1/C1/5056401).
- [Ma, identification-robust inference with high-dimensional covariates](https://arxiv.org/abs/2302.09756).
