# Randomized encouragement: implementation plan

Date: 2026-09-09

Status: the Python/R interfaces, calibration harnesses, Bayesian reduced-form
wrapper, and explicit richer-design reference estimators are implemented. The
optional structural milestone is delivered as its required identification and
benchmark assessment, with no structural sampler released. Calibration runs,
focused Python/R checks and the broad regression suite are complete.
Posterior calibration gates remain
unmet: configurations are explicitly unvalidated and no LPM/probit default is
promoted. Completion and validation evidence are recorded at the end of this file.

## 1. Objective and scope

Extend LongBet to analyze repeated outcomes after randomized encouragement when
actual adoption is voluntary. Expose three distinct quantities:

1. The effect of encouragement on outcomes, `itt_y`.
2. The effect of encouragement on adoption, `itt_d`.
3. Their ratio, `wald`, indexed by time since encouragement and eventually by
   prespecified baseline subgroup.

The first release supports a single encouragement wave and a permanent holdout,
with complete randomization at the individual-unit level: conditional on arm
sizes, every assignment of units to arms is equally likely. Equal-probability
independent Bernoulli assignment also has this property conditional on arm sizes.
Randomization is an assumption supplied by the analyst; observing a binary `z`
matrix does not establish it.

Initial data requirements:

- `z` and `d` are complete binary `(N, T)` panels. Both are absorbing.
- All encouraged units start in the same observed period. Both arms are present.
- `y` is a complete finite `(N, T)` panel, continuous or binary. Binary outcomes
  are analyzed on the probability scale by the design-based estimator.
- `t` is strictly increasing, with positive whole-unit gaps; uneven gaps are
  allowed. Defaults to `1..T`. Preserve calendar time and the existing
  `derive_exposure` index separately from the observation column index.
- Perfect compliance and encouragement beginning in the first period are valid.
  Baseline periods improve some analyses but are not required for randomized ITT
  identification. Always-takers and adoption before encouragement are allowed.
- If an arm has fewer than two units, retain point estimates but mark variance
  estimates and confidence sets unavailable. Do not invent zero standard errors.

Blocked, unequal-probability, cluster-randomized, staggered, sequentially
re-randomized, and incomplete panels are outside the first estimator's contract.
Reject detectable unsupported structures. State the randomization restrictions
in every inference entry point: they cannot all be detected from arrays.
Do not silently redefine repeated usage as ever-adoption; that changes the causal
treatment and may invalidate exclusion.

## 2. Estimands and assumptions

For a common post-encouragement period `h`, let `Y_i^h(z)` and `D_i^h(z)` be
potential outcomes and adoption stocks under assignment `z`. Target the same
finite population for both contrasts:

```
Delta_Y(h) = mean_i[Y_i^h(1) - Y_i^h(0)]
Delta_D(h) = mean_i[D_i^h(1) - D_i^h(0)]
W(h)       = Delta_Y(h) / Delta_D(h), when Delta_D(h) != 0.
```

Random assignment, consistency, and no interference identify the two ITTs.
Calling their ratio a CACE requires exclusion, monotonicity, a nonzero population
first stage, and a treatment definition that accounts for relevant history.
With duration-dependent outcomes, exclusion through an entire adoption path is
not enough to turn an adopted-by-`h` ratio into the effect of `h` treatment periods.
Encouragement can accelerate people who would adopt in either arm by `h`.
Complier membership can also change across horizons.

Observed adoption lags describe behavior; they cannot verify both potential
adoption times or identify individual compliance types. Do not return a
`timing_assumption_passed`, `monotonicity_passed`, or certified CACE label.
No model conditions on realized post-encouragement uptake as an ordinary
unconfounded treatment to estimate ITT. Future covariate adjustment must use
baseline or otherwise justified covariates, not encouragement-induced mediators.

## 3. Architecture

Keep statistics in Python. R delegates to Python and only converts or plots
results. Do not change sampler conditionals in the first delivery.

| Component | Files | Responsibility |
| --- | --- | --- |
| Design and summaries | `src/longbet/_encourage.py` | Validate the supported design, calculate adoption summaries, optional plot |
| Reference inference | `src/longbet/_encourage.py` | Joint outcome/adoption moments and analytic test inversion |
| Python exports | `src/longbet/__init__.py` | Public utility entry points |
| R interface | `R/encourage.R`, `NAMESPACE`, `man/` | Equivalent utilities and documentation |
| Contract | `contract/longbet-api.yaml`, `inst/contract/longbet-api.yaml` | Inference inputs, output columns, assumptions, defaults |
| Tests | `tests/test_encourage*.py`, `tests/testthat/test-encourage.R` | Validation, inferential properties, interface parity |
| Bayesian integration | `_encourage_model.py`, `_encourage_predict.py` | Joint ITTs, common targets, binary accumulation, ratios, diagnostics, persistence |
| Richer design inference | `_encourage_design.py` | Fixed targets, groups, blocks, clusters, probabilities, staggered contrasts and joint covariance |
| Calibration | `benchmarks/encouragement/` | Repeated experiments, convergence, recovery, coverage, research assessment |

Inference-only arguments such as `alpha` belong to utility signatures and an
`encouragement_api` contract section, not `LongBetConfig`: they do not configure
the sampler. Any future sampler option must follow the existing config/contract/R
mirroring rule.

## 4. Milestone 1 — validate and describe the experiment

Implement:

```python
validate_encouragement(z, d, t=None) -> dict
encouragement_summary(z, d, t=None, *, alpha=0.05) -> pandas.DataFrame
plot_encouragement(z, d, t=None, *, ax=None) -> matplotlib.axes.Axes
```

Validation returns serializable metadata: design name, panel dimensions, arm
counts, encouragement period and column, presence of baseline periods, and
perfect compliance. Internal normalized arrays remain an implementation detail.

The summary has one row per observed calendar period, including baseline periods.
Compare *assigned* arms throughout; all baseline `z` cells may be zero but the
future encouragement arm remains identifiable. Include:

- Calendar period, zero-based column index, model exposure horizon, and whether
  encouragement has begun.
- Counts of independent units in each arm, not pooled person-period counts.
- Arm take-up rates, their difference and its Neyman standard error and normal
  pointwise interval. Use sample variances with `ddof=1`.
- A negative-first-stage flag and whether its interval includes zero; neither
  constitutes a monotonicity test or proves instrument validity.
- Among encouraged units observed adopted by that period: number adopted,
  fraction adopting exactly at encouragement, median observed adoption lag, and
  fraction adopting before encouragement. Undefined fractions are missing.

Plot assigned-arm adoption curves against calendar time with the gap shaded and
encouragement time marked. Keep plotting imports optional. R returns a ggplot
from the Python-produced table.

Acceptance:

- Invalid dimensions, nonbinary/nonfinite panels, reversals (including boolean
  and unsigned arrays), multiple waves, invalid times, and missing arms fail
  with input-specific messages.
- Perfect compliance, first-period encouragement, negative sample first stages,
  pre-existing adoption, and delayed adoption work.
- Counts, moments, and lag summaries agree with small hand-computed panels.
- Exposure indexing agrees with `derive_exposure`, including uneven time gaps.
- No function claims to validate randomization, exclusion, or individual types.

## 5. Milestone 2 — reference ITTs and weak-instrument inference

Implement:

```python
encouragement_effects(y, d, z, t=None, *, alpha=0.05) -> pandas.DataFrame
```

Return one row per observed post-encouragement period. No smoothing, horizon
pooling, covariate adjustment, outcome imputation, or sampler is involved.
Every row uses both outcomes on the same units and preserves outcome/adoption
covariance. For arm `a`, form the sample covariance matrix of `(Y, D)` and set

```
V = sample_cov_arm_1 / n_1 + sample_cov_arm_0 / n_0.
```

This is the Neyman covariance estimate, conservative in the finite-population
sense because the covariance of individual assignment effects is unobserved.
It is not a finite-sample exact sampling covariance.

Report the two differences in means, their SEs and pointwise normal confidence
intervals, `V_YD`, the unmodified signed point ratio (missing if its denominator
is exactly zero), and a Fieller/Anderson–Rubin-style confidence set obtained by
inverting

```
(Delta_hat_Y - w * Delta_hat_D)^2
    <= z_(1-alpha/2)^2 * (V_YY - 2*w*V_YD + w^2*V_DD).
```

Solve the resulting quadratic inequality analytically. Preserve bounded,
disjoint, half-line, all-real, singleton, and empty sets; do not truncate to a
search grid or discard draws/values with negative denominators. This uses an
asymptotic normal critical value, **not** permutation critical values and **not**
an exact randomization interval. Its practical adequacy still requires suitable
arm sizes, moment conditions, and calibration. Inverting the test avoids dividing
by the estimated first stage to construct the set; it does not rescue invalid IV
assumptions or make an unidentified population ratio exist.

Use portable scalar columns for at most two closed interval components:
`wald_lower_1`, `wald_upper_1`, `wald_lower_2`, `wald_upper_2`, plus
`wald_set_type` and `wald_reason`. Infinity remains infinity, absent components
are `NaN`. An all-real set is informative evidence of inadequate identification,
not a computation failure. Label these frequentist confidence sets explicitly;
future posterior intervals occupy separate fields.

Acceptance:

- Hand-computed joint arm covariance, including the cross term.
- Set membership agrees with direct residualized-outcome tests across candidate
  effects; test bounded and unbounded outcomes and degenerate limits.
- Shifting outcomes leaves treatment effects unchanged; rescaling outcomes
  rescales the confidence set; negative rescaling reverses interval order.
- Binary outcomes give risk differences without a latent transformation.
- Perfect compliance reduces to the ordinary outcome difference interval.
- Zero or negative sample first stages remain representable.
- Tiny arms preserve point estimates and mark inference unavailable.
- Missing data fail before computing mismatched contrasts.
- R has the same numeric values, infinities, missing components, and defaults.

## 6. Milestone 3 — reproducible calibration harness

Create new benchmark scripts against the actual current sampler. Treat the old
branch's one-dataset comparisons as exploratory; do not carry its numerical
claims into tests as universal truths. Its configurations require explicit proper
innovation priors to run on current `main`.

Build a DGP returning both potential adoption paths, potential outcome paths,
assignment, baseline covariates, and finite-population ground truth. Include:

- Strong, weak, near-zero, and zero first stages over several sample sizes.
- One- and two-sided noncompliance, delayed adoption, adoption acceleration,
  negative and heterogeneous effects, baseline confounding of uptake and outcome.
- Correlated unit effects, serial errors, heteroskedasticity, binary outcomes.
- Exclusion violations and defiers as labeled assumption failures, not recovery
  cases that the estimator must somehow correct.
- Different seeds for data generation and posterior sampling.

First exercise reference inference over at least 500 replications for selected
settings; report Monte Carlo SEs for coverage and refusal/undefined rates.
Do not assert that all intervals from one seed cover truth. Record pointwise
coverage separately from simultaneous coverage; the first implementation supplies
only pointwise intervals. Preserve failed runs and exact environments.

Then compare proper-prior LPM, probability-scale probit, separate marginal fits,
and joint fits. Evaluate both finite-sample conditional standardization over
observed unit effects and population marginalization where that is the target;
do not mix their estimands. Record ESS, R-hat, MCSE, chain lengths, and sensitivity
to longer chains before diagnosing persistent errors as model bias.

Release gate: meaningful calibration failures must be explained or the affected
configuration remains unsupported. Passing convergence checks alone is not a
coverage claim. Passing average recovery alone is not convergence.

## 7. Milestone 4 — Bayesian encouragement wrapper

After milestones 1–3, implement `LongBetEncourage` around `LongBetMulti` for the
supported design. Initial model target: joint reduced forms of `Y` and `D` on
encouragement `Z`; both use the encouragement clock.

Tasks:

1. Choose and document proper variance priors on the standardized outcome scale.
   Keep LPM a candidate first stage until calibration supports the default.
2. Fit and predict the two ITTs on a common target population with equal unit
   weights. The finite-population reference averages over all study units;
   existing `get_att` averages over encouraged cells, so do **not** substitute
   its output silently. Implement common-population standardization or expose
   the distinction explicitly and align the reference before comparison.
3. Return `itt_y`, `itt_d`, posterior Wald medians/quantiles, and reference
   confidence sets with distinct names and metadata. No posterior mean of the
   ratio as the default summary.
4. Add natural-scale binary-outcome accumulation under `summary_only=True`, or
   explicitly restrict the wrapper to continuous outcomes until it exists.
5. Preserve actual joint draws. Verify covariance calibration at the independent
   unit level; shared draw indices alone are not evidence of correct uncertainty.
   Benchmark independent versus correlated unit intercept models before relying
   on posterior ratios. Correlated-intercept sampler work is a separate change
   with mandatory conditional/Geweke checks if needed.
6. Diagnose both ITTs and the reported ratio/subgroup summaries, including tails.
   Keep convergence, first-stage precision, practical relevance, and model/reference
   disagreement separate. A 5-percentage-point threshold is not instrument strength.
7. A disagreement flag compares the same target and accounts for uncertainty and
   covariance in the comparison; it does not declare any disagreement proof of
   misspecification.
8. Persist original design metadata and inference provenance, with archive
   versioning, R round trips, and replay of the common target.

Acceptance: reference-versus-model recovery and covariance calibration,
probability-scale tests, blockwise/full prediction equivalence, fixed-seed algebra,
save/load and Python/R equality, appropriate posterior diagnostics, repeated-DGP
coverage. No sampler default is promoted from an isolated successful experiment.

## 8. Milestone 5 — prespecified subgroups and richer randomization

First add fixed baseline groups to the single-wave analysis, using the same
within-group units and weights for both ITTs. Groups chosen after inspecting
effects require separate sample-splitting or selection-aware methodology. Report
group-specific arm support and first-stage precision. Small groups are retained
with unavailable uncertainty where appropriate.

Then extend the design contract in separate increments:

1. Blocked assignment: explicit block labels, arm counts, fixed target weights.
2. Cluster assignment: explicit randomization-unit identifiers, whole-cluster
   resampling/covariance, cluster count diagnostics, defined unit/cluster weighting.
3. Unequal probabilities: known assignment probabilities and explicit weighted
   estimators; no automatic equivalence to unweighted arm means.
4. Staggered encouragement: cohort/calendar contrasts against eligible controls,
   common support, no-anticipation assumption, declared cohort/horizon weights,
   and covariance for reused controls. Distinguish design support from model
   extrapolation after controls disappear.

For each increment, specify the estimand and covariance before generalizing the
API. A pooled person-period binomial SE is never the default longitudinal SE.

## 9. Milestone 6 — optional structural research

Keep a treatment-clock structural model and longitudinal principal stratification
out of the initial release. Before either is implemented, require a separate
identification argument, adoption-history model, target population, and benchmark
protocol. A binary instrument does not generally identify an unrestricted
treatment-duration response. Correlated Gaussian errors impose structure rather
than supply identification for free.

A structural model may reuse per-equation designs and covariance machinery, but
must encode absorption coherently, separate instrument columns from outcome
forests under exclusion, justify extrapolation beyond compliers, and undergo
prior sensitivity and misspecification studies. Disagreement with a reduced-form
Wald ratio can have several causes; it does not uniquely diagnose normality.

## 10. Documentation and delivery checks

Each delivery includes a Python example and matching R example, public docstrings,
R help, exports, mirrored contract, and tests appropriate to the changed layer.
No statistical code is reimplemented in R. Initial inference code needs no Geweke
test because no sampler conditional changes.

Run focused new tests, contract and rollout tests, then the repository's required
fast suite as resources permit. Run the R tests in the existing container workflow
when R is unavailable on the host. Record actual failures and unavailable checks;
do not describe a skipped runtime test as passed.

## 11. Sources

The distinction between ITT and complier effects follows the motivating
[randomized encouragement chapter](https://book.martinez.fyi/iv.html).
The statistical review of the
[earlier proposal](https://github.com/ignacio82/longbet-jax/blob/claude/longbet-randomized-encouragement-4rzej3/docs/encouragement-designs.md)
motivated the scope and validation changes above.

For test inversion with randomized encouragement and heterogeneous effects, see
[Aronow, Chang, and Lopatto, Randomization-based Confidence Sets for the Local
Average Treatment Effect](https://arxiv.org/html/2404.18786v2). Their work
distinguishes randomization critical values from normal approximations and exact
homogeneous-effect claims from asymptotic guarantees. Milestone 2 implements
normal-critical-value inversion only.

## 12. Implementation record

- Milestone 1: implemented in Python and R, including optional plots, input
  guards, assigned-arm baseline summaries, and observational lag summaries.
- Milestone 2: implemented in Python and R, including joint Neyman moments,
  normal-test inversion, unbounded confidence sets, and portable return tables.
- Validation: 126 focused Python tests pass (encouragement, contract, rollout).
  The broader non-slow Python suite also passed: 436 tests, 10 slow tests
  deselected, in 18m55s with four CPU cores available to the process. Two added
  input/numerical edge-case tests and the final error-message refinement are
  covered by the subsequent focused run. The slow sampler studies were not run
  because this delivery changes no sampler conditional.
  R encouragement, contract, and rollout tests pass in the existing
  `longbet-rtest` container; new R help files pass `tools::checkRd`.
  Return tables strip reticulate's live pandas-index handle and pass serialization
  checks. The contract is mirrored byte-for-byte in `inst/contract/`.
- Milestone 3: an initial fixed-population reference calibration harness is
  implemented. Both the initial 500-assignment run and its extension to 5,000
  assignments per regime are retained. See the
  [calibration report](../benchmarks/encouragement/README.md).
- Known inferential limitation: the normal first-stage interval undercovers in
  the tested zero-stage discrete population. Exact enumeration yields 93.1945%
  coverage for a nominal 95% interval. This is documented and not described as
  calibrated or exact; alternative interval methods and the full DGP matrix
  remain future work. Wald coverage is near nominal in the tested nonzero-stage
  populations, which does not establish general calibration.
- Subsequent implementation: the [user guide](encouragement-guide.md) documents
  the completed Python/R APIs, examples, assumptions and remaining release gates.
- Milestone 3: the expanded potential-history DGP covers 14 scenarios at
  N=80, 400 and 1,200, with 500 independent complete assignments per setting
  (21,000 reference runs). It retains exact finite-population truth, expected
  Neyman covariance, cross-horizon covariance, Monte Carlo SEs, pointwise versus
  simultaneous inclusion, undefined ratios, failures and environment metadata.
  Selected suspicious 500-run results were extended on the same populations to
  5,000 assignments, preserving the initial observations. The model harness
  compares separate/joint proper-prior LPM/probit configurations, both conditional
  and population targets, actual joint covariance, tail diagnostics and longer
  chains. See the [recorded results](../benchmarks/encouragement/calibration-report.md)
  and [reproduction instructions](../benchmarks/encouragement/README.md).
  All 132 model fits completed without runtime failure: 128 short runs over
  32 independently generated datasets (four scenarios, four model variants),
  plus four longer-chain sensitivity fits on the first strong dataset. At least
  one horizon fails the convergence checks for each reported quantity in every
  configuration. Coverage, covariance, quantile MCSE and refusal rates remain
  diagnostic evidence, not a passed calibration gate. Grouping in the report
  preserves sample sizes and DGP parameters. The matched-prior Gaussian
  intercept study uses 30 datasets in each of two settings and 4,000 independent
  prior-predictive conditional checks per variant; its largest standardized
  moment discrepancy is 1.45. A deliberately broken conditional fails the check.
- Milestone 4: `LongBetEncourage` fits both reduced forms on encouragement with
  proper innovation priors and an explicit first-stage choice. Prediction uses
  all original study units, common equal weights and fixed baseline groups.
  Natural-scale probit accumulation works in bounded memory. Ratios preserve
  negative and near-zero denominators, omit ratio means, and withhold summaries
  for undefined draws or failed joint diagnostics. Both ITTs and each ratio have
  bulk/tail ESS, rank R-hat and endpoint MCSE; practical relevance is separate.
  The covariance-aware comparison uses paired whole-unit bootstrap refits and
  withholds flags when support or diagnostics fail. Archives preserve original
  inputs, calendars, design and inference provenance; R models and predictions
  survive a fresh-process `saveRDS` round trip. Fractional calendar origins are
  preserved while checking that float32 model times retain the exposure clock.
- Milestone 4 release gate: **posterior interval calibration is not established**.
  Longer retained chains still show mixing failures in the measured cases, and
  the small independent-dataset matrix cannot establish reliable covariance or
  coverage. The production sampler retains independent unit-intercept priors;
  SUR couples innovations. A separate Gaussian correlated-intercept candidate
  was benchmarked with conditional/Geweke checks, but it is not equivalent to a
  calibrated correlated-intercept LongBet sampler. No production sampler default
  or conditional was changed. These limits remain visible in API metadata,
  documentation and results rather than being converted into a passed gate.
  For binary Y plus probit D, both innovation-loading rows are fixed to zero;
  with default separate forests the equation posteriors are independent despite
  `sur=True`. Actual coupling and shared-partition settings are recorded in model
  metadata and validated on load. Optional shared partitions do not add binary
  innovation or unit-intercept correlation.
- Milestone 5: `EncouragementDesign` and `design_encouragement_effects` implement
  fixed baseline subgroup targets, explicit blocked complete assignment, whole
  cluster covariance with unit/cluster weighting, known unequal Bernoulli
  probabilities, and independent categorical staggered-cohort assignment.
  Staggered comparisons declare eligible controls and fixed cohort weights and
  retain covariance from reused controls. Unsupported horizons/control support
  remain unavailable without renormalization or model extrapolation. Fixed
  cohort quotas and adaptive sequential randomization remain outside the stated
  contract. Generalized designs provide reference inference, while the Bayesian
  wrapper retains its single-wave scope. Exhaustive small-population tests verify
  HT unbiasedness and the exact expected Neyman covariance gap for cluster targets.
- Milestone 6: the [structural research assessment](../benchmarks/encouragement/structural_research.md)
  specifies an adoption-history likelihood, target/identification gates, exclusion
  restrictions, extrapolation limits and benchmark requirements. As explicitly
  required in section 9, an unrestricted treatment-clock structural model and
  longitudinal principal-stratification sampler are not released.
- Latest focused validation: 221 encouragement, benchmark/DGP, contract and
  rollout tests passed. The final binary/binary coupling-provenance addition and
  archive regression are covered by a subsequent 11-test wrapper run. Checks include
  including fractional-calendar and optional-covariate archive replay, refusal
  for unsupported reference uncertainty, signed ratios and bounded-memory
  probability contrasts. R encouragement, model, design and contract tests pass
  in `longbet-rtest`, including fresh-process replay and actual bootstrap refits;
  all encouragement R help files pass `tools::checkRd`.
- Broad regression: **524 Python tests passed**, 10 slow tests deselected,
  four existing generic-model warnings, in 21m10s on four CPU cores. This run
  began before the final calendar/support guards, binary-coupling provenance and
  calibration-harness review corrections; the later focused runs above cover
  those final changes. The production sampler's slow studies were not rerun
  because no production conditional changed. The new benchmark-only Gaussian
  candidate received its own prior-predictive conditional checks and a regression
  demonstrating detection of a deliberately broken conditional.
- Delivery checks: Python exports and signatures match contract version 2.3.0;
  R defaults are checked against the same contract. Both contract files are
  byte-identical, local guide links resolve, and `git diff --check` passes. Raw
  benchmark replay archives remain locally available and Git-ignored; compact
  result records and reproducible scripts accompany the implementation.

### Follow-up sampler repair experiment

The [repair protocol](encouragement-repair-protocol.md) and
[results](../benchmarks/encouragement/repair-report.md) implement the subsequent
controlled comparison requested after the original mixing failures. New research
modules provide a direct smooth Gaussian forest, a parametric analogue, exact
joint Gaussian coefficient/intercept updates, and randomization-based AR tests.

- The original and stump-restricted LongBet controls fail the effect checks in
  all four clean strong-stage runs at 2,000 retained draws per chain.
- Direct smooth leaves substantially improve effect mixing. Parameter diagnostics
  expose a remaining baseline/intercept dependence, motivating exact joint
  updates of all baseline and effect leaves with the unit intercepts. The three
  direct-forest update variants target the same posterior; changing from the
  original LongBet product model to direct smooth leaves changes the prior/model.
- The final forest passes every effect, variance, unit-intercept, function-probe
  and aggregate tree-usage check in all eight final runs across the clean and
  richer strong-stage scenarios, with four chains and 8,000 retained draws per
  chain. The parametric control passes all eight runs at 2,000 draws per chain.
  Both ITTs and the raw Wald ratio retain their ESS-growth and tail checks.
- All shorter-chain and unsuccessful-update records remain in the report. Exact
  continuation preserves every saved prefix and RNG state; no signs or small
  denominators are clipped. Forty run records include eight continuation records.
- **81 focused tests passed**, including independent Gaussian moments and
  cross-covariances, three forest-update Geweke checks, broken-kernel detection,
  residual preservation, exact archive continuation, DGP identities, and
  exhaustive randomization-AR size checks. Subsequent focused checks cover final
  report and continuation changes. No production sampler conditional changed.

These are benchmark modules, not a new `LongBetEncourage` backend. Two independent
datasets per scenario cannot establish calibration, and the final candidates
still use independent Gaussian/LPM equations. Correlated intercepts, identified
binary latent responses, absorbing adoption hazards, broader sample sizes and
frequentist coverage remain the next research stages. The public calibration
status therefore remains `not_established`.
