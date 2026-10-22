# Encouragement calibration and model comparison

## Implementation and metric correction, 2026-09-10

The [completed comparison report](comparison-report.md) recommends stopping the
IV-method superiority effort. Across 600 independent experiments and 300 matched
correlated-intercept fits, neither tested configuration improves on both adjusted
baselines while meeting the accuracy, coverage and diagnostic criteria. Complete
records and checksums are in [comparison-results](comparison-results/).

The [current comparison protocol](comparison-protocol.md) evaluates the repaired
direct-smoothing estimator against both linear and spline ANCOVA/Fieller
comparators on identical potential-outcome targets. It prespecifies strong,
weak, null, abrupt-effect and dependent-error settings, retains failed fits, and
requires coverage and diagnostics alongside estimation accuracy. Historical
recovery examples and short chains do not establish superiority over existing IV.

The legacy filename `sbc_calibration.py` now contains **fixed-parameter repeated-
sampling coverage checks**, not simulation-based calibration. The earlier
binary-bounds metric checked only whether a set was nonempty, and the earlier
hazard metric expanded each interval by 0.25 before checking containment. Neither
was coverage of the stated target by the reported interval. The corrected script
uses a prespecified binary threshold and known complier potential-outcome contrast,
reports binary bounds separately from continuous-effect intervals, checks actual
endpoints, and leaves all unevaluated AR grid regions and tails unknown.

`CoupledHazardIV` now declares its modular hazard/working-outcome distribution,
uses proper hazard priors and exact conditional outcome draws, and withholds
`causal_ci`. Its finite `beta_ci` describes an experimental working coefficient.
The residual regression is not generally a valid control function, and a zero
relevance draw leaves its two outcome columns collinear after time effects.
Spike-and-slab regularization does not eliminate weak-instrument bias or supply
identification. Optional working-interval containment checks are not validated
causal inference.

```bash
PYTHONPATH=src python benchmarks/encouragement/sbc_calibration.py --replications 100 --output coverage_report.json
```

The older benchmark programs and reports below remain separate evidence; the
corrections above do not reclassify their results as validation of a new model.

## Earlier reference and model studies

The reference estimator has been exercised in **21,000 complete-randomization
replications** across 14 scenarios and three sample sizes, followed by 5,000-run
extensions of two selected fixed populations. Model comparisons record actual
joint and independent LPM/probit fits, their transformed-effect diagnostics,
conditional and population targets, and sensitivity to longer chains. **No
posterior configuration has established frequentist calibration.** The Bayesian
API requires an explicit first-stage choice and exposes that limitation.

The subsequent [repair comparison](repair-report.md) evaluates direct smooth
encouragement forests and a parametric control against the current model. It
records both effect and parameter failures, exact joint Gaussian updates,
four-chain comparisons, and continuations from saved states. The
[protocol](../../docs/encouragement-repair-protocol.md) separates sampler repair
from the later calibration and binary-adoption work. These candidates remain
research modules; `LongBetEncourage` does not select them as a backend.

A separate [randomization AR reference](randomization-ar-guide.md) tests supplied
ratio values using complete assignment enumeration or Monte Carlo assignments.
It labels the sharp-null validity conditions and leaves unevaluated grid regions
and tails unknown.

The [binary encouragement benchmark](binary-report.md) implements the subsequent
single-horizon joint adoption/outcome model, optional monotone uptake, and
identification bounds with direct-effect sensitivity. It includes exact cell
controls, 600 independent simulated datasets, prior sensitivity, and 24 forest
fits. An additional marginal uptake update resolves the sampled diagnostic
failures in the eight monotone pilot datasets. The cell experiments also expose
substantial weak-stage bias from the monotone prior. This remains a research
module; these results do not establish forest coverage or add a public backend.

The original [three-population reference report](initial-reference-report.md),
including its 500- and 5,000-assignment results and exact discrete undercoverage
calculation, is preserved. The extended matrix does not erase that failure.

## Reproduction and artifacts

Run from this checkout's Python environment:

```bash
python benchmarks/encouragement/reference_matrix.py --replications 500
python benchmarks/encouragement/reference_matrix.py --scenarios strong --sizes 400 --seed 20261909 --replications 5000 --output benchmarks/encouragement/reference_strong_5000
python benchmarks/encouragement/reference_matrix.py --scenarios delayed --sizes 80 --seed 20275909 --replications 5000 --output benchmarks/encouragement/reference_delayed_5000
python benchmarks/encouragement/model_calibration.py --scenarios strong weak binary correlated_serial_heteroskedastic --datasets 8 --limit 128
python benchmarks/encouragement/model_calibration.py --scenarios strong --datasets 1 --chain-multiplier 4 --limit 4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/intercept_candidate.py --datasets 30
python benchmarks/encouragement/build_report.py
```

The model runner resumes completed runs and leaves failed inputs and tracebacks
intact; `--retry-failures` explicitly retries failed runs. `--limit` bounds newly
attempted fits. Data-generation, assignment, integration and sampler seeds are
recorded separately. Each completed model run stores its configuration, warnings,
chain lengths, actual synchronized timings, model archive, original potential
paths, and effect draws. JSON records include Python/package versions, platform,
Git revision and dirty status. Raw NPZ replay artifacts remain on disk and are
ignored by Git to keep large model traces out of version control. The compact
JSON results, scripts, seeds, and reports are versioned.
`source-manifest.json` records SHA256 hashes of the replay and engine source
snapshot at benchmark completion; later diagnostic refreshes leave posterior
draws and data unchanged.

## DGP and targets

`dgp.py` creates both potential adoption histories and both potential outcome
histories **before** assigning encouragement. All 14 scenarios use baseline
covariates, latent uptake propensity, heterogeneous effects, serial innovations,
and fixed compliance-type counts. Scenarios vary strong/weak/near-zero/zero
first stages; one-sided compliance; delays; acceleration; negative heterogeneous
effects; strong intercept and serial correlation with heteroskedasticity; binary
outcomes; and labeled exclusion or monotonicity failures. Sample sizes are 80,
400 and 1,200. A nominal 0.5% complier share rounds to zero at N=80; metadata and
truth retain that fact. Assumption-failure cases assess assignment contrasts and
their ratio, **not recovery of a CACE despite invalid assumptions**.

Reference inference targets the realized finite population and repeats assignment
while holding its potential paths fixed. Conditional model standardization holds
observed unit effects fixed and is compared with conditional expected potential
outcomes; the realized finite-population truth is also retained and evaluated.
For binary outcomes those truths differ because the probability contrast averages
outcome innovations. Population standardization integrates fresh latent units at
the empirical baseline-X distribution; the DGP integration uses 256 independent
draws and reports its own Monte Carlo error. It must not be substituted for the
finite-study-unit reference target. A fitted normal-intercept marginalization can
be misspecified for this DGP; that is a model assumption to evaluate, not a target
identity to impose.

## Reference evidence and limits

All 42 settings completed their 500 assignments without estimator failures.
The [matrix manifest](reference_matrix_results/manifest.json) links each result.
The [generated numerical report](calibration-report.md) gives coverage and
covariance comparisons. Rates include Monte Carlo SEs. Simultaneous coverage is
computed as inclusion at **all** horizons in each experiment; pointwise intervals
are not relabeled simultaneous intervals.

The matrix stores full cross-horizon, cross-outcome covariance. Under complete
randomization its exact target is `S1/n1 + S0/n0 - S_effect/N`, whereas the mean
Neyman covariance estimate targets `S1/n1 + S0/n0`. The tests enumerate all 70
assignments of an eight-unit population and verify both identities. Comparing an
estimated Neyman matrix with the smaller exact covariance as if equality were
required would incorrectly call deliberate conservatism a calibration bug.

Two low 500-run coverages prompted extensions of exactly the same populations and
assignment streams, preserving the initial results. Strong N=400 Wald coverage
was 92.6–93.8% initially and 94.92–95.44% over 5,000 assignments. Delayed N=80's
first-horizon coverage was 91.2% initially and 94.34% over 5,000; across horizons
the extended coverage was 94.24–94.80%. Extended Monte Carlo SEs are approximately
0.30–0.33 percentage points. These are population-specific results, not proof of
uniform 95% coverage.

The previously identified zero-stage discrete example still has **93.1945%
exact randomization coverage** for the nominal 95% normal first-stage interval.
Normal/Fieller inversion is an asymptotic method. Unbounded or disjoint sets and
undefined population ratios remain meaningful outputs; the implementation does
not impose artificial denominator bounds. An all-real set is not a software
failure, and a bounded set does not establish instrument validity.

## Model comparison and release decision

`model_calibration.py` compares joint triangular-SUR fits and the exact independent
likelihood obtained with `sur=False`, using independent forests (`num_shared_trees=0`)
and the same proper IG(2,1) innovation priors. Both variants retain **independent
unit-intercept priors across outcomes**. SUR models observation-level innovation
dependence; aligned draw indices or an estimated SUR covariance do not establish
the correct covariance of unit-level assignment effects.
When both outcome and uptake use probit, both SUR loading rows are fixed to zero
under the current scale identification. With shared trees disabled, the nominal
joint and separate binary/binary fits are the same uncoupled model, and their
identical results are expected. Mixed binary/continuous fits can have nonzero
SUR loadings. Actual coupling, rather than the requested switch alone, is
recorded in the run metadata.

The numerical report records coverage against both model targets, all ITT/ratio
ESS and R-hat values, tail diagnostics, mean MCSEs for the ITTs, quantile MCSEs
for both ITTs and the ratio's median and interval endpoints, and posterior
versus empirical error covariance across independently generated experiments.
Each short configuration has eight experiments, so its coverage estimates have
large Monte Carlo uncertainty and cannot validate a nominal 95% interval.
Unconditional inclusion is recorded even when an interval is unavailable; the
undefined-draw fraction is separate. Ratios retain negative and arbitrarily small
denominators. Their median and quantiles are summarized; the ratio's mean and
mean MCSE are omitted because its mathematical posterior mean need not exist.

Short runs use two chains with 400 burn-in and 200 saved draws, and five trees per
forest. The selected longer comparisons use 1,600 burn-in and 800 saved draws.
The stronger data do not eliminate mixing concerns: the longer first-dataset
runs still have maximum transformed-effect R-hat values around 1.07–1.16 and
insufficient tail ESS in several quantities. This is evidence that these runs
are inadequate for a calibration claim; it does not establish persistent model
bias or a general LPM-versus-probit ranking. No default is promoted from these
comparisons, and posterior calibration remains `not_established`.

## Correlated-intercept candidate

`intercept_candidate.py` implements a **benchmark-only Gaussian linear reduced
form** with either independent inverse-gamma intercept variances or a full
inverse-Wishart intercept covariance. Both candidates have matched IG(3,2)
marginal variance priors; the correlated prior is IW(7,4I). It uses known diagonal innovation
variances and an LPM working response for uptake. This is deliberately separate
from the production forests and does not enforce a probit or adoption hazard.

Its exact Gaussian/inverse-Wishart conditionals pass 4,000 independent
prior-to-data-to-Gibbs (Geweke) transitions per variant; the largest absolute
standardized discrepancy across checked moments is 1.45. The
[verification record](intercept_results/geweke.json) and
[30-dataset-per-setting comparison](intercept_results/summary.json) are retained.
Conditional correctness is not a claim that this working model fits the DGP.
The deliberately broken all-zero sampler fails the same verification in a
regression test, including the zero-MCSE/nonzero-discrepancy case. An earlier
comparison used unmatched IG(4,3) independent priors; it is preserved in
`intercept_results_unmatched_priors/` and is not the reported matched-prior comparison.

The correlated candidate permits nonzero posterior covariance of the ITTs, which
the independent prior/likelihood cannot represent in its exact posterior.
It does not in these comparisons reproduce every unit-level covariance under
serial correlation and heteroskedasticity. The small study and broad working
intervals do not establish a replacement model. Production correlated-intercept
work would need its own binary scale identification, conditional and Geweke
tests, archive/contract support, prior sensitivity and repeated-DGP calibration.

The [structural research assessment](structural_research.md) separately records
the identification gates for treatment-clock models and principal stratification.
Those models are outside this release, as specified in the implementation plan.
