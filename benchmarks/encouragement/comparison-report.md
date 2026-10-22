# Repaired LongBet IV: matched comparison

**Recommendation: stop developing this as a superior IV method on the present evidence.**
The repairs fix substantive statistical and interface errors, but neither evaluated
configuration demonstrates a reliable advantage over both adjusted IV comparators.
The evidence covers 600 independent simulated experiments and an additional 300
correlated-intercept fits on the same strong-stage datasets.

## What changed

- Direct smooth ITTs now average full-population assignment contrasts. Previously
  they saved assignment-weighted fitted contributions, which could nearly cancel.
- Correlated-intercept Gibbs updates now use the correct conditional prior mean
  and variance. Independent dense-Gaussian checks verify the corrected kernel.
- Python and R dispatch honor sampling controls, reject unsupported options and
  preserve the direct model and covariance traces in versioned archives. Legacy
  direct archives require refitting.
- Duration IV now uses rank-checked, unregularized panel 2SLS, unit-cluster CR1
  covariance with fixed-effect degrees of freedom, and t(N−1) critical values.
  It rejects deficient rank, nonzero ridge penalties and unsupported panels.
- Hazard relevance updates now use a coherent collapsed Gaussian block, with
  prior leaf draws under zero relevance and one fixed random basis across chains.
  This remains an experimental conditional-basis model.
- The coupled hazard/outcome class now targets an explicit modular distribution
  and always withholds causal intervals: `D - Dhat` is not generally a valid
  binary-adoption control function. Proper priors do not resolve its causal
  nonidentification at zero relevance.
- The old coverage script no longer counts nonempty bounds as coverage, pads
  coupled intervals, invents binary targets, or calls fixed-parameter simulations
  Bayesian simulation-based calibration.

Validation passed: **417 targeted Python tests and 65 Docker R assertions**.
The broader Python run excluded 14 tests marked slow; it was supplemented by the
15 repaired-hazard tests. Fresh R package installation and chapter rendering
exercise the matching R/Python interfaces.

## Matched design and continuation criterion

The [primary protocol](comparison-protocol.md) preceded inspection of comparative
evaluation results. Each scenario has 100 independently generated datasets, with
160 units, two baseline and four follow-up periods. Potential adoption and outcome
panels precede complete 80/80 assignment. Latent readiness affects both adoption
and untreated outcomes. Targets are the same finite-population ITTs and complier
effects for all methods. Effects are time-constant in the linear scenario but can
vary with covariates; the other strong scenarios have smooth or abrupt time effects.

The adjusted comparators fit outcome and adoption reduced forms separately by
assignment arm, using baseline outcome means and measured covariates. The spline
alternative adds fixed cubic bases. Both standardize to the same units and use
joint HC2 covariance with complete Fieller confidence sets. Independent algebraic
tests match a full interacted-regression sandwich and the intercept-only reference.

LongBet uses four baseline and four effect trees, four chains, 500 warmup and
1,000 retained draws per chain, with fixed default priors. Pilot seeds are separate.
The primary model has independent unit intercepts. A continuation case requires
at least a 10% CACE RMSE reduction and a paired-bootstrap upper ratio limit below
one against **both** adjusted alternatives, acceptable pointwise coverage, and
at least 95% of fits passing all reported diagnostics. Bootstrap intervals use
20,000 resamples of whole datasets, keeping horizons together. These are nominal,
scenario-specific intervals, without a simultaneous multiple-comparison guarantee.

## Primary strong-stage results

| Scenario | Unadjusted IV RMSE | Linear adjusted IV RMSE | Spline adjusted IV RMSE | LongBet direct RMSE | LongBet CACE coverage | LongBet diagnostic fits |
| --- | --- | --- | --- | --- | --- | --- |
| Linear | 1.195 | 0.347 | 0.371 | 0.356 | 97–99% | 90/100 |
| Smooth nonlinear | 1.341 | 0.464 | 0.409 | 0.425 | 98–100% | 94/100 |
| Abrupt nonlinear | 1.245 | 0.470 | 0.384 | 0.434 | 98–99% | 99/100 |

| Scenario | Comparator | RMSE ratio [95% paired interval] |
| --- | --- | --- |
| Linear | Linear adjusted IV | 1.023 [0.981, 1.069] |
| Linear | Spline adjusted IV | 0.957 [0.910, 1.007] |
| Smooth nonlinear | Linear adjusted IV | 0.914 [0.843, 0.994] |
| Smooth nonlinear | Spline adjusted IV | 1.038 [0.954, 1.126] |
| Abrupt nonlinear | Linear adjusted IV | 0.924 [0.868, 0.985] |
| Abrupt nonlinear | Spline adjusted IV | 1.130 [1.068, 1.197] |

No strong-stage scenario meets the accuracy criterion. The smooth and abrupt
settings improve slightly on linear adjustment, but spline adjustment performs
better. In the abrupt setting, LongBet has approximately 13% greater RMSE than
the spline comparator, with the paired interval entirely above one. Much of the
gain over raw Wald estimates is already obtained by inexpensive baseline adjustment.

## Calibration and stress cases

The following are minimum–maximum **pointwise** coverages over four horizons,
for nominal 95% intervals. Horizons are dependent observations of each experiment.
Raw intervals from diagnostic-failing fits remain in these calculations; their
coverage does not make those fits usable. Missing calculations are retained in
the coverage denominator, not silently discarded.

| Scenario | Outcome ITT | Adoption ITT | CACE |
| --- | --- | --- | --- |
| Linear | 91–97% | 77–78% | 97–99% |
| Smooth nonlinear | 88–97% | 86–87% | 98–100% |
| Abrupt nonlinear | 84–99% | 80–81% | 98–99% |
| Weak relevance | 82–89% | 74–76% | 95–99% |
| No relevance | 83–95% | 76–79% | Undefined target |
| Serial / heteroskedastic | 89–99% | 89–93% | 97–100% |

The primary adoption-ITT coverage is 77–87% even in the strong-stage scenarios.
For example, 77/100 has a 95% Wilson interval of [67.9%, 84.2%], so the discrepancy
is materially larger than Monte Carlo uncertainty. Near 95% coverage, the Monte
Carlo standard error with 100 experiments is about 2.2 percentage points.
High observed CACE coverage does not validate the joint reduced-form model.

Weak-stage ratio RMSE is excluded from the superiority criterion: nearly zero
sample denominators can make a regularized median look favorable against unstable
ratios without establishing reliable inference. At zero relevance, CACE is
undefined and is not assigned zero as a scoring target. Complete set components,
infinite endpoints, empty sets and undefined estimates are retained in the data.
Median widths are conditional on computable widths and are accompanied by
availability counts. Finite posterior quantiles do not establish identification.

## Correlated-intercept sensitivity

The [secondary protocol](correlated-sensitivity-protocol.md) was added after
inspection of the first 200 primary experiments. It enables correlated intercepts
on the same 100 datasets in each strong scenario, changing no other settings.
The proper inverse-Wishart prior keeps df 7 and scale 4. Seeds, truths and all
reference calculations match the primary runs. This is exploratory sensitivity
evidence; a promising result would require independent confirmation.

| Scenario | Unadjusted IV RMSE | Linear adjusted IV RMSE | Spline adjusted IV RMSE | LongBet direct RMSE | LongBet CACE coverage | LongBet diagnostic fits |
| --- | --- | --- | --- | --- | --- | --- |
| Linear | 1.195 | 0.347 | 0.371 | 0.362 | 97–99% | 86/100 |
| Smooth nonlinear | 1.341 | 0.464 | 0.409 | 0.440 | 95–99% | 95/100 |
| Abrupt nonlinear | 1.245 | 0.470 | 0.384 | 0.443 | 96–97% | 98/100 |

| Scenario | Comparator | RMSE ratio [95% paired interval] |
| --- | --- | --- |
| Linear | Linear adjusted IV | 1.043 [0.994, 1.094] |
| Linear | Spline adjusted IV | 0.976 [0.917, 1.037] |
| Smooth nonlinear | Linear adjusted IV | 0.947 [0.858, 1.045] |
| Smooth nonlinear | Spline adjusted IV | 1.075 [0.972, 1.184] |
| Abrupt nonlinear | Linear adjusted IV | 0.944 [0.880, 1.012] |
| Abrupt nonlinear | Spline adjusted IV | 1.154 [1.086, 1.228] |

| Scenario | Outcome ITT | Adoption ITT | CACE |
| --- | --- | --- | --- |
| Linear | 88–96% | 82–83% | 97–99% |
| Smooth nonlinear | 87–95% | 90–91% | 95–99% |
| Abrupt nonlinear | 82–97% | 85–88% | 96–97% |

The correlated option also fails the combined accuracy, coverage and diagnostic
criteria. The comparison therefore does not support adopting either evaluated
configuration as a better IV procedure.

## Interpretation and scope

Time smoothing trades variance against bias, and inexpensive covariate adjustment
already captures substantial precision gains in these designs. The Gaussian
adoption likelihood cannot fully represent dependence from absorbing binary
adoption. A correlated unit intercept supplies one dependence component; it does
not establish the full outcome/adoption error model. The observed coverage failure
is consistent with these working-likelihood limitations; this experiment does
not isolate a unique cause or evaluate every possible prior and model extension.

Duration 2SLS is useful when exposure history is the right treatment definition,
but independent full-dummy tests show that the repaired routine implements
conventional 2SLS. It adds no new identification guarantee. Likewise, randomization
AR and binary principal-stratum bounds are established methods with their own
assumptions; they are not evidence that the LongBet smoothing model dominates IV.

Flexible Bayesian IV already exists in [IVBART](https://arxiv.org/html/2102.01199v1).
[Generalized random forests](https://arxiv.org/html/1610.01271v4) cover heterogeneous
IV effects, and [DML](https://arxiv.org/html/1608.00060v7) supplies orthogonal scores
and cross-fitting for IV/LATE. [Dynamic exposure IV](https://arxiv.org/html/2501.01623v2)
has established stacked-2SLS alternatives under explicit assumptions.
[Weak-IV-robust inference](https://onlinelibrary.wiley.com/doi/10.1111/1468-0262.00438)
addresses weak identification under its stated conditions. These methods were
not directly benchmarked here and can target different parameters. The present
result is a failure to establish an advantage against the specified adjusted
baselines in the tested settings.

## Reproduction and complete records

The primary sources match commit `802cad26ac86f22ba59cd11860feffa556dfd896`;
the secondary driver matches `eb6e537819e3ef51b3727493017d4103018e5d4b`.
The launch manifests record the parent commits because evaluation began from
repaired working files before those commits were created. Source hashes verify
the evaluated code: see the [verification record](comparison-results/source-verification.json).
Only the separately repaired hazard and duration modules differ from the primary
launch-wide hashes; neither is called by this comparison. All 43 source hashes
match the recorded secondary revision.

Evaluation used Python 3.12.3, NumPy 2.5.2, SciPy 1.18.1, JAX 0.11.1, ArviZ 1.3.0,
bartz 0.12.1 and pandas 3.0.5. Runtime is illustrative on a shared CPU host:
each evaluation used six workers, with overlap between the two runs. Primary
strong-stage fits took about 22.5 seconds each, versus roughly 0.5–0.7 milliseconds
for adjusted reference calculations. The correlated runs incurred additional
covariance-sampling cost and resource contention; their elapsed times should not
be read as an isolated algorithmic speed ratio.

From a checkout of the secondary revision, with matching dependencies:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python benchmarks/encouragement/iv_comparison.py --output comparison-primary --replications 100 --workers 6
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python benchmarks/encouragement/iv_comparison.py --output comparison-correlated --replications 100 --scenarios linear smooth abrupt --workers 6 --correlated-intercepts
PYTHONPATH=src python benchmarks/encouragement/summarize_iv_comparison.py comparison-primary
PYTHONPATH=src python benchmarks/encouragement/summarize_iv_comparison.py comparison-correlated
```

The secondary driver's additional flag defaults to false, reproducing the primary
model's numerical configuration. Use fresh output directories when a source hash
or configuration changes. Full compressed run records, original-record hashes,
manifests, estimates, diagnostics and packet checksums are in repository directory
`benchmarks/encouragement/comparison-results/`. The [primary pointwise table](comparison-results/primary/pointwise.csv)
and [secondary pointwise table](comparison-results/correlated/pointwise.csv) retain
the complete coverage and Monte Carlo uncertainty results.
