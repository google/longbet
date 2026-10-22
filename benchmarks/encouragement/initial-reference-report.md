# Reference encouragement inference: initial calibration (preserved first delivery)

Date: 2026-09-09. This report covers `encouragement_effects`, the design-based
reference utility. It does not validate a Bayesian encouragement wrapper or
general staggered, blocked, or cluster-randomized designs.

## Reproduce

From an environment with this checkout installed:

```bash
python benchmarks/encouragement/reference_calibration.py --replications 500
python benchmarks/encouragement/reference_calibration.py --replications 5000
```

The first command's output is retained in `reference_calibration_500.json`.
The extended run is in `reference_calibration_results.json`. Both use 400 units,
seed 2026 (incremented by regime), and the same fixed populations and assignment
streams. The 500-assignment run is the prefix of the 5,000-assignment run, not an
independent replication. No estimator, DGP, or cutoff was tuned after inspecting
coverage. The extended run investigates the first run's low coverage estimates.

Environment: Python 3.12.3, NumPy 2.5.2, SciPy 1.18.1, pandas 3.0.5. No JAX
sampler is run. Timings are seconds to tens of seconds, depending on replication
count and other workloads.

## Design and target

For each regime, hold the finite population and potential outcomes fixed and
repeat complete randomization of half the units. Always-takers are 15% of the
population. Complier shares are 40%, 5%, and zero. Take-up is immediate for
encouraged compliers, and always-takers adopt before encouragement. Treatment
effects vary with a covariate and duration; baseline outcome levels are correlated
with compliance type. There are two post-encouragement observations, with all
panels complete. There are no defiers or exclusion violations.

Coverage is pointwise over random assignments. The two horizons are correlated
and are not counted as independent replications. Each reported Monte Carlo SE
uses 5,000 independent assignments. In the zero-first-stage population no true
Wald ratio exists, so Wald coverage is correctly left undefined.

## Results and limitations

Nominal confidence level: 95%.

| First stage | Outcome ITT coverage, h=1 / h=2 | Adoption ITT coverage | Wald set coverage, h=1 / h=2 |
| --- | --- | --- | --- |
| 40% | 94.72% / 94.84% | 98.16% / 98.16% | 94.44% / 94.56% |
| 5% | 95.06% / 95.26% | 95.82% / 95.82% | 95.02% / 95.06% |
| 0% | 95.04% / 95.12% | **92.76% / 92.76%** | Undefined population ratio |

Wald coverage Monte Carlo SEs are about 0.31–0.32 percentage points. The strong
case's deficit is about 1.4–1.7 Monte Carlo SEs; these measurements are useful
checks, not a universal 95% guarantee. The original 500-assignment run reported
92.6% Wald coverage in the strong case; both runs are preserved.

The null first-stage undercoverage is a real limitation of the discrete normal
approximation, not simply Monte Carlo noise. With 60 always-takers among 400
units and 200 encouraged units, the encouraged-adopter count `K` has distribution
`Hypergeometric(400, 60, 200)`. The nominal normal first-stage interval includes
zero only for `K=24,...,36`. Summing those exact randomization probabilities gives
**93.1945% coverage**, consistent with the measured 92.76%. This confirms the
implementation's formula while documenting that its nominal interval can
undercover even at this sample size. Do not describe these intervals as exact or
uniformly calibrated. Future permutation/score inference should be evaluated as
an explicit method extension rather than changing a cutoff to pass this seed.

For the weak case, across the two horizons and 5,000 assignments, 6,132 sets
were all-real, 1,278 were disjoint, and 2,590 were bounded. These are dependent
horizon counts, not an effective sample size of 10,000. The sample first stage
was sometimes negative despite monotonicity. The implementation retains both
facts rather than forcing bounded intervals or rejecting the data.

For the zero-stage case, some estimated confidence sets were bounded by chance.
Their shape does not establish the existence of a population CACE. Statistical
first-stage diagnostics cannot certify the IV identifying assumptions.

## What remains

This is one fixed population per regime, not the full validation matrix in the
[implementation plan](../../docs/encouragement-implementation-plan.md). Additional
sample sizes, populations, adoption delays and acceleration, binary outcomes,
heteroskedasticity, and alternative inferential methods remain benchmark work.
The first-stage undercoverage above remains an explicitly documented limitation
of `normal_ar` and the normal summary intervals.
