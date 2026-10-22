# Direct encouragement model repair comparison

These results evaluate mixing on a bounded experiment. They do not establish
frequentist calibration or authorize a production default change. See the
[protocol](../../docs/encouragement-repair-protocol.md).

The following completed settings pass every effect, parameter and ESS-growth
screen in all four runs:

- `direct_joint`, `clean_strong`, 8,000 retained draws per chain.
- `parametric`, `clean_strong`, 2,000 retained draws per chain.
- `direct_joint`, `strong`, 8,000 retained draws per chain.
- `parametric`, `strong`, 2,000 retained draws per chain.

Passing these screens supports further calibration work for the named
settings. The original forest, baseline-only block, and shorter-chain
failures remain below.

Every run uses four dispersed chains and 1,000 warmup sweeps. Retained
draw counts are shown per chain; extensions preserve their shorter prefixes.
There are two independent datasets and two sampler seeds per dataset in each
completed setting. Both ITTs and the raw Wald ratio are checked at all four
post-encouragement horizons in six-period panels with 80 study units. The
tail ESS refers to the 2.5% and 97.5% quantiles. Failed runs and
diagnostic failures are retained; repeated seeds do not count as new datasets.

| Scenario | Model | Draws | Effect passes | Failed / pending | Max R-hat | Min bulk ESS | Min tail ESS | Median fit seconds | Min bulk ESS/sec |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| clean_strong | current | 2000 | 0/4 | 0 / 0 | 1.3277 | 10 | 56 | 8.6 | 0.8 |
| clean_strong | current_stump | 2000 | 0/4 | 0 / 0 | 1.0378 | 98 | 518 | 8.1 | 12.3 |
| clean_strong | direct_blocked | 2000 | 4/4 | 0 / 0 | 1.0011 | 2721 | 4810 | 35.4 | 76.7 |
| clean_strong | direct_joint | 2000 | 4/4 | 0 / 0 | 1.0010 | 5539 | 6809 | 42.5 | 130.6 |
| clean_strong | direct_joint | 4000 | 4/4 | 0 / 0 | 1.0009 | 10941 | 13754 | 71.4 | 152.4 |
| clean_strong | direct_joint | 8000 | 4/4 | 0 / 0 | 1.0005 | 21366 | 28396 | 131.7 | 162.3 |
| clean_strong | direct_smooth | 2000 | 4/4 | 0 / 0 | 1.0023 | 2160 | 4831 | 27.4 | 70.3 |
| clean_strong | parametric | 2000 | 4/4 | 0 / 0 | 1.0006 | 7229 | 7016 | 1.4 | 5005.2 |
| strong | direct_joint | 8000 | 4/4 | 0 / 0 | 1.0005 | 21116 | 26968 | 137.2 | 153.9 |
| strong | parametric | 2000 | 4/4 | 0 / 0 | 1.0006 | 7221 | 7016 | 1.4 | 5009.5 |

## Other posterior diagnostics and longer windows

The first half of retained draws is compared with its enclosing full
window. The 1.5-fold ESS growth flag is a noisy screening criterion. Parameter
checks include every unit intercept, both equation variances, and, for the
new candidates, aggregate baseline/effect functions at three fixed original
covariate rows. Direct-forest checks also include counts of root trees and
each split feature across the forest, which do not depend on tree labels.
Individual labeled tree components are archived for audit;
exchangeability means their identities are not substantive estimands.

| Scenario | Model | Draws | Parameter passes | Growth passes | Min ESS growth | Median ESS growth |
|---|---|---:|---:|---:|---:|---:|
| clean_strong | current | 2000 | 0/4 | 0/4 | 0.18 | 1.67 |
| clean_strong | current_stump | 2000 | 0/4 | 0/4 | 0.75 | 1.85 |
| clean_strong | direct_blocked | 2000 | 0/4 | 4/4 | 1.67 | 1.98 |
| clean_strong | direct_joint | 2000 | 1/4 | 4/4 | 1.60 | 1.99 |
| clean_strong | direct_joint | 4000 | 3/4 | 4/4 | 1.71 | 2.01 |
| clean_strong | direct_joint | 8000 | 4/4 | 4/4 | 1.82 | 2.00 |
| clean_strong | direct_smooth | 2000 | 0/4 | 4/4 | 1.66 | 1.98 |
| clean_strong | parametric | 2000 | 4/4 | 4/4 | 1.81 | 2.00 |
| strong | direct_joint | 8000 | 4/4 | 4/4 | 1.72 | 2.00 |
| strong | parametric | 2000 | 4/4 | 4/4 | 1.81 | 2.00 |

`seed-agreement.json` compares independent sampler seeds on each fixed
dataset using ITT posterior means and Wald medians, with their combined
Monte Carlo SEs. These checks are meaningful only for adequately mixed
draws; they do not turn two datasets into four calibration replications.

## Interpretation and limitations

The current-stump ablation limits tree depth while retaining the current
product/coding model. New and old models share variance hyperparameters,
but their function priors and topology updates differ. These experiments
therefore compare complete modeling/sampling choices; they cannot isolate a
single causal explanation for a performance difference. New forests contain
sampled stumps with smooth time-vector leaves and cannot express general
covariate interactions. The parametric model has linear baseline covariates.

Matching average pointwise function variance does not match the prior
variance of an all-unit mean contrast, which also depends on covariance
across units. The unblocked and blocked direct forests do share exactly the
same statistical model: their differences are exact joint Gaussian updates
of baseline leaves (`direct_blocked`) or both forests' leaves (`direct_joint`)
and unit intercepts given the sampled topologies.

All models here use independent Gaussian equations, including an LPM working
likelihood for binary uptake. They do not model a correlated latent binary
system or an absorbing adoption hazard. Better mixing does not repair
likelihood misspecification, establish interval coverage, or identify a ratio
when its population denominator is zero. No probability-range restriction or
denominator clipping has been introduced.

Fit timings include JAX compilation when incurred; later fits in a process
may reuse cached compilation. Candidate fit timing includes integrated
effect/probe evaluation and construction of parameter arrays; old-model
timing excludes its separate standardization step. All fit timings exclude
NPZ serialization and post-fit diagnostic calculations. Total wall times
are also retained. Backend and bookkeeping differences limit comparisons. Raw
draws, full candidate states, current-model archives, source hashes, exact
seeds, configurations, truths and warnings accompany each result locally.
Compact JSON files and this report are reproducible from the scripts.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/repair_comparison.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/extend_repair.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/extend_repair.py --from-draws 4000 --to-draws 8000
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/repair_comparison.py --variants parametric --scenario strong
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/repair_comparison.py --variants direct_joint --scenario strong --draws 8000
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/build_repair_report.py
```
