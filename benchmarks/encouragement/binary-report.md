# Binary encouragement benchmark: implementation and results

The joint binary model is implemented and tested as a **research benchmark**. It captures dependence between the outcome and adoption ITTs, supports an explicit monotonicity assumption, and separates identified sets from point posteriors. A marginal uptake update substantially improves the monotone sampler in these experiments. The experiments also expose material weak-stage bias from the monotone prior. This is not evidence to make monotonicity automatic or to declare the model calibrated.

## Implemented behavior

- [binary_forest.py](binary_forest.py): six proper probit forests with sampled finite stumps, fixed unit utility variance, prior-dispersed chains, exact Gaussian leaf blocks, optional monotone adoption, optional marginal elliptical slice updates, and exact checkpoint continuation.
- [binary_cells.py](binary_cells.py): exact independent posterior draws for a joint finite-cell model, independent reduced-form control, ordered uptake control, and an explicitly null uptake model. It also implements sharp bounds under stated limits on direct encouragement risk differences.
- [binary_dgp.py](binary_dgp.py) and [binary_benchmark.py](binary_benchmark.py): potential responses generated before complete random assignment, separate seed streams, full replay archives, failures, ITT covariance, interval coverage, sensitivity summaries, and existing Fieller/randomization AR references.

The observed likelihood is `p(D|Z,X) p(Y|D,Z,X)`. Write `p_z=Pr(D=1|Z=z,X)` and `q_dz=Pr(Y=1|D=d,Z=z,X)`; the outcome reduced form is `r_z=p_z*q_1z+(1-p_z)*q_0z`. Both ITTs average over the original empirical X distribution. The regression on observed adoption is not itself causal. The forest supports baseline covariates but this calibration uses three prespecified baseline categories; deeper trees and longitudinal adoption histories are outside this experiment.

## Validation and scope

**112 focused tests passed**, including prior-preservation checks for both uptake parameterizations and both sampler variants; deliberately broken Gaussian and elliptical slice updates fail those checks. Tests cover extreme truncated-normal tails, the joint monotone augmentation, Gaussian block moments, exact checkpoint replay, empty outcome cells, residual invariants, analytic Beta covariance, and agreement with an independent linear program for every type of bound. A constant-X experiment independently compares the full forest posterior to exact Beta/ordered-Beta draws. Enumeration verifies the finite-population covariance identity.

The main study has **600 independent datasets: 100 in each of six scenarios, N=240**. Each cell fit uses four independent streams of 1,000 posterior draws. Two additional prior settings reuse the same 200 weak/zero datasets. The forest pilot reuses eight datasets (two strong, two weak, and one per remaining scenario), with four chains, 1,000 burn-in and 4,000 retained draws per chain. These eight datasets cannot establish forest coverage. No cell coverage rate is attributed to the forest.

All point intervals below have nominal 95% posterior/reference coverage. Rates show the Monte Carlo standard error and the number of available datasets; JSON also contains Wilson intervals. Conditional targets average potential-response probabilities at fixed empirical X. Separate finite-population targets retain realized potential outcomes. Estimates following failed fits are not available; successful-fit coverage must be read together with failure counts.

## Main cell results

| Scenario | Model | Completed/100 | Outcome ITT coverage | Adoption ITT coverage | Raw Wald coverage |
|---|---|---:|---:|---:|---:|
| strong | cells_joint | 100 | 94.0% ± 2.4 pp (100) | 98.0% ± 1.4 pp (100) | 96.0% ± 2.0 pp (100) |
| strong | cells_independent | 100 | 93.0% ± 2.6 pp (100) | 98.0% ± 1.4 pp (100) | 95.0% ± 2.2 pp (100) |
| strong | cells_monotone | 100 | 95.0% ± 2.2 pp (100) | 97.0% ± 1.7 pp (100) | 95.0% ± 2.2 pp (100) |
| strong | reference normal/Fieller | 100 | 94.0% ± 2.4 pp (100) | 97.0% ± 1.7 pp (100) | 94.0% ± 2.4 pp (100) |
| weak | cells_joint | 100 | 97.0% ± 1.7 pp (100) | 95.0% ± 2.2 pp (100) | 98.0% ± 1.4 pp (100) |
| weak | cells_independent | 100 | 97.0% ± 1.7 pp (100) | 95.0% ± 2.2 pp (100) | 99.0% ± 1.0 pp (100) |
| weak | cells_monotone | 99 | 96.0% ± 2.0 pp (99) | 66.7% ± 4.7 pp (99) | 96.0% ± 2.0 pp (99) |
| weak | reference normal/Fieller | 100 | 97.0% ± 1.7 pp (100) | 94.0% ± 2.4 pp (100) | 97.0% ± 1.7 pp (100) |
| zero | cells_joint | 100 | 97.0% ± 1.7 pp (100) | 95.0% ± 2.2 pp (100) | undefined |
| zero | cells_independent | 100 | 96.0% ± 2.0 pp (100) | 95.0% ± 2.2 pp (100) | undefined |
| zero | cells_monotone | 98 | 93.9% ± 2.4 pp (98) | 0.0% ± 0.0 pp (98) | undefined |
| zero | reference normal/Fieller | 100 | 96.0% ± 2.0 pp (100) | 93.0% ± 2.6 pp (100) | undefined |
| one_sided | cells_joint | 100 | 95.0% ± 2.2 pp (100) | 92.0% ± 2.7 pp (100) | 93.0% ± 2.6 pp (100) |
| one_sided | cells_independent | 100 | 95.0% ± 2.2 pp (100) | 92.0% ± 2.7 pp (100) | 93.0% ± 2.6 pp (100) |
| one_sided | cells_monotone | 100 | 95.0% ± 2.2 pp (100) | 93.0% ± 2.6 pp (100) | 93.0% ± 2.6 pp (100) |
| one_sided | reference normal/Fieller | 100 | 95.0% ± 2.2 pp (100) | 96.0% ± 2.0 pp (100) | 95.0% ± 2.2 pp (100) |
| direct_effect | cells_joint | 100 | 98.0% ± 1.4 pp (100) | 89.0% ± 3.1 pp (100) | 97.0% ± 1.7 pp (100) |
| direct_effect | cells_independent | 100 | 96.0% ± 2.0 pp (100) | 89.0% ± 3.1 pp (100) | 98.0% ± 1.4 pp (100) |
| direct_effect | cells_monotone | 100 | 97.0% ± 1.7 pp (100) | 89.0% ± 3.1 pp (100) | 96.0% ± 2.0 pp (100) |
| direct_effect | reference normal/Fieller | 100 | 96.0% ± 2.0 pp (100) | 94.0% ± 2.4 pp (100) | 94.0% ± 2.4 pp (100) |
| defiers | cells_joint | 100 | 96.0% ± 2.0 pp (100) | 96.0% ± 2.0 pp (100) | 97.0% ± 1.7 pp (100) |
| defiers | cells_independent | 100 | 97.0% ± 1.7 pp (100) | 96.0% ± 2.0 pp (100) | 98.0% ± 1.4 pp (100) |
| defiers | cells_monotone | 8 | 100.0% ± 0.0 pp (8) | 0.0% ± 0.0 pp (8) | 100.0% ± 0.0 pp (8) |
| defiers | reference normal/Fieller | 100 | 95.0% ± 2.2 pp (100) | 95.0% ± 2.2 pp (100) | 95.0% ± 2.2 pp (100) |

Wald coverage in the direct-effect and defier cases concerns the **ratio of randomized contrasts**, not an identified CACE. The null-stage model intentionally assumes equality of uptake across assignment; its adoption ITT is exactly zero and its Wald ratio remains undefined. That model has 100% zero-stage ITT coverage by construction and 0% in nonzero-stage scenarios.

### Covariance

A coherent binary joint posterior induces ITT covariance through shared uptake probabilities. In the unrestricted cell model, `Cov(r_z,p_z)=Var(p_z)*(E(q_1z)-E(q_0z))` exactly. Independent reduced-form posteriors set that covariance to zero. Each estimator is compared with its own empirical error covariance, because covariate adjustment changes variance.

| Scenario | Joint posterior covariance | Joint empirical error covariance ± MCSE | Independent posterior covariance |
|---|---:|---:|---:|
| strong | 0.000718 | 0.000285 ± 0.000256 | 0.000007 |
| weak | 0.000850 | 0.001124 ± 0.000313 | -0.000001 |
| zero | 0.000907 | 0.001103 ± 0.000265 | 0.000001 |
| one_sided | 0.000437 | 0.000117 ± 0.000227 | -0.000001 |
| direct_effect | 0.000888 | 0.000850 ± 0.000363 | 0.000011 |
| defiers | 0.000927 | 0.001284 ± 0.000358 | -0.000004 |

Finite-sample covariance estimates still have visible error, particularly in the one-sided scenario. The joint construction repairs a missing source of dependence; it does not force posterior covariance to equal repeated-sampling covariance. The independent and joint outcome priors also differ, so their interval differences are not a sampler comparison.

### Monotonicity and sensitivity to priors

Ordering uptake probabilities encodes a no-defiers assumption; it cannot establish that assumption from the data. Continuous ordered priors assign no mass to an exactly null first stage. Consequently their equal-tail credible intervals exclude zero even when zero is the truth. The weak-stage bias below is substantial and is not repaired by longer chains.

| Scenario | Symmetric Beta prior parameter | Monotone adoption-ITT bias | Adoption-ITT coverage | Rejection-budget failures |
|---|---:|---:|---:|---:|
| weak | 0.5 | 0.0602 | 62.6% ± 4.9 pp (99) | 1 |
| zero | 0.5 | 0.0786 | 0.0% ± 0.0 pp (97) | 3 |
| weak | 1 | 0.0593 | 66.7% ± 4.7 pp (99) | 1 |
| zero | 1 | 0.0780 | 0.0% ± 0.0 pp (98) | 2 |
| weak | 2 | 0.0575 | 68.7% ± 4.7 pp (99) | 1 |
| zero | 2 | 0.0770 | 0.0% ± 0.0 pp (100) | 0 |

The main monotone cell sampler exhausted its 500,000-candidate-per-cell budget in 1 weak, 2 zero, and 92 defier datasets. These are computational failures of exact rejection sampling; they are not a test proving defiers and do not mean that the constrained posterior does not exist. The defier scenario's eight successful constrained fits are highly selected. The forest product prior differs from the ordered-Beta prior, so the numerical cell bias cannot be transferred directly to the forest.

## Forest mixing and the added marginal update

The original monotone data augmentation mixes slowly when one latent probability approaches a boundary. The optional `marginal_uptake_refresh=True` step uses [elliptical slice sampling](https://proceedings.mlr.press/v9/murray10a.html) for uptake leaf coefficients with the observed adoption likelihood, integrating the auxiliary labels and utilities out. The next sweep regenerates those auxiliaries. Tree priors, probability parameterization and the posterior remain unchanged.

The table reports all sampled effect/function/leaf/topology checks. Required thresholds are R-hat ≤1.01 and bulk/tail ESS ≥400. Undefined sensitivity endpoints are labeled separately, not counted as successful diagnostics or discarded draws.

| Dataset | Original monotone failing quantities | Refreshed failing quantities | Adoption-ITT bulk ESS: original → refreshed | R-hat: original → refreshed |
|---|---:|---:|---:|---:|
| defiers / 0 | 5 | 0 | 298 → 2228 | 1.0149 → 1.0015 |
| direct_effect / 0 | 0 | 0 | 1501 → 5944 | 1.0034 → 1.0004 |
| one_sided / 0 | 4 | 0 | 1914 → 6682 | 1.0030 → 1.0005 |
| strong / 0 | 0 | 0 | 1811 → 6791 | 1.0007 → 1.0005 |
| strong / 1 | 0 | 0 | 1519 → 5956 | 1.0030 → 1.0006 |
| weak / 0 | 4 | 0 | 305 → 3342 | 1.0076 → 1.0024 |
| weak / 1 | 4 | 0 | 265 → 3283 | 1.0143 → 1.0018 |
| zero / 0 | 4 | 0 | 245 → 3799 | 1.0125 → 1.0010 |

Unrestricted forests pass all sampled checks in 8/8 pilot fits. With the marginal update, monotone forests pass in 8/8 fits. This is evidence about the sampler on these datasets, not a forest coverage study. Both original and refreshed records remain in the archive, including failures of the original diagnostics.

## Identification and direct encouragement effects

Under declared no-defiers, the code bounds three distinct contrasts: the joint encouragement/adoption effect among compliers (`Y(1,1)-Y(0,0)`), and treatment effects holding encouragement fixed (`Y(1,0)-Y(0,0)` and `Y(1,1)-Y(0,1)`). Bounds use complier masses before population aggregation, avoiding tiny cell-level denominators. A scalar delta, or separate always/never/complier deltas, limits absolute direct changes in conditional outcome probabilities. Delta=1 leaves direct effects unrestricted; delta=0 imposes exclusion and collapses compatible sets to the Wald contrast.

The unrestricted observed-data posterior is not an exclusion-restricted posterior. Some draws produce empty sets under stronger restrictions; the output records their fraction and withholds an unconditional endpoint summary whenever any draw is incompatible or has zero complier mass. No probability density is invented inside a set. Quantiles of bounds and their outer envelope are not labeled a 95% point-effect posterior or a confidence region with established coverage.

In this main study, no monotone cell dataset had every draw compatible with delta ≤0.1. All successfully sampled monotone fits have nonempty bounds for delta=1. The compatibility fraction is a property of an assumed-model transformation, not the posterior probability that an instrument is valid. Weak stages naturally yield broad bounds.

This follows the useful observed-mixture and monotone-probit ideas in the [StochTree IV vignette](https://stochtree.ai/vignettes/iv.html), while keeping unidentified effects as sets unless an additional sensitivity prior is explicitly supplied.

## Reference inference

The existing reference uses the joint Neyman covariance and analytical normal/Fieller inversion. The harness also evaluates randomization AR at the true finite and conditional Wald ratios with 499 Monte Carlo assignments and the plus-one p-value. It stores Monte Carlo intervals for the ideal enumeration tail probability. These DGPs have heterogeneous effects: finite-sample sharp-null validity does not make AR exact for their average ratios, and acceptance at supplied values is not an inverted full confidence set. Zero-stage truth leaves the ratio undefined. No causal CACE coverage claim is made in the direct-effect or defier cases.

## Reproduction

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/binary_benchmark.py --output /tmp/binary-cells --replicates 100 --draws 1000
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/encouragement/binary_benchmark.py --output /tmp/binary-forest --models forest_monotone --scenarios strong weak --replicates 2 --draws 4000 --burnin 1000 --marginal-uptake-refresh
pytest tests/test_encourage_binary.py tests/test_encourage_direct_smooth.py tests/test_encourage_randomization_ar.py tests/test_encourage.py -q
python benchmarks/encouragement/build_binary_report.py
```

Use a fresh output directory when source or configuration changes. Runners refuse to mix versions. `binary_results/*/manifest.json` records settings, package versions and source hashes. `binary_results/initial_source/` preserves the sources used before the optional marginal update. Raw NPZ files stay locally available and are ignored by Git; compact JSON records are versioned. Checkpoints can be continued with `fit_binary(..., burnin=0, resume=archive)` using unchanged data and configuration. `cells_null` is an explicit assumption, not automatic model selection.

## Practical assessment

**Benefits:** binary probabilities stay coherent; the model preserves ITT dependence; empirical-X standardization is explicit; there are independent exact sampler controls; the optional marginal step improves weak-stage sampling; and sensitivity analysis makes the unidentified part visible.

**Costs and limits:** monotonicity can bias weak first stages and exclude the exact null; exclusion and no-defiers remain substantive assumptions; some sensitivity sets are empty or broad; the exact ordered-cell sampler can be expensive; the forest is a finite-stump research model; and larger forest coverage studies plus longitudinal-history modeling remain necessary before considering a public backend.
