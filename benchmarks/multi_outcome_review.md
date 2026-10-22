# Multi-outcome chapter review — 2026-09-07

**Historical results:** all JAX runs below predate the dynamic leaf-precision
cache repair, including separate fits and isolated blocked-update trials.
They are retained as the investigation record, not evidence about the repaired
sampler. The C++ runs are unaffected by that bartz-specific bug. See the
[subsequent full-precision repair](full_sur_repair.md) for current semantics,
analytical tests and fresh same-data comparisons. Also see the
[mvbcf design review](mvbcf_design_review.md) for proposed shared-tree work.

The calendar/exposure-grid defect is fixed in JAX commit `1227b28`. The remaining
problem is not resolved by that correction: plausible point estimates still
come with serious chain disagreement. This review does not certify joint
posterior calibration or recommend restoring the chapter's worked example.

## Fixed comparison protocol

- JAX sampler source: `1227b28`. Subsequent edits in this review correct
  documentation and add benchmark tools; they do not change the sampler.
- Original C++ source: [ignacio82/longbet, commit
  4414406a7a4aafdc5772c37a203076342b0be746](https://github.com/ignacio82/longbet/tree/4414406a7a4aafdc5772c37a203076342b0be746).
  This was verified against the remote HEAD and compiled into an isolated R
  library. Neither the installed book image nor the reference source was
  modified.
- Data: the archived chapter simulation, with two continuous log outcomes and
  one binary complaint outcome, 750 sellers, 18 study periods, unchanged
  sign-changing effects and correlated shocks. Data SHA-256:
  `6e4904482c10530728705f0c4c1ce8bc943db6c2f0aa6aedfba5bb99284aaffd`.
- Target: all sellers launch in week 11; evaluate week 18, exposure 8. Never
  compare a never-treated unit's exposure-zero prediction with exposure-eight
  truth. Benefits: listings >= 80 and fulfilled (59 sellers). Harms: listings
  <= 25 and not fulfilled (138 sellers).
- Continuous metrics use `expm1` **on each log-effect draw**, then average.
  Binary metrics use the draw-wise risk difference, not `Phi(tau)`.
  Group intervals and diagnostics concern draw-wise group averages, not
  averages of sellers' interval endpoints or posterior probabilities.
- JAX diagnostics use ArviZ through `longbet`; the C++ benchmark uses R's
  `posterior` package. Both report rank-normalized split R-hat and bulk/tail
  ESS. Their ESS estimators can differ slightly; this does not change the
  failure of the prespecified diagnostic gates.
- JAX short budget: 2,000 burn-in, 250 retained draws/chain, thinning 2, four
  dispersed chains. Long budget: 8,000 burn-in, 500 retained draws/chain,
  thinning 4. Both use 20 trees per forest and `lambda_knl=2`.
- JAX multi seed 314159; scalar seeds 314160/314161/314162 for GMV, hours,
  complaint. Different public APIs do not use matched random streams merely
  because they receive the same seed.
- C++ budget: four independent fits with seeds 314159–314162, 750 **total**
  sweeps each, discarding 500, retaining 250. This is well above its documented
  default of 60 total/20 burn-in, but not a proof that still longer runs could
  not help. Use 20 trees/forest, `lambda_knl=2`, numeric one-hot covariates
  (`pcat=0`), random intercepts, SUR enabled and serial execution. Other options
  retain native defaults, including `a_scaling=TRUE`, `b_scaling=FALSE`.
- This is a same-data comparison, **not identical priors, equivalent sampler
  transitions, equal iterations, or equal compute**. The reference uses
  grow-from-root trees and different tree/leaf priors. Its independent fits
  also do not use JAX's dispersed initialization scheme. Timings under
  concurrent workloads are not a performance benchmark.

The C++ prediction method omits training-unit random intercepts. For a fair
in-sample binary estimand, the benchmark adds each seller's **matched posterior
intercept draw** to the untreated latent prediction before computing both
normal CDFs. Adding a posterior mean intercept instead would lose uncertainty.
Unadjusted public predictions are retained separately in the RDS archive.
No change to the C++ fitting algorithm is made.

## JAX repair experiments

The short fixed-coding result looked promising. The longer check did not
confirm it. Reporting only the short result would have hidden a failure.

| Complaint model | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| Separate, adaptive coding, short | .01310 | 1.953 / 5.6 | 1.364 / 9.3 |
| Separate, adaptive coding, long | .01271 | 2.128 / 5.3 | 2.010 / 5.5 |
| Separate, fixed coding, short | .01243 | 1.016 / 124.3 | 1.086 / 36.6 |
| Separate, fixed coding, long | .01088 | 1.272 / 11.3 | 1.359 / 9.3 |
| Coupled, adaptive coding, short | .01089 | 2.107 / 5.3 | 1.733 / 6.2 |
| Coupled, fixed coding, long | .00990 | 1.117 / 25.9 | 1.301 / 10.5 |
| Experimental blocked coding/GP, short | .01044 | 1.609 / 6.9 | 1.903 / 5.8 |
| Experimental blocked coding/GP, long | .01494 | 1.615 / 6.8 | 2.048 / 5.4 |

Fixed coding means `adaptive_coding=False` in Python / `b_scaling=FALSE` in R:
**both** coefficients become 1. It changes the prior/model; it is not an
equivalent faster transition. Do not change the global default based on this
single realization. The full three-outcome fixed-coding comparison also failed:

| Long fixed-coding model | GMV RMSE / max group R-hat | Hours RMSE / max group R-hat | Complaint RMSE / max group R-hat |
| --- | ---: | ---: | ---: |
| Separate | .03090 / 1.730 | .02946 / 1.419 | .01088 / 1.359 |
| Coupled | .03813 / 1.356 | .02471 / 1.278 | .00990 / 1.301 |

The experimental block marginalizes the Gaussian GP trajectory when proposing
the two coding coefficients, then redraws the GP conditionally. Random-walk and
prior-independence MH proposals were tried within the same block. Its collapsed
likelihood matched a dense multivariate-normal calculation, but neither tested
budget resolved the empirical chain disagreement. This remains an isolated
scratch experiment, **not installed package functionality**. These failures do
not establish that no other block or tree proposal could help.

## Original C++ results: four completed independent chains

The reference did **not** solve this chapter example at the tested budget.
Reference results below include all four prespecified seeds, not selected
chains. The fourth was run in a separate process; the three earlier results and
that fourth result were combined after checking the seed sequence exactly.

| Outcome | Reference RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| GMV proportional change | .08485 | 1.220 / 13.2 | 1.245 / 12.1 |
| Hours proportional change | .02792 | 1.110 / 27.3 | 1.168 / 17.6 |
| Complaint risk difference | .01504 | 1.018 / 172.3 | 1.037 / 156.7 |

For comparison, the corrected JAX coupled short/default-coding RMSEs were
.04244, .01940 and .01089. These same-data scores do not establish a universal
winner: configurations and priors differ, and neither run passes the required
diagnostic checks.

The reference's complaint benefit-region truth is **-2.814 percentage points**;
its mean estimate is -0.210 points, 95% interval [-2.969, 2.023], with 55.4% of
draws having the correct sign. The harm-region truth is **+2.186 points**;
its estimate is +0.417 points, interval [-1.438, 3.278], with 63.8% of draws
having the correct sign. These intervals and frequencies describe the sampled
output, not certified posterior uncertainty. The mean estimates substantially
attenuate the signed heterogeneity the example is supposed to demonstrate.

The full reference archive is `reference-complete.rds` in the scratch directory;
`reference-summary.csv` contains all subgroup scores and bulk/tail ESS;
`reference-draws.npz` contains the combined business-scale draws. Its smaller
R-hats on the binary outcome than JAX's do not establish accurate effect
recovery or calibrated joint inference. A longer C++ run, different priors or
another dataset could behave differently; this test does not rule those out.

## Statistical limitation and next repair

Both implementations use one-way recursive residual coupling. The reference's
[`mcmc_loop_multi`](https://github.com/ignacio82/longbet/blob/4414406a7a4aafdc5772c37a203076342b0be746/src/mcmc_loop.cpp)
constructs incoming offsets from earlier outcomes, skips them for binary
equations, and does not include downstream likelihood terms in earlier mean
updates. Therefore continuous outcomes cannot improve its binary fit through
this coupling. Pairing sweep numbers is not evidence of a calibrated joint
posterior. Matching that code cannot by itself repair this limitation.

A substantive repair needs mean and binary-latent conditionals using the full
error precision, separated from structural innovation-variance updates. With
`r=y_star-f`, `B=I-Gamma`, and `Omega=B' diag(1/v) B`, the conditional residual
for outcome `m` is `(Omega r)[m]/Omega[m,m]`, with conditional variance
`1/Omega[m,m]`. Missing cells require their own marginal/conditional treatment.
This change needs a new semantic/persistence version and analytical,
simulation-calibration, scalar-equivalence and convergence tests. The detailed
engine-development checklist is in `/home/ignacio/book/multi.md`, section 0.6.
It is a sampler redesign, not a correction that belongs inside chapter code.

Even an exact conditional update is not a guarantee of practical mixing or
better marginal RMSE on every outcome. Keep posterior correctness, chain
convergence, recovery and repeated-simulation calibration as separate checks.

## Reproduction and validation

The new scripts are deliberately independent of chapter fit caches:

```sh
# In an R environment with the JAX front door, reticulate and chapter libraries:
Rscript benchmarks/export_multi_chapter.R /home/ignacio/book audit-input.npz .

python benchmarks/bench_multi_chapter.py --data audit-input.npz \
  --output audit/joint-fixed-long --mode multi --coding fixed \
  --burnin 8000 --draws 500 --skip 4
python benchmarks/bench_multi_chapter.py --data audit-input.npz \
  --output audit/separate-fixed-long --mode scalar --coding fixed \
  --burnin 8000 --draws 500 --skip 4

# Install the PINNED reference checkout into a DIFFERENT R library first.
# RETICULATE_PYTHON must identify an interpreter with NumPy.
Rscript benchmarks/bench_multi_reference.R audit-input.npz \
  audit/cpp-reference /path/to/reference-r-library
```

The reference package and JAX R front door are both named `longbet`; never rely
on that name alone to identify the engine. The scripts check the namespace or
load an explicit source/library. Benchmark outputs refuse to overwrite prior
result files. Re-exporting produced exactly the original 13 arrays and the
same NPZ hash. The Python benchmark completed a short two-chain end-to-end
smoke run; its two summary regression tests passed. The new R scripts and
edited R documentation parsed successfully. The experimental transition was
not merged, and no full production test rerun is implied by these checks.

Scratch fits, compact draws and experimental scripts are retained under
`/tmp/longbet-reference-benchmark.9G6Fop`; original reviewed fits are under
`/tmp/longbet-multi-audit.xg5FSO`. These paths are audit artifacts, not runtime
dependencies. The synthetic input can be regenerated using the exporter.
