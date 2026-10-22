# Shared/private treatment forests: implementation review (2026-09-07)

## Status and scope

The mvbcf-inspired proposal is implemented as an **opt-in experimental LongBet
model**, retaining this package's panel, exposure GP and mixed probit/Gaussian
likelihood. It does not replace the engine with mvbcf or stochtree. Correctness
tests and chapter recovery checks are separate forms of evidence; neither a
passing unit test nor lower RMSE certifies usable posterior probabilities.

The default remains `num_shared_trees=0` and `full_precision_sur_v1`. The new
path uses `shared_private_sur_v1`. This review follows the
[full-precision/cache repair](full_sur_repair.md) and
[mvbcf design review](mvbcf_design_review.md); their unsuccessful experiments
remain part of the record.

## Public model and prior

Python's `LongBetConfig` and R's `longbet_multi()` expose:

```
num_trees_trt = 20              # TOTAL treatment trees per outcome
num_shared_trees = 10           # replaces 10 of the 20 private trees
shared_variance_fraction = 0.5  # fixed prior allocation, not a learned weight
```

Sharing is disabled by default. Counts must be integers with
`0 <= num_shared_trees < num_trees_trt`; the fraction must be finite and strictly
between zero and one. Scalar fitting rejects nonzero sharing. With sharing on,
let `Js` be the shared count, `Jp` the remaining private count and `f` the fixed
fraction. For each outcome:

```
nu_m(x) = sum_j private_leaf_jm(x) + sum_j shared_leaf_jm(x)
private_leaf ~ N(0, (1-f)/Jp)
shared_leaf_vector ~ N(0, (f/Js) I_M)
```

This preserves unit marginal prior variance of the total treatment forest
at a point, in internal standardized continuous/identified binary-latent units.
It does not equate shrinkage in dollars and probabilities, or fix the variance
of `b_Z beta_S nu_m`, whose other factors retain their own priors. The allocation
is not estimated from observed response correlations or simulated truth.

Shared partitions pool evidence about modifier locations, not signs. The same
leaf can increase GMV and decrease hours and complaints. Private trees permit
outcome-specific departures; **they do not guarantee freedom from negative
transfer**. Baseline forests, GPs, coding scales and unit intercepts remain
outcome-specific. `sur=False` removes residual coupling, but shared partitions
still make this a joint mean model. To recover independent mean models, disable
both SUR and sharing.

## Transition kernel

For raw residuals `r=y_star-f`, define `B=I-Gamma`. Under supported masks:

```
Omega_i = B' diag(observed[:,i] / innovation_variances) B
D_i = diag(b_m(Z_i) beta_m(S_i))
P_leaf = (Js/f) I + sum_i D_i' Omega_i D_i
h_leaf = sum_i D_i' Omega_i partial_raw_i
leaf | rest ~ N(solve(P_leaf,h_leaf), inverse(P_leaf))
log_integral = .5 (logdet(prior_precision) - logdet(P_leaf)
                  + h_leaf' solve(P_leaf,h_leaf))
```

The response-only quadratic cancels in topology comparisons. All outcome
weights, off-diagonal precision terms and determinant contributions are used.
Leaf sufficient statistics use indexed scatter reductions, not a dense array
over every observation and possible leaf. Conditional draws use Cholesky and
triangular solves, not an explicitly computed posterior inverse.

Each shared tree proposes grow or prune with probability one half, including
null boundary moves. Grow selects uniformly among actual non-bottom leaves,
then among globally split-enabled variables and that variable's global
cutpoints. Infeasible splits are **rejected, never redrawn**. Prune selects
uniformly among nodes whose two children are leaves. The acceptance ratio
includes the full depth/rule prior and forward/reverse selection probabilities.
Minimum child sizes count the fixed union of observed outcome cells. The
tree prior is conditioned on this fixed-design feasibility constraint.

This global rule prior differs from bartz's ancestor-restricted private-tree
rule prior, even though both use the configured depth parameters. Treatment
depth defaults `.25, 3` already matched the inspected mvbcf defaults; they were
not introduced as an additional repair. The new kernel was derived here, not
ported on the assumption that another package's acceptance calculation is
correct. Its tiny-tree posterior frequencies are tested independently.

A sweep updates outcome-local mean/latent/private blocks using the existing
full-precision conditionals, then backfits the shared trees jointly, then
updates SUR loadings and innovation variances. The beta–nu ridge move rescales
both ensembles for one outcome together. Its prior ratio includes both sets
of leaf values, and its Jacobian counts all **active** private/shared leaves.
The fitted likelihood and raw residuals remain invariant under this move.

The private-tree fit cache on the shared path is refreshed from stored leaves,
with deferred prunes respected. Recovering a small change by subtracting two
large weighted residuals caused measurable cancellation when GP/coding weights
were tiny. The direct leaf evaluation fixes this without changing the default
path's random schedule or transition behavior.

## Numerical precision and compatibility

An existing R fixture has two almost perfectly anticorrelated continuous
outcomes. Float32 matrix accumulation rounded away the regularizing prior and
caused NaN Cholesky factors. Float64 accumulation/factorization now preserves
that prior; there is no added arbitrary diagonal jitter or modified fixture.
Stored scalar states, vector leaves and traces remain float32. Non-finite
retained shared traces raise an explicit error rather than yielding an unusable
fitted object.

The shared MCMC driver holds JAX's X64 context through allocation, outer
tracing, batching, lowering and execution. Merely enabling it inside a matrix
helper was insufficient and produced dtype/lowering failures. Scalar bartz
initialization and local steps explicitly retain their float32 context; the
caller's global dtype setting is restored. Tested environment: **JAX 0.11.1,
bartz 0.12.1**, CPU, both Python and source-mounted R integration. Compatibility
with all older allowed dependency versions or accelerator backends has not
been empirically established. Distributed data meshes are explicitly rejected
on the shared path; single-device parallel chains are supported.

Shared multi archives use format 2, retain `precision_cache_version=1`, and
record counts, variance fraction, semantics and `private_then_shared` layout.
Each outcome's prediction trace concatenates private trees and its column of
shared leaves. Loading validates equal shared split/variable suffixes across
outcomes as well as axes, order and metadata. Extracted scalar children preserve
joint-fit provenance and make self-contained marginal predictions; they are
not independent fits and cannot reconstruct a joint posterior by being paired
arbitrarily. Old invalid precision-cache archives still require refitting.

## Correctness and regression evidence

The shared tests include independent dense observed-likelihood matrix and
integrated-density checks; actual vector-leaf Gaussian mean/full covariance;
nine-tree enumerated posterior frequencies; mixed probit/Gaussian observed-data
quadrature and a joint sign event; SUR-on/off mixed-panel residual identities;
prior variance accounting; dtype restoration; near-collinear precision;
parent/child persistence; malformed shared topology rejection; and an
independent ridge prior-ratio/Jacobian and likelihood-invariance check.

The original 17 shared tests passed in 134.04 seconds, and the added ridge
check passed in 29.68 seconds. The full R suite passed **268 assertions** in
263.0 seconds with no failures, warnings or skips, including fresh-process
shared parent/child persistence. All 42 chapter R chunks parse, and the setup,
engine guards and source-based cache fingerprint execute in a disposable,
source-mounted `book` container. The image was not rebuilt and the full
chapter was not rendered. A complete Python-suite run passed **216 tests**, with
one optional external-BCF comparison skipped, in 1,110.03 seconds. Four expected
warnings concern always-treated units whose random intercepts are inseparable
from treatment effects; there were no failures. The record is `pytest-full.xml`.
The existing high-signal, sign-changing two-continuous/one-binary statistical test
also passed with sharing on (68.18 seconds), preserving its data, seed, budget
and thresholds. With five shared and five private treatment trees, sign
accuracy was .988/.996/.988 and RMSE .146/.059/.086 for GMV/hours/complaints.
Every prespecified subgroup had sampled correct-sign probability 1.0. This
is a strong-signal regression fixture, **not the chapter panel or a calibrated
decision-probability demonstration**; it must not replace the harder benchmark.

Validation remains incomplete for repeated-panel interval coverage,
joint-decision calibration, null outcomes, deliberately different treatment
modifiers, weak/rare binary outcomes, outcome-order sensitivity and shared/private
prior sensitivity. The private component is a modeled escape from sharing,
not proof that these failure modes are resolved.

## Unchanged chapter comparison

The input is the exact 750-seller, 18-week panel used in the preceding repair:
two continuous log outcomes and one binary complaint outcome, with prespecified
benefit/harm regions and the same week-18/common-launch prediction scenario.
No outcomes, effects, seeds, group definitions or acceptance gates were changed.
Input SHA256:
`6e4904482c10530728705f0c4c1ce8bc943db6c2f0aa6aedfba5bb99284aaffd`.

Both joint models use four chains, 20 prognostic and **20 total treatment**
trees per outcome, adaptive coding, the existing GP/variance priors, and seed
314159. Sharing uses 10 shared + 10 private trees with fraction .5, declared
before the benchmark. The short budget is 2,000 burn-in, 250 saved draws and
thinning 2; the long budget is 8,000/500/4. Changing sharing changes the model,
topology prior and numerical path; this is not a claim of matched random draws
or matched wall-clock cost between the two model variants.

The freshly rerun default joint and separate-model short fits reproduced
**every saved effect draw bit for bit** from the preceding full-precision review
for all three outcomes. The previous longer default joint and separate-fit results are retained as
historical comparators, not described as fresh fits of this revision.

### Short budget (2,500 iterations per chain)

| Outcome/model | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| GMV, separate | .05056 | 1.250 / 13.0 | 1.149 / 19.1 |
| GMV, default joint | .05123 | 1.170 / 18.6 | 1.595 / 7.0 |
| GMV, shared/private | .02841 | 1.887 / 5.8 | 1.347 / 9.4 |
| Hours, separate | .02835 | 1.315 / 10.4 | 1.136 / 20.7 |
| Hours, default joint | .02306 | 1.087 / 35.5 | 1.216 / 13.9 |
| Hours, shared/private | .02192 | 1.044 / 77.3 | 1.155 / 17.9 |
| Complaint, separate | .01800 | 1.783 / 6.3 | 1.704 / 6.5 |
| Complaint, default joint | .01124 | 1.482 / 7.9 | 1.413 / 8.8 |
| Complaint, shared/private | .00735 | 1.047 / 61.9 | 1.143 / 19.9 |

Short shared complaint benefit/harm means are -.02230/.02390 against truths
-.02814/.02186. Their sampled 95% intervals exclude zero and contain the truths,
and sampled correct-sign probabilities are .997/.996. Those are encouraging
but **under-resolved MCMC estimates**, not certified probabilities: the
unchanged R-hat <= 1.01 and bulk/tail ESS >= 400 gates still fail. GMV's benefit
mixing is worse, demonstrating that the convergence issue is not binary-only.

Short fit wall times, including compilation and under concurrent workloads,
were 205.2 seconds for the default joint model and 358.6 seconds for sharing.
For complaint benefit/harm, bulk ESS per fit-second is .0387/.0431 versus
.1725/.0554. These are descriptive measurements from one nonconverged run,
not a controlled throughput or efficiency guarantee. Prediction time is excluded.

### Long budget and final disposition

The shared run completed at **10,000 iterations per chain**: 8,000 burn-in
and 500 retained draws thinned by four. The longer default joint and separate
comparators below are the completed, corrected runs from the preceding review;
their short-run equivalents were reproduced bit for bit with this revision.

| Outcome/model | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| GMV, separate (prior review) | .03490 | 1.058 / 49.5 | 1.056 / 49.5 |
| GMV, default joint (prior review) | .03459 | 1.032 / 120.7 | 1.332 / 9.7 |
| GMV, shared/private | .03560 | 1.407 / 8.5 | 1.305 / 10.3 |
| Hours, separate (prior review) | .02968 | 1.069 / 43.6 | 1.053 / 63.7 |
| Hours, default joint (prior review) | .02411 | 1.320 / 9.8 | 1.056 / 52.6 |
| Hours, shared/private | .02195 | 1.345 / 9.5 | 1.095 / 29.9 |
| Complaint, separate (prior review) | .01619 | 1.937 / 5.7 | 1.812 / 6.0 |
| Complaint, default joint (prior review) | .01143 | 1.530 / 7.3 | 1.484 / 7.9 |
| Complaint, shared/private | .00762 | 1.148 / 20.1 | 1.254 / 12.0 |

Complaint RMSE is approximately 33% lower than the default joint comparator
and 53% lower than the separate fit. That is point-recovery evidence from
this realization, **not a calibrated uncertainty claim**. Shared complaint
benefit/harm means are -.01651/.02231, versus truths -.02814/.02186. Their
sampled 95% intervals are [-.04009,.00130] and [.00003,.04869], with sampled
correct-sign probabilities .9595/.975. The benefit interval now includes zero;
the short-run sign probabilities and R-hats were more flattering. Bulk ESS
is only 20.1/12.0 and tail ESS 121.7/18.3. Those probabilities remain too
under-resolved for the chapter's decisions.

GMV's short-run RMSE improvement does not survive the longer comparison. Its
benefit-group mean is .30506 against truth .22503, and its interval
[.24721,.36864] misses that truth. GMV benefit/harm R-hats are 1.407/1.305;
hours R-hats are 1.345/1.095. Thus the remaining problem is **not binary-only**,
and a corrected binary augmentation cannot by itself certify the joint fit.

Long shared fitting took 1,062.1 seconds including compilation under concurrent
workloads. Complaint benefit/harm bulk ESS per fit-second is .0189/.0113,
lower than in the short run. This illustrates why a short-run efficiency
estimate from poorly mixing chains is not a stable performance claim.

**Disposition:** retain the implemented model as opt-in experimental
functionality; keep the default unchanged and the executable chapter example
withheld. All six shared subgroup summaries fail the R-hat <= 1.01 gate and
all have bulk ESS below 400. The separate/default joint comparators also fail
the original chapter gates. Do not promote a shorter run, lower RMSE, or the
high-signal regression fixture as a substitute for this result. Broad
negative-transfer and repeated-simulation calibration work remains outstanding.

An interruption stopped the first long-run and full-Python attempts before
they saved final results. Their partial execution is not counted as a pass.
Only the unfinished jobs were restarted, using unchanged source and settings.

## What the stochtree vignettes do and do not suggest

The [probit vignette](https://stochtree.ai/vignettes/probit.html) uses
Albert–Chib truncated-normal augmentation, fixed latent variance and an offset
for class imbalance. LongBet already has these ingredients, extended to
downstream-informed mixed-outcome conditionals. Its four-chain example does
not demonstrate convergence of LongBet's heterogeneous panel effects. There
is no missing elementary probit step to copy as an automatic R-hat repair.

The [ordinal vignette](https://stochtree.ai/vignettes/ordinal-outcome.html) uses
an asymmetric complementary-log-log model and stronger forest regularization.
That is a different likelihood, not a drop-in Gaussian-SUR probit augmentation.
Link and prior sensitivity could be worthwhile if scientifically justified,
but neither changing the link nor relabeling a binary observation as Gaussian
fixes the present sampler by itself. No stochtree benchmark was run and no
ordinal/cloglog feature is claimed here. The useful immediate lesson is to
keep identification explicit and check probability-scale quantities, rather
than infer causal uncertainty from a good latent-prediction plot.

## Reproduction and provenance

Working tree on engine HEAD `1227b28`; no commit, push or container rebuild was
performed. Python source fingerprint (filenames and bytes, as in the benchmark):
`3ea138472f0c38b3b387b985d7d730fdf01dda99c1007ee1d7910618121f96db`.
Python/R/contract fingerprint using the chapter plan's command:
`6430abd85f9e60a889ba1ce4f198803cca68d7d9c53119707d53d4f0b8b41d37`.
These identify source contents, not a released package version.

Audit artifacts are in `/tmp/longbet-shared-private.9CNbKH`, with per-fit
`report.json` and compact `draws.npz`, plus pytest XML. They are transient audit
files, not chapter runtime dependencies. Historical baselines are in
`/tmp/longbet-full-sur.n5aDoP/{scalar,joint}-indexed-{short,long}`.

```bash
taskset -c 0-3 .venv/bin/python benchmarks/bench_multi_chapter.py \
  --data /tmp/longbet-reference-benchmark.9G6Fop/reexported-input.npz \
  --output /tmp/longbet-shared-private.9CNbKH/shared-short \
  --mode multi --coding adaptive --shared-trees 10 \
  --burnin 2000 --draws 250 --skip 2
```

Choose a fresh output directory on rerun; the script refuses to overwrite a
completed result. For the long run use `--burnin 8000 --draws 500 --skip 4`; for
the default joint comparator omit `--shared-trees`. To regenerate the input,
use `export_multi_chapter.R` as described in the preceding repair report and
verify the byte hash before comparing results. A full chapter execution/render
with invalidated old caches remains separate work; no old scalar chapter
interval or economic valuation is certified by this multi-outcome audit.
