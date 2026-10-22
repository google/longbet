# Deeper multi-outcome repair (2026-09-07)

## What changed

The repaired working tree on `1227b28` uses `full_precision_sur_v1`, with
`precision_cache_version=1` for scalar and multi archives. The package version
has not changed. This repair fixes identified transition-kernel errors; it
does **not** claim that joint fitting always beats separate fits or that the
chapter's decision probabilities are ready to use.

The old reference-compatible update passed each downstream equation an
incoming residual offset but omitted the downstream likelihood when updating
earlier means and binary latent responses. The replacement conditions each
mean block on the full contemporaneous precision of the specified triangular
model. It retains separate outcome forests and priors.

### Full likelihood and supported missingness

Let `r = y_star-f`, `B=I-Gamma`, `e=B r`, and structural innovation variances
`v`. With complete observations the precision is
`Omega = B' diag(1/v) B`. For outcome `m`:

```
precision_m = Omega[m,m]
conditional_residual_m = (Omega @ r)[m] / precision_m
offset_m = r[m] - conditional_residual_m
```

The precision, not just the offset, enters both forest updates, the GP's
sufficient statistics, alpha, adaptive coding and unit-intercept updates.
Temporary chain-specific weights do not overwrite shared observation masks.
Raw residuals are restored at the state boundary.

Under the supported masks, every observed continuous equation has all its
predecessors observed. Unobserved downstream innovations integrate out. Use
cell-specific weights:

```
p_m = sum_k observed[k] * B[k,m]^2 / v[k]
q_m = sum_k observed[k] * B[k,m] * e[k] / v[k]
conditional_residual_m = q_m / p_m
```

Unobserved response cells get benign residual zero and precision one, and are
excluded from the likelihood. Binary-only masks need not be nested because
binary rows have zero incoming loadings. These formulas are **not** permission
to use arbitrary missing predecessors; the input validator still rejects them.

Binary latent responses follow their label-truncated conditional Gaussian,
with mean `current_latent-conditional_residual` and variance `1/p_m`.
That variance can be below one although the marginal latent variance remains
identified at one. Loading draws use current raw residuals. Continuous
variance draws use their structural innovations `e_m`, never the temporary
conditional residual. The binary identifying variances remain fixed.

### A second bug: dynamic precision versus bartz's cached leaf sums

Analytical tests of the actual forest kernel initially passed posterior means
but failed posterior variances and cross-outcome covariance. For example, a
single weighted constant leaf produced variance about .0844 instead of the
analytical .0976 in the initial prototype.

bartz 0.12.1 assumes observation weights remain fixed after initialization.
Its per-step precision reduction updates only nodes touched by a proposed
grow/prune; other leaves reuse `Forest.prec_tree`. LongBet changes
`w=b_Z*beta_S` between sweeps but had only rebuilt the observation-scale
arrays. Untouched leaves therefore used stale precisions. **This affects
the scalar sampler too**, and also affected the initial full-SUR prototype.

`_forest_cache.refresh_prec_tree` now recomputes all occupied-leaf precision
sums before every treatment-forest call and before dynamically weighted
prognostic-forest calls. It groups observations using pending-prune-adjusted
indices without changing the deferred-prune bookkeeping. Gaussian variance
and covariance oracles pass with this correction.

The first implementation of this refresh reused bartz's default precision
reduction, which builds dense membership arrays for just two affected leaves.
Refreshing *all* 1,024 default leaf slots that way caused excessive memory and
work: two large joint attempts terminated with exit 137. The final helper uses
an indexed reduction with at most 16 batches, retaining the configured
data-axis collective. Storage is proportional to observations plus batched
leaf slots, not their dense product. A regression test inspects intermediate
JAX array shapes and verifies sums; it would reject the dense implementation
even on a machine with enough RAM to run it. The probabilistic target is
unchanged by this reduction choice.

An additional small fix prevents the beta–nu rescaling move from changing beta
when `sample_beta=False`.

## Persistence and compatibility

Public array layouts and RNG splitting remain unchanged. Coupling-off and
zero-prior multi paths retain matched-key scalar equivalence. Full-SUR draws
are deliberately different from the old coupled sampler.

Old scalar and multi archives lack `precision_cache_version=1` and are
rejected with an instruction to refit from the original data. Multi archives
also require the new semantics tag. This is conservative for some fixed-weight
legacy configurations, but avoids certifying old draws with unknown sampler
provenance. Loading or changing tree cutpoints cannot upgrade posterior draws.
The earlier time-grid validation remains in force, including extracted scalar
children. Chapter cache keys must include actual engine source, not package
versions alone.

## Validation protocol

`tests/test_full_sur.py` exercises the actual weighted tree update:

- analytical Gaussian leaf mean and variance with unequal cell precisions;
- joint three-outcome Gaussian posterior mean and full covariance, with and
  without nested missingness;
- marginal observed-data precision versus independently inverted covariance
  submatrices, including non-nested binary-only masks;
- binary conditional truncation and moments;
- scalar treatment draws with changing weights, including zero weights;
- precision caches with pending prunes and untouched leaves;
- bounded-memory full-leaf reductions without dense membership arrays;
- fixed-beta behavior with the ridge option enabled;
- exact scheduled-key checks of precision-weighted alpha, GP, coding and
  unit-intercept sufficient statistics;
- a probit/Gaussian mean posterior and joint sign event compared with direct
  numerical integration of the **observed-data likelihood**, not a latent
  Gaussian stand-in.

Other regression tests cover current raw residuals against independently
evaluated forests, structural variance draws using exact scheduled keys,
off-coupling scalar equivalence, multi-chain mapping, malformed archive
rejection, and R/Python draw-order and serialization contracts.

These are correctness oracles and regression tests, not a large repeated-panel
coverage study of mixed-outcome treatment effects. That remains a separate
requirement before claims of calibrated rollout probabilities.

The final-engine full Python run returned **196 passed, 1 skipped, 1 failed**.
All 13 full-SUR/precision-cache oracle tests passed, as did the existing mixed
benefit/harm regression without weakened thresholds. The sole failure was the
existing Geweke mutation-sensitivity test: at its original 3,000-sweep budget,
the correct and deliberately corrupted samplers gave whitened second moments
.82246 and 1.44996. Its absolute-error checks passed, but the ratio of the two
errors missed the required factor of four. Increasing **both** runs to 12,000
sweeps with 2,000 warm-up, retaining seed 11, the same 50% mutation and all
original assertions, passed the targeted rerun. No engine change followed the
full run. Thus 197 distinct Python tests have passed across the full run and
that corrected-budget rerun, with one optional external-reference skip; a new
single full-suite run after this test-budget change was not performed.
The records are `pytest-indexed-final.xml` and `pytest-geweke-sensitivity.xml`
in the audit directory below. The failed short-budget record is retained.

The final indexed-cache R testthat run passed 245 assertions, with no failures,
warnings or skips, including the new legacy-prediction guard and a fresh-process
three-outcome serialization check. All chapter R chunks parse; the revised
setup, sampler guards and source-based cache fingerprint execute successfully
in a disposable `book` container with the current R/Python source mounted.
The persistent book image was not rebuilt and the full chapter was not rendered.

## Unchanged chapter-data comparison

Use `bench_multi_chapter.py` and the input with SHA-256
`6e4904482c10530728705f0c4c1ce8bc943db6c2f0aa6aedfba5bb99284aaffd`.
The exporter and full prior protocol are in
[the historical review](multi_outcome_review.md). Both corrected separate
and corrected joint models must be refitted; comparing against old scalar
draws would retain the precision-cache bug in the baseline.

The short protocol is four chains, 2,000 burn-in, 250 retained draws, thinning
2, 20 trees per forest, adaptive coding, and `lambda_knl=2`. Data, true
effects, groups, treatment scenario and seeds are unchanged. Group truth is
`.22503/-.17458` for GMV, `-.23114/.21804` for hours and
`-.02814/.02186` for complaint risk differences (benefit/harm groups).

Fresh audit directory: `/tmp/longbet-full-sur.n5aDoP`.
Re-exporting the final chapter's simulation with the usual native CPU setup
produced a byte-identical input archive (`reexported-after-full-repair-native.npz`).
A restricted-CPU export has only floating-point roundoff in `x` (maximum
`1.04e-16`, identical at float32); it is retained but is not used for fitting.
`chapter-short` contains an initial full-precision prototype **before the
cache fix** and is invalid as a posterior benchmark. It is retained only as a
failed development stage. Two subsequent joint attempts using dense full-cache
refresh terminated with exit 137 before writing results.
`scalar-cache-fixed-short` contains a complete but slow dense-refresh baseline.
The final same-source comparison uses `scalar-indexed-short` and
`joint-indexed-short`, followed by `scalar-indexed-long` and `joint-indexed-long`.
Timings under concurrent workloads or changed CPU affinity are not controlled
performance comparisons.

The final engine source fingerprint (Python/R/contract files, computed using
the command in the chapter plan) is
`8be30952cf36c9641bc4680d20bd9e98c554b4ef890c3db78ed971ab32bf6ee3`.
The benchmark's Python-only fingerprint, which hashes filenames and file bytes,
is `29bb0cbec29a12c0c7a230b074d6695e9971fda0ad32c6b2a41c4db489ad30f0`.

### Final indexed-cache short runs

All effect metrics below are on the business scale. R-hat and ESS refer to
the two prespecified **draw-wise group means**, not averaged individual diagnostics.

| Model | Outcome | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | --- | ---: | ---: | ---: |
| Separate | GMV | .05056 | 1.250 / 13.0 | 1.149 / 19.1 |
| Joint | GMV | .05123 | 1.170 / 18.6 | 1.595 / 7.0 |
| Separate | Hours | .02835 | 1.315 / 10.4 | 1.136 / 20.7 |
| Joint | Hours | .02306 | 1.087 / 35.5 | 1.216 / 13.9 |
| Separate | Complaint | .01800 | 1.783 / 6.3 | 1.704 / 6.5 |
| Joint | Complaint | .01124 | 1.482 / 7.9 | 1.413 / 8.8 |

The joint complaint group means are -.01138 and .01727, compared with truth
-.02814 and .02186. Their intervals are [-.03978,.00921] and
[-.00059,.04153]; empirical correct-sign probabilities are .755 and .951.
The separate complaint means are .00471 and .00812, with correct-sign
probabilities .430 and .799. **These probabilities are under-resolved by the
chains and must not be used for rollout decisions.** Lower RMSE does not
override failed convergence gates.

The longer check retains the same model settings and seeds, using 8,000 burn-in,
500 saved draws per chain and thinning 4 (10,000 iterations per chain). Both
separate and joint models use this budget; the old long runs are not reused.

### Final indexed-cache longer runs

Both jobs completed, without the dense-cache allocation failure. Their data
and Python source hashes match the short runs and the fingerprints above.

| Model | Outcome | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | --- | ---: | ---: | ---: |
| Separate | GMV | .03490 | 1.058 / 49.5 | 1.056 / 49.5 |
| Joint | GMV | .03459 | 1.032 / 120.7 | 1.332 / 9.7 |
| Separate | Hours | .02968 | 1.069 / 43.6 | 1.053 / 63.7 |
| Joint | Hours | .02411 | 1.320 / 9.8 | 1.056 / 52.6 |
| Separate | Complaint | .01619 | 1.937 / 5.7 | 1.812 / 6.0 |
| Joint | Complaint | .01143 | 1.530 / 7.3 | 1.484 / 7.9 |

For the binary outcome, in **percentage points**:

| Model | Region | Truth | Mean | 95% interval | Sampled probability of correct sign |
| --- | --- | ---: | ---: | --- | ---: |
| Separate | Benefit | -2.814 | .047 | [-3.158, 3.324] | .550 |
| Joint | Benefit | -2.814 | -.766 | [-3.858, 1.801] | .677 |
| Separate | Harm | 2.186 | 1.352 | [-.836, 3.485] | .807 |
| Joint | Harm | 2.186 | 1.042 | [-1.027, 4.267] | .779 |

These sampled probabilities describe the failed diagnostic run; they are not
certified posterior probabilities. The joint binary RMSE is about 29% lower
than the fresh separate-fit baseline in this realization, but neither recovers
both subgroup signs with the required uncertainty. Even the continuous group
summaries fail the chapter's R-hat <= 1.01 / bulk and tail ESS >= 400 gates.
The joint GMV benefit interval [.25255,.35419] also misses its .22503 truth in
this panel. One simulated interval is not a coverage study, but it is another
reason not to equate small overall RMSE with reliable subgroup uncertainty.
The short joint binary harm-sign probability .951 falls to .779 in the longer
run; selecting the short run would conceal this instability.

**Decision: keep the worked joint chapter example withheld.** The repair fixes
demonstrated conditional-update errors, but does not solve the chapter's
practical convergence and mixed-outcome uncertainty problem at either tested
budget. Separate fits are not a converged fallback at these settings either.
This does not prove that every prior, proposal or longer budget must fail.
The main scalar chapter results also need fresh execution and diagnostic
review because their sampler was affected; the benchmark is not a full render
of that 3,000-seller, 60-tree chapter analysis.

Reproduction, from the engine checkout with the environment installed:

```bash
.venv/bin/python benchmarks/bench_multi_chapter.py \
  --data /tmp/longbet-reference-benchmark.9G6Fop/reexported-input.npz \
  --output /tmp/longbet-full-sur-rerun-separate \
  --mode scalar --coding adaptive --burnin 8000 --draws 500 --skip 4
.venv/bin/python benchmarks/bench_multi_chapter.py \
  --data /tmp/longbet-reference-benchmark.9G6Fop/reexported-input.npz \
  --output /tmp/longbet-full-sur-rerun-joint \
  --mode multi --coding adaptive --burnin 8000 --draws 500 --skip 4
```

Use fresh output paths: the runner refuses overwrites. If temporary artifacts
are no longer present, use `export_multi_chapter.R` with the original book
simulation/archived example and verify the input hash before comparison.
Each completed output directory contains `report.json` and `draws.npz`;
the latter preserves chain-major effect draws, not a complete fitted model.

## What mvbcf adds to the design discussion

The [mvbcf source/design review](mvbcf_design_review.md) identifies a promising
next extension: shared treatment partitions with vector-valued leaves, and
regularized outcome-specific departures to protect against negative transfer.
This could pool evidence about effect modifiers rather than only residual
dependence. It is a model change, not part of this repair. mvbcf has not been
empirically validated here, and its inspected Gaussian-response kernel is not
a drop-in mixed-outcome longitudinal model.
