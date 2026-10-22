# Follow-up mixing repair and controlled DGP checks (2026-09-07)

**Outcome: no reliable repair was established.** All six short variants and
the selected 10,000-iteration-per-chain longer check failed the unchanged
chapter convergence gates. The new updates remain benchmark-only; public
defaults are unchanged and the worked joint example remains withheld. This
is evidence that the tested repairs are insufficient, **not proof that the
problem is mathematically impossible to solve**.

This review follows `shared_private_review.md`. Its purpose is to test repairs,
not to obtain a publishable-looking example by changing the simulated truth.
The chapter's original input, target, benefit/harm groups and acceptance gates
remain unchanged. Experiments are isolated behind benchmark-only flags; none
silently replaces the public Python/R sampler or its archive semantics.

## Same-posterior Gaussian blocks

`src/longbet/_joint_gaussian.py` implements an additional exact Gibbs transition,
conditional on current tree partitions/contrasts, coding coefficients, latent
responses, SUR loadings and variance parameters. It uses the full observed-data
precision, not just the outcome's own structural equation. It does **not**
pretend that GP coefficients and treatment-tree leaves are jointly Gaussian:
their product is nonlinear, and treatment leaves stay fixed within this block.

For each outcome, whiten the exposure GP by its kernel Cholesky factor and
the unit effects by their current prior SD. The prior precision in these
coordinates is identity. Accumulate the likelihood precision and score by
exposure and unit, integrate the block-diagonal unit effects with a Schur
complement, draw all GP coefficients jointly, and immediately redraw every
unit's outcome vector conditionally. That last redraw is essential: other
updates must not use old unit effects after a marginalized GP draw.

An optional larger block also includes one global location coordinate of each
prognostic forest. With J trees, L active leaves, physical leaf values theta,
and physical leaf precision p, define a = J sum(theta)/L. Subtracting a/J from
every active leaf holds the orthogonal contrasts fixed. Conditional on topology,
a has independent Gaussian prior variance J^2/(L p). Its likelihood loading is
the outcome's prognostic scale alpha. The larger block draws these coordinates
together with the GPs and unit effects and shifts the actual forest leaves,
cached fitted values and residuals consistently. This is **one global direction
per prognostic forest**, not a joint draw of every prognostic-tree coefficient.

The blocks are appended to an ordinary full sweep. This composition preserves
the same target; variances need not be redrawn again immediately to make the
composition valid. The calculations/factorizations run in float64, with stored
state in float32. Zero-design dummy coordinates preserve fixed-beta and disabled
random-intercept settings. Missing observations have zero likelihood precision.
Binary structural variances stay exactly one. No extra covariance jitter, prior
change, truncation approximation, or acceptance-threshold relaxation is used.

Correctness evidence: all **12** tests in `tests/test_joint_gaussian.py` passed.
They compare Schur means/covariances and 16,000 actual conditional draws with an
independently constructed dense observed-likelihood Gaussian posterior, with
SUR on/off, supported missingness, signed/zero weights, inactive coordinates,
prior-only units and the optional location coordinates. Integration checks cover
one/two chains, private/shared forests, fixed settings, raw residual identities,
actual shifted-forest predictions and the identifying binary variance. An
initial integration failure was a nonexistent prognostic-prior attribute name;
it was corrected before the expanded block's chapter benchmark was started.
These are conditional-kernel checks, not repeated-simulation calibration of
the complete nonlinear forest sampler.

### Additional unit-variance interweaving

Short expanded-block traces showed unit-variance R-hats 1.167 (complaint),
1.230 (GMV), and 1.074 (hours), while continuous innovation-variance R-hats were
1.011/1.001. This motivated an additional benchmark-only noncentered move in
`_unit_interweave.py`. It proposes gamma' = c gamma and v' = c^2 v, holding the
standardized unit effects fixed. With log(c) symmetric, the full Gaussian unit
prior, inverse-gamma variance prior and transformation/proposal Jacobian leave
the prior contribution `-2*a*log(c) + b/v*(1-c^-2)`. The likelihood contribution
uses the complete observed SUR quadratic, including cross-equation feedback.
Accept/reject updates the unit effects, their variance and raw residual together;
binary latent responses and their identifying variance are unchanged.

All five tests pass: comparison of the actual acceptance ratio with an
independently evaluated full joint density and Jacobian for every outcome,
accepted/rejected draws, chain-axis state/residual consistency, fixed settings,
and the noncentered-coordinate identity. This is a valid transition, not by
itself evidence that mixing improves. The short original-data fit composing it
with the expanded Gaussian block failed all six subgroup R-hat/bulk-ESS checks:

| Outcome | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| GMV proportional effect | .043016 | 1.5753 / 7.1 | 1.5684 / 7.1 |
| Hours proportional effect | .024269 | 1.4108 / 8.5 | 1.1043 / 28.0 |
| Complaint risk difference | .006183 | 1.1181 / 26.9 | 1.0256 / 144.3 |

Unit-variance R-hats became 1.179/1.668/1.015 (complaint/GMV/hours): this was
not a uniform improvement even for the targeted variance parameters. Complaint
RMSE was the best of these short variants, but that does not rescue inference.

## Distinct model restrictions

The benchmark also supports clearly labeled alternatives, without pretending
they preserve the default posterior:

- `--treatment-only` fixes b0=0 and b1=1, removing the untreated treatment-forest
  term. It differs from the public `adaptive_coding=False`, which fixes both
  coefficients to one. No data, covariates or simulated truth is changed.
- Adding `--fixed-beta` removes GP sampling and leaves beta equal to one; the
  treatment forest can still split on exposure and calendar time. This is not
  the full LongBet GP model and does not establish its extrapolation properties.
- `--no-calendar-trt` removes only calendar-time splits from treatment forests,
  retaining exposure splits and all prognostic calendar-time splits. The
  restriction is applied before forest/cache initialization, including shared
  trees. It matches the original simulated treatment mechanism, which depends
  on exposure rather than an additional calendar-time interaction. It is a
  model restriction requiring its own validation, not a general-purpose fix.

The preexisting public `split_time_trt=False` removes **exposure** splits, not
calendar-time splits. Do not confuse those two restrictions or infer a rank-one
treatment-effect requirement from the scalar GP multiplier.

## Controlled data checks

`make_mixing_controls.py` produces a paired simulation ladder with two continuous
outcomes, one genuinely binary outcome, and opposite-signed effects for observed
benefit/harm groups in every variant. The baseline has N=240, T=12, a common
randomized launch, simple exposure curves and independent residuals. Other
variants separately introduce rarity, residual correlation, unit intercepts,
or group-specific exposure shapes. They share covariates, assignment and the
underlying random numbers. All three generator/summary tests pass.

Completed fits use data seeds 20260907 and 20260908 and model seed 314159, four chains, 20
prognostic trees, 10 shared + 10 private treatment trees, fraction .5, depth 10,
and 2,000 warmup plus 1,000 subsequent iterations. First-seed baseline fits saved every
iteration; all remaining fits saved every second iteration. The paired rarity
comparison below also retains only iterations 2,4,...,1000 of each baseline
chain, preserving exact chain boundaries and matching the retained budget.

| Existing shared/private sampler | Complaint benefit R-hat / bulk ESS | Complaint harm R-hat / bulk ESS |
| --- | ---: | ---: |
| Seed 20260907, balanced-ish (44.2% positives) | 1.0054 / 673.3 | 1.0093 / 678.2 |
| Seed 20260907, rare (8.4% positives) | 1.0413 / 159.7 | 1.0071 / 546.4 |
| Seed 20260908, balanced-ish (42.4% positives) | 1.0030 / 481.5 | 1.0191 / 672.3 |
| Seed 20260908, rare (8.3% positives) | 1.0926 / 30.4 | 1.0539 / 60.8 |

The original GP/unit block on the balanced baseline, using all 1,000 saved
draws per chain, gives complaint R-hats 1.0067/1.0069 and bulk ESS 698.1/901.4.
The corresponding ordinary sampler gives 1.0073/1.0113 and 691.4/695.3. The block
is not a substantial universal improvement. Some continuous subgroup summaries
fail even on this easier fixture. Both paired seeds show worse benefit-group
binary mixing with rarity, but two pairs do not identify a sole cause or
quantify a general effect. The other generated ladder variants have not yet
been fitted here. This limited two-seed check is not a repeated-simulation
coverage study or simulation-based calibration.

## Original chapter acceptance test

Input: `/tmp/longbet-reference-benchmark.9G6Fop/reexported-input.npz`.
SHA256: `6e4904482c10530728705f0c4c1ce8bc943db6c2f0aa6aedfba5bb99284aaffd`.
Use the exact same 750 sellers, 18 weeks, common-launch scenario, week-18
evaluation, subgroup definitions, and response transformations as the preceding
review. In particular, transform **each** continuous draw with expm1 before
subgroup averaging; complaint effects are probability differences.

Short fits use four chains, 2,000 burn-in, 250 retained draws, thinning 2,
20 prognostic and 20 total treatment trees, depth 10 and model seed 314159.
The unchanged gates include R-hat <= 1.01 and bulk/tail ESS >= 400 for all six
subgroup effect summaries, plus recovery and uncertainty checks. Favorable
point RMSE or sign probabilities do not override these gates.

Completed short same-posterior GP/unit block:

| Outcome | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| GMV proportional effect | .039446 | 1.0586 / 49.7 | 1.5660 / 7.0 |
| Hours proportional effect | .023890 | 1.1951 / 14.4 | 1.1710 / 16.4 |
| Complaint risk difference | .006696 | 1.1801 / 15.7 | 1.2532 / 11.8 |

Completed short treatment-only model, with the GP retained:

| Outcome | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| GMV proportional effect | .033507 | 1.3447 / 9.8 | 1.1493 / 18.6 |
| Hours proportional effect | .025317 | 1.4602 / 8.0 | 1.0460 / 71.7 |
| Complaint risk difference | .009350 | 1.2378 / 14.1 | 1.1687 / 18.5 |

The remaining short comparisons also failed:

| Variant/outcome | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| GP/unit/prognostic-location block: GMV | .037879 | 1.0956 / 29.2 | 1.3075 / 10.5 |
| GP/unit/prognostic-location block: hours | .025172 | 1.1558 / 17.8 | 1.1128 / 26.4 |
| GP/unit/prognostic-location block: complaint | .008160 | 1.1140 / 26.8 | 1.1006 / 29.2 |
| Treatment-only, fixed beta: GMV | .040236 | 1.3512 / 9.5 | 1.4321 / 8.3 |
| Treatment-only, fixed beta: hours | .026181 | 1.0984 / 28.5 | 1.4994 / 7.6 |
| Treatment-only, fixed beta: complaint | .006759 | 1.1084 / 59.1 | 1.1933 / 15.3 |
| No calendar treatment splits: GMV | .036239 | 1.2034 / 15.0 | 1.3023 / 10.5 |
| No calendar treatment splits: hours | .022141 | 1.0761 / 36.3 | 1.2262 / 13.3 |
| No calendar treatment splits: complaint | .006847 | 1.1848 / 16.4 | 1.3721 / 9.6 |

All 30 subgroup summaries in these five short variants fail both R-hat and
bulk-ESS thresholds. Several complaint point estimates are encouraging, but
they do not validate the posterior decision probabilities. The expanded exact
block has the smallest worst-subgroup R-hat among these five variants and
preserves the original model; it was selected for an 8,000-burn-in, 500-saved,
thinning-4 longer check. That run uses the same seed, not an independent data
replicate. Including the subsequent variance-interweaving experiment, all 36
short subgroup summaries failed both thresholds.

### Longer expanded-block check: 10,000 iterations per chain

| Outcome | RMSE | Benefit R-hat / bulk ESS | Harm R-hat / bulk ESS |
| --- | ---: | ---: | ---: |
| GMV proportional effect | .032190 | 1.4081 / 8.5 | 1.0893 / 31.8 |
| Hours proportional effect | .024370 | 1.2219 / 13.3 | 1.1865 / 15.1 |
| Complaint risk difference | .006664 | 1.0870 / 36.5 | 1.1974 / 14.6 |

All six again failed R-hat and bulk-ESS thresholds. Complaint tail ESS was
534.5/89.0, so one adequate tail diagnostic did not establish convergence.
Its benefit/harm mean effects were -.01831/.02529, against truths -.02814/.02186;
the benefit interval [-.04133,.00035] included zero. Sampled sign probabilities
.9725/.9895 remain under-resolved. The GMV benefit interval [.22991,.36587]
missed truth .22503 and the corresponding chain means ranged .26869 to .32822.
Fit time was 1,242.8 seconds, including compilation under concurrent workloads;
it is not a controlled speed comparison. A longer run did not validate the
short-run impression of better mixing.

### Decision and remaining limitations

Do not enable these experiments as new public defaults, publish the worked
joint-decision example, or claim that an extra probit/variance update fixes the
problem. The controlled data show that binary outcomes are not inherently
unusable and that rarity can aggravate mixing. The original data also contain
complex exposure shapes, sparse target-exposure support in the benefit group,
and unobserved seller-level baseline variation. Gaussian/probit innovation
terms are compatible with the likelihood, but that does not make the entire
DGP an exact prior-generated draw of the fitted model; unmodeled seller drift
can matter. Neither a difficult DGP nor correct conditional kernels excuses
failed posterior diagnostics.

More substantial work on whole-tree/global mean moves or the model's exposure
and within-unit structure may still succeed. That is an untested direction,
not a promised fix or a reason to keep adding iterations indefinitely. It
would need independent correctness tests, properly specified priors for SBC,
repeated-seed recovery/coverage, and the same untouched chapter acceptance
test. Separate fits remain a starting comparison, but the original scalar
binary fit also failed diagnostics; separate fitting is not a blanket guarantee.

Final verification: 12 joint-block tests, 5 interweaving tests and 4 existing
sampler/benchmark regression checks passed (21 targeted tests). The earlier
full Python/R-suite results belong to the previous review and were not rerun
in full here. All 42 chapter R chunks parse and setup/cache guards execute
against source-mounted code. No full chapter render, persistent image rebuild,
commit or push was performed.

## Reproduction and audit record

All current results, compact effect draws, logs and test XML are under
`/tmp/longbet-mixing-repair.SQ7Mgz`. These are temporary audit artifacts, not
chapter dependencies. Each completed fit's JSON records its input hash, Python
source fingerprint, seed, config, experimental flags, chain means, diagnostic
summaries and fit time. Later runs also save compact parameter traces separately.
Those parameter arrays use **internal order: complaint, GMV, hours**; effect
archives and report outcome dictionaries use outcome names. Do not interpret
the parameter arrays in the caller's GMV/hours/complaint order.
Different fingerprints may include changes to the unused experimental module;
the public sampler path was not switched by these experiments.

Final Python/R/contract source fingerprint:
`0a20d5420c3826a3fb1b15b480aed8cb712fe19e8c8cb7d9f1bdb5e980a48069`.
The longer block-only run's Python-only fingerprint is
`ad6372f0779f4222c72a57069732ed4a111e6be7bdf5d8b0c72e27a4c1047678`;
that run preceded addition of the independently opt-in variance module.

The first ordinary baseline-control attempt was killed for exceeding available
memory while two other large fits ran concurrently. It produced no completed
report/draw archive. It was rerun successfully with less concurrency; the
original failed log and the retry log are both retained. This is a resource
failure, not a posterior diagnostic. More thinning was used only to limit
storage in the paired rarity check, with a matched retained baseline as above.

Example commands (choose fresh output directories):

```sh
# Re-export the unchanged input if the temporary audit archive is gone.
# Run in the book's R/JAX environment; this executes simulation chunks only.
Rscript benchmarks/export_multi_chapter.R /home/ignacio/book /tmp/chapter-input.npz /home/ignacio/longbet-jax
.venv/bin/python benchmarks/make_mixing_controls.py --output /tmp/mixing-controls
.venv/bin/python benchmarks/bench_multi_chapter.py \
  --data /tmp/chapter-input.npz \
  --output /tmp/joint-location-audit --shared-trees 10 \
  --joint-gaussian --joint-location --burnin 2000 --draws 250 --skip 2
.venv/bin/python -m pytest tests/test_joint_gaussian.py tests/test_unit_interweave.py tests/test_multi_step.py tests/test_multi_benchmark.py -q
```

Use `--burnin 8000 --draws 500 --skip 4` for the longer check and append
`--interweave-unit-variance` only for the separately labeled variance experiment.
The latter has only the short-budget empirical check in this review. Re-export
requires the preserved `multi-draft.md` generator chunks; its rejected narrative
and numerical claims must not be reused. The original archive hash records the
tested file; a re-export's ZIP metadata can change the file hash even when all
13 arrays are identical, so compare array contents as well as file provenance.
