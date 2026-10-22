# Localizing the multivariate-outcome mixing failure (2026-09-07 follow-up)

This review follows `mixing_repair_review.md`. It does two things:

1. It **localizes** the failure quantitatively. The previous reviews recorded
   that chains disagree; this one establishes *what kind* of disagreement it is
   and *which direction* of the posterior is responsible.
2. It implements and tests two exact, same-posterior Gibbs blocks aimed at the
   parts of that direction which can be characterized in closed form, and
   reports the measured effect on the unchanged chapter acceptance test.
3. It rules out two alternative explanations by experiment: that the DGP is
   unfittable (section 4.4) and that the chapter's estimand is the binding
   constraint (section 4.8).

**Status: the chapter convergence gates are still not met.** Both blocks are
benchmark-only; no public default changed, and the worked joint-decision example
stays withheld. What is new here is the diagnosis, two exact transitions
validated against dense oracles, a matched baseline replicated over three model
seeds, and the two ruled-out explanations.

## 1. Workspace verification

Verified before any fit, on the dirty worktree described by `bug.md`:

* Python/R/contract source fingerprint at handoff:
  `0a20d5420c3826a3fb1b15b480aed8cb712fe19e8c8cb7d9f1bdb5e980a48069` -- matches.
* Frozen acceptance input `/tmp/longbet-reference-benchmark.9G6Fop/reexported-input.npz`
  still exists, SHA256 `6e4904482c10530728705f0c4c1ce8bc943db6c2f0aa6aedfba5bb99284aaffd`
  -- matches. All 13 arrays present with the documented shapes and dtypes.
* Independently recomputed from that archive: `col = 11` (zero-based, week 18),
  complaint prevalence 6.2815%, benefit group n = 59, harm group n = 138, the
  evaluation rollout `ze` launching every seller in week 11 so that exposure at
  the target column is exactly S = 8, and all six business-scale subgroup truths
  to ten decimals (GMV +0.2250306663 / -0.1745837863, hours -0.2311425486 /
  +0.2180392275, complaint -0.0281442312 / +0.0218628292). Only five
  benefit-group sellers are observed at exposure 8 in week 18 factually.
* Python 3.12.3, JAX 0.11.1, bartz 0.12.1, numpy 2.5.2, CPU device only.
* Host: 16 logical CPUs, 28 GiB RAM, ~15 GiB available, swap fully committed.

**Practical note on the earlier exit-137.** `/tmp` on this host is a `tmpfs`
sized at 15 GiB, so every byte written to the previous audit directories was
resident memory. Large benchmark outputs written there compete with the fits
themselves. This audit therefore writes to disk-backed
`/home/ignacio/longbet-bug-audit`. This is an operational observation, not a
posterior diagnostic.

## 2. What kind of failure this is

### 2.1 The chains are internally healthy and mutually offset

From the preserved 10,000-iteration-per-chain run
`/tmp/longbet-mixing-repair.SQ7Mgz/chapter-location-long` (four chains, 500
retained each), GMV benefit-group draw-wise subgroup means:

| Chain | mean | within-chain SD | first half | second half | lag-1 autocorrelation |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | .32822 | .02272 | .32645 | .32999 | .092 |
| 2 | .28062 | .02987 | .29116 | .27008 | .369 |
| 3 | .31434 | .02765 | .30900 | .31969 | .172 |
| 4 | .26869 | .02769 | .27676 | .26062 | .338 |

Each chain is stationary-looking across its own halves and has *low* lag-1
autocorrelation, yet the chain means are separated by two to three within-chain
standard deviations. This is not slow local exploration. It is a set of chains
that each mix well in the directions they can see, sitting at different places
along a direction none of them traverses.

### 2.2 The slow direction is orthogonal to everything the sampler updates well

Same run, parameter traces (internal order complaint, GMV, hours):

| Quantity | R-hat | bulk ESS |
| --- | ---: | ---: |
| GMV innovation variance `sigma2` | 1.0077 | 724.9 |
| hours innovation variance `sigma2` | 1.0007 | 1603.2 |
| SUR loadings `Gamma[1,0]`, `Gamma[2,0]`, `Gamma[2,1]` | 1.0059-1.0068 | -- |
| GMV `b0` / `b1` | 2.2971 / 2.4039 | 5.1 / 5.0 |
| hours `b0` / `b1` | 1.9841 / 1.6428 | 5.6 / 6.9 |
| GMV exposure GP `beta`, per knot | 1.157-1.850 | -- |

The variance and loading blocks -- exact conditional draws with plenty of data
-- are converged. The coding coefficients and the exposure GP are not. Note
that `b`, `beta` and `nu` carry exact sign and scale symmetries that leave the
reported effect invariant, so their R-hat is not by itself a defect; the point
is the *contrast* with `sigma2`, which shows the failure is confined to the
multiplicative mean block rather than being a global sampler problem.

### 2.3 The autocorrelation gap: the slow mode is longer than the run

Using the fresh instrumented baseline of section 3, for the GMV benefit group,
per chain (250 retained draws each, thinning 2):

| Quantity | integrated autocorrelation time per chain | implied effective draws per chain |
| --- | --- | --- |
| effect `tau` | 3.5, 10.4, 7.3, 3.0 | 72, 24, 34, 82 |
| untreated surface `mu0` | 5.1, 31.3, 5.5, 24.8 | 49, 8, 46, 10 |
| contested split term `U` | 9.6, 14.6, --, -- | 26, 17, --, -- |

Within-chain autocorrelation predicts on the order of 200 effective draws in
total for `tau`. The benchmark's four-chain bulk ESS for the same quantity is
**5.8**. A roughly 35-fold gap between the within-chain and the between-chain
verdict is only possible if there is a mode whose autocorrelation time exceeds
the entire 2,500-iteration run. That is the defining fact about this failure:
**it is not slow mixing that more iterations fix at the observed rate; it is a
direction the current kernel barely moves at all.**

This also explains, without any further assumption, why the previously tested
repairs did not work. The joint GP/unit-intercept block, the global
prognostic-location coordinate, and the unit-variance interweaving move all
improve *conditional* Gaussian updates -- directions the chains already
explored adequately.

### 2.4 Which direction is slow

Every untreated cell has exposure `S = 0`. Writing the fitted mean as

```text
f_m = alpha_m*mu_m(x,t) + b_z,m * beta_m(S) * nu_m(x,S,t) + gamma_m(unit)
```

with `mu` a forest over `(x, t)` and `nu` a forest over `(x, t, s)`, the
treatment term on an untreated cell is `b0*beta(0)*nu(x,0,t)` -- a function of
`(x, t)`, the *same argument set* as the prognostic forest. Prediction makes
the consequence explicit (`_model.py`):

```text
mu0 = alpha*mu(x,t) + b0*beta(0)*nu(x,0,t) + gamma        (untreated surface)
tau = b1*beta(S)*nu(x,S,t) - b0*beta(0)*nu(x,0,t)         (effect)
```

Call the shared term `U = b0*beta(0)*nu(x,0,t)`. The data constrain the *sum*
`alpha*mu + U` on untreated cells; they constrain the *split* only through the
priors and through the treated cells. But `U` enters the reported effect with a
minus sign, and nothing cancels it.

This differs from cross-sectional BCF, where the analogous shift
`(b0, b1, mu) -> (b0 + d, b1 + d, mu - d*tau_forest)` leaves *both* potential
outcome surfaces and hence the effect unchanged. Here the treatment term is
evaluated at a *different exposure* for treated and untreated cells, so no
compensating shift of `b1` restores the treated surface. The direction is
therefore not a harmless symmetry: it is a near-flat but effect-relevant ridge,
pinned only by whatever treated cells share the relevant `(x, S, t)` region.

That is also why the benefit group is the worst subgroup: only five
benefit-group sellers are factually observed at the target exposure and week, so
the ridge is least constrained exactly where the chapter reads its answer.

### 2.5 Direct measurement of the ridge

`decompose_fit.py` (in the audit directory) recovers `U` from the public API
alone: a second prediction under an all-zero rollout returns
`tau_zero = (b1-b0)*beta(0)*nu(x,0,t)`, so `U = b0/(b1-b0) * tau_zero`. Four
chains, 2,000 burn-in, 250 retained, thinning 2, seed 314159, shared-trees 10.
Ratio column is between-chain SD over mean within-chain SD; near 1 means about
one effective draw per chain.

| Outcome / group | `tau` R-hat (b/w) | `mu0` R-hat (b/w) | `U` b/w |
| --- | --- | --- | ---: |
| GMV benefit | 1.8874 (1.85) | 1.6432 (1.38) | 3.08 |
| GMV harm | 1.3449 (0.94) | 1.3561 (1.03) | 4.97 |
| hours benefit | 1.0434 (0.31) | 1.2049 (0.71) | 3.90 |
| hours harm | 1.1581 (0.62) | 1.0704 (0.34) | 3.03 |
| complaint benefit | 1.0468 (0.11) | 1.0383 (0.25) | 1.22 |
| complaint harm | 1.1434 (0.43) | 1.2523 (0.68) | 1.66 |

Two facts follow. First, `U` has the worst between/within ratio of any
quantity measured, in every outcome and group: the chains genuinely disagree
about the split, as predicted. Second, the *consequences* differ by outcome. For
hours and complaint the effect largely survives the disagreement (`tau` ratios
0.11-0.62). For GMV it does not, and the damage reaches even `mu0` restricted to
never-treated sellers (R-hat 1.3354), whose untreated surface at week 18 is
directly observed rather than counterfactual. So for GMV the chains have not
equilibrated even in an identified direction, because that direction is strongly
correlated with the ridge.

These continuous-outcome diagnostics are computed on the model's own response
scale, not the business `expm1` scale used by the acceptance gates, so they are
internally comparable but not interchangeable with the tables in section 4.

### 2.6 What this does and does not establish

Established:

* The failure is a slow mode with autocorrelation time exceeding the run, not
  uniformly slow mixing (section 2.3).
* Variance and SUR-loading blocks are converged; the failure is confined to the
  multiplicative mean block (section 2.2).
* `alpha*mu` and `b0*beta(0)*nu(.,0,.)` are exactly confounded on untreated
  cells by construction, the effect subtracts one side of that split, and the
  split is empirically the least-converged quantity in the fit (2.4, 2.5).

Not established:

* That this ridge is the *only* slow direction. The archived `--treatment-only`
  variant fixes `b0 = 0`, removing `U` entirely, and still failed at GMV benefit
  R-hat 1.3447. Its improvement over the default is consistent with the
  diagnosis, but the residual failure shows at least one further slow direction
  -- most plausibly the competition between `mu(x,t)` and `nu(x,S,t)` on treated
  cells, together with tree-partition exploration for a target that requires
  extrapolating to `(S=8, week 18)` for sellers without factual support there.
* Any claim about misspecification. Poor recovery and poor mixing are separate
  failures and this section is about mixing only.

## 3. The repair: an exact joint prognostic-leaf / coding block

### 3.1 Mathematics

Conditional on the tree topologies, the treatment forest, the exposure GP, the
coding-free treatment factor `g = beta(S)*nu_fit`, the latent responses, the SUR
loadings and the variances, the prognostic leaf values and `(b0, b1)` enter the
mean **linearly**. Their joint conditional is therefore exactly Gaussian, and
drawing them together is an ordinary Gibbs block: **the target is unchanged and
this is the same model**, not a new one.

With `theta` the leaf values of one prognostic tree in data units, per cell
precision `Omega[m,m]` and score `(Omega r)[m]` taken from the raw full-model
residual (identical to `_sur.conditional_residual`, whose `precision` and
`score` are exactly these), the block precision and score are

```text
P = Lambda + D' diag(prec) D,     h = D' score - Lambda phi_cur
D = [ alpha * 1{cell in leaf l} | g*1{z=0} | g*1{z=1} ]
Lambda = diag(leaf_prior_cov_inv, ..., 1/sigma_b^2, 1/sigma_b^2)
```

and the draw is `phi_cur + N(P^-1 h, P^-1)`. Within one tree the leaves
partition the cells, so the leaf block of `P` is **diagonal** and `P` is an
arrowhead matrix. Sampling uses the exact conditional decomposition

```text
b     ~ N(S^-1 (h_b - A' D^-1 h_theta), S^-1),   S = B - A' D^-1 A
theta ~ N(D^-1 (h_theta - A b), D^-1)   given b
```

which costs `O(L)` per tree and factorizes only a 2x2 matrix. The block is
applied to each prognostic tree in turn, so `(b0, b1)` are re-drawn jointly with
a different tree's full flexibility 20 times per sweep, instead of once with all
20 trees held fixed. No approximation, jitter, prior change, truncation or
relaxed acceptance threshold is involved.

### 3.2 Invariants preserved

* Only *actual* leaves enter the block; unreachable slots are left untouched,
  matching `_ridge.py`. Cells always land on actual leaves.
* The forest `offset` is excluded by construction, because the block is
  parameterized by the change from the current value.
* `mu_fit` and the raw residual are updated together, so `resid = y - f` holds
  at the block boundary; `traverse_forest` is used for leaf membership, which is
  what `evaluate_forest` (and hence the residual identity test) uses.
* Full SUR precision and score, so earlier outcomes keep downstream feedback.
* Binary structural variances stay exactly one; nothing here rescales them.
* Unobserved cells carry zero precision and contribute nothing.
* `adaptive_coding=False` makes the block an exact no-op.
* Factorization in float64, stored state float32, chain axes preserved.
* `nu_fit` (which under sharing includes the shared ensemble) is read as the
  total treatment factor and is never modified. bartz's leaf-precision caches
  are refreshed each sweep by the existing `refresh_prec_tree` calls, so the
  changed coding weights are picked up.

### 3.3 Correctness evidence

`tests/test_joint_prognostic.py`, 7 tests, all passing:

* The actual draws of `_tree_block_draw` are compared with an **independently
  constructed dense Gaussian**: explicit design matrix, explicit `np.linalg.inv`,
  no reuse of the Schur derivation under test. 40,000 draws; every posterior
  mean matches within 4.5 Monte Carlo standard errors and the **full covariance
  matrix** matches to within 0.05 in correlation-scaled units. Run with and
  without any likelihood information; in the degenerate case the coding
  coefficients provably return to their exact prior.
* A separate test confirms the block is genuinely non-diagonal -- leaf/`b0`
  posterior correlation above 0.3 -- i.e. it moves the direction a separated
  scan cannot.
* Integration tests on a real multi-outcome state: the raw residual identity
  and the `mu_fit` identity both hold after the block, binary `sigma2` stays
  exactly 1.0, unreachable leaf slots are bit-identical, `adaptive_coding=False`
  is an exact no-op, and the chain-axis path gives distinct per-chain draws with
  the invariants holding within every chain.

These are conditional-kernel checks. They establish that the transition is the
claimed exact Gibbs block; they are **not** evidence that mixing improves. That
question is settled only by section 4.

Regression floor, run together with the new file:
`tests/test_full_sur.py tests/test_shared_forest.py tests/test_joint_gaussian.py
tests/test_unit_interweave.py tests/test_multi_step.py
tests/test_multi_benchmark.py tests/test_joint_prognostic.py` -- **59 passed**.

## 3b. The second ridge: a paired prognostic/treatment leaf block

Section 2.6 records that removing the first ridge cannot be the whole story:
the archived `--treatment-only` variant fixes `b0 = 0`, which deletes the
contested term outright, and still failed at GMV benefit R-hat 1.3447. The
reason is a *second* near-flat direction. On **treated** cells both forests are
active and both can represent the same function of `(x, t)`, so signal moves
between `alpha*mu(x,t)` and `b_z*beta(S)*nu(x,S,t)` nearly for free -- and the
reported effect again keeps only the `nu` side.

`src/longbet/_joint_forest_pair.py` blocks that direction. Conditional on the
topologies, the coding coefficients, the exposure GP and the variances, **both**
forests' leaf values are linear in the mean, hence jointly Gaussian. The block
pairs prognostic tree `j` with private treatment tree `j mod J_nu` and draws
their leaves together, cycling over all prognostic trees per sweep.

The coding coefficients are deliberately **excluded** here: `b_z*beta(S)*nu` is
a product of `b` and `nu`, so those two are *not* jointly Gaussian and must not
share a block. That is why this is a second, complementary transition rather
than one larger block -- `_joint_prognostic.py` moves `(mu leaves, b0, b1)` with
`nu` fixed, and this one moves `(mu leaves, nu leaves)` with `b` fixed.

Within each tree the leaves partition the cells, so both diagonal blocks are
diagonal and only the cross block is dense:

```text
P = [[Dmu, C], [C', Dnu]]
v ~ N(S^-1 (h_nu - C' Dmu^-1 h_mu), S^-1),  S = Dnu - C' Dmu^-1 C
u ~ N(Dmu^-1 (h_mu - C v), Dmu^-1)   given v
```

Only a `K x K` matrix is factorized. Active leaves are packed into `K = 48`
compact coordinates; a tree with more active leaves keeps its overflow
coordinates fixed, which is still a valid smaller Gibbs block because the
selection depends only on the topology the block conditions on. Under shared
treatment forests only the **private** ensemble moves; the shared contribution
stays inside `nu_fit` and is conditioned on.

Correctness evidence, `tests/test_joint_forest_pair.py`, 6 tests, all passing:
draws compared against an independently inverted dense Gaussian (mean within
4.5 MCSE, full covariance within 0.05 in correlation-scaled units) with both
partially overlapping and *identical* partitions; a check that the two forests
really are posterior-correlated (|corr| > 0.5, the ridge being crossed); and
integration checks that `mu_fit`, `nu_fit` and the raw residual identity all
survive, that `b0`/`b1` are untouched, that unreachable slots are bit-identical,
and that the per-chain path holds the invariants in every chain.

## 4. Measured effect on the unchanged chapter acceptance test

Input, target column, subgroup definitions, transformations and gates are
exactly those of `bug.md` section 5. Short settings: four dispersed chains,
2,000 burn-in, 250 retained per chain, thinning 2, 20 prognostic and 20 total
treatment trees (10 shared + 10 private, fraction .5), depth 10, `lambda_knl=2`,
adaptive coding, model seed 314159. The two arms differ **only** by
`--joint-prognostic`.

### 4.1 Matched pair at seed 314159

| Outcome / group | baseline R-hat / bulk ESS | with block R-hat / bulk ESS |
| --- | ---: | ---: |
| GMV benefit | 1.8871 / 5.8 | **1.2890 / 11.0** |
| GMV harm | 1.3467 / 9.4 | **1.0201 / 147.1** |
| hours benefit | 1.0440 / 77.3 | 1.1373 / 21.6 |
| hours harm | 1.1549 / 17.9 | 1.2567 / 11.6 |
| complaint benefit | 1.0468 / 61.9 | 1.2074 / 16.6 |
| complaint harm | 1.1434 / 19.9 | **1.0892 / 33.2** |
| **worst R-hat / min bulk ESS** | **1.8871 / 5.8** | **1.2890 / 11.0** |
| gates passed | 0 / 6 | 0 / 6 |

Fit time 476 s with the block versus 251 s without, both under concurrent
workloads; these are not controlled speed benchmarks.

The direction of the change matches the diagnosis. GMV is the outcome whose
*effect* section 2.5 showed to be damaged by the ridge, and both its subgroups
improve sharply -- the harm group by more than an order of magnitude in bulk
ESS. Hours and complaint, whose effects section 2.5 showed to be comparatively
insulated from the ridge, move slightly the other way.

**This single pair does not establish that the block helps.** With bulk ESS in
the tens, the sampling error of an R-hat estimate is itself large, and three of
six summaries moved the wrong way. Section 4.3 repeats the identical comparison
at further model seeds for that reason. Nothing here should be read as a
validated improvement, and no gate is met in either arm.

### 4.2 Coding coefficients are not the right report card

Under the block, GMV `b0`/`b1` still have R-hat 2.49 / 2.66 at bulk ESS ~5.
That is expected and is not evidence the block failed: `(b0, b1, nu)` carry an
exact sign symmetry and a scale symmetry that leave the reported effect
invariant, and the block holds `nu` fixed, so it cannot and should not traverse
them. The effect summaries above are the diagnostic that matters, exactly as
`bug.md` requires.

### 4.3 Joint decision event

`event_diagnostics.py` evaluates the chapter's three-outcome business event
(GMV >= +5%, hours <= -10%, complaint <= +1 percentage point) on **aligned joint
draws** at the common target seller/week -- never a product of marginal
probabilities. With the block, at seed 314159:

| Group | truth event share | posterior mean share | R-hat | bulk ESS |
| --- | ---: | ---: | ---: | ---: |
| benefit (n=59) | 1.0000 | .9918 | 1.0870 | 49.1 |
| harm (n=138) | 0.0000 | .0000 | -- | -- |
| all (n=750) | .3880 | .3602 | 1.0779 | 84.1 |

The harm-group indicator is **identically zero in all 4,000 draws**. As
`bug.md` warns, that is not evidence of good mixing; it means the statistic is
degenerate and the underlying effects must be inspected instead, which is what
section 4.1 does. The event probabilities are better behaved than the effect
means because thresholding at +5% / -10% is insensitive to shifts of the size
the ridge induces for most sellers -- but bulk ESS of 49 and 84 is still far
below the 400 gate, so these probabilities are **not** usable rollout guidance.

### 4.4 Is the DGP to blame? A model-generated control on the real design

The acceptance panel is hard, and `bug.md` records that it contains latent
seller drift a constant random intercept cannot represent. It is therefore fair
to ask whether the gates fail because the model cannot fit this DGP at all.
`make_model_control.py` settles that by changing exactly one thing.

It keeps the real covariates `x`, the real staggered treatment panel `z`, the
real evaluation rollout `ze`, the calendar grid, the target column, and the real
`listings`/`fulfilled` columns -- so the benefit and harm groups are the *same
sellers* with the same sparse support (five benefit-group sellers observed at
exposure 8 in week 18). Note that both subgroup variables are inside `x`:
`fulfilled` is exactly `x[:,3]`, and `x[:,2]` is `log(listings)` to machine
precision, so the model can represent the heterogeneity.

It replaces only the response-generating mechanism, with one drawn from the
fitted model class: sum-of-trees prognostic and treatment forests, a saturating
exposure curve, adaptive coding `b0 = .35`, `b1 = 1`, a Gaussian seller
intercept, and Gaussian/probit innovations with the original correlation
`[[1,.6,.45],[.6,1,.35],[.45,.35,1]]`. There is **no drift and no mean or
covariance misspecification**. Effect spread and prevalence were tuned to match:

| | real panel | model control |
| --- | ---: | ---: |
| complaint prevalence | 6.2815% | 6.2815% |
| GMV target-effect SD | .1350 | .1231 |
| hours target-effect SD | .1536 | .1542 |
| complaint risk-difference SD | .0168 | .0162 |

Same sampler, same settings, same seed 314159:

| Outcome / group | real panel R-hat / ESS | model control R-hat / ESS |
| --- | ---: | ---: |
| GMV benefit | 1.8871 / 5.8 | 1.7648 / 6.1 |
| GMV harm | 1.3467 / 9.4 | 1.2787 / 11.2 |
| hours benefit | 1.0440 / 77.3 | 1.6098 / 6.8 |
| hours harm | 1.1549 / 17.9 | 1.4608 / 8.1 |
| complaint benefit | 1.0468 / 61.9 | 1.0464 / 75.4 |
| complaint harm | 1.1434 / 19.9 | 1.0435 / 94.0 |
| **worst R-hat / min ESS** | **1.8871 / 5.8** | **1.7648 / 6.1** |
| gates passed | 0 / 6 | 0 / 6 |

Adding `--joint-prognostic` on the control improves it in the same direction
as on the real panel (worst R-hat 1.7648 -> 1.5502, min bulk ESS 6.1 -> 7.3),
so the block's mechanism is not specific to quirks of the original simulation.

**The failure reproduces in full on data the model generated itself.** The DGP
is therefore not the cause of the convergence failure: a dataset the model can
represent exactly, on the same design, produces the same worst-case R-hat and
the same minimum ESS. If anything the control is slightly harder for hours.

A second, independent argument points the same way. The block of section 3 is
an exact Gibbs kernel: same data, and provably the same target distribution. It
moved GMV harm from R-hat 1.3467 / ESS 9.4 to 1.0201 / ESS 147.1. A change that
leaves the posterior and the data identical cannot improve anything that was the
data's fault, so at least part of the failure is demonstrably algorithmic.

**What this does not settle.** It rules the DGP out as the explanation of the
*convergence* failure. It does not certify *recovery*: whether the intervals
would cover the truth at the right rate is unknowable until the chains actually
converge, and remains an open gate. On the short runs all twelve real-panel
subgroup intervals do contain the truth, but those intervals are not valid
posterior intervals while R-hat is 1.29-1.89, so that is reassurance, not
evidence. Two of the six control intervals miss, for the same reason.

### 4.5 Replicated seeds

The identical comparison at three model seeds. Data, target, subgroups, settings
and gates are unchanged; only the sampler's seed varies, and both arms always
share it.

| Seed | arm | worst R-hat | min bulk ESS | gates |
| --- | --- | ---: | ---: | ---: |
| 314159 | baseline | 1.8871 | 5.8 | 0/6 |
| 314159 | + block | **1.2890** | **11.0** | 0/6 |
| 271828 | baseline | 1.9135 | 5.7 | 0/6 |
| 271828 | + block | **1.1928** | **15.3** | 0/6 |
| 141421 | baseline | 1.7078 | 6.4 | 0/6 |
| 141421 | + block | **1.3810** | **9.1** | 0/6 |
| mean | baseline | 1.8361 | 6.0 | -- |
| mean | + block | **1.2876** | **11.8** | -- |

The worst-case R-hat improves in **3 of 3 paired seeds**, by .33 to .72, and the
minimum bulk ESS improves in 3 of 3, roughly doubling on average. That is now a
defensible claim rather than one lucky run. Total bulk ESS summed over the six
summaries is essentially unchanged (244.8 versus 227.1), so the block does not
buy general efficiency -- it specifically repairs the *worst* coordinate, which
is what a ridge-crossing move should do. **No run passes any gate.**

### 4.6 The longer run: hours un-freezes, GMV does not

`--joint-prognostic` at 8,000 burn-in, 500 retained, thinning 4 (10,000
iterations per chain, 2,049 s):

| Outcome / group | R-hat | bulk ESS | tail ESS |
| --- | ---: | ---: | ---: |
| GMV benefit | 1.3643 | 9.1 | 31.6 |
| GMV harm | 1.0912 | 29.1 | 306.0 |
| hours benefit | 1.0590 | 47.8 | 319.7 |
| hours harm | 1.0458 | 63.4 | 390.0 |
| complaint benefit | 1.1048 | 30.4 | 166.2 |
| complaint harm | 1.0883 | 36.2 | 556.9 |

Against the previous review's longer run, hours improves markedly (13.3 -> 47.8
and 15.1 -> 63.4 bulk ESS) and GMV does not (8.5 -> 9.1, 31.8 -> 29.1).

The diagnostic that matters is how ESS responds to a 4x increase in iterations,
comparing the short and long runs of the *same* configuration:

| Outcome / group | short ESS | long ESS | ratio |
| --- | ---: | ---: | ---: |
| GMV benefit | 11.0 | 9.1 | **0.82** |
| GMV harm | 147.1 | 29.1 | **0.20** |
| hours benefit | 21.6 | 47.8 | 2.21 |
| hours harm | 11.6 | 63.4 | 5.45 |
| complaint benefit | 16.6 | 30.4 | 1.83 |
| complaint harm | 33.2 | 36.2 | 1.09 |

For hours, ESS now grows with computation: the frozen mode is gone and only
ordinary slow mixing remains. For **GMV it does not grow at all** -- a
frozen mode survives the block. (The 147.1 in the short GMV harm run was a
favorable draw; at these ESS levels single numbers are noisy, which is why
section 4.5 replicates seeds.)

### 4.7 The second block did not help

Composing `--joint-forest-pair` on top of `--joint-prognostic` gave worst R-hat
**1.6905** (GMV harm, bulk ESS 6.4), against 1.2890 without it. The block is a
verified-correct transition -- section 3b's oracle tests are unambiguous -- it
simply does not improve mixing here, and in this matched run it hurt. It is
retained for the record and **not recommended**. Do not read this as evidence
that ridge B does not exist; it is evidence that pairing one prognostic tree
with one treatment tree is not an effective way to cross it.

### 4.8 The estimand is not the binding constraint

A natural hypothesis after section 4.6: GMV stays frozen because the chapter's
target is a whole-population counterfactual at `(S=8, week 18)`, where only the
week-11 cohort has factual support (five benefit-group sellers). If so, asking a
supported question instead would fix it. **It does not.** Three nested
restrictions were tested, each strictly better supported than the last.

*Restrict to factually treated sellers*, from the saved target-cell draws:

| | all sellers | factually treated only | never treated only |
| --- | ---: | ---: | ---: |
| GMV benefit, baseline | 1.887 / 5.8 | 1.865 / 5.9 | 1.900 / 5.8 |
| GMV harm, baseline | 1.347 / 9.4 | 1.334 / 9.7 | 1.370 / 9.0 |

Unit-level extrapolation is not the problem: dropping the 26 never-treated
benefit sellers changes nothing.

*Restrict to factually treated cells* (`estimand_test.py`, with the block):
estimand A is the chapter target; estimand B is the ATT over cells actually
observed treated, with no counterfactual extrapolation anywhere.

| Outcome / group | A: target cell | B: supported ATT | B cells |
| --- | ---: | ---: | ---: |
| GMV benefit | 1.289 / 11.0 | 1.294 / 10.9 | 397 |
| GMV harm | 1.020 / 147.1 | 1.101 / 28.1 | 1150 |
| hours benefit | 1.137 / 21.6 | 1.238 / 13.2 | 397 |
| hours harm | 1.257 / 11.6 | 1.242 / 12.4 | 1150 |
| complaint benefit | 1.207 / 16.6 | 1.209 / 16.3 | 397 |
| complaint harm | 1.089 / 33.2 | 1.186 / 16.8 | 1150 |

No improvement anywhere, and four of six are slightly worse.

*Restrict to one observed exposure* -- the plainest event-study quantity, an
average over 33 directly observed benefit-group cells. GMV at exposures 1 to 5
gives R-hat 1.61, 1.52, 1.52, 1.68, 1.58 at bulk ESS 7 to 8.

So the simplest fully supported quantity in the model fails as badly as the
extrapolated one. Combined with section 4.4, every knob other than the sampler
is inert: changing the data-generating mechanism to one the model can represent
exactly, restricting to supported units, restricting to supported cells, and
restricting to a single supported exposure all leave the failure intact, while
changing only the sampler -- with the posterior provably unchanged -- is what
moves it. **The binding constraint is posterior exploration, not the question
being asked.**

### 4.9 Where that leaves the diagnosis

Ridges A and B are linear-Gaussian directions and blocking A demonstrably works
for hours. Neither block touches the **discrete tree-partition space**, and that
is now the only major structure left unaddressed. The supporting evidence is
consistent: the archived `--treatment-only --fixed-beta` variant -- a plain BART
treatment forest with no coding ridge and no GP ridge at all -- still failed at
R-hat 1.3512, and section 4.8 shows the failure survives every reduction of the
estimand. This is a hypothesis with converging evidence, not an established
fact; section 8 of `bug.md` lists the cheap tests that would confirm it before
anyone builds whole-tree moves.

## 5. Reproduction

```sh
audit=/home/ignacio/longbet-bug-audit     # disk-backed; /tmp here is a tmpfs
data=/tmp/longbet-reference-benchmark.9G6Fop/reexported-input.npz

# matched baseline: identical settings, no experimental block
.venv/bin/python benchmarks/bench_multi_chapter.py --data "$data" \
  --output "$audit/chapter-baseline-short" --mode multi \
  --shared-trees 10 --shared-variance-fraction .5 \
  --burnin 2000 --draws 250 --skip 2 --chains 4 --seed 314159

# same, plus the exact joint prognostic-leaf / coding block
.venv/bin/python benchmarks/bench_multi_chapter.py --data "$data" \
  --output "$audit/chapter-jointprog-short" --mode multi \
  --shared-trees 10 --shared-variance-fraction .5 --joint-prognostic \
  --burnin 2000 --draws 250 --skip 2 --chains 4 --seed 314159

# the decomposition that measures the ridge directly
.venv/bin/python "$audit/decompose_fit.py" --data "$data" \
  --output "$audit/decompose-short" --burnin 2000 --draws 250 --skip 2 \
  --chains 4 --seed 314159

.venv/bin/python -m pytest tests/test_joint_prognostic.py -q
```

## 6. Test and provenance record

* Full Python suite in this environment: **247 passed, 1 skipped, 4 warnings**
  (835 s). The skip is the optional external comparison. No regressions.
* New tests: `tests/test_joint_prognostic.py` (7) and
  `tests/test_joint_forest_pair.py` (6), both validated against independently
  constructed dense Gaussian oracles.
* Final Python/R/contract source fingerprint: `bcf76115de639f38b8b79f3f0651dc957dfca49e397ab779e28549d063803337`
  (previous handoff: `0a20d5420c3826a3fb1b15b480aed8cb712fe19e8c8cb7d9f1bdb5e980a48069`;
  the difference is exactly the two added modules).
* Changed files this pass: `src/longbet/_joint_prognostic.py` (new),
  `src/longbet/_joint_forest_pair.py` (new), `tests/test_joint_prognostic.py`
  (new), `tests/test_joint_forest_pair.py` (new),
  `benchmarks/bench_multi_chapter.py` (two opt-in flags),
  `benchmarks/slow_mode_review.md` (this file), plus `bug.md` and
  `/home/ignacio/book/multi.md` as handoff/plan updates. The 28 previously
  modified tracked files are untouched.
* No commit, push, image rebuild, chapter render, or public default change was
  performed. No R suite was run in this pass; the previous review's 268
  assertions are historical.
