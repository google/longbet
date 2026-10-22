# Whole-tree refresh investigation (fifth pass, 2026-09-07)

This pass is complete without a reliable inference repair. Both whole-tree
trials fail all six gates. A subsequent mathematical audit proved that the
default chapter model has an improper posterior. Both explicit proper-prior
trials also fail all six gates. See the latest section of `bug.md`; outputs are under
`/home/ignacio/longbet-bug-audit/refresh-pass` and `proper-pass`.

## Mathematical target and implementation

The fourth-pass root proposals retain tree shape and only alter prognostic
rules. This experiment also permits a change in tree shape and updates private
treatment partitions. It does not presume topology is the sole slow coordinate:
the previous fixed-topology controls already refuted that claim.

For fixed coding, GP, alpha, covariance and augmented responses, write the
outcome's mean as `offset + A(T) x + U g`, with independent standard-normal
whitened leaves `x` and unit effects `g`. The full SUR conditional supplies
precision `Omega[m,m]` and score `(Omega r)[m]`. The existing exact Schur
complement in `_collapsed_exposure.py` integrates every prognostic, private
treatment and shared treatment coefficient for that outcome, plus its units.
Other outcome coefficients, and the shared partitions, are held fixed.

Select one private tree uniformly from a forest and independently draw its
complete shape and rules from the original geometric tree prior `p(T)`.
At each reachable node, use its actual `p_nonterminal`; choose among the
ancestor-eligible unblocked variables with their conditional `log_s` weights,
then uniformly among that variable's ancestor-restricted cutpoints. Terminate
when no geometric rule is available or maximum depth is reached. Sample counts
do not change the proposal's rule or terminal probabilities.

If the candidate violates the original count support, reject it once. For
supported candidates, the prior and independence-proposal terms cancel:

```
p(T') L(T') p(T) / [p(T) L(T) p(T')] = L(T') / L(T).
```

`L` is the ALL-forest/unit integrated likelihood, so changes in other forest
coefficients can compensate for a changed partition. Immediately redraw every
integrated coefficient after either acceptance or rejection. This proposal
preserves the existing unnormalized target; it does not change the model or
its priors. This invariance cannot give a finite normalizer to an improper
target. The same derivation applies when proper variance priors are specified.

The new module is `src/longbet/_refresh_trees.py`. The breadth-first sampler
only visits reachable nodes. When shape changes, the ENTIRE compact mapping,
including leaf indices, prior scales and treatment-column flags, switches to
the candidate before drawing coefficients. Both private forests' membership,
pending prunes, count and precision caches are synchronized; the selected
tree's affluence is recomputed. Shared coefficients and cached fits remain
aligned. Raw residuals account for the stored float32 coefficients.

The compact-system limit rejects proposed overflow; current overflow makes
this optional transition an identity. It never accepts or draws from an
integral missing an active coefficient. The ordinary sweep still has the full
tree support. The new step retains fixed beta/coding/alpha/covariance settings,
binary variance identification and supported missingness.

This is an independence Metropolis move, not Particle Gibbs. The
[Particle Gibbs BART paper](https://proceedings.mlr.press/v38/lakshminarayanan15.html)
motivates considering complete-tree proposals but does not establish the
correctness or efficiency of this implementation. Its conditional-SMC proposal
has not been implemented here.

## Prespecified comparisons and evidence

Opt in through `bench_multi_chapter.py --refresh-trees`. Every
`--collapsed-every` sweeps (default 10), the step refreshes one prognostic and
one private treatment tree for every outcome. The eight `refresh.npz` fields
are attempted, accepted, log likelihood ratio, current total leaves, candidate
total leaves, count-supported, accepted shape change and accepted root change.
Invalid candidates have log ratio negative infinity; this is an explicit
rejection diagnostic, not a nonfinite fitted state.

Initial independent tree-prior tests: seven passed in 7.36 s without warnings.
Initial end-to-end checks: four passed in 135.93 s. They enumerate 81 combined
private-tree configurations in two-outcome models, checking the actual
Gaussian posterior and an observed probit/Gaussian posterior via independently
constructed dense covariance matrices and numerical orthant probabilities.
Gaussian fitted means and full covariance also agree. Chained tests cover
sharing, missingness, zero/negative weights, changed shapes, fixed beta, raw
residuals and a subsequent ordinary sweep. These are correctness checks, not
original-data convergence evidence. The expanded suite, adding deliberate
nonzero offsets and the independent capacity-restricted Gaussian posterior,
passed **15 tests in 143.20 s**, with no failures, skips or warnings.
The broader regression floor passed 92 tests in 909.21 s, also without
failures, skips or warnings, for 107 distinct targeted tests. The two chained
tests additionally passed direct NumPy membership, count, precision and
affluence reconstruction (141.20 s; three unrelated cases deselected).

Exact commands for those 107 distinct checks and the two overlapping reruns:

```bash
taskset -c 0-3 .venv/bin/python -m pytest \
  tests/test_refresh_tree_prior.py tests/test_refresh_trees.py tests/test_multi_benchmark.py -q
taskset -c 8-11 .venv/bin/python -m pytest \
  tests/test_full_sur.py tests/test_shared_forest.py tests/test_joint_gaussian.py \
  tests/test_unit_interweave.py tests/test_joint_prognostic.py tests/test_joint_forest_pair.py \
  tests/test_multi_step.py tests/test_collapsed_exposure.py tests/test_coding_interweave.py \
  tests/test_change_prognostic.py tests/test_collapsed_tree.py tests/test_topology_diagnostic.py -q
taskset -c 0-3 .venv/bin/python -m pytest tests/test_refresh_trees.py -k chained -q
```

Their logs are `refresh-pass/expanded-tests.log`, `regression-tests.log` and
`cache-tests.log`. These are targeted checks, not a full Python/R suite.

The original-data smoke uses 30 sweeps, four chains and seed 314159. Full trials
are prespecified at seeds 314159 and 271828, four chains, 2000 burn-in, 250 retained
draws, skip 2, 20 prognostic and 20 total treatment trees including 10 shared.
The only addition to the matching fourth-pass coding candidates is whole-tree
refresh; joint-prognostic and collapsed coding remain enabled.

The first full trial completed in 806.11 s fit time (14:42.42 including
reporting/prediction), peaking at 4.25 GiB RSS. **None of the six gates passed**:
worst R-hat 1.2351, minimum bulk ESS 12.8 and minimum tail ESS 78.0. Retained
sampling accepted 470/1200 proposals, with 121 shape changes, 179 root changes
and 16 count-invalid proposals. No capacity overflow occurred. The GMV benefit
interval [.2393, .3623] misses its truth .2250.

The second trial completed in 783.74 s fit time (14:16.69 overall), peaking
at 4.23 GiB RSS. It also passes **0/6 gates**: worst R-hat 1.1626, minimum
bulk ESS 17.6, minimum tail ESS 30.7. Retained sampling accepted 496/1200
proposals, with 130 shape changes, 188 root changes and eight count-invalid
proposals; no capacity overflow occurred. Its GMV benefit interval
[.2308, .3329] also misses truth. Numerical improvement over the old coding
baselines (worst R-hat 1.6445 and 1.3534) cannot establish an efficiency
ranking at such low ESS.

`compare_refresh_runs.py --require-complete --validate-checkpoints` validates
all four matching archives. The new checkpoints reconstruct fitted means and
raw residuals from actual rules within maximum absolute error 5.82e-5, below
the 3e-4 tolerance, with binary variance exactly one. The independently
recomputed summaries are in `refresh-pass/comparison.json`; the trace figure
is `refresh-traces.png`/`.pdf`. It shows sustained between-chain differences.

Data and starting source fingerprints, preserved handoff and exact environment
are recorded in `bug.md`. No dependency or public API/default was changed.

## Variance-prior defect and explicit model trial

The later audit proved that the original chapter target has no finite
normalizing constant under `sigma_prior_a=sigma_prior_b=0`. The proof concerns
the observed probit/Gaussian likelihood, the actual data and original allowed
tree support. See [`variance_prior_review.md`](variance_prior_review.md) for
the derivation, independent certificate and reproduction commands. This does
not establish that the zero-variance boundary caused the failed finite-run
mixing diagnostics.

The next prespecified comparison uses the same data, sampler, two seeds and
budgets with explicit `--sigma-prior-a 2 --sigma-prior-b 1`. Proper IG(2,1)
continuous innovation variances remove the normalization defect. They change
the model and remain opt-in through existing config fields. Binary variance
remains exactly one. Results are under `proper-pass/proper-<SEED>`:

| Seed | Worst R-hat | Minimum bulk ESS | Minimum tail ESS | Gates |
| --- | ---: | ---: | ---: | ---: |
| 314159 | 1.6885 | 6.4 | 19.9 | 0/6 |
| 271828 | 1.3311 | 9.7 | 47.8 | 0/6 |

Fit times were 779.15/781.14 s; overall wall times were 14:11.76/14:13.50,
with peak RSS 4.24/4.26 GiB. Both checkpoints reconstruct within 5.81e-5
maximum residual error (3e-4 tolerance), and both retain binary variance one.
Retained accepted shape changes number 121/125, root changes 175/188; no
capacity overflow occurred. All six gates also fail using only the latter
half of each retained trace. The second GMV benefit interval [.2327, .3404]
misses truth .2250. Final diagnostics and labelled trace figures are in
`proper-pass/comparison.json` and `refresh-traces.png`/`.pdf`.

This repairs normalization for the explicit model but does not establish a
reliable inference procedure. The pass stops this candidate family. Neither
experimental change is promoted. Long-budget scaling, repeated-data coverage,
SBC, joint-event calibration, full Python/R tests, GPU verification, public
integration and chapter rendering were not performed for these failed
candidates. The chapter gates remain required.

Exact full-fit command, from the repository root (substitute an unused output
directory and one of the prespecified seeds):

```bash
/usr/bin/time -v taskset -c 0-7 .venv/bin/python benchmarks/bench_multi_chapter.py \
  --data /home/ignacio/longbet-bug-audit/topology-pass/chapter-input.npz \
  --output /tmp/longbet-proper-314159 --mode multi --shared-trees 10 \
  --joint-prognostic --collapsed-exposure --collapsed-proposal coding \
  --refresh-trees --save-topology --burnin 2000 --draws 250 --skip 2 \
  --chains 4 --seed 314159 --progress-every 250 \
  --sigma-prior-a 2 --sigma-prior-b 1
```

Omit only the two sigma-prior flags to reproduce an original-target refresh
trial. The recorded trials ran sequentially; each used about 4.3 GiB peak RSS.
Runtimes are contended and do not support fair efficiency comparisons.
Validate the existing completed audit layout using:

```bash
.venv/bin/python benchmarks/compare_refresh_runs.py --require-complete --validate-checkpoints
.venv/bin/python benchmarks/compare_refresh_runs.py \
  --root /home/ignacio/longbet-bug-audit/proper-pass --proper-prior \
  --require-complete --validate-checkpoints
.venv/bin/python benchmarks/plot_refresh_traces.py \
  --root /home/ignacio/longbet-bug-audit/proper-pass --runs proper-314159 proper-271828
```

The comparator verifies archive hashes, all outcomes, explicit chain axes,
paired settings and the stated change from zero hyperparameters to IG(2,1).
It marks that comparison as a model change. These scripts do not fit models.

## Conditional SUR coupling audit

An independent read-only analysis of both fourth-pass ensemble parameter
archives finds moderate residual correlations: approximately .480, .369 and
.600 for complaint/GMV, complaint/hours and GMV/hours. For complete observations
with every nonlinear/discrete coordinate and variance fixed, the normalized
cross-outcome Gaussian precision blocks have operator norms bounded by the
absolute noise partial correlations. The nonnegative Gauss-Seidel comparison
matrix therefore bounds the exact outcome-wise Gaussian mean-relaxation factor
by .438-.504 over the saved draws. The derivation and per-chain evidence are in
`refresh-pass/sur-coupling.json`; `analyze_sur_coupling.py` reproduces it.

This disfavors severe SUR conditioning alone under those fixed conditionals.
It says nothing about the mixing of local forest updates, binary augmentation,
nonlinear products, changing partitions, or full-chain treatment-effect ESS.
The underlying ensemble runs failed convergence, so the pooled parameter
ranges remain descriptive. Some correlation R-hats also narrowly miss 1.01.
