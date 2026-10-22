# Topology control and collapsed exposure experiment (2026-09-07)

**Status: no reliable repair established.** The candidate failed at both tested
proposal scales, including the prespecified two-seed adjustment. No new
original-posterior run passed any of the six subgroup convergence gates.
This follows `slow_mode_review.md`; all older work is preserved.

## Reproduction and topology measurement

The original input hash and the original Python/R/contract source hash matched
`bug.md`. Python 3.12.3, JAX 0.11.1, bartz 0.12.1, numpy 2.5.2, equinox 0.13.8;
CPU only. Output directory: `/home/ignacio/longbet-bug-audit/topology-pass`.
Fits use affinity 0-7; independent tests use 8-15. `/tmp` is RAM-backed here.
The audit directory now also contains a byte-identical durable `chapter-input.npz`
copy, SHA256 `6e4904482c10530728705f0c4c1ce8bc943db6c2f0aa6aedfba5bb99284aaffd`.

The matched instrumented run uses four dispersed chains, seed 314159,
2000 burn-in / 250 retained / thinning 2, 20 prognostic and 20 total treatment
trees with 10 shared, joint-prognostic enabled, all other benchmark defaults.
It took 529.3 s and exactly reproduced the archived joint-prognostic draws and
diagnostics. Its GMV benefit R-hat is 1.2890 / bulk ESS 11.0; no subgroup passes.

New instrumentation saves compressed rules for every retained tree and counters
from those retained sweeps, plus final topology for the fixed-tree experiment.
Previous benchmark archives did not save forests, despite the earlier handoff's
suggestion to analyze them without another fit. The new `topology.npz` is about
677 KiB compressed; final topology is about 12 KiB.

GMV prognostic root rules never change between saved draws in **19,18,18,19 of
20 trees** across the four chains. Only 9-41 root-rule changes occur in total
per chain; complaint has 400-556, with only 1-4 unchanged roots of 20. Grow and
prune acceptance rates for GMV prognostic forests are about .025-.042; GMV
treatment forests also have low acceptance and different rule histograms across
chains. These are observations, not evidence that a particular fixed root is
wrong. Between-save counts miss reversals between two retained draws.

`--fixed-topology` conditions on all prognostic, private treatment, and shared
treatment rules from reference chain 0. Only rules/membership are copied, not
the chain's leaves, GP, coding, units, or variances. Standard independent prior
dispersion is retained. The control vetoes topology proposals with process-local
hooks before tracing and keeps the existing conditional leaf draws. A mixed
two-chain integration test verifies fixed rules, nonzero moving leaves, dispersed
GP starts, fitted-value/residual consistency, and binary variance identification.
It samples a conditional posterior, not the full model posterior; even perfect
diagnostics would not certify the chapter estimator.

## A different move through the multiplicative ridge

The archived exposure GP traces stay in different *shapes*, not just different
global signs/scales. Previous Gaussian blocks hold treatment leaves fixed when
moving beta and hold beta fixed when moving leaves. The candidate instead
integrates all linear coefficients before proposing beta.

At fixed tree rules, b, alpha, residual covariance and augmented responses, take
one outcome's full-SUR conditional precision w and conditional response y.
Whiten all its active prognostic, private-treatment and shared-treatment leaves
to x ~ N(0,I), and all its enabled unit intercepts to g ~ N(0,I):

```
y = A(beta) x + unit_sd * g[unit] + epsilon
```

Prognostic design columns are alpha times leaf membership times leaf prior SD;
treatment columns use b(Z)*beta(S) times membership times the correct private
or shared leaf prior SD. Offsets are subtracted explicitly. Other outcomes'
shared leaf values remain fixed; their influence enters through the complete
SUR conditional score. Current shared leaf priors are diagonal across outcomes.

Let U be the whitened unit design and W = diag(w). The unit precision
D = I + U'WU is diagonal. Define

```
C = U'WA,    h_g = U'Wy
Q = I + A'WA - C'D^-1 C
h = A'Wy - C'D^-1 h_g
log L(beta) = -.5 [log|D| + log|Q| + y'Wy
                  - h_g'D^-1 h_g - h'Q^-1 h]
```

The omitted observation normalization is fixed during this transition. Propose
`beta_new = sqrt(1-rho^2)*beta + rho*K_chol*eta`, eta standard normal. This pCN
proposal is reversible with respect to the GP prior, so its MH acceptance is
`min(1, L(new)/L(old))`. The Gaussian prior/proposal cancellation can also be
checked directly from the symmetric joint Gaussian law of old/new beta; the
general approach is described by [Cotter et al. (2013), v3](https://arxiv.org/abs/1202.0709v3).
That source does not validate this LongBet adaptation.

Immediately after accept **or reject**, draw x using Q and h, then g conditional
on x. Restore physical leaf units, private/shared fits, gamma and raw residuals
together. Redrawing on rejection is necessary: the MH step targets beta's
marginal conditional, which no longer conditions on the old coefficients.
This is the same posterior and the same model, with no Gaussian approximation
to beta*nu, jitter, relaxed accept/reject thresholds, or new priors.

Implementation: `src/longbet/_collapsed_exposure.py`. All active leaves must fit
the packed system (benchmark capacity 128; module default 192). Overflow skips
the entire move, and selection depends only on topology, so it remains an exact
identity transition. Nothing is removed from the model. The benchmark records
every attempt, acceptance, marginal log ratio and active-leaf count in
`collapsed.npz`. The move is appended every 10 sweeps, after joint-prognostic,
with rho=.3. Rho=0 gives a useful all-leaves/unit-Gibbs-only control.

## Independent correctness evidence

Six tests passed in `tests/test_collapsed_exposure.py` (74.40 s for the first
five; 16.51 s for the sixth in a separate run):

* Log integral versus independently formed observation covariance V, using
  dense inversion/determinants; includes signed/collinear columns, absent
  likelihood cells, disabled units, dummy coordinates and a prior-only unit.
* 24,000 actual coefficient draws versus an independent dense posterior,
  checking means within 4.5 Monte Carlo standard errors and the full covariance
  to .05 in correlation-scaled units.
* Actual two-chain, three-outcome transitions with and without sharing,
  supported missingness, untouched latent values/variance identification and
  unreachable leaf slots, exact rules, and correct raw residual/fit identities.
* Fixed beta/intercepts, and an exact identity transition on capacity overflow.
* Actual mixed **observed probit/Gaussian posterior with sampled beta**, shared
  and private treatment leaves, versus independent numerical quadrature. The
  oracle integrates the GP prior by Gauss-Hermite quadrature and the two means
  over a dense grid, using the continuous density and conditional probit CDF.
  It checks mean/covariance, GP second moments, and a joint event against 18,000
  retained sampler draws. This tests the nonlinear product, not just a Gaussian
  latent-response stand-in. It is a small posterior oracle, not general SBC.

The topology-control test and three existing benchmark tests also pass (52.14 s).
An initial diagnostic-report attribute lookup error was corrected before the
full fit reached reporting. A 30-sweep full-data smoke test completed in 40.1 s
with finite output and no packing overflow; it is not inference evidence.

The handoff's regression floor plus the new topology-control test passed as a
single run: **66 passed in 368.78 s**, no skips/warnings. The full
Python/R suites have not been rerun in this pass; no public API has changed.

## Commands and pending results

Run from the engine repository, with the existing `.venv`. Output directories
must be new. The fixed run consumes the instrumented run's final topology.

```sh
audit=/home/ignacio/longbet-bug-audit/topology-pass
data="$audit/chapter-input.npz"
taskset -c 0-7 .venv/bin/python benchmarks/bench_multi_chapter.py \
  --data "$data" --output "$audit/instrumented-short" --mode multi \
  --shared-trees 10 --joint-prognostic --save-topology \
  --burnin 2000 --draws 250 --skip 2 --chains 4 --seed 314159
taskset -c 0-7 .venv/bin/python benchmarks/bench_multi_chapter.py \
  --data "$data" --output "$audit/fixed-short" --mode multi \
  --shared-trees 10 --joint-prognostic --save-topology \
  --fixed-topology "$audit/instrumented-short/final_topology.npz" \
  --burnin 2000 --draws 250 --skip 2 --chains 4 --seed 314159
```

The fixed-topology control completed in 490.4 s. Every saved split rule is fixed
across time and chains, as intended. Results (R-hat / bulk ESS / tail ESS):

| Outcome/group | Unrestricted joint-prognostic | Common fixed topology |
| --- | --- | --- |
| GMV benefit | 1.2890 / 11.0 / 42.8 | 1.2150 / 14.4 / 54.1 |
| GMV harm | 1.0201 / 147.1 / 285.5 | 1.0020 / 763.8 / 676.6 |
| Hours benefit | 1.1373 / 21.6 / 240.3 | 1.0339 / 205.6 / 738.9 |
| Hours harm | 1.2567 / 11.6 / 72.3 | 1.0009 / 915.5 / 924.7 |
| Complaint benefit | 1.2074 / 16.6 / 80.8 | 1.0158 / 148.9 / 430.6 |
| Complaint harm | 1.0892 / 33.2 / 288.5 | 1.0211 / 180.0 / 379.8 |

Tree topology is therefore **not the only remaining slow coordinate**. The
fixed-tree GMV benefit interval [.22853,.32717] also misses truth .22503, though
nonconvergence and conditioning on a single topology prohibit posterior recovery
claims from this control. Two conditional summaries passing does not establish
full-model convergence. These are different target distributions, not a fair
comparison of estimator efficiency or accuracy.

This motivated testing the collapsed nonlinear exposure move on the original
unrestricted posterior. `collapsed-short` matches the first command
above plus `--collapsed-exposure` (rho=.3, period 10, capacity 128). It completed
in about 680 s, with no capacity skips (maximum 115 leaves), and **0/6 gates**:

| Outcome/group | R-hat | Bulk ESS | Tail ESS |
| --- | ---: | ---: | ---: |
| GMV benefit | 1.0869 | 31.4 | 498.8 |
| GMV harm | 1.1321 | 21.0 | 147.2 |
| Hours benefit | 1.2016 | 14.8 | 100.7 |
| Hours harm | 1.4120 | 8.5 | 41.5 |
| Complaint benefit | 1.1933 | 15.7 | 71.8 |
| Complaint harm | 1.0940 | 37.3 | 222.9 |

At these effective sample sizes, the differences from the baseline do not
establish an efficiency ranking. The GMV benefit interval [.24221,.36033]
misses .22503; its apparently certain positive sign is not a convergence
certificate. All individual effects and summaries are preserved, including
failed results.

The fixed-topology/rho=.3 control also completed (669.5 s). Its GMV benefit
remains unconverged at 1.1383 / 20.9, with GMV harm 1.0055 / 691.1, hours benefit
1.0186 / 495.5, hours harm 1.0026 / 952.3, complaint benefit 1.0147 / 138.7, and
complaint harm 1.0080 / 218.9 (R-hat / bulk ESS). Only GMV harm and hours harm
pass all three conditional thresholds. Thus this new step also fails to
eliminate the continuous-coordinate problem when tree rules are held fixed.

GP acceptance in retained sampling is 30.5% for complaint, 3.5% for GMV, and
5.5% for hours. GMV's median marginal log likelihood ratio is -15.7 (10th/90th
percentiles -68.2/-4.0). Poor acceptance motivates one scale adjustment: rho=.05
with the same frequency and retained budget, **prespecified seeds 314159 and
271828**. The matching old joint-prognostic baselines already exist. Both
adjusted trials completed, in 736.3 / 729.5 s, and both pass **0/6 gates**:

| Outcome/group | Seed 314159: R-hat / bulk / tail ESS | Seed 271828: R-hat / bulk / tail ESS |
| --- | --- | --- |
| GMV benefit | 1.1793 / 16.7 / 84.9 | 1.1903 / 16.2 / 134.9 |
| GMV harm | 1.0710 / 37.4 / 585.4 | 1.4214 / 8.4 / 30.5 |
| Hours benefit | 1.3309 / 10.0 / 115.6 | 1.2128 / 14.4 / 201.0 |
| Hours harm | 1.0891 / 33.1 / 163.9 | 1.2143 / 13.7 / 113.9 |
| Complaint benefit | 1.1407 / 27.2 / 72.2 | 1.2046 / 14.7 / 32.6 |
| Complaint harm | 1.0901 / 35.7 / 257.9 | 1.1908 / 16.5 / 135.4 |

Retained-phase GP acceptance is now complaint 84.5% / 81.5%, GMV 63.0% / 47.5%,
and hours 71.0% / 68.5%. Better acceptance does not deliver reliable effect
draws. No packing overflow occurred in any full trial. The minimum bulk ESS
is 10.0 / 8.4. These runs do not establish an efficiency ranking at such low
ESS, and they provide no basis for recommending this kernel as a repair.

This bounded trial family was stopped here: no long-budget ESS-scaling test,
repeated-data coverage/SBC, joint-event calibration, full Python/R suite, or
chapter render was performed for this candidate. Those gates remain unmet,
rather than assumed to follow from the independent small posterior oracle.

## Reusable artifacts and next work

The complete six-run comparison, including every mean, interval, RMSE, sign
probability and failure, is `topology-pass/comparison.json`. Some fits overlapped
on separate CPU affinities, with about 4 GiB peak RAM each; the recorded runtimes
are contended, not fair speed comparisons. Source/data/environment provenance
is recorded in `bug.md` and the audit manifests. The final Python/R/contract
fingerprint is `62b26d4d1085d8ab52ba9291d804e3191501553bc55d7751b18ad4d48eed3735`;
Python-only is `c50736eb72e33fb95ca8bd535c71871ee5d5d8f0793fa408e004c3436f28983d`.

The two rho=.05 runs have an internal `final_state.eqx` checkpoint (about
22.7 MiB) in addition to retained rules and effect draws. Earlier runs predate
this addition. `benchmarks/load_multi_benchmark_state.py` loads the checkpoint
using the report settings and hash-checked original data, without running MCMC.
Its real-data load was verified: raw residual errors below 5e-5 and all four
binary variances exactly one. It is not a public fit archive or a complete
posterior trace. The diagnostic callback stores acceptance/log-ratio records
in float32; matrix calculations and acceptance decisions use float64.

Example for the adjusted runs (choose an unused output directory):

```sh
taskset -c 0-7 .venv/bin/python benchmarks/bench_multi_chapter.py \
  --data "$data" --output "$audit/collapsed-smallstep-314159" --mode multi \
  --shared-trees 10 --joint-prognostic --collapsed-exposure --collapsed-scale .05 \
  --save-topology --burnin 2000 --draws 250 --skip 2 --chains 4 --seed 314159
.venv/bin/python benchmarks/load_multi_benchmark_state.py \
  "$audit/collapsed-smallstep-314159" --data "$data"
```

The second adjusted run changes seed/output to 271828 and used affinity 8-15.
Do not repeat the exclusive-topology diagnosis. One concrete untested direction
is proposing coding coefficients under the collapsed likelihood: the new move
conditions on b0/b1, and the old coding block only joins one prognostic tree at
a time. All-leaf integration before moving coding could access a different
direction, but no such repair is established. The rho=0 all-leaves/unit-only
ablation was not fitted either. Use the saved states to inspect geometry before
building another expensive variant. Whole-tree proposals remain possible;
the present controls do not establish that they alone will solve the problem.

Public defaults, Python/R archive semantics, contracts, and the chapter warning
are unchanged.
