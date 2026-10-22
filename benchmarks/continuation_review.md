# Multiple-outcome continuation, 2026-09-08

The inference repair remains unresolved. This continuation fixes identifiable
implementation errors and prevents reuse of the documented improper joint
target. These changes do not establish that the chapter's effects or joint
probabilities are reliable. Work stays in `development/multiple-outcomes`, at
`/home/ignacio/longbet-chapter-review/development-engine`; `main` is unchanged.

## Correctness repairs

The inherited uncommitted changes were preserved in
`/home/ignacio/longbet-bug-audit/continuation-20260908/inherited.patch` before editing.

* Removed the deterministic reflection of beta and treatment leaves into
  `sum(beta) >= 0`. The declared parameter prior is an unrestricted Gaussian.
  Its Geweke test failed with a whitened-beta mean **10.30 MCSE** from zero.
  A symmetry-reduced representation would need an explicit parameter contract;
  preserving a treatment product alone does not validate the exposed beta law.
  This failure does not itself prove bias in every symmetry-invariant effect.
* Restored Gaussian prior starts for beta. The inherited
  `max(1 + 0.2 * prior_draw, 0.1)` initialization suppressed dispersion and
  sign-changing trajectories. A new test checks whitened means/covariances and
  both signs, rather than merely checking that chains differ.
* Restored documented fixed coding `b0=b1=1` and calendar-time visibility in the
  treatment forest. The inherited `b0=0` change was a different model and broke
  the independent observed-likelihood quadrature tests. `split_time_trt=False`
  again removes exposure splits while preserving calendar-time splits.
* Retained the reordered `mu -> alpha -> nu -> beta` sweep. It is a valid
  conditional ordering; beta must use the **updated** treatment forest. The
  independent sufficient-statistic test now reflects that schedule.
* Replaced residual differencing for forest fits with sums over cached leaf
  memberships. When a treatment weight is tiny, `R/w` is large and float32
  subtraction loses leaf updates. The existing binary/alpha residual test
  exposed a **1.8153e-4** discrepancy from independently traversed trees.
  New regressions cover weights `0`, `2e-6`, `-2e-6`, and `0.2`. Physical fit
  deltas update the full residual without traversing the trees again. This
  extends the stable calculation already used by the shared-tree path.

The last change addresses a numerical error, not a new prior or likelihood.
Its effect on large-panel mixing must be measured independently. No sole cause
of the remaining chain disagreement has been established.

## Proper target and persistence

Coupled fits containing continuous outcomes now require explicitly positive
`sigma_prior_a` and `sigma_prior_b`. The scalar defaults and uncoupled/binary-only
reductions are retained. An example for standardized outcomes is IG(2,1), with
prior mean innovation variance one. This is an explicit model choice, not a
validated universal default or a convergence repair.

Why this suffices for normalization: with binary innovation variances fixed at
one and zero incoming binary loadings, bound each observed continuous Gaussian
density by `(2*pi*v_m)^(-1/2)` and integrate the binary latent densities over
their sign regions. Thus the observed likelihood is bounded by
`product_m (2*pi*v_m)^(-n_m/2)`. All negative moments of IG(a,b) are finite when
`a,b > 0`, and the remaining priors are proper. The likelihood is positive for
finite parameters and positive variances. Its normalizer is finite and positive.
This argument addresses propriety; it gives no mixing or calibration guarantee.
The previous proof of impropriety under the chapter's `(0,0)` prior is in
`variance_prior_review.md`.

Validation runs before MCMC and at archive load. Extracted child archives retain
their parent's outcome types, so extracting a binary outcome cannot bypass the
joint-prior check. Older child archives lacking those types require a proper
prior when their config indicates coupling. Changing old archive metadata cannot
repair the original draws; refit from the data. Python/R documentation and both
contract copies describe the requirement.

## Original-data comparisons

Artifacts: `/home/ignacio/longbet-bug-audit/continuation-20260908`.
The frozen input is copied there as `chapter-input.npz`, SHA-256
`6e4904482c10530728705f0c4c1ce8bc943db6c2f0aa6aedfba5bb99284aaffd`.
All runs use the original 750 sellers, 18 periods, target week 18/exposure 8,
59 benefit and 138 harm sellers, seed 314159, four chains, 2,000 warmup,
250 retained draws thinned by two, and explicit IG(2,1) innovation priors.
Each outcome has 20 prognostic and 20 total treatment trees.

| Run | Change | Worst R-hat | Minimum bulk ESS | Gates passed |
| --- | --- | ---: | ---: | ---: |
| `inherited-shared` | Inherited positive initialization; 10 shared trees | 1.8396 | 5.91 | 0/6 |
| `corrected-shared` | Restored prior initialization; 10 shared trees | 1.7839 | 6.10 | 0/6 |
| `corrected-private` | Restored initialization; no shared trees, SUR retained | 1.6673 | 6.65 | 0/6 |
| `restricted-control` | Fixed beta=1, b0=0, b1=1; joint forest/unit draws every 10 sweeps | 1.5936 | 7.06 | 0/6 |
| `stable-private` | Final numerical fit-cache correction; full GP, adaptive coding, SUR, no shared trees | 1.9002 | 5.75 | 0/6 |

The restricted control is a **different model**, retaining exposure/calendar
splits, random intercepts and full SUR feedback. It is not a repair of the
original posterior or a validated GP-free forecast API. All its collapsed
updates and capacity diagnostics are saved: all 250 updates per chain/outcome
were attempted and accepted, with at most 115 leaves against capacity 256.
Its failure cannot be attributed to silently skipping those updates for lack
of capacity. It remains a benchmark experiment.

The first three fits took 265.6, 281.6, and 172.4 seconds respectively, excluding
prediction/persistence; the restricted control took 259.9 seconds and the final
numerically corrected private model took 208.6 seconds.
CPU affinity and concurrent tests differed, so these are not fair performance
comparisons. Reports contain each run's settings and Python source fingerprint.

`summarize_multi_runs.py` independently uses ArviZ on the saved aligned draws.
It diagnoses business-scale and log effects, the final half of each chain,
seller event shares and events on group-mean effects. Tail ESS is evaluated at
the **2.5% and 97.5% endpoints** and mean MCSE uses mean ESS. The original
benchmark used bulk ESS for mean MCSE and 5%/95% tails; it is corrected too.
Constant indicators remain undefined and fail checks. Near-constant indicators
can pass numerical thresholds despite failed underlying effects; they do not
certify the joint decision. All unrestricted runs fail the benefit seller-share
diagnostics, while harm events are constant zero. `final-comparison.json` and
`final-traces.png` contain the independent audit of all five runs.

The trace figure shows persistent chain separation. Neither accurate
signs, point estimates, interval containment, nor discarding half the retained
draws rescues these fits. No repeated-data coverage, SBC, joint-event calibration,
or chapter publication gate has been passed in this continuation.

## Reproduction

From the development checkout, use an interpreter with the repository's pinned
dependencies (the host interpreter used here is
`/home/ignacio/longbet-jax/.venv/bin/python`) and `PYTHONPATH=src`:

```sh
python benchmarks/bench_multi_chapter.py --data chapter-input.npz \
  --output audit/private --mode multi --sigma-prior-a 2 --sigma-prior-b 1 \
  --burnin 2000 --draws 250 --skip 2 --chains 4 --seed 314159 --save-model

# Shared model: add --shared-trees 10.
# Restricted diagnostic: add --treatment-only --fixed-beta \
#   --collapsed-exposure --collapsed-scale 0 --collapsed-capacity 256.

python benchmarks/summarize_multi_runs.py --data chapter-input.npz \
  --output comparison.json --plot traces.png audit/private audit/shared
python -m pytest -q
```

The saved `model.npz` archives for all three corrected unrestricted runs load in a
fresh Python process, pass archive validation, and match all **19** separately
saved parameter arrays exactly. `archive-checks.json` records their hashes.
The R suite was run in the existing `longbet-rtest` image with this checkout
mounted read-only at `/src` and `PYTHONPATH=/src/src`; it passed **269** expectations,
including a full rerun after the final numerical correction (`stable-r.log`).

The full Python run initially completed with 314 passed, 11 failed, and one
optional external-BCF skip. Eight failures were fixtures still using the now
rejected prior; one expected the old sweep schedule, one compared noisy
forecasts to the wrong common level, and one exposed the numerical fit drift
described above. The fixtures and independent forecast oracle were corrected;
their targeted reruns passed. The final sampler regression run passed **62**
tests in 412.79 seconds (`stable-fits.xml`), including both Geweke checks,
independent mixed-outcome quadrature, scalar equivalence, missing observations,
shared forests, dispersed starts, and the new tiny-weight regressions. This
supports the implementation repairs, not convergence on the chapter panel.

The final interface/diagnostics/forecast run passed **54** tests in 51.02 seconds
(`final-interface.xml`). Combining the full run with the subsequent targeted
reruns gives **330 distinct passing tests, one optional skip, and no unresolved
test failures**. `validation-ledger.json` lists the latest result and source XML
for each test. This was one full suite followed by targeted reruns, not a second
full-suite execution on the final source. Previously passing tests outside the
changed paths were not repeated. `final-source-manifest.json` pins the source,
which matches the final original-data run; `final-changes.patch` includes the
tracked edits and new files.

No commit, push, image rebuild or chapter edit was performed. The old handoff's
`/home/ignacio/book/multi.md` no longer exists in the current book checkout;
this report and `bug.md` carry the updated investigation record.
