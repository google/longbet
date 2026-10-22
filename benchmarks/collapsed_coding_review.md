# Collapsed coding and exposure experiments (fourth pass, 2026-09-07)

Status: investigation completed without a reliable repair. All seven four-chain
fits and final correctness/checkpoint checks are complete. This agent is
stopping after both longer replicated original-posterior trials failed every
subgroup gate and ESS failed to scale. This is not a proof that no repair is
possible. Read the fourth-pass section of `bug.md` for reproduction and the
handoff. All output is under `/home/ignacio/longbet-bug-audit/coding-pass`.

## Rationale and saved-state profiles

The third-pass GP proposal held b0/b1 fixed. Its joint leaf redraw did not fix
the conditional fixed-tree model, and raising GP acceptance did not fix the
original posterior. The older joint-prognostic block moves coding with one
tree at a time. This pass integrates all forests before updating coding.

`profile_collapsed_coding.py` loads each saved third-pass checkpoint, verifies
the new coding integral against the direct full-panel integral, and profiles
257 angles at the current coding radius. Both seed checkpoints were processed,
including all four chains and all outcomes. Maximum direct/cached integral
differences are about 1e-9. The GMV angle profiles have only 4-34 grid points
within two log units of their maximum. Current angles are within 0-1.51 log
units of those grid maxima, so this is not evidence of a large missing local
mode. These are conditional profiles, not posterior estimates or exact grid
updates. Different chains condition on different remaining parameters.

## Exact collapsed updates

The integrated likelihood and immediate all-leaf/unit redraw are inherited from
`_collapsed_exposure.py`; see `collapsed_exposure_review.md` for the underlying
Schur complement. The benchmark now accepts `--collapsed-proposal`:

* `pcn`: the original GP prior-reversible proposal.
* `coding`: elliptical slice update of (b0,b1), conditional on beta.
* `joint`: elliptical slice update of beta and coding together, respecting fixed
  settings and using their independent Gaussian prior blocks.

Elliptical slice sampling uses the original zero-mean Gaussian prior as its
reference. An independent prior draw supplies the ellipse direction. A uniform
slice threshold is drawn below the current collapsed likelihood; rejected
angles shrink a bracket containing angle zero. There is no iteration cap or
unaccepted fallback. After acceptance all marginalized coefficients are drawn
immediately. Method reference: https://arxiv.org/abs/1001.0175. This reference
does not establish that the LongBet adaptation mixes adequately.

Coding-only likelihood evaluations use a cached quadratic. Form H with the
unscaled leaf columns, response column, and treatment-offset column. Partition
rows by treatment. With C_z=U'W_z H and D=I+U'WU, cache

```
G00 = H'W0 H - C0'D^-1 C0
G11 = H'W1 H - C1'D^-1 C1
G01 = -C0'D^-1 C1
```

Each stratum has a column multiplier d_z: coding b_z on treatment leaf columns,
one on prognostic/response columns, and -b_z on the treatment-offset column.
Their scaled Gram sum gives the same leaf Schur system and integrated residual
quadratic as revisiting the full observation panel. Cross-stratum unit terms
must be included. The coefficient count/overflow policy and original priors are
unchanged. Diagnostics append likelihood evaluation count as a fifth field.

## Coding/GP interweaving

An additional experimental transition is implemented in `_coding_interweave.py`.
For the exposure convention z=0 iff S=0, define
`a=(b0*beta[0], b1*beta[1:])`. Holding a fixed preserves the likelihood and
within-grid potential-outcome contrasts. The conditional prior on coding in
these coordinates includes the determinant `|b0|^-1 |b1|^-(d-1)`.

For log magnitudes u=log|b| and group sizes n=(1,d-1), the target is

```
-.5 * [sum(exp(2u))/sigma_b^2 + ||L^-1(a / b[group])||^2]
+ sum((1-n)*u)
```

The final term includes BOTH the beta-to-a and coding-to-u Jacobians. Four
coding sign configurations are sampled from their exact conditional weights;
each log magnitude and their common shift are sliced using a Gaussian reference whose
log density is subtracted from the target. State restoration accounts for
float32 product roundoff in raw residuals. Unsupported exposure conventions or
fixed coding/beta disable the move. Benchmark flag: `--interweave-coding`.

## Validation and benchmark protocol

The first expanded test run passed 14 checks in 129.34 seconds. Three observed
probit/Gaussian quadrature cases passed in 42.32 seconds, two of them new, for
16 distinct tests at that point. These include sampled coding AND beta with
shared/private forests, full SUR dependence, posterior moments and a joint
event. They establish small-model correctness, not chapter convergence.

Fixed-tree coding and joint trials use the third-pass `instrumented-short`
common topology and fresh dispersed starts: four chains, 2000 burn-in, 250
retained, skip 2, 20 prognostic trees, 20 treatment trees including 10 shared,
joint-prognostic enabled, collapsed update every 10 sweeps, capacity 128. The
data, estimand, subgroups, priors and public defaults are unchanged. Runs and
small tests overlap on this CPU host; reported times are contended runtimes.

All original-posterior convergence/recovery, ESS-scaling, calibration and
integration gates in `bug.md` remain required. Conditional fixed-topology
success alone would not establish a repair.

## Root changes and ALL-forest topology integration

`_change_prognostic.py` implements three equally weighted component kernels for
each prognostic tree: replace its root rule from the root rule prior, move the
cutpoint by one in either direction, or swap the root rule with an internal
descendant. Tree shape stays fixed. The prior includes the rule probabilities
at every affected descendant, their ancestor-restricted cutpoint ranges, and
terminal probabilities that become one when no rule remains available. Count
constraints preserve the support reachable by the existing grow/prune steps.
All leaves of that tree are integrated for acceptance and drawn immediately.
Membership, pending-prune markers, count/precision caches and affluence are
updated together. Oversized trees skip this optional move. The benchmark flag
`--change-prognostic` runs the block every ten sweeps.

The additional `--collapsed-tree` option uses a stronger acceptance calculation:
all forests and unit effects for the outcome are integrated. It chooses one
prognostic tree uniformly per outcome and invokes the entire root MH kernel
with zero likelihood. This proposal Q is reversible with respect to the actual
tree prior. Consequently, for the collapsed target pi(T)*L(T),

```
pi(T') * L(T') * Q(T',T) / [pi(T) * L(T) * Q(T,T')] = L(T') / L(T).
```

This is a two-stage proposal, not an assumption that raw root rules have a
symmetric proposal. All integrated coefficients are redrawn after either
second-stage decision. Unlike one-tree leaf integration, this lets other
prognostic and treatment trees compensate for a changed partition. Shared
treatment partitions remain fixed during this outcome-specific move; their
leaf coefficients are included in the integral and redraw.

Correctness checks include independently enumerated root posteriors under
uniform and unequal variable priors, fitted means/full covariance, actual
binary-root probabilities from probit quadrature, and a two-outcome joint root
posterior from an independent dense SUR Gaussian covariance plus numerical
orthant probabilities. The latter checks all nine joint root configurations
with private/shared treatment coefficients integrated and observed binary data.

## Completed comparisons and stopping rationale

Seven full comparisons have completed; all five unrestricted trials failed
every subgroup gate. Times are contended wall times, not efficiency claims.

| Run | Target | Worst R-hat | Min bulk ESS | Gates |
| --- | --- | ---: | ---: | ---: |
| fixed-coding-314159 | Fixed topology | 1.0107 | 197.0 | 3/6 conditional |
| fixed-joint-314159 | Fixed topology | 1.0973 | 30.1 | 2/6 conditional |
| coding-314159 | Original | 1.6445 | 6.6 | 0/6 |
| coding-271828 | Original | 1.3534 | 9.4 | 0/6 |
| combined-314159 | Original, coding + interweaving + root changes | 1.3006 | 10.7 | 0/6 |
| ensemble-314159 | Combined + ALL-forest root proposals | 1.4007 | 8.5 | 0/6 |
| ensemble-271828 | Same final candidate, replicated seed | 1.4983 | 7.5 | 0/6 |

The fixed-coding GMV benefit has R-hat 1.0107, bulk ESS 198.6, tail ESS 547.8,
versus the old fixed-tree control's 1.2150 / 14.4 / 54.1. This conditional
improvement does not transfer to the original posterior. The combined run
reduces the number of GMV prognostic roots that never change during retained
sampling to 10-15 of 20 per chain, but still fails inference. Root movement
alone is not a convergence certificate. Every mean, interval, coverage failure,
RMSE, sign result and chain mean remains in the complete reports.

The final `ensemble-314159` and `ensemble-271828` trials add `--collapsed-tree`
and use 2000 burn-in / 500 retained / skip 4: 4000 sweeps per chain and a
fourfold longer retained window. Both completed, taking 1408.1 and 1622.0 seconds
including compilation and contention; the latter includes a pause to stagger
prediction memory. No fit or monitor remains running or paused. The completed
prefix comparison uses 125, 250 and 500 draws from each chain, plus the last 250
draws. Prefixes are dependent and cannot replace independent model seeds.

| Seed / summary | Bulk ESS at 125 | At 250 | At 500 | Growth for 4x window |
| --- | ---: | ---: | ---: | ---: |
| 314159 / GMV benefit | 24.36 | 17.12 | 20.14 | 0.83x |
| 314159 / GMV harm | 25.78 | 17.45 | 11.86 | 0.46x |
| 271828 / GMV benefit | 7.36 | 6.62 | 8.59 | 1.17x |
| 271828 / GMV harm | 9.55 | 8.00 | 7.46 | 0.78x |

Hours benefit stays around bulk ESS 8 in both trials. Dropping the first half
does not remove the failures. The final runs reduce unchanged GMV prognostic
roots to 3-10 and 5-11 of 20 per chain, respectively, without making inference
reliable. No active-leaf packing skips occurred in any fourth-pass full run.
The first final run's GMV harm interval misses truth; all estimates and coverage
results remain archived and must not be presented as validated inference.

These results justify stopping this investigation without claiming a repair.
Full subtree/Particle Gibbs proposals, broader treatment-partition integration
and other global mode-crossing methods remain untested possibilities, not
solutions established here. Public defaults remain unchanged. All public API/R
integration, generalization/calibration and chapter gates remain unmet.

The full-data 30-sweep combined and ensemble smoke tests completed in 57.4 and
57.6 seconds, respectively, with all effect arrays finite. These are execution
checks only. `coding-pass/compare.py` verifies complete archives and emits the
full comparison and prefix diagnostics. `bug.md` records detailed test runs,
including one original frequency assertion failure and its subsequent MCSE
and independent one-step detailed-balance checks; do not erase that history.

The final consolidated regression passed **95 tests in 504.51 seconds**. It
includes the existing regression floor plus independent dense Gaussian,
observed probit, prior-invariance, enumerated topology and joint SUR/root
posterior checks. Passing these does not establish convergence on the original
panel. The full Python/R suites were not rerun for the failed candidate.

All seven archives are complete, finite and correctly shaped. Both final
ensemble checkpoints were independently reconstructed from actual tree rules:
the maximum fitted-component discrepancy is below 4e-7 and maximum raw-residual
discrepancies are 6.448e-5 / 6.489e-5 after 4000 float32 sweeps, below the 3e-4
check tolerance. Every chain retains binary marginal variance exactly one.
The per-run `checkpoint-validation.json` files retain every component error.
Both final trace figures are saved; the second was visually inspected and
shows the same persistent chain separation as the diagnostics.
