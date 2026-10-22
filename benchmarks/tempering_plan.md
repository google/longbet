# Plan: parallel tempering for the LongBet forest sampler (2026-09-13)

## 1. What has to be fixed

On the chapter's onboarding panel (300 accounts, 8 weeks, one launch cohort)
four dispersed chains of the reduced-form LongBet fit disagree on the
population and unit-level intent-to-treat effects even with the exact CHANGE
and REGROW tree moves on (bulk ESS 58, R-hat 1.06 at 4 x 1,000 draws for the
week-8 population effect; the joint probit model is worse, ESS 6-16, R-hat
1.2-1.7). The per-chain effect surfaces show why: every chain finds the sharp
readiness/scale thresholds, but each represents the surface inside the
qualifying region with a different step function spread across the 20
treatment trees, and one chain sat for 4,000 sweeps with the scale threshold at
0.45 instead of -0.2. A single regrown tree cannot leave such a state because
the other trees, the unit intercepts and the exposure coefficients jointly
encode the current decomposition, so one-tree proposals are rejected. The
package's audits (`slow_mode_review.md` and follow-ups) reached the same dead
end with collapsed, whole-tree and joint-Gaussian blocks on a different panel.

Removing pieces of the model does not help: no GP, fixed or treatment-only
coding, no unit intercepts, no calendar-time treatment splits, larger minimum
leaf sizes and few deep trees all mixed the same or worse.

## 2. Options considered

**Whole-forest independence proposals (sequential data-driven forest, exact
MH).** Rejected as the first line. The acceptance ratio contains the reverse
density of the *current* forest under the sequential proposal. A converged
forest's trees are individually "unnatural" for the residual they would see in
a sequential rebuild (they stop where the rebuilt residual says to split), so
that density is astronomically small and the proposal never accepts. The
block-regrow experiment already showed this at two trees; conditional SMC with
a reference particle has the same failure in the form of an exploding
reference weight. Tempering the proposal helps only by making it prior-like,
which destroys its quality.

**Parallel tempering (PT).** Chosen. It is exact for any bridging family, it
reuses the whole existing sweep (all moves stay valid at every level), and it
is the only approach with a guarantee of eventually mixing over forest modes.
Its cost is a ladder of K tempered replicas per cold chain. A `conditional
precision` path already exists in the sweep (used by the SUR coupling): every
mean, latent and coefficient update takes a per-cell likelihood precision, so
a temperature beta is implemented exactly by passing `beta / sigma2` and then
drawing sigma2 from the tempered inverse-gamma conditional.

**Annealed burn-in.** Free by-product of the same machinery: ramp beta from
beta_min to 1 during burn-in on each chain, then sample at beta = 1. Discarded
burn-in draws leave the stationary target untouched. This finds the dominant
mode reliably but does not sample between modes; it is reported separately so
the two effects are not confused.

## 3. Ladder size from the measured heat capacity

Stage 0 measures the posterior standard deviation s of the log-likelihood at
beta = 1 (`runs/diag_ll.py`). With Var_beta(log L) ~ s^2 / beta, adjacent levels
exchange at a useful rate when their sqrt(beta) differ by about 1 / (2 s),
so a ladder equally spaced in sqrt(beta) needs about
K = 2 s (1 - sqrt(beta_min)) levels; geometric spacing was measured to
accept only 1-3 percent near beta = 1. With s ~ 14 (about 400 well-identified
coefficients: 300 unit intercepts plus leaves) and beta_min = 0.05 this is
K ~ 22. The hot end at beta = 0.05 turns a 40-unit split gain into 2 units, so
hot replicas are close to prior draws and carry no mode memory.

## 4. Implementation

1. `LongBetConfig`: `tempering_levels: int = 1`, `tempering_beta_min: float =
   0.05`, `anneal_burnin: bool = False`. Contract entries, R arguments, docs.
2. `LongBetState`: per-chain `temperature` field (chain axis), static
   `tempering_levels`. `init_longbet` builds `num_chains * levels` replicas;
   replica `c` has level `c % levels`, geometric ladder from 1 to beta_min.
3. `_step.py`: `longbet_step` passes `conditional_precision = temperature /
   sigma2` (masked) whenever the state is tempered or annealed, then draws
   sigma2 from IG(a + beta m / 2, b + beta SSR / 2) for continuous outcomes.
   The binary latent uses error scale 1 / sqrt(beta), which is the augmented
   tempered target N(z | f, 1/beta) 1{sign}.
4. `_loop.py`: after each sweep, a swap step per ladder: energy E = log L of
   the untempered likelihood (continuous: -(m/2) log sigma2 - SSR/(2 sigma2);
   binary/ordinal: -SSR/2), even-odd alternation of adjacent pairs, acceptance
   exp((beta_i - beta_j)(E_j - E_i)), and a gather of every chain-axis field
   except `temperature` along the accepted permutation. Only the cold replicas
   are written to the traces, so memory does not scale with K. Swap acceptance
   counts are returned for diagnostics.
5. Annealed burn-in: `temperature(i) = beta_min^(1 - i / n_ramp)` for
   `i < n_ramp = 0.6 * num_burnin`, then 1; no swaps.
6. Tests: (a) tempered sweep at beta = 1 is bit-for-bit the existing sweep;
   (b) with `tempering_levels = 3` on a small panel, the cold chains reproduce
   the plain sampler's posterior (Geweke-style moment comparison on the same
   data); (c) swap step preserves the residual identity and the per-replica
   temperatures; (d) annealed burn-in equals the plain sampler after burn-in.
7. Joint model (`_multi_step.py`): same precision scaling on the SUR
   conditional precision, tempered structural variance draws, swaps of the
   whole multi-state; energy = full SUR log-likelihood. Done only if the
   single-outcome result is positive.

## 5. Test protocol on the chapter panel

Same pilot (`generate_panel_data(300, seed = 2026)`), 4 cold chains, seed 7,
3,000 burn-in, 1,000 draws, moves on. Report for each arm: per-week population
ITT bulk ESS and R-hat, per-chain week-8 means, test RMSE on the 500-account
cohort, swap acceptance by level, and round trips per ladder.

Arms: (A) plain sampler with moves (baseline, done: ESS 58, R-hat 1.06);
(B) annealed burn-in only; (C) PT with K = 12 and K = 22; (D) B and C at 4x
the sweeps to measure ESS scaling. Pass criterion is the package gate (bulk
and tail ESS >= 400, R-hat <= 1.01). A partial result is reported as such.

## 6. Results (2026-09-14)

Implementation: `src/longbet/_tempering.py` (ladder, energies, exchange step,
replica selection, coupled-sampler variants), `_step.py`
(`_tempered_single_step`), `_loop.py` / `_multi_loop.py` (annealing schedule,
exchanges, posterior-replica traces), `_state.py` / `_multi_state.py`
(per-replica temperature, ladder assignment), `_multi_step.py` (tempered SUR
precision, loadings and innovation variances), config, contract, R front doors,
docs. Tests: `tests/test_tempering.py` (layout, exchange bookkeeping, joint
model, and a slow moment comparison of cold replicas against the plain
sampler, all passing); full non-slow suite 870 passed; targeted R suite 385
passed. Untempered runs keep their random stream bit for bit.

Calibration on the chapter panel (2 cold chains, 16 geometric levels, 500
sweeps): exchange acceptance 1-3 percent at the cold end and 39 percent at the
hot end, so the ladder was changed to equal steps in sqrt(beta).

Full run, reduced-form outcome model, 4 cold chains x 26 levels, beta_min =
0.05, 3,000 burn-in, 1,000 draws, moves on, seed 7 (fit 7,529 s on 12 CPU
cores, against 220 s for the plain sampler):

| arm | week-8 ITT chain means | ESS by week (4-8) | R-hat by week | log-lik by chain | RMSE |
| --- | --- | --- | --- | --- | ---: |
| plain sampler with moves | 2.82 2.66 2.71 2.71 | 58 at week 8 | 1.064 at week 8 | -4874 -4835 -4826 -4882 | 0.81 |
| annealed burn-in | 2.67 2.66 2.67 2.73 | 10 10 30 71 133 | 1.33 ... 1.03 | -4927 -4837 -4875 -4882 | -- |
| parallel tempering | 2.63 2.69 2.71 2.68 | 100 35 72 81 125 | 1.040 1.081 1.048 1.042 1.032 | -4813 -4830 -4842 -4862 | 0.89 |
| default sampler (frozen) | 2.67 2.64 2.74 2.81 | 38 at week 8 | 1.078 at week 8 | -4974 -4979 -4971 -4996 | 0.53 |

Exchange acceptance was 0.21-0.31 on the 15 coldest gaps and 0.10-0.15 on the
hottest. Every cold chain ends in the high-likelihood region (no chain 50-100
units below the others, as happens without tempering), and all four agree on
the shape of the effect surface, including the scale threshold that one
untempered chain misplaced for 4,000 sweeps. The remaining disagreement is
between forest modes of comparable likelihood that exchange states only every
few hundred sweeps: R-hat 1.03-1.08 and bulk ESS 35-125 against the gates of
1.01 and 400. ESS should grow linearly with sweeps under this kernel, so the
gates need roughly 4-10 times the budget: on this CPU 8-20 hours per fit, and
the joint model costs about twice that. Annealed burn-in alone does not remove
the disagreement: chains settle in different comparable-likelihood modes.

The accuracy ordering is the reverse of the likelihood ordering. The frozen
default chains sit 100-150 log-likelihood units below the tempered chains and
predict the held-out cohort best (RMSE 0.53 against 0.81-0.89). The exact
posterior of this model, sampled properly, resolves finer heterogeneity in the
treatment forest than the panel supports; that is a prior/model question, not
a sampler question, and it caps what better mixing can deliver on this panel.

## 7. What changed in the package (2026-09-14)

The finding in section 6 (the exact posterior of the reduced-form encouragement
model is less accurate than early-stopped chains) led to replacing that model
rather than its sampler: `LongBetEncourage` now fits the outcome on the observed
adoption clock and adoption as a discrete-time hazard on the offer clock, and
composes the two (frailty integrated outside the survival product). On the same
panel and budget as the plain runs above (4 chains, 2,000 + 1,000 sweeps, no
tempering, about six minutes) this gives held-out offer-effect RMSE 0.36-0.39,
bulk ESS 143-307 and R-hat 1.01-1.04 on the outcome side, and ESS 272-489 with
R-hat at most 1.02 on the adoption side. The CHANGE and REGROW moves run in
every sweep, the coding scales are fixed at 0 and 1, the treatment forest never
splits on the exposure clock, innovation priors are proper, and the retired
engines and experimental blocks were removed. Tempering remains available.
