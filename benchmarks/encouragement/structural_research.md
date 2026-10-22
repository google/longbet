# Structural treatment-clock research assessment

The implemented encouragement model estimates reduced forms on the encouragement
clock. A treatment-clock structural forest or longitudinal principal-stratification
sampler is **not part of this release**, as required by milestone 6 of the
[implementation plan](../../docs/encouragement-implementation-plan.md). This document
specifies what would justify a separate research implementation.

## What randomization identifies

For each observed horizon, random encouragement identifies the assignment
contrasts in outcomes and adoption for the declared target population. Calling
their ratio a complier effect additionally needs exclusion, monotonicity, a
nonzero population first stage, and a treatment whose relevant history is captured
by the estimand. An instrument's existence alone is insufficient to identify an
arbitrary population average treatment effect. See [Angrist and Imbens's original
identification argument](https://www.nber.org/papers/t0118) and the motivating
[randomized encouragement chapter](https://book.martinez.fyi/iv.html).

For longitudinal outcomes, distinguish adoption stock `D_h(z)` from the adoption
history `D_1(z), ..., D_h(z)`. Exclusion through the entire history permits an
effect of earlier adoption even when both assignment arms have adopted by `h`.
For example, suppose everyone adopts at period 1 when encouraged and period 3
otherwise, and the outcome effect is the square root of time exposed. At period
3, the adoption-stock first stage is zero, but the outcome ITT is `sqrt(3)-1`.
The stock Wald ratio is undefined although encouragement has an outcome effect.
`test_acceleration_changes_history_after_stock_first_stage_disappears` encodes
this example in the calibration DGP.

With mixtures of accelerated adopters, never-takers and induced adopters, a
positive stock first stage does not eliminate this problem: accelerated adopters
can contribute to the numerator without contributing to its denominator. This
is an algebraic consequence of the potential histories, not evidence that the
instrument or a normality assumption has failed.

## Identification gates for a future model

Before adding a likelihood or sampler, write down the following:

1. **Target and intervention.** Specify a dynamic adoption intervention, the
   outcome horizon, and whether the target is all study units, current stock
   compliers, adoption-time compliers, or new units. Do not identify an observed
   individual's latent principal stratum from a single realized adoption path.
2. **Adoption histories.** Use a discrete-time hazard conditional on not yet
   adopted, or an explicit distribution of potential adoption times. Generate an
   absorbing path from that representation; independently sampled binary period
   responses do not enforce absorption. Define always-adopted and never-adopted
   states, censoring, and any monotonicity restriction on the two potential times.
3. **Exclusion.** Under an exclusion-restricted structural outcome model,
   encouragement columns must be absent from outcome forests except through the
   declared adoption history. A direct encouragement effect belongs in an explicit
   sensitivity model with an identified or externally constrained parameter.
4. **Sufficient variation or justified restrictions.** A single binary instrument
   supplies one assignment contrast per outcome horizon, not arbitrary variation
   in every treatment duration and every baseline subgroup. Identify which
   duration parameters the observed moments recover. A homogeneous additive
   distributed-lag response may be recoverable under rank and stability conditions
   on a matrix of lagged first stages; an unrestricted duration-by-person forest
   generally is not. State those restrictions before fitting.
5. **Latent dependence.** Correlated Gaussian adoption/outcome errors and random
   intercepts are distributional assumptions. They can change extrapolation and
   posterior regularization; they do not create randomization or identify an
   otherwise unrestricted response. Separate within-unit serial dependence,
   between-outcome innovations and correlated persistent unit effects.
6. **Extrapolation.** Effects for never-takers, always-takers or other populations
   require assumptions beyond the local IV argument. Specify these explicitly;
   do not label a prior-driven extrapolation a randomized comparison. The limits
   of extending local effects are discussed in [Angrist's treatment-heterogeneity
   paper](https://www.nber.org/papers/w9708).

Principal stratification does not remove these gates. It organizes missing
potential treatment responses into latent strata. The observed likelihood can
leave contrasts unidentified without additional assumptions or prior information;
the [clustered-encouragement principal-stratification study by Forastiere, Mealli
and VanderWeele](https://pmc.ncbi.nlm.nih.gov/articles/PMC5166715/) gives explicit
examples of required cross-world extrapolation and mechanism assumptions.

## Benchmark protocol required before a release

Begin with an analytically identifiable linear or finite-stratum special case,
derive every conditional, and verify prior-predictive marginal/conditional
(Geweke) identities. Prove observational equivalence or rank failure in at least
one intentionally unidentified case: a sampler that concentrates there without
declared identifying prior restrictions is a failure, not successful recovery.

Use the existing potential-history DGP as one source of misspecification cases,
and add correctly specified hazard and structural-outcome generators. Cross
sample sizes, weak and zero first stages, acceleration, delayed adoption,
heterogeneous duration responses, exclusion violations, defiers, latent
non-Gaussian dependence and prior strength. Independently seed populations,
assignment and chains. Preserve failed datasets and prior/posterior draws.

Evaluate the **declared** structural target, both ITTs, and the observable joint
distribution; report pointwise and simultaneous coverage separately. Compare
conditional study-unit standardization with a separately defined marginal target.
Diagnose all transformed effects and tail quantiles, retain longer-chain
sensitivity, and separate nonconvergence from persistent model error. Weak-IV
test inversion remains a reference for the assignment ratio; its validity does
not automatically transfer to a different structural duration target. See
[Aronow, Chang and Lopatto](https://arxiv.org/html/2404.18786v2) for the distinction
between randomization-based sets and normal approximations.

The benchmark-only `intercept_candidate.py` is a preliminary covariance
diagnostic with exact Gaussian conditionals. It uses linear reduced forms and a
working LPM, and does not enforce the binary adoption likelihood or identify
structural treatment effects. Its results therefore cannot justify silently
replacing LongBet's unit-intercept prior or claiming a general IV extension.
