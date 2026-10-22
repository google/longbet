# Original chapter posterior is improper under the default variance prior

Status, 2026-09-07: the normalization defect below is proved for the original
chapter configuration and independently checked against the current public
engine's support. This does **not** establish that the boundary causes the
observed ESS failures, or that proper priors repair mixing. Both explicit
`IG(2,1)` variance-prior benchmarks fail all six convergence gates.
No public default or inference
recommendation has been promoted. See [`bug.md`](../bug.md) for subsequent
benchmark status and all remaining validation gates.

The intended continuous-parameter chapter posterior with
`sigma_prior_a=sigma_prior_b=0` has no finite normalizing constant. The proof
constructs an allowed positive-prior topology and finite means for which the
**observed** probit/Gaussian likelihood has a positive limit as the last
continuous innovation variance goes to zero. The construction uses the actual
frozen design and observed responses. Topology construction uses covariates
only; neither it nor the feasibility check uses true treatment effects.

The audit is preserved under
`/home/ignacio/longbet-bug-audit/refresh-pass/propriety*`. This review contains
the argument independently of those external files; the files provide its
reproducible computational certificate.

## The variance prior and observed SUR boundary

The public defaults in [`_config.py`](../src/longbet/_config.py) are both zero.
The variance conditionals in [`_step.py`](../src/longbet/_step.py) and
[`_multi_step.py`](../src/longbet/_multi_step.py) use the inverse-gamma convention

```
p(v) proportional to v^(-a-1) exp(-b/v),  v=sigma2>0.
```

Thus the default is `dv/v`. Leaf, GP, coding, unit-effect and loading priors
do not shrink with this innovation variance. For example, the loading
conditional adds `I/sur_prior_var` to `X'X/v`; its Gaussian prior variance
is independent of `v`. Binary innovation variance is fixed at one.

The actual internal outcome order is complaint, GMV, hours, indexed 0, 1, 2.
Write `r_m=y_m-f_m`, with augmented binary response `y_0=z`. The triangular
model implemented in [`_sur.py`](../src/longbet/_sur.py) is

```
r0 = e0,                     e0 ~ N(0,1)
r1 = G10*r0 + e1,            e1 ~ N(0,v1)
r2 = G20*r0 + G21*r1 + e2,   e2 ~ N(0,v2).
```

Hold `G20!=0` and `v1>0`. At `v2=0`, the map from `(e0,e1)` to the observed
continuous residuals `(r1,r2)` has determinant `-G20`, so remains nonsingular.
Put `r0*=(r2-G21*r1)/G20`. The limiting observed-data density per cell is

```
phi(r0*) * phi_v1(r1-G10*r0*) / abs(G20)
  * 1{ (2*b-1)*(f0+r0*) > 0 },
```

where `b` is the observed binary response. It is strictly positive when the
binary inequality is strict. This integrates the binary latent value; holding
that value fixed in a Gaussian surrogate would miss the mechanism.

## A real design constraint rules out the simpler proposed proof

There are 750 sellers and 18 weeks, giving 13,500 cells. The binned covariates
have only 749 distinct seller profiles. Sellers 339 and 726, using zero-based
indices, have identical six-dimensional covariate bins. Seller 339 is never
treated; seller 726 starts exposure one at internal week 6. Their complete
predictors, coding and exposure weights therefore coincide in weeks 0 through 5.

All forests give identical pair contributions during those six weeks,
regardless of topology or coefficient values. Unit effects supply only a
constant pair difference. Hence these five independent contrasts annihilate
**every** mean design:

```
(seller339, week t) - (seller726, week t)
 - (seller339, week0) + (seller726, week0),  t=1,...,5.
```

Mean-design rank is at most 13,495. A fixed-latent proof assuming full row rank
therefore does not apply to this chapter. Counting potential coefficients or
checking structural rank would not suffice. The stronger observed-SUR proof
only requires feasible binary inequalities, which the next construction supplies.

## Allowed topology and exact combined-rank certificate

Choose `G10=0`, `G20=1`, `G21=1/100`, and define, with this sign convention,

```
q = G20*f_complaint + G21*f_gmv - f_hours.
```

The nonzero `G21` is necessary for the GMV forests to contribute. Set every
GP component and each coding coefficient to one, with the original fixed
prognostic multiplier one. Gaussian priors have positive density in a
neighborhood of this point. All resulting forest columns are ordinary leaf
indicators multiplied by nonzero constants.

The combined budget is exactly that of the original shared-forest model:
60 prognostic partitions (20 per outcome), 30 private treatment partitions
(10 per outcome), and 10 shared treatment partitions. The shared partitions
are counted once; their nonsingular vector-leaf Gaussian prior lets each
shared leaf have an arbitrary combined `q` coefficient.

The constructor uses 99 of these trees for profile functions. Each separates
calendar weeks and then partitions sellers by their binned covariates. All
terminal depths are at most 9, corresponding to the original `max_depth=10`
storage convention. Every leaf contains at least 10 cells, and all rules lie
inside their ancestor-restricted cut ranges. A final covariate-only repair
replaces one redundant week-3 subtree with a binary-covariate split.

For **each** week, the 750-row profile-indicator matrix has exact rank 749
over GF(2). A nonzero GF(2) minor has an odd integer determinant, giving real
rank at least 749. The duplicate profile bounds it above by 749. Since every
leaf is confined to one week, the direct-sum real rank is `18*749=13,482`.
This certificate does not use a floating-point singular-value threshold.

The remaining private hours-treatment tree separates weeks, then splits at
exposure one in weeks 6 through 17. It has 30 leaves, maximum depth 6 and
minimum leaf count 291. It adds one independent pair contrast in each of
those twelve weeks. A seller-339 unit-effect column adds one remaining
direction for the constant pair difference in the six early weeks. Thus

```
combined real rank = 13,482 + 12 + 1 = 13,495.
```

The five displayed contrasts are the entire left null space. Every constructed
topology has positive original prior mass: allowed split probabilities and
rule probabilities are positive, terminal probabilities are positive or
forced to one, and count constraints hold. The shared-tree global rule prior
also gives positive mass to these geometrically valid partitions.

The main verifier reconstructs memberships from actual heap rules, validates
support, and recomputes the ranks. An independent audit reconstructed the
trees and recomputed GF(2) rank with the opposite pivot ordering. It also
verified that the exported design matches a fresh checkpoint load.

## Strict binary feasibility and finite coefficients

Let the actual internal standardized responses define
`k=y_hours-(1/100)*y_gmv`. At the boundary, `z*=k+q`. Outside the six early
duplicate-pair weeks choose `q=(2*b-1)-k`, giving the desired latent sign
with magnitude one.

For the early pair let `d_t=k_339,t-k_726,t`. Their observed bits are

```
seller339: 0 1 0 0 0 0
seller726: 0 0 0 0 0 1.
```

The stored float32 responses are exact dyadic rationals. At `G21=1/100`,

```
d_1-d_5 = 135108327 / 1677721600 > 0.
```

Choose the constant pair difference in `q` to be `u=-(d_1+d_5)/2`. At weeks
1 and 5 center the two latent values around zero. Their signs have the required
opposite ordering, each with magnitude

```
(d_1-d_5)/4 = 135108327 / 6710886400
            = approximately 0.0201327095925808.
```

At the other four early weeks shift both latents to at most minus one, preserving
the same pair difference in `q`. All 13,500 binary inequalities are therefore
strict. The five null constraints hold exactly by construction; the exported
floating check has maximum roundoff `5.6e-17`.

The exact rank and null-space certificate imply that this `q` has a finite
combined leaf/unit coefficient realization. Fixed forest offsets only add
a fixed combined constant here; the profile space contains that constant,
so subtracting it causes no obstruction. No approximate least-squares
coefficient fit is needed to establish existence. Proper Gaussian leaf and
unit priors give positive density near every finite realization, conditional
on positive unit variances.

## Divergence and scope of the conclusion

Strict margins persist on an open neighborhood of the finite coefficients,
coding, GP, loadings and remaining positive variances. This neighborhood has
positive prior measure and can be chosen independently of `v_hours`. Select
a compact sub-neighborhood away from zero in `abs(G20)` and `v_gmv` and away
from all binary sign boundaries. The observed likelihood has a positive
continuous limit as `v_hours` tends to zero, uniformly on that compact set.
It is therefore bounded below by a positive constant for sufficiently small
`v_hours`. Its independent default prior contributes

```
integral from 0 to epsilon of dv_hours/v_hours = infinity.
```

The topology prior mass can be extremely small but remains positive; multiplying
a divergent integral by that constant does not make it finite. The original
chapter posterior is consequently improper. Proper individual full conditionals
and fixed-variance posterior tests cannot establish a missing joint normalizer.

This is a statement about the intended continuous statistical target. Floating
storage, underflow or numerical guards do not define a documented proper
replacement prior. The sampler's gamma-variate floor prevents division by zero;
it is not an explicit positive lower bound on innovation variance.

The proof does **not** show that the chains approach this boundary, or that it
causes their separated effect traces. The two original-prior whole-tree trials
retained hours innovation variances from 0.2497 to 0.2741 and GMV variances
from 0.07179 to 0.08100. Their
[archived ranges](/home/ignacio/longbet-bug-audit/refresh-pass/original-variance-ranges.json)
describe failed chains under the improper target; they are neither certified
posterior intervals nor evidence against the divergent integral at zero.

## Why positive inverse-gamma hyperparameters give a proper chapter posterior

For this complete-observation chapter, independent `IG(a_m,b_m)` innovation
priors with `a_m>0`, `b_m>0` for every continuous outcome, together with the
remaining proper model priors, suffice for a finite positive observed-data
normalizer. This is stronger than only removing the constructed boundary.

Condition on all means, loadings and other parameters. Let `C` index the
continuous outcomes and `B` the binary outcomes, with unit-variance independent
binary innovations. Their marginal continuous residual covariance is

```
Sigma_CC = B_CC^-1 (D_C + G_CB G_CB') B_CC^-T,
D_C = diag(v_m : m in C),     det(B_CC)=1.
```

Here `B_CC` is the unit lower-triangular continuous block of `I-Gamma`.
Since `G_CB G_CB'` is positive semidefinite,
`det(Sigma_CC) >= det(D_C) = product_m v_m`. The observed binary-event
probability conditional on the continuous responses is at most one.
Consequently, for `n=13,500` complete cells and `M_C=2`, the joint observed
likelihood satisfies the uniform bound

```
L(binary, y_continuous | parameters)
 <= (2*pi)^(-n*M_C/2) * product_m v_m^(-n/2).
```

This bound is independent of all means, tree topologies and loadings. The
inverse moments required to integrate it are finite:

```
E_IG(a,b)[v^(-p)] = b^(-p) * Gamma(a+p) / Gamma(a),  p=n/2.
```

All other proper priors integrate to one, establishing a finite normalizer.
The likelihood is positive for positive innovation variances on positive-prior
sets, so the normalizer is also positive. This argument integrates the actual
observed binary events and allows arbitrary permitted continuous loadings;
it does not require a covariance bound estimated from sampled traces.

The `IG(2,1)` experiment is a disclosed model change using existing
public `LongBetConfig(sigma_prior_a=2, sigma_prior_b=1)` fields on the same
standardized response scales. In this convention it has mean one and infinite
variance. Posterior propriety provides no guarantee of efficient sampling,
recovery or calibrated decisions. Those checks remain required, and no defaults
are changed by this review.

The two complete proper-prior fits (seeds 314159/271828, four chains, 2000
burn-in, 250 retained, skip 2) each pass 0/6 subgroup gates. Worst R-hat is
1.6885/1.3311, minimum bulk ESS 6.4/9.7, minimum tail ESS 19.9/47.8.
Independent checkpoint reconstruction passes, with maximum residual error
4.94e-5/5.81e-5. These are empirical failures of the tested inference procedure
despite a now proper posterior. Final comparisons and validation limits are
in `bug.md`; the pass is complete without a reliable inference repair.

## Separate public-model edge case

An elementary edge case is a one-time-point continuous panel with
random unit intercepts. Then `U=I_N`; integrating their proper Gaussian prior
adds `tau2*I_N` to the observation covariance. For fixed `tau2>0`, the likelihood
stays positive as innovation variance tends to zero, so the default `dv/v`
prior also makes that model improper. This case does not rely on the chapter
tree budget, SUR or binary augmentation.

## Checks against current engine support

No conflicting support assumption was found. The audit checks the following
against initialization and the current source, without relying on experimental
kernel capacity or on the sampler reaching the constructed topology:

* [`_state.py`](../src/longbet/_state.py) initializes the exact tree budgets,
  positive leaf priors independent of innovation variance, Gaussian coding/GP
  support, fixed alpha one and random unit effects with proper variance priors.
* [`_shared_forest.py`](../src/longbet/_shared_forest.py) allows the witness
  shared partitions and gives their vector leaves a positive diagonal prior.
* [`_multi_model.py`](../src/longbet/_multi_model.py) supplies the same design
  and allowed split grids to every outcome. All chapter responses are observed,
  so the proof needs no extension to unsupported missingness.
* [`_multi_step.py`](../src/longbet/_multi_step.py) permits the nonzero direct
  hours-from-binary loading and the nonzero hours-from-GMV loading, with an
  innovation-variance-independent Gaussian prior. Binary identification imposes
  no constraint removing either loading.
* `max_depth=10` permits decision depths 0 through 8 and leaves at depth 9.
  The min-leaf requirement is 10 and no extra decision-node threshold applies.
  Every witness tree can be built by allowed grow moves because its final
  supported leaves make every intermediate ancestor subtree supported too.
* The experimental 128-leaf packing cap is an identity/overflow policy for an
  optional transition. It does not truncate the public model's topology prior
  or change the witness's positive support under ordinary tree updates.

## Reproduction, artifacts and hashes

Run these from `/home/ignacio/longbet-jax`; they perform no MCMC fit:

```sh
taskset -c 8-11 .venv/bin/python /home/ignacio/longbet-bug-audit/refresh-pass/propriety_export.py
taskset -c 8-11 .venv/bin/python /home/ignacio/longbet-bug-audit/refresh-pass/propriety_rank_witness.py
taskset -c 8-11 .venv/bin/python /home/ignacio/longbet-bug-audit/refresh-pass/propriety_verify.py
```

The [exporter](/home/ignacio/longbet-bug-audit/refresh-pass/propriety_export.py)
loads `refresh-pass/smoke/final_state.eqx` with
[`load_multi_benchmark_state.py`](load_multi_benchmark_state.py). It checks the
frozen input hash and asserts the actual budgets, priors and support settings.
Its exported design/responses are byte-identical to the independently reviewed
inputs. The [constructor](/home/ignacio/longbet-bug-audit/refresh-pass/propriety_rank_witness.py)
uses fixed covariate-only construction seed `975318642`. The
[verifier](/home/ignacio/longbet-bug-audit/refresh-pass/propriety_verify.py)
checks actual rules, exact rank, the additional exposure/unit directions and
the strict observed binary target.

The [certificate](/home/ignacio/longbet-bug-audit/refresh-pass/propriety-witness.json)
records support counts, exact rank, the rational margin and artifact hashes.
The [independent review](/home/ignacio/longbet-bug-audit/refresh-pass/propriety-independent-review.json)
records separate support/rank/sign checks. The NPZ files preserve all heap
rules, memberships, target means and annihilator contrasts. Logs are
`propriety-export.log`, `propriety-rank.log` and `propriety-verify.log` in the
same audit directory.

| Artifact | SHA256 |
| --- | --- |
| `topology-pass/chapter-input.npz` | `6e4904482c10530728705f0c4c1ce8bc943db6c2f0aa6aedfba5bb99284aaffd` |
| `refresh-pass/smoke/final_state.eqx` | `ab6c1a24616d0ecb666bfb94fd861904546711a8efe17b44c204fcb78500bd0d` |
| `propriety-design.npz` | `8e17bc754745e7549a4aaf32aceeda94a5b78c5a5d67f667953b18f6fc7171f5` |
| `propriety-responses.npz` | `7e84161c4fdcf283d514e41fa53947d9b21a9384982ebc014881f9cf6c33b89b` |
| `propriety-profile-topologies.npz` | `af323e03f9398f25a1f2658bb7d5ed3175e76e3f8f95ad211325ec111d52e5ba` |
| `propriety-exposure-tree.npz` | `3528616e1f6307800439e606e110789d8270f1e8ebcc057272636660491b2e53` |
| `propriety-boundary-target.npz` | `197ca7c0cb21a5281a4c7b5f854b7c86498c7d6cd2eb46d2e0b70351824bbdff` |
| `propriety_export.py` | `19e90fb3fb34f6c42fb46a5928d2b504f9c09ed7b5c4168dee8b5a704908680a` |
| `propriety_rank_witness.py` | `e34d089582c699fae46a304e96d5bc1c9183b4ebf8c04c5eabe22b2895139850` |
| `propriety_verify.py` | `d4cbf1bea607fdfb7dc87bd9462dca234e6aa2c5e21b90fb21d4dd78c954f5f5` |
| `propriety-witness.json` | `e34861403f8d14a152dc09369aa2dcc09faf7ad8b518507568cfc35828dcd195` |
| `propriety-independent-review.json` | `53fa2ba12415ae2e3c7b8493f64d4a18aee378ed46c0c717e02dc459327a91ac` |

Paths in the first two rows are relative to
`/home/ignacio/longbet-bug-audit`; the remaining artifact names are relative
to its `refresh-pass` directory. These checks establish the original-model
normalization defect, not a successful inference repair.
