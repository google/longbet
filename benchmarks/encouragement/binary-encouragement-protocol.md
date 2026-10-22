# Single-horizon binary encouragement benchmark

This research experiment follows the review of the
[StochTree IV vignette](https://stochtree.ai/vignettes/iv.html). It adds no public
LongBet backend. The existing longitudinal sampler and calibration records stay
available while a binary observed-data model is evaluated independently.

## Model and target

For binary adoption D, outcome Y and randomized encouragement Z, fit
`p_z(x)=Pr(D=1|Z=z,X=x)` and `q_dz(x)=Pr(Y=1|D=d,Z=z,X=x)`. Their product defines
one joint distribution of (Y,D) given (Z,X). Derive each outcome reduced form as
`r_z=p_z*q_1z+(1-p_z)*q_0z`, and average contrasts over the same original empirical
baseline distribution. The outcome regression is observational; conditioning on
adoption does not make its coefficient causal.

Compare an exact finite-cell posterior, a probit forest, and an independent
reduced-form cell baseline. The finite-cell model has independent Beta(1,1)
priors on the conditional probabilities. Its monotone version explicitly
conditions the uptake prior on `p_1>=p_0`. A bounded exact rejection sampler must
report failure if its computational budget is exhausted. The forest instead uses
the monotone-probit construction `p_1=Phi(f)`, `p_0=Phi(f)*Phi(g)`; these are
different priors and must not be described as a sampler-only comparison.
An explicit null-stage cell model has `p_1=p_0` and an undefined Wald ratio.

## Identification

Monotonicity is a declared no-defiers assumption, not something proved by ordered
observed uptake rates. For that model, estimate always-taker, complier and
never-taker probabilities. Bound complier outcome probabilities using the
observed mixtures. Aggregate bounds using complier masses and original-unit
weights, avoiding division by a near-zero cell-level complier share.

Report the effect of encouragement plus adoption among compliers separately from
the treatment effect at fixed encouragement. A sensitivity parameter delta
limits the absolute direct encouragement effect on each stratum's conditional
outcome **probability**. It is not a restriction on each binary potential-outcome
difference. Delta=1 leaves that direct effect unrestricted; delta=0 imposes
exclusion. Empty sets remain empty. An unconstrained observed-data posterior is
not silently filtered into an exclusion-restricted posterior, and its raw Wald
ratio is not relabeled a proper posterior on a bounded binary CACE.

## Validation and experiment

Before fitting simulations, check Gaussian leaf conditionals, truncated-normal
tails, the joint monotone augmentation, residual invariants, and Geweke prior
preservation. Deliberately broken updates must fail. Compare exact cell posterior
moments and covariance with analytic Beta identities; compare bounds with an
independent linear-program formulation.

Generate all potential outcomes and treatment responses before complete balanced
random assignment. Use prespecified baseline categories, heterogeneous compliance
and outcomes, and strong, weak, zero, one-sided, exclusion-violation and defier
scenarios. Record both realized finite-population truths and expected truths at
the original empirical X distribution. Sampler and assignment seeds are separate.

Use larger repeated-population studies for the fast cell controls, and a bounded
four-chain forest study before expansion. Record full ITT covariance, bias,
coverage with Monte Carlo uncertainty, interval widths, weak-stage probabilities,
and identification-bound uncertainty. Compare each estimator's posterior/estimated
covariance with its own empirical estimation-error covariance; covariate adjustment
does not require equality with the unadjusted reference variance. Reference normal
Fieller/AR and randomization AR results retain their respective validity conditions.
The randomization AR grid is not relabeled a complete confidence set.

Raw draws, seeds, data, failures and configurations are archived. Diagnose both
ITTs, the raw ratio, bound endpoints, first-stage probabilities and forest
functions. Passing R-hat/ESS checks does not establish coverage. Priors differ
between the finite-cell and forest controls; neither is promoted by a small
simulation study. Longitudinal adoption histories remain a subsequent stage.

## Completed experiment and sampler amendment

The [implementation and numerical report](../benchmarks/encouragement/binary-report.md)
records 112 passing focused tests, 600 independent datasets, two additional
prior settings on the same weak/zero datasets, and 24 forest fits. Full original
and modified sampler records, failed cell fits, source hashes and replay arrays
are retained in `benchmarks/encouragement/binary_results/`.

The initial forest pilot revealed slow mixing for monotone uptake near a
probability boundary. In response, an optional marginal elliptical slice update
was added for uptake leaf coefficients. It integrates the auxiliary binary
labels and probit utilities out of that step and regenerates them next sweep.
This amendment changes the transition, leaving the model and priors fixed.
The original runs remain available; the extra eight runs compare the update on
the same datasets and chain lengths. Both Geweke checks and an independent exact
posterior comparison validate the additional transition.

All sampled diagnostics pass in the eight modified monotone fits. That result
does not resolve the cell study's weak-stage bias under monotonicity, the lack
of exact-null mass in continuous ordered priors, or the remaining need for a
larger forest coverage study.

## Using the research model

From a checkout of this repository:

```python
from benchmarks.encouragement.binary_dgp import make_binary_data
from benchmarks.encouragement.binary_forest import BinaryConfig, fit_binary
from benchmarks.encouragement.binary_cells import effects, identification_bounds

data = make_binary_data("strong", n=240, population_seed=123, assignment_seed=456)
draws = fit_binary(
    data["y"], data["d"], data["z"], data["x"], seed=789,
    chains=4, burnin=1000, draws=4000,
    first_stage="monotone",  # explicitly assumes no defiers
    config=BinaryConfig(marginal_uptake_refresh=True),
)
itt = effects(draws)  # raw Wald ratio is not automatically a causal CACE
bounds = identification_bounds(draws["p"], draws["q"], draws["weights"], delta=.1)
available = bounds["available"]  # retain empty and undefined draws
```

`first_stage="unrestricted"` fits two unconstrained uptake forests and is the
default. Principal-stratum bounds require a declared no-defiers model.
`delta` can be a scalar or three values for always-takers, never-takers and
compliers; these constrain conditional mean risk differences. Inspect
`compatible`, `defined` and `available` before interpreting endpoints. The raw
archive includes all probabilities, leaf coefficients, rules, RNG state and
configuration; `binary_benchmark.summarize` produces the full diagnostics used
in this report. The research modules are not installed public APIs.
