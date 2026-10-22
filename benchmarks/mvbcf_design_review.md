# Ideas from mvbcf for LongBet (2026-09-07)

## Scope and evidence

Inspected Nathan-McJames/mvbcf at commit
`fc3b89b0a78ce8a31ae75c43a6ec75f1945ca0c8`, including the R front door,
vignettes and C++ leaf/tree updates. This is a source/design review, **not an
empirical validation of that package**. It was not fitted to the chapter data.

Primary sources:

- [Model paper, Section 3 and Algorithm 1](https://doi.org/10.1093/jrsssa/qnae049).
- [Pinned C++ implementation](https://github.com/Nathan-McJames/mvbcf/blob/fc3b89b0a78ce8a31ae75c43a6ec75f1945ca0c8/src/mvbcf_source_code.cpp).
- [Pinned R interface and defaults](https://github.com/Nathan-McJames/mvbcf/blob/fc3b89b0a78ce8a31ae75c43a6ec75f1945ca0c8/R/mvbcf_function.R).
- [Author's full guide](https://nathan-mcjames.github.io/mvbcf/articles/full-guide.html).

The initial repair was full-precision SUR with separate forests. The user
subsequently authorized the shared/private extension below, now implemented
as an **opt-in experimental model**, not a replacement default. See the
[implementation and validation report](shared_private_review.md) for its
actual behavior and benchmark results. This document records the design
rationale; the list of required comparisons is not a claim they all passed.

## Most useful distinction: sharing heterogeneity versus sharing noise

mvbcf has one prognostic ensemble and one treatment ensemble whose trees have
vector-valued leaves. A split is selected using a multivariate likelihood;
each outcome gets its own value in the resulting leaf. LongBet's default has
one pair of forests per outcome, coupled only through the residual likelihood;
the new opt-in model adds shared treatment partitions with private departures.

For the chapter, a common partition could pool evidence about whether catalog
size and fulfillment modify treatment. Within that partition, GMV, hours and
complaints may have different magnitudes and opposite signs. **Sharing a split
does not require sharing a sign, a value, or a benefit/harm label.** The weaker
binary response could help choose splits and benefit from splits supported by
the continuous responses. This is a plausible benefit, not a measured result.

Residual correlation alone need not produce much extra precision when each
outcome has similar predictors and its own unconstrained mean. Sharing
partitions supplies an additional modeling assumption about common structure.
This also explains why a correct SUR implementation is not a guarantee that
joint fits beat separate models in effect recovery.

## Separate the different covariance assumptions

mvbcf exposes separate fixed prior covariance matrices `sigma_mu` and
`sigma_tau` for leaf vectors, and samples the residual covariance `Sigma`.
Its default leaf covariance matrices are diagonal; the interface does **not**
automatically learn a treatment-effect correlation matrix. Shared topology can
pool evidence about splits even with diagonal leaf priors.

For LongBet, keep prognostic pooling and treatment-effect pooling distinct.
Baseline correlation need not imply shared treatment heterogeneity. If signed
cross-outcome leaf shrinkage is added, specify its prior separately from the
error covariance. Do not estimate it by plugging in the sample correlation of
observed outcomes or by using the known simulated treatment effects.

The inspected defaults also put a much stronger shallow-tree prior on treatment
trees (`alpha_tau=.25`, `beta_tau=3`) than prognostic trees (`alpha=.95`,
`beta=2`). **These treatment-depth defaults already match LongBet's defaults**;
copying them is not a new repair. They motivate a declared prior-sensitivity
comparison with a flexible baseline and simpler heterogeneity. LongBet's
GP and coding multipliers also affect the induced effect prior. Apply any such
change to the separate-fit baseline too, and check coverage as well as mixing.
No new tree-depth or GP prior defaults were selected in this repair.

## Extension design (now implemented; validation remains a separate question)

1. Keep the repaired separate-forest SUR model as the baseline. Do not silently
   change the meaning of existing fits or `sur=False`.
2. Prototype **shared treatment partitions with outcome-specific leaf values**,
   initially leaving prognostic forests outcome-specific. Retain LongBet's
   calendar/exposure design and outcome-specific GP trajectories. Define one
   shared predictor layout, cutpoint grid and split prior explicitly; reject
   incompatible per-outcome tree settings instead of silently using the first
   outcome's configuration. Preserve the documented external outcome order
   when mapping binary-first internal arrays and joint posterior draws.
3. Add a regularized outcome-specific treatment component, or another explicit
   escape from sharing, before adopting shared partitions broadly. This is our
   proposed protection against negative transfer, not a feature claimed for
   the inspected mvbcf package. Put declared priors on the common/private
   variance allocation so adding ensembles does not accidentally double the
   prior effect variance. Check posterior identification and mixing.
4. Introduce sharing through an explicit opt-in setting. With sharing off and
   SUR off, require matched-key scalar equivalence. With shared trees on,
   `sur=False` means independent residual errors, **not independent models**.
5. Version trace layouts, semantics, R serialization and cache keys; shared
   topology cannot be represented honestly as independently fitted forests
   without preserving the joint draw provenance.

### Correct weighted vector-leaf updates

The following is a proposed LongBet adaptation, not a verbatim mvbcf formula.
For one shared tree leaf with outcome-vector parameter `theta`, prior
`N(0,V)`, partial residual `r_i` excluding that tree, outcome precision `Omega_i`
and diagonal multiplier `D_i`:

```
P = inverse(V) + sum_i D_i' Omega_i D_i
h = sum_i D_i' Omega_i r_i
theta | rest ~ Normal(solve(P,h), inverse(P))
```

For a treatment tree, the diagonal entries are the **outcome-specific** current
`b_m(Z_i) * beta_m(S_i)`. For a prognostic tree they are `alpha_m`.
Specify `V` on the internal standardized continuous/identified binary-latent
scale, and document the induced prior after conversion to each response scale.
An equal numeric variance in dollars and probabilities is not equal shrinkage.
Use Cholesky solves and triangular Gaussian sampling, not repeated dense
inverses. Use the corresponding integrated vector-leaf likelihood in every
tree acceptance calculation, including determinant terms and the complete
forward/reverse proposal probabilities. Do not reuse scalar leaf precision
by merely broadcasting across outcomes. Heterogeneous multipliers mean the
multivariate leaf precision generally differs by leaf.

The new dynamic-precision cache tests must extend to matrix-valued leaf sums,
changing GP/coding weights, changed covariance draws, and deferred prunes.
Masked cells require the observed-data likelihood; never replace missing
responses by zero and count them as observed.

### Binary and longitudinal outcomes are not a drop-in extension

The inspected mvbcf interface/kernel uses Gaussian response matrices and an
unconstrained inverse-Wishart residual covariance update. Its binary `Z` is
the **treatment**, not a supported binary response likelihood. Passing 0/1
complaint observations as an ordinary Gaussian outcome would not implement
LongBet's probit model.

Keep latent probit augmentation, downstream conditional means/variances, and
the identification constraints. Do not replace the current covariance block
by unconstrained inverse-Wishart sampling when some marginal latent variances
must equal one. mvbcf also does not provide LongBet's exposure GP, absorbing
rollout semantics or panel unit-intercept machinery; these must be retained
and included in the joint likelihood, not dropped in a benchmark conversion.

The source is useful for design, not a trusted replacement transition kernel.
In particular, its grow/prune/change/swap acceptance code warrants a separate
proposal-balance audit before any port: the inspected acceptance expression
is `exp(lnew-lold)`, while its tree score visibly includes integrated likelihood
and depth-prior terms. This review has not independently established that its
proposal probabilities cancel. Do not copy it on the assumption that matching
an existing implementation proves posterior correctness.

## Required comparison

Compare corrected separate fits, corrected separate-forest SUR, and the
proposed shared/private model, with matched declared computational budgets and
also ESS per second. Keep the unchanged chapter panel and prespecified
benefit/harm groups. Add repeated simulation seeds with:

- shared sign-changing modifiers, including opposite signs across outcomes;
- weak or zero residual correlation but shared treatment modifiers;
- correlated residuals but different treatment modifiers;
- a null-effect outcome, to detect induced false positives;
- weak/rare binary outcomes and supported missingness patterns;
- outcome-order permutations and prior sensitivity.

Check actual vector-leaf Gaussian posterior moments, tiny-tree enumerated
posteriors/detailed balance, mixed latent-response conditionals, interval
coverage and joint-event calibration. Diagnose risk differences and group
averages with chain boundaries retained. A lower RMSE or narrower interval
without correct coverage and convergence is not an improvement.

Recommendation: this provides a mechanism to **learn shared heterogeneity**,
particularly for the weaker complaint outcome. Implementation alone must not
be described as demonstrated superiority, nor used to justify restoring the
chapter example prematurely. Consult the implementation report for evidence.
