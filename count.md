# Count outcomes: implementation specification

Status: proposed; no count backend is implemented by this document. Reviewed
against the working tree on 2026-09-12. Implement the numbered phases in section
10 and record evidence for their gates. Section 4 contains a **required sampler
feasibility gate**, not an assertion that an exact general-shape JAX sampler
already exists.

This specification is self-contained. Use [README.md](README.md) and the source
links below for existing behavior; do not depend on the removed `ordinal.md`.
Preserve unrelated working-tree changes during implementation.

## 1. Deliverable and scope

Implement scalar `LongBet(outcome="negbin")` in Python and R, with observation
exposure, inferred or fixed dispersion, count-scale posterior means and ATTs,
incidence rate ratios (IRRs), bounded-memory prediction, diagnostics, and archive
round trips. Preserve existing continuous, binary, ordinal, and multi-outcome
behavior, including their random-key schedules.

Decisions for the first release:

| Area | Required behavior |
|---|---|
| Likelihood | Negative binomial with log **mean** link and one global dispersion per chain. |
| Augmentation | Pólya–Gamma (PG), subject to section 4's feasibility gate. Gaussian parameter conditionals do not by themselves establish an exact sampler. |
| Observation exposure | `exposure=` on `fit` and `predict`; positive amount at risk, distinct from treatment duration `s` / `exposure_idx`. |
| Fixed offset | `log(exposure)` has coefficient exactly one, outside the prognostic forest and its `alpha` scaling. No separate public `offset=` alias in this release. |
| Estimands | Absolute expected-count ATT, arithmetic mean cell IRR, and ratio of aggregate expected counts; define each separately. |
| Hurdle, zero-inflated NB, joint count outcomes | Deferred to section 11. Reject these unsupported combinations explicitly; do not register an unfinished `"hurdle"` outcome. |
| Poisson | A likelihood limit at fixed mean, not a separately implemented backend. |
| Existing ordinal alternative | For a known finite support `0..K-1`, document `outcome="ordinal", num_categories=K` and `att_expected_score(weights=np.arange(K))`. It imposes ordinal rather than NB assumptions and cannot extrapolate beyond the declared categories. No arbitrary cutoff such as `Y <= 8`. |

Keep tree topology Metropolis–Hastings moves, Gaussian leaf integration, GP,
adaptive coding, and unit-intercept updates. Do not describe tree topology or
variance updates as Gaussian. Do not automatically switch likelihoods based on
the fraction of zeros, or transform the count response to `log(y + 1)`.

## 2. Statistical contract

### 2.1 Predictor and negative-binomial parameterization

Let `n = N*T` cells be flattened in unit-major order. Let `C` be chains, `L`
retained sweeps per chain, and `D = C*L` flattened draws in chain-major order.
Use `A_it` for treatment and `e_it > 0` for observation exposure. Obtain treatment
duration `S_it` from the existing `derive_exposure(z, t)`; untreated cells have
`S_it=0`. Do not create a second adoption-index implementation.

Define the forest predictor without the observation offset:

$$
f_{it}=\alpha m_{it}+b_{A_{it}}\beta_{S_{it}}\nu_{it}+\gamma_i,
\qquad m_{it}=m_0+f_\mu(X_i,t,X^{\mathrm{tv}}_{it}).
$$

Here `m_it` is exactly `state.mu_fit`, including bartz's scalar forest offset
`m0`. The treatment forest uses the current treatment design. Keep existing
priors for forests, `alpha`, coding, GP, and unit intercepts; their scales are
now log-rate units. Defaults currently give forest standard deviations `1.5`
and `1.0` before multiplication by the scale parameters. Record these priors in
the benchmark; do not infer their scale from the spread of raw counts.

Set

$$
o_{it}=\log e_{it},\quad \eta_{it}=o_{it}+f_{it},\quad
\lambda_{it}=\exp(\eta_{it}),\quad \psi_{it}=\eta_{it}-\log r.
$$

The NB probability mass function is

$$
p(y\mid\lambda,r)=\frac{\Gamma(y+r)}{\Gamma(r)\Gamma(y+1)}
\left(\frac{r}{r+\lambda}\right)^r
\left(\frac{\lambda}{r+\lambda}\right)^y,
\qquad r>0.
$$

Thus `E[Y]=lambda` and `Var[Y]=lambda+lambda**2/r`. Holding `eta` fixed while
varying `r` leaves the mean fixed. This removes the earlier draft's unintended
`r*exp(eta)` mean. The Poisson limit is `r -> infinity` at fixed `lambda`.
For NumPy/SciPy NB calls the success probability is `r/(r+lambda)`, whereas
`expit(psi)=lambda/(r+lambda)` is its complement. Test this convention against
[SciPy's NB definition](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html).

Keep current coding semantics: adaptive coding starts at `b0=-0.5, b1=0.5`;
turning it off fixes **both** to one. Neither `b0` nor `beta[0]` is fixed to zero.
The control contribution is part of the counterfactual model, not a standalone
pre-trend test.

Causal interpretation requires the existing treatment-assignment, overlap,
consistency, and missingness assumptions. Both arms below use the same specified
`e_it`. If treatment changes the amount at risk, this is a contrast holding
exposure fixed; it does not estimate the total effect through changing exposure.

### 2.2 PG working likelihood and the residual invariant

With `h_it=y_it+r` and `kappa_it=(y_it-r)/2`, the likelihood kernel in `psi` is
`exp(kappa*psi)/(2*cosh(psi/2))**h`. PG augmentation gives

$$
\omega_{it}\mid y,f,r\sim\mathrm{PG}(y_{it}+r,\psi_{it}),\qquad
\log p(f\mid\omega,r,y)=
-\tfrac12\sum_{\mathrm{obs}}\omega_{it}(q_{it}-f_{it})^2+\log p(f)+\text{const},
$$

$$
q_{it}=\frac{\kappa_{it}}{\omega_{it}}-o_{it}+\log r.
$$

This is a conditional likelihood kernel; `q` is a working response, not a newly
observed Gaussian outcome. The mixing identity uses `PG(h,0)`; the conditional
draw uses `PG(h,psi)`. See [Polson, Scott, and Windle (2013)](https://arxiv.org/abs/1205.0310).

For the count path, reuse inherited `state.z` for **offset-adjusted** `q`:

```python
state.y       # original counts, finite zero placeholders for missing cells
state.z       # q, NOT treatment z_vec and NOT kappa / omega alone
state.resid   # where(obs_mask, state.z - f, 0)
state.sigma2  # exactly 1; conditional precision supplies the observation scale
```

These equalities must hold after initialization, chain overdispersion, every
sweep, and every enabled move. A change to the cell offset or dispersion requires
rebuilding `q` before using any Gaussian update. Never learn or multiply the
cell offset through `alpha` or a forest leaf.

### 2.3 Dispersion move and partially collapsed ordering

Use the shape/rate prior `r ~ Gamma(a_r, b_r)`. For a log-dispersion proposal
`u_new = u + proposal_sigma * Normal(0,1)`, `u=log(r)`, evaluate the **observed-data**
log target, integrating out PG variables and holding `f` and `o` fixed:

$$
\begin{aligned}
\ell(u)&=\sum_{\mathrm{obs}}\{\log\Gamma(y+e^u)-\log\Gamma(e^u)
-\log\Gamma(y+1)\\
&\quad-y\,\operatorname{softplus}[-(\eta-u)]
-e^u\,\operatorname{softplus}(\eta-u)\}
+a_r u-b_r e^u.
\end{aligned}
$$

Accept if `log(U) < min(0, ell(u_new)-ell(u))`. The `a_r*u` term includes the
log-transform Jacobian. Recompute `psi=eta-u` for each candidate; dropping the
`y` term because it was constant in the old parameterization is incorrect.
A nonrepresentable proposal is rejected, never clipped into the parameter space;
a nonfinite current target is a surfaced numerical error.

Required sweep order:

1. Reconstruct current `f` from cached fits and parameters, then `eta=o+f`.
2. Update `r | y, eta` using the marginalized MH move, or keep its configured
   fixed value. Old PG variables are discarded for this block.
3. Immediately draw fresh `omega | y, eta, r`, then rebuild `q` and `resid`.
4. Run one complete weighted Gaussian mean-parameter sweep with `r` fixed.
5. Apply existing likelihood-preserving ridge moves as enabled; return state.

No update may condition on stale PG variables between steps 2 and 3. Do not
interpret the marginal NB target as `p(r | omega, eta, y)`. The prescribed order
combines a marginal `r` transition, regeneration of the augmentation, and
conditional mean updates. Verify the joint stationary distribution using the
low-dimensional posterior oracle in section 9.

## 3. Public inputs, configuration, and initialization

### 3.1 Configuration

Extend [LongBetConfig](src/longbet/_config.py) with:

```python
outcome: Literal["continuous", "binary", "ordinal", "negbin"] = "continuous"
dispersion_prior_shape: float = 2.0
dispersion_prior_rate: float = 0.5
dispersion_proposal_sigma: float = 0.15
sample_dispersion: bool = True
initial_dispersion: float | None = None
```

Require finite, strictly positive real values for the three numeric dispersion
options and for a non-`None` initial value; reject booleans masquerading as
numbers. Require a boolean `sample_dispersion`. When `sample_dispersion=False`,
require an explicit `initial_dispersion` and keep it identical across chains.
When sampling and no initial value is supplied, initialize at the prior mean
`shape/rate`. Proposal tuning is fixed in this release; no adaptation during
retained sampling.

Require `num_categories=None` for NB. Leave the public `standardize=True`
default compatible with other outcomes but document that it applies only to
continuous responses. For NB always use `meany=0`, `sdy=1`, and no outcome
standardization. Do not mutate the frozen config to pretend the user set False.
Error-variance prior options are inactive for NB, just as its working variance
is fixed at one. Dispersion options are inactive for non-count outcomes.

### 3.2 `fit` and `predict` exposure contract

Append keyword-only `exposure=None` **after existing parameters** in both
[LongBet.fit and LongBet.predict](src/longbet/_model.py), forwarding it through
their internal implementations. Do not insert it between existing positional
arguments such as `t` and `x_trt`.

- Count `y` accepts shape `(N,T)` or `(N*T,)` in unit-major order. Validate the
  original dtype before coercion: numeric integers or real floats only; reject
  boolean, string, object, and complex arrays. Only `NaN` denotes missingness.
  Reject infinities, negatives, fractions, and an entirely missing panel.
  Accept all-zero and other constant count panels.
- Validate before float32 conversion so rounding cannot turn a fraction into
  an integer. The initial float32 engine supports integer counts through
  `2**24`; reject larger values with an explicit precision-limit error. This
  bound does not imply that all large-count panels are computationally cheap.
- `exposure` accepts `(N,T)` or `(N*T,)`, never implicit row/column broadcasting
  or a scalar. Require numeric, finite, strictly positive values at **every**
  cell, including cells with missing `y`, because prediction covers them too.
  Values below one are valid: negative log offsets must be accepted.
- At fit, `None` means a panel of ones. Store the normalized host float64
  `log_exposure_fit_` and a boolean `requires_exposure_` indicating whether the
  user supplied exposure explicitly. Compute logarithms before device casting.
- At prediction, require exposure if `requires_exposure_` is true, even for a
  panel with the fitted shape. Otherwise `None` means ones. Explicit exposure
  may be supplied to any count prediction, including a unit-exposure fit.
  Shape alone cannot establish panel identity, so do not reuse fitted offsets
  implicitly. Preserve this rule across save/load.
- Reject non-`None` exposure for other outcomes and unsupported multi/count or
  encouragement/count combinations. Scalar `sur` retains its current inert
  behavior; shared trees remain unsupported for scalar fits.

Exposure is fixed response metadata, not a binned predictor. Do not add it to
[Design](src/longbet/_design.py), allow trees to split on it automatically, or
confuse it with the existing raw design source `"s"`.

### 3.3 Initialization

On observed cells in host float64, set

$$
m_0=\log\left(\frac{\sum y+0.5}{\sum e}\right).
$$

Use a stable log-sum calculation if the sums overflow. The half-count makes
all-zero initialization finite; it is an empirical initialization/prior-centering
choice, not an added observation in the likelihood. Do not use the variance of
`y/e` as a method-of-moments estimator of dispersion under unequal exposures.
Initialize `r` as specified above, `mu_fit=m0`, and `nu_fit=0`; retain current
starting values for the remaining parameters.

Construct both bartz forests through their Gaussian path using a finite
log-scale initialization response, not raw counts as the Gaussian outcome.
Retain the original counts in `LongBetState.y`, set `binary_indices=None`,
`num_categories=0`, and an empty cutpoint array. Ensure both forest views have
fixed observation variance; bartz must neither draw binary latents nor resample
a Gaussian error variance for this outcome.

After all parameters (including chain-specific intercepts) are initialized,
compute actual `f`, draw PG variables, and construct `q` and the residual. For
multiple chains, preserve existing parameter overdispersion; additionally use
`r_c=r_initial*exp(0.3*Normal(0,1))` **only when dispersion is sampled**. Set
`sigma2=1` in every chain. Regenerate PG variables after these changes rather
than broadcasting one chain's working response. Initialize the single-chain
path too; it does not pass through `_overdisperse_chains`.

## 4. Pólya–Gamma implementation feasibility gate

Create `src/longbet/_count.py` with pure functions for
validation/preparation, stable NB log mass, PG draws, working-response
construction, and dispersion MH. The following is the intended PG interface:

```python
def sample_polya_gamma(key, h, tilt):
    """Return draws and per-cell failure flags, with h/tilt's broadcast shape."""
```

Before wiring a public backend, produce `benchmarks/count_sampler_review.md`
with the selected algorithm, primary reference and algorithm location, any
ported code's revision/license, supported parameter domain, exactness claim,
numerical failure policy, and measured CPU/GPU behavior.

The Gamma dispersion prior assigns positive probability to `0<r<1`; an observed
zero then requires `PG(h,tilt)` with `0<h<1`. Fractional and large `h` are ordinary
cases, not optional extensions. Devroye's unit-shape sampler plus integer
additivity does not cover this domain. A truncated gamma series and a normal
approximation change the augmentation distribution. A moment-matched normal is
not a saddlepoint rejection sampler. Windle et al.'s alternate construction
covers `h>=1` subject to its stated envelope conjecture; its large-shape section
is explicitly approximate. See [Windle, Polson, and Scott (2014), alternate and
approximate samplers](https://arxiv.org/html/1405.0506).

**Gate:** identify and implement a justified sampler covering all positive real
shapes admitted by the model, and pass section 9's distribution tests. Distinguish
mathematical exactness from finite-precision error. A dependency's name or a
hybrid default is not evidence that its draws are exact. If that cannot be
established, report the unresolved domain and keep public NB support unreleased.
Do not silently restrict the Gamma prior, round `r`, clip PG draws, or substitute
an approximation. An explicitly approximate backend requires a revised design
with approximation controls and posterior sensitivity evidence; it is not
implicitly authorized by this specification.

JAX implementation requirements:

- Work under `jit`, single-chain execution, and `vmap` over chains. Use JAX keys
  exclusively in the production transition; no host RNG callbacks in the sweep.
- Use independent keys for cells/proposals and masked rejection loops which
  preserve accepted cells. Finite iteration guards must return failure, never
  a moment approximation or the last rejected proposal.
- Keep the existing float32 forest engine. Use scoped float64 for PG or NB
  arithmetic where required, without changing global JAX precision or legacy
  kernels. The current sweep is decorated with `enable_x64(False)`; explicitly
  test nested precision behavior rather than assuming requested float64 survives.
- Accepted observed PG values and working responses must be positive/finite
  and finite, respectively. Underflow, overflow, or non-convergence must surface
  as an error at the public fit boundary after synchronization. Do not clamp
  `omega` to `1e-12` or silently discard failing observed cells.
- Handle zero tilt with analytic limits, negative tilt by symmetry, and large
  absolute tilt with stable expressions rather than overflowing `sinh/cosh`.
  Stream series/proposal work; avoid an unbounded `(chains,cells,terms)` array.

## 5. State and sweep integration

Relevant files are [state construction](src/longbet/_state.py),
[weighted sweeps](src/longbet/_step.py), [scale attributes](src/longbet/_scales.py),
[precision caches](src/longbet/_forest_cache.py), and [GP updates](src/longbet/_gp.py).

Add optional fields with non-count defaults so existing state constructors
continue to work:

| Field | Single-chain shape / role | Axis handling |
|---|---|---|
| `dispersion` | scalar `r`; `None` for other outcomes | `chains=CHAIN_AXIS` |
| `omega` | `(n,)`; `None` for other outcomes | `chains=CHAIN_AXIS, data=-1` |
| `log_exposure` | `(n,)`; `None` for other outcomes | shared, `data=-1` |
| Existing `z` | count working response `(n,)` | inherited chain metadata; do not redeclare |
| Dispersion options | prior, proposal scale, sampling flag | static |
| `count_failed` | sticky scalar boolean, initially false | per chain |
| `dispersion_mh_accepted`, `dispersion_mh_proposed` | scalar cumulative integer counters | per chain |

Keep `y`, the design, observation mask, panel indices, and offsets shared.
Update `chain_filter_spec` tests to verify that only parameter/working arrays
acquire a chain axis. Do not put observation-sized arrays in `LongBetConfig`.

For missing outcomes, sanitize inputs *before* divisions, PG calls, logarithms,
and likelihood evaluation. Use harmless PG arguments `(h,tilt)=(1,0)` on masked
cells; their retained `omega=1`, `q=0`, `resid=0`. Effective information precision
is `P=where(obs_mask,omega,0)`. Likelihoods and sufficient statistics must exclude
these cells; multiplying an existing NaN by zero is not sufficient.

At the top of `longbet_single_step`, branch on `outcome_type_str=="negbin"`,
perform section 2.3's dispersion/PG block, and set local
`conditional_precision=omega` with `obs_mask` still authoritative. The count
branch must reject externally supplied SUR precision and `shared_treatment=True`.
Construct `conditional_attrs` **after** refreshing PG weights. Reuse the current
Gaussian update body rather than recursively calling the same count branch.

Audit every conditional, not just the two forests:

| Update | Required observation information |
|---|---|
| Prognostic leaves | `P*alpha**2`, accounting for bartz residual and leaf units. |
| `alpha` | `sum(P*mu_fit**2)` plus its existing prior precision. |
| Treatment leaves | `P*(b_z*beta[S])**2`; preserve the existing treatment-multiplier exclusion threshold independently of `P`. |
| GP | Group weighted squares and weighted residual products by `exposure_idx`. |
| `b0`, `b1` | Group weighted squares/products of `beta[S]*nu_fit` by treatment arm. |
| `gamma_i` | `sum_t(P_it)` plus `1/sigma_gamma2`, and the corresponding weighted residual sum; not unweighted `unit_counts`. |
| `sigma_gamma2` | Existing inverse-gamma draw from sampled intercepts. |
| `sigma2` | Always one; never the Gaussian residual-based draw. |

Refresh bartz leaf precision caches whenever PG, coding, or trajectory weights
change. Preserve shared baseline observation-scale fields when returning from
the conditional path; do not leak a chain's temporary precisions into shared
state. Retain `z=None` on the temporary forest views to prevent a second latent
update. Preserve the ridge product and the residual invariant after ridge moves.

Use count-only fold-in streams (reserve tags `8201` for initial PG, `8202` for
initial dispersion spread, `8211` for sweep dispersion, `8212` for sweep PG,
after checking for collisions). Split within those streams as needed. Preserve
existing `random.split(key,10)` and five-key chain-initialization schedules for
all other outcomes.

## 6. Prediction and estimands

### 6.1 Paired arms

Within each existing prediction cell block, evaluate the treatment forest at
factual `S` and control `S=0` using the same draw. Preserve fitted-grid clipping
for the forest and GP projection beyond the fitted horizon. Compute

$$
\begin{aligned}
f_0 &= \alpha m+b_0\beta_0\nu(X,0,t)+\gamma_i,\\
f_1 &= \alpha m+b_1\beta_S\nu(X,S,t)+\gamma_i,\\
\lambda_0 &= e\exp(f_0),\quad \lambda_1=e\exp(f_1),\\
\rho_a &= \lambda_a/e=\exp(f_a),\\
\delta &= \lambda_1-\lambda_0,\qquad I=\exp(f_1-f_0).
\end{aligned}
$$

`lambda` is expected **count**; `rho` is expected count **per unit exposure**.
Use `lambda_factual=where(A==1,lambda1,lambda0)`. These are posterior draws of
conditional means, not noisy posterior predictive observations. They require
no Gaussian smearing factor and no multiplication by `r`.

Maintain current row-based unit-intercept matching: same fitted units must be
in the same order. If row count differs, retain the existing warning and
`gamma=0` behavior and describe the result as conditional on zero intercept.
Do not call it a population-marginal mean: integrating a normal intercept would
change a log-link mean. New-unit identification/marginalization is outside scope.

For `T_s = {it: A_it=1, S_it=s}` in the **prediction panel**, compute per draw:

$$
\mathrm{ATT}_s=\operatorname{mean}_{T_s}(\delta),\quad
\mathrm{IRR}_s=\operatorname{mean}_{T_s}(I),\quad
\mathrm{AggregateRatio}_s=\frac{\sum_{T_s}\lambda_1}{\sum_{T_s}\lambda_0},\quad
\mathrm{Lift}_s=\mathrm{AggregateRatio}_s-1.
$$

The first IRR is an arithmetic average of cell ratios. It is generally unequal
to the ratio of sums or `exp(mean(log(I)))`. Lift is a fraction; multiply by 100
only for percentage presentation. Aggregate draws before posterior summaries.
Missing training outcomes do not remove cells from this prediction target.
Return `NaN` draws and zero support for unsupported horizons, following the
current `max(S_max_test,1)` convention for an entirely untreated panel.

### 6.2 Public prediction contract

Extend [LongBetPrediction](src/longbet/_model.py) with optional count-only fields;
existing constructors and other outcomes must remain compatible.

| Field | Shape / NB meaning |
|---|---|
| Existing `yhats`, `muhats0`, `tauhats` | `(N,T,D)` factual expected count, control expected count, and count difference; `None` in summary mode. |
| Existing `y_summary`, `mu0_summary`, `tau_summary` | Exact per-cell summaries on those same count scales, `(N,T)`. |
| Existing `att_full` | `(S_out,D)` absolute count ATT, float64 for NB. |
| `att_irr_full` | `(S_out,D)` arithmetic mean cell IRR. |
| `att_aggregate_ratio_full` | `(S_out,D)` ratio of aggregate expected counts. |
| `dispersion_samples` | `(D,)`, flattened in exactly the forest draw order. |
| `exposure` | Normalized prediction exposure `(N,T)`. |
| `rate_summary`, `rate_control_summary` | `(N,T)` summaries of factual and control `rho`, not `lambda`. |

Methods:

- `pred.att(alpha=0.05)` and `pred.catt()` retain their existing result schema
  and now report counts for NB; `pred.stability()` diagnoses count ATT draws.
- `pred.irr(alpha=0.05, aggregation="mean")` returns `{"irr": (S_out,),
  "intervals": (2,S_out), "irr_full": (S_out,D), "exposure": (S_out,),
  "aggregation": str}`. Accept only `"mean"` and `"aggregate"`, selecting the
  corresponding stored draw array. Reject non-count predictions.
- `pred.rate_counterfactuals()` returns `{"factual": rate_summary,
  "control": rate_control_summary}`. Its bounds use the `alpha` supplied to
  `predict`, like `catt`; it does not reconstruct draws from summary moments.

No separate lift draw array is needed: aggregate ratio draws minus one give it
exactly. Keep ordinal probability methods guarded against NB.

### 6.3 Numerical and memory requirements

Promote block-level link transformations and count/ratio outputs to host
float64. Evaluate small contrasts using `expm1(f1-f0)` where safe; use stable
log-sum accumulation for aggregate ratios. Do not exponentiate posterior mean
log rates or clip transformed draws to hide overflow. A nonfinite or
nonrepresentable transformed result must raise an informative error rather
than yield a misleading finite ATT; rates and mean counts are mathematically
positive, while differences may be negative or zero.

Reuse [BlockAccumulator](src/longbet/_summary.py) with float64 outputs. All draws
for one cell are present together, so per-cell means, standard deviations, and
quantiles remain exact. Running sums and squares alone cannot recover quantiles.
Accumulate horizon sums separately for `delta`, cell IRRs, and each count arm.

Add a count-specific block budget that accounts for simultaneously live float64
buffers and quantile workspace (target approximately 32 MiB, at least one cell).
Clamp an oversized user `block_size` to this bound in summary mode. Without the
explicit forest cache, working memory scales as `O(block_size*D + S_out*D)` plus
required `O(N*T)` summaries/design and model traces, never `O(N*T*D)`.
`cache_forest_evaluations=True` retains its documented memory opt-in; cached
forest evaluations must be transformed again using the current prediction
exposure. Test two calls with the same design and different exposure.

## 7. Traces and persistence

Extend [LongBetTrace and the loop carry](src/longbet/_loop.py) with optional
`dispersion` traces, allocated, saved, and assembled alongside `beta`, coding,
and forests at the same retained iteration. Shapes are `(L,)` without a chain
axis or `(C,L)` with chains. Save fixed dispersion too. Do not save `omega`, `q`,
or cell residuals at every sweep.

Count acceptance counters must count actual MH proposals, including skipped
iterations; `n_skip` must not turn a last-sweep boolean into an acceptance-rate
estimate. Record accepted/proposed totals per chain, distinguish warm-up and
retained-sampling iterations, and expose them through `model.count_diagnostics_`.
Record PG numerical failure status and synchronize it before `fit` succeeds.
For fixed dispersion there are zero proposals and an unavailable acceptance
rate. Retain `r` draws for diagnostics; diagnose scientific effects as well.

[Scalar NPZ I/O](src/longbet/_io.py) currently has precision-cache and ordinal
schema metadata, **not** the multi-model `format_version` sequence. Add
`COUNT_SCHEMA_VERSION=1` for count archives; do not invent scalar "format 4"
or change [multi archive versions](src/longbet/_multi_io.py).

Count archives must contain:

- Existing config, design, times, dimensions, forest offsets/traces, scaling,
  and `precision_cache_version` metadata.
- `dispersion` aligned with retained parameter draws and `log_exposure_fit`
  shaped `(N,T)` as numeric NPZ arrays; use an optional I/O argument for the
  latter, rather than a large JSON list.
- Metadata `count_schema_version=1`, `count_parameterization="nb2_log_mean"`,
  `count_sampler_semantics="pg_marginal_r_then_omega_v1"`,
  `requires_exposure`, and count diagnostics. Also record the verified PG
  algorithm/version selected in phase 0.

Validate these fields on write and load: schema and parameterization match,
finite positive dispersion, trace/dimension/chain consistency, finite offsets,
boolean exposure flag, unit working variance and `meany=0,sdy=1`. Fixed-dispersion
traces must equal the configured value. Reject missing or corrupted count data;
do not guess a default `r` or reinterpret an earlier parameterization. Preserve
existing validation and loading behavior for non-count archives.

Round trips must reproduce every prediction field and summary with the same
prediction inputs/key, including forecast draws. Loading restores a model for
prediction, as today; it does not promise MCMC continuation or reconstruct the
last augmentation state.

## 8. R front door and shared contract

Update [R/longbet.R](R/longbet.R), [R/predict.longbet.R](R/predict.longbet.R),
[R/att.R](R/att.R), and [R/rehydrate.R](R/rehydrate.R) as needed:

- Append new named options without changing existing positional argument
  meanings. Expose all five dispersion options with Python defaults, plus
  `exposure=NULL` at fit and prediction; add `"negbin"` to scalar outcomes only.
- Validate numeric count labels and exposures before lossy coercion. Forward
  `(N,T)` matrices through reticulate; do not use R's column-major `as.vector`
  as the Python unit-major flattened input. Keep singleton dimensions intact.
- Existing count `get_att`, `get_catt`, and prediction fields mirror Python.
  Add exported `att_irr(object, alpha=0.05, aggregation="mean")`, delegating to
  `py_pred$irr`. Return a data frame with `exposure`, `irr`, `lower`, `upper`,
  `aggregation`, and `n_treated`, one row per horizon. Expose rate summaries and
  dispersion draws in the prediction object without duplicating statistics in R.
- Persist exposure flags and count data through raw NPZ/RDS storage. Test
  fresh-process rehydration of both models and prediction objects.
- Update `NAMESPACE`, relevant `man/*.Rd`, Python docstrings, README examples,
  and both [contract YAML](contract/longbet-api.yaml) and its
  [installed copy](inst/contract/longbet-api.yaml). The YAML copies must remain
  byte-identical; their bidirectional config checks catch undocumented options.

Example of the intended Python usage after all release gates pass:

```python
model = LongBet(outcome="negbin", initial_dispersion=4.0).fit(
    y, x, z, t=t, exposure=person_days,
)
pred = model.predict(x, z, t=t, exposure=person_days, summary_only=True)
count_att = pred.att()
mean_cell_irr = pred.irr()
aggregate_irr = pred.irr(aggregation="aggregate")
```

## 9. Verification and acceptance criteria

### 9.1 Deterministic and fast integration checks

Use independent expected calculations; avoid tests that simply call the same
helper twice. Add focused tests in `tests/test_count_inputs.py`,
`test_count_primitives.py`, `test_count_model.py`, `test_count_predict.py`, and
`test_count_io.py`.

1. **Inputs:** dtype, shape, pre-cast fractions, infinities, zeros, constants,
   missingness, exposure below one, precision limits, fixed dispersion, and
   unsupported wrapper combinations. Flattened/matrix equivalence and fitted
   explicit-exposure omission errors must survive loading.
2. **Algebra:** compare NB log masses and log-r MH ratios with SciPy float64,
   including `r<1`, large counts, zero counts, and candidate-dependent `psi`.
   Check the Gaussian working kernel's log-density differences against the
   expanded PG quadratic with a nonconstant offset. Compare NB and Poisson PMFs
   at a fixed mean for increasing `r`; do not require sample-path equality.
3. **Sweeps:** check residual/working-response invariants at initialization and
   after successive sweeps, one/multiple chains, missing rows and columns,
   fixed/sampled `r`, `sample_alpha`, adaptive coding, random intercept, GP, and
   ridge toggles. Use direct weighted Gaussian linear algebra for GP/coding/
   intercept conditionals and explicit leaf sums for precision-cache checks.
   Missing-only units must draw intercepts from their prior conditional.
4. **Trace alignment:** thinning and batched/unbatched loops save dispersion
   with the same parameters/forests, preserve chain axes, and count all MH
   attempts. Exercise failure propagation through `jit` and chain `vmap`.
5. **Predictions:** use hand-built posterior traces with nonzero `b0`, unequal
   `nu(S)`/`nu(0)`, varying exposure, nonunit `alpha`, and nonzero intercepts.
   At fixed draws, multiplying exposure by `k>0` multiplies counts/ATT by `k`
   but leaves rates and IRRs unchanged; changing only `r` leaves means unchanged.
   Test factual versus treated arms on controls, heterogeneous mean-versus-
   aggregate IRRs, no treated cells, sparse horizons, and projected horizons.
6. **Summaries/memory:** full versus summary parity, singleton axes, ragged
   blocks, at least three block sizes, exact quantiles, and cached predictions
   with changed exposure. Instrument allocations/transforms with cache disabled
   to show that no full panel-by-draw buffer appears; output-shape checks alone
   do not establish bounded memory.
7. **Persistence:** scalar and R fresh-process round trips; deliberately missing
   arrays, bad schemas, nonpositive dispersion, wrong draw axes, wrong offsets,
   inconsistent fixed dispersion, and nonunit working variance must fail.
8. **Regression:** run the existing contract, scalar, ordinal, chain-axis,
   precision-cache, prediction, and multi scalar-equivalence tests. Adding an
   outcome must not reroute NB through a probit or joint-SUR path.

### 9.2 Distribution and posterior tests (mark expensive tests `slow`)

For independent PG draws with `h>0` and `c=abs(tilt)`, the moment and Laplace
checks are

$$
E[\omega]=\frac{h}{2c}\tanh(c/2),\quad
\operatorname{Var}(\omega)=\frac{h}{2c^3}\frac{\sinh(c)-c}{\cosh(c)+1},
$$

with zero-tilt limits `h/4` and `h/24`, and

$$
E[e^{-t\omega}]=\left[\frac{\cosh(c/2)}
{\cosh(\sqrt{c^2+2t}/2)}\right]^h,\quad t\ge0.
$$

Use stable evaluations/series near zero. These distribution identities follow
from the [PG Laplace transform](https://arxiv.org/html/1405.0506).

- Cover `h={0.1,0.5,1,1.7,5,20,100}` and
  `tilt={0,1e-6,1.5,-2,20,100}`, plus every algorithm-switch boundary.
  Compare means, variances, and several Laplace values using Monte Carlo
  standard errors (variance errors need fourth moments or independent batches).
  Use a stated simultaneous error budget, such as Bonferroni-adjusted 99%
  intervals; fixed three-standard-error thresholds across a large grid are
  prone to false failures. Include a separately evaluated CDF/quantile oracle
  where feasible. Moments alone do not establish the distribution or exactness.
- With `eta` fixed, compare a dispersion chain's distribution with normalized
  one-dimensional numerical quadrature under the specified prior. Recovery
  means agreement with this finite-data posterior, not equality to generating r.
- For a tiny intercept-only NB model with a proper normal log-rate prior and
  nonconstant exposure, compare joint `(f,r)` samples from the complete
  partially collapsed transition with independent two-dimensional quadrature.
  Check marginals and covariance using effective Monte Carlo errors. This is
  the required oracle for the dispersion/augmentation ordering.

### 9.3 Statistical benchmark and evidence report

Create `benchmarks/count_validation.py` and `benchmarks/count_validation_report.md`.
The main reproducible panel uses `N=400`, `T=12`, `t=1..12`, equal cohorts adopting
at periods `4,6,8,never`, and independent standard-normal `X1,X2`:

```text
baseline_i = 0.5*sin(X1_i) + 0.3*X2_i
unit_intercept_i ~ Normal(0, 0.4**2)
e_it ~ Uniform(0.8, 1.5)
beta[0] = 0
beta[s] = 1 + 0.8*(beta[s-1] - 1) for s >= 1  # deterministic trajectory
log_rate0 = baseline_i + unit_intercept_i
log_rate1 = log_rate0 + beta[S_it]*(0.3 + 0.2*X1_i)
lambda_a = e_it*exp(log_rate_a)
r_true = 3
Y ~ NB(r_true, success_probability=r_true/(r_true+lambda_factual))
```

Truth is obtained by applying section 6's formulas to the finite panel, with
exactly the same horizon membership as prediction. Add an intercept-free base
case, a 10% independently missing-outcome case, a high-zero/low-rate case, and
a high-count case. Randomize cohort assignment independently of covariates for
the primary benchmark so effect recovery is not confounded by an unspecified
assignment mechanism. Do not label the deterministic trajectory a random AR(1)
draw.

Start with ten fixed seeds `20260912..20260921`, fit seeds offset by `1000`, four
chains, `num_burnin=2000`, `num_sweeps=1000`, and `n_skip=1`, with other options
recorded explicitly. Report count ATT and both IRRs: bias/RMSE, interval widths,
coverage per horizon, treated counts, R-hat, bulk/tail ESS, and dispersion MH
acceptance. Use R-hat <= 1.04 and bulk ESS >= 250 as initial diagnostic goals;
report failures and investigate them before claiming successful validation.
Record seed-level and horizon-level results rather than hiding disagreement
in one average. Changes to run lengths/settings must be explicit in the report.

Ten datasets are a smoke benchmark, not evidence for a narrow 91–98% coverage
requirement: each horizon's observed coverage moves in 10-point increments.
Report binomial uncertainty across independent datasets for each horizon;
correlated horizons are not independent replications. Claims of calibrated
coverage require a larger, predeclared replication study, with uncertainty
reported. A `log(y+1)` comparison is optional and descriptive; NB is not required
to outperform another estimator on every seed.

Artifacts must include config, seeds, package/device versions, source digest,
DGP counts/exposures/true counterfactuals, draws needed to reproduce estimands,
and metrics. Synchronize JAX with `block_until_ready` before timing, separate
compilation from sampling, and report peak memory. Stale artifacts do not count
as evidence for changed sampler code. Document prior sensitivity of dispersion
and log-rate scales, especially in weakly informed/all-zero panels.

## 10. Implementation phases and definition of done

| Phase | Files / deliverables | Gate before advancing |
|---|---|---|
| 0. PG feasibility | `_count.py` prototype; `benchmarks/count_sampler_review.md` | Algorithm and full real-shape domain justified; independent PG tests pass; exactness/numerical limitations explicit. |
| 1. Inputs and primitives | `_config.py`, `_count.py`, `_model.py` preparation; contract entries | NB/log-r algebra and validation tests pass; configuration is documented. |
| 2. Scalar sampler | `_state.py`, `_step.py`, `_loop.py`, cache/scale paths as needed | Weighted conditional oracles, residual invariants, joint posterior oracle, chain/thinning checks, and legacy regressions pass. |
| 3. Prediction and archives | `_model.py`, `_summary.py`, `_io.py` | Estimand algebra, memory instrumentation, full/summary parity, corruption checks and round trips pass. |
| 4. R and public docs | R wrappers, `NAMESPACE`, `man/`, README, both contracts | Python/R parity and fresh-process RDS rehydration pass; `R CMD check --no-manual` completes. |
| 5. Statistical validation | Benchmark, report, distribution/validation tests | All prescribed runs/artifacts inspected; numerical and diagnostic failures resolved or explicitly leave the release incomplete. |

Suggested checks once the corresponding tests exist:

```bash
.venv/bin/python -m pytest tests/test_count_inputs.py tests/test_count_primitives.py tests/test_count_model.py tests/test_count_predict.py tests/test_count_io.py tests/test_contract.py -m 'not slow' -q
.venv/bin/python -m pytest tests/test_step.py tests/test_chain_axis.py tests/test_forest_cache.py tests/test_predict_options.py tests/test_ordinal_model.py tests/test_ordinal_predict.py tests/test_multi_scalar_equivalence.py -m 'not slow' -q
.venv/bin/python -m pytest tests/test_count_primitives.py tests/test_count_model.py -m slow -q
```

Record actual completed commands, results, skipped checks, and unresolved gates
in `benchmarks/count_implementation_status.md`. A fast test pass or existing
output file is not proof that the statistical benchmark ran successfully.
Completion means all first-release Python/R functionality and validation above
are delivered; section 11 is future work and does not enter that completion gate.

## 11. Deferred hurdle and joint-count work

The earlier proposal's binary-probit plus log-normal positive component models
a semicontinuous positive outcome, not an integer-valued positive count. A true
count hurdle requires a positive-count distribution (for example zero-truncated
NB), its normalization, and a separate inference design; the ordinary NB PG
kernel cannot be reused unchanged after truncation.

Also, correlated probit/log-normal errors invalidate the independent-product
mean in the earlier draft. To see the issue, suppose
`D=1[eta_D+epsilon_D>0]`, `log(Y*)=eta_Y+epsilon_Y`, jointly normal errors,
`Var(epsilon_D)=1`, `Var(epsilon_Y)=v_Y`, and
`Cov(epsilon_D,epsilon_Y)=c`. Exponential tilting of that bivariate normal gives

$$
E[D Y^*]=\exp(\eta_Y+v_Y/2)\,\Phi(\eta_D+c).
$$

The product `Phi(eta_D)*exp(eta_Y+v_Y/2)` follows only when `c=0`. This is a
derivation for that selection model, not a claim that existing `LongBetMulti`
implements a validated count hurdle. Selection, positive-outcome support,
conditional versus structural variances, identifiability, and missing-data
semantics need their own specification and tests.

Do not add NB to multi-outcome discrete-first probit validation, or multiply SUR
precision by PG precision, to simulate joint support. A joint count likelihood
requires an explicit dependence model. Keep the current allowlists rejecting
NB in `LongBetMulti` and encouragement wrappers until that work is designed.
