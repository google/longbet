# Randomized encouragement with a continuous treatment intensity

`LongBetEncourage(treatment="intensity")` extends the adoption-clock
encouragement model ([encouragement guide](encouragement-guide.md)) to an
encouraged behaviour that is **continuous and non-absorbing**: hours of
exercise, study or screen time, quarterly training-months. Such a dose rises,
varies fractionally and switches off, so there is no adoption time and no
exposure clock. The model replaces the clock with a **distributed-lag kernel on
the calendar clock** and keeps everything else the encouragement layer already
had: a design-based reference that needs no model, a Bayesian two-stage model
whose offer effects follow by composition, unit-level effects and the decision
layer, archives, and the R front door.

## Data

A single randomized offer wave: `z[i, t] = Z_i * 1{t >= t0}` is binary and
absorbing with a permanent control arm (as in the binary model). The dose
`d[i, t] >= 0` is a complete finite panel; the outcome `y[i, t]` is continuous.
Periods are consecutive and unit spaced (the kernel is indexed by lag). Pass
`strata` (sites) to `fit()` when randomization was blocked; the reference and
the dip screen then analyse the experiment as it was randomized.

## The model

Stage 1 is the dose on the offer clock and Stage 2 the outcome on the dose
history:

```math
D_{it} = \mu_D(X_i) + \eta^D_{it} + \beta^D_{S^Z_{it}}\,\nu_D(X_i)\,Z_{it} + \gamma^D_i + v^\perp_{it},
\qquad
Y_{it} = \mu_Y(X_i) + \eta^Y_{it} + \nu_Y(X_i)\sum_{\ell=0}^{L}\kappa_\ell\,D_{i,t-\ell} + \mathrm{CF}_{it} + \gamma^Y_i + \epsilon_{it},
```

with $\kappa = (\kappa_0,\dots,\kappa_L) \sim \mathcal N(0, s^2 K_L)$ a
Gaussian process over the lag index, both moderators unit-anchored
($N^{-1}\sum_i \nu(X_i) = 1$), unit intercepts and AR(1) errors in both
equations, and $L = T - t_0$ by default (one lag per post-offer period; pre-panel
lags are zero, or set `demean_pre=True` to subtract each unit's pre-offer mean
dose first). $\kappa_\ell$ is the effect on $Y_{it}$ of one more unit of dose
received $\ell$ periods ago; the GP prior allows sign changes, so a negative
contemporaneous lock-in ($\kappa_0 < 0$) followed by positive returns
($\kappa_{\ell\ge1} > 0$) is representable. The binary model is the special
case $\kappa_\ell = \tau(\ell+1) - \tau(\ell)$ of an absorbing dose (Abel
summation).

Offer effects are composed draw by draw, pairing draw $b$ of Stage 1 (stream
`key_d`) with draw $b$ of Stage 2 (stream `key_y`):

```math
\mathrm{CITT}_Y^{(b)}(x,t) = \nu_Y^{(b)}(x)\sum_{\ell=0}^{t-t_0}\kappa^{(b)}_\ell\,\mathrm{CITT}_D^{(b)}(x,t-\ell).
```

Pre-offer doses do not respond to the offer (no anticipation) and the
control-function terms have mean zero in both arms given $(X, Z, t)$, so
neither enters the offer contrast.

## The three identification tiers

`identification` chooses what Stage 2 learns the kernel from. All three tiers
share Stage 1 and the Stage-2 engine; they differ in the lag matrix and in the
presence of the control function.

| `identification` | Tier | Stage-2 lag regressor | Extra conjugate columns | Robust to | Variance |
|---|---|---|---|---|---|
| `"deconvolution"` | T0 | $(Z_i - \hat e_i)\,\widehat{\mathrm{CITT}}_D(X_i,t-\ell)\,\mathbb{1}\{t-\ell \ge t_0\}$, the propensity-centered Stage-1 offer-induced dose increment | $[\widehat D^{(0)}_{it},\ \sum_{\ell=1}^{L}\widehat D^{(0)}_{i,t-\ell}]\otimes[1,\tilde X_i]$ (unconstrained organic-dose nuisance controls; no CF on $v_{it}$) | any within-arm dose endogeneity (static, time-varying, dynamic feedback $Y_{t-1}\to D_t$) and classical measurement error in the dose | highest; the GP prior regularizes the Toeplitz deconvolution |
| `"control_function"` (default) | T1 | observed $D_{i,t-\ell}$ | $[\hat v_{it},\ \hat v^{\mathrm{lag}}_{it},\ \bar{\hat v}_i]\otimes[1,\tilde X_i]$ | static selection and linear contemporaneous / cumulative-lag selection | intermediate: $\kappa_0$ and the level of $\kappa_{1:L}$ anchored by the offer, the lag shape sharpened by within-arm variation |
| `"unconfounded"` | T2 | observed $D_{i,t-\ell}$ | none | neither | lowest |

The T1 control function is built from the Stage-1 posterior-mean residuals
$\hat v_{it} = D_{it} - \hat D_{it}$: the contemporaneous residual, the
**cumulative lagged residual** $\hat v^{\mathrm{lag}}_{it} = \sum_{\ell=1}^{L}\hat v_{i,t-\ell}$,
and the unit's post-offer mean residual, each interacted with the standardized
covariates. All three terms are needed. Without $\hat v^{\mathrm{lag}}$ the
lags $1..L$ stay confounded (even static selection enters them through an
active window that grows with calendar time) and the free offer-driven column
at lag 0 absorbs the whole gap, producing a spurious contemporaneous effect and
wrong-signed later lags; `tests/test_intensity_reference.py` replicates the
mechanism. The columns are fixed during Stage-2 sampling and their coefficients
are drawn in the same conjugate block as the trend, the intercepts and the
kernel, which keeps the two MCMC streams independent and the composition exact
for the product of the two posteriors. (A draw-by-draw "cut" posterior would
need one Stage-2 chain per Stage-1 draw.)

### Which tier when

* Start with **T0 as the benchmark**: it is the only tier whose kernel does not
  depend on how units choose their dose within an arm. Report it always.
* Use **T1** when the dose is chosen with a time-invariant or linearly
  selected component and the first stage is informative: it is several times
  more precise at lags $\ge 1$ and agrees with T0 when its working assumption
  holds. Its linear control function is exact under linear conditional-mean
  dependence of the outcome shock on $(\hat v_{it}, \hat v^{\mathrm{lag}}_{it}, \bar{\hat v}_i)$
  and a working projection under hurdle/zero-mass censoring; the discrepancy
  statistic below audits it.
* Use **T2** only when within-arm dose variation is credibly unconfounded given
  covariates and the unit intercept (a designed dose, a mechanical schedule).
  It is the most precise and the most fragile.
* Under **dynamic feedback** ($Y_{i,t-1}$ drives $D_{it}$) only T0 is
  consistent; under **measurement error** in the dose only T0 is unbiased.

## The design-based reference

`intensity_reference(y, d, z, t=None, x=None, strata=None, lag_horizon=None,
demean_pre=False, alpha=0.05)` needs no sampler and is what `predict()` attaches
to every fit. It reports:

* `itt`: stratum-weighted $\widehat{\mathrm{ITT}}_D(t)$, $\widehat{\mathrm{ITT}}_Y(t)$
  for every period (pre-offer rows are placebo contrasts) with Neyman standard
  errors, the per-period ratio `omega` with its Fieller/Anderson–Rubin set, and
  for post-offer periods the **cumulative return per unit of cumulative dose**
  `omega_cum(t) = sum_{r<=t} ITT_Y(r) / sum_{r<=t} ITT_D(r)` with an exact
  Fieller/AR set (`omega_cum_set_type`, `omega_cum_lower_1`, ...). The
  per-period ratio divides by a first stage that may have returned to zero;
  the cumulative return does not.
* `covariance`: the joint Neyman covariance of $(\mathrm{ITT}_Y, \mathrm{ITT}_D)$
  over post-offer periods.
* `kernel`: the unregularized Toeplitz deconvolution
  $\hat\kappa = T_{\widehat{\mathrm{ITT}}_D}^{-1}\widehat{\mathrm{ITT}}_Y$ with
  delta-method standard errors and the cumulative kernel $\bar\kappa(h)$;
  `condition_number` is $\operatorname{cond}(T_d)$: near 1–3 for a
  front-loaded compliance pulse, exponentially large for a gradual ramp
  ($d_1 > d_0$), in which case the Bayesian tiers' GP prior does the
  regularizing and the reference's lag-level SEs are honest about it.
* `anderson_rubin(kappa)`: the joint $\chi^2$ test of a hypothesized kernel
  (one degree of freedom per post-offer period), valid under a weak first
  stage; `anderson_rubin_projection()` projects the region onto each lag by
  Monte Carlo (an inner approximation; `hits_envelope` flags directions in
  which the set should be regarded as unbounded).
* `ecdf`, `fosd`: arm-wise ECDFs of cumulative post-offer dose on the pooled
  support, their difference $F_0 - F_1$ (positive everywhere under monotone
  compliance), the extensive-margin shift at zero dose and a one-sided Smirnov
  check of first-order stochastic dominance.
* `dip`: the **dose-history-conditioned pre-treatment dip screen**: the
  coefficient of $\Delta Y_{i,t-1}$ in a regression of $\Delta D_{it}$ on it,
  the dose lags $D_{i,t-1},\dots,D_{i,t-L-2}$, covariates, the offer and period
  fixed effects, with a unit-clustered standard error. Conditioning on the dose
  history blocks the mechanical path from past doses to the lagged outcome
  innovation, so the screen has nominal size when past doses have real effects
  and power when units react to their own outcome innovations (dynamic
  feedback), the case in which T1 and T2 fail.
* `discrepancy(kappa_draws)`: given posterior draws of a kernel,
  $z_{\mathrm{tier}-\mathrm{T0}}(h)$ compares the cumulative kernels, dividing by
  $\sqrt{\mathrm{SE}_{T0}(h)^2 + \mathrm{SD}_{\mathrm{post}}(h)^2}$ (the two
  estimators are positively dependent, so this is alarm-averse).

What the reference identifies: under randomization, exclusion and no
anticipation, linear calendar-homogeneous lag effects and time-invariant
conditional complier weights, the population deconvolution recovers the
**compliance-weighted** kernel
$\kappa^{\mathrm{acr}}_\ell = \mathbb E[\nu_D\nu_Y]/\mathbb E[\nu_D]\cdot\kappa_\ell$,
which equals the unweighted population kernel when compliance and benefit are
uncorrelated across covariates. `predict().kernel()` therefore reports both the
population kernel (`kappa`) and the compliance-weighted kernel (`kappa_acr`)
and compares the latter with the reference.

## Diagnostics and what to do with them

`pred.diagnostics()` collects, per group: `condition_number`, `z_dip`,
`fosd_pvalue` / `fosd_sup_violation`, and `z_discrepancy` (the headline
$z_{\mathrm{T1-T0}}$ or $z_{\mathrm{T2-T0}}$ at the last horizon; the full
trajectory is `pred.discrepancy()`).

* **$|z_{\mathrm{T2-T0}}|$ large, $|z_{\mathrm{T1-T0}}|$ small**: within-arm
  dose selection that the control function absorbs; report T1 with T0.
* **Both large**: selection the linear control function does not capture
  (dynamic feedback, measurement error, nonlinear selection); report T0.
* **$z_{\mathrm{dip}}$ large**: units react to their own outcome innovations;
  only T0 is consistent.
* **FOSD violated**: defiers or reactance; the Wald/ACR summaries carry signed
  weights and the model's monotone-compliance reading does not hold.
* **Large condition number**: the first stage ramps up slowly; prefer the
  cumulative return `omega_cum` and the cumulative kernel over lag-level
  estimates, and expect the GP prior to matter in the Bayesian tiers.

`pred.cumulative_discrepancy()` and `fit.bootstrap_comparison()` compare the
model's ITTs with the reference exactly as in the binary model.

## Example

```python
import numpy as np
from longbet import LongBetEncourage, intensity_reference

# y, d, z: (N, T) outcome, dose and offer panels; x: (N, P); site: (N,) strata.
ref = intensity_reference(y, d, z, x=x, strata=site)
print(ref.kernel)                       # Toeplitz kernel with delta-method SEs
print(ref.itt[["period", "itt_d", "itt_y", "omega_cum", "omega_cum_lower_1", "omega_cum_upper_1"]])
print(ref.condition_number, ref.dip["z_stat"], ref.fosd["pvalue"])

fit = LongBetEncourage(treatment="intensity", identification="control_function",
                       num_chains=4, num_burnin=2000, num_sweeps=1000)
fit.fit(y, d, z, x, t=t, strata=site)

pred = fit.predict(groups=segment)
print(pred.effects())                   # itt_y, itt_d, wald, omega_cum by group and horizon
print(pred.kernel())                    # posterior kappa, kappa_acr beside the reference kernel
print(pred.discrepancy())               # z_{T1-T0} by horizon
print(pred.diagnostics())
print(pred.cumulative_return())         # omega_cum with the exact Fieller/AR set

cond = fit.predict_conditional(x=new_cohort)
cond.citt_y.mean, cond.citt_d.mean, cond.kernel.mean, cond.omega_cum.median
cond.breakeven_probability(cost=20.0)
cond.knapsack_policy(cost=20.0, capacity=80, hurdle=0.5)

fit.save("fit.npz"); again = LongBetEncourage.load("fit.npz")
```

For a sensitivity analysis fit the three tiers with the same `key` or
`random_seed`; they share Stage 1 exactly. `stage2_backend="linear"` replaces
Stage 2 with a closed-form two-way fixed-effects regression (no covariate
moderation, independent errors) for smoke tests and quick checks; it is not the
model of the paper.

In R: `longbet_encourage(y, d, z, x, treatment = "intensity", identification =
"control_function", strata = site, config = list(...))`, then `predict()`
(adds `kernel`, `discrepancy`, `design_diagnostics`, `cumulative_return`),
`predict_conditional()` (adds `kernel`, `kernel_cum`, `omega_cum`) and
`intensity_reference()`.

## What is fixed and what is free

As in the binary model, the components are fixed: horseshoe linear term in
both moderators, learned trajectory and kernel scales, interweaving, AR(1)
errors, unit intercepts, unit-anchored moderators. The lag kernel's prior is
configurable through `LongBetConfig` (`lag_sig_knl`, `lag_lambda_knl`,
`lag_kernel_type`, `lag_sigma_m`, `lag_gp_constant_mean`, `lag_learn_scale`,
`lag_scale_prior`) and the control-function prior through
`extra_cols_prior_scale`. `lag_horizon` shortens the kernel; `demean_pre`
handles behaviours with a pre-study steady state. A binary outcome is not
supported with `treatment="intensity"`.

Archives (`save`/`load`) nest both stages, record `treatment`,
`identification`, `lag_horizon`, `demean_pre`, the Stage-2 backend, the
control-function column names, `component_streams: independent` and the
engine's version and git commit. Posterior calibration status is
`not_established`, as for every encouragement fit.
