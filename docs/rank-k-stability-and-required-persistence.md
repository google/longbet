# Is the effect stable? In-sample persistence diagnostics and required-persistence decisions

A short longitudinal experiment identifies the cumulative effect over its observed window far better than the
*slope* of that effect, yet any return-on-investment question that reaches past the end of the panel hinges on the
slope. Version 0.8.0 adds two layers that keep those two questions apart:

* **In-sample diagnostics** (`longbet._stability`): how much the observed window can say about persistence, and how
  much it cannot, without refitting anything.
* **A decision layer** (`longbet._rank_k_plots`): post-panel questions are inverted into the *required persistence*
  `X` (how long the effect must last for an investment to pay back), reported next to explicitly labelled scenario
  brackets and a placebo null benchmark.

> **Package rule.** Following `att_stability` ("report the numbers, withhold an uncalibrated verdict"), no function
> certifies a multi-year extrapolation. `TrajectoryStabilityReport.verdict` is always `None`; the retention rule is
> the analyst's pre-registered rule, evaluated literally and under robustness variants.

## What a short panel can and cannot say

| Periods | Status | What the numbers rest on |
|---|---|---|
| Observed cells of the treated-vs-control contrast | Identified by design | Randomization (conditional on strata) |
| Unobserved cells inside the horizon of a masked fit | Model-based | The fit's covariates; check `rank_k_transport_audit` and its `unobserved_effect_share` |
| Periods after the panel ends | **Not identified** | A persistence scenario, the plateau / pulse split, a survival curve and constant monetization |

Three failure modes make naive slope extrapolation look more precise than it is:

1. **Level versus slope signal-to-noise.** A panel that is well powered for the cumulative effect can have a
   minimum detectable wear-out larger than the effect itself. `wearout_ancova` prints the minimum detectable
   difference (MDD) next to every estimate, so "no detectable wear-out" is never read as "no wear-out".
2. **Boundary pile-up.** Positivity constraints (`catt_floor = 0`, half-normal priors, outcomes with a hard cap) push
   retention draws onto `rho = 1`; a slope fitter then reads truncation as precision.
   `estimate_post_pulse_retention` reports the share of draws on each bound.
3. **Correlated "replicates".** Retention estimates of several outcomes measured on the same units are not
   independent studies. Their Cochran `Q` falls in the lower tail and the pooled standard error is understated (by
   `sqrt(1 + (M - 1) r)` for `M` outcomes with equal variances, a common correlation `r` and `tau^2-hat = 0`, or
   `sqrt(1 + (M - 1) r v / (v + tau^2))` when `tau^2 > 0`).
   `pool_multi_outcome_retention` re-evaluates the pooled retention under correlation-adjusted, Q-moment,
   leave-one-outcome-out and equal-weight variants.

## Workflow

1. **Trajectory draws.** From a two-arm fit, `event_time_trajectory(res, weighting="unit")` gives the `(S, H)`
   treated-average trajectory. From a nested three-arm design (C control, B base treatment, A base plus an
   increment), `reconcile_nested_3arm_draws(res_bc, res_ab, res_ac)` gives coherent `BC`
   (persistent plateau), `AB` (increment, often a fading pulse) and `AC = BC + AB` draws. Every function below
   accepts the reconciliation result directly.
2. **In-sample stability** (`longbet._stability`):
   * `wearout_ancova`: design-based early-versus-late ANCOVA on the balanced subpanel (strata fixed effects, HC1).
   * `estimate_post_pulse_retention`: per-draw monthly retention `rho` over a post-pulse window, weighted least
     squares on levels (robust to exact zeros), with bound pile-up shares.
   * `unit_slope_correlation`: within-arm correlation across outcomes of per-unit slopes, the proxy for the
     correlation of the outcomes' retention estimators.
   * `pool_multi_outcome_retention`: DerSimonian-Laird pooling of `lambda = -ln(rho)`, draw-preserving Morris
     shrinkage of a focal outcome, robustness variants and the break-even correlation.
   * `RetentionRuleSpec` and `evaluate_retention_rule`: the pre-registered gate, with a literal and a robust reading.
   * `diagnose_trajectory_stability`: all of the above in one `TrajectoryStabilityReport`.
3. **Size a hindcast before running it.** `pseudo_horizon_separation` computes, from existing draws, how far apart
   the no-decay and decay scenarios are over a pseudo-horizon, the power to tell them apart and how much shorter the
   pseudo-horizon lever is than the real tail. A leak-free refit is worth its cost only when that power is high.
4. **Beyond the panel: scenarios, not forecasts.** `extrapolate_rank_k_decay(..., pulse_decay_tau=..., survival=...)`
   extends the plateau with `rho^h`, lets the pulse fade with `exp(-h / tau)` and weights both by a survival curve
   (`weibull_survival_weights`); `discount_weights` gives present-value weights.
5. **Decision inversion.** `compute_required_persistence` reports the required persistence `X`, the break-even
   retention, payback CDFs and cumulative value under Bound U (`rho = 1`), Bound L (default `rho = 0.85`) and the
   no-persistence floor (`rho = 0`), optional extra regimes (for instance the pooled posterior retention, labelled
   *diagnostic* unless the rule passed) and a placebo gate. `plot_required_persistence_and_payback` draws the
   three-panel figure.

## Example (synthetic draws, runs in seconds)

```python
import numpy as np
import longbet

rng = np.random.default_rng(2026)
S, L, H = 400, 12, 60
months = np.arange(1, L + 1)

# Posterior draws of a nested three-arm design: a persistent plateau (B vs C) with a slow drift and noise,
# and an increment (A vs B) that fades.
level = np.abs(rng.normal(0.20, 0.08, size=(S, 1)))
drift = rng.normal(0.0, 0.01, size=(S, 1))
plateau = np.clip(level * (1.0 + drift * (months - 6)) + rng.normal(0.0, 0.01, size=(S, L)), 0.0, None)
pulse = np.abs(rng.normal(0.15, 0.05, size=(S, 1))) * np.exp(-(months - 1) / 3.0)
traj = {"BC": plateau, "AB": pulse, "AC": plateau + pulse}

# Three secondary outcomes on the same units, sharing the plateau's slope.
outcomes = {
    name: np.clip(plateau * scale + rng.normal(0.0, 0.01, size=(S, L)), 0.0, None)
    for name, scale in (("minutes", 3.0), ("active_days", 1.5), ("sessions", 2.0))
}

# 1. In-sample stability: retention per outcome, pooled retention and the pre-registered rule.
report = longbet.diagnose_trajectory_stability(
    {**outcomes, "value": traj},           # any (S, L) draws, fits, reconciliations or {contrast: draws}
    window=(7, 12),                        # post-pulse window of the retention fits
    focal="value",                         # shrunk toward the pooled outcomes and excluded from the pool
    correlation=0.5,                       # or a SlopeCorrelationResult from unit_slope_correlation
    rule=longbet.RetentionRuleSpec(p10_threshold=0.90),
    degenerate_policy="resample",          # keeps pooled_rho_draws finite for extra_regimes when a window is all zero
)
print(report.summary_df)
print(report.rule.reasons, report.headline_regimes)

# 2. Decision inversion under explicit scenarios.
result = longbet.compute_required_persistence(
    traj,
    cost_grid=(2.0, 4.0),
    alpha_grid=(0.5, 1.0),                 # incrementality of the alpha-scaled share of the effect
    horizon=H,
    pulse_decay_tau=3.0,                   # AB fades; BC persists at rho (U: 1.0, L: 0.85)
    survival=longbet.weibull_survival_weights(H, L, median_periods=60.0),
    discount=longbet.discount_weights(H, annual_rate=0.10),
    # Pooled posterior retention of the focal outcome, paired with the draws by index (literal pooling; the draws at
    # the supplied correlation are in report.pooled.pooled_rho_draws_adjusted). Pass degenerate_policy="resample"
    # above when exact-zero post-pulse windows are possible (e.g. catt_floor = 0), since extra_regimes must be finite.
    extra_regimes={"P": report.pooled.pooled_rho_draws["value"]},
    headline_regimes=report.headline_regimes,              # P stays diagnostic unless the rule passed
)
print(result.required_persistence_df)
print(result.decision_df.head())
```

To draw the figure (matplotlib required):

```python
axes = longbet.plot_required_persistence_and_payback(result, cost_highlight=4.0, alpha_highlight=1.0, unit="$")
```

With real fits, replace `traj` by `longbet.reconcile_nested_3arm_draws(res_bc, res_ab, res_ac)`
(three arms) or `longbet.event_time_trajectory(res, weighting="unit")` (two arms, no pulse split), and replace the
secondary outcomes by their own trajectories. For the model-free check, pass the outcome by event time to
`longbet.wearout_ancova(y_event, treated, baseline, strata=strata, early=(1, 6), late=(7, 12))` and hand the result
to `diagnose_trajectory_stability(..., wearout=...)`.

## Definitions

**Tail model.** With the split on, for `h = 1..H-L` and survival weights `sigma_h = S(L+h) / S(L)`:

$$\tau_{BC}(L+h) = \tau_{BC}(L)\,\rho^h\,\sigma_h,\qquad \tau_{AB}(L+h) = \tau_{AB}(L)\,e^{-h/\tau_p}\,\sigma_h,\qquad
\tau_{AC}(L+h) = \tau_{BC}(L+h) + \tau_{AB}(L+h).$$

Without the split every arm follows `tau(L) rho^h sigma_h`. Observed periods are never altered.

**Value.** `v(t) = tau(t) m_t(alpha) delta_t` with `m_t(alpha) = alpha s_t + (1 - s_t)`, where `s_t` (`alpha_share`)
is the share of the effect whose incrementality is uncertain and `delta_t` the discount weight. `alpha_share = 1`
reproduces `roi_break_even_table` when `survival` and `discount` are omitted. The payback month is the first `m` with
cumulative value `V(m) > C`.

**Required persistence.** Along Bound U (`rho = 1` on the plateau, the pulse fading), `X = 0` if `V(L) > C`,
`X = inf` if `V(H) <= C`, and otherwise, with first payback month `m*`,

$$X = (m^* - 1 - L) + \frac{C - V(m^*-1)}{v(m^*)}.$$

When every arm is non-negative, `Pr(X <= k) = Pr(payback under Bound U <= L + k)` for every integer `k <= H - L`
(checkpoints `k > H - L` are reported as `NaN` with a warning), and the payback CDFs are ordered U >= L >= `rho = 0`;
both are asserted (otherwise `cdf_identity_checked` is `False` and a warning is recorded).

**Break-even retention.** The smallest `rho` in `[0, 1]` with `V_rho(H) > C`: `0` when the observed window plus the
fading pulse already pay back, `inf` when even `rho = 1` does not.

**Placebo gate.** With `placebo_draws` (a placebo test run through the same machinery), a cell is
`decision_immune` when the placebo's `Pr(ROI > 0)` stays below 0.5. A decision should not rest on cells that are
not immune.

**Pooling.** See the docstring of `pool_multi_outcome_retention` for the DerSimonian-Laird, correlation-adjusted and
Morris formulas.

## Reporting checklist

* Report the wear-out estimate with its MDD, and say when the test cannot detect even a full wear-out.
* Report the bound pile-up shares of every retention fit.
* Report the literal and the robust reading of the pre-registered rule. When the robust reading fails, the headline
  shows Bounds U and L only and the posterior retention stays a labelled diagnostic.
* Report `X` with Bounds U and L, the `rho = 0` floor and the placebo benchmark, and state the survival curve, the
  discount rate and `alpha_share`.
* Do not present the hierarchical `rho ~ Beta` band of `extrapolate_rank_k_decay` as a credible interval: it is an
  assumed scenario band. Since 0.9.0, `rho_sd` defaults to `None` (no band); pass a positive value explicitly to
  compute the assumed scenario band.
