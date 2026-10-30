# Masked event-time Rank-K fits: population effects under administrative censoring

*Added in longbet-jax 0.7.0.*

## When to use this

You ran a (randomized) rollout, every unit adopted at the same moment *in its
own clock*, and the panel was cut on a fixed calendar date. Units that adopted
early are observed for many post-adoption periods, late adopters for few: unit
`i` is observed on the event-time prefix `s = 1, ..., L_i` of a horizon
`H = max_i L_i`, and nothing is observed beyond `L_i`. You want the effect of a
decision horizon `h` (say 12 months) **for the whole population**, not for the
subset that happened to be observed through `h`, and you do not want to impute
the unobserved cells.

`fit_rank_k_longbet` now takes such a panel directly: set the unobserved cells
of `y` to `NaN` (or pass `mask=`), and every unit stays in the fit with the
cells it has.

```python
import numpy as np
from longbet import (
    RankKLongBetConfig, fit_rank_k_longbet,
    cumulative_effect, event_time_trajectory, rank_k_transport_audit,
)

# y: (N, T) in event time -- t0 pre-adoption columns, then s = 1..H; NaN beyond L_i
# z: (N, T) with z[i, t0:] = 1 for treated units; x: (N, P) pre-treatment covariates
cfg = RankKLongBetConfig(
    rank_k=2,                 # level + saturation ramp (K=3 adds a linear tilt; report both)
    att_weighting="unit",     # every treated unit counts once: the population estimand
    mode_coupling="propagated",  # higher-mode posterior draws propagated into Mode 1 (default)
    n_forecast_draws=8, draws_per_component=25, num_chains=2,  # 400 draws
)
res = fit_rank_k_longbet(y, z, x, cfg=cfg)          # or features=..., horizon=12, mask=...

ate_12 = cumulative_effect(res, 12)                # (S,) draws of (1/N_trt) sum_i sum_{s<=12} tau_i(s)
traj = event_time_trajectory(res)                  # (S, H) unit-weighted tau-bar(s)
audit = rank_k_transport_audit(
    res, strata=np.where(res.post_window_length >= 12, "L>=12", np.where(res.post_window_length <= 7, "L<=7", None)),
    overlap_horizon=5, cumulative_horizon=12,
)
print(audit.summary())
```

## Estimand

With `tau_i(s) = E[Y_i(1, s) - Y_i(0, s) | x_i]` the cell effect at event time
`s`, the `h`-period population effect is

    ATE_h = (1 / N_trt) * sum_i sum_{s <= h} tau_i(s),

every treated unit counting once (`att_weighting="unit"`). `res.att_draws` is
the per-period level averaged over units (`ATE_H / H` under unit weighting),
`cumulative_effect(res, h)` the sum over the first `h` event times, and
`event_time_trajectory(res)` the mean profile `tau-bar(s)`.

Cells `s > L_i` are never observed for unit `i`; they enter only through the
posterior predictive `tau_1(x_i) + sum_{k >= 2} tau_k(x_i) psi_k(s)`.
`res.catt_observed` marks the observed treated cells, and the transport audit
reports the share of the cumulative effect that rests on unobserved cells.

## What the fit does

1. **Common basis.** One orthonormal exposure basis `Psi (H, K)` on `s = 1..H`
   (`res.bases = {t0: Psi}`, `res.horizon = H`): mode 1 is the constant, so
   `tau_1(x_i)` is the time-averaged effect over the *full* horizon for every
   unit, whatever its window; the higher modes are zero-sum dynamic contrasts
   comparable across units. (The paper's per-cohort bases make mode 1 the
   average over each cohort's *own* window, which is not comparable across
   windows when effects decay.)
2. **NaN-aware demeaning.** Each unit's pre-period anchor is the mean of its
   observed pre cells (at least one is required); control period means use the
   observed never-treated units of that period (one is required wherever a
   treated outcome is observed).
3. **Two-pass projection** (`projection="two_pass"`, default; see
   `longbet._rank_k_event_time`):
   * higher modes `k >= 2` on the units whose window supports them
     (`min_support_modes`, default `(1, 6, 8)`: mode 2 needs `L_i >= 6`, mode 3
     `L_i >= 8`): `theta_hat_ik = Cov_{O_i}(psi_k, r_i) / Var_{O_i}(psi_k)`, the
     per-unit OLS slope of the observed residual on `{1, psi_k}` (on a partial
     window the constant and `psi_k` are not orthogonal, so the level has to be
     partialled out; on a full window this is the paper's `<psi_k, r_i> / H`),
     with precision weight `sum_{O_i} (psi_k - mean psi_k)^2`; units below
     support get weight 0 and their `tau_k(x_i)` comes from the cross-sectional
     posterior (the forest's prediction at `x_i`);
   * mode 1 on every unit as the observed-window mean of the residual after the
     fitted treatment dynamics are removed,
     `theta_hat_i1 = (1/L_i) sum_{s <= L_i} [r_i(s) - z_i sum_{k>=2} tau~_ik psi_k(s)]`,
     with precision weight `L_i` (plus the pre-period anchor term). Because
     `sum_{s=1}^H psi_k(s) = 0`, the correction is exactly zero on fully observed
     windows, where the fit reduces bit for bit to the released estimator.
   `mode_coupling` picks `tau~_ik`: `"propagated"` (default) runs Mode 1 on
   `n_forecast_draws` tracks, each with one posterior draw of the higher modes
   (and one forecast draw when there are forecast draws), exactly like forecast
   draws are propagated; `"plug_in"` uses the posterior means; `"none"` skips the
   correction (a window-mean estimator; use only as an ablation).
   `projection="restricted_wls"` fits every unit's supported modes jointly,
   `(Psi_O' Psi_O)^{-1} Psi_O' r_O`, with precision weights from the restricted
   Gram matrix (a sensitivity option: it extrapolates a unit's own ramp from its
   few points).
4. **Covariates.** The observed window length `L_i` is appended to the forest
   covariates (`support_covariate=True`), so the forests can learn
   cohort-specific baselines and effects where the data support them and shrink
   elsewhere. Add the activation month (and anything else that drives both
   censoring and dynamics) to `x` yourself.
5. **Weights.** `att_weighting="cell"` keeps the paper's `L_i / sum_j L_j`
   (observed cells); `"unit"` is `1 / N_trt`. `res.unit_weights` stores them and
   `res.arm_metadata` carries the treated-average trajectory draws and the arm
   residual covariances that `reconcile_nested_3arm_draws` and
   `reconcile_multi_arm_draws` read (they never consume `att_draws`).
6. **Cell draws.** `res.catt_draws[:, i, t0 + s - 1]` is filled for every post
   cell `s <= H` of every treated unit, observed or not; the Stage-1 options
   `catt_floor` / `tau_clip` act on the full horizon as before (the multi-arm
   reconciliation applies no projection of its own).

## Assumptions

* **Administrative censoring.** `L_i` is set by the activation date and the
  panel cut, not by the treatment or the outcome (no differential attrition).
  The validator only allows trailing censoring of the post window and missing
  pre cells with at least one observed pre cell per unit.
* **No cohort-by-event-time effect modification given `x`.** The effect
  profile of a unit observed for 5 months at event times 8-12 is the profile of
  comparable (same `x`, same `L_i` where supported) units that were observed
  there. `rank_k_transport_audit` tests this on the event times both strata
  observe (default `s <= 5`): the posterior of the difference of the two
  strata's unit-weighted mean profiles (`P5 / P30 / P50 / P70 / P95`,
  `prob_positive`) and a `flag` when its 90% interval excludes 0 (or an absolute
  `threshold`). Also report `unobserved_effect_share`.
* **Basis shape.** The decay form between the observed and the predicted cells
  is the span of the basis; report `rank_k=2` against `rank_k=3` and, if in
  doubt, a basis family sensitivity (`ExposureBasisSpec`). Under a mask the ATT
  depends on the fitted dynamics through the correction, so Proposition 2's
  `K`-invariance holds only on fully observed panels; the identity
  `sum_{s <= H} tau_i(s) = H tau_i1` still holds draw by draw.
* **Uncertainty.** The plug-in coupling understates Mode-1 uncertainty for
  short-window units (the higher-mode posterior is treated as known); the
  propagated coupling carries it through the tracks. Unit-level intervals remain
  approximate (see the paper).

## Forecast features and heterogeneous pre coverage

The counterfactual forecast must be complete (`features.point_forecast` or
`baseline_forecast` is `(N, T)` and finite). When the pre-period is complete,
build `features` as usual on the same `(y, z)`. When pre coverage is
heterogeneous, either build the features on the common observed pre window
(the last `c = min_i #observed pre cells` pre columns, passing `t_0_panel=c` to
the fit on the shortened panel) or run without a forecast
(`use_forecast=False`, `baseline_forecast=None`); both paths are tested.
Forecast *draws* come from `baseline_forecast_draws` or from `features`
(drawn on demand); the internal Normal-inverse-gamma generator needs complete
control outcomes, so on a masked panel with a bare `baseline_forecast` the fit
conditions on the point forecast and says so (`RankKConfigWarning`).

## Limits of this release

* Masked panels need a single adoption period (align the panel in event time;
  the natural way to handle staggered rollouts here) and the JAX engine; the
  plug-in `stochtree` engine, `fit_rank_k_hybrid` and `rank_k_placebo_audit`
  raise a clear error on masked panels.
* `horizon` is at most the post-adoption length of the panel; to predict beyond
  it, append all-`NaN` columns to `y`.
* Out-of-sample prediction (`predict_catt` for units not in the fit) is
  deferred: every unit of the population is a fitted unit here.
