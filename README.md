# LongBet

Bayesian causal forests for longitudinal experiments and panel data, implemented in JAX with Python and R interfaces calling the same compiled engine. Four estimators are provided:

- **`fit_rank_k_longbet`** (Forecast-Augmented Rank-$K$ Exposure-Basis BCF): unit-level heterogeneous *and* dynamic treatment effects in randomized or covariate-stratified staggered rollouts. One estimator with two uncertainty variants: "propagated" (the default: `uncertainty_propagation=True`, engine `"jax"`), which jointly propagates forecast-draw, control-group slope, pre-anchor factor-loading, and spatial/cluster uncertainty under the control-only cut, and "plug-in" (`uncertainty_propagation=False`, engine `"stochtree"`), a fast plug-in approximation for unit-level screening and targeting under unit-level randomization without spatial clustering. The recommended default for longitudinal experiments.
- **`LongBet`**: cell-level model $\tau_{it} = \beta_{S_{it}}\nu(X_i)$ for scalar continuous, binary, or ordinal outcomes under staggered or simultaneous adoption, with a conjugate counterfactual-forecast channel.
- **`LongBetMulti`**: joint seemingly unrelated regression (SUR) across multiple continuous, binary, and ordinal panel outcomes.
- **`LongBetEncourage`**: longitudinal encouragement designs with an assigned offer wave and endogenous, absorbing treatment adoption.

All of them can consume the leak-free pre-treatment features of **`build_causal_features`** (counterfactual baseline forecast, trajectory summaries, selection score), computed strictly from data before the first adoption date.

## Installation

From a source checkout using Python 3.12 on Linux or macOS:

```sh
git clone https://github.com/ignacio82/longbet-jax.git
cd longbet-jax
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -e .
```

The Python distribution is named `longbet-jax` and imports as `longbet`. The base install runs every estimator, including both variants of `fit_rank_k_longbet`: `stochtree`, which the fast "plug-in" variant uses, is a core dependency since 0.5.0, and the `[rank-k]` extra is kept (empty) so that existing `".[rank-k]"` installs keep working. For NVIDIA CUDA acceleration of the JAX engine, optional extras `.[cuda12]` and `.[cuda13]` are provided.

To install the R wrapper from the repository root:

```sh
Rscript -e 'install.packages(c("reticulate", "rlang"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

Before loading `longbet` in R, bind `reticulate` to the virtual environment containing the Python engine:

```r
reticulate::use_virtualenv("/absolute/path/to/longbet-jax/.venv", required = TRUE)
library(longbet)
```

---

## What to use when

The recommendations below are the practitioner guidelines of the accompanying paper (*Forecast-Augmented Rank-K Bayesian Causal Forests for Longitudinal Experiments*), whose evidence comes from a pre-registered real-data semi-synthetic benchmark ($R = 64$ paired replications per cell, $N = 250$ units, NOAA GHCN-M station temperature and BLS LAUS county employment), synthetic diagnostics, and a 50-state application. Numbers quoted here are from that paper.

| Your situation | Use | Why |
|---|---|---|
| **Randomized or covariate-stratified rollout**, longitudinal outcomes, effects may differ across units in size, sign, or shape; $N$ of a few hundred; at least $L_0 \ge 12$ pre-treatment periods (history + in-panel) | `fit_rank_k_longbet(y, z, x, features=build_causal_features(...), coords=latlon)` with the defaults: the "propagated" variant (`uncertainty_propagation=True`, engine `"jax"`) with `rank_k = 2` (pass `coords=latlon` or `cluster_ids=` whenever units can share regional shocks). For fast unit-level screening or targeting under unit-level randomization without spatial clustering, the fast "plug-in" variant: `cfg=RankKLongBetConfig(uncertainty_propagation=False)` | In the paper's RCT study the propagated default passed the pre-registered ATT coverage rule on every development panel (93.8% on GHCN-M temperature and 90.6% on BLS county employment) and every hold-out set of the second frozen confirmation (pooled 95% coverage 94.4%: 96.1%, 93.0% and 94.1% across Sets E–G, SD ratios 0.93–1.02; rules H1′–H7 passed), and remained calibrated in synthetic panels with spatially clustered shocks, loading-correlated assignment or uncertain forecast slopes where the plug-in variant under-covered. The fast plug-in variant missed the pooled coverage rule on one Stage-1 hold-out set (89.8% against 90%), though both variants passed in the second confirmation (93.6% plug-in vs 94.4% propagated, McNemar $p = 0.53$, with propagated intervals 24% longer). Propagating the control-only baseline posterior carries a modest finite-sample variance price on unit-level point recovery (propagated unit PEHE was 0.010 ± 0.004 higher than plug-in on covariate-stratified cells, failing non-inferiority rule H8; unit sign accuracy 71.7% vs 83.4% on GHCN-M and 61.0% vs 70.0% on BLS in the sign-flip regime). Unit-level (CATE) intervals are approximate and below nominal (95% unit-level intervals covered 72–93% with the propagated and 88–93% with the plug-in variant): treat them as descriptive. By construction the ATT and unit-averaged effects do not change when $K$ is raised; report $K = 1$ alongside $K = 2$ (identical ATT, different cell profiles). |
| Same design, but you need the **shape of the dynamic response** beyond a level + ramp (novelty wear-out vs. delayed ramp) | `rank_k = 2` first; raise to 3 only when every treated wave has $L_g \ge K$ post-adoption periods | In the paper's sign-reversing-dynamics benchmark the terminal-effect sign accuracy was 79.6% for the propagated default (87.2% for plug-in) at $K = 2$, against $\le 48\%$ for $K = 1$ and cohort-average estimators; higher modes add parameters per wave, and only $K \le 3$ was validated. |
| **Short or uninformative history**: `feats.gate_open` is `False` ($L_0 < 12$ or out-of-fold non-polynomial $R^2 < 0.05$) | Default `LongBet` with no forecast features (`x=feats.augment_x(x)`, `ps=feats.recommended_ps`, `baseline_forecast=feats.baseline_forecast` do exactly that automatically), or `fit_rank_k_longbet` with `use_forecast=False` in a randomized design | On smooth polynomial baselines or short contexts the trajectory summaries alone add variance without removing confounding; the gate-aware accessors revert identically to the default model. |
| **Small observational panel** (tens of units, selection on observables), and you want a cell-level exposure profile $\beta_S$ | Cell-level `LongBet` with the gate-aware features (`augment_x`, `recommended_ps`, `baseline_forecast`), plus `n_forecast_draws=8` pooled across fits when the series are serially persistent, **and report DR-DiD alongside** | The conjugate channel cut CATT PEHE by about 75% on the seasonal synthetic designs and removed the seasonal confounding of default `LongBet` on raw temperature series, but its plug-in ATT intervals failed the pre-registered rule on both panels (88.3% / 65.6% pooled coverage); the cut posterior raised BLS null-cell coverage from 56% / 73% to 72% / 94% while also moving the centres. DR-DiD was the only comparator within Monte Carlo error of nominal in every nominal cell (94%–100%), with intervals about three times wider. On the 50-state application ($N = 50$, 25 treated) the unit-grain Rank-$K$ estimator is shrunk and the cell-level channel is preferable. |
| **Selection on the recent untreated trajectory** (Ashenfelter-dip-type adoption) | No estimator here is unbiased; treat the forecaster's late-context kernels as assumptions, re-run under `late_context="none"` and moved kernels (`LateContextSpec`), inspect `late_context_evidence`, report $\hat R_{\max}$ and ATT-level ESS | Under trajectory-based selection every estimator in the benchmark was biased and on BLS the raw-input channel was worse than the default; poor mixing accompanied the failure. No validated feasible statistic identifies this selection from a single panel. |
| Several outcomes on the same panel | `LongBetMulti` | Shared design and adoption schedule, outcome-specific forests, SUR-coupled innovations, joint counterfactual probabilities. |
| Randomized **encouragement** (offer) with endogenous adoption | `LongBetEncourage` | Outcome equation on the adoption clock composed with a discrete-time adoption hazard on the offer clock. |
| **Administrative censoring in event time** (every unit adopts at `t0` of its own clock, the panel is cut on a calendar date, so unit `i` is observed for `L_i` post periods) and you want the `h`-period effect **for the whole population**, not for the units observed through `h` | `fit_rank_k_longbet(y, z, x, cfg=RankKLongBetConfig(att_weighting="unit"))` with the unobserved cells of `y` set to `NaN`; then `cumulative_effect(res, h)`, `event_time_trajectory(res)` and `rank_k_transport_audit(res, strata)`. See [docs/rank-k-masked-event-time.md](docs/rank-k-masked-event-time.md). | Every unit stays in the fit with the cells it has: one common exposure basis on `s = 1..H`, a two-pass projection whose Mode 1 is the window mean after the fitted dynamics are removed (exactly the paper's estimator on complete panels, bit for bit), zero weight for modes a short window cannot support, and posterior-predictive cell effects on the unobserved cells. The transport audit compares the early-month profiles of long- and short-window units, the assumption that carries the prediction. |
| **Multi-arm experiments** ($A \ge 3$ arms in a nested chain or partial order, such as a dose ladder or $2 \times 2$ factorial) across **one or more outcomes and composites** | `reconcile_multi_arm_draws(results, arms, ...)` for a single outcome, or `reconcile_multi_arm_multi_outcome_draws(results_by_outcome, arms, ..., composite_outcomes=..., composite_mode="coherent_total")` across $M$ outcomes; then `multi_arm_decision_table(...)` and `roi_break_even_table(..., contrast=...)`. See [docs/rank-k-multi-arm-multi-outcome.md](docs/rank-k-multi-arm-multi-outcome.md). | Independent pairwise Rank-$K$ BCF fits ignore shared reference-arm sampling covariance $\operatorname{Cov}(\hat{\tau}_{(a,b)}, \hat{\tau}_{(b,c)}) = -v_b$, violate cycle additivity $\tau_{(a,c)} = \tau_{(a,b)} + \tau_{(b,c)}$, allow quantile reversals under positivity constraints, and inflate null composite means when $K$ truncated component outcomes are summed. Cut-space projection, NNLS arm-variance recovery (with automatic bipartite detection), Dirichlet chain / order-polytope poset bundle priors, Iman-Conover copula coupling, and single-truncation coherent composite rescaling enforce exact draw-wise cycle and composite additivity ($\max|\text{gap}| < 10^{-12}$, $\text{AVR} = 0$). |
| You need to know whether the effect is **stable or wearing out**, or a return-on-investment question reaches **past the end of the panel** (payback over several years from a one-year experiment) | `diagnose_trajectory_stability(...)` for the in-sample evidence (`wearout_ancova`, `estimate_post_pulse_retention`, `pool_multi_outcome_retention` with a pre-registered `RetentionRuleSpec`), then `compute_required_persistence(...)` with Bounds U / L, the `rho = 0` floor and a placebo; `extrapolate_rank_k_decay(..., pulse_decay_tau=..., survival=...)` for scenario trajectories. See [docs/rank-k-stability-and-required-persistence.md](docs/rank-k-stability-and-required-persistence.md). | A panel powered for the cumulative effect is usually not powered for its slope: the minimum detectable wear-out can exceed the effect itself, positivity constraints pile retention draws up at `rho = 1`, and outcomes measured on the same units are not independent replicates. Periods after the panel are not identified, so the decision is inverted into the persistence it requires and reported next to labelled scenario brackets rather than as a forecast. No function certifies a multi-year extrapolation. |

Things you should never do, because each one was measured to hurt:

- **Never pass a forecast through `x_tv`** (the time-varying tree covariates). It disables the collapsed unit intercepts and discretizes a smooth trajectory; naive injection cost 1.7–4.5× in minimum ESS. Forecasts enter only through `baseline_forecast` (linearly, doubly centred, inside the conjugate trend block) or through the Rank-$K$ CUPED slope.
- **Never build a forecast on a rolling origin or on post-adoption data.** `build_causal_features` anchors every unit at the first adoption date and masks every post-anchor cell of ever-treated units; a rolling origin attenuated the dynamic-confounding ATT toward zero.
- **Never pass a time-varying propensity `(N, T)` to `LongBet`.** `feats.ps` is a unit-constant selection score for that reason.

### Defaults and why

| Option | Default | Why this default |
|---|---|---|
| `build_causal_features(forecaster=...)` | `None` = `patch_ridge_forecaster` | The guarded-core spectral forecaster with cross-series patch-ridge refinement is the backbone used for every forecast arm of the paper. The refinement validates itself out-of-fold and switches off (reducing to the spectral backbone) when it does not predict; `"spectral"` selects the backbone alone. A pretrained time-series foundation model (TimesFM) was evaluated and removed: it did not dominate the built-in backbone on the causal targets, and it required a GPU and a 330M-parameter checkpoint. |
| `late_context` | `"fixed"` | Slope-break hinge 10 periods before the first adoption and a Gaussian pre-adoption dip (width 1.8, centred 1.5 periods before the first post-anchor period), amplitudes projected on the cross-unit harmonic factor. These are analyst-set identification restrictions; a data-adaptive search never dominated them and selected on pre-treatment noise (see §1 below). |
| gate (`feats.gate_open`) | $L_0 \ge 12$ and `beyond_poly_r2 >= 0.05` | Below either threshold the forecast column is redundant with `LongBet`'s polynomial + DCT-II trend block and the gate-aware accessors revert to the default model. The Rank-$K$ estimator does not use the gate: its default propagated variant samples the slope $\theta_k$ from never-treated units under a $\mathcal{N}(1, 0.5^2)$ prior truncated to $[0, 1.5]$, and the plug-in variant clips a ridge slope of the forecast fitted on never-treated units to the same interval. |
| `n_forecast_draws` | `0` | Forecast draws in `build_causal_features` cost extra bootstrap work when only `LongBet` is used. Use `8` for the cell-level channel on serially persistent series (within-unit AR(1) $\rho_1 \approx 0.76$ on monthly county employment) to pool a modular posterior. `fit_rank_k_longbet` does not need them from `build_causal_features`: its default propagated variant draws its own posterior forecast panels on demand (`RankKLongBetConfig.n_forecast_draws`, default 8), and the plug-in variant uses none. |
| `LongBetConfig.baseline_cells` | `"all"` | See §5: fitting the cell-level baseline on untreated cells only raised ATT bias on short pre-adoption windows. |
| `RankKLongBetConfig` | `rank_k=2`; the propagated variant (`uncertainty_propagation=True`, engine `"jax"`): 2 chains, 100 burn-in sweeps, `n_forecast_draws=8` without adaptive expansion (`max_forecast_draws=None`) and 25 draws per chain and forecast draw ($2 \times 8 \times 25 = 400$ draws), 50 / 30 / 15 trees, spatial field and factor loadings ($q = 3$) on, `baseline_cells="untreated"`, `theta_f_clip=(0, 1.5)`. With `uncertainty_propagation=False` (engine `"stochtree"`, the fast plug-in variant): one chain per mode with 10 grow-from-root sweeps, 200 burn-in sweeps and 600 retained draws, 75 / 35 / 20 trees | Each variant gets the budget behind its results in the paper: for the propagated default the benchmark budget of the development benchmark and second frozen confirmation, and for the plug-in variant the configuration of every plug-in result. `engine`, `num_chains`, `num_burnin` and the tree counts default to `None` and are filled in for the selected variant. For small panels use the paper's application budget, `RankKLongBetConfig(num_chains=4, num_burnin=300, draws_per_component=100)`. `max_forecast_draws > n_forecast_draws` (adaptive expansion) and `rank_k > 3` were not validated. Mode 1 is seeded independently of $K$, so ATT and unit-averaged effects are invariant to $K$. |

---

## 1. `fit_rank_k_longbet`: Forecast-Augmented Rank-$K$ Exposure-Basis BCF

### Model

For a unit with first adoption period $G_i$ and post-adoption horizon $L_g = T - G_i$, the cell effects are expanded on an orthonormal exposure basis $\{\psi_{L_g,k}(s)\}_{k=1}^{K}$ on $s = 1, \dots, L_g$ (constant first mode, then an orthonormalized concave ramp $1 - e^{-s/3}$, a linear tilt, and cosine harmonics),

```math
\tau_{it} = \sum_{k=1}^{K} \psi_{G_i,k}(S_{it})\, \nu_k(X_i, u_i),
```

so that mode 1 is exactly the unit's time-averaged post-adoption effect and every higher mode integrates to zero over the unit's post-adoption window. The estimator proceeds in five steps:

1. **Leak-free counterfactual forecast** $F_{it}$ from `build_causal_features` (built-in backbone anchored at $t_0 = \min_i G_i - 1$, refined by a $K$-fold cross-fitted control-group adapter).
2. **Pre-anchor level centring and control-group calendar demeaning** of $Y$ and $F$, which remove unit intercepts and arbitrary common calendar shocks across staggered waves.
3. **Within-unit projection** of the centred post-adoption outcomes and forecasts onto the exposure basis, giving one scalar response $y_{ik}$ and one scalar forecast $f_{ik}$ per unit and mode; never-treated units receive the wave-weighted average over the treated waves' windows.
4. **One cross-sectional BCF per mode** with the forecast as a control variate. In the default **"propagated"** variant (`uncertainty_propagation=True`, engine `"jax"`) every baseline component is learned from never-treated units under the control-only cut: a prognostic forest $\mu_k(X_i, u_i, \tilde f_{ik}, \hat\pi_i)$; the forecast control-variate slope $\theta_k$ with prior $\mathcal{N}(1, 0.5^2)$ truncated to $[0, 1.5]$; pre-anchor factor loadings $\lambda_i' f_k$ ($q = 3$ loadings of the pre-anchor forecast residuals); and a spatial random field $\zeta_{ik}$, a Matérn-3/2 Gaussian process on great-circle distances between unit coordinates (`coords`) or a cluster random effect (`cluster_ids`), kriged from never-treated to treated units. They are sampled together with the treatment forest $\tau_k(X_i, u_i)$, fitted on the treated units' remainder, and each unit gets a precision weight $w_{ik}$ inversely proportional to the variance of its projection, so units with short post-adoption windows weigh less. Mode 1 is refitted on each of $M = 8$ posterior forecast draws and the draws are pooled into a modular (cut) posterior (higher modes condition on the point forecast). The sampler is an approximate partially collapsed Gibbs sampler (the forest steps use working likelihoods), not an exact posterior sampler.
5. **Reconstruction** of cell draws $\tau_{it}$, unit-averaged effects, and the treated-cell ATT from the mode draws.

The fast **"plug-in"** variant (`uncertainty_propagation=False`, engine `"stochtree"`) replaces step 4 by a plug-in ridge regression slope $\hat\theta_k$ of $y_{ik}$ on $f_{ik}$ among never-treated units, clipped to $[0, 1.5]$, and fits one `stochtree.BCFModel` per mode to $y_{ik} - \hat\theta_k f_{ik}$ with covariates $(X_i, u_i, \tilde f_{ik})$ and the selection score $\hat\pi_i$ in the propensity slot, without forecast draws, factor loadings or spatial field.

**Which variant.** Use the propagated variant (the default) for full uncertainty propagation across forecast draws, the control-group slope, pre-anchor factor loadings, and spatial/cluster correlation (pass `coords` or `cluster_ids` whenever units can share regional shocks). Use the fast plug-in variant (`uncertainty_propagation=False`) for rapid unit-level screening and targeting under unit-level randomization without spatial clustering (see [What to use when](#what-to-use-when)). In both variants the 95% ATT interval is the 2.5% / 97.5% posterior percentiles of the ATT draws. Unit-level (CATE) intervals are approximate and below nominal in the paper's benchmarks, and the cell-level intervals (2.5% / 97.5% posterior percentiles) are approximate as well; treat both as descriptive.

Because mode 1 is fitted with the same seed whatever $K$ is, and all higher modes telescope to zero within each unit's window, the posterior of the ATT and of every unit's time-averaged effect is identical for $K = 1, 2, 3, \dots$ (Proposition 2 of the paper); raising $K$ only refines the within-unit dynamic profile.

### Python

```python
import numpy as np
from longbet import RankKLongBetConfig, build_causal_features, fit_rank_k_longbet, rank_k_placebo_audit

# y: (N, T) outcomes, z: (N, T) absorbing 0/1 treatment, x: (N, P) covariates,
# y_hist: (N, T_hist) strictly pre-panel history (optional but recommended),
# latlon: (N, 2) unit (latitude, longitude) in degrees, or region: (N,) cluster labels
feats = build_causal_features(y=y, z=z, x=x, y_pre_history=y_hist, seed=1)
res = fit_rank_k_longbet(y, z, x, features=feats, coords=latlon)    # propagated default, K = 2
# without coordinates: fit_rank_k_longbet(y, z, x, features=feats, cluster_ids=region)

print(res.att_hat, res.att_ci_low, res.att_ci_high)                 # treated-cell ATT
treated = z.any(axis=1)
unit_avg = res.mode_tau_draws[0].mean(axis=1)[treated]              # unit-averaged effects
res_k1 = fit_rank_k_longbet(y, z, x, features=feats, coords=latlon, cfg=RankKLongBetConfig(rank_k=1))
np.testing.assert_allclose(res.att_draws, res_k1.att_draws)        # ATT invariant to K
audit = rank_k_placebo_audit(y, z, x, y_pre_history=y_hist, coords=latlon)
print(audit.passes_audit)                                           # 95% placebo interval covers 0

# Fast plug-in approximation (stochtree engine; no forecast draws, loadings, or spatial GP)
cfg_plugin = RankKLongBetConfig(uncertainty_propagation=False)
res_plugin = fit_rank_k_longbet(y, z, x, features=feats, cfg=cfg_plugin)

# Domain-informed Mode-1 ATT prior (e.g., non-negativity constraint ATT >= 0 and informative scale)
cfg_informed = RankKLongBetConfig(
    tau_prior_mean=0.0,
    tau_prior_sd=1.5,
    tau_clip=(0.0, float("inf")),
    catt_floor=0.0,  # optional post-adoption period floor enforcing tau_{it} >= 0
)
res_informed = fit_rank_k_longbet(y, z, x, features=feats, coords=latlon, cfg=cfg_informed)

# Coherent 3-arm nested reconciliation (Group C = Control, B = Base, A = Base + Incremental):
# accounts for shared Group B covariance Cov(tau_BC, tau_AB) = -v_B < 0 and enforces tau_AC = tau_BC + tau_AB
from longbet import reconcile_nested_3arm_draws
rec = reconcile_nested_3arm_draws(res_bc, res_ab, res_ac, prior_scale=0.3)  # explicit scale in outcome units
trt_catt_bc, trt_catt_ab, trt_catt_ac = rec                                # unpacks as the three trajectories

# Is the effect stable, and how long must it last to pay back a per-unit cost? Scenario brackets, not a forecast
# (docs/rank-k-stability-and-required-persistence.md)
from longbet import compute_required_persistence, discount_weights, weibull_survival_weights
decision = compute_required_persistence(
    rec, cost_grid=(2.0,), alpha_grid=(0.5, 1.0), horizon=60, pulse_decay_tau=3.0,
    survival=weibull_survival_weights(60, rec.trt_catt_ac.shape[1]), discount=discount_weights(60),
)
print(decision.required_persistence_df)                            # X quantiles and Pr(X <= 12 / 24 / 48)
```

`res.catt_draws` has shape `(S, N, T)` (posterior draws of $\tau_{it}$, zero on untreated cells), `res.catt_hat` / `res.catt_ci_low` / `res.catt_ci_high` are `(N, T)` summaries, `res.att_draws` is `(S,)`, `res.mode_tau_draws` holds one `(N, S)` array per mode (row $i$ of the first is unit $i$'s time-averaged effect, defined for treated units), `res.theta_f_modes` holds the forecast slopes, `res.rhat_att` / `res.ess_att` are convergence diagnostics (`rhat_att` needs two chains, so it is NaN for the plug-in variant, which runs one chain per mode), and `res.bases` holds the exposure basis per wave. Domain-informed priors on the per-period Mode-1 treatment effect intercept are supported via `tau_prior_mean`, `tau_prior_sd` (in outcome units; default `None` uses `10.0` on the standardized mode scale), `tau_clip` (e.g., `(0.0, float("inf"))` for non-negative average treatment effects), and `catt_floor` (optional post-adoption cell floor). For nested 3-arm designs (`B vs. C` base treatment, `A vs. B` incremental treatment, `A vs. C` combined treatment), `reconcile_nested_3arm_draws` reconciles the three pairwise `RankKLongBetResult` objects under the shared Group B sampling covariance and returns coherent `(S, T_post)` trajectory draws satisfying `trt_catt_ac = trt_catt_bc + trt_catt_ab` with `Combined >= Incremental >= 0` and `Combined >= Base >= 0` at every posterior quantile (`P30`, `P50`, `P70`). At the defaults $S = 2 \text{ chains} \times 8 \text{ forecast draws} \times 25 = 400$ for the propagated variant and $S$ = `num_mcmc` = 600 for the plug-in variant. Pass `cohorts=` and `t_0_panel=` explicitly if they should not be inferred from `z`. Invalid inputs raise `ValueError` before any computation, and so do the contradictory settings `engine="jax"` with `uncertainty_propagation=False` and `engine="stochtree"` with `uncertainty_propagation=True` (`engine` is kept for backward compatibility; on its own it selects its variant). A `longbet.RankKConfigWarning` flags explicit inputs that the fit would ignore or that contradict each other (for example `coords` with the plug-in variant). `dataclasses.replace` keeps a configuration's resolved engine and budget: replacing `uncertainty_propagation` alone raises `ValueError`, and replacing it together with `engine=None` keeps the other variant's budget (the fit then warns), so build a new `RankKLongBetConfig` to switch variants.

### R

```r
feats <- build_causal_features(y = y, z = z, x = x, y_pre_history = y_hist, seed = 1)
res   <- fit_rank_k_longbet(y, z, x, features = feats, coords = latlon)   # propagated default, rank_k = 2
c(res$att_hat, res$att_ci_low, res$att_ci_high)
treated  <- rowSums(z) > 0
unit_avg <- rowMeans(res$mode_tau_draws[[1]])[treated]                    # unit-averaged effects
res_k1   <- fit_rank_k_longbet(y, z, x, features = feats, coords = latlon, rank_k = 1)
audit    <- rank_k_placebo_audit(y, z, x, y_pre_history = y_hist, coords = latlon)
# Fast plug-in approximation (stochtree engine)
res_plugin <- fit_rank_k_longbet(y, z, x, features = feats,
                                 uncertainty_propagation = FALSE)
```

The R defaults are those of `RankKLongBetConfig()`: the variant-dependent arguments (`engine`, `uncertainty_propagation`, `max_forecast_draws`, `num_chains`, `num_burnin` and the three tree counts) default to `NULL` and are resolved in Python for the selected variant. `cohorts` is the zero-based first treated period (`ncol(y)` for never-treated units) and `t_0_panel` a single integer, as in Python. Counts must be single whole numbers and flags single `TRUE` / `FALSE` values. The arguments added in 0.5.0 (`pre_anchor_loadings`, `cohorts`, `max_forecast_draws` and `draws_per_component`) come after `num_threads`, so 0.4.0 positional calls keep their meaning; pass them by name. A `RankKConfigWarning` is issued through Python's `warnings` module, not as an R condition.

---

## 2. `LongBet`: Heterogeneous Dynamic Treatment Effects in Panel Data

### Model

For unit $i \in \{1, \dots, N\}$ in period $t \in \{1, \dots, T\}$, let $Z_{it} \in \{0, 1\}$ be an absorbing treatment indicator and $S_{it} \ge 0$ be the exposure clock ($S_{it} = 0$ prior to adoption, $S_{it} = 1$ in the adoption period, and incrementing by one each subsequent period). For a continuous outcome,

```math
Y_{it} = \mu(X_i, X^{\mathrm{tv}}_{it}) + \eta_{it}
       + \beta_{S_{it}}\, \nu(X_i, X^{\mathrm{trt,tv}}_{it})\, \mathbf{1}\{S_{it} \ge 1\}
       + \gamma_i + \epsilon_{it},
\qquad \epsilon_{it} \sim \mathcal{N}(0, \sigma^2).
```

- **Prognostic and moderating forests ($\mu, \nu$).** $\mu$ is a BART sum of trees over baseline covariates $X_i$ and optional time-varying covariates $X^{\mathrm{tv}}_{it}$; $\nu$ is a moderating BART sum of trees over $X_i$ (or a distinct moderation matrix `x_trt`) and optional time-varying covariates $X^{\mathrm{trt,tv}}_{it}$. Following Bayesian Causal Forests (Hahn, Murray, and Carvalho, 2020), a node at depth $d$ splits with probability $\alpha/(1+d)^\beta$ ($0.95/(1+d)^2$ for $\mu$, $0.25/(1+d)^3$ for $\nu$), and leaf parameters have Gaussian priors $\mathcal{N}(0, \tau^2)$ on the standardized response scale with $\tau_\mu = 1.5/\sqrt{m_\mu}$ ($m_\mu = 20$ default trees) and $\tau_\nu = 0.5/\sqrt{m_\nu}$ ($m_\nu = 60$ default trees). Neither $\mu$ nor $\nu$ splits directly on raw calendar time $t$ or exposure $S_{it}$, preserving clean separation between the exposure profile $\beta_s$, covariate moderation $\nu(X_i)$, and the untreated baseline trend $\eta_{it}$.
- **Conjugate baseline trend block ($\eta_{it}$, `use_trend_block = True`).** Common period shocks, optional doubly-centered counterfactual baseline forecasts $\widetilde{F}_{it} \in \mathbb{R}^M$ (`baseline_forecast`), and covariate-dependent secular trends are captured by
  ```math
  \eta_{it} = \Phi^{\mathrm{time}}_t c_0 + \widetilde{F}_{it}^\top c_{\mathrm{forecast}} + \sum_{k=1}^{d} \phi_k(t)\, b(X_i)^\top c_k,
  ```
  where $\Phi^{\mathrm{time}} \in \mathbb{R}^{T \times (T-1)}$ (`trend_period_effects = True`) is an orthonormal centered DCT-II basis with Gaussian prior scale `trend_period_scale = 2.0`, $\widetilde{F}_{it}$ is a unit- and period-centered, RMS-standardized counterfactual forecast block with Gaussian prior scale `trend_forecast_scale = 1.0`, and $\phi(t)$ is an orthogonal polynomial basis of degree `trend_degree = 2`. Entering $\widetilde{F}_{it}$ linearly inside the conjugate trend block (rather than inside the tree ensemble $\mu$) preserves smooth slope shrinkage, joint block-Cholesky sampling of $(c, \gamma, \beta)$, and analytical marginalization of the unit intercepts $\gamma_i$ (`collapse_gamma_mu = True`). The unit feature map $b(X_i)$ ($p_0 + F$ columns, with $p_0 = \min(p, 24)$ and $F = \text{`trend\_num\_features`} = 64$ default nonlinear columns) concatenates standardized linear covariates $\widetilde{X}_i$ with a standardized coordinate-based dictionary of 1D piecewise-linear hinges $\{|\widetilde{X}_{ij}-1|, |\widetilde{X}_{ij}|\}_{1 \le j \le \min(p_0, 8)}$, 2D bilinear products $\{\widetilde{X}_{ij}\widetilde{X}_{ij'}\}_{1 \le j < j' \le \min(p_0, 8)}$, and directed 2D hinge interactions $\{\widetilde{X}_{ij}|\widetilde{X}_{ij'}-1|\}_{1 \le j \ne j' \le \min(p_0, 6)}$, truncated or padded with isotropic random Fourier features $\psi(\widetilde{X}_i)$ to $F$ columns, blended elementwise with $+0.05\,\psi(\widetilde{X}_i)$, and centered and RMS-normalized across units. Fourier features use the fixed seed `20260918`; blending does not guarantee independent columns, and proper coefficient priors regularize rank-deficient designs. By default (`trend_horseshoe = True`), the interaction coefficients $c_k$ carry a global–local horseshoe prior $c_j \mid \lambda_j, \tau \sim \mathcal{N}(0, \tau^2 \lambda_j^2)$ that adapts between parallel trends ($\tau \to 0$) and diverging covariate-specific slopes. The inverse-Gamma scale updates use `float64` arithmetic with `float32` Gamma draws and stored scales; numerical bounds $[10^{-30},10^{30}]$ approximate the unbounded conditionals if reached.
- **Exposure trajectory ($\beta_S$) and $\mathbb{Z}_2$ orientation.** The dynamic shape $\beta = (\beta_1, \dots, \beta_{\bar S})^\top$ follows a Gaussian process over exposure time,
  ```math
  \beta \sim \mathcal{N}\!\big(0,\, K + \sigma_m^2 \mathbf{1}\mathbf{1}^\top\big),
  ```
  with a squared-exponential (`kernel_type = "se"`, default), Matérn-3/2 (`"matern32"`), Matérn-5/2 (`"matern52"`), or AR(1) (`"ar1"`) kernel $K$ of standard deviation `sig_knl = 1.0`, lengthscale `lambda_knl = 1.0`, and common-level standard deviation `sigma_m = 1.0`. After every trajectory draw, the idempotent two-to-one $\mathbb{Z}_2$ quotient projection $(\beta, \nu) \mapsto \mathrm{sign}(\sum_{s \ge 1}\beta_s)\,(\beta, \nu)$ maps the factors onto the canonical half-space $\sum_{s \ge 1}\beta_s \ge 0$, eliminating mirror-mode splitting across chains while leaving every cell-level treatment effect $\tau_{it} = \beta_{S_{it}}\nu(X_i)$ strictly invariant.
- **Unit random intercepts and MCMC transitions.** $\gamma_i \sim \mathcal{N}(0, \sigma_\gamma^2)$ models persistent unobserved unit heterogeneity (`random_intercept = True`), with $\sigma_\gamma^2 \sim \mathrm{IG}(1, 0.1)$ and $\sigma^2 \sim \mathrm{IG}(2, 1)$ on the standardized outcome scale. Whenever a forest's splits are unit-constant (`collapse_gamma_mu` / `collapse_gamma_nu`), the unit intercepts $\gamma_i$ are analytically marginalized out during forest updates via the Sherman–Morrison unit-sum projection (with exact unit-specific precision weighting $\sqrt{\kappa_{\mu, i}}$ on unbalanced panels and an intermediate $\gamma_i$ refresh after $\mu$ whenever `collapse_gamma_mu = True` and `collapse_gamma_nu = False` so the partially collapsed Gibbs transition remains target-invariant). Each sweep combines BART GROW/PRUNE, CHANGE, internal-CHANGE, and likelihood-guided REGROW proposals, a multiplicative scale-ridge Metropolis move on $(\beta, \nu)$, inter-ensemble subspace transfers (`use_inter_ensemble_move = True`), and an exact block Cholesky Gaussian draw of $(c, \gamma, \beta)$.
- **Which cells inform the untreated baseline (`baseline_cells`).** Under the default `baseline_cells = "all"`, the baseline components $(\mu, \eta, \gamma, \sigma^2)$ are updated from every observed cell (with treated cells entering through $Y_{it} - \beta_{S_{it}}\nu(X_i)$) and `trend_block_step` draws $(c, \gamma, \beta)$ from their exact joint Gaussian conditional given the trees. Setting `baseline_cells = "untreated"` fits the baseline as a separate module from untreated cells only ($Z_{it} = 0$) and conditions $(\beta, \nu)$ on it, targeting the modular (*cut*) posterior $p(\text{baseline} \mid Y_{Z=0})\,p(\beta, \nu \mid \text{baseline}, Y_{Z=1})$ (Plummer, 2015; cf. Borusyak, Jaravel, and Spiess, 2024; Liu, Wang, and Xu, 2024). Because treated cohorts contribute only their pre-adoption periods to $\eta$ and $\gamma_i$, their post-adoption baseline under `"untreated"` is an out-of-sample polynomial extrapolation over $t \ge G_i$; on short panels with covariate-driven selection or unit-specific seasonality without a forecast channel, shrinkage on the pre-adoption slope estimate increases rather than reduces finite-sample ATT bias (see Section 5 below for the empirical trade-off). Requires `outcome = "continuous"` and `tempering_levels = 1`, and is not available in `LongBetMulti`.
- **Binary and ordinal outcomes.** Setting `outcome = "binary"` fits a probit model $Y_{it} = \mathbf{1}\{Y^*_{it} > 0\}$ with unit latent innovation variance. Setting `outcome = "ordinal"` with `num_categories = K` fits an ordered probit model with marginalized log-gap Metropolis updates on the interior cutpoints.

### Python

```python
import numpy as np
from longbet import LongBet, LongBetConfig, get_att

rng = np.random.default_rng(0)
N, T = 200, 10
x = rng.normal(size=(N, 4))
z = np.zeros((N, T)); z[:100, 4:] = 1
t = np.arange(1, T + 1)
s = np.where(z == 1, np.maximum(t[None, :] - 4, 0), 0)
y = 0.5 * x[:, [0]] + 0.1 * t + 0.8 * np.sqrt(s) + rng.normal(0, 0.3, (N, T))

model = LongBet(LongBetConfig(num_chains=4, num_burnin=2000, num_sweeps=1000, n_skip=1))
model.fit(y=y, x=x, z=z, t=t)

pred = model.predict(x=x, z=z, t=t, summary_only=True)
att = get_att(pred)
diag = pred.stability()
print(att["exposure"], att["att"], diag.summary)
model.save("longbet.npz"); model = LongBet.load("longbet.npz")
```

### R

```r
library(longbet)
set.seed(0)
N <- 200; Tn <- 10
x <- matrix(rnorm(N * 4), N, 4)
z <- matrix(0, N, Tn); z[1:100, 5:Tn] <- 1
s <- ifelse(z == 1, pmax(col(z) - 4, 0), 0)
y <- 0.5 * x[, 1] + 0.1 * col(z) + 0.8 * sqrt(s) + matrix(rnorm(N * Tn, 0, 0.3), N, Tn)

fit  <- longbet(y = y, x = x, z = z, t = 1:Tn, num_chains = 4, num_burnin = 2000, num_sweeps = 1000)
pred <- predict(fit, x = x, z = z, t = 1:Tn, summary_only = TRUE)
att  <- get_att(pred)
diag <- att_stability(pred)
saveRDS(fit, "fit.rds"); fit <- readRDS("fit.rds")
```

### Counterfactual Forecast Augmentation (`baseline_forecast` & `build_causal_features`)

When a pre-treatment context history ($L_{\mathrm{ctx}} = T_{\mathrm{hist}} + t_0 \ge 12$) exhibits multi-frequency seasonality, local trend breaks, or pre-adoption dips that are not spanned by a low-degree polynomial in calendar time, `build_causal_features` constructs leak-free counterfactual forecasts and trajectory summaries anchored strictly at the common pre-treatment cutoff $t_0 = \min_{i \in \mathcal{I}_{\mathrm{trt}}} G_i$:

1. **Zero-shot counterfactual trajectory $\hat{y}^{\mathrm{ZS}}_{it}(0)$**: forecast from $y_{i, 1:t_0}$ by the built-in backbone (`forecaster=None`, i.e. `patch_ridge_forecaster`: guarded-core multi-harmonic spectral fit with late-context hinge/dip kernels, plus a self-validating cross-series patch-ridge refinement), by the spectral backbone alone (`forecaster="spectral"`), or by any callable `fn(context, horizon) -> (N, horizon)` you supply. Whatever the forecaster, it only ever sees the strictly pre-anchor context.
2. **Out-of-fold control-group trajectory adapter $\hat{y}^{\mathrm{CF}}_{it}(0)$**: $K$-fold cross-fitted across never-treated units $\mathcal{I}_0$ with out-of-fold variance-reduction shrinkage so control trajectories never use their own post-$t_0$ realizations and noisy corrections shrink automatically toward zero. The result is `feats.point_forecast` (always available).
3. **Trajectory summaries ($u_i$) and selection score ($\hat{\pi}_i$)**: a 2-column summary $u_i \in \mathbb{R}^2$ of unit-specific pre-anchor drift/break amplitude and post-anchor cyclical dispersion (`feats.unit_features`) and a cross-fitted score $\hat{\pi}_i = \widehat{\Pr}(E_i = 1 \mid X_i, u_i)$ (`feats.ps`) for the prognostic forest $\mu(X_i, u_i, \hat{\pi}_i)$ (Hahn, Murray, and Carvalho, 2020), keeping $\mu$ inputs strictly unit-constant (`collapse_gamma_mu = True`) and passing only baseline covariates `x_trt = x` to the moderation forest $\nu(X_i)$.
4. **Automatic pre-treatment gating (`feats.gate_open`, `feats.augment_x(x)`, `feats.recommended_ps`, `feats.baseline_forecast`)**: `feats.beyond_poly_r2` measures the out-of-fold fraction of non-polynomial post-$t_0$ control-trajectory variance explained by the forecast. The gate opens if and only if $L_{\mathrm{ctx}} \ge 12$ and `beyond_poly_r2 >= 0.05`. Using the gate-aware accessors `x_mu = feats.augment_x(x)`, `ps = feats.recommended_ps`, and `baseline_forecast = feats.baseline_forecast` activates all three channels together when the gate is open and **reverts identically to default `LongBet`** (`x_mu = x`, `ps = None`, `baseline_forecast = None`) when the gate is closed.

```python
from longbet import LongBet, LongBetConfig, build_causal_features

feats = build_causal_features(y=y, z=z, x=x, y_pre_history=y_hist)
x_mu = feats.augment_x(x)  # appends feats.unit_features when feats.gate_open else returns x

model = LongBet(LongBetConfig(num_chains=4, num_burnin=2000, num_sweeps=1000))
model.fit(
    y=y, x=x_mu, x_trt=x, z=z, t=t,
    ps=feats.recommended_ps,
    baseline_forecast=feats.baseline_forecast,
)
pred = model.predict(
    x=x_mu, x_trt=x, z=z, t=t,
    ps=feats.recommended_ps,
    baseline_forecast=feats.baseline_forecast,
)
```

```r
feats <- build_causal_features(y = y, z = z, x = x, y_pre_history = y_hist)
x_mu  <- feats$augment_x(x)

fit  <- longbet(y = y, x = x_mu, x_trt = x, z = z, t = 1:Tn,
                ps = feats$recommended_ps,
                baseline_forecast = feats$baseline_forecast)
pred <- predict(fit, x = x_mu, x_trt = x, z = z, t = 1:Tn,
                ps = feats$recommended_ps,
                baseline_forecast = feats$baseline_forecast)
```

#### Defaults of the Built-in Forecaster (`late_context`) and Propagating Forecast Uncertainty (`n_forecast_draws`)

- **`late_context`** (`"fixed"` default, `"none"`, or a `LateContextSpec`): controls how the built-in forecaster models structural trend breaks and pre-adoption dips (Ashenfelter's dip) inside the late pre-treatment window.
  - `"fixed"` (**default**): excludes the last `guard = 4` pre-anchor periods from the secular/harmonic fit, then adds a structural-break hinge kinked `break_kink = 10` periods before $t_0$ and a Gaussian pre-adoption dip kernel centred `transient_centre = -1.5` periods before the first post-anchor period with width `transient_width = 1.8`, with per-unit amplitudes projected onto the dominant cross-sectional harmonic factor. The kernel *locations* are an explicit, analyst-set identification restriction (analogous to an event-study anticipation window), and the post-anchor recovery of the dip follows the kernel (an explicit mean-reversion assumption). They are parameters of `LateContextSpec`, so the sensitivity of an estimate to that restriction is reported by re-running with the kink or dip moved, e.g. `late_context=LateContextSpec(break_kink=7, transient_centre=-2.5, transient_width=3.0)`, or with a component switched off (`break_mode="none"` / `transient_mode="none"`). The defaults reproduce the forecaster of the accompanying paper bit for bit.
  - `"none"`: secular trend and multi-harmonic seasonal backbone only (no hinge, no dip). Use when domain knowledge rules out pre-adoption anticipation and late breaks, or as a sensitivity bound.
  - **Why there is no data-adaptive mode.** An earlier release (0.1.x) estimated the kink date and the dip centre/width by joint profile least squares with a Mallows $C_p$ decision. In the paired sweep of the accompanying paper (R = 64 panels per cell; archived in the paper repository as `repro/out/archive/v2_engine_23a4680/`) it lowered the treated forecast RMSE but did not remove the treated-minus-control post-anchor forecast-error differential (the quantity the forecast channel turns into ATT bias). On the design whose break and dip sit exactly where the fixed kernels assume, the search made that differential *worse* than using no kernels at all (+0.18 vs +0.13, against +0.07 for the fixed kernels), because a 0.02-SD-per-period slope break and a 0.4-SD dip overlap inside the last ten pre-anchor periods and cannot be separated at typical panel sizes, so the search selects on pre-treatment noise -- noise that co-moves with selection into treatment in exactly the designs that motivate a dip kernel. When the kernels were misplaced it narrowed the forecast-level differential but did not restore ATT coverage in full MCMC re-runs (47-72% vs 59-91% for the fixed kernels) and raised the ATT bias in three of five misplacements. An option that is never preferable is not offered. The profile machinery survives only as the descriptive diagnostic `late_context_evidence(y, z, y_pre_history=...)`, which returns the treated-minus-control residual in the guard window (a pre-adoption-dip check in the spirit of a pre-trend plot), what the fixed kernels attribute to treated units, and the $C_p$ curve/surface over kink dates and dip placements -- evidence to report, not a rule the forecaster acts on.
- **`n_forecast_draws = M`** (default `0`): when set to $M \ge 2$ (e.g., $M = 8$), returns `feats.baseline_forecast_draws` of shape `(M, N, T)` drawn from the forecaster's predictive distribution (joint moving-block bootstrap of pre-anchor context residuals across units + Bayesian-ridge posterior draws of the control-adapter weights), automatically calibrated against never-treated units' double-demeaned post-anchor forecast error (`feats.forecast_draws["calibration"]["kappa"]`). Fitting `LongBet` across the $M$ forecast draws and pooling the resulting posterior draws yields a modular posterior that propagates both estimation and control-observable misspecification uncertainty of the forecast into ATT and CATT credible intervals. This describes the default `draw_method="bootstrap"`; `draw_method="posterior"` returns conjugate Normal-inverse-gamma posterior draws of the patch-ridge and control-adapter weights without the `kappa` inflation. `fit_rank_k_longbet` needs no draws from the builder: its default propagated variant draws its own posterior panels on demand, and the plug-in variant uses none.

---

## 3. `LongBetMulti`: Seemingly Unrelated Panel Outcomes

### Model

Multiple outcomes observed on the same panel share the covariate design, adoption schedule, and calendar-time basis while maintaining outcome-specific forests, trend coefficients, exposure trajectories, and unit intercepts. With `sur = True` (default), a triangular Seemingly Unrelated Regression (SUR) factorization couples the standardized/latent innovations:

```math
Y^{*(m)}_{it} = f^{(m)}_{it}
             + \sum_{j < m} \Gamma_{mj}\big(Y^{*(j)}_{it}-f^{(j)}_{it}\big) + \epsilon^{(m)}_{it},
\qquad \epsilon^{(m)}_{it} \sim \mathcal{N}(0, \sigma_m^2),
```

with $\Gamma_{mj} \sim \mathcal{N}(0, \text{`sur\_prior\_var`})$ updated by analytically marginalizing out the unit random intercepts $\gamma^{(m)}_i$ before drawing $(\Gamma_m, \sigma_m^2)$, followed by an exact Gaussian refreshment of $\gamma^{(m)}_i$. Posterior draws are aligned across outcomes, supporting joint counterfactual probability queries via `effect_draws()` and `joint_prob()`.

### Python

```python
from longbet import LongBetMulti, LongBetConfig, effect_draws, joint_prob, outcome_correlation

y = {"revenue": revenue_panel, "churn": churn_panel}
model = LongBetMulti(LongBetConfig(num_chains=4, num_burnin=2000, num_sweeps=1000))
model.fit(y=y, x=x, z=z, t=t, outcome={"revenue": "continuous", "churn": "binary"})

pred = model.predict(x=x, z=z, t=t, summary_only=False)
tau_revenue = effect_draws(pred, "revenue")
tau_churn = effect_draws(pred, "churn")
p_win_win = joint_prob(pred, {"revenue": lambda e: e > 0, "churn": lambda e: e < 0})
R = outcome_correlation(model)
model.save("multi.npz"); model = LongBetMulti.load("multi.npz")
```

### R

```r
fit <- longbet_multi(y = list(revenue = revenue_mat, churn = churn_mat), x = x, z = z, t = 1:Tn,
                     outcome = c(revenue = "continuous", churn = "binary"),
                     num_chains = 4, num_burnin = 2000, num_sweeps = 1000)
pred <- predict(fit, x = x, z = z, t = 1:Tn, summary_only = FALSE)
eff_revenue <- effect_draws(pred, outcome = "revenue")
eff_churn   <- effect_draws(pred, outcome = "churn")
p_win_win   <- joint_prob(pred, list(revenue = function(e) e > 0, churn = function(e) e < 0))
R <- outcome_correlation(fit)
```

---

## 4. `LongBetEncourage`: Longitudinal Encouragement Designs

### Model

Given an individually randomized encouragement offer $Z_{it}$ and an endogenous absorbing adoption indicator $D_{it}$, `LongBetEncourage` fits two LongBet equations and composes their posterior draws:

1. **Outcome equation on the adoption clock ($S^D_{it}$):**
   ```math
   Y_{it} = \mu_Y(X_i) + \eta^Y_{it} + \beta^Y_{S^D_{it}}\, \nu_Y(X_i)\, \mathbf{1}\{S^D_{it} \ge 1\} + \gamma^Y_i + \epsilon^Y_{it}.
   ```
2. **Discrete-time adoption hazard on the offer clock ($S^Z_{it}$):** on the at-risk set $\{D_{i,t-1} = 0\}$,
   ```math
   \lambda_{it} = \Pr(D_{it} = 1 \mid D_{i,t-1} = 0)
   = \Phi\!\Big(\mu_D(X_i) + \eta^D_{it} + \beta^D_{S^Z_{it}}\, \nu_D(X_i)\, \mathbf{1}\{S^Z_{it} \ge 1\}\Big).
   ```

Each moderator is a horseshoe linear term plus trees, $\nu(X_i) = \widetilde{X}_i'\delta + \sum_k g_k(X_i)$, the exposure trajectory's scale is learned ($\beta \sim \mathcal{GP}(0, s^2K)$, $s \sim C^+(0,1)$), a continuous outcome has stationary AR(1) errors, and the hazard is population-averaged (no unit frailty); these are fixed parts of the model, and draws are thinned by 2 unless `n_skip` is set. `LongBetEncourage.model_config()` shows the configuration a fit will use. See [the encouragement guide](docs/encouragement-guide.md).

Composing the induced adoption-timing distribution $\pi_a(x, z)$ with the dynamic treatment response $\tau_t(x, s)$ yields intention-to-treat outcome and adoption trajectories ($\mathrm{CITT}_Y$, $\mathrm{CITT}_D$), dynamic Wald/CACE summaries, principal-stratum weights, breakeven probabilities, and budget-constrained knapsack allocation rules.

### Python

```python
from longbet import LongBetEncourage

fit = LongBetEncourage(num_chains=4, num_burnin=2000, num_sweeps=1000)
fit.fit(y, d, z, x, t=t)

pred = fit.predict(groups=segment)
print(pred.table)
print(pred.reference)
print(pred.stability())

cond = fit.predict_conditional(x=new_cohort)
cond.citt_y.mean, cond.citt_d.mean, cond.cace.median
cond.principal_strata()
cond.breakeven_probability(cost=20.0)
cond.knapsack_policy(cost=20.0, capacity=80, budget=1600, hurdle=0.5)
```

### Continuous treatment intensity (`treatment="intensity"`)

When the encouraged behaviour is a continuous, non-absorbing dose $D_{it} \ge 0$ (hours; quarterly training-months that start, vary and stop), there is no adoption clock. `LongBetEncourage(treatment="intensity")` replaces it with a **distributed-lag kernel on the calendar clock**, fitted in two stages with independent PRNG streams:

```math
D_{it} = \mu_D(X_i) + \eta^D_{it} + \beta^D_{S^Z_{it}}\,\nu_D(X_i)\,Z_{it} + \gamma^D_i + v^\perp_{it},
\qquad
Y_{it} = \mu_Y(X_i) + \eta^Y_{it} + \nu_Y(X_i)\sum_{\ell=0}^{L}\kappa_\ell\,D_{i,t-\ell} + \mathrm{CF}_{it} + \gamma^Y_i + \epsilon_{it},
```

with $\kappa \sim \mathcal{GP}(0, s^2 K_L)$ over the lag index ($L = T - t_0$ by default), so sign-changing kernels (in-treatment lock-in $\kappa_0 < 0$, later returns $\kappa_{\ell \ge 1} > 0$) are representable; the binary model is the special case $\kappa_\ell = \tau(\ell+1) - \tau(\ell)$. Offer effects are composed draw by draw, $\mathrm{CITT}_Y(x,t) = \nu_Y(x)\sum_{\ell=0}^{t-t_0}\kappa_\ell\,\mathrm{CITT}_D(x,t-\ell)$. `identification` selects what Stage 2 learns the kernel from: `"deconvolution"` (T0: the propensity-centered Stage-1 offer-induced dose increment $(Z_i - \hat e_i)\,\widehat{\mathrm{CITT}}_D(X_i,t-\ell)\,\mathbb{1}\{t-\ell \ge t_0\}$ with unconstrained organic-dose nuisance controls, immune to within-arm dose endogeneity and measurement error), `"control_function"` (T1, default: the observed dose plus the 3-term plug-in control function $[\hat v_{it}, \sum_{\ell\ge1}\hat v_{i,t-\ell}, \bar{\hat v}_i] \otimes [1, \tilde X_i]$ of the Stage-1 residuals, drawn in the conjugate block) or `"unconfounded"` (T2: the observed dose alone). `intensity_reference()` is the model-free T0 benchmark (stratified ITTs, Toeplitz deconvolution with delta-method SEs and $\operatorname{cond}(T_d)$, exact Fieller/AR sets for the cumulative return per unit dose, joint AR region, ECDF/FOSD check, dose-history-conditioned $z_{\mathrm{dip}}$), and every prediction reports $z_{\mathrm{T1-T0}}$ / $z_{\mathrm{T2-T0}}$ discrepancies against it. See [the intensity guide](docs/encouragement-intensity-guide.md).

```python
from longbet import LongBetEncourage, intensity_reference

ref = intensity_reference(y, d, z, x=x, strata=site)          # design-based T0, no sampler
ref.kernel, ref.itt[["period", "itt_d", "itt_y", "omega_cum", "omega_cum_set_type"]], ref.condition_number

fit = LongBetEncourage(treatment="intensity", identification="control_function",
                       num_chains=4, num_burnin=2000, num_sweeps=1000)
fit.fit(y, d, z, x, t=t, strata=site)
pred = fit.predict(groups=segment)
pred.effects(), pred.kernel(), pred.discrepancy(), pred.diagnostics(), pred.cumulative_return()
cond = fit.predict_conditional(x=new_cohort)
cond.citt_y.mean, cond.kernel.mean, cond.omega_cum.median, cond.knapsack_policy(cost=20.0, capacity=80)
```

---

## 5. Key Configuration Options (`LongBetConfig`)

| Option | Default | Description |
|---|---|---|
| `num_chains`, `num_burnin`, `num_sweeps` | `4`, `2000`, `1000` | Number of parallel MCMC chains, warm-up sweeps, and retained posterior sweeps per chain. |
| `num_trees_pr`, `num_trees_trt` | `20`, `60` | Number of trees in the prognostic ($\mu$) and treatment moderation ($\nu$) ensembles. |
| `use_trend_block` | `True` | Include the conjugate calendar-time basis $\eta_{it}$ (centered DCT-II period effects + polynomial $\times$ covariate interactions) drawn jointly with $(\gamma, \beta)$. |
| `trend_degree`, `trend_num_features` | `2`, `64` | Polynomial degree in calendar time and number of standardized nonlinear hinge, bilinear, hinge-interaction, and random Fourier features in $b(X_i)$. |
| `trend_horseshoe` | `True` | Place a global–local horseshoe shrinkage prior on the covariate-by-time trend interactions $c_k$. |
| `baseline_cells` | `"all"` | Cells that inform the untreated baseline $(\mu, \eta, \gamma, \sigma^2)$: `"all"` (single-likelihood Gibbs; treated cells enter through $Y_{it} - \beta_{S_{it}}\nu(X_i)$) or `"untreated"` (two-module *cut* posterior: baseline fitted from $Z_{it}=0$ cells only, $(\beta,\nu)$ conditioned on it; continuous outcomes, `tempering_levels = 1`, `LongBet` only). |
| `kernel_type`, `sig_knl`, `lambda_knl`, `sigma_m` | `"se"`, `1.0`, `1.0`, `1.0` | Exposure GP covariance kernel (`"se"`, `"matern32"`, `"matern52"`, or `"ar1"`), kernel SD, lengthscale, and common-level SD. |
| `random_intercept` | `True` | Include unit random intercepts $\gamma_i \sim \mathcal{N}(0, \sigma_\gamma^2)$, analytically marginalized during unit-constant forest updates and jointly sampled with $(c, \beta)$. |
| `outcome`, `num_categories` | `"continuous"`, `None` | Outcome distribution (`"continuous"`, `"binary"`, or `"ordinal"` with `num_categories = K`). |

### What the Defaults Do and Why (`baseline_cells` and `trend_num_features`)

In the paired $R = 64$ baseline-module factorial of the accompanying paper (`repro/src/run_baseline_cells_factorial.py`), the two baseline-module options (`baseline_cells` $\in \{\text{`"all"`}, \text{`"untreated"`}\}$ and `trend_num_features` $F \in \{64, 0\}$) were evaluated across polynomial-trend (`S1`), unit-specific seasonal (`S2`), and dynamic-confounding (`S3`) designs under both non-null and $\tau \equiv 0$ effects:

- **Why `baseline_cells = "all"` is the default.** Under `baseline_cells = "all"`, `trend_block_step` samples the trend coefficients $c$, unit intercepts $\gamma_i$, and exposure trajectory $\beta_S$ from their exact joint Gaussian conditional across all observed cells. Under `baseline_cells = "untreated"`, the baseline module $(\mu, \eta, \gamma, \sigma^2)$ is restricted to untreated cells ($Z_{it} = 0$), so treated cohorts contribute only their pre-adoption periods ($t < G_i$) to $\eta$ and $\gamma_i$ and their post-adoption baseline is an out-of-sample polynomial extrapolation over $t \ge G_i$. When adoption is selective on $X_i$ and pre-adoption windows are short ($3\text{--}5$ periods), horseshoe shrinkage on the pre-adoption slope estimate pulls the extrapolated post-adoption baseline toward the grand trend, **increasing** ATT bias on `S1` from $+0.086$ to $+0.132$ under non-null effects (and from $+0.024$ with 97% coverage to $+0.127$ with 62% coverage under $\tau \equiv 0$) and on `S2` without the forecast channel from $-0.148$ to $+0.334$. Use `baseline_cells = "untreated"` only as a modular sensitivity check when the in-panel pre-adoption window is long or when an external `baseline_forecast` already carries the post-adoption trajectory.
- **When to use `trend_num_features = 0` vs `64` (default).** On `S1` (where true differential trends are linear in $X_i$ and adopters are concentrated in the tails of $X_i$), the $F = 64$ nonlinear dictionary columns in $b(X_i)$ compete with $\mu(X_i)$ and $\beta_S \nu(X_i)$ for post-adoption variation under horseshoe shrinkage. Setting `trend_num_features = 0` restricts $b(X_i)$ to the $p_0$ standardized linear covariates, cutting `S1` ATT bias by more than half ($+0.086 \to +0.031$ for default `LongBet`, $+0.097 \to +0.046$ for `Conjugate+PS`), restoring `S1` ATT 95% coverage from $81\text{--}83\%$ to $92\text{--}97\%$ ($100\%$ under $\tau \equiv 0$), and increasing median minimum exposure ESS by $2\times\text{--}3\times$, while remaining non-inferior across `S2` and `S3`. Conversely, when true differential trends load nonlinearly on $X_i$ (`trend_nonlin`), omitting the nonlinear dictionary introduces bias unless those trends are captured by `baseline_forecast`. Keep the default `trend_num_features = 64` when covariate-by-time trend loadings may be nonlinear, and compare against `trend_num_features = 0` when linear covariate slopes are plausible and faster mixing or lower tail-selection shrinkage bias is desired.

---

## Changes in 0.10.0 (in progress on `feat/v0.10-methodology`)

Methodology v2 for multi-arm, multi-outcome reconciliation ([docs/rank-k-multi-arm-multi-outcome.md](docs/rank-k-multi-arm-multi-outcome.md), section 5). Every changed default has an explicit opt-out that reproduces 0.9.0 bit for bit (`tests/test_rank_k_v010_parity.py` checks it against draws stored from 0.9.0). Changes relative to 0.9.0:

- **`composite_mode` defaults to `"coherent_total"`** in `reconcile_multi_arm_multi_outcome_draws` and `reconcile_nested_3arm_multi_draws` (`composite_mode=None` resolves to `"coherent_total"` when the composites have disjoint components and the prior truncates at zero, else to `"bottom_up"`). `composite_mode="bottom_up"` restores 0.9.0 and now emits a `RankKConfigWarning` with the expected null-floor inflation $\sum_m \sigma_m / \sqrt{\boldsymbol{\sigma}^\top\mathbf{R}\boldsymbol{\sigma}}$ whenever the unconstrained composite posterior on the key contrast has $\Pr(\tau \le 0) > 0.10$; both results gain `composite_floor_diagnostics`.
- **Bipartite contrast graphs** (`reconcile_multi_arm_draws`, `reconcile_multi_arm_multi_outcome_draws`): the warning is now a `RankKDesignWarning` (new, exported, a `UserWarning` subclass) naming the missing odd cycle and the one extra fit that closes it; new `bipartite_policy` in `{"allocation" (default, 0.9.0), "conservative", "marginalize"}` resolves the unidentified arm-variance direction (results gain `bipartite_policy`, `null_shift_interval`, `null_shift`); `prior_mode=None` is the new default sentinel (resolves to `"bundle_halfnormal"`); and a chain order fitted without all adjacent edges (e.g. a vs-control star) is accepted with `topology="chain"`.
- **Slice sampler auto-continuation for $d \ge 4$** (`auto_continue=True` default, with `target_rhat=1.01`, `target_ess=400`, `max_iter_multiplier=8`, in `reconcile_multi_arm_draws` and `reconcile_multi_arm_multi_outcome_draws`): blocks are appended until split-$\hat{R}$ and bulk ESS meet their targets, then a `RankKConfigWarning` if the cap is hit; `SamplerDiagnostics` gains `n_iterations` and `converged`. `auto_continue=False` restores 0.9.0.

## Changes in 0.9.0

General $A$-arm ($A \ge 2$) and $M$-outcome coherent reconciliation for Rank-$K$ BCF ([docs/rank-k-multi-arm-multi-outcome.md](docs/rank-k-multi-arm-multi-outcome.md)) and completion of the `extrapolate_rank_k_decay` `rho_sd=None` default flip announced in 0.8.0. Changes relative to 0.8.0:

- **New module `longbet._rank_k_multi_arm`** (exported; general $A$-arm Mode-1 and trajectory reconciliation across chain and poset topologies): `reconcile_multi_arm_draws`, `MultiArmReconciliationResult`, `ContrastReconciliationResult`, and `SamplerDiagnostics`. Projects $P \ge A - 1$ pairwise Mode-1 posterior means onto the $(A - 1)$-dimensional cut space $\hat{\boldsymbol{\mu}} = (\mathbf{B}^\top\mathbf{B})^{-1}\mathbf{B}^\top\hat{\mathbf{y}}$; recovers arm sampling variances $\hat{\mathbf{v}} \in \mathbb{R}_+^A$ via bounded NNLS on non-bipartite contrast graphs (reproducing `_backsolve_group_variances` in closed form for the $A = 3$ triangle) or sample-size-informed ridge-regularized NNLS with a `UserWarning` when the contrast graph is bipartite; supports chain `Dirichlet(alpha)` bundle priors, unique-top poset order-polytope bundle priors via uniform linear extensions, independent truncated Normal/Half-Normal priors, and flat priors; samples the posterior via 50/50 defensive-mixture SIR for $d = A - 1 \le 3$ (with automatic fallback when $\text{ESS} < \max(100, 0.02 \, n_{\text{proposals}})$) or a vectorized multi-chain radial-cone / polyhedron hit-and-run slice sampler for $d \ge 4$; and projects spanning-tree edge trajectories so cycle identities hold to $< 10^{-12}$ across all contrasts and periods. At $A = 3$ with `topology="chain"`, default `split_alpha`, `monotone=True`, and `sampler in ("auto", "sir")`, `reconcile_multi_arm_draws` delegates directly to `reconcile_nested_3arm_draws` for bit-for-bit equality with 0.8.0.
- **Extended `longbet._rank_k_multi`** (exported; multi-arm multi-outcome coupling and composite synthesis): `reconcile_multi_arm_multi_outcome_draws` and `MultiArmMultiOutcomeReconciliationResult` generalize `reconcile_nested_3arm_multi_draws` to arbitrary $A \ge 2$ arms and $P \ge A - 1$ contrasts, coupling $M$ structural outcomes via `iman_conover_couple` on `key_contrast` (with per-arm Ledoit-Wolf ANCOVA residual correlation from `estimate_cross_outcome_correlation`) and synthesizing `composite_outcomes` under `"bottom_up"` or `"coherent_total"` with exact draw-wise additivity (`additivity_violation_rate == 0.0`, $\max|\text{gap}| < 10^{-12}$) across all contrasts and periods. `ge_revenue_draws` also accepts `MultiArmMultiOutcomeReconciliationResult` with any `(treated, ref)` contrast pair.
- **Extended decision and ROI tables in `longbet._rank_k_plots`**: `roi_break_even_table` accepts `MultiArmReconciliationResult` and gains an optional `contrast` argument to filter to a specific `(treated, ref)` pair; `extrapolate_rank_k_decay` accepts `MultiArmReconciliationResult`; and new `multi_arm_decision_table` (accepting `costs` or `arm_costs`, `step_costs` or `incremental_costs`, and `alpha` or `outcome_weights`) computes per-arm probabilities of being best (`pr_best_arm` / `p_best_net`, `pr_best_treated_arm`, `p_best_gross`), positive net value (`pr_positive_net_value` / `p_beats_control`), spanning-tree incremental step probabilities (`pr_incremental_step_worthwhile` / `p_step_roi_positive`, `p_positive_step`), `expected_net_value`, and net value quantiles (`net_value_p30`, `net_value_p50`, `net_value_p70`).
- **Completed deprecation (`extrapolate_rank_k_decay`)**: As announced in 0.8.0, the default of `rho_sd` in `extrapolate_rank_k_decay` is now `None` (no assumed `rho ~ Beta` scenario band unless a positive `rho_sd` is passed explicitly), and the `FutureWarning` when omitting `rho_sd` has been removed.

## Changes in 0.8.0

In-sample persistence diagnostics and required-persistence decisions ([docs/rank-k-stability-and-required-persistence.md](docs/rank-k-stability-and-required-persistence.md)). No estimator, prior or sampler changed: every fit, reconciliation and aggregate of 0.7.0 is reproduced bit for bit. Changes relative to 0.7.0:

- **New module `longbet._rank_k_multi`** (exported; multi-outcome nested 3-arm reconciliation and copula coupling): `estimate_cross_outcome_correlation` (within-arm ANCOVA residual correlation across outcomes with analytical Ledoit-Wolf shrinkage and eigenvalue flooring), `iman_conover_couple` (marginal-preserving Iman-Conover rank-copula coupling across outcomes), `reconcile_nested_3arm_multi_draws` (`Nested3ArmMultiReconciliationResult`: per-outcome 2-D SIR nested 3-arm reconciliation coupled via Iman-Conover with `"bottom_up"` and `"coherent_total"` composite synthesis satisfying exact draw-wise additivity across contrasts and components), `ge_revenue_draws` (coupled general-equilibrium weighted linear combinations across outcomes), `check_rank_k_convergence` (`RankKConvergenceCheck` per-mode split R-hat and bulk ESS gates), and `paid_winsor_cap` (strictly pre-adoption positive-part quantile cap for sparse/zero-inflated outcomes); `roi_break_even_table` in `longbet._rank_k_plots` also gains `component_cum_draws` and `component_weights` for vector-`\alpha` multi-outcome break-even evaluation.
- **New module `longbet._stability`** (exported; pure numpy / scipy / pandas on posterior draws and panels, no refits): `wearout_ancova` (design-based early-versus-late ANCOVA on the balanced subpanel with strata fixed effects, HC1 standard errors and the minimum detectable difference next to every estimate), `estimate_post_pulse_retention` (per-draw monthly retention by weighted least squares on levels, robust to exact zeros, with the share of draws piled up on each `rho` bound and a degenerate-window policy), `unit_slope_correlation` (within-arm correlation across outcomes of per-unit slopes, Fisher-z pooled), `pool_multi_outcome_retention` (DerSimonian-Laird pooling of `lambda = -ln(rho)`, draw-preserving Morris shrinkage of a focal outcome, correlation-adjusted, Q-moment, leave-one-outcome-out and equal-weight variants and the break-even correlation), `RetentionRuleSpec` / `evaluate_retention_rule` (a pre-registerable gate with a literal and a robust reading), `diagnose_trajectory_stability` (`TrajectoryStabilityReport`, whose `verdict` is always `None`) and `pseudo_horizon_separation` (sizes a leak-free hindcast from existing draws), plus `StabilityWarning`.
- **New decision layer** in `longbet._rank_k_plots`: `compute_required_persistence` (`RequiredPersistenceResult`: the required persistence `X` along Bound U, break-even retention, payback CDFs and cumulative value under Bound U, Bound L and the `rho = 0` floor, optional `extra_regimes` with `headline` / `diagnostic` roles, split incrementality multipliers through `alpha_share`, survival and discount weights, a placebo decision-immunity gate and a data-health certificate; the identity `Pr(X <= k) = Pr(payback under U <= L + k)` and the ordering of the regimes are asserted whenever every arm is non-negative), `plot_required_persistence_and_payback`, `weibull_survival_weights` and `discount_weights`.
- **`extrapolate_rank_k_decay`** gains `pulse_decay_tau`, `plateau_arm`, `pulse_arm`, `total_arm` and `survival`, and accepts a `Nested3ArmReconciliationResult`: the plateau decays at `rho^h`, the pulse fades at `exp(-h / tau)`, both are weighted by `S(L + h) / S(L)`, and the total arm's tail is their sum, so `AC = BC + AB` holds draw by draw. Without these arguments the 0.7.0 draws and tables are reproduced bit for bit (plus the new `rho_band_label` column). `DecayExtrapolationResult` gains `scenario_note`, `pulse_decay_tau`, `survival` and `arm_roles`, and `cumulative_summary_df` gains `rho_band_label`.
- **Deprecation.** The hierarchical `rho ~ Beta` band of `extrapolate_rank_k_decay` is an assumed scenario, not a posterior interval, and is now labelled as such in its tables and in `plot_rank_k_decay_extrapolation`. Calling `extrapolate_rank_k_decay` without `rho_sd` emits a `FutureWarning`; the default changes from `0.07` to `None` (no band) in 0.9.0. Pass `rho_sd` explicitly to keep the band, or `rho_sd=None` to drop it.
- README: the three-arm example passes the explicit `prior_scale` that `reconcile_nested_3arm_draws` requires.

## Changes in 0.7.0

Masked event-time Rank-$K$ fits: population effects under administrative censoring ([docs/rank-k-masked-event-time.md](docs/rank-k-masked-event-time.md)). Changes relative to 0.6.0:

- **`fit_rank_k_longbet` accepts unobserved cells.** `y` may hold `NaN` (or pass `mask=`, an `(N, T)` boolean that must agree with the `NaN` pattern): each unit's observed post-adoption cells must be a contiguous prefix `s = 1..L_i` of the event-time horizon (administrative censoring) and at least one pre-period cell must be observed; masked panels need a single adoption period and the JAX engine (the plug-in `stochtree` engine, `fit_rank_k_hybrid` and `rank_k_placebo_audit` raise a clear error). New `horizon=` argument (default: the longest observed window; append all-`NaN` columns to predict beyond the panel).
- **Common event-time basis** (`RankKLongBetConfig(common_basis=None | True | False)`; `None` means `True` for masked and single-wave panels, the paper's per-cohort bases otherwise): one orthonormal basis on `s = 1..H` for every unit, `res.bases = {t0: Psi}`, `res.horizon = H`. Staggered panels without a mask keep the per-cohort bases unless `common_basis=True` (then the later cohorts' Mode 1 is corrected to the common horizon as well).
- **Two-pass projection** (`projection="two_pass"`, default; `"restricted_wls"` for the joint per-unit fit with a condition-number guard): higher modes on the units whose window supports them (`min_support_modes=(1, 6, 8)`; the per-unit OLS slope on `{1, psi_k}`, so the level does not leak into a partial window's projection; unsupported units get precision weight 0 and a cross-sectional posterior prediction), Mode 1 on every unit as the observed-window mean of the residual after the fitted treatment dynamics are removed. `mode_coupling="propagated"` (default) runs Mode 1 on `n_forecast_draws` tracks with one higher-mode posterior draw each (paired with the forecast draws when present; `n_components` and `att_var_between` then count these components), `"plug_in"` uses the posterior means, `"none"` skips the correction. Per-unit precision weights from the observed windows (`1 / n_pre_i` pre-anchor term). The observed window length `L_i` is appended to the forest covariates (`support_covariate=True`).
- **Population weighting**: `att_weighting="cell"` (default, the paper's `L_i / |O_trt|`) or `"unit"` (`1 / N_trt`); `att_draws`, `att_hat`, `att_se`, `att_ci_*` honour it, both engines. `RankKLongBetResult` gains `unit_weights`, `att_weighting`, `cohorts`, `catt_observed`, `horizon`, `support_mask`, `post_window_length` and the `treated_units` property; `catt_draws` holds the posterior-predictive cell effect on every post cell `s <= H` of every treated unit. `reconcile_nested_3arm_draws` reads the stored cohorts and unit weights and follows `res.horizon`.
- **New** `cumulative_effect(res, horizon, weighting="unit")`, `event_time_trajectory(res, weighting="unit")` and `rank_k_transport_audit(res, strata, overlap_horizon=5)` (`TransportAuditResult`: posterior difference of the two strata's early event-time profiles with P5/P30/P50/P70/P95, a flag, and the share of the cumulative effect carried by unobserved cells).
- **Sampler.** `sample_rank_k_mode_jax` takes `dy_k_draws` (per-track Mode-1 targets); units with precision weight 0 no longer count in the shape parameter of the inverse-gamma variance updates (`N_ctrl_dof`, `N_trt_dof`, `N_all_dof`) nor in the outcome scale; on complete panels nothing changes.
- **Bit-for-bit compatibility.** On fully observed panels the correction is exactly zero and the common-basis path reproduces 0.6.0 bit for bit (`tests/test_rank_k_masked.py`, reference `tests/data/masked_reference_6f714d5.npz` generated from the unmodified 0.6.0 sources with `tests/data/make_masked_reference.py`); the per-cohort path is unchanged. Under a mask the ATT depends on the fitted dynamics (`rank_k`, basis family) through the correction: report a `rank_k` sensitivity.

## Changes in 0.5.0

Release hardening of the Rank-$K$ estimator (`fit_rank_k_longbet`, `RankKLongBetConfig`, `rank_k_placebo_audit`) in Python and R. Changes relative to 0.4.0:

- **The default remains the "propagated" variant** (`uncertainty_propagation=True`, engine `"jax"`), now pinned to the companion paper's validated benchmark budget: 2 chains instead of 4, 100 burn-in sweeps instead of 200, `n_forecast_draws=8` instead of 16 without adaptive expansion (`max_forecast_draws` 32 → `None`, meaning `n_forecast_draws`), and 50 / 30 / 15 trees instead of 75 / 35 / 20, that is $2 \times 8 \times 25 = 400$ posterior draws instead of $4 \times 16 \times 25 = 1600$. Select the fast "plug-in" variant with `uncertainty_propagation=False` (or, as before, `engine="stochtree"`), which keeps the 0.4.0 plug-in budget (one chain per mode, 10 grow-from-root sweeps, 200 burn-in sweeps, 600 retained draws, 75 / 35 / 20 trees). `stochtree>=0.4.5` is now a core dependency; the `[rank-k]` extra is kept, empty.
- **Variant-dependent defaults.** `engine` and `uncertainty_propagation` default to `None`: `uncertainty_propagation` selects the variant (`None` means `False` only with `engine="stochtree"`, and `True` otherwise) and the engine follows from it. `num_chains`, `num_burnin` and the three tree counts also default to `None` and are resolved at construction to the budget validated for the selected variant, so a constructed configuration holds concrete values. For small panels use the application budget `RankKLongBetConfig(num_chains=4, num_burnin=300, draws_per_component=100)`. `dataclasses.replace` keeps the resolved engine and budget: replacing `uncertainty_propagation` alone raises `ValueError`, and replacing it together with `engine=None` keeps the other variant's budget, which triggers a `RankKConfigWarning`; build a new configuration to switch variants. `max_forecast_draws=None` stays `None`, so `replace(cfg, n_forecast_draws=...)` never turns the expansion on.
- **Contradictory variant settings raise `ValueError`:** `engine="jax"` with `uncertainty_propagation=False` (the JAX sampler without forecast draws is only a synthetic-study analogue of the plug-in variant, and its intervals under-covered, 60–84%, under spatially correlated shocks) and `engine="stochtree"` with `uncertainty_propagation=True` (not implemented; 0.4.0 ignored `uncertainty_propagation` with `engine="stochtree"`). The computation of the former remains available as `uncertainty_propagation=True, n_forecast_draws=1`, which triggers a `RankKConfigWarning`.
- **`warm_burnin` is deprecated** (default 50 → `None`). It was never a warm start: when forecast draws were propagated, every mode-1 track ran `max(num_burnin, warm_burnin)` burn-in sweeps. Setting it now triggers a `RankKConfigWarning`; set `num_burnin` instead. Fits with `num_burnin < 50` that relied on the old default now run `num_burnin` mode-1 burn-in sweeps instead of 50.
- **Removed:** `RankKLongBetConfig.use_forests` (a verification switch that replaced both forests by constants; passing it now raises `TypeError`) and the `moran_i_residual` / `moran_z_residual` fields of `PlaceboAuditResult` (never computed, always 0).
- **Validation.** 0.4.0 did not validate `RankKLongBetConfig`: an unknown `engine` or `baseline_cells` failed only when fitting, the samplers cast counts with `int()` (which truncates `2.5` to 2), and `num_chains < 1`, `draws_per_component < 2` and `max_forecast_draws < n_forecast_draws` were silently raised to their minimum. `RankKLongBetConfig` now raises `ValueError` at construction for an invalid type or value: counts must be whole numbers in range (for example `rank_k < 1`, a negative `n_forecast_draws`, `max_forecast_draws < n_forecast_draws`, `num_chains < 1`, `draws_per_component < 2`, `num_mcmc < 2`, tree counts below 1, `num_threads` 0 or below -1, a negative `random_seed`); flags must be `True` or `False` (0 and 1 are accepted; `uncertainty_propagation` may also be `None`); `theta_prior_mean`, `theta_prior_sd`, `ridge_eps` and `between_var_threshold` must be finite numbers, with `theta_prior_sd > 0`, `ridge_eps >= 0` and `between_var_threshold` in $[0, 1]$; `theta_f_clip` must be a `(low, high)` pair with `low <= high`, `lengthscale_grid_km` a non-empty sequence of positive, finite ranges, and `engine` and `baseline_cells` known values. `None` is accepted only by the fields that default to `None`, and strings and booleans are rejected in numeric fields. Whole-number floats and NumPy scalars are stored as Python `int` / `float` and sequences as tuples, so a configuration stays hashable. `fit_rank_k_longbet` and `rank_k_placebo_audit` check every input before any computation and raise `ValueError` for non-finite or misshapen arrays, a `z` that is not 0/1 or not absorbing, a panel without treated or never-treated units or with units treated in period 0, `cohorts` that disagree with `z` (never-treated units need `cohorts[i] = T`; values above `T` are now rejected), `t_0_panel` outside `[1, first adoption period]`, `features` built on another panel, `ps` outside $[0, 1]$, `coords` outside the latitude/longitude ranges, and a `pseudo_anchor_offset` that is not a whole number (0.4.0 truncated it).
- **New `longbet.RankKConfigWarning`**, emitted only when an explicit input or a non-default option would be silently ignored or contradicts another one, for example: with the plug-in variant, `coords`, `cluster_ids`, `pre_anchor_loadings`, `baseline_forecast_draws` or a non-default propagated-only option (forecast-draw, chain, slope-prior, factor-loading, spatial-field and `baseline_cells` options); with the propagated variant, a non-default `num_gfr`, `ridge_eps` or `num_threads`, `num_mcmc` while forecast draws are propagated, `uncertainty_propagation=True` with `n_forecast_draws` 0 or 1, fewer `baseline_forecast_draws` than it uses, or a non-default `n_factor_loadings` that differs from the column count of the supplied loadings; both `coords` and `cluster_ids`; slope options without a forecast; a configuration that holds the default budget of the other variant (as `dataclasses.replace(cfg, engine=None, uncertainty_propagation=...)` leaves it). It is also emitted when the deprecated `warm_burnin` is set and when a cohort has fewer post-adoption periods than `rank_k`. Choosing a supported variant or option never warns. `rank_k_placebo_audit` also warns that it ignores `features` (it always rebuilds them on the placebo window, now with the configured `n_factor_loadings`; 0.4.0 always used 3).
- **Robust spatial kriging.** When the Cholesky factorisation of the never-treated units' spatial kernel fails (long ranges on spread-out, for example global, unit sets, where the Matérn-3/2 correlation on great-circle distances can be indefinite), the kriging operator falls back to an eigenvalue-floored pseudo-inverse instead of raising; results are bit-identical whenever the factorisation succeeds.
- **R wrapper.** The defaults of `fit_rank_k_longbet()` now mirror `RankKLongBetConfig()`, with `NULL` (resolved in Python) for the variant-dependent arguments `engine`, `uncertainty_propagation`, `max_forecast_draws`, `num_chains`, `num_burnin` and the three tree counts, so the R default is the propagated variant with its validated benchmark budget. 0.4.0 used `engine = "jax"`, `uncertainty_propagation = TRUE`, `n_forecast_draws = 8L`, 4 chains, 150 burn-in sweeps, `num_mcmc = 100` and 50 / 30 / 15 trees, and left `max_forecast_draws` at the then Python default of 32, so the adaptive expansion was on. New arguments `max_forecast_draws`, `draws_per_component`, `pre_anchor_loadings` and `cohorts`, added after `num_threads` so that the 0.4.0 argument positions are kept (pass the 0.5.0 arguments by name); `t_0_panel` is a single integer (it was documented and passed as a history panel); logical flags are validated (`NA` and vectors are errors), and counts must be single whole numbers (0.4.0 passed them through `as.integer()`, which truncates `2.5` to 2); a `features` list restored with `readRDS()` supplies only the stored arrays that the resolved variant uses (its forecast draws only when draws are propagated, that is with a forecast and `n_forecast_draws >= 2`); the result now contains `mode_tau_draws`, `theta_f_draws`, `s_zeta_modes`, `ell_modes`, `n_components`, `rhat_att` and `ess_att`, and `att_within_var` / `att_between_var` read the right attributes (`spatial_var_share_modes` is gone). `rank_k_placebo_audit()` passed a `placebo_horizon` argument that the Python function does not accept; it now takes `pseudo_anchor_offset` and `cohorts`, forwards configuration arguments through `...` (so `rank_k` defaults to 2 instead of a hard-coded 1; the placebo ATT is the same for every `rank_k`) and returns `flag_parallel_trends_violation` instead of Moran fields. `build_causal_features()` defaults to `n_forecast_draws = 0L` (was `8L`; the propagated variant draws its own forecast panels and the plug-in variant uses none), gains `draw_method`, `forecast_draw_block_length` and `n_factor_loadings`, and returns `pre_anchor_loadings`.

## Changes in 0.3.0

- `build_causal_timesfm_features` → **`build_causal_features`** and `CausalTimesFMFeatures` → **`CausalFeatures`** (the old names remain as deprecated aliases for one release). The result gained a public `point_forecast` (the ungated forecast consumed by `fit_rank_k_longbet`).
- The pretrained-TimesFM loader (`load_timesfm_forecaster`), the `[timesfm]` extra, and the TimesFM-object forecaster protocol were removed. The only external-forecaster protocol is a callable `fn(context, horizon) -> (N, horizon)`.
- The default forecaster is now **`patch_ridge_forecaster`** (the backbone of every forecast arm in the paper); `forecaster="spectral"` gives the previous default.
- New: **`fit_rank_k_longbet`**, `RankKLongBetConfig`, `RankKLongBetResult`, `orthonormal_exposure_basis` (extra `[rank-k]`), in Python and R.
