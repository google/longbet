# Coherent Multi-Arm and Multi-Outcome Reconciliation for Rank-K BCF

When a longitudinal experiment compares $A \ge 3$ arms (a graded dose ladder $C \prec B_1 \prec \dots \prec A_{\text{top}}$, or a $2 \times 2$ factorial partial order) across $M \ge 1$ outcomes and linear composites, fitting each pairwise contrast $(a, b)$ and each outcome $m$ independently with `fit_rank_k_longbet` leaves four coherence gaps:

1. **Shared reference-arm covariance.** Contrasts that share an arm $b$ have non-zero sampling covariance, $\operatorname{Cov}(\hat{\tau}_{(a,b)}, \hat{\tau}_{(b,c)}) = -v_b < 0$ and $\operatorname{Cov}(\hat{\tau}_{(j,0)}, \hat{\tau}_{(\ell,0)}) = +v_0 > 0$. Treating the pairwise posteriors as independent misstates the variance of increments and sums.
2. **Cycle additivity violations.** Over any cycle of the contrast graph, independent pairwise draws violate $\tau_{(a,c)} = \tau_{(a,b)} + \tau_{(b,c)}$ on every draw and period.
3. **Order reversals and compounded truncation.** Independent priors across nested arms allow non-monotone posterior quantiles ($Q_q(\tau_{\text{Combined}}) < Q_q(\tau_{\text{Base}})$), and summing $K$ independently truncated component outcomes inflates the null composite mean.
4. **Cross-outcome correlation.** Outcomes measured on the same units covary, so monetization ($\sum_m \alpha_m \mathrm{Cum}_m$) and composite outcomes need a joint posterior, not independent marginals.

Since 0.11.0 the library closes all four gaps with **one joint hierarchical Bayesian model** fitted to the Stage-1 posterior summaries (a two-stage, meta-analytic formulation): every unknown that affects inference has a prior and is integrated out by an exact sampler. There is no Type-II maximum likelihood, no plug-in tuning constant, no IRLS, no copula reshuffling, and no post-hoc clipping, flooring or rescaling of draws. Cycle identities, composite additivity and trajectory coherence hold **exactly on every draw** because contrasts and composites are linear maps of the arm-level parameters (`max_cycle_discrepancy` $< 10^{-13}$, `additivity_violation_rate == 0.0`).

The public surface is two functions with no tuning knobs:

```python
reconcile_multi_arm_draws(results, *, arms=None, order=None, poset=None, monotone=True,
                          prior_scale=None, horizon=None, n_modes=2, seed=2026)
reconcile_multi_arm_multi_outcome_draws(results_by_outcome, *, composite_outcomes=None,
                          composite_results=None, cross_outcome_residuals=None, arms=None,
                          order=None, poset=None, monotone=True, prior_scale=None,
                          horizon=None, n_modes=2, seed=2026)
```

`reconcile_nested_3arm_draws(res_bc, res_ab, res_ac)` and `reconcile_nested_3arm_multi_draws(results_bc, results_ab, results_ac, ...)` are thin convenience wrappers for the nested three-arm design $C \preceq B \preceq A$ (`arms=("C", "B", "A")`, `order=("C", "B", "A")`); they run the same engine and return the same draws as the general functions. The zero-keyword call *is* the shipped configuration: convergence control (4 chains, split-$\hat{R} \le 1.01$, bulk ESS $\ge 400$, automatic continuation), the robust Student-$t_4$ forest error, per-period positivity, the SUR cross-outcome block and dynamic shrinkage are parts of the model, not options. `tests/test_api_hygiene_v011.py` checks that every option removed in 0.11.0 raises and that the zero-keyword call resolves to the shipped model.

---

## 1. Stage-1 inputs

Each entry of `results` is the `RankKLongBetResult` of one pairwise fit `(arm_high, arm_low)` from `fit_rank_k_longbet` on the two arms' units (the lower arm plays the role of the never-treated group). The fitted pairs must form a connected graph over the arms, use one orientation per pair, and share the post-adoption window. The reconciliation reads three things from every fit:

- **Treated-average trajectory draws** $C^{(s)}_{m,p,t}$ on event time $t = 1, \dots, T_{\text{post}}$, the treated-unit average of `catt_draws` (the 15-month `"cell"`-weighted `att_draws` are never consumed; the window Mode 1 below is a parameter of the model, not a by-product).
- **Arm metadata** (`RankKLongBetResult.arm_metadata`, attached automatically by both engines since 0.11.0, about 100 KB per fit): per-arm unit counts `n_treated` / `n_control`, the post-window residual covariance over event time of each arm (`post_residual_cov_treated` / `post_residual_cov_control`, $T_{\text{post}} \times T_{\text{post}}$, after projecting out the within-arm prognostic predictors), the pre-period residual variances, and the trajectory draws above. The arm covariances are the shared-arm sampling noise $V_{m,a}$ of §2.2, so bipartite contrast graphs (vs-control stars, trees, even cycles) are identified without an odd cycle. A synthetic result without `arm_metadata` falls back to an isotropic allocation of the fits' own variances across arms, and `RankKDesignWarning` names the one extra fit that would close an odd cycle.
- **Treated-arm excess variance.** `RankKLongBetConfig.treated_extra_variance` defaults to `None`, which resolves to `True` under the default `baseline_cells="untreated"` (and to `False` under `"all"`; `True` with `"all"` raises `ValueError`): the JAX engine learns a treated-arm excess residual variance $s_\eta^2$ with a half-$t$ prior (draws in `diagnostics["s_eta2_draws"]`), so heteroskedastic treated arms carry calibrated Stage-1 uncertainty into the reconciliation.

Stage-1 defaults are unprojected (`tau_clip=(-inf, inf)`, `catt_floor=None`), so no truncation enters Stage 2 twice.

---

## 2. The joint hierarchical Bayesian model

### 2.1 Window and exposure basis

Let $H$ be the estimand window (`horizon`; default the full post-adoption window) and $\Psi = $ `orthonormal_exposure_basis(H, n_modes)` $\in \mathbb{R}^{H \times K'}$ the Euclidean-orthonormal basis whose first column is $\mathbf{1}_H / \sqrt{H}$ and whose second (default `n_modes=2`) is the orthonormalized concave ramp $1 - e^{-t/3}$. Each fitted pair $p = (a, b)$ of outcome $m$ yields the measurement

$$\hat{\theta}_{m,p} = \Psi^\top \bar{C}_{m,p} \in \mathbb{R}^{K'}, \qquad S_{m,p} = \Psi^\top \widehat{\operatorname{Cov}}(C_{m,p}) \Psi \in \mathbb{R}^{K' \times K'},$$

the posterior mean and covariance of the projected trajectory draws. The window Mode 1 (the $H$-period average effect) is $\hat{\theta}_{m,p,1} / \sqrt{H}$.

The parameters are arm-level coefficients $\Theta_m \in \mathbb{R}^{(A-1) \times K'}$ (control row fixed at 0). Every contrast is a linear map, $\theta_{m,p} = b_p^\top \Theta_m$ with $b_p \in \{-1, 0, 1\}^{A-1}$ the signed cut-space row of pair $p$, and every composite is a sum, $\Theta_c = \sum_{m \in c} \Theta_m$. Cycle identities and composite additivity therefore hold exactly for every draw, with no projection step.

### 2.2 Likelihood: SUR with shared-arm noise, robust forest error, composites as measurements

$$\hat{\theta}_{m,p} = b_p^\top \Theta_m + \xi_{m,a} - \xi_{m,b} + \eta_{m,p}, \qquad \hat{\theta}_{c,p} = b_p^\top \sum_{m \in c} \Theta_m + \sum_{m \in c} (\xi_{m,a} - \xi_{m,b}) + \eta_{c,p}.$$

- **Shared-arm sampling noise** $\xi_{m,a} \sim \mathcal{N}(0, V_{m,a})$, with $V_{m,a} = \Psi^\top \Sigma_{m,a}[1\!:\!H, 1\!:\!H]\, \Psi$ from the arm's Stage-1 residual covariance (averaged over the fits that contain arm $a$). Across outcomes, $\operatorname{Cov}(\xi_{m,a}, \xi_{\ell,a}) = V_{m,a}^{1/2} R_{m\ell} V_{\ell,a}^{1/2}$ with $R$ the within-unit residual correlation. The arm noise of a composite fit is *exactly* the sum of its components' arm noises (same units), so a composite fit adds information only through its own forest error; its marginal value is measured by the composite-fit ablation of the companion paper, and the composite is coherent with or without it.
- **Fit-specific (forest) error** $\eta_{m,p} \mid \omega_{m,p} \sim \mathcal{N}(0, F_{m,p} / \omega_{m,p})$. The forest covariance is identified from Stage 1 as $F_{m,p} = S_{m,p} - \rho_m (V_{m,a} + V_{m,b})$, floored at a share $\phi_m$ of $S_{m,p}$ in the simultaneous eigenbasis of $(S_{m,p}, V_{m,a} + V_{m,b})$, with hierarchical variance-decomposition parameters $\rho_m \sim \operatorname{Beta}(5, 1.5)$ and $\phi_m \sim \operatorname{Beta}(2, 10)$ (no cycle-based $\hat{\lambda}$). Robustness: $\omega_{m,p} \sim \operatorname{Gamma}(2, 2)$, that is a Student-$t_4$ forest error per fitted pair. $\mathbb{E}[\omega_{m,p}^{-1} \mid \text{data}]$ is reported as `contrast_scale_inflation` (it replaces the `flagged_pairs` of 0.10.0): on trees and bipartite graphs it equals the prior, since no cycle can contradict a fit, and on a triangle suspicion is shared across the three edges by their precisions.
- **Cross-outcome correlation** $R$. With `cross_outcome_residuals` (an `(N_units, M)` matrix of pre-period unit residuals, or the `CrossOutcomeCorrelationResult` of `estimate_cross_outcome_correlation`, which carries the covariate-adjusted within-unit residual scatter matrix and its degrees of freedom), $R$ is sampled every iteration from the conjugate inverse-Wishart posterior of the residual covariance, rescaled to a correlation matrix (`posterior_correlation_draws`). Without it, $R = I_M$. Zellner equivalence fails under cone priors and unequal shrinkage, so the joint posterior sharpens truncated outcomes using pure noise correlation, with no mediation or persistence assumption. When `composite_results` are supplied they already carry the cross-outcome noise cancellation of the composite, so the SUR block should not be credited with composite precision in that case; its measurable contribution is to the null-outcome CRPS and to composite intervals when no composite fit exists.

### 2.3 Priors: the poset cone as support, independent per-outcome scales

- **Positivity and monotonicity as prior support** (`monotone=True`). Let $\mathcal{E}_{\text{cover}}$ be the cover relations of the order (adjacent steps of a chain; the transitive reduction of `poset`). The support of $\Theta_m$ is the polyhedral cone $\mathcal{K} = \{\Theta_m : c^\top \mu_m \ge 0 \text{ and } c^\top \Theta_m \Psi_t^\top \ge 0 \ \forall c \in \mathcal{E}_{\text{cover}},\ t = 1..H\}$, where $\mu_m = \Theta_m \Psi^\top \mathbf{1}_H / H$ is the window Mode 1: every cover increment is non-negative on average *and* in every period of the window. This replaces the `tau_clip` / `catt_floor` draw projections of earlier versions and removes their amplitude crushing at the source (the constraint set is never empty, and the window mean is a parameter). Non-cover contrasts inherit monotonicity by transitivity. `monotone=False` (per outcome, as a mapping, in the multi-outcome function) drops the cone.
- **Cover-increment structure.** $\Theta_m \mid \sigma_m, \kappa_m \sim \mathcal{N}(0, \Sigma_0(\sigma_m, \kappa_m))\, \mathbf{1}_{\mathcal{K}}$ with Mode-1 precision $\sum_{c \in \mathcal{E}_{\text{cover}}} w_c\, c\, c^\top$ scaled so that the top arm has prior variance $\sigma_m^2$: independent half-normal increments on chains and trees, a prior symmetric under the $10 \leftrightarrow 01$ relabelling of a $2 \times 2$ diamond, and no radial Jacobian. The dynamic coefficients ($\Theta_{m, \cdot, 2:K'}$) have their own scale $\kappa_m$ and shrink toward a flat trajectory at low signal-to-noise.
- **Independent per-outcome generalized-horseshoe scale.** With $s_{0,m}$ the pre-treatment outcome SD (the median pre-period residual variance across the arms of outcome $m$; $2 \times$ the top contrast's Mode-1 SE on synthetic inputs without metadata), the standardized scale $\tau_m = \sigma_m / s_{0,m}$ has prior $\tau_m^2 \sim F(1/2, 1/2) = \operatorname{BetaPrime}(1/4, 1/4)$, that is a shrinkage factor $1 / (1 + \tau_m^2) \sim \operatorname{Beta}(1/4, 1/4)$: one smooth family with an integrable pole at $\tau_m = 0$ (null outcomes shrink to zero) and a $\tau_m^{-3/2}$ tail (active outcomes stay unshrunk), median $\tau_m = 1$. The dynamic-to-level ratio $\gamma_m = \kappa_m / \sigma_m$ has a half-$t_3(0.5)$ prior. Outcomes are independent in prior scale; cross-outcome borrowing arises only through $R$ and the composite measurements. Because $\mathcal{K}$ is a cone, the truncation normalizer $\Pr_{\mathcal{N}(0, \sigma_m^2 \Sigma_0)}(\mathcal{K})$ does not depend on $\sigma_m$, so the scale conditionals are exact one-dimensional updates (Huang–Wand representation for the half-$t$; slice updates on the log scale for $\tau_m$). `prior_scale=<float>` (one value, or a per-outcome mapping) fixes $\sigma_m$ to an externally calibrated value in outcome units; strings such as the former `"empirical_bayes"` / `"in_sample"` raise `ValueError`.

### 2.4 Inference: exact, tuning-free

Blocked Gibbs with the shared-arm noises $\xi$ marginalized analytically when updating $\Theta$ (so the chains mix even when the forest error is small relative to the arm noise):

1. $\Theta \mid \sigma, \kappa, \omega, \rho, \phi, R$: a multivariate Gaussian truncated to the polyhedral cone, sampled by **exact Hamiltonian Monte Carlo for truncated Gaussians** (Pakman & Paninski, 2014: exact harmonic trajectories with wall reflections, no step size) in posterior-whitened coordinates, started at the cone-constrained posterior mode (a small QP);
2. $\omega_{m,p}$: exact Gamma conditionals given the latent $\xi$ (drawn from its closed-form Gaussian conditional);
3. $\sigma_m, \kappa_m, \rho_m, \phi_m$: one-dimensional conditionals (inverse-gamma for the half-$t$ scale mixtures, slice updates otherwise) with an ancillarity–sufficiency interweaving step on the scales;
4. $R$: inverse-Wishart.

Four chains start overdispersed across the scale prior and run 200 warm-up and 400 retained sweeps each; blocks of the same length are appended (at most four times) until every parameter block meets split-$\hat{R} \le 1.01$ and bulk ESS $\ge 400$. `SamplerDiagnostics` reports `r_hat`, `ess`, `ess_tail`, `n_chains`, `n_draws_per_chain`, `n_warmup`, `n_continuation_rounds`, `converged`, and the per-block extrema `max_rhat_by_block` / `min_ess_by_block`. The sampler is deterministic given `seed`.

Simulation-based calibration over the full prior (every block: $\Theta$, $\sigma$, $\kappa$, $\rho$, $\phi$, $\omega$, $R$; 1,000 deterministic trials per configuration across a 3-arm chain, a $2 \times 2$ diamond, a 5-arm chain and a two-outcome design with $R$) passes the $\chi^2_{19}$ rank-uniformity test on every block with $\max \hat{R} \le 1.01$ and $\min \text{ESS} \ge 400$; the unit tests compare the truncated-Gaussian HMC against rejection sampling in $d \le 3$, the conjugate conditionals against numerical quadrature, and the diamond prior against its $10 \leftrightarrow 01$ symmetry.

### 2.5 Model checking shipped with the result

- `ppc_cycle_pvalue`: posterior predictive $p$-value of the cycle discrepancy $N^\top \hat{y}$ of the redundant (non-tree) fits; the cycle-space identity is a *check*, not an estimator.
- `ppc_composite_pvalue`: composite-fit-versus-sum discrepancy when `composite_results` are supplied.
- `ppc_basis_pvalue`: basis misfit, the residual trajectory energy outside $\Psi$ (raise `n_modes` if it rejects).
- `contrast_scale_inflation` ($\mathbb{E}[\omega_{m,p}^{-1}]$ per fitted pair) and `posterior_corruption_prob`.
- `RankKDesignWarning` for a bipartite contrast graph without arm metadata, naming the one fit that would add an odd cycle (information, not a policy).

---

## 3. Outputs

`MultiArmReconciliationResult` (single outcome) holds the design (`arms`, `order`, `cover_edges`, `topology`, `tree_edges`, `fitted_pairs`, `redundant_pairs`, `all_pairs`, `horizon`, `n_modes`, `basis`), the parameter draws (`theta_draws` of shape `(S, A - 1, K')`, `sigma_draws`, `kappa_draws`, `tau_draws`, `rho_draws`, `phi_draws`, `omega_draws`, `local_scale_draws`), derived draws for every pair of arms (`arm_trajectory_draws`, `tau_mu_draws`, `tau_increment_draws`, `contrast_mode1_draws`, `contrast_trajectory_draws` of shape `(S, H)`, `contrast_cumulative_draws`, `contrast_catt_draws`), the checks of §2.5, `max_cycle_discrepancy`, `s0_scale`, `prior_scale` and `sampler_diagnostics`. `get_contrast_result(treated, ref)` returns a `ContrastReconciliationResult` (`mode1_draws`, `trajectory_draws`, `cumulative_draws`, `mean`, `sd`, `ci_95`) for any ordered pair, fitted or implied; `trt_catt_bc` / `trt_catt_ab` / `trt_catt_ac` are the nested three-arm trajectories.

`MultiArmMultiOutcomeReconciliationResult` holds one `MultiArmReconciliationResult` per structural outcome and per composite in `by_outcome`, the joint draws keyed by outcome, `posterior_correlation_draws` (and `posterior_correlation_mean`), `additivity_violation_rate` and `max_abs_composite_discrepancy` per composite (always `0.0` and $< 10^{-12}$), the pooled checks of §2.5 and `sampler_diagnostics`. `Nested3ArmReconciliationResult` (from `reconcile_nested_3arm_draws`) unpacks as `(trt_catt_bc, trt_catt_ab, trt_catt_ac)`, carries the three `post_att_*` Mode-1 draws and the full `multi_arm_result`.

`multi_arm_decision_table`, `roi_break_even_table`, `extrapolate_rank_k_decay`, `compute_required_persistence` and `ge_revenue_draws` accept these results directly.

---

## 4. Python example

```python
import longbet

# 1. Single outcome, 4-arm graded chain C < B1 < B2 < A with 6 pairwise Rank-K BCF fits
#    (keys are (arm_high, arm_low); any connected set of pairs works, redundant pairs
#    feed the cycle posterior predictive check)
arms = ("C", "B1", "B2", "A")
rec = longbet.reconcile_multi_arm_draws(
    {
        ("B1", "C"): fit_b1_c,
        ("B2", "B1"): fit_b2_b1,
        ("A", "B2"): fit_a_b2,
        ("B2", "C"): fit_b2_c,
        ("A", "B1"): fit_a_b1,
        ("A", "C"): fit_a_c,
    },
    arms=arms,
    order=arms,
    horizon=12,                                        # 12-period estimand window
)
print(rec.sampler_diagnostics.converged, rec.sampler_diagnostics.r_hat, rec.sampler_diagnostics.ess)
print(rec.max_cycle_discrepancy, rec.ppc_cycle_pvalue, rec.ppc_basis_pvalue)
print(rec.contrast_scale_inflation)                    # E[1 / omega_p] per fitted pair
ac = rec.get_contrast_result("A", "C")                 # any (treated, ref) pair, fitted or implied
print(ac.mean, ac.ci_95)
traj_b2_b1 = rec.contrast_trajectory_draws[("B2", "B1")]   # (S, 12) coherent cover-increment draws

# 2. A 2 x 2 factorial as a poset of cover relations (arm_high, arm_low); "00" is control
rec_2x2 = longbet.reconcile_multi_arm_draws(
    fits_2x2, arms=("00", "10", "01", "11"),
    poset=(("11", "10"), ("11", "01"), ("10", "00"), ("01", "00")),
)

# 3. Decision and ROI break-even tables
dec_df = longbet.multi_arm_decision_table(
    rec,
    arm_costs={"C": 0.0, "B1": 1.0, "B2": 2.0, "A": 3.5},
    incremental_costs={("B1", "C"): 1.0, ("B2", "B1"): 1.0, ("A", "B2"): 1.5},
)
roi_df = longbet.roi_break_even_table(rec, cost_grid=(1.0, 2.0, 3.5), alpha_grid=(0.5, 1.0), contrast=("A", "C"))

# 4. Two outcomes, their composite (a fitted composite is an extra measurement, not a
#    requirement), and the cross-outcome correlation from pre-period residuals
corr = longbet.estimate_cross_outcome_correlation(
    y_post={"core_rev": y_core, "addon_rev": y_addon}, y_pre={"core_rev": y_core_pre, "addon_rev": y_addon_pre},
    x_covs=x, arm_labels=arm_of_unit,
)
multi_rec = longbet.reconcile_multi_arm_multi_outcome_draws(
    {"core_rev": fits_core, "addon_rev": fits_addon},
    composite_outcomes={"total_rev": ("core_rev", "addon_rev")},
    composite_results={"total_rev": fits_total},       # optional
    cross_outcome_residuals=corr,                      # or an (N_units, M) residual matrix; None -> R = I
    arms=arms,
    order=arms,
    horizon=12,
)
print(multi_rec.additivity_violation_rate)             # {'total_rev': 0.0}
print(multi_rec.max_abs_composite_discrepancy)         # {'total_rev': < 1e-12}
print(multi_rec.posterior_correlation_mean)            # posterior mean of R
total = multi_rec.by_outcome["total_rev"].get_contrast_result("A", "C")

# 5. Nested three-arm design (C = control, B = base, A = base + incremental): a thin wrapper
bc, ab, ac = longbet.reconcile_nested_3arm_draws(res_bc, res_ab, res_ac)   # ac == bc + ab on every draw
```

---

## 5. Changes from 0.10.0 and migration

Everything removed below was a knob, a plug-in or a projection that the joint model makes unnecessary; passing any of them raises `TypeError` (unknown keyword) or `ValueError` (`prior_scale` strings).

| Removed in 0.11.0 | Replaced by |
|---|---|
| `prior_scale="empirical_bayes"` / `"in_sample"`, `prior_scale_floor` | the per-outcome generalized-horseshoe scale prior (§2.3); `prior_scale=<float>` still fixes the scale |
| `prior_mode` (`bundle_*`, `independent_*`, `flat`, `bundle_spike_halfnormal`), `spike_prob`, `split_alpha`, `prior_mean` | the poset-cone prior with cover-increment structure (§2.3); a null outcome is handled by the scale prior's pole at zero, not by a point mass |
| `projection="huber"`, `projection_weights`, `flagged_pairs` | the Student-$t_4$ forest error with latent precisions; `contrast_scale_inflation` |
| `sampler`, `n_proposals`, `n_draws`, `n_burnin`, `n_chains`, `auto_continue`, `target_rhat`, `target_ess`, `max_iter_multiplier` | one exact sampler with unconditional convergence control (§2.4) |
| `tau_clip`, `catt_floor` and the damped-trajectory projection inside the reconciliation | positivity as prior support (§2.3); the Stage-1 `RankKLongBetConfig` options keep their single-fit meaning |
| `composite_mode` (`"bottom_up"` / `"coherent_total"`), `key_contrast`, `composite_floor_diagnostics` | composite = sum of component parameters; a composite fit is a measurement (§2.2) |
| `iman_conover_couple`, `ImanConoverResult`, `cross_outcome_corr=` | the joint SUR posterior with $R$ sampled from its inverse-Wishart conditional (`cross_outcome_residuals=`) |
| `arm_variance_source`, `arm_outcome_variances`, `bipartite_policy`, `bipartite_ridge*`, `min_group_var_share`, `null_shift_interval`, `sample_sizes`, `num_burnin`, `return_prior_draws`, `topology=` | arm noise covariances from Stage-1 `arm_metadata` (§1); the topology follows from `order` / `poset` |
| the separate bodies of `reconcile_nested_3arm_draws` / `reconcile_nested_3arm_multi_draws` (`share_mean_logit`, `share_sd_logit`, `prior_sd_*`, `prior_mean_*`) and the `A = 3` delegation path | thin wrappers over the A-arm model |

Migration from a 0.10.0 call: drop every keyword except `arms`, `horizon`, `n_modes`, `seed` and (optionally) a float `prior_scale`; replace `topology="chain"` by `order=arms` and `topology="poset"` by the cover relations in `poset=`; replace `cross_outcome_corr=estimate_cross_outcome_correlation(...)` by `cross_outcome_residuals=estimate_cross_outcome_correlation(...)`; pass the composite's own fits, if any, as `composite_results`. Re-run the Stage-1 fits with 0.11.0 when they were produced by an earlier version, so that `arm_metadata` is present (fits without it fall back to the isotropic allocation and, on bipartite graphs, warn).
