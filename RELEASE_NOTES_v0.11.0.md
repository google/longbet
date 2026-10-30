# Release Notes — `longbet-jax` v0.11.0

`longbet-jax` `v0.11.0` introduces a unified Stage-2 joint hierarchical Bayesian reconciliation engine for multi-arm and multi-outcome longitudinal experiments (`reconcile_multi_arm_draws` and `reconcile_multi_arm_multi_outcome_draws`), stage-1 arm residual covariance propagation, and a streamlined zero-knob public reconciliation API.

## Highlights

### 1. Stage-1 Arm Residual Metadata & Default Heteroskedastic Treated Variance
- **`RankKLongBetResult.arm_metadata` & `RankKLongBetResult.rank_k`**: Stage-1 fits now attach lightweight arm-level sufficient statistics (`n_treated`, `n_control`, `post_residual_cov_treated`, `post_residual_cov_control`, `pre_residual_var_treated`, `pre_residual_var_control`, and event-time `trt_trajectory_draws`) together with the fitted exposure-basis rank `rank_k: int`. Stage-2 reconciliation projects $V_{m,a}$ onto any post-adoption window $1..H$ and exposure basis $\Psi$ without requiring raw unit-level panels, and automatically builds partial-mode observation selectors whenever pairwise contrasts have heterogeneous basis ranks $K_p \le K'$.
- **Default `treated_extra_variance` on `"untreated"` baselines**: `RankKLongBetConfig.treated_extra_variance` now defaults to `None` (resolving to `True` when `baseline_cells == "untreated"` and `False` otherwise), propagating learned excess treated residual variance $s_\eta^2$ into `catt_draws` for calibrated Stage-1 uncertainty under heteroskedastic treatment effects.

### 2. Joint Hierarchical Bayesian Multi-Arm & Multi-Outcome Reconciliation
- **Exact Structural Coherence by Construction**: Parameters $\Theta_m \in \mathbb{R}^{(A-1)\times K'}$ live directly on non-control arm basis coordinates over the Euclidean-orthonormal exposure basis $\Psi \in \mathbb{R}^{H\times K'}$ ($\Psi_{\cdot,1} = \mathbf{1}_H / \sqrt{H}$). Every pairwise contrast $\theta_{m,p} = b_p^\top \Theta_m$ and accounting composite $\Theta_c = \sum_{m \in c} \Theta_m$ is an exact linear map at every posterior draw (`AVR = 0`, cycle/composite discrepancy `< 1e-13`).
- **Collapsed Seemingly Unrelated Regressions (SUR) + Robust Student-$t_4$ Forest Likelihood**: Shared arm sampling shocks $\xi_{\cdot, a}$ (correlated across outcomes via $R$) are analytically marginalized in the collapsed Gaussian likelihood $\Omega(\rho, \phi, \omega, R)$ across structural and composite Stage-1 summaries, paired with Matheron conditional updates for latent contrast precision weights $\omega_q \sim \mathrm{Gamma}(2, 2)$ and conjugate pre-period Inverse-Wishart correlation updates for $R$.
- **Poset-Cone Support & Generalized-Horseshoe Per-Outcome Scale Prior**: Under `monotone=True` (default), draws are supported on the polyhedral cone enforcing cover-edge window-mean monotonicity and per-period trajectory non-negativity, sampled via exact Pakman–Paninski (2014) wall-reflecting HMC in posterior-whitened coordinates. Independent per-outcome standardized scales $\tau_m = \sigma_m / s_{0,m}$ follow $(\tau_m)^2 \sim \mathrm{BetaPrime}(1/4, 1/4)$ anchored strictly to pre-period outcome SD $s_{0,m}$ (`prior_scale=None`), with dynamic decay scale $\kappa_m / \sigma_m \sim \mathrm{Half}\text{-}t_3(0.50)$ and Yu–Meng (2011) ASIS interweaving across 4 auto-continued chains (`split-Rhat <= 1.01`, `ESS >= 400`).

### 3. Posterior Predictive Diagnostics & Pre-Treatment Correlation Estimation
- **`MultiArmReconciliationResult` & `MultiArmMultiOutcomeReconciliationResult` diagnostics**:
  - `ppc_cycle_pvalue`: Cycle-space discrepancy posterior predictive tail probability (`None` on trees/bipartite graphs without cycles).
  - `ppc_composite_pvalue`: Composite-vs-component-sum measurement discrepancy p-value on composite outcomes (`None` on structural outcomes).
  - `ppc_basis_pvalue`: Out-of-basis trajectory residual energy p-value (`None` when Stage-1 trajectory draws already lie strictly inside $\operatorname{col}(\Psi)$).
  - `contrast_scale_inflation`: Posterior mean $E[1/\omega_p \mid \text{data}]$ per contrast pair.
- **`estimate_cross_outcome_correlation`**: Accepts `y` as its primary first parameter (with `y_post` retained as a backward-compatible keyword alias) and handles staggered `NaN` rows without empty-slice warnings.

### 4. API Streamlining (`v0.11.0`)
- Removed deprecated/ad-hoc plug-in reconciliation options (`prior_scale="empirical_bayes"` / `"in_sample"`, radial `prior_mode` variants, Huber cut-space projection, SIR/slice post-hoc samplers, post-sampling `tau_clip`/`catt_floor` reconciliation clamping, Iman–Conover copula coupling, synthetic composite rescaling, and manual bipartite ridge heuristics).
- `reconcile_nested_3arm_draws` and `reconcile_nested_3arm_multi_draws` are retained as thin convenience wrappers delegating directly to `reconcile_multi_arm_draws` and `reconcile_multi_arm_multi_outcome_draws`.
