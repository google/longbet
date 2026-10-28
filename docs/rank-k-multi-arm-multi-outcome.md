# Coherent Multi-Arm and Multi-Outcome Reconciliation for Rank-K BCF

When a longitudinal experiment compares $A \ge 3$ arms (such as a graded dose ladder $C \prec B_1 \prec B_2 \prec \dots \prec A_{\text{top}}$ or a $2 \times 2$ factorial partial order) across $M \ge 1$ outcomes and linear composites, fitting each pairwise contrast $(j, \ell)$ and each outcome $m$ independently with `fit_rank_k_longbet` introduces four coherence gaps:

1. **Shared reference-arm covariance.** Contrasts that share a reference or intermediate arm $b$ have non-zero sampling covariance $\operatorname{Cov}(\hat{\tau}_{(a,b)}, \hat{\tau}_{(b,c)}) = -v_b < 0$ and $\operatorname{Cov}(\hat{\tau}_{(j,0)}, \hat{\tau}_{(\ell,0)}) = +v_0 > 0$. Treating pairwise posteriors as independent misstates the variance of increments and sums.
2. **Cycle additivity violations.** Over any cycle in the contrast graph, independent pairwise draws violate $\tau_{(a,c)} = \tau_{(a,b)} + \tau_{(b,c)}$ on every posterior draw and period.
3. **Prior order-reversal and compounded truncation.** Independent priors across nested arms allow non-monotone posterior quantiles ($Q_q(\tau_{\text{Combined}}) < Q_q(\tau_{\text{Base}})$), while summing $K$ independently truncated component outcomes inflates the null composite mean proportionally to $K$.
4. **Cross-outcome correlation.** Multiple outcomes measured on the same panel units covary across draws, so multi-outcome monetization ($\sum_m \alpha_m \text{Cum}_m^{(s)}$) and composite outcomes require rank-copula coupling that preserves each outcome's marginal distribution and within-outcome cycle identities.

Version 0.9.0 generalizes the 3-arm nested reconciliation of 0.8.0 (`reconcile_nested_3arm_draws` and `reconcile_nested_3arm_multi_draws`) to arbitrary $A \ge 2$ arms and $P \ge A - 1$ fitted contrasts via `reconcile_multi_arm_draws` (`longbet._rank_k_multi_arm`) and `reconcile_multi_arm_multi_outcome_draws` (`longbet._rank_k_multi`), while delegating the default $A = 3$ chain case to `reconcile_nested_3arm_draws` for bit-for-bit backward compatibility.

---

## 1. Cut-Space Likelihood and Arm-Variance Identification

Let $\mathcal{V} = \{0, 1, \dots, A - 1\}$ denote the $A$ arms with $0$ as the control arm, and let $\mathcal{E}_{\text{fit}} = \{(j_p, \ell_p)\}_{p=1}^P$ be the $P \ge A - 1$ fitted contrasts forming a connected graph $\mathcal{G} = (\mathcal{V}, \mathcal{E}_{\text{fit}})$. Every coherent vector of pairwise Mode-1 effects lies in the $(A - 1)$-dimensional cut space parameterized by the vs-control effects $\boldsymbol{\mu} = (\mu_1, \dots, \mu_{A-1})^\top \in \mathbb{R}^{A-1}$ (with $\mu_0 \equiv 0$):

$$\tau_{(j,\ell)} = \mu_j - \mu_\ell = \mathbf{b}_p^\top \boldsymbol{\mu}, \qquad \mathbf{B} \in \{-1, 0, 1\}^{P \times (A-1)}.$$

Because $\mathcal{G}$ is connected, $\operatorname{rank}(\mathbf{B}) = A - 1$. Furthermore, under the arm-additive sampling covariance $\mathbf{V}_{\text{arm}} = \operatorname{diag}(v_0, v_1, \dots, v_{A-1})$, the covariance of the $P$ pairwise Mode-1 summaries is $\mathbf{E}\mathbf{V}_{\text{arm}}\mathbf{E}^\top$, whose column space is invariant to $\mathbf{V}_{\text{arm}}$ on the cut space. Consequently, the generalized least squares (GLS) sufficient statistics on $\boldsymbol{\mu}$ are:

$$\hat{\boldsymbol{\mu}} = (\mathbf{B}^\top \mathbf{B})^{-1} \mathbf{B}^\top \hat{\mathbf{y}}, \qquad \boldsymbol{\Sigma}_{\text{lik}} = \operatorname{diag}(\hat{v}_1, \dots, \hat{v}_{A-1}) + \hat{v}_0 \mathbf{1}_{A-1}\mathbf{1}_{A-1}^\top,$$

with cycle discrepancy diagnostic $r_{\text{cycle}} = \|\hat{\mathbf{y}} - \mathbf{B}\hat{\boldsymbol{\mu}}\|_2 / \max(\|\hat{\mathbf{y}}\|_2, 10^{-12})$.

Since 0.10.0, `projection="huber"` (opt-in) replaces the Euclidean projection by iteratively reweighted least squares $\hat{\boldsymbol{\mu}} = (\mathbf{B}^\top\mathbf{W}\mathbf{B})^{-1}\mathbf{B}^\top\mathbf{W}\hat{\mathbf{y}}$ with Huber weights $w_p = \min(1, 1.345 / |r_p / \hat{s}|)$ on the cycle residual $\mathbf{r} = \hat{\mathbf{y}} - \mathbf{B}\hat{\boldsymbol{\mu}}$, $\hat{s} = 1.4826\,\operatorname{MAD}(\mathbf{r})$. Under the exact shared-arm model $\hat{\mathbf{y}} \in \operatorname{col}(\mathbf{B})$, all weights are one and Proposition 1 is unchanged; a single corrupted fit is down-weighted and reported in `flagged_pairs` (with `projection_weights`). $\boldsymbol{\Sigma}_{\text{lik}}$ is unchanged, since any weighted projection of a column-space vector has the same sampling covariance.

### Arm-Variance Recovery: Non-Bipartite vs. Bipartite Contrast Graphs

Each pairwise fit $p = (j_p, \ell_p)$ provides a Mode-1 posterior variance $s_p^2 \approx v_{j_p} + v_{\ell_p}$, defining the unsigned incidence system $\mathbf{A}_{\text{var}} \mathbf{v} \approx \mathbf{s}^2$ with $\mathbf{A}_{\text{var}} \in \{0, 1\}^{P \times A}$:

* **Non-bipartite contrast graphs (`is_bipartite == False`, e.g., any graph containing a 3-cycle or odd cycle):** $\operatorname{rank}(\mathbf{A}_{\text{var}}) = A$, so all $A$ arm variances are point-identified. `reconcile_multi_arm_draws` solves the bounded non-negative least squares (NNLS) problem
  $$\hat{\mathbf{v}} = \arg\min_{\mathbf{v} \ge \mathbf{v}_{\min}} \|\mathbf{A}_{\text{var}}\mathbf{v} - \mathbf{s}^2\|_2^2, \qquad v_{\min, a} = \eta \min_{p \ni a} s_p^2 \quad (\eta = 0.15),$$
  which reduces in closed form to `_backsolve_group_variances` for the $A = 3$ triangle.
* **Bipartite contrast graphs (`is_bipartite == True`, e.g., trees, star graphs of vs-control fits, or even cycles):** $\operatorname{rank}(\mathbf{A}_{\text{var}}) = A - 1$ with 1-dimensional null vector $z_a = (-1)^{\text{color}(a)}$, so shifting $v_a \mapsto v_a + \xi (-1)^{\text{color}(a)}$ leaves every pairwise sum $v_{j_p} + v_{\ell_p}$ unchanged while altering off-diagonal covariances (in particular the shared-control covariance $v_0$). `reconcile_multi_arm_draws` detects bipartiteness via 2-coloring and emits a `RankKDesignWarning` (a `UserWarning` subclass, new in 0.10.0) that names the missing odd cycle and the one extra pairwise fit that would close it (for a vs-control star on $C \prec B_1 \prec B_2 \prec B_3$: fit $(B_2, B_1)$ to close $B_1 \to C \to B_2 \to B_1$). The unidentified shift $\xi$ is resolved by `bipartite_policy` over the feasible interval $[\xi_{\text{lo}}, \xi_{\text{hi}}]$ that keeps every $v_a$ above its floor (reported as `null_shift_interval` / `null_shift`):
  * `"allocation"` (default, 0.9.0 behavior): sample-size-informed ridge-regularized NNLS toward $v_{\text{target}, a} \propto 1 / n_a$ ($\xi = 0$).
  * `"conservative"`: the $\xi$ that maximizes the Gaussian-approximate posterior variance of the top contrast (prior second moments from prior draws).
  * `"marginalize"`: integrates $\xi$ under a flat prior on the feasible interval, i.e. replaces the Gaussian likelihood by the 9-point scale mixture $\frac{1}{9}\sum_k \mathcal{N}(\hat{\boldsymbol{\mu}} \mid \boldsymbol{\mu}, \boldsymbol{\Sigma}(\xi_k))$.

* **Direct arm variances (`arm_variance_source="direct"`, new in 0.10.0, opt-in):** instead of back-solving $\mathbf{A}_{\text{var}}\mathbf{v} \approx \mathbf{s}^2$, set $v_a = c\,\hat{\sigma}_a^2 / n_a$ from the arms' own residual outcome variances (`arm_outcome_variances`, e.g. within-arm pre-period residual variances) and unit counts, with one scale $c = \sum_p s_p^2 w_p / \sum_p w_p^2$, $w_p = \hat{\sigma}_{j_p}^2/n_{j_p} + \hat{\sigma}_{\ell_p}^2/n_{\ell_p}$. The allocation across arms then comes from data, so bipartite designs are identified and odd cycles serve only as a diagnostic.

  Since 0.10.0 a chain order fitted without all of its adjacent edges (e.g. only the vs-control star) is accepted with `topology="chain"`: it is reconciled as the poset whose cover relations are the chain edges plus the fitted pairs, so the spanning tree comes from the fitted pairs and the chain order remains the monotonicity constraint (with the default uniform split this prior equals the chain's $\operatorname{Dirichlet}(\mathbf{1})$ bundle prior).

---

## 2. Topologies, Priors, and Dimension-Adaptive Sampling

Given a spanning tree $\mathcal{T}$ of $\mathcal{G}$ with $A - 1$ directed edges $e_a = (\text{child}_a, \text{parent}_a)$ and incidence matrix $\mathbf{T}_{\text{tree}} \in \{-1, 0, 1\}^{(A-1) \times (A-1)}$, the edge increments are $\boldsymbol{\theta} = \mathbf{T}_{\text{tree}} \boldsymbol{\mu}$ and $\boldsymbol{\mu} = \mathbf{M}_{\text{tree}} \boldsymbol{\theta}$ where $\mathbf{M}_{\text{tree}} = \mathbf{T}_{\text{tree}}^{-1} \in \{0, 1\}^{(A-1) \times (A-1)}$ is the path-indicator matrix.

### Prior Families (`prior_mode`)

1. **Nested chain bundle prior (`topology="chain"`, `prior_mode` in `{"bundle_halfnormal", "bundle_halfstudentt", "bundle_gamma_spike"}`):**
   Places the bundle prior on the top-vs-control effect $\tau_{\text{top}} = \mu_{A-1} \sim p_{\text{bundle}}(\cdot; \sigma)$ and splits it across the $A - 1$ adjacent increments via $\mathbf{u} \sim \operatorname{Dirichlet}(\boldsymbol{\alpha})$ (`split_alpha`, default $\mathbf{1}_{A-1}$), setting $\theta_a = u_a \tau_{\text{top}} \ge 0$ and $\mu_j = \sum_{a=1}^j \theta_a$.
2. **Unique-top poset order-polytope prior (`topology="poset"`, `prior_mode` in `{"bundle_halfnormal", ...}`):**
   For a partial order with a unique maximal arm `top` (such as a $2 \times 2$ factorial design $C \prec \{B_1, B_2\} \prec A_{12}$), draws $\mu_{\text{top}} = \tau_{\text{top}} \sim p_{\text{bundle}}(\cdot; \sigma)$ and sets $\mu_j = \tau_{\text{top}} w_j$ for interior arms $j \notin \{0, \text{top}}$, where $\mathbf{w}$ is drawn uniformly from the order polytope $\mathcal{O}(\mathcal{P}_{\text{int}}) = \{\mathbf{w} \in [0, 1]^{A-2} : w_\ell \le w_j \text{ whenever } \ell \preceq j\}$ via Stanley's triangulation into uniform linear extensions.
3. **Independent truncated prior (`prior_mode` in `{"independent", "independent_halfnormal", "independent_normal"}`):**
   Places independent Normal or Half-Normal priors on the coordinates of $\boldsymbol{\mu}$ (or on $\boldsymbol{\theta}$ for `"independent_halfnormal"` on chains), truncated to `tau_clip` and the poset monotonicity cone. Required when a poset has multiple maximal arms and no single top arm.
4. **Flat prior (`prior_mode="flat"`):**
   Uniform prior on `tau_clip` (and the monotonicity cone if `monotone=True`). When `monotone=False` and `tau_clip=(-np.inf, np.inf)`, the posterior is the exact closed-form Gaussian $\mathcal{N}(\hat{\boldsymbol{\mu}}, \boldsymbol{\Sigma}_{\text{lik}})$.
5. **Spike-and-slab bundle prior (`prior_mode="bundle_spike_halfnormal"`, new in 0.10.0, opt-in):**
   $\tau_{\text{top}} \sim \pi\,\delta_0 + (1 - \pi)\,\operatorname{HalfNormal}(\sigma)$ with $\pi$ = `spike_prob` (default $0.5$, which is what $\pi \sim \operatorname{Beta}(1, 1)$ integrates to for a single reconciliation); on the slab the increments keep the bundle prior's split, on the spike $\boldsymbol{\mu} = \mathbf{0}$. The posterior spike weight is exact given the slab evidence, $\Pr(\text{spike} \mid \hat{\boldsymbol{\mu}}) = \pi\,\mathcal{N}(\hat{\boldsymbol{\mu}}; \mathbf{0}, \boldsymbol{\Sigma}_{\text{lik}}) / [\pi\,\mathcal{N}(\hat{\boldsymbol{\mu}}; \mathbf{0}, \boldsymbol{\Sigma}_{\text{lik}}) + (1 - \pi) Z_{\text{slab}}]$, where the spike term is closed form and $Z_{\text{slab}}$ is the deterministic-mixture importance-sampling evidence of the defensive SIR proposals; draws on the spike are exactly zero on every contrast and period (`spike_posterior_prob` on the result). Under a `coherent_total` composite, a spiked composite draw sets all its components to zero on that draw.

### Prior Scale (`prior_scale`)

`"empirical_bayes"` (default since 0.10.0) is the Type-II maximum-likelihood scale under the Gaussian cut-space likelihood of the top contrast, $\hat{\sigma} = \arg\max_{\sigma \ge \sigma_{\text{floor}}} \log \int p_{\text{bundle}}(\tau; \sigma)\, \mathcal{N}(\hat{\tau}_{\text{top}}; \tau, s_{\text{top}}^2)\, d\tau$; for the half-normal family the marginal is the closed-form skew-normal $2\,\phi(\hat{\tau}; 0, \sigma^2 + s^2)\,\Phi(\hat{\tau}\sigma / (s\sqrt{\sigma^2 + s^2}))$, otherwise a 1-D quadrature (the spike is mixed in exactly for the spike-and-slab prior; the independent family uses each coordinate's truncated-normal marginal). The floor $\sigma_{\text{floor}}$ (`prior_scale_floor`) defaults to the in-sample noise scale. `"in_sample"` (0.9.0 behavior, opt-out) sets the bundle scale $\sigma$ to the in-sample posterior SD of the top contrast (its noise scale), which shrinks a top effect several posterior SDs from zero substantially toward zero.

### Dimension-Adaptive Posterior Sampler (`sampler`)

* **`"sir"` (default under `"auto"` for $d = A - 1 \le 3$):** 50/50 defensive-mixture Sampling Importance Resampling mixing prior proposals on the cone with Gaussian likelihood proposals $\mathcal{N}(\hat{\boldsymbol{\mu}}, \boldsymbol{\Sigma}_{\text{lik}})$. If `sampler="auto"` and `ess < max(100, 0.02 * n_proposals)`, automatically falls back to `"slice"`.
* **`"slice"` (default under `"auto"` for $d = A - 1 \ge 4$):** Vectorized multi-chain analytical-polyhedron hit-and-run slice sampler with radial-cone parameterization $\mathbf{z} = (\mathbf{w}, \tau_{\text{top}})$ that eliminates the Cartesian $\tau_{\text{top}}^{-(d-1)}$ Jacobian singularity at the origin, reporting split-$\hat{R}$ and bulk ESS across chains in `sampler_diagnostics`. Since 0.10.0 (`auto_continue=True`) the chains run in blocks: after the default `num_burnin + ceil(S / num_chains) * thin` sweeps the sampler checks split-$\hat{R} \le 1.01$ (`target_rhat`) and bulk ESS $\ge \min(400, 0.75 S)$ (`target_ess`) on every coordinate and, if either fails, appends blocks of the same length (discarding the first half of the blocks as extra warm-up and thinning back to $S$ draws) up to `max_iter_multiplier = 8` times the default, then emits a `RankKConfigWarning`. A run that already meets the targets is bit-identical to 0.9.0; `auto_continue=False` always is. Sweeps used are reported in `sampler_diagnostics.n_iterations` (with `converged`).

### Spanning-Tree Trajectory Projection

Each spanning-tree edge increment $\theta_e^{(s)}$ is projected onto its post-adoption trajectory $\mathbf{c}_e^{(s)} \in \mathbb{R}^{T_{\text{post}}}$ using that edge's fitted `RankKLongBetResult` via `_project_cohort_damped_trajectory` (preserving higher-mode exposure dynamics and `catt_floor` damping). Every fitted or implied contrast $(j, \ell)$ is then assembled as the exact signed path sum along $\mathcal{T}$:

$$\mathbf{c}_{(j,\ell)}^{(s)} = \sum_{e \in \text{path}_{\mathcal{T}}(\ell \to j)} \operatorname{sign}(e)\, \mathbf{c}_e^{(s)},$$

guaranteeing that cycle identities hold to machine precision ($\max |\text{cycle gap}| < 10^{-12}$) across all contrasts and periods. For redundant (non-tree) fitted contrasts, `own_vs_implied_max_gap` (and `trajectory_additivity_gaps`) records the maximum absolute discrepancy between the contrast's own projected trajectory and its spanning-tree-implied trajectory.

---

## 3. Multi-Outcome Coupling and Coherent Composites

`reconcile_multi_arm_multi_outcome_draws` reconciles $M$ structural outcomes across $A$ arms in three steps:

1. **Per-outcome multi-arm reconciliation:** Runs `reconcile_multi_arm_draws` for each structural outcome $m \in \{1, \dots, M\}$.
2. **Marginal-preserving Iman-Conover coupling on `key_contrast`:** Permutes the draw indices $s \in \{0, \dots, S - 1\}$ of each outcome $m$ via `iman_conover_couple` on `key_contrast` (default `(arms[-1], arms[0])`) to match the target cross-outcome correlation matrix (`estimate_cross_outcome_correlation`). Because the same permutation $\pi_m$ is applied to all contrasts and periods of outcome $m$, every marginal distribution and within-outcome cycle identity is preserved bit for bit.
3. **Composite outcome synthesis (`composite_outcomes`):** For each requested composite $C = \sum_{m \in \mathcal{M}_C} m$:
   * `composite_mode="coherent_total"` (default since 0.10.0 when the composites have disjoint components and the prior truncates at zero): samples the composite's Mode-1 vector $\boldsymbol{\mu}_C^{(s)}$ under a single bundle prior truncation (preventing $K$-fold truncation inflation under the null), rank-couples $\boldsymbol{\mu}_C$ to the component sum on `key_contrast`, and proportionally rescales the component draws by $r^{(s)} = \tau_{C, \text{key}}^{(s)} / \sum_{m \in \mathcal{M}_C} \tau_{m, \text{key}}^{(s)}$ so that exact draw-wise additivity (`additivity_violation_rate == 0.0`, `max_abs_composite_discrepancy < 1e-12`) holds simultaneously across all $P$ contrasts and all post-adoption periods. Without `composite_results`, the composite likelihood is formed from the components' GLS moments and the cross-outcome correlation (the synthetic-composite route), so no extra composite fit is needed.
   * `composite_mode="bottom_up"`: sums the coupled component draws directly across all contrasts and periods (the 0.9.0 default; bit-identical to 0.9.0 when passed explicitly). Near the null this inflates the composite floor by about $\sum_m \sigma_m / \sqrt{\boldsymbol{\sigma}^\top \mathbf{R} \boldsymbol{\sigma}}$; a `RankKConfigWarning` reports that factor whenever the unconstrained composite posterior on `key_contrast` has $\Pr(\tau \le 0) > 0.10$. The numbers are stored in `composite_floor_diagnostics`.

---

## 4. Python Example

```python
import longbet

# 1. Single-outcome 4-arm graded chain (C < B1 < B2 < A) with 6 pairwise Rank-K BCF fits
arms = ["C", "B1", "B2", "A"]
rec = longbet.reconcile_multi_arm_draws(
    results={
        ("B1", "C"): fit_b1_c,
        ("B2", "B1"): fit_b2_b1,
        ("A", "B2"): fit_a_b2,
        ("B2", "C"): fit_b2_c,
        ("A", "B1"): fit_a_b1,
        ("A", "C"): fit_a_c,
    },
    arms=arms,
    topology="chain",
    prior_mode="bundle_halfnormal",
    prior_scale=0.35,
)
print(rec.arm_variances.to_dict())                # recovered per-arm sampling variances
print(rec.cycle_discrepancy_ratio)                # cut-space projection residual ratio
print(rec.sampler_diagnostics)                    # SamplerDiagnostics(method='sir', ess=..., r_hat=...)
ac = rec.get_contrast_result("A", "C")            # ContrastReconciliationResult for any (treated, ref)
print(ac.mean, ac.ci_low, ac.ci_high)

# 2. Multi-arm decision and ROI break-even tables
# Returns canonical columns (pr_best_arm, pr_positive_net_value, pr_incremental_step_worthwhile,
# expected_net_value, net_value_p30/p50/p70) alongside shorthand aliases (p_best_gross, p_best_net,
# p_beats_control, p_positive_step, p_step_roi_positive)
dec_df = longbet.multi_arm_decision_table(
    rec,
    arm_costs={"C": 0.0, "B1": 1.0, "B2": 2.0, "A": 3.5},
    incremental_costs={("B1", "C"): 1.0, ("B2", "B1"): 1.0, ("A", "B2"): 1.5},
    alpha=1.0,
)
roi_df = longbet.roi_break_even_table(
    rec,
    cost_grid=(1.0, 2.0, 3.5),
    alpha_grid=(0.5, 1.0),
    contrast=("A", "C"),
)

# 3. Multi-outcome multi-arm reconciliation with coherent composite synthesis
multi_rec = longbet.reconcile_multi_arm_multi_outcome_draws(
    results_by_outcome={
        "core_rev": fits_core,
        "addon_rev": fits_addon,
        "retention_days": fits_ret,
    },
    arms=arms,
    topology="chain",
    prior_scales={"core_rev": 0.30, "addon_rev": 0.15, "retention_days": 0.25},
    cross_outcome_corr=corr_res,
    composite_outcomes={"total_rev": ("core_rev", "addon_rev")},
    composite_mode="coherent_total",
)
print(multi_rec.max_abs_composite_discrepancy)    # {'total_rev': < 1e-12}
print(multi_rec.additivity_violation_rate)        # {'total_rev': 0.0}
```
