# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Multi-outcome Bayesian Rank-K reconciliation, pre-period IW correlation, and revenue aggregation."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Any

import numpy as np

from longbet._rank_k import RankKLongBetResult
from longbet._rank_k_multi_arm import (
    MultiArmReconciliationResult,
    SamplerDiagnostics,
    _psd_clip,
    _reconcile_joint_bayesian_core,
    _split_rhat_and_ess_1d,
    _validate_prior_scale_scalar,
)

__all__ = [
    "CrossOutcomeCorrelationResult",
    "MultiArmMultiOutcomeReconciliationResult",
    "Nested3ArmMultiReconciliationResult",
    "RankKConvergenceCheck",
    "check_rank_k_convergence",
    "estimate_cross_outcome_correlation",
    "ge_revenue_draws",
    "paid_winsor_cap",
    "reconcile_multi_arm_multi_outcome_draws",
    "reconcile_nested_3arm_multi_draws",
]


@dataclass(frozen=True)
class CrossOutcomeCorrelationResult:
    """Sufficient statistics of within-unit cross-outcome residuals for the conjugate IW prior on ``R``."""

    outcomes: tuple[str, ...]
    correlation: np.ndarray
    scatter_matrix: np.ndarray
    df: int
    n_units: int
    per_arm_correlation: dict[str, np.ndarray]
    min_eigenvalue_after_floor: float


@dataclass(frozen=True)
class RankKConvergenceCheck:
    """Mode-by-mode split-Rhat and bulk ESS check on a two-arm :class:`RankKLongBetResult`."""

    converged: bool
    rhat_modes: tuple[float, ...]
    ess_modes: tuple[float, ...]
    failed_modes: tuple[int, ...]
    rhat_max: float
    ess_min: float


@dataclass
class MultiArmMultiOutcomeReconciliationResult:
    """Posterior draws and diagnostics for a jointly reconciled multi-arm multi-outcome study."""

    by_outcome: dict[str, MultiArmReconciliationResult]
    arms: tuple[str, ...]
    order: tuple[str, ...] | None
    cover_edges: tuple[tuple[str, str], ...]
    topology: str
    tree_edges: tuple[tuple[str, str], ...]
    fitted_pairs: tuple[tuple[str, str], ...]
    redundant_pairs: tuple[tuple[str, str], ...]
    outcomes: tuple[str, ...]
    structural_outcomes: tuple[str, ...]
    composite_outcomes: dict[str, tuple[str, ...]]
    horizon: int
    n_modes: int
    basis: np.ndarray
    arm_trajectory_draws: dict[str, dict[str, np.ndarray]]
    tau_mu_draws: dict[str, np.ndarray]
    tau_increment_draws: dict[str, np.ndarray]
    contrast_mode1_draws: dict[str, dict[tuple[str, str], np.ndarray]]
    contrast_trajectory_draws: dict[str, dict[tuple[str, str], np.ndarray]]
    contrast_catt_draws: dict[str, dict[tuple[str, str], np.ndarray]]
    contrast_cumulative_draws: dict[str, dict[tuple[str, str], np.ndarray]]
    posterior_correlation_draws: np.ndarray
    ppc_cycle_pvalue: dict[str, float | None]
    ppc_composite_pvalue: dict[str, float | None]
    ppc_basis_pvalue: dict[str, float | None]
    contrast_scale_inflation: dict[str, dict[tuple[str, str], float]]
    posterior_corruption_prob: dict[str, dict[tuple[str, str], float]]
    additivity_violation_rate: dict[str, float]
    max_abs_composite_discrepancy: dict[str, float]
    sampler_diagnostics: SamplerDiagnostics

    @property
    def diagnostics(self) -> SamplerDiagnostics:
        """Alias for :attr:`sampler_diagnostics`."""
        return self.sampler_diagnostics

    @property
    def posterior_correlation_mean(self) -> np.ndarray:
        """Posterior mean cross-outcome correlation matrix ``E[R | data]`` of shape ``(M, M)``."""
        return np.mean(self.posterior_correlation_draws, axis=0)


Nested3ArmMultiReconciliationResult = MultiArmMultiOutcomeReconciliationResult


# ---- Utility functions: winsorization, convergence check, GE revenue, IW correlation ------------


def paid_winsor_cap(y_pre: np.ndarray, *, q: float = 0.995) -> float:
    """Compute positive-cell pre-period quantile cap on ``y_pre`` (0.0 when no positive cells exist)."""
    q_val = float(q)
    if not (0.0 < q_val < 1.0):
        raise ValueError(f"q must lie strictly in (0, 1); got {q!r}")
    arr = np.asarray(y_pre, dtype=np.float64)
    pos = arr[np.isfinite(arr) & (arr > 0.0)]
    if pos.size == 0:
        return 0.0
    return float(np.quantile(pos, q_val))


def check_rank_k_convergence(
    res: RankKLongBetResult,
    *,
    num_chains: int = 4,
    rhat_max: float = 1.05,
    ess_min: float = 100.0,
) -> RankKConvergenceCheck:
    """Check split-Rhat and bulk ESS across all exposure modes of a Stage-1 :class:`RankKLongBetResult`."""
    n_c = max(int(num_chains), 1)
    rhat_lim = float(rhat_max)
    ess_lim = float(ess_min)

    cohorts = getattr(res, "cohorts", None)
    n_units = int(res.catt_draws.shape[1])
    n_periods = int(res.catt_draws.shape[2])
    if cohorts is not None:
        ever = np.asarray(cohorts, dtype=np.int32) < n_periods
    else:
        ever = np.any(np.abs(res.catt_draws) > 0.0, axis=(0, 2))
    if not np.any(ever):
        ever = np.ones(n_units, dtype=bool)

    rhats: list[float] = []
    esses: list[float] = []
    failed: list[int] = []
    for k_idx, m_draws in enumerate(res.mode_tau_draws):
        arr = np.asarray(m_draws, dtype=np.float64)
        trace_1d = np.mean(arr[ever], axis=0)
        s_len = int(trace_1d.size)
        per_c = s_len // n_c
        if n_c >= 2 and per_c >= 8:
            chains_2d = trace_1d[: n_c * per_c].reshape(n_c, per_c)
            rh, ess_b, _ = _split_rhat_and_ess_1d(chains_2d)
        else:
            rh = float(res.rhat_att) if (k_idx == 0 and res.rhat_att is not None) else 1.0
            ess_b = float(res.ess_att) if (k_idx == 0 and res.ess_att is not None) else float(s_len)
        rhats.append(float(rh))
        esses.append(float(ess_b))
        if (not math.isfinite(rh)) or rh > rhat_lim or ess_b < ess_lim:
            failed.append(k_idx)

    return RankKConvergenceCheck(
        converged=len(failed) == 0,
        rhat_modes=tuple(rhats),
        ess_modes=tuple(esses),
        failed_modes=tuple(failed),
        rhat_max=rhat_lim,
        ess_min=ess_lim,
    )


def ge_revenue_draws(
    draws_by_outcome: Mapping[str, np.ndarray],
    *,
    weights: Mapping[str, float] | None = None,
    horizon: int | None = None,
    discount_rate: float = 0.0,
    periods_per_year: int = 12,
) -> np.ndarray:
    """Compute weighted (and optionally discounted) general-equilibrium revenue draws across outcomes."""
    if not isinstance(draws_by_outcome, Mapping) or not draws_by_outcome:
        raise ValueError("draws_by_outcome must be a non-empty mapping of {outcome: array}")
    rate = float(discount_rate)
    if not (math.isfinite(rate) and rate > -1.0):
        raise ValueError(f"discount_rate must be finite and > -1; got {discount_rate!r}")
    per_yr = int(periods_per_year)
    if per_yr < 1:
        raise ValueError(f"periods_per_year must be >= 1; got {periods_per_year!r}")

    total: np.ndarray | None = None
    for m_name, raw in draws_by_outcome.items():
        arr = np.asarray(raw, dtype=np.float64)
        w_m = 1.0 if weights is None else float(weights.get(m_name, 1.0))
        if arr.ndim == 1:
            contrib = w_m * arr
        elif arr.ndim == 2:
            h_use = arr.shape[1] if horizon is None else int(horizon)
            if not (1 <= h_use <= arr.shape[1]):
                raise ValueError(f"horizon={horizon!r} out of range [1, {arr.shape[1]}] for {m_name!r}")
            disc = (1.0 + rate) ** (-np.arange(h_use, dtype=np.float64) / float(per_yr))
            contrib = w_m * (arr[:, :h_use] @ disc)
        else:
            raise ValueError(f"draws_by_outcome[{m_name!r}] must be 1-D or 2-D; got shape {arr.shape}")
        if total is None:
            total = contrib.copy()
        else:
            if contrib.shape != total.shape:
                raise ValueError(f"Shape mismatch across outcomes: {contrib.shape} vs {total.shape}")
            total = total + contrib
    assert total is not None
    return total


def estimate_cross_outcome_correlation(
    y_post: Mapping[str, np.ndarray] | np.ndarray,
    y_pre: Mapping[str, np.ndarray] | np.ndarray | None = None,
    *,
    x_covs: np.ndarray | None = None,
    arm_labels: Sequence[Any] | np.ndarray | None = None,
    outcomes: Sequence[str] | None = None,
    min_eigenvalue: float = 1e-4,
) -> CrossOutcomeCorrelationResult:
    """Extract covariate-adjusted within-unit residual scatter and conjugate IW correlation summary."""
    min_eig = float(min_eigenvalue)
    if not (math.isfinite(min_eig) and 0.0 < min_eig < 1.0):
        raise ValueError(f"min_eigenvalue must lie in (0, 1); got {min_eigenvalue!r}")

    if isinstance(y_post, Mapping):
        if not y_post:
            raise ValueError("y_post mapping must be non-empty")
        out_names = tuple(str(k) for k in (outcomes if outcomes is not None else y_post.keys()))
        cols_post = []
        for k in out_names:
            v = np.asarray(y_post[k], dtype=np.float64)
            cols_post.append(np.nanmean(v, axis=1) if v.ndim == 2 else v.reshape(-1))
        ymat = np.column_stack(cols_post)
    else:
        arr = np.asarray(y_post, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] < 2:
            raise ValueError(f"y_post array must have shape (N, M) with M >= 2; got {arr.shape}")
        ymat = arr
        out_names = (
            tuple(str(k) for k in outcomes)
            if outcomes is not None
            else tuple(f"outcome_{i}" for i in range(ymat.shape[1]))
        )
        if len(out_names) != ymat.shape[1]:
            raise ValueError(f"outcomes length {len(out_names)} != {ymat.shape[1]}")

    n_units, m_dim = ymat.shape
    if n_units < m_dim + 2:
        raise ValueError(f"Need at least M + 2 = {m_dim + 2} units; got {n_units}")

    design_cols = [np.ones((n_units, 1), dtype=np.float64)]
    if x_covs is not None:
        x_arr = np.asarray(x_covs, dtype=np.float64)
        if x_arr.ndim == 1:
            x_arr = x_arr[:, None]
        if x_arr.shape[0] != n_units:
            raise ValueError(f"x_covs has {x_arr.shape[0]} rows; expected {n_units}")
        design_cols.append(x_arr)
    if arm_labels is not None:
        arms_arr = np.asarray(arm_labels)
        if arms_arr.shape[0] != n_units:
            raise ValueError(f"arm_labels length {arms_arr.shape[0]} != {n_units}")
        uniq = list(dict.fromkeys(arms_arr.tolist()))
        if len(uniq) > 1:
            dummies = np.column_stack([(arms_arr == u).astype(np.float64) for u in uniq[1:]])
            design_cols.append(dummies)
    else:
        arms_arr = None

    pre_mat: np.ndarray | None = None
    if y_pre is not None:
        if isinstance(y_pre, Mapping):
            pcols = []
            for k in out_names:
                pv = np.asarray(y_pre[k], dtype=np.float64)
                pcols.append(np.nanmean(pv, axis=1) if pv.ndim == 2 else pv.reshape(-1))
            pre_mat = np.column_stack(pcols)
        else:
            pre_mat = np.asarray(y_pre, dtype=np.float64)
            if pre_mat.shape != ymat.shape:
                raise ValueError(f"y_pre shape {pre_mat.shape} != y_post shape {ymat.shape}")

    resid = np.zeros_like(ymat)
    base_x = np.hstack(design_cols)
    p_rank = 1
    for m_i in range(m_dim):
        xm = np.column_stack([base_x, pre_mat[:, m_i]]) if pre_mat is not None else base_x
        p_rank = max(p_rank, int(np.linalg.matrix_rank(xm)))
        coef, *_ = np.linalg.lstsq(xm, ymat[:, m_i], rcond=None)
        resid[:, m_i] = ymat[:, m_i] - xm @ coef

    df_eff = max(int(n_units - p_rank), 1)

    def _corr_from_rows(rmat: np.ndarray, df_r: int) -> tuple[np.ndarray, np.ndarray]:
        centered = rmat - np.mean(rmat, axis=0, keepdims=True)
        sds = np.std(centered, axis=0, ddof=1) if rmat.shape[0] > 1 else np.ones(m_dim)
        sds = np.where(sds > 1e-12, sds, 1.0)
        std_r = centered / sds[None, :]
        scatter = std_r.T @ std_r
        # Conjugate Barnard-McCulloch-Meng IW(M + 1 + df, I_M + scatter) posterior mean correlation
        post_scale = np.eye(m_dim, dtype=np.float64) + scatter
        d_inv = 1.0 / np.sqrt(np.maximum(np.diag(post_scale), 1e-24))
        corr = post_scale * (d_inv[:, None] * d_inv[None, :])
        corr = _psd_clip(corr, floor=min_eig)
        d_inv2 = 1.0 / np.sqrt(np.maximum(np.diag(corr), 1e-24))
        corr = corr * (d_inv2[:, None] * d_inv2[None, :])
        np.fill_diagonal(corr, 1.0)
        corr = _psd_clip(corr, floor=min_eig)
        return corr, scatter

    corr_all, scatter_all = _corr_from_rows(resid, df_eff)
    per_arm: dict[str, np.ndarray] = {}
    if arms_arr is not None:
        for u in dict.fromkeys(arms_arr.tolist()):
            mask_u = arms_arr == u
            if int(np.sum(mask_u)) >= 3:
                cu, _ = _corr_from_rows(resid[mask_u], max(int(np.sum(mask_u)) - 1, 1))
                per_arm[str(u)] = cu
            else:
                per_arm[str(u)] = np.eye(m_dim, dtype=np.float64)

    min_ev = float(np.min(np.linalg.eigvalsh(corr_all)))
    return CrossOutcomeCorrelationResult(
        outcomes=out_names,
        correlation=corr_all,
        scatter_matrix=scatter_all,
        df=df_eff,
        n_units=n_units,
        per_arm_correlation=per_arm,
        min_eigenvalue_after_floor=min_ev,
    )


def _resolve_iw_spec(
    cross_outcome_residuals: np.ndarray | CrossOutcomeCorrelationResult | None,
    structural_outcomes: tuple[str, ...],
) -> dict[str, Any] | None:
    """Build Inverse-Wishart ``(df_post, scale_post, r_init)`` from ``cross_outcome_residuals``."""
    m_dim = len(structural_outcomes)
    if cross_outcome_residuals is None or m_dim < 2:
        return None

    if isinstance(cross_outcome_residuals, CrossOutcomeCorrelationResult):
        res_c = cross_outcome_residuals
        if set(structural_outcomes) <= set(res_c.outcomes):
            idx = [res_c.outcomes.index(m) for m in structural_outcomes]
            scatter = np.asarray(res_c.scatter_matrix, dtype=np.float64)[np.ix_(idx, idx)]
            r_init = np.asarray(res_c.correlation, dtype=np.float64)[np.ix_(idx, idx)]
        elif res_c.scatter_matrix.shape == (m_dim, m_dim):
            scatter = np.asarray(res_c.scatter_matrix, dtype=np.float64)
            r_init = np.asarray(res_c.correlation, dtype=np.float64)
        else:
            raise ValueError(
                f"CrossOutcomeCorrelationResult outcomes {res_c.outcomes} do not match {structural_outcomes}"
            )
        df_pre = max(int(res_c.df), 1)
        scale_post = np.eye(m_dim, dtype=np.float64) + _psd_clip(scatter, floor=0.0)
        return {
            "df_post": float(m_dim + 1 + df_pre),
            "scale_post": scale_post,
            "r_init": _psd_clip(r_init, floor=1e-6),
        }

    arr = np.asarray(cross_outcome_residuals, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != m_dim:
        raise ValueError(
            f"cross_outcome_residuals must be a 2-D (N_units, {m_dim}) residual matrix or "
            f"CrossOutcomeCorrelationResult; got shape {arr.shape}"
        )
    valid_rows = np.all(np.isfinite(arr), axis=1)
    clean = arr[valid_rows]
    if clean.shape[0] < 3:
        raise ValueError(f"cross_outcome_residuals needs at least 3 finite rows; got {clean.shape[0]}")
    centered = clean - np.mean(clean, axis=0, keepdims=True)
    sds = np.std(centered, axis=0, ddof=1)
    sds = np.where(sds > 1e-12, sds, 1.0)
    std_r = centered / sds[None, :]
    df_pre = int(clean.shape[0] - 1)
    scatter = std_r.T @ std_r
    scale_post = np.eye(m_dim, dtype=np.float64) + scatter
    d_inv = 1.0 / np.sqrt(np.maximum(np.diag(scale_post), 1e-24))
    r_init = scale_post * (d_inv[:, None] * d_inv[None, :])
    np.fill_diagonal(r_init, 1.0)
    return {
        "df_post": float(m_dim + 1 + df_pre),
        "scale_post": scale_post,
        "r_init": _psd_clip(r_init, floor=1e-6),
    }


def reconcile_multi_arm_multi_outcome_draws(
    results_by_outcome: Mapping[str, Mapping[tuple[str, str], RankKLongBetResult]],
    *,
    composite_outcomes: Mapping[str, Sequence[str]] | None = None,
    composite_results: Mapping[str, Mapping[tuple[str, str], RankKLongBetResult]] | None = None,
    cross_outcome_residuals: np.ndarray | CrossOutcomeCorrelationResult | None = None,
    arms: Sequence[str] | None = None,
    order: Sequence[str] | None = None,
    poset: Sequence[tuple[str, str]] | None = None,
    monotone: bool | Mapping[str, bool] = True,
    prior_scale: float | Mapping[str, float | None] | None = None,
    horizon: int | None = None,
    n_modes: int = 2,
    seed: int = 2026,
    _model_options: Mapping[str, Any] | None = None,
) -> MultiArmMultiOutcomeReconciliationResult:
    """Jointly reconcile pairwise Stage-1 Rank-K fits across arms, outcomes, and composites."""
    if not isinstance(results_by_outcome, Mapping) or not results_by_outcome:
        raise ValueError("`results_by_outcome` must be a non-empty mapping of {outcome: {(u, v): res}}")
    structural_outcomes = tuple(str(m) for m in results_by_outcome.keys())

    comp_specs: dict[str, tuple[str, ...]] = {}
    if composite_outcomes is not None:
        if not isinstance(composite_outcomes, Mapping):
            raise ValueError("`composite_outcomes` must be a mapping of {composite: sequence_of_outcomes}")
        for raw_c, raw_members in composite_outcomes.items():
            c_name = str(raw_c)
            if c_name in structural_outcomes:
                raise ValueError(f"Composite outcome {c_name!r} collides with a structural outcome name")
            members = tuple(str(m) for m in raw_members)
            if not members or any(m not in structural_outcomes for m in members):
                raise ValueError(
                    f"Composite {c_name!r} members {members} must be non-empty subset of {structural_outcomes}"
                )
            comp_specs[c_name] = members

    comp_res_clean: dict[str, Mapping[tuple[str, str], RankKLongBetResult]] = {}
    if composite_results is not None:
        if not isinstance(composite_results, Mapping):
            raise ValueError("`composite_results` must be a mapping or None")
        for raw_c, cmap in composite_results.items():
            c_name = str(raw_c)
            if c_name not in comp_specs:
                raise ValueError(f"`composite_results` key {c_name!r} must be declared in `composite_outcomes`")
            if not isinstance(cmap, Mapping) or not cmap:
                raise ValueError(f"`composite_results[{c_name!r}]` must be a non-empty mapping")
            comp_res_clean[c_name] = cmap

    if isinstance(monotone, Mapping):
        unknown_m = set(str(k) for k in monotone) - set(structural_outcomes)
        if unknown_m:
            raise ValueError(f"`monotone` mapping contains unknown outcomes {sorted(unknown_m)}")
        monotone_map = {m: bool(monotone.get(m, True)) for m in structural_outcomes}
    elif isinstance(monotone, (bool, np.bool_)):
        monotone_map = {m: bool(monotone) for m in structural_outcomes}
    else:
        raise ValueError(f"`monotone` must be a bool or mapping; got {monotone!r}")

    if isinstance(prior_scale, Mapping):
        unknown_ps = set(str(k) for k in prior_scale) - set(structural_outcomes)
        if unknown_ps:
            raise ValueError(f"`prior_scale` mapping contains unknown outcomes {sorted(unknown_ps)}")
        ps_map = {
            m: _validate_prior_scale_scalar(prior_scale.get(m), f"prior_scale[{m!r}]")
            for m in structural_outcomes
        }
    else:
        ps_scalar = _validate_prior_scale_scalar(prior_scale, "prior_scale")
        ps_map = {m: ps_scalar for m in structural_outcomes}

    iw_spec = _resolve_iw_spec(cross_outcome_residuals, structural_outcomes)

    core = _reconcile_joint_bayesian_core(
        structural_outcomes=structural_outcomes,
        results_by_outcome=results_by_outcome,
        composite_outcomes=comp_specs,
        composite_results=comp_res_clean,
        arms=arms,
        order=order,
        poset=poset,
        monotone_by_outcome=monotone_map,
        prior_scale_by_outcome=ps_map,
        iw_spec=iw_spec,
        horizon=horizon,
        n_modes=n_modes,
        seed=seed,
        model_options=_model_options,
    )

    by_out: dict[str, MultiArmReconciliationResult] = core["by_outcome"]
    topo = core["topology"]
    first_rec = by_out[structural_outcomes[0]]
    all_outcomes = (*structural_outcomes, *tuple(comp_specs.keys()))

    avr_map: dict[str, float] = {}
    disc_map: dict[str, float] = {}
    for c_name, members in comp_specs.items():
        max_gap = 0.0
        viol = np.zeros(first_rec.theta_draws.shape[0], dtype=bool)
        for p in first_rec.all_pairs:
            comp_tr = by_out[c_name].contrast_trajectory_draws[p]
            sum_tr = sum(by_out[m].contrast_trajectory_draws[p] for m in members)
            diff_p = np.max(np.abs(comp_tr - sum_tr), axis=1)
            max_gap = max(max_gap, float(np.max(diff_p)))
            viol |= diff_p > 1e-12
        avr_map[c_name] = float(np.mean(viol))
        disc_map[c_name] = max_gap

    return MultiArmMultiOutcomeReconciliationResult(
        by_outcome=by_out,
        arms=topo["arms"],
        order=topo["order"],
        cover_edges=topo["cover_edges"],
        topology=topo["topology"],
        tree_edges=topo["tree_edges"],
        fitted_pairs=topo["fitted_pairs"],
        redundant_pairs=topo["redundant_pairs"],
        outcomes=all_outcomes,
        structural_outcomes=structural_outcomes,
        composite_outcomes=comp_specs,
        horizon=first_rec.horizon,
        n_modes=first_rec.n_modes,
        basis=first_rec.basis,
        arm_trajectory_draws={m: by_out[m].arm_trajectory_draws for m in all_outcomes},
        tau_mu_draws={m: by_out[m].tau_mu_draws for m in all_outcomes},
        tau_increment_draws={m: by_out[m].tau_increment_draws for m in all_outcomes},
        contrast_mode1_draws={m: by_out[m].contrast_mode1_draws for m in all_outcomes},
        contrast_trajectory_draws={m: by_out[m].contrast_trajectory_draws for m in all_outcomes},
        contrast_catt_draws={m: by_out[m].contrast_catt_draws for m in all_outcomes},
        contrast_cumulative_draws={m: by_out[m].contrast_cumulative_draws for m in all_outcomes},
        posterior_correlation_draws=core["corr_draws"],
        ppc_cycle_pvalue={m: by_out[m].ppc_cycle_pvalue for m in all_outcomes},
        ppc_composite_pvalue={m: by_out[m].ppc_composite_pvalue for m in all_outcomes},
        ppc_basis_pvalue={m: by_out[m].ppc_basis_pvalue for m in all_outcomes},
        contrast_scale_inflation={m: by_out[m].contrast_scale_inflation for m in all_outcomes},
        posterior_corruption_prob={m: by_out[m].posterior_corruption_prob for m in all_outcomes},
        additivity_violation_rate=avr_map,
        max_abs_composite_discrepancy=disc_map,
        sampler_diagnostics=core["diagnostics"],
    )


def reconcile_nested_3arm_multi_draws(
    results_bc: Mapping[str, RankKLongBetResult],
    results_ab: Mapping[str, RankKLongBetResult],
    results_ac: Mapping[str, RankKLongBetResult],
    *,
    composite_outcomes: Mapping[str, Sequence[str]] | None = None,
    composite_results: Mapping[str, Mapping[Any, RankKLongBetResult]] | None = None,
    cross_outcome_residuals: np.ndarray | CrossOutcomeCorrelationResult | None = None,
    monotone: bool | Mapping[str, bool] = True,
    prior_scale: float | Mapping[str, float | None] | None = None,
    horizon: int | None = None,
    n_modes: int = 2,
    seed: int = 2026,
    _model_options: Mapping[str, Any] | None = None,
) -> MultiArmMultiOutcomeReconciliationResult:
    """Convenience wrapper reconciling nested 3-arm multi-outcome fits ``(C <= B <= A)``."""
    outcomes = tuple(results_bc.keys())
    if set(results_ab.keys()) != set(outcomes) or set(results_ac.keys()) != set(outcomes):
        raise ValueError("results_bc, results_ab, and results_ac must share identical outcome keys")
    by_out = {m: {("B", "C"): results_bc[m], ("A", "B"): results_ab[m], ("A", "C"): results_ac[m]} for m in outcomes}
    key_map = {"BC": ("B", "C"), "AB": ("A", "B"), "AC": ("A", "C")}
    comp_res = (
        {c: {key_map.get(k, k): v for k, v in cmap.items()} for c, cmap in composite_results.items()}
        if composite_results is not None
        else None
    )
    return reconcile_multi_arm_multi_outcome_draws(
        by_out,
        composite_outcomes=composite_outcomes,
        composite_results=comp_res,
        cross_outcome_residuals=cross_outcome_residuals,
        arms=("C", "B", "A"),
        order=("C", "B", "A"),
        monotone=monotone,
        prior_scale=prior_scale,
        horizon=horizon,
        n_modes=n_modes,
        seed=seed,
        _model_options=_model_options,
    )
