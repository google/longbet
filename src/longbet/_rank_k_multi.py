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

"""Multi-outcome nested 3-arm reconciliation, Iman-Conover copula coupling, and composite synthesis."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Any

import numpy as np
from scipy import special

from longbet._rank_k import (
    _BUNDLE_PRIOR_MODES,
    _bundle_prior_logpdf,
    _draw_bundle_prior,
    _get_rhat_ess,
    _infer_panel_masks_from_result,
    _project_cohort_damped_trajectory,
    _sir_posterior_2d,
    _warn_rank_k,
    Nested3ArmReconciliationResult,
    RankKLongBetResult,
    reconcile_nested_3arm_draws,
)

__all__ = [
    "CrossOutcomeCorrelationResult",
    "ImanConoverResult",
    "Nested3ArmMultiReconciliationResult",
    "RankKConvergenceCheck",
    "check_rank_k_convergence",
    "estimate_cross_outcome_correlation",
    "ge_revenue_draws",
    "iman_conover_couple",
    "paid_winsor_cap",
    "reconcile_nested_3arm_multi_draws",
]

_COMPOSITE_MODES = ("bottom_up", "coherent_total")
_CORRELATION_TYPES = ("spearman", "pearson", "score")


@dataclass(frozen=True)
class CrossOutcomeCorrelationResult:
    """Within-arm ANCOVA residual correlation across ``M`` outcomes with Ledoit-Wolf shrinkage.

    Attributes:
      outcomes: Tuple of ``M`` outcome labels.
      correlation: ``(M, M)`` pooled positive-definite correlation matrix across arms.
      per_arm_correlation: ``{arm: (M, M)}`` shrunk positive-definite correlation matrix per arm.
      shrinkage_intensities: ``{arm: float}`` analytical Ledoit-Wolf shrinkage intensity in ``[0, 1]``.
      effective_sample_sizes: ``{arm: int}`` residual degrees of freedom per arm.
      min_eigenvalue_before_floor: Smallest eigenvalue across arm correlation matrices before flooring.
      min_eigenvalue_after_floor: Smallest eigenvalue across arm and pooled matrices after flooring.
    """

    outcomes: tuple[str, ...]
    correlation: np.ndarray
    per_arm_correlation: dict[str, np.ndarray]
    shrinkage_intensities: dict[str, float]
    effective_sample_sizes: dict[str, int]
    min_eigenvalue_before_floor: float
    min_eigenvalue_after_floor: float

    @property
    def per_arm(self) -> dict[str, np.ndarray]:
        """Alias for ``per_arm_correlation``."""
        return self.per_arm_correlation

    @property
    def pooled(self) -> np.ndarray:
        """Alias for ``correlation``."""
        return self.correlation


@dataclass(frozen=True)
class ImanConoverResult:
    """Output of :func:`iman_conover_couple`.

    Attributes:
      outcomes: Tuple of ``M`` outcome labels.
      coupled_draws: ``{outcome: np.ndarray}`` draw arrays reordered along axis 0 so that
        each outcome's exact marginal draws are preserved while inducing ``target_correlation``.
      permutations: ``{outcome: (S,) int64 array}`` permutation indices applied to each outcome.
      target_correlation: ``(M, M)`` SPD target correlation matrix.
      score_correlation: ``(M, M)`` SPD normal-score correlation matrix used in Cholesky coloring.
      realized_spearman: ``(M, M)`` realized Spearman rank correlation of the coupled scalar keys.
      realized_pearson: ``(M, M)`` realized Pearson correlation of the coupled scalar keys.
    """

    outcomes: tuple[str, ...]
    coupled_draws: dict[str, np.ndarray]
    permutations: dict[str, np.ndarray]
    target_correlation: np.ndarray
    score_correlation: np.ndarray
    realized_spearman: np.ndarray
    realized_pearson: np.ndarray

    @property
    def draws(self) -> dict[str, np.ndarray]:
        """Alias for ``coupled_draws``."""
        return self.coupled_draws


@dataclass(frozen=True)
class Nested3ArmMultiReconciliationResult:
    """Output of :func:`reconcile_nested_3arm_multi_draws`.

    Attributes:
      outcomes: Tuple of all outcome labels (structural outcomes followed by composite outcomes).
      structural_outcomes: Tuple of structural outcome labels.
      composite_outcomes: ``{composite_name: tuple_of_component_names}``.
      composite_mode: ``"bottom_up"`` or ``"coherent_total"``.
      by_outcome: ``{outcome: Nested3ArmReconciliationResult}`` with coupled draw order.
      post_att_bc: ``{outcome: (S,)}`` posterior Mode-1 ``B - C`` ATT draws.
      post_att_ab: ``{outcome: (S,)}`` posterior Mode-1 ``A - B`` ATT draws.
      post_att_ac: ``{outcome: (S,)}`` posterior Mode-1 ``A - C`` ATT draws.
      trt_catt_bc: ``{outcome: (S, T_post)}`` posterior ``B - C`` trajectory draws.
      trt_catt_ab: ``{outcome: (S, T_post)}`` posterior ``A - B`` trajectory draws.
      trt_catt_ac: ``{outcome: (S, T_post)}`` posterior ``A - C`` trajectory draws.
      cum_att_bc: ``{outcome: (S,)}`` cumulative ``B - C`` effect draws over ``horizon``.
      cum_att_ab: ``{outcome: (S,)}`` cumulative ``A - B`` effect draws over ``horizon``.
      cum_att_ac: ``{outcome: (S,)}`` cumulative ``A - C`` effect draws over ``horizon``.
      permutations: ``{outcome: (S,)}`` Iman-Conover permutation indices for structural outcomes.
      target_correlation: ``(M, M)`` target correlation matrix across structural outcomes.
      realized_spearman: ``(M, M)`` realized Spearman correlation across structural outcomes' ``post_att_ac``.
      max_abs_composite_discrepancy: ``{composite_name: float}`` max draw-wise absolute discrepancy
        between the composite and the sum of its components across all contrasts and periods.
      additivity_violation_rate: ``{composite_name: float}`` fraction of draws where
        ``|composite - sum(components)| > 1e-12``.
    """

    outcomes: tuple[str, ...]
    structural_outcomes: tuple[str, ...]
    composite_outcomes: dict[str, tuple[str, ...]]
    composite_mode: str
    by_outcome: dict[str, Nested3ArmReconciliationResult]
    post_att_bc: dict[str, np.ndarray]
    post_att_ab: dict[str, np.ndarray]
    post_att_ac: dict[str, np.ndarray]
    trt_catt_bc: dict[str, np.ndarray]
    trt_catt_ab: dict[str, np.ndarray]
    trt_catt_ac: dict[str, np.ndarray]
    cum_att_bc: dict[str, np.ndarray]
    cum_att_ab: dict[str, np.ndarray]
    cum_att_ac: dict[str, np.ndarray]
    permutations: dict[str, np.ndarray]
    target_correlation: np.ndarray
    realized_spearman: np.ndarray
    max_abs_composite_discrepancy: dict[str, float]
    additivity_violation_rate: dict[str, float]

    def __getitem__(self, outcome: str) -> Nested3ArmReconciliationResult:
        return self.by_outcome[outcome]


@dataclass(frozen=True)
class RankKConvergenceCheck:
    """Convergence diagnostics gate for a :class:`RankKLongBetResult`.

    Attributes:
      rhat_att: Split R-hat of the overall Mode-1 ATT draws.
      ess_att: Effective sample size of the overall Mode-1 ATT draws.
      rhat_modes: Tuple of per-mode treated-mean split R-hat values ``(mode_0, ..., mode_{K-1})``.
      ess_modes: Tuple of per-mode treated-mean bulk ESS values ``(mode_0, ..., mode_{K-1})``.
      converged: ``True`` iff ATT and all modes satisfy ``rhat <= rhat_max`` and ``ess >= ess_min``.
      failed_modes: 0-based mode indices ``k`` that failed the R-hat or ESS gate.
    """

    rhat_att: float
    ess_att: float
    rhat_modes: tuple[float, ...]
    ess_modes: tuple[float, ...]
    converged: bool
    failed_modes: tuple[int, ...]


def _floor_correlation_eigenvalues(
    matrix: np.ndarray, min_eigenvalue: float
) -> tuple[np.ndarray, float, float]:
    """Project a symmetric matrix to a unit-diagonal SPD correlation matrix with eigenvalue floor."""
    sym = 0.5 * (np.asarray(matrix, dtype=np.float64) + np.asarray(matrix, dtype=np.float64).T)
    np.fill_diagonal(sym, 1.0)
    eigvals, eigvecs = np.linalg.eigh(sym)
    min_before = float(np.min(eigvals))
    if min_before < min_eigenvalue:
        floored = np.maximum(eigvals, float(min_eigenvalue))
        rebuilt = (eigvecs * floored[None, :]) @ eigvecs.T
        scale = np.sqrt(np.maximum(np.diag(rebuilt), 1e-12))
        sym = rebuilt / np.outer(scale, scale)
        sym = 0.5 * (sym + sym.T)
        np.fill_diagonal(sym, 1.0)
        eigvals_after = np.linalg.eigvalsh(sym)
        min_after = float(np.min(eigvals_after))
    else:
        min_after = min_before
    return sym, min_before, min_after


def _ledoit_wolf_correlation_shrinkage(
    standardized_residuals: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Analytical Ledoit-Wolf shrinkage of a sample correlation matrix toward ``I_M``."""
    n_obs, n_outcomes = standardized_residuals.shape
    if n_outcomes == 1:
        return np.ones((1, 1), dtype=np.float64), 0.0
    sample_corr = (standardized_residuals.T @ standardized_residuals) / float(n_obs)
    sample_corr = 0.5 * (sample_corr + sample_corr.T)
    np.fill_diagonal(sample_corr, 1.0)
    off_diag_sq = float(np.sum((sample_corr - np.eye(n_outcomes)) ** 2))
    if off_diag_sq <= 1e-16:
        return np.eye(n_outcomes, dtype=np.float64), 1.0
    z_sq = standardized_residuals**2
    second_moment_sum = float(np.sum((z_sq.T @ z_sq) / float(n_obs) - sample_corr**2))
    shrinkage = float(np.clip(second_moment_sum / (float(n_obs) * off_diag_sq), 0.0, 1.0))
    shrunk = (1.0 - shrinkage) * sample_corr + shrinkage * np.eye(n_outcomes, dtype=np.float64)
    np.fill_diagonal(shrunk, 1.0)
    return 0.5 * (shrunk + shrunk.T), shrinkage


def _coerce_outcome_matrix(
    data: Mapping[str, np.ndarray] | np.ndarray,
    outcomes: Sequence[str] | None,
    name: str,
) -> tuple[tuple[str, ...], np.ndarray]:
    """Convert ``{outcome: (N,) or (N, T)}`` or ``(N, M)`` into ``(outcomes, (N, M) array)``."""
    if isinstance(data, Mapping):
        if not data:
            raise ValueError(f"{name} must not be empty")
        keys = tuple(str(k) for k in data.keys()) if outcomes is None else tuple(str(k) for k in outcomes)
        cols = []
        for k in keys:
            if k not in data:
                raise ValueError(f"{name} is missing outcome {k!r}")
            arr = np.asarray(data[k], dtype=np.float64)
            if arr.ndim == 2:
                arr = np.nanmean(arr, axis=1)
            elif arr.ndim != 1:
                raise ValueError(f"{name}[{k!r}] must be 1-D (N,) or 2-D (N, T); got shape {arr.shape}")
            cols.append(arr)
        lengths = {c.shape[0] for c in cols}
        if len(lengths) != 1:
            raise ValueError(f"{name} outcomes must share the same unit count N; got {sorted(lengths)}")
        mat = np.column_stack(cols)
    else:
        mat = np.asarray(data, dtype=np.float64)
        if mat.ndim != 2 or mat.shape[1] < 1:
            raise ValueError(f"{name} array must be 2-D (N, M); got shape {mat.shape}")
        if outcomes is None:
            keys = tuple(f"outcome_{i}" for i in range(mat.shape[1]))
        else:
            keys = tuple(str(k) for k in outcomes)
            if len(keys) != mat.shape[1]:
                raise ValueError(f"outcomes length ({len(keys)}) does not match {name} columns ({mat.shape[1]})")
    if not np.all(np.isfinite(mat)):
        raise ValueError(f"{name} must contain only finite values")
    return keys, mat


def estimate_cross_outcome_correlation(
    y_post_means: Mapping[str, np.ndarray] | np.ndarray,
    y_pre_means: Mapping[str, np.ndarray] | np.ndarray | None = None,
    x_covs: np.ndarray | None = None,
    arm_labels: np.ndarray | Sequence[Any] | None = None,
    *,
    outcomes: Sequence[str] | None = None,
    min_eigenvalue: float = 1e-4,
) -> CrossOutcomeCorrelationResult:
    """Estimate within-arm ANCOVA residual correlation across ``M`` outcomes with Ledoit-Wolf shrinkage.

    Within each arm ``a``, each outcome's post-period unit mean ``y_post_means[:, m]`` is
    residualized by OLS against ``[1, y_pre_means[:, m], x_covs]``. The ``(n_a, M)`` residual
    matrix is standardized per outcome, shrunk toward ``I_M`` via analytical Ledoit-Wolf
    shrinkage, and floored at ``min_eigenvalue``. Arm correlation matrices are pooled across
    arms weighted by residual degrees of freedom.

    Args:
      y_post_means: ``{outcome: (N,) or (N, T_post)}`` mapping or ``(N, M)`` array of post-period means.
      y_pre_means: Optional ``{outcome: (N,) or (N, T_pre)}`` mapping or ``(N, M)`` array of pre-period means.
      x_covs: Optional ``(N, P)`` baseline covariate matrix.
      arm_labels: Optional length-``N`` arm indicators (for example ``"C"``, ``"B"``, ``"A"``).
        When ``None``, all units are treated as a single arm ``"all"``.
      outcomes: Optional explicit outcome ordering.
      min_eigenvalue: Minimum eigenvalue floor in ``(0, 1)`` (default ``1e-4``).

    Returns:
      :class:`CrossOutcomeCorrelationResult`.
    """
    min_eig = float(min_eigenvalue)
    if not (0.0 < min_eig < 1.0):
        raise ValueError(f"min_eigenvalue must lie in (0, 1); got {min_eigenvalue!r}")

    outcome_names, post_mat = _coerce_outcome_matrix(y_post_means, outcomes, "y_post_means")
    n_units, n_outcomes = post_mat.shape
    if n_units < 4:
        raise ValueError(f"estimate_cross_outcome_correlation requires at least 4 units; got N={n_units}")

    pre_mat: np.ndarray | None = None
    if y_pre_means is not None:
        _, pre_mat = _coerce_outcome_matrix(y_pre_means, outcome_names, "y_pre_means")
        if pre_mat.shape != post_mat.shape:
            raise ValueError(f"y_pre_means shape {pre_mat.shape} does not match y_post_means shape {post_mat.shape}")

    cov_mat: np.ndarray | None = None
    if x_covs is not None:
        cov_mat = np.asarray(x_covs, dtype=np.float64)
        if cov_mat.ndim == 1:
            cov_mat = cov_mat[:, None]
        if cov_mat.ndim != 2 or cov_mat.shape[0] != n_units:
            raise ValueError(f"x_covs must have shape ({n_units}, P); got {cov_mat.shape}")
        if not np.all(np.isfinite(cov_mat)):
            raise ValueError("x_covs must contain only finite values")

    if arm_labels is None:
        arms_arr = np.full(n_units, "all", dtype=object)
        unique_arms = ["all"]
    else:
        arms_arr = np.asarray(arm_labels)
        if arms_arr.shape != (n_units,):
            raise ValueError(f"arm_labels must have length {n_units}; got shape {arms_arr.shape}")
        seen: list[str] = []
        for label in arms_arr:
            s_label = str(label)
            if s_label not in seen:
                seen.append(s_label)
        unique_arms = seen
        arms_arr = np.asarray([str(v) for v in arms_arr], dtype=object)

    per_arm_corr: dict[str, np.ndarray] = {}
    shrinkage_by_arm: dict[str, float] = {}
    dof_by_arm: dict[str, int] = {}
    min_before_all = 1.0
    min_after_all = 1.0

    for arm in unique_arms:
        mask = arms_arr == arm
        n_a = int(np.sum(mask))
        if n_a < 3:
            raise ValueError(f"arm {arm!r} has only {n_a} units; need at least 3")
        residuals = np.zeros((n_a, n_outcomes), dtype=np.float64)
        n_predictors = 1
        for m in range(n_outcomes):
            y_m = post_mat[mask, m]
            design_cols = [np.ones(n_a, dtype=np.float64)]
            if pre_mat is not None:
                design_cols.append(pre_mat[mask, m])
            if cov_mat is not None and n_a > len(design_cols) + cov_mat.shape[1] + 2:
                design_cols.append(cov_mat[mask])
            design = np.column_stack(design_cols)
            n_predictors = design.shape[1]
            if n_a <= n_predictors + 1:
                design = np.ones((n_a, 1), dtype=np.float64)
                n_predictors = 1
            coef, _, _, _ = np.linalg.lstsq(design, y_m, rcond=None)
            res_m = y_m - design @ coef
            res_m -= float(np.mean(res_m))
            sd_m = float(np.std(res_m, ddof=1))
            if sd_m <= 1e-12:
                raise ValueError(f"outcome {outcome_names[m]!r} has zero ANCOVA residual variance in arm {arm!r}")
            residuals[:, m] = res_m / sd_m

        shrunk_a, lam_a = _ledoit_wolf_correlation_shrinkage(residuals)
        floored_a, min_bef, min_aft = _floor_correlation_eigenvalues(shrunk_a, min_eig)
        per_arm_corr[arm] = floored_a
        shrinkage_by_arm[arm] = lam_a
        dof_by_arm[arm] = max(n_a - n_predictors, 1)
        min_before_all = min(min_before_all, min_bef)
        min_after_all = min(min_after_all, min_aft)

    total_dof = float(sum(dof_by_arm.values()))
    pooled_raw = sum(
        (float(dof_by_arm[arm]) / total_dof) * per_arm_corr[arm] for arm in unique_arms
    )
    pooled_corr, min_bef_pool, min_aft_pool = _floor_correlation_eigenvalues(pooled_raw, min_eig)
    min_before_all = min(min_before_all, min_bef_pool)
    min_after_all = min(min_after_all, min_aft_pool)

    return CrossOutcomeCorrelationResult(
        outcomes=outcome_names,
        correlation=pooled_corr,
        per_arm_correlation=per_arm_corr,
        shrinkage_intensities=shrinkage_by_arm,
        effective_sample_sizes=dof_by_arm,
        min_eigenvalue_before_floor=min_before_all,
        min_eigenvalue_after_floor=min_after_all,
    )


def _ranks(values: np.ndarray) -> np.ndarray:
    """0-based stable ranks ``0..S-1`` of a 1-D vector."""
    order = np.argsort(values, kind="stable")
    ranks = np.empty_like(order)
    ranks[order] = np.arange(values.size, dtype=order.dtype)
    return ranks


def _extract_coupling_key(
    arr: np.ndarray,
    outcome: str,
    key: str | int | Mapping[str, np.ndarray | Callable[[np.ndarray], np.ndarray]] | Callable[[np.ndarray], np.ndarray],
) -> np.ndarray:
    """Extract the 1-D ``(S,)`` scalar key used to rank-order an outcome's draws."""
    n_draws = arr.shape[0]
    if isinstance(key, Mapping):
        if outcome not in key:
            raise ValueError(f"key mapping is missing outcome {outcome!r}")
        spec = key[outcome]
        if callable(spec):
            vec = np.asarray(spec(arr), dtype=np.float64)
        else:
            vec = np.asarray(spec, dtype=np.float64)
    elif callable(key):
        vec = np.asarray(key(arr), dtype=np.float64)
    elif isinstance(key, (int, np.integer)):
        vec = arr if arr.ndim == 1 else arr[:, int(key)]
    elif isinstance(key, str):
        if arr.ndim == 1 or key == "sum":
            vec = arr if arr.ndim == 1 else arr.reshape(n_draws, -1).sum(axis=1)
        elif key.startswith("col") and key[3:].isdigit():
            col_idx = int(key[3:])
            vec = arr[:, col_idx]
        else:
            raise ValueError(f"unknown string key {key!r}; expected 'sum' or 'col<idx>'")
    else:
        raise ValueError(f"unsupported key specification {key!r}")
    vec = np.asarray(vec, dtype=np.float64)
    if vec.shape != (n_draws,) or not np.all(np.isfinite(vec)):
        raise ValueError(f"coupling key for outcome {outcome!r} must be a finite ({n_draws},) vector; got {vec.shape}")
    return vec


def iman_conover_couple(
    draws_by_outcome: Mapping[str, np.ndarray],
    target_corr: np.ndarray | CrossOutcomeCorrelationResult,
    *,
    key: str | int | Mapping[str, np.ndarray | Callable[[np.ndarray], np.ndarray]] | Callable[[np.ndarray], np.ndarray] = "sum",
    seed: int = 0,
    min_eigenvalue: float = 1e-6,
    correlation_type: str = "spearman",
) -> ImanConoverResult:
    """Couple independent per-outcome posterior draws via the Iman-Conover rank copula.

    Permutes the draw index ``s = 0..S-1`` of each outcome so that the scalar summaries
    specified by ``key`` exhibit the target rank correlation ``target_corr`` while preserving
    each outcome's marginal distribution (mean, variance, quantiles, empirical CDF) exactly.

    Args:
      draws_by_outcome: ``{outcome: (S, ...)}`` mapping of posterior draw arrays.
      target_corr: ``(M, M)`` correlation matrix or :class:`CrossOutcomeCorrelationResult`.
      key: Scalar key per draw used to rank each outcome: ``"sum"`` (sum over trailing axes,
        e.g. ``theta_BC + theta_AB = tau_AC``), ``"col0"`` / integer column index, a callable,
        or a ``{outcome: 1D_array_or_callable}`` mapping.
      seed: Random seed for the reference Van der Waerden normal scores.
      min_eigenvalue: Minimum eigenvalue floor for Cholesky factorization.
      correlation_type: ``"spearman"`` (applies ``2 * sin(pi * R / 6)`` normal-copula map)
        or ``"pearson"`` / ``"score"``.

    Returns:
      :class:`ImanConoverResult`.
    """
    if not isinstance(draws_by_outcome, Mapping) or not draws_by_outcome:
        raise ValueError("draws_by_outcome must be a non-empty mapping of {outcome: draws}")
    if correlation_type not in _CORRELATION_TYPES:
        raise ValueError(f"correlation_type must be one of {_CORRELATION_TYPES}; got {correlation_type!r}")
    min_eig = float(min_eigenvalue)
    if not (0.0 < min_eig < 1.0):
        raise ValueError(f"min_eigenvalue must lie in (0, 1); got {min_eigenvalue!r}")

    names = tuple(str(k) for k in draws_by_outcome.keys())
    n_outcomes = len(names)
    arrays: dict[str, np.ndarray] = {}
    keys: dict[str, np.ndarray] = {}
    n_draws_set: set[int] = set()
    for name, raw in zip(names, draws_by_outcome.values(), strict=True):
        arr = np.asarray(raw, dtype=np.float64).copy()
        if arr.ndim < 1 or arr.shape[0] < 4:
            raise ValueError(f"draws_by_outcome[{name!r}] must have at least 4 draws along axis 0; got {arr.shape}")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"draws_by_outcome[{name!r}] contains NaN or Inf values")
        arrays[name] = arr
        n_draws_set.add(int(arr.shape[0]))
        keys[name] = _extract_coupling_key(arr, name, key)

    if len(n_draws_set) != 1:
        raise ValueError(f"all outcomes in draws_by_outcome must share the same draw count S; got {sorted(n_draws_set)}")
    (n_draws,) = n_draws_set

    if isinstance(target_corr, CrossOutcomeCorrelationResult):
        if set(target_corr.outcomes) != set(names):
            raise ValueError(
                f"target_corr outcomes {target_corr.outcomes} do not match draws_by_outcome {names}"
            )
        order = [target_corr.outcomes.index(name) for name in names]
        raw_corr = target_corr.correlation[np.ix_(order, order)]
    else:
        raw_corr = np.asarray(target_corr, dtype=np.float64)

    if raw_corr.shape != (n_outcomes, n_outcomes) or not np.all(np.isfinite(raw_corr)):
        raise ValueError(f"target_corr must be a finite ({n_outcomes}, {n_outcomes}) matrix; got {raw_corr.shape}")
    if not np.allclose(raw_corr, raw_corr.T, atol=1e-8) or not np.allclose(np.diag(raw_corr), 1.0, atol=1e-6):
        raise ValueError("target_corr must be symmetric with unit diagonal")
    if np.any(np.abs(raw_corr) > 1.0 + 1e-8):
        raise ValueError("target_corr entries must lie in [-1, 1]")

    if n_outcomes == 1:
        only = names[0]
        perm = np.arange(n_draws, dtype=np.int64)
        one = np.ones((1, 1), dtype=np.float64)
        return ImanConoverResult(
            outcomes=names,
            coupled_draws={only: arrays[only].copy()},
            permutations={only: perm},
            target_correlation=one,
            score_correlation=one.copy(),
            realized_spearman=one.copy(),
            realized_pearson=one.copy(),
        )

    target_spd, _, _ = _floor_correlation_eigenvalues(np.clip(raw_corr, -1.0, 1.0), min_eig)
    if correlation_type == "spearman":
        sin_mapped = 2.0 * np.sin((math.pi / 6.0) * np.clip(target_spd, -1.0, 1.0))
        np.fill_diagonal(sin_mapped, 1.0)
        score_target, _, _ = _floor_correlation_eigenvalues(sin_mapped, min_eig)
    else:
        score_target = target_spd.copy()

    rng = np.random.default_rng(int(seed))
    vdw_scores = special.ndtri(np.arange(1, n_draws + 1, dtype=np.float64) / float(n_draws + 1))
    score_matrix = np.column_stack([rng.permutation(vdw_scores) for _ in names])
    emp_corr, _, _ = _floor_correlation_eigenvalues(np.corrcoef(score_matrix, rowvar=False), min_eig)
    chol_emp = np.linalg.cholesky(emp_corr)
    chol_target = np.linalg.cholesky(score_target)
    coupled_scores = score_matrix @ np.linalg.solve(chol_emp.T, chol_target.T)

    permutations: dict[str, np.ndarray] = {}
    for col_idx, name in enumerate(names):
        key_order = np.argsort(keys[name], kind="stable")
        permutations[name] = key_order[_ranks(coupled_scores[:, col_idx])].astype(np.int64)

    # Anchor outcome 0 to its original index order 0..S-1 so outcome 0 is unpermuted
    anchor_inv = np.argsort(permutations[names[0]], kind="stable")
    permutations = {name: perm[anchor_inv] for name, perm in permutations.items()}

    coupled_draws = {name: arrays[name][permutations[name]] for name in names}
    coupled_keys = np.column_stack([keys[name][permutations[name]] for name in names])
    coupled_ranks = np.column_stack([_ranks(coupled_keys[:, i]) for i in range(n_outcomes)])
    realized_spearman = np.corrcoef(coupled_ranks, rowvar=False)
    realized_pearson = np.corrcoef(coupled_keys, rowvar=False)

    return ImanConoverResult(
        outcomes=names,
        coupled_draws=coupled_draws,
        permutations=permutations,
        target_correlation=target_spd,
        score_correlation=score_target,
        realized_spearman=realized_spearman,
        realized_pearson=realized_pearson,
    )


def _permute_reconciliation_result(
    rec: Nested3ArmReconciliationResult, perm: np.ndarray
) -> Nested3ArmReconciliationResult:
    """Reorder all draw arrays inside a :class:`Nested3ArmReconciliationResult` by ``perm``."""
    return Nested3ArmReconciliationResult(
        trt_catt_bc=rec.trt_catt_bc[perm],
        trt_catt_ab=rec.trt_catt_ab[perm],
        trt_catt_ac=rec.trt_catt_ac[perm],
        post_att_bc=rec.post_att_bc[perm],
        post_att_ab=rec.post_att_ab[perm],
        post_att_ac=rec.post_att_ac[perm],
        prior_att_bc=None if rec.prior_att_bc is None else rec.prior_att_bc[perm],
        prior_att_ab=None if rec.prior_att_ab is None else rec.prior_att_ab[perm],
        prior_att_ac=None if rec.prior_att_ac is None else rec.prior_att_ac[perm],
        prior_catt_bc=None if rec.prior_catt_bc is None else rec.prior_catt_bc[perm],
        prior_catt_ab=None if rec.prior_catt_ab is None else rec.prior_catt_ab[perm],
        prior_catt_ac=None if rec.prior_catt_ac is None else rec.prior_catt_ac[perm],
        lik_mean=None if rec.lik_mean is None else rec.lik_mean.copy(),
        lik_cov=None if rec.lik_cov is None else rec.lik_cov.copy(),
        prior_mode=rec.prior_mode,
        prior_scale=rec.prior_scale,
        sir_ess=rec.sir_ess,
    )


def _backsolve_group_variances(
    res_bc: RankKLongBetResult,
    res_ab: RankKLongBetResult,
    res_ac: RankKLongBetResult,
    min_group_var_share: float = 0.15,
) -> tuple[float, float, float]:
    """Back-solve ``(v_A, v_B, v_C)`` from the three pairwise fits' ``att_draws`` variances."""
    v_bc = max(float(np.var(np.asarray(res_bc.att_draws, dtype=np.float64), ddof=1)), 1e-12)
    v_ab = max(float(np.var(np.asarray(res_ab.att_draws, dtype=np.float64), ddof=1)), 1e-12)
    v_ac = max(float(np.var(np.asarray(res_ac.att_draws, dtype=np.float64), ddof=1)), 1e-12)
    share = float(min_group_var_share)
    v_b = max(0.5 * (v_bc + v_ab - v_ac), share * min(v_bc, v_ab))
    v_c = max(v_bc - v_b, share * v_bc)
    v_a = max(v_ab - v_b, share * v_ab)
    return v_a, v_b, v_c


def _resolve_target_correlation(
    cross_outcome_corr: np.ndarray | CrossOutcomeCorrelationResult | None,
    structural_outcomes: tuple[str, ...],
    group_vars: Mapping[str, tuple[float, float, float]],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Resolve the ``(M, M)`` AC target correlation and ``{arm: (M, M)}`` arm correlations."""
    m_count = len(structural_outcomes)
    if cross_outcome_corr is None:
        eye = np.eye(m_count, dtype=np.float64)
        return eye, {"A": eye.copy(), "B": eye.copy(), "C": eye.copy()}

    if isinstance(cross_outcome_corr, CrossOutcomeCorrelationResult):
        if set(cross_outcome_corr.outcomes) != set(structural_outcomes):
            raise ValueError(
                f"cross_outcome_corr outcomes {cross_outcome_corr.outcomes} do not match "
                f"structural_outcomes {structural_outcomes}"
            )
        order = [cross_outcome_corr.outcomes.index(name) for name in structural_outcomes]
        arm_corrs: dict[str, np.ndarray] = {}
        for arm in ("A", "B", "C"):
            if arm in cross_outcome_corr.per_arm_correlation:
                arm_corrs[arm] = cross_outcome_corr.per_arm_correlation[arm][np.ix_(order, order)]
            else:
                arm_corrs[arm] = cross_outcome_corr.correlation[np.ix_(order, order)]
        # Implied AC sampling-error correlation from S_A + S_C
        sd_a = np.array([math.sqrt(group_vars[m][0]) for m in structural_outcomes], dtype=np.float64)
        sd_c = np.array([math.sqrt(group_vars[m][2]) for m in structural_outcomes], dtype=np.float64)
        cov_ac = np.outer(sd_a, sd_a) * arm_corrs["A"] + np.outer(sd_c, sd_c) * arm_corrs["C"]
        scale_ac = np.sqrt(np.maximum(np.diag(cov_ac), 1e-12))
        corr_ac = cov_ac / np.outer(scale_ac, scale_ac)
        target_spd, _, _ = _floor_correlation_eigenvalues(corr_ac, 1e-6)
        return target_spd, arm_corrs

    mat = np.asarray(cross_outcome_corr, dtype=np.float64)
    if mat.shape != (m_count, m_count) or not np.all(np.isfinite(mat)):
        raise ValueError(f"cross_outcome_corr must have shape ({m_count}, {m_count}); got {mat.shape}")
    if not np.allclose(mat, mat.T, atol=1e-8) or not np.allclose(np.diag(mat), 1.0, atol=1e-6):
        raise ValueError("cross_outcome_corr must be symmetric with unit diagonal")
    target_spd, _, _ = _floor_correlation_eigenvalues(np.clip(mat, -1.0, 1.0), 1e-6)
    return target_spd, {"A": target_spd.copy(), "B": target_spd.copy(), "C": target_spd.copy()}


def _sample_synthetic_composite_2d(
    comp_members: Sequence[str],
    structural_outcomes: tuple[str, ...],
    results_bc: Mapping[str, RankKLongBetResult],
    results_ab: Mapping[str, RankKLongBetResult],
    results_ac: Mapping[str, RankKLongBetResult],
    raw_by_outcome: Mapping[str, Nested3ArmReconciliationResult],
    group_vars: Mapping[str, tuple[float, float, float]],
    arm_corrs: Mapping[str, np.ndarray],
    target_corr: np.ndarray,
    *,
    comp_prior_scale: float | str | None,
    prior_mode: str,
    tau_clip: tuple[float, float],
    n_draws: int,
    n_proposals: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Sample single-truncation ``(theta_BC, theta_AB)`` for a composite outcome from component moments."""
    idx_list = [structural_outcomes.index(m) for m in comp_members]
    m_bc = sum(float(np.mean(results_bc[m].att_draws)) for m in comp_members)
    m_ab = sum(float(np.mean(results_ab[m].att_draws)) for m in comp_members)
    m_ac = sum(float(np.mean(results_ac[m].att_draws)) for m in comp_members)

    sub_a = arm_corrs["A"][np.ix_(idx_list, idx_list)]
    sub_b = arm_corrs["B"][np.ix_(idx_list, idx_list)]
    sub_c = arm_corrs["C"][np.ix_(idx_list, idx_list)]
    sd_a = np.array([math.sqrt(group_vars[m][0]) for m in comp_members], dtype=np.float64)
    sd_b = np.array([math.sqrt(group_vars[m][1]) for m in comp_members], dtype=np.float64)
    sd_c = np.array([math.sqrt(group_vars[m][2]) for m in comp_members], dtype=np.float64)

    v_a_comp = max(float(sd_a @ sub_a @ sd_a), 1e-12)
    v_b_comp = max(float(sd_b @ sub_b @ sd_b), 1e-12)
    v_c_comp = max(float(sd_c @ sub_c @ sd_c), 1e-12)
    v_ac_comp = v_a_comp + v_c_comp

    sigma_lik = np.array(
        [[v_b_comp + v_c_comp, -v_b_comp], [-v_b_comp, v_a_comp + v_b_comp]],
        dtype=np.float64,
    )
    delta = (m_ac - (m_bc + m_ab)) / 3.0
    mu_lik = np.array([m_bc + delta, m_ab + delta], dtype=np.float64)

    if comp_prior_scale is None:
        comp_scales = []
        for m in comp_members:
            sc = raw_by_outcome[m].prior_scale
            if isinstance(sc, tuple):
                comp_scales.append(math.hypot(float(sc[0]), float(sc[1])))
            elif sc is not None and math.isfinite(float(sc)):
                comp_scales.append(float(sc))
            else:
                comp_scales.append(math.sqrt(group_vars[m][0] + group_vars[m][2]))
        sc_vec = np.asarray(comp_scales, dtype=np.float64)
        sub_target = target_corr[np.ix_(idx_list, idx_list)]
        sd_ac = math.sqrt(max(float(sc_vec @ sub_target @ sc_vec), 1e-12))
    elif isinstance(comp_prior_scale, str):
        if comp_prior_scale != "in_sample":
            raise ValueError(f"composite prior scale must be > 0 or 'in_sample'; got {comp_prior_scale!r}")
        sd_ac = math.sqrt(max(v_ac_comp, 1e-12))
    else:
        sd_ac = float(comp_prior_scale)
        if not sd_ac > 0.0:
            raise ValueError(f"composite prior scale must be > 0; got {comp_prior_scale!r}")

    rng = np.random.default_rng(int(seed))
    low_clip, high_clip = float(tau_clip[0]), float(tau_clip[1])

    if prior_mode in _BUNDLE_PRIOR_MODES:

        def _prior_draw(n: int) -> np.ndarray:
            tot = _draw_bundle_prior(rng, prior_mode, sd_ac, n)
            u = rng.uniform(0.0, 1.0, size=n)
            return np.column_stack([u * tot, (1.0 - u) * tot])

        def _prior_logpdf(pts: np.ndarray) -> np.ndarray:
            t = pts[:, 0] + pts[:, 1]
            with np.errstate(divide="ignore"):
                return _bundle_prior_logpdf(prior_mode, sd_ac, t) - np.log(np.maximum(t, 1e-300))

    else:
        sd_pair = sd_ac / math.sqrt(2.0)

        def _prior_draw(n: int) -> np.ndarray:
            return sd_pair * rng.standard_normal((n, 2))

        def _prior_logpdf(pts: np.ndarray) -> np.ndarray:
            return -0.5 * np.sum((pts / sd_pair) ** 2, axis=1)

    post_pts, _, sir_ess = _sir_posterior_2d(
        rng,
        prior_draw=_prior_draw,
        prior_logpdf=_prior_logpdf,
        mu_lik=mu_lik,
        sigma_lik=sigma_lik,
        low_clip=low_clip,
        high_clip=high_clip,
        n_draws=n_draws,
        n_prop=max(int(n_proposals), 10000),
    )
    return post_pts[:, 0], post_pts[:, 1], mu_lik, sigma_lik, sd_ac, sir_ess


def _per_outcome_arg(arg: Any, outcome: str, arg_name: str) -> Any:
    """Extracts a per-outcome argument when ``arg`` is a mapping, else returns ``arg``."""
    if isinstance(arg, Mapping):
        if outcome not in arg:
            raise ValueError(f"{arg_name} mapping is missing outcome {outcome!r}")
        return arg[outcome]
    return arg


def reconcile_nested_3arm_multi_draws(
    results_bc: Mapping[str, RankKLongBetResult],
    results_ab: Mapping[str, RankKLongBetResult],
    results_ac: Mapping[str, RankKLongBetResult],
    *,
    prior_scales: Mapping[str, float | str] | float | str = "in_sample",
    prior_mode: str = "bundle_halfnormal",
    share_mean_logit: float = 0.0,
    share_sd_logit: float = 1.5,
    cross_outcome_corr: np.ndarray | CrossOutcomeCorrelationResult | None = None,
    composite_outcomes: Mapping[str, Sequence[str]] | None = None,
    composite_mode: str = "bottom_up",
    composite_results: Mapping[str, Mapping[str, RankKLongBetResult]] | None = None,
    composite_prior_scales: Mapping[str, float | str] | None = None,
    cohorts_bc: np.ndarray | Mapping[str, np.ndarray] | None = None,
    ever_bc: np.ndarray | Mapping[str, np.ndarray] | None = None,
    cohorts_ab: np.ndarray | Mapping[str, np.ndarray] | None = None,
    ever_ab: np.ndarray | Mapping[str, np.ndarray] | None = None,
    tau_clip: tuple[float, float] | None = (0.0, np.inf),
    catt_floor: float | None = 0.0,
    horizon: int | None = None,
    n_proposals: int = 250_000,
    seed: int = 0,
) -> Nested3ArmMultiReconciliationResult:
    """Reconcile nested 3-arm Rank-K BCF posteriors across ``M`` outcomes with copula coupling and composites.

    Implements the ``M1`` multi-outcome reconciliation architecture:
    1. Each structural outcome ``m`` is reconciled across ``(BC, AB, AC)`` via
       :func:`reconcile_nested_3arm_draws` (singular-GLS shared-B-arm reduction plus
       2-D defensive-mixture SIR under the bundle prior and zero-sum dynamic damping).
    2. Draw indices ``s = 0..S-1`` are coupled across the ``M`` structural outcomes via
       :func:`iman_conover_couple` on ``post_att_ac = post_att_bc + post_att_ab``,
       preserving every outcome's marginal posterior distribution while inducing the
       cross-outcome sampling-error correlation ``cross_outcome_corr``.
    3. Any requested ``composite_outcomes`` (e.g. ``{"total": ("outcome_a", "outcome_b")}``) are
       synthesized under ``composite_mode``:
       - ``"bottom_up"``: composite draws and trajectories are the exact draw-wise sums
         of the coupled component outcomes (``AVR == 0``, ``max |gap| < 1e-12``,
         component marginals 100% unchanged).
       - ``"coherent_total"``: reconciles the single-truncation composite total draw once
         (from ``composite_results`` if provided, or from the composite GLS moments and
         combined prior scale), couples it comonotonically with the sum of the coupled
         components, and rescales the non-negative component draws on each draw ``s``
         (``scale_s = total_ac_s / sum_m comp_ac_m_s``) so that both exact draw-wise
         additivity (``AVR == 0``, ``max |gap| < 1e-12``) and single-truncation composite
         null-floor preservation hold simultaneously.

    Args:
      results_bc: ``{outcome: RankKLongBetResult}`` for the ``B - C`` fits.
      results_ab: ``{outcome: RankKLongBetResult}`` for the ``A - B`` fits.
      results_ac: ``{outcome: RankKLongBetResult}`` for the ``A - C`` fits.
      prior_scales: Positive float, ``"in_sample"``, or ``{outcome: scale}`` mapping.
      prior_mode: Bundle prior family (``"bundle_halfnormal"``, ``"bundle_halfstudentt"``,
        ``"bundle_gamma_spike"``, or ``"independent_halfnormal"``).
      share_mean_logit: Prior mean of attribution share logit (reserved for API symmetry).
      share_sd_logit: Prior SD of attribution share logit (reserved for API symmetry).
      cross_outcome_corr: ``None`` (identity), an ``(M, M)`` correlation matrix, or a
        :class:`CrossOutcomeCorrelationResult`.
      composite_outcomes: Optional ``{composite_name: (component_1, component_2, ...)}``.
      composite_mode: ``"bottom_up"`` (default) or ``"coherent_total"``.
      composite_results: Optional ``{composite_name: {"BC": fit_bc, "AB": fit_ab, "AC": fit_ac}}``
        used when ``composite_mode="coherent_total"``.
      composite_prior_scales: Optional ``{composite_name: scale}`` for ``"coherent_total"``.
      cohorts_bc: Optional array or ``{outcome: array}`` of ``B - C`` adoption cohorts.
      ever_bc: Optional array or ``{outcome: array}`` of ``B - C`` ever-treated masks.
      cohorts_ab: Optional array or ``{outcome: array}`` of ``A - B`` adoption cohorts.
      ever_ab: Optional array or ``{outcome: array}`` of ``A - B`` ever-treated masks.
      tau_clip: ``(low, high)`` Mode-1 clip bounds (default ``(0.0, inf)``; ``None`` for ``(-inf, inf)``).
      catt_floor: Optional per-period trajectory floor for zero-sum dynamic damping (default ``0.0``).
      horizon: Optional number of post-adoption periods summed into ``cum_att_*`` (default: all ``T_post``).
      n_proposals: Importance-sampling proposal count per outcome (default ``250_000``).
      seed: Random seed.

    Returns:
      :class:`Nested3ArmMultiReconciliationResult`.
    """
    _ = (share_mean_logit, share_sd_logit)
    if not isinstance(results_bc, Mapping) or not results_bc:
        raise ValueError("results_bc must be a non-empty mapping of {outcome: RankKLongBetResult}")
    structural_outcomes = tuple(str(k) for k in results_bc.keys())
    if set(results_ab.keys()) != set(structural_outcomes) or set(results_ac.keys()) != set(structural_outcomes):
        raise ValueError(
            "results_bc, results_ab, and results_ac must share the exact same outcome keys; "
            f"got {list(results_bc)}, {list(results_ab)}, {list(results_ac)}"
        )
    if composite_mode not in _COMPOSITE_MODES:
        raise ValueError(f"composite_mode must be one of {_COMPOSITE_MODES}; got {composite_mode!r}")

    resolved_tau_clip = (-float("inf"), float("inf")) if tau_clip is None else (float(tau_clip[0]), float(tau_clip[1]))

    # Validate composite_outcomes specification up front
    comp_spec: dict[str, tuple[str, ...]] = {}
    seen_coherent_components: set[str] = set()
    if composite_outcomes is not None:
        for comp_name, members in composite_outcomes.items():
            c_str = str(comp_name)
            if c_str in structural_outcomes:
                raise ValueError(f"composite outcome {c_str!r} collides with a structural outcome name")
            member_tuple = tuple(str(m) for m in members)
            if len(member_tuple) < 2 or len(set(member_tuple)) != len(member_tuple):
                raise ValueError(
                    f"composite outcome {c_str!r} must sum at least 2 distinct components; got {member_tuple}"
                )
            missing = [m for m in member_tuple if m not in structural_outcomes]
            if missing:
                raise ValueError(
                    f"composite outcome {c_str!r} references unknown component(s) {missing}; "
                    f"structural outcomes are {structural_outcomes}"
                )
            if composite_mode == "coherent_total":
                overlap = seen_coherent_components & set(member_tuple)
                if overlap:
                    raise ValueError(
                        f"composite_mode='coherent_total' requires disjoint composite component sets, "
                        f"but composite {c_str!r} shares component(s) {sorted(overlap)} with an earlier composite"
                    )
                seen_coherent_components.update(member_tuple)
            comp_spec[c_str] = member_tuple

    # Step 1: Per-outcome 3-arm reconciliation
    raw_by_outcome: dict[str, Nested3ArmReconciliationResult] = {}
    group_vars: dict[str, tuple[float, float, float]] = {}
    for idx, m in enumerate(structural_outcomes):
        scale_m = _per_outcome_arg(prior_scales, m, "prior_scales")
        c_bc_m = _per_outcome_arg(cohorts_bc, m, "cohorts_bc")
        e_bc_m = _per_outcome_arg(ever_bc, m, "ever_bc")
        c_ab_m = _per_outcome_arg(cohorts_ab, m, "cohorts_ab")
        e_ab_m = _per_outcome_arg(ever_ab, m, "ever_ab")
        raw_by_outcome[m] = reconcile_nested_3arm_draws(
            results_bc[m],
            results_ab[m],
            results_ac[m],
            prior_scale=scale_m,
            prior_mode=prior_mode,
            cohorts_bc=c_bc_m,
            ever_bc=e_bc_m,
            cohorts_ab=c_ab_m,
            ever_ab=e_ab_m,
            tau_clip=resolved_tau_clip,
            catt_floor=catt_floor,
            n_proposals=int(n_proposals),
            seed=int(seed) + 1009 * idx,
            return_prior_draws=True,
        )
        group_vars[m] = _backsolve_group_variances(results_bc[m], results_ab[m], results_ac[m])

    draw_counts = {int(raw_by_outcome[m].post_att_ac.shape[0]) for m in structural_outcomes}
    if len(draw_counts) != 1:
        raise ValueError(f"structural outcomes produced different draw counts S: {sorted(draw_counts)}")
    (n_draws,) = draw_counts

    post_lengths = {int(raw_by_outcome[m].trt_catt_ac.shape[1]) for m in structural_outcomes}
    if len(post_lengths) != 1:
        raise ValueError(f"structural outcomes produced different T_post lengths: {sorted(post_lengths)}")
    (t_post,) = post_lengths

    if horizon is None:
        h_eff = t_post
    else:
        h_eff = int(horizon)
        if not (1 <= h_eff <= t_post):
            raise ValueError(f"horizon must lie in [1, {t_post}]; got {horizon!r}")

    # Step 2: Cross-outcome Iman-Conover coupling on post_att_ac = post_att_bc + post_att_ab
    target_corr, arm_corrs = _resolve_target_correlation(cross_outcome_corr, structural_outcomes, group_vars)
    theta_by_outcome = {
        m: np.column_stack([raw_by_outcome[m].post_att_bc, raw_by_outcome[m].post_att_ab])
        for m in structural_outcomes
    }
    ic_res = iman_conover_couple(
        theta_by_outcome,
        target_corr,
        key="sum",
        seed=int(seed) + 7919,
    )

    by_outcome: dict[str, Nested3ArmReconciliationResult] = {
        m: _permute_reconciliation_result(raw_by_outcome[m], ic_res.permutations[m])
        for m in structural_outcomes
    }

    # Step 3: Synthesize composite outcomes
    max_abs_disc: dict[str, float] = {}
    avr_by_comp: dict[str, float] = {}

    for comp_idx, (comp_name, comp_members) in enumerate(comp_spec.items()):
        if composite_mode == "coherent_total":
            comp_scale_arg = (
                None
                if composite_prior_scales is None
                else composite_prior_scales.get(comp_name)
            )
            if composite_results is not None and comp_name in composite_results:
                c_fits = composite_results[comp_name]
                fit_bc = c_fits.get("BC") or c_fits.get("bc")
                fit_ab = c_fits.get("AB") or c_fits.get("ab")
                fit_ac = c_fits.get("AC") or c_fits.get("ac")
                if fit_bc is None or fit_ab is None or fit_ac is None:
                    raise ValueError(
                        f"composite_results[{comp_name!r}] must contain 'BC', 'AB', and 'AC' fits"
                    )
                if comp_scale_arg is None:
                    comp_scale_arg = (
                        prior_scales
                        if isinstance(prior_scales, (int, float, str))
                        else "in_sample"
                    )
                rec_comp_direct = reconcile_nested_3arm_draws(
                    fit_bc,
                    fit_ab,
                    fit_ac,
                    prior_scale=comp_scale_arg,
                    prior_mode=prior_mode,
                    tau_clip=resolved_tau_clip,
                    catt_floor=catt_floor,
                    n_proposals=int(n_proposals),
                    seed=int(seed) + 20011 + 997 * comp_idx,
                    return_prior_draws=True,
                )
                comp_bc_raw = rec_comp_direct.post_att_bc
                comp_ab_raw = rec_comp_direct.post_att_ab
                comp_lik_mean = rec_comp_direct.lik_mean
                comp_lik_cov = rec_comp_direct.lik_cov
                comp_scale_used = rec_comp_direct.prior_scale
                comp_ess = rec_comp_direct.sir_ess
            else:
                (
                    comp_bc_raw,
                    comp_ab_raw,
                    comp_lik_mean,
                    comp_lik_cov,
                    comp_scale_used,
                    comp_ess,
                ) = _sample_synthetic_composite_2d(
                    comp_members,
                    structural_outcomes,
                    results_bc,
                    results_ab,
                    results_ac,
                    raw_by_outcome,
                    group_vars,
                    arm_corrs,
                    target_corr,
                    comp_prior_scale=comp_scale_arg,
                    prior_mode=prior_mode,
                    tau_clip=resolved_tau_clip,
                    n_draws=n_draws,
                    n_proposals=int(n_proposals),
                    seed=int(seed) + 20011 + 997 * comp_idx,
                )

            comp_ac_raw = comp_bc_raw + comp_ab_raw
            sum_comp_ac = sum(by_outcome[m].post_att_ac for m in comp_members)
            # Couple composite total comonotonically with the sum of coupled components
            comp_order = np.argsort(comp_ac_raw, kind="stable")
            sum_ranks = _ranks(sum_comp_ac)
            perm_comp = comp_order[sum_ranks]
            total_ac_coupled = comp_ac_raw[perm_comp]
            total_bc_coupled = comp_bc_raw[perm_comp]
            total_ab_coupled = comp_ab_raw[perm_comp]

            # Rescale component draws on each draw s so sum_m comp_m^(s) == total^(s)
            positive_sum = sum_comp_ac > 1e-12
            scale_s = np.where(
                positive_sum,
                total_ac_coupled / np.maximum(sum_comp_ac, 1e-12),
                0.0,
            )
            k_members = float(len(comp_members))
            for m in comp_members:
                rec_m = by_outcome[m]
                new_post_bc = np.where(
                    positive_sum,
                    rec_m.post_att_bc * scale_s,
                    total_bc_coupled / k_members,
                )
                new_post_ab = np.where(
                    positive_sum,
                    rec_m.post_att_ab * scale_s,
                    total_ab_coupled / k_members,
                )
                new_post_ac = new_post_bc + new_post_ab
                new_trt_bc = np.where(
                    positive_sum[:, None],
                    rec_m.trt_catt_bc * scale_s[:, None],
                    (total_bc_coupled / k_members)[:, None],
                )
                new_trt_ab = np.where(
                    positive_sum[:, None],
                    rec_m.trt_catt_ab * scale_s[:, None],
                    (total_ab_coupled / k_members)[:, None],
                )
                new_trt_ac = new_trt_bc + new_trt_ab
                by_outcome[m] = Nested3ArmReconciliationResult(
                    trt_catt_bc=new_trt_bc,
                    trt_catt_ab=new_trt_ab,
                    trt_catt_ac=new_trt_ac,
                    post_att_bc=new_post_bc,
                    post_att_ab=new_post_ab,
                    post_att_ac=new_post_ac,
                    prior_att_bc=rec_m.prior_att_bc,
                    prior_att_ab=rec_m.prior_att_ab,
                    prior_att_ac=rec_m.prior_att_ac,
                    prior_catt_bc=rec_m.prior_catt_bc,
                    prior_catt_ab=rec_m.prior_catt_ab,
                    prior_catt_ac=rec_m.prior_catt_ac,
                    lik_mean=rec_m.lik_mean,
                    lik_cov=rec_m.lik_cov,
                    prior_mode=rec_m.prior_mode,
                    prior_scale=rec_m.prior_scale,
                    sir_ess=rec_m.sir_ess,
                )
        else:
            comp_lik_mean = sum(by_outcome[m].lik_mean for m in comp_members if by_outcome[m].lik_mean is not None)
            comp_lik_cov = sum(by_outcome[m].lik_cov for m in comp_members if by_outcome[m].lik_cov is not None)
            comp_scale_used = None
            comp_ess = float(min(by_outcome[m].sir_ess for m in comp_members))

        comp_post_bc = sum(by_outcome[m].post_att_bc for m in comp_members)
        comp_post_ab = sum(by_outcome[m].post_att_ab for m in comp_members)
        comp_post_ac = comp_post_bc + comp_post_ab
        comp_trt_bc = sum(by_outcome[m].trt_catt_bc for m in comp_members)
        comp_trt_ab = sum(by_outcome[m].trt_catt_ab for m in comp_members)
        comp_trt_ac = comp_trt_bc + comp_trt_ab

        comp_pri_att_bc = (
            sum(by_outcome[m].prior_att_bc for m in comp_members)
            if all(by_outcome[m].prior_att_bc is not None for m in comp_members)
            else None
        )
        comp_pri_att_ab = (
            sum(by_outcome[m].prior_att_ab for m in comp_members)
            if all(by_outcome[m].prior_att_ab is not None for m in comp_members)
            else None
        )
        comp_pri_att_ac = (
            comp_pri_att_bc + comp_pri_att_ab
            if comp_pri_att_bc is not None and comp_pri_att_ab is not None
            else None
        )
        comp_pri_catt_bc = (
            sum(by_outcome[m].prior_catt_bc for m in comp_members)
            if all(by_outcome[m].prior_catt_bc is not None for m in comp_members)
            else None
        )
        comp_pri_catt_ab = (
            sum(by_outcome[m].prior_catt_ab for m in comp_members)
            if all(by_outcome[m].prior_catt_ab is not None for m in comp_members)
            else None
        )
        comp_pri_catt_ac = (
            comp_pri_catt_bc + comp_pri_catt_ab
            if comp_pri_catt_bc is not None and comp_pri_catt_ab is not None
            else None
        )

        by_outcome[comp_name] = Nested3ArmReconciliationResult(
            trt_catt_bc=comp_trt_bc,
            trt_catt_ab=comp_trt_ab,
            trt_catt_ac=comp_trt_ac,
            post_att_bc=comp_post_bc,
            post_att_ab=comp_post_ab,
            post_att_ac=comp_post_ac,
            prior_att_bc=comp_pri_att_bc,
            prior_att_ab=comp_pri_att_ab,
            prior_att_ac=comp_pri_att_ac,
            prior_catt_bc=comp_pri_catt_bc,
            prior_catt_ab=comp_pri_catt_ab,
            prior_catt_ac=comp_pri_catt_ac,
            lik_mean=comp_lik_mean if isinstance(comp_lik_mean, np.ndarray) else None,
            lik_cov=comp_lik_cov if isinstance(comp_lik_cov, np.ndarray) else None,
            prior_mode=prior_mode,
            prior_scale=comp_scale_used,
            sir_ess=comp_ess,
        )

        sum_member_trt_ac = sum(by_outcome[m].trt_catt_ac for m in comp_members)
        sum_member_post_ac = sum(by_outcome[m].post_att_ac for m in comp_members)
        draw_gaps = np.maximum(
            np.max(np.abs(comp_trt_ac - sum_member_trt_ac), axis=1),
            np.abs(comp_post_ac - sum_member_post_ac),
        )
        max_abs_disc[comp_name] = float(np.max(draw_gaps))
        avr_by_comp[comp_name] = float(np.mean(draw_gaps > 1e-12))

    all_outcomes = (*structural_outcomes, *comp_spec.keys())
    post_att_bc = {name: by_outcome[name].post_att_bc for name in all_outcomes}
    post_att_ab = {name: by_outcome[name].post_att_ab for name in all_outcomes}
    post_att_ac = {name: by_outcome[name].post_att_ac for name in all_outcomes}
    trt_catt_bc = {name: by_outcome[name].trt_catt_bc for name in all_outcomes}
    trt_catt_ab = {name: by_outcome[name].trt_catt_ab for name in all_outcomes}
    trt_catt_ac = {name: by_outcome[name].trt_catt_ac for name in all_outcomes}
    cum_att_bc = {name: by_outcome[name].trt_catt_bc[:, :h_eff].sum(axis=1) for name in all_outcomes}
    cum_att_ab = {name: by_outcome[name].trt_catt_ab[:, :h_eff].sum(axis=1) for name in all_outcomes}
    cum_att_ac = {name: by_outcome[name].trt_catt_ac[:, :h_eff].sum(axis=1) for name in all_outcomes}

    return Nested3ArmMultiReconciliationResult(
        outcomes=all_outcomes,
        structural_outcomes=structural_outcomes,
        composite_outcomes=comp_spec,
        composite_mode=composite_mode,
        by_outcome=by_outcome,
        post_att_bc=post_att_bc,
        post_att_ab=post_att_ab,
        post_att_ac=post_att_ac,
        trt_catt_bc=trt_catt_bc,
        trt_catt_ab=trt_catt_ab,
        trt_catt_ac=trt_catt_ac,
        cum_att_bc=cum_att_bc,
        cum_att_ab=cum_att_ab,
        cum_att_ac=cum_att_ac,
        permutations=ic_res.permutations,
        target_correlation=ic_res.target_correlation,
        realized_spearman=ic_res.realized_spearman,
        max_abs_composite_discrepancy=max_abs_disc,
        additivity_violation_rate=avr_by_comp,
    )


def ge_revenue_draws(
    reconciled_or_draws: Nested3ArmMultiReconciliationResult | Mapping[str, np.ndarray],
    *,
    weights: Mapping[str, float],
    contrast: str = "AC",
    horizon: int | None = None,
    discount_rate: float = 0.0,
    periods_per_year: int = 12,
) -> np.ndarray:
    """Compute coupled general-equilibrium (GE) revenue draws ``sum_m alpha_m * Cum_m^(s)``.

    Args:
      reconciled_or_draws: A :class:`Nested3ArmMultiReconciliationResult` or a mapping
        ``{outcome: (S,) cumulative draws or (S, T) monthly trajectory draws}``.
      weights: ``{outcome: alpha_m}`` incrementality weights (for example
        ``{"outcome_a": alpha_a, "outcome_b": 1.0}``).
      contrast: Contrast used when ``reconciled_or_draws`` is a
        :class:`Nested3ArmMultiReconciliationResult` (``"BC"``, ``"AB"``, or ``"AC"``).
      horizon: Optional number of periods to sum when inputs are 2-D trajectories.
      discount_rate: Annual discount rate ``r > -1`` applied to 2-D trajectories (default ``0.0``).
      periods_per_year: Periods per year for discounting (default ``12``).

    Returns:
      ``(S,)`` float64 array of GE revenue draws.
    """
    if not isinstance(weights, Mapping) or not weights:
        raise ValueError("weights must be a non-empty mapping of {outcome: weight}")
    rate = float(discount_rate)
    if not (math.isfinite(rate) and rate > -1.0):
        raise ValueError(f"discount_rate must be a finite number > -1; got {discount_rate!r}")
    if int(periods_per_year) < 1:
        raise ValueError(f"periods_per_year must be >= 1; got {periods_per_year!r}")

    if isinstance(reconciled_or_draws, Nested3ArmMultiReconciliationResult):
        c_upper = contrast.upper()
        if c_upper == "AC":
            source: Mapping[str, np.ndarray] = reconciled_or_draws.trt_catt_ac
        elif c_upper == "BC":
            source = reconciled_or_draws.trt_catt_bc
        elif c_upper == "AB":
            source = reconciled_or_draws.trt_catt_ab
        else:
            raise ValueError(f"contrast must be 'BC', 'AB', or 'AC'; got {contrast!r}")
    elif isinstance(reconciled_or_draws, Mapping):
        source = reconciled_or_draws
    else:
        raise ValueError(
            "reconciled_or_draws must be a Nested3ArmMultiReconciliationResult or Mapping[str, np.ndarray]"
        )

    total: np.ndarray | None = None
    for outcome, w_raw in weights.items():
        w = float(w_raw)
        if not math.isfinite(w):
            raise ValueError(f"weights[{outcome!r}] must be finite; got {w_raw!r}")
        if outcome not in source:
            raise ValueError(f"outcome {outcome!r} in weights is not present in reconciled_or_draws")
        arr = np.asarray(source[outcome], dtype=np.float64)
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"draws for outcome {outcome!r} contain NaN or Inf values")
        if arr.ndim == 1:
            if horizon is not None:
                raise ValueError(
                    f"horizon={horizon!r} cannot be applied to 1-D already-cumulative draws for {outcome!r}"
                )
            cum_m = arr
        elif arr.ndim == 2:
            t_avail = arr.shape[1]
            h_eff = t_avail if horizon is None else int(horizon)
            if not (1 <= h_eff <= t_avail):
                raise ValueError(f"horizon must lie in [1, {t_avail}]; got {horizon!r}")
            sliced = arr[:, :h_eff]
            if rate == 0.0:
                cum_m = sliced.sum(axis=1)
            else:
                disc = ((1.0 + rate) ** (-1.0 / float(periods_per_year))) ** np.arange(h_eff, dtype=np.float64)
                cum_m = sliced @ disc
        else:
            raise ValueError(f"draws for {outcome!r} must be 1-D (S,) or 2-D (S, T); got shape {arr.shape}")

        contrib = w * cum_m
        if total is None:
            total = contrib
        else:
            if total.shape != contrib.shape:
                raise ValueError(f"draw count mismatch across outcomes: {total.shape} vs {contrib.shape}")
            total = total + contrib

    assert total is not None
    return total


def check_rank_k_convergence(
    result: RankKLongBetResult,
    *,
    treated_mask: np.ndarray | None = None,
    num_chains: int = 4,
    rhat_max: float = 1.05,
    ess_min: float = 100.0,
) -> RankKConvergenceCheck:
    """Evaluate split R-hat and bulk ESS convergence gates across ATT and all ``K`` modes.

    Args:
      result: Fitted :class:`RankKLongBetResult`.
      treated_mask: Optional boolean ``(N,)`` mask of treated units (defaults to ``result.treated_units``).
      num_chains: Number of chains used to split the retained draws (default ``4``).
      rhat_max: Maximum acceptable split R-hat (default ``1.05``).
      ess_min: Minimum acceptable bulk ESS (default ``100.0``).

    Returns:
      :class:`RankKConvergenceCheck`.
    """
    compute_rhat, compute_ess = _get_rhat_ess()
    n_chains_req = max(int(num_chains), 1)

    def _chains_view(trace: np.ndarray) -> np.ndarray:
        s_len = int(trace.size)
        if n_chains_req >= 2 and s_len >= 2 * n_chains_req and s_len % n_chains_req == 0:
            return trace.reshape(n_chains_req, s_len // n_chains_req)
        if s_len >= 8 and s_len % 2 == 0:
            return trace.reshape(2, s_len // 2)
        return trace.reshape(1, s_len)

    att_draws = np.asarray(result.att_draws, dtype=np.float64)
    att_chains = _chains_view(att_draws)
    rhat_att = (
        float(result.rhat_att)
        if (np.isfinite(result.rhat_att) and result.ess_att > 0.0)
        else float(compute_rhat(att_chains))
    )
    ess_att = (
        float(result.ess_att)
        if (np.isfinite(result.ess_att) and result.ess_att > 0.0)
        else float(compute_ess(att_chains))
    )

    if treated_mask is None:
        mask = np.asarray(result.treated_units, dtype=bool)
    else:
        mask = np.asarray(treated_mask, dtype=bool)
    if mask.ndim != 1 or not np.any(mask):
        mask = np.ones(result.mode_tau_draws[0].shape[0], dtype=bool)

    rhat_modes: list[float] = []
    ess_modes: list[float] = []
    failed_modes: list[int] = []
    for k_idx, mode_draws in enumerate(result.mode_tau_draws):
        arr_k = np.asarray(mode_draws, dtype=np.float64)
        trace_k = arr_k[mask].mean(axis=0)
        chains_k = _chains_view(trace_k)
        rk = float(compute_rhat(chains_k))
        ek = float(compute_ess(chains_k))
        rhat_modes.append(rk)
        ess_modes.append(ek)
        if not (math.isfinite(rk) and rk <= float(rhat_max) and math.isfinite(ek) and ek >= float(ess_min)):
            failed_modes.append(k_idx)

    att_ok = math.isfinite(rhat_att) and rhat_att <= float(rhat_max) and math.isfinite(ess_att) and ess_att >= float(ess_min)
    converged = bool(att_ok and not failed_modes)
    return RankKConvergenceCheck(
        rhat_att=rhat_att,
        ess_att=ess_att,
        rhat_modes=tuple(rhat_modes),
        ess_modes=tuple(ess_modes),
        converged=converged,
        failed_modes=tuple(failed_modes),
    )


def paid_winsor_cap(y_pre: np.ndarray, *, q: float = 0.99) -> float:
    """Compute the pre-treatment positive-part quantile cap for sparse/zero-inflated outcomes.

    Uses strictly pre-adoption outcomes ``y_pre`` only (no post-treatment outcome leakage).
    Returns the ``q``-quantile of strictly positive finite entries ``y_pre[y_pre > 0]``,
    or ``0.0`` when no strictly positive entries exist.

    Args:
      y_pre: Pre-adoption outcome array (1-D or 2-D).
      q: Quantile in ``(0, 1]`` among strictly positive pre-period values (default ``0.99``).

    Returns:
      Non-negative float cap.
    """
    q_val = float(q)
    if not (0.0 < q_val <= 1.0):
        raise ValueError(f"q must lie in (0, 1]; got {q!r}")
    arr = np.asarray(y_pre, dtype=np.float64)
    if arr.size == 0 or np.any(np.isinf(arr)):
        raise ValueError("y_pre must be a non-empty array without Inf values")
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError("y_pre contains no finite values")
    pos = finite[finite > 0.0]
    if pos.size == 0:
        return 0.0
    return float(np.quantile(pos, q_val))
