# Copyright 2026 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""In-sample trajectory stability and a pre-registerable retention rule for Rank-K analyses.

A short longitudinal experiment identifies the cumulative effect over its observed window far
better than the slope of that effect, yet a long-horizon extrapolation hinges on the slope. This
module measures how much the observed window can say about persistence -- and how much it
cannot -- without refitting anything. Following :func:`~longbet.att_stability` it reports the
numbers and withholds a verdict: :class:`TrajectoryStabilityReport` carries ``verdict = None``
by design, and no function here certifies a long-horizon extrapolation.

* :func:`wearout_ancova` is the design-based (model-free) wear-out test on the balanced
  subpanel: an ANCOVA with strata fixed effects and HC1 standard errors for the early and late
  levels and their change, with the minimum detectable difference (MDD) next to every estimate.
* :func:`estimate_post_pulse_retention` fits a per-period retention ``rho`` to every posterior draw
  of a treated-average trajectory over a post-pulse window (weighted least squares on levels,
  amplitude profiled out, grid search) and reports how many draws pile up on the ``rho`` bounds:
  positivity constraints (``catt_floor = 0``, half-normal priors, capped outcomes) push draws onto
  ``rho = 1``, and a slope fitter then reads truncation as precision.
* :func:`unit_slope_correlation` estimates the within-arm correlation across outcomes of per-unit
  post-pulse slopes, a proxy for the correlation of the outcomes' retention estimators.
* :func:`pool_multi_outcome_retention` pools ``lambda = -ln(rho)`` across outcomes
  (DerSimonian-Laird), shrinks a focal outcome toward the pooled mean (Morris 1983,
  draw-preserving) and re-evaluates the pooled ``rho`` P10 under correlation-adjusted, Q-moment,
  leave-one-outcome-out and equal-weight variants, with the break-even correlation.
* :class:`RetentionRuleSpec` and :func:`evaluate_retention_rule` turn a pre-registered retention
  rule into a machine-checked gate with a literal and a robust reading.
* :func:`diagnose_trajectory_stability` orchestrates the above into a
  :class:`TrajectoryStabilityReport`.
* :func:`pseudo_horizon_separation` sizes a leak-free hindcast before it is run: how far apart
  the no-decay (U) and decay (L) scenarios are over a pseudo-horizon, the power to tell them
  apart, and how much smaller the pseudo-horizon lever is than the real tail.

Every function is pure numpy / scipy / pandas on posterior draws and panels (no MCMC),
deterministic given its seed, raises ``ValueError`` on malformed or non-finite inputs and
returns a frozen dataclass with a ``summary_df`` and, where data are summarized, a
``data_health_certificate_df`` (row counts, NaN / Inf counts and Min / P25 / P50 / P95 / Max /
Mean). Internal invariants that cannot fail on valid inputs raise ``AssertionError``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Any
import warnings

import numpy as np
import pandas as pd
from scipy import optimize as sp_optimize
from scipy import stats as sp_stats

from longbet._aggregates import TransportAuditResult, event_time_trajectory
from longbet._rank_k import Nested3ArmReconciliationResult, RankKLongBetResult
from longbet._rank_k_plots import (
    _as_draw_matrix,
    _check_quantiles,
    _positive_scalar,
    _quantile_label,
    _summarize,
)

#: Probabilities of the retention summaries; P10 is the statistic of the pre-registered retention rule.
RETENTION_QUANTILES = (0.10, 0.50, 0.90)
#: Persistence regimes of the value layer: posterior ``rho`` (P), no decay (U) and the decay bound (L).
PERSISTENCE_REGIMES = ("P", "U", "L")
#: Regimes left in the headline when the retention rule does not pass robustly.
BOUND_REGIMES = ("U", "L")
#: Treatments of draws whose post-pulse window has no positive projection (``rho`` unidentified).
DEGENERATE_POLICIES = ("nan", "resample")
#: Per-unit slope estimators of :func:`unit_slope_correlation`.
SLOPE_METHODS = ("ols", "delta")
#: Correlation types of :func:`unit_slope_correlation`.
CORRELATION_TYPES = ("pearson", "spearman")

_HEALTH_STATS = ("Min", "P25", "P50", "P95", "Max", "Mean")
_HEALTH_COLUMNS = ("kind", "quantity", "n", "n_nan", "n_inf", *_HEALTH_STATS)
_RULE_QUANTILE = 0.10
_FISHER_Z_CLIP = 1e-12
_BLOCK_CELLS = 1 << 22  # draws x grid points per block of the grid search (32 MB of float64)
_PSD_TOL = 1e-8
_CORRELATION_VARIANTS = ("correlation_matrix", "correlation_common")
_ABOVE_ONE_TOL = 1e-9


class StabilityWarning(UserWarning):
    """Diagnostic warning of the trajectory-stability layer (pile-up, under-dispersion, fallback)."""


# ---- Shared validation and data-health helpers -------------------------------------------------


def _as_int(value: Any, name: str, minimum: int | None = None) -> int:
    """Return ``value`` as an ``int`` (integral floats allowed), raising ``ValueError`` otherwise."""
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer; got {value!r}")
    try:
        out = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be an integer; got {value!r}") from exc
    if out != value:
        raise ValueError(f"{name} must be an integer; got {value!r}")
    if minimum is not None and out < minimum:
        raise ValueError(f"{name} must be >= {minimum}; got {value!r}")
    return out


def _as_float(value: Any, name: str) -> float:
    """Return ``value`` as a finite float, raising ``ValueError`` otherwise."""
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number; got {value!r}") from exc
    if not math.isfinite(out):
        raise ValueError(f"{name} must be a finite number; got {value!r}")
    return out


def _check_open_unit(value: Any, name: str) -> float:
    """Return ``value`` as a float strictly inside ``(0, 1)``."""
    out = _as_float(value, name)
    if not 0.0 < out < 1.0:
        raise ValueError(f"{name} must lie strictly inside (0, 1); got {value!r}")
    return out


def _check_window(window: Any, n_periods: int, name: str, min_length: int = 1) -> tuple[int, int]:
    """Validate a 1-based inclusive ``(first, last)`` window inside ``1..n_periods``."""
    try:
        first_raw, last_raw = window
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a (first, last) pair of 1-based periods; got {window!r}") from exc
    first = _as_int(first_raw, f"{name}[0]", 1)
    last = _as_int(last_raw, f"{name}[1]", 1)
    if last - first + 1 < min_length or last > n_periods:
        raise ValueError(
            f"{name} must satisfy 1 <= first <= last <= {n_periods} and span at least {min_length} "
            f"period(s); got {window!r}"
        )
    return first, last


def _as_vector(value: Any, name: str, length: int | None = None, *, allow_nan: bool = False) -> np.ndarray:
    """Coerce ``value`` to a non-empty 1-D float64 array, raising ``ValueError`` on bad shapes or entries."""
    try:
        arr = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric 1-D array; got {type(value).__name__}") from exc
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a non-empty 1-D array; got shape {arr.shape}")
    if length is not None and arr.shape[0] != length:
        raise ValueError(f"{name} must have length {length}; got {arr.shape[0]}")
    bad = np.isinf(arr) if allow_nan else ~np.isfinite(arr)
    if np.any(bad):
        kind = "infinite" if allow_nan else "NaN / Inf"
        raise ValueError(f"{name} contains {int(np.sum(bad))} {kind} entries")
    return arr


def _as_labels(value: Any, length: int, name: str) -> np.ndarray:
    """Validate a 1-D array of hashable labels with one entry per unit and no missing labels."""
    labels = np.asarray(value)
    if labels.ndim != 1 or labels.shape[0] != length:
        raise ValueError(f"{name} must be 1-D with one label per unit ({length}); got shape {labels.shape}")
    if labels.dtype.kind == "f" and np.any(np.isnan(labels)):
        raise ValueError(f"{name} contains {int(np.sum(np.isnan(labels)))} NaN labels")
    if labels.dtype.kind == "O":
        n_missing = sum(1 for v in labels if v is None or (isinstance(v, float) and math.isnan(v)))
        if n_missing:
            raise ValueError(f"{name} contains {n_missing} missing (None / NaN) labels")
    return labels


def _mean_off_diagonal(matrix: np.ndarray) -> float:
    """Mean of the strictly upper-triangular entries of a square matrix."""
    return float(np.mean(matrix[np.triu_indices(matrix.shape[0], k=1)]))


def _summary_or_nan(values: Any, quantiles: Sequence[float], prefix: str) -> dict[str, float]:
    """:func:`_summarize` over the finite entries, or NaN summaries when there are none."""
    flat = np.asarray(values, dtype=np.float64).ravel()
    finite = flat[np.isfinite(flat)]
    if finite.size:
        return _summarize(finite, quantiles, prefix)
    keys = [f"{prefix}{_quantile_label(q)}" for q in quantiles] + [f"{prefix}Mean", f"{prefix}SD"]
    return dict.fromkeys(keys, float("nan"))


def _distribution_row(quantity: str, values: Any) -> dict[str, Any]:
    """Data-health row: NaN / Inf counts and Min / P25 / P50 / P95 / Max / Mean of the finite values."""
    flat = np.asarray(values, dtype=np.float64).ravel()
    finite = flat[np.isfinite(flat)]
    row: dict[str, Any] = {
        "kind": "distribution",
        "quantity": quantity,
        "n": int(flat.size),
        "n_nan": int(np.sum(np.isnan(flat))),
        "n_inf": int(np.sum(np.isinf(flat))),
    }
    if finite.size:
        row.update(
            Min=float(np.min(finite)),
            P25=float(np.quantile(finite, 0.25)),
            P50=float(np.quantile(finite, 0.50)),
            P95=float(np.quantile(finite, 0.95)),
            Max=float(np.max(finite)),
            Mean=float(np.mean(finite)),
        )
    else:
        row.update(dict.fromkeys(_HEALTH_STATS, float("nan")))
    return row


def _count_row(quantity: str, count: int) -> dict[str, Any]:
    """Row-count audit row of a data-health certificate (the distribution columns are NaN)."""
    row: dict[str, Any] = {"kind": "row_count", "quantity": quantity, "n": int(count), "n_nan": 0, "n_inf": 0}
    row.update(dict.fromkeys(_HEALTH_STATS, float("nan")))
    return row


def _health_frame(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    """Data-health certificate with the fixed columns ``kind, quantity, n, n_nan, n_inf, Min..Mean``."""
    return pd.DataFrame([dict(row) for row in rows], columns=list(_HEALTH_COLUMNS))


# ---- Design-based wear-out test ------------------------------------------------------------


@dataclass(frozen=True, eq=False)
class WearoutAncovaResult:
    """Output of :func:`wearout_ancova`: model-free wear-out test on the balanced subpanel.

    Attributes:
      early: ``(first, last)`` event times of the early level.
      late: ``(first, last)`` event times of the late level.
      late_split: The two late sub-windows of the late-split change, or ``None``.
      min_post_periods: Minimum observed window ``L_i`` of the balanced subpanel.
      alpha: Two-sided level of the test and of the MDD.
      power: Power of the MDD.
      mdd_multiplier: ``z_{1 - alpha / 2} + z_power`` (2.80 at the defaults).
      n_units: Units in the balanced subpanel.
      n_treated: Treated units in the balanced subpanel.
      n_control: Control units in the balanced subpanel.
      early_estimate: ANCOVA estimate of the early-level effect.
      change_estimate: ANCOVA estimate of the change ``mean(late) - mean(early)``.
      change_mdd: MDD of that change.
      mdd_pct_of_early: ``100 * change_mdd / |early_estimate|`` (NaN when the early estimate is 0).
      powered_for_half_wearout: Whether ``change_mdd <= 0.5 * max(early_estimate, 0)`` with a
        positive standard error; when it is ``False`` a null change is not evidence of persistence.
      summary_df: One row per measure (``level_early``, ``level_late``, ``change_late_minus_early``
        and, with ``late_split``, ``change_late_split``) with ``measure``, ``label``, ``n``,
        ``n_trt``, ``n_ctrl``, ``est``, ``se_hc1``, ``t_stat``, ``p_value`` (two-sided, normal),
        ``mdd``, ``mdd_pct_of_early`` and ``powered_for_half_wearout`` (change rows only).
      data_health_certificate_df: Row counts (input, balanced, by arm and stratum) and the
        distribution of the outcome window, the baseline and every measure.
      warnings: Warnings emitted by the call (for instance when the treated share of the balanced
        subpanel differs from the input panel).
    """

    early: tuple[int, int]
    late: tuple[int, int]
    late_split: tuple[tuple[int, int], tuple[int, int]] | None
    min_post_periods: int
    alpha: float
    power: float
    mdd_multiplier: float
    n_units: int
    n_treated: int
    n_control: int
    early_estimate: float
    change_estimate: float
    change_mdd: float
    mdd_pct_of_early: float
    powered_for_half_wearout: bool
    summary_df: pd.DataFrame
    data_health_certificate_df: pd.DataFrame
    warnings: tuple[str, ...] = ()


def _as_treated(value: Any, length: int) -> np.ndarray:
    """Boolean ``(N,)`` treatment indicator from booleans or 0 / 1 values."""
    arr = np.asarray(value)
    if arr.ndim != 1 or arr.shape[0] != length:
        raise ValueError(f"treated must be 1-D with one entry per unit ({length}); got shape {arr.shape}")
    if arr.dtype == bool:
        return arr.copy()
    try:
        numeric = arr.astype(np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("treated must hold booleans or 0 / 1 values") from exc
    if not np.all(np.isin(numeric, (0.0, 1.0))):
        raise ValueError("treated must hold booleans or 0 / 1 values")
    return numeric == 1.0


def _as_covariates(value: Any, length: int) -> np.ndarray | None:
    """``(N, P)`` float covariates (a 1-D input is one column) without infinite entries, or ``None``."""
    if value is None:
        return None
    try:
        covariates = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"x must be a numeric (N, P) array; got {type(value).__name__}") from exc
    covariates = covariates[:, None] if covariates.ndim == 1 else covariates
    if covariates.ndim != 2 or covariates.shape[0] != length or np.any(np.isinf(covariates)):
        raise ValueError(f"x must be an (N, P) array with N = {length} and no infinite entries; got {covariates.shape}")
    return covariates


def _finite_prefix_lengths(values: np.ndarray) -> np.ndarray:
    """Length of the finite prefix of every row (observed periods come first under administrative censoring)."""
    finite = np.isfinite(values)
    return np.where(np.all(finite, axis=1), values.shape[1], np.argmin(finite, axis=1)).astype(np.float64)


def _strata_dummies(strata: np.ndarray) -> list[np.ndarray]:
    """Fixed-effect dummies of every stratum present in the sample but the first (sorted) one."""
    try:
        levels = list(np.unique(strata))
    except TypeError:
        levels = sorted(set(strata.tolist()), key=repr)
    columns = []
    for level in levels[1:]:
        column = (strata == level).astype(np.float64)
        if np.std(column) > 1e-8:
            columns.append(column)
    return columns


def _ancova_hc1(outcome: np.ndarray, design: np.ndarray) -> tuple[float, float]:
    """OLS coefficient of design column 1 (the treatment) and its HC1 standard error."""
    n_obs, n_params = design.shape
    beta, _, _, _ = np.linalg.lstsq(design, outcome, rcond=None)
    resid = outcome - design @ beta
    xtx_inv = np.linalg.inv(design.T @ design)
    meat = (design * (resid[:, None] ** 2)).T @ design
    cov_hc1 = (n_obs / float(n_obs - n_params)) * (xtx_inv @ meat @ xtx_inv)
    return float(beta[1]), float(np.sqrt(cov_hc1[1, 1]))


def _late_split_windows(late_split: Any, n_periods: int) -> tuple[tuple[int, int], tuple[int, int]] | None:
    """Validate the two late sub-windows of the late-split change (``None`` skips that measure)."""
    if late_split is None:
        return None
    try:
        first_half, second_half = late_split
    except (TypeError, ValueError) as exc:
        raise ValueError(f"late_split must be a pair of (first, last) windows or None; got {late_split!r}") from exc
    return (
        _check_window(first_half, n_periods, "late_split[0]"),
        _check_window(second_half, n_periods, "late_split[1]"),
    )


def _wearout_design(
    is_treated: np.ndarray, baseline: np.ndarray, strata: np.ndarray | None, covariates: np.ndarray | None
) -> np.ndarray:
    """ANCOVA design ``[1, Z, b - mean(b), strata dummies, X]`` of the balanced subpanel."""
    columns = [np.ones(is_treated.shape[0]), is_treated.astype(np.float64), baseline - float(np.mean(baseline))]
    if strata is not None:
        columns.extend(_strata_dummies(strata))
    if covariates is not None:
        columns.extend(covariates.T)
    design = np.column_stack(columns)
    if design.shape[0] <= design.shape[1]:
        raise ValueError(
            f"the balanced subpanel has {design.shape[0]} units for {design.shape[1]} ANCOVA parameters; "
            "lower min_post_periods or drop covariates"
        )
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("the ANCOVA design is rank deficient (collinear baseline, strata or covariates)")
    return design


def _wearout_measures(
    y_balanced: np.ndarray,
    early: tuple[int, int],
    late: tuple[int, int],
    late_split: tuple[tuple[int, int], tuple[int, int]] | None,
) -> list[tuple[str, str, np.ndarray]]:
    """``[(key, label, (n,) values)]`` of the early / late levels and the late-minus-early and late-split changes."""
    means = {}
    for window in (early, late) + (() if late_split is None else late_split):
        means[window] = y_balanced[:, window[0] - 1 : window[1]].mean(axis=1)
    label = {window: f"periods {window[0]}-{window[1]}" for window in means}
    measures = [
        ("level_early", f"Level, {label[early]}", means[early]),
        ("level_late", f"Level, {label[late]}", means[late]),
        ("change_late_minus_early", f"Change: mean({label[late]}) - mean({label[early]})", means[late] - means[early]),
    ]
    if late_split is not None:
        first, second = late_split
        measures.append(
            ("change_late_split", f"Change: mean({label[second]}) - mean({label[first]})", means[second] - means[first])
        )
    for key, measure_label, values in measures:
        if float(np.ptp(values)) == 0.0:
            if key == "change_late_split":
                remedy = "pass late_split=None or other late_split windows"
            else:
                remedy = "check the outcome and change the early / late windows"
            raise ValueError(
                f"the measure {measure_label!r} is constant on the balanced subpanel, so its ANCOVA is degenerate; "
                f"{remedy}"
            )
    return measures


def _wearout_rows(
    fits: Mapping[str, tuple[str, float, float]], multiplier: float, n_units: int, n_treated: int
) -> list[dict[str, Any]]:
    """One summary row per measure; power is claimed only with a positive standard error."""
    early_est = fits["level_early"][1]
    rows = []
    for key, (label, est, se) in fits.items():
        mdd = multiplier * se
        t_stat = est / se if se > 0.0 else float("nan")
        powered = None
        if key.startswith("change"):
            powered = bool(se > 0.0 and mdd <= 0.5 * max(early_est, 0.0))
        rows.append({
            "measure": key,
            "label": label,
            "n": n_units,
            "n_trt": n_treated,
            "n_ctrl": n_units - n_treated,
            "est": est,
            "se_hc1": se,
            "t_stat": t_stat,
            "p_value": float(2.0 * sp_stats.norm.sf(abs(t_stat))) if math.isfinite(t_stat) else float("nan"),
            "mdd": mdd,
            "mdd_pct_of_early": 100.0 * mdd / abs(early_est) if early_est != 0.0 else float("nan"),
            "powered_for_half_wearout": powered,
        })
    return rows


def _attrition_messages(
    n_units: int, n_trt_in: int, n_bal: int, n_trt_bal: int, min_periods: int, alpha: float
) -> list[str]:
    """Warn when selection into the balanced subpanel (``L_i >= min_periods``) differs by treatment arm."""
    n_drop = n_units - n_bal
    if n_drop <= 0:
        return []
    share_in = n_trt_in / float(n_units)
    share_bal = n_trt_bal / float(n_bal)
    share_drop = (n_trt_in - n_trt_bal) / float(n_drop)
    pooled_var = share_in * (1.0 - share_in) * (1.0 / n_bal + 1.0 / n_drop)
    if pooled_var <= 0.0:
        return []
    z_stat = abs(share_bal - share_drop) / math.sqrt(pooled_var)
    p_val = float(2.0 * sp_stats.norm.sf(z_stat))
    if p_val >= alpha:
        return []
    return [
        f"the treated share in the balanced subpanel (L_i >= {min_periods}: {n_trt_bal}/{n_bal} = "
        f"{share_bal:.1%}) differs from the input panel ({n_trt_in}/{n_units} = {share_in:.1%}; dropped units "
        f"{share_drop:.1%}, p = {p_val:.3g} < {alpha:g}): L_i must reflect administrative censoring "
        "independent of treatment, not post-treatment attrition."
    ]


def wearout_ancova(
    y_event: Any,
    treated: Any,
    baseline: Any,
    *,
    strata: Any = None,
    x: Any = None,
    post_window_length: Any = None,
    min_post_periods: int | None = None,
    early: tuple[int, int] = (1, 6),
    late: tuple[int, int] = (7, 12),
    late_split: tuple[tuple[int, int], tuple[int, int]] | None = ((7, 9), (10, 12)),
    alpha: float = 0.05,
    power: float = 0.80,
) -> WearoutAncovaResult:
    """Design-based wear-out test: ANCOVA of early and late levels and their change on the balanced subpanel.

    The model-free counterpart of the slope a decay extrapolation relies on. The balanced subpanel
    holds the units with ``L_i >= min_post_periods``, so every measure averages the same observed
    periods for every unit. For each measure ``Y_i`` (early level, late level, the change
    ``mean(late) - mean(early)`` and the late-split change
    ``mean(late_split[1]) - mean(late_split[0])``)::

        Y_i = alpha_{g(i)} + beta Z_i + gamma (b_i - b_bar) + delta' X_i + e_i,

    with stratum fixed effects ``alpha_g``, the baseline centred on the subpanel mean, HC1
    covariance, a two-sided normal p-value and the minimum detectable difference
    ``MDD = (z_{1 - alpha / 2} + z_power) SE`` (2.80 SE at the defaults). A level panel can be well
    powered for the cumulative effect and still have an MDD of the change larger than the early
    effect itself: then no detectable wear-out is not evidence of persistence.

    Args:
      y_event: ``(N, H)`` outcome by event time ``s = 1..H`` (NaN = unobserved).
      treated: ``(N,)`` treatment indicator (booleans or 0 / 1).
      baseline: ``(N,)`` pre-period rate (centred inside on the balanced subpanel); it must be
        finite on the balanced subpanel.
      strata: Optional ``(N,)`` stratum labels (for instance the randomization block or the
        enrolment wave) entering as fixed effects.
      x: Optional ``(N, P)`` (or ``(N,)``) covariates, finite on the balanced subpanel.
      post_window_length: Optional ``(N,)`` observed window lengths ``L_i``; default: the length of
        each row's finite prefix of ``y_event``. ``L_i`` must reflect administrative censoring
        independent of treatment (not post-treatment attrition); a :class:`StabilityWarning` is
        emitted when the treated share of the balanced subpanel differs from the input panel.
      min_post_periods: Minimum ``L_i`` of the balanced subpanel; default (and minimum allowed) the
        last event time any measure uses.
      early: ``(first, last)`` event times of the early level.
      late: ``(first, last)`` event times of the late level.
      late_split: Two ``(first, last)`` sub-windows of the late period for the late-split change, or
        ``None`` to skip that measure.
      alpha: Two-sided level of the test and the MDD, in ``(0, 1)``.
      power: Power of the MDD, in ``(0, 1)``.

    Returns:
      :class:`WearoutAncovaResult`.

    Raises:
      ValueError: On malformed inputs, non-finite outcomes, baseline or covariates inside the
        balanced subpanel, fewer than 2 treated or 2 control units, a measure that is constant on
        the balanced subpanel, or a rank-deficient design.
    """
    try:
        y = np.asarray(y_event, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"y_event must be a numeric (N, H) array; got {type(y_event).__name__}") from exc
    if y.ndim != 2 or y.shape[0] < 4 or y.shape[1] < 1:
        raise ValueError(f"y_event must be an (N, H) array with N >= 4 units; got shape {y.shape}")
    if np.any(np.isinf(y)):
        raise ValueError(f"y_event contains {int(np.sum(np.isinf(y)))} infinite entries")
    n_units, n_periods = y.shape
    is_treated = _as_treated(treated, n_units)
    base = _as_vector(baseline, "baseline", n_units, allow_nan=True)
    strata_arr = None if strata is None else _as_labels(strata, n_units, "strata")
    covariates = _as_covariates(x, n_units)
    early_win = _check_window(early, n_periods, "early")
    late_win = _check_window(late, n_periods, "late")
    split = _late_split_windows(late_split, n_periods)
    last_needed = max([early_win[1], late_win[1]] + ([] if split is None else [split[0][1], split[1][1]]))
    min_periods = last_needed if min_post_periods is None else _as_int(min_post_periods, "min_post_periods", 1)
    if not last_needed <= min_periods <= n_periods:
        raise ValueError(
            f"min_post_periods must lie in [{last_needed}, {n_periods}] (the measures use event times up "
            f"to {last_needed}); got {min_post_periods!r}"
        )
    alpha_val = _check_open_unit(alpha, "alpha")
    power_val = _check_open_unit(power, "power")
    if post_window_length is None:
        lengths = _finite_prefix_lengths(y)
    else:
        lengths = _as_vector(post_window_length, "post_window_length", n_units)
    balanced = lengths >= min_periods
    y_bal = y[balanced, :last_needed]
    for name, values in (("outcomes", y_bal), ("baseline", base[balanced]),
                         ("covariates", None if covariates is None else covariates[balanced])):
        if values is not None and not np.all(np.isfinite(values)):
            raise ValueError(
                f"{int(np.sum(~np.isfinite(values)))} {name} of the balanced subpanel (event times 1.."
                f"{last_needed}) are NaN; check post_window_length, the baseline and x"
            )
    trt = is_treated[balanced]
    n_bal, n_trt = int(balanced.sum()), int(trt.sum())
    if n_trt < 2 or n_bal - n_trt < 2:
        raise ValueError(
            f"the balanced subpanel (L_i >= {min_periods}) has {n_trt} treated and {n_bal - n_trt} control "
            "units; at least 2 of each are needed"
        )
    messages = _attrition_messages(n_units, int(is_treated.sum()), n_bal, n_trt, min_periods, alpha_val)
    for message in messages:
        warnings.warn(message, StabilityWarning, stacklevel=2)
    design = _wearout_design(
        trt,
        base[balanced],
        None if strata_arr is None else strata_arr[balanced],
        None if covariates is None else covariates[balanced],
    )
    measures = _wearout_measures(y_bal, early_win, late_win, split)
    multiplier = float(sp_stats.norm.ppf(1.0 - alpha_val / 2.0) + sp_stats.norm.ppf(power_val))
    fits = {key: (label, *_ancova_hc1(values, design)) for key, label, values in measures}
    summary = pd.DataFrame(_wearout_rows(fits, multiplier, n_bal, n_trt))
    change = summary.loc[summary["measure"] == "change_late_minus_early"].iloc[0]

    health = [
        _count_row("units (input)", n_units),
        _count_row("treated units (input)", int(is_treated.sum())),
        _count_row(f"units with L_i >= {min_periods} (balanced subpanel)", n_bal),
        _count_row("treated units (balanced)", n_trt),
        _count_row("control units (balanced)", n_bal - n_trt),
    ]
    if strata_arr is not None:
        labels, counts = np.unique(strata_arr[balanced].astype(str), return_counts=True)
        for stratum, count in zip(labels, counts, strict=True):
            health.append(_count_row(f"stratum {stratum} (balanced)", count))
    health.append(_distribution_row(f"y_event periods 1..{last_needed} (balanced)", y_bal))
    health.append(_distribution_row("baseline (balanced)", base[balanced]))
    health.extend(_distribution_row(f"measure: {label}", values) for _, label, values in measures)
    return WearoutAncovaResult(
        early=early_win,
        late=late_win,
        late_split=split,
        min_post_periods=min_periods,
        alpha=alpha_val,
        power=power_val,
        mdd_multiplier=multiplier,
        n_units=n_bal,
        n_treated=n_trt,
        n_control=n_bal - n_trt,
        early_estimate=float(fits["level_early"][1]),
        change_estimate=float(change["est"]),
        change_mdd=float(change["mdd"]),
        mdd_pct_of_early=float(change["mdd_pct_of_early"]),
        powered_for_half_wearout=bool(change["powered_for_half_wearout"]),
        summary_df=summary,
        data_health_certificate_df=_health_frame(health),
        warnings=tuple(messages),
    )


# ---- Post-pulse retention per posterior draw -----------------------------------------------


@dataclass(frozen=True, eq=False)
class PostPulseRetentionResult:
    """Output of :func:`estimate_post_pulse_retention`: per-draw retention per period of one trajectory.

    Attributes:
      label: Outcome label used by the pooling and report tables.
      window: ``(s0, s1)`` 1-based event times of the post-pulse window (inclusive).
      weights: ``(s1 - s0 + 1,)`` support weights of the window periods used by the fit.
      support_source: ``"user"``, ``"post_window_length"`` (units whose observed window reaches each
        event time) or ``"uniform"``.
      trajectory_source: Where the trajectory draws came from.
      rho_bounds: ``(rho_lo, rho_hi)`` of the constrained grid.
      unconstrained_upper: Upper end of the unconstrained diagnostic grid.
      grid_step: Spacing of both grids.
      degenerate_policy: ``"nan"`` or ``"resample"``.
      seed: Seed of the ``"resample"`` policy.
      rho_draws: ``(S,)`` constrained grid fits (NaN, or resampled from the identified draws, for
        degenerate draws).
      rho_unconstrained_draws: ``(S,)`` fits on ``[rho_lo, unconstrained_upper]`` (NaN for degenerate
        draws under either policy).
      lambda_draws: ``(S,)`` ``-ln(rho_draws)``.
      degenerate: ``(S,)`` draws whose window has no positive projection on the grid (for instance
        an identically zero or an entirely non-positive window), so that the best amplitude is 0
        for every ``rho`` and ``rho`` is unidentified.
      share_at_upper: Share of all draws whose identified constrained fit sits on ``rho_hi``.
      share_at_lower: Share of all draws whose identified constrained fit sits on ``rho_lo``.
      share_degenerate: Share of degenerate draws.
      share_unconstrained_above_one: Share of all draws whose identified unconstrained fit exceeds 1.
      pileup_threshold: Threshold of ``boundary_pileup_flag``.
      boundary_pileup_flag: ``share_at_upper >= pileup_threshold``.
      summary_df: One row with the window, draw counts, quantiles / mean / SD (finite draws) of
        ``rho``, the unconstrained ``rho`` and ``lambda``, the boundary shares and the flag.
      support_df: One row per event time of the trajectory with ``event_time``, ``support`` (NaN
        where unknown), ``in_window``, ``weight`` and ``weight_share`` (NaN outside the window).
      data_health_certificate_df: Row counts and the distribution of every window period and fit.
    """

    label: str
    window: tuple[int, int]
    weights: np.ndarray
    support_source: str
    trajectory_source: str
    rho_bounds: tuple[float, float]
    unconstrained_upper: float
    grid_step: float
    degenerate_policy: str
    seed: int
    rho_draws: np.ndarray
    rho_unconstrained_draws: np.ndarray
    lambda_draws: np.ndarray
    degenerate: np.ndarray
    share_at_upper: float
    share_at_lower: float
    share_degenerate: float
    share_unconstrained_above_one: float
    pileup_threshold: float
    boundary_pileup_flag: bool
    summary_df: pd.DataFrame
    support_df: pd.DataFrame
    data_health_certificate_df: pd.DataFrame


def _support_counts(fit: RankKLongBetResult, n_periods: int) -> np.ndarray | None:
    """Units (treated and control) whose observed window reaches event time ``s = 1..n_periods``."""
    if fit.post_window_length is None:
        return None
    lengths = np.asarray(fit.post_window_length, dtype=np.float64).reshape(-1)
    periods = np.arange(1, n_periods + 1, dtype=np.float64)
    return np.sum(lengths[None, :] >= periods[:, None], axis=1).astype(np.float64)


def _resolve_trajectory(
    traj: Any, *, contrast: str, weighting: str, name: str
) -> tuple[np.ndarray, np.ndarray | None, str]:
    """``(draws (S, L), support counts (L,) or None, source)`` from draws, a mapping, a fit or a reconciliation."""
    if isinstance(traj, RankKLongBetResult):
        draws = _as_draw_matrix(event_time_trajectory(traj, weighting=weighting), f"event_time_trajectory({name})")
        return draws, _support_counts(traj, draws.shape[1]), f"event_time_trajectory(weighting={weighting!r})"
    if isinstance(traj, Nested3ArmReconciliationResult):
        by_contrast = {"BC": traj.trt_catt_bc, "AB": traj.trt_catt_ab, "AC": traj.trt_catt_ac}
        if contrast not in by_contrast:
            raise ValueError(
                f"contrast must be 'BC', 'AB' or 'AC' for a Nested3ArmReconciliationResult; got {contrast!r}"
            )
        source = f"reconciliation.trt_catt_{contrast.lower()}"
        return _as_draw_matrix(by_contrast[contrast], f"{name}.trt_catt_{contrast.lower()}"), None, source
    if isinstance(traj, Mapping):
        if contrast not in traj:
            raise ValueError(f"{name} is a mapping without the contrast key {contrast!r}; got {sorted(map(str, traj))}")
        return _as_draw_matrix(traj[contrast], f"{name}[{contrast!r}]"), None, f"mapping[{contrast!r}]"
    return _as_draw_matrix(traj, name), None, "draws"


def _resolve_support(
    support: Any, counts: np.ndarray | None, n_periods: int, window: tuple[int, int]
) -> tuple[np.ndarray, np.ndarray, str]:
    """``(window weights, full-length support with NaN where unknown, source)`` of the retention fit."""
    first, last = window
    n_window = last - first + 1
    if support is not None:
        values = _as_vector(support, "support")
        if values.shape[0] == n_periods:
            full = values.copy()
        elif values.shape[0] == n_window:
            full = np.full(n_periods, np.nan)
            full[first - 1 : last] = values
        else:
            raise ValueError(
                f"support must hold one weight per event time ({n_periods}) or per window period ({n_window}); "
                f"got {values.shape[0]}"
            )
        source = "user"
    elif counts is not None:
        full = np.asarray(counts, dtype=np.float64).copy()
        source = "post_window_length"
    else:
        full = np.full(n_periods, np.nan)
        full[first - 1 : last] = 1.0
        source = "uniform"
    weights = full[first - 1 : last].copy()
    if np.any(weights < 0.0) or not np.sum(weights) > 0.0:
        raise ValueError(f"support weights over the window must be >= 0 with a positive sum; got {weights.tolist()}")
    return weights, full, source


def _profiled_retention_fit(
    tau_window: np.ndarray, weights: np.ndarray, grid: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Per-draw grid maximizer of the profiled gain of ``tau(m) = a rho^(m - s0)`` with ``a >= 0``.

    With ``r_m = rho^(m - s0)`` the weighted SSE profiled over ``a >= 0`` is
    ``sum w tau^2 - max(sum w tau r, 0)^2 / sum w r^2``, so the fit maximizes the gain
    ``g(rho) = max(sum w tau r, 0)^2 / sum w r^2`` on the grid, vectorized over blocks of draws.
    A draw whose gain is zero on the whole grid has no positive projection: ``a* = 0`` for every
    ``rho`` and ``rho`` is unidentified.

    Returns:
      ``(rho, identified)``: ``(S,)`` maximizers (the first grid point when unidentified) and flags.
    """
    exponents = np.arange(tau_window.shape[1], dtype=np.float64)
    basis = grid[:, None] ** exponents[None, :]
    norm = (basis**2) @ weights
    n_draws = tau_window.shape[0]
    rho = np.empty(n_draws, dtype=np.float64)
    identified = np.empty(n_draws, dtype=bool)
    block = max(1, _BLOCK_CELLS // grid.size)
    for start in range(0, n_draws, block):
        stop = min(start + block, n_draws)
        cross = (tau_window[start:stop] * weights[None, :]) @ basis.T
        gain = np.where(cross > 0.0, cross**2 / norm[None, :], 0.0)
        rho[start:stop] = grid[np.argmax(gain, axis=1)]
        identified[start:stop] = np.any(gain > 0.0, axis=1)
    return rho, identified


def _check_rho_bounds(rho_bounds: Any, name: str = "rho_bounds") -> tuple[float, float]:
    """Validate ``(rho_lo, rho_hi)`` with ``0 < rho_lo < rho_hi``."""
    try:
        low_raw, high_raw = rho_bounds
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a (low, high) pair; got {rho_bounds!r}") from exc
    low, high = _as_float(low_raw, f"{name}[0]"), _as_float(high_raw, f"{name}[1]")
    if not 0.0 < low < high:
        raise ValueError(f"{name} must satisfy 0 < low < high; got {rho_bounds!r}")
    return low, high


def estimate_post_pulse_retention(
    traj: Any,
    *,
    window: tuple[int, int],
    support: Any = None,
    contrast: str = "AC",
    weighting: str = "unit",
    rho_bounds: tuple[float, float] = (0.50, 1.00),
    unconstrained_upper: float = 1.50,
    grid_points: int = 1001,
    degenerate_policy: str = "nan",
    pileup_threshold: float = 0.25,
    seed: int = 2026,
    quantiles: Sequence[float] = RETENTION_QUANTILES,
    label: str = "effect",
) -> PostPulseRetentionResult:
    """Per-draw per-period retention ``rho`` of a treated-average trajectory over a post-pulse window.

    For every posterior draw ``d`` the estimator is weighted least squares on levels::

        min_{a >= 0, rho in [rho_lo, rho_hi]}  sum_{m = s0..s1} w_m (tau_m^(d) - a rho^(m - s0))^2,

    with the amplitude profiled out in closed form and ``rho`` maximized on a grid. Levels, not
    logs, so exact zeros from ``catt_floor = 0`` are admissible. An unconstrained fit on
    ``[rho_lo, unconstrained_upper]`` measures how much of the mass at ``rho_hi`` is boundary
    pile-up rather than evidence: positivity constraints and capped outcomes push draws onto
    ``rho = 1``, which a slope fitter reads as precision.

    A draw whose window has no positive projection on the grid (``sum_m w_m tau_m rho^(m - s0) <= 0``
    for every grid ``rho``; for instance an identically zero or an entirely non-positive window)
    has a zero gain for every ``rho``, so the best amplitude is 0 and ``rho`` is unidentified; such
    draws are *degenerate*. The rule is a property of the gain, hence invariant to the scale of the
    outcome. ``degenerate_policy="nan"`` returns NaN for them (and the summaries ignore them);
    ``"resample"`` replaces each by a draw sampled with replacement from the identified draws.
    Either way the share is reported.

    Args:
      traj: ``(S, L)`` trajectory draws; a :class:`~longbet.RankKLongBetResult` (uses
        :func:`~longbet.event_time_trajectory` with ``weighting`` and, as support, the number of
        units whose ``post_window_length`` reaches each event time); a
        :class:`~longbet.Nested3ArmReconciliationResult` or a ``{contrast: (S, L)}`` mapping (uses
        ``contrast``).
      window: ``(s0, s1)`` 1-based inclusive post-pulse event times (at least two periods).
      support: Optional non-negative weights, one per event time ``(L,)`` or per window period;
        default: the fit's support counts, else uniform.
      contrast: Contrast of a reconciliation or mapping input (``"AC"`` by default).
      weighting: ``"unit"`` (population trajectory, default) or ``"cell"`` for a fit input.
      rho_bounds: ``(rho_lo, rho_hi)`` of the constrained grid, ``0 < rho_lo < rho_hi``.
      unconstrained_upper: Upper end (``>= rho_hi``) of the unconstrained diagnostic grid, which
        keeps the constrained grid's step.
      grid_points: Points of the constrained grid (``>= 3``; 1001 gives a 0.0005 step on
        ``[0.50, 1.00]``).
      degenerate_policy: ``"nan"`` (default) or ``"resample"``.
      pileup_threshold: ``boundary_pileup_flag`` fires when ``share_at_upper >= pileup_threshold``.
      seed: Seed of the ``"resample"`` policy.
      quantiles: Probabilities of the summary (odd length; default P10 / P50 / P90).
      label: Outcome label carried into the summary tables.

    Returns:
      :class:`PostPulseRetentionResult`.

    Raises:
      ValueError: On malformed or non-finite inputs, or when ``"resample"`` finds no identified draw.
    """
    draws, counts, trajectory_source = _resolve_trajectory(traj, contrast=contrast, weighting=weighting, name="traj")
    n_draws, n_periods = draws.shape
    first, last = _check_window(window, n_periods, "window", min_length=2)
    qs = _check_quantiles(quantiles)
    low, high = _check_rho_bounds(rho_bounds)
    upper = _as_float(unconstrained_upper, "unconstrained_upper")
    if upper < high:
        raise ValueError(f"unconstrained_upper must be >= rho_bounds[1] = {high}; got {unconstrained_upper!r}")
    n_grid = _as_int(grid_points, "grid_points", 3)
    if degenerate_policy not in DEGENERATE_POLICIES:
        raise ValueError(f"degenerate_policy must be one of {DEGENERATE_POLICIES}; got {degenerate_policy!r}")
    threshold = _as_float(pileup_threshold, "pileup_threshold")
    if not 0.0 < threshold <= 1.0:
        raise ValueError(f"pileup_threshold must lie in (0, 1]; got {pileup_threshold!r}")
    seed_value = _as_int(seed, "seed", 0)
    weights, support_full, support_source = _resolve_support(support, counts, n_periods, (first, last))

    tau_window = draws[:, first - 1 : last]
    step = (high - low) / (n_grid - 1)
    grid = np.linspace(low, high, n_grid)
    n_unconstrained = max(n_grid, int(math.ceil((upper - low) / step - 1e-9)) + 1)
    grid_unconstrained = np.linspace(low, low + step * (n_unconstrained - 1), n_unconstrained)
    rho_raw, identified = _profiled_retention_fit(tau_window, weights, grid)
    rho_unconstrained_raw, _ = _profiled_retention_fit(tau_window, weights, grid_unconstrained)
    degenerate = ~identified
    n_degenerate = int(degenerate.sum())
    if degenerate_policy == "resample" and n_degenerate:
        if n_degenerate == n_draws:
            raise ValueError("every draw's window is degenerate (rho unidentified); there is nothing to resample from")
        rng = np.random.default_rng(seed_value)
        rho = rho_raw.copy()
        rho[degenerate] = rng.choice(rho_raw[identified], size=n_degenerate, replace=True)
    elif degenerate_policy == "resample":
        rho = rho_raw.copy()
    else:
        rho = np.where(identified, rho_raw, np.nan)
    rho_unconstrained = np.where(identified, rho_unconstrained_raw, np.nan)
    lam = 0.0 - np.log(rho)
    fitted = rho[np.isfinite(rho)]
    if fitted.size and (float(np.min(fitted)) < low or float(np.max(fitted)) > high):
        raise AssertionError("constrained rho draws fell outside rho_bounds")

    share_at_upper = float(np.mean(identified & (rho_raw >= high)))
    share_at_lower = float(np.mean(identified & (rho_raw <= low)))
    share_above_one = float(np.mean(identified & (rho_unconstrained_raw > 1.0 + _ABOVE_ONE_TOL)))
    pileup = bool(share_at_upper >= threshold)
    row: dict[str, Any] = {
        "outcome": label,
        "window_first": first,
        "window_last": last,
        "n_draws": n_draws,
        "n_identified": int(identified.sum()),
        "degenerate_policy": degenerate_policy,
    }
    row.update(_summary_or_nan(rho, qs, "rho_"))
    row.update(_summary_or_nan(rho_unconstrained, qs, "rho_unconstrained_"))
    row.update(_summary_or_nan(lam, qs, "lambda_"))
    row.update(
        share_at_upper=share_at_upper,
        share_at_lower=share_at_lower,
        share_degenerate=n_degenerate / n_draws,
        share_unconstrained_above_one=share_above_one,
        boundary_pileup_flag=pileup,
    )
    event_time = np.arange(1, n_periods + 1)
    in_window = (event_time >= first) & (event_time <= last)
    weight_full = np.full(n_periods, np.nan)
    weight_full[in_window] = weights
    support_df = pd.DataFrame({
        "event_time": event_time,
        "support": support_full,
        "in_window": in_window,
        "weight": weight_full,
        "weight_share": weight_full / float(np.sum(weights)),
    })
    health = [
        _count_row("posterior draws", n_draws),
        _count_row("event times in the trajectory", n_periods),
        _count_row("window periods", last - first + 1),
        _count_row("identified draws", int(identified.sum())),
        _count_row("degenerate draws", n_degenerate),
    ]
    health.extend(_distribution_row(f"tau(s={s})", draws[:, s - 1]) for s in range(first, last + 1))
    health.extend([
        _distribution_row("rho", rho),
        _distribution_row("rho_unconstrained", rho_unconstrained),
        _distribution_row("lambda", lam),
    ])
    return PostPulseRetentionResult(
        label=str(label),
        window=(first, last),
        weights=weights,
        support_source=support_source,
        trajectory_source=trajectory_source,
        rho_bounds=(low, high),
        unconstrained_upper=float(grid_unconstrained[-1]),
        grid_step=float(step),
        degenerate_policy=degenerate_policy,
        seed=seed_value,
        rho_draws=rho,
        rho_unconstrained_draws=rho_unconstrained,
        lambda_draws=lam,
        degenerate=degenerate,
        share_at_upper=share_at_upper,
        share_at_lower=share_at_lower,
        share_degenerate=n_degenerate / n_draws,
        share_unconstrained_above_one=share_above_one,
        pileup_threshold=threshold,
        boundary_pileup_flag=pileup,
        summary_df=pd.DataFrame([row]),
        support_df=support_df,
        data_health_certificate_df=_health_frame(health),
    )


# ---- Within-arm correlation of per-unit post-pulse slopes ----------------------------------


@dataclass(frozen=True, eq=False)
class SlopeCorrelationResult:
    """Output of :func:`unit_slope_correlation`.

    Attributes:
      outcomes: Outcome labels, in the order of ``correlation_matrix``.
      arms_used: Arm labels whose units enter (each correlated within its own arm).
      slope_method: ``"ols"`` or ``"delta"``.
      correlation_type: ``"pearson"`` or ``"spearman"``.
      min_obs: Minimum observed window periods per outcome of an eligible unit.
      correlation_matrix: ``(M, M)`` within-arm matrices pooled by Fisher-z averaging with weights
        ``n_arm - 3`` (unit diagonal). For ``M >= 3`` an entrywise Fisher-z average of positive
        semi-definite matrices need not be positive semi-definite; :func:`pool_multi_outcome_retention`
        rejects such a matrix.
      r_hat: Mean off-diagonal entry of ``correlation_matrix``.
      per_arm_matrix: ``{arm: (M, M)}`` within-arm correlation matrices.
      per_arm_r_hat: ``{arm: mean off-diagonal}``.
      n_units: Units in the input.
      n_in_arms: Units in ``arms_used``.
      n_eligible: Units used (``>= min_obs`` observed periods in every outcome, finite slopes).
      n_eligible_by_arm: ``{arm: units used}``.
      n_dropped_nonfinite: Units with enough periods but a non-finite slope (for instance a
        ``"delta"`` half without observations), dropped.
      eligible: ``(N,)`` mask of the units used.
      slopes: ``(n_eligible, M)`` per-unit slopes.
      summary_df: One row per outcome pair with the pooled and per-arm correlations.
      data_health_certificate_df: Row counts and the distribution of the window values and slopes.
    """

    outcomes: tuple[str, ...]
    arms_used: tuple[Any, ...]
    slope_method: str
    correlation_type: str
    min_obs: int
    correlation_matrix: np.ndarray
    r_hat: float
    per_arm_matrix: dict[Any, np.ndarray]
    per_arm_r_hat: dict[Any, float]
    n_units: int
    n_in_arms: int
    n_eligible: int
    n_eligible_by_arm: dict[Any, int]
    n_dropped_nonfinite: int
    eligible: np.ndarray
    slopes: np.ndarray
    summary_df: pd.DataFrame
    data_health_certificate_df: pd.DataFrame


def _unit_slopes(window: np.ndarray, periods: np.ndarray, method: str, early_periods: int) -> np.ndarray:
    """Per-unit slope over the observed window periods (NaN-aware).

    ``"ols"`` is the OLS slope of the outcome on the period over the observed periods; ``"delta"``
    the mean of the observed late periods (positions ``>= early_periods``) minus the mean of the
    observed early periods.
    """
    observed = np.isfinite(window)
    values = np.where(observed, window, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        if method == "ols":
            counts = observed.sum(axis=1, keepdims=True)
            period_obs = np.where(observed, periods[None, :], 0.0)
            period_dev = np.where(observed, periods[None, :] - period_obs.sum(1, keepdims=True) / counts, 0.0)
            value_dev = np.where(observed, window - values.sum(1, keepdims=True) / counts, 0.0)
            return (period_dev * value_dev).sum(axis=1) / (period_dev**2).sum(axis=1)
        late = np.arange(window.shape[1]) >= early_periods
        late_mean = values[:, late].sum(1) / observed[:, late].sum(1)
        early_mean = values[:, ~late].sum(1) / observed[:, ~late].sum(1)
        return late_mean - early_mean


def _within_arm_matrix(arm_slopes: np.ndarray, correlation: str) -> np.ndarray:
    """Pearson or Spearman (Pearson of average ranks) correlation matrix of the columns."""
    if correlation == "spearman":
        arm_slopes = sp_stats.rankdata(arm_slopes, axis=0)
    return np.atleast_2d(np.corrcoef(arm_slopes, rowvar=False))


def unit_slope_correlation(
    y_event_by_outcome: Mapping[str, Any],
    arms: Any,
    *,
    arms_used: Sequence[Any] = ("A", "C"),
    periods: Sequence[float] | None = None,
    min_obs: int = 5,
    slope_method: str = "ols",
    correlation: str = "pearson",
    delta_early_periods: int | None = None,
) -> SlopeCorrelationResult:
    """Within-arm correlation across outcomes of per-unit post-pulse slopes.

    The retention estimates of several outcomes measured on the same units are arm-mean
    differences of the same units, so they inherit the within-arm correlation of the units' own
    slopes; this is the proxy :func:`pool_multi_outcome_retention` needs to stop
    treating the outcomes as independent replicates. Every unit of ``arms_used`` with at least
    ``min_obs`` observed window periods in every outcome gets one slope per outcome; the
    correlation matrix is computed within each arm (so the treatment effect does not inflate it)
    and pooled across arms by Fisher-z averaging with weights ``n_arm - 3``; ``r_hat`` is the mean
    off-diagonal entry.

    Args:
      y_event_by_outcome: ``{outcome: (N, M)}`` outcome levels over the window periods (NaN =
        unobserved), at least two outcomes with a common shape.
      arms: ``(N,)`` arm labels.
      arms_used: Arm labels whose units enter (each needs at least 4 eligible units).
      periods: ``(M,)`` strictly increasing event times of the window columns (default ``1..M``);
        only their differences matter for ``"ols"``.
      min_obs: Minimum observed window periods (in every outcome) of an eligible unit, ``>= 2``.
      slope_method: ``"ols"`` (slope of the outcome on the event time) or ``"delta"`` (late minus
        early mean).
      correlation: ``"pearson"`` or ``"spearman"``.
      delta_early_periods: Leading window columns forming the early half of ``"delta"`` (default
        ``M // 2``).

    Returns:
      :class:`SlopeCorrelationResult`.

    Raises:
      ValueError: On malformed inputs, an arm with fewer than 4 eligible units, or an outcome whose
        slopes are constant within an arm.
    """
    if not isinstance(y_event_by_outcome, Mapping) or len(y_event_by_outcome) < 2:
        raise ValueError("y_event_by_outcome must map at least two outcome labels to (N, M) window arrays")
    if slope_method not in SLOPE_METHODS:
        raise ValueError(f"slope_method must be one of {SLOPE_METHODS}; got {slope_method!r}")
    if correlation not in CORRELATION_TYPES:
        raise ValueError(f"correlation must be one of {CORRELATION_TYPES}; got {correlation!r}")
    windows: dict[str, np.ndarray] = {}
    for key, value in y_event_by_outcome.items():
        try:
            window = np.asarray(value, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"y_event_by_outcome[{key!r}] must be a numeric (N, M) array") from exc
        if window.ndim != 2:
            raise ValueError(f"y_event_by_outcome[{key!r}] must be 2-D (units, periods); got shape {window.shape}")
        if np.any(np.isinf(window)):
            raise ValueError(f"y_event_by_outcome[{key!r}] contains infinite entries")
        windows[str(key)] = window
    if len(windows) != len(y_event_by_outcome):
        raise ValueError("y_event_by_outcome keys must stay distinct after conversion to str")
    shapes = {window.shape for window in windows.values()}
    if len(shapes) != 1:
        raise ValueError(f"every outcome window must have the same (N, M) shape; got {sorted(shapes)}")
    n_units, n_window = next(iter(shapes))
    if n_window < 2:
        raise ValueError(f"the window needs at least 2 periods; got {n_window}")
    arm_labels = _as_labels(arms, n_units, "arms")
    used = tuple(arms_used)
    if not used or len(set(used)) != len(used):
        raise ValueError(f"arms_used must hold distinct arm labels; got {arms_used!r}")
    present = set(arm_labels.tolist())
    missing = [arm for arm in used if arm not in present]
    if missing:
        raise ValueError(f"arms_used labels {missing} do not occur in arms")
    period_grid = np.arange(1.0, n_window + 1.0) if periods is None else _as_vector(periods, "periods", n_window)
    if np.any(np.diff(period_grid) <= 0.0):
        raise ValueError("periods must be strictly increasing")
    min_count = _as_int(min_obs, "min_obs", 2)
    if min_count > n_window:
        raise ValueError(f"min_obs ({min_count}) exceeds the number of window periods ({n_window})")
    if delta_early_periods is None:
        early_periods = n_window // 2
    else:
        early_periods = _as_int(delta_early_periods, "delta_early_periods")
    if not 1 <= early_periods <= n_window - 1:
        raise ValueError(f"delta_early_periods must lie in [1, {n_window - 1}]; got {delta_early_periods!r}")

    in_arms = np.isin(arm_labels, list(used))
    observed_periods = np.min(np.stack([np.isfinite(window).sum(axis=1) for window in windows.values()]), axis=0)
    candidate = in_arms & (observed_periods >= min_count)
    slopes_all = np.column_stack(
        [_unit_slopes(window[candidate], period_grid, slope_method, early_periods) for window in windows.values()]
    )
    finite_rows = np.all(np.isfinite(slopes_all), axis=1)
    eligible = candidate.copy()
    eligible[np.nonzero(candidate)[0][~finite_rows]] = False
    slopes = slopes_all[finite_rows]
    eligible_arms = arm_labels[eligible]
    outcomes = tuple(windows)
    n_outcomes = len(outcomes)
    z_sum = np.zeros((n_outcomes, n_outcomes))
    weight_sum = 0.0
    per_arm_matrix: dict[Any, np.ndarray] = {}
    per_arm_r: dict[Any, float] = {}
    n_by_arm: dict[Any, int] = {}
    for arm in used:
        arm_slopes = slopes[eligible_arms == arm]
        n_by_arm[arm] = int(arm_slopes.shape[0])
        if arm_slopes.shape[0] < 4:
            raise ValueError(f"arm {arm!r} has {arm_slopes.shape[0]} eligible units; Fisher-z pooling needs >= 4")
        constant = [outcomes[j] for j in np.nonzero(np.ptp(arm_slopes, axis=0) == 0.0)[0]]
        if constant:
            raise ValueError(f"slopes of {constant} are constant within arm {arm!r}; the correlation is undefined")
        matrix = _within_arm_matrix(arm_slopes, correlation)
        weight = float(arm_slopes.shape[0] - 3)
        z_sum += weight * np.arctanh(np.clip(matrix, -1.0 + _FISHER_Z_CLIP, 1.0 - _FISHER_Z_CLIP))
        weight_sum += weight
        per_arm_matrix[arm] = matrix
        per_arm_r[arm] = _mean_off_diagonal(matrix)
    pooled = np.tanh(z_sum / weight_sum)
    np.fill_diagonal(pooled, 1.0)
    if not (np.allclose(pooled, pooled.T) and np.all(np.abs(pooled) <= 1.0) and np.allclose(np.diag(pooled), 1.0)):
        raise AssertionError("the pooled slope correlation matrix is not symmetric with a unit diagonal and |r| <= 1")

    pair_rows = []
    for i in range(n_outcomes):
        for j in range(i + 1, n_outcomes):
            pair = {"outcome_a": outcomes[i], "outcome_b": outcomes[j], "r_pooled": float(pooled[i, j])}
            pair.update({f"r_{arm}": float(per_arm_matrix[arm][i, j]) for arm in used})
            pair_rows.append(pair)
    health = [
        _count_row("units (input)", n_units),
        _count_row("units in arms_used", int(in_arms.sum())),
        _count_row(f"units with >= {min_count} observed periods in every outcome", int(candidate.sum())),
        _count_row("units with finite slopes (eligible)", int(eligible.sum())),
    ]
    health.extend(_count_row(f"eligible units in arm {arm}", count) for arm, count in n_by_arm.items())
    for column, (name, window) in enumerate(windows.items()):
        health.append(_distribution_row(f"window values: {name} (eligible)", window[eligible]))
        health.append(_distribution_row(f"slope ({slope_method}): {name}", slopes[:, column]))
    return SlopeCorrelationResult(
        outcomes=outcomes,
        arms_used=used,
        slope_method=slope_method,
        correlation_type=correlation,
        min_obs=min_count,
        correlation_matrix=pooled,
        r_hat=_mean_off_diagonal(pooled),
        per_arm_matrix=per_arm_matrix,
        per_arm_r_hat=per_arm_r,
        n_units=n_units,
        n_in_arms=int(in_arms.sum()),
        n_eligible=int(eligible.sum()),
        n_eligible_by_arm=n_by_arm,
        n_dropped_nonfinite=int(np.sum(~finite_rows)),
        eligible=eligible,
        slopes=slopes,
        summary_df=pd.DataFrame(pair_rows),
        data_health_certificate_df=_health_frame(health),
    )


# ---- Multi-outcome pooling of lambda = -ln(rho) ---------------------------------------------


@dataclass(frozen=True, eq=False)
class PooledRetentionResult:
    """Output of :func:`pool_multi_outcome_retention`.

    Moments, quantiles and every P10 use each outcome's finite draws only: NaN draws (degenerate
    windows under ``degenerate_policy="nan"``) are excluded, and ``degenerate_shares`` reports how
    many there were.

    Attributes:
      outcomes: Pooling outcomes, in input order.
      focal_label: Label of the focal outcome (shrunk, not pooled), or ``None``.
      rule_outcomes: Outcomes whose pooled ``rho`` P10 forms the rule statistic: the focal outcome,
        or every pooling outcome (minimum over them) when there is none.
      lambda_means: Posterior mean of ``lambda`` per outcome (pooling outcomes and focal).
      lambda_variances: Posterior variance (ddof 1) of ``lambda`` per outcome.
      q_statistic: Cochran's ``Q`` of the pooling outcomes.
      q_dof: Its degrees of freedom ``M - 1``.
      q_p_value_upper: ``1 - F_{chi2(M - 1)}(Q)`` (heterogeneity).
      q_p_value_lower: ``F_{chi2(M - 1)}(Q)``; small values flag under-dispersion, which positively
        correlated outcomes produce.
      tau2: DerSimonian-Laird between-outcome variance ``tau^2-hat``.
      i2: ``I^2 = max(0, (Q - df) / Q)``.
      weight_shares: Normalized random-effects weights of the pooling outcomes.
      dominant_outcome: Pooling outcome with the largest weight.
      dominant_weight: Its weight share.
      mu_hat: Pooled mean of ``lambda``.
      var_mu_independent: ``Var(mu-hat) = 1 / sum w`` (outcomes treated as independent).
      rho_at_mu: ``exp(-mu_hat)``.
      correlation_matrix: The ``(M, M)`` correlation used, or ``None``.
      r_common: Common correlation used (the mean off-diagonal entry of a matrix), or ``None``.
      var_mu_correlated: ``Var_R(mu-hat)`` at ``correlation_matrix``, or ``None``.
      var_mu_common: ``Var_R(mu-hat)`` at ``r_common``, or ``None``.
      q_moment_r: Moment correlation ``max(0, 1 - Q / (M - 1))``.
      shrinkage_b: Morris shrinkage factor ``B = v / (v + tau^2)`` per outcome (literal pooling).
      kappa: Spread factor of the draw-preserving shrinkage per outcome (literal pooling).
      clip_share: Share of pooled ``lambda`` draws below 0 before clipping, per outcome.
      pooled_lambda_draws: Literal (``r = 0``) pooled ``lambda`` draws per outcome.
      pooled_rho_draws: ``exp(-pooled_lambda_draws)`` per outcome.
      adjustment: ``"matrix"`` or ``"common"``: the correlation behind the adjusted draws, or ``None``.
      pooled_rho_draws_adjusted: Pooled ``rho`` draws at that correlation, or ``None``.
      literal_p10: Rule statistic (pooled ``rho`` P10) of the literal pooling.
      adjusted_p10: Rule statistic at the supplied correlation -- the matrix when one was given
        (directly or through a :class:`SlopeCorrelationResult`), else the common ``r`` -- or
        ``None`` without a correlation.
      p10_threshold: Threshold of the ``passed`` column and of ``breakeven_r``.
      breakeven_r: Common correlation in ``[0, 1]`` at which the rule statistic reaches the
        threshold; ``None`` when the literal statistic does not exceed it or the statistic still
        exceeds it at ``r = 1``.
      rule_variants_df: One row per variant: ``variant`` (``literal``, ``correlation_matrix``,
        ``correlation_common``, ``q_moment``, ``leave_one_out``, ``equal_weights``,
        ``equal_weights_correlated``), ``label``, ``left_out``, ``correlation``, ``tau2``,
        ``mu_hat``, ``var_mu``, ``se_mu``, ``rho_at_mu``, ``rule_outcome``, ``rho_P10`` and
        ``passed`` (``rho_P10 > p10_threshold``).
      cap_hit_shares: ``share_at_upper`` per outcome when the inputs were
        :class:`PostPulseRetentionResult` objects (only those outcomes), else ``None``.
      degenerate_shares: Share of draws per outcome (pooling and focal) whose ``rho`` is
        unidentified: ``share_degenerate`` of a :class:`PostPulseRetentionResult` input (under
        either policy), else the share of NaN draws.
      with_mu_uncertainty: Whether the Morris ``Var(mu-hat)`` term widens the pooled draws.
      clip_lambda_at_zero: Whether pooled ``lambda`` draws were clipped at 0 (``rho <= 1``).
      warnings: Warnings raised by the pooling (``tau^2-hat = 0``, under-dispersed ``Q``).
      summary_df: One row per outcome with its role, moments, weight, ``B``, ``kappa``, clip,
        cap-hit and degenerate shares and the unpooled, pooled and (if any) adjusted ``rho``
        quantiles.
      data_health_certificate_df: Row counts and the distribution of every input and pooled series.
    """

    outcomes: tuple[str, ...]
    focal_label: str | None
    rule_outcomes: tuple[str, ...]
    lambda_means: dict[str, float]
    lambda_variances: dict[str, float]
    q_statistic: float
    q_dof: int
    q_p_value_upper: float
    q_p_value_lower: float
    tau2: float
    i2: float
    weight_shares: dict[str, float]
    dominant_outcome: str
    dominant_weight: float
    mu_hat: float
    var_mu_independent: float
    rho_at_mu: float
    correlation_matrix: np.ndarray | None
    r_common: float | None
    var_mu_correlated: float | None
    var_mu_common: float | None
    q_moment_r: float
    shrinkage_b: dict[str, float]
    kappa: dict[str, float]
    clip_share: dict[str, float]
    pooled_lambda_draws: dict[str, np.ndarray]
    pooled_rho_draws: dict[str, np.ndarray]
    adjustment: str | None
    pooled_rho_draws_adjusted: dict[str, np.ndarray] | None
    literal_p10: float
    adjusted_p10: float | None
    p10_threshold: float
    breakeven_r: float | None
    rule_variants_df: pd.DataFrame
    cap_hit_shares: dict[str, float] | None
    degenerate_shares: dict[str, float]
    with_mu_uncertainty: bool
    clip_lambda_at_zero: bool
    warnings: tuple[str, ...]
    summary_df: pd.DataFrame
    data_health_certificate_df: pd.DataFrame


@dataclass(frozen=True, eq=False)
class _Pooling:
    """DerSimonian-Laird pooling of per-outcome lambda draws (internal)."""

    means: np.ndarray
    variances: np.ndarray
    q_stat: float
    dof: int
    tau2: float
    i2: float
    weights: np.ndarray
    mu: float
    var_mu: float
    p_upper: float
    p_lower: float


@dataclass(frozen=True, eq=False)
class _ShrinkTarget:
    """One outcome whose pooled rho P10 enters the rule statistic (internal)."""

    label: str
    draws: np.ndarray
    mean: float
    variance: float


@dataclass(frozen=True)
class _PoolCentre:
    """What a variant shrinks toward: ``tau^2``, ``mu-hat`` and ``Var(mu-hat)`` (internal)."""

    tau2: float
    mu: float
    var_mu: float


@dataclass(frozen=True)
class _ShrinkOptions:
    """Shrinkage switches and the P10 threshold shared by every variant (internal)."""

    with_mu_uncertainty: bool
    clip: bool
    threshold: float


@dataclass(frozen=True, eq=False)
class _CorrelationInputs:
    """The supplied correlation and its ``Var_R(mu-hat)`` values; ``None`` fields are absent (internal)."""

    matrix: np.ndarray | None
    r_common: float | None
    var_matrix: float | None
    var_common: float | None


def _dersimonian_laird(samples: Sequence[np.ndarray]) -> _Pooling:
    """DerSimonian-Laird random-effects pooling of the per-sample means and variances.

    A single sample is its own pool (``Q = 0``, ``tau^2 = 0``, NaN p-values), which the
    leave-one-out variant of two outcomes needs.
    """
    means = np.array([float(np.mean(sample)) for sample in samples])
    variances = np.array([float(np.var(sample, ddof=1)) for sample in samples])
    dof = len(samples) - 1
    if dof == 0:
        q_stat, tau2, i2, p_upper, p_lower = 0.0, 0.0, 0.0, float("nan"), float("nan")
    else:
        fixed_weights = 1.0 / variances
        total_weight = float(np.sum(fixed_weights))
        mu_fixed = float(np.sum(fixed_weights * means) / total_weight)
        q_stat = float(np.sum(fixed_weights * (means - mu_fixed) ** 2))
        scale = total_weight - float(np.sum(fixed_weights**2)) / total_weight
        tau2 = max(0.0, (q_stat - dof) / scale) if scale > 0.0 else 0.0
        i2 = max(0.0, (q_stat - dof) / q_stat) if q_stat > 0.0 else 0.0
        p_upper = float(sp_stats.chi2.sf(q_stat, dof))
        p_lower = float(sp_stats.chi2.cdf(q_stat, dof))
    random_weights = 1.0 / (variances + tau2)
    return _Pooling(
        means=means,
        variances=variances,
        q_stat=q_stat,
        dof=dof,
        tau2=tau2,
        i2=i2,
        weights=random_weights,
        mu=float(np.sum(random_weights * means) / np.sum(random_weights)),
        var_mu=float(1.0 / np.sum(random_weights)),
        p_upper=p_upper,
        p_lower=p_lower,
    )


def _var_mu_common(variances: np.ndarray, tau2: float, correlation: float) -> float:
    """``Var_R(mu-hat)`` of the random-effects mean under a common correlation ``r`` (``r = 0``: ``1 / sum w``)."""
    weights = 1.0 / (variances + tau2)
    scaled_sd = weights * np.sqrt(variances)
    diagonal = float(np.sum(weights**2 * (variances + tau2)))
    off_diagonal = float(np.sum(scaled_sd) ** 2 - np.sum(scaled_sd**2))
    return (diagonal + correlation * off_diagonal) / float(np.sum(weights)) ** 2


def _var_mu_weighted(weights: np.ndarray, variances: np.ndarray, tau2: float, matrix: np.ndarray) -> float:
    """``[sum w^2 (v + tau^2) + sum_{j != k} R_jk w_j w_k sqrt(v_j v_k)] / (sum w)^2`` for any weights."""
    scaled_sd = weights * np.sqrt(variances)
    diagonal = float(np.sum(weights**2 * (variances + tau2)))
    off_diagonal = float(scaled_sd @ (matrix - np.diag(np.diag(matrix))) @ scaled_sd)
    return (diagonal + off_diagonal) / float(np.sum(weights)) ** 2


def _shrink_toward_pool(
    lam: np.ndarray, mean: float, variance: float, centre: _PoolCentre, options: _ShrinkOptions
) -> tuple[np.ndarray, float, float, float]:
    """Draw-preserving shrinkage of one outcome's lambda draws toward ``mu-hat`` (Morris 1983).

    ``B = v / (v + tau^2)`` and ``lambda' = (1 - B) m + B mu + kappa (lambda - m)`` with the Morris
    spread ``kappa = sqrt(((1 - B) v + B^2 Var(mu-hat)) / v)`` or ``sqrt(1 - B)``; NaN draws stay
    NaN.

    Args:
      lam: ``(S,)`` lambda draws (NaN allowed).
      mean: Their posterior mean ``m``.
      variance: Their posterior variance ``v``.
      centre: ``tau^2``, ``mu-hat`` and ``Var(mu-hat)`` of the pool.
      options: Whether to use the Morris spread and to clip the pooled draws at 0.

    Returns:
      ``(pooled draws, B, kappa, share of finite pooled draws below 0 before clipping)``.
    """
    shrink = variance / (variance + centre.tau2)
    if options.with_mu_uncertainty:
        spread_var = (1.0 - shrink) * variance + shrink**2 * centre.var_mu
        spread = float(np.sqrt(spread_var / variance))
    else:
        spread = float(np.sqrt(1.0 - shrink))
    pooled = (1.0 - shrink) * mean + shrink * centre.mu + spread * (lam - mean)
    finite = pooled[np.isfinite(pooled)]
    clip_share = float(np.mean(finite < 0.0)) if finite.size else float("nan")
    if options.clip:
        pooled = np.maximum(pooled, 0.0)
    return pooled, float(shrink), spread, clip_share


def _rule_p10(targets: Sequence[_ShrinkTarget], centre: _PoolCentre, options: _ShrinkOptions) -> tuple[float, str]:
    """Rule statistic: the smallest pooled ``rho`` P10 (finite draws) over the targets, and its target."""
    best_p10, best_label = math.inf, ""
    for target in targets:
        pooled, _, _, _ = _shrink_toward_pool(target.draws, target.mean, target.variance, centre, options)
        p10 = float(np.quantile(np.exp(-pooled[np.isfinite(pooled)]), _RULE_QUANTILE))
        if p10 < best_p10:
            best_p10, best_label = p10, target.label
    return best_p10, best_label


def _lambda_input(value: Any, name: str) -> tuple[np.ndarray, float | None, float]:
    """``(lambda draws with NaN allowed, cap-hit share or None, degenerate share)`` from draws or a retention result."""
    if isinstance(value, PostPulseRetentionResult):
        lam = np.asarray(value.lambda_draws, dtype=np.float64)
        return lam, float(value.share_at_upper), float(value.share_degenerate)
    lam = _as_vector(value, name, allow_nan=True)
    return lam, None, float(np.mean(np.isnan(lam)))


def _finite_moments(lam: np.ndarray, name: str) -> tuple[np.ndarray, float, float]:
    """``(finite draws, mean, variance ddof 1)``; at least two finite, non-constant draws."""
    finite = lam[np.isfinite(lam)]
    if finite.size < 2:
        raise ValueError(f"{name} needs at least 2 finite lambda draws; got {finite.size}")
    variance = float(np.var(finite, ddof=1))
    if float(np.ptp(finite)) == 0.0 or not variance > 0.0:
        raise ValueError(f"{name} has constant lambda draws; pooling needs a positive posterior variance")
    return finite, float(np.mean(finite)), variance


def _check_correlation_matrix(matrix: np.ndarray, n_outcomes: int) -> np.ndarray:
    """Validate a finite, symmetric, positive semi-definite ``(M, M)`` correlation with a unit diagonal."""
    if matrix.shape != (n_outcomes, n_outcomes):
        raise ValueError(
            f"correlation must be ({n_outcomes}, {n_outcomes}), one row per pooling outcome; got {matrix.shape}"
        )
    if not np.all(np.isfinite(matrix)):
        raise ValueError("correlation contains NaN / Inf entries")
    if not np.allclose(matrix, matrix.T, atol=1e-10) or not np.allclose(np.diag(matrix), 1.0, atol=1e-10):
        raise ValueError("correlation must be symmetric with a unit diagonal")
    if np.any(np.abs(matrix) > 1.0 + 1e-12):
        raise ValueError("correlation entries must lie in [-1, 1]")
    min_eigenvalue = float(np.linalg.eigvalsh(matrix).min())
    if min_eigenvalue < -_PSD_TOL:
        raise ValueError(
            f"correlation must be positive semi-definite; its smallest eigenvalue is {min_eigenvalue:.3g}"
        )
    return matrix


def _resolve_correlation(correlation: Any, outcomes: Sequence[str]) -> tuple[np.ndarray | None, float | None]:
    """``(matrix or None, common r or None)`` from a float, an ``(M, M)`` matrix or a slope correlation."""
    if correlation is None:
        return None, None
    if isinstance(correlation, SlopeCorrelationResult):
        if set(correlation.outcomes) != set(outcomes):
            raise ValueError(
                f"the slope correlation covers {list(correlation.outcomes)} but the pooling outcomes are "
                f"{list(outcomes)}; use the same labels"
            )
        order = [correlation.outcomes.index(name) for name in outcomes]
        reordered = np.asarray(correlation.correlation_matrix, dtype=np.float64)[np.ix_(order, order)]
        matrix = _check_correlation_matrix(reordered, len(outcomes))
        return matrix, _mean_off_diagonal(matrix)
    try:
        values = np.asarray(correlation, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"correlation must be a float, an (M, M) matrix or a SlopeCorrelationResult; got {correlation!r}"
        ) from exc
    if values.ndim == 0:
        common = float(values)
        floor = -1.0 / (len(outcomes) - 1)
        if not (math.isfinite(common) and floor <= common <= 1.0):
            raise ValueError(
                f"a common correlation of {len(outcomes)} outcomes must lie in [{floor:.4g}, 1] (positive "
                f"semi-definite); got {correlation!r}"
            )
        return None, common
    matrix = _check_correlation_matrix(values, len(outcomes))
    return matrix, _mean_off_diagonal(matrix)


def _variant_row(
    variant: str,
    label: str,
    centre: _PoolCentre,
    targets: Sequence[_ShrinkTarget],
    options: _ShrinkOptions,
    *,
    left_out: str | None = None,
    correlation: float = float("nan"),
) -> dict[str, Any]:
    """One row of the rule-variant table: pooled mean, its SE and the rule statistic."""
    if not centre.var_mu > 0.0:
        raise ValueError(
            f"Var(mu-hat) of the {variant!r} variant is not positive ({centre.var_mu:.3g}); check the correlation"
        )
    p10, rule_outcome = _rule_p10(targets, centre, options)
    return {
        "variant": variant,
        "label": label,
        "left_out": left_out,
        "correlation": correlation,
        "tau2": centre.tau2,
        "mu_hat": centre.mu,
        "var_mu": centre.var_mu,
        "se_mu": math.sqrt(centre.var_mu),
        "rho_at_mu": math.exp(-centre.mu),
        "rule_outcome": rule_outcome,
        "rho_P10": p10,
        "passed": bool(p10 > options.threshold),
    }


def _breakeven_correlation(
    targets: Sequence[_ShrinkTarget], pooling: _Pooling, options: _ShrinkOptions
) -> float | None:
    """Common ``r`` in ``[0, 1]`` where the rule statistic crosses the threshold (Brent's method)."""

    def gap(correlation: float) -> float:
        centre = _PoolCentre(pooling.tau2, pooling.mu, _var_mu_common(pooling.variances, pooling.tau2, correlation))
        return _rule_p10(targets, centre, options)[0] - options.threshold

    if not gap(0.0) > 0.0 > gap(1.0):
        return None
    return float(sp_optimize.brentq(gap, 0.0, 1.0))


def _q_moment_correlation(pooling: _Pooling) -> float:
    """Equal-variance moment estimate of a common correlation, ``max(0, 1 - Q / (M - 1))``."""
    return max(0.0, 1.0 - pooling.q_stat / pooling.dof)


def _rule_variant_rows(
    pooling: _Pooling,
    samples: Mapping[str, np.ndarray],
    targets: Sequence[_ShrinkTarget],
    options: _ShrinkOptions,
    correlation: _CorrelationInputs,
) -> list[dict[str, Any]]:
    """Rule-variant rows: literal, correlation-adjusted, Q-moment, leave-one-out and equal outcome weights.

    Args:
      pooling: DerSimonian-Laird pooling of every pooling outcome.
      samples: Finite lambda draws per pooling outcome, in pooling order.
      targets: Outcomes whose pooled P10 forms the rule statistic.
      options: Shrinkage switches and the P10 threshold.
      correlation: The supplied correlation and its ``Var_R(mu-hat)`` values.

    Returns:
      One dict per variant (see ``PooledRetentionResult.rule_variants_df``).
    """
    outcomes = tuple(samples)
    tau2, mu, r_common = pooling.tau2, pooling.mu, correlation.r_common
    rows = [_variant_row(
        "literal", "Pre-registered (literal): outcomes treated as independent", _PoolCentre(tau2, mu, pooling.var_mu),
        targets, options,
    )]
    if correlation.var_matrix is not None:
        rows.append(_variant_row(
            "correlation_matrix", f"Correlation-adjusted at the R matrix (mean off-diagonal {r_common:.3f})",
            _PoolCentre(tau2, mu, correlation.var_matrix), targets, options, correlation=float(r_common),
        ))
    if correlation.var_common is not None:
        rows.append(_variant_row(
            "correlation_common", f"Correlation-adjusted at a common r = {r_common:.3f}",
            _PoolCentre(tau2, mu, correlation.var_common), targets, options, correlation=float(r_common),
        ))
    q_moment_r = _q_moment_correlation(pooling)
    rows.append(_variant_row(
        "q_moment", f"Correlation-adjusted at the Q-moment r = {q_moment_r:.3f} (r = 1 - Q / df)",
        _PoolCentre(tau2, mu, _var_mu_common(pooling.variances, tau2, q_moment_r)), targets, options,
        correlation=q_moment_r,
    ))
    for left_out in outcomes:
        loo = _dersimonian_laird([samples[name] for name in outcomes if name != left_out])
        rows.append(_variant_row(
            "leave_one_out", f"Leave one out: without {left_out}", _PoolCentre(loo.tau2, loo.mu, loo.var_mu), targets,
            options, left_out=left_out,
        ))
    mu_equal = float(np.mean(pooling.means))
    var_equal = float(np.sum(pooling.variances + tau2) / pooling.variances.size**2)
    rows.append(_variant_row(
        "equal_weights", "Equal outcome weights (R = I)", _PoolCentre(tau2, mu_equal, var_equal), targets, options,
        correlation=0.0,
    ))
    if r_common is not None:
        matrix = correlation.matrix
        if matrix is None:
            matrix = np.full((len(outcomes), len(outcomes)), r_common)
            np.fill_diagonal(matrix, 1.0)
        var_equal_correlated = _var_mu_weighted(np.ones(len(outcomes)), pooling.variances, tau2, matrix)
        rows.append(_variant_row(
            "equal_weights_correlated", "Equal outcome weights at the supplied correlation",
            _PoolCentre(tau2, mu_equal, var_equal_correlated), targets, options, correlation=float(r_common),
        ))
    return rows


@dataclass(frozen=True, eq=False)
class _PoolingInputs:
    """Validated lambda draws of the pooling outcomes, then the optional focal outcome (internal)."""

    raw: dict[str, np.ndarray]
    outcomes: tuple[str, ...]
    focal_name: str | None
    cap_hit: dict[str, float]
    degenerate: dict[str, float]


@dataclass(frozen=True, eq=False)
class _Shrunk:
    """Draw-preserving shrinkage of every outcome toward one pool centre (internal)."""

    lambda_draws: dict[str, np.ndarray]
    rho_draws: dict[str, np.ndarray]
    shrinkage_b: dict[str, float]
    kappa: dict[str, float]
    clip_share: dict[str, float]


def _pooling_inputs(lambda_draws: Any, focal: Any, focal_label: str | None) -> _PoolingInputs:
    """Validate the pooling and focal inputs and collect their cap-hit and degenerate shares."""
    if not isinstance(lambda_draws, Mapping) or len(lambda_draws) < 2:
        raise ValueError("lambda_draws must map at least two pooling outcomes to (S,) lambda draws")
    raw: dict[str, np.ndarray] = {}
    cap_hit: dict[str, float] = {}
    degenerate: dict[str, float] = {}
    for key, value in lambda_draws.items():
        lam, share, degenerate_share = _lambda_input(value, f"lambda_draws[{key!r}]")
        raw[str(key)] = lam
        degenerate[str(key)] = degenerate_share
        if share is not None:
            cap_hit[str(key)] = share
    if len(raw) != len(lambda_draws):
        raise ValueError("lambda_draws keys must stay distinct after conversion to str")
    outcomes = tuple(raw)
    focal_name: str | None = None
    if focal is not None:
        lam, share, degenerate_share = _lambda_input(focal, "focal")
        default_label = focal.label if isinstance(focal, PostPulseRetentionResult) else "focal"
        focal_name = str(default_label if focal_label is None else focal_label)
        if focal_name in raw:
            raise ValueError(
                f"focal_label {focal_name!r} collides with a pooling outcome; the focal outcome is not pooled"
            )
        raw[focal_name] = lam
        degenerate[focal_name] = degenerate_share
        if share is not None:
            cap_hit[focal_name] = share
    elif focal_label is not None:
        raise ValueError("focal_label was given without focal draws")
    return _PoolingInputs(raw, outcomes, focal_name, cap_hit, degenerate)


def _shrink_every_outcome(
    raw: Mapping[str, np.ndarray],
    moments: Mapping[str, tuple[np.ndarray, float, float]],
    centre: _PoolCentre,
    options: _ShrinkOptions,
) -> _Shrunk:
    """:func:`_shrink_toward_pool` applied to every outcome (pooling and focal)."""
    lambda_draws: dict[str, np.ndarray] = {}
    shrinkage_b: dict[str, float] = {}
    kappa: dict[str, float] = {}
    clip_share: dict[str, float] = {}
    for name, lam in raw.items():
        lambda_draws[name], shrinkage_b[name], kappa[name], clip_share[name] = _shrink_toward_pool(
            lam, moments[name][1], moments[name][2], centre, options
        )
    rho_draws = {name: np.exp(-draws) for name, draws in lambda_draws.items()}
    return _Shrunk(lambda_draws, rho_draws, shrinkage_b, kappa, clip_share)


def _pooling_messages(pooling: _Pooling, focal_name: str | None) -> list[str]:
    """Warnings of the pooling: ``tau^2-hat = 0`` and an under-dispersed Cochran's ``Q``."""
    messages = []
    if pooling.tau2 == 0.0:
        who = f"the focal outcome {focal_name!r}" if focal_name is not None else "every outcome"
        messages.append(
            f"tau^2-hat = 0, so B = 1 for {who}: it carries no weight of its own and its pooled rho is the "
            "pooling outcomes' mean exp(-mu-hat), spread only through Var(mu-hat)."
        )
    if math.isfinite(pooling.p_lower) and pooling.p_lower < 0.05:
        messages.append(
            f"Cochran's Q = {pooling.q_stat:.3g} on {pooling.dof} df lies in the lower chi2 tail (p = "
            f"{pooling.p_lower:.3g}): the outcomes spread less than independent sampling implies, as positively "
            "correlated outcomes do, so Var(mu-hat) under independence is understated."
        )
    return messages


def _pooling_tables(
    inputs: _PoolingInputs,
    moments: Mapping[str, tuple[np.ndarray, float, float]],
    weight_shares: Mapping[str, float],
    literal: _Shrunk,
    adjusted_rho: Mapping[str, np.ndarray] | None,
    quantiles: Sequence[float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``(summary_df, data_health_certificate_df)`` of the pooling, one summary row per outcome."""
    rows = []
    health = [
        _count_row("pooling outcomes", len(inputs.outcomes)),
        _count_row("focal outcomes", int(inputs.focal_name is not None)),
    ]
    for name, lam in inputs.raw.items():
        row: dict[str, Any] = {
            "outcome": name,
            "role": "focal" if name == inputs.focal_name else "pooled",
            "n_draws": int(lam.size),
            "n_finite": int(np.sum(np.isfinite(lam))),
            "lambda_mean": moments[name][1],
            "lambda_variance": moments[name][2],
            "weight_share": weight_shares.get(name, float("nan")),
            "shrinkage_b": literal.shrinkage_b[name],
            "kappa": literal.kappa[name],
            "clip_share": literal.clip_share[name],
            "cap_hit_share": inputs.cap_hit.get(name, float("nan")),
            "degenerate_share": inputs.degenerate[name],
        }
        row.update(_summary_or_nan(np.exp(-lam), quantiles, "unpooled_rho_"))
        row.update(_summary_or_nan(literal.rho_draws[name], quantiles, "pooled_rho_"))
        if adjusted_rho is not None:
            row.update(_summary_or_nan(adjusted_rho[name], quantiles, "adjusted_rho_"))
        rows.append(row)
        health.append(_distribution_row(f"lambda: {name}", lam))
        health.append(_distribution_row(f"pooled lambda (literal): {name}", literal.lambda_draws[name]))
    return pd.DataFrame(rows), _health_frame(health)


def pool_multi_outcome_retention(
    lambda_draws: Mapping[str, Any],
    *,
    focal: Any = None,
    focal_label: str | None = None,
    correlation: Any = None,
    with_mu_uncertainty: bool = True,
    clip_lambda_at_zero: bool = True,
    p10_threshold: float = 0.90,
    quantiles: Sequence[float] = RETENTION_QUANTILES,
) -> PooledRetentionResult:
    """Pool ``lambda = -ln(rho)`` across outcomes and test the pooled ``rho`` P10 under robustness variants.

    **DerSimonian-Laird.** On the per-outcome posterior means ``m_j`` and variances ``v_j`` of
    ``lambda``: ``w_j^FE = 1 / v_j``, ``Q = sum w_j^FE (m_j - mu_FE)^2``,
    ``tau^2 = max(0, (Q - (M - 1)) / (sum w^FE - sum (w^FE)^2 / sum w^FE))``,
    ``w_j = 1 / (v_j + tau^2)``, ``mu-hat = sum w_j m_j / sum w_j`` and
    ``Var_ind(mu-hat) = 1 / sum w_j``. Outcomes measured on the same units are not independent
    studies: positively correlated estimates spread less than their sampling variances imply, so
    ``Q`` falls in the lower ``chi2`` tail and ``Var_ind`` is understated. With a correlation ``R``,
    ``Var_R(mu-hat) = [sum w_j^2 (v_j + tau^2) + sum_{j != k} R_jk w_j w_k sqrt(v_j v_k)] /
    (sum w_j)^2``, which equals ``Var_ind`` at ``R = I`` and inflates it by ``1 + (M - 1) r`` for
    equal variances, a common ``r`` and ``tau^2-hat = 0`` (or ``1 + (M - 1) r v / (v + tau^2)`` when
    ``tau^2 > 0``).

    **Draw-preserving shrinkage (Morris 1983).** The focal outcome (the outcome a decision rests on,
    excluded from the pooling set) is shrunk toward ``mu-hat``: ``B = v_f / (v_f + tau^2)`` and
    ``lambda' = (1 - B) m_f + B mu-hat + kappa (lambda_f - m_f)`` with
    ``kappa = sqrt(((1 - B) v_f + B^2 Var(mu-hat)) / v_f)`` (``sqrt(1 - B)`` without the Morris
    term), clipped at ``lambda' >= 0``. ``tau^2-hat = 0`` gives ``B = 1``: the focal outcome then
    carries no weight of its own and its pooled ``rho`` is the pooling outcomes' mean. Without a
    focal outcome every pooling outcome is shrunk and the rule statistic is the smallest of their
    P10s (in the ``leave_one_out`` variant the left-out outcome stays in ``rule_outcomes`` and is
    shrunk toward the pool of the remaining ``M - 1`` outcomes).

    **Rule variants** (``rule_variants_df``): the literal pre-registered computation (``r = 0``);
    ``Var_R`` at the correlation matrix and at its common (mean off-diagonal) ``r``; the
    Q-moment ``r_Q = max(0, 1 - Q / (M - 1))``; DerSimonian-Laird recomputed without each outcome
    in turn; equal outcome weights (``w_j = 1 / M``, at ``R = I`` and, with a correlation, at
    ``R``); plus ``breakeven_r``, the common ``r`` at which the rule statistic reaches the
    threshold.

    Args:
      lambda_draws: ``{outcome: (S,) lambda draws}`` or ``{outcome: PostPulseRetentionResult}`` for
        the ``M >= 2`` pooling outcomes. NaN draws (degenerate windows) are excluded from that
        outcome's moments and P10s and stay NaN in its pooled draws; ``degenerate_shares`` reports
        their share (for a retention result, its ``share_degenerate`` under either policy).
      focal: Optional focal outcome's lambda draws (or retention result), shrunk toward the pool.
      focal_label: Label of the focal outcome (default: the retention result's label or ``"focal"``).
      correlation: ``None``; a common correlation in ``[-1 / (M - 1), 1]``; a positive
        semi-definite ``(M, M)`` matrix ordered like ``lambda_draws``; or a
        :class:`SlopeCorrelationResult` over the same outcome labels (reordered to them). The
        adjusted draws and ``adjusted_p10`` use the matrix when there is one, else the common ``r``.
      with_mu_uncertainty: Use the Morris ``Var(mu-hat)`` term (default) or ``sqrt(1 - B)``.
      clip_lambda_at_zero: Clip pooled ``lambda`` at 0, i.e. ``rho <= 1`` (default).
      p10_threshold: Threshold of the ``passed`` column and of ``breakeven_r``, in ``(0, 1]``.
      quantiles: Probabilities of the per-outcome summary (odd length; default P10 / P50 / P90).

    Returns:
      :class:`PooledRetentionResult`.

    Raises:
      ValueError: On malformed inputs, fewer than two finite draws or zero variance for an outcome,
        a focal label that collides with a pooling outcome, a correlation that is not positive
        semi-definite, or a correlation under which ``Var(mu-hat)`` is not positive.
    """
    threshold = _as_float(p10_threshold, "p10_threshold")
    if not 0.0 < threshold <= 1.0:
        raise ValueError(f"p10_threshold must lie in (0, 1]; got {p10_threshold!r}")
    qs = _check_quantiles(quantiles)
    inputs = _pooling_inputs(lambda_draws, focal, focal_label)
    raw, outcomes, focal_name = inputs.raw, inputs.outcomes, inputs.focal_name
    moments = {name: _finite_moments(lam, f"lambda draws of {name!r}") for name, lam in raw.items()}
    pooling = _dersimonian_laird([moments[name][0] for name in outcomes])
    rule_names = (focal_name,) if focal_name is not None else outcomes
    targets = [_ShrinkTarget(name, raw[name], moments[name][1], moments[name][2]) for name in rule_names]
    options = _ShrinkOptions(bool(with_mu_uncertainty), bool(clip_lambda_at_zero), threshold)
    literal = _shrink_every_outcome(raw, moments, _PoolCentre(pooling.tau2, pooling.mu, pooling.var_mu), options)

    matrix, r_common = _resolve_correlation(correlation, outcomes)
    var_matrix = None if matrix is None else _var_mu_weighted(pooling.weights, pooling.variances, pooling.tau2, matrix)
    var_common = None if r_common is None else _var_mu_common(pooling.variances, pooling.tau2, r_common)
    samples = {name: moments[name][0] for name in outcomes}
    variants = pd.DataFrame(
        _rule_variant_rows(
            pooling, samples, targets, options, _CorrelationInputs(matrix, r_common, var_matrix, var_common)
        )
    )
    adjustment = None if r_common is None else ("matrix" if matrix is not None else "common")
    var_adjusted = var_matrix if matrix is not None else var_common
    pooled_rho_adjusted, adjusted_p10 = None, None
    if var_adjusted is not None:
        adjusted = _PoolCentre(pooling.tau2, pooling.mu, var_adjusted)
        pooled_rho_adjusted = _shrink_every_outcome(raw, moments, adjusted, options).rho_draws
        adjusted_p10 = _rule_p10(targets, adjusted, options)[0]

    messages = _pooling_messages(pooling, focal_name)
    for message in messages:
        warnings.warn(message, StabilityWarning, stacklevel=2)
    shares = pooling.weights / float(np.sum(pooling.weights))
    weight_shares = {name: float(share) for name, share in zip(outcomes, shares, strict=True)}
    dominant = max(weight_shares, key=weight_shares.__getitem__)
    summary, health = _pooling_tables(inputs, moments, weight_shares, literal, pooled_rho_adjusted, qs)
    return PooledRetentionResult(
        outcomes=outcomes,
        focal_label=focal_name,
        rule_outcomes=tuple(rule_names),
        lambda_means={name: moments[name][1] for name in raw},
        lambda_variances={name: moments[name][2] for name in raw},
        q_statistic=pooling.q_stat,
        q_dof=pooling.dof,
        q_p_value_upper=pooling.p_upper,
        q_p_value_lower=pooling.p_lower,
        tau2=pooling.tau2,
        i2=pooling.i2,
        weight_shares=weight_shares,
        dominant_outcome=dominant,
        dominant_weight=weight_shares[dominant],
        mu_hat=pooling.mu,
        var_mu_independent=pooling.var_mu,
        rho_at_mu=math.exp(-pooling.mu),
        correlation_matrix=matrix,
        r_common=r_common,
        var_mu_correlated=var_matrix,
        var_mu_common=var_common,
        q_moment_r=_q_moment_correlation(pooling),
        shrinkage_b=literal.shrinkage_b,
        kappa=literal.kappa,
        clip_share=literal.clip_share,
        pooled_lambda_draws=literal.lambda_draws,
        pooled_rho_draws=literal.rho_draws,
        adjustment=adjustment,
        pooled_rho_draws_adjusted=pooled_rho_adjusted,
        literal_p10=float(variants["rho_P10"].iloc[0]),
        adjusted_p10=adjusted_p10,
        p10_threshold=threshold,
        breakeven_r=_breakeven_correlation(targets, pooling, options),
        rule_variants_df=variants,
        cap_hit_shares=inputs.cap_hit or None,
        degenerate_shares=inputs.degenerate,
        with_mu_uncertainty=options.with_mu_uncertainty,
        clip_lambda_at_zero=options.clip,
        warnings=tuple(messages),
        summary_df=summary,
        data_health_certificate_df=health,
    )


# ---- Pre-registerable retention rule -------------------------------------------------------


@dataclass(frozen=True, kw_only=True)
class RetentionRuleSpec:
    """A pre-registerable retention rule: when may the posterior-``rho`` regime P enter the headline?

    The literal reading only asks that the pooled ``rho`` P10 exceed ``p10_threshold``; the robust
    reading also asks every enabled clause below. Declare the spec before seeing the data; its
    fields are keyword-only.

    Attributes:
      p10_threshold: The pooled ``rho`` P10 must exceed this value (strictly).
      require_corr_adjusted: Also at the estimated between-outcome correlation (every
        ``correlation_matrix`` / ``correlation_common`` variant); fails when no correlation was
        supplied to :func:`pool_multi_outcome_retention`.
      require_leave_one_out: Also after dropping any single pooling outcome.
      max_single_outcome_weight: Largest admissible DerSimonian-Laird weight share of one outcome
        (``None`` disables the clause). With ``M`` pooling outcomes the largest share is at least
        ``1 / M``, so the clause cannot pass when ``M < 1 / max_single_outcome_weight``.
      max_cap_hit_share: Largest admissible share of draws at the upper ``rho`` bound, per outcome
        (pooling and focal), from :class:`PostPulseRetentionResult` inputs; fails when those shares
        are unavailable (``None`` disables the clause).
      max_degenerate_share: Largest admissible share of draws whose ``rho`` is unidentified
        (degenerate windows), per outcome (pooling and focal); above it the moments and P10s
        describe a selected subset of the draws (``None`` disables the clause).
      max_q_lower_tail_p: The clause fails when Cochran's ``Q`` lower-tail p-value falls below this
        level, i.e. when the outcomes are under-dispersed as correlated outcomes are (``None``
        disables the under-dispersion check).
    """

    p10_threshold: float = 0.90
    require_corr_adjusted: bool = True
    require_leave_one_out: bool = True
    max_single_outcome_weight: float | None = 0.40
    max_cap_hit_share: float | None = 0.25
    max_degenerate_share: float | None = 0.05
    max_q_lower_tail_p: float | None = 0.05

    def __post_init__(self) -> None:
        """Validates the thresholds."""
        threshold = _as_float(self.p10_threshold, "p10_threshold")
        if not 0.0 < threshold <= 1.0:
            raise ValueError(f"p10_threshold must lie in (0, 1]; got {self.p10_threshold!r}")
        if not isinstance(self.require_corr_adjusted, bool):
            raise ValueError(f"require_corr_adjusted must be a bool; got {self.require_corr_adjusted!r}")
        if not isinstance(self.require_leave_one_out, bool):
            raise ValueError(f"require_leave_one_out must be a bool; got {self.require_leave_one_out!r}")
        if self.max_single_outcome_weight is not None:
            weight = _as_float(self.max_single_outcome_weight, "max_single_outcome_weight")
            if not 0.0 < weight <= 1.0:
                raise ValueError(f"max_single_outcome_weight must lie in (0, 1]; got {weight!r}")
        if self.max_cap_hit_share is not None:
            share = _as_float(self.max_cap_hit_share, "max_cap_hit_share")
            if not 0.0 <= share <= 1.0:
                raise ValueError(f"max_cap_hit_share must lie in [0, 1]; got {self.max_cap_hit_share!r}")
        if self.max_degenerate_share is not None:
            share = _as_float(self.max_degenerate_share, "max_degenerate_share")
            if not 0.0 <= share <= 1.0:
                raise ValueError(f"max_degenerate_share must lie in [0, 1]; got {self.max_degenerate_share!r}")
        if self.max_q_lower_tail_p is not None:
            _check_open_unit(self.max_q_lower_tail_p, "max_q_lower_tail_p")


@dataclass(frozen=True, eq=False)
class RetentionRuleEvaluation:
    """Output of :func:`evaluate_retention_rule`.

    Attributes:
      spec: The evaluated :class:`RetentionRuleSpec`.
      rule_outcomes: Outcomes whose pooled ``rho`` P10 forms the rule statistic.
      literal_p10: Rule statistic of the literal (``r = 0``) pooling.
      literal_passed: Whether ``literal_p10 > spec.p10_threshold``.
      robust_passed: Whether the literal clause and every enabled robustness clause pass.
      fallback_recommended: ``not robust_passed``: report the bound regimes U and L only.
      reasons: One message per failed clause.
      headline_regimes: ``("P", "U", "L")`` if ``robust_passed`` else ``("U", "L")``.
      warnings: Configuration warnings: a spec threshold that differs from the pooling's, or a
        weight clause that cannot pass with this many pooling outcomes.
      summary_df: One row per clause with ``clause``, ``enabled``, ``value``, ``threshold``,
        ``passed`` and ``detail``.
    """

    spec: RetentionRuleSpec
    rule_outcomes: tuple[str, ...]
    literal_p10: float
    literal_passed: bool
    robust_passed: bool
    fallback_recommended: bool
    reasons: tuple[str, ...]
    headline_regimes: tuple[str, ...]
    warnings: tuple[str, ...]
    summary_df: pd.DataFrame


def _clause(
    clause: str, enabled: bool, value: float, threshold: float, passed: bool | None, detail: str
) -> dict[str, Any]:
    """One row of the rule-evaluation table (``passed`` is ``None`` for disabled clauses)."""
    return {
        "clause": clause,
        "enabled": enabled,
        "value": value,
        "threshold": threshold,
        "passed": passed if enabled else None,
        "detail": detail,
    }


def _minimum_variant(variants: pd.DataFrame, names: Sequence[str]) -> pd.Series | None:
    """Row with the smallest ``rho_P10`` among the given variants, or ``None`` if there is none."""
    rows = variants.loc[variants["variant"].isin(list(names))]
    return None if rows.empty else rows.loc[rows["rho_P10"].idxmin()]


def _weight_clause_unreachable(n_pooling: int, weight_limit: float) -> bool:
    """Whether no weighting of ``n_pooling`` outcomes can keep the largest share within ``weight_limit``."""
    return 1.0 / n_pooling > weight_limit


def _rule_clauses(pooled: PooledRetentionResult, spec: RetentionRuleSpec) -> list[dict[str, Any]]:
    """Clause rows: literal P10, correlation-adjusted, leave-one-out, weight, cap-hit, degenerate and Q tail."""
    threshold = float(spec.p10_threshold)
    variants = pooled.rule_variants_df
    who = " / ".join(pooled.rule_outcomes) if len(pooled.rule_outcomes) <= 3 else "min over outcomes"
    clauses = [_clause(
        "literal_p10", True, pooled.literal_p10, threshold, pooled.literal_p10 > threshold,
        f"literal pooled rho P10 ({who}) = {pooled.literal_p10:.4f}; must exceed {threshold:g}",
    )]
    adjusted = _minimum_variant(variants, _CORRELATION_VARIANTS)
    if adjusted is None:
        clauses.append(_clause(
            "correlation_adjusted_p10", spec.require_corr_adjusted, float("nan"), threshold, False,
            "no between-outcome correlation was supplied to pool_multi_outcome_retention, so robustness to "
            "correlated outcomes cannot be verified",
        ))
    else:
        clauses.append(_clause(
            "correlation_adjusted_p10", spec.require_corr_adjusted, float(adjusted["rho_P10"]), threshold,
            bool(adjusted["rho_P10"] > threshold),
            f"pooled rho P10 at the estimated correlation ({adjusted['variant']}, r = {adjusted['correlation']:.3f}) "
            f"= {adjusted['rho_P10']:.4f}; must exceed {threshold:g}",
        ))
    loo = _minimum_variant(variants, ("leave_one_out",))
    clauses.append(_clause(
        "leave_one_out_p10", spec.require_leave_one_out, float(loo["rho_P10"]), threshold,
        bool(loo["rho_P10"] > threshold),
        f"smallest leave-one-out pooled rho P10 = {loo['rho_P10']:.4f} (without {loo['left_out']!r}); must exceed "
        f"{threshold:g}",
    ))
    weight_limit = spec.max_single_outcome_weight
    weight_detail = f"{pooled.dominant_outcome!r} carries {pooled.dominant_weight:.1%} of the DerSimonian-Laird weight"
    if weight_limit is not None:
        weight_detail += f"; at most {weight_limit:.0%} allowed"
        if _weight_clause_unreachable(len(pooled.outcomes), weight_limit):
            weight_detail += (
                f" (unreachable: with {len(pooled.outcomes)} pooling outcomes the largest share is at least "
                f"{1.0 / len(pooled.outcomes):.0%})"
            )
    clauses.append(_clause(
        "max_single_outcome_weight", weight_limit is not None, pooled.dominant_weight,
        float("nan") if weight_limit is None else float(weight_limit),
        weight_limit is not None and pooled.dominant_weight <= weight_limit, weight_detail,
    ))
    cap_limit = spec.max_cap_hit_share
    needed = list(pooled.outcomes) + ([] if pooled.focal_label is None else [pooled.focal_label])
    shares = pooled.cap_hit_shares or {}
    missing = [name for name in needed if name not in shares]
    if missing:
        clauses.append(_clause(
            "max_cap_hit_share", cap_limit is not None, float("nan"),
            float("nan") if cap_limit is None else float(cap_limit), False,
            f"cap-hit shares unavailable for {missing}; pass PostPulseRetentionResult objects to "
            "pool_multi_outcome_retention",
        ))
    else:
        worst = max(needed, key=shares.__getitem__)
        clauses.append(_clause(
            "max_cap_hit_share", cap_limit is not None, shares[worst],
            float("nan") if cap_limit is None else float(cap_limit),
            cap_limit is not None and shares[worst] <= cap_limit,
            f"{shares[worst]:.1%} of the {worst!r} draws sit at the upper rho bound"
            + ("" if cap_limit is None else f"; at most {cap_limit:.0%} allowed"),
        ))
    degenerate_limit = spec.max_degenerate_share
    degenerate = pooled.degenerate_shares
    worst_degenerate = max(needed, key=degenerate.__getitem__)
    clauses.append(_clause(
        "max_degenerate_share", degenerate_limit is not None, degenerate[worst_degenerate],
        float("nan") if degenerate_limit is None else float(degenerate_limit),
        degenerate_limit is not None and degenerate[worst_degenerate] <= degenerate_limit,
        f"{degenerate[worst_degenerate]:.1%} of the {worst_degenerate!r} draws have an unidentified rho (degenerate "
        "window); the moments and P10s use the identified draws only"
        + ("" if degenerate_limit is None else f"; at most {degenerate_limit:.0%} allowed"),
    ))
    q_limit = spec.max_q_lower_tail_p
    p_lower = pooled.q_p_value_lower
    clauses.append(_clause(
        "q_lower_tail", q_limit is not None, p_lower, float("nan") if q_limit is None else float(q_limit),
        q_limit is not None and not p_lower < q_limit,
        f"Cochran's Q = {pooled.q_statistic:.3g} on {pooled.q_dof} df, lower-tail p = {p_lower:.3g}"
        + ("" if q_limit is None else f"; below {q_limit:g} the outcomes are not independent replicates"),
    ))
    return clauses


def evaluate_retention_rule(
    pooled: PooledRetentionResult, spec: RetentionRuleSpec | None = None
) -> RetentionRuleEvaluation:
    """Evaluate a pre-registered :class:`RetentionRuleSpec` on a pooled retention result.

    ``literal_passed`` is the rule as written (``literal_p10 > p10_threshold``); ``robust_passed``
    also requires every enabled clause of the spec, and ``fallback_recommended = not
    robust_passed``. The rule is the user's declared decision rule, evaluated literally and under
    the robustness variants; no calibration of the P10 threshold is claimed. A spec threshold that
    differs from the one the pooling used, or a weight clause that cannot pass with this many
    pooling outcomes, is reported as a :class:`StabilityWarning`.

    Args:
      pooled: Output of :func:`pool_multi_outcome_retention`.
      spec: The pre-registered rule (default: :class:`RetentionRuleSpec` defaults).

    Returns:
      :class:`RetentionRuleEvaluation`.

    Raises:
      ValueError: If ``pooled`` or ``spec`` has the wrong type.
    """
    if not isinstance(pooled, PooledRetentionResult):
        raise ValueError(f"pooled must be a PooledRetentionResult; got {type(pooled).__name__}")
    rule = RetentionRuleSpec() if spec is None else spec
    if not isinstance(rule, RetentionRuleSpec):
        raise ValueError(f"spec must be a RetentionRuleSpec or None; got {type(rule).__name__}")
    clauses = _rule_clauses(pooled, rule)
    literal_passed = bool(clauses[0]["passed"])
    robust_passed = all(bool(clause["passed"]) for clause in clauses if clause["enabled"])
    failed = [clause for clause in clauses if clause["enabled"] and not clause["passed"]]
    reasons = tuple(f"{clause['clause']}: {clause['detail']}" for clause in failed)
    messages = []
    if float(rule.p10_threshold) != pooled.p10_threshold:
        messages.append(
            f"spec.p10_threshold = {float(rule.p10_threshold):g} differs from the pooling's p10_threshold = "
            f"{pooled.p10_threshold:g}: the clauses use the spec's, but breakeven_r and the passed column of "
            "rule_variants_df use the pooling's; pass p10_threshold=spec.p10_threshold to pool_multi_outcome_retention."
        )
    weight_limit = rule.max_single_outcome_weight
    if weight_limit is not None and _weight_clause_unreachable(len(pooled.outcomes), weight_limit):
        messages.append(
            f"the max_single_outcome_weight clause cannot pass: with {len(pooled.outcomes)} pooling outcomes the "
            f"largest weight share is at least {1.0 / len(pooled.outcomes):.0%}, above the {weight_limit:.0%} allowed."
        )
    for message in messages:
        warnings.warn(message, StabilityWarning, stacklevel=2)
    return RetentionRuleEvaluation(
        spec=rule,
        rule_outcomes=pooled.rule_outcomes,
        literal_p10=pooled.literal_p10,
        literal_passed=literal_passed,
        robust_passed=robust_passed,
        fallback_recommended=not robust_passed,
        reasons=reasons,
        headline_regimes=PERSISTENCE_REGIMES if robust_passed else BOUND_REGIMES,
        warnings=tuple(messages),
        summary_df=pd.DataFrame(clauses),
    )


# ---- Orchestration: trajectory stability report --------------------------------------------


@dataclass(frozen=True, eq=False)
class TrajectoryStabilityReport:
    """Output of :func:`diagnose_trajectory_stability`; ``verdict`` is ``None`` by design.

    The report states what the observed window can and cannot say about persistence; it does not
    certify a long-horizon extrapolation (see the module docstring and :func:`~longbet.att_stability`).

    Attributes:
      outcomes: Every outcome analysed (pooling outcomes, then the focal outcome).
      focal: Focal outcome label, or ``None``.
      window: ``(s0, s1)`` post-pulse window of the retention fits.
      retention: ``{outcome: PostPulseRetentionResult}``.
      pooled: The :class:`PooledRetentionResult`.
      rule: The :class:`RetentionRuleEvaluation` of the pre-registered rule.
      wearout: ``{name: WearoutAncovaResult}`` supplied by the caller (may be empty).
      slope_correlation: The :class:`SlopeCorrelationResult` behind the correlation, or ``None``.
      unobserved_effect_share: ``TransportAuditResult.unobserved_effect_share`` when a transport
        audit was supplied (share of the cumulative effect carried by unobserved cells), else ``None``.
      unobserved_cell_share: ``TransportAuditResult.unobserved_cell_share``, or ``None``.
      headline_regimes: ``rule.headline_regimes``.
      warnings: Wear-out tests that cannot speak to persistence, boundary pile-up, degenerate
        windows, a dominant outcome, the pooling and rule warnings and the rule fallback.
      summary_df: One row per outcome (pooling summary plus the retention boundary shares).
      data_health_certificate_df: The components' certificates stacked, with a ``component`` column.
      verdict: Always ``None``; any other value raises ``ValueError``.
    """

    outcomes: tuple[str, ...]
    focal: str | None
    window: tuple[int, int]
    retention: dict[str, PostPulseRetentionResult]
    pooled: PooledRetentionResult
    rule: RetentionRuleEvaluation
    wearout: dict[str, WearoutAncovaResult]
    slope_correlation: SlopeCorrelationResult | None
    unobserved_effect_share: dict[Any, float] | None
    unobserved_cell_share: dict[Any, float] | None
    headline_regimes: tuple[str, ...]
    warnings: tuple[str, ...]
    summary_df: pd.DataFrame
    data_health_certificate_df: pd.DataFrame
    verdict: None = None

    def __post_init__(self) -> None:
        """Rejects a verdict: the report withholds one by design."""
        if self.verdict is not None:
            raise ValueError(f"TrajectoryStabilityReport withholds a verdict by design; got verdict={self.verdict!r}")


def _wearout_mapping(wearout: Any) -> dict[str, WearoutAncovaResult]:
    """Normalize ``None`` / one result / a mapping of results to ``{name: WearoutAncovaResult}``."""
    if wearout is None:
        return {}
    if isinstance(wearout, WearoutAncovaResult):
        return {"wearout": wearout}
    if isinstance(wearout, Mapping) and all(isinstance(value, WearoutAncovaResult) for value in wearout.values()):
        return {str(key): value for key, value in wearout.items()}
    raise ValueError("wearout must be None, a WearoutAncovaResult or a mapping of them")


def _wearout_message(name: str, result: WearoutAncovaResult) -> str | None:
    """Why a wear-out test cannot speak to persistence, or ``None`` when it is powered for half wear-out."""
    if result.powered_for_half_wearout:
        return None
    if result.early_estimate <= 0.0:
        return (
            f"wear-out test {name!r} has a non-positive early-level estimate ({result.early_estimate:.4g}): there is "
            f"no positive effect to wear out, so its change (MDD {result.change_mdd:.4g}) says nothing about "
            "persistence."
        )
    if not result.change_mdd > 0.0:
        return (
            f"wear-out test {name!r} has a zero standard error for the late-minus-early change (the ANCOVA fits it "
            "exactly), so its power is not assessed."
        )
    return (
        f"wear-out test {name!r} is underpowered: the MDD of the late-minus-early change ({result.change_mdd:.4g}) "
        f"is {result.mdd_pct_of_early:.0f}% of the early-level effect ({result.early_estimate:.4g}), above the 50% "
        "needed to detect half wear-out; a null change is not evidence of persistence."
    )


def _stability_messages(
    retention: Mapping[str, PostPulseRetentionResult],
    pooled: PooledRetentionResult,
    rule: RetentionRuleEvaluation,
    wearouts: Mapping[str, WearoutAncovaResult],
) -> list[str]:
    """Report warnings in a fixed order: wear-out, pile-up, degenerate windows, weight, pooling, rule."""
    messages = []
    for name, result in wearouts.items():
        messages.extend(f"wear-out test {name!r}: {warning}" for warning in result.warnings)
        message = _wearout_message(name, result)
        if message is not None:
            messages.append(message)
    for name, result in retention.items():
        if result.boundary_pileup_flag:
            messages.append(
                f"boundary pile-up: {result.share_at_upper:.1%} of the {name!r} retention draws sit at the upper "
                f"bound rho = {result.rho_bounds[1]:g} (threshold {result.pileup_threshold:.0%}); bounded or floored "
                "effects read as precision there."
            )
        if result.share_degenerate > 0.0:
            messages.append(
                f"{result.share_degenerate:.1%} of the {name!r} draws have a degenerate window (rho unidentified; "
                f"policy {result.degenerate_policy!r})."
            )
    weight_limit = rule.spec.max_single_outcome_weight
    if weight_limit is not None and pooled.dominant_weight > weight_limit:
        messages.append(
            f"dominant outcome: {pooled.dominant_outcome!r} carries {pooled.dominant_weight:.1%} of the "
            f"DerSimonian-Laird weight (more than {weight_limit:.0%}); the pooled rho is mostly that outcome's."
        )
    messages.extend(pooled.warnings)
    messages.extend(rule.warnings)
    if rule.fallback_recommended:
        messages.append(
            "retention rule fallback: report the bound regimes U and L only (" + "; ".join(rule.reasons) + ")."
        )
    return messages


def diagnose_trajectory_stability(
    traj_by_outcome: Mapping[str, Any],
    *,
    window: tuple[int, int],
    focal: str | None = None,
    pooling_outcomes: Sequence[str] | None = None,
    support: Any = None,
    contrast: str = "AC",
    weighting: str = "unit",
    correlation: Any = None,
    wearout: Any = None,
    transport_audit: TransportAuditResult | None = None,
    rule: RetentionRuleSpec | None = None,
    rho_bounds: tuple[float, float] = (0.50, 1.00),
    unconstrained_upper: float = 1.50,
    grid_points: int = 1001,
    degenerate_policy: str = "nan",
    pileup_threshold: float = 0.25,
    seed: int = 2026,
    with_mu_uncertainty: bool = True,
    clip_lambda_at_zero: bool = True,
    quantiles: Sequence[float] = RETENTION_QUANTILES,
) -> TrajectoryStabilityReport:
    """In-sample stability report of a 2- or 3-arm Rank-K analysis; ``verdict`` is ``None`` by design.

    Runs :func:`estimate_post_pulse_retention` on every outcome, pools the pooling outcomes with
    :func:`pool_multi_outcome_retention` (shrinking ``focal`` toward them), evaluates the
    pre-registered ``rule`` with :func:`evaluate_retention_rule`, and collects the panel-level
    evidence the caller supplies: :func:`wearout_ancova` results, the
    :func:`unit_slope_correlation` behind ``correlation`` and a
    :func:`~longbet.rank_k_transport_audit` whose ``unobserved_effect_share`` is carried into the
    report. Every warning is also emitted as a :class:`StabilityWarning`.

    Args:
      traj_by_outcome: ``{outcome: trajectory}``; each trajectory is anything
        :func:`estimate_post_pulse_retention` accepts (``(S, L)`` draws, a fit, a 3-arm
        reconciliation or a ``{contrast: (S, L)}`` mapping).
      window: ``(s0, s1)`` post-pulse window of every retention fit.
      focal: Optional outcome shrunk toward the pool and excluded from it (the outcome a decision
        rests on).
      pooling_outcomes: Outcomes pooled (default: every outcome but ``focal``; at least two).
      support: ``None``, one weight vector shared by every outcome, or ``{outcome: weights}`` keyed by
        ``traj_by_outcome`` labels (an outcome without an entry uses the default support; an unknown
        key raises).
      contrast: Contrast of reconciliation / mapping inputs.
      weighting: ``"unit"`` or ``"cell"`` for fit inputs.
      correlation: Between-outcome correlation of the pooling outcomes: ``None``, a float, an
        ``(M, M)`` matrix or a :class:`SlopeCorrelationResult`.
      wearout: ``None``, a :class:`WearoutAncovaResult` or ``{name: WearoutAncovaResult}`` (for
        instance one per contrast of a 3-arm design).
      transport_audit: Optional :class:`~longbet.TransportAuditResult` of a masked fit.
      rule: The pre-registered :class:`RetentionRuleSpec` (default: its defaults).
      rho_bounds: Passed to :func:`estimate_post_pulse_retention`.
      unconstrained_upper: Passed to :func:`estimate_post_pulse_retention`.
      grid_points: Passed to :func:`estimate_post_pulse_retention`.
      degenerate_policy: Passed to :func:`estimate_post_pulse_retention`.
      pileup_threshold: Passed to :func:`estimate_post_pulse_retention`.
      seed: Passed to :func:`estimate_post_pulse_retention`.
      with_mu_uncertainty: Passed to :func:`pool_multi_outcome_retention`.
      clip_lambda_at_zero: Passed to :func:`pool_multi_outcome_retention`.
      quantiles: Probabilities of the summaries (odd length; default P10 / P50 / P90).

    Returns:
      :class:`TrajectoryStabilityReport`.

    Raises:
      ValueError: On malformed inputs (see the component functions).
    """
    if not isinstance(traj_by_outcome, Mapping) or not traj_by_outcome:
        raise ValueError("traj_by_outcome must be a non-empty mapping {outcome: trajectory}")
    by_name = {str(key): value for key, value in traj_by_outcome.items()}
    if len(by_name) != len(traj_by_outcome):
        raise ValueError("traj_by_outcome keys must stay distinct after conversion to str")
    focal_name = None if focal is None else str(focal)
    if focal_name is not None and focal_name not in by_name:
        raise ValueError(f"focal {focal!r} is not a key of traj_by_outcome")
    if pooling_outcomes is None:
        pool_names = [name for name in by_name if name != focal_name]
    else:
        pool_names = [str(name) for name in pooling_outcomes]
    unknown = [name for name in pool_names if name not in by_name]
    if unknown or len(set(pool_names)) != len(pool_names) or focal_name in pool_names or len(pool_names) < 2:
        raise ValueError(
            f"pooling_outcomes must be at least two distinct keys of traj_by_outcome other than focal; got "
            f"{pool_names} (unknown: {unknown})"
        )
    spec = RetentionRuleSpec() if rule is None else rule
    if not isinstance(spec, RetentionRuleSpec):
        raise ValueError(f"rule must be a RetentionRuleSpec or None; got {type(spec).__name__}")
    wearouts = _wearout_mapping(wearout)
    if transport_audit is not None and not isinstance(transport_audit, TransportAuditResult):
        raise ValueError(f"transport_audit must be a TransportAuditResult; got {type(transport_audit).__name__}")
    names = pool_names + ([] if focal_name is None else [focal_name])
    if isinstance(support, Mapping):
        supports = {str(key): value for key, value in support.items()}
        unknown_support = sorted(key for key in supports if key not in by_name)
        if unknown_support or len(supports) != len(support):
            raise ValueError(
                f"support keys must be distinct keys of traj_by_outcome {sorted(by_name)}; got unknown keys "
                f"{unknown_support}"
            )
    else:
        supports = dict.fromkeys(names, support)
    retention: dict[str, PostPulseRetentionResult] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", StabilityWarning)
        for name in names:
            retention[name] = estimate_post_pulse_retention(
                by_name[name], window=window, support=supports.get(name), contrast=contrast, weighting=weighting,
                rho_bounds=rho_bounds, unconstrained_upper=unconstrained_upper, grid_points=grid_points,
                degenerate_policy=degenerate_policy, pileup_threshold=pileup_threshold, seed=seed,
                quantiles=quantiles, label=name,
            )
        pooled = pool_multi_outcome_retention(
            {name: retention[name] for name in pool_names},
            focal=None if focal_name is None else retention[focal_name],
            focal_label=focal_name,
            correlation=correlation,
            with_mu_uncertainty=with_mu_uncertainty,
            clip_lambda_at_zero=clip_lambda_at_zero,
            p10_threshold=spec.p10_threshold,
            quantiles=quantiles,
        )
        evaluation = evaluate_retention_rule(pooled, spec)
    messages = list(dict.fromkeys(_stability_messages(retention, pooled, evaluation, wearouts)))
    for message in messages:
        warnings.warn(message, StabilityWarning, stacklevel=2)

    boundary = pd.DataFrame([
        {
            "outcome": name,
            "share_at_upper": result.share_at_upper,
            "share_at_lower": result.share_at_lower,
            "share_degenerate": result.share_degenerate,
            "share_unconstrained_above_one": result.share_unconstrained_above_one,
            "boundary_pileup_flag": result.boundary_pileup_flag,
        }
        for name, result in retention.items()
    ])
    summary = pooled.summary_df.merge(boundary, on="outcome", how="left")
    parts = []
    for name, result in retention.items():
        parts.append(result.data_health_certificate_df.assign(component=f"retention[{name}]"))
    parts.append(pooled.data_health_certificate_df.assign(component="pooling"))
    for name, result in wearouts.items():
        parts.append(result.data_health_certificate_df.assign(component=f"wearout[{name}]"))
    slope = correlation if isinstance(correlation, SlopeCorrelationResult) else None
    if slope is not None:
        parts.append(slope.data_health_certificate_df.assign(component="slope_correlation"))
    health = pd.concat(parts, ignore_index=True)[["component", *_HEALTH_COLUMNS]]
    return TrajectoryStabilityReport(
        outcomes=tuple(names),
        focal=focal_name,
        window=retention[names[0]].window,
        retention=retention,
        pooled=pooled,
        rule=evaluation,
        wearout=wearouts,
        slope_correlation=slope,
        unobserved_effect_share=None if transport_audit is None else dict(transport_audit.unobserved_effect_share),
        unobserved_cell_share=None if transport_audit is None else dict(transport_audit.unobserved_cell_share),
        headline_regimes=evaluation.headline_regimes,
        warnings=tuple(messages),
        summary_df=summary,
        data_health_certificate_df=health,
        verdict=None,
    )


# ---- Sizing a leak-free hindcast before running it -----------------------------------------


@dataclass(frozen=True, eq=False)
class PseudoHorizonSeparationResult:
    """Output of :func:`pseudo_horizon_separation`.

    Attributes:
      pseudo_horizon: Anchor event time ``L0`` of the pseudo-horizon forecasts.
      observed_horizon: Observed horizon ``L``; the pseudo window is ``L0 + 1..L``.
      horizon: Extrapolation horizon ``H`` of the real tail ``L + 1..H``.
      rho_bounds: ``(rho_L, rho_U)``: Bound L and Bound U retention.
      pulse_decay_tau: Pulse decay time of the plateau / pulse split, or ``None`` (whole level).
      alpha: One-sided level of the power calculation.
      power: Power of the minimum detectable separation.
      forecast_sum_draws: ``{"U": (S,), "L": (S,)}`` forecast sums over the pseudo window.
      gap_draws: ``(S,)`` per-draw gap ``U - L`` between the cumulative forecasts.
      delta: ``mean(U sum) - mean(L sum)``, the separation the hindcast must detect.
      observed_window_sum_draws: ``(S,)`` observed cumulative over the pseudo window, or ``None``
        when the draws stop before ``L``.
      observed_window_sd: Its posterior SD, or ``None``.
      min_detectable_gap: Gap threshold of ``share_gap_above_min`` (default: the first target's MDD).
      share_gap_above_min: Share of draws whose gap exceeds ``min_detectable_gap``.
      lever_ratio: ``sum_{h <= L - L0} w_h / sum_{h <= H - L} w_h`` with the scenario-gap weights
        ``w_h = (rho_U^h - rho_L^h) sigma_h`` (``sigma_h = 1`` in the observed pseudo window and
        ``survival[L + h - 1]`` in the tail).
      tail_lever_multiple: ``1 / lever_ratio``: how many times larger the real tail lever is.
      lever_df: One row per step with ``segment`` (``pseudo`` / ``tail``), ``h``, ``gap_weight``,
        ``survival`` and ``weight``.
      power_df: One row per target with its SD (and mean), the Bound L forecast SD, the combined
        SD, ``delta``, ``z_if_persistence_true``, ``power`` and ``mdd``, plus the target's z-scores
        against both forecasts when its mean is known.
      summary_df: One row with the horizons, bounds, ``delta``, gap quantiles, lever and power range.
      data_health_certificate_df: Row counts and the distribution of the anchors, forecasts and gap.
    """

    pseudo_horizon: int
    observed_horizon: int
    horizon: int
    rho_bounds: tuple[float, float]
    pulse_decay_tau: float | None
    alpha: float
    power: float
    forecast_sum_draws: dict[str, np.ndarray]
    gap_draws: np.ndarray
    delta: float
    observed_window_sum_draws: np.ndarray | None
    observed_window_sd: float | None
    min_detectable_gap: float
    share_gap_above_min: float
    lever_ratio: float
    tail_lever_multiple: float
    lever_df: pd.DataFrame
    power_df: pd.DataFrame
    summary_df: pd.DataFrame
    data_health_certificate_df: pd.DataFrame


def _resolve_arms(traj_draws: Any, *, split: bool, arm_names: tuple[str, str, str]) -> dict[str, np.ndarray]:
    """``{"total": (S, T)}``, plus ``"plateau"`` and ``"pulse"`` under the split, from any trajectory input."""
    plateau_arm, pulse_arm, total_arm = arm_names
    if isinstance(traj_draws, Nested3ArmReconciliationResult):
        source: Mapping[str, Any] = {
            "BC": traj_draws.trt_catt_bc, "AB": traj_draws.trt_catt_ab, "AC": traj_draws.trt_catt_ac,
        }
    elif isinstance(traj_draws, Mapping):
        source = traj_draws
    elif split:
        raise ValueError(
            f"the plateau / pulse split (pulse_decay_tau) needs {{{plateau_arm!r}, {pulse_arm!r}}} trajectories: pass "
            "a mapping or a Nested3ArmReconciliationResult"
        )
    else:
        return {"total": _as_draw_matrix(traj_draws, "traj_draws")}
    if not split:
        if total_arm not in source:
            raise ValueError(f"traj_draws has no total arm {total_arm!r}; got {sorted(map(str, source))}")
        return {"total": _as_draw_matrix(source[total_arm], f"traj_draws[{total_arm!r}]")}
    missing = [arm for arm in (plateau_arm, pulse_arm) if arm not in source]
    if missing:
        raise ValueError(f"traj_draws lacks the arm(s) {missing} of the plateau / pulse split")
    arms = {
        "plateau": _as_draw_matrix(source[plateau_arm], f"traj_draws[{plateau_arm!r}]"),
        "pulse": _as_draw_matrix(source[pulse_arm], f"traj_draws[{pulse_arm!r}]"),
    }
    if total_arm in source:
        arms["total"] = _as_draw_matrix(source[total_arm], f"traj_draws[{total_arm!r}]")
    else:
        arms["total"] = arms["plateau"] + arms["pulse"]
    if len({arr.shape for arr in arms.values()}) != 1:
        raise ValueError(f"the arms must share their (S, T) shape; got {[arr.shape for arr in arms.values()]}")
    scale = max(1.0, float(np.max(np.abs(arms["total"]))))
    gap = float(np.max(np.abs(arms["total"] - (arms["plateau"] + arms["pulse"]))))
    if gap > 1e-8 * scale:
        raise ValueError(
            f"traj_draws[{total_arm!r}] must equal traj_draws[{plateau_arm!r}] + traj_draws[{pulse_arm!r}] under the "
            f"plateau / pulse split (max |total - (plateau + pulse)| = {gap:.3e})"
        )
    return arms


def _parse_targets(targets: Mapping[str, Any] | None) -> list[tuple[str, float, float]]:
    """``[(label, mean or NaN, sd)]`` from ``{label: sd}`` or ``{label: (mean, sd)}``."""
    if targets is None:
        return []
    if not isinstance(targets, Mapping):
        raise ValueError("targets must be a mapping {label: sd} or {label: (mean, sd)}")
    parsed = []
    for label, spec in targets.items():
        if isinstance(spec, (tuple, list)):
            if len(spec) != 2:
                raise ValueError(f"targets[{label!r}] must be an SD or a (mean, sd) pair; got {spec!r}")
            mean, sd = _as_float(spec[0], f"targets[{label!r}] mean"), _as_float(spec[1], f"targets[{label!r}] sd")
        else:
            mean, sd = float("nan"), _as_float(spec, f"targets[{label!r}]")
        if sd < 0.0:
            raise ValueError(f"targets[{label!r}] SD must be >= 0; got {sd}")
        parsed.append((str(label), mean, sd))
    return parsed


def _survival_weights(survival: Any, horizon: int, observed_horizon: int) -> np.ndarray:
    """Validate ``(H,)`` survival weights: finite, in ``[0, 1]``, exactly 1 on ``t <= L`` (``None``: all ones)."""
    if survival is None:
        return np.ones(horizon)
    weights = _as_vector(survival, "survival", horizon)
    if np.any(weights < 0.0) or np.any(weights > 1.0):
        raise ValueError("survival weights must lie in [0, 1]")
    if np.any(weights[:observed_horizon] != 1.0):
        raise ValueError(
            "survival must be 1 on t <= observed_horizon and S(t) / S(L) beyond (the trajectory-layer convention)"
        )
    return weights


def _check_scenario_bounds(rho_bounds: Any) -> tuple[float, float]:
    """Validate scenario retentions ``(rho_L, rho_U)`` with ``0 <= rho_L < rho_U <= 1``."""
    try:
        low_raw, high_raw = rho_bounds
    except (TypeError, ValueError) as exc:
        raise ValueError(f"rho_bounds must be a (rho_L, rho_U) pair; got {rho_bounds!r}") from exc
    low, high = _as_float(low_raw, "rho_bounds[0]"), _as_float(high_raw, "rho_bounds[1]")
    if not 0.0 <= low < high <= 1.0:
        raise ValueError(f"rho_bounds must satisfy 0 <= rho_L < rho_U <= 1; got {rho_bounds!r}")
    return low, high


def _power_rows(
    sums: Mapping[str, np.ndarray],
    delta: float,
    targets: Sequence[tuple[str, float, float]],
    z_values: tuple[float, float],
) -> list[dict[str, Any]]:
    """Power to separate U from L per target, and each target's z-score against both forecasts."""
    z_alpha, z_power = z_values
    var_fast = float(np.var(sums["L"], ddof=1))
    rows = []
    for label, mean, sd in targets:
        combined_sd = float(np.sqrt(var_fast + sd**2))
        z_expected = delta / combined_sd if combined_sd > 0.0 else float("inf")
        row: dict[str, Any] = {
            "target": label,
            "target_mean": mean,
            "target_sd": sd,
            "forecast_sd_bound_l": math.sqrt(var_fast),
            "combined_sd": combined_sd,
            "delta": delta,
            "z_if_persistence_true": z_expected,
            "power": float(sp_stats.norm.cdf(z_expected - z_alpha)),
            "mdd": (z_alpha + z_power) * combined_sd,
        }
        for regime in ("U", "L"):
            regime_sd = float(np.sqrt(np.var(sums[regime], ddof=1) + sd**2))
            z_target = float("nan")
            if math.isfinite(mean) and regime_sd > 0.0:
                z_target = (mean - float(np.mean(sums[regime]))) / regime_sd
            row[f"z_target_minus_{regime}"] = z_target
        rows.append(row)
    return rows


def pseudo_horizon_separation(
    traj_draws: Any,
    *,
    pseudo_horizon: int,
    observed_horizon: int,
    horizon: int,
    rho_bounds: tuple[float, float] = (0.85, 1.00),
    pulse_decay_tau: float | None = None,
    plateau_arm: str = "BC",
    pulse_arm: str = "AB",
    total_arm: str = "AC",
    survival: Any = None,
    targets: Mapping[str, Any] | None = None,
    min_detectable_gap: float | None = None,
    alpha: float = 0.05,
    power: float = 0.80,
    quantiles: Sequence[float] = RETENTION_QUANTILES,
) -> PseudoHorizonSeparationResult:
    """Can a leak-free hindcast through ``L0 < L`` tell no decay from decay? Sized from draws, before any refit.

    Anchored at event time ``L0`` the scenario forecasts over the pseudo window ``j = 1..L - L0``
    are ``tau_rho(L0 + j) = tau_plateau(L0) rho^j + tau_pulse(L0) exp(-j / pulse_decay_tau)`` with
    the plateau / pulse split, or ``tau(L0) rho^j`` for the whole level, for Bound U (``rho_U``) and
    Bound L (``rho_L``); the pulse cancels in their gap. With ``delta = mean(U sum) - mean(L sum)``
    and, per target, ``combined_sd = sqrt(Var(L sum) + target_sd^2)``, a one-sided level-``alpha``
    hindcast has ``power = Phi(delta / combined_sd - z_{1 - alpha})`` and detects
    ``mdd = (z_{1 - alpha} + z_power) combined_sd``. Targets are the observed cumulative over
    ``L0 + 1..L`` (its posterior SD from the same draws, when they reach ``L``) and any external
    target (a model-free ANCOVA SE, another posterior). The lever ratio compares the scenario-gap
    exposure ``w_h = (rho_U^h - rho_L^h) sigma_h`` of the pseudo window with that of the real tail
    ``L + 1..H`` (``sigma_h`` from ``survival``): a hindcast that cannot separate the scenarios
    over the short lever says little about the long one.

    Args:
      traj_draws: ``(S, T)`` total-effect draws, a ``{arm: (S, T)}`` mapping or a
        :class:`~longbet.Nested3ArmReconciliationResult`; ``T >= L0`` (``T >= L`` for the
        same-draws target).
      pseudo_horizon: Anchor event time ``L0 >= 1``.
      observed_horizon: Observed horizon ``L > L0``.
      horizon: Extrapolation horizon ``H > L`` of the tail lever.
      rho_bounds: ``(rho_L, rho_U)`` with ``0 <= rho_L < rho_U <= 1``.
      pulse_decay_tau: Pulse decay time (periods) of the plateau / pulse split, or ``None``.
      plateau_arm: Plateau arm key of the split.
      pulse_arm: Pulse arm key of the split.
      total_arm: Total arm key (``plateau + pulse`` when absent under the split).
      survival: Optional ``(H,)`` survival weights, 1 on ``t <= L`` and ``S(t) / S(L)`` beyond, not
        all zero on the tail; used by the lever only (the pseudo window is observed).
      targets: Optional external targets ``{label: sd}`` or ``{label: (mean, sd)}`` of the observed
        cumulative over ``L0 + 1..L``.
      min_detectable_gap: Gap threshold of ``share_gap_above_min`` (default: the first target's MDD).
      alpha: One-sided level, in ``(0, 1)``.
      power: Power of the MDD, in ``(0, 1)``.
      quantiles: Probabilities of the gap summary (odd length; default P10 / P50 / P90).

    Returns:
      :class:`PseudoHorizonSeparationResult`.

    Raises:
      ValueError: On malformed inputs, a survival that is zero over the whole tail, or when no
        target is available.
    """
    split = pulse_decay_tau is not None
    tau_pulse = None if pulse_decay_tau is None else _positive_scalar(pulse_decay_tau, "pulse_decay_tau")
    arms = _resolve_arms(traj_draws, split=split, arm_names=(plateau_arm, pulse_arm, total_arm))
    n_draws, n_periods = arms["total"].shape
    anchor = _as_int(pseudo_horizon, "pseudo_horizon", 1)
    observed = _as_int(observed_horizon, "observed_horizon", 1)
    extrapolated = _as_int(horizon, "horizon", 1)
    if not anchor < observed < extrapolated:
        raise ValueError(
            f"need pseudo_horizon < observed_horizon < horizon; got {anchor}, {observed}, {extrapolated}"
        )
    if anchor > n_periods:
        raise ValueError(f"pseudo_horizon ({anchor}) exceeds the {n_periods} periods of traj_draws")
    rho_low, rho_high = _check_scenario_bounds(rho_bounds)
    alpha_val, power_val = _check_open_unit(alpha, "alpha"), _check_open_unit(power, "power")
    qs = _check_quantiles(quantiles)
    survival_weights = _survival_weights(survival, extrapolated, observed)

    steps = np.arange(1, observed - anchor + 1, dtype=np.float64)
    base = arms["plateau" if split else "total"][:, [anchor - 1]]
    forecasts = {}
    for regime, rho in (("U", rho_high), ("L", rho_low)):
        forecast = base * rho ** steps[None, :]
        if split:
            forecast = forecast + arms["pulse"][:, [anchor - 1]] * np.exp(-steps / tau_pulse)[None, :]
        forecasts[regime] = forecast
    sums = {regime: forecast.sum(axis=1) for regime, forecast in forecasts.items()}
    gap = sums["U"] - sums["L"]
    delta = float(np.mean(sums["U"]) - np.mean(sums["L"]))
    window_sum, window_sd = None, None
    parsed = []
    if n_periods >= observed:
        window_sum = arms["total"][:, anchor:observed].sum(axis=1)
        window_sd = float(np.std(window_sum, ddof=1))
        parsed.append(("posterior SD of the observed cumulative (same draws)", float("nan"), window_sd))
    parsed.extend(_parse_targets(targets))
    if not parsed:
        raise ValueError(
            f"no target: traj_draws stops at {n_periods} < observed_horizon = {observed}; pass targets={{label: sd}}"
        )
    z_values = (float(sp_stats.norm.ppf(1.0 - alpha_val)), float(sp_stats.norm.ppf(power_val)))
    power_df = pd.DataFrame(_power_rows(sums, delta, parsed, z_values))
    if min_detectable_gap is None:
        threshold = float(power_df["mdd"].iloc[0])
    else:
        threshold = _as_float(min_detectable_gap, "min_detectable_gap")
    pseudo_steps = np.arange(1, observed - anchor + 1, dtype=np.float64)
    tail_steps = np.arange(1, extrapolated - observed + 1, dtype=np.float64)
    pseudo_gap = rho_high**pseudo_steps - rho_low**pseudo_steps
    tail_gap = rho_high**tail_steps - rho_low**tail_steps
    tail_survival = survival_weights[observed:extrapolated]
    tail_total = float(np.sum(tail_gap * tail_survival))
    if not tail_total > 0.0:
        raise ValueError("survival is zero over the whole tail L + 1..H, so the tail lever is zero")
    lever_ratio = float(np.sum(pseudo_gap)) / tail_total
    lever_df = pd.concat([
        pd.DataFrame({"segment": "pseudo", "h": pseudo_steps.astype(int), "gap_weight": pseudo_gap, "survival": 1.0,
                      "weight": pseudo_gap}),
        pd.DataFrame({"segment": "tail", "h": tail_steps.astype(int), "gap_weight": tail_gap, "survival": tail_survival,
                      "weight": tail_gap * tail_survival}),
    ], ignore_index=True)
    summary: dict[str, Any] = {
        "pseudo_horizon": anchor,
        "observed_horizon": observed,
        "horizon": extrapolated,
        "rho_L": rho_low,
        "rho_U": rho_high,
        "split": split,
        "delta": delta,
        "observed_window_sd": float("nan") if window_sd is None else window_sd,
        "power_min": float(power_df["power"].min()),
        "power_max": float(power_df["power"].max()),
        "min_detectable_gap": threshold,
        "share_gap_above_min": float(np.mean(gap > threshold)),
        "lever_ratio": lever_ratio,
        "tail_lever_multiple": 1.0 / lever_ratio,
    }
    summary.update(_summarize(gap, qs, "gap_"))
    health = [_count_row("posterior draws", n_draws), _count_row("periods in traj_draws", n_periods)]
    health.extend(_distribution_row(f"anchor tau_{arm}(L0={anchor})", arr[:, anchor - 1]) for arm, arr in arms.items())
    health.extend(_distribution_row(f"forecast sum {regime} (L0+1..L)", values) for regime, values in sums.items())
    health.append(_distribution_row("gap U - L", gap))
    if window_sum is not None:
        health.append(_distribution_row("observed cumulative (L0+1..L)", window_sum))
    return PseudoHorizonSeparationResult(
        pseudo_horizon=anchor,
        observed_horizon=observed,
        horizon=extrapolated,
        rho_bounds=(rho_low, rho_high),
        pulse_decay_tau=tau_pulse,
        alpha=alpha_val,
        power=power_val,
        forecast_sum_draws=sums,
        gap_draws=gap,
        delta=delta,
        observed_window_sum_draws=window_sum,
        observed_window_sd=window_sd,
        min_detectable_gap=threshold,
        share_gap_above_min=float(np.mean(gap > threshold)),
        lever_ratio=lever_ratio,
        tail_lever_multiple=1.0 / lever_ratio,
        lever_df=lever_df,
        power_df=power_df,
        summary_df=pd.DataFrame([summary]),
        data_health_certificate_df=_health_frame(health),
    )


__all__ = [
    "BOUND_REGIMES",
    "CORRELATION_TYPES",
    "DEGENERATE_POLICIES",
    "PERSISTENCE_REGIMES",
    "PooledRetentionResult",
    "PostPulseRetentionResult",
    "PseudoHorizonSeparationResult",
    "RETENTION_QUANTILES",
    "RetentionRuleEvaluation",
    "RetentionRuleSpec",
    "SLOPE_METHODS",
    "SlopeCorrelationResult",
    "StabilityWarning",
    "TrajectoryStabilityReport",
    "WearoutAncovaResult",
    "diagnose_trajectory_stability",
    "estimate_post_pulse_retention",
    "evaluate_retention_rule",
    "pool_multi_outcome_retention",
    "pseudo_horizon_separation",
    "unit_slope_correlation",
    "wearout_ancova",
]
