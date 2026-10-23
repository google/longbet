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

"""Population aggregates and the transport audit of Rank-K event-time fits.

All functions read the posterior-predictive cell effects ``res.catt_draws``
``(S, N, T)`` of a :class:`~longbet.RankKLongBetResult` in *event time*: for a
treated unit adopting in period ``g`` the event-time cell ``s = 1, ..., H`` is
panel column ``g + s - 1``. On masked panels (:func:`~longbet.fit_rank_k_longbet`
with ``NaN`` in ``y``) the cells beyond a unit's observed window are model
predictions, so the cumulative effects below are *population* quantities --
every treated unit counts (``weighting="unit"``) whether or not it was observed
through the horizon -- and :func:`rank_k_transport_audit` checks the assumption
that carries the prediction: no cohort-by-event-time effect modification given
the covariates, tested on the event times that both strata observe.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from longbet._rank_k import RankKLongBetResult, _infer_panel_masks_from_result

_WEIGHTINGS = ("unit", "cell")


def _event_time_cells(res: RankKLongBetResult, horizon: int | None) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """``(tau (S, N_trt, H), observed (N_trt, H), treated index, H)`` of the treated units in event time."""
    catt = np.asarray(res.catt_draws, dtype=np.float64)
    S, N, T = catt.shape
    cohorts, ever, _ = _infer_panel_masks_from_result(res, None, None, None)
    t0 = int(min(res.bases)) if res.bases else 1
    h_fit = int(res.horizon) if getattr(res, "horizon", None) is not None else T - t0
    H = h_fit if horizon is None else int(horizon)
    if not 1 <= H <= h_fit:
        raise ValueError(f"horizon must be an integer in [1, {h_fit}] (the fit's event-time horizon); got {horizon!r}")
    trt = np.nonzero(ever)[0]
    tau = np.zeros((S, trt.size, H), dtype=np.float64)
    observed = np.zeros((trt.size, H), dtype=bool)
    catt_obs = getattr(res, "catt_observed", None)
    for j, i in enumerate(trt):
        g = int(cohorts[i])
        width = min(H, T - g)
        tau[:, j, :width] = catt[:, i, g : g + width]
        if catt_obs is not None:
            observed[j, :width] = np.asarray(catt_obs, dtype=bool)[i, g : g + width]
        else:
            observed[j, :width] = True
    return tau, observed, trt, H


def _unit_weights(res: RankKLongBetResult, trt: np.ndarray, weighting: str, units: np.ndarray | None) -> np.ndarray:
    """Normalized weights over the treated units (``(N_trt,)``), restricted to ``units`` when given."""
    if weighting not in _WEIGHTINGS:
        raise ValueError(f"weighting must be one of {_WEIGHTINGS}; got {weighting!r}")
    if weighting == "unit":
        w = np.ones(trt.size, dtype=np.float64)
    else:
        stored = getattr(res, "post_window_length", None)
        if stored is not None:
            w = np.asarray(stored, dtype=np.float64)[trt]
        else:
            uw = getattr(res, "unit_weights", None)
            w = np.asarray(uw, dtype=np.float64)[trt] if uw is not None else np.ones(trt.size)
    if units is not None:
        sel = np.asarray(units)
        if sel.dtype == bool:
            if sel.shape[0] != int(res.catt_draws.shape[1]):
                raise ValueError("a boolean `units` mask must have one entry per unit of the panel")
            keep = sel[trt]
        else:
            keep = np.isin(trt, sel.astype(np.int64))
        w = np.where(keep, w, 0.0)
    total = float(np.sum(w))
    if total <= 0.0:
        raise ValueError("no treated unit selected (or all selected units have zero weight)")
    return w / total


def event_time_trajectory(
    res: RankKLongBetResult,
    *,
    weighting: str = "unit",
    units: np.ndarray | None = None,
    horizon: int | None = None,
) -> np.ndarray:
    """Posterior draws of the weighted mean cell effect ``tau-bar(s)`` by event time.

    Args:
      res: A :class:`~longbet.RankKLongBetResult`.
      weighting: ``"unit"`` (every treated unit counts once; the population
        trajectory) or ``"cell"`` (units weighted by their observed window length
        ``L_i``, the paper's ATT convention).
      units: Optional subset of units -- a boolean ``(N,)`` mask or integer indices;
        never-treated units are ignored.
      horizon: Number of event times to return (default: the fit's horizon ``H``).

    Returns:
      ``(S, H)`` draws of ``sum_i w_i tau_i(s)`` for ``s = 1..H``.
    """
    tau, _, trt, H = _event_time_cells(res, horizon)
    w = _unit_weights(res, trt, weighting, units)
    return np.einsum("snh,n->sh", tau, w)


def cumulative_effect(
    res: RankKLongBetResult,
    horizon: int,
    *,
    weighting: str = "unit",
    units: np.ndarray | None = None,
) -> np.ndarray:
    """Posterior draws of the ``horizon``-period cumulative effect per treated unit.

    ``sum_i w_i sum_{s <= horizon} tau_i(s)``: with ``weighting="unit"`` the
    population estimand ``ATE_h`` of the methodology plan (every treated unit counts
    once, observed or not through ``h``). On a masked fit the cells beyond a unit's
    window are posterior-predictive; see :func:`rank_k_transport_audit` for the
    share they carry.

    Args:
      res: A :class:`~longbet.RankKLongBetResult`.
      horizon: Number of post-adoption periods ``h`` to accumulate (``1 <= h <= H``).
      weighting: ``"unit"`` or ``"cell"`` (see :func:`event_time_trajectory`).
      units: Optional subset of units (boolean mask or indices).

    Returns:
      ``(S,)`` posterior draws of the cumulative effect.
    """
    traj = event_time_trajectory(res, weighting=weighting, units=units, horizon=int(horizon))
    return traj.sum(axis=1)


@dataclass
class TransportAuditResult:
    """Output of :func:`rank_k_transport_audit`.

    Attributes:
      labels: The two stratum labels compared, ``(a, b)``; differences are ``b - a``.
      overlap_horizon: Number of event times compared (both strata observe them).
      n_units: Treated units per stratum, ``{label: count}``.
      profile_draws: ``{label: (S, overlap_horizon)}`` unit-weighted mean effect draws.
      diff_draws: ``(S, overlap_horizon)`` draws of ``profile_b - profile_a``.
      mean_diff_draws: ``(S,)`` draws of the difference averaged over the overlap window.
      quantiles: ``{"p5", "p30", "p50", "p70", "p95"}`` of ``mean_diff_draws``.
      prob_positive: Posterior probability that the overlap-window difference is > 0.
      flag: Whether the 90% interval (P5, P95) of the overlap-window difference
        excludes 0 (or excludes ``[-threshold, threshold]`` when a threshold is set).
      threshold: The absolute tolerance used by ``flag`` (``None``: zero).
      unobserved_cell_share: Per stratum (and ``"all"``), share of the
        ``cumulative_horizon`` treated unit-cells that were not observed.
      unobserved_effect_share: Per stratum (and ``"all"``), share of the posterior
        mean ``cumulative_horizon``-period cumulative effect carried by unobserved
        cells (NaN when the cumulative effect is ~0).
      cumulative_horizon: Horizon of the unobserved-share diagnostics.
    """

    labels: tuple[Any, Any]
    overlap_horizon: int
    n_units: dict[Any, int]
    profile_draws: dict[Any, np.ndarray]
    diff_draws: np.ndarray
    mean_diff_draws: np.ndarray
    quantiles: dict[str, float]
    prob_positive: float
    flag: bool
    threshold: float | None
    unobserved_cell_share: dict[Any, float]
    unobserved_effect_share: dict[Any, float]
    cumulative_horizon: int

    def summary(self) -> str:
        """One-paragraph human-readable summary."""
        q = self.quantiles
        a, b = self.labels
        return (
            f"Transport audit on event times 1..{self.overlap_horizon}: mean effect difference "
            f"({b!r} minus {a!r}) P50 = {q['p50']:.4g} [P5 {q['p5']:.4g}, P95 {q['p95']:.4g}; "
            f"P30 {q['p30']:.4g}, P70 {q['p70']:.4g}], Pr(> 0) = {self.prob_positive:.2f}; "
            f"flag = {self.flag}. Unobserved share of the {self.cumulative_horizon}-period cumulative: "
            + ", ".join(
                f"{k!r}: {self.unobserved_cell_share[k]:.1%} of cells / "
                f"{self.unobserved_effect_share[k]:.1%} of effect"
                for k in self.unobserved_cell_share
            )
        )


def rank_k_transport_audit(
    res: RankKLongBetResult,
    strata: np.ndarray,
    *,
    overlap_horizon: int = 5,
    pair: tuple[Any, Any] | None = None,
    threshold: float | None = None,
    cumulative_horizon: int | None = None,
) -> TransportAuditResult:
    """Compare the early event-time effect profiles of two strata of treated units.

    Masked event-time fits predict the late cells of short-window units from the
    long-window units through the covariates. This audit compares the
    unit-weighted mean effect ``tau-bar(s)`` of two strata (for instance
    ``L_i >= 12`` vs ``L_i <= 7``) on the first ``overlap_horizon`` event times,
    which both strata observe: a difference there is direct evidence of a
    cohort-by-event-time effect modification that the covariates do not capture,
    and hence of a transport risk for the predicted cells. The audit also reports
    how much of the ``cumulative_horizon``-period cumulative effect rests on
    unobserved cells, per stratum.

    Args:
      res: A :class:`~longbet.RankKLongBetResult` (typically a masked fit).
      strata: ``(N,)`` stratum labels (any hashable dtype); ``None`` / ``NaN`` and
        never-treated units are ignored.
      overlap_horizon: Event times ``1..overlap_horizon`` compared (both strata
        should observe them; the audit does not check this).
      pair: The two labels ``(a, b)`` to compare; default: the two distinct labels
        present (sorted), error if there are not exactly two.
      threshold: Optional absolute tolerance: the flag fires when the 90% interval
        of the overlap-window difference lies entirely outside
        ``[-threshold, threshold]``; ``None`` tests against 0.
      cumulative_horizon: Horizon of the unobserved-share diagnostics (default: the
        fit's horizon).

    Returns:
      A :class:`TransportAuditResult`.
    """
    tau, observed, trt, H = _event_time_cells(res, None)
    h_ov = int(overlap_horizon)
    if not 1 <= h_ov <= H:
        raise ValueError(f"overlap_horizon must be in [1, {H}]; got {overlap_horizon!r}")
    h_cum = H if cumulative_horizon is None else int(cumulative_horizon)
    if not 1 <= h_cum <= H:
        raise ValueError(f"cumulative_horizon must be in [1, {H}]; got {cumulative_horizon!r}")
    labels_all = np.asarray(strata, dtype=object)
    if labels_all.shape[0] != int(res.catt_draws.shape[1]):
        raise ValueError("strata must hold one label per unit of the panel")
    lab = labels_all[trt]
    valid = np.array([not (v is None or (isinstance(v, float) and np.isnan(v))) for v in lab], dtype=bool)
    present = sorted({v for v in lab[valid]}, key=lambda v: (str(type(v)), v))
    if pair is None:
        if len(present) != 2:
            raise ValueError(f"strata must hold exactly two labels among the treated units (got {present}); pass pair=(a, b)")
        a, b = present
    else:
        a, b = pair
        if a not in present or b not in present:
            raise ValueError(f"pair {pair!r} not found among the treated units' labels {present}")
    sel = {a: valid & (lab == a), b: valid & (lab == b)}
    for key, mask in sel.items():
        if int(mask.sum()) < 2:
            raise ValueError(f"stratum {key!r} has fewer than 2 treated units")
    profiles = {}
    for key, mask in sel.items():
        w = mask.astype(np.float64) / float(mask.sum())
        profiles[key] = np.einsum("snh,n->sh", tau[:, :, :h_ov], w)
    diff = profiles[b] - profiles[a]
    mean_diff = diff.mean(axis=1)
    p5, p30, p50, p70, p95 = np.percentile(mean_diff, [5, 30, 50, 70, 95])
    tol = 0.0 if threshold is None else abs(float(threshold))
    flag = bool(p5 > tol or p95 < -tol)
    cum_tau_hat = tau[:, :, :h_cum].mean(axis=0)  # (N_trt, h_cum) posterior means
    unobs = ~observed[:, :h_cum]
    cell_share: dict[Any, float] = {}
    effect_share: dict[Any, float] = {}
    for key, mask in list(sel.items()) + [("all", np.ones(trt.size, dtype=bool))]:
        cells = unobs[mask]
        cell_share[key] = float(cells.mean()) if cells.size else float("nan")
        total = float(np.sum(cum_tau_hat[mask]))
        effect_share[key] = float(np.sum(cum_tau_hat[mask][cells]) / total) if abs(total) > 1e-12 else float("nan")
    return TransportAuditResult(
        labels=(a, b),
        overlap_horizon=h_ov,
        n_units={a: int(sel[a].sum()), b: int(sel[b].sum())},
        profile_draws=profiles,
        diff_draws=diff,
        mean_diff_draws=mean_diff,
        quantiles={"p5": float(p5), "p30": float(p30), "p50": float(p50), "p70": float(p70), "p95": float(p95)},
        prob_positive=float(np.mean(mean_diff > 0.0)),
        flag=flag,
        threshold=threshold,
        unobserved_cell_share=cell_share,
        unobserved_effect_share=effect_share,
        cumulative_horizon=h_cum,
    )
