"""Single-wave randomized encouragement: descriptions and reference inference.

Inference assumes complete randomization across individual units and a permanent
holdout. Array validation cannot establish randomization, exclusion, monotonicity,
or assumptions about potential adoption histories. The ratio is called Wald,
not automatically CACE. No sampler or latent outcome scale is involved here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy.special import ndtri_exp

from longbet._model import _check_absorbing, _check_time_vector, derive_exposure


@dataclass(frozen=True)
class _Panel:
    z: np.ndarray
    d: np.ndarray
    t: np.ndarray
    assigned: np.ndarray
    start: int
    exposure: np.ndarray

    def metadata(self) -> dict[str, Any]:
        return {
            "design": "single_wave_complete",
            "n_units": self.z.shape[0],
            "n_periods": self.z.shape[1],
            "n_encouraged": int(self.assigned.sum()),
            "n_control": int((~self.assigned).sum()),
            "encouragement_period": float(self.t[self.start]),
            "encouragement_index": self.start,
            "has_pre_periods": self.start > 0,
            "perfect_compliance": bool(np.array_equal(self.z, self.d)),
        }


def _binary_panel(value: Any, name: str) -> np.ndarray:
    a = np.asarray(value)
    if a.dtype.kind not in "biuf":
        raise ValueError(f"{name} must contain finite binary numeric indicators (0 or 1).")
    # Convert BEFORE differencing: boolean differences and unsigned underflow
    # otherwise hide a 1 -> 0 transition from the shared absorbing check.
    a = a.astype(np.float64)
    try:
        _check_absorbing(a)
    except ValueError as exc:
        detail = str(exc).replace("z must", name + " must")
        if "switches back" in detail:
            # Dropping units because of post-assignment uptake would change the
            # randomized comparison. Do not inherit the scalar fitting remedy.
            detail = detail.partition(" Drop or reshape")[0]
            detail += " Use an analysis that supports nonabsorbing panels."
        raise ValueError(f"{name}: {detail}") from exc
    return a


def _validate(z: Any, d: Any, t: Any) -> _Panel:
    z = _binary_panel(z, "z")
    d = _binary_panel(d, "d")
    if d.shape != z.shape:
        raise ValueError(f"d has shape {d.shape}; expected z's shape {z.shape}.")
    t = (np.arange(1, z.shape[1] + 1, dtype=np.float64) if t is None
         else np.atleast_1d(np.asarray(t, dtype=np.float64)))
    if t.ndim != 1 or t.size != z.shape[1]:
        raise ValueError(f"t must be a 1-D vector of length {z.shape[1]}.")
    if np.all(np.isfinite(t)) and not np.allclose(
        np.diff(t), np.rint(np.diff(t)), rtol=0, atol=1e-8,
    ):
        raise ValueError("t must have whole-unit gaps for encouragement exposure indexing.")
    _check_time_vector(t)
    assigned = z[:, -1] == 1
    if not assigned.any() or assigned.all():
        raise ValueError("z must include both an encouraged arm and a permanent control arm.")
    starts = np.argmax(z[assigned] == 1, axis=1)
    if np.unique(starts).size != 1:
        raise ValueError(
            "Only a single encouragement wave is supported; z has multiple start periods. "
            "Staggered encouragement needs cohort-specific contrasts and weights."
        )
    start = int(starts[0])
    exposure = derive_exposure(z, t)[np.flatnonzero(assigned)[0]]
    return _Panel(z, d, t, assigned, start, exposure)


def validate_encouragement(z: Any, d: Any, t: Any = None) -> dict[str, Any]:
    """Validate a single-wave encouragement panel and return design metadata.

    ``z`` (assignment from the encouragement period onward) and ``d`` (actual
    adoption) must be complete, binary, absorbing ``(N, T)`` arrays. All assigned
    units must start encouragement together, with a permanent control arm.
    ``t`` defaults to ``1..T``; strictly increasing whole-unit gaps are supported.

    Perfect compliance, adoption before encouragement, negative observed first
    stages, and experiments without pre-periods are accepted. Adoption need not
    be immediate. Restrictions on absorption are this implementation's scope,
    not general requirements of instrumental variables.

    Returns a dictionary of dimensions, arm counts, encouragement period and
    zero-based column index, ``has_pre_periods``, and ``perfect_compliance``.
    ``design='single_wave_complete'`` names the assumed randomization scheme;
    it does not certify that the experiment was randomized. Individual-unit
    complete randomization is required for the inference functions. Blocked,
    unequal-probability, and cluster assignment need different estimators.
    """
    return _validate(z, d, t).metadata()


def _critical(alpha: float) -> float:
    if isinstance(alpha, (bool, np.bool_)) or not np.isscalar(alpha):
        raise ValueError("alpha must be a finite number strictly between 0 and 1.")
    try:
        alpha = float(alpha)
    except (TypeError, ValueError) as exc:
        raise ValueError("alpha must be a finite number strictly between 0 and 1.") from exc
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be a finite number strictly between 0 and 1.")
    # Evaluate the tail in log space so alpha/2 cannot underflow to zero.
    return float(-ndtri_exp(np.log(alpha) - np.log(2)))


def _difference(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    delta = float(a.mean() - b.mean())
    variance = (float(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size)
                if min(a.size, b.size) >= 2 else np.nan)
    return delta, variance


def _period(panel: _Panel, j: int) -> dict[str, Any]:
    return {
        "period": float(panel.t[j]),
        "period_index": j,
        "horizon": int(panel.exposure[j]),
        "post_encouragement": j >= panel.start,
        "n_encouraged": int(panel.assigned.sum()),
        "n_control": int((~panel.assigned).sum()),
    }


def encouragement_summary(
    z: Any, d: Any, t: Any = None, *, alpha: float = 0.05,
) -> pd.DataFrame:
    """Describe adoption in the assigned arms, one row per calendar period.

    Inputs follow :func:`validate_encouragement`. Baseline periods are included;
    arm membership always refers to eventual randomized assignment. ``horizon``
    uses LongBet's elapsed-time exposure index (zero before encouragement), and
    ``period_index`` is a zero-based observed column, so uneven gaps are explicit.

    ``first_stage`` is the take-up-rate difference. Its SE uses separate arm
    sample variances divided by unit counts. Bounds are pointwise normal
    confidence intervals under individual-unit complete randomization, not
    Bayesian or finite-sample exact intervals. Inference is missing if either
    arm has fewer than two units. Intervals are not clipped to [-1, 1].

    Lag columns describe encouraged units *observed adopted by that period*:
    ``share_immediate`` is the fraction first adopting at encouragement,
    ``median_observed_lag`` is calendar adoption time minus encouragement time,
    and ``share_pre_encouragement_adoption`` is the fraction adopting beforehand.
    These fractions are missing when no such adopters exist. They do not identify
    compliers or test counterfactual timing assumptions. Negative first stages
    and intervals including zero are retained, not treated as invalid inputs.
    """
    panel = _validate(z, d, t)
    critical = _critical(alpha)
    takeup = panel.d[panel.assigned]
    first = np.argmax(takeup == 1, axis=1)
    rows = []
    for j in range(panel.t.size):
        a, b = takeup[:, j], panel.d[~panel.assigned, j]
        delta, variance = _difference(a, b)
        se = float(np.sqrt(variance))
        lower, upper = delta - critical * se, delta + critical * se
        adopted = a == 1
        lags = panel.t[first[adopted]] - panel.t[panel.start]
        rows.append({
            **_period(panel, j),
            "takeup_encouraged": float(a.mean()),
            "takeup_control": float(b.mean()),
            "first_stage": delta,
            "first_stage_se": se,
            "first_stage_lower": lower,
            "first_stage_upper": upper,
            "first_stage_includes_zero": (bool(lower <= 0 <= upper)
                                           if np.isfinite(se) else np.nan),
            "negative_first_stage": delta < 0,
            "n_encouraged_adopters": int(adopted.sum()),
            "share_immediate": float(np.mean(lags == 0)) if lags.size else np.nan,
            "median_observed_lag": float(np.median(lags)) if lags.size else np.nan,
            "share_pre_encouragement_adoption": (float(np.mean(lags < 0))
                                                 if lags.size else np.nan),
        })
    return pd.DataFrame(rows)


def _quadratic_set(a: float, b: float, c: float) -> tuple[str, list[tuple[float, float]]]:
    """Solve a*x**2 + b*x + c <= 0, preserving unbounded components."""
    if a == 0:
        if b == 0:
            return ("all_real", [(-np.inf, np.inf)]) if c <= 0 else ("empty", [])
        root = -c / b
        return "half_line", [(-np.inf, root) if b > 0 else (root, np.inf)]
    discriminant = b * b - 4 * a * c
    # An exactly tangent quadratic can acquire a tiny negative discriminant
    # through cancellation. Scale tolerance to the arithmetic, not effect units.
    tolerance = 16 * np.finfo(float).eps * (b * b + abs(4 * a * c))
    if -tolerance <= discriminant < 0:
        discriminant = 0.0
    if discriminant < 0:
        return ("all_real", [(-np.inf, np.inf)]) if a < 0 else ("empty", [])
    if discriminant == 0:
        root = -b / (2 * a)
        return (("all_real", [(-np.inf, np.inf)]) if a < 0
                else ("singleton", [(root, root)]))
    # Stable quadratic formula: avoid cancellation at the smaller root.
    q = -0.5 * (b + np.copysign(np.sqrt(discriminant), b))
    lo, hi = sorted((q / a, c / q))
    if a > 0:
        return "bounded", [(lo, hi)]
    return "disjoint", [(-np.inf, lo), (hi, np.inf)]


def _wald_set(dy: float, dd: float, vy: float, vd: float, cyd: float,
              critical: float) -> dict[str, Any]:
    result = {
        "wald_set_type": "unavailable", "wald_reason": "insufficient_arm_size",
        "wald_lower_1": np.nan, "wald_upper_1": np.nan,
        "wald_lower_2": np.nan, "wald_upper_2": np.nan,
    }
    if not np.all(np.isfinite([vy, vd, cyd])):
        return result
    # Normalize the two outcome scales independently before forming coefficients.
    # A change from dollars to cents must not alter the set's topology.
    sy = max(abs(dy), np.sqrt(vy)) or 1.0
    sd = max(abs(dd), np.sqrt(vd)) or 1.0
    y, d = dy / sy, dd / sd
    q = critical**2
    a = d**2 - q * (vd / sd / sd)
    b = 2 * (q * (cyd / sy / sd) - y * d)
    c = y**2 - q * (vy / sy / sy)
    kind, intervals = _quadratic_set(a, b, c)
    result.update(wald_set_type=kind, wald_reason="")
    for k, (lo, hi) in enumerate(intervals, 1):
        result[f"wald_lower_{k}"] = float(lo * sy / sd)
        result[f"wald_upper_{k}"] = float(hi * sy / sd)
    return result


def encouragement_effects(
    y: Any, d: Any, z: Any, t: Any = None, *, alpha: float = 0.05,
) -> pd.DataFrame:
    """Estimate encouragement ITTs and Wald confidence sets without fitting a model.

    Requires individual-unit complete randomization, one encouragement wave, a
    permanent control arm, and complete panels. See :func:`validate_encouragement`
    for ``z``, ``d``, and ``t``. ``y`` must be finite numeric ``(N, T)``; binary
    outcomes give risk differences directly. No missing cells are dropped.

    Returns one row per observed post-encouragement period. ``itt_y`` and
    ``itt_d`` are assigned-arm differences in means, targeting the same study
    population. Their SEs and ``itt_y_d_cov`` use the joint Neyman covariance
    estimate, the sum of within-arm sample covariance matrices divided by arm
    sizes. Each row has one observation per independent unit; intervals are
    pointwise, not simultaneous across horizons.

    ``wald`` is the signed ratio, missing only for an exactly zero denominator.
    ``wald_set_type`` and up to two pairs ``wald_lower_1``/``wald_upper_1`` and
    ``wald_lower_2``/``wald_upper_2`` describe the closed confidence set. Infinite
    endpoints are preserved; unused components are NaN. Sets can be bounded,
    disjoint, half_line, all_real, singleton, empty, or unavailable. The last
    case retains point estimates but sets ``wald_reason='insufficient_arm_size'``.

    Inference analytically inverts the studentized contrast of ``Y - w*D`` using
    a normal critical value (Fieller/Anderson–Rubin-style inference). It does not
    divide by the first stage to form a confidence set, but remains an asymptotic
    approximation requiring adequate arm sizes and moment conditions. It is not
    an exact permutation procedure or a Bayesian credible interval. No 5% first
    stage threshold or truncation of unbounded sets is used.

    Randomization identifies ITTs, not automatically a treatment effect from
    their ratio. A CACE interpretation additionally needs exclusion, monotonicity,
    relevance, and appropriate adoption-history assumptions. The function neither
    verifies these assumptions nor assigns latent compliance types.
    """
    panel = _validate(z, d, t)
    critical = _critical(alpha)
    y = np.asarray(y)
    if y.shape != panel.z.shape or y.dtype.kind not in "biuf":
        raise ValueError(f"y must be a finite numeric array with shape {panel.z.shape}.")
    y = y.astype(np.float64)
    if not np.all(np.isfinite(y)):
        raise ValueError("y must be complete and finite; missing cells are not silently dropped.")
    rows = []
    for j in range(panel.start, panel.t.size):
        ya, yb = y[panel.assigned, j], y[~panel.assigned, j]
        da, db = panel.d[panel.assigned, j], panel.d[~panel.assigned, j]
        dy, vy = _difference(ya, yb)
        dd, vd = _difference(da, db)
        cyd = (float(np.cov(ya, da, ddof=1)[0, 1] / ya.size
                     + np.cov(yb, db, ddof=1)[0, 1] / yb.size)
               if min(ya.size, yb.size) >= 2 else np.nan)
        if min(ya.size, yb.size) >= 2 and not np.all(np.isfinite([dy, vy, dd, vd, cyd])):
            raise ValueError("y produces nonfinite moments; rescale outcomes before inference.")
        row = {**_period(panel, j), "itt_y_d_cov": cyd,
               "wald": dy / dd if dd != 0 else np.nan,
               "inference": "normal_ar"}
        for name, delta, variance in (("itt_y", dy, vy), ("itt_d", dd, vd)):
            se = float(np.sqrt(variance))
            row.update({name: delta, f"{name}_se": se,
                        f"{name}_lower": delta - critical * se,
                        f"{name}_upper": delta + critical * se})
        row.update(_wald_set(dy, dd, vy, vd, cyd, critical))
        rows.append(row)
    return pd.DataFrame(rows)


def plot_encouragement(z: Any, d: Any, t: Any = None, *, ax: Any = None) -> Any:
    """Plot adoption rates in assigned arms; matplotlib is an optional dependency.

    Curves include baseline periods, the signed gap is shaded, and the common
    encouragement time is marked. No compliance classification is inferred.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError("plot_encouragement() needs matplotlib; install longbet-jax[plot].") from exc
    frame = encouragement_summary(z, d, t)
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(frame.period, frame.takeup_encouraged, label="Assigned encouragement")
    ax.plot(frame.period, frame.takeup_control, label="Assigned control")
    ax.fill_between(frame.period, frame.takeup_encouraged, frame.takeup_control, alpha=0.15)
    start = frame.loc[frame.post_encouragement, "period"].iloc[0]
    ax.axvline(start, color="grey", linestyle="--", label="Encouragement begins")
    ax.set(xlabel="Period", ylabel="Adoption rate", ylim=(-0.02, 1.02))
    ax.legend()
    return ax
