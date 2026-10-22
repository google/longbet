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

"""Prior diagnostics, prior-sensitivity tiers, decay extrapolation and figures for the Rank-K estimator.

The decision layer of a Rank-K BCF analysis asks the same questions every time:
how much did the data move the prior, how sensitive is the horizon effect to the
prior family and scale, and what does the observed post-adoption window imply
once it is extrapolated with a geometric decay. This module turns the hand-built
notebook versions of those answers into library calls, split -- as in
:mod:`longbet._rollout` -- into a pure numpy / pandas data layer and thin
``plot_*`` functions that import matplotlib only when called:

* :func:`sample_rank_k_prior_draws` draws the mode-1 ATT prior of a two-arm fit
  and projects it onto post-adoption trajectories with the same zero-sum
  higher-mode damping as the posterior (:func:`_project_cohort_damped_trajectory`).
* :func:`prior_vs_posterior_summary` tabulates prior and posterior quantiles of
  the horizon statistic (the sum of the trajectory over its periods), the lift
  against a control baseline total, the interval-width change and an overlap flag;
  :func:`assert_nested_additivity` checks ``AC = BC + AB``.
* :func:`evaluate_prior_sensitivity` reruns :func:`reconcile_nested_3arm_draws`
  across :class:`PriorTier` tiers (:func:`ytlr_prior_tiers`), tabulating the
  standardized prior families analytically, the observed- and extrapolated-horizon
  posteriors and the commercial break-even ``Pr(ROI > 0 | alpha, c)``.
* :func:`extrapolate_rank_k_decay` extends observed trajectory draws with
  ``tau(L + h) = rho^h tau(L)`` ("LongBet Regime 4"), with an optional hierarchical
  ``rho ~ Beta`` band, a ``rho`` grid, milestone tables and a data-health certificate.

Every data-layer function raises ``ValueError`` on malformed or non-finite inputs
instead of propagating NaNs, and every random step takes an explicit seed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Any
import warnings

import numpy as np
import pandas as pd
from scipy import special as sp_special
from scipy import stats as sp_stats

from longbet._rank_k import (
    Nested3ArmReconciliationResult,
    RankKConfigWarning,
    RankKLongBetResult,
    _infer_panel_masks_from_result,
    _project_cohort_damped_trajectory,
    reconcile_nested_3arm_draws,
)

#: Mode-1 ATT prior families accepted by :func:`sample_rank_k_prior_draws`.
PRIOR_FAMILIES = ("normal", "halfnormal", "halfstudentt", "gamma_spike")
#: Bundle prior modes of :func:`reconcile_nested_3arm_draws` accepted in a :class:`PriorTier`.
BUNDLE_PRIOR_MODES = ("bundle_halfnormal", "bundle_halfstudentt", "bundle_gamma_spike")
#: Contrast labels of a nested three-arm design, in the order the reconciliation returns them.
NESTED_CONTRASTS = ("BC", "AB", "AC")
#: Default central quantiles of every summary table (P30 / P50 / P70).
DEFAULT_QUANTILES = (0.30, 0.50, 0.70)
#: Default months at which the decay extrapolation is tabulated.
DEFAULT_MILESTONES = (1, 3, 5, 12, 24, 36, 48, 60)

_FAMILY_OF_MODE = {
    "bundle_halfnormal": "halfnormal",
    "bundle_halfstudentt": "halfstudentt",
    "bundle_gamma_spike": "gamma_spike",
}
_FAMILY_DISPLAY = {
    "normal": "Normal",
    "halfnormal": "Half-Normal",
    "halfstudentt": "Half-Student-t (df 2.5, scale 0.5 s)",
    "gamma_spike": "Gamma spike (shape 0.40, scale 1.60 s)",
}
# Parameterization shared with ``_draw_bundle_prior`` / ``tau_prior_family`` of ``_rank_k``.
_HALF_T_DF = 2.5
_HALF_T_SCALE_FACTOR = 0.50
_GAMMA_SHAPE = 0.40
_GAMMA_SCALE_FACTOR = 1.60
_ARM_COLORS = ("#3b6ea5", "#d1495b", "#2e8b57", "#8e6c8a", "#c98b2a")


# ---- Shared validation and summary helpers -------------------------------------------------------


def _as_draw_matrix(value: Any, name: str) -> np.ndarray:
    """Coerce ``value`` to a finite float64 ``(S, T)`` matrix of draws, raising ``ValueError`` otherwise."""
    try:
        arr = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric array of draws; got {type(value).__name__}") from exc
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2-D (draws, periods) array; got shape {arr.shape}")
    if arr.shape[0] < 2 or arr.shape[1] < 1:
        raise ValueError(f"{name} needs at least 2 draws and 1 period; got shape {arr.shape}")
    if not np.all(np.isfinite(arr)):
        n_bad = int(np.sum(~np.isfinite(arr)))
        raise ValueError(f"{name} contains {n_bad} NaN / Inf entries; the draws must be finite")
    return arr


def _as_traj_dict(value: Any, name: str, default_label: str = "effect") -> dict[str, np.ndarray]:
    """Normalize a single ``(S, T)`` array or a ``{label: (S, T)}`` mapping to an ordered dict of matrices."""
    if isinstance(value, Mapping):
        if not value:
            raise ValueError(f"{name} is an empty mapping; pass at least one (S, T) array")
        return {str(key): _as_draw_matrix(arr, f"{name}[{key!r}]") for key, arr in value.items()}
    return {default_label: _as_draw_matrix(value, name)}


def _check_quantiles(quantiles: Sequence[float]) -> tuple[float, ...]:
    """Validate a strictly increasing tuple of probabilities in ``(0, 1)`` with an odd length."""
    qs = tuple(float(q) for q in quantiles)
    if len(qs) < 1 or len(qs) % 2 == 0:
        raise ValueError(f"quantiles must have an odd length (low, ..., central, ..., high); got {qs}")
    if any(not 0.0 < q < 1.0 for q in qs) or any(b <= a for a, b in zip(qs, qs[1:])):
        raise ValueError(f"quantiles must be strictly increasing probabilities in (0, 1); got {qs}")
    return qs


def _quantile_label(q: float) -> str:
    """``0.30 -> 'P30'``; probabilities that are not whole percents keep one decimal (``'P2.5'``)."""
    pct = 100.0 * q
    return f"P{int(round(pct))}" if abs(pct - round(pct)) < 1e-9 else f"P{pct:g}"


def _summarize(values: np.ndarray, quantiles: Sequence[float], prefix: str = "") -> dict[str, float]:
    """Quantiles, mean and SD (ddof 1) of a 1-D sample as ``{prefix + 'P30': ..., prefix + 'Mean': ...}``."""
    values = np.asarray(values, dtype=np.float64)
    out = {f"{prefix}{_quantile_label(q)}": float(np.quantile(values, q)) for q in quantiles}
    out[f"{prefix}Mean"] = float(np.mean(values))
    out[f"{prefix}SD"] = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
    return out


def _positive_scalar(value: Any, name: str) -> float:
    """Return ``value`` as a finite positive float, raising ``ValueError`` otherwise."""
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive number; got {value!r}") from exc
    if not (math.isfinite(out) and out > 0.0):
        raise ValueError(f"{name} must be a finite positive number; got {value!r}")
    return out


def _check_tau_clip(tau_clip: Sequence[float]) -> tuple[float, float]:
    """Validate ``(low, high)`` with ``low < high`` (endpoints may be infinite)."""
    if len(tau_clip) != 2:
        raise ValueError(f"tau_clip must be a (low, high) pair; got {tau_clip!r}")
    low, high = float(tau_clip[0]), float(tau_clip[1])
    if math.isnan(low) or math.isnan(high) or not low < high:
        raise ValueError(f"tau_clip must satisfy low < high; got {tau_clip!r}")
    return low, high


# ---- Prior families ------------------------------------------------------------------------------


class _HalfStudentT:
    """Half-Student-t ``|scale * T_df|`` with the frozen-scipy method subset this module needs.

    scipy has no half-t family; the quantile, distribution and density functions
    follow from the symmetric ``t`` (``F(x) = 2 F_t(x / s) - 1`` on ``x >= 0``), and
    the first two moments are the standard closed forms (finite for ``df > 2``).
    """

    def __init__(self, df: float, scale: float) -> None:
        self.df = float(df)
        self.scale = float(scale)
        self._t = sp_stats.t(df=self.df)

    def pdf(self, x: Any) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        return np.where(x >= 0.0, 2.0 * self._t.pdf(x / self.scale) / self.scale, 0.0)

    def cdf(self, x: Any) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        return np.where(x >= 0.0, 2.0 * self._t.cdf(x / self.scale) - 1.0, 0.0)

    def ppf(self, p: Any) -> np.ndarray:
        p = np.asarray(p, dtype=np.float64)
        return self.scale * self._t.ppf(0.5 * (1.0 + p))

    def mean(self) -> float:
        if self.df <= 1.0:
            return float("inf")
        log_ratio = sp_special.gammaln(0.5 * (self.df + 1.0)) - sp_special.gammaln(0.5 * self.df)
        numer = 2.0 * math.sqrt(self.df) * math.exp(log_ratio)
        return float(self.scale * numer / (math.sqrt(math.pi) * (self.df - 1.0)))

    def var(self) -> float:
        if self.df <= 2.0:
            return float("inf")
        return float(self.scale**2 * self.df / (self.df - 2.0) - self.mean() ** 2)

    def std(self) -> float:
        return float(math.sqrt(self.var()))


def _prior_family_distribution(family: str, scale: float, mean: float = 0.0) -> Any:
    """Frozen distribution of a mode-1 ATT prior family, parameterized exactly as the sampler does.

    ``"halfnormal"`` is ``N^+(0, scale^2)``, ``"halfstudentt"`` is ``T^+_{2.5}(0, (0.5 scale)^2)``,
    ``"gamma_spike"`` is ``Gamma(0.40, 1.60 scale)`` and ``"normal"`` is ``N(mean, scale^2)``
    (truncation to ``tau_clip`` is applied by the caller). The returned object exposes
    ``pdf``, ``cdf``, ``ppf``, ``mean()`` and ``std()``.
    """
    if family not in PRIOR_FAMILIES:
        raise ValueError(f"prior_family must be one of {PRIOR_FAMILIES}; got {family!r}")
    scale = _positive_scalar(scale, "prior_scale")
    if family == "halfnormal":
        return sp_stats.halfnorm(loc=0.0, scale=scale)
    if family == "halfstudentt":
        return _HalfStudentT(df=_HALF_T_DF, scale=_HALF_T_SCALE_FACTOR * scale)
    if family == "gamma_spike":
        return sp_stats.gamma(a=_GAMMA_SHAPE, scale=_GAMMA_SCALE_FACTOR * scale)
    mean = float(mean)
    if not math.isfinite(mean):
        raise ValueError(f"prior_mean must be finite; got {mean!r}")
    return sp_stats.norm(loc=mean, scale=scale)


def _family_of_prior_mode(prior_mode: str) -> str:
    """Map a bundle ``prior_mode`` of the reconciliation onto its ATT prior family."""
    if prior_mode not in BUNDLE_PRIOR_MODES:
        raise ValueError(f"prior_mode must be one of {BUNDLE_PRIOR_MODES}; got {prior_mode!r}")
    return _FAMILY_OF_MODE[prior_mode]


@dataclass(frozen=True)
class PriorTier:
    """One prior-sensitivity tier of :func:`evaluate_prior_sensitivity`.

    Attributes:
      name: Display name (unique within a tier list).
      prior_mode: Bundle prior mode passed to :func:`reconcile_nested_3arm_draws`; one of
        :data:`BUNDLE_PRIOR_MODES`.
      prior_scale: Positive scale in outcome units, or ``"in_sample"`` to use the fit's own
        posterior SD (the reconciliation warns about that choice; the sensitivity driver
        silences the warning for this tier only).
    """

    name: str
    prior_mode: str
    prior_scale: float | str

    def __post_init__(self) -> None:
        if not str(self.name):
            raise ValueError("PriorTier.name must be a non-empty string")
        _family_of_prior_mode(self.prior_mode)
        if isinstance(self.prior_scale, str):
            if self.prior_scale != "in_sample":
                raise ValueError(
                    f"PriorTier.prior_scale must be a positive number or 'in_sample'; got {self.prior_scale!r}"
                )
        else:
            _positive_scalar(self.prior_scale, "PriorTier.prior_scale")

    @property
    def family(self) -> str:
        """ATT prior family implied by ``prior_mode`` (``"halfnormal"``, ``"halfstudentt"`` or ``"gamma_spike"``)."""
        return _family_of_prior_mode(self.prior_mode)

    @property
    def is_in_sample(self) -> bool:
        """Whether the scale is the data-dependent ``"in_sample"`` posterior SD."""
        return isinstance(self.prior_scale, str)


def ytlr_prior_tiers(sigma_hist: float) -> tuple[PriorTier, ...]:
    """The four YTLR prior-sensitivity tiers for a historical benchmark scale ``sigma_hist``.

    Returns ``("Out-of-sample Half-Normal", "bundle_halfnormal", sigma_hist)``,
    ``("In-sample SE Half-Normal", "bundle_halfnormal", "in_sample")``,
    ``("Half-Student-t", "bundle_halfstudentt", sigma_hist)`` and
    ``("Gamma spike", "bundle_gamma_spike", sigma_hist)``.

    Args:
      sigma_hist: Positive out-of-sample calibrated scale of the total bundle effect, in
        per-period outcome units.
    """
    scale = _positive_scalar(sigma_hist, "sigma_hist")
    return (
        PriorTier("Out-of-sample Half-Normal", "bundle_halfnormal", scale),
        PriorTier("In-sample SE Half-Normal", "bundle_halfnormal", "in_sample"),
        PriorTier("Half-Student-t", "bundle_halfstudentt", scale),
        PriorTier("Gamma spike", "bundle_gamma_spike", scale),
    )


# ---- Prior draws of a two-arm fit ----------------------------------------------------------------


@dataclass(frozen=True)
class RankKPriorDraws:
    """Prior draws of the mode-1 ATT of a two-arm Rank-K fit, projected onto trajectories.

    Attributes:
      att_draws: ``(S,)`` prior draws of the mode-1 ATT, restricted to ``tau_clip``.
      traj_draws: ``(S, T_post)`` prior trajectory draws: ``att_draws`` plus the fit's
        cohort-share-weighted higher-mode dynamics, damped with ``catt_floor`` exactly as
        the posterior is (:func:`_project_cohort_damped_trajectory`).
      posterior_att_draws: ``(S_post,)`` posterior ATT draws of the fit (``res.att_draws``).
      posterior_traj_draws: ``(S_post, T_post)`` posterior trajectory draws under the same
        projection, so prior and posterior are compared on the same footing.
      family: Prior family (one of :data:`PRIOR_FAMILIES`).
      scale: Prior scale in per-period outcome units.
      mean: Location of the ``"normal"`` family (ignored by the zero-mode families).
      tau_clip: Truncation interval applied to the ATT draws.
      catt_floor: Floor used by the zero-sum damping, or ``None``.
      seed: Seed of the prior sampler.
    """

    att_draws: np.ndarray
    traj_draws: np.ndarray
    posterior_att_draws: np.ndarray
    posterior_traj_draws: np.ndarray
    family: str
    scale: float
    mean: float
    tau_clip: tuple[float, float]
    catt_floor: float | None
    seed: int


def _truncated_family_draws(dist: Any, low: float, high: float, n_draws: int, rng: np.random.Generator) -> np.ndarray:
    """Inverse-CDF draws of ``dist`` restricted to ``[low, high]`` (exact, no rejection loop)."""
    cdf_low = float(dist.cdf(low)) if math.isfinite(low) else (0.0 if low < 0 else 1.0)
    cdf_high = float(dist.cdf(high)) if math.isfinite(high) else (1.0 if high > 0 else 0.0)
    if not cdf_high - cdf_low > 1e-12:
        raise ValueError(f"the prior puts no mass inside tau_clip=({low}, {high}); widen the clip or the prior")
    u = rng.uniform(cdf_low, cdf_high, size=n_draws)
    draws = np.asarray(dist.ppf(u), dtype=np.float64)
    return np.clip(draws, low, high)


def sample_rank_k_prior_draws(
    res: RankKLongBetResult,
    *,
    prior_family: str = "halfnormal",
    prior_mean: float = 0.0,
    prior_scale: float,
    tau_clip: tuple[float, float] = (0.0, float("inf")),
    catt_floor: float | None = 0.0,
    n_draws: int | None = None,
    seed: int = 2026,
    cohorts: np.ndarray | None = None,
    ever: np.ndarray | None = None,
) -> RankKPriorDraws:
    """Draw the mode-1 ATT prior of a two-arm fit and project it onto post-adoption trajectories.

    The prior on the time-averaged (mode-1) effect is sampled exactly (inverse CDF,
    truncated to ``tau_clip``) from ``"normal"`` (``N(prior_mean, prior_scale^2)``),
    ``"halfnormal"`` (``N^+(0, prior_scale^2)``), ``"halfstudentt"``
    (``T^+_{2.5}(0, (0.5 prior_scale)^2)``) or ``"gamma_spike"``
    (``Gamma(0.40, 1.60 prior_scale)``) -- the parameterization of
    ``RankKLongBetConfig.tau_prior_family`` -- and projected with the fit's own
    higher-mode dynamics and the zero-sum damping of ``catt_floor``, which is what the
    posterior trajectory goes through, so a prior-vs-posterior comparison isolates what
    the data did to the level.

    Args:
      res: Fitted two-arm result (needs ``att_draws``, ``mode_tau_draws`` and ``bases``).
      prior_family: One of :data:`PRIOR_FAMILIES`.
      prior_mean: Mean of the ``"normal"`` family before truncation (ignored otherwise).
      prior_scale: Positive scale in per-period outcome units.
      tau_clip: ``(low, high)`` truncation of the ATT draws (default non-negative).
      catt_floor: Floor of the zero-sum higher-mode damping, or ``None`` for none.
      n_draws: Number of prior draws (default: the number of posterior draws).
      seed: Seed of the prior sampler.
      cohorts: Optional ``(N,)`` first-treated period per unit (``T`` for never treated);
        inferred from ``res`` when omitted, which assigns every treated unit to the first wave.
      ever: Optional ``(N,)`` ever-treated mask; inferred from ``res`` when omitted.

    Returns:
      :class:`RankKPriorDraws`.
    """
    low, high = _check_tau_clip(tau_clip)
    scale = _positive_scalar(prior_scale, "prior_scale")
    if prior_family != "normal" and low < 0.0:
        raise ValueError(f"prior_family={prior_family!r} lives on [0, inf); tau_clip low must be >= 0, got {low}")
    if catt_floor is not None and not math.isfinite(float(catt_floor)):
        raise ValueError(f"catt_floor must be finite or None; got {catt_floor!r}")
    posterior_att = np.asarray(res.att_draws, dtype=np.float64)
    if posterior_att.ndim != 1 or posterior_att.shape[0] < 2:
        raise ValueError(f"res.att_draws must be a 1-D array of at least 2 draws; got shape {posterior_att.shape}")
    if not res.bases:
        raise ValueError("res.bases is empty; the result holds no treated cohort to project onto")
    n_out = int(posterior_att.shape[0]) if n_draws is None else int(n_draws)
    if n_out < 2:
        raise ValueError(f"n_draws must be >= 2; got {n_draws!r}")

    dist = _prior_family_distribution(prior_family, scale, prior_mean)
    rng = np.random.default_rng(int(seed))
    att = _truncated_family_draws(dist, low, high, n_out, rng)

    cohorts_arr, ever_arr, _ = _infer_panel_masks_from_result(res, cohorts, ever, None)
    floor = None if catt_floor is None else float(catt_floor)
    traj = _project_cohort_damped_trajectory(att, res, ever_arr, cohorts_arr, (low, high), floor)
    post_traj = _project_cohort_damped_trajectory(posterior_att, res, ever_arr, cohorts_arr, (low, high), floor)
    return RankKPriorDraws(
        att_draws=att,
        traj_draws=np.asarray(traj, dtype=np.float64),
        posterior_att_draws=posterior_att,
        posterior_traj_draws=np.asarray(post_traj, dtype=np.float64),
        family=prior_family,
        scale=scale,
        mean=float(prior_mean),
        tau_clip=(low, high),
        catt_floor=floor,
        seed=int(seed),
    )


# ---- Prior-vs-posterior tables -------------------------------------------------------------------


def _resolve_labels(keys: Sequence[str], labels: Any, name: str = "labels") -> dict[str, str]:
    """Display names per key: ``None`` keeps the keys, a ``str`` names a single key, a mapping renames."""
    if labels is None:
        return {key: key for key in keys}
    if isinstance(labels, str):
        if len(keys) != 1:
            raise ValueError(f"{name} is a single string but there are {len(keys)} contrasts; pass a mapping")
        return {keys[0]: labels}
    if isinstance(labels, Mapping):
        return {key: str(labels.get(key, key)) for key in keys}
    raise ValueError(f"{name} must be None, a string or a mapping {{key: display name}}; got {type(labels).__name__}")


def _resolve_per_key(value: Any, keys: Sequence[str], name: str) -> dict[str, float | None]:
    """A scalar (shared by every key), a mapping (per key, missing -> None) or ``None`` -> ``{key: float | None}``."""
    if value is None:
        return {key: None for key in keys}
    if isinstance(value, Mapping):
        out: dict[str, float | None] = {}
        for key in keys:
            raw = value.get(key)
            out[key] = None if raw is None else _positive_scalar(raw, f"{name}[{key!r}]")
        return out
    shared = _positive_scalar(value, name)
    return {key: shared for key in keys}


def prior_vs_posterior_summary(
    prior_traj: Any,
    post_traj: Any,
    *,
    labels: Any = None,
    quantiles: Sequence[float] = DEFAULT_QUANTILES,
    ref_total: Any = None,
    unit: str = "",
    assert_overlap: bool = False,
) -> pd.DataFrame:
    """Tabulate prior and posterior quantiles of the horizon statistic per contrast.

    The horizon statistic of a trajectory is the sum of its periods (the cumulative
    effect over the observed post-adoption window), computed draw by draw; prior and
    posterior may hold different numbers of draws.

    Args:
      prior_traj: Prior trajectory draws, one ``(S, T_post)`` array or ``{contrast: array}``.
      post_traj: Posterior trajectory draws with the same keys (or a single array).
      labels: Display names: ``None``, one string (single-array inputs) or ``{key: name}``.
      quantiles: Strictly increasing probabilities of odd length; the first and last bound
        the reported interval and the middle one is the point summary (default P30/P50/P70).
      ref_total: Control baseline total in outcome units -- a scalar shared by every
        contrast or ``{key: value}`` -- against which the lift is expressed, or ``None``.
      unit: Outcome unit written into the ``unit`` column (display only).
      assert_overlap: Raise ``AssertionError`` when any contrast's prior and posterior
        central intervals do not overlap (a prior-data conflict).

    Returns:
      pandas.DataFrame with one row per contrast and the columns

      ``contrast``
          Display name of the contrast.
      ``horizon_periods``
          Number of post-adoption periods summed into the horizon statistic.
      ``unit``
          The ``unit`` string.
      ``prior_P30`` / ``prior_P50`` / ``prior_P70`` (one per entry of ``quantiles``), ``prior_Mean``, ``prior_SD``
          Prior summaries of the horizon statistic.
      ``post_P30`` / ``post_P50`` / ``post_P70``, ``post_Mean``, ``post_SD``
          Posterior summaries of the horizon statistic.
      ``ref_total``
          Baseline total used for the lift (NaN when not supplied).
      ``prior_lift_pct`` / ``post_lift_pct``
          ``100 * central quantile / ref_total`` (NaN when ``ref_total`` is not supplied).
      ``width_change_pct``
          Percent change of the central-interval width from prior to posterior
          (``100 * (post_width - prior_width) / prior_width``; negative means the data narrowed it).
      ``interval_overlap``
          ``"PASS"`` when the prior and posterior central intervals intersect, else ``"FAIL"``.
    """
    qs = _check_quantiles(quantiles)
    prior = _as_traj_dict(prior_traj, "prior_traj")
    post = _as_traj_dict(post_traj, "post_traj")
    if isinstance(prior_traj, Mapping) != isinstance(post_traj, Mapping):
        raise ValueError("prior_traj and post_traj must both be mappings or both be single arrays")
    if set(prior) != set(post):
        raise ValueError(
            f"prior_traj and post_traj must share their contrast keys; got {sorted(prior)} vs {sorted(post)}"
        )
    keys = list(post)
    names = _resolve_labels(keys, labels)
    refs = _resolve_per_key(ref_total, keys, "ref_total")
    low_name, high_name = _quantile_label(qs[0]), _quantile_label(qs[-1])
    mid_name = _quantile_label(qs[len(qs) // 2])

    rows = []
    failures = []
    for key in keys:
        if prior[key].shape[1] != post[key].shape[1]:
            raise ValueError(
                f"prior and posterior trajectories of {key!r} span different horizons: "
                f"{prior[key].shape[1]} vs {post[key].shape[1]} periods"
            )
        row: dict[str, Any] = {"contrast": names[key], "horizon_periods": int(post[key].shape[1])}
        row["unit"] = unit
        row.update(_summarize(prior[key].sum(axis=1), qs, "prior_"))
        row.update(_summarize(post[key].sum(axis=1), qs, "post_"))
        ref = refs[key]
        row["ref_total"] = float("nan") if ref is None else ref
        row["prior_lift_pct"] = float("nan") if ref is None else 100.0 * row[f"prior_{mid_name}"] / ref
        row["post_lift_pct"] = float("nan") if ref is None else 100.0 * row[f"post_{mid_name}"] / ref
        prior_width = row[f"prior_{high_name}"] - row[f"prior_{low_name}"]
        post_width = row[f"post_{high_name}"] - row[f"post_{low_name}"]
        if prior_width <= 0.0:
            raise ValueError(f"the prior central interval of {key!r} has zero width; the prior draws are degenerate")
        row["width_change_pct"] = 100.0 * (post_width - prior_width) / prior_width
        upper = min(row[f"prior_{high_name}"], row[f"post_{high_name}"])
        lower = max(row[f"prior_{low_name}"], row[f"post_{low_name}"])
        overlap = upper >= lower
        row["interval_overlap"] = "PASS" if overlap else "FAIL"
        if not overlap:
            failures.append(names[key])
        rows.append(row)
    if assert_overlap and failures:
        raise AssertionError(
            f"prior and posterior [{low_name}, {high_name}] intervals do not overlap for {failures}; "
            "the prior and the data conflict -- report the unconstrained fit alongside"
        )
    return pd.DataFrame(rows)


def assert_nested_additivity(
    bc: Any,
    ab: Any,
    ac: Any,
    atol: float = 1e-8,
    *,
    elementwise: bool = True,
) -> None:
    """Assert the structural identity ``AC = BC + AB`` of a nested three-arm design.

    Checks the means of the horizon statistic (sum over periods of each draw) and, when
    ``elementwise`` is set, every draw and period, so a reconciled posterior whose AC
    draws were not built as ``BC + AB`` is caught. Pass ``elementwise=False`` for draws
    that are only coherent in distribution (for example independently sampled priors).

    Args:
      bc: ``(S, T_post)`` trajectory draws of the base-vs-control contrast.
      ab: ``(S, T_post)`` trajectory draws of the incremental-vs-base contrast.
      ac: ``(S, T_post)`` trajectory draws of the combined-vs-control contrast.
      atol: Absolute tolerance of the comparison, in outcome units.
      elementwise: Also compare draw by draw (requires equal shapes).

    Raises:
      AssertionError: With the largest discrepancy when the identity fails.
      ValueError: On malformed or non-finite inputs.
    """
    bc_arr, ab_arr, ac_arr = (_as_draw_matrix(v, n) for v, n in ((bc, "bc"), (ab, "ab"), (ac, "ac")))
    atol = float(atol)
    if not atol >= 0.0:
        raise ValueError(f"atol must be >= 0; got {atol!r}")
    if not bc_arr.shape[1] == ab_arr.shape[1] == ac_arr.shape[1]:
        raise ValueError(
            f"bc, ab and ac must span the same horizon; got {bc_arr.shape}, {ab_arr.shape}, {ac_arr.shape}"
        )
    mean_ac = float(np.mean(ac_arr.sum(axis=1)))
    mean_gap = abs(mean_ac - float(np.mean(bc_arr.sum(axis=1)) + np.mean(ab_arr.sum(axis=1))))
    if mean_gap > atol:
        raise AssertionError(f"mean(AC) differs from mean(BC) + mean(AB) by {mean_gap:.3e} (atol {atol:.1e})")
    if elementwise:
        if not bc_arr.shape == ab_arr.shape == ac_arr.shape:
            raise ValueError(
                f"elementwise additivity needs equal shapes; got {bc_arr.shape}, {ab_arr.shape}, {ac_arr.shape}"
            )
        gap = float(np.max(np.abs(ac_arr - (bc_arr + ab_arr))))
        if gap > atol:
            raise AssertionError(f"AC != BC + AB draw by draw: max |AC - BC - AB| = {gap:.3e} (atol {atol:.1e})")


def nested_3arm_trajectories(
    rec: Nested3ArmReconciliationResult,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Split a reconciliation result into ``(prior, posterior)`` trajectory dicts keyed by contrast.

    Both dicts use the keys of :data:`NESTED_CONTRASTS` (``"BC"``, ``"AB"``, ``"AC"``) and
    feed :func:`prior_vs_posterior_summary` / :func:`plot_rank_k_prior_vs_posterior`
    directly. Requires the result of ``reconcile_nested_3arm_draws(...,
    return_prior_draws=True)``.
    """
    if rec.prior_catt_bc is None or rec.prior_catt_ab is None or rec.prior_catt_ac is None:
        raise ValueError(
            "the reconciliation holds no prior draws; rerun reconcile_nested_3arm_draws(..., return_prior_draws=True)"
        )
    prior = {"BC": rec.prior_catt_bc, "AB": rec.prior_catt_ab, "AC": rec.prior_catt_ac}
    post = {"BC": rec.trt_catt_bc, "AB": rec.trt_catt_ab, "AC": rec.trt_catt_ac}
    return (
        {k: _as_draw_matrix(v, f"prior_catt_{k.lower()}") for k, v in prior.items()},
        {k: _as_draw_matrix(v, f"trt_catt_{k.lower()}") for k, v in post.items()},
    )


# ---- Decay extrapolation ("LongBet Regime 4") ----------------------------------------------------


def _extend_with_decay(observed: np.ndarray, rho: Any, horizon: int) -> np.ndarray:
    """``(S, L)`` observed draws -> ``(S, H)``: observed months, then ``tau(L + h) = rho^h tau(L)``.

    ``rho`` is a scalar or an ``(S,)`` vector (one decay per draw). ``H == L`` returns a copy.
    """
    n_draws, n_obs = observed.shape
    n_tail = int(horizon) - n_obs
    if n_tail < 0:
        raise ValueError(f"horizon ({horizon}) must be >= the number of observed periods ({n_obs})")
    if n_tail == 0:
        return np.array(observed, dtype=np.float64, copy=True)
    rho_arr = np.broadcast_to(np.asarray(rho, dtype=np.float64), (n_draws,))
    steps = np.arange(1, n_tail + 1, dtype=np.float64)
    tail = observed[:, -1][:, None] * np.power(rho_arr[:, None], steps[None, :])
    return np.concatenate([observed, tail], axis=1)


def _beta_from_moments(mean: float, sd: float) -> tuple[float, float]:
    """Method-of-moments ``Beta(a, b)`` with the given mean and SD (requires ``sd^2 < mean (1 - mean)``)."""
    if not 0.0 < mean < 1.0:
        raise ValueError(f"rho must lie strictly inside (0, 1) for a Beta decay prior; got {mean}")
    var = sd * sd
    if not var < mean * (1.0 - mean):
        raise ValueError(
            f"rho_sd={sd} is too large for a Beta with mean {mean}: "
            f"need rho_sd^2 < rho (1 - rho) = {mean * (1 - mean):.4g}"
        )
    kappa = mean * (1.0 - mean) / var - 1.0
    return mean * kappa, (1.0 - mean) * kappa


def _resolve_milestones(milestones: Sequence[int], horizon: int) -> tuple[int, ...]:
    """Positive integer months clipped to ``horizon`` (deduplicated, ascending); ``horizon`` itself is always kept."""
    months = set()
    for m in milestones:
        month = int(m)
        if month != m or month < 1:
            raise ValueError(f"milestones must be positive integers; got {m!r}")
        months.add(min(month, int(horizon)))
    months.add(int(horizon))
    return tuple(sorted(months))


@dataclass(frozen=True)
class DecayExtrapolationResult:
    """Output of :func:`extrapolate_rank_k_decay`; every per-arm mapping is keyed by arm label.

    Attributes:
      arms: Arm labels in input order (a single-array input is the arm ``"effect"``).
      n_observed: Number of observed post-adoption periods ``L``.
      horizon: Extrapolation horizon ``H`` in periods.
      rho: Fixed geometric decay of the point extrapolation.
      rho_sd: SD of the hierarchical ``Beta`` decay, or ``None`` when it was not requested.
      beta_params: ``(a, b)`` of the method-of-moments ``Beta`` decay prior, or ``None``.
      monthly_draws: ``{arm: (S, H)}`` monthly-rate draws (observed months then the decayed tail).
      cumulative_draws: ``{arm: (S, H)}`` running sums of ``monthly_draws``.
      monthly_summary_df: Long table with the columns ``arm``, ``month`` (1-based), ``segment``
        (``"observed"`` / ``"extrapolated"``), ``P30``, ``P50``, ``P70``, ``Mean`` of the monthly rate.
      cumulative_summary_df: Same layout for the cumulative effect, plus ``rho_band_P30`` /
        ``rho_band_P50`` / ``rho_band_P70`` from the hierarchical-``rho`` draws (NaN when
        ``rho_sd`` is ``None``).
      milestone_table_df: One row per arm and milestone month with ``arm``, ``month``, ``segment``,
        ``rate_P30``, ``rate_P50``, ``rate_P70``, ``rate_Mean``, ``cum_P30``, ``cum_P50``, ``cum_P70``,
        ``cum_Mean``.
      data_health_certificate_df: One row per arm and milestone month with ``arm``, ``month``,
        ``n_draws``, ``n_nan``, ``n_inf`` (both asserted to be 0), ``Min``, ``P25``, ``P50``, ``P95``,
        ``Max``, ``Mean`` of the cumulative effect.
      rho_grid_cumulative: ``{rho: {arm: (S, H)}}`` cumulative draws for each decay on ``rho_grid``.
      random_rho_cumulative_draws: ``{arm: (S * num_mc_per_draw, H)}`` cumulative draws under
        ``rho ~ Beta`` (the same decay draws are shared across arms, so additivity survives), or
        ``None``.
    """

    arms: tuple[str, ...]
    n_observed: int
    horizon: int
    rho: float
    rho_sd: float | None
    beta_params: tuple[float, float] | None
    monthly_draws: dict[str, np.ndarray]
    cumulative_draws: dict[str, np.ndarray]
    monthly_summary_df: pd.DataFrame
    cumulative_summary_df: pd.DataFrame
    milestone_table_df: pd.DataFrame
    data_health_certificate_df: pd.DataFrame
    rho_grid_cumulative: dict[float, dict[str, np.ndarray]]
    random_rho_cumulative_draws: dict[str, np.ndarray] | None


def _long_summary(draws: dict[str, np.ndarray], n_observed: int, quantiles: Sequence[float]) -> pd.DataFrame:
    """Long ``arm x month`` table of quantiles and mean of ``(S, H)`` draws, flagging the observed window."""
    rows = []
    for arm, arr in draws.items():
        for j in range(arr.shape[1]):
            segment = "observed" if j < n_observed else "extrapolated"
            row: dict[str, Any] = {"arm": arm, "month": j + 1, "segment": segment}
            summary = _summarize(arr[:, j], quantiles)
            summary.pop("SD")
            row.update(summary)
            rows.append(row)
    return pd.DataFrame(rows)


def extrapolate_rank_k_decay(
    traj_draws: Any,
    *,
    horizon: int = 60,
    rho: float = 0.85,
    rho_sd: float | None = 0.07,
    rho_grid: Sequence[float] = (0.75, 0.80, 0.85, 0.90, 0.95),
    num_mc_per_draw: int = 25,
    seed: int = 2026,
    milestones: Sequence[int] = DEFAULT_MILESTONES,
    quantiles: Sequence[float] = DEFAULT_QUANTILES,
) -> DecayExtrapolationResult:
    """Extend observed trajectory draws to a longer horizon with a geometric decay of the last observed rate.

    For every draw ``s`` the observed monthly rates ``tau^(s)(1..L)`` are kept and the tail
    is ``tau^(s)(L + h) = rho^h tau^(s)(L)``; the cumulative effect is the running sum.
    The hierarchical variant draws ``rho ~ Beta(a, b)`` with mean ``rho`` and SD ``rho_sd``
    (method of moments), ``num_mc_per_draw`` decays per posterior draw, shared across arms
    so ``AC = BC + AB`` survives the extrapolation. The ``rho_grid`` sensitivity repeats the
    point extrapolation for each decay in the grid.

    Two invariants are asserted after the computation: the cumulative effect is
    non-decreasing wherever the monthly rates are non-negative, and the horizon-``L``
    cumulative equals the observed one to ``1e-10``.

    Args:
      traj_draws: One ``(S, L)`` array (arm label ``"effect"``) or ``{arm: (S, L)}`` with a
        common ``S`` and ``L``.
      horizon: Extrapolation horizon ``H >= L`` in periods.
      rho: Fixed decay in ``[0, 1]`` (strictly inside ``(0, 1)`` when ``rho_sd`` is set).
      rho_sd: SD of the ``Beta`` decay prior, or ``None`` / ``0`` to skip the hierarchical band.
      rho_grid: Decays of the sensitivity grid (each in ``[0, 1]``).
      num_mc_per_draw: Decay draws per posterior draw in the hierarchical variant.
      seed: Seed of the decay draws.
      milestones: Months to tabulate; clipped to ``horizon`` and always including it.
      quantiles: Probabilities of the summary tables (odd length, default P30/P50/P70).

    Returns:
      :class:`DecayExtrapolationResult`.
    """
    draws = _as_traj_dict(traj_draws, "traj_draws")
    qs = _check_quantiles(quantiles)
    shapes = {arr.shape for arr in draws.values()}
    if len(shapes) != 1:
        raise ValueError(f"every arm must have the same (S, L) shape; got {sorted(shapes)}")
    n_draws, n_obs = next(iter(shapes))
    horizon = int(horizon)
    if horizon < n_obs:
        raise ValueError(f"horizon ({horizon}) must be >= the number of observed periods ({n_obs})")
    rho = float(rho)
    if not 0.0 <= rho <= 1.0:
        raise ValueError(f"rho must lie in [0, 1]; got {rho}")
    grid = tuple(float(r) for r in rho_grid)
    if any(not 0.0 <= r <= 1.0 for r in grid):
        raise ValueError(f"rho_grid entries must lie in [0, 1]; got {grid}")
    n_mc = int(num_mc_per_draw)
    if n_mc < 1:
        raise ValueError(f"num_mc_per_draw must be >= 1; got {num_mc_per_draw!r}")
    months = _resolve_milestones(milestones, horizon)

    monthly = {arm: _extend_with_decay(arr, rho, horizon) for arm, arr in draws.items()}
    cumulative = {arm: np.cumsum(arr, axis=1) for arm, arr in monthly.items()}
    grid_cumulative = {
        r: {arm: np.cumsum(_extend_with_decay(arr, r, horizon), axis=1) for arm, arr in draws.items()} for r in grid
    }

    sd = None if rho_sd is None or float(rho_sd) == 0.0 else float(rho_sd)
    beta_params = None
    random_cumulative = None
    if sd is not None:
        if sd < 0.0:
            raise ValueError(f"rho_sd must be >= 0 or None; got {rho_sd!r}")
        beta_params = _beta_from_moments(rho, sd)
        rng = np.random.default_rng(int(seed))
        rho_draws = rng.beta(beta_params[0], beta_params[1], size=n_draws * n_mc)
        random_cumulative = {
            arm: np.cumsum(_extend_with_decay(np.repeat(arr, n_mc, axis=0), rho_draws, horizon), axis=1)
            for arm, arr in draws.items()
        }

    for arm in draws:
        obs_cum = np.cumsum(draws[arm], axis=1)[:, -1]
        if np.max(np.abs(cumulative[arm][:, n_obs - 1] - obs_cum)) > 1e-10:
            raise AssertionError(f"horizon-L cumulative of {arm!r} differs from the observed cumulative")
        if np.all(monthly[arm] >= 0.0) and np.any(np.diff(cumulative[arm], axis=1) < -1e-12):
            raise AssertionError(f"cumulative effect of {arm!r} decreases although every monthly rate is >= 0")

    monthly_df = _long_summary(monthly, n_obs, qs)
    cumulative_df = _long_summary(cumulative, n_obs, qs)
    for q in qs:
        label = _quantile_label(q)
        cumulative_df[f"rho_band_{label}"] = float("nan")
        if random_cumulative is not None:
            for arm in draws:
                mask = cumulative_df["arm"] == arm
                cumulative_df.loc[mask, f"rho_band_{label}"] = np.quantile(random_cumulative[arm], q, axis=0)

    milestone_rows = []
    health_rows = []
    for arm in draws:
        for month in months:
            rate = monthly[arm][:, month - 1]
            cum = cumulative[arm][:, month - 1]
            segment = "observed" if month <= n_obs else "extrapolated"
            row: dict[str, Any] = {"arm": arm, "month": month, "segment": segment}
            row.update({f"rate_{k}": v for k, v in _summarize(rate, qs).items() if k != "SD"})
            row.update({f"cum_{k}": v for k, v in _summarize(cum, qs).items() if k != "SD"})
            milestone_rows.append(row)
            n_nan, n_inf = int(np.sum(np.isnan(cum))), int(np.sum(np.isinf(cum)))
            if n_nan or n_inf:
                raise AssertionError(
                    f"cumulative effect of {arm!r} at month {month} has {n_nan} NaN and {n_inf} Inf draws"
                )
            health_rows.append(
                {
                    "arm": arm, "month": month, "n_draws": int(cum.shape[0]), "n_nan": n_nan, "n_inf": n_inf,
                    "Min": float(np.min(cum)), "P25": float(np.quantile(cum, 0.25)),
                    "P50": float(np.quantile(cum, 0.50)), "P95": float(np.quantile(cum, 0.95)),
                    "Max": float(np.max(cum)), "Mean": float(np.mean(cum)),
                }
            )
    return DecayExtrapolationResult(
        arms=tuple(draws), n_observed=int(n_obs), horizon=horizon, rho=rho, rho_sd=sd, beta_params=beta_params,
        monthly_draws=monthly, cumulative_draws=cumulative, monthly_summary_df=monthly_df,
        cumulative_summary_df=cumulative_df, milestone_table_df=pd.DataFrame(milestone_rows),
        data_health_certificate_df=pd.DataFrame(health_rows), rho_grid_cumulative=grid_cumulative,
        random_rho_cumulative_draws=random_cumulative,
    )


# ---- Prior sensitivity across tiers --------------------------------------------------------------


def standardized_prior_family_table(
    tiers: Sequence[PriorTier],
    *,
    quantiles: Sequence[float] = DEFAULT_QUANTILES,
    thresholds: Sequence[float] = (0.25, 0.50),
    realized_scales: Mapping[str, float] | None = None,
) -> pd.DataFrame:
    """Analytic summaries of each tier's prior family in units of its scale.

    Every family is a scale family, so its quantiles, moments and near-zero mass are
    computed once at unit scale with scipy: in units of ``sigma_hist`` for the numeric
    tiers and in units of the fit's own posterior SD for an ``"in_sample"`` tier.

    Args:
      tiers: The tiers to tabulate (unique names).
      quantiles: Probabilities to report (odd length; default P30/P50/P70).
      thresholds: Multiples of the scale whose prior mass ``Pr(tau < c * scale)`` is reported.
      realized_scales: Optional ``{tier name: scale in outcome units}`` (the reconciliation's
        realized scale, which resolves ``"in_sample"``) written into ``scale_value``.

    Returns:
      pandas.DataFrame with one row per tier and the columns ``tier``, ``prior_mode``,
      ``family``, ``scale_unit`` (``"sigma_hist"`` or ``"in-sample posterior SD"``),
      ``scale_value`` (outcome units; NaN when unknown), one column per quantile
      (``P30``, ``P50``, ``P70``), ``Mean``, ``SD``, ``Var`` and one
      ``pr_below_<c>_scale`` column per threshold (for example ``pr_below_0.25_scale``).
    """
    qs = _check_quantiles(quantiles)
    names = [t.name for t in tiers]
    if not names:
        raise ValueError("tiers must hold at least one PriorTier")
    if len(set(names)) != len(names):
        raise ValueError(f"tier names must be unique; got {names}")
    cuts = tuple(float(c) for c in thresholds)
    if any(not c > 0.0 for c in cuts):
        raise ValueError(f"thresholds must be positive multiples of the scale; got {cuts}")
    rows = []
    for tier in tiers:
        dist = _prior_family_distribution(tier.family, 1.0)
        row: dict[str, Any] = {
            "tier": tier.name,
            "prior_mode": tier.prior_mode,
            "family": _FAMILY_DISPLAY[tier.family],
            "scale_unit": "in-sample posterior SD" if tier.is_in_sample else "sigma_hist",
            "scale_value": float("nan"),
        }
        if realized_scales is not None and tier.name in realized_scales:
            row["scale_value"] = float(realized_scales[tier.name])
        elif not tier.is_in_sample:
            row["scale_value"] = float(tier.prior_scale)
        for q in qs:
            row[_quantile_label(q)] = float(dist.ppf(q))
        row["Mean"] = float(dist.mean())
        row["SD"] = float(dist.std())
        row["Var"] = float(dist.std()) ** 2
        for c in cuts:
            row[f"pr_below_{c:g}_scale"] = float(dist.cdf(c))
        rows.append(row)
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class PriorSensitivityResult:
    """Output of :func:`evaluate_prior_sensitivity`.

    Attributes:
      tiers: The evaluated tiers, in order.
      rho: Geometric decay used for the extrapolated horizon.
      extrap_horizon: Extrapolated horizon ``H`` in periods.
      n_observed: Observed post-adoption periods ``L``.
      cost_grid: Per-unit costs of the break-even table (outcome units).
      alpha_grid: Incrementality shares of the break-even table.
      standardized_families_df: :func:`standardized_prior_family_table` of the tiers.
      posterior_comparison_df: Long table with one row per tier, contrast and horizon:
        ``tier``, ``prior_mode``, ``prior_scale`` (realized, outcome units), ``sir_ess``,
        ``contrast`` (``"BC"`` / ``"AB"`` / ``"AC"``), ``horizon`` (``"observed"`` /
        ``"extrapolated"``), ``horizon_periods``, ``P30``, ``P50``, ``P70``, ``Mean``, ``SD`` of
        the cumulative effect.
      roi_break_even_df: One row per tier, contrast, cost and alpha: ``tier``, ``contrast``,
        ``horizon_periods``, ``cost``, ``alpha``, ``pr_roi_positive``
        (``Pr(alpha * cumulative_H > cost)`` over the draws) and ``break_even_alpha_50``
        (``cost / median(cumulative_H)``: the incrementality at which ``Pr(ROI > 0) = 0.5``).
      draws_by_tier: ``{tier name: Nested3ArmReconciliationResult}`` with prior draws attached.
      observed_cumulative_by_tier: ``{tier name: {contrast: (S,)}}`` sums over the observed window.
      extrapolated_cumulative_by_tier: ``{tier name: {contrast: (S,)}}`` cumulative at ``H``.
    """

    tiers: tuple[PriorTier, ...]
    rho: float
    extrap_horizon: int
    n_observed: int
    cost_grid: tuple[float, ...]
    alpha_grid: tuple[float, ...]
    standardized_families_df: pd.DataFrame
    posterior_comparison_df: pd.DataFrame
    roi_break_even_df: pd.DataFrame
    draws_by_tier: dict[str, Nested3ArmReconciliationResult]
    observed_cumulative_by_tier: dict[str, dict[str, np.ndarray]]
    extrapolated_cumulative_by_tier: dict[str, dict[str, np.ndarray]]


def roi_break_even_table(
    cumulative_draws: Mapping[str, np.ndarray],
    *,
    cost_grid: Sequence[float],
    alpha_grid: Sequence[float],
    horizon_periods: int,
    tier: str = "",
) -> pd.DataFrame:
    """Commercial break-even of a cumulative effect: ``ROI > 0`` iff ``alpha * effect > cost``.

    Args:
      cumulative_draws: ``{contrast: (S,)}`` cumulative-effect draws at the decision horizon.
      cost_grid: Per-unit costs (outcome units, ``> 0``).
      alpha_grid: Incrementality shares in ``[0, 1]`` (the fraction of the effect that is causal revenue).
      horizon_periods: Horizon of the cumulative effect, written into the table.
      tier: Tier name written into the ``tier`` column.

    Returns:
      pandas.DataFrame with the columns ``tier``, ``contrast``, ``horizon_periods``, ``cost``,
      ``alpha``, ``pr_roi_positive`` and ``break_even_alpha_50`` (see
      :class:`PriorSensitivityResult`). ``break_even_alpha_50`` is ``inf`` when the median
      cumulative effect is not positive.
    """
    costs = tuple(_positive_scalar(c, "cost_grid entry") for c in cost_grid)
    alphas = tuple(float(a) for a in alpha_grid)
    if not costs or not alphas:
        raise ValueError("cost_grid and alpha_grid must be non-empty")
    if any(not 0.0 <= a <= 1.0 for a in alphas):
        raise ValueError(f"alpha_grid entries must lie in [0, 1]; got {alphas}")
    rows = []
    for contrast, draws in cumulative_draws.items():
        arr = np.asarray(draws, dtype=np.float64)
        if arr.ndim != 1 or arr.size < 2 or not np.all(np.isfinite(arr)):
            raise ValueError(f"cumulative_draws[{contrast!r}] must be a finite 1-D array of draws")
        median = float(np.median(arr))
        for cost in costs:
            alpha_50 = cost / median if median > 0.0 else float("inf")
            for alpha in alphas:
                rows.append(
                    {
                        "tier": tier, "contrast": contrast, "horizon_periods": int(horizon_periods), "cost": cost,
                        "alpha": alpha, "pr_roi_positive": float(np.mean(alpha * arr > cost)),
                        "break_even_alpha_50": alpha_50,
                    }
                )
    return pd.DataFrame(rows)


def _reconcile_tier(
    tier: PriorTier,
    res_bc: RankKLongBetResult,
    res_ab: RankKLongBetResult,
    res_ac: RankKLongBetResult,
    *,
    seed: int,
    n_proposals: int | None,
    reconcile_kwargs: Mapping[str, Any],
) -> Nested3ArmReconciliationResult:
    """Run the reconciliation for one tier, silencing only the ``in_sample`` scale warning of that tier."""
    with warnings.catch_warnings():
        if tier.is_in_sample:
            warnings.filterwarnings("ignore", message=r".*prior_scale='in_sample'.*", category=RankKConfigWarning)
        return reconcile_nested_3arm_draws(
            res_bc, res_ab, res_ac, prior_mode=tier.prior_mode, prior_scale=tier.prior_scale, seed=seed,
            n_proposals=n_proposals, return_prior_draws=True, **dict(reconcile_kwargs),
        )


def evaluate_prior_sensitivity(
    res_bc: RankKLongBetResult,
    res_ab: RankKLongBetResult,
    res_ac: RankKLongBetResult,
    *,
    tiers: Sequence[PriorTier],
    rho: float = 0.85,
    extrap_horizon: int = 12,
    cost_grid: Sequence[float] = (1.10, 2.85),
    alpha_grid: Sequence[float] = (0.50, 1.00),
    seed: int = 2026,
    n_proposals: int | None = None,
    quantiles: Sequence[float] = DEFAULT_QUANTILES,
    reconcile_kwargs: Mapping[str, Any] | None = None,
) -> PriorSensitivityResult:
    """Rerun the nested three-arm reconciliation under several prior tiers and tabulate the decision quantities.

    For each tier the reconciliation is run with ``prior_mode=tier.prior_mode`` and
    ``prior_scale=tier.prior_scale`` (prior draws attached); the ``in_sample`` warning of
    the reconciliation is silenced for that tier's call only. The observed-horizon
    statistic is the sum of each posterior trajectory over its ``L`` periods and the
    extrapolated-horizon statistic is the cumulative effect at ``extrap_horizon`` under the
    geometric decay ``rho`` of :func:`extrapolate_rank_k_decay`.

    Args:
      res_bc, res_ab, res_ac: The three pairwise fits (base vs control, incremental vs base,
        combined vs control).
      tiers: Prior tiers, for example :func:`ytlr_prior_tiers`.
      rho: Geometric decay of the extrapolated horizon, in ``[0, 1]``.
      extrap_horizon: Decision horizon ``H >= L`` in periods.
      cost_grid: Per-unit costs of the break-even table.
      alpha_grid: Incrementality shares of the break-even table.
      seed: Seed shared by every tier's reconciliation.
      n_proposals: Importance-sampling proposals per tier (``None`` for the reconciliation default).
      quantiles: Probabilities of the summary tables (odd length).
      reconcile_kwargs: Extra keyword arguments forwarded to :func:`reconcile_nested_3arm_draws`
        (for example ``tau_clip``, ``catt_floor`` or explicit cohorts).

    Returns:
      :class:`PriorSensitivityResult`.
    """
    tiers = tuple(tiers)
    qs = _check_quantiles(quantiles)
    families_df = standardized_prior_family_table(tiers, quantiles=qs)  # validates names and modes
    extra = dict(reconcile_kwargs or {})
    for forbidden in ("prior_mode", "prior_scale", "seed", "n_proposals", "return_prior_draws"):
        if forbidden in extra:
            raise ValueError(f"reconcile_kwargs must not set {forbidden!r}; it is controlled by the tiers / arguments")
    rho = float(rho)
    if not 0.0 <= rho <= 1.0:
        raise ValueError(f"rho must lie in [0, 1]; got {rho}")
    horizon = int(extrap_horizon)

    draws_by_tier: dict[str, Nested3ArmReconciliationResult] = {}
    observed: dict[str, dict[str, np.ndarray]] = {}
    extrapolated: dict[str, dict[str, np.ndarray]] = {}
    realized_scales: dict[str, float] = {}
    comparison_rows = []
    roi_frames = []
    n_obs = None
    for tier in tiers:
        rec = _reconcile_tier(tier, res_bc, res_ab, res_ac, seed=seed, n_proposals=n_proposals, reconcile_kwargs=extra)
        _, post = nested_3arm_trajectories(rec)
        n_obs = int(post["AC"].shape[1])
        if horizon < n_obs:
            raise ValueError(f"extrap_horizon ({horizon}) must be >= the observed post-adoption length ({n_obs})")
        draws_by_tier[tier.name] = rec
        scale = rec.prior_scale if isinstance(rec.prior_scale, (int, float)) else float("nan")
        realized_scales[tier.name] = float(scale)
        observed[tier.name] = {k: v.sum(axis=1) for k, v in post.items()}
        extrapolated[tier.name] = {k: _extend_with_decay(v, rho, horizon).sum(axis=1) for k, v in post.items()}
        for contrast in NESTED_CONTRASTS:
            for kind, stat, periods in (("observed", observed, n_obs), ("extrapolated", extrapolated, horizon)):
                row: dict[str, Any] = {
                    "tier": tier.name, "prior_mode": tier.prior_mode, "prior_scale": float(scale),
                    "sir_ess": float(rec.sir_ess), "contrast": contrast, "horizon": kind, "horizon_periods": periods,
                }
                row.update(_summarize(stat[tier.name][contrast], qs))
                comparison_rows.append(row)
        roi_frames.append(
            roi_break_even_table(
                extrapolated[tier.name], cost_grid=cost_grid, alpha_grid=alpha_grid, horizon_periods=horizon,
                tier=tier.name,
            )
        )
    families_df = standardized_prior_family_table(tiers, quantiles=qs, realized_scales=realized_scales)
    return PriorSensitivityResult(
        tiers=tiers, rho=rho, extrap_horizon=horizon, n_observed=int(n_obs or 0),
        cost_grid=tuple(float(c) for c in cost_grid), alpha_grid=tuple(float(a) for a in alpha_grid),
        standardized_families_df=families_df, posterior_comparison_df=pd.DataFrame(comparison_rows),
        roi_break_even_df=pd.concat(roi_frames, ignore_index=True), draws_by_tier=draws_by_tier,
        observed_cumulative_by_tier=observed, extrapolated_cumulative_by_tier=extrapolated,
    )


# ---- Figures (matplotlib imported lazily) --------------------------------------------------------


def _pyplot(func_name: str) -> Any:
    """Import ``matplotlib.pyplot`` on demand so matplotlib stays an optional dependency."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            f"{func_name}() needs matplotlib. Install it with 'pip install matplotlib' (or longbet-jax[plot]), "
            "or use the data-layer function and plot the table with whatever you already have."
        ) from exc
    return plt


def _axes_array(plt: Any, axes: Any, n_panels: int, figsize: tuple[float, float]) -> np.ndarray:
    """Return ``n_panels`` Axes as a flat object array, creating a one-row figure when ``axes`` is ``None``."""
    if axes is None:
        _, created = plt.subplots(1, n_panels, figsize=figsize)
        return np.atleast_1d(created).ravel()
    arr = np.atleast_1d(np.asarray(axes, dtype=object)).ravel()
    if arr.size != n_panels:
        raise ValueError(f"expected {n_panels} Axes; got {arr.size}")
    return arr


def _density_curve(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Gaussian KDE of ``values`` on ``grid``; degenerate samples fall back to a normalized histogram."""
    try:
        return np.asarray(sp_stats.gaussian_kde(values)(grid), dtype=np.float64)
    except (np.linalg.LinAlgError, ValueError):
        n_bins = min(30, max(5, values.size // 4))
        counts, edges = np.histogram(values, bins=n_bins, range=(grid[0], grid[-1]), density=True)
        idx = np.clip(np.searchsorted(edges, grid, side="right") - 1, 0, counts.size - 1)
        return counts[idx].astype(np.float64)


def _unit_suffix(unit: str) -> str:
    return f" ({unit})" if unit else ""


def _draw_forest(
    ax: Any, summary: pd.DataFrame, low: str, mid: str, high: str, ref_totals: Mapping[str, float | None], unit: str
) -> None:
    """Prior (grey) and posterior (blue) central intervals per contrast, one pair of rows per contrast."""
    n = len(summary)
    y_post = np.arange(n, dtype=np.float64) * 2.0 + 0.3
    y_prior = y_post - 0.6
    for i, row in enumerate(summary.itertuples(index=False)):
        rec = row._asdict()
        prior_span = [rec[f"prior_{low}"], rec[f"prior_{high}"]]
        post_span = [rec[f"post_{low}"], rec[f"post_{high}"]]
        ax.plot(prior_span, [y_prior[i]] * 2, color="#9e9e9e", lw=4, solid_capstyle="butt")
        ax.plot(rec[f"prior_{mid}"], y_prior[i], "o", color="#616161", ms=6)
        ax.plot(post_span, [y_post[i]] * 2, color=_ARM_COLORS[0], lw=4, solid_capstyle="butt")
        ax.plot(rec[f"post_{mid}"], y_post[i], "o", color=_ARM_COLORS[0], ms=6)
    ticks, labels = [], []
    for i, row in enumerate(summary.itertuples(index=False)):
        rec = row._asdict()
        ref = ref_totals.get(rec["contrast"])
        lift = ""
        if ref is not None and math.isfinite(rec["post_lift_pct"]):
            lift = f"\n(+{rec['post_lift_pct']:.1f}% of baseline)"
        ticks.append(0.5 * (y_post[i] + y_prior[i]))
        labels.append(f"{rec['contrast']}{lift}")
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_ylim(-1.0, 2.0 * n)
    ax.axvline(0.0, color="black", lw=0.8, alpha=0.5)
    ax.set_xlabel(f"Horizon effect{_unit_suffix(unit)}")
    ax.plot([], [], color="#9e9e9e", lw=4, label=f"Prior [{low}, {high}]")
    ax.plot([], [], color=_ARM_COLORS[0], lw=4, label=f"Posterior [{low}, {high}]")
    ax.legend(frameon=False, loc="best")


def plot_rank_k_prior_vs_posterior(
    prior_traj: Any,
    post_traj: Any,
    *,
    labels: Any = None,
    quantiles: Sequence[float] = DEFAULT_QUANTILES,
    unit: str = "",
    ref_totals: Any = None,
    kind: str = "forest_and_density",
    ax: Any = None,
    title: str | None = None,
) -> Any:
    """Draw prior and posterior horizon effects per contrast as a forest plot and / or densities.

    ``matplotlib`` is imported here rather than at module load. The numbers come from
    :func:`prior_vs_posterior_summary`.

    Args:
      prior_traj, post_traj, labels, quantiles, unit: As in :func:`prior_vs_posterior_summary`.
      ref_totals: Control baseline totals (scalar or ``{key: value}``) used to annotate the
        posterior lift in the forest tick labels, or ``None``.
      kind: ``"forest"`` (one Axes), ``"density"`` (one Axes) or ``"forest_and_density"`` (two Axes).
      ax: Axes to draw on -- one Axes for the single-panel kinds, a sequence of two for
        ``"forest_and_density"`` -- or ``None`` to create a figure.
      title: Figure-level title (suptitle for two panels), or ``None``.

    Returns:
      A single ``matplotlib.axes.Axes`` for the one-panel kinds, else a length-2 object array of Axes.
    """
    plt = _pyplot("plot_rank_k_prior_vs_posterior")
    if kind not in ("forest", "density", "forest_and_density"):
        raise ValueError(f"kind must be 'forest', 'density' or 'forest_and_density'; got {kind!r}")
    qs = _check_quantiles(quantiles)
    summary = prior_vs_posterior_summary(
        prior_traj, post_traj, labels=labels, quantiles=qs, ref_total=ref_totals, unit=unit
    )
    prior = _as_traj_dict(prior_traj, "prior_traj")
    post = _as_traj_dict(post_traj, "post_traj")
    names = _resolve_labels(list(post), labels)
    refs = {names[k]: v for k, v in _resolve_per_key(ref_totals, list(post), "ref_totals").items()}
    low, mid, high = _quantile_label(qs[0]), _quantile_label(qs[len(qs) // 2]), _quantile_label(qs[-1])

    n_panels = 2 if kind == "forest_and_density" else 1
    axes = _axes_array(plt, ax, n_panels, (6.0 * n_panels + 1.0, 4.2))
    panel = 0
    if kind in ("forest", "forest_and_density"):
        _draw_forest(axes[panel], summary, low, mid, high, refs, unit)
        axes[panel].set_title("Central intervals of the horizon effect")
        panel += 1
    if kind in ("density", "forest_and_density"):
        dax = axes[panel]
        totals = {k: (prior[k].sum(axis=1), post[k].sum(axis=1)) for k in post}
        lo = min(float(min(p.min(), q.min())) for p, q in totals.values())
        hi = max(float(max(p.max(), q.max())) for p, q in totals.values())
        pad = 0.05 * (hi - lo if hi > lo else 1.0)
        grid = np.linspace(lo - pad, hi + pad, 400)
        for i, key in enumerate(post):
            color = _ARM_COLORS[i % len(_ARM_COLORS)]
            dax.plot(grid, _density_curve(totals[key][0], grid), color=color, ls="--", lw=1.2, alpha=0.8)
            dax.fill_between(grid, _density_curve(totals[key][1], grid), color=color, alpha=0.25, label=names[key])
        dax.plot([], [], color="black", ls="--", lw=1.2, label="Prior (dashed)")
        dax.set_xlabel(f"Horizon effect{_unit_suffix(unit)}")
        dax.set_ylabel("Density")
        dax.set_title("Prior (dashed) vs posterior (filled)")
        dax.legend(frameon=False)
    if title:
        if n_panels == 2:
            axes[0].figure.suptitle(title)
        else:
            axes[0].set_title(title)
    return axes if n_panels == 2 else axes[0]


def _two_halfnormal_sum_pdf(x: np.ndarray, scale: float) -> np.ndarray:
    """Density of ``X + Y`` with ``X, Y`` iid ``N^+(0, s^2)``: ``2 / (sqrt(pi) s) exp(-x^2 / 4 s^2) erf(x / 2 s)``."""
    x = np.asarray(x, dtype=np.float64)
    val = 2.0 / (math.sqrt(math.pi) * scale) * np.exp(-(x**2) / (4.0 * scale**2))
    val = val * sp_special.erf(x / (2.0 * scale))
    return np.where(x >= 0.0, val, 0.0)


def plot_rank_k_prior_calibration(
    prior_scale: float,
    *,
    hist_mean: float | None = None,
    hist_se: float | None = None,
    prior_mode: str = "bundle_halfnormal",
    split_share: float = 0.5,
    horizon_periods: int = 1,
    unit: str = "",
    show_rayleigh_comparison: bool = True,
    ax: Any = None,
    title: str | None = None,
) -> Any:
    """Show how the zero-mode bundle prior on the total effect relates to a historical benchmark.

    Panel A draws the prior density of the total bundle effect over ``horizon_periods``
    (the family of ``prior_mode`` with scale ``horizon_periods * prior_scale``, since the
    prior is on the per-period level), shades ``[P30, P70]``, marks P50 and the implied
    single-component P50 (``split_share * P50``), and places the historical benchmark
    ``hist_mean +/- hist_se`` with the percent discount of the prior P50 relative to it.
    Panel B (``show_rayleigh_comparison``) contrasts the bundle construction
    ``tau_AC ~ N^+(0, s^2)`` split by ``u ~ U(0, 1)`` -- whose total is the half-normal
    itself -- with the Rayleigh-shaped density of a sum of two independent ``N^+(0, s^2)``
    components (the total implied by independent half-normal priors on ``BC`` and ``AB``,
    which vanishes at zero and roughly doubles the prior mean).

    Args:
      prior_scale: Per-period bundle prior scale ``s`` in outcome units.
      hist_mean: Historical benchmark point estimate of the horizon total, or ``None``.
      hist_se: Its standard error (needs ``hist_mean``), or ``None``.
      prior_mode: One of :data:`BUNDLE_PRIOR_MODES`.
      split_share: Share of the bundle attributed to one component for the implied P50 marker.
      horizon_periods: Periods summed into the horizon total (``>= 1``).
      unit: Outcome unit for the axis label.
      show_rayleigh_comparison: Add panel B.
      ax: One Axes (panel A only) or a sequence of two Axes, or ``None`` to create a figure.
      title: Figure-level title, or ``None``.

    Returns:
      A single Axes without panel B, else a length-2 object array of Axes.
    """
    plt = _pyplot("plot_rank_k_prior_calibration")
    scale = _positive_scalar(prior_scale, "prior_scale")
    family = _family_of_prior_mode(prior_mode)
    periods = int(horizon_periods)
    if periods < 1:
        raise ValueError(f"horizon_periods must be >= 1; got {horizon_periods!r}")
    if not 0.0 <= float(split_share) <= 1.0:
        raise ValueError(f"split_share must lie in [0, 1]; got {split_share!r}")
    if hist_se is not None and hist_mean is None:
        raise ValueError("hist_se needs hist_mean")
    total_scale = periods * scale
    dist = _prior_family_distribution(family, total_scale)
    p30, p50, p70 = (float(dist.ppf(q)) for q in (0.30, 0.50, 0.70))
    x_max = float(dist.ppf(0.995))
    if hist_mean is not None:
        x_max = max(x_max, float(hist_mean) + 3.0 * float(hist_se or 0.0))
    grid = np.linspace(0.0, 1.05 * x_max, 500)

    n_panels = 2 if show_rayleigh_comparison else 1
    axes = _axes_array(plt, ax, n_panels, (6.0 * n_panels + 1.0, 4.2))
    pax = axes[0]
    density = np.asarray(dist.pdf(grid), dtype=np.float64)
    finite = np.isfinite(density)
    pax.plot(
        grid[finite], density[finite], color=_ARM_COLORS[0], lw=1.8,
        label=f"{_FAMILY_DISPLAY[family]} prior, s = {scale:g}",
    )
    band = (grid >= p30) & (grid <= p70) & finite
    pax.fill_between(
        grid[band], density[band], color=_ARM_COLORS[0], alpha=0.25, label=f"[P30, P70] = [{p30:.3g}, {p70:.3g}]"
    )
    pax.axvline(p50, color=_ARM_COLORS[0], ls="-", lw=1.2, label=f"P50 = {p50:.3g}")
    pax.axvline(
        float(split_share) * p50, color=_ARM_COLORS[2], ls=":", lw=1.4,
        label=f"Implied single-component P50 = {split_share * p50:.3g}",
    )
    if hist_mean is not None:
        h_mean = float(hist_mean)
        discount = 100.0 * (p50 - h_mean) / h_mean if h_mean != 0.0 else float("nan")
        y_mark = 0.5 * float(np.nanmax(density[finite])) if finite.any() else 1.0
        se_text = "" if hist_se is None else f" +/- {float(hist_se):.2g}"
        pax.errorbar(
            [h_mean], [y_mark], xerr=None if hist_se is None else [float(hist_se)], fmt="D", color=_ARM_COLORS[1],
            capsize=4, label=f"Historical benchmark {h_mean:.3g}{se_text}",
        )
        pax.annotate(
            f"prior P50 vs benchmark: {discount:+.0f}%", xy=(h_mean, y_mark), xytext=(0, 12),
            textcoords="offset points", ha="center", fontsize=9, color=_ARM_COLORS[1],
        )
    pax.set_xlabel(f"Total bundle effect over {periods} period(s){_unit_suffix(unit)}")
    pax.set_ylabel("Prior density")
    pax.set_title("A. Zero-mode bundle prior on the total effect")
    pax.legend(frameon=False, fontsize=8)
    if show_rayleigh_comparison:
        rax = axes[1]
        hn = sp_stats.halfnorm(scale=total_scale)
        rax.plot(
            grid, hn.pdf(grid), color=_ARM_COLORS[0], lw=1.8, label="Bundle: tau_AC ~ N+(0, s^2), u ~ U(0, 1) split"
        )
        rax.plot(
            grid, _two_halfnormal_sum_pdf(grid, total_scale), color=_ARM_COLORS[1], lw=1.8, ls="--",
            label="Sum of two independent N+(0, s^2) (Rayleigh-shaped)",
        )
        rax.axvline(float(hn.ppf(0.5)), color=_ARM_COLORS[0], lw=1.0, ls=":")
        rax.set_xlabel(f"Total effect{_unit_suffix(unit)}")
        rax.set_ylabel("Prior density")
        rax.set_title("B. Bundle prior vs independent-component prior")
        rax.legend(frameon=False, fontsize=8)
    if title:
        if n_panels == 2:
            axes[0].figure.suptitle(title)
        else:
            axes[0].set_title(title)
    return axes if n_panels == 2 else axes[0]


def plot_rank_k_prior_sensitivity(
    sens_result: PriorSensitivityResult,
    *,
    metric_label: str = "",
    unit: str = "",
    cost_highlight: float = 1.10,
    contrast: str = "AC",
    axes: Any = None,
    title: str | None = None,
) -> np.ndarray:
    """Three-panel prior-sensitivity figure: standardized priors, posterior comparison, break-even curves.

    Panel A draws each distinct prior family at unit scale with its near-zero mass
    annotations (``Pr(tau < 0.25 s)``, ``Pr(tau < 0.50 s)``); panel B is a forest of the
    posterior P50 with ``[P30, P70]`` per tier for the observed and extrapolated horizons of
    ``contrast``; panel C draws ``Pr(ROI > 0 | alpha, cost_highlight)`` against the
    incrementality ``alpha`` for each tier, with the break-even ``alpha_50`` markers where
    the curves cross one half.

    Args:
      sens_result: Output of :func:`evaluate_prior_sensitivity`.
      metric_label: Name of the outcome for the axis labels.
      unit: Outcome unit for the axis labels.
      cost_highlight: Per-unit cost of the break-even curves (outcome units, ``> 0``).
      contrast: Contrast shown in panels B and C (``"BC"``, ``"AB"`` or ``"AC"``).
      axes: Sequence of three Axes, or ``None`` to create a one-row figure.
      title: Figure-level title, or ``None``.

    Returns:
      Length-3 object array of Axes.
    """
    plt = _pyplot("plot_rank_k_prior_sensitivity")
    cost = _positive_scalar(cost_highlight, "cost_highlight")
    if contrast not in NESTED_CONTRASTS:
        raise ValueError(f"contrast must be one of {NESTED_CONTRASTS}; got {contrast!r}")
    axs = _axes_array(plt, axes, 3, (18.0, 4.6))
    label = metric_label or "effect"

    fam_df = sens_result.standardized_families_df
    pr_cols = [c for c in fam_df.columns if c.startswith("pr_below_")]
    grid = np.linspace(1e-3, 3.0, 400)
    seen: dict[str, list[str]] = {}
    for row in fam_df.itertuples(index=False):
        seen.setdefault(row.prior_mode, []).append(row.tier)
    for i, (mode, tier_names) in enumerate(seen.items()):
        dist = _prior_family_distribution(_FAMILY_OF_MODE[mode], 1.0)
        rec = fam_df.loc[fam_df["prior_mode"] == mode].iloc[0]
        mass = ", ".join(f"Pr(<{c.split('_')[2]}s)={rec[c]:.2f}" for c in pr_cols)
        color = _ARM_COLORS[i % len(_ARM_COLORS)]
        axs[0].plot(grid, dist.pdf(grid), color=color, lw=1.8, label=f"{', '.join(tier_names)}\n{mass}")
    axs[0].set_ylim(0.0, min(float(axs[0].get_ylim()[1]), 3.0))
    axs[0].set_xlabel("tau in units of the prior scale s")
    axs[0].set_ylabel("Prior density")
    axs[0].set_title("A. Standardized prior families")
    axs[0].legend(frameon=False, fontsize=7)

    comp = sens_result.posterior_comparison_df
    comp = comp.loc[comp["contrast"] == contrast]
    tiers = [t.name for t in sens_result.tiers]
    horizons = (("observed", _ARM_COLORS[0], 0.15), ("extrapolated", _ARM_COLORS[1], -0.15))
    for kind, color, offset in horizons:
        sub = comp.loc[comp["horizon"] == kind].set_index("tier").loc[tiers]
        y = np.arange(len(tiers), dtype=np.float64) + offset
        periods = int(sub["horizon_periods"].iloc[0])
        axs[1].errorbar(
            sub["P50"].to_numpy(), y, xerr=[sub["P50"] - sub["P30"], sub["P70"] - sub["P50"]], fmt="o", color=color,
            capsize=3, label=f"{kind.capitalize()} horizon ({periods} periods): P50 [P30, P70]",
        )
    axs[1].set_yticks(np.arange(len(tiers)))
    axs[1].set_yticklabels(tiers, fontsize=8)
    axs[1].invert_yaxis()
    axs[1].set_xlabel(f"Cumulative {label} of {contrast}{_unit_suffix(unit)}")
    axs[1].set_title("B. Posterior across prior tiers")
    axs[1].legend(frameon=False, fontsize=8)

    alpha_fine = np.linspace(0.0, 1.0, 201)
    for i, tier_name in enumerate(tiers):
        draws = np.asarray(sens_result.extrapolated_cumulative_by_tier[tier_name][contrast], dtype=np.float64)
        pr = np.array([float(np.mean(a * draws > cost)) for a in alpha_fine])
        color = _ARM_COLORS[i % len(_ARM_COLORS)]
        axs[2].plot(alpha_fine, pr, color=color, lw=1.6, label=tier_name)
        median = float(np.median(draws))
        if median > 0.0 and cost / median <= 1.0:
            axs[2].plot([cost / median], [0.5], marker="v", color=color, ms=8, ls="none")
    axs[2].axhline(0.5, color="grey", lw=0.8, ls="--")
    axs[2].set_xlabel("Incrementality alpha")
    axs[2].set_ylabel(f"Pr(ROI > 0 | alpha, cost = {cost:g}{_unit_suffix(unit)})")
    axs[2].set_ylim(-0.02, 1.02)
    axs[2].set_title(f"C. Break-even of {contrast} over {sens_result.extrap_horizon} periods (markers: alpha_50)")
    axs[2].legend(frameon=False, fontsize=8)
    if title:
        axs[0].figure.suptitle(title)
    return axs


def plot_rank_k_decay_extrapolation(
    decay_result: DecayExtrapolationResult,
    *,
    prior_decay_result: DecayExtrapolationResult | None = None,
    unit: str = "",
    milestones: Sequence[int] = DEFAULT_MILESTONES,
    axes: Any = None,
    title: str | None = None,
) -> np.ndarray:
    """Two-panel decay figure: monthly-rate ribbons by arm and the cumulative fan with the ``rho`` bands.

    Panel A draws the P30-P70 ribbon and P50 line of the monthly rate per arm (solid over
    the observed window, which is shaded; dashed over the extrapolated tail), with the
    matching prior ribbon in grey when ``prior_decay_result`` is given. Panel B draws the
    cumulative P30-P70 fan per arm under the fixed ``rho``, the wider hierarchical
    ``Beta``-``rho`` band when available, thin dotted P50 lines for each decay of the ``rho``
    grid and light vertical markers at the milestones.

    Args:
      decay_result: Output of :func:`extrapolate_rank_k_decay` (posterior draws).
      prior_decay_result: Optional extrapolation of the prior draws with the same arms and horizon.
      unit: Outcome unit for the axis labels.
      milestones: Months marked in panel B (clipped to the horizon).
      axes: Sequence of two Axes, or ``None`` to create a one-row figure.
      title: Figure-level title, or ``None``.

    Returns:
      Length-2 object array of Axes.
    """
    plt = _pyplot("plot_rank_k_decay_extrapolation")
    if prior_decay_result is not None and (
        prior_decay_result.horizon != decay_result.horizon or set(prior_decay_result.arms) != set(decay_result.arms)
    ):
        raise ValueError("prior_decay_result must share the arms and horizon of decay_result")
    months_all = np.arange(1, decay_result.horizon + 1)
    n_obs = decay_result.n_observed
    marks = _resolve_milestones(milestones, decay_result.horizon)
    axs = _axes_array(plt, axes, 2, (13.0, 4.6))

    rate_ax, cum_ax = axs[0], axs[1]
    rho_label = f"rho = {decay_result.rho:g}"
    grid_rhos = sorted(decay_result.rho_grid_cumulative)
    grid_label = f"rho grid {grid_rhos[0]:g}-{grid_rhos[-1]:g}" if grid_rhos else ""

    def _arm_table(frame: pd.DataFrame, arm: str) -> pd.DataFrame:
        return frame.loc[frame["arm"] == arm].sort_values("month")

    for i, arm in enumerate(decay_result.arms):
        color = _ARM_COLORS[i % len(_ARM_COLORS)]
        monthly = _arm_table(decay_result.monthly_summary_df, arm)
        p30, p50, p70 = (monthly[c].to_numpy() for c in ("P30", "P50", "P70"))
        rate_ax.fill_between(months_all, p30, p70, color=color, alpha=0.22)
        rate_ax.plot(months_all[:n_obs], p50[:n_obs], color=color, lw=1.8, label=f"{arm} P50 (observed)")
        rate_ax.plot(
            months_all[n_obs - 1 :], p50[n_obs - 1 :], color=color, lw=1.4, ls="--",
            label=f"{arm} P50 ({rho_label} tail)",
        )
        if prior_decay_result is not None:
            prior = _arm_table(prior_decay_result.monthly_summary_df, arm)
            rate_ax.fill_between(
                months_all, prior["P30"].to_numpy(), prior["P70"].to_numpy(), color="#9e9e9e", alpha=0.18,
                label=f"{arm} prior P30-P70" if i == 0 else None,
            )
        cum = _arm_table(decay_result.cumulative_summary_df, arm)
        cum_ax.fill_between(months_all, cum["P30"].to_numpy(), cum["P70"].to_numpy(), color=color, alpha=0.25)
        cum_ax.plot(months_all, cum["P50"].to_numpy(), color=color, lw=1.8, label=f"{arm} P50 [P30, P70], {rho_label}")
        if decay_result.random_rho_cumulative_draws is not None:
            cum_ax.fill_between(
                months_all, cum["rho_band_P30"].to_numpy(), cum["rho_band_P70"].to_numpy(), color=color, alpha=0.10,
                label=f"{arm} rho ~ Beta band (sd {decay_result.rho_sd:g})",
            )
        for k, rho_k in enumerate(grid_rhos):
            cum_ax.plot(
                months_all, np.median(decay_result.rho_grid_cumulative[rho_k][arm], axis=0), color=color, lw=0.8,
                ls=":", alpha=0.7, label=f"{arm} {grid_label}" if k == 0 else None,
            )
    rate_ax.axvspan(0.5, n_obs + 0.5, color="#f0f0f0", zorder=0, label="Observed window")
    rate_ax.axhline(0.0, color="black", lw=0.8, alpha=0.5)
    rate_ax.set_xlabel("Month since adoption")
    rate_ax.set_ylabel(f"Monthly effect{_unit_suffix(unit)}")
    rate_ax.set_title("A. Monthly rate: observed window and geometric-decay tail")
    rate_ax.legend(frameon=False, fontsize=7)
    for m in marks:
        cum_ax.axvline(m, color="grey", lw=0.6, ls=":", alpha=0.6)
    cum_ax.set_xlabel("Month since adoption")
    cum_ax.set_ylabel(f"Cumulative effect{_unit_suffix(unit)}")
    cum_ax.set_title("B. Cumulative fan with rho bands (milestones dotted)")
    cum_ax.legend(frameon=False, fontsize=7)
    if title:
        axs[0].figure.suptitle(title)
    return axs


__all__ = [
    "BUNDLE_PRIOR_MODES",
    "DEFAULT_MILESTONES",
    "DEFAULT_QUANTILES",
    "DecayExtrapolationResult",
    "NESTED_CONTRASTS",
    "PRIOR_FAMILIES",
    "PriorSensitivityResult",
    "PriorTier",
    "RankKPriorDraws",
    "assert_nested_additivity",
    "evaluate_prior_sensitivity",
    "extrapolate_rank_k_decay",
    "nested_3arm_trajectories",
    "plot_rank_k_decay_extrapolation",
    "plot_rank_k_prior_calibration",
    "plot_rank_k_prior_sensitivity",
    "plot_rank_k_prior_vs_posterior",
    "prior_vs_posterior_summary",
    "roi_break_even_table",
    "sample_rank_k_prior_draws",
    "standardized_prior_family_table",
    "ytlr_prior_tiers",
]
