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

"""Prior diagnostics, prior-sensitivity tiers, decay extrapolation, required persistence and Rank-K figures.

The decision layer of a Rank-K BCF analysis asks the same questions every time:
how much did the data move the prior, how sensitive is the horizon effect to the
prior family and scale, what does the observed post-adoption window imply once
it is extrapolated with an assumed decay, and how long would the effect have to
persist to pay back a cost. This module turns the hand-built notebook versions of
those answers into library calls, split -- as in :mod:`longbet._rollout` -- into a
pure numpy / pandas data layer and thin ``plot_*`` functions that import
matplotlib only when called:

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
  ``tau(L + h) = rho^h tau(L)`` ("LongBet Regime 4"), optionally split into a
  plateau that decays as ``rho^h`` and a pulse that decays as ``exp(-h / tau_p)``
  (the total arm is recomputed as their sum, so ``AC = BC + AB`` holds at every
  period) and weighted by a survival curve (:func:`weibull_survival_weights`),
  with a ``rho`` grid, an assumed ``Beta(rho, sd)`` scenario band that is labelled
  as not identified by the data, milestone tables and a data-health certificate.
* :func:`compute_required_persistence` is the value layer on top of it: it prices
  the trajectories (incrementality share, :func:`discount_weights`) and inverts the
  post-horizon ROI question into the required persistence ``X`` (periods of full
  plateau persistence after ``L`` needed to pay back a cost), the break-even
  retention, payback CDFs under the pre-registered scenarios Bound U, Bound L and
  the ``rho = 0`` floor (plus labelled extra retentions, each reported in a
  headline or diagnostic role), and a placebo decision gate;
  :func:`plot_required_persistence_and_payback` draws it.

Every data-layer function raises ``ValueError`` on malformed or non-finite inputs
instead of propagating NaNs, and every random step takes an explicit seed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import functools
import itertools
import math
from typing import Any
import warnings

import numpy as np
import pandas as pd
from scipy import optimize as sp_optimize
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
#: Default periods after the observed window ``L`` at which ``Pr(X <= k)`` of the required persistence is reported.
DEFAULT_PERSISTENCE_CHECKPOINTS = (12, 24, 48)

# Label of the hierarchical ``rho ~ Beta`` band of :func:`extrapolate_rank_k_decay` in tables and legends.
_RHO_BAND_LABEL = "assumed Beta(rho, sd) scenario band — not identified by the data"
_RHO_SD_FUTURE_WARNING = (
    "extrapolate_rank_k_decay(): rho_sd was not passed, so the 0.07 default adds an assumed Beta(rho, sd) "
    "scenario band, which is not identified by the data. The default becomes None (no band) in longbet 0.9.0; "
    "pass rho_sd explicitly (a float to keep the band, None to drop it) to silence this warning."
)
# Caveat carried by the labels of Bound U and Bound L in :func:`compute_required_persistence`.
_SCENARIO_CAVEAT = "pre-registered scenario, not a statistical bound; faster decay is not excluded"
# Regimes of the value layer: Bound U, Bound L and the rho = 0 floor (plateau retention per period).
_REGIMES = ("U", "L", "rho0")
# Regime label of the decision rows at the end of the observed window, where every regime coincides.
_OBSERVED_REGIME = "observed"
# (scope, role) of the placebo immunity counts by role: "observed" and the headline regimes, the diagnostic
# regimes and the rho = 0 reference floor.
_PLACEBO_ROLE_SCOPES = (("observed+headline", "headline"), ("diagnostic", "diagnostic"), ("reference", "reference"))
# Labels an extra regime may not take: the built-in regimes and the placebo-summary scopes.
_RESERVED_REGIME_LABELS = (*_REGIMES, _OBSERVED_REGIME, *(scope for scope, _ in _PLACEBO_ROLE_SCOPES), "all")
# Relative tolerance of the observed-window identity ``total = plateau + pulse`` under the split.
_ADDITIVITY_RTOL = 1e-8
# Grid of the break-even plateau retention: step 5e-5 on [0, 1].
_BREAKEVEN_GRID_POINTS = 20001
# Elements per block when the plateau value is tabulated on that grid (about 32 MB of float64).
_GRID_BLOCK_ELEMENTS = 4_000_000
# A placebo decision cell is immune when Pr(ROI > 0) under the placebo draws is below this share.
_PLACEBO_DECISION_THRESHOLD = 0.5

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
# (color, line style) of the value-layer regimes in the figures.
_REGIME_STYLES = {"U": ("#188038", "-"), "L": ("#d93025", "--"), "rho0": ("#5f6368", ":")}


class _ImplicitRhoSd(float):
    """The implicit ``rho_sd`` default (0.07) of :func:`extrapolate_rank_k_decay`.

    The default is this distinct float instance, so an explicit ``rho_sd=0.07`` is
    not mistaken for it and only the implicit use emits the ``FutureWarning`` of the
    deprecation path (the default becomes ``None`` in 0.9.0).
    """


_IMPLICIT_RHO_SD = _ImplicitRhoSd(0.07)


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


def _positive_int(value: Any, name: str) -> int:
    """Return ``value`` as an ``int >= 1``, raising ``ValueError`` unless it is a whole positive number."""
    try:
        out = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a positive integer; got {value!r}") from exc
    if out != value or out < 1:
        raise ValueError(f"{name} must be a positive integer; got {value!r}")
    return out


def _unit_interval(value: Any, name: str) -> float:
    """Return ``value`` as a float in ``[0, 1]``, raising ``ValueError`` otherwise."""
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a number in [0, 1]; got {value!r}") from exc
    if not 0.0 <= out <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]; got {value!r}")
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
    """The four standard bundle prior-sensitivity tiers for a historical benchmark scale ``sigma_hist``.

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


def _check_periods(horizon: Any, n_observed: Any) -> tuple[int, int]:
    """Validate integer ``n_observed >= 1`` and ``horizon >= n_observed``; return ``(horizon, n_observed)``."""
    n_obs = _positive_int(n_observed, "n_observed")
    horizon_int = _positive_int(horizon, "horizon")
    if horizon_int < n_obs:
        raise ValueError(f"horizon ({horizon_int}) must be >= the number of observed periods ({n_obs})")
    return horizon_int, n_obs


def weibull_survival_weights(
    horizon: int,
    n_observed: int,
    *,
    median_periods: float = 60.0,
    shape: float = 2.5,
) -> np.ndarray:
    """Survival weights of a Weibull life curve, conditional on surviving the observed window.

    The survival curve is ``S(t) = exp(-(t / lambda)^k)`` with the scale
    ``lambda = median / (ln 2)^(1 / k)``, so ``S(median) = 1 / 2``; the defaults (median 60
    periods, ``k = 2.5``) give ``lambda = 69.47``, ``S(12) = 0.9877`` and ``S(60) = 0.50``. The
    weights are 1 for ``t <= L`` and ``S(t) / S(L)`` for ``t > L``: the share of the units alive
    at the end of the observed window that still carry the effect at ``t``. They are the
    ``survival`` argument of :func:`extrapolate_rank_k_decay` and
    :func:`compute_required_persistence`; the life curve is an assumption, not estimated here.

    Args:
      horizon: Horizon ``H >= L`` in periods.
      n_observed: Observed post-adoption periods ``L >= 1``.
      median_periods: Median life in periods (``> 0``).
      shape: Weibull shape ``k > 0`` (``k > 1``: the hazard increases with age).

    Returns:
      ``(H,)`` float64 weights in ``[0, 1]``, exactly 1 on ``t <= L`` and non-increasing.
    """
    horizon, n_obs = _check_periods(horizon, n_observed)
    median = _positive_scalar(median_periods, "median_periods")
    shape = _positive_scalar(shape, "shape")
    periods = np.arange(1, horizon + 1, dtype=np.float64)
    scale = median / np.log(2.0) ** (1.0 / shape)
    survival = np.exp(-((periods / scale) ** shape))
    if not survival[n_obs - 1] > 0.0:
        raise ValueError(
            f"S(L) underflows to 0 at L = {n_obs} with median_periods={median:g} and shape={shape:g}; "
            "the median life must not be tiny relative to the observed window"
        )
    weights = np.ones(horizon, dtype=np.float64)
    weights[n_obs:] = survival[n_obs:] / survival[n_obs - 1]
    return weights


def discount_weights(horizon: int, *, annual_rate: float = 0.10, periods_per_year: int = 12) -> np.ndarray:
    """Present-value weights ``delta^(t - 1)`` with ``delta = (1 + r)^(-1 / p)`` for ``t = 1..H``.

    The first period is undiscounted and period ``1 + p`` is worth ``1 / (1 + r)``.

    Args:
      horizon: Number of periods ``H >= 1``.
      annual_rate: Annual discount rate ``r > -1`` (default 10%).
      periods_per_year: Periods per year ``p >= 1`` (12 for monthly data).

    Returns:
      ``(H,)`` float64 weights, 1 at ``t = 1``.
    """
    horizon = _positive_int(horizon, "horizon")
    per_year = _positive_int(periods_per_year, "periods_per_year")
    try:
        rate = float(annual_rate)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"annual_rate must be a number > -1; got {annual_rate!r}") from exc
    if not (math.isfinite(rate) and rate > -1.0):
        raise ValueError(f"annual_rate must be a finite number > -1; got {annual_rate!r}")
    period_factor = (1.0 + rate) ** (-1.0 / per_year)
    return period_factor ** np.arange(horizon, dtype=np.float64)


def _extend_with_decay(
    observed: np.ndarray, rho: Any, horizon: int, survival_tail: np.ndarray | None = None
) -> np.ndarray:
    """``(S, L)`` observed draws -> ``(S, H)``: observed months, then ``tau(L + h) = rho^h tau(L) sigma_h``.

    ``rho`` is a scalar or an ``(S,)`` vector (one decay per draw); ``survival_tail`` holds the
    ``(H - L,)`` factors ``sigma_h = S(L + h) / S(L)`` (none when ``None``). ``H == L`` returns a copy.
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
    if survival_tail is not None:
        tail = tail * survival_tail[None, :]
    return np.concatenate([observed, tail], axis=1)


def _extend_with_pulse(
    observed: np.ndarray, pulse_decay_tau: float, horizon: int, survival_tail: np.ndarray | None
) -> np.ndarray:
    """``(S, L)`` observed draws -> ``(S, H)``: observed months, then ``tau(L) exp(-h / tau_p) sigma_h``."""
    n_tail = int(horizon) - observed.shape[1]
    steps = np.arange(1, n_tail + 1, dtype=np.float64)
    tail = observed[:, -1][:, None] * np.exp(-steps / pulse_decay_tau)[None, :]
    if survival_tail is not None:
        tail = tail * survival_tail[None, :]
    return np.concatenate([observed, tail], axis=1)


def _survival_tail(survival: Any, horizon: int, n_obs: int) -> np.ndarray | None:
    """Validate ``(H,)`` survival weights and return the tail factors ``S(L + h) / S(L)`` (``None`` if absent)."""
    if survival is None:
        return None
    try:
        weights = np.asarray(survival, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"survival must be a numeric (H,) array; got {type(survival).__name__}") from exc
    if weights.shape != (horizon,):
        raise ValueError(f"survival must have shape ({horizon},), one weight per period; got {weights.shape}")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0.0) or np.any(weights > 1.0):
        raise ValueError("survival weights must be finite and lie in [0, 1]")
    if np.any(weights[:n_obs] != 1.0):
        raise ValueError(
            f"survival must be exactly 1 over the observed periods t <= {n_obs}: the weights condition on "
            "surviving the observed window (S(t) / S(L)); see weibull_survival_weights"
        )
    return np.array(weights[n_obs:], dtype=np.float64, copy=True)


def _decay_arms(traj_draws: Any, roles: tuple[str, str, str], name: str) -> dict[str, np.ndarray]:
    """``traj_draws`` as ``{arm: (S, L)}``; a reconciliation result maps BC / AB / AC onto ``roles``."""
    if isinstance(traj_draws, Nested3ArmReconciliationResult):
        if len(set(roles)) != 3:
            raise ValueError(f"plateau_arm, pulse_arm and total_arm must be three distinct labels; got {roles}")
        contrasts = (traj_draws.trt_catt_bc, traj_draws.trt_catt_ab, traj_draws.trt_catt_ac)
        return {
            arm: _as_draw_matrix(arr, f"{name}.trt_catt_{key}")
            for arm, arr, key in zip(roles, contrasts, ("bc", "ab", "ac"), strict=True)
        }
    return _as_traj_dict(traj_draws, name)


def _split_roles(
    draws: Mapping[str, np.ndarray], pulse_decay_tau: float | None, roles: tuple[str, str, str]
) -> tuple[str, str, str] | None:
    """Validate the ``(plateau, pulse, total)`` arms of the split; ``None`` when ``pulse_decay_tau`` is ``None``.

    The split needs exactly the three arms, and the observed total must equal plateau + pulse to a
    relative ``_ADDITIVITY_RTOL`` (reconciled draws satisfy it to machine precision).
    """
    if pulse_decay_tau is None:
        return None
    if len(set(roles)) != 3:
        raise ValueError(f"plateau_arm, pulse_arm and total_arm must be three distinct labels; got {roles}")
    missing = [arm for arm in roles if arm not in draws]
    if missing:
        raise ValueError(
            f"the plateau / pulse split (pulse_decay_tau set) needs the arms {list(roles)}; missing {missing}. "
            "Pass {plateau: (S, L), pulse: (S, L), total: (S, L)} or a Nested3ArmReconciliationResult"
        )
    extra = [arm for arm in draws if arm not in roles]
    if extra:
        raise ValueError(f"the plateau / pulse split takes exactly the arms {list(roles)}; got extra arms {extra}")
    plateau, pulse, total = (draws[arm] for arm in roles)
    gap = float(np.max(np.abs(total - (plateau + pulse))))
    tolerance = _ADDITIVITY_RTOL * max(1.0, float(np.max(np.abs(total))))
    if gap > tolerance:
        raise ValueError(
            f"the total arm {roles[2]!r} must equal {roles[0]!r} + {roles[1]!r} over the observed window "
            f"(max |gap| = {gap:.3e} > {tolerance:.1e}); reconcile the draws first (reconcile_nested_3arm_draws)"
        )
    return roles


def _extend_arms(
    draws: Mapping[str, np.ndarray],
    rho: Any,
    horizon: int,
    *,
    roles: tuple[str, str, str] | None,
    pulse_decay_tau: float | None,
    survival_tail: np.ndarray | None,
) -> dict[str, np.ndarray]:
    """``{arm: (S, H)}`` monthly draws of every arm under the plateau decay ``rho`` (scalar or one per draw).

    Without ``roles`` every arm follows ``tau(L) rho^h sigma_h``. With ``roles`` the plateau follows
    ``tau(L) rho^h sigma_h``, the pulse ``tau(L) exp(-h / tau_p) sigma_h``, and the total keeps its
    observed periods while its tail is the plateau tail plus the pulse tail (never extrapolated on
    its own), so the split is additive draw by draw.
    """
    if roles is None:
        return {arm: _extend_with_decay(arr, rho, horizon, survival_tail) for arm, arr in draws.items()}
    plateau_arm, pulse_arm, total_arm = roles
    n_obs = draws[total_arm].shape[1]
    plateau = _extend_with_decay(draws[plateau_arm], rho, horizon, survival_tail)
    pulse = _extend_with_pulse(draws[pulse_arm], pulse_decay_tau, horizon, survival_tail)
    total = np.concatenate([draws[total_arm], plateau[:, n_obs:] + pulse[:, n_obs:]], axis=1)
    extended = {plateau_arm: plateau, pulse_arm: pulse, total_arm: total}
    return {arm: extended[arm] for arm in draws}


def _cumulative_arms(monthly: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Running sums over periods of every ``(S, H)`` array."""
    return {arm: np.cumsum(arr, axis=1) for arm, arr in monthly.items()}


def _assert_split_additivity(per_arm: Mapping[str, np.ndarray], roles: tuple[str, str, str], what: str) -> None:
    """Raise ``AssertionError`` unless ``total = plateau + pulse`` draw by draw (relative ``_ADDITIVITY_RTOL``)."""
    plateau, pulse, total = (per_arm[arm] for arm in roles)
    gap = float(np.max(np.abs(total - plateau - pulse)))
    tolerance = _ADDITIVITY_RTOL * total.shape[1] * max(1.0, float(np.max(np.abs(total))))
    if gap > tolerance:
        raise AssertionError(
            f"{what}: {roles[2]!r} != {roles[0]!r} + {roles[1]!r} draw by draw (max gap {gap:.3e} > {tolerance:.1e})"
        )


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
      rho: Fixed geometric decay of the point extrapolation (of the plateau under the split).
      rho_sd: SD of the assumed ``Beta`` decay scenario, or ``None`` when it was not requested.
      beta_params: ``(a, b)`` of the method-of-moments ``Beta`` decay scenario, or ``None``.
      monthly_draws: ``{arm: (S, H)}`` monthly-rate draws (observed months then the decayed tail).
      cumulative_draws: ``{arm: (S, H)}`` running sums of ``monthly_draws``.
      monthly_summary_df: Long table with the columns ``arm``, ``month`` (1-based), ``segment``
        (``"observed"`` / ``"extrapolated"``), ``P30``, ``P50``, ``P70``, ``Mean`` of the monthly rate.
      cumulative_summary_df: Same layout for the cumulative effect, plus ``rho_band_P30`` /
        ``rho_band_P50`` / ``rho_band_P70`` from the ``rho ~ Beta`` draws (NaN when ``rho_sd`` is
        ``None``) and ``rho_band_label`` (``"assumed Beta(rho, sd) scenario band — not identified
        by the data"``, empty without the band).
      milestone_table_df: One row per arm and milestone month with ``arm``, ``month``, ``segment``,
        ``rate_P30``, ``rate_P50``, ``rate_P70``, ``rate_Mean``, ``cum_P30``, ``cum_P50``, ``cum_P70``,
        ``cum_Mean``.
      data_health_certificate_df: One row per arm and milestone month with ``arm``, ``month``,
        ``n_draws``, ``n_nan``, ``n_inf`` (both asserted to be 0), ``Min``, ``P25``, ``P50``, ``P95``,
        ``Max``, ``Mean`` of the cumulative effect.
      rho_grid_cumulative: ``{rho: {arm: (S, H)}}`` cumulative draws for each decay on ``rho_grid``
        (applied to the plateau only under the split).
      random_rho_cumulative_draws: ``{arm: (S * num_mc_per_draw, H)}`` cumulative draws under the
        assumed ``rho ~ Beta`` scenario (the same decay draws are shared across arms, so additivity
        survives; plateau only under the split), or ``None``.
      scenario_note: Plain-language statement of what is observed and what is an assumed scenario.
      pulse_decay_tau: Time constant of the pulse decay ``exp(-h / tau_p)``, or ``None`` (no split).
      survival: ``(H,)`` survival weights applied to the tail, or ``None``.
      arm_roles: ``{arm: role}``: ``"plateau"``, ``"pulse"`` and ``"total"`` under the split, else
        ``"geometric"`` for every arm.
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
    scenario_note: str = ""
    pulse_decay_tau: float | None = None
    survival: np.ndarray | None = None
    arm_roles: dict[str, str] | None = None


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


def _check_decay_invariants(
    draws: Mapping[str, np.ndarray],
    monthly: Mapping[str, np.ndarray],
    cumulative: Mapping[str, np.ndarray],
    roles: tuple[str, str, str] | None,
) -> None:
    """Assert finiteness, the horizon-``L`` identity, monotone cumulatives for rates >= 0 and split additivity."""
    for arm, observed in draws.items():
        if not (np.all(np.isfinite(monthly[arm])) and np.all(np.isfinite(cumulative[arm]))):
            raise AssertionError(f"the extrapolation of {arm!r} has NaN / Inf values")
        obs_cum = np.cumsum(observed, axis=1)[:, -1]
        if np.max(np.abs(cumulative[arm][:, observed.shape[1] - 1] - obs_cum)) > 1e-10:
            raise AssertionError(f"horizon-L cumulative of {arm!r} differs from the observed cumulative")
        if np.all(monthly[arm] >= 0.0) and np.any(np.diff(cumulative[arm], axis=1) < -1e-12):
            raise AssertionError(f"cumulative effect of {arm!r} decreases although every monthly rate is >= 0")
    if roles is not None:
        _assert_split_additivity(monthly, roles, "monthly draws")


def _decay_milestone_tables(
    monthly: Mapping[str, np.ndarray],
    cumulative: Mapping[str, np.ndarray],
    months: Sequence[int],
    n_obs: int,
    quantiles: Sequence[float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Milestone table and data-health certificate (cumulative effect) of :func:`extrapolate_rank_k_decay`."""
    milestone_rows = []
    health_rows = []
    for arm in monthly:
        for month in months:
            rate = monthly[arm][:, month - 1]
            cum = cumulative[arm][:, month - 1]
            segment = "observed" if month <= n_obs else "extrapolated"
            row: dict[str, Any] = {"arm": arm, "month": month, "segment": segment}
            row.update({f"rate_{k}": v for k, v in _summarize(rate, quantiles).items() if k != "SD"})
            row.update({f"cum_{k}": v for k, v in _summarize(cum, quantiles).items() if k != "SD"})
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
    return pd.DataFrame(milestone_rows), pd.DataFrame(health_rows)


def _decay_scenario_note(
    n_obs: int,
    rho: float,
    rho_sd: float | None,
    roles: tuple[str, str, str] | None,
    pulse_decay_tau: float | None,
    has_survival: bool,
) -> str:
    """Plain-language statement of what :func:`extrapolate_rank_k_decay` observed and what it assumed."""
    if roles is None:
        tail = f"every arm continues as tau(L) * rho^h with rho = {rho:g}"
    else:
        tail = (
            f"the plateau {roles[0]!r} continues as tau(L) * rho^h with rho = {rho:g}, the pulse {roles[1]!r} as "
            f"tau(L) * exp(-h / {pulse_decay_tau:g}) and the total {roles[2]!r} as their sum"
        )
    survival_text = ", times the survival weights S(L + h) / S(L)" if has_survival else ""
    note = (
        f"Periods 1-{n_obs} are posterior draws; after period {n_obs} {tail}{survival_text}. "
        "The decay is an assumed scenario, not estimated from the data"
    )
    if rho_sd is not None:
        note += f"; the rho ~ Beta band (mean {rho:g}, sd {rho_sd:g}) is an {_RHO_BAND_LABEL}"
    return note + "."


def extrapolate_rank_k_decay(
    traj_draws: Any,
    *,
    horizon: int = 60,
    rho: float = 0.85,
    rho_sd: float | None = _IMPLICIT_RHO_SD,
    rho_grid: Sequence[float] = (0.75, 0.80, 0.85, 0.90, 0.95),
    num_mc_per_draw: int = 25,
    seed: int = 2026,
    milestones: Sequence[int] = DEFAULT_MILESTONES,
    quantiles: Sequence[float] = DEFAULT_QUANTILES,
    pulse_decay_tau: float | None = None,
    plateau_arm: str = "BC",
    pulse_arm: str = "AB",
    total_arm: str = "AC",
    survival: np.ndarray | None = None,
) -> DecayExtrapolationResult:
    """Extend observed trajectory draws to a longer horizon with a geometric decay of the last observed rate.

    For every draw ``s`` the observed monthly rates ``tau^(s)(1..L)`` are kept and the tail
    is ``tau^(s)(L + h) = rho^h tau^(s)(L) sigma_h`` for ``h = 1..H - L``, where
    ``sigma_h = S(L + h) / S(L)`` are the survival weights of ``survival`` (1 when it
    is ``None``); the cumulative effect is the running sum. The ``rho_grid`` sensitivity
    repeats the point extrapolation for each decay in the grid.

    With ``pulse_decay_tau`` set, the tail of a nested three-arm design is split into a
    plateau and a pulse: the plateau arm (base vs control) decays as ``rho^h``, the pulse
    arm (incremental vs base) as ``exp(-h / pulse_decay_tau)``, and the total arm keeps its
    observed periods while its tail is recomputed as plateau + pulse -- never extrapolated
    on its own -- so ``AC = BC + AB`` holds at every period (the observed totals must equal
    plateau + pulse to a relative ``1e-8``, as reconciled draws do). The ``rho`` grid and
    the ``Beta`` band then act on the plateau only.

    The hierarchical variant draws ``rho ~ Beta(a, b)`` with mean ``rho`` and SD ``rho_sd``
    (method of moments), ``num_mc_per_draw`` decays per posterior draw, shared across arms
    so ``AC = BC + AB`` survives the extrapolation. It is an assumed scenario band, not
    identified by the data, and is labelled as such in the tables and figures. Calling
    without ``rho_sd`` keeps the 0.07 default but emits a ``FutureWarning``: the default
    becomes ``None`` (no band) in longbet 0.9.0.

    Invariants asserted after the computation: the horizon-``L`` cumulative equals the
    observed one to ``1e-10``, the cumulative effect is non-decreasing wherever the monthly
    rates are non-negative, every output is finite and, under the split, the total equals
    plateau + pulse draw by draw in the point, grid and band extrapolations.

    Args:
      traj_draws: One ``(S, L)`` array (arm label ``"effect"``), ``{arm: (S, L)}`` with a
        common ``S`` and ``L``, or a :class:`Nested3ArmReconciliationResult` (its BC / AB /
        AC draws become the arms ``plateau_arm`` / ``pulse_arm`` / ``total_arm``).
      horizon: Extrapolation horizon ``H >= L`` in periods.
      rho: Fixed decay in ``[0, 1]`` (strictly inside ``(0, 1)`` when ``rho_sd`` is set).
      rho_sd: SD of the assumed ``Beta`` decay scenario, or ``None`` / ``0`` to skip the band.
        Defaults to 0.07 with a ``FutureWarning`` when not passed (``None`` from 0.9.0).
      rho_grid: Decays of the sensitivity grid (each in ``[0, 1]``).
      num_mc_per_draw: Decay draws per posterior draw in the hierarchical variant.
      seed: Seed of the decay draws.
      milestones: Months to tabulate; clipped to ``horizon`` and always including it.
      quantiles: Probabilities of the summary tables (odd length, default P30/P50/P70).
      pulse_decay_tau: Time constant (periods, ``> 0``) of the pulse decay, or ``None`` for no
        split (every arm decays as ``rho^h``).
      plateau_arm: Label of the plateau arm under the split.
      pulse_arm: Label of the pulse arm under the split.
      total_arm: Label of the total arm under the split.
      survival: Optional ``(H,)`` survival weights, finite, in ``[0, 1]`` and exactly 1
        on ``t <= L`` (for example :func:`weibull_survival_weights`); they multiply every tail.

    Returns:
      :class:`DecayExtrapolationResult`.
    """
    roles_requested = (str(plateau_arm), str(pulse_arm), str(total_arm))
    draws = _decay_arms(traj_draws, roles_requested, "traj_draws")
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
    tau_p = None if pulse_decay_tau is None else _positive_scalar(pulse_decay_tau, "pulse_decay_tau")
    roles = _split_roles(draws, tau_p, roles_requested)
    survival_tail = _survival_tail(survival, horizon, n_obs)
    if rho_sd is _IMPLICIT_RHO_SD:
        warnings.warn(_RHO_SD_FUTURE_WARNING, FutureWarning, stacklevel=2)
    sd = None if rho_sd is None or float(rho_sd) == 0.0 else float(rho_sd)
    if sd is not None and sd < 0.0:
        raise ValueError(f"rho_sd must be >= 0 or None; got {rho_sd!r}")

    extend = functools.partial(
        _extend_arms, horizon=horizon, roles=roles, pulse_decay_tau=tau_p, survival_tail=survival_tail
    )
    monthly = extend(draws, rho)
    cumulative = _cumulative_arms(monthly)
    grid_cumulative = {r: _cumulative_arms(extend(draws, r)) for r in grid}
    beta_params = None
    random_cumulative = None
    if sd is not None:
        beta_params = _beta_from_moments(rho, sd)
        rng = np.random.default_rng(int(seed))
        rho_draws = rng.beta(beta_params[0], beta_params[1], size=n_draws * n_mc)
        repeated = {arm: np.repeat(arr, n_mc, axis=0) for arm, arr in draws.items()}
        random_cumulative = _cumulative_arms(extend(repeated, rho_draws))
    _check_decay_invariants(draws, monthly, cumulative, roles)
    if roles is not None:
        for grid_rho, per_arm in grid_cumulative.items():
            _assert_split_additivity(per_arm, roles, f"rho grid {grid_rho:g}")
        if random_cumulative is not None:
            _assert_split_additivity(random_cumulative, roles, "rho ~ Beta band")

    monthly_df = _long_summary(monthly, n_obs, qs)
    cumulative_df = _long_summary(cumulative, n_obs, qs)
    for q in qs:
        label = _quantile_label(q)
        cumulative_df[f"rho_band_{label}"] = float("nan")
        if random_cumulative is not None:
            for arm in draws:
                mask = cumulative_df["arm"] == arm
                cumulative_df.loc[mask, f"rho_band_{label}"] = np.quantile(random_cumulative[arm], q, axis=0)
    cumulative_df["rho_band_label"] = _RHO_BAND_LABEL if random_cumulative is not None else ""
    milestone_df, health_df = _decay_milestone_tables(monthly, cumulative, months, n_obs, qs)
    role_of = {} if roles is None else dict(zip(roles, ("plateau", "pulse", "total"), strict=True))
    return DecayExtrapolationResult(
        arms=tuple(draws), n_observed=int(n_obs), horizon=horizon, rho=rho, rho_sd=sd, beta_params=beta_params,
        monthly_draws=monthly, cumulative_draws=cumulative, monthly_summary_df=monthly_df,
        cumulative_summary_df=cumulative_df, milestone_table_df=milestone_df, data_health_certificate_df=health_df,
        rho_grid_cumulative=grid_cumulative, random_rho_cumulative_draws=random_cumulative,
        scenario_note=_decay_scenario_note(n_obs, rho, sd, roles, tau_p, survival_tail is not None),
        pulse_decay_tau=tau_p,
        survival=None if survival is None else np.array(survival, dtype=np.float64, copy=True),
        arm_roles={arm: role_of.get(arm, "geometric") for arm in draws},
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


def _solve_draw_wise_break_even(
    focal_draws: np.ndarray,
    fixed_offset: np.ndarray,
    cost: float,
    quantile: float = 0.50,
    alpha_upper: float = 100.0,
) -> float:
    """Smallest ``alpha >= 0`` such that ``Q_q(alpha * focal_draws + fixed_offset) >= cost``."""

    def _gap(alpha: float) -> float:
        return float(np.quantile(alpha * focal_draws + fixed_offset, quantile)) - cost

    if _gap(0.0) >= 0.0:
        return 0.0
    fine = np.linspace(0.0, min(2.0, alpha_upper), 201)
    grid = (
        fine
        if alpha_upper <= 2.0
        else np.concatenate([fine, np.geomspace(2.0, alpha_upper, 100)[1:]])
    )
    vals = (
        np.quantile(grid[:, None] * focal_draws[None, :] + fixed_offset[None, :], quantile, axis=1)
        - cost
    )
    hit = np.nonzero(vals >= 0.0)[0]
    if hit.size == 0:
        return float("inf")
    idx = int(hit[0])
    return float(sp_optimize.brentq(_gap, float(grid[idx - 1]), float(grid[idx]), xtol=1e-12))


def roi_break_even_table(
    cumulative_draws: Mapping[str, np.ndarray] | None = None,
    *,
    cost_grid: Sequence[float] = (1.10, 2.85),
    alpha_grid: Sequence[float] = (0.50, 1.00),
    horizon_periods: int = 12,
    tier: str = "",
    component_cum_draws: Mapping[str, np.ndarray] | Mapping[str, Mapping[str, np.ndarray]] | None = None,
    component_weights: Mapping[str, float] | Sequence[Mapping[str, float]] | None = None,
) -> pd.DataFrame:
    """Commercial break-even of a cumulative effect: ``ROI > 0`` iff ``alpha * effect > cost``.

    Supports both single-outcome / scalar-``alpha`` evaluation via ``cumulative_draws`` and
    coupled multi-outcome vector-``alpha`` evaluation via ``component_cum_draws`` and
    ``component_weights`` (evaluating ``sum_m alpha_m * Cum_m^(s) > cost`` draw by draw and
    solving the median break-even ``alpha`` of the focal first component holding the other
    components' weights fixed).

    Args:
      cumulative_draws: ``{contrast: (S,)}`` cumulative-effect draws at the decision horizon.
      cost_grid: Per-unit costs (outcome units, ``> 0``).
      alpha_grid: Incrementality shares in ``[0, 1]`` (the fraction of the effect that is causal revenue).
      horizon_periods: Horizon of the cumulative effect, written into the table.
      tier: Tier name written into the ``tier`` column.
      component_cum_draws: Optional coupled per-outcome cumulative draws ``{outcome: (S,)}``
        (or ``{contrast: {outcome: (S,)}}``).
      component_weights: Optional ``{outcome: weight}`` mapping or sequence of mappings for
        multi-outcome vector-``alpha`` evaluation. When ``None``, the first component varies
        over ``alpha_grid`` and remaining components have weight ``1.0``.

    Returns:
      pandas.DataFrame with the columns ``tier``, ``contrast``, ``horizon_periods``, ``cost``,
      ``alpha``, ``pr_roi_positive`` and ``break_even_alpha_50`` (see
      :class:`PriorSensitivityResult`). ``break_even_alpha_50`` is ``inf`` when no non-negative
      ``alpha`` reaches break-even at the median.
    """
    costs = tuple(_positive_scalar(c, "cost_grid entry") for c in cost_grid)
    alphas = tuple(float(a) for a in alpha_grid)
    if not costs or not alphas:
        raise ValueError("cost_grid and alpha_grid must be non-empty")
    if any(not 0.0 <= a <= 1.0 for a in alphas):
        raise ValueError(f"alpha_grid entries must lie in [0, 1]; got {alphas}")

    if component_cum_draws is None:
        if cumulative_draws is None or not cumulative_draws:
            raise ValueError("pass cumulative_draws or component_cum_draws")
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

    if not isinstance(component_cum_draws, Mapping) or not component_cum_draws:
        raise ValueError("component_cum_draws must be a non-empty mapping")
    first_val = next(iter(component_cum_draws.values()))
    by_contrast: dict[str, dict[str, np.ndarray]] = {}
    if isinstance(first_val, Mapping):
        outer_keys = tuple(str(k) for k in component_cum_draws.keys())
        if set(outer_keys) <= set(NESTED_CONTRASTS):
            for c_key, inner in component_cum_draws.items():
                by_contrast[str(c_key)] = {str(m): np.asarray(v, dtype=np.float64) for m, v in inner.items()}
        else:
            inner_keys = tuple(str(k) for k in first_val.keys())
            for c_key in inner_keys:
                by_contrast[c_key] = {
                    str(m): np.asarray(component_cum_draws[m][c_key], dtype=np.float64)
                    for m in component_cum_draws
                }
    else:
        by_contrast["AC"] = {str(m): np.asarray(v, dtype=np.float64) for m, v in component_cum_draws.items()}

    rows = []
    for contrast, comp_map in by_contrast.items():
        if not comp_map:
            raise ValueError(f"component_cum_draws for contrast {contrast!r} must not be empty")
        comp_names = tuple(comp_map.keys())
        draw_sizes: set[int] = set()
        for m_name, arr in comp_map.items():
            if arr.ndim != 1 or arr.size < 2 or not np.all(np.isfinite(arr)):
                raise ValueError(f"component_cum_draws[{m_name!r}] must be a finite 1-D array with >= 2 draws")
            draw_sizes.add(int(arr.size))
        if len(draw_sizes) != 1:
            raise ValueError(f"all outcomes in component_cum_draws must share the same draw count S; got {sorted(draw_sizes)}")

        focal = comp_names[0]
        if component_weights is None:
            base_non_focal = {m: 1.0 for m in comp_names[1:]}
            weight_specs = [{focal: a, **base_non_focal} for a in alphas]
        elif isinstance(component_weights, Mapping):
            base_non_focal = {m: float(component_weights.get(m, 1.0)) for m in comp_names[1:]}
            if focal in component_weights and float(component_weights[focal]) not in alphas:
                weight_specs = [{focal: float(component_weights[focal]), **base_non_focal}]
            else:
                weight_specs = [{focal: a, **base_non_focal} for a in alphas]
        else:
            weight_specs = [
                {m: float(w_map.get(m, 1.0)) for m in comp_names}
                for w_map in component_weights
            ]
            if not weight_specs:
                raise ValueError("component_weights sequence must not be empty")

        focal_arr = comp_map[focal]
        alpha50_cache: dict[tuple[float, tuple[float, ...]], float] = {}
        for cost in costs:
            for w_spec in weight_specs:
                non_focal_key = (cost, tuple(float(w_spec[m]) for m in comp_names[1:]))
                if non_focal_key not in alpha50_cache:
                    fixed_offset = sum(
                        (w_spec[m] * comp_map[m] for m in comp_names[1:]),
                        start=np.zeros_like(focal_arr),
                    )
                    alpha50_cache[non_focal_key] = _solve_draw_wise_break_even(
                        focal_arr, fixed_offset, cost, quantile=0.50
                    )
                alpha_50 = alpha50_cache[non_focal_key]
                ge_rev = sum(w_spec[m] * comp_map[m] for m in comp_names)
                rows.append(
                    {
                        "tier": tier,
                        "contrast": contrast,
                        "horizon_periods": int(horizon_periods),
                        "cost": cost,
                        "alpha": float(w_spec[focal]),
                        "pr_roi_positive": float(np.mean(ge_rev > cost)),
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


# ---- Required persistence and payback (value / decision layer) -----------------------------------


@dataclass(frozen=True)
class RequiredPersistenceResult:
    """Output of :func:`compute_required_persistence`.

    Probabilities are shares of posterior draws in ``[0, 1]``. The built-in regimes are ``"U"``
    (Bound U: plateau retention ``rho_bounds[1]`` per period, 1 by default, i.e. full persistence),
    ``"L"`` (Bound L: ``rho_bounds[0]``) and ``"rho0"`` (the floor: nothing of the plateau after
    ``L``, only the decaying pulse); ``extra_regimes`` add labelled plateau retentions (a scalar, or
    one per draw paired by index); ``"observed"`` labels the decision rows at ``L``, where every
    regime coincides. Bounds U and L are pre-registered scenarios, not statistical bounds, and
    faster decay is not excluded. Every table that carries regimes has a ``role`` column right after
    ``regime``: ``"headline"`` (``"observed"`` and the ``headline_regimes``), ``"diagnostic"`` (the
    other regimes) or ``"reference"`` (the ``rho0`` floor).

    Attributes:
      n_observed: Observed periods ``L``.
      horizon: Decision horizon ``H``.
      cost_grid: Costs per unit, in outcome units.
      alpha_grid: Incrementality shares ``alpha``.
      rho_bounds: ``(Bound L, Bound U)`` plateau retentions per period.
      pulse_decay_tau: Time constant of the pulse decay, or ``None`` for the whole-level upper
        envelope.
      arm_roles: ``{arm: role}`` of the arms used: ``"plateau"``, ``"pulse"`` and ``"total"`` under
        the split, else one ``"total"`` arm whose whole level persists.
      alpha_share: ``(H,)`` shares ``s_t``; the value multiplier is ``m_t(alpha) = alpha s_t + 1 - s_t``.
      survival: ``(H,)`` survival weights used (ones when ``survival`` was ``None``).
      discount: ``(H,)`` discount weights used (ones when undiscounted).
      checkpoints: Periods ``k`` after ``L`` of the ``pr_X_le_<k>`` columns (NaN with a warning when
        ``k > H - L``).
      report_months: Months of the ``pr_payback_le_<m>`` columns and decision horizons: ``L``,
        ``L + k`` per checkpoint (capped at ``H``) and ``H``.
      quantiles: Probabilities of the quantile columns.
      regime_rhos: ``{"U": rho_U, "L": rho_L, "rho0": 0.0}`` plus each extra regime's retention (a
        float, or the ``(S,)`` per-draw array).
      regime_labels: Display label per regime and for ``"observed"``; Bounds U and L carry
        ``"pre-registered scenario, not a statistical bound; faster decay is not excluded"``.
      regime_roles: ``{regime: role}`` for ``"observed"`` and every regime (see above).
      placebo_regime_rhos: Plateau retention applied to the placebo draws per gated regime (the
        built-in ones plus each extra regime of ``placebo_extra_regimes``), or ``None`` without
        placebo draws.
      monthly_draws: ``{regime: (S, H)}`` monthly effect of the total arm (survival applied,
        undiscounted, before the value multipliers).
      cumulative_value_draws: ``{(regime, alpha): (S, H)}`` cumulative value ``V(m)``.
      required_persistence_draws: ``{(cost, alpha): (S,)}`` required persistence ``X`` in periods
        after ``L``, always along Bound U (``inf``: not paid back by ``H`` even under Bound U).
      breakeven_rho_draws: ``{(cost, alpha): (S,)}`` break-even plateau retention on a grid of
        step ``5e-5`` (0: the observed window and the decaying pulse already pay back by ``H``;
        ``inf``: not even ``rho = 1`` does).
      payback_month_draws: ``{(regime, cost, alpha): (S,)}`` first period with ``V > cost``
        (``H + 1``: never).
      required_persistence_df: One row per ``(cost, alpha)``: ``cost``, ``alpha``,
        ``persisting_level`` (what persists along Bound U), ``observed_value_P50`` (``V(L)``),
        ``persisting_level_value_P50`` (``m_{L+1}(alpha)`` times the persisting level at ``L``),
        ``X_P30`` / ``X_P50`` / ``X_P70`` (inverted-CDF quantiles, ``inf`` when never) with their
        ``*_label`` (``">H-L"`` when never), ``pr_X_le_<k>`` per checkpoint (NaN when ``k > H - L``,
        with a warning), ``pr_X_never`` and ``pr_paid_back_in_window`` (the share with ``X = 0``).
      breakeven_rho_df: One row per ``(cost, alpha)``: ``rho_be_P30`` / ``rho_be_P50`` /
        ``rho_be_P70`` (inverted CDF, ``inf`` when never) with labels, ``pr_rho_be_zero``,
        ``pr_rho_be_never`` and, per extra regime, ``pr_rho_ge_breakeven_<label>`` (the share of
        draws whose retention, paired by index, is at least the break-even retention).
      payback_df: One row per regime, cost and alpha: ``regime``, ``role``, ``rho`` (the median for a
        per-draw retention), ``cost``, ``alpha``, ``payback_P30`` / ``payback_P50`` / ``payback_P70``
        (inverted CDF; ``H + 1`` = never) with labels, ``pr_payback_le_<m>`` per report month and
        ``pr_payback_never``.
      payback_cdf_df: Long table ``regime``, ``role``, ``cost``, ``alpha``, ``month``,
        ``months_after_observed``, ``pr_paid_back`` (``Pr(payback <= month)``) for every month.
      monthly_summary_df: Long table ``regime``, ``role``, ``month``, ``segment``, one column per
        quantile and ``Mean`` of the monthly effect.
      cumulative_value_df: Long table ``regime``, ``role``, ``alpha``, ``month``, ``segment``, one
        column per quantile and ``Mean`` of ``V(month)``.
      milestone_table_df: Scenario table at the milestones: ``regime``, ``role``, ``cost``,
        ``alpha``, ``month``, ``segment``, ``value_P30`` / ``value_P50`` / ``value_P70`` /
        ``value_Mean`` of ``V(month)``, ``pr_roi_positive`` (``Pr(V(month) > cost)``) and
        ``pr_paid_back`` (``Pr(payback <= month)``).
      decision_df: One row per regime, decision horizon (``report_months``), cost and alpha with
        ``regime``, ``role``, ``horizon``, ``cost``, ``alpha`` and ``pr_roi_positive``
        (``Pr(V(horizon) > cost)``); with placebo draws also ``pr_roi_positive_placebo``,
        ``decision_immune`` (nullable boolean: the placebo share is below 0.5) and
        ``immunity_margin`` (0.5 minus the placebo share), missing for an extra regime without
        ``placebo_extra_regimes``.
      placebo_summary_df: Immune-cell counts of the placebo gate, one row per ``scope`` with its
        ``kind`` and ``role``: each regime (``"regime"``, with that regime's role), the roles
        ``"observed+headline"``, ``"diagnostic"`` and ``"reference"`` (``"role"``) and every cell
        (``"all"``, role missing), with ``n_cells`` / ``n_immune`` / ``n_failing`` over the cells with a
        placebo value, ``n_cells_ungated`` (cells of the scope without one: extra regimes without
        ``placebo_extra_regimes``), ``max_pr_roi_positive_placebo`` and ``max_at`` (missing when no cell
        is gated); scopes without any decision row are omitted. ``None`` without placebo draws.
      data_health_certificate_df: One row per ``(cost, alpha)`` and quantity (``required_persistence_X``,
        ``breakeven_rho``, ``payback_month_<regime>``, ``value_at_H_<regime>`` and, for every regime with
        placebo retentions, ``placebo_value_at_H_<regime>``): ``cost``, ``alpha``, ``quantity``,
        ``regime`` and ``role`` (missing for the regime-free ``X`` and break-even rows), ``never_code``
        (the value standing for "never", NaN when not applicable), ``n_draws``, ``n_nan``, ``n_inf``
        (both asserted 0), ``n_never``, ``Min``, ``P25``, ``P50``, ``P95``, ``Max``, ``Mean``.
      cdf_identity_checked: Whether every monthly effect was ``>= 0``, so the identity
        ``Pr(X <= k) = Pr(payback_U <= L + k)`` and the ordering ``U >= L >= rho0`` of the payback CDFs
        (and of every extra regime above ``rho0``) were asserted.
      warnings: Warnings emitted by the call.
      scenario_note: Plain-language statement of what is observed and what is assumed.
    """

    n_observed: int
    horizon: int
    cost_grid: tuple[float, ...]
    alpha_grid: tuple[float, ...]
    rho_bounds: tuple[float, float]
    pulse_decay_tau: float | None
    arm_roles: dict[str, str]
    alpha_share: np.ndarray
    survival: np.ndarray
    discount: np.ndarray
    checkpoints: tuple[int, ...]
    report_months: tuple[int, ...]
    quantiles: tuple[float, ...]
    regime_rhos: dict[str, float | np.ndarray]
    regime_labels: dict[str, str]
    regime_roles: dict[str, str]
    placebo_regime_rhos: dict[str, float | np.ndarray] | None
    monthly_draws: dict[str, np.ndarray]
    cumulative_value_draws: dict[tuple[str, float], np.ndarray]
    required_persistence_draws: dict[tuple[float, float], np.ndarray]
    breakeven_rho_draws: dict[tuple[float, float], np.ndarray]
    payback_month_draws: dict[tuple[str, float, float], np.ndarray]
    required_persistence_df: pd.DataFrame
    breakeven_rho_df: pd.DataFrame
    payback_df: pd.DataFrame
    payback_cdf_df: pd.DataFrame
    monthly_summary_df: pd.DataFrame
    cumulative_value_df: pd.DataFrame
    milestone_table_df: pd.DataFrame
    decision_df: pd.DataFrame
    placebo_summary_df: pd.DataFrame | None
    data_health_certificate_df: pd.DataFrame
    cdf_identity_checked: bool
    warnings: tuple[str, ...]
    scenario_note: str


def _value_arms(
    traj_draws: Any, pulse_decay_tau: float | None, roles: tuple[str, str, str], name: str
) -> tuple[dict[str, np.ndarray], tuple[str, str, str] | None, str]:
    """``(draws, split roles or None, total arm)`` of the value layer.

    Under the split the three role arms are required; otherwise the total arm (or the only arm) is
    used and its whole level persists.
    """
    draws = _decay_arms(traj_draws, roles, name)
    shapes = {arr.shape for arr in draws.values()}
    if len(shapes) != 1:
        raise ValueError(f"every arm of {name} must have the same (S, L) shape; got {sorted(shapes)}")
    split = _split_roles(draws, pulse_decay_tau, roles)
    if split is not None:
        return draws, split, roles[2]
    if len(draws) == 1:
        (only_arm,) = draws
        return draws, None, only_arm
    if roles[2] not in draws:
        raise ValueError(
            f"{name} has the arms {list(draws)}: without the plateau / pulse split (pulse_decay_tau=None) the "
            f"value layer uses the total arm {roles[2]!r} or a single arm"
        )
    return {roles[2]: draws[roles[2]]}, None, roles[2]


def _decision_grid(values: Sequence[float], name: str, *, positive: bool) -> tuple[float, ...]:
    """Non-empty tuple of distinct costs (``> 0``) or incrementality shares (in ``[0, 1]``)."""
    try:
        entries = tuple(values)
    except TypeError as exc:
        raise ValueError(f"{name} must be a sequence of numbers; got {values!r}") from exc
    if positive:
        grid = tuple(_positive_scalar(value, f"{name} entry") for value in entries)
    else:
        grid = tuple(_unit_interval(value, f"{name} entry") for value in entries)
    if not grid:
        raise ValueError(f"{name} must be non-empty")
    if len(set(grid)) != len(grid):
        raise ValueError(f"{name} entries must be distinct; got {grid}")
    return grid


def _alpha_share_vector(alpha_share: Any, horizon: int) -> np.ndarray:
    """``(H,)`` shares ``s_t`` in ``[0, 1]`` from a scalar or an ``(H,)`` array."""
    try:
        shares = np.asarray(alpha_share, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"alpha_share must be a number or an (H,) array; got {type(alpha_share).__name__}") from exc
    if shares.ndim == 0:
        shares = np.full(horizon, float(shares))
    elif shares.shape != (horizon,):
        raise ValueError(f"alpha_share must be a scalar or have shape ({horizon},); got {shares.shape}")
    if not np.all(np.isfinite(shares)) or np.any(shares < 0.0) or np.any(shares > 1.0):
        raise ValueError("alpha_share entries must be finite and lie in [0, 1]")
    return np.array(shares, dtype=np.float64, copy=True)


def _discount_vector(discount: Any, discount_rate: float | None, horizon: int) -> np.ndarray:
    """``(H,)`` positive discount weights from ``discount`` or ``discount_rate`` (ones when neither is given)."""
    if discount is not None and discount_rate is not None:
        raise ValueError("pass discount (an (H,) array of weights) or discount_rate, not both")
    if discount_rate is not None:
        return discount_weights(horizon, annual_rate=discount_rate)
    if discount is None:
        return np.ones(horizon, dtype=np.float64)
    try:
        weights = np.asarray(discount, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"discount must be a numeric (H,) array; got {type(discount).__name__}") from exc
    if weights.shape != (horizon,):
        raise ValueError(f"discount must have shape ({horizon},), one weight per period; got {weights.shape}")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("discount weights must be finite and positive")
    return np.array(weights, dtype=np.float64, copy=True)


def _check_rho_bounds(rho_bounds: Sequence[float]) -> tuple[float, float]:
    """Validate ``(lower, upper)`` plateau retentions with ``0 <= lower <= upper <= 1``."""
    try:
        lower, upper = (float(value) for value in rho_bounds)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"rho_bounds must be a (lower, upper) pair of numbers; got {rho_bounds!r}") from exc
    if not 0.0 <= lower <= upper <= 1.0:
        raise ValueError(f"rho_bounds must satisfy 0 <= lower <= upper <= 1; got {rho_bounds!r}")
    return lower, upper


def _check_checkpoints(checkpoints: Sequence[int]) -> tuple[int, ...]:
    """Sorted distinct positive integers ``k`` (periods after ``L``) of the ``Pr(X <= k)`` columns."""
    try:
        entries = tuple(checkpoints)
    except TypeError as exc:
        raise ValueError(f"checkpoints must be a sequence of positive integers; got {checkpoints!r}") from exc
    ks = tuple(sorted({_positive_int(k, "checkpoints entry") for k in entries}))
    if not ks:
        raise ValueError("checkpoints must hold at least one positive integer")
    return ks


def _extra_regime_rhos(extra_regimes: Any, n_draws: int, name: str) -> dict[str, float | np.ndarray]:
    """Validated ``{label: rho}`` of extra regimes: a float, or an ``(S,)`` array paired with the draws by index."""
    if extra_regimes is None:
        return {}
    if not isinstance(extra_regimes, Mapping):
        raise ValueError(f"{name} must be a mapping {{label: rho}}; got {type(extra_regimes).__name__}")
    rhos: dict[str, float | np.ndarray] = {}
    for label, value in extra_regimes.items():
        if not isinstance(label, str) or not label:
            raise ValueError(f"{name} labels must be non-empty strings; got {label!r}")
        if label in _RESERVED_REGIME_LABELS:
            raise ValueError(
                f"{name} label {label!r} is reserved: {list(_RESERVED_REGIME_LABELS)} name the built-in regimes "
                "and the placebo-summary scopes"
            )
        try:
            rho = np.asarray(value, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{name}[{label!r}] must be a number or an ({n_draws},) array; got {type(value).__name__}"
            ) from exc
        if rho.ndim != 0 and rho.shape != (n_draws,):
            raise ValueError(
                f"{name}[{label!r}] must be a scalar or have shape ({n_draws},), one retention per draw paired by "
                f"index; got {rho.shape}"
            )
        if not np.all(np.isfinite(rho)) or np.any(rho < 0.0) or np.any(rho > 1.0):
            raise ValueError(f"{name}[{label!r}] must be finite and lie in [0, 1]")
        rhos[label] = float(rho) if rho.ndim == 0 else np.array(rho, dtype=np.float64, copy=True)
    return rhos


def _resolve_headline_regimes(headline_regimes: Any, extra_labels: Sequence[str]) -> tuple[str, ...]:
    """Distinct headline regimes among ``"U"``, ``"L"`` and the extra labels (one label may be passed as a string)."""
    entries = (headline_regimes,) if isinstance(headline_regimes, str) else headline_regimes
    try:
        labels = tuple(entries)
    except TypeError as exc:
        raise ValueError(f"headline_regimes must be a sequence of regime labels; got {headline_regimes!r}") from exc
    allowed = ("U", "L", *extra_labels)
    for label in labels:
        if label not in allowed:
            raise ValueError(
                f"headline_regimes entries must be among {list(allowed)}; got {label!r} ('observed' is always a "
                "headline row and 'rho0' is the reference floor)"
            )
    if len(set(labels)) != len(labels):
        raise ValueError(f"headline_regimes entries must be distinct; got {labels}")
    return labels


def _regime_roles(extra_labels: Sequence[str], headline: Sequence[str]) -> dict[str, str]:
    """``{regime: role}``: ``"observed"`` and the headline regimes are ``"headline"``, ``rho0`` the ``"reference"``."""
    roles = {_OBSERVED_REGIME: "headline"}
    for regime in (*_REGIMES, *extra_labels):
        if regime == "rho0":
            roles[regime] = "reference"
        else:
            roles[regime] = "headline" if regime in headline else "diagnostic"
    return roles


def _rho_display(rho: float | np.ndarray) -> float:
    """A regime's retention as one number: the value itself, or the median of a per-draw retention."""
    return float(np.median(rho)) if isinstance(rho, np.ndarray) else float(rho)


def _cumulative_values(
    monthly: Mapping[str, np.ndarray], value_weights: Mapping[float, np.ndarray]
) -> dict[tuple[str, float], np.ndarray]:
    """``{(regime, alpha): (S, H)}`` running sums of the priced monthly effects ``tau(t) m_t(alpha) delta_t``."""
    out = {}
    for regime, arr in monthly.items():
        for alpha, weights in value_weights.items():
            out[(regime, alpha)] = np.cumsum(arr * weights[None, :], axis=1)
    return out


def _payback_months(cumulative: np.ndarray, cost: float) -> np.ndarray:
    """``(S,)`` first period with cumulative value ``> cost`` (strict), ``H + 1`` when never."""
    paid = cumulative > cost
    first = np.argmax(paid, axis=1) + 1
    return np.where(paid.any(axis=1), first, cumulative.shape[1] + 1)


def _required_persistence(value: np.ndarray, cumulative: np.ndarray, cost: float, n_obs: int) -> np.ndarray:
    """``(S,)`` periods after ``L`` until the Bound-U value exceeds ``cost``, interpolated within the period.

    With ``m*`` the first period with ``V(m*) > cost``, ``X = max(0, m* - 1 - L + f)`` and
    ``f = (cost - V(m* - 1)) / v(m*)`` in ``[0, 1)``; ``inf`` when ``V`` never exceeds the cost. When
    ``m* > L`` and ``f`` is 0 (an exact tie: the value equals the cost at the end of period ``m* - 1``)
    or below half an ulp of ``m* - 1 - L``, the sum would land on that integer although the cost is
    only exceeded one period later, so ``X`` is moved to the next float above it: ``Pr(X <= k) =
    Pr(payback <= L + k)`` keeps the strict payback rule.
    """
    paid = cumulative > cost
    first = np.argmax(paid, axis=1)
    rows = np.arange(cumulative.shape[0])
    before = np.where(first > 0, cumulative[rows, np.maximum(first - 1, 0)], 0.0)
    out = np.full(cumulative.shape[0], np.inf)
    hit = np.flatnonzero(paid.any(axis=1))
    fraction = (cost - before[hit]) / value[hit, first[hit]]
    whole = (first[hit] - n_obs).astype(np.float64)
    months = np.where(whole >= 0.0, whole + fraction, 0.0)
    lost = (whole >= 0.0) & (months <= whole)
    months[lost] = np.nextafter(whole[lost], np.inf)
    out[hit] = months
    return out


def _plateau_value_on_grid(grid: np.ndarray, tail_weights: np.ndarray) -> np.ndarray:
    """``sum_h rho^h w_h`` for every ``rho`` on ``grid``, tabulated in row blocks to bound the memory."""
    steps = np.arange(1, tail_weights.size + 1, dtype=np.float64)
    block = max(1, _GRID_BLOCK_ELEMENTS // max(1, tail_weights.size))
    return np.concatenate(
        [(grid[start : start + block, None] ** steps[None, :]) @ tail_weights for start in range(0, grid.size, block)]
    )


def _breakeven_parts(
    observed_total: np.ndarray,
    pulse_level: np.ndarray | None,
    weights: np.ndarray,
    survival_tail: np.ndarray | None,
    pulse_decay_tau: float | None,
    grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cost-free pieces of ``V_rho(H)`` for one ``alpha``: ``(observed value, pulse value, plateau value)``.

    ``V_rho(H) = V_obs + pulse value + tau_plateau(L) sum_h rho^h w_h`` with ``w_h`` the priced survival
    weights of period ``L + h``; the last term is returned per unit plateau level on ``grid``.
    """
    n_obs = observed_total.shape[1]
    tail_weights = weights[n_obs:] if survival_tail is None else weights[n_obs:] * survival_tail
    observed_value = (observed_total * weights[None, :n_obs]).sum(axis=1)
    pulse_value = np.zeros(observed_total.shape[0])
    if pulse_level is not None:
        steps = np.arange(1, tail_weights.size + 1, dtype=np.float64)
        pulse_value = pulse_level * float(np.sum(tail_weights * np.exp(-steps / pulse_decay_tau)))
    return observed_value, pulse_value, _plateau_value_on_grid(grid, tail_weights)


def _breakeven_rho(
    cost: float,
    observed_value: np.ndarray,
    pulse_value: np.ndarray,
    plateau_level: np.ndarray,
    plateau_value: np.ndarray,
    grid: np.ndarray,
) -> np.ndarray:
    """``(S,)`` smallest grid ``rho`` with ``V_rho(H) > cost`` (0: ``rho = 0`` suffices; ``inf``: ``rho = 1`` fails).

    ``V_rho(H)`` is affine in the plateau level and non-decreasing in ``rho``, so each draw is solved by
    ``searchsorted`` on the tabulated plateau value. The rule is strict, like the payback month: at an
    exact tie the cost is not yet recovered.
    """
    residual = (cost - observed_value) - pulse_value
    rho_be = np.full(residual.shape[0], np.inf)
    rho_be[residual < 0.0] = 0.0
    feasible = (residual >= 0.0) & (plateau_level * plateau_value[-1] > residual)
    ratio = residual[feasible] / plateau_level[feasible]
    index = np.minimum(np.searchsorted(plateau_value, ratio, side="right"), grid.size - 1)
    rho_be[feasible] = grid[index]
    return rho_be


def _payback_cdf(months_drawn: np.ndarray, horizon: int) -> np.ndarray:
    """``(H,)`` ``Pr(payback <= m)`` for ``m = 1..H``."""
    return np.mean(months_drawn[:, None] <= np.arange(1, horizon + 1)[None, :], axis=0)


def _all_nonnegative(extended: Mapping[str, Mapping[str, np.ndarray]]) -> bool:
    """Whether every monthly draw of every arm in every regime is ``>= 0``."""
    for arms in extended.values():
        for arr in arms.values():
            if np.any(arr < 0.0):
                return False
    return True


def _check_value_invariants(
    draws: Mapping[str, np.ndarray],
    total_key: str,
    monthly: Mapping[str, np.ndarray],
    cumulative: Mapping[tuple[str, float], np.ndarray],
    required: Mapping[tuple[float, float], np.ndarray],
    breakeven: Mapping[tuple[float, float], np.ndarray],
) -> None:
    """Hard assertions (i) and (v): one observed cumulative across regimes, and no NaN / Inf."""
    n_obs = draws[total_key].shape[1]
    observed_sum = draws[total_key].sum(axis=1)
    tolerance = 1e-10 * max(1.0, float(np.max(np.abs(observed_sum))))
    for regime, arr in monthly.items():
        if not np.all(np.isfinite(arr)):
            raise AssertionError(f"the monthly effect of regime {regime!r} has NaN / Inf values")
        gap = float(np.max(np.abs(np.cumsum(arr, axis=1)[:, n_obs - 1] - observed_sum)))
        if gap > tolerance:
            raise AssertionError(
                f"cumulative effect at L of regime {regime!r} differs from the observed sum by {gap:.3e}"
            )
    for (regime, alpha), arr in cumulative.items():
        if not np.all(np.isfinite(arr)):
            raise AssertionError(f"the cumulative value of ({regime!r}, alpha {alpha:g}) has NaN / Inf values")
        if np.max(np.abs(arr[:, n_obs - 1] - cumulative[("U", alpha)][:, n_obs - 1])) > tolerance:
            raise AssertionError(f"the cumulative value at L of regime {regime!r} differs from Bound U")
    for key, values in required.items():
        if np.any(np.isnan(values)) or np.any(values < 0.0):
            raise AssertionError(f"required persistence at (cost, alpha) = {key} has NaN or negative draws")
    for key, values in breakeven.items():
        finite = values[np.isfinite(values)]
        if np.any(np.isnan(values)) or np.any(finite < 0.0) or np.any(finite > 1.0):
            raise AssertionError(f"break-even rho at (cost, alpha) = {key} is NaN or outside [0, 1]")


def _check_payback_identities(
    required: Mapping[tuple[float, float], np.ndarray],
    payback: Mapping[tuple[str, float, float], np.ndarray],
    regimes: Sequence[str],
    n_obs: int,
    horizon: int,
) -> None:
    """Hard assertions (iii) and (iv), valid when every monthly effect is ``>= 0``."""
    steps = np.arange(horizon - n_obs + 1)
    for (cost, alpha), persistence in required.items():
        lhs = np.mean(persistence[:, None] <= steps[None, :], axis=0)
        rhs = np.mean(payback[("U", cost, alpha)][:, None] <= (n_obs + steps)[None, :], axis=0)
        gap = float(np.max(np.abs(lhs - rhs)))
        if gap > 1e-12:
            raise AssertionError(
                f"Pr(X <= k) differs from Pr(payback_U <= L + k) by {gap:.3e} at (cost, alpha) = ({cost:g}, {alpha:g})"
            )
        cdf = {regime: _payback_cdf(payback[(regime, cost, alpha)], horizon) for regime in regimes}
        if np.any(cdf["U"] < cdf["L"] - 1e-12) or np.any(cdf["L"] < cdf["rho0"] - 1e-12):
            raise AssertionError(
                f"payback CDFs are not ordered U >= L >= rho0 at (cost, alpha) = ({cost:g}, {alpha:g})"
            )
        for regime in regimes:
            if regime not in _REGIMES and np.any(cdf[regime] < cdf["rho0"] - 1e-12):
                raise AssertionError(
                    f"the payback CDF of regime {regime!r} is below the rho0 floor at (cost, alpha) = "
                    f"({cost:g}, {alpha:g})"
                )


def _discrete_quantile(values: np.ndarray, q: float) -> float:
    """Inverted-CDF quantile: always an observed draw, so ``inf`` (never) draws are never interpolated."""
    return float(np.quantile(values, q, method="inverted_cdf"))


def _periods_label(value: float, cap: int) -> str:
    """One-decimal label of a period count; values beyond ``cap`` (never) show as ``'>cap'``."""
    return f">{cap}" if not math.isfinite(value) or value > cap else f"{value:.1f}"


def _breakeven_label(value: float, horizon: int, split: bool) -> str:
    """Label of a break-even retention."""
    if value <= 0.0:
        return "0 (observed window + pulse suffice)" if split else "0 (observed window suffices)"
    if not math.isfinite(value) or value > 1.0:
        return f">1 (never by period {horizon})"
    return f"{value:.3f}"


def _required_persistence_table(
    required: Mapping[tuple[float, float], np.ndarray],
    cumulative: Mapping[tuple[str, float], np.ndarray],
    shares: np.ndarray,
    plateau_level: np.ndarray,
    persisting_level: str,
    checkpoints: Sequence[int],
    quantiles: Sequence[float],
    n_obs: int,
) -> pd.DataFrame:
    """One row per ``(cost, alpha)`` of the required persistence ``X``."""
    horizon = shares.size
    rows = []
    for cost, alpha in required:
        persistence = required[(cost, alpha)]
        row: dict[str, Any] = {"cost": cost, "alpha": alpha, "persisting_level": persisting_level}
        row["observed_value_P50"] = float(np.quantile(cumulative[("U", alpha)][:, n_obs - 1], 0.5))
        row["persisting_level_value_P50"] = float("nan")
        if horizon > n_obs:
            multiplier = alpha * shares[n_obs] + (1.0 - shares[n_obs])
            row["persisting_level_value_P50"] = float(np.quantile(multiplier * plateau_level, 0.5))
        for q in quantiles:
            label = _quantile_label(q)
            value = _discrete_quantile(persistence, q)
            row[f"X_{label}"] = value
            row[f"X_{label}_label"] = _periods_label(value, horizon - n_obs)
        max_k = horizon - n_obs
        for k in checkpoints:
            row[f"pr_X_le_{k}"] = float("nan") if k > max_k else float(np.mean(persistence <= k))
        row["pr_X_never"] = float(np.mean(np.isinf(persistence)))
        row["pr_paid_back_in_window"] = float(np.mean(persistence == 0.0))
        rows.append(row)
    return pd.DataFrame(rows)


def _breakeven_table(
    breakeven: Mapping[tuple[float, float], np.ndarray],
    extra_rhos: Mapping[str, float | np.ndarray],
    quantiles: Sequence[float],
    horizon: int,
    split: bool,
) -> pd.DataFrame:
    """One row per ``(cost, alpha)`` of the break-even plateau retention."""
    rows = []
    for (cost, alpha), rho_be in breakeven.items():
        row: dict[str, Any] = {"cost": cost, "alpha": alpha}
        for q in quantiles:
            label = _quantile_label(q)
            value = _discrete_quantile(rho_be, q)
            row[f"rho_be_{label}"] = value
            row[f"rho_be_{label}_label"] = _breakeven_label(value, horizon, split)
        row["pr_rho_be_zero"] = float(np.mean(rho_be <= 0.0))
        row["pr_rho_be_never"] = float(np.mean(np.isinf(rho_be)))
        for regime, rho in extra_rhos.items():  # a per-draw retention is paired with rho_be by index
            row[f"pr_rho_ge_breakeven_{regime}"] = float(np.mean(rho >= rho_be))
        rows.append(row)
    return pd.DataFrame(rows)


def _payback_tables(
    payback: Mapping[tuple[str, float, float], np.ndarray],
    regime_rhos: Mapping[str, float | np.ndarray],
    roles: Mapping[str, str],
    report_months: Sequence[int],
    quantiles: Sequence[float],
    n_obs: int,
    horizon: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Payback summary (one row per regime, cost and alpha) and the long payback CDF over every month."""
    rows = []
    frames = []
    months = np.arange(1, horizon + 1)
    for (regime, cost, alpha), months_drawn in payback.items():
        row: dict[str, Any] = {
            "regime": regime, "role": roles[regime], "rho": _rho_display(regime_rhos[regime]), "cost": cost,
            "alpha": alpha,
        }
        for q in quantiles:
            label = _quantile_label(q)
            value = int(np.quantile(months_drawn, q, method="inverted_cdf"))
            row[f"payback_{label}"] = value
            row[f"payback_{label}_label"] = f">{horizon}" if value > horizon else str(value)
        for month in report_months:
            row[f"pr_payback_le_{month}"] = float(np.mean(months_drawn <= month))
        row["pr_payback_never"] = float(np.mean(months_drawn > horizon))
        rows.append(row)
        frames.append(
            pd.DataFrame(
                {
                    "regime": regime, "role": roles[regime], "cost": cost, "alpha": alpha, "month": months,
                    "months_after_observed": months - n_obs, "pr_paid_back": _payback_cdf(months_drawn, horizon),
                }
            )
        )
    return pd.DataFrame(rows), pd.concat(frames, ignore_index=True)


def _per_period_table(draws: np.ndarray, quantiles: Sequence[float], n_obs: int, **keys: Any) -> pd.DataFrame:
    """Per-period quantiles and mean of ``(S, H)`` draws as a long table led by the ``keys`` columns."""
    months = np.arange(1, draws.shape[1] + 1)
    table: dict[str, Any] = dict(keys)
    table["month"] = months
    table["segment"] = np.where(months <= n_obs, "observed", "extrapolated")
    for q in quantiles:
        table[_quantile_label(q)] = np.quantile(draws, q, axis=0)
    table["Mean"] = draws.mean(axis=0)
    return pd.DataFrame(table)


def _value_milestone_table(
    cumulative: Mapping[tuple[str, float], np.ndarray],
    payback: Mapping[tuple[str, float, float], np.ndarray],
    roles: Mapping[str, str],
    milestones: Sequence[int],
    quantiles: Sequence[float],
    n_obs: int,
) -> pd.DataFrame:
    """Scenario table: cumulative value quantiles, ``Pr(ROI > 0)`` and the payback CDF at the milestones."""
    rows = []
    for (regime, cost, alpha), months_drawn in payback.items():
        values_by_month = cumulative[(regime, alpha)]
        for month in milestones:
            values = values_by_month[:, month - 1]
            row: dict[str, Any] = {
                "regime": regime, "role": roles[regime], "cost": cost, "alpha": alpha, "month": month,
            }
            row["segment"] = "observed" if month <= n_obs else "extrapolated"
            row.update({f"value_{k}": v for k, v in _summarize(values, quantiles).items() if k != "SD"})
            row["pr_roi_positive"] = float(np.mean(values > cost))
            row["pr_paid_back"] = float(np.mean(months_drawn <= month))
            rows.append(row)
    return pd.DataFrame(rows)


def _decision_table(
    cumulative: Mapping[tuple[str, float], np.ndarray],
    placebo_cumulative: Mapping[tuple[str, float], np.ndarray] | None,
    regimes: Sequence[str],
    roles: Mapping[str, str],
    costs: Sequence[float],
    alphas: Sequence[float],
    horizons: Sequence[int],
    n_obs: int,
) -> pd.DataFrame:
    """``Pr(V(horizon) > cost)`` per regime, horizon, cost and alpha, with the placebo gate when available.

    A regime absent from ``placebo_cumulative`` (an extra regime without placebo retentions) gets
    missing placebo columns.
    """
    rows = []
    for month in horizons:
        month_regimes = (_OBSERVED_REGIME,) if month == n_obs else regimes
        for regime, cost, alpha in itertools.product(month_regimes, costs, alphas):
            source = "U" if regime == _OBSERVED_REGIME else regime
            row: dict[str, Any] = {
                "regime": regime, "role": roles[regime], "horizon": month, "cost": cost, "alpha": alpha,
            }
            row["pr_roi_positive"] = float(np.mean(cumulative[(source, alpha)][:, month - 1] > cost))
            if placebo_cumulative is not None:
                placebo_values = placebo_cumulative.get((source, alpha))
                placebo_pr = math.nan
                if placebo_values is not None:
                    placebo_pr = float(np.mean(placebo_values[:, month - 1] > cost))
                row["pr_roi_positive_placebo"] = placebo_pr
                row["decision_immune"] = pd.NA if placebo_values is None else placebo_pr < _PLACEBO_DECISION_THRESHOLD
                row["immunity_margin"] = _PLACEBO_DECISION_THRESHOLD - placebo_pr
            rows.append(row)
    decision = pd.DataFrame(rows)
    if placebo_cumulative is not None:
        decision["decision_immune"] = decision["decision_immune"].astype("boolean")
    return decision


def _placebo_summary(decision: pd.DataFrame) -> pd.DataFrame:
    """Immune-cell counts of the placebo gate per regime, per role and overall, with the ungated cells apart."""
    gated = decision["pr_roi_positive_placebo"].notna()
    scopes = []
    for regime in dict.fromkeys(decision["regime"]):
        members = decision["regime"] == regime
        scopes.append((regime, "regime", decision.loc[members, "role"].iloc[0], members))
    scopes.extend((scope, "role", role, decision["role"] == role) for scope, role in _PLACEBO_ROLE_SCOPES)
    scopes.append(("all", "all", None, pd.Series(True, index=decision.index)))
    rows = []
    for scope, kind, role, members in scopes:
        if not members.any():
            continue
        cells = decision.loc[members & gated]
        n_immune = int(cells["decision_immune"].sum())
        row: dict[str, Any] = {
            "scope": scope, "kind": kind, "role": role, "n_cells": int(len(cells)), "n_immune": n_immune,
            "n_failing": int(len(cells)) - n_immune, "n_cells_ungated": int((members & ~gated).sum()),
            "max_pr_roi_positive_placebo": math.nan, "max_at": None,
        }
        if not cells.empty:
            worst = cells.loc[cells["pr_roi_positive_placebo"].idxmax()]
            row["max_pr_roi_positive_placebo"] = float(worst["pr_roi_positive_placebo"])
            row["max_at"] = (
                f"{worst['regime']}, horizon {int(worst['horizon'])}, cost {float(worst['cost']):g}, "
                f"alpha {float(worst['alpha']):g}"
            )
        rows.append(row)
    return pd.DataFrame(rows)


def _certificate_row(
    cost: float,
    alpha: float,
    quantity: str,
    coded: np.ndarray,
    n_never: int,
    never_code: float,
    regime: str | None = None,
    role: str | None = None,
) -> dict[str, Any]:
    """One data-health row: counts, Null / NaN / Inf audit and the five-number summary plus the mean."""
    coded = np.asarray(coded, dtype=np.float64)
    return {
        "cost": cost, "alpha": alpha, "quantity": quantity, "regime": regime, "role": role, "never_code": never_code,
        "n_draws": int(coded.size), "n_nan": int(np.sum(np.isnan(coded))), "n_inf": int(np.sum(np.isinf(coded))),
        "n_never": int(n_never), "Min": float(np.min(coded)), "P25": float(np.quantile(coded, 0.25)),
        "P50": float(np.quantile(coded, 0.50)), "P95": float(np.quantile(coded, 0.95)), "Max": float(np.max(coded)),
        "Mean": float(np.mean(coded)),
    }


def _value_certificate(
    required: Mapping[tuple[float, float], np.ndarray],
    breakeven: Mapping[tuple[float, float], np.ndarray],
    payback: Mapping[tuple[str, float, float], np.ndarray],
    cumulative: Mapping[tuple[str, float], np.ndarray],
    placebo_cumulative: Mapping[tuple[str, float], np.ndarray] | None,
    roles: Mapping[str, str],
    regimes: Sequence[str],
    n_obs: int,
    horizon: int,
) -> pd.DataFrame:
    """Data-health certificate per ``(cost, alpha)``; "never" is coded explicitly and counted in ``n_never``."""
    rows = []
    never_x = float(horizon - n_obs + 1)
    for cost, alpha in required:
        persistence = required[(cost, alpha)]
        rows.append(
            _certificate_row(
                cost, alpha, "required_persistence_X", np.where(np.isinf(persistence), never_x, persistence),
                int(np.sum(np.isinf(persistence))), never_x,
            )
        )
        rho_be = breakeven[(cost, alpha)]
        rows.append(
            _certificate_row(
                cost, alpha, "breakeven_rho", np.where(np.isinf(rho_be), 2.0, rho_be),
                int(np.sum(np.isinf(rho_be))), 2.0,
            )
        )
        for regime in regimes:
            months_drawn = payback[(regime, cost, alpha)]
            rows.append(
                _certificate_row(
                    cost, alpha, f"payback_month_{regime}", months_drawn, int(np.sum(months_drawn > horizon)),
                    float(horizon + 1), regime, roles[regime],
                )
            )
        for regime in regimes:
            at_horizon = cumulative[(regime, alpha)][:, -1]
            rows.append(
                _certificate_row(
                    cost, alpha, f"value_at_H_{regime}", at_horizon, 0, math.nan, regime, roles[regime]
                )
            )
            if placebo_cumulative is not None and (regime, alpha) in placebo_cumulative:
                at_horizon = placebo_cumulative[(regime, alpha)][:, -1]
                rows.append(
                    _certificate_row(
                        cost, alpha, f"placebo_value_at_H_{regime}", at_horizon, 0, math.nan, regime, roles[regime]
                    )
                )
    certificate = pd.DataFrame(rows)
    if (certificate["n_nan"] > 0).any() or (certificate["n_inf"] > 0).any():
        raise AssertionError("the data-health certificate of compute_required_persistence has NaN / Inf entries")
    return certificate


def _regime_labels(
    rho_lower: float,
    rho_upper: float,
    split: bool,
    extra_rhos: Mapping[str, float | np.ndarray],
    roles: Mapping[str, str],
) -> dict[str, str]:
    """Display labels; Bounds U and L carry the pre-registered-scenario caveat, extra regimes their role."""
    level = "plateau" if split else "whole level"
    floor = "observed window and decaying pulse only" if split else "observed window only"
    labels = {
        _OBSERVED_REGIME: "Observed window (identical across regimes)",
        "U": f"Bound U: {level} rho = {rho_upper:g} ({_SCENARIO_CAVEAT})",
        "L": f"Bound L: {level} rho = {rho_lower:g} ({_SCENARIO_CAVEAT})",
        "rho0": f"Floor: {level} rho = 0 ({floor}; a reference, not a bound)",
    }
    for regime, rho in extra_rhos.items():
        value = f"per draw (median {_rho_display(rho):.4g})" if isinstance(rho, np.ndarray) else f"= {rho:g}"
        labels[regime] = f"{regime}: {level} rho {value} ({roles[regime]} regime)"
    return labels


def _value_scenario_note(
    n_obs: int,
    horizon: int,
    level_name: str,
    pulse_text: str,
    rho_bounds: tuple[float, float],
    extra_labels: Sequence[str],
    has_survival: bool,
    discounted: bool,
) -> str:
    """Plain-language statement of what :func:`compute_required_persistence` observed and what it assumed."""
    rho_lower, rho_upper = rho_bounds
    if horizon == n_obs:
        tail = "there are no periods after the observed window"
    else:
        survival_text = ", times the survival weights S(L + h) / S(L)," if has_survival else ","
        regimes_text = f"{rho_upper:g} (Bound U), {rho_lower:g} (Bound L) and 0 (floor)"
        if extra_labels:
            regimes_text = (
                f"{rho_upper:g} (Bound U), {rho_lower:g} (Bound L), 0 (floor) and the given rho of each extra regime "
                f"({', '.join(extra_labels)}; fixed or one per draw)"
            )
        tail = (
            f"periods {n_obs + 1}-{horizon} are scenarios, not data: {level_name} keeps a share rho of its level "
            f"per period{pulse_text}{survival_text} with rho = {regimes_text}; Bounds U and L are pre-registered "
            "scenarios, not statistical bounds, and faster decay is not excluded"
        )
    persistence = "full persistence" if rho_upper == 1.0 else f"persistence at rho = {rho_upper:g} (Bound U)"
    discount_text = "discounted" if discounted else "undiscounted"
    return (
        f"Periods 1-{n_obs} are posterior draws; {tail}. Values are effect x (alpha s_t + 1 - s_t), "
        f"{discount_text}. X counts the periods of {persistence} after period {n_obs} needed for the cumulative "
        "value to exceed the cost."
    )


def compute_required_persistence(
    traj_draws: Any,
    *,
    cost_grid: Sequence[float],
    alpha_grid: Sequence[float],
    horizon: int = 60,
    alpha_share: float | Sequence[float] | np.ndarray = 1.0,
    rho_bounds: tuple[float, float] = (0.85, 1.00),
    extra_regimes: Mapping[str, float | np.ndarray] | None = None,
    headline_regimes: Sequence[str] = ("U", "L"),
    pulse_decay_tau: float | None = None,
    plateau_arm: str = "BC",
    pulse_arm: str = "AB",
    total_arm: str = "AC",
    survival: np.ndarray | None = None,
    discount: np.ndarray | None = None,
    discount_rate: float | None = None,
    checkpoints: Sequence[int] = DEFAULT_PERSISTENCE_CHECKPOINTS,
    milestones: Sequence[int] = DEFAULT_MILESTONES,
    placebo_draws: Any = None,
    placebo_extra_regimes: Mapping[str, float | np.ndarray] | None = None,
    quantiles: Sequence[float] = DEFAULT_QUANTILES,
) -> RequiredPersistenceResult:
    """Invert a post-horizon ROI question into required persistence, break-even retention and payback.

    The value layer on top of :func:`extrapolate_rank_k_decay`: every regime extends the observed
    draws with the same tail model -- the plateau keeps ``rho`` of its level per period (Bound U
    ``rho_bounds[1]``, Bound L ``rho_bounds[0]``, floor 0 and each retention of ``extra_regimes``), the
    pulse decays as ``exp(-h / pulse_decay_tau)`` and the survival weights multiply both -- and
    prices each period as ``v(t) = tau(t) * m_t(alpha) * delta_t`` with
    ``m_t(alpha) = alpha * s_t + (1 - s_t)``, where ``s_t`` (``alpha_share``) is the share of the effect
    whose incrementality is ``alpha`` (the rest counts in full) and ``delta_t`` the discount weight.
    ``V(m)`` is the running sum and the payback month is the first ``m`` with ``V(m) > cost`` (``H + 1``
    when never).

    * Required persistence ``X``: periods after ``L`` that the persisting level must last along
      Bound U (full persistence with the default ``rho_U = 1``; persistence at ``rho_U`` otherwise, as
      the labels state) for ``V`` to exceed the cost, with linear interpolation inside the payback
      period; 0 when the observed window already pays back, ``inf`` when not even ``H`` does. With
      ``pulse_decay_tau=None`` the whole level ``tau(L)`` persists (the upper envelope).
    * Break-even retention: the smallest plateau ``rho`` in ``[0, 1]`` (grid step ``5e-5``) with
      ``V_rho(H) > cost``, solved in closed form because ``V_rho(H)`` is affine in the plateau level,
      and per extra regime the share of draws whose retention reaches it.
    * Scenario tables: payback CDFs and cumulative-value quantiles for U, L, ``rho = 0`` and the extra
      regimes, and ``Pr(ROI > 0) = Pr(V(m) > cost)`` at ``L``, ``L + k`` and ``H``. Every regime has a
      role: ``"headline"`` (``"observed"`` and ``headline_regimes``), ``"diagnostic"`` (the other
      regimes) or ``"reference"`` (``rho0``).
    * Placebo gate: ``placebo_draws`` (same structure, with its own plateau / pulse split) go through
      the identical machinery, the extra regimes with ``placebo_extra_regimes``; a cell is
      decision-immune when its placebo ``Pr(ROI > 0)`` is below 0.5, and the immune cells are counted per
      regime, for ``"observed"`` with the headline regimes, for the diagnostic regimes and overall. Pass
      the placebo draws as ``traj_draws`` of a separate call for their own ``X`` table.

    Hard assertions: (i) the cumulative effect at ``L`` is equal across regimes and to the observed
    sum; (ii) the split is additive draw by draw; when every monthly effect is ``>= 0``, (iii)
    ``Pr(X <= k) = Pr(payback_U <= L + k)`` for every integer ``k <= H - L`` and (iv)
    ``Pr(payback_U <= m) >= Pr(payback_L <= m) >= Pr(payback_rho0 <= m)``, and every extra regime's
    payback CDF at least the floor's, for every ``m`` -- otherwise the cumulative value is not
    monotone, payback is a first crossing, ``cdf_identity_checked`` is ``False`` and a warning is
    emitted; (v) no NaN / Inf; (vi) a data-health certificate for every ``(cost, alpha)``. With
    ``alpha_share = 1`` and neither discounting nor survival, ``Pr(ROI > 0)`` at ``H`` reproduces
    :func:`roi_break_even_table` on the cumulative effect of the same regime.

    Payback, ``X`` and the break-even retention all use the strict rule ``V > cost``: at an exact tie
    the cost is not yet recovered, and ``X`` is then moved to the next float above its integer value
    so that ``Pr(X <= k)`` agrees with the payback CDF. There is no random step and hence no seed: the
    assumed ``Beta`` band of :func:`extrapolate_rank_k_decay` is deliberately not carried into the
    value layer, whose scenarios are the labelled Bounds U and L, the floor and the extra regimes.

    Args:
      traj_draws: Observed trajectory draws: one ``(S, L)`` array, ``{arm: (S, L)}`` or a
        :class:`Nested3ArmReconciliationResult`. The split needs the plateau, pulse and total arms;
        without it the total arm (or the only arm) is used.
      cost_grid: Costs per unit (outcome units, ``> 0``, distinct).
      alpha_grid: Incrementality shares in ``[0, 1]`` (distinct).
      horizon: Decision horizon ``H >= L`` in periods.
      alpha_share: Share ``s_t`` in ``[0, 1]`` of each period's effect that ``alpha`` applies to: a
        scalar or an ``(H,)`` array (for example ``[0.9] * 6 + [0.5] * (H - 6)``). 1 prices the whole
        effect at ``alpha``.
      rho_bounds: ``(Bound L, Bound U)`` plateau retentions with ``0 <= L <= U <= 1``.
      extra_regimes: Optional ``{label: rho}`` of further plateau retentions, each run through the
        trajectory, payback, value, decision and placebo machinery of Bounds U and L: a scalar, or an
        ``(S,)`` array with one retention per draw, paired with ``traj_draws`` by index; finite and in
        ``[0, 1]``. Labels are non-empty strings other than ``"U"``, ``"L"``, ``"rho0"``,
        ``"observed"`` and the placebo-summary scopes (``"observed+headline"``, ``"diagnostic"``,
        ``"reference"``, ``"all"``). ``X`` stays on Bound U.
      headline_regimes: Regimes with the ``"headline"`` role besides ``"observed"``: distinct labels
        among ``"U"``, ``"L"`` and the extra regimes (a single label may be passed as a string). The
        other regimes are ``"diagnostic"`` and ``rho0`` is the ``"reference"``.
      pulse_decay_tau: Time constant (periods, ``> 0``) of the pulse decay, or ``None`` for the
        whole-level upper envelope.
      plateau_arm: Label of the plateau arm under the split.
      pulse_arm: Label of the pulse arm under the split.
      total_arm: Label of the total arm.
      survival: Optional ``(H,)`` survival weights, finite, in ``[0, 1]`` and exactly 1 on
        ``t <= L`` (see :func:`weibull_survival_weights`).
      discount: Optional ``(H,)`` positive discount weights (see :func:`discount_weights`).
      discount_rate: Convenience annual rate for ``discount_weights(horizon, annual_rate=...)`` with
        12 periods per year; exclusive with ``discount``. Undiscounted when both are ``None``.
      checkpoints: Periods ``k`` after ``L`` at which ``Pr(X <= k)`` is reported (positive integers;
        any ``k > horizon - L`` yields NaN in ``pr_X_le_<k>`` and emits a warning).
      milestones: Months of the scenario table; clipped to ``horizon`` and always including it.
      placebo_draws: Optional null-placebo draws with the structure of ``traj_draws`` and the same
        ``L``, or ``None`` to skip the gate.
      placebo_extra_regimes: Retentions of the extra regimes for the placebo draws, with the keys of
        ``extra_regimes``: a scalar, or one per placebo draw paired by index (for example the study's
        own per-draw retention, to gate the decision under the same scenario). Requires
        ``placebo_draws``. When it is omitted, the extra regimes skip the placebo gate -- missing
        placebo columns, left out of the counts -- and a warning says so.
      quantiles: Probabilities of the quantile columns (odd length, default P30/P50/P70).

    Returns:
      :class:`RequiredPersistenceResult`.
    """
    roles_requested = (str(plateau_arm), str(pulse_arm), str(total_arm))
    tau_p = None if pulse_decay_tau is None else _positive_scalar(pulse_decay_tau, "pulse_decay_tau")
    draws, roles, total_key = _value_arms(traj_draws, tau_p, roles_requested, "traj_draws")
    horizon, n_obs = _check_periods(horizon, draws[total_key].shape[1])
    costs = _decision_grid(cost_grid, "cost_grid", positive=True)
    alphas = _decision_grid(alpha_grid, "alpha_grid", positive=False)
    shares = _alpha_share_vector(alpha_share, horizon)
    rho_lower, rho_upper = _check_rho_bounds(rho_bounds)
    extra_rhos = _extra_regime_rhos(extra_regimes, draws[total_key].shape[0], "extra_regimes")
    extra_labels = tuple(extra_rhos)
    regimes = (*_REGIMES, *extra_labels)
    regime_roles = _regime_roles(extra_labels, _resolve_headline_regimes(headline_regimes, extra_labels))
    survival_tail = _survival_tail(survival, horizon, n_obs)
    discount_vec = _discount_vector(discount, discount_rate, horizon)
    ks = _check_checkpoints(checkpoints)
    milestone_months = _resolve_milestones(milestones, horizon)
    qs = _check_quantiles(quantiles)
    placebo = None
    if placebo_draws is not None:
        placebo = _value_arms(placebo_draws, tau_p, roles_requested, "placebo_draws")
        placebo_periods = placebo[0][placebo[2]].shape[1]
        if placebo_periods != n_obs:
            raise ValueError(
                f"placebo_draws must span the {n_obs} observed periods of traj_draws; got {placebo_periods}"
            )
    placebo_extra: dict[str, float | np.ndarray] = {}
    if placebo_extra_regimes is not None:
        if placebo is None:
            raise ValueError("placebo_extra_regimes holds the extra regimes' retentions for placebo_draws; pass both")
        placebo_draw_count = placebo[0][placebo[2]].shape[0]
        placebo_extra = _extra_regime_rhos(placebo_extra_regimes, placebo_draw_count, "placebo_extra_regimes")
        if set(placebo_extra) != set(extra_labels):
            raise ValueError(
                f"placebo_extra_regimes must have the keys of extra_regimes {sorted(extra_labels)}; "
                f"got {sorted(placebo_extra)}"
            )
    emitted = []
    beyond = [k for k in ks if k > horizon - n_obs]
    if beyond:
        message = (
            f"compute_required_persistence(): checkpoints {beyond} exceed the extrapolated window "
            f"horizon - L = {horizon - n_obs} (horizon = {horizon}, L = {n_obs}), so X is censored at "
            f"{horizon - n_obs} and their pr_X_le_<k> columns are NaN; pass horizon >= L + max(checkpoints) "
            f"({n_obs + max(beyond)}) or checkpoints <= {horizon - n_obs}"
        )
        warnings.warn(message, UserWarning, stacklevel=2)
        emitted.append(message)
    if placebo is not None and extra_labels and placebo_extra_regimes is None:
        message = (
            "compute_required_persistence(): placebo_extra_regimes was not given, so the extra regimes "
            f"{list(extra_labels)} skip the placebo gate: their placebo columns are missing and they are left out "
            "of the placebo counts; pass placebo_extra_regimes with the same keys to gate them"
        )
        warnings.warn(message, UserWarning, stacklevel=2)
        emitted.append(message)

    regime_rhos: dict[str, float | np.ndarray] = {"U": rho_upper, "L": rho_lower, "rho0": 0.0, **extra_rhos}
    extend = functools.partial(
        _extend_arms, horizon=horizon, roles=roles, pulse_decay_tau=tau_p, survival_tail=survival_tail
    )
    extended = {regime: extend(draws, rho_regime) for regime, rho_regime in regime_rhos.items()}
    monthly = {regime: arms[total_key] for regime, arms in extended.items()}
    value_weights = {alpha: (alpha * shares + (1.0 - shares)) * discount_vec for alpha in alphas}
    cumulative = _cumulative_values(monthly, value_weights)
    payback = {
        (regime, cost, alpha): _payback_months(cumulative[(regime, alpha)], cost)
        for regime, cost, alpha in itertools.product(regimes, costs, alphas)
    }
    plateau_level = draws[total_key if roles is None else roles[0]][:, -1]
    pulse_level = None if roles is None else draws[roles[1]][:, -1]
    grid = np.linspace(0.0, 1.0, _BREAKEVEN_GRID_POINTS)
    upper_values = {alpha: monthly["U"] * value_weights[alpha][None, :] for alpha in alphas}
    breakeven_parts = {
        alpha: _breakeven_parts(draws[total_key], pulse_level, value_weights[alpha], survival_tail, tau_p, grid)
        for alpha in alphas
    }
    required = {}
    breakeven = {}
    for cost, alpha in itertools.product(costs, alphas):
        required[(cost, alpha)] = _required_persistence(upper_values[alpha], cumulative[("U", alpha)], cost, n_obs)
        observed_value, pulse_value, plateau_value = breakeven_parts[alpha]
        breakeven[(cost, alpha)] = _breakeven_rho(
            cost, observed_value, pulse_value, plateau_level, plateau_value, grid
        )
    placebo_cumulative = None
    placebo_rhos: dict[str, float | np.ndarray] | None = None
    if placebo is not None:
        placebo_arms, _, placebo_total = placebo
        placebo_rhos = {"U": rho_upper, "L": rho_lower, "rho0": 0.0}
        placebo_rhos.update((regime, placebo_extra[regime]) for regime in extra_labels if regime in placebo_extra)
        placebo_extended = {regime: extend(placebo_arms, rho_regime) for regime, rho_regime in placebo_rhos.items()}
        if roles is not None:
            for regime, arms in placebo_extended.items():
                _assert_split_additivity(arms, roles, f"placebo regime {regime}")
        placebo_monthly = {regime: arms[placebo_total] for regime, arms in placebo_extended.items()}
        placebo_cumulative = _cumulative_values(placebo_monthly, value_weights)

    if roles is not None:
        for regime, arms in extended.items():
            _assert_split_additivity(arms, roles, f"regime {regime}")
    _check_value_invariants(draws, total_key, monthly, cumulative, required, breakeven)
    nonnegative = _all_nonnegative(extended)
    if nonnegative:
        _check_payback_identities(required, payback, regimes, n_obs, horizon)
    else:
        message = (
            "compute_required_persistence(): some monthly effects are negative, so the cumulative value is not "
            "monotone and payback months and X are first-crossing times; the identity Pr(X <= k) = "
            "Pr(payback_U <= L + k) and the ordering U >= L >= rho0 were not checked (cdf_identity_checked=False)"
        )
        warnings.warn(message, UserWarning, stacklevel=2)
        emitted.append(message)

    split = roles is not None
    upper_text = "in full"
    if rho_upper != 1.0:
        upper_text = f"at rho = {rho_upper:g} per period (Bound U, not full persistence)"
    if split:
        pulse_decay = f"pulse {roles[1]!r} decays as exp(-h / {tau_p:g})"
        level_name = f"the plateau {roles[0]!r}"
        pulse_text = f" while the {pulse_decay}"
        persisting_level = f"plateau {roles[0]!r}, persisting {upper_text}; {pulse_decay}"
    else:
        level_name = f"the whole level of {total_key!r}"
        pulse_text = ""
        persisting_level = f"whole level {total_key!r}, persisting {upper_text} (upper envelope)"
    report_months = tuple(sorted({n_obs, horizon} | {min(n_obs + k, horizon) for k in ks}))
    payback_df, payback_cdf_df = _payback_tables(payback, regime_rhos, regime_roles, report_months, qs, n_obs, horizon)
    decision_df = _decision_table(
        cumulative, placebo_cumulative, regimes, regime_roles, costs, alphas, report_months, n_obs
    )
    monthly_frames = [
        _per_period_table(monthly[regime], qs, n_obs, regime=regime, role=regime_roles[regime]) for regime in regimes
    ]
    value_frames = [
        _per_period_table(arr, qs, n_obs, regime=regime, role=regime_roles[regime], alpha=alpha)
        for (regime, alpha), arr in cumulative.items()
    ]
    certificate = _value_certificate(
        required, breakeven, payback, cumulative, placebo_cumulative, regime_roles, regimes, n_obs, horizon
    )
    if set(zip(certificate["cost"], certificate["alpha"])) != set(required):
        raise AssertionError("the data-health certificate does not cover every (cost, alpha) cell")
    role_names = ("plateau", "pulse", "total") if split else ("total",)
    return RequiredPersistenceResult(
        n_observed=n_obs, horizon=horizon, cost_grid=costs, alpha_grid=alphas, rho_bounds=(rho_lower, rho_upper),
        pulse_decay_tau=tau_p,
        arm_roles=dict(zip(roles if split else (total_key,), role_names, strict=True)),
        alpha_share=shares,
        survival=np.ones(horizon) if survival_tail is None else np.concatenate([np.ones(n_obs), survival_tail]),
        discount=discount_vec, checkpoints=ks, report_months=report_months, quantiles=qs, regime_rhos=regime_rhos,
        regime_labels=_regime_labels(rho_lower, rho_upper, split, extra_rhos, regime_roles),
        regime_roles=regime_roles, placebo_regime_rhos=placebo_rhos, monthly_draws=monthly,
        cumulative_value_draws=cumulative,
        required_persistence_draws=required, breakeven_rho_draws=breakeven, payback_month_draws=payback,
        required_persistence_df=_required_persistence_table(
            required, cumulative, shares, plateau_level, persisting_level, ks, qs, n_obs
        ),
        breakeven_rho_df=_breakeven_table(breakeven, extra_rhos, qs, horizon, split), payback_df=payback_df,
        payback_cdf_df=payback_cdf_df, monthly_summary_df=pd.concat(monthly_frames, ignore_index=True),
        cumulative_value_df=pd.concat(value_frames, ignore_index=True),
        milestone_table_df=_value_milestone_table(cumulative, payback, regime_roles, milestone_months, qs, n_obs),
        decision_df=decision_df,
        placebo_summary_df=None if placebo_cumulative is None else _placebo_summary(decision_df),
        data_health_certificate_df=certificate, cdf_identity_checked=nonnegative, warnings=tuple(emitted),
        scenario_note=_value_scenario_note(
            n_obs, horizon, level_name, pulse_text, (rho_lower, rho_upper), extra_labels, survival_tail is not None,
            bool(np.any(discount_vec != 1.0)),
        ),
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


def _decay_tail_label(role: str, decay_result: DecayExtrapolationResult) -> str:
    """Legend text of an arm's tail in :func:`plot_rank_k_decay_extrapolation`."""
    if role == "pulse":
        return f"exp(-h / {decay_result.pulse_decay_tau:g}) tail"
    if role == "total":
        return f"plateau rho = {decay_result.rho:g} + pulse tail"
    return f"rho = {decay_result.rho:g} tail"


def plot_rank_k_decay_extrapolation(
    decay_result: DecayExtrapolationResult,
    *,
    prior_decay_result: DecayExtrapolationResult | None = None,
    unit: str = "",
    milestones: Sequence[int] = DEFAULT_MILESTONES,
    axes: Any = None,
    title: str | None = None,
) -> np.ndarray:
    """Two-panel decay figure: monthly-rate ribbons by arm and the cumulative fan with the ``rho`` scenarios.

    Panel A draws the P30-P70 ribbon and P50 line of the monthly rate per arm (solid over
    the observed window, which is shaded; dashed over the extrapolated tail, labelled by its
    tail model under the plateau / pulse split), with the matching prior ribbon in grey when
    ``prior_decay_result`` is given. Panel B draws the cumulative P30-P70 fan per arm under the
    fixed ``rho``, the wider ``Beta``-``rho`` band when available -- labelled as an "assumed
    Beta(rho, sd) scenario band — not identified by the data" -- thin dotted P50 lines for
    each decay of the ``rho`` grid and light vertical markers at the milestones.

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
    roles = decay_result.arm_roles or {}
    grid_rhos = sorted(decay_result.rho_grid_cumulative)
    grid_label = f"rho grid {grid_rhos[0]:g}-{grid_rhos[-1]:g}" if grid_rhos else ""
    band_label = None
    if decay_result.random_rho_cumulative_draws is not None:
        band_label = f"{_RHO_BAND_LABEL} (rho {decay_result.rho:g}, sd {decay_result.rho_sd:g}; one per arm)"

    def _arm_table(frame: pd.DataFrame, arm: str) -> pd.DataFrame:
        return frame.loc[frame["arm"] == arm].sort_values("month")

    band_drawn = False
    for i, arm in enumerate(decay_result.arms):
        color = _ARM_COLORS[i % len(_ARM_COLORS)]
        role = roles.get(arm, "geometric")
        tail_label = _decay_tail_label(role, decay_result)
        monthly = _arm_table(decay_result.monthly_summary_df, arm)
        p30, p50, p70 = (monthly[c].to_numpy() for c in ("P30", "P50", "P70"))
        rate_ax.fill_between(months_all, p30, p70, color=color, alpha=0.22)
        rate_ax.plot(months_all[:n_obs], p50[:n_obs], color=color, lw=1.8, label=f"{arm} P50 (observed)")
        rate_ax.plot(
            months_all[n_obs - 1 :], p50[n_obs - 1 :], color=color, lw=1.4, ls="--", label=f"{arm} P50 ({tail_label})",
        )
        if prior_decay_result is not None:
            prior = _arm_table(prior_decay_result.monthly_summary_df, arm)
            rate_ax.fill_between(
                months_all, prior["P30"].to_numpy(), prior["P70"].to_numpy(), color="#9e9e9e", alpha=0.18,
                label=f"{arm} prior P30-P70" if i == 0 else None,
            )
        cum = _arm_table(decay_result.cumulative_summary_df, arm)
        cum_ax.fill_between(months_all, cum["P30"].to_numpy(), cum["P70"].to_numpy(), color=color, alpha=0.25)
        cum_ax.plot(
            months_all, cum["P50"].to_numpy(), color=color, lw=1.8, label=f"{arm} P50 [P30, P70], {tail_label}"
        )
        if band_label is not None and role != "pulse":
            cum_ax.fill_between(
                months_all, cum["rho_band_P30"].to_numpy(), cum["rho_band_P70"].to_numpy(), color=color, alpha=0.10,
                label=band_label if not band_drawn else None,
            )
            band_drawn = True
        if role != "pulse":
            for k, rho_k in enumerate(grid_rhos):
                cum_ax.plot(
                    months_all, np.median(decay_result.rho_grid_cumulative[rho_k][arm], axis=0), color=color, lw=0.8,
                    ls=":", alpha=0.7, label=f"{arm} {grid_label}" if k == 0 else None,
                )
    rate_ax.axvspan(0.5, n_obs + 0.5, color="#f0f0f0", zorder=0, label="Observed window")
    rate_ax.axhline(0.0, color="black", lw=0.8, alpha=0.5)
    rate_ax.set_xlabel("Month since adoption")
    rate_ax.set_ylabel(f"Monthly effect{_unit_suffix(unit)}")
    rate_ax.set_title("A. Monthly rate: observed window and assumed decay tail")
    rate_ax.legend(frameon=False, fontsize=7)
    for m in marks:
        cum_ax.axvline(m, color="grey", lw=0.6, ls=":", alpha=0.6)
    cum_ax.set_xlabel("Month since adoption")
    cum_ax.set_ylabel(f"Cumulative effect{_unit_suffix(unit)}")
    cum_ax.set_title("B. Cumulative fan with assumed rho scenarios (milestones dotted)")
    cum_ax.legend(frameon=False, fontsize=7)
    if title:
        axs[0].figure.suptitle(title)
    return axs


def _grid_member(value: Any, grid: Sequence[float], name: str) -> float:
    """The entry of ``grid`` equal to ``value`` (to ``1e-12``), raising ``ValueError`` when there is none."""
    try:
        target = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a number on the grid {tuple(grid)}; got {value!r}") from exc
    for entry in grid:
        if math.isclose(entry, target, rel_tol=1e-12, abs_tol=1e-12):
            return entry
    raise ValueError(f"{name}={value!r} is not on the result's grid {tuple(grid)}")


def _band_columns(quantiles: Sequence[float]) -> tuple[str, str, str]:
    """``(low, central, high)`` quantile column labels."""
    return (
        _quantile_label(quantiles[0]),
        _quantile_label(quantiles[len(quantiles) // 2]),
        _quantile_label(quantiles[-1]),
    )


def _draw_regime_effects(ax: Any, result: RequiredPersistenceResult, unit: str) -> None:
    """Panel A of :func:`plot_required_persistence_and_payback`: observed band and the regime tails."""
    table = result.monthly_summary_df
    n_obs = result.n_observed
    low, mid, high = _band_columns(result.quantiles)
    ax.axvspan(0.5, n_obs + 0.5, color="#f0f0f0", zorder=0, label="Observed window")
    observed = table.loc[(table["regime"] == "U") & (table["month"] <= n_obs)]
    ax.fill_between(observed["month"], observed[low], observed[high], color="#202124", alpha=0.15, lw=0)
    ax.plot(observed["month"], observed[mid], color="#202124", lw=1.8, label=f"Observed {mid} [{low}, {high}]")
    for regime in _REGIMES:
        color, style = _REGIME_STYLES[regime]
        tail = table.loc[(table["regime"] == regime) & (table["month"] >= n_obs)]
        ax.plot(tail["month"], tail[mid], color=color, ls=style, lw=1.6, label=result.regime_labels[regime])
    ax.axhline(0.0, color="black", lw=0.8, alpha=0.5)
    ax.set_xlim(0.5, result.horizon + 0.5)
    ax.set_xlabel("Period since adoption")
    ax.set_ylabel(f"Monthly effect{_unit_suffix(unit)}")
    ax.set_title(f"A. Monthly effect {mid} by regime (tails are scenarios)")
    ax.legend(frameon=False, fontsize=6.5, loc="best")


def _draw_payback_cdfs(ax: Any, result: RequiredPersistenceResult, cost: float, alpha: float) -> None:
    """Panel B: payback CDFs against the periods after ``L``, with the interpolated ``X`` and its quantiles."""
    cdf = result.payback_cdf_df
    n_obs = result.n_observed
    level = "plateau" if result.pulse_decay_tau is not None else "whole-level"
    cell = cdf.loc[(cdf["cost"] == cost) & (cdf["alpha"] == alpha) & (cdf["month"] >= n_obs)]
    for regime in _REGIMES:
        color, style = _REGIME_STYLES[regime]
        rows = cell.loc[cell["regime"] == regime]
        ax.step(
            rows["months_after_observed"], rows["pr_paid_back"], where="post", color=color, ls=style, lw=1.8,
            label=f"{regime}: {level} rho = {result.regime_rhos[regime]:g} (scenario)",
        )
    # Pr(X <= x) for fractional x: it meets the Bound-U payback CDF at every integer and carries the markers.
    upper_color = _REGIME_STYLES["U"][0]
    persistence = result.required_persistence_draws[(cost, alpha)]
    finite = np.sort(persistence[np.isfinite(persistence)])
    tail_end = max(1.0, float(result.horizon - n_obs))
    shares = np.arange(1, finite.size + 1) / persistence.size
    ax.step(
        np.r_[0.0, finite, tail_end], np.r_[0.0, shares, finite.size / persistence.size], where="post",
        color=upper_color, lw=0.9, alpha=0.7, label="X along U (interpolated): Pr(X <= x)",
    )
    summary = result.required_persistence_df
    row = summary.loc[(summary["cost"] == cost) & (summary["alpha"] == alpha)].iloc[0]
    for q in result.quantiles:
        label = _quantile_label(q)
        value = float(row[f"X_{label}"])
        if math.isfinite(value):
            ax.plot([value], [q], marker="o", ms=6, color=upper_color, mfc="white", ls="none")
            ax.annotate(
                f"X {label} = {value:.1f}", xy=(value, q), xytext=(6, -10), textcoords="offset points", fontsize=7,
                color=upper_color,
            )
    ax.set_xlim(0.0, tail_end)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(f"Periods after L = {n_obs}")
    ax.set_ylabel("Pr(paid back by then)")
    ax.set_title(f"B. Payback CDF at cost {cost:g}, alpha {alpha:g} (X read off U)")
    ax.legend(frameon=False, fontsize=7, loc="lower right")


def _draw_cumulative_value(ax: Any, result: RequiredPersistenceResult, cost: float, alpha: float, unit: str) -> None:
    """Panel C: cumulative value band per regime against the cost lines."""
    table = result.cumulative_value_df
    low, mid, high = _band_columns(result.quantiles)
    for regime in _REGIMES:
        color, style = _REGIME_STYLES[regime]
        rows = table.loc[(table["regime"] == regime) & (table["alpha"] == alpha)]
        ax.fill_between(rows["month"], rows[low], rows[high], color=color, alpha=0.12, lw=0)
        ax.plot(rows["month"], rows[mid], color=color, ls=style, lw=1.8, label=f"{regime} {mid} [{low}, {high}]")
    for entry in result.cost_grid:
        highlighted = entry == cost
        ax.axhline(
            entry, color="black", lw=1.4 if highlighted else 0.8, ls="-." if highlighted else ":",
            label=f"cost {entry:g}" + (" (highlighted)" if highlighted else ""),
        )
    ax.axvline(result.n_observed + 0.5, color="grey", lw=0.8, ls="--")
    discounted = "discounted " if bool(np.any(result.discount != 1.0)) else ""
    ax.set_xlim(0.5, result.horizon + 0.5)
    ax.set_xlabel("Period since adoption")
    ax.set_ylabel(f"Cumulative {discounted}value{_unit_suffix(unit)}")
    ax.set_title(f"C. Cumulative {discounted}value at alpha {alpha:g}")
    ax.legend(frameon=False, fontsize=7, loc="upper left")


def plot_required_persistence_and_payback(
    result: RequiredPersistenceResult,
    *,
    cost_highlight: float,
    alpha_highlight: float,
    unit: str = "",
    axes: Any = None,
    title: str | None = None,
) -> np.ndarray:
    """Three-panel decision figure: regime trajectories, payback CDFs with ``X``, cumulative value vs cost.

    Panel A draws the monthly-effect P50 of the observed window (with its P30-P70 band, the window
    shaded) and of each regime's tail -- Bound U, Bound L and the dotted ``rho = 0`` floor -- with
    the scenario caveats in the legend. Panel B draws the payback CDFs of the three regimes against
    the periods after ``L`` at ``(cost_highlight, alpha_highlight)`` and marks the quantiles of the
    required persistence ``X`` on the Bound-U curve (``Pr(X <= k)`` is that curve at ``k``). Panel C
    draws the cumulative (discounted) value P50 with its P30-P70 band per regime against the cost
    lines, the highlighted cost in bold. ``matplotlib`` is imported here rather than at module load.

    Args:
      result: Output of :func:`compute_required_persistence`.
      cost_highlight: Cost of panels B and C; must be on ``result.cost_grid``.
      alpha_highlight: Incrementality of panels B and C; must be on ``result.alpha_grid``.
      unit: Outcome unit for the axis labels.
      axes: Sequence of three Axes, or ``None`` to create a one-row figure.
      title: Figure-level title, or ``None``.

    Returns:
      Length-3 object array of Axes.
    """
    plt = _pyplot("plot_required_persistence_and_payback")
    cost = _grid_member(cost_highlight, result.cost_grid, "cost_highlight")
    alpha = _grid_member(alpha_highlight, result.alpha_grid, "alpha_highlight")
    axs = _axes_array(plt, axes, 3, (19.0, 5.0))
    _draw_regime_effects(axs[0], result, unit)
    _draw_payback_cdfs(axs[1], result, cost, alpha)
    _draw_cumulative_value(axs[2], result, cost, alpha, unit)
    if title:
        axs[0].figure.suptitle(title)
    return axs


__all__ = [
    "BUNDLE_PRIOR_MODES",
    "DEFAULT_MILESTONES",
    "DEFAULT_PERSISTENCE_CHECKPOINTS",
    "DEFAULT_QUANTILES",
    "DecayExtrapolationResult",
    "NESTED_CONTRASTS",
    "PRIOR_FAMILIES",
    "PriorSensitivityResult",
    "PriorTier",
    "RankKPriorDraws",
    "RequiredPersistenceResult",
    "assert_nested_additivity",
    "compute_required_persistence",
    "discount_weights",
    "evaluate_prior_sensitivity",
    "extrapolate_rank_k_decay",
    "nested_3arm_trajectories",
    "plot_rank_k_decay_extrapolation",
    "plot_rank_k_prior_calibration",
    "plot_rank_k_prior_sensitivity",
    "plot_rank_k_prior_vs_posterior",
    "plot_required_persistence_and_payback",
    "prior_vs_posterior_summary",
    "roi_break_even_table",
    "sample_rank_k_prior_draws",
    "standardized_prior_family_table",
    "weibull_survival_weights",
    "ytlr_prior_tiers",
]
