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

"""Randomized encouragement with absorbing adoption: the adoption-clock model.

Two LongBet equations are fitted and composed:

* **Outcome on the adoption clock.** ``y`` is regressed on the *observed*
  adoption ``d`` with LongBet's exposure clock counting periods since
  adoption, so the outcome equation contains the jump that adoption causes,
  its exposure profile ``beta_S`` and its covariate heterogeneity ``nu(x)``,
  plus a unit intercept and, for a continuous outcome, stationary AR(1)
  errors. Both arms contribute adoption events.
* **Adoption as a discrete-time hazard.** A binary LongBet is fitted to the
  first-adoption indicator on the at-risk set (cells after adoption are
  masked out), with the randomized encouragement ``z`` as the treatment and
  periods since the offer as its clock. The hazard is population-averaged
  (no unit frailty): its survival product is the marginal adoption-time
  distribution, which is what the composition needs.

In both equations the moderator is a horseshoe linear term plus trees and the
exposure trajectory's scale is learned (:data:`MODEL_COMPONENTS`).

The offer's effects follow by composition: the hazard equation gives each
unit's adoption-time distribution with and without the offer, the outcome
equation gives the effect of having adopted ``s`` periods ago, and the offer
effect on the outcome is the difference of the expected exposure responses.
For a binary outcome the unit intercept is integrated out exactly,
``E[Phi(a + sigma Z)] = Phi(a / sqrt(1 + sigma^2))``.

This uses the exclusion restriction (the offer moves outcomes only through
adoption) and the assumption that, given covariates and the unit intercept,
the timing of adoption is not confounded with the outcome innovations. The
model-free design-based reference in ``reference``/``encouragement_effects``
needs neither and is the check on the population offer effects.

A ``reduced-form`` alternative (regressing ``y`` on ``z`` directly) was
retired: an adoption jump at a random period that the outcome equation does
not model is unit-specific, time-structured residual noise, and flexible
forests fit it, so the exact posterior of that model was less accurate out of
sample than early-stopped chains. See ``benchmarks/tempering_plan.md``.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import tempfile
import warnings
from pathlib import Path
from typing import Any

import arviz as az
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.special import ndtr

from longbet._config import LongBetConfig, resolve_model_defaults
from longbet._diagnostics import compute_ess, compute_rhat
from longbet._encourage import _critical, _validate, encouragement_effects
from longbet._encourage_intensity import (
    IDENTIFICATION_TIERS,
    STAGE2_BACKENDS,
    TIER_LABELS,
    IntensityReference,
    LinearLagStage2,
    LongBetLagStage2,
    compose_offer_effect,
    control_function_columns,
    cumulative_return_draws,
    default_lag_horizon,
    intensity_reference,
    lag_matrix,
    offer_induced_dose,
    organic_dose_columns,
    pre_offer_columns,
    pre_offer_demean,
    standardized_covariates,
    validate_intensity,
)
from longbet._model import LongBet, derive_exposure

ARCHIVE_VERSION = 2
INFERENCE_VERSION = "adoption_clock_v2"
# Provenance of the two equations' random streams, recorded by every fit: the
# composition pairs draw b of the two fits, which samples the product of the
# two posteriors only if their chains are independent (engine releases before
# this marker seeded both equations identically).
COMPONENT_STREAMS = "independent"
ENGINE = "adoption_clock"
#: The continuous-intensity model (``treatment="intensity"``): Stage 1 is a
#: continuous LongBet for the dose on the offer clock, Stage 2 the
#: distributed-lag kernel for the outcome; see ``_encourage_intensity``.
ENGINE_INTENSITY = "distributed_lag"
INFERENCE_VERSION_INTENSITY = "distributed_lag_v1"
ARCHIVE_VERSION_INTENSITY = 1
TREATMENTS: tuple[str, ...] = ("absorbing", "intensity")
DEFAULT_IDENTIFICATION = "control_function"
#: The encouragement model's components, fixed (LongBet-IV paper, Section 3;
#: chosen on development simulations and confirmed on held-out replicates):
#:
#: * ``linear_trt``: a horseshoe linear term in both moderators; the tree prior
#:   alone makes realistic heterogeneity a priori negligible.
#: * ``learn_beta_scale``: the trajectory's scale is learned, so the effect can
#:   shrink under a null and grow under a strong signal.
#: * ``interweave_gamma``: interweaving of the intercept variance (a mixing aid;
#:   same posterior).
#: * ``ar1_errors`` for a continuous outcome: serially correlated shocks are the
#:   rule in panels, and independent errors overstate a unit's information.
#:
#: The adoption hazard has no unit frailty: with one adoption event per unit a
#: frailty's variance is weakly identified, mixes poorly and flattens the
#: offer's effect on adoption across covariates.
MODEL_COMPONENTS: dict[str, str] = {
    "linear_trt": "the horseshoe linear term in the moderators",
    "learn_beta_scale": "the learned trajectory scale",
    "interweave_gamma": "the interweaving update of the intercept variance",
}
#: Thinning when ``n_skip`` is left unset (every validated fit thins by 2).
DEFAULT_N_SKIP = 2


def _finite_matrix(value: Any, name: str, n: int, shape: tuple | None = None) -> np.ndarray:
    a = np.asarray(value)
    if (a.dtype.kind not in "biuf" or a.ndim != 2 or a.shape[0] != n
            or (shape is not None and a.shape != shape)):
        raise ValueError(f"{name} must be a finite numeric matrix with "
                         + (f"shape {shape}." if shape else f"{n} rows."))
    a = np.array(a, dtype=np.float64, copy=True)
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be complete and finite.")
    return a


@dataclasses.dataclass(frozen=True)
class EncouragementDraws:
    """Aligned offer-effect draws and the population over which they average.

    ``draws[name]`` has axes ``(group, horizon, chain, retained_draw)``.
    ``weights[g]`` sum to one over the study units in group ``g``.
    """

    draws: dict[str, np.ndarray]
    group_labels: tuple[str, ...]
    group_counts: np.ndarray
    weights: np.ndarray
    periods: np.ndarray
    period_indices: np.ndarray
    horizons: np.ndarray
    standardization: str
    provenance: str

    @property
    def effects(self) -> np.ndarray:
        """Equivalent array with axes ``(chain, draw, group, horizon, outcome)``."""
        return np.stack(tuple(self.draws.values()), axis=-1).transpose(2, 3, 0, 1, 4)


def _target_weights(
    n_units: int, groups: Any, weights: Any,
) -> tuple[tuple[str, ...], np.ndarray, np.ndarray]:
    """Normalize fixed, nonnegative unit weights within each baseline group."""
    if weights is None:
        unit_weights = np.ones(n_units, dtype=np.float64)
    else:
        unit_weights = np.asarray(weights, dtype=np.float64)
        if (unit_weights.shape != (n_units,) or
                not np.all(np.isfinite(unit_weights)) or np.any(unit_weights < 0)):
            raise ValueError("weights must be a finite nonnegative vector of length N.")
        if unit_weights.max(initial=0) > 0:
            unit_weights = unit_weights / unit_weights.max()
    if unit_weights.sum() <= 0:
        raise ValueError("weights must have a positive total.")

    masks = [np.ones(n_units, dtype=bool)]
    labels = ["all"]
    if groups is not None:
        raw_groups = np.asarray(groups)
        if raw_groups.shape != (n_units,) or np.any(pd.isna(raw_groups)):
            raise ValueError("groups must be a length-N vector of nonmissing baseline labels.")
        try:
            values = list(pd.unique(raw_groups))
        except TypeError as exc:
            raise ValueError("groups must contain scalar baseline labels.") from exc
        for value in values:
            label = str(value)
            if label in labels:
                raise ValueError("group labels must be unique as strings and cannot be 'all'.")
            labels.append(label)
            masks.append(raw_groups == value)

    counts = np.asarray([mask.sum() for mask in masks], dtype=np.int64)
    target = np.stack([unit_weights * mask for mask in masks])
    totals = target.sum(axis=1)
    if np.any(totals <= 0):
        raise ValueError("every group must have positive total weight.")
    return tuple(labels), counts, target / totals[:, None]


def _reference_groups(data: dict, standardized: EncouragementDraws, alpha: float) -> pd.DataFrame:
    """Conditional-on-subgroup-arm-count reference, on the identical unit target."""
    frames = []
    for g, label in enumerate(standardized.group_labels):
        mask = standardized.weights[g] > 0
        z = data["z"][mask]
        counts = [int((z[:, -1] == a).sum()) for a in (0, 1)]
        if min(counts) > 0:
            frame = encouragement_effects(data["y"][mask], data["d"][mask], z,
                                          data["t"], alpha=alpha)
            n_post = len(frame)
            assigned = z[:, -1] == 1
            if min(counts) >= 2 and n_post > 0:
                for q_name, arr in (("itt_y", data["y"][mask][:, -n_post:]),
                                    ("itt_d", data["d"][mask][:, -n_post:])):
                    s_i = np.asarray(arr, dtype=np.float64).sum(axis=1)
                    v_cum = (np.var(s_i[assigned], ddof=1) / counts[1]
                             + np.var(s_i[~assigned], ddof=1) / counts[0])
                    frame[f"{q_name}_cum_se"] = float(np.sqrt(v_cum))
            else:
                frame["itt_y_cum_se"] = np.nan
                frame["itt_d_cum_se"] = np.nan
        else:
            # Retain unsupported groups. There is no observed arm mean to fill in.
            frame = encouragement_effects(data["y"], data["d"], data["z"],
                                          data["t"], alpha=alpha)
            for name in frame:
                if name.startswith(("itt_", "wald")):
                    frame[name] = np.nan
            frame["itt_y_cum_se"] = np.nan
            frame["itt_d_cum_se"] = np.nan
            frame["wald_set_type"] = "unavailable"
            frame["wald_reason"] = "missing_group_arm"
            frame["n_control"], frame["n_encouraged"] = counts
        frame.insert(0, "group", label)
        frame["n_units"] = int(mask.sum())
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Composition of the two equations
# ---------------------------------------------------------------------------

def _flat_draws(model: LongBet, name: str) -> np.ndarray:
    """A traced parameter with chains flattened chain-major into one draw axis."""
    a = np.asarray(getattr(model.trace, name))
    return a.reshape(-1, *a.shape[2:]) if a.ndim > 1 and model._chained else a


def _predict_quiet(model: LongBet, x: np.ndarray, z: np.ndarray, t: np.ndarray,
                   x_trt: np.ndarray | None = None) -> Any:
    saved_n = model.N_
    try:
        model.N_ = -1
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=r"predict\(\) received .* units", category=UserWarning)
            return model.predict(x, z, t=t, x_trt=x_trt, summary_only=False)
    finally:
        model.N_ = saved_n


def _without_matched_intercepts(model: LongBet, values: np.ndarray, n_units: int) -> np.ndarray:
    """No-op because ``_predict_quiet`` already evaluates without unit intercepts."""
    del model, n_units
    return values


def adoption_distribution(
    adoption_model: LongBet, x: np.ndarray, z: np.ndarray, t: np.ndarray, *,
    x_trt: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """First-adoption probabilities and adoption stock under an offer schedule.

    Returns ``(p_adopt, stock)``, both ``(N, T, draws)``: the probability of
    adopting first in each period and of having adopted by each period, from
    the survival product of the population-averaged hazard.
    """
    n_units = z.shape[0]
    if adoption_model.config.random_intercept:
        raise ValueError("the adoption hazard is population-averaged; this model has unit intercepts.")
    eta = np.asarray(_predict_quiet(adoption_model, x, z, t, x_trt).yhats, dtype=np.float64)
    lam = ndtr(eta)
    surv = np.cumprod(1.0 - lam, axis=1)
    ones = np.ones((n_units, 1, eta.shape[-1]))
    p_adopt = lam * np.concatenate([ones, surv[:, :-1, :]], axis=1)
    return p_adopt, 1.0 - surv


def _is_shift_invariant(outcome_model: LongBet, t: np.ndarray, n_periods: int) -> bool:
    """Whether exposure responses at adoption period ``a`` equal period-0 responses shifted by ``a``."""
    if outcome_model.config.split_calendar_trt or outcome_model.config.split_exposure_trt:
        return False
    if n_periods <= 1:
        return True
    s0 = derive_exposure(np.ones((1, n_periods), dtype=np.float32), t)[0]
    for a in range(1, n_periods):
        sched = np.zeros((1, n_periods), dtype=np.float32)
        sched[0, a:] = 1.0
        if not np.array_equal(derive_exposure(sched, t)[0, a:], s0[: n_periods - a]):
            return False
    return True


def offer_effect_on_outcome(
    outcome_model: LongBet, x: np.ndarray, t: np.ndarray,
    p_adopt_offer: np.ndarray, p_adopt_control: np.ndarray, *,
    x_trt: np.ndarray | None = None, fitted_intercepts: np.ndarray | None = None,
) -> np.ndarray:
    """Offer effect on the outcome, ``(N, T, draws)``, by composition.

    For every candidate adoption period ``a`` the outcome equation is evaluated
    under adoption from ``a`` onward; the effect at calendar period ``t`` is
    then the exposure response weighted by the difference in the probability of
    adopting at ``a`` with and without the offer. For a binary outcome the
    response is on the probability scale, with the outcome's unit intercept
    integrated out, or set to the study units' own draws when
    ``fitted_intercepts`` (``(draws, N)``) is given. The integral is exact,
    ``E[Phi(a + sigma Z)] = Phi(a / sqrt(1 + sigma^2))``.
    """
    n_units, n_periods = p_adopt_offer.shape[:2]
    binary = outcome_model.config.outcome == "binary"
    marginal = fitted_intercepts is None
    diff = p_adopt_offer - p_adopt_control
    effect = np.zeros_like(diff)
    if binary:
        sd = (np.sqrt(_flat_draws(outcome_model, "sigma_gamma2"))
              if marginal and outcome_model.config.random_intercept else np.zeros(diff.shape[-1]))
        scale = np.sqrt(1.0 + np.square(sd))[None, None, :]
    if _is_shift_invariant(outcome_model, t, n_periods):
        schedule0 = np.ones((n_units, n_periods), dtype=np.float32)
        pred0 = _predict_quiet(outcome_model, x, schedule0, t, x_trt)
        tau0 = np.asarray(pred0.tauhats, dtype=np.float64)
        if not binary:
            for a in range(n_periods):
                effect[:, a:, :] += diff[:, a, None, :] * tau0[:, : n_periods - a, :]
            return effect
        mu0 = _without_matched_intercepts(outcome_model, np.asarray(pred0.muhats0, dtype=np.float64), n_units)
        if not marginal:
            mu0 = mu0 + np.asarray(fitted_intercepts, dtype=np.float64).T[:, None, :]
        p0_base = ndtr(mu0 / scale)
        for a in range(n_periods):
            tau_a = tau0[:, : n_periods - a, :]
            effect[:, a:, :] += diff[:, a, None, :] * (
                ndtr((mu0[:, a:, :] + tau_a) / scale) - p0_base[:, a:, :]
            )
        return effect
    for a in range(n_periods):
        schedule = np.zeros((n_units, n_periods), dtype=np.float32)
        schedule[:, a:] = 1.0
        pred = _predict_quiet(outcome_model, x, schedule, t, x_trt)
        tau = np.asarray(pred.tauhats, dtype=np.float64)
        weight_a = diff[:, a, None, :]
        if not binary:
            effect += weight_a * tau
            continue
        mu0 = _without_matched_intercepts(outcome_model, np.asarray(pred.muhats0, dtype=np.float64), n_units)
        if not marginal:
            mu0 = mu0 + np.asarray(fitted_intercepts, dtype=np.float64).T[:, None, :]
        effect += weight_a * (ndtr((mu0 + tau) / scale) - ndtr(mu0 / scale))
    return effect


@dataclasses.dataclass(frozen=True)
class EncouragementComparison:
    """Paired whole-unit bootstrap comparison of model and reference ITTs.

    ``replicates`` has axes (replicate, group, horizon, quantity, estimator),
    where quantities are Y,D and estimators are model,reference. Failures remain
    NaN and are listed in ``failures``; incomplete runs never issue flags.
    """

    table: pd.DataFrame
    replicates: np.ndarray
    failures: pd.DataFrame
    metadata: dict[str, Any]

    def cumulative_discrepancy(self) -> pd.DataFrame:
        """Paired bootstrap test of cumulative model vs. reference ITT across horizons.

        Aggregates per-horizon model and reference estimates across all post-offer
        horizons within each group and quantity (``itt_y``, ``itt_d``), returning
        ``cumulative_model``, ``cumulative_reference``, ``cumulative_difference``,
        ``cumulative_se``, ``z_stat``, and ``disagreement``.
        """
        alpha = float(self.metadata.get("alpha", 0.05))
        critical = _critical(alpha)
        replicates_n = int(self.metadata.get("replicates", self.replicates.shape[0]))
        complete = (
            self.failures.empty
            and bool(self.metadata.get("original_diagnostics_passed", False))
            and bool(self.metadata.get("original_reference_inference_available", False))
            and replicates_n >= 100
        )
        groups = list(pd.unique(self.table["group"]))
        rows = []
        for g, group in enumerate(groups):
            for q, quantity in enumerate(("itt_y", "itt_d")):
                sub = self.table.loc[(self.table.group == group) & (self.table.quantity == quantity)]
                cum_model = float(sub["model"].sum())
                cum_ref = float(sub["reference"].sum())
                cum_diff = float(sub["difference"].sum()) if sub["difference"].notna().all() else np.nan
                pairs_cum = self.replicates[:, g, :, q, :].sum(axis=1)  # (replicates, 2)
                if np.isfinite(pairs_cum).all() and pairs_cum.shape[0] >= 2:
                    se = float(np.std(pairs_cum[:, 0] - pairs_cum[:, 1], ddof=1))
                else:
                    se = np.nan
                z_stat = float(cum_diff / se) if (np.isfinite(cum_diff) and np.isfinite(se) and se > 0) else np.nan
                disagree = bool(abs(z_stat) > critical) if (complete and np.isfinite(z_stat)) else pd.NA
                rows.append(dict(
                    group=group,
                    quantity=quantity,
                    cumulative_model=cum_model,
                    cumulative_reference=cum_ref,
                    cumulative_difference=cum_diff,
                    cumulative_se=se,
                    z_stat=z_stat,
                    disagreement=disagree,
                ))
        return pd.DataFrame(rows)


class EncouragementPrediction:
    """Common-population posterior ITTs, ratios, and separate reference inference.

    ``draws`` maps ``itt_y``, ``itt_d``, and ``wald`` to arrays with axes
    (group, horizon, chain, retained draw). No draws are removed or clipped.
    Exactly zero denominators produce NaN and make that ratio's posterior
    interval unavailable. Near-zero and negative denominators remain intact.
    Quantiles of a Wald ratio are not a design-based weak-IV confidence set.
    """

    #: Quantities that are ratios of ITTs: no posterior mean or mean MCSE is reported for them.
    _ratio_quantities: tuple[str, ...] = ("wald",)

    def __init__(self, standardized: EncouragementDraws, reference: pd.DataFrame,
                 metadata: dict[str, Any], *, alpha: float = 0.05):
        _critical(alpha)
        self.standardized = standardized
        self.reference = reference.copy()
        self.metadata = dict(metadata, standardization=standardized.standardization,
                             calibration_status="not_established",
                             inference_version=INFERENCE_VERSION,
                             posterior_inference="model_posterior",
                             reference_inference="normal_ar",
                             target_aligned=standardized.standardization == "conditional",
                             ratio_interpretation="Wald; CACE requires additional assumptions")
        self.alpha = float(alpha)
        dy, dd = standardized.draws["outcome"], standardized.draws["takeup"]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            ratio = np.divide(dy, dd, out=np.full_like(dy, np.nan), where=dd != 0)
        self.draws = {"itt_y": dy, "itt_d": dd, "wald": ratio}
        self.group_labels = standardized.group_labels
        self.periods = standardized.periods
        self.horizons = standardized.horizons

    def stability(self, *, min_ess: float = 400, max_rhat: float = 1.01) -> pd.DataFrame:
        """Diagnose each group/horizon ITT and ratio, including interval endpoints.

        Report bulk/tail ESS, rank R-hat, and quantile MCSE. Mean MCSE is only
        reported for the ITTs: a ratio's posterior mean need not exist. Missing,
        constant, short, or single-chain traces cannot pass the full check.
        Threshold checks never certify interval coverage or IV identification.
        """
        if not np.isfinite(min_ess) or min_ess <= 0:
            raise ValueError("min_ess must be positive and finite.")
        if not np.isfinite(max_rhat) or max_rhat < 1:
            raise ValueError("max_rhat must be finite and at least 1.")
        rows = []
        for quantity, values in self.draws.items():
            for g, group in enumerate(self.group_labels):
                for h, horizon in enumerate(self.horizons):
                    a = values[g, h]
                    row = dict(group=group, period=float(self.periods[h]), horizon=int(horizon),
                               quantity=quantity, ess_bulk=np.nan, ess_tail=np.nan,
                               rhat=np.nan, mcse_mean=np.nan, mcse_median=np.nan,
                               mcse_lower=np.nan, mcse_upper=np.nan)
                    valid = np.all(np.isfinite(a)) and a.shape[1] >= 4 and np.ptp(a) > 0
                    if valid:
                        row["ess_bulk"] = float(compute_ess(a, method="bulk"))
                        row["ess_tail"] = float(compute_ess(
                            a, method="tail", prob=(self.alpha / 2, 1 - self.alpha / 2)))
                        row["rhat"] = float(compute_rhat(a))
                        for field, p in (("median", .5), ("lower", self.alpha / 2),
                                         ("upper", 1 - self.alpha / 2)):
                            row[f"mcse_{field}"] = float(az.mcse(a, method="quantile", prob=p))
                        if quantity not in self._ratio_quantities:
                            row["mcse_mean"] = float(az.mcse(a, method="mean"))
                    row["ess_ok"] = bool(np.isfinite(row["ess_bulk"]) and
                                         np.isfinite(row["ess_tail"]) and
                                         min(row["ess_bulk"], row["ess_tail"]) >= min_ess)
                    row["rhat_ok"] = bool(np.isfinite(row["rhat"]) and row["rhat"] <= max_rhat)
                    row["diagnostics_passed"] = row["ess_ok"] and row["rhat_ok"]
                    row["calibration_status"] = "not_established"
                    rows.append(row)
        return pd.DataFrame(rows)

    def effects(self, *, min_ess: float = 400, max_rhat: float = 1.01) -> pd.DataFrame:
        """Tidy posterior summaries with diagnostics; ratio means are absent.

        Quantiles use the empirical inverse CDF, preserving infinite values.
        Numeric summaries are retained when diagnostics fail, with an explicit
        status. Use ``wald()`` for a table that withholds unconverged ratios.
        """
        rows = []
        for quantity, values in self.draws.items():
            for g, group in enumerate(self.group_labels):
                for h, horizon in enumerate(self.horizons):
                    a = values[g, h].ravel()
                    missing = float(np.mean(np.isnan(a)))
                    q = (np.full(3, np.nan) if missing else
                         np.quantile(a, [self.alpha / 2, .5, 1 - self.alpha / 2],
                                     method="inverted_cdf"))
                    rows.append(dict(group=group, period=float(self.periods[h]), horizon=int(horizon),
                                     quantity=quantity, posterior_lower=q[0], posterior_median=q[1],
                                     posterior_upper=q[2], posterior_mean=(float(a.mean())
                                         if quantity not in self._ratio_quantities else np.nan),
                                     undefined_fraction=missing,
                                     infinite_fraction=float(np.mean(np.isinf(a)))))
        frame = pd.DataFrame(rows).merge(self.stability(min_ess=min_ess, max_rhat=max_rhat),
            on=["group", "period", "horizon", "quantity"], validate="one_to_one")
        frame["posterior_status"] = np.where(frame.undefined_fraction > 0, "undefined_draws",
            np.where(frame.diagnostics_passed, "diagnostics_passed_calibration_unestablished",
                     "diagnostics_failed"))
        support = self.reference[["group", "horizon", "n_units", "n_encouraged", "n_control"]]
        frame = frame.merge(support, on=["group", "horizon"], validate="many_to_one")
        frame["model_extrapolation"] = (frame.n_encouraged == 0) | (frame.n_control == 0)
        return frame

    @property
    def table(self) -> pd.DataFrame:
        return self.effects()

    def itt_covariance(self) -> pd.DataFrame:
        """Empirical joint posterior covariance across groups, horizons and ITTs.

        All entries come from the aligned joint draws. This is a posterior
        covariance under the fitted model; design-based sampling covariance is
        a different quantity. Ratio covariance is omitted because its moments
        need not exist even when its quantiles are well defined.
        """
        values = np.stack([self.draws["itt_y"], self.draws["itt_d"]], axis=2)
        values = values.reshape(-1, np.prod(values.shape[-2:]))
        index = pd.MultiIndex.from_product([self.group_labels, self.horizons, ["itt_y", "itt_d"]],
                                           names=["group", "horizon", "quantity"])
        cov = (np.cov(values, ddof=1) if values.shape[1] >= 2 else
               np.full((len(index), len(index)), np.nan))
        return pd.DataFrame(cov, index=index, columns=index)

    def wald(self, *, require_convergence: bool = True, min_ess: float = 400,
             max_rhat: float = 1.01) -> pd.DataFrame:
        """Posterior ratio medians/intervals, withheld on failed diagnostics by default.

        ``require_convergence=False`` permits inspection of unconverged empirical
        summaries. The raw draws are always available. Passing diagnostics does
        not turn these posterior intervals into calibrated weak-IV inference;
        the design-based confidence sets are in ``reference``.
        """
        all_effects = self.effects(min_ess=min_ess, max_rhat=max_rhat)
        passed = all_effects.groupby(["group", "horizon"]).diagnostics_passed.all()
        frame = all_effects.loc[all_effects.quantity == "wald"].copy().reset_index(drop=True)
        frame["joint_diagnostics_passed"] = [bool(passed.loc[(r.group, r.horizon)])
                                               for r in frame.itertuples()]
        if require_convergence:
            frame.loc[~frame.joint_diagnostics_passed,
                      ["posterior_lower", "posterior_median", "posterior_upper"]] = np.nan
        return frame

    def first_stage(self, *, practical_threshold: float | None = None,
                    min_ess: float = 400, max_rhat: float = 1.01) -> pd.DataFrame:
        """Precision/sign and optional practical relevance, without a strength cutoff.

        ``practical_threshold`` is an analyst-chosen positive adoption risk
        difference; its posterior exceedance probability is separate from
        first-stage precision and MCMC diagnostics.
        """
        if practical_threshold is not None and (
            not np.isscalar(practical_threshold) or not np.isfinite(practical_threshold)
            or practical_threshold <= 0
        ):
            raise ValueError("practical_threshold must be positive and finite.")
        frame = self.effects(min_ess=min_ess, max_rhat=max_rhat).query(
            "quantity == 'itt_d'").copy().reset_index(drop=True)
        a = self.draws["itt_d"]
        frame["posterior_probability_positive"] = (a > 0).mean(axis=(-2, -1)).ravel()
        frame["posterior_interval_includes_zero"] = (
            (frame.posterior_lower <= 0) & (frame.posterior_upper >= 0))
        frame["practical_threshold"] = practical_threshold
        frame["posterior_probability_above_threshold"] = (
            np.nan if practical_threshold is None else
            (a > practical_threshold).mean(axis=(-2, -1)).ravel())
        return frame

    def comparison(self) -> pd.DataFrame:
        """Descriptive model/reference differences on the matching finite-unit target.

        No disagreement flag is fabricated from independent standard errors of
        two estimates from the same data. ``bootstrap_comparison`` on the fitted
        wrapper obtains the required paired sampling covariance by refitting.
        Population-marginal predictions do not match the finite-unit reference.
        """
        rows = []
        aligned = self.metadata["target_aligned"]
        for g, group in enumerate(self.group_labels):
            ref = self.reference.loc[self.reference.group == group].reset_index(drop=True)
            for h, horizon in enumerate(self.horizons):
                for quantity in ("itt_y", "itt_d"):
                    model = float(self.draws[quantity][g, h].mean())
                    reference = float(ref.loc[h, quantity])
                    rows.append(dict(group=group, period=float(self.periods[h]), horizon=int(horizon),
                        quantity=quantity, model=model, reference=reference,
                        difference=model - reference if aligned else np.nan,
                        difference_se=np.nan, disagreement=pd.NA,
                        comparison_status="paired_sampling_covariance_required" if aligned
                                          else "different_standardization_targets"))
        return pd.DataFrame(rows)

    def cumulative_discrepancy(self) -> pd.DataFrame:
        """Compare cumulative model and design-based reference ITTs across all post-offer horizons."""
        rows = []
        aligned = self.metadata["target_aligned"]
        for g, group in enumerate(self.group_labels):
            ref = self.reference.loc[self.reference.group == group].reset_index(drop=True)
            for quantity in ("itt_y", "itt_d"):
                cum_draws = self.draws[quantity][g].sum(axis=0).ravel()
                cum_model = float(cum_draws.mean())
                cum_model_sd = float(cum_draws.std(ddof=1)) if cum_draws.size >= 2 else np.nan
                cum_ref = float(ref[quantity].sum())
                if f"{quantity}_cum_se" in ref.columns and len(ref) > 0:
                    cum_ref_se = float(ref[f"{quantity}_cum_se"].iloc[0])
                else:
                    ref_se_sq = float((ref[f"{quantity}_se"] ** 2).sum())
                    cum_ref_se = float(np.sqrt(ref_se_sq)) if np.isfinite(ref_se_sq) else np.nan
                cum_diff = cum_model - cum_ref if aligned else np.nan
                cum_se = cum_ref_se if aligned else np.nan
                z_stat = (float(cum_diff / cum_se)
                          if (aligned and np.isfinite(cum_diff) and np.isfinite(cum_se) and cum_se > 0)
                          else np.nan)
                rows.append(dict(
                    group=group,
                    quantity=quantity,
                    n_periods=int(len(self.horizons)),
                    cumulative_model=cum_model,
                    cumulative_model_sd=cum_model_sd,
                    cumulative_reference=cum_ref,
                    cumulative_reference_se=cum_ref_se,
                    cumulative_difference=cum_diff,
                    cumulative_se=cum_se,
                    z_stat=z_stat,
                ))
        return pd.DataFrame(rows)


@dataclasses.dataclass(frozen=True)
class ConditionalEffect:
    """Conditional treatment effect summary across units and post-treatment horizons.

    Attributes
    ----------
    draws : np.ndarray
        Array of shape (N, H, D) containing MCMC draws.
    mean : np.ndarray
        Posterior mean of shape (N, H).
    median : np.ndarray
        Posterior median of shape (N, H).
    lower : np.ndarray
        Posterior lower quantile (level alpha/2) of shape (N, H).
    upper : np.ndarray
        Posterior upper quantile (level 1 - alpha/2) of shape (N, H).
    sd : np.ndarray
        Posterior standard deviation of shape (N, H).
    """

    draws: np.ndarray
    mean: np.ndarray
    median: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    sd: np.ndarray


class ConditionalEncouragementPrediction:
    """Unit-level conditional ITTs, CACE, and risk-aware policy targeting.

    Attributes
    ----------
    citt_y : ConditionalEffect
        Conditional Intent-to-Treat on the outcome (N, H, D).
    citt_d : ConditionalEffect
        Conditional Intent-to-Treat on adoption / compliance (N, H, D).
    cace : ConditionalEffect
        Conditional Complier Average Causal Effect (N, H, D).
    periods : np.ndarray
        Calendar periods corresponding to post-treatment horizons (H,).
    horizons : np.ndarray
        Exposure duration indices (H,).
    alpha : float
        Significance level for credible intervals.
    """

    def __init__(
        self,
        citt_y_draws: np.ndarray,
        citt_d_draws: np.ndarray,
        periods: np.ndarray,
        horizons: np.ndarray,
        *,
        p0_draws: np.ndarray | None = None,
        alpha: float = 0.05,
        cace_stabilization: float = 0.02,
        monotonic_first_stage: bool = True,
    ):
        self.alpha = float(alpha)
        self.periods = np.asarray(periods)
        self.horizons = np.asarray(horizons)
        if not (np.isfinite(cace_stabilization) and cace_stabilization > 0):
            raise ValueError("cace_stabilization must be positive and finite.")
        self.cace_stabilization = float(cace_stabilization)
        self.monotonic_first_stage = bool(monotonic_first_stage)

        # 1. CITT_Y
        self.citt_y = self._build_effect(citt_y_draws, alpha)

        # 2. CITT_D with optional monotonicity enforcement
        if self.monotonic_first_stage:
            # Encouragement can never discourage adoption (eliminating defiers):
            citt_d_draws = np.maximum(0.0, citt_d_draws)
        self.citt_d = self._build_effect(citt_d_draws, alpha)

        # 2b. Principal Strata (Compliers, Always-Takers, Never-Takers) under Monotonicity
        self.p0_draws = np.asarray(p0_draws) if p0_draws is not None else None
        self.compliers = self.citt_d
        if self.p0_draws is not None:
            p_at = np.clip(self.p0_draws, 0.0, 1.0)
            p_c = np.clip(self.citt_d.draws, 0.0, 1.0)
            p_nt = np.clip(1.0 - (p_at + p_c), 0.0, 1.0)
            self.always_takers = self._build_effect(p_at, alpha)
            self.never_takers = self._build_effect(p_nt, alpha)
        else:
            self.always_takers = None
            self.never_takers = None

        # 3. CACE: the ratio of the two offer effects, draw by draw, with a
        # small floor on the adoption contrast so a near-zero first stage
        # cannot produce arbitrarily large ratios.
        cace_draws_eff = citt_y_draws / (self.citt_d.draws + self.cace_stabilization)
        self.cace = self._build_effect(cace_draws_eff, alpha)

    @staticmethod
    def _build_effect(draws: np.ndarray, alpha: float) -> ConditionalEffect:
        mean = np.mean(draws, axis=-1)
        median = np.median(draws, axis=-1)
        lower = np.quantile(draws, alpha / 2, axis=-1)
        upper = np.quantile(draws, 1.0 - alpha / 2, axis=-1)
        sd = np.std(draws, axis=-1, ddof=1) if draws.shape[-1] > 1 else np.zeros_like(mean)
        return ConditionalEffect(
            draws=draws, mean=mean, median=median, lower=lower, upper=upper, sd=sd
        )

    def _resolve_horizon_idx(self, horizon: int) -> int:
        h_int = int(horizon)
        matches = np.flatnonzero(np.asarray(self.horizons, dtype=int) == h_int)
        if len(matches) > 0:
            return int(matches[0])
        if 0 <= h_int < len(self.horizons):
            return h_int
        raise ValueError(
            f"horizon={horizon} not found in horizons {list(self.horizons)} "
            f"or valid 0-based indices 0..{len(self.horizons) - 1}."
        )

    def principal_strata(
        self,
        *,
        horizon: int | None = None,
        as_counts: bool = True,
    ) -> pd.DataFrame:
        """Posterior distribution of principal strata proportions and cohort counts.

        Under the monotonicity assumption (no defiers: D(1) >= D(0)), each unit at each
        post-treatment horizon belongs to one of three latent principal strata:
        - Compliers: adopt if and only if encouraged (D(1) = 1, D(0) = 0)
        - Always-Takers: adopt even in the absence of encouragement (D(1) = 1, D(0) = 1)
        - Never-Takers: refuse adoption even when encouraged (D(1) = 0, D(0) = 0)

        Parameters
        ----------
        horizon : int, optional
            Specific post-treatment exposure horizon (matching ``self.horizons``, e.g. 1..H,
            or 0 for the first post-treatment horizon). If None, evaluates all horizons.
        as_counts : bool, default True
            If True, reports both population proportions and expected cohort counts (N * prob).

        Returns
        -------
        pd.DataFrame
            DataFrame summarizing posterior mean, median, lower/upper credible bounds,
            and expected counts for each principal stratum.
        """
        if self.always_takers is None or self.never_takers is None:
            strata_dict = {"complier": self.compliers.draws}
        else:
            strata_dict = {
                "complier": self.compliers.draws,
                "always_taker": self.always_takers.draws,
                "never_taker": self.never_takers.draws,
            }

        n_units = self.citt_d.draws.shape[0]
        h_indices = (
            range(len(self.horizons))
            if horizon is None
            else [self._resolve_horizon_idx(horizon)]
        )

        rows = []
        for h_idx in h_indices:
            p_val = float(self.periods[h_idx])
            h_val = int(self.horizons[h_idx])

            for name, draws in strata_dict.items():
                unit_draws = draws[:, h_idx, :]  # (N, D)
                pop_draws = np.mean(unit_draws, axis=0)  # (D,)
                count_draws = np.sum(unit_draws, axis=0)  # (D,)

                row = {
                    "period": p_val,
                    "horizon": h_val,
                    "stratum": name,
                    "prob_mean": float(np.mean(pop_draws)),
                    "prob_median": float(np.median(pop_draws)),
                    "prob_lower": float(np.quantile(pop_draws, self.alpha / 2)),
                    "prob_upper": float(np.quantile(pop_draws, 1.0 - self.alpha / 2)),
                    "prob_sd": float(np.std(pop_draws, ddof=1)) if len(pop_draws) > 1 else 0.0,
                }
                if as_counts:
                    row.update({
                        "count_mean": float(np.mean(count_draws)),
                        "count_median": float(np.median(count_draws)),
                        "count_lower": float(np.quantile(count_draws, self.alpha / 2)),
                        "count_upper": float(np.quantile(count_draws, 1.0 - self.alpha / 2)),
                        "n_total": n_units,
                    })
                rows.append(row)

        return pd.DataFrame(rows)

    def strata_summary(self, horizon: int | None = None) -> str:
        """Human-readable statement of principal strata estimates matching BIVA conventions."""
        h_idx = 0 if horizon is None else self._resolve_horizon_idx(horizon)
        df = self.principal_strata(horizon=int(self.horizons[h_idx]), as_counts=True)
        p_val = self.periods[h_idx]
        h_val = self.horizons[h_idx]
        n_units = self.citt_d.draws.shape[0]

        lines = [
            f"Principal Strata Estimates at Horizon {h_val} (Period {p_val:.2f}, N = {n_units}):"
        ]
        for _, r in df.iterrows():
            st = r["stratum"].replace("_", "-").title()
            pct_mean = f"{r['prob_mean'] * 100:.1f}%"
            pct_ci = f"[{r['prob_lower'] * 100:.1f}%, {r['prob_upper'] * 100:.1f}%]"
            cnt_mean = f"{r['count_mean']:.1f}"
            cnt_ci = f"[{r['count_lower']:.1f}, {r['count_upper']:.1f}]"
            lines.append(f"  • {st:<13}: {pct_mean:>6} (95% CI: {pct_ci})  —  ~{cnt_mean} units (95% CI: {cnt_ci})")

        c_rows = df[df["stratum"] == "complier"]
        if not c_rows.empty:
            c_mean = c_rows.iloc[0]["prob_mean"]
            if c_mean < 0.10:
                lines.append("  [WARNING: Weak Instrument detected - complier share < 10%]")
            else:
                lines.append("  [Robust Instrument - complier share exceeds 10% threshold]")

        return "\n".join(lines)

    def cumulative_lift(self, start_h: int = 0, end_h: int | None = None) -> np.ndarray:
        """Draws of cumulative outcome lift: sum_{t} CITT_Y(X_i, t). Shape: (N, D)."""
        sl = slice(start_h, end_h)
        return np.sum(self.citt_y.draws[:, sl, :], axis=1)

    def cumulative_lift_summary(self, start_h: int = 0, end_h: int | None = None) -> pd.DataFrame:
        """Summary table of cumulative outcome lift across units."""
        cum = self.cumulative_lift(start_h, end_h)
        return pd.DataFrame({
            "unit": np.arange(cum.shape[0]),
            "cum_mean": np.mean(cum, axis=-1),
            "cum_median": np.median(cum, axis=-1),
            "cum_lower": np.quantile(cum, self.alpha / 2, axis=-1),
            "cum_upper": np.quantile(cum, 1.0 - self.alpha / 2, axis=-1),
            "cum_sd": np.std(cum, axis=-1, ddof=1) if cum.shape[-1] > 1 else 0.0,
        })

    def breakeven_probability(
        self, cost: float, start_h: int = 0, end_h: int | None = None
    ) -> np.ndarray:
        """Posterior probability that cumulative lift exceeds cost: Pr(sum CITT_Y > cost | Data)."""
        cum = self.cumulative_lift(start_h, end_h)
        return np.mean(cum > float(cost), axis=-1)

    def optimal_policy(
        self,
        cost: float,
        hurdle: float = 0.50,
        start_h: int = 0,
        end_h: int | None = None,
        *,
        require_positive_expected_net: bool = False,
    ) -> np.ndarray:
        """Risk-aware encouragement rule: 1{Pr(sum CITT_Y > cost | Data) >= hurdle}."""
        p = self.breakeven_probability(cost, start_h, end_h)
        eligible = p >= float(hurdle)
        if require_positive_expected_net:
            cum_draws = self.cumulative_lift(start_h, end_h)
            expected_net = np.mean(cum_draws, axis=-1) - float(cost)
            eligible = eligible & (expected_net > 0)
        return eligible

    def knapsack_policy(
        self,
        cost: float,
        *,
        budget: float | None = None,
        capacity: int | None = None,
        hurdle: float = 0.50,
        ranking_metric: str = "expected_net_value",
        start_h: int = 0,
        end_h: int | None = None,
        require_positive_expected_net: bool = False,
    ) -> np.ndarray:
        """Compute optimal targeted policy under budget and/or unit capacity constraints.

        Parameters
        ----------
        cost : float
            Cost per encouraged account.
        budget : float, optional
            Maximum allowable total expenditure (cost * sum(decision) <= budget).
        capacity : int, optional
            Maximum allowable count of encouraged accounts (sum(decision) <= capacity).
        hurdle : float
            Minimum posterior breakeven probability required for eligibility.
        ranking_metric : str
            Priority criterion when constraints bind:
            - 'expected_net_value': E[Lift_i - cost]
            - 'certainty_adjusted': E[Lift_i - cost] / (SD(Lift_i) + 1e-6)
            - 'breakeven_probability': Pr(Lift_i > cost)
        start_h, end_h : int, optional
            Horizon slice for cumulative lift.
        require_positive_expected_net : bool, default False
            If True, also requires E[Lift_i - cost] > 0 in addition to Pr(Lift_i > cost) >= hurdle.

        Returns
        -------
        np.ndarray of shape (N,)
            Boolean indicator of encouraged units satisfying all constraints.
        """
        cost = float(cost)
        if cost <= 0:
            raise ValueError("cost must be positive.")
        if ranking_metric not in ("expected_net_value", "certainty_adjusted", "breakeven_probability"):
            raise ValueError(
                f"Unknown ranking_metric '{ranking_metric}'. Available: "
                "'expected_net_value', 'certainty_adjusted', 'breakeven_probability'."
            )
        n_units = self.citt_y.draws.shape[0]

        # 1. Eligibility check (consistent with optimal_policy when unconstrained)
        p_breakeven = self.breakeven_probability(cost, start_h, end_h)
        cum_draws = self.cumulative_lift(start_h, end_h)
        expected_lift = np.mean(cum_draws, axis=-1)
        expected_net = expected_lift - cost
        sd_lift = np.std(cum_draws, axis=-1, ddof=1) if cum_draws.shape[-1] > 1 else np.ones(n_units)

        eligible = p_breakeven >= float(hurdle)
        if require_positive_expected_net:
            eligible = eligible & (expected_net > 0)
        eligible_indices = np.flatnonzero(eligible)

        # 2. Determine max allowed units under constraints
        max_units = n_units
        if capacity is not None:
            max_units = min(max_units, int(capacity))
        if budget is not None:
            max_units_by_budget = int(np.floor(float(budget) / cost))
            max_units = min(max_units, max_units_by_budget)

        if max_units <= 0:
            return np.zeros(n_units, dtype=bool)

        if len(eligible_indices) <= max_units:
            decision = np.zeros(n_units, dtype=bool)
            decision[eligible_indices] = True
            return decision

        # 3. Ranking metric calculation
        if ranking_metric == "expected_net_value":
            scores = expected_net
        elif ranking_metric == "certainty_adjusted":
            scores = expected_net / (sd_lift + 1e-6)
        elif ranking_metric == "breakeven_probability":
            scores = p_breakeven
        else:
            raise ValueError(
                f"Unknown ranking_metric '{ranking_metric}'. Available: "
                "'expected_net_value', 'certainty_adjusted', 'breakeven_probability'."
            )

        # 4. Rank eligible units and select top max_units
        eligible_scores = scores[eligible_indices]
        top_order = np.argsort(-eligible_scores)[:max_units]
        selected_indices = eligible_indices[top_order]

        decision = np.zeros(n_units, dtype=bool)
        decision[selected_indices] = True
        return decision

    def policy_value(
        self,
        cost: float,
        hurdle: float = 0.50,
        true_cumulative_lift: np.ndarray | None = None,
        budget: float | None = None,
        capacity: int | None = None,
        ranking_metric: str = "expected_net_value",
        start_h: int = 0,
        end_h: int | None = None,
    ) -> float:
        """Net policy value per account: mean_i 1{g(X_i)=1} * (Lift_i - cost)."""
        if budget is not None or capacity is not None:
            dec = self.knapsack_policy(
                cost,
                budget=budget,
                capacity=capacity,
                hurdle=hurdle,
                ranking_metric=ranking_metric,
                start_h=start_h,
                end_h=end_h,
            )
        else:
            dec = self.optimal_policy(cost, hurdle, start_h, end_h)
        lift = (
            np.mean(self.cumulative_lift(start_h, end_h), axis=-1)
            if true_cumulative_lift is None
            else np.asarray(true_cumulative_lift)
        )
        return float(np.mean(np.where(dec, lift - float(cost), 0.0)))


# ---------------------------------------------------------------------------
# Continuous intensity: posterior containers
# ---------------------------------------------------------------------------

class IntensityEncouragementPrediction(EncouragementPrediction):
    """Population offer effects of the continuous-intensity model, with the design-based T0 reference.

    Beside ``itt_y``, ``itt_d`` and the per-period ratio ``wald`` (which
    divides by a first stage that may have returned to zero), ``draws`` carries
    ``omega_cum``, the cumulative return per unit of cumulative dose
    ``sum_{r<=t} ITT_Y(r) / sum_{r<=t} ITT_D(r)``. ``kernel_draws[group]`` holds
    the population lag kernel ``kappa_l * mean_i nu_Y(x_i)`` and
    ``kernel_acr_draws[group]`` the compliance-weighted kernel the design-based
    deconvolution targets, both with axes ``(lag, chain, draw)``;
    ``references[group]`` is the matching :class:`IntensityReference`.
    """

    _ratio_quantities = ("wald", "omega_cum")

    def __init__(self, standardized: EncouragementDraws, reference: pd.DataFrame,
                 metadata: dict[str, Any], *, kernel_draws: dict[str, np.ndarray],
                 kernel_acr_draws: dict[str, np.ndarray], references: dict[str, IntensityReference | None],
                 alpha: float = 0.05):
        super().__init__(standardized, reference, metadata, alpha=alpha)
        dy, dd = standardized.draws["outcome"], standardized.draws["takeup"]
        self.draws["omega_cum"] = cumulative_return_draws(dy, dd)
        self.kernel_draws = kernel_draws
        self.kernel_acr_draws = kernel_acr_draws
        self.references = references
        self.metadata.update(
            inference_version=INFERENCE_VERSION_INTENSITY, engine=ENGINE_INTENSITY,
            ratio_interpretation="omega_cum: cumulative return per unit of offer-induced cumulative dose; "
                                 "wald: per-period ratio (unstable when the first stage fades)")

    def kernel(self) -> pd.DataFrame:
        """Posterior lag kernel per group beside the design-based Toeplitz kernel.

        ``kappa`` is the population-average kernel, ``kappa_acr`` the
        compliance-weighted kernel that the model-free T0 deconvolution
        (``reference_kappa``) identifies; compare like with like.
        """
        rows = []
        for group in self.group_labels:
            ref = self.references.get(group)
            for name, store in (("kappa", self.kernel_draws), ("kappa_acr", self.kernel_acr_draws)):
                draws = store[group]
                cum = np.cumsum(draws, axis=0)
                for lag in range(draws.shape[0]):
                    a, c = draws[lag].ravel(), cum[lag].ravel()
                    row = dict(group=group, quantity=name, lag=lag,
                               posterior_mean=float(a.mean()), posterior_sd=float(a.std(ddof=1)) if a.size > 1 else 0.0,
                               posterior_lower=float(np.quantile(a, self.alpha / 2)),
                               posterior_median=float(np.median(a)),
                               posterior_upper=float(np.quantile(a, 1 - self.alpha / 2)),
                               cumulative_mean=float(c.mean()),
                               cumulative_lower=float(np.quantile(c, self.alpha / 2)),
                               cumulative_upper=float(np.quantile(c, 1 - self.alpha / 2)))
                    if ref is not None and lag < len(ref.kernel):
                        k = ref.kernel.iloc[lag]
                        row.update(reference_kappa=float(k.kappa), reference_kappa_se=float(k.kappa_se),
                                   reference_kappa_cum=float(k.kappa_cum),
                                   reference_kappa_cum_se=float(k.kappa_cum_se))
                    else:
                        row.update(reference_kappa=np.nan, reference_kappa_se=np.nan,
                                   reference_kappa_cum=np.nan, reference_kappa_cum_se=np.nan)
                    rows.append(row)
        return pd.DataFrame(rows)

    def discrepancy(self) -> pd.DataFrame:
        """``z_{tier - T0}`` per group and horizon: posterior cumulative compliance-weighted kernel against the T0 reference."""
        label = TIER_LABELS.get(self.metadata.get("identification", ""), "model")
        frames = []
        for group in self.group_labels:
            ref = self.references.get(group)
            if ref is None:
                continue
            draws = self.kernel_acr_draws[group]
            flat = draws.reshape(draws.shape[0], -1).T
            frame = ref.discrepancy(flat, label=label)
            frame.insert(0, "group", group)
            frames.append(frame)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def cumulative_return(self, *, min_ess: float = 400, max_rhat: float = 1.01) -> pd.DataFrame:
        """Posterior ``omega_cum`` by horizon beside the reference's exact Fieller/AR confidence set."""
        frame = self.effects(min_ess=min_ess, max_rhat=max_rhat)
        frame = frame.loc[frame.quantity == "omega_cum"].reset_index(drop=True)
        columns = [c for c in self.reference.columns if c.startswith("omega_cum")]
        return frame.merge(self.reference[["group", "horizon"] + columns], on=["group", "horizon"],
                           validate="one_to_one")

    def diagnostics(self) -> pd.DataFrame:
        """Design-based diagnostics per group: Toeplitz conditioning, dip screen, FOSD check and the headline discrepancy."""
        rows = []
        disc = self.discrepancy()
        for group in self.group_labels:
            ref = self.references.get(group)
            if ref is None:
                rows.append(dict(group=group, condition_number=np.nan, z_dip=np.nan, dip_pvalue=np.nan,
                                 fosd_pvalue=np.nan, fosd_sup_violation=np.nan, z_discrepancy=np.nan))
                continue
            z_disc = disc.loc[disc.group == group, "z_stat"]
            rows.append(dict(group=group, condition_number=ref.condition_number,
                             z_dip=ref.dip["z_stat"], dip_pvalue=ref.dip["pvalue"],
                             fosd_pvalue=ref.fosd["pvalue"], fosd_sup_violation=ref.fosd["sup_violation"],
                             z_discrepancy=float(z_disc.iloc[-1]) if len(z_disc) else np.nan))
        return pd.DataFrame(rows)


class IntensityConditionalPrediction(ConditionalEncouragementPrediction):
    """Unit-level offer effects of the continuous-intensity model and the decision layer.

    ``citt_y`` and ``citt_d`` are the conditional offer effects on the outcome
    and on the dose (``(N, H, draws)``); ``kernel`` is the conditional lag
    kernel ``kappa(x, l) = nu_Y(x) kappa_l`` (``(N, L+1, draws)``),
    ``kernel_cum`` its running sum, and ``omega_cum`` the conditional cumulative
    return per unit of cumulative dose (exact zero denominators give NaN; no
    stabilization floor). ``cace`` is an alias of ``omega_cum``. The lift,
    break-even, optimal and knapsack policies are inherited unchanged.
    Principal strata are undefined for a continuous dose; use the arm-wise ECDF
    shift in the design-based reference instead.
    """

    def __init__(self, citt_y_draws: np.ndarray, citt_d_draws: np.ndarray, kernel_draws: np.ndarray,
                 periods: np.ndarray, horizons: np.ndarray, *, alpha: float = 0.05):
        self.alpha = float(alpha)
        self.periods = np.asarray(periods)
        self.horizons = np.asarray(horizons)
        self.cace_stabilization = 0.0
        self.monotonic_first_stage = False
        self.citt_y = self._build_effect(np.asarray(citt_y_draws, dtype=np.float64), self.alpha)
        self.citt_d = self._build_effect(np.asarray(citt_d_draws, dtype=np.float64), self.alpha)
        self.kernel = self._build_effect(np.asarray(kernel_draws, dtype=np.float64), self.alpha)
        self.kernel_cum = self._build_effect(np.cumsum(np.asarray(kernel_draws, dtype=np.float64), axis=1),
                                             self.alpha)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            omega = cumulative_return_draws(self.citt_y.draws, self.citt_d.draws)
            self.omega_cum = ConditionalEffect(
                draws=omega, mean=np.nanmean(omega, axis=-1), median=np.nanmedian(omega, axis=-1),
                lower=np.nanquantile(omega, self.alpha / 2, axis=-1),
                upper=np.nanquantile(omega, 1 - self.alpha / 2, axis=-1),
                sd=np.nanstd(omega, axis=-1, ddof=1) if omega.shape[-1] > 1 else np.zeros(omega.shape[:-1]))
        self.cace = self.omega_cum
        self.compliers = self.citt_d
        self.p0_draws = None
        self.always_takers = None
        self.never_takers = None

    def principal_strata(self, *, horizon: int | None = None, as_counts: bool = True) -> pd.DataFrame:
        raise ValueError("principal strata are not defined for a continuous dose; see the arm-wise "
                         "ECDF shift of cumulative dose in the design-based reference.")

    def strata_summary(self, horizon: int | None = None) -> str:
        raise ValueError("principal strata are not defined for a continuous dose; see the arm-wise "
                         "ECDF shift of cumulative dose in the design-based reference.")

    def kernel_summary(self) -> pd.DataFrame:
        """Posterior mean and interval of the conditional kernel averaged over the evaluated units."""
        draws = self.kernel.draws.mean(axis=0)          # (L+1, draws)
        cum = np.cumsum(draws, axis=0)
        return pd.DataFrame({
            "lag": np.arange(draws.shape[0]),
            "kappa_mean": draws.mean(axis=-1), "kappa_lower": np.quantile(draws, self.alpha / 2, axis=-1),
            "kappa_upper": np.quantile(draws, 1 - self.alpha / 2, axis=-1),
            "kappa_cum_mean": cum.mean(axis=-1), "kappa_cum_lower": np.quantile(cum, self.alpha / 2, axis=-1),
            "kappa_cum_upper": np.quantile(cum, 1 - self.alpha / 2, axis=-1),
        })


class LongBetEncourage:
    """The encouragement model: adoption clock (binary) or distributed lags (continuous intensity).

    ``treatment="absorbing"`` (default) is the adoption-clock model of the
    module docstring: ``config`` (or keyword overrides) applies to both
    equations; the adoption equation is always fitted as a binary hazard.
    ``outcome`` is ``'continuous'`` or ``'binary'``. Fit a complete single-wave
    individually randomized panel with absorbing offer and adoption indicators.

    ``treatment="intensity"`` fits a continuous, non-absorbing encouraged
    behaviour ``d >= 0`` (hours, training-months) with two stages: a continuous
    LongBet for the dose on the offer clock (unit intercept, AR(1) errors,
    PRNG stream ``key_d``) and the distributed-lag kernel for the outcome
    ``nu_Y(x) sum_l kappa_l D_{i,t-l}`` with a GP prior over lags (stream
    ``key_y``). ``identification`` chooses what Stage 2 learns the kernel from:

    * ``"deconvolution"`` (T0) -- the Stage-1 fitted offer-induced dose
      ``CITT_D(X_i, t-l) Z_{i,t-l}``: no within-arm dose variation, immune to
      dose endogeneity and measurement error, highest variance.
    * ``"control_function"`` (T1, default) -- the observed dose plus the 3-term
      plug-in control function ``[v_it, sum_{l>=1} v_{i,t-l}, vbar_i] (x) [1, x]``
      built from the Stage-1 posterior-mean residuals ``v_it``.
    * ``"unconfounded"`` (T2) -- the observed dose alone.

    ``lag_horizon`` defaults to ``T - t0`` (one lag per post-offer period) and
    ``demean_pre`` subtracts each unit's pre-offer mean dose before building
    lags. ``stage2_backend`` is ``"longbet"`` (the sampler) or ``"linear"`` (a
    closed-form two-way fixed-effects regression without covariate moderation,
    for smoke tests and fallbacks). The offer effects follow by composition,
    draw by draw: ``CITT_Y(x,t) = nu_Y(x) sum_l kappa_l CITT_D(x,t-l)``.

    The model's components are fixed (:data:`MODEL_COMPONENTS`): a horseshoe
    linear term in both moderators, a learned trajectory scale, interweaving
    of the intercept variance, AR(1) errors for a continuous outcome, and a
    population-averaged adoption hazard. The other options (sampler budget,
    forests, priors) are free; ``n_skip`` left unset thins by 2.
    """

    def __init__(self, config: LongBetConfig | None = None, *, outcome: str = "continuous",
                 treatment: str = "absorbing", identification: str | None = None,
                 lag_horizon: int | None = None, demean_pre: bool = False,
                 stage2_backend: str = "longbet", **config_kwargs: Any):
        if outcome not in ("continuous", "binary"):
            raise ValueError("outcome must be 'continuous' or 'binary'.")
        if treatment not in TREATMENTS:
            raise ValueError("treatment must be 'absorbing' or 'intensity'.")
        if treatment == "absorbing":
            if identification is not None:
                raise ValueError("identification applies to treatment='intensity' only.")
            if lag_horizon is not None or demean_pre:
                raise ValueError("lag_horizon and demean_pre apply to treatment='intensity' only.")
        else:
            if outcome != "continuous":
                raise ValueError("treatment='intensity' supports a continuous outcome only.")
            identification = DEFAULT_IDENTIFICATION if identification is None else identification
            if identification not in IDENTIFICATION_TIERS:
                raise ValueError(f"identification must be one of {IDENTIFICATION_TIERS}.")
            if lag_horizon is not None and (isinstance(lag_horizon, bool)
                                            or not isinstance(lag_horizon, (int, np.integer))
                                            or lag_horizon < 0):
                raise ValueError("lag_horizon must be a nonnegative integer or None.")
        if stage2_backend not in STAGE2_BACKENDS:
            raise ValueError(f"stage2_backend must be one of {STAGE2_BACKENDS}.")
        if config is None:
            config = LongBetConfig(**config_kwargs)
        elif isinstance(config, LongBetConfig):
            config = dataclasses.replace(config, **config_kwargs) if config_kwargs else config
        else:
            raise TypeError("config must be a LongBetConfig or None.")
        if config.outcome != "continuous" or config.num_categories is not None:
            raise ValueError("pass the outcome type as outcome=..., not in the config.")
        self.config = self._model_config(config, outcome)
        self.outcome = outcome
        self.treatment = treatment
        self.identification = identification
        self.lag_horizon = None if lag_horizon is None else int(lag_horizon)
        self.demean_pre = bool(demean_pre)
        self.stage2_backend = stage2_backend
        self.outcome_model: LongBet | None = None
        self.adoption_model: LongBet | None = None
        #: Stage 1 of the intensity model (the continuous dose equation).
        self.dose_model: LongBet | None = None
        #: Stage 2 of the intensity model behind the backend protocol.
        self.stage2: Any = None
        self._data: dict[str, np.ndarray] = {}
        self.metadata: dict[str, Any] = {}

    @staticmethod
    def _model_config(config: LongBetConfig, outcome: str) -> LongBetConfig:
        """``config`` with the model's components filled in; refuses a config that
        switches one off or rules one out."""
        for name, what in MODEL_COMPONENTS.items():
            if getattr(config, name) is False:
                raise ValueError(f"LongBetEncourage always uses {what}; leave {name} unset.")
        if not config.random_intercept:
            raise ValueError("LongBetEncourage has unit intercepts in the outcome equation; "
                             "leave random_intercept=True.")
        ar1 = outcome == "continuous"
        if config.ar1_errors is not None and bool(config.ar1_errors) != ar1:
            raise ValueError("LongBetEncourage models AR(1) errors for a continuous outcome and "
                             "not for a binary one; leave ar1_errors unset.")
        return resolve_model_defaults(config, {**{name: True for name in MODEL_COMPONENTS},
                                               "ar1_errors": ar1, "n_skip": DEFAULT_N_SKIP})

    @classmethod
    def model_config(cls, outcome: str = "continuous", **overrides: Any) -> LongBetConfig:
        """The configuration a fit with these ``overrides`` would use."""
        return cls._model_config(LongBetConfig(**overrides), outcome)

    # -- inputs ---------------------------------------------------------------
    def _inputs(self, y: Any, d: Any, z: Any, x: Any, t: Any, x_trt: Any) -> dict:
        panel = _validate(z, d, t)
        y = _finite_matrix(y, "y", len(panel.z), panel.z.shape)
        if self.outcome == "binary" and not np.isin(y, [0, 1]).all():
            raise ValueError("binary outcome y must contain only 0 and 1.")
        data = dict(y=y, d=panel.d.copy(), z=panel.z.copy(), t=panel.t.copy(),
                    x=_finite_matrix(x, "x", len(panel.z)))
        if x_trt is not None:
            data["x_trt"] = _finite_matrix(x_trt, "x_trt", len(panel.z))
        model_t = data["t"].astype(np.float32)
        for name in ("z", "d"):
            if not np.all(np.diff(model_t) > 0) or not np.array_equal(
                    derive_exposure(data[name], model_t), derive_exposure(data[name], data["t"])):
                raise ValueError("t loses clock precision in the float32 model; shift the "
                                 "calendar origin before fitting.")
        return data

    @staticmethod
    def _first_adoption_events(d: np.ndarray) -> np.ndarray:
        """First-adoption indicator on the at-risk set; NaN once a unit has adopted."""
        previous = np.concatenate([np.zeros((len(d), 1)), d[:, :-1]], axis=1)
        return np.where(previous == 1, np.nan, d)

    def _equation_configs(self) -> tuple[LongBetConfig, LongBetConfig]:
        outcome_cfg = dataclasses.replace(
            self.config, outcome=self.outcome, num_categories=None,
            anchor_moderator=(self.outcome == "continuous"))
        # The adoption hazard is a binary, population-averaged equation: no
        # AR(1) innovations and no unit frailty.
        adoption_cfg = dataclasses.replace(
            self.config, outcome="binary", num_categories=None, ar1_errors=False,
            random_intercept=False, anchor_moderator=False)
        return outcome_cfg, adoption_cfg

    # -- continuous intensity: inputs and configs ------------------------------
    def _inputs_intensity(self, y: Any, d: Any, z: Any, x: Any, t: Any, x_trt: Any,
                          strata: Any) -> dict:
        panel = validate_intensity(z, d, t)
        n = len(panel.z)
        data = dict(y=_finite_matrix(y, "y", n, panel.z.shape), d=panel.d.copy(), z=panel.z.copy(),
                    t=panel.t.copy(), x=_finite_matrix(x, "x", n))
        if x_trt is not None:
            data["x_trt"] = _finite_matrix(x_trt, "x_trt", n)
        if strata is not None:
            s = np.asarray(strata)
            if s.shape != (n,) or np.any(pd.isna(s)):
                raise ValueError(f"strata must be a complete length-{n} vector of baseline labels.")
            data["strata"] = s.astype(str)
        model_t = data["t"].astype(np.float32)
        if not np.all(np.diff(model_t) > 0) or not np.array_equal(
                derive_exposure(data["z"], model_t), derive_exposure(data["z"], data["t"])):
            raise ValueError("t loses clock precision in the float32 model; shift the "
                             "calendar origin before fitting.")
        if self.lag_horizon is not None and self.lag_horizon > default_lag_horizon(panel.z.shape[1], panel.start):
            raise ValueError(f"lag_horizon cannot exceed T - t0 = "
                             f"{default_lag_horizon(panel.z.shape[1], panel.start)} for this panel.")
        return data

    def _stage_configs(self) -> tuple[LongBetConfig, LongBetConfig]:
        """Stage 1 (dose) and Stage 2 (outcome) configurations of the intensity model.

        Both are continuous equations with a unit intercept, AR(1) errors, the
        horseshoe linear moderator term, a learned trajectory scale and the
        unit-anchored moderator; Stage 1 keeps the absorbing offer clock, Stage 2
        swaps the exposure clock for the distributed-lag kernel.
        """
        dose_cfg = dataclasses.replace(self.config, outcome="continuous", num_categories=None,
                                       anchor_moderator=True)
        outcome_cfg = dataclasses.replace(self.config, outcome="continuous", num_categories=None,
                                          anchor_moderator=True)
        return dose_cfg, outcome_cfg

    def _resolved_lag_horizon(self) -> int:
        panel = validate_intensity(self._data["z"], self._data["d"], self._data["t"])
        return (default_lag_horizon(panel.z.shape[1], panel.start) if self.lag_horizon is None
                else int(self.lag_horizon))

    # -- fitting --------------------------------------------------------------
    def fit(self, y: Any, d: Any, z: Any, x: Any, t: Any = None,
            x_trt: Any = None, key: Any = None, strata: Any = None) -> LongBetEncourage:
        """Fit the two equations: outcome on adoption and adoption hazard on the offer (absorbing),
        or dose on the offer and outcome on the dose lags (intensity).

        ``strata`` (intensity only) are baseline randomization strata such as
        sites: the design-based reference and the dip screen then analyse the
        experiment as it was randomized. They do not enter the sampler; put
        site indicators in ``x`` for that.
        """
        if self.treatment == "intensity":
            return self._fit_intensity(y, d, z, x, t, x_trt, key, strata)
        if strata is not None:
            raise ValueError("strata is used by treatment='intensity' only.")
        data = self._inputs(y, d, z, x, t, x_trt)
        outcome_cfg, adoption_cfg = self._equation_configs()
        import jax
        if key is not None:
            key_y, key_d = jax.random.split(key)
        else:
            # The two equations must draw independent random streams: composed
            # draws pair draw b of one fit with draw b of the other, which
            # samples the product of the two posteriors only if the chains are
            # independent. Seeding both from ``random_seed`` made chain c of one
            # equation reuse chain c's random numbers in the other. The outcome
            # equation keeps the stream ``random_seed`` gives it; the adoption
            # equation's is folded off that key.
            key_y = jax.random.key(outcome_cfg.random_seed)
            key_d = jax.random.fold_in(key_y, 1)
        with warnings.catch_warnings():
            # Units adopting in the first period are 'always treated' for the
            # outcome equation; their exposure clock is right (adoption at
            # period 1), so the scalar model's advice to drop them does not apply.
            warnings.filterwarnings("ignore", message=r"\d+ unit\(s\) are treated in every period,.*",
                                    category=UserWarning)
            self.outcome_model = LongBet(outcome_cfg).fit(
                data["y"], data["x"], data["d"], t=data["t"], x_trt=data.get("x_trt"), key=key_y)
            self.adoption_model = LongBet(adoption_cfg).fit(
                self._first_adoption_events(data["d"]), data["x"], data["z"], t=data["t"],
                x_trt=data.get("x_trt"), key=key_d)
        self._data = data
        for value in self._data.values():
            value.flags.writeable = False
        self.metadata = dict(_validate(z, d, t).metadata(),
                             target="all_original_units_equal_weight", engine=ENGINE,
                             outcome=self.outcome, inference_version=INFERENCE_VERSION,
                             unit_intercept_covariance="independent_prior",
                             calibration_status="not_established",
                             component_streams=COMPONENT_STREAMS)
        return self

    # -- continuous intensity: the two stages ----------------------------------
    @staticmethod
    def _split_keys(key: Any, seed: int) -> tuple[Any, Any]:
        """``(key_y, key_d)``: independent streams for the two stages (see ``fit``)."""
        import jax
        if key is not None:
            return tuple(jax.random.split(key))
        key_y = jax.random.key(seed)
        return key_y, jax.random.fold_in(key_y, 1)

    def _stage1_summaries(self, data: dict) -> tuple[np.ndarray, np.ndarray]:
        """Posterior-mean conditional dose ``m_D(X_i, Z_i, t)`` and offer-induced dose ``CITT_D``."""
        panel = validate_intensity(data["z"], data["d"], data["t"])
        offered = np.zeros_like(data["z"], dtype=np.float32)
        offered[:, panel.start:] = 1.0
        saved_n = self.dose_model.N_
        try:
            # First-stage innovations v_{it} = D_{it} - m_D(X_i, Z_i, t) are defined relative to the
            # covariate-and-arm conditional expectation (without the unobserved unit intercept gamma_i).
            self.dose_model.N_ = -1
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message=r"predict\(\) received .* units", category=UserWarning)
                offer = self.dose_model.predict(data["x"], offered, t=data["t"], x_trt=data.get("x_trt"),
                                                summary_only=True)
        finally:
            self.dose_model.N_ = saved_n
        citt_d_mean = np.asarray(offer.tau_summary.mean, dtype=np.float64)
        citt_d_mean[:, : panel.start] = 0.0
        d0_hat = np.asarray(offer.mu0_summary.mean, dtype=np.float64)
        self._d0_hat = d0_hat
        d_hat = d0_hat + (data["z"] == 1) * citt_d_mean
        arm = (data["z"] == 1).any(axis=1)
        strata = data.get("strata")
        cells = arm if strata is None else np.char.add(np.asarray(strata, dtype=str), np.where(arm, ":1", ":0"))
        for cell in np.unique(cells):
            idx = cells == cell
            if idx.any():
                d_hat[idx] += (data["d"][idx] - d_hat[idx]).mean(axis=0, keepdims=True)
        return d_hat, citt_d_mean

    def _stage2_design(self, data: dict, d_hat: np.ndarray, citt_d_mean: np.ndarray,
                       d0_hat: np.ndarray | None = None,
                       ) -> tuple[np.ndarray, np.ndarray | None, list[str]]:
        """``(lag matrix (N,T,L+1), extra columns (N,T,Q) or None, extra column names)`` for the tier."""
        panel = validate_intensity(data["z"], data["d"], data["t"])
        lag_horizon = (default_lag_horizon(panel.z.shape[1], panel.start) if self.lag_horizon is None
                       else int(self.lag_horizon))
        x_tilde = standardized_covariates(data["x"])
        if self.identification == "deconvolution":
            source = offer_induced_dose(citt_d_mean, data["z"], strata=data.get("strata"))
            lags = lag_matrix(source, lag_horizon)
            d0_surface = d0_hat if d0_hat is not None else getattr(self, "_d0_hat", None)
            if d0_surface is None or d0_surface.shape != d_hat.shape:
                d0_surface = d_hat - (data["z"] == 1) * citt_d_mean
            extra, names = organic_dose_columns(d0_surface, lag_horizon, x_tilde)
            return lags, extra, names
        dose = pre_offer_demean(data["d"], panel.start) if self.demean_pre else data["d"]
        lags = lag_matrix(dose, lag_horizon)
        if self.identification == "control_function":
            v_hat = data["d"] - d_hat
            extra, names = control_function_columns(v_hat, lag_horizon, x_tilde, panel.start)
        else:
            extra, names = None, []
        if self.demean_pre and panel.start > 0:
            pre_cols, pre_names = pre_offer_columns(data["d"], panel.start, x_tilde)
            extra = pre_cols if extra is None else np.concatenate([extra, pre_cols], axis=2)
            names = list(names) + pre_names
        return lags, extra, names

    def _fit_stage2(self, data: dict, lags: np.ndarray, extra: np.ndarray | None, names: list[str],
                    outcome_cfg: LongBetConfig, key_y: Any) -> Any:
        if self.stage2_backend == "linear":
            return LinearLagStage2().fit(
                data["y"], lags, extra, num_draws=outcome_cfg.num_sweeps, num_chains=outcome_cfg.num_chains,
                random_seed=int(outcome_cfg.random_seed), extra_names=names)
        keep_extra = None
        extra_fit = None
        if extra is not None and extra.shape[2] > 0:
            rms = np.sqrt(np.mean(np.square(np.asarray(extra, dtype=np.float64)), axis=(0, 1)))
            keep_extra = rms > 1e-12
            if keep_extra.any():
                extra_fit = extra[:, :, keep_extra].astype(np.float32)
        model = LongBet(outcome_cfg).fit(
            data["y"], data["x"], None, t=data["t"], x_trt=data.get("x_trt"), key=key_y,
            exposure_lags=lags.astype(np.float32), extra_cols=extra_fit)
        return LongBetLagStage2(model, keep_extra=keep_extra)

    def _fit_intensity(self, y: Any, d: Any, z: Any, x: Any, t: Any, x_trt: Any, key: Any,
                       strata: Any) -> LongBetEncourage:
        data = self._inputs_intensity(y, d, z, x, t, x_trt, strata)
        dose_cfg, outcome_cfg = self._stage_configs()
        key_y, key_d = self._split_keys(key, outcome_cfg.random_seed)
        with warnings.catch_warnings():
            # Without pre-offer periods the offered units are 'treated in every
            # period' of the dose equation; their offer clock is right.
            warnings.filterwarnings("ignore", message=r"\d+ unit\(s\) are treated in every period,.*",
                                    category=UserWarning)
            self.dose_model = LongBet(dose_cfg).fit(
                data["d"], data["x"], data["z"], t=data["t"], x_trt=data.get("x_trt"), key=key_d)
        d_hat, citt_d_mean = self._stage1_summaries(data)
        lags, extra, names = self._stage2_design(data, d_hat, citt_d_mean, getattr(self, "_d0_hat", None))
        self.stage2 = self._fit_stage2(data, lags, extra, names, outcome_cfg, key_y)
        self.outcome_model = getattr(self.stage2, "model", None)
        self.adoption_model = self.dose_model
        self._data = data
        for value in self._data.values():
            value.flags.writeable = False
        panel = validate_intensity(data["z"], data["d"], data["t"])
        self.metadata = dict(panel.metadata(), target="all_original_units_equal_weight",
                             engine=ENGINE_INTENSITY, outcome=self.outcome, treatment=self.treatment,
                             identification=self.identification, tier=TIER_LABELS[self.identification],
                             lag_horizon=int(lags.shape[2] - 1), demean_pre=self.demean_pre,
                             stage2_backend=self.stage2_backend, extra_columns=list(names),
                             n_strata=int(np.unique(data["strata"]).size) if "strata" in data else 1,
                             inference_version=INFERENCE_VERSION_INTENSITY,
                             unit_intercept_covariance="independent_prior",
                             calibration_status="not_established", component_streams=COMPONENT_STREAMS,
                             engine_sha=_engine_sha(), engine_version=_engine_version())
        return self

    def _require_fit(self) -> None:
        if self.treatment == "intensity":
            if self.dose_model is None or self.stage2 is None:
                raise RuntimeError("LongBetEncourage must be fitted before prediction.")
            return
        if self.outcome_model is None or self.adoption_model is None:
            raise RuntimeError("LongBetEncourage must be fitted before prediction.")

    # -- population effects ---------------------------------------------------
    def predict(self, *, summary_only: bool = True, groups: Any = None, weights: Any = None,
                block_size: int | None = None, alpha: float = .05) -> EncouragementPrediction:
        """Offer effects standardized over the study units, with the design-based reference.

        Every study unit is evaluated under both schedules (offered from the
        common start onward, never offered) with its own fitted intercepts, so
        the target is the finite study population the reference describes.
        ``groups`` are prespecified baseline labels; ``weights`` are fixed
        nonnegative unit weights normalized within each group.
        """
        self._require_fit()
        _critical(alpha)
        if self.treatment == "intensity":
            return self._predict_intensity(groups=groups, weights=weights, block_size=block_size, alpha=alpha)
        del summary_only  # the standardized draws are small; kept for API symmetry
        data = self._data
        panel = _validate(data["z"], data["d"], data["t"])
        n_units, n_periods = data["z"].shape
        labels, counts, target = _target_weights(n_units, groups, weights)
        post_idx = np.arange(panel.start, n_periods)
        num_chains = self.config.num_chains
        n_draws = self.config.num_sweeps * num_chains
        offered = np.zeros((n_units, n_periods), dtype=np.float32)
        offered[:, panel.start:] = 1.0
        never = np.zeros_like(offered)
        if block_size is None:
            block_size = max(1, int(2.5e8 // (8 * n_periods * n_draws)))
        elif not isinstance(block_size, (int, np.integer)) or block_size < 1:
            raise ValueError("block_size must be a positive integer.")
        sums = {name: np.zeros((len(labels), len(post_idx), n_draws)) for name in ("outcome", "takeup")}
        t_vec = data["t"]
        gamma_y = (_flat_draws(self.outcome_model, "gamma")
                   if self.config.random_intercept and self.outcome == "binary"
                   else np.zeros((n_draws, n_units)))
        x_trt_all = data.get("x_trt")
        for lo in range(0, n_units, block_size):
            sl = slice(lo, min(lo + block_size, n_units))
            x_blk = data["x"][sl]
            x_trt_blk = None if x_trt_all is None else x_trt_all[sl]
            p1, stock1 = adoption_distribution(self.adoption_model, x_blk, offered[sl], t_vec,
                                               x_trt=x_trt_blk)
            p0, stock0 = adoption_distribution(self.adoption_model, x_blk, never[sl], t_vec,
                                               x_trt=x_trt_blk)
            citt_y = offer_effect_on_outcome(self.outcome_model, x_blk, t_vec, p1, p0,
                                             x_trt=x_trt_blk, fitted_intercepts=gamma_y[:, sl])
            citt_d = stock1 - stock0
            for name, values in (("outcome", citt_y), ("takeup", citt_d)):
                sums[name] += np.einsum("gn,nhd->ghd", target[:, sl], values[:, post_idx, :])
        draws = {name: sums[name].reshape(len(labels), len(post_idx), num_chains, -1)
                 for name in sums}
        standardized = EncouragementDraws(
            draws=draws, group_labels=labels, group_counts=counts, weights=target,
            periods=panel.t[post_idx].copy(), period_indices=post_idx,
            horizons=panel.exposure[post_idx].copy(), standardization="conditional",
            provenance=INFERENCE_VERSION)
        reference = _reference_groups(data, standardized, alpha)
        return EncouragementPrediction(standardized, reference, self.metadata, alpha=alpha)

    # -- continuous intensity: composition -------------------------------------
    def _offer_schedule(self, n_units: int, n_periods: int, start: int) -> np.ndarray:
        offered = np.zeros((n_units, n_periods), dtype=np.float32)
        offered[:, start:] = 1.0
        return offered

    def _compose_block(self, x_blk: np.ndarray, x_trt_blk: np.ndarray | None, z_blk: np.ndarray,
                       t_vec: np.ndarray, start: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """``(citt_y, citt_d, nu, kernel)`` for a block of units: ``(N_b,T,B)``, ``(N_b,T,B)``, ``(B,N_b)``, ``(N_b,B,L+1)``.

        Draw ``b`` of Stage 1 (``CITT_D``) is paired with draw ``b`` of Stage 2
        (``nu_Y``, ``kappa``); the streams are independent, so the pairs are
        draws from the product of the two posteriors.
        """
        citt_d = np.asarray(_predict_quiet(self.dose_model, x_blk, z_blk, t_vec, x_trt_blk).tauhats,
                            dtype=np.float64)
        citt_d[:, :start, :] = 0.0
        kappa = self.stage2.kappa_draws                                   # (B, L+1)
        nu = self.stage2.nu_draws(x_blk, x_trt_blk)                       # (B, N_b)
        if nu.shape != (kappa.shape[0], len(x_blk)) or citt_d.shape[-1] != kappa.shape[0]:
            raise RuntimeError("the two stages returned different numbers of draws; they must share "
                               "num_chains and num_sweeps.")
        kernel = kappa[None, :, :] * nu.T[:, :, None]                     # (N_b, B, L+1)
        citt_y = compose_offer_effect(kernel, citt_d, start)
        return citt_y, citt_d, nu, kernel

    def _predict_intensity(self, *, groups: Any, weights: Any, block_size: int | None,
                           alpha: float) -> IntensityEncouragementPrediction:
        data = self._data
        panel = validate_intensity(data["z"], data["d"], data["t"])
        n_units, n_periods = data["z"].shape
        labels, counts, target = _target_weights(n_units, groups, weights)
        post_idx = np.arange(panel.start, n_periods)
        num_chains = self.config.num_chains
        n_draws = self.config.num_sweeps * num_chains
        n_lags = self._resolved_lag_horizon() + 1
        offered = self._offer_schedule(n_units, n_periods, panel.start)
        if block_size is None:
            block_size = max(1, int(2.5e8 // (8 * n_periods * n_draws)))
        elif isinstance(block_size, (bool, np.bool_)) or not isinstance(block_size, (int, np.integer)) or block_size < 1:
            raise ValueError("block_size must be a positive integer.")
        sums = {name: np.zeros((len(labels), len(post_idx), n_draws)) for name in ("outcome", "takeup")}
        nu_bar = np.zeros((len(labels), n_draws))
        acr_num = np.zeros((len(labels), n_draws))
        acr_den = np.zeros((len(labels), n_draws))
        x_trt_all = data.get("x_trt")
        for lo in range(0, n_units, block_size):
            sl = slice(lo, min(lo + block_size, n_units))
            x_trt_blk = None if x_trt_all is None else x_trt_all[sl]
            citt_y, citt_d, nu, _ = self._compose_block(data["x"][sl], x_trt_blk, offered[sl], data["t"],
                                                        panel.start)
            for name, values in (("outcome", citt_y), ("takeup", citt_d)):
                sums[name] += np.einsum("gn,nhd->ghd", target[:, sl], values[:, post_idx, :])
            compliance = citt_d[:, post_idx, :].sum(axis=1)               # (N_b, B): cumulative CITT_D
            nu_bar += target[:, sl] @ nu.T
            acr_num += target[:, sl] @ (compliance * nu.T)
            acr_den += target[:, sl] @ compliance
        draws = {name: sums[name].reshape(len(labels), len(post_idx), num_chains, -1) for name in sums}
        kappa = self.stage2.kappa_draws                                   # (B, L+1)
        with np.errstate(divide="ignore", invalid="ignore"):
            acr_weight = np.where(acr_den != 0, acr_num / acr_den, np.nan)
        kernel_draws = {label: (kappa.T * nu_bar[g][None, :]).reshape(n_lags, num_chains, -1)
                        for g, label in enumerate(labels)}
        kernel_acr_draws = {label: (kappa.T * acr_weight[g][None, :]).reshape(n_lags, num_chains, -1)
                            for g, label in enumerate(labels)}
        standardized = EncouragementDraws(
            draws=draws, group_labels=labels, group_counts=counts, weights=target,
            periods=panel.t[post_idx].copy(), period_indices=post_idx,
            horizons=panel.exposure[post_idx].copy(), standardization="conditional",
            provenance=INFERENCE_VERSION_INTENSITY)
        reference, references = _intensity_reference_groups(
            data, standardized, alpha, lag_horizon=n_lags - 1, demean_pre=self.demean_pre)
        return IntensityEncouragementPrediction(
            standardized, reference, self.metadata, kernel_draws=kernel_draws,
            kernel_acr_draws=kernel_acr_draws, references=references, alpha=alpha)

    def _predict_conditional_intensity(self, x_np: np.ndarray, x_trt_np: np.ndarray | None,
                                       z_mat: np.ndarray, t_vec: np.ndarray, start: int,
                                       block_size: int, alpha: float) -> IntensityConditionalPrediction:
        n_units, n_periods = z_mat.shape
        blocks_y, blocks_d, blocks_k = [], [], []
        for lo in range(0, n_units, block_size):
            sl = slice(lo, min(lo + block_size, n_units))
            x_trt_blk = None if x_trt_np is None else x_trt_np[sl]
            citt_y, citt_d, _, kernel = self._compose_block(x_np[sl], x_trt_blk, z_mat[sl], t_vec, start)
            blocks_y.append(citt_y)
            blocks_d.append(citt_d)
            blocks_k.append(kernel.transpose(0, 2, 1))                   # (N_b, L+1, B)
        citt_y = np.concatenate(blocks_y, axis=0)
        citt_d = np.concatenate(blocks_d, axis=0)
        kernel = np.concatenate(blocks_k, axis=0)
        post_idx = np.arange(start, n_periods)
        return IntensityConditionalPrediction(
            citt_y_draws=citt_y[:, post_idx, :], citt_d_draws=citt_d[:, post_idx, :], kernel_draws=kernel,
            periods=np.asarray(t_vec)[post_idx].copy(),
            horizons=np.asarray(t_vec)[post_idx] - np.asarray(t_vec)[start] + 1, alpha=alpha)

    # -- unit-level effects ---------------------------------------------------
    def predict_conditional(
        self, x: Any = None, z: Any = None, t: Any = None, *, x_trt: Any = None,
        alpha: float = 0.05, cace_stabilization: float = 0.02, monotonic_first_stage: bool = True,
        block_size: int | None = None,
    ) -> ConditionalEncouragementPrediction:
        """Unit-level offer effects, principal strata, lift and policy quantities.

        Describes a unit with covariate profile ``x`` drawn from the fitted
        population: both unit intercepts are integrated out, so the same rule
        applies to study units and to a new cohort. ``x`` defaults to the study
        units, ``z`` to an offer from the study's launch period onward, ``t``
        to the fitted calendar. The comparison schedule is never offering.
        Adoption contrasts are clipped at zero when ``monotonic_first_stage``
        (no defiers); the CACE is the draw-wise ratio of the two offer effects
        with ``cace_stabilization`` added to the adoption contrast.

        For ``treatment="intensity"`` the result is an
        :class:`IntensityConditionalPrediction`: the dose contrast is never
        clipped, ``cace_stabilization`` and ``monotonic_first_stage`` are
        ignored, and the ratio reported is the cumulative return per unit of
        cumulative dose.
        """
        self._require_fit()
        _critical(alpha)
        x_np = self._data["x"] if x is None else np.asarray(x, dtype=np.float32)
        if x_np.ndim != 2 or x_np.shape[1] != self._data["x"].shape[1]:
            raise ValueError(f"x must be a 2-D array with {self._data['x'].shape[1]} columns.")
        fitted_x_trt = self._data.get("x_trt")
        if x is None:
            x_trt_np = fitted_x_trt
        elif fitted_x_trt is None:
            if x_trt is not None:
                raise ValueError("the model was fitted without x_trt.")
            x_trt_np = None
        else:
            if x_trt is None:
                raise ValueError("the model was fitted with x_trt, so new units need x_trt too.")
            x_trt_np = _finite_matrix(x_trt, "x_trt", len(x_np), (len(x_np), fitted_x_trt.shape[1]))
        t_vec = self._data["t"] if t is None else np.asarray(t, dtype=np.float64)
        n_units, n_periods = len(x_np), len(t_vec)
        panel = self._panel()
        if z is None:
            if t is None:
                start = panel.start
            else:
                launch_t = float(self._data["t"][panel.start])
                cand = np.flatnonzero(t_vec >= launch_t)
                start = int(cand[0]) if len(cand) > 0 else 0
            z_mat = np.zeros((n_units, n_periods), dtype=np.float32)
            z_mat[:, start:] = 1.0
        else:
            z_mat = np.asarray(z, dtype=np.float32)
            if z_mat.shape != (n_units, n_periods):
                raise ValueError(f"z must have shape ({n_units}, {n_periods}).")
            pos_cols = np.flatnonzero((z_mat > 0).any(axis=0))
            start = int(pos_cols[0]) if len(pos_cols) > 0 else min(panel.start, n_periods - 1)
        never = np.zeros_like(z_mat)
        if block_size is None:
            block_size = n_units
        elif isinstance(block_size, (bool, np.bool_)) or not isinstance(block_size, (int, np.integer)) or block_size < 1:
            raise ValueError("block_size must be a positive integer.")
        block_size = min(int(block_size), n_units)
        if self.treatment == "intensity":
            return self._predict_conditional_intensity(x_np, x_trt_np, z_mat, t_vec, start, block_size, alpha)
        if block_size == n_units:
            p1, stock1 = adoption_distribution(self.adoption_model, x_np, z_mat, t_vec, x_trt=x_trt_np)
            p0, stock0 = adoption_distribution(self.adoption_model, x_np, never, t_vec, x_trt=x_trt_np)
            citt_y = offer_effect_on_outcome(self.outcome_model, x_np, t_vec, p1, p0, x_trt=x_trt_np)
        else:
            citt_y_blocks, stock1_blocks, stock0_blocks = [], [], []
            for lo in range(0, n_units, block_size):
                sl = slice(lo, min(lo + block_size, n_units))
                x_blk = x_np[sl]
                x_trt_blk = None if x_trt_np is None else x_trt_np[sl]
                p1_b, s1_b = adoption_distribution(self.adoption_model, x_blk, z_mat[sl], t_vec, x_trt=x_trt_blk)
                p0_b, s0_b = adoption_distribution(self.adoption_model, x_blk, never[sl], t_vec, x_trt=x_trt_blk)
                cy_b = offer_effect_on_outcome(self.outcome_model, x_blk, t_vec, p1_b, p0_b, x_trt=x_trt_blk)
                citt_y_blocks.append(cy_b)
                stock1_blocks.append(s1_b)
                stock0_blocks.append(s0_b)
            citt_y = np.concatenate(citt_y_blocks, axis=0)
            stock1 = np.concatenate(stock1_blocks, axis=0)
            stock0 = np.concatenate(stock0_blocks, axis=0)
        post_idx = np.arange(start, n_periods)
        return ConditionalEncouragementPrediction(
            citt_y_draws=citt_y[:, post_idx, :],
            citt_d_draws=(stock1 - stock0)[:, post_idx, :],
            periods=np.asarray(t_vec)[post_idx].copy(),
            horizons=np.asarray(t_vec)[post_idx] - np.asarray(t_vec)[start] + 1,
            p0_draws=stock0[:, post_idx, :],
            alpha=alpha, cace_stabilization=cace_stabilization,
            monotonic_first_stage=monotonic_first_stage,
        )

    # -- paired bootstrap comparison -----------------------------------------
    def bootstrap_comparison(self, *, replicates: int = 200, groups: Any = None,
                             random_seed: int = 0, alpha: float = .05,
                             min_ess: float = 400, max_rhat: float = 1.01,
                             callback: Any = None) -> EncouragementComparison:
        """Refit paired whole-unit bootstrap samples to compare model and reference ITTs.

        Resample entire unit records within assigned arm and prespecified
        group, holding their counts fixed; refit both equations with new
        sampler seeds; recompute both estimators on the same resampled
        population, so the comparison variance includes their sampling
        covariance. Only the ITTs are compared. This requires ``replicates``
        complete fits. A disagreement flag is withheld if the original fit or
        any replicate fails diagnostics, if reference inference is unavailable,
        or if fewer than 100 replications are requested.
        """
        critical = _critical(alpha)
        if isinstance(replicates, (bool, np.bool_)) or not isinstance(replicates, (int, np.integer)) or replicates < 2:
            raise ValueError("replicates must be an integer of at least 2.")
        if callback is not None and not callable(callback):
            raise ValueError("callback must be callable or None.")
        original = self.predict(groups=groups, alpha=alpha)
        original_checks = original.stability(min_ess=min_ess, max_rhat=max_rhat)
        itts = ["itt_y", "itt_d"]
        original_ok = original_checks.query("quantity in @itts").diagnostics_passed.all()
        reference_columns = ["itt_y_se", "itt_d_se", "itt_y_d_cov"]
        original_reference_ok = bool(np.isfinite(original.reference[reference_columns]).all().all())
        n = len(self._data["z"])
        group_values = np.zeros(n) if groups is None else np.asarray(groups)
        strata = []
        for label in pd.unique(group_values):
            for arm in (0, 1):
                idx = np.flatnonzero((group_values == label) & (self._data["z"][:, -1] == arm))
                if idx.size:
                    strata.append(idx)
        rng = np.random.default_rng(random_seed)
        paired = np.full((replicates, len(original.group_labels), len(original.horizons), 2, 2), np.nan)
        failures = []
        for b in range(replicates):
            indices = np.concatenate([rng.choice(idx, len(idx), replace=True) for idx in strata])
            seed = int(rng.integers(0, 2**31))
            try:
                fit = self._clone(random_seed=seed)
                fit.fit(self._data["y"][indices], self._data["d"][indices], self._data["z"][indices],
                        self._data["x"][indices], self._data["t"],
                        x_trt=None if "x_trt" not in self._data else self._data["x_trt"][indices],
                        strata=None if "strata" not in self._data else self._data["strata"][indices])
                pred = fit.predict(groups=None if groups is None else group_values[indices], alpha=alpha)
                for g, label in enumerate(original.group_labels):
                    pg = pred.group_labels.index(label)
                    ref = pred.reference.loc[pred.reference.group == label]
                    for q, quantity in enumerate(("itt_y", "itt_d")):
                        paired[b, g, :, q, 0] = pred.draws[quantity][pg].mean(axis=(-2, -1))
                        paired[b, g, :, q, 1] = ref[quantity].to_numpy()
                checks = pred.stability(min_ess=min_ess, max_rhat=max_rhat)
                if not checks.query("quantity in @itts").diagnostics_passed.all():
                    failures.append(dict(replicate=b, sampler_seed=seed, reason="diagnostics_failed"))
                if not np.isfinite(paired[b]).all():
                    failures.append(dict(replicate=b, sampler_seed=seed, reason="unavailable_estimate"))
                if not np.isfinite(pred.reference[reference_columns]).all().all():
                    failures.append(dict(replicate=b, sampler_seed=seed,
                                         reason="unavailable_reference_inference"))
            except (ValueError, RuntimeError, FloatingPointError) as exc:
                failures.append(dict(replicate=b, sampler_seed=seed,
                                     reason=f"{type(exc).__name__}: {exc}"))
            if callback is not None:
                callback(b + 1, replicates)
        frame = original.comparison()
        frame["model_sampling_variance"] = np.nan
        frame["reference_sampling_variance"] = np.nan
        frame["model_reference_covariance"] = np.nan
        complete = not failures and bool(original_ok) and original_reference_ok and replicates >= 100
        for r, row in frame.iterrows():
            g = original.group_labels.index(row.group)
            h = int(np.flatnonzero(original.horizons == row.horizon)[0])
            q = ("itt_y", "itt_d").index(row.quantity)
            pairs = paired[:, g, h, q]
            if np.isfinite(pairs).all():
                cov = np.cov(pairs, rowvar=False, ddof=1)
                se = float(np.std(pairs[:, 0] - pairs[:, 1], ddof=1))
                frame.loc[r, ["model_sampling_variance", "reference_sampling_variance",
                              "model_reference_covariance", "difference_se"]] = (
                                  cov[0, 0], cov[1, 1], cov[0, 1], se)
                if complete and np.isfinite(row.difference) and np.isfinite(se) and se > 0:
                    frame.loc[r, "disagreement"] = bool(abs(row.difference) > critical * se)
                    frame.loc[r, "comparison_status"] = "paired_unit_bootstrap_normal_pointwise"
                else:
                    frame.loc[r, "comparison_status"] = "bootstrap_checks_incomplete"
            else:
                frame.loc[r, "comparison_status"] = "bootstrap_failed_estimates"
        return EncouragementComparison(frame, paired,
            pd.DataFrame(failures, columns=["replicate", "sampler_seed", "reason"]),
            dict(method="paired_unit_bootstrap", replicates=int(replicates), random_seed=random_seed,
                 alpha=alpha, original_diagnostics_passed=bool(original_ok),
                 original_reference_inference_available=original_reference_ok,
                 resampling_unit="entire_unit_record", strata="assigned_arm_and_baseline_group",
                 target="same_resampled_units_for_both_estimators", exact_randomization=False,
                 calibration_status="not_established"))

    # -- helpers shared by both treatments -------------------------------------
    def _clone(self, **config_overrides: Any) -> LongBetEncourage:
        """An unfitted model with this one's treatment, tiers and options."""
        kwargs: dict[str, Any] = dict(outcome=self.outcome, treatment=self.treatment)
        if self.treatment == "intensity":
            kwargs.update(identification=self.identification, lag_horizon=self.lag_horizon,
                          demean_pre=self.demean_pre, stage2_backend=self.stage2_backend)
        return LongBetEncourage(self.config, **kwargs, **config_overrides)

    def _panel(self) -> Any:
        """The validated fitted panel (``.start`` is the offer index under either treatment)."""
        if self.treatment == "intensity":
            return validate_intensity(self._data["z"], self._data["d"], self._data["t"])
        return _validate(self._data["z"], self._data["d"], self._data["t"])

    # -- persistence ----------------------------------------------------------
    def save(self, path: str | Path) -> None:
        """Save a versioned, pickle-free archive with both equations and the study data."""
        self._require_fit()
        if self.treatment == "intensity":
            return self._save_intensity(path)
        arrays = dict(self._data)
        with tempfile.TemporaryDirectory() as tmp:
            for name, model in (("outcome_model", self.outcome_model),
                                ("adoption_model", self.adoption_model)):
                nested = Path(tmp) / f"{name}.npz"
                model.save(nested)
                arrays[name] = np.frombuffer(nested.read_bytes(), dtype=np.uint8)
        meta = dict(self.metadata, kind="LongBetEncourage", archive_version=ARCHIVE_VERSION,
                    config=self.config.to_dict(),
                    input_hashes={name: _digest(value) for name, value in arrays.items()})
        arrays["metadata"] = np.asarray(json.dumps(meta, allow_nan=False))
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as stream:
            np.savez_compressed(stream, **arrays)

    def _save_intensity(self, path: str | Path) -> None:
        arrays = dict(self._data)
        with tempfile.TemporaryDirectory() as tmp:
            nested = Path(tmp) / "dose_model.npz"
            self.dose_model.save(nested)
            arrays["dose_model"] = np.frombuffer(nested.read_bytes(), dtype=np.uint8)
            if self.stage2_backend == "longbet":
                nested = Path(tmp) / "outcome_model.npz"
                self.stage2.model.save(nested)
                arrays["outcome_model"] = np.frombuffer(nested.read_bytes(), dtype=np.uint8)
            else:
                arrays.update(self.stage2.to_arrays())
        keep_extra = getattr(self.stage2, "keep_extra", None)
        options = dict(treatment=self.treatment, identification=self.identification,
                       lag_horizon=self.lag_horizon, demean_pre=self.demean_pre,
                       stage2_backend=self.stage2_backend,
                       extra_keep=None if keep_extra is None else np.asarray(keep_extra, dtype=bool).tolist())
        meta = dict(self.metadata, kind="LongBetEncourage", archive_version=ARCHIVE_VERSION_INTENSITY,
                    config=self.config.to_dict(), options=options,
                    input_hashes={name: _digest(value) for name, value in arrays.items()})
        arrays["metadata"] = np.asarray(json.dumps(meta, allow_nan=False))
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as stream:
            np.savez_compressed(stream, **arrays)

    @classmethod
    def _load_intensity(cls, meta: dict, arrays: dict) -> LongBetEncourage:
        if (meta.get("archive_version") != ARCHIVE_VERSION_INTENSITY
                or meta.get("inference_version") != INFERENCE_VERSION_INTENSITY):
            raise ValueError("unsupported kind or version")
        options = dict(meta["options"])
        backend = options.get("stage2_backend", "longbet")
        required = {"y", "d", "z", "x", "t", "dose_model"}
        optional = {"x_trt", "strata"}
        if backend == "longbet":
            required |= {"outcome_model"}
        else:
            required |= {"linear_kappa", "linear_sigma2", "linear_num_chains"}
            optional |= {"linear_extra_coef"}
        if not required <= arrays.keys() or arrays.keys() - required - optional:
            raise ValueError("unexpected or missing archive arrays")
        if meta["input_hashes"] != {name: _digest(value) for name, value in arrays.items()}:
            raise ValueError("input or model digest mismatch")
        config = LongBetConfig.from_dict(meta["config"])
        obj = cls(config, outcome=meta["outcome"], treatment="intensity",
                  identification=options["identification"], lag_horizon=options.get("lag_horizon"),
                  demean_pre=bool(options.get("demean_pre", False)), stage2_backend=backend)
        with tempfile.TemporaryDirectory() as tmp:
            blob = arrays.pop("dose_model")
            if blob.dtype != np.uint8 or blob.ndim != 1:
                raise ValueError("invalid nested model bytes")
            nested = Path(tmp) / "dose_model.npz"
            nested.write_bytes(blob.tobytes())
            dose_model = LongBet.load(nested)
            outcome_model = None
            if backend == "longbet":
                blob = arrays.pop("outcome_model")
                if blob.dtype != np.uint8 or blob.ndim != 1:
                    raise ValueError("invalid nested model bytes")
                nested = Path(tmp) / "outcome_model.npz"
                nested.write_bytes(blob.tobytes())
                outcome_model = LongBet.load(nested)
        linear = {name: arrays.pop(name) for name in list(arrays) if name.startswith("linear_")}
        strata = arrays.pop("strata", None)
        obj._data = obj._inputs_intensity(**{**arrays, "x_trt": arrays.get("x_trt"), "strata": strata})
        panel = validate_intensity(arrays["z"], arrays["d"], arrays["t"])
        dose_cfg, outcome_cfg = obj._stage_configs()
        if (dose_model.config != dose_cfg or dose_model.N_ != len(arrays["y"])
                or dose_model.T_ != arrays["y"].shape[1]
                or not np.array_equal(dose_model.t_fit_, arrays["t"].astype(np.float32))):
            raise ValueError("model/design provenance mismatch")
        if backend == "longbet":
            if (outcome_model.config != outcome_cfg or outcome_model.N_ != len(arrays["y"])
                    or outcome_model.T_ != arrays["y"].shape[1]):
                raise ValueError("model/design provenance mismatch")
            keep_extra = options.get("extra_keep")
            obj.stage2 = LongBetLagStage2(outcome_model,
                                           keep_extra=None if keep_extra is None else np.asarray(keep_extra, dtype=bool))
        else:
            obj.stage2 = LinearLagStage2.from_arrays(linear, list(meta.get("extra_columns", [])))
        if (any(meta.get(k) != v for k, v in panel.metadata().items())
                or meta.get("target") != "all_original_units_equal_weight"
                or meta.get("calibration_status") != "not_established"):
            raise ValueError("model/design provenance mismatch")
        for value in obj._data.values():
            value.flags.writeable = False
        obj.dose_model = dose_model
        obj.adoption_model = dose_model
        obj.outcome_model = outcome_model
        obj.metadata = {k: v for k, v in meta.items()
                        if k not in ("kind", "archive_version", "input_hashes", "config", "options")}
        return obj

    @classmethod
    def load(cls, path: str | Path) -> LongBetEncourage:
        """Validate and replay an encouragement archive."""
        try:
            with np.load(path, allow_pickle=False) as archive:
                meta = json.loads(str(archive["metadata"]))
                arrays = {name: archive[name] for name in archive.files if name != "metadata"}
            if isinstance(meta, dict) and meta.get("inference_version") == "adoption_clock_v1":
                raise ValueError("the archive was made by an earlier release of the encouragement model "
                                 "(inference adoption_clock_v1); load it with that release")
            if (isinstance(meta, dict) and meta.get("kind") == "LongBetEncourage"
                    and meta.get("engine") == ENGINE_INTENSITY):
                return cls._load_intensity(meta, arrays)
            if (not isinstance(meta, dict) or meta.get("kind") != "LongBetEncourage" or
                    meta.get("archive_version") != ARCHIVE_VERSION or
                    meta.get("inference_version") != INFERENCE_VERSION or
                    meta.get("engine") != ENGINE):
                raise ValueError("unsupported kind or version")
            required = {"y", "d", "z", "x", "t", "outcome_model", "adoption_model"}
            if not required <= arrays.keys() or arrays.keys() - required - {"x_trt"}:
                raise ValueError("unexpected or missing archive arrays")
            if meta["input_hashes"] != {name: _digest(value) for name, value in arrays.items()}:
                raise ValueError("input or model digest mismatch")
            config = LongBetConfig.from_dict(meta["config"])
            obj = cls(config, outcome=meta["outcome"])
            models = {}
            with tempfile.TemporaryDirectory() as tmp:
                for name in ("outcome_model", "adoption_model"):
                    blob = arrays.pop(name)
                    if blob.dtype != np.uint8 or blob.ndim != 1:
                        raise ValueError("invalid nested model bytes")
                    nested = Path(tmp) / f"{name}.npz"
                    nested.write_bytes(blob.tobytes())
                    models[name] = LongBet.load(nested)
            obj._data = obj._inputs(**{**arrays, "x_trt": arrays.get("x_trt")})
            panel = _validate(arrays["z"], arrays["d"], arrays["t"])
            outcome_cfg, adoption_cfg = obj._equation_configs()
            for name, expected in (("outcome_model", outcome_cfg), ("adoption_model", adoption_cfg)):
                model = models[name]
                if (model.config != expected or model.N_ != len(arrays["y"])
                        or model.T_ != arrays["y"].shape[1]
                        or not np.array_equal(model.t_fit_, arrays["t"].astype(np.float32))):
                    raise ValueError("model/design provenance mismatch")
            if (any(meta.get(k) != v for k, v in panel.metadata().items())
                    or meta.get("target") != "all_original_units_equal_weight"
                    or meta.get("calibration_status") != "not_established"):
                raise ValueError("model/design provenance mismatch")
            for value in obj._data.values():
                value.flags.writeable = False
            obj.outcome_model, obj.adoption_model = models["outcome_model"], models["adoption_model"]
            obj.metadata = {k: v for k, v in meta.items()
                            if k not in ("kind", "archive_version", "input_hashes", "config")}
            if meta.get("component_streams") != COMPONENT_STREAMS:
                warnings.warn(
                    "this archive was fitted by an engine release that drew both equations from the "
                    "same random streams, so its paired draws are correlated and composed quantities "
                    "are not exact draws from the product of the two posteriors; refit it to obtain "
                    "independent streams", UserWarning, stacklevel=2)
            return obj
        except (KeyError, TypeError, ValueError, OSError) as exc:
            raise ValueError(f"Invalid encouragement archive: {exc}") from exc


def _digest(a: np.ndarray) -> str:
    digest = hashlib.sha256()
    digest.update(str(a.dtype).encode())
    digest.update(str(a.shape).encode())
    digest.update(np.ascontiguousarray(a).tobytes())
    return digest.hexdigest()


def _engine_version() -> str:
    try:
        from longbet import __version__
        return str(__version__)
    except Exception:  # pragma: no cover - import cycle guard
        return "unknown"


def _engine_sha() -> str:
    """The engine's git commit when running from a checkout; ``"unknown"`` from an installed wheel."""
    import subprocess
    try:
        out = subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=5, check=False)
        sha = out.stdout.strip()
        return sha if out.returncode == 0 and len(sha) >= 7 else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


_INTENSITY_REFERENCE_COLUMNS = (
    "itt_y", "itt_y_se", "itt_y_lower", "itt_y_upper", "itt_d", "itt_d_se", "itt_d_lower", "itt_d_upper",
    "itt_y_d_cov", "wald", "wald_set_type", "wald_reason", "wald_lower_1", "wald_upper_1", "wald_lower_2",
    "wald_upper_2", "cumulative_itt_y", "cumulative_itt_d", "cumulative_itt_y_se", "cumulative_itt_d_se",
    "omega_cum", "omega_cum_se_delta", "omega_cum_set_type", "omega_cum_reason", "omega_cum_lower_1",
    "omega_cum_upper_1", "omega_cum_lower_2", "omega_cum_upper_2",
)


def _intensity_reference_groups(data: dict, standardized: EncouragementDraws, alpha: float, *,
                                lag_horizon: int, demean_pre: bool,
                                ) -> tuple[pd.DataFrame, dict[str, IntensityReference | None]]:
    """Design-based (T0) reference per baseline group on the identical unit target.

    Strata that lose an arm inside a group fall back to the unstratified
    contrast for that group; groups missing an arm keep NaN rows flagged
    ``missing_group_arm`` like the binary model.
    """
    frames: list[pd.DataFrame] = []
    references: dict[str, IntensityReference | None] = {}
    post = standardized.period_indices
    for g, label in enumerate(standardized.group_labels):
        mask = standardized.weights[g] > 0
        z = data["z"][mask]
        counts = [int((z[:, -1] == a).sum()) for a in (0, 1)]
        ref = None
        if min(counts) > 0:
            strata = data["strata"][mask] if "strata" in data else None
            try:
                ref = intensity_reference(data["y"][mask], data["d"][mask], z, data["t"],
                                          strata=strata, lag_horizon=lag_horizon, demean_pre=demean_pre,
                                          alpha=alpha)
            except ValueError:
                if strata is None:
                    raise
                ref = intensity_reference(data["y"][mask], data["d"][mask], z, data["t"], strata=None,
                                          lag_horizon=lag_horizon, demean_pre=demean_pre, alpha=alpha)
        references[label] = ref
        if ref is not None:
            table = ref.itt.loc[ref.itt.post_encouragement].reset_index(drop=True)
            frame = pd.DataFrame({"period": table.period, "period_index": table.period_index,
                                  "horizon": table.horizon, "post_encouragement": True,
                                  "n_encouraged": counts[1], "n_control": counts[0]})
            for name in _INTENSITY_REFERENCE_COLUMNS:
                source = name.replace("wald", "omega") if name.startswith("wald") else name
                frame[name] = table[source].to_numpy()
            frame["inference"] = "normal_ar"
            frame["itt_y_cum_se"] = float(table["cumulative_itt_y_se"].iloc[-1])
            frame["itt_d_cum_se"] = float(table["cumulative_itt_d_se"].iloc[-1])
            frame["condition_number"] = ref.condition_number
            frame["z_dip"] = ref.dip["z_stat"]
            frame["fosd_pvalue"] = ref.fosd["pvalue"]
        else:
            frame = pd.DataFrame({"period": data["t"][post], "period_index": post,
                                  "horizon": standardized.horizons, "post_encouragement": True,
                                  "n_encouraged": counts[1], "n_control": counts[0]})
            for name in _INTENSITY_REFERENCE_COLUMNS:
                frame[name] = np.nan
            frame["wald_set_type"] = "unavailable"
            frame["wald_reason"] = "missing_group_arm"
            frame["omega_cum_set_type"] = "unavailable"
            frame["omega_cum_reason"] = "missing_group_arm"
            frame["inference"] = "normal_ar"
            frame["itt_y_cum_se"] = np.nan
            frame["itt_d_cum_se"] = np.nan
            frame["condition_number"] = np.nan
            frame["z_dip"] = np.nan
            frame["fosd_pvalue"] = np.nan
        frame.insert(0, "group", label)
        frame["n_units"] = int(mask.sum())
        frames.append(frame)
    return pd.concat(frames, ignore_index=True), references
