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

"""Randomized encouragement with a continuous, non-absorbing treatment intensity.

This module holds the parts of the encouragement layer that do not need the
sampler: panel validation for a continuous dose, the **design-based
reference** (Tier T0 of the methodology), the lag-matrix and control-function
constructions that Stage 2 of :class:`longbet.LongBetEncourage` consumes when
``treatment="intensity"``, and the draw-by-draw composition of the two stages.

Setting
-------
A single randomized offer wave ``Z_i`` at period ``t0`` (``z[i, t] = Z_i *
1{t >= t0}`` is absorbing), a continuous encouraged behaviour ``d[i, t] >= 0``
that may rise, vary and switch off, and a continuous panel outcome ``y``. The
structural model is the distributed-lag kernel on the calendar clock,

    Y_it = mu_Y(X_i) + eta_it + nu_Y(X_i) * sum_{l=0}^{L} kappa_l D_{i,t-l} + gamma_i + eps_it,

whose offer-arm contrast is the triangular convolution

    ITT_Y(t) = sum_{l=0}^{t-t0} kappa_l ITT_D(t-l),       t >= t0,

because pre-offer doses do not respond to the offer (no anticipation) and the
unit intercept and innovations are balanced across arms by randomization.

Design-based reference (T0)
---------------------------
:func:`intensity_reference` estimates the arm contrasts ``ITT_D(t)``,
``ITT_Y(t)`` (stratum-weighted differences in means with the joint Neyman
covariance across periods and quantities), inverts the lower-triangular
Toeplitz system ``T_{ITT_D} kappa = ITT_Y`` without regularization, and reports
delta-method standard errors, the condition number of ``T_{ITT_D}``, the
cumulative kernel, the cumulative Wald return per unit of cumulative dose with
exact Fieller/Anderson--Rubin confidence sets, the joint Anderson--Rubin
``chi^2`` test for the whole kernel, the arm-wise ECDF shift of cumulative dose
with a first-order stochastic dominance check, and the dose-history-conditioned
pre-treatment dip screen ``z_dip``. Nothing here uses within-arm dose variation
to learn the kernel, so the reference is immune to within-arm dose endogeneity
and to classical measurement error in the dose; it is also the benchmark for
the discrepancy statistics ``z_{T2-T0}`` and ``z_{T1-T0}``.

What the reference identifies
-----------------------------
Under randomization, exclusion and no anticipation, linear calendar-homogeneous
lag effects, and time-invariant conditional complier weights, the population
deconvolution identifies the **compliance-weighted** kernel
``kappa^acr_l = E[nu_D(X) nu_Y(X)] / E[nu_D(X)] * kappa_l``; it equals the
unweighted population kernel when compliance and benefit are uncorrelated
across covariates. The conditional (model-based) tiers recover both.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy import stats

from longbet._encourage import _binary_panel, _critical, _wald_set
from longbet._model import _check_time_vector, derive_exposure

#: Kinds of the three identification tiers, in increasing order of what they
#: assume about within-arm dose variation.
IDENTIFICATION_TIERS: tuple[str, ...] = ("deconvolution", "control_function", "unconfounded")
TIER_LABELS: dict[str, str] = {"deconvolution": "T0", "control_function": "T1", "unconfounded": "T2"}


# ---------------------------------------------------------------------------
# Panel validation
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class IntensityPanel:
    """A validated single-wave encouragement panel with a continuous dose."""

    z: np.ndarray
    d: np.ndarray
    t: np.ndarray
    assigned: np.ndarray
    start: int
    exposure: np.ndarray

    @property
    def n_post(self) -> int:
        return int(self.z.shape[1] - self.start)

    def metadata(self) -> dict[str, Any]:
        return {
            "design": "single_wave_complete",
            "treatment": "intensity",
            "n_units": int(self.z.shape[0]),
            "n_periods": int(self.z.shape[1]),
            "n_encouraged": int(self.assigned.sum()),
            "n_control": int((~self.assigned).sum()),
            "encouragement_period": float(self.t[self.start]),
            "encouragement_index": int(self.start),
            "has_pre_periods": bool(self.start > 0),
            "dose_min": float(self.d.min()),
            "dose_max": float(self.d.max()),
            "dose_zero_share": float(np.mean(self.d == 0)),
            "dose_nonabsorbing": bool(np.any(np.diff(self.d, axis=1) < 0)),
        }


def validate_intensity(z: Any, d: Any, t: Any = None) -> IntensityPanel:
    """Validate a single-wave offer panel ``z`` and a continuous dose panel ``d``.

    ``z`` must be a complete, binary, absorbing ``(N, T)`` offer indicator with
    one encouragement wave and a permanent control arm (as in
    :func:`longbet.validate_encouragement`). ``d`` is a complete finite
    ``(N, T)`` dose panel; it need not be binary, monotone or nonnegative
    (a pre-offer-demeaned dose is negative in places), although the model's
    zero-dose boundary condition is stated for raw nonnegative doses.
    """
    z = _binary_panel(z, "z")
    d = np.asarray(d)
    if d.dtype.kind not in "biuf" or d.shape != z.shape:
        raise ValueError(f"d must be a finite numeric dose panel with shape {z.shape}.")
    d = d.astype(np.float64)
    if not np.all(np.isfinite(d)):
        raise ValueError("d must be complete and finite; missing cells are not silently dropped.")
    t = (np.arange(1, z.shape[1] + 1, dtype=np.float64) if t is None
         else np.atleast_1d(np.asarray(t, dtype=np.float64)))
    if t.ndim != 1 or t.size != z.shape[1]:
        raise ValueError(f"t must be a 1-D vector of length {z.shape[1]}.")
    if np.all(np.isfinite(t)) and not np.allclose(np.diff(t), np.rint(np.diff(t)), rtol=0, atol=1e-8):
        raise ValueError("t must have whole-unit gaps for lag indexing.")
    _check_time_vector(t)
    if z.shape[1] > 1 and not np.allclose(np.diff(t), 1.0, rtol=0, atol=1e-8):
        raise ValueError("the distributed-lag model needs consecutive unit-spaced periods; "
                         "reindex t before fitting a continuous-intensity panel.")
    assigned = z[:, -1] == 1
    if not assigned.any() or assigned.all():
        raise ValueError("z must include both an encouraged arm and a permanent control arm.")
    starts = np.argmax(z[assigned] == 1, axis=1)
    if np.unique(starts).size != 1:
        raise ValueError("Only a single encouragement wave is supported; z has multiple start periods.")
    start = int(starts[0])
    exposure = derive_exposure(z, t)[np.flatnonzero(assigned)[0]]
    return IntensityPanel(z, d, t, assigned, start, exposure)


def _finite_outcome(y: Any, shape: tuple[int, int]) -> np.ndarray:
    y = np.asarray(y)
    if y.shape != shape or y.dtype.kind not in "biuf":
        raise ValueError(f"y must be a finite numeric array with shape {shape}.")
    y = y.astype(np.float64)
    if not np.all(np.isfinite(y)):
        raise ValueError("y must be complete and finite; missing cells are not silently dropped.")
    return y


def _strata_labels(strata: Any, n: int) -> np.ndarray:
    if strata is None:
        return np.zeros(n, dtype=np.int64)
    s = np.asarray(strata)
    if s.shape != (n,) or np.any(pd.isna(s)):
        raise ValueError(f"strata must be a complete length-{n} vector of baseline labels.")
    return pd.factorize(s.astype(str) if s.dtype.kind in "OUS" else s, sort=True)[0]


def _covariates(x: Any, n: int) -> np.ndarray | None:
    if x is None:
        return None
    a = np.asarray(x, dtype=np.float64)
    if a.ndim == 1:
        a = a[:, None]
    if a.ndim != 2 or a.shape[0] != n or not np.all(np.isfinite(a)):
        raise ValueError(f"x must be a finite numeric array with {n} rows.")
    return a


# ---------------------------------------------------------------------------
# Lag matrices, control-function columns and composition
# ---------------------------------------------------------------------------

def default_lag_horizon(n_periods: int, start: int) -> int:
    """``L = T - t0`` in one-based period labels: one lag per post-offer step."""
    return max(int(n_periods - start - 1), 0)


def pre_offer_demean(d: np.ndarray, start: int) -> np.ndarray:
    """Subtract each unit's strictly pre-offer mean dose; the identity without pre-periods.

    The pre-offer mean is a function of pre-treatment data only, so it cannot
    leak offer-induced dose into the pre-period, and its constant lag sum is
    absorbed by the outcome's unit intercept.
    """
    d = np.asarray(d, dtype=np.float64)
    if start <= 0:
        return d.copy()
    return d - d[:, :start].mean(axis=1, keepdims=True)


def lag_matrix(d: np.ndarray, lag_horizon: int) -> np.ndarray:
    """Stack ``D_{i,t-l}`` for ``l = 0..L`` into ``(N, T, L+1)`` with zero pre-panel padding.

    Cells before the panel starts carry dose zero. That is exact when the
    behaviour cannot have occurred before study entry (post-randomization
    training); otherwise demean by the pre-offer baseline first
    (:func:`pre_offer_demean`) so the padded steady state is zero.
    """
    d = np.asarray(d, dtype=np.float64)
    if d.ndim != 2:
        raise ValueError("d must be an (N, T) dose panel.")
    if not isinstance(lag_horizon, (int, np.integer)) or isinstance(lag_horizon, bool) or lag_horizon < 0:
        raise ValueError("lag_horizon must be a nonnegative integer.")
    n, nt = d.shape
    out = np.zeros((n, nt, lag_horizon + 1), dtype=np.float64)
    for lag in range(lag_horizon + 1):
        if lag < nt:
            out[:, lag:, lag] = d[:, : nt - lag]
    return out


def offer_induced_dose(citt_d_mean: np.ndarray, z: np.ndarray, *, strata: Any = None,
                       center: bool = True) -> np.ndarray:
    """Tier T0's Robinson/FWL regressor ``(Z_i - e_i) * CITT_D(X_i, t) * 1{t >= t0}``.

    Propensity-centering the offer indicator ``Z_i`` by the sample or stratum
    offer rate ``e_i`` makes the regressor mean-zero given ``(X_i, t)`` and
    therefore orthogonal to any control-arm organic dose profile
    ``nu_Y(X_i) sum_l kappa_l E[D_{i,t-l}(0) | X_i, t]``, while identifying
    ``kappa_l`` and ``nu_Y(X_i)`` purely from randomized offer variation.
    Pass ``center=False`` for the uncentered increment ``Z_it * CITT_D(X_i, t)``.
    """
    citt = np.asarray(citt_d_mean, dtype=np.float64)
    z_arr = (np.asarray(z, dtype=np.float64) == 1).astype(np.float64)
    if not center:
        return citt * z_arr
    assigned = z_arr.any(axis=1).astype(np.float64)
    active = z_arr.any(axis=0).astype(np.float64)
    if strata is None:
        e_hat = np.full_like(assigned, assigned.mean())
    else:
        labels = _strata_labels(strata, len(assigned))
        e_hat = np.empty_like(assigned)
        for s in range(int(labels.max()) + 1):
            in_s = labels == s
            e_hat[in_s] = assigned[in_s].mean()
    return ((assigned - e_hat)[:, None] * active[None, :]) * citt


def _interact_terms(terms: list[tuple[str, np.ndarray]],
                    x_tilde: np.ndarray | None) -> tuple[np.ndarray, list[str]]:
    cols, names = [], []
    for name, term in terms:
        cols.append(term)
        names.append(name)
        if x_tilde is not None:
            for p in range(x_tilde.shape[1]):
                cols.append(term * x_tilde[:, p][:, None])
                names.append(f"{name}:x{p}")
    return np.stack(cols, axis=2), names


def organic_dose_columns(d0_hat: np.ndarray, lag_horizon: int,
                         x_tilde: np.ndarray | None) -> tuple[np.ndarray, list[str]]:
    """Tier T0's unconstrained nuisance controls ``[d0_it, d0^lag_it] (x) [1, x_tilde]``.

    ``d0_hat`` is the Stage-1 posterior-mean control-arm dose surface
    ``E[D_it(0) | X_i, t]`` (a function of ``(X_i, t)`` only, excluding the
    unit intercept and AR(1) innovations), and ``d0^lag_it = sum_{l=1}^{L} d0_{i,t-l}``
    with zero pre-panel padding. Entering them as unconstrained conjugate
    columns absorbs the growing-lag-window baseline profile of organic doses
    without constraining its slope to equal ``kappa_l``.
    """
    d0 = np.asarray(d0_hat, dtype=np.float64)
    lagged = lag_matrix(d0, lag_horizon)[:, :, 1:].sum(axis=2) if lag_horizon >= 1 else np.zeros_like(d0)
    return _interact_terms([("d0", d0), ("d0_lag", lagged)], x_tilde)


def pre_offer_columns(d: np.ndarray, start: int,
                      x_tilde: np.ndarray | None) -> tuple[np.ndarray, list[str]]:
    """Mundlak controls ``d_pre_i (x) [1, x_tilde]`` for the absorbed unit constant under ``demean_pre=True``."""
    d_arr = np.asarray(d, dtype=np.float64)
    n, nt = d_arr.shape
    if start <= 0:
        return np.zeros((n, nt, 0), dtype=np.float64), []
    d_pre = np.repeat(d_arr[:, :start].mean(axis=1, keepdims=True), nt, axis=1)
    return _interact_terms([("d_pre", d_pre)], x_tilde)


def control_function_columns(v_hat: np.ndarray, lag_horizon: int, x_tilde: np.ndarray | None,
                             start: int) -> tuple[np.ndarray, list[str]]:
    """The 3-term plug-in control-function basis of Tier T1, ``(N, T, 3 * (1 + P))``.

    Columns are ``[v_it, v^lag_it, vbar_i] (x) [1, x_tilde]`` with
    ``v^lag_it = sum_{l=1}^{L} v_{i,t-l}`` (zero pre-panel padding) and
    ``vbar_i`` the unit mean of the post-offer residuals. All three terms are
    mandatory: the contemporaneous residual deconfounds lag 0, the cumulative
    lagged residual deconfounds the level of the post-dose kernel (whose
    active window grows with calendar time even under static selection), and
    the unit mean absorbs static Mundlak selection. A two-term basis without
    ``v^lag`` forces ``kappa_0`` to absorb the whole multi-lag IV-vs-OLS gap.
    """
    v = np.asarray(v_hat, dtype=np.float64)
    n, nt = v.shape
    lagged = lag_matrix(v, lag_horizon)[:, :, 1:].sum(axis=2) if lag_horizon >= 1 else np.zeros_like(v)
    post = slice(start, nt)
    v_bar = np.repeat(v[:, post].mean(axis=1, keepdims=True), nt, axis=1)
    return _interact_terms([("v", v), ("v_lag", lagged), ("v_bar", v_bar)], x_tilde)


def standardized_covariates(x: np.ndarray) -> np.ndarray:
    """Centered, unit-scale covariates for the control-function interactions; constants drop."""
    x = np.asarray(x, dtype=np.float64)
    sd = x.std(axis=0)
    keep = sd > 0
    return (x[:, keep] - x[:, keep].mean(axis=0)) / sd[keep]


def compose_offer_effect(kappa: np.ndarray, citt_d: np.ndarray, start: int) -> np.ndarray:
    """Draw-by-draw composition ``CITT_Y^(b)(x,t) = sum_{l=0}^{t-t0} kappa_l^(b) CITT_D^(b)(x,t-l)``.

    ``kappa`` has shape ``(B, L+1)`` and already includes the moderator
    ``nu_Y^(b)(x)`` when the caller wants the conditional kernel (pass
    ``kappa[None] * nu[:, :, None]`` of shape ``(N, B, L+1)``); ``citt_d`` is
    ``(N, T, B)`` and must be zero before ``start``. Control-function terms are
    absent by construction: they have mean zero in both arms given ``(X, Z, t)``,
    so they cancel from the offer contrast.
    """
    citt_d = np.asarray(citt_d, dtype=np.float64)
    kappa = np.asarray(kappa, dtype=np.float64)
    n, nt, n_draws = citt_d.shape
    if kappa.ndim == 2:
        kappa = np.broadcast_to(kappa[None], (n,) + kappa.shape)
    if kappa.shape[:2] != (n, n_draws):
        raise ValueError("kappa must have shape (B, L+1) or (N, B, L+1).")
    out = np.zeros_like(citt_d)
    n_lags = kappa.shape[2]
    for lag in range(min(n_lags, nt - start)):
        # out[:, t, b] += kappa[:, b, lag] * citt_d[:, t - lag, b] for t >= start + lag.
        out[:, start + lag:, :] += kappa[:, None, :, lag] * citt_d[:, start: nt - lag, :]
    return out


# ---------------------------------------------------------------------------
# Toeplitz deconvolution and Anderson--Rubin inference
# ---------------------------------------------------------------------------

def toeplitz_lower(first_column: np.ndarray, size: int | None = None) -> np.ndarray:
    """Lower-triangular Toeplitz matrix with the given first column (zero-padded)."""
    c = np.asarray(first_column, dtype=np.float64)
    size = c.size if size is None else int(size)
    out = np.zeros((size, size))
    for k in range(min(c.size, size)):
        out[np.arange(k, size), np.arange(0, size - k)] = c[k]
    return out


def deconvolve(itt_y: np.ndarray, itt_d: np.ndarray, cov: np.ndarray, lag_horizon: int | None = None,
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unregularized Toeplitz deconvolution with the delta-method covariance.

    ``itt_y`` and ``itt_d`` are the ``K`` post-offer contrasts and ``cov`` is the
    ``(2K, 2K)`` joint covariance of ``(itt_y, itt_d)``. With the default
    ``lag_horizon = K - 1`` the system is square and ``kappa = T_d^{-1} y``;
    with fewer lags it is solved by least squares over all ``K`` equations.
    Returns ``(kappa, covariance, jacobian)`` where the Jacobian is with respect
    to the stacked ``(itt_y, itt_d)``.
    """
    y = np.asarray(itt_y, dtype=np.float64)
    d = np.asarray(itt_d, dtype=np.float64)
    k = y.size
    n_lags = k if lag_horizon is None else int(lag_horizon) + 1
    if not 1 <= n_lags <= k:
        raise ValueError("lag_horizon must lie between 0 and the number of post-offer periods minus one.")
    a = toeplitz_lower(d)[:, :n_lags]
    if n_lags == k:
        if d[0] == 0:
            raise ValueError("the first-stage contrast at the offer period is exactly zero; "
                             "the Toeplitz system is singular.")
        kappa = np.linalg.solve(a, y)
        resid = np.zeros(k)
        a_inv = np.linalg.inv(a)
        j_y = a_inv
        j_d = -a_inv @ toeplitz_lower(kappa, k)
    else:
        gram = a.T @ a
        if np.linalg.matrix_rank(gram) < n_lags:
            raise ValueError("the first-stage contrasts do not identify the requested kernel.")
        gram_inv = np.linalg.inv(gram)
        kappa = gram_inv @ a.T @ y
        resid = y - a @ kappa
        hankel = np.zeros((n_lags, k))
        for lag in range(n_lags):
            hankel[lag, : k - lag] = resid[lag:]
        j_y = gram_inv @ a.T
        j_d = gram_inv @ (hankel - a.T @ toeplitz_lower(kappa, k))
    jac = np.concatenate([j_y, j_d], axis=1)
    cov_kappa = jac @ cov @ jac.T if np.all(np.isfinite(cov)) else np.full((n_lags, n_lags), np.nan)
    return kappa, cov_kappa, jac


def anderson_rubin_statistic(kappa: np.ndarray, itt_y: np.ndarray, itt_d: np.ndarray, cov: np.ndarray,
                             ) -> float:
    """``m' Var(m)^{-1} m`` for ``m(kappa) = itt_y - T_kappa itt_d``, asymptotically ``chi^2_K`` under H0.

    The moment vector is linear in the two estimated trajectories, so its
    covariance follows exactly from the joint Neyman covariance for every
    hypothesized kernel; no division by the first stage is involved, which is
    what makes the set valid under a weak or vanishing first stage.
    """
    y = np.asarray(itt_y, dtype=np.float64)
    d = np.asarray(itt_d, dtype=np.float64)
    k = y.size
    t_kappa = toeplitz_lower(np.asarray(kappa, dtype=np.float64), k)
    a = np.concatenate([np.eye(k), -t_kappa], axis=1)
    m = y - t_kappa @ d
    v = a @ cov @ a.T
    try:
        return float(m @ np.linalg.solve(v, m))
    except np.linalg.LinAlgError:
        return float(m @ np.linalg.pinv(v) @ m)


# ---------------------------------------------------------------------------
# Stratified arm contrasts
# ---------------------------------------------------------------------------

def _stratified_contrast(values: np.ndarray, assigned: np.ndarray, strata: np.ndarray,
                         ) -> tuple[np.ndarray, np.ndarray, bool]:
    """Stratum-share-weighted difference in arm means of the unit vectors ``values`` ``(N, Q)``.

    Returns ``(estimate, neyman_covariance, inference_available)``; the
    covariance is ``sum_s w_s^2 (S_s1 / n_s1 + S_s0 / n_s0)`` with within-arm
    sample covariance matrices ``S``, which is conservative for the
    randomization variance by the finite-population treatment-effect
    heterogeneity term.
    """
    n, q = values.shape
    estimate = np.zeros(q)
    cov = np.zeros((q, q))
    available = True
    for s in np.unique(strata):
        in_s = strata == s
        w = in_s.mean()
        a = values[in_s & assigned]
        b = values[in_s & ~assigned]
        if a.shape[0] == 0 or b.shape[0] == 0:
            raise ValueError("every stratum must contain units in both arms.")
        estimate += w * (a.mean(axis=0) - b.mean(axis=0))
        if min(a.shape[0], b.shape[0]) >= 2:
            cov += w * w * (np.atleast_2d(np.cov(a, rowvar=False, ddof=1)) / a.shape[0]
                            + np.atleast_2d(np.cov(b, rowvar=False, ddof=1)) / b.shape[0])
        else:
            available = False
    if not available:
        cov = np.full((q, q), np.nan)
    return estimate, cov, available


def _weighted_ecdf(values: np.ndarray, weights: np.ndarray, grid: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="stable")
    v, w = values[order], weights[order]
    cum = np.cumsum(w) / w.sum()
    idx = np.searchsorted(v, grid, side="right")
    return np.where(idx > 0, cum[np.maximum(idx - 1, 0)], 0.0)


def _cluster_ols(y: np.ndarray, design: np.ndarray, clusters: np.ndarray, coefficient: int,
                 ) -> tuple[float, float, int]:
    """OLS coefficient with a unit-clustered sandwich SE (CR1 small-sample factor)."""
    n, k = design.shape
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ coef
    gram_inv = np.linalg.pinv(design.T @ design)
    meat = np.zeros((k, k))
    labels = np.unique(clusters)
    for c in labels:
        rows = clusters == c
        score = design[rows].T @ resid[rows]
        meat += np.outer(score, score)
    g = labels.size
    rank = int(np.linalg.matrix_rank(design))
    factor = (g / max(g - 1, 1)) * ((n - 1) / max(n - rank, 1))
    cov = factor * gram_inv @ meat @ gram_inv
    se = float(np.sqrt(max(cov[coefficient, coefficient], 0.0)))
    return float(coef[coefficient]), se, g


# ---------------------------------------------------------------------------
# The design-based reference
# ---------------------------------------------------------------------------

@dataclass
class IntensityReference:
    """Design-based (Tier T0) estimates and diagnostics for a continuous-intensity encouragement panel.

    Attributes
    ----------
    itt
        One row per calendar period: arm contrasts ``itt_d``, ``itt_y`` with
        SEs and pointwise intervals, the per-period Wald ratio ``omega`` with
        its Fieller/Anderson--Rubin set, and, for post-offer periods, the
        cumulative return per unit of cumulative dose ``omega_cum`` with its
        exact set (``omega_cum_set_type``, ``omega_cum_lower_1`` ...).
    covariance
        Joint Neyman covariance of ``(itt_y, itt_d)`` over the post-offer
        periods, indexed by ``(quantity, period_index)``.
    kernel
        One row per lag ``l = 0..L``: the Toeplitz kernel ``kappa`` with
        delta-method SE and interval, and the cumulative kernel ``kappa_cum``.
    kernel_covariance
        Delta-method covariance of ``kappa``, shape ``(L+1, L+1)``.
    condition_number
        ``cond(T_{ITT_D})`` in the 2-norm: near 1 for a front-loaded compliance
        pulse, exponentially large for a gradual ramp.
    ecdf, fosd
        Arm-wise (stratum-weighted) ECDFs of cumulative post-offer dose on the
        pooled support and the first-order-stochastic-dominance check implied
        by monotone compliance (``F_1 <= F_0`` everywhere).
    dip
        The dose-history-conditioned pre-treatment dip screen: the coefficient
        of ``Delta Y_{i,t-1}`` in a regression of ``Delta D_it`` on it, lagged
        doses ``D_{i,t-1..t-L-1}``, covariates, the offer and period fixed
        effects, with a unit-clustered SE.
    metadata
        Design, dimensions, horizon and assumptions.
    """

    itt: pd.DataFrame
    covariance: pd.DataFrame
    kernel: pd.DataFrame
    kernel_covariance: np.ndarray
    condition_number: float
    ecdf: pd.DataFrame
    fosd: dict[str, Any]
    dip: dict[str, Any]
    metadata: dict[str, Any]
    _post: dict[str, np.ndarray] = field(default_factory=dict, repr=False)

    # -- joint Anderson--Rubin region ------------------------------------------
    def anderson_rubin(self, kappa: Any) -> dict[str, Any]:
        """Joint AR test of ``H0: kernel = kappa`` (``chi^2`` with one degree of freedom per post period)."""
        k = int(self.metadata["n_post_periods"])
        kappa = np.asarray(kappa, dtype=np.float64).ravel()
        if kappa.size > k:
            raise ValueError(f"kappa may have at most {k} lags.")
        statistic = anderson_rubin_statistic(kappa, self._post["itt_y"], self._post["itt_d"],
                                             self._post["covariance"])
        critical = float(stats.chi2.ppf(1 - self.metadata["alpha"], df=k))
        return {"statistic": statistic, "df": k, "pvalue": float(stats.chi2.sf(statistic, df=k)),
                "critical": critical, "in_region": bool(statistic <= critical)}

    def anderson_rubin_projection(self, *, n_candidates: int = 4000, inflation: float = 3.0,
                                  random_seed: int = 0) -> pd.DataFrame:
        """Monte Carlo projection of the joint AR region onto each lag and cumulative lag.

        Candidates are drawn from the delta-method normal approximation inflated
        by ``inflation`` (plus axis-aligned rays); the per-coordinate range of
        the accepted candidates is an **inner** approximation of the exact
        projection, so it understates an unbounded set. ``hits_envelope`` flags
        coordinates whose accepted range reaches the sampling envelope, in
        which case the region should be regarded as unbounded in that direction.
        """
        if not np.isfinite(self.kernel_covariance).all():
            raise ValueError("the reference covariance is unavailable (fewer than two units per stratum-arm).")
        kappa_hat = self.kernel["kappa"].to_numpy()
        n_lags = kappa_hat.size
        rng = np.random.default_rng(random_seed)
        try:
            chol = np.linalg.cholesky(self.kernel_covariance + 1e-12 * np.eye(n_lags))
        except np.linalg.LinAlgError:
            chol = np.diag(np.sqrt(np.maximum(np.diag(self.kernel_covariance), 1e-12)))
        normal = rng.normal(size=(int(n_candidates), n_lags))
        rays = np.concatenate([np.eye(n_lags), -np.eye(n_lags)]) * inflation * 2.0
        offsets = np.concatenate([normal * inflation, rays])
        candidates = kappa_hat + offsets @ chol.T
        candidates = np.concatenate([kappa_hat[None], candidates])
        stat = np.array([anderson_rubin_statistic(c, self._post["itt_y"], self._post["itt_d"],
                                                  self._post["covariance"]) for c in candidates])
        critical = float(stats.chi2.ppf(1 - self.metadata["alpha"], df=int(self.metadata["n_post_periods"])))
        accepted = candidates[stat <= critical]
        cum = np.cumsum(accepted, axis=1)
        env_lo = candidates.min(axis=0)
        env_hi = candidates.max(axis=0)
        rows = []
        for lag in range(n_lags):
            rows.append(dict(lag=lag, kappa_lower=float(accepted[:, lag].min()),
                             kappa_upper=float(accepted[:, lag].max()),
                             kappa_cum_lower=float(cum[:, lag].min()),
                             kappa_cum_upper=float(cum[:, lag].max()),
                             hits_envelope=bool(np.isclose(accepted[:, lag].min(), env_lo[lag])
                                                or np.isclose(accepted[:, lag].max(), env_hi[lag])),
                             n_accepted=int(accepted.shape[0]), n_candidates=int(candidates.shape[0])))
        return pd.DataFrame(rows)

    # -- discrepancy statistics ------------------------------------------------
    def discrepancy(self, kappa_draws: Any, *, label: str = "model", cumulative: bool = False,
                    ) -> pd.DataFrame:
        """``z_{model - T0}``: posterior cumulative kernel against the design-based cumulative kernel.

        ``kappa_draws`` holds posterior draws with lags on the last axis,
        ``(B, L+1)``; pass ``cumulative=True`` when they are already cumulated.
        The z statistic divides the difference at each horizon ``h`` by
        ``sqrt(SE_T0(h)^2 + SD_post(h)^2)``; the two estimators share the data
        and are positively dependent, so this is a conservative (alarm-averse)
        standardization. The row at the last horizon is the headline
        ``z_{T2-T0}`` / ``z_{T1-T0}`` statistic.
        """
        draws = np.asarray(kappa_draws, dtype=np.float64)
        if draws.ndim != 2:
            raise ValueError("kappa_draws must be a (draws, lags) array.")
        n_lags = int(len(self.kernel))
        if draws.shape[1] < n_lags:
            raise ValueError(f"kappa_draws must cover at least {n_lags} lags.")
        cum = draws[:, :n_lags] if cumulative else np.cumsum(draws[:, :n_lags], axis=1)
        critical = _critical(self.metadata["alpha"])
        rows = []
        for h in range(n_lags):
            model_mean = float(cum[:, h].mean())
            model_sd = float(cum[:, h].std(ddof=1)) if cum.shape[0] >= 2 else np.nan
            ref = float(self.kernel.loc[h, "kappa_cum"])
            ref_se = float(self.kernel.loc[h, "kappa_cum_se"])
            diff = model_mean - ref
            se = float(np.sqrt(ref_se ** 2 + model_sd ** 2))
            z = diff / se if np.isfinite(se) and se > 0 else np.nan
            rows.append(dict(label=label, horizon=h, model_mean=model_mean, model_sd=model_sd,
                             reference=ref, reference_se=ref_se, difference=diff, se=se, z_stat=z,
                             flag=bool(abs(z) > critical) if np.isfinite(z) else None))
        return pd.DataFrame(rows)


def intensity_reference(
    y: Any, d: Any, z: Any, t: Any = None, *, x: Any = None, strata: Any = None,
    lag_horizon: int | None = None, demean_pre: bool = False, alpha: float = 0.05,
) -> IntensityReference:
    """Design-based (Tier T0) analysis of a continuous-intensity encouragement panel without a model.

    Parameters
    ----------
    y, d, z
        Outcome, dose and offer panels, ``(N, T)``; see :func:`validate_intensity`.
    t
        Calendar vector with unit gaps; defaults to ``1..T``.
    x
        Optional baseline covariates ``(N, P)``; used only by the dip screen.
    strata
        Optional length-``N`` baseline randomization strata (sites). Contrasts
        are stratum-share-weighted differences in arm means, so a site-blocked
        design is analysed as it was randomized.
    lag_horizon
        Number of lags beyond zero in the kernel; defaults to ``T - t0`` (one
        unknown per post-offer equation, the exactly identified square system).
    demean_pre
        Demean the dose by each unit's strictly pre-offer mean before building
        lags (for behaviours with a pre-study steady state).
    alpha
        Tail probability for intervals and confidence sets.

    Returns
    -------
    IntensityReference
        See the class documentation. All intervals are pointwise, asymptotic
        normal/chi-square approximations under individual-unit (stratified)
        complete randomization; none is a test of exclusion or linearity.
    """
    panel = validate_intensity(z, d, t)
    critical = _critical(alpha)
    n, nt = panel.z.shape
    y = _finite_outcome(y, panel.z.shape)
    strata_idx = _strata_labels(strata, n)
    x_arr = _covariates(x, n)
    dose = pre_offer_demean(panel.d, panel.start) if demean_pre else panel.d
    k = panel.n_post
    n_lags = default_lag_horizon(nt, panel.start) + 1 if lag_horizon is None else int(lag_horizon) + 1
    if not 1 <= n_lags <= k:
        raise ValueError(f"lag_horizon must lie between 0 and {k - 1} for this panel.")

    # Arm contrasts for every period, with the joint covariance of (y, d).
    stacked = np.concatenate([y, dose], axis=1)                    # (N, 2T): y then d
    estimate, cov, available = _stratified_contrast(stacked, panel.assigned, strata_idx)
    itt_y_all, itt_d_all = estimate[:nt], estimate[nt:]
    post = np.arange(panel.start, nt)
    idx_post = np.concatenate([post, nt + post])
    cov_post = cov[np.ix_(idx_post, idx_post)]
    itt_y, itt_d = itt_y_all[post], itt_d_all[post]

    rows = []
    for j in range(nt):
        vy, vd, cyd = cov[j, j], cov[nt + j, nt + j], cov[j, nt + j]
        row = {
            "period": float(panel.t[j]), "period_index": j, "horizon": int(panel.exposure[j]),
            "post_encouragement": bool(j >= panel.start),
            "n_encouraged": int(panel.assigned.sum()), "n_control": int((~panel.assigned).sum()),
            "itt_y_d_cov": float(cyd), "inference": "normal_ar",
        }
        for name, delta, variance in (("itt_y", itt_y_all[j], vy), ("itt_d", itt_d_all[j], vd)):
            se = float(np.sqrt(variance)) if np.isfinite(variance) else np.nan
            row.update({name: float(delta), f"{name}_se": se, f"{name}_lower": delta - critical * se,
                        f"{name}_upper": delta + critical * se})
        row["omega"] = float(itt_y_all[j] / itt_d_all[j]) if itt_d_all[j] != 0 else np.nan
        row.update({key.replace("wald", "omega"): value for key, value in
                    _wald_set(itt_y_all[j], itt_d_all[j], vy, vd, cyd, critical).items()})
        if j >= panel.start:
            c = np.zeros(2 * k)
            m = j - panel.start + 1
            c[:m] = 1.0                                   # cumulative y
            cy = c @ np.concatenate([itt_y, np.zeros(k)])
            cd = itt_d[:m].sum()
            sel_y = np.zeros(2 * k); sel_y[:m] = 1.0
            sel_d = np.zeros(2 * k); sel_d[k:k + m] = 1.0
            v_cy = float(sel_y @ cov_post @ sel_y)
            v_cd = float(sel_d @ cov_post @ sel_d)
            c_cyd = float(sel_y @ cov_post @ sel_d)
            row["cumulative_itt_y"] = float(cy)
            row["cumulative_itt_d"] = float(cd)
            row["cumulative_itt_y_se"] = float(np.sqrt(v_cy)) if np.isfinite(v_cy) else np.nan
            row["cumulative_itt_d_se"] = float(np.sqrt(v_cd)) if np.isfinite(v_cd) else np.nan
            row["omega_cum"] = float(cy / cd) if cd != 0 else np.nan
            if cd != 0 and np.isfinite([v_cy, v_cd, c_cyd]).all():
                grad = np.array([1.0 / cd, -cy / cd ** 2])
                v_omega = grad @ np.array([[v_cy, c_cyd], [c_cyd, v_cd]]) @ grad
                row["omega_cum_se_delta"] = float(np.sqrt(max(v_omega, 0.0)))
            else:
                row["omega_cum_se_delta"] = np.nan
            row.update({key.replace("wald", "omega_cum"): value for key, value in
                        _wald_set(cy, cd, v_cy, v_cd, c_cyd, critical).items()})
        else:
            for name in ("cumulative_itt_y", "cumulative_itt_d", "cumulative_itt_y_se",
                         "cumulative_itt_d_se", "omega_cum", "omega_cum_se_delta",
                         "omega_cum_lower_1", "omega_cum_upper_1", "omega_cum_lower_2", "omega_cum_upper_2"):
                row[name] = np.nan
            row["omega_cum_set_type"] = "pre_encouragement"
            row["omega_cum_reason"] = ""
        rows.append(row)
    itt_frame = pd.DataFrame(rows)

    # Toeplitz deconvolution.
    kappa, cov_kappa, _ = deconvolve(itt_y, itt_d, cov_post, n_lags - 1)
    t_d = toeplitz_lower(itt_d)[:, :n_lags]
    condition = float(np.linalg.cond(t_d))
    cum_mat = np.tril(np.ones((n_lags, n_lags)))
    kappa_cum = cum_mat @ kappa
    cov_cum = cum_mat @ cov_kappa @ cum_mat.T
    se_k = np.sqrt(np.maximum(np.diag(cov_kappa), 0.0)) if np.isfinite(cov_kappa).all() else np.full(n_lags, np.nan)
    se_c = np.sqrt(np.maximum(np.diag(cov_cum), 0.0)) if np.isfinite(cov_cum).all() else np.full(n_lags, np.nan)
    kernel = pd.DataFrame({
        "lag": np.arange(n_lags), "kappa": kappa, "kappa_se": se_k,
        "kappa_lower": kappa - critical * se_k, "kappa_upper": kappa + critical * se_k,
        "kappa_cum": kappa_cum, "kappa_cum_se": se_c,
        "kappa_cum_lower": kappa_cum - critical * se_c, "kappa_cum_upper": kappa_cum + critical * se_c,
    })

    # Arm-wise ECDF shift of cumulative post-offer dose and the FOSD check.
    cum_dose = dose[:, post].sum(axis=1)
    weights = np.zeros(n)
    for s in np.unique(strata_idx):
        in_s = strata_idx == s
        for arm in (True, False):
            cell = in_s & (panel.assigned == arm)
            weights[cell] = in_s.mean() / cell.sum()
    grid = np.unique(cum_dose)
    f1 = _weighted_ecdf(cum_dose[panel.assigned], weights[panel.assigned], grid)
    f0 = _weighted_ecdf(cum_dose[~panel.assigned], weights[~panel.assigned], grid)
    ecdf = pd.DataFrame({"dose": grid, "F_control": f0, "F_encouraged": f1, "difference": f0 - f1})
    n1, n0 = int(panel.assigned.sum()), int((~panel.assigned).sum())
    sup_violation = float(np.max(f1 - f0))
    scale = np.sqrt(n1 * n0 / (n1 + n0))
    fosd_stat = max(sup_violation, 0.0) * scale
    fosd = {
        "max_control_minus_encouraged": float(np.max(f0 - f1)),
        "sup_violation": max(sup_violation, 0.0),
        "statistic": float(fosd_stat),
        "pvalue": float(np.exp(-2.0 * fosd_stat ** 2)),
        "holds_at_level": bool(np.exp(-2.0 * fosd_stat ** 2) > alpha),
        "extensive_margin_shift": float((f0 - f1)[grid <= 0].max()) if np.any(grid <= 0) else np.nan,
        "mean_shift": float(itt_frame.loc[itt_frame.post_encouragement, "itt_d"].sum()),
        "method": "one_sided_smirnov_asymptotic",
    }

    # Dose-history-conditioned dip screen.
    dip = _dip_screen(y, dose, panel, x_arr, strata_idx, n_lags - 1)

    index = pd.MultiIndex.from_arrays(
        [np.repeat(["itt_y", "itt_d"], k), np.concatenate([post, post])], names=["quantity", "period_index"])
    metadata = dict(panel.metadata(), estimator="stratified_difference_in_means",
                    covariance="stratified_neyman", n_strata=int(np.unique(strata_idx).size),
                    inference_available=bool(available), n_post_periods=int(k), lag_horizon=int(n_lags - 1),
                    demean_pre=bool(demean_pre), alpha=float(alpha), identification="deconvolution",
                    kernel_estimand="compliance_weighted_population_kernel",
                    assumptions=[
                        "Individual-unit (stratified) complete randomization of the offer.",
                        "Exclusion and no anticipation of the offer or the dose.",
                        "Linear, additive, calendar-homogeneous lag effects of the dose on the outcome.",
                        "Time-invariant conditional complier weights across post-offer periods.",
                    ],
                    limitations=[
                        "Intervals are pointwise asymptotic approximations, not simultaneous bands.",
                        "The unregularized Toeplitz inverse amplifies noise when compliance ramps up gradually.",
                    ])
    return IntensityReference(
        itt=itt_frame, covariance=pd.DataFrame(cov_post, index=index, columns=index), kernel=kernel,
        kernel_covariance=cov_kappa, condition_number=condition, ecdf=ecdf, fosd=fosd, dip=dip,
        metadata=metadata, _post={"itt_y": itt_y, "itt_d": itt_d, "covariance": cov_post})


def _dip_screen(y: np.ndarray, dose: np.ndarray, panel: IntensityPanel, x: np.ndarray | None,
                strata: np.ndarray, lag_horizon: int) -> dict[str, Any]:
    """Regress ``Delta D_it`` on ``Delta Y_{i,t-1}`` given dose history, covariates, offer and period FE.

    ``Delta Y_{i,t-1} = Y_{i,t-1} - Y_{i,t-2}`` is moved by doses back to
    ``D_{i,t-(L+2)}`` under the lag kernel, so the regression conditions on
    ``D_{i,t-1}, ..., D_{i,t-L-2}`` (zero before the panel), which blocks the
    mechanical path from past doses to the lagged outcome innovation.
    """
    n, nt = dose.shape
    n_history = lag_horizon + 2
    if nt < 3:
        return {"beta_dip": np.nan, "se": np.nan, "z_stat": np.nan, "pvalue": np.nan, "n_obs": 0,
                "n_units": 0, "periods": [], "method": "unavailable_fewer_than_three_periods"}
    lags = lag_matrix(dose, n_history)                            # D_{i,t-l}, l = 0..L+2
    periods = np.arange(2, nt)
    rows_y, rows_x, rows_c = [], [], []
    for j in periods:
        delta_d = dose[:, j] - dose[:, j - 1]
        delta_y_prev = y[:, j - 1] - y[:, j - 2]
        cols = [delta_y_prev] + [lags[:, j, lag] for lag in range(1, n_history + 1)]
        cols.append(panel.z[:, j])
        if x is not None:
            cols.extend(x.T)
        period_dummies = (periods[None, :] == j).astype(np.float64).repeat(n, axis=0)
        block = np.column_stack(cols + [period_dummies])
        rows_y.append(delta_d)
        rows_x.append(block)
        rows_c.append(np.arange(n))
    design = np.concatenate(rows_x)
    uniq = np.unique(strata)
    if uniq.size > 1:
        s_dummies = (np.tile(strata, periods.size)[:, None] == uniq[None, 1:]).astype(np.float64)
        design = np.concatenate([design, s_dummies], axis=1)
    response = np.concatenate(rows_y)
    clusters = np.concatenate(rows_c)
    beta, se, n_clusters = _cluster_ols(response, design, clusters, 0)
    z_stat = beta / se if se > 0 else np.nan
    return {"beta_dip": beta, "se": se, "z_stat": float(z_stat),
            "pvalue": float(2 * stats.norm.sf(abs(z_stat))) if np.isfinite(z_stat) else np.nan,
            "n_obs": int(response.size), "n_units": int(n_clusters),
            "periods": [float(panel.t[j]) for j in periods],
            "controls": f"D_{{t-1..t-{n_history}}}, Z_it, period FE"
                        + (", X" if x is not None else "") + (", strata FE" if uniq.size > 1 else ""),
            "method": "ols_unit_clustered_cr1"}


# ---------------------------------------------------------------------------
# Stage-2 backends
# ---------------------------------------------------------------------------
#
# Stage 2 regresses the outcome on the lag matrix (plus, for Tier T1, the fixed
# control-function columns). Two backends share one small protocol so the
# encouragement layer can be exercised end to end without the sampler:
#
# * ``"longbet"`` -- the LongBet sampler with the distributed-lag kernel
#   (``exposure_lags``/``extra_cols`` in ``LongBet.fit``): forests for
#   ``mu_Y(x)`` and ``nu_Y(x)``, the horseshoe linear moderator term, the GP
#   prior over lags, the unit intercept, AR(1) errors and the conjugate trend
#   block. This is the model of the paper.
# * ``"linear"`` -- a two-way fixed-effects Gaussian distributed-lag regression
#   with a flat prior, whose posterior is available in closed form. It has no
#   covariate moderation (``nu_Y(x) = 1``) and independent errors; it exists
#   for fast tests, smoke runs and as a transparent fallback.

STAGE2_BACKENDS: tuple[str, ...] = ("longbet", "linear")


def _flatten_chain_major(a: np.ndarray, num_chains: int) -> np.ndarray:
    a = np.asarray(a)
    if a.ndim > 1 and num_chains > 1:
        return a.reshape(-1, *a.shape[2:])
    return a


class LinearLagStage2:
    """Closed-form two-way fixed-effects distributed-lag regression (the ``"linear"`` backend).

    ``y_it = alpha_t + gamma_i + d_it' kappa + c_it' rho + e_it`` with
    ``e_it ~ N(0, sigma^2)`` i.i.d., unit and period effects swept out by
    within transformations, a flat prior on ``(kappa, rho)`` and the Jeffreys
    prior on ``sigma^2``, so the posterior is Normal-Inverse-Gamma and draws are
    exact. The moderator is identically one.
    """

    backend = "linear"

    def __init__(self) -> None:
        self.kappa_draws: np.ndarray | None = None
        self.extra_coef_draws: np.ndarray | None = None
        self.sigma2_draws: np.ndarray | None = None
        self.coef_mean: np.ndarray | None = None
        self.coef_cov_scale: np.ndarray | None = None
        self.n_lags: int = 0
        self.extra_names: list[str] = []
        self.num_chains: int = 1

    def fit(self, y: np.ndarray, lag_mat: np.ndarray, extra_cols: np.ndarray | None, *,
            num_draws: int, num_chains: int, random_seed: int,
            extra_names: list[str] | None = None) -> LinearLagStage2:
        y = np.asarray(y, dtype=np.float64)
        n, nt, n_lags = lag_mat.shape
        cols = [lag_mat.reshape(n * nt, n_lags)]
        if extra_cols is not None and extra_cols.shape[2] > 0:
            cols.append(np.asarray(extra_cols, dtype=np.float64).reshape(n * nt, -1))
        design = np.concatenate(cols, axis=1)
        resp = y.reshape(n * nt)
        unit = np.repeat(np.arange(n), nt)
        period = np.tile(np.arange(nt), n)

        def within(a: np.ndarray) -> np.ndarray:
            # Two-way demeaning of a balanced panel: unit means, then period means.
            a = a - np.asarray(pd.DataFrame(a).groupby(unit).transform("mean"))
            return a - np.asarray(pd.DataFrame(a).groupby(period).transform("mean"))

        design_w = within(design)
        resp_w = within(resp[:, None])[:, 0]
        keep = design_w.std(axis=0) > 1e-12          # time-invariant columns are absorbed
        gram = design_w[:, keep].T @ design_w[:, keep]
        gram_inv = np.linalg.pinv(gram)
        coef_k = gram_inv @ design_w[:, keep].T @ resp_w
        resid = resp_w - design_w[:, keep] @ coef_k
        dof = max(n * nt - n - nt + 1 - int(keep.sum()), 1)
        rss = float(resid @ resid)
        rng = np.random.default_rng(random_seed)
        total = int(num_draws) * int(num_chains)
        sigma2 = rss / rng.chisquare(dof, size=total)
        chol = np.linalg.cholesky(gram_inv + 1e-12 * np.eye(gram_inv.shape[0]))
        normal = rng.normal(size=(total, gram_inv.shape[0]))
        coef = coef_k[None, :] + np.sqrt(sigma2)[:, None] * (normal @ chol.T)
        full = np.zeros((total, design.shape[1]))
        full[:, keep] = coef
        self.kappa_draws = full[:, :n_lags]
        self.extra_coef_draws = full[:, n_lags:] if design.shape[1] > n_lags else None
        self.sigma2_draws = sigma2
        full_mean = np.zeros(design.shape[1]); full_mean[keep] = coef_k
        self.coef_mean = full_mean
        self.coef_cov_scale = gram_inv
        self.n_lags = n_lags
        self.extra_names = list(extra_names or [])
        self.num_chains = int(num_chains)
        self.absorbed_columns = [i for i, k in enumerate(keep) if not k]
        return self

    def nu_draws(self, x: np.ndarray, x_trt: np.ndarray | None = None) -> np.ndarray:
        del x_trt
        return np.ones((self.kappa_draws.shape[0], len(x)))

    def to_arrays(self) -> dict[str, np.ndarray]:
        out = {"linear_kappa": self.kappa_draws, "linear_sigma2": self.sigma2_draws,
               "linear_num_chains": np.asarray(self.num_chains)}
        if self.extra_coef_draws is not None:
            out["linear_extra_coef"] = self.extra_coef_draws
        return out

    @classmethod
    def from_arrays(cls, arrays: dict[str, np.ndarray], extra_names: list[str]) -> LinearLagStage2:
        obj = cls()
        obj.kappa_draws = np.asarray(arrays["linear_kappa"], dtype=np.float64)
        obj.sigma2_draws = np.asarray(arrays["linear_sigma2"], dtype=np.float64)
        obj.extra_coef_draws = (np.asarray(arrays["linear_extra_coef"], dtype=np.float64)
                                if "linear_extra_coef" in arrays else None)
        obj.num_chains = int(arrays["linear_num_chains"])
        obj.n_lags = obj.kappa_draws.shape[1]
        obj.extra_names = list(extra_names)
        return obj


class LongBetLagStage2:
    """Adapter over :class:`longbet.LongBet` fitted with the distributed-lag kernel (the ``"longbet"`` backend).

    Thin on purpose: it normalizes what the composition needs -- kernel draws on
    the outcome's original scale ``(B, L+1)``, the extra-column coefficients
    ``(B, Q)`` and moderator draws ``nu_Y(x)`` ``(B, N)`` -- so that the
    encouragement layer never reaches into the sampler's trace itself.
    """

    backend = "longbet"

    def __init__(self, model: Any, *, keep_extra: np.ndarray | None = None) -> None:
        self.model = model
        self.keep_extra = None if keep_extra is None else np.asarray(keep_extra, dtype=bool)
        self._extra_coef_override: np.ndarray | None = None

    # The sampler owns these names; keep every access in one place.
    @property
    def num_chains(self) -> int:
        return int(self.model.config.num_chains)

    @property
    def kappa_draws(self) -> np.ndarray:
        draws = getattr(self.model, "kappa_draws_", None)
        if draws is None:
            trace = self.model.trace
            raw = getattr(trace, "kappa", None)
            if raw is None:
                raise RuntimeError("the fitted LongBet model carries no lag-kernel draws (kappa).")
            draws = _flatten_chain_major(np.asarray(raw), self.num_chains) * float(self.model.sdy)
        return np.asarray(draws, dtype=np.float64)

    @property
    def extra_coef_draws(self) -> np.ndarray | None:
        if self._extra_coef_override is not None:
            return self._extra_coef_override
        draws = getattr(self.model, "extra_coef_draws_", None)
        if draws is None:
            raw = getattr(self.model.trace, "extra_coef", None)
            if raw is not None:
                draws = _flatten_chain_major(np.asarray(raw), self.num_chains) * float(self.model.sdy)
        if self.keep_extra is None:
            return None if draws is None else np.asarray(draws, dtype=np.float64)
        n_draws = self.kappa_draws.shape[0]
        full = np.zeros((n_draws, self.keep_extra.size), dtype=np.float64)
        if draws is not None and self.keep_extra.any():
            full[:, self.keep_extra] = np.asarray(draws, dtype=np.float64)
        return full

    @extra_coef_draws.setter
    def extra_coef_draws(self, value: np.ndarray | None) -> None:
        self._extra_coef_override = None if value is None else np.asarray(value, dtype=np.float64)

    def nu_draws(self, x: np.ndarray, x_trt: np.ndarray | None = None) -> np.ndarray:
        """Moderator draws ``nu_Y(x)``, ``(B, N)``: trees plus the horseshoe linear term, anchored."""
        fn = getattr(self.model, "moderator_draws", None)
        if fn is None:
            raise RuntimeError("the fitted LongBet model does not expose moderator draws.")
        return np.asarray(fn(np.asarray(x, dtype=np.float32),
                             None if x_trt is None else np.asarray(x_trt, dtype=np.float32)),
                          dtype=np.float64)


# ---------------------------------------------------------------------------
# Posterior containers
# ---------------------------------------------------------------------------

def _summaries(draws: np.ndarray, alpha: float, axis: int = -1) -> dict[str, np.ndarray]:
    with np.errstate(invalid="ignore"):
        return {
            "mean": np.nanmean(draws, axis=axis), "median": np.nanmedian(draws, axis=axis),
            "lower": np.nanquantile(draws, alpha / 2, axis=axis),
            "upper": np.nanquantile(draws, 1 - alpha / 2, axis=axis),
            "sd": np.nanstd(draws, axis=axis, ddof=1) if draws.shape[axis] > 1 else np.zeros(
                np.delete(draws.shape, axis)),
        }


def cumulative_return_draws(citt_y: np.ndarray, citt_d: np.ndarray, *, floor: float = 0.0) -> np.ndarray:
    """``omega_cum(t) = sum_{r<=t} CITT_Y / sum_{r<=t} CITT_D`` along the horizon axis; exact zeros give NaN."""
    horizon_axis = 1 if np.ndim(citt_y) >= 3 else -2
    num = np.cumsum(citt_y, axis=horizon_axis)
    den = np.cumsum(citt_d, axis=horizon_axis)
    if floor > 0:
        den = np.where(np.abs(den) < floor, np.sign(den) * floor + (den == 0) * floor, den)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(num, den, out=np.full_like(num, np.nan), where=den != 0)
