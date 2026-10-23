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

"""Event-time windows, the common exposure basis and the two-pass projection of masked Rank-K panels.

A *masked* panel observes unit ``i`` on a contiguous event-time prefix
``O_i = {1, ..., L_i}`` of the post-adoption horizon ``s = 1, ..., H``
(administrative censoring); the cells ``s > L_i`` are ``NaN`` in ``y`` and leave
the likelihood. All units share one orthonormal exposure basis
``Psi (H, K)`` on ``s = 1..H`` (:func:`longbet.orthonormal_exposure_basis`), so
mode 1 is the time-averaged effect over the *full* horizon for every unit and
the higher modes are comparable across units whatever their window.

Per-unit mode coordinates are estimated by the two-pass ("mode-coupling")
projection of the methodology plan:

1. higher modes ``k >= 2`` on units with ``L_i >= min_support_modes[k]``:
   ``theta_hat_ik = Cov_{O_i}(psi_k, r_i) / Var_{O_i}(psi_k)`` -- the per-unit OLS
   slope of the observed residual on ``{1, psi_k}`` (on a partial window the
   constant and ``psi_k`` are not orthogonal, so the level must be partialled
   out; on a full window this is the paper's ``<psi_k, r_i> / H``) -- with
   precision weight proportional to ``sum_{O_i} (psi_k - mean psi_k)^2`` (units
   below support get weight 0, so their ``tau_k(x_i)`` comes from the
   cross-sectional posterior only);
2. mode 1 on every unit: the observed-window mean of the residual after the
   fitted treatment dynamics are removed,
   ``theta_hat_i1 = (1 / L_i) sum_{s in O_i} [r_i(s) - z_i sum_{k >= 2} tau~_ik psi_k(s)]``,
   with precision weight proportional to ``L_i``. Because
   ``sum_{s = 1}^{H} psi_k(s) = 0`` the correction vanishes on fully observed
   windows, where the projection reduces *exactly* (in floating point as well)
   to the paper's per-cohort projection.

The alternative ``projection="restricted_wls"`` fits every unit's supported
modes jointly, ``(Psi_O' Psi_O)^{-1} Psi_O' r_O``, with precision weights from the
restricted Gram matrix and a condition-number guard.

Everything here is plain NumPy and independent of the samplers; the arithmetic on
fully observed windows mirrors :func:`longbet._rank_k._project_panel_to_mode` and
:func:`longbet._rank_k._compute_mode_precision_weights` operation by operation so
that complete panels reproduce the released estimator bit for bit.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class EventTimeLayout:
    """Observation structure of a (possibly masked) panel on a common event-time basis.

    Attributes:
      t0: Pre-period length (periods ``0 .. t0 - 1`` are pre-adoption for every unit).
      horizon: Event-time horizon ``H`` of the common basis.
      n_periods: Number of panel periods ``T``.
      cohorts: ``(N,)`` first-adoption period of every unit (``T`` for never-treated).
      ever: ``(N,)`` ever-treated indicator.
      obs: ``(N, T)`` observed-cell indicator.
      pre_count: ``(N,)`` number of observed pre-period cells.
      waves: Sorted first-adoption periods of the treated units.
      wave_weights: Wave ``g`` -> share of the observed treated cells held by wave ``g``.
      wave_horizon: Wave ``g`` -> ``H_g = min(H, T - g)``, the rows of its basis.
      window_len: Wave ``g`` -> ``(N,)`` observed prefix length of every unit's
        window ``[g, g + H_g)`` (treated units of other waves hold 0).
      post_len: ``(N,)`` ``L_i``: a treated unit's own window length; a never-treated
        unit's window length for the first wave.
      basis: ``(H, K_eff)`` common orthonormal exposure basis ``Psi``.
      bases: Wave ``g`` -> ``Psi[:H_g]`` (the per-wave view consumed by the result).
      min_support: Per-mode minimum window length (index ``k``, 0-based).
      masked: Whether any cell of the panel is unobserved.
    """

    t0: int
    horizon: int
    n_periods: int
    cohorts: np.ndarray
    ever: np.ndarray
    obs: np.ndarray
    pre_count: np.ndarray
    waves: list[int]
    wave_weights: dict[int, float]
    wave_horizon: dict[int, int]
    window_len: dict[int, np.ndarray]
    post_len: np.ndarray
    basis: np.ndarray
    bases: dict[int, np.ndarray]
    min_support: tuple[int, ...]
    masked: bool
    _cumsum: dict[int, np.ndarray] = dc_field(default_factory=dict, repr=False, compare=False)

    @property
    def n_units(self) -> int:
        return int(self.cohorts.shape[0])

    @property
    def k_eff(self) -> int:
        return int(self.basis.shape[1])

    @property
    def ctrl(self) -> np.ndarray:
        return ~self.ever

    def support_for(self, k: int) -> int:
        """Minimum window length required to project mode index ``k`` (0-based)."""
        return int(self.min_support[min(k, len(self.min_support) - 1)])

    def active_waves(self, k: int) -> list[int]:
        """Waves whose basis holds mode ``k``."""
        return [g for g in self.waves if k < self.bases[g].shape[1]]

    def support_mask(self, k: int) -> np.ndarray:
        """``(N,)`` units whose window supports mode ``k`` (any wave for never-treated units)."""
        need = self.support_for(k)
        out = np.zeros(self.n_units, dtype=bool)
        for g in self.active_waves(k):
            out |= self.window_len[g] >= need
        return out

    def support_matrix(self) -> np.ndarray:
        """``(N, K_eff)`` support of every mode."""
        return np.column_stack([self.support_mask(k) for k in range(self.k_eff)])

    def has_partial_windows(self) -> bool:
        """Whether any unit observes less than the common horizon ``H`` (the correction is then live)."""
        return any(bool(np.any((self.window_len[g] > 0) & (self.window_len[g] < self.horizon))) for g in self.waves)

    def psi_norm2(self, g: int, k: int, length: int) -> float:
        """Squared norm of mode ``k`` on the first ``length`` rows of the basis (exactly ``H`` on a full window).

        Orthonormality under ``(1 / H) sum`` makes the full-window norm ``H`` exactly;
        using the integer avoids the rounding of ``np.sum(psi ** 2)`` and reproduces
        the paper's ``/ L_g``. On a partial window the higher modes (``k >= 1``) are
        centered on the window first (see :func:`project_event_time_mode`), so their
        norm is ``sum_{s <= L} (psi_k(s) - mean_L psi_k)^2``.
        """
        if length == self.horizon:
            return float(length)
        psi = self.bases[g][:length, k]
        if k >= 1:
            psi = psi - psi.mean()
        return float(np.sum(psi ** 2))

    def partial_mean(self, g: int, k: int, length: int) -> float:
        """``(1 / length) sum_{s <= length} psi_k(s)``; exactly 0 on a full window (``length == H``)."""
        if length == self.horizon or length <= 0:
            return 0.0
        if g not in self._cumsum:
            self._cumsum[g] = np.cumsum(self.bases[g], axis=0)
        return float(self._cumsum[g][length - 1, k]) / float(length)

    def treated_cells_observed(self) -> np.ndarray:
        """``(N, T)`` treated post cells inside the horizon that are observed."""
        out = np.zeros(self.obs.shape, dtype=bool)
        for g in self.waves:
            rows = np.nonzero(self.cohorts == g)[0]
            for i in rows:
                L = int(self.window_len[g][i])
                out[i, g : g + L] = True
        return out

    def treated_cells_modeled(self) -> np.ndarray:
        """``(N, T)`` treated post cells inside the horizon (observed or predicted)."""
        out = np.zeros(self.obs.shape, dtype=bool)
        for g in self.waves:
            rows = self.cohorts == g
            out[rows, g : g + self.wave_horizon[g]] = True
        return out


def observed_prefix_lengths(obs_window: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-row count of observed cells and whether the observed cells form a prefix.

    Args:
      obs_window: ``(n, L)`` boolean observation indicator of a window.

    Returns:
      ``(lengths, is_prefix)``: ``lengths[i]`` is the number of ``True`` cells of row
      ``i`` and ``is_prefix[i]`` whether they are exactly the first ``lengths[i]`` cells.
    """
    obs_window = np.asarray(obs_window, dtype=bool)
    lengths = obs_window.sum(axis=1).astype(np.int64)
    n, L = obs_window.shape
    prefix = np.arange(L)[None, :] < lengths[:, None]
    is_prefix = np.all(obs_window == prefix, axis=1)
    return lengths, is_prefix


def check_observation_structure(
    obs: np.ndarray,
    cohorts: np.ndarray,
    t0: int,
    *,
    horizon: int | None,
    min_support: tuple[int, ...],
) -> dict[str, Any]:
    """Validate the administrative-censoring structure of a masked panel.

    Rules (``ValueError`` otherwise): every unit observes at least one pre-period
    cell; within its window ``[g, g + H)`` (``g = t0`` for never-treated units) the
    observed cells of every unit form a contiguous prefix of length
    ``L_i >= min_support[0]``; and every post period in which a treated unit is
    observed has at least one observed never-treated unit (the control period
    mean must exist).

    Args:
      obs: ``(N, T)`` observed-cell indicator.
      cohorts: ``(N,)`` first-adoption periods (``T`` for never-treated units).
      t0: Pre-period length.
      horizon: Event-time horizon ``H`` or ``None`` (the longest treated window).
      min_support: Per-mode minimum window lengths.

    Returns:
      ``{"horizon": H, "post_len": (N,) L_i, "pre_count": (N,)}``.
    """
    obs = np.asarray(obs, dtype=bool)
    N, T = obs.shape
    cohorts = np.asarray(cohorts, dtype=np.int64)
    ever = cohorts < T
    pre_count = obs[:, :t0].sum(axis=1)
    if np.any(pre_count < 1):
        bad = np.nonzero(pre_count < 1)[0]
        raise ValueError(
            f"{bad.size} unit(s) have no observed pre-period cell (t < {t0}), e.g. unit {int(bad[0])}: every unit "
            "needs at least one observed pre-adoption outcome for its pre-period anchor; drop these units"
        )
    h_max = T - int(np.min(cohorts[ever]))
    if horizon is None:
        horizon_eff = 1
        for g in sorted({int(v) for v in cohorts[ever]}):
            rows = cohorts == g
            lengths, _ = observed_prefix_lengths(obs[rows, g:])
            horizon_eff = max(horizon_eff, int(lengths.max()) if lengths.size else 1)
    else:
        horizon_eff = int(horizon)
        if not 1 <= horizon_eff <= h_max:
            raise ValueError(
                f"horizon must be an integer in [1, {h_max}] (at most the post-adoption length of the first "
                f"wave); got {horizon!r}. To predict beyond the last panel period, append all-NaN columns to y."
            )
    post_len = np.zeros(N, dtype=np.int64)
    starts = np.where(ever, cohorts, int(np.min(cohorts[ever]))).astype(np.int64)
    for g in sorted({int(v) for v in starts}):
        rows = np.nonzero(starts == g)[0]
        h_g = min(horizon_eff, T - g)
        lengths, is_prefix = observed_prefix_lengths(obs[rows, g : g + h_g])
        if not np.all(is_prefix):
            i = int(rows[np.argmax(~is_prefix)])
            win = obs[i, g : g + h_g]
            gap = int(np.argmin(win))
            raise ValueError(
                f"unit {i} has an unobserved cell inside its post-adoption window (first gap at t = {g + gap}, "
                f"window [{g}, {g + h_g})): only trailing (administrative) censoring is supported, i.e. the "
                "observed post-adoption cells of a unit must be a contiguous prefix s = 1, ..., L_i"
            )
        post_len[rows] = lengths
    need = int(min_support[0])
    short = np.nonzero(post_len < need)[0]
    if short.size:
        raise ValueError(
            f"{short.size} unit(s) observe fewer than min_support_modes[0] = {need} post-adoption period(s) inside "
            f"the horizon H = {horizon_eff} (e.g. unit {int(short[0])} with L_i = {int(post_len[short[0]])}); drop them "
            "or lower min_support_modes[0]"
        )
    # Control period means must exist wherever a treated outcome is observed.
    ctrl = ~ever
    for g in sorted({int(v) for v in cohorts[ever]}):
        h_g = min(horizon_eff, T - g)
        for t in range(g, g + h_g):
            if np.any(obs[(cohorts == g), t]) and not np.any(obs[ctrl, t]):
                raise ValueError(
                    f"period t = {t} has observed treated outcomes but no observed never-treated outcome, so the "
                    "control period mean of the double demeaning is undefined; mask those treated cells or extend "
                    "the never-treated units' observation"
                )
    return {"horizon": horizon_eff, "post_len": post_len, "pre_count": pre_count.astype(np.int64)}


def build_event_time_layout(
    obs: np.ndarray,
    cohorts: np.ndarray,
    t0: int,
    *,
    horizon: int | None,
    rank_k: int,
    min_support: tuple[int, ...],
    basis_factory: Any,
) -> EventTimeLayout:
    """Build the :class:`EventTimeLayout` of a panel (after :func:`check_observation_structure`).

    Args:
      obs: ``(N, T)`` observed-cell indicator.
      cohorts: ``(N,)`` first-adoption periods.
      t0: Pre-period length.
      horizon: Horizon ``H`` or ``None`` (default of :func:`check_observation_structure`).
      rank_k: Requested number of modes ``K``.
      min_support: Per-mode minimum window lengths.
      basis_factory: ``basis_factory(H, K) -> (H, min(K, H))`` orthonormal basis.

    Returns:
      The layout; ``bases[g]`` is ``Psi[:H_g]`` for every wave.
    """
    info = check_observation_structure(obs, cohorts, t0, horizon=horizon, min_support=min_support)
    obs = np.asarray(obs, dtype=bool)
    N, T = obs.shape
    cohorts = np.asarray(cohorts, dtype=np.int32)
    ever = cohorts < T
    H = int(info["horizon"])
    waves = sorted(int(g) for g in np.unique(cohorts[ever]))
    basis = np.asarray(basis_factory(H, int(rank_k)), dtype=np.float64)
    window_len: dict[int, np.ndarray] = {}
    wave_horizon: dict[int, int] = {}
    bases: dict[int, np.ndarray] = {}
    for g in waves:
        h_g = min(H, T - g)
        wave_horizon[g] = h_g
        bases[g] = basis[:h_g]
        rows = (cohorts == g) | ~ever
        lengths, _ = observed_prefix_lengths(obs[rows, g : g + h_g])
        wl = np.zeros(N, dtype=np.int64)
        wl[rows] = lengths
        window_len[g] = wl
    total = float(sum(int(np.sum(window_len[g][cohorts == g])) for g in waves))
    wave_weights = {g: float(np.sum(window_len[g][cohorts == g])) / max(total, 1.0) for g in waves}
    return EventTimeLayout(
        t0=int(t0),
        horizon=H,
        n_periods=int(T),
        cohorts=cohorts,
        ever=ever,
        obs=obs,
        pre_count=np.asarray(info["pre_count"], dtype=np.int64),
        waves=waves,
        wave_weights=wave_weights,
        wave_horizon=wave_horizon,
        window_len=window_len,
        post_len=np.asarray(info["post_len"], dtype=np.int64),
        basis=basis,
        bases=bases,
        min_support=tuple(int(v) for v in min_support),
        masked=bool(not np.all(obs)),
    )


def _project_rows(
    panel: np.ndarray, rows: np.ndarray, start: int, length: int, psi: np.ndarray, *, center: bool
) -> np.ndarray:
    """``<block, psi[:length]>`` per row of ``panel[rows, start:start + length]`` (2-D panel or 3-D draw stack).

    With ``center`` both the rows and ``psi`` are centered on the window first
    (the slope of the per-unit OLS of the row on ``{1, psi}``).
    """
    psi_w = psi[:length]
    if panel.ndim == 2:
        block = panel[rows, start : start + length]
        if center:
            block = block - block.mean(axis=1, keepdims=True)
            psi_w = psi_w - psi_w.mean()
        return block @ psi_w
    block = panel[:, rows, start : start + length]
    if center:
        block = block - block.mean(axis=2, keepdims=True)
        psi_w = psi_w - psi_w.mean()
    return np.einsum("mit,t->mi", block, psi_w)


def project_event_time_mode(
    panel: np.ndarray,
    layout: EventTimeLayout,
    k: int,
    *,
    supported_only: bool = True,
) -> np.ndarray:
    """Single-mode projection of every unit's observed window onto mode ``k``.

    On a full window this is ``<psi_k, r_i> / H``, the paper's projection. On a
    partial window the constant and a higher mode are no longer orthogonal, so the
    level would leak into ``<psi_k, r_i>_{O_i}``; higher modes (``k >= 1``) are
    therefore projected after the window constant is partialled out,
    ``Cov_{O_i}(psi_k, r_i) / Var_{O_i}(psi_k)`` (the per-unit OLS slope of ``r_i``
    on ``{1, psi_k}``), and mode 1 is the plain window mean (its dynamic correction
    is :func:`mode1_dynamic_correction`).

    Treated units are projected on their own wave window; never-treated units on
    every active wave window, combined with the wave weights (the paper's
    convention). Units whose window does not support mode ``k`` return 0 when
    ``supported_only`` is set (they receive precision weight 0 as well).

    Args:
      panel: ``(N, T)`` double-demeaned panel or ``(M, N, T)`` draw stack; ``NaN``
        is allowed outside the observed windows only.
      layout: Panel layout.
      k: Mode index (0-based).
      supported_only: Zero the projection of unsupported units.

    Returns:
      ``(N,)`` or ``(M, N)`` projections.
    """
    panel = np.asarray(panel, dtype=np.float64)
    active = layout.active_waves(k)
    w_sum = max(sum(layout.wave_weights[g] for g in active), 1e-12)
    need = layout.support_for(k) if supported_only else 1
    N = layout.n_units
    out = np.zeros(N if panel.ndim == 2 else (panel.shape[0], N), dtype=np.float64)
    ctrl_rows = np.nonzero(layout.ctrl)[0]
    for g in active:
        psi = layout.bases[g][:, k]
        w_norm = layout.wave_weights[g] / w_sum
        lengths = layout.window_len[g]
        trt_rows = np.nonzero(layout.cohorts == g)[0]
        for rows, is_ctrl in ((trt_rows, False), (ctrl_rows, True)):
            if rows.size == 0:
                continue
            for L in np.unique(lengths[rows]):
                L = int(L)
                if L < need:
                    continue
                sel = rows[lengths[rows] == L]
                center = bool(k >= 1 and L < layout.horizon)
                proj = _project_rows(panel, sel, g, L, psi, center=center) / layout.psi_norm2(g, k, L)
                if is_ctrl:
                    out[..., sel] += w_norm * proj
                else:
                    out[..., sel] = proj
    return out


def mode1_dynamic_correction(layout: EventTimeLayout, tau_tilde: list[np.ndarray]) -> np.ndarray:
    """``z_i sum_{k >= 2} tau~_ik (1 / L_i) sum_{s <= L_i} psi_k(s)``: the partial-window mean of the fitted dynamics.

    Exactly 0 for never-treated units and for fully observed windows.

    Args:
      layout: Panel layout.
      tau_tilde: ``K_eff - 1`` arrays of shape ``(N,)`` or ``(M, N)`` with the fitted
        treatment coordinates of modes ``2 .. K_eff`` (posterior means or draws).

    Returns:
      ``(N,)`` or ``(M, N)`` correction to subtract from the mode-1 window means.
    """
    if not tau_tilde:
        return np.zeros(layout.n_units, dtype=np.float64)
    lead = tau_tilde[0].shape[:-1]
    out = np.zeros((*lead, layout.n_units), dtype=np.float64)
    for g in layout.waves:
        lengths = layout.window_len[g]
        rows = np.nonzero(layout.cohorts == g)[0]
        for i in rows:
            L = int(lengths[i])
            if L == layout.horizon or L <= 0:
                continue
            total = 0.0
            for j, tt in enumerate(tau_tilde, start=1):
                if j >= layout.bases[g].shape[1]:
                    continue
                total = total + tt[..., i] * layout.partial_mean(g, j, L)
            out[..., i] = total
    return out


def event_time_precision_weights(layout: EventTimeLayout, k: int) -> np.ndarray:
    """Relative precision weights ``w_ik`` of the single-mode projections (control mean 1).

    Mirrors :func:`longbet._rank_k._compute_mode_precision_weights` operation by
    operation on fully observed windows: a treated unit's variance factor is
    ``1 / ||psi_k||^2_{O_i}`` (plus ``1 / n_pre_i`` for mode 1), a never-treated
    unit's the squared norm of its wave-weighted time-weight vector; unsupported
    units get weight 0.

    Args:
      layout: Panel layout.
      k: Mode index (0-based).

    Returns:
      ``(N,)`` weights.
    """
    N, T = layout.n_units, layout.n_periods
    active = layout.active_waves(k)
    w_sum = max(sum(layout.wave_weights[g] for g in active), 1e-12)
    need = layout.support_for(k)
    pre_var = (1.0 / np.maximum(layout.pre_count, 1).astype(np.float64)) if k == 0 else np.zeros(N)
    var_unit = np.full(N, np.inf, dtype=np.float64)
    supported = np.zeros(N, dtype=bool)
    # Never-treated units: time-weight vector over their (wave-weighted) windows.
    ctrl_rows = np.nonzero(layout.ctrl)[0]
    for i in ctrl_rows:
        c = np.zeros(T, dtype=np.float64)
        any_support = False
        for g in active:
            L = int(layout.window_len[g][i])
            if L < need:
                continue
            any_support = True
            psi = layout.bases[g][:, k]
            w_norm = layout.wave_weights[g] / w_sum
            psi_w = psi[:L]
            if k >= 1 and L < layout.horizon:
                psi_w = psi_w - psi_w.mean()
            c[g : g + L] += w_norm * (psi_w / layout.psi_norm2(g, k, L))
        if any_support:
            var_unit[i] = max(float(np.sum(c ** 2)) + pre_var[i], 1e-8)
            supported[i] = True
    for g in active:
        rows = np.nonzero(layout.cohorts == g)[0]
        for i in rows:
            L = int(layout.window_len[g][i])
            if L < need:
                continue
            var_unit[i] = max((1.0 / layout.psi_norm2(g, k, L)) + pre_var[i], 1e-8)
            supported[i] = True
    inv_var = np.where(supported, 1.0 / var_unit, 0.0)
    ctrl_sup = supported & layout.ctrl
    if not np.any(ctrl_sup):
        raise ValueError(
            f"no never-treated unit supports mode {k + 1} (window length >= {need}); lower min_support_modes "
            "or rank_k"
        )
    ctrl_mean = float(np.mean(inv_var[ctrl_sup]))
    return inv_var / max(ctrl_mean, 1e-8)


def restricted_wls_projection(
    panel: np.ndarray,
    layout: EventTimeLayout,
    tau_tilde: list[np.ndarray] | None,
    *,
    cond_max: float = 1e8,
) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray]:
    """Joint per-unit WLS of the supported modes, ``(Psi_O' Psi_O)^{-1} Psi_O' r_O``.

    For unit ``i`` the supported set is ``S_i = {k : L_i >= min_support[k]}`` (nested
    by construction). The omitted modes bias the joint estimate by
    ``(Psi_S' Psi_S)^{-1} Psi_S' psi_k tau_ik``; for treated units this is removed with
    the fitted ``tau~_ik`` (``tau_tilde``, as in the two-pass projection). Units
    whose restricted Gram matrix is ill-conditioned (``cond > cond_max``) fall back
    to the single-mode (two-pass) projection and its weights. On fully observed
    windows the Gram matrix is ``H I`` and the estimates equal the single-mode
    projections exactly.

    Args:
      panel: ``(N, T)`` double-demeaned panel.
      layout: Panel layout.
      tau_tilde: Fitted higher-mode treatment coordinates ``(N,)`` for modes ``2..K``
        (posterior means), or ``None`` for no omitted-mode correction.
      cond_max: Condition-number guard.

    Returns:
      ``(theta, var_factor, fell_back)``: per mode ``k`` the ``(N,)`` estimates and the
      ``(N,)`` variance factors ``[(Psi_S' Psi_S)^{-1}]_kk`` (``inf`` when unsupported),
      and the ``(N,)`` fallback indicator.
    """
    if len(layout.waves) != 1:
        raise ValueError("projection='restricted_wls' needs a single adoption period (event-time panel)")
    panel = np.asarray(panel, dtype=np.float64)
    N = layout.n_units
    K = layout.k_eff
    g = layout.waves[0]
    Psi_g = layout.bases[g]
    lengths = layout.window_len[g]
    theta = [np.zeros(N, dtype=np.float64) for _ in range(K)]
    var_factor = [np.full(N, np.inf, dtype=np.float64) for _ in range(K)]
    fell_back = np.zeros(N, dtype=bool)
    single = [project_event_time_mode(panel, layout, k) for k in range(K)]
    for i in range(N):
        L = int(lengths[i])
        S = [k for k in range(K) if L >= layout.support_for(k)]
        if not S:
            continue
        if L == layout.horizon:
            for k in S:
                theta[k][i] = single[k][i]
                var_factor[k][i] = 1.0 / float(L)
            continue
        Psi_S = Psi_g[:L, S]
        gram = Psi_S.T @ Psi_S
        if not np.all(np.isfinite(gram)) or np.linalg.cond(gram) > cond_max:
            fell_back[i] = True
            for k in S:
                theta[k][i] = single[k][i]
                var_factor[k][i] = 1.0 / layout.psi_norm2(g, k, L)
            continue
        gram_inv = np.linalg.inv(gram)
        est = gram_inv @ (Psi_S.T @ panel[i, g : g + L])
        if tau_tilde and layout.ever[i]:
            for k in range(1, K):
                if k in S or k - 1 >= len(tau_tilde):
                    continue
                b_k = gram_inv @ (Psi_S.T @ Psi_g[:L, k])
                est = est - b_k * float(tau_tilde[k - 1][i])
        for j, k in enumerate(S):
            theta[k][i] = est[j]
            var_factor[k][i] = gram_inv[j, j]
    return theta, var_factor, fell_back


def restricted_wls_precision_weights(layout: EventTimeLayout, var_factor: list[np.ndarray], k: int) -> np.ndarray:
    """Relative precision weights from the restricted-Gram variance factors (control mean 1)."""
    N = layout.n_units
    pre_var = (1.0 / np.maximum(layout.pre_count, 1).astype(np.float64)) if k == 0 else np.zeros(N)
    vf = np.asarray(var_factor[k], dtype=np.float64)
    supported = np.isfinite(vf)
    inv_var = np.where(supported, 1.0 / np.maximum(vf + pre_var, 1e-8), 0.0)
    ctrl_sup = supported & layout.ctrl
    if not np.any(ctrl_sup):
        raise ValueError(f"no never-treated unit supports mode {k + 1}; lower min_support_modes or rank_k")
    return inv_var / max(float(np.mean(inv_var[ctrl_sup])), 1e-8)
