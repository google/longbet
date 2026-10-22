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

"""Late-context structure of the built-in counterfactual forecaster.

The spectral zero-shot forecaster (:func:`longbet._forecast._spectral_damped_zero_shot`)
fits a secular trend and unit-specific harmonics on a *guarded* core window
(the last ``guard = 4`` pre-anchor periods are held out of that fit) and then
adds two late-context components whose *shapes* are fixed:

* a slope-break hinge kinked ``break_kink = 10`` periods before the first
  adoption, and
* a Gaussian pre-adoption dip of width ``transient_width = 1.8`` centred
  ``transient_centre = -1.5`` periods before the first post-anchor period,

with per-unit amplitudes projected onto the dominant cross-sectional harmonic
factor.  The kernel *locations* are an identification restriction the analyst
sets -- the analogue of the anticipation window of an event study -- and
:class:`LateContextSpec` exposes them as parameters so that the sensitivity of
an estimate to that restriction can be reported (move the kink, move or widen
the dip, or switch a component off with ``break_mode="none"`` /
``transient_mode="none"``).  The defaults reproduce the original forecaster
bit for bit.

**Why there is no data-adaptive mode.**  An earlier version of this module
(longbet 0.1.x, commit 23a4680) estimated the kink date and the dip centre/width
by joint profile least squares with empirical-Bayes amplitude pooling and a
Mallows ``C_p`` decision.  In the paired sweep of the companion paper (R = 64
panels per cell, N = 80, sigma = 0.35, ten late-context periods; archived as
``repro/out/archive/v2_engine_23a4680/forecast_bias_diagnostic.json`` in the
paper repository) that estimator lowered the treated forecast RMSE but did
**not** remove the treated-minus-control post-anchor forecast-error
differential -- the quantity the plug-in forecast channel turns into ATT bias.
On the panel whose break and dip sit exactly where the fixed kernels assume, the
search made the differential *worse* than using no kernels at all (+0.18 vs
+0.13, against +0.07 for the fixed kernels): a transient dip and a late slope
change are collinear on the short pre-anchor window, the ``C_p`` search prefers
a break-only hinge about 15 periods back (which also fits curvature of the
secular trend) in most replications whatever the true structure, and its dip
centre, although unbiased on average, varies by about one period from panel to
panel -- selection on pre-treatment noise that in that design co-moves with
selection into treatment.  When the kernels were misplaced the search narrowed
the forecast-level differential (worst case 0.32 vs 0.58) without restoring
ATT coverage in the full MCMC re-runs (47-72% vs 59-91% for the fixed
kernels), raising the ATT bias in three of five misplacements.  Estimating the
fixed shapes jointly with the secular fit did not dominate the guarded two-stage
fit either.  An option that is never preferable is not offered; the machinery
survives only as the descriptive diagnostic :func:`late_context_evidence`, which
reports the profile ``C_p`` curve over kink dates, the ``C_p`` surface over dip
centres and widths, and the treated-minus-control residual in the guard window,
so that an analyst can see what the extrapolation rests on without the
forecaster acting on it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np

__all__ = [
    "LateContextSpec",
    "CoreFit",
    "fit_core",
    "apply_fixed_kernels",
    "apply_fixed_kernels_full",
    "resolve_late_context",
    "assemble_context",
    "late_context_evidence",
]

LateContextMode = Literal["fixed", "none"]
ComponentMode = Literal["fixed", "none"]

#: Hinge scale of the original forecaster (``max(0, t - kink) / 16``).
_HINGE_SCALE = 16.0


@dataclass(frozen=True)
class LateContextSpec:
    """How the forecaster treats the last pre-anchor periods.

    Parameters
    ----------
    mode
        Preset.  ``"fixed"`` (default) -- the hinge and Gaussian dip kernels at
        the locations below; ``"none"`` -- secular trend and harmonics only.
        The component fields override the preset when given.
    break_mode, transient_mode
        ``"fixed"`` or ``"none"`` for the slope-break hinge and the pre-adoption
        dip individually (``None`` = follow ``mode``).
    guard
        Pre-anchor periods held out of the secular/harmonic fit (only when the
        context has at least 12 periods; shorter contexts use no guard and no
        kernels).  The dip amplitude is estimated on this guard window.
    break_kink
        Location of the hinge kink in periods *before the first adoption*
        (default 10).  Must be positive.
    transient_centre
        Centre of the Gaussian dip in periods relative to the first post-anchor
        period (negative = before it; default -1.5).
    transient_width
        Width (SD) of the Gaussian dip in periods (default 1.8).  Must be
        positive.

    Notes
    -----
    The post-anchor recovery of the dip follows the fitted Gaussian -- an
    explicit mean-reversion assumption, the only one under which a transient
    concentrated in treated units can be extrapolated at all.  Report results
    under at least one alternative placement; see :func:`late_context_evidence`
    for the data's (weak) evidence on the placement.
    """

    mode: LateContextMode = "fixed"
    break_mode: ComponentMode | None = None
    transient_mode: ComponentMode | None = None
    guard: int = 4
    break_kink: float = 10.0
    transient_centre: float = -1.5
    transient_width: float = 1.8

    def __post_init__(self) -> None:
        if self.mode not in ("fixed", "none"):
            raise ValueError(
                f"late_context mode must be 'fixed' or 'none', got {self.mode!r} "
                "(the data-adaptive mode was removed: it never dominated the fixed kernels; "
                "use late_context_evidence() for the profile diagnostics and vary break_kink / "
                "transient_centre / transient_width for a sensitivity analysis)"
            )
        if int(self.guard) < 0:
            raise ValueError("guard must be non-negative")
        for name in ("break_mode", "transient_mode"):
            val = getattr(self, name)
            if val is not None and val not in ("fixed", "none"):
                raise ValueError(f"{name} must be 'fixed', 'none' or None, got {val!r}")
        if not np.isfinite(self.break_kink) or float(self.break_kink) <= 0.0:
            raise ValueError("break_kink must be a positive number of periods before the first adoption")
        if not np.isfinite(self.transient_centre):
            raise ValueError("transient_centre must be finite")
        if not np.isfinite(self.transient_width) or float(self.transient_width) <= 0.0:
            raise ValueError("transient_width must be positive")

    # -- resolved components ---------------------------------------------------------------
    @property
    def resolved_break(self) -> str:
        return self.break_mode if self.break_mode is not None else self.mode

    @property
    def resolved_transient(self) -> str:
        return self.transient_mode if self.transient_mode is not None else self.mode

    @property
    def uses_original_path(self) -> bool:
        """True when the spec is the original fixed-kernel forecaster (default locations)."""
        return (
            self.resolved_break == "fixed"
            and self.resolved_transient == "fixed"
            and int(self.guard) == 4
            and float(self.break_kink) == 10.0
            and float(self.transient_centre) == -1.5
            and float(self.transient_width) == 1.8
        )


def resolve_late_context(late_context: LateContextSpec | str | None) -> LateContextSpec:
    """Accept a spec, a mode string or ``None`` (``"fixed"``)."""
    if late_context is None:
        return LateContextSpec()
    if isinstance(late_context, LateContextSpec):
        return late_context
    if isinstance(late_context, str):
        return LateContextSpec(mode=late_context)  # type: ignore[arg-type]
    raise TypeError("late_context must be a LateContextSpec, a mode string or None")


@dataclass
class CoreFit:
    """Secular + harmonic fit on the guarded core of the pre-anchor context.

    All arrays are in the forecaster's internal time coordinates
    (``t_ctx = index - core_end``; ``t_panel`` continues them across the
    analysis panel, so the first post-anchor period sits at ``t = guard``).
    ``ctx_resid`` is the residual of the cross-sectionally demeaned context
    after the core fit, including the held-out guard periods.
    """

    N: int
    L: int
    T_panel: int
    anchor_t: int
    guard: int
    core_end: int
    grand_mean: float
    ctx_dm: np.ndarray
    t_ctx: np.ndarray
    t_panel: np.ndarray
    s_t: float
    Phi_ctx: np.ndarray
    Phi_pan: np.ndarray
    W_Phi: np.ndarray
    gram_Phi: np.ndarray
    coef: np.ndarray
    sec_pred: np.ndarray
    ctx_resid: np.ndarray
    sigma2_core: float
    u1: np.ndarray
    U_design: np.ndarray | None
    P_u: np.ndarray | None
    kernels_allowed: bool


def fit_core(context: np.ndarray, T_panel: int, anchor_t: int, guard: int = 4) -> CoreFit:
    """Fit unit-specific secular trend and harmonics on the guarded core window.

    This is the first stage of :func:`longbet._forecast._spectral_damped_zero_shot`;
    the arithmetic is kept identical so the default spec reproduces the original
    forecaster exactly.
    """
    context = np.asarray(context, dtype=np.float64)
    N, L = context.shape
    guard = int(guard) if L >= 12 else 0
    core_end = max(L - guard, 1)
    grand_mean = float(np.mean(context[:, :core_end]))
    ctx_dm = context - np.mean(context, axis=0, keepdims=True)

    # Anchor-invariant time coordinates centered at core_end
    t_ctx = np.arange(L, dtype=np.float64) - float(core_end)
    t_panel = np.arange(T_panel, dtype=np.float64) + float(guard - anchor_t)
    t_c = float(np.mean(t_ctx[:core_end]))
    s_t = max(float(np.std(t_ctx[:core_end])), 1.0)
    t_lin_ctx = (t_ctx - t_c) / s_t
    t_lin_pan = (t_panel - t_c) / s_t

    basis_ctx_cols = [np.ones(L, dtype=np.float64), t_lin_ctx]
    basis_pan_cols = [np.ones(T_panel, dtype=np.float64), t_lin_pan]
    reg_diag = [1e-6, 1e-4 * float(L)]

    # Require at least 2 full cycles in core history when L >= 12 to avoid half-cycle aliasing
    periods = [
        p
        for p in (4.0, 6.0, 8.0, 12.0, 16.0, 24.0)
        if (core_end >= max(int(2.0 * p), 6) if L >= 12 else L >= max(int(1.25 * p), 4))
    ]
    for p in periods:
        omega = 2.0 * np.pi / p
        basis_ctx_cols.extend([np.cos(omega * t_ctx), np.sin(omega * t_ctx)])
        basis_pan_cols.extend([np.cos(omega * t_panel), np.sin(omega * t_panel)])
        reg_diag.extend([1e-3 * float(L), 1e-3 * float(L)])

    Phi_ctx = np.column_stack(basis_ctx_cols)  # (L, B)
    Phi_pan = np.column_stack(basis_pan_cols)  # (T_panel, B)

    weights = np.ones(L, dtype=np.float64)
    if guard > 0:
        weights[core_end:] = 0.0
    W_Phi = Phi_ctx * weights[:, None]
    gram_Phi = Phi_ctx.T @ W_Phi + np.diag(reg_diag)
    coef = np.linalg.solve(gram_Phi, W_Phi.T @ ctx_dm.T)  # (B, N)
    sec_pred = (Phi_pan @ coef).T  # (N, T_panel)
    ctx_resid = ctx_dm - (Phi_ctx @ coef).T  # (N, L)

    u1 = np.zeros(N, dtype=np.float64)
    U_design: np.ndarray | None = None
    P_u: np.ndarray | None = None
    sigma2_core = float(np.mean(ctx_resid[:, :core_end] ** 2))
    kernels_allowed = bool(guard > 0 and L >= 12 and coef.shape[0] > 2 and N >= 6)
    if kernels_allowed:
        U_harm = coef[2:, :].T
        U_harm_c = U_harm - U_harm.mean(axis=0, keepdims=True)
        u_svd, s_svd, _ = np.linalg.svd(U_harm_c, full_matrices=False)
        harm_power = float(s_svd[0] ** 2) / float(N)
        if harm_power > 2.0 * (sigma2_core / float(core_end)):
            u1 = u_svd[:, 0] * np.sqrt(float(N))
            top_col = int(np.argmax(np.var(U_harm_c, axis=0)))
            if np.std(u1) > 1e-8 and np.std(U_harm_c[:, top_col]) > 1e-8:
                if np.corrcoef(u1, U_harm_c[:, top_col])[0, 1] < 0:
                    u1 = -u1
            u_pos = np.maximum(0.0, u1) - np.mean(np.maximum(0.0, u1))
            U_design = np.column_stack([u1, u_pos / max(float(np.std(u_pos)), 1e-6)])
            k_u = U_design.shape[1]
            P_u = U_design @ np.linalg.solve(U_design.T @ U_design + 0.2 * np.eye(k_u), U_design.T)

    return CoreFit(
        N=N, L=L, T_panel=int(T_panel), anchor_t=int(anchor_t), guard=guard, core_end=core_end,
        grand_mean=grand_mean, ctx_dm=ctx_dm, t_ctx=t_ctx, t_panel=t_panel, s_t=s_t,
        Phi_ctx=Phi_ctx, Phi_pan=Phi_pan, W_Phi=W_Phi, gram_Phi=gram_Phi, coef=coef,
        sec_pred=sec_pred, ctx_resid=ctx_resid, sigma2_core=sigma2_core, u1=u1,
        U_design=U_design, P_u=P_u, kernels_allowed=kernels_allowed,
    )


def apply_fixed_kernels(core: CoreFit, spec: LateContextSpec | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The hinge + Gaussian dip kernels at the spec's locations (default: the original kernels).

    Returns ``(panel_pred, u1, dip_amp, brk_amp)`` exactly as the pre-refactor
    forecaster did for the default spec.
    """
    res = apply_fixed_kernels_full(core, spec)
    return res["panel_pred"], res["u1"], res["dip_amp"], res["brk_amp"]


def apply_fixed_kernels_full(core: CoreFit, spec: LateContextSpec | None = None) -> dict[str, np.ndarray]:
    """:func:`apply_fixed_kernels` plus the in-context residual after both kernels.

    The hinge is kinked ``spec.break_kink`` periods before the first adoption
    (``t_ctx = guard - break_kink``; ``-6`` for the defaults) and orthogonalised
    to the core basis; the dip is the Gaussian ``exp(-0.5 ((t - (guard +
    transient_centre)) / transient_width)^2)`` whose amplitude is estimated on
    the guard window.  Both amplitudes are projected onto the harmonic-factor
    subspace ``P_u``; a component set to ``"none"`` is skipped.  ``ctx_resid``
    (``(N, L)``, cross-sectionally demeaned scale) is what remains of the context
    after secular trend, harmonics and the kept kernels; the forecast-draw
    bootstrap of :func:`longbet._forecast.build_causal_features`
    resamples it.
    """
    spec = LateContextSpec() if spec is None else spec
    N = core.N
    sec_pred = core.sec_pred
    ctx_resid = core.ctx_resid
    t_ctx, t_panel, core_end, guard = core.t_ctx, core.t_panel, core.core_end, core.guard
    dip_amp = np.zeros(N, dtype=np.float64)
    brk_amp = np.zeros(N, dtype=np.float64)
    use_break = spec.resolved_break == "fixed"
    use_dip = spec.resolved_transient == "fixed"
    if core.kernels_allowed and core.P_u is not None and (use_break or use_dip):
        P_u = core.P_u
        if use_break:
            # Late-context local structural break hinge orthogonalized to Phi_ctx
            kink = float(guard) - float(spec.break_kink)
            h_brk_ctx = np.maximum(0.0, t_ctx - kink) / _HINGE_SCALE
            h_brk_pan = np.maximum(0.0, t_panel - kink) / _HINGE_SCALE
            beta_h = np.linalg.solve(core.gram_Phi, core.W_Phi.T @ h_brk_ctx)
            h_orth_ctx = h_brk_ctx - core.Phi_ctx @ beta_h
            h_orth_pan = h_brk_pan - core.Phi_pan @ beta_h
            h_norm2 = float(np.sum(h_orth_ctx[:core_end] ** 2)) + 1e-8
            if h_norm2 > 0.02:
                brk_raw = (ctx_resid[:, :core_end] @ h_orth_ctx[:core_end]) / h_norm2
                brk_amp = P_u @ brk_raw
                sec_pred = sec_pred + brk_amp[:, None] * h_orth_pan[None, :]
                ctx_resid = ctx_resid - brk_amp[:, None] * h_orth_ctx[None, :]

        if use_dip:
            # Mean-reverting pre-adoption dip kernel in the guard window
            t_dip_center = float(guard) + float(spec.transient_centre)
            width = float(spec.transient_width)
            k_dip_ctx = np.exp(-0.5 * ((t_ctx - t_dip_center) / width) ** 2)
            k_dip_pan = np.exp(-0.5 * ((t_panel - t_dip_center) / width) ** 2)
            k_win = k_dip_ctx[core_end:]
            k_norm2 = float(np.sum(k_win**2)) + 1e-8
            dip_raw = (ctx_resid[:, core_end:] @ k_win) / k_norm2
            dip_amp = P_u @ dip_raw
            panel_pred = sec_pred + dip_amp[:, None] * k_dip_pan[None, :] + core.grand_mean
            ctx_resid = ctx_resid - dip_amp[:, None] * k_dip_ctx[None, :]
        else:
            panel_pred = sec_pred + core.grand_mean
    else:
        panel_pred = sec_pred + core.grand_mean
    return {"panel_pred": panel_pred, "u1": core.u1, "dip_amp": dip_amp, "brk_amp": brk_amp,
            "ctx_resid": np.ascontiguousarray(ctx_resid)}


# ----------------------------------------------------------------------------
# Context assembly shared with the feature builder
# ----------------------------------------------------------------------------


def assemble_context(
    y: np.ndarray, z: np.ndarray, y_pre_history: np.ndarray | None = None
) -> tuple[np.ndarray, int, np.ndarray, np.ndarray]:
    """Strictly pre-treatment context, common anchor, control indices and the firewalled panel.

    Returns ``(context (N, L_ctx), anchor_t, ctrl_idx, y_safe)`` where the
    context is ``[y_pre_history | y[:, :anchor_t]]`` (missing in-panel cells
    filled with the unit's pre-anchor mean), ``anchor_t = min_i G_i`` is the
    first period in which any unit is treated, ``ctrl_idx`` indexes the
    never-treated units and ``y_safe`` is ``y`` with every post-anchor cell of
    an ever-treated unit (and every treated cell) set to ``NaN``.
    """
    import warnings

    y_arr = np.asarray(y, dtype=np.float64)
    z_arr = np.asarray(z)
    if y_arr.ndim != 2 or z_arr.ndim != 2 or y_arr.shape != z_arr.shape:
        raise ValueError(
            f"y and z must be 2D arrays of identical shape (N, T), got {y_arr.shape} and {z_arr.shape}"
        )
    N, T = y_arr.shape
    if not np.all(np.isin(z_arr, (0, 1))):
        raise ValueError("z must contain only binary indicators (0 or 1)")
    if T > 1 and np.any(np.diff(z_arr, axis=1) < 0):
        raise ValueError("z must represent an absorbing treatment (no 1 -> 0 transitions)")

    treated_any_unit = np.any(z_arr == 1, axis=0)
    anchor_t = int(np.argmax(treated_any_unit)) if np.any(treated_any_unit) else T

    ctx_parts: list[np.ndarray] = []
    if y_pre_history is not None:
        hist = np.asarray(y_pre_history, dtype=np.float64)
        if hist.ndim != 2 or hist.shape[0] != N:
            raise ValueError(f"y_pre_history must have shape (N, T_hist) with N={N}, got {hist.shape}")
        if not np.all(np.isfinite(hist)):
            raise ValueError("y_pre_history must contain only finite values")
        ctx_parts.append(hist)
    if anchor_t > 0:
        ctx_in_panel = y_arr[:, :anchor_t].copy()
        if not np.all(np.isfinite(ctx_in_panel)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                row_means = np.nanmean(ctx_in_panel, axis=1, keepdims=True)
            row_means = np.where(np.isfinite(row_means), row_means, 0.0)
            ctx_in_panel = np.where(np.isfinite(ctx_in_panel), ctx_in_panel, row_means)
        ctx_parts.append(ctx_in_panel)
    if not ctx_parts:
        raise ValueError(
            "No pre-treatment context periods available (anchor_t=0 and y_pre_history=None). "
            "Provide y_pre_history or ensure at least one pre-treatment period exists."
        )
    context = np.hstack(ctx_parts)

    ever_treated_mask = np.any(z_arr == 1, axis=1)
    ctrl_idx = np.flatnonzero(~ever_treated_mask)
    y_safe = y_arr.copy()
    y_safe[ever_treated_mask, anchor_t:] = np.nan
    y_safe[z_arr == 1] = np.nan
    return context, anchor_t, ctrl_idx, y_safe


# ----------------------------------------------------------------------------
# Descriptive evidence on the late-context placement (no forecasts)
# ----------------------------------------------------------------------------


def _factor_design(core: CoreFit) -> np.ndarray:
    """Intercept + harmonic-factor columns used as the empirical-Bayes prior mean regression."""
    cols = [np.ones(core.N, dtype=np.float64)]
    if core.U_design is not None:
        cols.append(core.U_design)
    return np.column_stack(cols)


def _eb_shrink(raw: np.ndarray, noise_var: float, design: np.ndarray) -> tuple[np.ndarray, float, float, float]:
    """Shrink per-unit estimates toward their ridge prediction from ``design``.

    ``raw_i ~ N(theta_i, noise_var)``, ``theta_i ~ N(design_i' c, s2)``: the
    posterior mean is ``m_i + w (raw_i - m_i)`` with ``w = s2 / (s2 + noise_var)``
    and ``s2`` the method-of-moments excess variance of ``raw - m``.
    Returns ``(shrunk, w, s2, r2_design)``.
    """
    n, k = design.shape
    gram = design.T @ design + 1e-6 * np.eye(k)
    c = np.linalg.solve(gram, design.T @ raw)
    m = design @ c
    resid = raw - m
    dof = max(n - k, 1)
    s2 = max(float(np.sum(resid ** 2)) / dof - float(noise_var), 0.0)
    w = s2 / (s2 + float(noise_var)) if (s2 + noise_var) > 0 else 0.0
    tot = float(np.sum((raw - raw.mean()) ** 2))
    r2 = float(np.clip(1.0 - float(np.sum(resid ** 2)) / max(tot, 1e-12), 0.0, 1.0)) if tot > 0 else 0.0
    return m + w * resid, float(w), float(s2), r2


def _noise_variance(core: CoreFit) -> float:
    """Core residual variance per cell, with the per-unit basis size charged as degrees of freedom."""
    N, core_end = core.N, core.core_end
    ssr0 = float(np.sum(core.ctx_resid[:, :core_end] ** 2))
    dof = max(N * core_end - N * int(core.coef.shape[0]), N)
    return ssr0 / dof


def _gaussian_kernel(t: np.ndarray, centre: float, width: float) -> np.ndarray:
    return np.exp(-0.5 * ((t - float(centre)) / float(width)) ** 2)


def _profile_least_squares(
    core: CoreFit,
    design: np.ndarray,
    sigma2: float,
    kinks: list[float],
    transients: list[tuple[float, float]],
    grid_allowance: float,
) -> dict[str, Any]:
    """Break and transient by joint profile least squares on the whole context (diagnostic).

    For every candidate (kink ``b``, centre ``c``, width ``w``) the design
    ``[Phi | h_b | k_cw]`` -- secular trend, harmonics, hinge and Gaussian
    transient -- is fitted to every unit on all ``L`` context periods (no
    guard), the per-unit hinge and transient amplitudes are shrunk (empirical
    Bayes, :func:`_eb_shrink`) toward their prediction from the harmonic
    factor, the secular coefficients are re-estimated given the shrunk
    amplitudes, and Mallows' ``C_p = SSR + 2 sigma^2 df`` is formed with the
    shrinkage degrees of freedom.  The families {none, break, transient, both}
    compete after adding the grid-search allowance ``grid_allowance * sigma^2 *
    log(#grid)`` for every searched dimension.
    """
    N, L, B = core.N, core.L, int(core.coef.shape[0])
    Phi, Phi_pan = core.Phi_ctx, core.Phi_pan
    reg = np.diag(core.gram_Phi - Phi.T @ core.W_Phi).copy()  # the core fit's ridge diagonal
    Y = core.ctx_dm.T  # (L, N)
    G_phi = Phi.T @ Phi + np.diag(reg)
    k_design = int(design.shape[1])
    anchor_ctx = float(core.guard)
    min_energy = 0.25  # a component needs this much squared norm inside the context to be identifiable

    kinks = [b for b in kinks if float(np.sum((np.maximum(0.0, core.t_ctx - b) / _HINGE_SCALE) ** 2)) > 0.02]
    transients = [(c, w) for (c, w) in transients
                  if float(np.sum(_gaussian_kernel(core.t_ctx, anchor_ctx + c, w) ** 2)) >= min_energy]
    cands_b: list[float | None] = [None] + kinks
    cands_t: list[tuple[float, float] | None] = [None] + transients
    allow_b = float(grid_allowance) * sigma2 * float(np.log(max(len(kinks), 2))) if kinks else 0.0
    allow_t = float(grid_allowance) * sigma2 * float(np.log(max(len(transients), 2))) if transients else 0.0

    def _family(b: float | None, t: tuple[float, float] | None) -> str:
        return {(False, False): "none", (True, False): "break", (False, True): "transient", (True, True): "both"}[
            (b is not None, t is not None)
        ]

    def _fit(b: float | None, t: tuple[float, float] | None) -> dict[str, Any]:
        cols_ctx: list[np.ndarray] = []
        names: list[str] = []
        if b is not None:
            cols_ctx.append(np.maximum(0.0, core.t_ctx - b) / _HINGE_SCALE)
            names.append("break")
        if t is not None:
            cols_ctx.append(_gaussian_kernel(core.t_ctx, anchor_ctx + t[0], t[1]))
            names.append("transient")
        m = len(cols_ctx)
        if m == 0:
            coef_sec = np.linalg.solve(G_phi, Phi.T @ Y)
            resid = Y - Phi @ coef_sec
            ssr = float(np.sum(resid ** 2))
            return {"b": b, "t": t, "family": "none", "ssr": ssr, "cp": ssr, "score": ssr, "A_hat": None,
                    "names": names, "shrink_w": {}, "factor_r2": {}}
        E = np.column_stack(cols_ctx)  # (L, m)
        X = np.hstack([Phi, E])
        G = X.T @ X + np.diag(np.concatenate([reg, np.full(m, 1e-6 * float(L))]))
        G_inv = np.linalg.inv(G)
        C = G_inv @ (X.T @ Y)  # (B + m, N)
        A_raw = C[B:, :]  # (m, N)
        A_hat = np.empty_like(A_raw)   # empirical-Bayes amplitudes
        df_eb = 0.0
        shrink_w: dict[str, float] = {}
        factor_r2: dict[str, float] = {}
        for j in range(m):
            noise_var = sigma2 * float(G_inv[B + j, B + j])
            a_hat, w_eb, _s2, r2 = _eb_shrink(A_raw[j], noise_var, design)
            A_hat[j] = a_hat
            df_eb += k_design + w_eb * max(N - k_design, 0)
            shrink_w[names[j]] = float(w_eb)
            factor_r2[names[j]] = float(r2)
        coef_sec = np.linalg.solve(G_phi, Phi.T @ (Y - E @ A_hat))
        resid = Y - Phi @ coef_sec - E @ A_hat
        ssr_eb = float(np.sum(resid ** 2))
        cp_eb = ssr_eb + 2.0 * sigma2 * df_eb
        score = cp_eb + (allow_b if b is not None else 0.0) + (allow_t if t is not None else 0.0)
        return {"b": b, "t": t, "family": _family(b, t), "ssr": ssr_eb, "cp": cp_eb, "score": score,
                "A_hat": A_hat, "names": names, "shrink_w": shrink_w, "factor_r2": factor_r2}

    del Phi_pan
    fits: dict[tuple[Any, Any], dict[str, Any]] = {}
    best: dict[str, Any] | None = None
    best_by_family: dict[str, dict[str, Any]] = {}
    for b in cands_b:
        for t in cands_t:
            f = _fit(b, t)
            fits[(b, t)] = f
            fam = f["family"]
            if fam not in best_by_family or f["score"] < best_by_family[fam]["score"]:
                best_by_family[fam] = f
            if best is None or f["score"] < best["score"]:
                best = f
    assert best is not None
    b_star, t_star = best["b"], best["t"]
    curve = [{"periods_before_anchor": float(core.guard - b), "ssr": fits[(b, t_star)]["ssr"], "cp": fits[(b, t_star)]["cp"]}
             for b in cands_b if b is not None and (b, t_star) in fits]
    surface = [{"centre_rel_anchor": float(t[0]), "width": float(t[1]), "ssr": fits[(b_star, t)]["ssr"], "cp": fits[(b_star, t)]["cp"]}
               for t in cands_t if t is not None and (b_star, t) in fits]
    scores = {fam: {"ssr": f["ssr"], "cp": f["cp"], "score": f["score"],
                    "kink_periods_before_anchor": (float(core.guard - f["b"]) if f["b"] is not None else None),
                    "centre_rel_anchor": (float(f["t"][0]) if f["t"] is not None else None),
                    "width": (float(f["t"][1]) if f["t"] is not None else None)}
              for fam, f in best_by_family.items()}
    none_score = best_by_family["none"]["score"] if "none" in best_by_family else float("nan")
    return {"best": best, "curve": curve, "surface": surface, "scores": scores, "n_grid_break": len(kinks),
            "n_grid_transient": len(transients), "allowance_break": allow_b, "allowance_transient": allow_t,
            "gain_vs_none": float(none_score - best["score"]) if np.isfinite(none_score) else None}


def late_context_evidence(
    y: np.ndarray,
    z: np.ndarray,
    *,
    y_pre_history: np.ndarray | None = None,
    late_context: LateContextSpec | str | None = None,
    kink_grid: tuple[float, ...] | None = None,
    transient_centre_grid: tuple[float, ...] = (-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5),
    transient_width_grid: tuple[float, ...] = (1.0, 1.4, 1.8, 2.4, 3.0),
    break_min_pre: int = 6,
    break_min_post: int = 3,
    grid_allowance: float = 2.0,
) -> dict[str, Any]:
    """Descriptive evidence on the late-context identification restriction (no forecast is produced).

    Builds the same strictly pre-treatment context as
    :func:`longbet.build_causal_features` (so nothing a treated unit
    does after the anchor enters), fits the guarded core, and reports

    * ``guard_window``: the treated-minus-control mean of the core residual in
      each held-out guard period (a pre-adoption-dip check in the spirit of a
      pre-trend plot), its overall mean and a crude t-statistic;
    * ``fixed_kernels``: the per-unit break and dip amplitudes the spec's fixed
      kernels would attribute, summarised as treated-minus-control means;
    * ``profile``: the joint profile-least-squares ``C_p`` curve over kink
      dates (``kink_grid``; default every admissible core period) and the ``C_p``
      surface over dip centres and widths, with the best family under the grid
      allowance.  This is the data-adaptive estimator the package used to offer
      as ``late_context="adaptive"``; at typical panel sizes it cannot separate a
      late slope break from a late dip (see the module docstring), which is why
      it is reported and not acted upon.

    Returns a JSON-serialisable dictionary.
    """
    spec = resolve_late_context(late_context)
    context, anchor_t, ctrl_idx, _y_safe = assemble_context(y, z, y_pre_history)
    N, L = context.shape
    T = int(np.asarray(y).shape[1])
    core = fit_core(context, T, anchor_t, guard=spec.guard)
    ctrl = np.zeros(N, dtype=bool)
    ctrl[ctrl_idx] = True
    treated = ~ctrl
    out: dict[str, Any] = {
        "anchor_t": int(anchor_t), "context_length": int(L), "guard": int(core.guard),
        "n_controls": int(ctrl.sum()), "n_treated": int(treated.sum()),
        "kernels_allowed": bool(core.kernels_allowed and core.P_u is not None),
        "spec": {"mode": spec.mode, "break_mode": spec.resolved_break, "transient_mode": spec.resolved_transient,
                 "guard": int(spec.guard), "break_kink": float(spec.break_kink),
                 "transient_centre": float(spec.transient_centre), "transient_width": float(spec.transient_width)},
    }
    if core.guard == 0 or not core.kernels_allowed or core.P_u is None:
        out["skipped"] = "context shorter than 12 periods or too few units: no late-context structure is estimated"
        return out
    both = bool(treated.any() and ctrl.any())

    # 1. guard-window residual differential (what the dip kernel acts on)
    R_guard = core.ctx_resid[:, core.core_end:]  # (N, guard)
    diff_t = (R_guard[treated].mean(axis=0) - R_guard[ctrl].mean(axis=0)) if both else np.full(R_guard.shape[1], np.nan)
    d_unit = R_guard.mean(axis=1)
    se = (np.sqrt(d_unit[treated].var(ddof=1) / max(treated.sum(), 2) + d_unit[ctrl].var(ddof=1) / max(ctrl.sum(), 2))
          if both and treated.sum() > 1 and ctrl.sum() > 1 else float("nan"))
    out["guard_window"] = {
        "periods_before_first_adoption": [float(core.guard - j) for j in range(R_guard.shape[1])],
        "treated_minus_control_residual": [float(v) for v in diff_t],
        "mean": float(np.mean(diff_t)) if both else None,
        "t_stat": (float(np.mean(diff_t) / se) if both and np.isfinite(se) and se > 0 else None),
    }

    # 2. what the spec's fixed kernels attribute
    fk = apply_fixed_kernels_full(core, spec)
    out["fixed_kernels"] = {
        "break_amp_treated_minus_control": float(fk["brk_amp"][treated].mean() - fk["brk_amp"][ctrl].mean()) if both else None,
        "dip_amp_treated_minus_control": float(fk["dip_amp"][treated].mean() - fk["dip_amp"][ctrl].mean()) if both else None,
        "break_amp_sd": float(np.std(fk["brk_amp"])),
        "dip_amp_sd": float(np.std(fk["dip_amp"])),
    }

    # 3. profile least squares over the grids (descriptive)
    design = _factor_design(core)
    sigma2 = _noise_variance(core)
    if kink_grid is None:
        kinks = [float(core.t_ctx[j]) for j in range(int(break_min_pre), core.core_end - int(break_min_post))]
    else:
        kinks = [float(core.guard) - float(k) for k in kink_grid]
    trans = [(float(c), float(w)) for c in transient_centre_grid for w in transient_width_grid]
    pr = _profile_least_squares(core, design, sigma2, kinks, trans, grid_allowance)
    best = pr["best"]
    out["profile"] = {
        "sigma2": float(sigma2),
        "best_family": best["family"],
        "best_kink_periods_before_anchor": (float(core.guard - best["b"]) if best["b"] is not None else None),
        "best_centre_rel_anchor": (float(best["t"][0]) if best["t"] is not None else None),
        "best_width": (float(best["t"][1]) if best["t"] is not None else None),
        "gain_vs_none": pr["gain_vs_none"],
        "shrink_w": best["shrink_w"],
        "factor_r2": best["factor_r2"],
        "family_scores": pr["scores"],
        "kink_curve": pr["curve"],
        "transient_surface": pr["surface"],
        "n_grid_break": pr["n_grid_break"],
        "n_grid_transient": pr["n_grid_transient"],
        "allowance_break": pr["allowance_break"],
        "allowance_transient": pr["allowance_transient"],
        "note": "descriptive: profile least squares cannot reliably separate a late slope break from a "
                "late dip at typical panel sizes; the forecaster does not act on this",
    }
    return out
