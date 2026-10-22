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

"""Leak-free counterfactual-forecast features for longitudinal causal forests.

A long pre-treatment history carries information about each unit's untreated
trajectory -- seasonality, local trend breaks, pre-adoption dips -- that a
low-degree polynomial in calendar time does not span.  Using it inside a panel
causal forest (:class:`longbet.LongBet`, :func:`longbet.fit_rank_k_longbet`)
requires three pieces of causal discipline, and this module enforces all of
them by construction:

1. **No post-treatment contamination.**  Every unit's forecast context ends at
   the common anchor ``anchor_t = min_i G_i`` (the first period in which any
   unit is treated).  A forecast built on a rolling window ``Y_{i,1:t-1}`` for
   ``t >= G_i`` would absorb the treatment effect and attenuate ``tau`` toward
   zero.
2. **No horizon asymmetry.**  Treated and untreated units are forecast from the
   same origin over the same horizon, so forecast damping and error variance are
   not functions of exposure ``S_it``.
3. **Unit-constant inputs to the trees.**  The trajectory summaries and the
   selection score are per-unit scalars, so the prognostic forest keeps its
   analytically collapsed unit intercepts (``collapse_gamma_mu = True``) and no
   tree has to discretize a smooth temporal trajectory.

:func:`build_causal_features` returns

* ``point_forecast`` -- an ``(N, T)`` counterfactual baseline forecast from the
  built-in guarded-core spectral forecaster with cross-series patch-ridge
  refinement (:func:`patch_ridge_forecaster`, the backbone validated in the
  companion paper), refined by a K-fold out-of-fold ridge adapter fitted on
  never-treated units only;
* ``baseline_forecast`` -- the same forecast, or ``None`` when the pre-treatment
  gate (context length and out-of-fold non-polynomial R^2) says the cell-level
  conjugate trend block would not benefit from it;
* ``unit_features`` -- ``(N, 2)`` unit-constant trajectory summaries;
* ``ps`` -- an ``(N,)`` cross-fitted selection score from ``(X_i, u_i)``;
* optionally ``baseline_forecast_draws`` -- ``M`` draws from the forecaster's
  predictive law for a modular ("cut") posterior.

Which estimator consumes which output is documented in the package README
("What to use when").
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
import warnings
from typing import Any, Callable

import numpy as np

from longbet._latecontext import (
    LateContextSpec,
    apply_fixed_kernels_full,
    assemble_context,
    fit_core,
    resolve_late_context,
)

#: Minimum pre-treatment context length (in periods) and minimum out-of-fold
#: non-polynomial R^2 required for ``build_causal_features`` to activate
#: ``baseline_forecast``. Below these thresholds, ``LongBet``'s built-in
#: polynomial + DCT-II period basis already spans the signal without adding
#: harmonic estimation variance.
MIN_CONTEXT_LENGTH: int = 12
MIN_BEYOND_POLY_R2: float = 0.05

#: Type of a user-supplied forecaster: ``fn(context, horizon) -> (N, horizon)``
#: array (an optional ``late_context`` keyword receives the resolved
#: :class:`~longbet.LateContextSpec`).
Forecaster = Callable[..., np.ndarray]


@dataclass(frozen=True)
class CausalFeatures:
    r"""Outputs of :func:`build_causal_features`.

    Attributes
    ----------
    baseline_forecast
        Gated counterfactual baseline forecast panel of shape ``(N, T)``,
        float32, for the conjugate trend block of :meth:`longbet.LongBet.fit`:
        equal to :attr:`point_forecast` when the pre-treatment context window
        has length ``L_ctx >= 12`` and ``beyond_poly_r2 >= 0.05``; otherwise
        ``None`` so that :meth:`longbet.LongBet.fit` and
        :meth:`longbet.LongBet.predict` fall back to ``LongBet``'s default
        polynomial + DCT-II trend block on short or purely polynomial panels.
    unit_features
        Unit-constant trajectory summaries of shape ``(N, 2)``, float32:
        column 0 is standardized predicted baseline drift, column 1 is
        standardized predicted seasonal/cyclical amplitude. Can be concatenated
        to ``x`` while keeping ``collapse_gamma_mu = True``.
    ps
        Out-of-fold selection score of shape ``(N,)``, float32, estimating the
        treatment adoption probability from static covariates ``x`` and the
        pre-anchor trajectory summaries.
    anchor_t
        Read-only integer index ``0 <= anchor_t <= T`` on the analysis panel
        marking the common pre-treatment cutoff $t_0 = \min_{i \in \mathcal{I}_{\mathrm{trt}}} G_i$.
    beyond_poly_r2
        Out-of-fold $R^2 \in [0, 1]$ of the counterfactual forecast on
        post-anchor untreated control residuals after double-centering and
        projecting out a degree-2 polynomial trend in calendar time (``0.0``
        when ``L_ctx < 12``).
    point_forecast
        The ungated ``(N, T)`` float32 counterfactual baseline forecast (always
        present).  This is the forecast :func:`longbet.fit_rank_k_longbet`
        consumes: its mode-wise CUPED slope is estimated on never-treated units
        and clipped, so an uninformative forecast is simply weighted toward
        zero there and no gate is needed.
    baseline_forecast_draws
        ``None`` unless ``n_forecast_draws > 0`` was requested *and* the gate is
        open: then ``(M, N, T)`` float32 draws from the forecaster's predictive
        law (see :func:`build_causal_features`), centred exactly on
        ``baseline_forecast``.  Fit one :class:`longbet.LongBet` per draw and
        pool the posterior draws for a modular ("cut") posterior that carries
        forecast uncertainty into the treatment-effect intervals.
    forecast_draws
        Diagnostics of the draws (``n_draws``, ``block_length``, mean
        post-anchor / pre-anchor draw SD, the share of the post-anchor draw
        variance that comes from the control adapter, the control-group
        calibration factor ``kappa``) or ``None``.
    """

    baseline_forecast: np.ndarray | None
    unit_features: np.ndarray
    ps: np.ndarray
    anchor_t: int
    beyond_poly_r2: float = 0.0
    point_forecast: np.ndarray | None = dc_field(default=None, repr=False)
    baseline_forecast_draws: np.ndarray | None = dc_field(default=None, repr=False)
    forecast_draws: dict[str, Any] | None = dc_field(default=None, repr=False)
    #: Ungated counterpart of :attr:`baseline_forecast_draws` (private, mirrors
    #: :attr:`point_forecast`); ``None`` unless ``n_forecast_draws > 0``.
    _unfiltered_forecast_draws: np.ndarray | None = dc_field(default=None, repr=False)
    #: Pre-anchor factor loadings ``lambda_i`` (left singular vectors of the
    #: all-unit pre-anchor residual matrix), shape ``(N, q)``, float32.
    pre_anchor_loadings: np.ndarray | None = dc_field(default=None, repr=False)
    #: Cached leak-free pre-anchor context and control-adapter fit for on-demand
    #: posterior forecast draws inside :func:`longbet.fit_rank_k_longbet`.
    _posterior_draw_state: dict[str, Any] | None = dc_field(
        default=None, repr=False, compare=False
    )

    def draw_posterior_forecasts(
        self, n_draws: int, seed: int | None = None
    ) -> tuple[np.ndarray, dict[str, Any]] | None:
        """Draw ``n_draws`` exact NIG posterior forecast panels on demand."""
        st = self._posterior_draw_state
        if st is None or int(n_draws) < 2:
            return None
        s_eff = int(st["seed"]) + 9973 if seed is None else int(seed)
        draws, diag = _forecast_draws_posterior(
            n_draws=int(n_draws),
            context=st["context"],
            unfiltered_forecast=st["unfiltered_forecast"],
            adapter_fit=st["adapter_fit"],
            T=int(st["T"]),
            anchor_t=int(st["anchor_t"]),
            lc_spec=st["lc_spec"],
            forecaster=st["forecaster"],
            y_safe=st["y_safe"],
            ctrl_idx=st["ctrl_idx"],
            rng=np.random.default_rng(s_eff),
        )
        diag["gate_open"] = bool(self.gate_open)
        return draws.astype(np.float32), diag

    @property
    def point_forecast_draws(self) -> np.ndarray | None:
        """Ungated ``(M, N, T)`` float32 draws centred on :attr:`point_forecast`."""
        return self._unfiltered_forecast_draws

    @property
    def _unfiltered_forecast(self) -> np.ndarray | None:
        """Private alias of :attr:`point_forecast` kept for existing scripts."""
        return self.point_forecast

    @property
    def gate_open(self) -> bool:
        """True when the pre-treatment gate activated the forecast channel."""
        return self.baseline_forecast is not None

    @property
    def recommended_forecast(self) -> np.ndarray | None:
        """Alias for :attr:`baseline_forecast` (which is automatically gated by default)."""
        return self.baseline_forecast

    @property
    def recommended_unit_features(self) -> np.ndarray | None:
        """:attr:`unit_features` when the gate is open, else ``None``.

        The history summaries are estimated from the same pre-anchor context as
        the forecast; when that context is too short or carries no signal beyond
        a quadratic trend, feeding them to the prognostic forest adds estimation
        noise to a model that already spans the signal.  Use this (and
        :attr:`recommended_ps`) to revert fully to default ``LongBet`` on gated-off
        panels.
        """
        return self.unit_features if self.gate_open else None

    @property
    def recommended_ps(self) -> np.ndarray | None:
        """:attr:`ps` when the gate is open, else ``None`` (see :attr:`recommended_unit_features`)."""
        return self.ps if self.gate_open else None

    def augment_x(self, x: np.ndarray | None = None) -> np.ndarray | None:
        """Append :attr:`unit_features` to ``x`` when the gate is open, else return ``x`` unchanged.

        Passing ``x=feats.augment_x(x), x_trt=x, ps=feats.recommended_ps,
        baseline_forecast=feats.baseline_forecast`` activates the full dual-channel
        augmentation when :attr:`gate_open` is ``True`` and reverts identically to
        default ``LongBet`` when :attr:`gate_open` is ``False``.
        """
        if not self.gate_open:
            return x
        if x is None:
            return self.unit_features.copy()
        x_arr = np.asarray(x)
        if x_arr.ndim == 1:
            x_arr = x_arr[:, None]
        if x_arr.shape[0] != self.unit_features.shape[0]:
            raise ValueError(
                f"x must have N={self.unit_features.shape[0]} rows, got {x_arr.shape}"
            )
        return np.hstack([x_arr, self.unit_features])


def _invoke_forecaster(
    forecaster: Forecaster,
    context: np.ndarray,
    horizon: int,
    late_context: LateContextSpec | None = None,
) -> np.ndarray:
    """Call a ``fn(context_2d, horizon) -> (N, horizon)`` forecaster and validate its output.

    If the callable accepts a ``late_context`` keyword, the resolved
    :class:`~longbet.LateContextSpec` is forwarded so that a backbone built on
    :func:`_spectral_damped_zero_shot` applies the same kernels as the context
    stage (this is how :func:`patch_ridge_forecaster` receives it).
    """
    N = int(context.shape[0])
    h = int(horizon)
    if h <= 0:
        return np.zeros((N, 0), dtype=np.float64)
    if not callable(forecaster):
        raise TypeError(
            "forecaster must be None, 'patch_ridge', 'spectral', or a callable "
            "fn(context_2d, horizon) -> np.ndarray of shape (N, horizon)"
        )
    kwargs: dict[str, Any] = {}
    try:
        import inspect

        sig_params = inspect.signature(forecaster).parameters
        if late_context is not None and "late_context" in sig_params:
            kwargs["late_context"] = late_context
    except (TypeError, ValueError):
        pass
    arr = np.asarray(forecaster(context, h, **kwargs), dtype=np.float64)
    if arr.shape != (N, h):
        raise ValueError(f"forecaster returned shape {arr.shape}, expected ({N}, {h})")
    if not np.all(np.isfinite(arr)):
        raise ValueError("forecaster returned non-finite values")
    return arr


def _spectral_damped_zero_shot(
    context: np.ndarray,
    T_panel: int,
    anchor_t: int,
    late_context: LateContextSpec | str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic multi-harmonic + spectral subspace dip/break zero-shot forecaster.

    Fits unit-specific Fourier harmonics and secular trend on the guarded core
    pre-treatment context window ``context[:, :core_end]`` of shape
    ``(N, L_ctx)`` (holding out a boundary guard window ``guard = 4`` when
    ``L_ctx >= 12`` so pre-adoption dips do not contaminate the secular slope),
    and forecasts smoothly across all ``T_panel`` analysis periods.  With
    ``late_context="fixed"`` (default) it then projects a late-context
    structural break (hinge kinked ``break_kink = 10`` periods before the first
    adoption) and a mean-reverting Gaussian pre-adoption dip (width 1.8, centred
    1.5 periods before it) onto the dominant cross-sectional harmonic factor
    subspace -- the original behaviour.  A :class:`~longbet.LateContextSpec`
    moves either kernel or switches it off; ``"none"`` returns the
    secular/harmonic fit only.

    Returns
    -------
    panel_pred : np.ndarray
        Denoised/forecast trajectory over the ``T_panel`` analysis periods,
        shape ``(N, T_panel)``.
    ctx_features : np.ndarray
        Per-unit spectral and factor summary features extracted strictly from
        ``context``, shape ``(N, K_ctx)``: core mean, core SD, harmonic factor
        score ``u1``, break amplitude, and the secular/harmonic coefficients.
    """
    spec = resolve_late_context(late_context)
    panel_pred, ctx_features, _core, _resid = _spectral_damped_zero_shot_full(context, T_panel, anchor_t, spec)
    return panel_pred, ctx_features


def _spectral_damped_zero_shot_full(
    context: np.ndarray,
    T_panel: int,
    anchor_t: int,
    spec: LateContextSpec,
) -> tuple[np.ndarray, np.ndarray, Any, np.ndarray]:
    """:func:`_spectral_damped_zero_shot` plus the :class:`CoreFit` and the in-context residual.

    The residual ``(N, L)`` is what remains of the (cross-sectionally demeaned)
    context after every component the backbone fitted -- secular trend,
    harmonics and the kept hinge / dip kernels -- and is the quantity the
    forecast-draw bootstrap resamples.
    """
    core = fit_core(np.asarray(context, dtype=np.float64), int(T_panel), int(anchor_t), guard=spec.guard)
    res = apply_fixed_kernels_full(core, spec)
    panel_pred, u1, brk_amp, ctx_resid = res["panel_pred"], res["u1"], res["brk_amp"], res["ctx_resid"]

    core_end = core.core_end
    ctx_mean = np.mean(np.asarray(context, dtype=np.float64)[:, :core_end], axis=1, keepdims=True)
    ctx_sd = np.std(core.ctx_dm[:, :core_end], axis=1, keepdims=True)
    ctx_features = np.hstack(
        [
            ctx_mean,
            ctx_sd,
            u1[:, None],
            brk_amp[:, None],
            core.coef[1:, :].T,
        ]
    )
    return panel_pred, ctx_features, core, ctx_resid


def patch_ridge_forecaster(
    context: np.ndarray,
    horizon: int,
    patch_len: int = 4,
    late_context: LateContextSpec | str | None = None,
) -> np.ndarray:
    """The default backbone: guarded-core spectral forecaster + cross-series patch-ridge refinement.

    Fits :func:`_spectral_damped_zero_shot` (unit-specific secular trend and
    harmonics on the guarded core window, the late-context hinge and dip
    kernels on the harmonic-factor subspace) and then refines its first
    ``patch_len`` post-anchor steps with a cross-series ridge regression of
    multi-resolution patch summaries (last-patch level, patch-to-patch change,
    two-patch dispersion) onto the backbone's own in-context residual at the
    end of the core window.  The refinement is validated by a two-half
    out-of-fold ``R^2`` on the core window and dropped when that ``R^2`` is
    below ``0.05``, when the core is shorter than ``4 * patch_len`` periods, or
    when ``N < 12`` -- in which case the output is exactly the spectral
    backbone.  Its contribution decays geometrically over the horizon.

    This is the forecaster used for every forecast arm of the companion paper
    (synthetic and real-data studies and the application).  It is deterministic,
    dependency-free, and runs in milliseconds; pass ``forecaster="spectral"`` to
    :func:`build_causal_features` for the backbone alone.

    Parameters
    ----------
    context
        Pre-anchor context matrix of shape ``(N, L_ctx)``.
    horizon
        Number of post-anchor periods to forecast.
    patch_len
        Patch length (periods) of the residual refinement; default ``4``.
    late_context
        Forwarded to the spectral backbone (``None`` = ``"fixed"``).

    Returns
    -------
    np.ndarray
        Forecast of shape ``(N, horizon)``, float32.
    """
    ctx = np.asarray(context, dtype=np.float64)
    N, L = ctx.shape
    horizon = int(horizon)
    if horizon <= 0:
        return np.zeros((N, 0), dtype=np.float32)

    def _backbone(c: np.ndarray, h: int) -> np.ndarray:
        pred, _ = _spectral_damped_zero_shot(c, h, 0, late_context=late_context)
        return np.asarray(pred, dtype=np.float64)

    base_pred = _backbone(ctx, horizon)

    guard = 4 if L >= 12 else 0
    L_core = max(L - guard, 1)
    if L_core >= 4 * patch_len and N >= 12:
        L_in = L_core - patch_len
        in_pred = _backbone(ctx[:, :L_in], patch_len)
        resid_target = ctx[:, L_in:L_core] - in_pred
        resid_target = resid_target - np.mean(resid_target, axis=0, keepdims=True)

        def _tokens(mat: np.ndarray) -> np.ndarray:
            raw = np.column_stack([
                np.mean(mat[:, -2 * patch_len : -patch_len], axis=1),
                np.mean(mat[:, -patch_len:], axis=1)
                - np.mean(mat[:, -2 * patch_len : -patch_len], axis=1),
                np.std(mat[:, -2 * patch_len :], axis=1),
            ])
            return raw - np.mean(raw, axis=0, keepdims=True)

        tok_in = _tokens(ctx[:, :L_in])
        pen = 8.0 * float(N) * np.eye(tok_in.shape[1], dtype=np.float64)
        half = N // 2
        W_1 = np.linalg.solve(tok_in[:half].T @ tok_in[:half] + 0.5 * pen, tok_in[:half].T @ resid_target[:half])
        W_2 = np.linalg.solve(tok_in[half:].T @ tok_in[half:] + 0.5 * pen, tok_in[half:].T @ resid_target[half:])
        oof_pred = np.vstack([tok_in[:half] @ W_2, tok_in[half:] @ W_1])
        null_mse = float(np.mean(resid_target**2))
        oof_mse = float(np.mean((resid_target - oof_pred) ** 2))
        r2_patch = float(np.clip(1.0 - oof_mse / max(null_mse, 1e-12), 0.0, 1.0))
        if r2_patch >= 0.05:
            W_patch = np.linalg.solve(tok_in.T @ tok_in + pen, tok_in.T @ resid_target)
            tok_eval = _tokens(ctx[:, :L_core])
            step_adj = tok_eval @ W_patch
            decay = np.exp(-np.arange(horizon, dtype=np.float64) / float(patch_len))
            for h in range(horizon):
                base_pred[:, h] += r2_patch * decay[h] * step_adj[:, min(h, patch_len - 1)]

    return base_pred.astype(np.float32)


def _patch_ridge_fit_details(
    context: np.ndarray,
    horizon: int,
    patch_len: int = 4,
    late_context: LateContextSpec | str | None = None,
) -> dict[str, Any] | None:
    """Return exact conjugate NIG sufficient statistics for :func:`patch_ridge_forecaster`."""
    ctx = np.asarray(context, dtype=np.float64)
    N, L = ctx.shape
    horizon = int(horizon)
    guard = 4 if L >= 12 else 0
    L_core = max(L - guard, 1)
    if horizon <= 0 or L_core < 4 * patch_len or N < 12:
        return None

    def _backbone(c: np.ndarray, h: int) -> np.ndarray:
        pred, _ = _spectral_damped_zero_shot(c, h, 0, late_context=late_context)
        return np.asarray(pred, dtype=np.float64)

    L_in = L_core - patch_len
    in_pred = _backbone(ctx[:, :L_in], patch_len)
    resid_target = ctx[:, L_in:L_core] - in_pred
    resid_target = resid_target - np.mean(resid_target, axis=0, keepdims=True)

    def _tokens(mat: np.ndarray) -> np.ndarray:
        raw = np.column_stack([
            np.mean(mat[:, -2 * patch_len : -patch_len], axis=1),
            np.mean(mat[:, -patch_len:], axis=1)
            - np.mean(mat[:, -2 * patch_len : -patch_len], axis=1),
            np.std(mat[:, -2 * patch_len :], axis=1),
        ])
        return raw - np.mean(raw, axis=0, keepdims=True)

    tok_in = _tokens(ctx[:, :L_in])
    p_tok = tok_in.shape[1]
    pen = 8.0 * float(N) * np.eye(p_tok, dtype=np.float64)
    half = N // 2
    W_1 = np.linalg.solve(tok_in[:half].T @ tok_in[:half] + 0.5 * pen, tok_in[:half].T @ resid_target[:half])
    W_2 = np.linalg.solve(tok_in[half:].T @ tok_in[half:] + 0.5 * pen, tok_in[half:].T @ resid_target[half:])
    oof_pred = np.vstack([tok_in[:half] @ W_2, tok_in[half:] @ W_1])
    null_mse = float(np.mean(resid_target**2))
    oof_mse = float(np.mean((resid_target - oof_pred) ** 2))
    r2_patch = float(np.clip(1.0 - oof_mse / max(null_mse, 1e-12), 0.0, 1.0))
    gram_patch = tok_in.T @ tok_in + pen
    W_patch = np.linalg.solve(gram_patch, tok_in.T @ resid_target)
    tok_eval = _tokens(ctx[:, :L_core])
    resid_patch = resid_target - tok_in @ W_patch
    dof_unit_patch = max(N - p_tok, 1)
    dof_patch = max(dof_unit_patch * patch_len, 1)
    ssr_patch = float(np.sum(resid_patch**2))
    sigma2_patch = ssr_patch / float(dof_patch)
    sigma_patch_mat = (resid_patch.T @ resid_patch) / float(dof_unit_patch)
    sigma_patch_mat = 0.5 * (sigma_patch_mat + sigma_patch_mat.T)
    decay = np.exp(-np.arange(horizon, dtype=np.float64) / float(patch_len))
    return {
        "tok_eval": tok_eval,
        "W_patch": W_patch,
        "gram_patch": gram_patch,
        "gram_inv": np.linalg.inv(gram_patch),
        "r2_patch": float(r2_patch),
        "sigma2_patch": float(sigma2_patch),
        "sigma_patch_mat": sigma_patch_mat,
        "dof_patch": int(dof_patch),
        "ssr_patch": float(ssr_patch),
        "patch_len": int(patch_len),
        "horizon": int(horizon),
        "decay": decay,
    }


def extract_pre_anchor_loadings(ctx_resid: np.ndarray, q: int = 3) -> np.ndarray:
    """Extract standardized left singular vectors ``lambda_i`` of the pre-anchor residual matrix.

    Uses only the strictly pre-anchor context residual matrix ``ctx_resid`` of
    shape ``(N, L_ctx)``, double-centered across units and context periods and
    row-standardized with a soft scale floor so heteroskedastic outlier units do
    not dominate the singular vectors, with deterministic column sign so the output
    is invariant across BLAS backends and satisfies Proposition 1 (exact invariance
    to post-anchor outcomes).
    """
    R = np.asarray(ctx_resid, dtype=np.float64)
    N, L = R.shape
    q_int = int(max(q, 0))
    if q_int == 0 or L < 2 or N < 2:
        return np.zeros((N, 0), dtype=np.float32)
    R_dc = R - R.mean(axis=1, keepdims=True) - R.mean(axis=0, keepdims=True) + R.mean()
    row_sd = np.std(R_dc, axis=1, keepdims=True)
    med_sd = max(float(np.median(row_sd)), 1e-6)
    R_norm = R_dc / np.maximum(row_sd, 0.5 * med_sd)
    R_norm = R_norm - R_norm.mean(axis=0, keepdims=True)
    u_svd, s_svd, _ = np.linalg.svd(R_norm, full_matrices=False)
    q_eff = min(q_int, u_svd.shape[1], max(N - 1, 1))
    cols: list[np.ndarray] = []
    for j in range(q_eff):
        if s_svd[j] <= 1e-10:
            cols.append(np.zeros(N, dtype=np.float64))
            continue
        col = u_svd[:, j].copy()
        idx_max = int(np.argmax(np.abs(col)))
        if col[idx_max] < 0.0:
            col = -col
        col = col - np.mean(col)
        sd = float(np.std(col))
        if sd > 1e-8:
            z_col = col / sd
            z_clip = 2.5 * np.tanh(z_col / 2.5)
            z_clip = z_clip - np.mean(z_clip)
            sd_clip = float(np.std(z_clip))
            cols.append(z_clip / sd_clip if sd_clip > 1e-8 else np.zeros(N, dtype=np.float64))
        else:
            cols.append(np.zeros(N, dtype=np.float64))
    return np.column_stack(cols).astype(np.float32) if cols else np.zeros((N, 0), dtype=np.float32)


def _resolve_forecaster(forecaster: Forecaster | str | None) -> Forecaster | None:
    """Map the ``forecaster`` argument of :func:`build_causal_features` to a callable.

    ``None`` / ``"patch_ridge"`` -> :func:`patch_ridge_forecaster` (the default);
    ``"spectral"`` -> ``None`` internally, i.e. the spectral backbone alone;
    any callable -> itself.
    """
    if forecaster is None or forecaster == "patch_ridge":
        return patch_ridge_forecaster
    if forecaster == "spectral":
        return None
    if callable(forecaster):
        return forecaster
    raise TypeError(
        "forecaster must be None, 'patch_ridge', 'spectral', or a callable "
        f"fn(context_2d, horizon) -> np.ndarray; got {forecaster!r}"
    )


def _kfold_indices(n: int, n_folds: int, rng: np.random.Generator) -> list[np.ndarray]:
    """Partition ``0 .. n-1`` into ``n_folds`` approximately equal random folds."""
    k = max(2, min(int(n_folds), n))
    perm = rng.permutation(n)
    return [perm[i::k] for i in range(k)]


def _cross_fit_control_adapter(
    y_ctrl: np.ndarray,
    ctx_mean_ctrl: np.ndarray,
    H_all: np.ndarray,
    ctrl_idx: np.ndarray,
    zero_shot_all: np.ndarray,
    n_folds: int,
    rng: np.random.Generator,
    return_fit: bool = False,
) -> np.ndarray | tuple[np.ndarray, dict[str, Any] | None]:
    """Refine zero-shot forecasts via K-fold cross-fitting on never-treated controls.

    Uses **only** never-treated control units ``ctrl_idx`` (where ``Z_iT = 0``) to
    learn a multi-output ridge mapping from pre-treatment embeddings ``H_i`` to
    double-centered untreated panel trajectory residuals
    ``Y_{i, 1:T} - zero_shot_{i, 1:T}``.  With ``return_fit=True`` the full
    control fit (standardised design, ridge weights, Gram matrix, residual
    variance, shrink weight) is returned too, so that
    :func:`build_causal_features` can draw adapter weights from their
    Bayesian-ridge posterior for the forecast draws.
    """
    del ctx_mean_ctrl
    n_ctrl = len(ctrl_idx)
    if n_ctrl < 6:
        return (zero_shot_all.copy(), None) if return_fit else zero_shot_all.copy()

    # Standardize design features H using control-group statistics
    H_ctrl = H_all[ctrl_idx]
    h_mean = H_ctrl.mean(axis=0, keepdims=True)
    h_scale = H_ctrl.std(axis=0, keepdims=True)
    h_scale = np.where(h_scale < 1e-8, 1.0, h_scale)
    Hs_all = (H_all - h_mean) / h_scale
    Hs_ctrl = Hs_all[ctrl_idx]

    # Double-center residual target so the adapter refines within-unit trajectories
    # without perturbing unit means or common period shocks
    target_raw = y_ctrl - zero_shot_all[ctrl_idx]  # (n_ctrl, T)
    target_ctrl = (
        target_raw
        - target_raw.mean(axis=1, keepdims=True)
        - target_raw.mean(axis=0, keepdims=True)
        + target_raw.mean()
    )
    p_dim = Hs_ctrl.shape[1]

    folds = _kfold_indices(n_ctrl, n_folds, rng)
    oof_null = np.zeros_like(target_ctrl)
    for val_local in folds:
        train_mask = np.ones(n_ctrl, dtype=bool)
        train_mask[val_local] = False
        oof_null[val_local] = target_ctrl[train_mask].mean(axis=0, keepdims=True)
    null_mse = float(np.mean((target_ctrl - oof_null) ** 2))

    candidate_alphas = (5.0, 20.0, 80.0, 250.0, 1000.0)
    best_alpha = 80.0
    best_mse = np.inf
    best_oof_dev = np.zeros_like(target_ctrl)
    eye = np.eye(p_dim, dtype=np.float64)

    for alpha in candidate_alphas:
        oof_dev = np.zeros_like(target_ctrl)
        for val_local in folds:
            train_mask = np.ones(n_ctrl, dtype=bool)
            train_mask[val_local] = False
            x_tr_mean = Hs_ctrl[train_mask].mean(axis=0, keepdims=True)
            y_tr_mean = target_ctrl[train_mask].mean(axis=0, keepdims=True)
            X_tr = Hs_ctrl[train_mask] - x_tr_mean
            Y_tr = target_ctrl[train_mask] - y_tr_mean
            alpha_fold = alpha * (X_tr.shape[0] / float(n_ctrl))
            gram = X_tr.T @ X_tr + alpha_fold * eye
            W = np.linalg.solve(gram, X_tr.T @ Y_tr)
            oof_dev[val_local] = (Hs_ctrl[val_local] - x_tr_mean) @ W
        mse = float(np.mean(((target_ctrl - oof_null) - oof_dev) ** 2))
        if mse < best_mse:
            best_mse = mse
            best_alpha = alpha
            best_oof_dev = oof_dev

    r2_ctrl = float(np.clip(1.0 - best_mse / max(null_mse, 1e-12), 0.0, 1.0))
    shrink_weight = r2_ctrl if r2_ctrl >= 0.05 else 0.0

    # Full control-group centered fit for treated units (out-of-sample by construction)
    Y_ctrl_centered = target_ctrl - target_ctrl.mean(axis=0, keepdims=True)
    gram_full = Hs_ctrl.T @ Hs_ctrl + best_alpha * eye
    W_full = np.linalg.solve(gram_full, Hs_ctrl.T @ Y_ctrl_centered)
    correction_all = Hs_all @ W_full
    correction_all[ctrl_idx] = best_oof_dev
    correction_all = correction_all - correction_all.mean(axis=1, keepdims=True)

    refined = zero_shot_all + shrink_weight * correction_all
    if not return_fit:
        return refined
    resid_full = Y_ctrl_centered - Hs_ctrl @ W_full
    dof = max(n_ctrl * Y_ctrl_centered.shape[1] - p_dim * Y_ctrl_centered.shape[1], 1)
    fit = {
        "Hs_all": Hs_all, "W_full": W_full, "gram_full": gram_full, "alpha": float(best_alpha),
        "shrink_weight": float(shrink_weight), "r2_ctrl": r2_ctrl,
        "sigma2": float(np.sum(resid_full ** 2)) / dof,
        "resid_full": resid_full,
        "n_ctrl": int(n_ctrl),
    }
    return refined, fit


def _cross_fit_propensity(
    features: np.ndarray,
    z: np.ndarray,
    n_folds: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Estimate out-of-fold scalar propensity/timing score from (x, unit_features).

    Combines static covariates ``X_i`` and the counterfactual trajectory
    summaries ``u_i`` to predict unit adoption propensity/intensity via
    K-fold cross-fitted L2-penalized regression without fold-intercept jackknife bias.
    """
    N, _ = z.shape
    ever_treated = np.any(z == 1, axis=1).astype(np.float64)
    exposure_share = np.mean(z == 1, axis=1).astype(np.float64)
    y_ps = 0.5 * ever_treated + 0.5 * exposure_share

    mean_rate = float(np.clip(np.mean(ever_treated), 0.05, 0.95))
    if N < 6 or np.all(ever_treated == ever_treated[0]):
        return np.full(N, mean_rate, dtype=np.float32)

    f_mean = features.mean(axis=0, keepdims=True)
    f_scale = features.std(axis=0, keepdims=True)
    f_scale = np.where(f_scale < 1e-8, 1.0, f_scale)
    X_std = (features - f_mean) / f_scale
    p0 = min(X_std.shape[1], 6)
    quad_cols = [X_std[:, j] ** 2 - 1.0 for j in range(p0)]
    X_design = (
        np.hstack([X_std, np.column_stack(quad_cols)]) if quad_cols else X_std
    )

    p_dim = X_design.shape[1]
    eye = np.eye(p_dim, dtype=np.float64)
    folds = _kfold_indices(N, n_folds, rng)

    oof_null = np.zeros(N, dtype=np.float64)
    for val_idx in folds:
        tr_mask = np.ones(N, dtype=bool)
        tr_mask[val_idx] = False
        oof_null[val_idx] = float(np.mean(y_ps[tr_mask]))
    null_mse = float(np.mean((y_ps - oof_null) ** 2))

    best_mse = np.inf
    best_oof_dev = np.zeros(N, dtype=np.float64)

    for alpha in (2.0, 10.0, 30.0, 100.0):
        oof_dev = np.zeros(N, dtype=np.float64)
        for val_idx in folds:
            tr_mask = np.ones(N, dtype=bool)
            tr_mask[val_idx] = False
            x_tr_mean = X_design[tr_mask].mean(axis=0, keepdims=True)
            y_tr_mean = float(np.mean(y_ps[tr_mask]))
            X_tr = X_design[tr_mask] - x_tr_mean
            y_tr = y_ps[tr_mask] - y_tr_mean
            alpha_fold = alpha * (X_tr.shape[0] / float(N))
            w = np.linalg.solve(X_tr.T @ X_tr + alpha_fold * eye, X_tr.T @ y_tr)
            oof_dev[val_idx] = (X_design[val_idx] - x_tr_mean) @ w
        mse = float(np.mean(((y_ps - oof_null) - oof_dev) ** 2))
        if mse < best_mse:
            best_mse = mse
            best_oof_dev = oof_dev

    r2_oof = float(np.clip(1.0 - best_mse / max(null_mse, 1e-12), 0.0, 1.0))
    ps_sd = float(np.std(best_oof_dev))
    if ps_sd < 1e-6 or r2_oof <= 0.0:
        return np.full(N, mean_rate, dtype=np.float32)
    z_score = (best_oof_dev - np.mean(best_oof_dev)) / ps_sd
    logit_base = np.log(mean_rate / (1.0 - mean_rate))
    signal_scale = 0.85 * np.sqrt(r2_oof)
    ps_prob = 1.0 / (
        1.0 + np.exp(-np.clip(logit_base + signal_scale * z_score, -4.0, 4.0))
    )
    return np.clip(ps_prob, 0.02, 0.98).astype(np.float32)


def _zero_shot_stage(
    context: np.ndarray,
    T: int,
    anchor_t: int,
    lc_spec: LateContextSpec,
    forecaster: Forecaster | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Backbone forecast for one context matrix.

    Returns ``(zero_shot (N, T), ctx_features (N, K), ctx_resid (N, L))`` where
    ``ctx_resid`` is the in-context residual after every component the
    spectral backbone fitted.  ``forecaster`` is the *resolved* forecaster
    (see :func:`_resolve_forecaster`): a callable populates the post-anchor
    horizon and the spectral fit is kept on the pre-anchor periods; ``None``
    keeps the spectral fit everywhere.
    """
    spectral_pred, ctx_features, _core, ctx_resid = _spectral_damped_zero_shot_full(context, T, anchor_t, lc_spec)
    if forecaster is not None:
        horizon = T - anchor_t
        ext_post = _invoke_forecaster(forecaster, context, horizon, late_context=lc_spec)
        # Use denoised spectral fit on [:anchor_t] and the forecaster on [anchor_t:]
        zero_shot = spectral_pred.copy()
        if horizon > 0:
            zero_shot[:, anchor_t:] = ext_post
    else:
        zero_shot = spectral_pred
    return zero_shot, ctx_features, ctx_resid


def _moving_block_indices(L: int, block_length: int, rng: np.random.Generator) -> np.ndarray:
    """Time indices of a moving-block bootstrap replicate of length ``L``."""
    ell = int(max(1, min(block_length, L)))
    n_blocks = int(np.ceil(L / ell))
    starts = rng.integers(0, L - ell + 1, size=n_blocks)
    idx = np.concatenate([np.arange(s, s + ell) for s in starts])[:L]
    return idx


def _forecast_draws(
    *,
    n_draws: int,
    block_length: int | None,
    context: np.ndarray,
    ctx_resid: np.ndarray,
    zero_shot: np.ndarray,
    unfiltered_forecast: np.ndarray,
    adapter_fit: dict[str, Any] | None,
    T: int,
    anchor_t: int,
    lc_spec: LateContextSpec,
    forecaster: Forecaster | None,
    y_safe: np.ndarray,
    ctrl_idx: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, dict[str, Any]]:
    r"""Draws from the forecaster's predictive law, centred on the point forecast.

    Two sources of estimation uncertainty are propagated and then calibrated on
    never-treated control units:

    1. **Backbone (spectral fit, kernel amplitudes, cross-unit factor):** a
       moving-block bootstrap of the in-context residual matrix that resamples
       *time blocks jointly across units* (so cross-sectional correlation and
       common-timing shocks are preserved), added to the fitted context, and the
       whole forecaster re-run on the pseudo-context.  The deviation
       ``zero_shot* - zero_shot`` is the backbone perturbation.
    2. **Control adapter:** one draw of the ridge weights from their
       Bayesian-ridge posterior ``N(W, sigma^2 (H'H + alpha I)^{-1})``, mapped
       through the standardised design of *all* units and scaled by the
       adapter's shrink weight.
    3. **Control-group error calibration:** the centred deviations are inflated
       by a factor ``kappa >= 1`` chosen so that, on never-treated units'
       post-anchor cells, sample draw variance (``ddof=1``) + noise variance
       matches the observed squared forecast error (all three double-demeaned,
       as LongBet's period effects and unit intercepts would).  Nothing a
       treated unit does after the anchor enters any draw: the bootstrap touches
       only the pre-anchor context, and the adapter and calibration touch only
       never-treated units.
    """
    M = int(n_draws)
    N, L = context.shape
    ell = int(block_length) if block_length is not None else int(max(2, round(L ** (1.0 / 3.0)) + 1))
    resid = np.asarray(ctx_resid, dtype=np.float64)
    resid = resid - resid.mean(axis=1, keepdims=True)  # centre per unit so E*[pseudo-context] = fitted context
    fitted_ctx = np.asarray(context, dtype=np.float64) - resid
    backbone_dev = np.zeros((M, N, T), dtype=np.float64)
    for m in range(M):
        idx = _moving_block_indices(L, ell, rng)
        ctx_star = fitted_ctx + resid[:, idx]
        zs_star, _feat, _res = _zero_shot_stage(ctx_star, T, anchor_t, lc_spec, forecaster)
        backbone_dev[m] = zs_star - zero_shot

    adapter_dev = np.zeros((M, N, T), dtype=np.float64)
    if adapter_fit is not None and adapter_fit["shrink_weight"] > 0.0:
        cov = adapter_fit["sigma2"] * np.linalg.inv(adapter_fit["gram_full"])
        cov = 0.5 * (cov + cov.T)
        try:
            chol = np.linalg.cholesky(cov + 1e-12 * np.eye(cov.shape[0]))
        except np.linalg.LinAlgError:
            w_, v_ = np.linalg.eigh(cov)
            chol = v_ * np.sqrt(np.maximum(w_, 0.0))[None, :]
        p_dim, T_out = adapter_fit["W_full"].shape
        for m in range(M):
            dW = chol @ rng.standard_normal((p_dim, T_out))
            corr = adapter_fit["Hs_all"] @ dW
            corr = corr - corr.mean(axis=1, keepdims=True)  # the point correction is row-centred too
            adapter_dev[m] = adapter_fit["shrink_weight"] * corr

    dev = backbone_dev + adapter_dev
    if M > 1:
        dev = dev - dev.mean(axis=0, keepdims=True)
    post = slice(anchor_t, T) if anchor_t < T else slice(0, T)
    var_backbone = float(np.mean(np.var(backbone_dev[:, :, post], axis=0, ddof=1))) if M > 1 else 0.0
    var_adapter = float(np.mean(np.var(adapter_dev[:, :, post], axis=0, ddof=1))) if M > 1 else 0.0

    # Leak-free calibration on never-treated units: the bootstrap carries the forecaster's
    # *estimation* uncertainty but not its misspecification.  Compare the observable
    # post-anchor control error (y - F), double-demeaned as LongBet's period effects and unit
    # intercepts would, with the draws' double-demeaned dispersion plus the noise variance
    # estimated from the pre-anchor core residuals, and inflate the deviations so that the two
    # agree.  Treated post-anchor cells are never touched.
    kappa = 1.0
    calib: dict[str, Any] = {"applied": False}
    T_post = T - anchor_t
    ctrl = np.asarray(ctrl_idx, dtype=int)
    if M > 1 and T_post >= 2 and ctrl.size >= 6:
        y_ctl = np.asarray(y_safe, dtype=np.float64)[ctrl][:, post]
        if np.all(np.isfinite(y_ctl)):
            def _dd(a: np.ndarray) -> np.ndarray:
                return a - a.mean(axis=0, keepdims=True) - a.mean(axis=1, keepdims=True) + a.mean()

            n_c = int(ctrl.size)
            v_tot = float(np.mean(_dd(y_ctl - unfiltered_forecast[ctrl][:, post]) ** 2))
            v_est = float(np.mean(np.var(np.stack([_dd(dev[m][ctrl][:, post]) for m in range(M)]), axis=0, ddof=1)))
            core_ = fit_core(np.asarray(context, dtype=np.float64), T, anchor_t, guard=lc_spec.guard)
            ssr0 = float(np.sum(core_.ctx_resid[:, : core_.core_end] ** 2))
            dof0 = max(core_.N * core_.core_end - core_.N * int(core_.coef.shape[0]), core_.N)
            sigma2_dd = (ssr0 / dof0) * (1.0 - 1.0 / n_c) * (1.0 - 1.0 / T_post)
            excess = v_tot - sigma2_dd
            if v_est > 0.0:
                kappa = float(np.sqrt(np.clip(excess / v_est, 1.0, 9.0)))
            calib = {"applied": True, "kappa": kappa, "control_error_var_dd": v_tot, "noise_var_dd": sigma2_dd,
                     "draw_var_dd": v_est, "n_controls": n_c}
            dev = dev * kappa
    draws = unfiltered_forecast[None, :, :] + dev
    diag = {
        "draw_method": "bootstrap",
        "n_draws": M,
        "block_length": ell,
        "calibration": calib,
        "post_anchor_sd_mean": float(np.mean(np.std(dev[:, :, post], axis=0, ddof=1))) if M > 1 else 0.0,
        "pre_anchor_sd_mean": (float(np.mean(np.std(dev[:, :, :anchor_t], axis=0, ddof=1))) if (M > 1 and anchor_t > 0) else 0.0),
        "adapter_variance_share": (var_adapter / (var_backbone + var_adapter)) if (var_backbone + var_adapter) > 0 else 0.0,
        "adapter_sigma2": (float(adapter_fit["sigma2"]) if adapter_fit is not None else None),
    }
    return draws, diag


def _forecast_draws_posterior(
    *,
    n_draws: int,
    context: np.ndarray,
    unfiltered_forecast: np.ndarray,
    adapter_fit: dict[str, Any] | None,
    T: int,
    anchor_t: int,
    lc_spec: LateContextSpec,
    forecaster: Forecaster | None,
    y_safe: np.ndarray,
    ctrl_idx: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Exact conjugate Matrix-Normal-Inverse-Gamma posterior draws of Module 1.

    Propagates parameter uncertainty from the two cross-series ridge regressions
    in Module 1 -- the pre-anchor patch-ridge refinement ``W_patch`` and the
    never-treated control adapter ``W_full`` -- preserving empirical cross-period
    residual covariance without bootstrap resampling or ``kappa`` inflation.
    A posterior-predictive check on never-treated post-anchor cells is recorded
    in ``diag["calibration"]`` as a diagnostic only.
    """
    M = int(n_draws)
    N, _L = context.shape
    T_post = max(int(T - anchor_t), 0)
    post = slice(anchor_t, T) if anchor_t < T else slice(0, T)

    backbone_dev = np.zeros((M, N, T), dtype=np.float64)
    analytical_var_cell = np.zeros((N, T), dtype=np.float64)
    analytical_cov_post_mean = np.zeros((N, N), dtype=np.float64)

    pr_fit = (
        _patch_ridge_fit_details(context, T_post, late_context=lc_spec)
        if (forecaster is patch_ridge_forecaster and T_post > 0)
        else None
    )
    if pr_fit is not None:
        tok_eval = pr_fit["tok_eval"]
        gram_inv_p = 0.5 * (pr_fit["gram_inv"] + pr_fit["gram_inv"].T)
        sigma_p_mat = pr_fit["sigma_patch_mat"]
        dof_p = int(pr_fit["dof_patch"])
        patch_len = int(pr_fit["patch_len"])
        w_p = max(float(pr_fit["r2_patch"]), 0.50)
        decay = pr_fit["decay"]
        p_tok = tok_eval.shape[1]

        try:
            chol_p = np.linalg.cholesky(gram_inv_p + 1e-12 * np.eye(p_tok))
        except np.linalg.LinAlgError:
            w_, v_ = np.linalg.eigh(gram_inv_p)
            chol_p = v_ * np.sqrt(np.maximum(w_, 0.0))[None, :]

        try:
            chol_tp = np.linalg.cholesky(sigma_p_mat + 1e-12 * np.eye(patch_len))
        except np.linalg.LinAlgError:
            wt_, vt_ = np.linalg.eigh(sigma_p_mat)
            chol_tp = vt_ * np.sqrt(np.maximum(wt_, 0.0))[None, :]

        # Scalar NIG inverse-gamma scale multiplier with exact mean E[eta^2] = 1
        a_p = max(0.5 * float(dof_p), 2.5)
        b_p = a_p - 1.0
        for m in range(M):
            g_draw = max(float(rng.gamma(a_p, 1.0)), 1e-30)
            eta_m = np.sqrt(b_p / g_draw)
            Z_p = rng.standard_normal((p_tok, patch_len))
            dW_p = eta_m * (chol_p @ Z_p @ chol_tp.T)
            step_adj = tok_eval @ dW_p
            for h in range(T_post):
                backbone_dev[m, :, anchor_t + h] = (
                    w_p * decay[h] * step_adj[:, min(h, patch_len - 1)]
                )

        cov_unit_p = tok_eval @ gram_inv_p @ tok_eval.T
        diag_unit_p = np.diag(cov_unit_p)
        diag_time_p = np.diag(sigma_p_mat)
        c_weights = np.zeros(patch_len, dtype=np.float64)
        for h in range(T_post):
            scale_h = w_p * decay[h]
            k_h = min(h, patch_len - 1)
            analytical_var_cell[:, anchor_t + h] += (scale_h ** 2) * diag_time_p[k_h] * diag_unit_p
            c_weights[k_h] += scale_h / float(T_post)
        time_quad_p = float(c_weights @ sigma_p_mat @ c_weights)
        analytical_cov_post_mean = analytical_cov_post_mean + time_quad_p * cov_unit_p

    adapter_dev = np.zeros((M, N, T), dtype=np.float64)
    if adapter_fit is not None:
        Hs_all = adapter_fit["Hs_all"]
        gram_inv_a = np.linalg.inv(adapter_fit["gram_full"])
        gram_inv_a = 0.5 * (gram_inv_a + gram_inv_a.T)
        w_a = max(float(adapter_fit["shrink_weight"]), 0.50)
        p_dim, T_out = adapter_fit["W_full"].shape
        n_ctrl = int(len(ctrl_idx))
        dof_unit_a = max(n_ctrl - p_dim, 1)
        dof_a = max(dof_unit_a * T_out, 6)

        if "resid_full" in adapter_fit:
            R_a = np.asarray(adapter_fit["resid_full"], dtype=np.float64)
            R_a = R_a - R_a.mean(axis=1, keepdims=True)
            sigma_t_mat = (R_a.T @ R_a) / float(dof_unit_a)
            sigma_t_mat = 0.5 * (sigma_t_mat + sigma_t_mat.T)
        else:
            J_t = np.eye(T_out, dtype=np.float64) - np.ones((T_out, T_out), dtype=np.float64) / float(max(T_out, 1))
            sigma_t_mat = float(adapter_fit["sigma2"]) * J_t

        try:
            chol_a = np.linalg.cholesky(gram_inv_a + 1e-12 * np.eye(p_dim))
        except np.linalg.LinAlgError:
            w_, v_ = np.linalg.eigh(gram_inv_a)
            chol_a = v_ * np.sqrt(np.maximum(w_, 0.0))[None, :]

        wt_a, vt_a = np.linalg.eigh(sigma_t_mat)
        chol_ta = vt_a * np.sqrt(np.maximum(wt_a, 0.0))[None, :]

        a_a = max(0.5 * float(dof_a), 2.5)
        b_a = a_a - 1.0
        for m in range(M):
            g_draw = max(float(rng.gamma(a_a, 1.0)), 1e-30)
            eta_m = np.sqrt(b_a / g_draw)
            Z_a = rng.standard_normal((p_dim, T_out))
            dW_a = eta_m * (chol_a @ Z_a @ chol_ta.T)
            corr = Hs_all @ dW_a
            corr = corr - corr.mean(axis=1, keepdims=True)
            adapter_dev[m] = w_a * corr

        cov_unit_a = (w_a ** 2) * (Hs_all @ gram_inv_a @ Hs_all.T)
        diag_unit_a = np.diag(cov_unit_a)
        diag_time_a = np.diag(sigma_t_mat)
        analytical_var_cell += np.outer(diag_unit_a, diag_time_a)
        if T_post > 0:
            v_post = np.zeros(T_out, dtype=np.float64)
            v_post[anchor_t:T] = 1.0 / float(T_post)
            post_mean_factor = float(v_post @ sigma_t_mat @ v_post)
            analytical_cov_post_mean = analytical_cov_post_mean + post_mean_factor * cov_unit_a

    dev = backbone_dev + adapter_dev
    if M > 1:
        dev = dev - dev.mean(axis=0, keepdims=True)

    var_backbone = float(np.mean(np.var(backbone_dev[:, :, post], axis=0, ddof=1))) if M > 1 else 0.0
    var_adapter = float(np.mean(np.var(adapter_dev[:, :, post], axis=0, ddof=1))) if M > 1 else 0.0

    # Diagnostic-only posterior-predictive check on never-treated post-anchor cells (no inflation)
    calib: dict[str, Any] = {"applied": False, "kappa": 1.0}
    ctrl = np.asarray(ctrl_idx, dtype=int)
    if M > 1 and T_post >= 2 and ctrl.size >= 6:
        y_ctl = np.asarray(y_safe, dtype=np.float64)[ctrl][:, post]
        if np.all(np.isfinite(y_ctl)):
            def _dd(a: np.ndarray) -> np.ndarray:
                return a - a.mean(axis=0, keepdims=True) - a.mean(axis=1, keepdims=True) + a.mean()

            n_c = int(ctrl.size)
            v_tot = float(np.mean(_dd(y_ctl - unfiltered_forecast[ctrl][:, post]) ** 2))
            v_est = float(np.mean(np.var(np.stack([_dd(dev[m][ctrl][:, post]) for m in range(M)]), axis=0, ddof=1)))
            core_ = fit_core(np.asarray(context, dtype=np.float64), T, anchor_t, guard=lc_spec.guard)
            ssr0 = float(np.sum(core_.ctx_resid[:, : core_.core_end] ** 2))
            dof0 = max(core_.N * core_.core_end - core_.N * int(core_.coef.shape[0]), core_.N)
            sigma2_dd = (ssr0 / dof0) * (1.0 - 1.0 / n_c) * (1.0 - 1.0 / T_post)
            calib = {
                "applied": False,
                "kappa": 1.0,
                "control_error_var_dd": v_tot,
                "noise_var_dd": sigma2_dd,
                "draw_var_dd": v_est,
                "n_controls": n_c,
            }

    draws = unfiltered_forecast[None, :, :] + dev
    emp_post_var_mean = float(np.mean(np.var(dev[:, :, post], axis=0, ddof=1))) if M > 1 else 0.0
    ana_post_var_mean = float(np.mean(analytical_var_cell[:, post]))
    diag = {
        "draw_method": "posterior",
        "n_draws": M,
        "block_length": 0,
        "calibration": calib,
        "post_anchor_sd_mean": float(np.mean(np.std(dev[:, :, post], axis=0, ddof=1))) if M > 1 else 0.0,
        "pre_anchor_sd_mean": (float(np.mean(np.std(dev[:, :, :anchor_t], axis=0, ddof=1))) if (M > 1 and anchor_t > 0) else 0.0),
        "adapter_variance_share": (var_adapter / (var_backbone + var_adapter)) if (var_backbone + var_adapter) > 0 else 0.0,
        "adapter_sigma2": (float(adapter_fit["sigma2"]) if adapter_fit is not None else None),
        "analytical_post_anchor_var_mean": ana_post_var_mean,
        "empirical_post_anchor_var_mean": emp_post_var_mean,
        "analytical_cov_post_mean": analytical_cov_post_mean,
    }
    return draws, diag


def build_causal_features(
    y: np.ndarray,
    z: np.ndarray,
    x: np.ndarray | None = None,
    *,
    y_pre_history: np.ndarray | None = None,
    forecaster: Forecaster | str | None = None,
    n_folds: int = 5,
    seed: int = 0,
    late_context: LateContextSpec | str | None = None,
    n_forecast_draws: int = 0,
    forecast_draw_block_length: int | None = None,
    draw_method: str = "bootstrap",
    n_factor_loadings: int = 3,
) -> CausalFeatures:
    """Build leak-free counterfactual baseline forecasts, trajectory summaries, and a selection score.

    Enforces three hard causal and statistical invariants by construction:

    1. **Strict Pre-Rollout Context:** Every unit's temporal context window ends
       at the common pre-treatment cutoff ``anchor_t = min G_i`` (the earliest
       period in which any unit is treated), where ``np.all(z[:, :anchor_t] == 0)``.
    2. **Post-Treatment Masking & Out-of-Fold Control Cross-Learning:** Before
       calling ``forecaster`` or fitting the $K$-fold control-group trajectory
       adapter, all post-``anchor_t`` entries of ever-treated units are masked
       with ``NaN``, guaranteeing exact invariance to ``y[z == 1]``.
    3. **Automatic Pre-Treatment Gating:** ``baseline_forecast``,
       ``recommended_unit_features``, and ``recommended_ps`` are active only
       when the pre-treatment context has at least ``L_ctx >= 12`` periods and
       ``beyond_poly_r2 >= 0.05``; otherwise ``feats.gate_open`` is ``False``
       and passing ``x=feats.augment_x(x), x_trt=x, ps=feats.recommended_ps,
       baseline_forecast=feats.baseline_forecast`` automatically reverts to
       default ``LongBet``.

    Parameters
    ----------
    y
        Outcome panel of shape ``(N, T)``.
    z
        Binary absorbing treatment panel of shape ``(N, T)``.
    x
        Optional static baseline covariates of shape ``(N, P)``.
    y_pre_history
        Optional strictly pre-treatment historical outcome panel of shape
        ``(N, T_hist)`` observed prior to period 0 of ``y``. A longer
        pre-rollout history (``L_ctx = T_hist + anchor_t >= 12``, ideally
        ``>= 24``) is what makes seasonal and momentum identification possible.
    forecaster
        ``None`` or ``"patch_ridge"`` (default): the built-in guarded-core
        spectral forecaster with cross-series patch-ridge refinement,
        :func:`patch_ridge_forecaster` -- the backbone used throughout the
        companion paper.  ``"spectral"``: the spectral backbone alone.  Or a
        callable ``fn(context_2d, horizon) -> np.ndarray`` of shape
        ``(N, horizon)`` (an optional ``late_context`` keyword receives the
        resolved :class:`~longbet.LateContextSpec`); the callable sees only
        the strictly pre-anchor context, so any forecaster can be plugged in
        without weakening the firewall.
    n_folds
        Number of folds for out-of-fold control-group cross-learning and
        propensity estimation.
    seed
        Random seed for deterministic fold assignment.
    late_context
        How the forecaster treats the last pre-anchor periods: ``"fixed"``
        (default; a slope-break hinge kinked 10 periods before the first
        adoption and a Gaussian pre-adoption dip of width 1.8 centred 1.5
        periods before the first post-anchor period, an analyst-set
        identification restriction) or ``"none"``.  A
        :class:`~longbet.LateContextSpec` moves either kernel
        (``break_kink``, ``transient_centre``, ``transient_width``) or switches
        it off, which is how the sensitivity of an estimate to that restriction
        is reported; :func:`longbet.late_context_evidence` gives the data's
        (descriptive) evidence on the placement.  There is no data-adaptive
        mode: see :mod:`longbet._latecontext` for why.
    n_forecast_draws
        If ``> 0`` (it must then be ``>= 2``), also return
        ``baseline_forecast_draws`` and ``point_forecast_draws`` -- ``M`` draws
        from the forecaster's posterior/predictive law, centred on
        ``point_forecast``.  Default ``0`` leaves every other output
        bit-identical.  Not needed for :func:`~longbet.fit_rank_k_longbet`: the
        Rank-K fit draws its ``max_forecast_draws`` posterior forecast panels on
        demand (:meth:`CausalFeatures.draw_posterior_forecasts`) and reuses draws
        built here only when ``draw_method="posterior"`` and ``M`` is at least
        ``max_forecast_draws``.
    forecast_draw_block_length
        Block length of the moving-block bootstrap when ``draw_method="bootstrap"``;
        default ``round(L_ctx ** (1/3)) + 1``.
    draw_method
        ``"bootstrap"`` (default; legacy moving-block bootstrap whose
        cell-level channel is inflated by a control-error factor ``kappa`` in
        ``[1, 3]``) or ``"posterior"`` (conjugate Normal-inverse-gamma
        posterior draws of the patch-ridge and control-adapter weights, with no
        ``kappa`` inflation; the law the Rank-K fit uses).  Only used when
        ``n_forecast_draws > 0``.
    n_factor_loadings
        Number of pre-anchor residual singular vectors ``lambda_i`` to extract
        into :attr:`CausalFeatures.pre_anchor_loadings` (default ``3``).

    Returns
    -------
    CausalFeatures
        Frozen dataclass with ``(baseline_forecast, point_forecast,
        unit_features, ps, anchor_t, beyond_poly_r2, pre_anchor_loadings)``,
        gate-aware properties ``(gate_open, recommended_unit_features,
        recommended_ps, augment_x, point_forecast_draws)``, and optional
        ``baseline_forecast_draws`` / ``forecast_draws`` diagnostics.
    """
    n_draws_int = int(n_forecast_draws)
    if n_draws_int != 0 and n_draws_int < 2:
        raise ValueError(f"n_forecast_draws must be 0 or >= 2, got {n_forecast_draws!r}")
    if draw_method not in ("bootstrap", "posterior"):
        raise ValueError(f"draw_method must be 'bootstrap' or 'posterior', got {draw_method!r}")
    forecaster_fn = _resolve_forecaster(forecaster)

    # Strictly pre-treatment context, common anchor, control indices and the hard causal
    # firewall (every post-anchor cell of an ever-treated unit masked with NaN)
    context, anchor_t, ctrl_idx, y_safe = assemble_context(y, z, y_pre_history)
    N, T = y_safe.shape
    z_arr = np.asarray(z)

    # Step 1: Zero-shot temporal forecast across the panel window
    lc_spec = resolve_late_context(late_context)
    zero_shot, ctx_features, ctx_resid = _zero_shot_stage(context, T, anchor_t, lc_spec, forecaster_fn)
    pre_anchor_loadings = extract_pre_anchor_loadings(ctx_resid, q=int(n_factor_loadings))

    # Step 2: Control-group out-of-fold cross-learning adapter
    rng = np.random.default_rng(seed)
    if x is not None:
        x_arr = np.asarray(x, dtype=np.float64)
        if x_arr.ndim == 1:
            x_arr = x_arr[:, None]
        if x_arr.shape[0] != N:
            raise ValueError(f"x must have N={N} rows, got {x_arr.shape}")
        H_embed = np.hstack([x_arr, ctx_features])
    else:
        x_arr = np.zeros((N, 0), dtype=np.float64)
        H_embed = ctx_features

    adapter_fit: dict[str, Any] | None = None
    if len(ctrl_idx) >= 6:
        y_ctrl = y_safe[ctrl_idx].copy()
        if not np.all(np.isfinite(y_ctrl)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                col_means = np.nanmean(y_ctrl, axis=0, keepdims=True)
            col_means = np.where(np.isfinite(col_means), col_means, 0.0)
            y_ctrl = np.where(np.isfinite(y_ctrl), y_ctrl, col_means)
        ctx_mean_ctrl = np.mean(context[ctrl_idx], axis=1, keepdims=True)
        unfiltered_forecast, adapter_fit = _cross_fit_control_adapter(
            y_ctrl=y_ctrl,
            ctx_mean_ctrl=ctx_mean_ctrl,
            H_all=H_embed,
            ctrl_idx=ctrl_idx,
            zero_shot_all=zero_shot,
            n_folds=n_folds,
            rng=rng,
            return_fit=True,
        )
    else:
        unfiltered_forecast = zero_shot.copy()

    # Step 3: Extract unit-constant trajectory summaries u_i = [drift_or_u1, seasonal_sd]
    guard = 4 if context.shape[1] >= 12 else 0
    core_end = max(context.shape[1] - guard, 1)
    ctx_mean_all = np.mean(context[:, :core_end], axis=1)
    post_slice = unfiltered_forecast[:, anchor_t:] if anchor_t < T else unfiltered_forecast
    pred_drift = np.mean(post_slice, axis=1) - ctx_mean_all
    # Center each period across units before measuring unit-specific cyclical amplitude
    post_centered = post_slice - np.mean(post_slice, axis=0, keepdims=True)
    pred_amp = np.std(post_centered, axis=1)
    u1_feat = ctx_features[:, 2]
    primary_u = u1_feat if np.std(u1_feat) > 1e-6 else pred_drift

    raw_u = np.column_stack([primary_u, pred_amp])
    u_mean = raw_u.mean(axis=0, keepdims=True)
    u_scale = raw_u.std(axis=0, keepdims=True)
    u_scale = np.where(u_scale < 1e-8, 1.0, u_scale)
    unit_features = ((raw_u - u_mean) / u_scale).astype(np.float32)

    # Step 4: Cross-fitted scalar propensity score combining x and unit_features
    ps_features = np.hstack([x_arr, unit_features.astype(np.float64)])
    ps_hat = _cross_fit_propensity(ps_features, z_arr, n_folds=n_folds, rng=rng)

    # Step 5: Out-of-fold non-polynomial validation R^2 & automatic pre-treatment gating
    r2_rng = np.random.default_rng(int(seed) + 7919)
    beyond_poly_r2 = _compute_beyond_poly_r2(
        y_safe=y_safe,
        ctrl_idx=ctrl_idx,
        baseline_forecast=unfiltered_forecast,
        context=context,
        anchor_t=anchor_t,
        n_folds=n_folds,
        rng=r2_rng,
    )
    unfiltered_f32 = unfiltered_forecast.astype(np.float32)
    is_active = (
        context.shape[1] >= MIN_CONTEXT_LENGTH
        and beyond_poly_r2 >= MIN_BEYOND_POLY_R2
    )
    gated_forecast = unfiltered_f32 if is_active else None

    # Step 6 (optional): draws from the forecaster's predictive/posterior law for a cut posterior
    draws_f32: np.ndarray | None = None
    unfiltered_draws_f32: np.ndarray | None = None
    draws_diag: dict[str, Any] | None = None
    if n_draws_int >= 2:
        draw_rng = np.random.default_rng(int(seed) + 9973)
        if draw_method == "posterior":
            draws, draws_diag = _forecast_draws_posterior(
                n_draws=n_draws_int,
                context=context,
                unfiltered_forecast=unfiltered_forecast,
                adapter_fit=adapter_fit,
                T=T,
                anchor_t=anchor_t,
                lc_spec=lc_spec,
                forecaster=forecaster_fn,
                y_safe=y_safe,
                ctrl_idx=ctrl_idx,
                rng=draw_rng,
            )
        else:
            draws, draws_diag = _forecast_draws(
                n_draws=n_draws_int,
                block_length=forecast_draw_block_length,
                context=context,
                ctx_resid=ctx_resid,
                zero_shot=zero_shot,
                unfiltered_forecast=unfiltered_forecast,
                adapter_fit=adapter_fit,
                T=T,
                anchor_t=anchor_t,
                lc_spec=lc_spec,
                forecaster=forecaster_fn,
                y_safe=y_safe,
                ctrl_idx=ctrl_idx,
                rng=draw_rng,
            )
        draws_diag["gate_open"] = bool(is_active)
        unfiltered_draws_f32 = draws.astype(np.float32)
        if is_active:
            draws_f32 = unfiltered_draws_f32

    post_draw_state = {
        "context": context,
        "unfiltered_forecast": unfiltered_forecast,
        "adapter_fit": adapter_fit,
        "T": int(T),
        "anchor_t": int(anchor_t),
        "lc_spec": lc_spec,
        "forecaster": forecaster_fn,
        "y_safe": y_safe,
        "ctrl_idx": ctrl_idx,
        "seed": int(seed),
    }

    return CausalFeatures(
        baseline_forecast=gated_forecast,
        unit_features=unit_features,
        ps=ps_hat,
        anchor_t=int(anchor_t),
        beyond_poly_r2=beyond_poly_r2,
        point_forecast=unfiltered_f32,
        baseline_forecast_draws=draws_f32,
        forecast_draws=draws_diag,
        _unfiltered_forecast_draws=unfiltered_draws_f32,
        pre_anchor_loadings=pre_anchor_loadings,
        _posterior_draw_state=post_draw_state,
    )


def _compute_beyond_poly_r2(
    y_safe: np.ndarray,
    ctrl_idx: np.ndarray,
    baseline_forecast: np.ndarray,
    context: np.ndarray,
    anchor_t: int,
    degree: int = 2,
    n_folds: int = 5,
    rng: np.random.Generator | None = None,
) -> float:
    """Out-of-fold $R^2$ of ``baseline_forecast`` orthogonal to polynomial trends.

    Returns ``0.0`` when the pre-treatment context window has fewer than
    ``MIN_CONTEXT_LENGTH`` (12) periods. Otherwise measures the fraction of
    post-anchor untreated variance explained by ``baseline_forecast`` after
    double-centering and projecting out a degree-``degree`` polynomial in
    calendar time for every unit.
    """
    if context.shape[1] < MIN_CONTEXT_LENGTH:
        return 0.0
    if rng is None:
        rng = np.random.default_rng(0)

    T = y_safe.shape[1]
    if len(ctrl_idx) >= 4 and T >= 4:
        y_mat = y_safe[ctrl_idx].copy()
        if not np.all(np.isfinite(y_mat)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                col_means = np.nanmean(y_mat, axis=0, keepdims=True)
            col_means = np.where(np.isfinite(col_means), col_means, 0.0)
            y_mat = np.where(np.isfinite(y_mat), y_mat, col_means)
        f_mat = baseline_forecast[ctrl_idx].astype(np.float64)
        post_start = anchor_t if anchor_t < T else 0
    elif context.shape[1] >= 6:
        L = context.shape[1]
        val_h = min(6, max(2, L // 4))
        sub_pred, _ = _spectral_damped_zero_shot(context[:, : L - val_h], L, L - val_h)
        y_mat = context.copy()
        f_mat = sub_pred
        post_start = L - val_h
    else:
        return 0.0

    n_u, n_t = y_mat.shape
    y_dc = (
        y_mat
        - y_mat.mean(axis=1, keepdims=True)
        - y_mat.mean(axis=0, keepdims=True)
        + y_mat.mean()
    )
    f_dc = (
        f_mat
        - f_mat.mean(axis=1, keepdims=True)
        - f_mat.mean(axis=0, keepdims=True)
        + f_mat.mean()
    )

    deg_eff = min(int(degree), max(n_t - 2, 0))
    t_grid = np.linspace(-1.0, 1.0, n_t, dtype=np.float64)
    V = np.vander(t_grid, N=deg_eff + 1, increasing=True)
    Q, _ = np.linalg.qr(V)

    y_perp = y_dc - (y_dc @ Q) @ Q.T
    f_perp = f_dc - (f_dc @ Q) @ Q.T

    y_post = y_perp[:, post_start:]
    f_post = f_perp[:, post_start:]
    if y_post.size == 0:
        return 0.0

    k_eff = min(max(int(n_folds), 2), n_u)
    perm = rng.permutation(n_u)
    folds = np.array_split(perm, k_eff)
    pred_oof = np.zeros_like(y_post)
    for val_idx in folds:
        tr_mask = np.ones(n_u, dtype=bool)
        tr_mask[val_idx] = False
        num = float(np.sum(f_post[tr_mask] * y_post[tr_mask]))
        den = float(np.sum(f_post[tr_mask] ** 2)) + 1e-8
        slope = float(np.clip(num / den, 0.0, 1.5))
        pred_oof[val_idx] = slope * f_post[val_idx]

    null_mse = float(np.mean(y_post ** 2))
    oof_mse = float(np.mean((y_post - pred_oof) ** 2))
    return float(np.clip(1.0 - oof_mse / max(null_mse, 1e-12), 0.0, 1.0))


__all__ = [
    "CausalFeatures",
    "Forecaster",
    "build_causal_features",
    "patch_ridge_forecaster",
]

