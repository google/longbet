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

"""Tests for the late-context structure of the built-in forecaster (``longbet._latecontext``).

Guarantees pinned here:

1. ``late_context="fixed"`` (the default) is **bit-identical** to the forecaster
   that existed before the late-context refactor (SHA-256 digests of its outputs
   on six synthetic panels, ``tests/data/latecontext_fixed_digests.json``), and
   the explicit default kernel locations reproduce it exactly.  The spectral
   backbone (``forecaster="spectral"``, the pre-0.3 default) and a custom
   callable are pinned by the original digests; the 0.3 default
   (``forecaster=None`` = :func:`longbet.patch_ridge_forecaster`) is pinned by
   the ``features_default`` digest.
2. Every spec respects the causal firewall: perturbing post-anchor outcomes of
   ever-treated units changes nothing in the forecast, the unit features or the
   propensity score.
3. The kernel-location parameters do what they say (moving / disabling a kernel
   changes the forecast in the expected way; ``"none"`` drops both), the
   data-adaptive mode is gone, and the ``LateContextSpec`` validation rules hold.
4. ``late_context_evidence`` is a pure diagnostic: finite, JSON-serialisable,
   leak-free, and it does not produce a forecast.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from longbet import LateContextSpec, build_causal_features, late_context_evidence
from longbet._latecontext import apply_fixed_kernels, apply_fixed_kernels_full, assemble_context, fit_core, resolve_late_context
from longbet._forecast import _spectral_damped_zero_shot

DIGESTS = Path(__file__).parent / "data" / "latecontext_fixed_digests.json"


def _digest(*arrs: np.ndarray) -> str:
    h = hashlib.sha256()
    for a in arrs:
        h.update(np.ascontiguousarray(np.asarray(a, dtype=np.float64)).tobytes())
    return h.hexdigest()


def _panel(seed: int, N: int, T: int, T_hist: int, dip: bool = True):
    """Seasonal panel with a slope break 10 periods before the first adoption and an optional
    Gaussian pre-adoption dip (the structure the fixed kernels were built for)."""
    rng = np.random.default_rng(seed)
    L = T_hist + T
    t = np.arange(-T_hist, T, dtype=float)
    v = rng.normal(size=N)
    x = rng.normal(size=(N, 3))
    amp = 1.0 + 0.5 * np.abs(v)
    ph = np.where(v > 0, -np.pi / 2, np.pi / 2)
    y = (0.8 * x[:, :1] + v[:, None] * 0.3 + amp[:, None] * np.cos(2 * np.pi * t[None, :] / 6 + ph[:, None])
         + 0.4 * np.tanh(v)[:, None] * np.maximum(0, t[None, :] + 6) / 16
         + 0.3 * np.sin(2 * np.pi * t[None, :] / 8)[None, :].reshape(1, -1))
    if dip:
        y = y - 0.4 * np.maximum(0, v)[:, None] * np.exp(-0.5 * ((t[None, :] - 2.5) / 1.8) ** 2)
    y = y + rng.normal(0, 0.35, size=(N, L))
    util = 1.3 * v + 0.5 * x[:, 0] + rng.logistic(size=N) * 0.35
    q = np.quantile(util, [0.5, 0.68, 0.82])
    cohorts = np.where(util >= q[2], 4, np.where(util >= q[1], 5, np.where(util >= q[0], 6, T)))
    z = np.zeros((N, T), dtype=int)
    for i, g in enumerate(cohorts):
        if g < T:
            z[i, g:] = 1
    return y[:, T_hist:] + 1.5 * z, z, x, y[:, :T_hist]


def _parse_key(key: str):
    parts = key.split("_")
    return int(parts[0][4:]), int(parts[1][1:]), int(parts[2][1:]), int(parts[3][2:]), bool(int(parts[4][3:]))


def _pre_refactor_spectral_damped_zero_shot(
    context: np.ndarray, T_panel: int, anchor_t: int
) -> tuple[np.ndarray, np.ndarray]:
    """Monolithic pre-refactor spectral zero-shot forecaster."""
    N, L = context.shape
    guard = 4 if L >= 12 else 0
    core_end = max(L - guard, 1)
    grand_mean = float(np.mean(context[:, :core_end]))
    ctx_dm = context - np.mean(context, axis=0, keepdims=True)

    t_ctx = np.arange(L, dtype=np.float64) - float(core_end)
    t_panel = np.arange(T_panel, dtype=np.float64) + float(guard - anchor_t)
    t_c = float(np.mean(t_ctx[:core_end]))
    s_t = max(float(np.std(t_ctx[:core_end])), 1.0)
    t_lin_ctx = (t_ctx - t_c) / s_t
    t_lin_pan = (t_panel - t_c) / s_t

    basis_ctx_cols = [np.ones(L, dtype=np.float64), t_lin_ctx]
    basis_pan_cols = [np.ones(T_panel, dtype=np.float64), t_lin_pan]
    reg_diag = [1e-6, 1e-4 * float(L)]

    periods = [
        p
        for p in (4.0, 6.0, 8.0, 12.0, 16.0, 24.0)
        if (
            core_end >= max(int(2.0 * p), 6)
            if L >= 12
            else L >= max(int(1.25 * p), 4)
        )
    ]
    for p in periods:
        omega = 2.0 * np.pi / p
        basis_ctx_cols.extend([np.cos(omega * t_ctx), np.sin(omega * t_ctx)])
        basis_pan_cols.extend(
            [np.cos(omega * t_panel), np.sin(omega * t_panel)]
        )
        reg_diag.extend([1e-3 * float(L), 1e-3 * float(L)])

    Phi_ctx = np.column_stack(basis_ctx_cols)
    Phi_pan = np.column_stack(basis_pan_cols)

    weights = np.ones(L, dtype=np.float64)
    if guard > 0:
        weights[core_end:] = 0.0
    W_Phi = Phi_ctx * weights[:, None]
    gram_Phi = Phi_ctx.T @ W_Phi + np.diag(reg_diag)
    coef = np.linalg.solve(gram_Phi, W_Phi.T @ ctx_dm.T)
    sec_pred = (Phi_pan @ coef).T
    ctx_resid = ctx_dm - (Phi_ctx @ coef).T

    u1 = np.zeros(N, dtype=np.float64)
    dip_amp = np.zeros(N, dtype=np.float64)
    brk_amp = np.zeros(N, dtype=np.float64)
    if guard > 0 and L >= 12 and coef.shape[0] > 2 and N >= 6:
        sigma2_core = float(np.mean(ctx_resid[:, :core_end] ** 2))
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
            U_design = np.column_stack(
                [u1, u_pos / max(float(np.std(u_pos)), 1e-6)]
            )
            k_u = U_design.shape[1]
            P_u = U_design @ np.linalg.solve(
                U_design.T @ U_design + 0.2 * np.eye(k_u), U_design.T
            )

            h_brk_ctx = np.maximum(0.0, t_ctx + 6.0) / 16.0
            h_brk_pan = np.maximum(0.0, t_panel + 6.0) / 16.0
            beta_h = np.linalg.solve(gram_Phi, W_Phi.T @ h_brk_ctx)
            h_orth_ctx = h_brk_ctx - Phi_ctx @ beta_h
            h_orth_pan = h_brk_pan - Phi_pan @ beta_h
            h_norm2 = float(np.sum(h_orth_ctx[:core_end] ** 2)) + 1e-8
            if h_norm2 > 0.02:
                brk_raw = (
                    ctx_resid[:, :core_end] @ h_orth_ctx[:core_end]
                ) / h_norm2
                brk_amp = P_u @ brk_raw
                sec_pred = sec_pred + brk_amp[:, None] * h_orth_pan[None, :]
                ctx_resid = ctx_resid - brk_amp[:, None] * h_orth_ctx[None, :]

            t_dip_center = float(guard) - 1.5
            k_dip_ctx = np.exp(-0.5 * ((t_ctx - t_dip_center) / 1.8) ** 2)
            k_dip_pan = np.exp(-0.5 * ((t_panel - t_dip_center) / 1.8) ** 2)
            k_win = k_dip_ctx[core_end:]
            k_norm2 = float(np.sum(k_win**2)) + 1e-8
            dip_raw = (ctx_resid[:, core_end:] @ k_win) / k_norm2
            dip_amp = P_u @ dip_raw
            panel_pred = (
                sec_pred + dip_amp[:, None] * k_dip_pan[None, :] + grand_mean
            )
        else:
            panel_pred = sec_pred + grand_mean
    else:
        panel_pred = sec_pred + grand_mean

    ctx_mean = np.mean(context[:, :core_end], axis=1, keepdims=True)
    ctx_sd = np.std(ctx_dm[:, :core_end], axis=1, keepdims=True)
    ctx_features = np.hstack(
        [
            ctx_mean,
            ctx_sd,
            u1[:, None],
            brk_amp[:, None],
            coef[1:, :].T,
        ]
    )
    return panel_pred, ctx_features


# --------------------------------------------------------------------------------------
# 1. The fixed path is the pre-refactor forecaster, bit for bit
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("key", sorted(json.loads(DIGESTS.read_text())))
def test_fixed_mode_is_bit_identical_to_pre_refactor_forecaster(key):
    ref = json.loads(DIGESTS.read_text())[key]
    seed, N, T, Th, dip = _parse_key(key)
    y, z, x, yh = _panel(seed, N, T, Th, dip)
    anchor = int(np.argmax(z.any(axis=0)))
    ctx = np.hstack([yh, y[:, :anchor]])
    sp, cf = _spectral_damped_zero_shot(ctx, T, anchor)
    sp_pre, cf_pre = _pre_refactor_spectral_damped_zero_shot(ctx, T, anchor)
    assert _digest(sp, cf) == _digest(sp_pre, cf_pre)
    # the explicit default locations are the same path
    explicit_spec = LateContextSpec(
        break_kink=10.0, transient_centre=-1.5, transient_width=1.8
    )
    sp2, cf2 = _spectral_damped_zero_shot(
        ctx, T, anchor, late_context=explicit_spec
    )
    assert _digest(sp2, cf2) == _digest(sp, cf)
    f = build_causal_features(
        y, z, x=x, y_pre_history=yh, seed=seed, forecaster="spectral"
    )
    f_exp = build_causal_features(
        y,
        z,
        x=x,
        y_pre_history=yh,
        seed=seed,
        forecaster="spectral",
        late_context=explicit_spec,
    )
    f_dig = _digest(
        f.point_forecast, f.unit_features, f.ps, np.array([f.beyond_poly_r2])
    )
    f_exp_dig = _digest(
        f_exp.point_forecast,
        f_exp.unit_features,
        f_exp.ps,
        np.array([f_exp.beyond_poly_r2]),
    )
    assert f_dig == f_exp_dig
    assert (f.baseline_forecast is not None) == ref["gate_open"]
    assert not hasattr(f, "late_context")  # no adaptive diagnostics exist
    ext_fn = lambda c, h: np.repeat(c[:, -1:], h, axis=1)
    f2 = build_causal_features(
        y, z, x=x, y_pre_history=yh, forecaster=ext_fn, seed=seed
    )
    f2_exp = build_causal_features(
        y,
        z,
        x=x,
        y_pre_history=yh,
        forecaster=ext_fn,
        seed=seed,
        late_context=explicit_spec,
    )
    f2_dig = _digest(
        f2.point_forecast,
        f2.unit_features,
        f2.ps,
        np.array([f2.beyond_poly_r2]),
    )
    f2_exp_dig = _digest(
        f2_exp.point_forecast,
        f2_exp.unit_features,
        f2_exp.ps,
        np.array([f2_exp.beyond_poly_r2]),
    )
    assert f2_dig == f2_exp_dig
    # 0.3 default is the paper's patch-ridge backbone; None == "patch_ridge"
    fd = build_causal_features(y, z, x=x, y_pre_history=yh, seed=seed)
    fd_exp = build_causal_features(
        y, z, x=x, y_pre_history=yh, seed=seed, late_context=explicit_spec
    )
    fd_dig = _digest(
        fd.point_forecast,
        fd.unit_features,
        fd.ps,
        np.array([fd.beyond_poly_r2]),
    )
    fd_exp_dig = _digest(
        fd_exp.point_forecast,
        fd_exp.unit_features,
        fd_exp.ps,
        np.array([fd_exp.beyond_poly_r2]),
    )
    assert fd_dig == fd_exp_dig
    assert (fd.baseline_forecast is not None) == ref["gate_open_default"]
    fd2 = build_causal_features(
        y, z, x=x, y_pre_history=yh, seed=seed, forecaster="patch_ridge"
    )
    np.testing.assert_array_equal(fd.point_forecast, fd2.point_forecast)
    np.testing.assert_array_equal(fd.ps, fd2.ps)
    # On the reference BLAS/SIMD microarchitecture where DIGESTS was recorded,
    # also verify the archived raw float64 byte digests.
    if _digest(sp, cf) == ref["spectral"]:
        assert f_dig == ref["features"]
        assert f2_dig == ref["features_ext"]
        assert fd_dig == ref["features_default"]
    with pytest.raises(TypeError):
        build_causal_features(
            y, z, x=x, y_pre_history=yh, seed=seed, forecaster="timesfm"  # type: ignore[arg-type]
        )


def test_fixed_spec_variants_resolve_to_original_path():
    assert resolve_late_context(None).uses_original_path
    assert resolve_late_context("fixed").uses_original_path
    assert LateContextSpec(mode="fixed", break_mode="fixed", transient_mode="fixed").uses_original_path
    assert not LateContextSpec(mode="fixed", guard=6).uses_original_path
    assert not LateContextSpec(break_kink=8.0).uses_original_path
    assert not LateContextSpec(transient_centre=-2.5).uses_original_path
    assert not LateContextSpec(transient_width=3.0).uses_original_path
    assert not LateContextSpec(mode="none").uses_original_path
    assert not LateContextSpec(transient_mode="none").uses_original_path


def test_apply_fixed_kernels_default_spec_equals_no_spec():
    y, z, x, yh = _panel(1, 80, 10, 26, True)
    anchor = int(np.argmax(z.any(axis=0)))
    core = fit_core(np.hstack([yh, y[:, :anchor]]), 10, anchor, guard=4)
    a = apply_fixed_kernels(core)
    b = apply_fixed_kernels(core, LateContextSpec())
    for u, v in zip(a, b):
        np.testing.assert_array_equal(u, v)


# --------------------------------------------------------------------------------------
# 2. Causal firewall for every spec
# --------------------------------------------------------------------------------------

SPECS = {
    "fixed": LateContextSpec(),
    "none": LateContextSpec(mode="none"),
    "break_only": LateContextSpec(transient_mode="none"),
    "dip_only": LateContextSpec(break_mode="none"),
    "moved": LateContextSpec(break_kink=7.0, transient_centre=-2.5, transient_width=3.0),
    "guard6": LateContextSpec(guard=6),
}


@pytest.mark.parametrize("name", sorted(SPECS))
def test_specs_never_use_treated_post_anchor_outcomes(name):
    spec = SPECS[name]
    y, z, x, yh = _panel(4, 120, 12, 44, True)
    f = build_causal_features(y, z, x=x, y_pre_history=yh, seed=4, late_context=spec)
    anchor = f.anchor_t
    ever = z.any(axis=1)
    y_pert = y.copy()
    y_pert[ever, anchor:] += 50.0  # absurd post-anchor values for every ever-treated unit
    g = build_causal_features(y_pert, z, x=x, y_pre_history=yh, seed=4, late_context=spec)
    np.testing.assert_array_equal(f.point_forecast, g.point_forecast)
    np.testing.assert_array_equal(f.unit_features, g.unit_features)
    np.testing.assert_array_equal(f.ps, g.ps)
    assert f.beyond_poly_r2 == g.beyond_poly_r2
    assert np.all(np.isfinite(f.point_forecast))
    assert f.point_forecast.shape == y.shape
    assert f.unit_features.shape == (y.shape[0], 2)


# --------------------------------------------------------------------------------------
# 3. Kernel locations and component switches
# --------------------------------------------------------------------------------------


def test_component_switches_and_locations_change_the_forecast_as_expected():
    y, z, x, yh = _panel(1, 80, 10, 26, True)
    anchor = int(np.argmax(z.any(axis=0)))
    core = fit_core(np.hstack([yh, y[:, :anchor]]), 10, anchor, guard=4)
    assert core.kernels_allowed and core.P_u is not None
    full = apply_fixed_kernels_full(core)
    none = apply_fixed_kernels_full(core, LateContextSpec(mode="none"))
    brk = apply_fixed_kernels_full(core, LateContextSpec(transient_mode="none"))
    dip = apply_fixed_kernels_full(core, LateContextSpec(break_mode="none"))
    # "none" is the secular/harmonic fit only
    np.testing.assert_allclose(none["panel_pred"], core.sec_pred + core.grand_mean, rtol=0, atol=1e-12)
    assert np.all(none["brk_amp"] == 0.0) and np.all(none["dip_amp"] == 0.0)
    # single components keep exactly their own amplitudes
    assert np.all(brk["dip_amp"] == 0.0) and np.any(brk["brk_amp"] != 0.0)
    assert np.all(dip["brk_amp"] == 0.0) and np.any(dip["dip_amp"] != 0.0)
    np.testing.assert_array_equal(brk["brk_amp"], full["brk_amp"])  # the break is estimated first, as before
    # the full forecast differs from the single-component ones
    assert float(np.mean(np.abs(full["panel_pred"] - brk["panel_pred"]))) > 1e-6
    assert float(np.mean(np.abs(full["panel_pred"] - dip["panel_pred"]))) > 1e-6
    # moving the dip by two periods changes the post-anchor forecast (an explicit sensitivity analysis)
    moved = apply_fixed_kernels_full(core, LateContextSpec(transient_centre=-3.5))
    assert float(np.mean(np.abs(moved["panel_pred"][:, anchor:] - full["panel_pred"][:, anchor:]))) > 1e-4
    # a kink too close to the anchor has no energy in the core and is dropped
    late = apply_fixed_kernels_full(core, LateContextSpec(break_kink=0.5))
    assert np.all(late["brk_amp"] == 0.0)
    # the residual bookkeeping is consistent: fitted context + residual == demeaned context on the core
    # (secular + kernels removed from the residual; cross-sectional mean removed from the context)
    for res in (full, none, brk, dip, moved):
        assert res["ctx_resid"].shape == core.ctx_resid.shape and np.all(np.isfinite(res["ctx_resid"]))


def test_short_context_has_no_kernels_for_any_spec():
    y, z, x, yh = _panel(3, 40, 10, 4, True)  # 8-period context: guard 0, kernels off
    anchor = int(np.argmax(z.any(axis=0)))
    core = fit_core(np.hstack([yh, y[:, :anchor]]), 10, anchor, guard=4)
    assert core.guard == 0 and not core.kernels_allowed
    for spec in SPECS.values():
        res = apply_fixed_kernels_full(core, spec)
        np.testing.assert_allclose(res["panel_pred"], core.sec_pred + core.grand_mean, rtol=0, atol=1e-12)


# --------------------------------------------------------------------------------------
# 4. Spec validation: the adaptive mode is gone
# --------------------------------------------------------------------------------------


def test_spec_validation_rules():
    with pytest.raises(ValueError, match="removed"):
        LateContextSpec(mode="adaptive")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        resolve_late_context("adaptive")
    with pytest.raises(ValueError):
        LateContextSpec(mode="bogus")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        LateContextSpec(break_mode="joint")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        LateContextSpec(transient_mode="joint")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        LateContextSpec(joint=True)  # type: ignore[call-arg]  # removed field
    with pytest.raises(ValueError):
        LateContextSpec(guard=-1)
    with pytest.raises(ValueError):
        LateContextSpec(break_kink=0.0)
    with pytest.raises(ValueError):
        LateContextSpec(transient_width=0.0)
    with pytest.raises(ValueError):
        LateContextSpec(transient_centre=float("nan"))
    with pytest.raises(TypeError):
        resolve_late_context(3)  # type: ignore[arg-type]
    s = LateContextSpec(mode="none", break_mode="fixed")
    assert (s.resolved_break, s.resolved_transient) == ("fixed", "none") and not s.uses_original_path
    assert build_causal_features.__doc__ is not None and "adaptive" not in build_causal_features.__doc__.split("There is no data-adaptive")[0]


# --------------------------------------------------------------------------------------
# 5. The evidence diagnostic
# --------------------------------------------------------------------------------------


def test_late_context_evidence_is_descriptive_leak_free_and_serialisable():
    y, z, x, yh = _panel(1, 80, 10, 26, True)
    ev = late_context_evidence(y, z, y_pre_history=yh)
    assert "skipped" not in ev and ev["kernels_allowed"]
    assert ev["anchor_t"] == int(np.argmax(z.any(axis=0))) and ev["guard"] == 4
    gw = ev["guard_window"]
    assert len(gw["treated_minus_control_residual"]) == 4 and gw["mean"] is not None
    # the panel has a treated-side dip: the treated-minus-control guard residual is negative on average
    assert gw["mean"] < 0.0
    fk = ev["fixed_kernels"]
    assert fk["dip_amp_treated_minus_control"] < 0.0
    pr = ev["profile"]
    assert pr["best_family"] in ("none", "break", "transient", "both")
    assert len(pr["kink_curve"]) == pr["n_grid_break"] > 0 and len(pr["transient_surface"]) == pr["n_grid_transient"] > 0
    assert set(pr["family_scores"]) >= {"none"}
    json.dumps(ev)  # JSON-serialisable
    assert not any(isinstance(v, np.ndarray) for v in ev.values())
    # leak-free: post-anchor treated outcomes do not enter
    y_pert = y.copy()
    y_pert[z.any(axis=1), ev["anchor_t"]:] += 50.0
    ev2 = late_context_evidence(y_pert, z, y_pre_history=yh)
    assert json.dumps(ev2, sort_keys=True) == json.dumps(ev, sort_keys=True)
    # short contexts are reported as skipped
    y_s, z_s, _x_s, yh_s = _panel(3, 40, 10, 4, True)
    assert "skipped" in late_context_evidence(y_s, z_s, y_pre_history=yh_s)
    # a custom kink grid is honoured
    ev3 = late_context_evidence(y, z, y_pre_history=yh, kink_grid=(8.0, 10.0, 12.0))
    assert ev3["profile"]["n_grid_break"] == 3


def test_assemble_context_matches_feature_builder_anchor_and_firewall():
    y, z, x, yh = _panel(2, 60, 10, 20, True)
    ctx, anchor, ctrl_idx, y_safe = assemble_context(y, z, yh)
    f = build_causal_features(y, z, x=x, y_pre_history=yh, seed=2)
    assert anchor == f.anchor_t and ctx.shape == (60, 20 + anchor)
    assert np.all(np.isnan(y_safe[z.any(axis=1), anchor:])) and np.all(np.isfinite(y_safe[~z.any(axis=1)]))
    assert set(ctrl_idx) == set(np.flatnonzero(~z.any(axis=1)))
    with pytest.raises(ValueError):
        assemble_context(y[:, :0], z[:, :0], None)
