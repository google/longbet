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

"""The coupled model must carry the calendar-time trend block, not just the scalar one.

Why this file exists
--------------------
``LongBetMulti`` shares ``LongBetConfig`` with ``LongBet``. When the trend block
was added, :attr:`longbet.LongBetConfig.split_calendar_mu` became ``None`` by
default and resolves to ``not (use_trend_block and trend_period_effects)``, i.e.
to ``False`` with the shipped defaults: the prognostic forest is forbidden from
splitting on calendar time *because* the additive block now carries it.

The two halves are a package deal. A coupled fit that inherited the blinded
forest without also getting the block would have nothing at all able to
represent calendar time -- the measured failure mode, and a silent one, since
such a fit still runs and still returns numbers. Every test below fails loudly
on that state of the world:

* the block is actually constructed and reaches every child equation;
* it is the *same* block for every equation, since it depends only on ``(x, t)``;
* its coefficients are traced, so ``child.predict`` can use them;
* save/load round-trips it exactly;
* the sampler's own residual bookkeeping accounts for it.
"""

import numpy as np
import pytest

from longbet import LongBetConfig, LongBetMulti
from longbet._trend_basis import (
    MAX_RAW_COLS,
    build_period_basis,
    build_time_basis,
    build_trend_design,
    build_unit_features,
)

N, T, P = 14, 6, 3


def _panel(seed: int = 20260918):
    """A small staggered panel with a genuine calendar-time trend.

    The trend is common to both outcomes and is *not* a step function, so a
    model whose only time-varying prognostic component is a piecewise-constant
    forest cannot reproduce it. That is deliberate: it is the configuration the
    block exists for.
    """
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, P)).astype(np.float32)
    t = np.arange(1, T + 1)
    z = np.zeros((N, T), dtype=np.float32)
    z[:6, 3:] = 1.0
    z[6:10, 4:] = 1.0
    ramp = (t - 1.0) / (T - 1.0)
    base = x[:, 0][:, None] + 2.0 * ramp[None, :] * x[:, 1][:, None]
    tau = np.where(z == 1.0, 0.8, 0.0)
    y1 = (base + tau + 0.2 * rng.normal(size=(N, T))).astype(np.float32)
    y2 = (0.5 * base - tau + 0.2 * rng.normal(size=(N, T))).astype(np.float32)
    return x, z, {"y1": y1, "y2": y2}


def _config(**over):
    base = dict(
        sigma_prior_a=2.0,
        sigma_prior_b=1.0,
        num_burnin=3,
        num_sweeps=4,
        num_chains=2,
        num_trees_pr=3,
        num_trees_trt=2,
        max_depth_pr=3,
        max_depth_trt=3,
        min_points_per_leaf_pr=2,
        min_points_per_leaf_trt=2,
        random_seed=11,
        device="cpu",
    )
    base.update(over)
    return LongBetConfig(**base)


def _expected_q_tot(cfg):
    """Column count of the assembled design, derived from the config alone."""
    n_period = (T - 1) if cfg.trend_period_effects else 0
    q_feat = min(P, MAX_RAW_COLS) + cfg.trend_num_features
    return n_period + build_time_basis(T, cfg.trend_degree).shape[1] * q_feat


@pytest.fixture(scope="module")
def fitted():
    x, z, y = _panel()
    cfg = _config()
    return x, z, y, cfg, LongBetMulti(cfg).fit(y=y, x=x, z=z)


# -- construction ----------------------------------------------------------


def test_default_config_blinds_the_forest_and_therefore_needs_the_block():
    """The premise of every other test here: guard it explicitly."""
    cfg = _config()
    assert cfg.use_trend_block is True
    assert cfg.trend_period_effects is True
    # Resolved, not supplied: with the block on, the forest must not also be
    # allowed to split on calendar time.
    assert cfg.split_calendar_mu is False


def test_every_child_state_carries_the_trend_block(fitted):
    _, _, _, cfg, model = fitted
    assert model.state is not None
    assert len(model.state.states) == 2
    for m, st in enumerate(model.state.states):
        assert st.use_trend_block == cfg.use_trend_block, f"outcome {m}"
        assert st.trend_design.shape[1] > 0, f"outcome {m} has an empty trend design"
        assert st.trend_design.shape[1] == _expected_q_tot(cfg), f"outcome {m}"
        assert st.trend_num_period == T - 1, f"outcome {m}"


def test_trend_design_is_shared_across_outcomes(fitted):
    """Same ``x`` and same ``t`` must mean the same design, cell for cell.

    The panel is fully observed, so the per-outcome observed-cell masking that
    ``init_longbet`` applies is a no-op here and the arrays must be bit-equal.
    A per-outcome design would mean the coefficients of different equations
    were not comparable.
    """
    _, _, y, _, model = fitted
    assert all(np.isfinite(v).all() for v in y.values()), "fixture must be fully observed"
    first = np.asarray(model.state.states[0].trend_design)
    for m, st in enumerate(model.state.states[1:], start=1):
        np.testing.assert_array_equal(
            np.asarray(st.trend_design), first, err_msg=f"outcome {m}"
        )


def test_trend_design_matches_an_independent_reconstruction(fitted):
    """The model must use the documented basis, not merely *a* basis."""
    x, _, _, cfg, model = fitted
    unit_idx = np.repeat(np.arange(N, dtype=np.int32), T)
    time_idx = np.tile(np.arange(T, dtype=np.int32), N)
    phi = build_time_basis(T, cfg.trend_degree)
    period = build_period_basis(T)
    features, _ = build_unit_features(x, cfg.trend_num_features)
    expected = np.asarray(
        build_trend_design(features, phi, unit_idx, time_idx, period=period)
    )
    np.testing.assert_allclose(
        np.asarray(model.state.states[0].trend_design), expected, rtol=0, atol=0
    )
    # The fitted model must expose the same bases it sampled against.
    np.testing.assert_array_equal(np.asarray(model.trend_time_basis_), phi)
    np.testing.assert_array_equal(np.asarray(model.trend_period_basis_), period)
    for fit in model.fits:
        assert fit.trend_time_basis_ is model.trend_time_basis_
        assert fit.trend_period_basis_ is model.trend_period_basis_
        assert fit.trend_feature_spec_ is model.trend_feature_spec_


# -- tracing ---------------------------------------------------------------


def test_trend_coefficients_are_traced_per_outcome(fitted):
    _, _, _, cfg, model = fitted
    expected = (cfg.num_chains, cfg.num_sweeps, _expected_q_tot(cfg))
    for m in range(2):
        coef = np.asarray(model.trace.traces[m].trend_coef)
        assert coef.shape == expected, f"outcome {m}: {coef.shape}"
        assert np.isfinite(coef).all(), f"outcome {m} has non-finite trend draws"
        # An all-zero trace would mean the draws were allocated but never
        # written -- the bug this file guards against, one level down.
        assert np.any(coef != 0.0), f"outcome {m} trend trace was never written"
    # Each child fit sees its own equation's coefficients.
    for u, name in enumerate(model.outcome_names):
        internal = model.inverse_order[u]
        np.testing.assert_array_equal(
            np.asarray(model.fits[name].trace.trend_coef),
            np.asarray(model.trace.traces[internal].trend_coef),
        )


def test_block_off_gives_a_zero_width_trace_not_none():
    """Turning the block off must stay on the same code path, only narrower.

    ``q_tot == 0`` rather than ``None`` is what lets predict(), save/load and
    the diagnostics avoid a multi-specific special case.
    """
    x, z, y = _panel()
    cfg = _config(use_trend_block=False, split_calendar_mu=True)
    model = LongBetMulti(cfg).fit(y=y, x=x, z=z)
    assert model.trend_time_basis_ is None
    for m, st in enumerate(model.state.states):
        assert st.trend_design.shape[1] == 0, f"outcome {m}"
        assert st.use_trend_block is False, f"outcome {m}"
        coef = model.trace.traces[m].trend_coef
        assert coef is not None, f"outcome {m}"
        assert np.asarray(coef).shape == (cfg.num_chains, cfg.num_sweeps, 0)
    # And prediction still runs, without reconstructing a trend.
    model.predict(x=x, z=z)


# -- persistence -----------------------------------------------------------


def test_save_load_preserves_the_block_and_the_predictions(fitted, tmp_path):
    x, z, _, cfg, model = fitted
    before = model.predict(x=x, z=z)

    path = tmp_path / "multi_trend.npz"
    model.save(path)
    loaded = LongBetMulti.load(path)

    np.testing.assert_allclose(
        np.asarray(loaded.trend_time_basis_), np.asarray(model.trend_time_basis_)
    )
    np.testing.assert_allclose(
        np.asarray(loaded.trend_period_basis_), np.asarray(model.trend_period_basis_)
    )
    spec_a, spec_b = model.trend_feature_spec_, loaded.trend_feature_spec_
    assert spec_b is not None
    assert spec_b.num_raw == spec_a.num_raw
    for field in ("mean", "scale", "omega", "phase", "centre", "col_scale"):
        np.testing.assert_allclose(
            np.asarray(getattr(spec_b, field)),
            np.asarray(getattr(spec_a, field)),
            err_msg=field,
        )
    for m in range(2):
        np.testing.assert_array_equal(
            np.asarray(loaded.trace.traces[m].trend_coef),
            np.asarray(model.trace.traces[m].trend_coef),
        )

    after = loaded.predict(x=x, z=z)
    for name in model.outcome_names:
        np.testing.assert_array_equal(
            np.asarray(after[name].tauhats), np.asarray(before[name].tauhats)
        )
        np.testing.assert_array_equal(
            np.asarray(after[name].muhats0), np.asarray(before[name].muhats0)
        )


# -- the sampler's own bookkeeping ----------------------------------------


def test_residual_invariant_holds_for_every_child_state(fitted):
    """``state.resid`` must equal ``y - fitted`` with the trend block included.

    This is the check that would catch a block that is constructed and traced
    but not actually subtracted from the residual: the forests would then be
    fitting a surface the rest of the model had already explained. Recomputing
    the mean from scratch, in NumPy, is the only way to see it.

    Under SUR coupling the children store *raw* residuals at state boundaries
    (``longbet._multi_step``), so the invariant is per equation, exactly as in
    the scalar model.
    """
    _, _, _, _, model = fitted
    for m, st in enumerate(model.state.states):
        mask = np.asarray(st.obs_mask)
        unit_idx = np.asarray(st.unit_idx)
        exposure_idx = np.asarray(st.exposure_idx)
        z_vec = np.asarray(st.z_vec)

        b_z = np.where(z_vec[None, :] == 1.0,
                       np.asarray(st.b1)[:, None],
                       np.asarray(st.b0)[:, None])
        fitted_mean = (
            np.asarray(st.alpha)[:, None] * np.asarray(st.mu_fit)
            + np.asarray(st.trend_coef) @ np.asarray(st.trend_design).T
            + b_z * np.asarray(st.beta)[:, exposure_idx] * np.asarray(st.nu_fit)
            + np.asarray(st.gamma)[:, unit_idx]
        )
        response = np.asarray(st.y)
        recon = np.where(mask[None, :], response[None, :] - fitted_mean, 0.0)
        err = float(np.abs(recon - np.asarray(st.resid)).max())
        scale = float(np.abs(response[mask]).max())
        # float32 accumulation over q_tot + trees terms; tie the tolerance to
        # the data scale rather than hard-coding an absolute number.
        assert err < 1e-4 * max(scale, 1.0), f"outcome {m}: max residual error {err:.3e}"


def test_predictions_track_the_calendar_trend(fitted):
    """End to end: the fit must explain the panel it was given.

    With the forest blind to calendar time and no block, the fitted mean cannot
    move with ``t`` at all and this comparison is not close.
    """
    x, z, y, _, model = fitted
    pred = model.predict(x=x, z=z)
    for name in model.outcome_names:
        yhat = np.asarray(pred[name].y_summary.mean)
        resid = np.abs(np.asarray(y[name]) - yhat).mean()
        assert resid < 0.8 * float(np.asarray(y[name]).std()), (
            f"{name}: mean |y - yhat| = {resid:.3f}"
        )
