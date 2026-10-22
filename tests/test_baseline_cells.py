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

"""``baseline_cells='untreated'``: the two-module (cut) baseline.

The property that defines the option is that **no treated outcome can reach
the baseline module**. That is checked exactly, not statistically: with
standardization off, adding an arbitrary shift to every treated cell must
leave every baseline draw -- the prognostic forest's fit, the trend
coefficients, the unit intercepts and both variances -- bit-for-bit unchanged,
while the same experiment under ``'all'`` moves them.
"""

from __future__ import annotations

import dataclasses

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from longbet import LongBet, LongBetConfig
from longbet._trend_basis import trend_block_step


def _null_panel(seed=3, N=36, T=8, P=3, noise=0.3):
    """Staggered adoption, no treatment effect, a covariate-dependent trend."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, P)).astype(np.float32)
    t = np.arange(T, dtype=np.float32)
    cohort = np.where(rng.uniform(size=N) < 0.5, rng.integers(3, T - 1, size=N), 10**6)
    z = (t[None, :] >= cohort[:, None]).astype(np.float32)
    gamma = rng.normal(0.0, 0.5, size=N).astype(np.float32)
    y = (
        0.6 * x[:, :1]
        + 0.15 * t[None, :] * (1.0 + 0.5 * x[:, 1:2])
        + gamma[:, None]
        + rng.normal(0.0, noise, size=(N, T))
    ).astype(np.float32)
    return x, y, z


def _small_config(cells, **kw):
    base = dict(
        num_burnin=15, num_sweeps=10, num_chains=2, num_trees_pr=6, num_trees_trt=8,
        trend_num_features=8, random_seed=11, standardize=False, baseline_cells=cells,
    )
    base.update(kw)
    return LongBetConfig(**base)


# --------------------------------------------------------------------------
# configuration surface
# --------------------------------------------------------------------------


def test_config_validates_baseline_cells():
    assert LongBetConfig().baseline_cells == "all"
    with pytest.raises(ValueError, match="baseline_cells"):
        LongBetConfig(baseline_cells="treated")
    with pytest.raises(ValueError, match="continuous"):
        LongBetConfig(baseline_cells="untreated", outcome="binary")
    with pytest.raises(ValueError, match="tempering"):
        LongBetConfig(baseline_cells="untreated", tempering_levels=2)
    with pytest.raises(ValueError, match="ar1_errors"):
        LongBetConfig(baseline_cells="untreated", ar1_errors=True)
    cfg = LongBetConfig(baseline_cells="untreated")
    assert LongBetConfig.from_dict(cfg.to_dict()).baseline_cells == "untreated"
    # Archives written before the option existed load with the old behaviour.
    old = {k: v for k, v in cfg.to_dict().items() if k != "baseline_cells"}
    assert LongBetConfig.from_dict(old).baseline_cells == "all"


def test_multi_outcome_rejects_the_cut():
    from longbet import LongBetMulti

    x, y, z = _null_panel()
    model = LongBetMulti(_small_config("untreated"))
    with pytest.raises(NotImplementedError, match="LongBetMulti"):
        model.fit(y={"a": y, "b": y + 1.0}, x=x, z=z)


# --------------------------------------------------------------------------
# the defining property: treated outcomes never reach the baseline module
# --------------------------------------------------------------------------


def _baseline_draws(model, x, z):
    pred = model.predict(x=x, z=z)
    tr = model.trace
    return {
        "mu0": np.asarray(pred.muhats0),
        "gamma": np.asarray(tr.gamma),
        "trend_coef": np.asarray(tr.trend_coef),
        "sigma2": np.asarray(tr.sigma2),
        "sigma_gamma2": np.asarray(tr.sigma_gamma2),
    }, np.asarray(pred.att_full)


@pytest.mark.parametrize("cells", ["untreated", "all"])
def test_treated_cells_are_invisible_to_the_baseline_iff_cut(cells):
    x, y, z = _null_panel()
    shift = 5.0
    y_shifted = (y + shift * z).astype(np.float32)

    m0 = LongBet(_small_config(cells)).fit(x=x, y=y, z=z)
    m1 = LongBet(_small_config(cells)).fit(x=x, y=y_shifted, z=z)
    d0, att0 = _baseline_draws(m0, x, z)
    d1, att1 = _baseline_draws(m1, x, z)

    if cells == "untreated":
        # Exact invariance: the baseline module's random streams and
        # sufficient statistics do not involve a single treated outcome.
        for name in d0:
            np.testing.assert_array_equal(d0[name], d1[name], err_msg=name)
        assert m1.state.baseline_untreated_only
        assert not m1.state.use_inter_ensemble_move
        base = np.asarray(m1.state.base_mask)
        assert base.sum() == int(np.sum(z == 0))
        # The treatment module did see the shift: the ATT moves by about it.
        gap = att1.mean() - att0.mean()
        assert 0.6 * shift < gap < 1.4 * shift, gap
    else:
        # Under the joint posterior the shift reaches the baseline.
        assert not np.array_equal(d0["mu0"], d1["mu0"])
        assert not np.array_equal(d0["gamma"], d1["gamma"])


def test_cut_state_keeps_the_residual_invariant():
    """``resid == y - fitted`` on every observed cell after the cut sweep."""
    x, y, z = _null_panel()
    y = y.copy()
    y[2, 1] = np.nan  # one missing cell, so obs_mask != all-true
    m = LongBet(_small_config("untreated")).fit(x=x, y=y, z=z)
    st = m.state
    y_vec = np.asarray(st.y)
    obs = np.asarray(st.obs_mask)
    b_z = np.where(np.asarray(st.z_vec) == 1.0, np.asarray(st.b1)[:, None], np.asarray(st.b0)[:, None])
    beta = np.asarray(st.beta)[:, np.asarray(st.exposure_idx)]
    fitted = (
        np.asarray(st.alpha)[:, None] * np.asarray(st.mu_fit)
        + np.asarray(st.trend_coef) @ np.asarray(st.trend_design).T
        + b_z * beta * np.asarray(st.nu_fit)
        + np.asarray(st.gamma)[:, np.asarray(st.unit_idx)]
    )
    resid = np.asarray(st.resid)
    np.testing.assert_allclose(resid[:, obs], (y_vec[None, :] - fitted)[:, obs], atol=2e-4, rtol=0)
    assert np.all(resid[:, ~obs] == 0.0)
    assert np.asarray(st.base_mask)[~obs].sum() == 0


# --------------------------------------------------------------------------
# module 1 is the missing-data posterior
# --------------------------------------------------------------------------


@pytest.mark.slow
def test_cut_baseline_matches_fit_with_treated_cells_missing():
    """``untreated`` on the full panel == ``all`` on the panel with treated
    cells set to missing, for every baseline quantity, up to Monte Carlo error.

    The two samplers consume their random streams differently (the cut draws
    beta separately from the trend block), so the comparison is on posterior
    means with a tolerance of a few Monte Carlo standard errors.
    """
    x, y, z = _null_panel(seed=5, N=40, T=8, noise=0.25)
    y_masked = np.where(z == 1, np.nan, y).astype(np.float32)
    kw = dict(num_burnin=150, num_sweeps=150, num_chains=2, num_trees_pr=10,
              num_trees_trt=10, trend_num_features=8, random_seed=4, standardize=False)
    m_cut = LongBet(LongBetConfig(baseline_cells="untreated", **kw)).fit(x=x, y=y, z=z)
    m_mis = LongBet(LongBetConfig(baseline_cells="all", **kw)).fit(x=x, y=y_masked, z=z)
    p_cut = m_cut.predict(x=x, z=z)
    p_mis = m_mis.predict(x=x, z=z)

    mu_cut = np.asarray(p_cut.muhats0)
    mu_mis = np.asarray(p_mis.muhats0)
    # Posterior means of the untreated surface on *treated* cells (the
    # imputations) and on untreated cells.
    for mask, label in ((z == 1, "treated cells"), (z == 0, "untreated cells")):
        a, b = mu_cut[mask].mean(axis=-1), mu_mis[mask].mean(axis=-1)
        sd = 0.5 * (mu_cut[mask].std(axis=-1) + mu_mis[mask].std(axis=-1))
        # draws are autocorrelated; allow ~1 posterior SD on the RMS gap
        rms_gap = np.sqrt(np.mean((a - b) ** 2))
        assert rms_gap < 0.6 * np.mean(sd) + 0.05, (label, rms_gap, np.mean(sd))
    s_cut, s_mis = np.asarray(m_cut.trace.sigma2), np.asarray(m_mis.trace.sigma2)
    assert abs(s_cut.mean() - s_mis.mean()) < 0.25 * max(s_cut.std(), s_mis.std()) + 0.01


# --------------------------------------------------------------------------
# trend_block_step: resid_mask
# --------------------------------------------------------------------------


def test_trend_block_step_resid_mask_default_is_identity_and_superset_updates_extra_cells():
    rng = np.random.default_rng(0)
    N, T, q, S1 = 12, 6, 5, 4
    n = N * T
    unit_idx = jnp.asarray(np.repeat(np.arange(N), T), jnp.int32)
    exposure_idx = jnp.asarray(rng.integers(0, S1, size=n), jnp.int32)
    H = jnp.asarray(rng.normal(size=(n, q)), jnp.float32)
    R = jnp.asarray(rng.normal(size=n), jnp.float32)
    d_vec = jnp.asarray(rng.normal(size=n), jnp.float32)
    obs = jnp.asarray(rng.uniform(size=n) < 0.6)
    all_cells = jnp.ones(n, dtype=bool)
    L = jnp.asarray(np.linalg.cholesky(np.eye(S1) + 0.3), jnp.float32)
    args = dict(
        trend_coef=jnp.asarray(rng.normal(size=q), jnp.float32),
        gamma=jnp.asarray(rng.normal(size=N), jnp.float32),
        beta=jnp.asarray(rng.normal(size=S1), jnp.float32),
        R=R, H=H, d_vec=d_vec, unit_idx=unit_idx, N_units=N, exposure_idx=exposure_idx,
        sigma2=jnp.float32(0.5), sigma_gamma2=jnp.float32(0.4), K_chol=L,
        trend_col_scale=jnp.ones(q, jnp.float32), random_intercept=True, sample_beta=False,
    )
    key = jax.random.key(1)
    ref = trend_block_step(key, obs_mask=obs, **args)
    same = trend_block_step(key, obs_mask=obs, resid_mask=obs, **args)
    for a, b in zip(ref, same):
        np.testing.assert_array_equal(np.asarray(a), np.asarray(b))

    sup = trend_block_step(key, obs_mask=obs, resid_mask=all_cells, **args)
    # Same (c, gamma) draw: the extra cells did not inform it ...
    np.testing.assert_array_equal(np.asarray(ref[0]), np.asarray(sup[0]))
    np.testing.assert_array_equal(np.asarray(ref[1]), np.asarray(sup[1]))
    # ... and the residual agrees on the informing cells while the extra
    # cells carry the partial residual minus the new fit.
    obs_np = np.asarray(obs)
    np.testing.assert_array_equal(np.asarray(ref[3])[obs_np], np.asarray(sup[3])[obs_np])
    c_new, g_new = np.asarray(sup[0], np.float64), np.asarray(sup[1], np.float64)
    R_tilde = (np.asarray(R, np.float64) + np.asarray(H, np.float64) @ np.asarray(args["trend_coef"], np.float64)
               + np.asarray(args["gamma"], np.float64)[np.asarray(unit_idx)])
    expected = R_tilde - np.asarray(H, np.float64) @ c_new - g_new[np.asarray(unit_idx)]
    np.testing.assert_allclose(np.asarray(sup[3])[~obs_np], expected[~obs_np], atol=1e-5)


# --------------------------------------------------------------------------
# effect recovery under the cut (coarse)
# --------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("tau", [0.0, 3.0])
def test_cut_recovers_a_constant_effect(tau):
    x, y, z = _null_panel(seed=9, N=60, T=8, noise=0.3)
    y = (y + tau * z).astype(np.float32)
    cfg = LongBetConfig(baseline_cells="untreated", num_burnin=200, num_sweeps=150,
                        num_chains=2, num_trees_pr=10, num_trees_trt=20,
                        trend_num_features=8, random_seed=2)
    pred = LongBet(cfg).fit(x=x, y=y, z=z).predict(x=x, z=z)
    att = np.asarray(pred.att_full)  # (S_max, draws)
    draws = att.mean(axis=0)
    assert abs(draws.mean() - tau) < 3.0 * draws.std() + 0.15, (draws.mean(), draws.std())
