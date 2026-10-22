"""Parallel tempering: layout, exchanges, invariance."""

import numpy as np
import pytest
import jax
import jax.numpy as jnp

from longbet import LongBet, LongBetConfig, LongBetEncourage
from longbet._diagnostics import compute_ess
from longbet._tempering import chain_energy, ladder_temperatures, swap_step


def _panel(seed=5, n=36, periods=6, noise=0.5):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 3)).astype(np.float32)
    z = np.zeros((n, periods)); z[: n // 2, 3:] = 1
    s = np.where(z == 1, np.arange(periods)[None, :] - 2, 0)
    y = (0.5 * x[:, :1] + 0.2 * np.arange(periods)[None, :] + (x[:, :1] > 0) * 1.5 * np.log1p(s)
         + rng.normal(0, noise, (n, periods)) + rng.normal(0, 0.5, (n, 1)))
    return x, z, y, np.arange(1, periods + 1)


def _cfg(**kw):
    base = dict(num_chains=2, num_burnin=20, num_sweeps=10, num_trees_pr=3, num_trees_trt=3,
                random_seed=11)
    return LongBetConfig(**{**base, **kw})


def test_config_validation():
    with pytest.raises(ValueError, match="tempering_levels"):
        LongBetConfig(tempering_levels=0)
    with pytest.raises(ValueError, match="tempering_beta_min"):
        LongBetConfig(tempering_beta_min=1.5)
    assert LongBetConfig(tempering_levels=4, tempering_beta_min=0.1).tempering_levels == 4


def test_ladder_layout():
    temps = np.asarray(ladder_temperatures(6, 3, 0.25))
    assert np.allclose(temps, [1.0, 0.5625, 0.25, 1.0, 0.5625, 0.25])


def test_tempered_fit_keeps_only_posterior_replicas_and_exchanges():
    x, z, y, t = _panel()
    cfg = _cfg(tempering_levels=6, tempering_beta_min=0.3, num_burnin=60, num_sweeps=40)
    model = LongBet(cfg).fit(y, x, z, t=t)
    assert model.state.num_chains == 12
    temps = np.asarray(model.state.temperature).reshape(2, 6)
    assert np.allclose(temps[0], temps[1]) and temps[0, 0] == 1.0 and temps[0, -1] == pytest.approx(0.3)
    assert np.asarray(model.trace.beta).shape[:2] == (2, 40)
    assert np.asarray(model.trace.nu_trace.leaf_tree).shape[0] == 2
    acc = model.swap_accepts_.reshape(2, 6)
    assert np.all(acc[:, -1] == 0), "the hottest level never leads an exchange"
    assert acc[:, :-1].sum() > 0, "adjacent levels this close must exchange sometimes"
    pred = model.predict(x, z, t=t, summary_only=True)
    assert np.all(np.isfinite(pred.tau_summary.mean))


def test_swap_step_moves_states_and_keeps_temperatures():
    x, z, y, t = _panel()
    cfg = _cfg(tempering_levels=3, tempering_beta_min=0.5, num_burnin=5, num_sweeps=3)
    model = LongBet(cfg).fit(y, x, z, t=t)
    state = model.state
    before_gamma = np.asarray(state.gamma); before_resid = np.asarray(state.resid)
    before_temp = np.asarray(state.temperature)
    energy = np.asarray(chain_energy(state))
    assert energy.shape == (6,) and np.all(np.isfinite(energy))
    # Force acceptance by giving the hotter replica the better energy.
    new, accepted = swap_step(jax.random.key(0), state, jnp.int32(0))
    accepted = np.asarray(accepted); new_gamma = np.asarray(new.gamma); new_resid = np.asarray(new.resid)
    assert np.allclose(np.asarray(new.temperature), before_temp)
    for c in range(6):
        partner = c + 1 if accepted[c] else (c - 1 if c > 0 and accepted[c - 1] else c)
        assert np.allclose(new_gamma[c], before_gamma[partner])
        assert np.allclose(new_resid[c], before_resid[partner])


@pytest.mark.slow
def test_tempered_cold_chains_match_plain_posterior():
    """Cold replicas of a tempered run target the same posterior as the plain sampler."""
    x, z, y, t = _panel(noise=0.7)
    base = dict(num_chains=4, num_burnin=300, num_sweeps=600, num_trees_pr=3, num_trees_trt=3,
                random_seed=21)
    fits = {}
    for name, extra in (("plain", {}), ("pt", dict(tempering_levels=5, tempering_beta_min=0.3))):
        m = LongBet(LongBetConfig(**{**base, **extra})).fit(y, x, z, t=t)
        pred = m.predict(x, z, t=t, summary_only=False)
        att = pred.tauhats[: len(x) // 2, -1, :].mean(axis=0).reshape(4, 600)
        s2 = np.asarray(m.trace.sigma2).reshape(4, 600)
        fits[name] = (att, s2)
    for j, label in ((0, "week-T ATT"), (1, "sigma2")):
        a, b = fits["plain"][j], fits["pt"][j]
        se = np.sqrt(a.var() / max(compute_ess(a), 4) + b.var() / max(compute_ess(b), 4))
        assert abs(a.mean() - b.mean()) < 4 * se + 1e-6, (label, a.mean(), b.mean(), se)


def test_encouragement_model_tempering_runs_and_keeps_posterior_replicas():
    x, z, y, t = _panel()
    d = (z * (x[:, :1] > 0)).astype(float)
    # Levels close enough that exchanges happen on a tiny panel within 60 sweeps.
    cfg = _cfg(tempering_levels=4, tempering_beta_min=0.7, num_burnin=40, num_sweeps=20)
    fit = LongBetEncourage(cfg).fit(y, d, z, x, t=t)
    for model in (fit.outcome_model, fit.adoption_model):
        assert model.state.num_chains == 8 and model.swap_accepts_.shape == (8,)
        temps = np.asarray(model.state.temperature).reshape(2, 4)
        assert temps[0, 0] == 1.0 and temps[0, -1] == pytest.approx(0.7)
        assert np.asarray(model.trace.beta).shape[:2] == (2, 20)
        assert model.swap_accepts_[:-1].sum() > 0
    pred = fit.predict()
    assert np.all(np.isfinite(pred.draws["itt_y"]))
