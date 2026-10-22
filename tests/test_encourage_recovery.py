"""The adoption-clock model recovers offer effects on a panel with exact expected truth."""

import numpy as np
import pytest
from scipy.stats import norm

from longbet import LongBetEncourage


def _arm_stats(h, g, tt):
    """P(adopted by t) and E[g(S_t)] for one unit under one hazard path."""
    surv = np.cumprod(1 - h)
    p_adopt = h * np.concatenate([[1.0], surv[:-1]])
    cdf = np.cumsum(p_adopt)
    eg = np.array([sum(p_adopt[a] * g(t - a + 1) for a in range(t + 1)) for t in range(tt)])
    return cdf, eg


def onboarding_panel(n, seed, tt=6, t_enc=3, sd_eps=1.5):
    """Absorbing adoption with organic and offer-driven hazards; effect tau(x) * log1p(exposure)."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 3))
    u = rng.binomial(1, 0.5, n)
    z_unit = rng.binomial(1, 0.5, n)
    z = np.zeros((n, tt)); z[z_unit == 1, t_enc - 1:] = 1
    p0 = np.where(u == 1, 0.20, 0.05)
    p_comp = 0.20 + 0.65 / (1 + np.exp(-3 * x[:, 0]))
    qualifying = (x[:, 0] > 0) & (x[:, 1] > -0.2)
    tau = 1.2 + 7.8 * qualifying * (0.5 + 0.5 * norm.cdf(x[:, 1]))
    g = np.log1p
    v = rng.uniform(size=(n, tt))
    h0 = np.tile(p0[:, None], (1, tt))
    h1 = h0.copy(); h1[:, t_enc - 1:] = p_comp[:, None]
    def first(h):
        hit = v < h
        out = np.full(n, np.inf)
        any_hit = hit.any(axis=1)
        out[any_hit] = np.argmax(hit[any_hit], axis=1) + 1
        return out
    a0, a1 = first(h0), first(h1)
    def exposure(a):
        s = np.maximum(0, np.arange(1, tt + 1)[None, :] - a[:, None] + 1)
        return np.where(np.isfinite(a)[:, None], s, 0)
    s0, s1 = exposure(a0), exposure(a1)
    d = np.where(z == 1, s1 > 0, s0 > 0).astype(float)
    alpha_i = rng.normal(10, 1.5, n)
    mu = alpha_i[:, None] + 0.5 * np.arange(1, tt + 1)[None, :] - 2.5 * u[:, None] + 0.8 * x[:, :1] - 0.5 * x[:, 1:2]
    eps = rng.normal(0, sd_eps, (n, tt))
    y = mu + tau[:, None] * g(np.where(z == 1, s1, s0)) + eps
    # Exact expected offer effects, averaging over the latent type u.
    citt_y = np.zeros((n, tt)); citt_d = np.zeros((n, tt))
    for i in range(n):
        for level in (0.05, 0.20):
            hz0 = np.full(tt, level); hz1 = hz0.copy(); hz1[t_enc - 1:] = p_comp[i]
            c0, e0 = _arm_stats(hz0, g, tt); c1, e1 = _arm_stats(hz1, g, tt)
            citt_d[i] += (c1 - c0) / 2; citt_y[i] += tau[i] * (e1 - e0) / 2
    return dict(y=y, d=d, z=z, x=x, t=np.arange(1, tt + 1, dtype=float), citt_y=citt_y, citt_d=citt_d,
                start=t_enc - 1)


@pytest.mark.slow
def test_offer_effects_are_recovered_on_a_simulated_pilot():
    train = onboarding_panel(240, seed=2026)
    test = onboarding_panel(200, seed=9999)
    fit = LongBetEncourage(num_chains=2, num_burnin=600, num_sweeps=400, num_trees_pr=10,
                           num_trees_trt=10, random_seed=7).fit(train["y"], train["d"], train["z"],
                                                                 train["x"], t=train["t"])
    post = slice(train["start"], None)
    cond = fit.predict_conditional(x=test["x"])
    rmse_y = float(np.sqrt(np.mean((cond.citt_y.mean - test["citt_y"][:, post]) ** 2)))
    rmse_d = float(np.sqrt(np.mean((cond.citt_d.mean - test["citt_d"][:, post]) ** 2)))
    # A pooled offer-effect estimate ignores covariates entirely; the model must beat it clearly.
    pooled = np.tile(train["citt_y"][:, post].mean(axis=0), (200, 1))
    rmse_pooled = float(np.sqrt(np.mean((pooled - test["citt_y"][:, post]) ** 2)))
    assert rmse_y < 0.6 * rmse_pooled, (rmse_y, rmse_pooled)
    assert rmse_d < 0.15, rmse_d
    pred = fit.predict()
    draws = pred.draws["itt_y"][0, -1]                      # week T, all units
    truth = train["citt_y"][:, -1].mean()
    assert abs(draws.mean() - truth) < 3 * draws.std() + 0.3, (draws.mean(), truth, draws.std())
    # The model's population effect and the design-based reference must agree at week T.
    ref = pred.reference.query("group == 'all'").iloc[-1]
    assert abs(draws.mean() - ref.itt_y) < 3 * ref.itt_y_se + 0.3, (draws.mean(), ref.itt_y, ref.itt_y_se)
