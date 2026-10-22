"""Potential-outcome binary encouragement experiments with explicit truths."""
from __future__ import annotations

import numpy as np

SCENARIOS = ("strong", "weak", "zero", "one_sided", "direct_effect", "defiers")
PATTERNS = np.array([[1, 1], [0, 0], [0, 1], [1, 0]])  # always, never, complier, defier


def make_binary_data(scenario, *, n=240, population_seed=1, assignment_seed=2):
    if scenario not in SCENARIOS:
        raise ValueError("Unknown binary scenario.")
    if isinstance(n, bool) or not isinstance(n, int) or n < 12 or n % 6:
        raise ValueError("n must be a multiple of six, at least twelve.")
    x = np.tile(np.arange(3), n // 3)[:, None]
    index = x[:, 0]
    pi = np.zeros((3, 4))
    pi[:, 0] = [.15, .2, .25]
    pi[:, 2] = [.35, .45, .55]
    if scenario == "weak":
        pi[:, 2] = [.02, .04, .06]
    elif scenario == "zero":
        pi[:, 2] = 0
    elif scenario == "one_sided":
        pi[:, 0] = 0
    elif scenario == "defiers":
        pi[:, 2], pi[:, 3] = .08, .35
    pi[:, 1] = 1 - pi.sum(axis=1)
    gamma = np.empty((3, 4, 2, 2))
    for g in range(3):
        for s in range(4):
            for dd in (0, 1):
                for zz in (0, 1):
                    gamma[g, s, dd, zz] = (.15 + .08 * g + [.20, 0, .08, .16][s]
                                           + (.1 + .06 * g) * dd
                                           + (zz * [.08, .04, .12, .1][s] if scenario == "direct_effect" else 0))
    assert np.all((gamma >= 0) & (gamma <= 1)) and np.all(pi >= 0)
    p = pi @ PATTERNS
    q = np.empty((3, 2, 2))
    for dd in (0, 1):
        for zz in (0, 1):
            share = np.sum(pi * (PATTERNS[:, zz] == dd), axis=-1)
            mass = np.sum(pi * (PATTERNS[:, zz] == dd) * gamma[:, :, dd, zz], axis=-1)
            # An empty observed D,Z cell has no defined conditional outcome
            # probability. This arbitrary completion receives exactly zero mass.
            q[:, dd, zz] = np.divide(mass, share, out=np.full(3, .5), where=share > 0)
    r = p * q[:, 1] + (1 - p) * q[:, 0]
    weights = np.full(3, 1 / 3)
    conditional_itt = np.array([weights @ (r[:, 1] - r[:, 0]), weights @ (p[:, 1] - p[:, 0])])
    if scenario == "zero":
        conditional_itt[:] = 0.  # algebraically exact null, not numerical division noise
    conditional_complier = np.full(3, np.nan)
    mass = weights @ pi[:, 2]
    if mass > 0:
        conditional_complier = np.array([weights @ (pi[:, 2] * contrast) / mass for contrast in
            (gamma[:, 2, 1, 1] - gamma[:, 2, 0, 0], gamma[:, 2, 1, 0] - gamma[:, 2, 0, 0],
             gamma[:, 2, 1, 1] - gamma[:, 2, 0, 1])])
    rng = np.random.default_rng(population_seed)
    strata = np.array([rng.choice(4, p=pi[g]) for g in index])
    potential_d = PATTERNS[strata]
    # Every potential response exists before assignment. A common latent uniform
    # enforces exclusion exactly when gamma_s,d,0=gamma_s,d,1.
    potential_y = (rng.random(n)[:, None, None] < gamma[index, strata]).astype(int)
    assignment_rng = np.random.default_rng(assignment_seed)
    z = assignment_rng.permutation(np.arange(n) % 2)
    d = potential_d[np.arange(n), z]
    y = potential_y[np.arange(n), d, z]
    potential_observed_y = np.stack([potential_y[np.arange(n), potential_d[:, zz], zz] for zz in (0, 1)], axis=-1)
    finite_itt = np.array([(potential_observed_y[:, 1] - potential_observed_y[:, 0]).mean(),
                           (potential_d[:, 1] - potential_d[:, 0]).mean()])
    complier = strata == 2
    finite_complier = np.full(3, np.nan)
    if complier.any():
        cy = potential_y[complier]
        finite_complier = np.array([(cy[:, 1, 1] - cy[:, 0, 0]).mean(),
                                   (cy[:, 1, 0] - cy[:, 0, 0]).mean(),
                                   (cy[:, 1, 1] - cy[:, 0, 1]).mean()])
    v = np.stack([potential_observed_y, potential_d], axis=-1)  # N,Z,(Y,D)
    neyman = (np.cov(v[:, 0].T, ddof=1) + np.cov(v[:, 1].T, ddof=1)) / (n // 2)
    gap = np.cov((v[:, 1] - v[:, 0]).T, ddof=1) / n
    return dict(y=y, d=d, z=z, x=x, strata=strata, potential_y=potential_y, potential_d=potential_d,
                p_truth=p, q_truth=q, r_truth=r, weights=weights, pi_truth=pi, gamma_truth=gamma,
                conditional_itt=conditional_itt, finite_itt=finite_itt,
                conditional_complier=conditional_complier, finite_complier=finite_complier,
                assignment_covariance=neyman-gap, expected_neyman_covariance=neyman,
                population_seed=np.asarray(population_seed), assignment_seed=np.asarray(assignment_seed),
                scenario=np.asarray(scenario), n=np.asarray(n))


def targets(data):
    result = {}
    for label in ("conditional", "finite"):
        y, d = data[label + "_itt"]
        result[label] = dict(itt_y=float(y), itt_d=float(d), wald=float(y / d) if d != 0 else np.nan)
        for key, value in zip(("complier_encouragement", "treatment_z0", "treatment_z1"), data[label + "_complier"]):
            result[label][key] = float(value)
    return result
