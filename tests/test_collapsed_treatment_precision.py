# Copyright 2026 Google LLC
# SPDX-License-Identifier: Apache-2.0

"""Gaussian covariance oracles for the collapsed treatment precision repair."""

import jax
import jax.numpy as jnp
import equinox as eqx
import numpy as np
import pytest

from longbet._step import WEIGHT_EPS, _stabilize_collapsed_nu_stats
from longbet._scales import compute_treatment_scale_attrs


def _forest_statistics(w, u, mask, units, sigma2, sigma_gamma2, cell_precision):
    """Exercise the production repair and the forest's physical scale encoding."""
    num_units = int(np.max(units)) + 1
    w, u = jnp.asarray(w, jnp.float32), jnp.asarray(u, jnp.float32)
    mask, units = jnp.asarray(mask), jnp.asarray(units)
    w_obs = jnp.where(mask, w, 0.0)
    u_obs = jnp.where(mask, u, 0.0)
    p = jnp.where(mask, 1.0 if cell_precision is None else jnp.asarray(cell_precision), 0.0)
    p = p.astype(jnp.float32)

    def total(a):
        return jnp.zeros(num_units, jnp.float32).at[units].add(a)

    c = sigma2 / sigma_gamma2 if cell_precision is None else 1.0 / sigma_gamma2
    second = total(p * w_obs**2)
    d = total(p) + c
    kappa = jnp.minimum(1.0 - total(p * w_obs)**2 / (jnp.maximum(second, 1e-12) * jnp.maximum(d, 1e-6)), 1.0)
    score = total(p * w_obs * (u_obs - (total(p * u_obs) / jnp.maximum(d, 1e-6))[units]))
    target = score / jnp.maximum(kappa * second, 1e-12)

    @jax.jit
    def repaired(kappa, target):
        return _stabilize_collapsed_nu_stats(
            kappa, target, second, d, u, w, mask, units, num_units, c,
            None if cell_precision is None else p,
        )

    stable_kappa, stable_target, repair = repaired(kappa, target)
    precision_scale = 1.0 if cell_precision is None else p * sigma2
    effective_w = w * jnp.sqrt(stable_kappa[units] * precision_scale)
    attrs, active = compute_treatment_scale_attrs(
        effective_w, mask, eps=jnp.where(repair[units], 0.0, WEIGHT_EPS)
    )
    physical_precision = attrs.prec_scale * attrs.inv_sdev_unit**2 / sigma2
    Q = np.asarray(total(physical_precision))
    h = np.asarray(total(physical_precision * stable_target[units]))
    return Q, h, np.asarray(repair), np.asarray(active), (kappa, target), (stable_kappa, stable_target)


@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("case", ["constant", "missing", "tiny", "zero"])
def test_collapsed_treatment_matches_dense_gaussian_oracle(weighted, case):
    n = 100
    units = np.repeat(np.arange(3), n).astype(np.int32)
    mask = np.ones(3 * n, dtype=bool)
    w = np.ones(3 * n, dtype=np.float32)
    u = np.linspace(-2.0, 3.0, 3 * n).astype(np.float32)
    if case == "missing":
        mask[::7] = False
        mask[2 * n:] = False
        # A missing response must not contaminate sufficient statistics.
        u[~mask] = np.nan
    elif case == "tiny":
        w *= 1e-8
        w[n:2*n] *= np.linspace(-1.0, 1.0, n)
    elif case == "zero":
        w[n:2*n] = 0.0
        mask[2*n:] = False
    sigma2, sigma_gamma2 = 1.0, 1e4
    p = np.linspace(0.25, 2.0, 3*n).astype(np.float32) if weighted else None
    Q, h, repair, active, _, _ = _forest_statistics(w, u, mask, units, sigma2, sigma_gamma2, p)
    np.testing.assert_array_equal(
        repair, [np.any(mask[units == unit] & (w[units == unit] != 0.0)) for unit in range(3)]
    )
    for unit in range(3):
        selected = (units == unit) & mask
        wi, ui = w[selected].astype(float), u[selected].astype(float)
        pi = np.ones(wi.size) / sigma2 if p is None else p[selected].astype(float)
        covariance = np.diag(1.0 / pi) + sigma_gamma2 * np.ones((wi.size, wi.size))
        oracle_Q = wi @ np.linalg.solve(covariance, wi)
        oracle_h = wi @ np.linalg.solve(covariance, ui)
        np.testing.assert_allclose(Q[unit], oracle_Q, rtol=2e-5, atol=1e-28)
        np.testing.assert_allclose(h[unit], oracle_h, rtol=2e-5, atol=1e-20)
        # Independent N(0,1) leaf prior gives this exact posterior.
        np.testing.assert_allclose(1.0/(1.0+Q[unit]), 1.0/(1.0+oracle_Q), rtol=1e-7)
        np.testing.assert_allclose(h[unit]/(1.0+Q[unit]), oracle_h/(1.0+oracle_Q), rtol=2e-5, atol=1e-20)
    np.testing.assert_array_equal(active, mask & (w != 0.0))
    if case == "constant" and not weighted:
        # Before this repair the floor yielded Q=0.01, about 100 times too large.
        np.testing.assert_allclose(Q, n/(1.0+n*sigma_gamma2), rtol=1e-6)


@pytest.mark.parametrize("weighted", [False, True])
def test_ordinary_collapsed_statistics_are_preserved(weighted):
    units = np.repeat(np.arange(2), 12).astype(np.int32)
    w = np.tile([0.0]*6 + [0.2, 0.4, 0.5, 0.7, 0.8, 1.0], 2)
    u = np.linspace(-2.0, 1.0, 24)
    p = np.linspace(0.5, 1.5, 24) if weighted else None
    _, _, repair, _, old, new = _forest_statistics(w, u, np.ones(24, bool), units, 1.0, 2.0, p)
    assert not np.any(repair)
    for before, after in zip(old, new):
        np.testing.assert_array_equal(before, after)


@pytest.mark.parametrize("weighted", [False, True])
def test_extreme_collapsed_sampler_step_compiles_and_preserves_residual(weighted):
    from longbet._config import LongBetConfig
    from longbet._state import init_longbet
    from longbet._step import longbet_single_step

    n, periods = 200, 100
    units = jnp.repeat(jnp.arange(2), periods)
    y = jnp.linspace(-1.0, 2.0, n, dtype=jnp.float32)
    state = init_longbet(
        X_unified=jnp.zeros((1, n), jnp.uint8), y=y, unit_idx=units,
        time_idx=jnp.tile(jnp.arange(periods), 2),
        exposure_idx=jnp.ones(n, jnp.int32), z_vec=jnp.ones(n, jnp.float32),
        obs_mask=jnp.ones(n, bool), max_split_mu=jnp.zeros(1, jnp.uint8),
        max_split_nu=jnp.zeros(1, jnp.uint8),
        config=LongBetConfig(num_trees_pr=1, num_trees_trt=1, sample_beta=False),
    )
    state = eqx.tree_at(lambda s: (s.sigma2, s.sigma_gamma2), state,
                       (jnp.float32(1.0), jnp.float32(1e4)))
    assert state.collapse_gamma_mu and state.collapse_gamma_nu
    precision = jnp.linspace(0.25, 2.0, n, dtype=jnp.float32) if weighted else None
    result = longbet_single_step(jax.random.key(281), state, precision)
    weights = jnp.where(result.z_vec == 1, result.b1, result.b0) * result.beta[result.exposure_idx]
    expected = (y - result.alpha*result.mu_fit - weights*result.nu_fit
                - result.gamma[units] - result.trend_design @ result.trend_coef)
    assert np.all(np.isfinite(np.asarray(result.resid)))
    np.testing.assert_allclose(result.resid, expected, rtol=2e-5, atol=2e-5)
