"""Unit tests for ``longbet._rank_k``."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from longbet import (
    RankKConfigWarning,
    RankKLongBetConfig,
    fit_rank_k_longbet,
    orthonormal_exposure_basis,
    reconcile_nested_3arm_draws,
)


def test_orthonormal_exposure_basis_gram_and_zero_sum() -> None:
    for L in (4, 6, 7, 8, 12):
        psi = orthonormal_exposure_basis(L, 4)
        assert psi.shape == (L, 4)
        gram = (psi.T @ psi) / float(L)
        np.testing.assert_allclose(gram, np.eye(4), atol=1e-12)
        np.testing.assert_allclose(psi[:, 0], np.ones(L), atol=1e-12)
        for k in range(1, 4):
            np.testing.assert_allclose(np.sum(psi[:, k]), 0.0, atol=1e-12)


@pytest.mark.parametrize("uncertainty_propagation", [False, True], ids=["plug-in", "propagated"])
def test_exact_mode1_decoupling_theorem(uncertainty_propagation: bool) -> None:
    if not uncertainty_propagation:
        pytest.importorskip("stochtree")
    rng = np.random.default_rng(2026)
    N, T = 80, 10
    x = rng.normal(size=(N, 3))
    cohorts = np.full(N, T, dtype=np.int32)
    cohorts[:20] = 4
    cohorts[20:40] = 6
    z = (np.arange(T)[None, :] >= cohorts[:, None]).astype(np.int32)
    y = rng.normal(size=(N, T)) + z * (1.2 + 0.5 * x[:, :1])
    f = 0.8 * y + rng.normal(scale=0.3, size=(N, T))
    u = rng.normal(size=(N, 2))
    ps = np.full(N, 0.5)

    # Short runs; the propagated variant keeps its default draws (num_mcmc is unused there).
    budget = {"num_burnin": 50} if uncertainty_propagation else {"num_burnin": 50, "num_mcmc": 100}
    results = {}
    with warnings.catch_warnings():
        warnings.simplefilter("error", RankKConfigWarning)
        for K in (1, 2, 3):
            cfg = RankKLongBetConfig(
                rank_k=K, uncertainty_propagation=uncertainty_propagation, random_seed=77, **budget
            )
            results[K] = fit_rank_k_longbet(y, z, x, unit_features=u, ps=ps, baseline_forecast=f, cfg=cfg)
    res1, res2, res3 = results[1], results[2], results[3]

    np.testing.assert_allclose(res1.att_draws, res2.att_draws, atol=1e-12)
    np.testing.assert_allclose(res1.att_draws, res3.att_draws, atol=1e-12)
    for i in range(40):
        g = cohorts[i]
        u1 = res1.catt_hat[i, g:].mean()
        u2 = res2.catt_hat[i, g:].mean()
        u3 = res3.catt_hat[i, g:].mean()
        np.testing.assert_allclose(u1, u2, atol=1e-12)
        np.testing.assert_allclose(u1, u3, atol=1e-12)


def test_reconcile_nested_3arm_draws_additivity_and_dominance() -> None:
    pytest.importorskip("stochtree")
    rng = np.random.default_rng(2026)
    n_c, n_b, n_a, t_periods = 30, 30, 30, 6
    x_c = rng.normal(size=(n_c, 3))
    x_b = rng.normal(size=(n_b, 3))
    x_a = rng.normal(size=(n_a, 3))

    base_c = rng.normal(size=(n_c, t_periods))
    base_b = rng.normal(size=(n_b, t_periods))
    base_a = rng.normal(size=(n_a, t_periods))

    z_bc = np.zeros((n_b + n_c, t_periods), dtype=np.int32)
    z_bc[:n_b, 1:] = 1
    y_bc = np.vstack([base_b + 0.8 * z_bc[:n_b], base_c])
    x_bc = np.vstack([x_b, x_c])

    z_ab = np.zeros((n_a + n_b, t_periods), dtype=np.int32)
    z_ab[:n_a, 1:] = 1
    y_ab = np.vstack([base_a + 1.9 * z_ab[:n_a], base_b + 0.8 * z_bc[:n_b]])
    x_ab = np.vstack([x_a, x_b])

    z_ac = np.zeros((n_a + n_c, t_periods), dtype=np.int32)
    z_ac[:n_a, 1:] = 1
    y_ac = np.vstack([base_a + 1.9 * z_ac[:n_a], base_c])
    x_ac = np.vstack([x_a, x_c])

    cfg = RankKLongBetConfig(
        rank_k=2,
        uncertainty_propagation=False,
        num_burnin=40,
        num_mcmc=80,
        random_seed=101,
    )
    res_bc = fit_rank_k_longbet(y_bc, z_bc, x_bc, cfg=cfg)
    res_ab = fit_rank_k_longbet(y_ab, z_ab, x_ab, cfg=cfg)
    res_ac = fit_rank_k_longbet(y_ac, z_ac, x_ac, cfg=cfg)

    trt_bc, trt_ab, trt_ac = reconcile_nested_3arm_draws(res_bc, res_ab, res_ac, prior_scale=2.0, seed=2026)
    assert trt_bc.shape == (1600, t_periods - 1)
    assert trt_ab.shape == (1600, t_periods - 1)
    assert trt_ac.shape == (1600, t_periods - 1)

    cum_bc = trt_bc.sum(axis=1)
    cum_ab = trt_ab.sum(axis=1)
    cum_ac = trt_ac.sum(axis=1)

    # Exact structural additivity of posterior means and draw-by-draw sums
    np.testing.assert_allclose(cum_ac, cum_bc + cum_ab, atol=1e-12)
    np.testing.assert_allclose(float(np.mean(cum_ac)), float(np.mean(cum_bc) + np.mean(cum_ab)), atol=1e-12)

    # First-order stochastic dominance / quantile monotonicity at P30, P50, P70
    for q in (30.0, 50.0, 70.0):
        q_bc = float(np.percentile(cum_bc, q))
        q_ab = float(np.percentile(cum_ab, q))
        q_ac = float(np.percentile(cum_ac, q))
        assert q_bc >= 0.0
        assert q_ab >= 0.0
        assert q_ac >= q_ab >= 0.0
        assert q_ac >= q_bc >= 0.0

