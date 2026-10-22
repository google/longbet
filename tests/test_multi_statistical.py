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

"""Statistical properties, shock recovery, conditional variance, and full-SUR distinctions."""

import numpy as np
import pytest
import scipy.stats as stats

from longbet import LongBetConfig, LongBetMulti, outcome_correlation


@pytest.mark.slow
def test_mixed_three_outcome_sign_changing_recovery():
    """Learn benefits AND harms from observed covariates in a randomized panel.

    This high-signal regression test is not the chapter DGP or a calibration
    study. All outcomes have positive and negative conditional effects.
    """
    from longbet import effect_draws, joint_prob

    rng = np.random.default_rng(20260907)
    n, t = 240, 8
    x = rng.normal(size=(n, 2)).astype(np.float32)
    sign = np.where(x[:, 0] >= 0, 1.0, -1.0)
    z = np.zeros((n, t), dtype=np.float32)
    treated = rng.permutation(n)[: n // 2]
    z[treated, 3:] = 1
    u, v, w = rng.normal(size=(3, n, t))
    baseline = 0.25 * x[:, 1, None] + 0.04 * np.arange(t)
    intercept = rng.normal(0, 0.15, (n, 1))
    mu_q = -0.5 + 0.2 * x[:, 1, None] + 0.04 * np.arange(t) + intercept
    truth = {
        "gmv": 0.8 * sign,
        "hours": -0.6 * sign,
        "complaint": stats.norm.cdf(mu_q[:, -1] - sign) - stats.norm.cdf(mu_q[:, -1]),
    }
    y = {
        "gmv": baseline + intercept + 0.8 * sign[:, None] * z + 0.3 * u,
        "hours": baseline + intercept - 0.6 * sign[:, None] * z + 0.3 * (0.6 * u + 0.8 * v),
        "complaint": (mu_q - sign[:, None] * z + 0.45 * u + 0.1 * v + np.sqrt(0.7875) * w > 0).astype(float),
    }
    cfg = LongBetConfig(
        sigma_prior_a=2, sigma_prior_b=1,
        num_burnin=600, num_sweeps=120, n_skip=2, num_chains=2,
        num_trees_pr=10, num_trees_trt=10, max_depth_pr=5, max_depth_trt=5,
        min_points_per_leaf_pr=5, min_points_per_leaf_trt=5,
        lambda_knl=2, random_seed=314159, device="cpu",
    )
    model = LongBetMulti(cfg).fit(
        y, x, z, outcome={"gmv": "continuous", "hours": "continuous", "complaint": "binary"},
    )
    z_eval = np.zeros_like(z)
    z_eval[:, 3:] = 1
    pred = model.predict(x, z_eval)
    for name in y:
        draws = effect_draws(pred, name)[:, -1, :]
        estimate = draws.mean(axis=-1)
        accuracy = np.mean(np.sign(estimate) == np.sign(truth[name]))
        rmse = np.sqrt(np.mean((estimate - truth[name]) ** 2))
        print(f"{name}: sign accuracy={accuracy:.3f}, RMSE={rmse:.3f}")
        assert accuracy >= 0.85, (name, accuracy)
        assert rmse < (0.20 if name == "complaint" else 0.35), (name, rmse)
        for group_sign in (1, -1):
            group = sign == group_sign
            group_draws = draws[group].mean(axis=0)
            correct = np.sign(truth[name][group].mean())
            prob = np.mean(group_draws * correct > 0)
            print(f"  group {group_sign:+d}: P(correct sign)={prob:.3f}")
            assert prob >= 0.90, (name, group_sign, prob)

    conditions = {"gmv": lambda a: a > 0, "hours": lambda a: a < 0,
                  "complaint": lambda a: a < 0}
    joint = joint_prob(pred, conditions)
    manual = np.mean(
        (effect_draws(pred, "gmv") > 0)
        & (effect_draws(pred, "hours") < 0)
        & (effect_draws(pred, "complaint") < 0), axis=-1,
    )
    np.testing.assert_allclose(joint, manual)


def test_full_sur_downstream_precision_derivation():
    """Every conditional includes the likelihood of downstream equations.

    For e0=y0-f0 and e1=y1-f1-g*e0, the f0 likelihood precision is
    1/v0 + g**2/v1. The old recursive update used only 1/v0. The final equation's
    conditional, in contrast, is valid given the earlier residual realization.
    """
    import jax
    import jax.numpy as jnp

    v0, v1, g = 1.2, 0.8, 0.65
    y = jnp.array([0.4, -0.2])

    def nll(f):
        r = y - f
        return 0.5 * (r[0] ** 2 / v0 + (r[1] - g * r[0]) ** 2 / v1)

    hessian = np.asarray(jax.hessian(nll)(jnp.zeros(2)))
    expected = np.array([[1 / v0 + g**2 / v1, -g / v1],
                         [-g / v1, 1 / v1]])
    np.testing.assert_allclose(hessian, expected, rtol=1e-6)
    assert hessian[0, 0] > 1 / v0
    np.testing.assert_allclose(hessian[1, 1], 1 / v1)
    from longbet._sur import conditional_residual
    G = jnp.array([[0., 0.], [g, 0.]])
    for m in range(2):
        r, precision = conditional_residual(
            y[:, None], jnp.ones((2, 1), bool), G, jnp.array([v0, v1]), m)
        np.testing.assert_allclose(precision[0], hessian[m, m], rtol=1e-6)
        np.testing.assert_allclose(r[0]*precision[0], (hessian@y)[m], rtol=1e-6)


def test_joint_prob_gaussian_calibration():
    """Verify joint_prob against scipy.stats.multivariate_normal CDF on known Gaussian draws."""
    N, T, D = 1, 1, 100000
    rho = 0.5
    cov = np.array([[1.0, rho], [rho, 1.0]])

    rng = np.random.default_rng(42)
    # Generate joint draws of treatment effects
    draws = rng.multivariate_normal([0.0, 0.0], cov, size=D).T  # (2, D)
    eff1 = draws[0].reshape(N, T, D)
    eff2 = draws[1].reshape(N, T, D)

    # Event: eff1 > 0.5 AND eff2 > 0.2
    c1 = 0.5
    c2 = 0.2

    from longbet._multi_model import reduce_joint_masks
    mask1 = eff1 > c1
    mask2 = eff2 > c2
    emp_p = float(reduce_joint_masks([mask1, mask2])[0, 0])

    # Analytic probability from bivariate normal
    mvn = stats.multivariate_normal(mean=[0.0, 0.0], cov=cov)
    # P(X1 > c1, X2 > c2) = 1 - P(X1 <= c1) - P(X2 <= c2) + P(X1 <= c1, X2 <= c2)
    p_analytic = 1.0 - stats.norm.cdf(c1) - stats.norm.cdf(c2) + mvn.cdf([c1, c2])

    # Monte Carlo standard error: sqrt(p * (1 - p) / D)
    se = np.sqrt(p_analytic * (1.0 - p_analytic) / D)
    np.testing.assert_allclose(emp_p, p_analytic, atol=3.5 * se)


def test_all_binary_no_coupling():
    """Verify all-binary models run with Gamma held strictly at zero."""
    N, T = 6, 4
    rng = np.random.default_rng(888)
    y_b1 = rng.choice([0.0, 1.0], size=(N, T)).astype(np.float32)
    y_b2 = rng.choice([0.0, 1.0], size=(N, T)).astype(np.float32)
    X = rng.standard_normal((N, 2)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)

    cfg = LongBetConfig(
        num_sweeps=4,
        num_burnin=2,
        num_chains=1,
        num_trees_pr=3,
        num_trees_trt=3,
        random_seed=15,
        sur=True,
    )

    model = LongBetMulti(cfg).fit(
        y={"b1": y_b1, "b2": y_b2},
        x=X,
        z=z,
        outcome="binary",
    )

    # Both equations are binary -> Gamma loadings must be strictly 0
    np.testing.assert_allclose(model.trace.gamma_loadings, 0.0)
    R = outcome_correlation(model)
    np.testing.assert_allclose(R, np.eye(2))
