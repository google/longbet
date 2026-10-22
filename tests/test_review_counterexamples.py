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

"""Regression tests for all mathematical counterchecks in review-2026-09-20."""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np

from longbet._config import LongBetConfig
from longbet._state import init_longbet
from longbet._step import compute_treatment_scale_attrs
from longbet._trend_basis import build_unit_features, init_horseshoe, sample_horseshoe


def test_horseshoe_unclipped_exact_conditional():
    """Finding 1: global_aux ~ IG(1, 1/A^2 + 1/tau^2) is not clipped at 25."""
    q_inter = 138
    global_scale = 2.0 / np.sqrt(q_inter)
    A2 = global_scale**2  # 4 / 138
    hs = init_horseshoe(q_inter, global_scale)
    coef = jnp.zeros(q_inter, dtype=jnp.float32)

    keys = jax.random.split(jax.random.key(42), 2000)

    def one_draw(k):
        _, _, _, g_aux = sample_horseshoe(
            k,
            coef,
            jnp.asarray(hs["local2"]),
            jnp.asarray(hs["local_aux"]),
            jnp.asarray(hs["global2"]),
            jnp.asarray(hs["global_aux"]),
            global_scale,
        )
        return g_aux

    draws = np.asarray(jax.vmap(one_draw)(keys))
    assert np.max(draws) > 25.0
    assert not np.any(draws == 25.0)
    # Since b = 1/A^2 + 1/global2_new >= 1/A^2 = 34.5, P(xi > 25) >= 1 - exp(-34.5 / 25) = 0.7484
    assert np.mean(draws > 25.0) > 0.70


def test_unbalanced_collapsed_mu_precision_oracle():
    """Finding 2: two-unit unbalanced panel [1, 9] matches exact Gaussian oracle."""
    unit_idx = jnp.array([0] + [1] * 9, dtype=jnp.int32)
    obs_mask = jnp.ones(10, dtype=jnp.bool_)
    y = jnp.array([1.0] + [0.0] * 9, dtype=jnp.float32)
    X_unified = jnp.zeros((1, 10), dtype=jnp.uint8)
    cfg = LongBetConfig(num_trees_pr=1, num_trees_trt=1, random_intercept=True)
    st = init_longbet(
        X_unified=X_unified,
        y=y,
        unit_idx=unit_idx,
        time_idx=jnp.arange(10, dtype=jnp.int32),
        exposure_idx=jnp.zeros(10, dtype=jnp.int32),
        z_vec=jnp.zeros(10, dtype=jnp.float32),
        obs_mask=obs_mask,
        max_split_mu=jnp.zeros(1, dtype=jnp.uint8),
        max_split_nu=jnp.zeros(1, dtype=jnp.uint8),
        config=cfg,
    )
    assert st.has_unbalanced_units is True
    assert st.collapse_gamma_mu is True

    sigma2 = 1.0
    sigma_gamma2 = 1.0
    kappa_mu = jnp.reciprocal(1.0 + st.unit_counts * (sigma_gamma2 / sigma2))
    attrs_mu, _ = compute_treatment_scale_attrs(
        jnp.sqrt(kappa_mu[st.unit_idx]), obs_mask, eps=0.0
    )
    u_mu_sum = jnp.zeros(2, dtype=jnp.float32).at[st.unit_idx].add(y)
    u_mu_bar = u_mu_sum / st.unit_counts
    resid_mu = u_mu_bar[st.unit_idx]

    phys_prec = attrs_mu.prec_scale * jnp.square(attrs_mu.inv_sdev_unit)
    data_prec = float(jnp.sum(phys_prec))
    score = float(jnp.sum(phys_prec * resid_mu))
    prior_prec = 1.0
    post_var = 1.0 / (prior_prec + data_prec)
    post_mean = post_var * score

    np.testing.assert_allclose(data_prec, 1.4, rtol=1e-6)
    np.testing.assert_allclose(post_var, 5.0 / 12.0, rtol=1e-6)
    np.testing.assert_allclose(post_mean, 0.5 / 2.4, rtol=1e-6)


def test_time_varying_predictors_disable_collapsed_leaf_likelihood():
    """Finding 3: time-varying predictors or calendar/exposure splits disable collapsing."""
    unit_idx = jnp.array([0, 0, 1, 1], dtype=jnp.int32)
    time_idx = jnp.array([0, 1, 0, 1], dtype=jnp.int32)
    obs_mask = jnp.ones(4, dtype=jnp.bool_)
    y = jnp.zeros(4, dtype=jnp.float32)

    # Time-varying X_unified within unit 0
    X_tv = jnp.array([[0, 1, 0, 0]], dtype=jnp.uint8)
    st_tv = init_longbet(
        X_unified=X_tv,
        y=y,
        unit_idx=unit_idx,
        time_idx=time_idx,
        exposure_idx=jnp.zeros(4, dtype=jnp.int32),
        z_vec=jnp.zeros(4, dtype=jnp.float32),
        obs_mask=obs_mask,
        max_split_mu=jnp.array([1], dtype=jnp.uint8),
        max_split_nu=jnp.array([1], dtype=jnp.uint8),
        config=LongBetConfig(random_intercept=True),
    )
    assert st_tv.collapse_gamma_mu is False
    assert st_tv.collapse_gamma_nu is False

    # Static X_unified within every unit, but split_calendar_mu=True / split_calendar_trt=True
    X_static = jnp.array([[0, 0, 1, 1]], dtype=jnp.uint8)
    st_cal = init_longbet(
        X_unified=X_static,
        y=y,
        unit_idx=unit_idx,
        time_idx=time_idx,
        exposure_idx=jnp.zeros(4, dtype=jnp.int32),
        z_vec=jnp.zeros(4, dtype=jnp.float32),
        obs_mask=obs_mask,
        max_split_mu=jnp.array([1], dtype=jnp.uint8),
        max_split_nu=jnp.array([1], dtype=jnp.uint8),
        config=LongBetConfig(
            random_intercept=True, split_calendar_mu=True, split_calendar_trt=True
        ),
    )
    assert st_cal.collapse_gamma_mu is False
    assert st_cal.collapse_gamma_nu is False


def test_unit_features_location_scale_invariance():
    """Finding 4: build_unit_features is strictly invariant to affine rescaling/shifting."""
    rng = np.random.default_rng(123)
    x = rng.normal(size=(80, 5))
    scales = np.array([2.5, 0.4, 10.0, 1.7, 0.1])
    shifts = np.array([-3.0, 5.0, 12.0, -0.8, 100.0])
    x_transformed = x * scales[None, :] + shifts[None, :]

    B1, _ = build_unit_features(x, num_random_features=64)
    B2, _ = build_unit_features(x_transformed, num_random_features=64)
    assert B1.shape == (80, 69)
    np.testing.assert_allclose(B1, B2, rtol=1e-5, atol=1e-5)


def test_mixed_collapse_gamma_mu_true_nu_false_preserves_posterior():
    """Verify PCG target invariance when collapse_gamma_mu=True and collapse_gamma_nu=False."""
    # 1. Exact Gaussian covariance propagation check for (mu, nu, gamma):
    #    y = mu + nu + gamma + eps = 0 with independent N(0, 1) priors & noise.
    #    Target posterior covariance is Sigma = I - 11'/4 (diag = 0.75).
    Sigma = np.eye(3) - 0.25 * np.ones((3, 3))
    # Step 1: mu' ~ N(-nu/3, 2/3) (gamma marginalized)
    A1 = np.array([[0.0, -1.0 / 3.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    Q1 = np.diag([2.0 / 3.0, 0.0, 0.0])
    # Step 1b (intermediate gamma refresh): gamma_pre ~ N(-(mu' + nu)/2, 1/2)
    A1b = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-0.5, -0.5, 0.0]])
    Q1b = np.diag([0.0, 0.0, 0.5])
    # Step 3: nu' ~ N(-(mu' + gamma_pre)/2, 1/2)
    A3 = np.array([[1.0, 0.0, 0.0], [-0.5, 0.0, -0.5], [0.0, 0.0, 1.0]])
    Q3 = np.diag([0.0, 0.5, 0.0])
    # Step 3b (post-nu gamma refresh): gamma' ~ N(-(mu' + nu')/2, 1/2)
    A3b = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-0.5, -0.5, 0.0]])
    Q3b = np.diag([0.0, 0.0, 0.5])

    S = A1 @ Sigma @ A1.T + Q1
    S = A1b @ S @ A1b.T + Q1b
    np.testing.assert_allclose(S, Sigma, rtol=1e-12, atol=1e-12)
    S = A3 @ S @ A3.T + Q3
    np.testing.assert_allclose(S, Sigma, rtol=1e-12, atol=1e-12)
    S = A3b @ S @ A3b.T + Q3b
    np.testing.assert_allclose(S, Sigma, rtol=1e-12, atol=1e-12)

    # 2. Verify init_longbet with static mu predictor (row 0) and time-varying nu
    #    predictor (row 1) sets collapse_gamma_mu=True, collapse_gamma_nu=False,
    #    and that longbet_step refreshes gamma and maintains exact residual identity.
    from longbet._step import longbet_step

    unit_idx = jnp.array([0, 0, 1, 1], dtype=jnp.int32)
    time_idx = jnp.array([0, 1, 0, 1], dtype=jnp.int32)
    exposure_idx = jnp.array([0, 1, 0, 1], dtype=jnp.int32)
    z_vec = jnp.array([0.0, 1.0, 0.0, 1.0], dtype=jnp.float32)
    obs_mask = jnp.ones(4, dtype=jnp.bool_)
    y = jnp.array([0.5, -0.2, 0.3, 0.8], dtype=jnp.float32)
    # Row 0: unit-constant [0, 0, 1, 1]; Row 1: time-varying [0, 1, 0, 1]
    X_mixed = jnp.array([[0, 0, 1, 1], [0, 1, 0, 1]], dtype=jnp.uint8)
    st_mixed = init_longbet(
        X_unified=X_mixed,
        y=y,
        unit_idx=unit_idx,
        time_idx=time_idx,
        exposure_idx=exposure_idx,
        z_vec=z_vec,
        obs_mask=obs_mask,
        max_split_mu=jnp.array([1, 0], dtype=jnp.uint8),
        max_split_nu=jnp.array([1, 1], dtype=jnp.uint8),
        config=LongBetConfig(num_trees_pr=2, num_trees_trt=2, random_intercept=True),
    )
    assert st_mixed.collapse_gamma_mu is True
    assert st_mixed.collapse_gamma_nu is False

    st_next = longbet_step(jax.random.key(123), st_mixed)
    b_z = jnp.where(st_next.z_vec == 1.0, st_next.b1, st_next.b0)
    w = b_z * st_next.beta[st_next.exposure_idx]
    expected_resid = (
        y
        - st_next.alpha * st_next.mu_fit
        - (st_next.trend_design @ st_next.trend_coef)
        - w * st_next.nu_fit
        - st_next.gamma[unit_idx]
    )
    np.testing.assert_allclose(np.asarray(st_next.resid), np.asarray(expected_resid), rtol=1e-5, atol=1e-5)

