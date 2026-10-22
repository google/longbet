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

"""Design-based reference for continuous intensity: exact recovery, immunity to
within-arm confounding, delta-method calibration, weak-IV sets, FOSD and dip screens."""

import numpy as np
import pandas as pd
import pytest

from longbet import intensity_reference, lag_matrix, validate_intensity
from longbet._encourage_intensity import (
    LinearLagStage2,
    anderson_rubin_statistic,
    compose_offer_effect,
    control_function_columns,
    deconvolve,
    offer_induced_dose,
    pre_offer_demean,
    standardized_covariates,
    toeplitz_lower,
)

TRUTH = np.array([-0.3, 0.5, 0.4, 0.2, 0.1, 0.05, 0.0])
PULSE = np.array([1.0, 0.6, 0.3, 0.1, 0.0, 0.0, 0.0])
RAMP = np.array([0.1, 0.3, 0.6, 1.0, 1.0, 1.0, 1.0])


def panel(rng, *, n=400, T=10, t0=4, response=PULSE, kappa=TRUTH, noise=0.0, confound=0.0,
          feedback=0.0, defiers=0.0, strata=False):
    """Encouragement panel with a known lag kernel; ``t0`` is one-based.

    ``confound`` adds a static unit confounder to dose and outcome, ``feedback``
    lets the dose respond to the previous outcome innovation (dynamic feedback),
    and ``defiers`` makes a share of the encouraged arm *reduce* its dose.
    """
    z = np.zeros((n, T))
    perm = rng.permutation(n)
    enc = perm[: n // 2]
    z[enc, t0 - 1:] = 1
    u = rng.normal(size=(n, 1))
    base = 0.5 + confound * u + noise * 0.2 * rng.normal(size=(n, T))
    d = base.copy()
    d[enc, t0 - 1:] += response[: T - t0 + 1]
    if defiers > 0:
        bad = enc[: int(defiers * len(enc))]
        d[bad, t0 - 1:] -= 2.5 * response[: T - t0 + 1]
    eps = noise * rng.normal(size=(n, T))
    y = np.zeros((n, T))
    if feedback == 0:
        lags = lag_matrix(d, len(kappa) - 1)
        y = lags @ kappa + 2.0 * confound * u + eps
    else:
        for j in range(T):
            if j >= 2:
                d[:, j] += feedback * (y[:, j - 1] - y[:, j - 2])
            lags_j = np.stack([d[:, j - lag] if j - lag >= 0 else np.zeros(n) for lag in range(len(kappa))], axis=1)
            y[:, j] = lags_j @ kappa + 2.0 * confound * u[:, 0] + eps[:, j]
    s = (np.arange(n) % 3) if strata else None
    return y, d, z, s


def test_noiseless_pulse_and_ramp_recover_the_kernel_exactly():
    rng = np.random.default_rng(0)
    for response, well_conditioned in ((PULSE, True), (RAMP, False)):
        y, d, z, _ = panel(rng, response=response)
        ref = intensity_reference(y, d, z)
        np.testing.assert_allclose(ref.kernel.kappa.to_numpy(), TRUTH, atol=1e-9)
        np.testing.assert_allclose(ref.kernel.kappa_cum.to_numpy(), np.cumsum(TRUTH), atol=1e-9)
        post = ref.itt.loc[ref.itt.post_encouragement]
        np.testing.assert_allclose(post.itt_d.to_numpy(), response, atol=1e-9)
        np.testing.assert_allclose(post.itt_y.to_numpy(), toeplitz_lower(response) @ TRUTH, atol=1e-9)
        # omega_cum(t) = sum_l kappa_l C_D(t-l) / C_D(t) with C_D the cumulative first stage.
        c_d = np.cumsum(response)
        expected = np.array([sum(TRUTH[lag] * c_d[t - lag] for lag in range(t + 1)) / c_d[t] for t in range(7)])
        np.testing.assert_allclose(post.omega_cum.to_numpy(), expected, atol=1e-9)
        assert (ref.condition_number < 5) is well_conditioned
    assert ref.metadata["lag_horizon"] == 6 and ref.metadata["n_post_periods"] == 7


def test_inverse_toeplitz_growth_follows_the_generating_function_root():
    # The inverse of a lower-triangular Toeplitz matrix is Toeplitz with first column equal to the
    # power-series coefficients of 1 / D(x), D(x) = sum_k d_k x^k. For a two-period first stage
    # d = (d_0, d_1, 0, ...) that is (T_d^{-1})_{l,0} = (-d_1/d_0)^l / d_0: geometric growth iff
    # d_1 > d_0 (a root of D inside the unit disk).
    d0, alpha = 0.2, 1.8
    inv = np.linalg.inv(toeplitz_lower(np.array([d0, alpha * d0, 0, 0, 0, 0])))
    np.testing.assert_allclose(inv[:, 0], (-alpha) ** np.arange(6) / d0)
    # A geometric ramp d_k = d_0 alpha^k has D(x) = d_0 / (1 - alpha x), so its inverse is the
    # banded (1, -alpha) / d_0: no growth at all. The ramp is ill-conditioned through cond(T_d),
    # not through this entry.
    inv = np.linalg.inv(toeplitz_lower(d0 * alpha ** np.arange(6)))
    np.testing.assert_allclose(inv[:, 0], np.r_[1.0, -alpha, 0, 0, 0, 0] / d0, atol=1e-12)
    # A single-period pulse (d = d0 e_0) makes the inverse diagonal: kappa = y / d0.
    np.testing.assert_allclose(np.linalg.inv(toeplitz_lower(np.array([d0, 0, 0, 0]))), np.eye(4) / d0)


def test_t0_is_unbiased_under_strong_within_arm_confounding_where_ols_is_not():
    rng = np.random.default_rng(1)
    estimates, ses, ols = [], [], []
    for _ in range(24):
        y, d, z, _ = panel(rng, n=800, noise=1.0, confound=1.0)
        ref = intensity_reference(y, d, z)
        estimates.append(ref.kernel.kappa.to_numpy())
        ses.append(ref.kernel.kappa_se.to_numpy())
        design = np.column_stack([np.ones(y.size), lag_matrix(d, 6).reshape(-1, 7)])
        ols.append(np.linalg.lstsq(design, y.ravel(), rcond=None)[0][1:])
    estimates, ses, ols = map(np.asarray, (estimates, ses, ols))
    mcse = estimates.std(axis=0, ddof=1) / np.sqrt(len(estimates))
    assert np.all(np.abs(estimates.mean(axis=0) - TRUTH) < 3.5 * mcse)
    # Delta-method SEs track the Monte Carlo spread.
    ratio = ses.mean(axis=0) / estimates.std(axis=0, ddof=1)
    assert np.all((ratio > 0.6) & (ratio < 1.6))
    # Naive within-arm OLS on the lag matrix is far off at lag 0 (the confounded contemporaneous dose).
    assert abs(ols.mean(axis=0)[0] - TRUTH[0]) > 1.0


def test_measurement_error_in_the_dose_leaves_t0_unbiased():
    rng = np.random.default_rng(7)
    estimates = []
    for _ in range(16):
        y, d, z, _ = panel(rng, n=800, noise=1.0)
        noisy = d + 0.8 * rng.normal(size=d.shape)         # classical measurement error
        estimates.append(intensity_reference(y, noisy, z).kernel.kappa.to_numpy())
    estimates = np.asarray(estimates)
    mcse = estimates.std(axis=0, ddof=1) / np.sqrt(len(estimates))
    assert np.all(np.abs(estimates.mean(axis=0) - TRUTH) < 3.5 * mcse)


def test_stratified_contrast_equals_blocked_difference_in_means():
    rng = np.random.default_rng(2)
    y, d, z, s = panel(rng, n=300, noise=1.0, confound=0.5, strata=True)
    ref = intensity_reference(y, d, z, strata=s)
    assigned = z[:, -1] == 1
    expected = sum(np.mean(s == k) * (y[(s == k) & assigned].mean(axis=0) - y[(s == k) & ~assigned].mean(axis=0))
                   for k in np.unique(s))
    np.testing.assert_allclose(ref.itt.itt_y.to_numpy(), expected)
    plain = intensity_reference(y, d, z)
    np.testing.assert_allclose(plain.itt.itt_y.to_numpy(), y[assigned].mean(axis=0) - y[~assigned].mean(axis=0))
    assert ref.metadata["n_strata"] == 3 and ref.metadata["inference_available"]
    with pytest.raises(ValueError, match="both arms"):
        intensity_reference(y, d, z, strata=np.where(assigned, "a", "b"))


def test_fieller_sets_cover_the_cumulative_return_and_flag_weak_denominators():
    rng = np.random.default_rng(3)
    covered, kinds = [], []
    truth_last = (TRUTH @ np.cumsum(PULSE)[::-1]) / np.cumsum(PULSE)[-1]
    for _ in range(30):
        y, d, z, _ = panel(rng, n=600, noise=1.0, confound=0.5)
        row = intensity_reference(y, d, z).itt.query("post_encouragement").iloc[-1]
        kinds.append(row.omega_cum_set_type)
        covered.append(row.omega_cum_lower_1 <= truth_last <= row.omega_cum_upper_1)
    assert set(kinds) == {"bounded"} and np.mean(covered) >= 0.8
    # A null first stage makes the per-period set unbounded (all_real / disjoint / half_line), never a fake interval.
    y, d, z, _ = panel(rng, n=400, response=np.zeros(7), kappa=np.zeros(7), noise=1.0)
    ref = intensity_reference(y, d, z)
    assert set(ref.itt.query("post_encouragement").omega_set_type) <= {"all_real", "disjoint", "half_line"}


def test_anderson_rubin_region_is_exact_at_the_estimate_and_rejects_distant_kernels():
    rng = np.random.default_rng(4)
    y, d, z, _ = panel(rng, n=800, noise=1.0, confound=0.5)
    ref = intensity_reference(y, d, z)
    at_estimate = ref.anderson_rubin(ref.kernel.kappa)
    assert at_estimate["statistic"] == pytest.approx(0.0, abs=1e-8) and at_estimate["in_region"]
    assert ref.anderson_rubin(TRUTH)["pvalue"] > 1e-3
    far = ref.anderson_rubin(TRUTH + 3.0)
    assert not far["in_region"] and far["pvalue"] < 1e-6
    # The statistic is the quadratic form of the linear moment with its exact covariance.
    k = ref.metadata["n_post_periods"]
    stat = anderson_rubin_statistic(TRUTH, ref._post["itt_y"], ref._post["itt_d"], ref._post["covariance"])
    assert stat == pytest.approx(ref.anderson_rubin(TRUTH)["statistic"]) and at_estimate["df"] == k
    proj = ref.anderson_rubin_projection(n_candidates=400, random_seed=0)
    assert list(proj.columns[:3]) == ["lag", "kappa_lower", "kappa_upper"]
    assert np.all(proj.kappa_lower <= ref.kernel.kappa + 1e-12) and np.all(proj.kappa_upper >= ref.kernel.kappa - 1e-12)


def test_truncated_lag_horizon_uses_least_squares_with_a_consistent_jacobian():
    rng = np.random.default_rng(5)
    kappa_short = np.array([-0.3, 0.5, 0.4, 0.0, 0.0, 0.0, 0.0])
    y, d, z, _ = panel(rng, kappa=kappa_short)
    ref = intensity_reference(y, d, z, lag_horizon=2)
    np.testing.assert_allclose(ref.kernel.kappa.to_numpy(), kappa_short[:3], atol=1e-9)
    # Finite differences of the deconvolution agree with the analytic Jacobian.
    itt_y, itt_d, cov = ref._post["itt_y"], ref._post["itt_d"], np.eye(14)
    kappa, _, jac = deconvolve(itt_y, itt_d, cov, lag_horizon=2)
    numeric = np.zeros_like(jac)
    for j in range(14):
        step = np.zeros(14); step[j] = 1e-6
        plus = deconvolve(itt_y + step[:7], itt_d + step[7:], cov, lag_horizon=2)[0]
        minus = deconvolve(itt_y - step[:7], itt_d - step[7:], cov, lag_horizon=2)[0]
        numeric[:, j] = (plus - minus) / 2e-6
    np.testing.assert_allclose(jac, numeric, atol=1e-5)
    with pytest.raises(ValueError, match="lag_horizon"):
        intensity_reference(y, d, z, lag_horizon=7)


def test_ecdf_shift_and_fosd_check():
    rng = np.random.default_rng(6)
    y, d, z, s = panel(rng, n=600, noise=1.0, strata=True)
    ref = intensity_reference(y, d, z, strata=s)
    assert (ref.ecdf.difference >= -1e-12).all() and ref.fosd["sup_violation"] == pytest.approx(0.0, abs=1e-12)
    assert ref.fosd["holds_at_level"] and ref.fosd["mean_shift"] == pytest.approx(PULSE.sum(), abs=0.2)
    y, d, z, _ = panel(rng, n=600, noise=1.0, defiers=0.4)
    bad = intensity_reference(y, d, z)
    assert bad.fosd["sup_violation"] > 0.1 and bad.fosd["pvalue"] < 0.01 and not bad.fosd["holds_at_level"]


def test_dip_screen_has_size_under_no_feedback_and_power_under_feedback():
    rng = np.random.default_rng(8)
    null, alt = [], []
    for _ in range(12):
        y, d, z, _ = panel(rng, n=400, noise=1.0, confound=0.5)
        null.append(intensity_reference(y, d, z, x=rng.normal(size=(400, 2))).dip["z_stat"])
        y, d, z, _ = panel(rng, n=400, noise=1.0, confound=0.5, feedback=0.6)
        alt.append(intensity_reference(y, d, z).dip["z_stat"])
    assert np.mean(np.abs(null) > 1.96) <= 0.25 and abs(np.mean(null)) < 1.0
    assert np.mean(np.abs(alt) > 1.96) >= 0.9
    ref = intensity_reference(y, d, z)
    assert ref.dip["n_units"] == 400 and ref.dip["method"] == "ols_unit_clustered_cr1"
    assert "D_{t-1..t-8}" in ref.dip["controls"]


def test_dip_screen_does_not_fire_on_past_dose_effects_alone():
    # Autocorrelated, non-absorbing dose with a strong causal lag effect and no feedback: the
    # dose-history conditioning blocks the mechanical path D_{t-1} -> Delta Y_{t-1}.
    rng = np.random.default_rng(9)
    zs = []
    for _ in range(12):
        n, T = 400, 10
        z = np.zeros((n, T)); z[: n // 2, 3:] = 1
        d = np.zeros((n, T)); shock = rng.normal(size=(n, T))
        for j in range(T):
            d[:, j] = 0.7 * (d[:, j - 1] if j else 0) + shock[:, j] + 0.5 * z[:, j]
        y = lag_matrix(d, 6) @ TRUTH + rng.normal(size=(n, T))
        zs.append(intensity_reference(y, d, z).dip["z_stat"])
    assert np.mean(np.abs(zs) > 1.96) <= 0.25


def test_discrepancy_statistics_flag_biased_model_draws_and_pass_correct_ones():
    rng = np.random.default_rng(10)
    y, d, z, _ = panel(rng, n=800, noise=1.0)
    ref = intensity_reference(y, d, z)
    good = TRUTH[None] + 0.01 * rng.normal(size=(300, 7))
    bad = TRUTH[None] + np.array([2.0, 0, 0, 0, 0, 0, 0]) + 0.01 * rng.normal(size=(300, 7))
    ok = ref.discrepancy(good, label="T2")
    flagged = ref.discrepancy(bad, label="T1")
    assert set(ok.label) == {"T2"} and len(ok) == 7
    assert np.all(np.abs(ok.z_stat) < 3.5) and np.all(np.abs(flagged.z_stat) > 4)
    assert flagged.flag.all() and not ok.flag.iloc[-1]
    np.testing.assert_allclose(ref.discrepancy(np.cumsum(good, axis=1), cumulative=True).z_stat, ok.z_stat)


def test_validation_and_metadata():
    rng = np.random.default_rng(11)
    y, d, z, _ = panel(rng, n=40, noise=1.0)
    p = validate_intensity(z, d)
    assert p.start == 3 and p.n_post == 7 and p.metadata()["treatment"] == "intensity"
    assert p.metadata()["dose_nonabsorbing"]
    with pytest.raises(ValueError, match="d must be a finite numeric dose panel"):
        validate_intensity(z, d[:, :-1])
    with pytest.raises(ValueError, match="complete and finite"):
        validate_intensity(z, np.where(np.arange(d.size).reshape(d.shape) == 5, np.nan, d))
    with pytest.raises(ValueError, match="encouraged arm and a permanent control arm"):
        validate_intensity(np.ones_like(z), d)
    with pytest.raises(ValueError, match="consecutive unit-spaced"):
        validate_intensity(z, d, t=np.arange(0, 20, 2))
    ref = intensity_reference(y, d, z)
    assert ref.metadata["kernel_estimand"] == "compliance_weighted_population_kernel"
    assert ref.metadata["identification"] == "deconvolution"
    assert list(ref.covariance.index.names) == ["quantity", "period_index"]
    assert ref.covariance.shape == (14, 14)


def test_lag_matrix_demeaning_and_control_function_columns():
    d = np.array([[1.0, 2.0, 3.0, 4.0], [0.0, 0.0, 1.0, 0.0]])
    lags = lag_matrix(d, 2)
    assert lags.shape == (2, 4, 3)
    np.testing.assert_array_equal(lags[0, :, 0], d[0])
    np.testing.assert_array_equal(lags[0, :, 1], [0, 1, 2, 3])
    np.testing.assert_array_equal(lags[0, :, 2], [0, 0, 1, 2])
    np.testing.assert_array_equal(pre_offer_demean(d, 2)[0], d[0] - 1.5)
    np.testing.assert_array_equal(pre_offer_demean(d, 0), d)
    z = np.array([[0, 0, 1, 1], [0, 0, 0, 0]])
    np.testing.assert_array_equal(offer_induced_dose(np.full((2, 4), 0.5), z, center=False), 0.5 * z)
    np.testing.assert_allclose(
        offer_induced_dose(np.full((2, 4), 0.5), z),
        np.array([[0.0, 0.0, 0.25, 0.25], [0.0, 0.0, -0.25, -0.25]]),
    )
    v = np.array([[1.0, -1.0, 2.0, 0.0], [0.5, 0.5, -0.5, -0.5]])
    x_tilde = standardized_covariates(np.array([[1.0, 3.0], [2.0, 3.0]]))
    assert x_tilde.shape == (2, 1)                           # the constant column drops
    cols, names = control_function_columns(v, 2, x_tilde, start=2)
    assert cols.shape == (2, 4, 6) and names == ["v", "v:x0", "v_lag", "v_lag:x0", "v_bar", "v_bar:x0"]
    np.testing.assert_array_equal(cols[0, :, 0], v[0])
    np.testing.assert_array_equal(cols[0, :, 2], [0, 1, 0, 1])              # v_{t-1} + v_{t-2}
    np.testing.assert_allclose(cols[0, :, 4], np.full(4, 1.0))               # post-offer mean of v
    np.testing.assert_allclose(cols[1, :, 4], np.full(4, -0.5))
    np.testing.assert_allclose(cols[:, :, 1], cols[:, :, 0] * x_tilde[:, 0][:, None])
    with pytest.raises(ValueError, match="lag_horizon"):
        lag_matrix(d, -1)


def test_composition_matches_the_toeplitz_convolution_and_abel_summation():
    rng = np.random.default_rng(12)
    kappa = rng.normal(size=(5, 4))                          # 5 draws, L = 3
    citt_d = np.zeros((3, 9, 5))
    citt_d[:, 3:, :] = rng.normal(size=(3, 6, 5))
    out = compose_offer_effect(kappa, citt_d, 3)
    for b in range(5):
        for i in range(3):
            np.testing.assert_allclose(out[i, 3:, b], toeplitz_lower(np.r_[kappa[b], 0, 0])[:6, :6] @ citt_d[i, 3:, b])
    assert np.all(out[:, :3, :] == 0)
    # Abel summation: an absorbing adoption stock F(t) = sum_{a<=t} pi_a with kappa_l = tau(l+1) - tau(l)
    # composes to sum_a pi_a tau(t - a + 1), the binary paper's composition.
    tau = np.array([0.0, 0.4, 0.7, 0.9, 1.0, 1.05, 1.08])
    kappa_bin = np.diff(tau)[None]                           # tau(l+1) - tau(l), l = 0..5
    pi = np.array([0.3, 0.2, 0.1, 0.05, 0.05, 0.0])          # adoption at a = t0..t0+5
    stock = np.zeros((1, 9, 1)); stock[0, 3:, 0] = np.cumsum(pi)
    composed = compose_offer_effect(kappa_bin, stock, 3)[0, 3:, 0]
    expected = np.array([sum(pi[a] * tau[t - a + 1] for a in range(t + 1)) for t in range(6)])
    np.testing.assert_allclose(composed, expected)
    # Per-unit kernels (N, B, L+1) broadcast the moderator nu(x).
    nu = np.array([1.0, 2.0, 0.5])
    scaled = compose_offer_effect(kappa[None] * nu[:, None, None], citt_d, 3)
    np.testing.assert_allclose(scaled, out * nu[:, None, None])


def test_three_term_control_function_recovers_the_kernel_where_two_terms_fail():
    """Replicates the mechanism of the JTPA audit (Plan v3 Design C10): outcome shocks load on the
    contemporaneous dose innovation and on the innovations of the active lag window, whose length
    grows with calendar time. Controlling only for the contemporaneous residual and its unit mean
    leaves lags 1..L confounded and pushes the gap into kappa_0; the cumulative lagged residual
    ``v_lag`` restores the kernel. Static selection alone is absorbed by the unit effects."""
    rng = np.random.default_rng(13)
    n, T, t0 = 3000, 10, 4
    kappa = np.array([-0.3, 0.5, 0.4, 0.2, 0.1, 0.05, 0.0])
    offer = np.zeros((n, T)); offer[: n // 2, t0 - 1:] = PULSE[: T - t0 + 1]
    x = rng.normal(size=(n, 2))
    e = rng.normal(size=(n, T))                                 # dose innovations
    d = 0.5 + 0.5 * e + offer
    e_lag = lag_matrix(e, 6)[:, :, 1:].sum(axis=2)
    y = lag_matrix(d, 6) @ kappa + 0.8 * e + 0.6 * e_lag + 0.3 * rng.normal(size=(n, T))
    v = d - (0.5 + offer)                                       # the first stage is known here
    lags = lag_matrix(d, 6)
    three, names = control_function_columns(v, 6, standardized_covariates(x), start=t0 - 1)
    two = three[:, :, [k for k, name in enumerate(names) if not name.startswith("v_lag")]]
    fit3 = LinearLagStage2().fit(y, lags, three, num_draws=200, num_chains=1, random_seed=0, extra_names=names)
    fit2 = LinearLagStage2().fit(y, lags, two, num_draws=200, num_chains=1, random_seed=0)
    fit0 = LinearLagStage2().fit(y, lags, None, num_draws=200, num_chains=1, random_seed=0)
    err3 = np.abs(fit3.kappa_draws.mean(axis=0) - kappa)
    err2 = np.abs(fit2.kappa_draws.mean(axis=0) - kappa)
    err0 = np.abs(fit0.kappa_draws.mean(axis=0) - kappa)
    assert err3.max() < 0.05, err3
    assert err2[0] > 0.3 and err2[1:].min() > 0.3, err2        # kappa_0 and every lag stay confounded
    assert err0[0] > 0.3, err0
    # The unit-mean term is time invariant and absorbed by the unit fixed effects of the backend.
    assert set(fit3.absorbed_columns) >= {7 + names.index("v_bar"), 7 + names.index("v_bar:x0")}
