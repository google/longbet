"""Reference IV inference checked against directly residualized outcome tests."""

import numpy as np
import pytest
from scipy.stats import norm

from longbet import encouragement_effects
from longbet._encourage import _quadratic_set


def _contains(row, candidates):
    candidates = np.asarray(candidates)
    included = np.zeros(candidates.shape, dtype=bool)
    for k in (1, 2):
        lo, hi = row[f"wald_lower_{k}"], row[f"wald_upper_{k}"]
        if not np.isnan(lo):
            included |= (candidates >= lo) & (candidates <= hi)
    return included


def _experiment(seed=1, n=100):
    rng = np.random.default_rng(seed)
    z = np.zeros((n, 3))
    z[:n // 2, 1:] = 1
    # Vary the empirical first stage, including negative and zero contrasts.
    p = (0.15, 0.4, 0.8)[seed % 3]
    uptake = rng.random(n) < np.where(z[:, -1], p, 0.4)
    d = np.zeros_like(z)
    d[:, 1:] = uptake[:, None]
    y = rng.normal(size=z.shape) + 2 * d + (seed % 4) * z
    return y, d, z


def test_joint_neyman_moments_against_hand_calculation():
    z = np.array([[1], [1], [1], [0], [0], [0]])
    d = np.array([[1], [1], [0], [1], [0], [0]])
    y = np.array([[1], [3], [5], [2], [4], [8]])
    row = encouragement_effects(y, d, z).iloc[0]
    assert row.itt_y == pytest.approx(-5 / 3)
    assert row.itt_d == pytest.approx(1 / 3)
    assert row.itt_y_se**2 == pytest.approx(40 / 9)
    assert row.itt_d_se**2 == pytest.approx(2 / 9)
    assert row.itt_y_d_cov == pytest.approx(-7 / 9)
    assert row.wald == pytest.approx(-5)
    assert row.inference == "normal_ar"


@pytest.mark.parametrize("seed", range(15))
def test_confidence_set_matches_direct_studentized_residual_contrast(seed):
    y, d, z = _experiment(seed)
    frame = encouragement_effects(y, d, z, alpha=0.1)
    candidates = np.r_[np.linspace(-30, 30, 601), -1e8, 1e8]
    assigned = z[:, -1] == 1
    critical = norm.isf(0.05)
    for row in frame.to_dict("records"):
        j = row["period_index"]
        residual = y[:, j, None] - d[:, j, None] * candidates[None, :]
        a, b = residual[assigned], residual[~assigned]
        difference = a.mean(0) - b.mean(0)
        variance = a.var(0, ddof=1) / len(a) + b.var(0, ddof=1) / len(b)
        direct = difference**2 <= critical**2 * variance
        np.testing.assert_array_equal(_contains(row, candidates), direct)


@pytest.mark.parametrize("coefficients,kind,intervals", [
    ((1, 0, -1), "bounded", [(-1, 1)]),
    ((-1, 0, 1), "disjoint", [(-np.inf, -1), (1, np.inf)]),
    ((-1, 0, -1), "all_real", [(-np.inf, np.inf)]),
    ((1, 0, 1), "empty", []),
    ((0, 2, -4), "half_line", [(-np.inf, 2)]),
    ((0, -2, 4), "half_line", [(2, np.inf)]),
    ((0, 0, -1), "all_real", [(-np.inf, np.inf)]),
    ((0, 0, 0), "all_real", [(-np.inf, np.inf)]),
    ((0, 0, 1), "empty", []),
    ((1, -4, 4), "singleton", [(2, 2)]),
    ((-1, 4, -4), "all_real", [(-np.inf, np.inf)]),
])
def test_quadratic_boundary_and_unbounded_cases(coefficients, kind, intervals):
    actual_kind, actual_intervals = _quadratic_set(*coefficients)
    assert actual_kind == kind
    np.testing.assert_allclose(actual_intervals, intervals)


def test_near_linear_quadratic_keeps_far_endpoint():
    kind, intervals = _quadratic_set(1e-16, 1, -1)
    assert kind == "bounded"
    lo, hi = intervals[0]
    assert lo == pytest.approx(-1e16)
    assert hi == pytest.approx(1)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4, 5])
@pytest.mark.parametrize("scale", [1e-12, -1e-12, 1e12, -1e12])
def test_confidence_sets_respect_outcome_units_and_sign(seed, scale):
    y, d, z = _experiment(seed)
    original = encouragement_effects(y, d, z)
    changed = encouragement_effects(y * scale, d, z)
    np.testing.assert_array_equal(original.wald_set_type, changed.wald_set_type)
    for before, after in zip(original.to_dict("records"), changed.to_dict("records")):
        candidates = np.linspace(-100, 100, 401)
        np.testing.assert_array_equal(_contains(before, candidates),
                                      _contains(after, candidates * scale))
        assert after["wald"] / scale == pytest.approx(before["wald"], nan_ok=True)


def test_outcome_shift_does_not_change_the_estimates_or_confidence_set():
    y, d, z = _experiment(2)
    before = encouragement_effects(y, d, z)
    after = encouragement_effects(y + 100, d, z)
    for field in ("itt_y", "itt_d", "itt_y_d_cov", "wald", "wald_lower_1", "wald_upper_1"):
        np.testing.assert_allclose(before[field], after[field], equal_nan=True)


def test_perfect_compliance_reduces_to_outcome_difference_interval():
    y, _, z = _experiment(3)
    result = encouragement_effects(y, z, z)
    np.testing.assert_allclose(result.itt_d, 1)
    np.testing.assert_allclose(result.itt_d_se, 0)
    np.testing.assert_allclose(result.wald, result.itt_y)
    np.testing.assert_allclose(result.wald_lower_1, result.itt_y_lower)
    np.testing.assert_allclose(result.wald_upper_1, result.itt_y_upper)
    assert result.wald_lower_2.isna().all()


def test_binary_outcome_is_on_risk_difference_scale():
    z = np.array([[1], [1], [1], [0], [0], [0]])
    y = np.array([[1], [1], [0], [1], [0], [0]])
    row = encouragement_effects(y, z, z).iloc[0]
    assert row.itt_y == pytest.approx(1 / 3)
    assert row.wald == pytest.approx(1 / 3)


def test_zero_first_stage_retains_all_real_or_empty_set():
    z = np.array([[1], [1], [0], [0]])
    d = np.zeros_like(z)
    row = encouragement_effects(np.array([[1], [2], [1], [2]]), d, z).iloc[0]
    assert np.isnan(row.wald)
    assert row.wald_set_type == "all_real"
    assert row.wald_lower_1 == -np.inf and row.wald_upper_1 == np.inf
    row = encouragement_effects(10 * z, d, z).iloc[0]
    assert row.wald_set_type == "empty"
    assert np.isnan(row.wald_lower_1)


def test_negative_first_stage_is_not_masked():
    y, _, z = _experiment(1)
    d = np.tile((1 - z[:, -1])[:, None], (1, z.shape[1]))
    result = encouragement_effects(y, d, z)
    assert result.itt_d.eq(-1).all()
    np.testing.assert_allclose(result.wald, -result.itt_y)
    assert result.wald_set_type.eq("bounded").all()


def test_small_arm_retains_points_but_not_uncertainty():
    z = np.array([[1], [0], [0]])
    result = encouragement_effects(3 * z, z, z)
    assert result.iloc[0].wald == 3
    assert result.iloc[0].wald_set_type == "unavailable"
    assert result.iloc[0].wald_reason == "insufficient_arm_size"
    assert result.itt_y_se.isna().all()


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_missing_outcomes_are_not_silently_dropped(invalid):
    y, d, z = _experiment()
    y[0, 1] = invalid
    with pytest.raises(ValueError, match="complete and finite"):
        encouragement_effects(y, d, z)


def test_outcome_shape_and_confidence_level_validation():
    y, d, z = _experiment()
    with pytest.raises(ValueError, match="y must be.*shape"):
        encouragement_effects(y[:1], d, z)
    with pytest.raises(ValueError, match="alpha"):
        encouragement_effects(y, d, z, alpha=False)


def test_only_observed_post_periods_are_returned_with_model_exposure():
    y, d, z = _experiment()
    result = encouragement_effects(y, d, z, t=[2, 5, 9])
    assert result.period.tolist() == [5, 9]
    assert result.period_index.tolist() == [1, 2]
    assert result.horizon.tolist() == [3, 7]
