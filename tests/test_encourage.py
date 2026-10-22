"""Design validation and observable encouragement diagnostics."""

import numpy as np
import pytest

from longbet import (
    derive_exposure,
    encouragement_summary,
    plot_encouragement,
    validate_encouragement,
)


def panel():
    z = np.zeros((6, 4))
    z[:3, 1:] = 1
    d = np.array([[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 0, 1],
                  [1, 1, 1, 1], [0, 0, 0, 0], [0, 0, 1, 1]])
    return z, d


def test_validation_metadata_and_perfect_compliance_are_supported():
    z, _ = panel()
    meta = validate_encouragement(z, z)
    assert meta == {
        "design": "single_wave_complete", "n_units": 6, "n_periods": 4,
        "n_encouraged": 3, "n_control": 3, "encouragement_period": 2.0,
        "encouragement_index": 1, "has_pre_periods": True,
        "perfect_compliance": True,
    }


def test_single_period_randomization_needs_no_baseline():
    z = np.array([[1], [1], [0], [0]])
    meta = validate_encouragement(z, z, t=np.array(7.0))
    assert meta["has_pre_periods"] is False
    assert meta["encouragement_period"] == 7
    assert encouragement_summary(z, z, 7).horizon.tolist() == [1]


@pytest.mark.parametrize("dtype", [bool, np.uint8, np.uint64, float])
@pytest.mark.parametrize("name", ["z", "d"])
def test_reversals_are_rejected_for_boolean_and_unsigned_panels(dtype, name):
    z, d = panel()
    arrays = {"z": z.astype(dtype), "d": d.astype(dtype)}
    arrays[name][0, :] = [0, 1, 0, 1]
    with pytest.raises(ValueError, match=rf"{name}: .*switches back") as exc:
        validate_encouragement(**arrays)
    assert "Drop or reshape" not in str(exc.value)


@pytest.mark.parametrize("name", ["z", "d"])
@pytest.mark.parametrize("value", [0.5, np.nan, np.inf, -1])
def test_bad_binary_values_name_the_input(name, value):
    z, d = panel()
    arrays = {"z": z.astype(float), "d": d.astype(float)}
    arrays[name][0, 0] = value
    with pytest.raises(ValueError, match=rf"{name}.*binary"):
        validate_encouragement(**arrays)


@pytest.mark.parametrize("bad", [np.zeros(3), np.zeros((0, 4)),
                                  np.zeros((4, 0)), np.zeros((2, 2, 2))])
def test_bad_panel_shapes(bad):
    with pytest.raises(ValueError, match="nonempty"):
        validate_encouragement(bad, bad)


def test_mismatched_shapes_and_nonnumeric_panels():
    z, d = panel()
    with pytest.raises(ValueError, match="d has shape"):
        validate_encouragement(z, d[:2])
    with pytest.raises(ValueError, match="d.*numeric"):
        validate_encouragement(z, d.astype(str))


@pytest.mark.parametrize("t", [[1, 2], [[1, 2, 3, 4]], [1, 2, np.nan, 4],
                                [1, 2, 2, 4], [4, 3, 2, 1], [1, 1.5, 2, 3]])
def test_bad_times(t):
    with pytest.raises(ValueError, match="t|time"):
        validate_encouragement(*panel(), t=t)


def test_missing_arms_and_multiple_waves_are_explicit_errors():
    z, d = panel()
    for arm in (np.zeros_like(z), np.ones_like(z)):
        with pytest.raises(ValueError, match="both an encouraged arm"):
            validate_encouragement(arm, d)
    z[0, 1] = 0
    with pytest.raises(ValueError, match="single encouragement wave"):
        validate_encouragement(z, d)


def test_summary_arm_membership_moments_and_observed_lags():
    z, d = panel()
    frame = encouragement_summary(z, d)
    assert frame.post_encouragement.tolist() == [False, True, True, True]
    # Assigned arms remain meaningful before encouragement; both have one adopter.
    assert frame.iloc[0].takeup_encouraged == pytest.approx(1 / 3)
    assert frame.iloc[0].takeup_control == pytest.approx(1 / 3)
    row = frame.iloc[1]
    assert row.first_stage == pytest.approx(1 / 3)
    assert row.first_stage_se == pytest.approx(np.sqrt(2 / 9))
    assert row.n_encouraged_adopters == 2
    assert row.share_immediate == 0.5
    assert row.share_pre_encouragement_adoption == 0.5
    assert row.median_observed_lag == -0.5
    row = frame.iloc[-1]
    assert row.share_immediate == pytest.approx(1 / 3)
    assert row.share_pre_encouragement_adoption == pytest.approx(1 / 3)
    assert row.median_observed_lag == 0
    assert not any("complier" in name or "passed" in name for name in frame.columns)


def test_negative_sample_first_stage_is_retained():
    z, d = panel()
    d = d[[3, 4, 5, 0, 1, 2]]
    frame = encouragement_summary(z, d)
    assert frame.iloc[1].first_stage < 0
    assert frame.iloc[1].negative_first_stage


def test_empty_lag_denominator_is_missing_and_no_variance_is_invented():
    z = np.array([[0, 1], [0, 0], [0, 0]])
    frame = encouragement_summary(z, np.zeros_like(z))
    assert frame.n_encouraged_adopters.eq(0).all()
    assert frame.share_immediate.isna().all()
    assert frame.median_observed_lag.isna().all()
    assert frame.first_stage.eq(0).all()
    assert frame.first_stage_se.isna().all()
    assert frame.first_stage_includes_zero.isna().all()


def test_exposure_uses_model_clock_and_lags_use_calendar_time():
    z, d = panel()
    t = np.array([2, 5, 9, 12])
    frame = encouragement_summary(z, d, t)
    np.testing.assert_array_equal(frame.horizon, derive_exposure(z, t)[0])
    np.testing.assert_array_equal(frame.horizon, [0, 3, 7, 10])
    assert frame.iloc[1].median_observed_lag == -1.5


def test_fractional_calendar_origin_preserves_whole_unit_gaps():
    z, d = panel()
    frame = encouragement_summary(z, d, t=[0.2, 1.2, 2.2, 3.2])
    assert frame.horizon.tolist() == [0, 1, 2, 3]


def test_extreme_valid_alpha_has_finite_normal_bounds():
    frame = encouragement_summary(*panel(), alpha=np.nextafter(0.0, 1.0))
    assert np.isfinite(frame.first_stage_lower).all()
    assert np.isfinite(frame.first_stage_upper).all()


@pytest.mark.parametrize("alpha", [0, 1, -0.5, np.nan, np.inf, True, [0.05], "bad"])
def test_invalid_confidence_level(alpha):
    with pytest.raises(ValueError, match="alpha"):
        encouragement_summary(*panel(), alpha=alpha)


def test_plot_has_assigned_arm_curves_and_encouragement_marker():
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    z, d = panel()
    ax = plot_encouragement(z, d)
    np.testing.assert_allclose(ax.lines[0].get_ydata(), [1/3, 2/3, 2/3, 1])
    np.testing.assert_allclose(ax.lines[1].get_ydata(), [1/3, 1/3, 2/3, 2/3])
    np.testing.assert_allclose(ax.lines[2].get_xdata(), [2, 2])
    assert len(ax.collections) == 1
    plt.close(ax.figure)
