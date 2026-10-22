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

"""Describing a staggered rollout before any model is fitted."""

import numpy as np
import pandas as pd
import pytest

from longbet import (
    LongBet,
    adoption_cohorts,
    derive_exposure,
    plot_rollout,
    rollout_summary,
)


def _waves(n=200, n_periods=14, launch=(4, 6, 8, 10), holdout_share=0.4, seed=0):
    """The canonical design: several launch waves plus a never-treated holdout."""
    rng = np.random.default_rng(seed)
    p = [(1 - holdout_share) / len(launch)] * len(launch) + [holdout_share]
    assign = rng.choice(list(launch) + [0], size=n, p=p)
    z = np.zeros((n, n_periods))
    for wk in launch:
        z[assign == wk, wk - 1:] = 1.0
    return z, np.arange(1, n_periods + 1), assign


def test_cohorts_are_derived_from_z():
    z, t, assign = _waves()
    first, times = adoption_cohorts(z, t)
    np.testing.assert_array_equal(times, [4, 6, 8, 10])
    assert np.isnan(first[assign == 0]).all()
    for wk in (4, 6, 8, 10):
        assert np.allclose(first[assign == wk], wk)


def test_summary_shape_and_status_grid():
    z, t, assign = _waves()
    df = rollout_summary(z, t, labels={4: "W1", 6: "W2", 8: "W3", 10: "W4"},
                         never_treated_label="Holdout")

    assert list(df["cohort"].cat.categories) == ["W1", "W2", "W3", "W4", "Holdout"]
    assert len(df) == 5 * len(t)

    # Each cohort is untreated up to its launch week and treated from then on.
    for label, wk in (("W1", 4), ("W2", 6), ("W3", 8), ("W4", 10)):
        rows = df[df["cohort"] == label].sort_values("period")
        treated = (rows["status"] == "Treated").to_numpy()
        np.testing.assert_array_equal(treated, rows["period"].to_numpy() >= wk)

    holdout = df[df["cohort"] == "Holdout"]
    assert (holdout["status"] == "Never treated").all()
    assert (holdout["exposure"] == 0).all()

    counts = df.groupby("cohort", observed=True)["n_units"].first()
    assert counts.sum() == z.shape[0]
    assert counts["Holdout"] == int((assign == 0).sum())


def test_exposure_matches_the_model_definition():
    """The plot's event-time axis must be the sampler's, not a re-derivation."""
    z, t, _ = _waves(seed=3)
    df = rollout_summary(z, t)
    exposure = derive_exposure(z, t)

    first, _ = adoption_cohorts(z, t)
    for cohort in df["cohort"].cat.categories:
        rows = df[df["cohort"] == cohort].sort_values("period_index")
        if str(cohort) == "Never treated":
            member = int(np.argmax(np.isnan(first)))
        else:
            member = int(np.argmax(np.isclose(first, float(str(cohort)))))
        np.testing.assert_array_equal(rows["exposure"].to_numpy(), exposure[member])


def test_labels_are_optional_and_partial():
    z, t, _ = _waves()
    df = rollout_summary(z, t, labels={4: "First"})
    cats = list(df["cohort"].cat.categories)
    assert cats[0] == "First"
    assert cats[1:] == ["6", "8", "10", "Never treated"]


def test_all_treated_and_none_treated_panels():
    t = np.arange(1, 6)
    everyone = np.ones((10, 5))
    df = rollout_summary(everyone, t)
    assert list(df["cohort"].cat.categories) == ["1"]
    assert (df["status"] == "Treated").all()

    nobody = np.zeros((10, 5))
    df = rollout_summary(nobody, t)
    assert list(df["cohort"].cat.categories) == ["Never treated"]
    assert (df["status"] == "Never treated").all()
    assert df["n_units"].iloc[0] == 10


def test_uneven_calendar_time_is_carried_through():
    """Exposure is elapsed time, not a count of periods.

    With t = (2, 5, 9) and adoption at t = 5, the period before adoption is
    t = 2, so exposure is 3 at adoption and 7 two periods later -- not 1 and 2.
    Pinning this here keeps the rollout picture on the same axis as the ATT.
    """
    z = np.zeros((4, 3))
    z[:2, 1:] = 1.0
    df = rollout_summary(z, t=np.array([2.0, 5.0, 9.0]))
    assert sorted(df["period"].unique()) == [2.0, 5.0, 9.0]

    cohort = df[df["cohort"] == "5"].sort_values("period")
    np.testing.assert_array_equal(cohort["exposure"].to_numpy(), [0, 3, 7])
    np.testing.assert_array_equal(
        cohort["status"].astype(str).to_numpy(),
        ["Not yet treated", "Treated", "Treated"],
    )


def test_estimator_exposes_it_without_fitting():
    z, t, _ = _waves(n=40)
    direct = rollout_summary(z, t)
    via_class = LongBet.rollout_summary(z, t)
    pd.testing.assert_frame_equal(direct, via_class)


def test_single_period_panel():
    """A length-1 t must not index like a scalar (it arrives 0-d from R)."""
    z = np.array([[1.0], [0.0], [1.0]])
    df = rollout_summary(z, t=np.array(7.0))
    assert list(df["cohort"].cat.categories) == ["7", "Never treated"]
    assert int(df.loc[df["cohort"] == "7", "exposure"].iloc[0]) == 1
    assert int(df.loc[df["cohort"] == "7", "n_units"].iloc[0]) == 2


def test_bad_input_is_refused():
    with pytest.raises(ValueError, match=r"\(N, T\)"):
        rollout_summary(np.zeros(5))
    with pytest.raises(ValueError, match=r"\(N, T\)"):
        rollout_summary(np.zeros((2, 3, 4)))
    with pytest.raises(ValueError, match="empty"):
        rollout_summary(np.zeros((0, 3)))
    with pytest.raises(ValueError, match="length"):
        rollout_summary(np.zeros((3, 4)), t=np.arange(3))


def test_plot_draws_one_row_per_cohort():
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")

    z, t, _ = _waves(n=60)
    ax = plot_rollout(z, t, labels={4: "W1", 6: "W2", 8: "W3", 10: "W4"},
                      never_treated_label="Holdout")

    # Read the axis the way a viewer does: descending y is top to bottom.
    rows = sorted(zip(ax.get_yticks(), ax.get_yticklabels()), reverse=True)
    ordered = [lab.get_text() for _, lab in rows]
    assert len(ordered) == 5
    assert [name.split(" (")[0] for name in ordered] == [
        "W1", "W2", "W3", "W4", "Holdout"
    ], "cohorts should read top to bottom in adoption order, holdout last"
    assert all("(n=" in name for name in ordered)

    # One tile per cohort-period.
    assert len(ax.patches) == 5 * len(t)
