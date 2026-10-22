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

"""Describing a staggered rollout, before any model is fitted.

Every staggered-adoption analysis starts with the same question -- who was
treated, and when -- and answering it by hand means reconstructing the adoption
pattern that :func:`longbet.derive_exposure` already computes.  Getting that
reconstruction subtly wrong (an off-by-one in the exposure index, a cohort
boundary that does not match the one the model uses) is easy and quiet.

Two functions, deliberately separated:

* :func:`rollout_summary` returns a tidy table and pulls in nothing beyond
  pandas.  It is what you want for a plot in any library, a printed table, or a
  sanity check in a script.
* :func:`plot_rollout` draws the usual tile chart on top of it.  ``matplotlib``
  is imported only when it is called, so it stays an optional dependency.

Cohorts are **derived** from ``z`` rather than supplied: a cohort is the set of
units sharing a first-treated period.  When a design has named launch waves that
is exactly the waves, and it still works when it does not.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

#: Status of a cohort in a period.  Ordered as they should appear in a legend.
STATUSES = ("Not yet treated", "Treated", "Never treated")


def _validate(z: Any, t: Any) -> tuple[np.ndarray, np.ndarray]:
    """Coerce and check ``z`` and ``t``; shared so both entry points agree."""
    z = np.asarray(z)
    if z.ndim != 2:
        raise ValueError(f"z must be an (N, T) array, got shape {z.shape}")
    if z.size == 0:
        raise ValueError("z is empty")
    n_periods = z.shape[1]
    # atleast_1d: a single period arrives from R as a 0-d array, which indexes
    # like a scalar and would fail deep inside instead of here.
    t = (
        np.arange(1, n_periods + 1, dtype=np.float64)
        if t is None
        else np.atleast_1d(np.asarray(t, dtype=np.float64))
    )
    if t.size != n_periods:
        raise ValueError(f"t has length {t.size}, expected {n_periods}")
    return z, t


def adoption_cohorts(
    z: np.ndarray, t: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Assign each unit to an adoption cohort.

    Parameters
    ----------
    z
        Treatment indicator, shape ``(N, T)``. Must be absorbing.
    t
        Calendar time, length ``T``. Defaults to ``1..T``.

    Returns
    -------
    first_treated : np.ndarray
        The time at which each unit is first treated, ``NaN`` for units that
        never are. Shape ``(N,)``.
    cohort_times : np.ndarray
        The distinct adoption times, ascending. Never-treated units are not
        represented here; they form their own cohort in
        :func:`rollout_summary`.
    """
    z, t = _validate(z, t)

    treated = z == 1
    ever = treated.any(axis=1)
    first_idx = np.argmax(treated, axis=1)
    first_treated = np.where(ever, t[first_idx], np.nan)
    cohort_times = np.unique(first_treated[ever]) if ever.any() else np.array([])
    return first_treated, cohort_times


def rollout_summary(
    z: np.ndarray,
    t: np.ndarray | None = None,
    labels: dict[Any, str] | None = None,
    never_treated_label: str = "Never treated",
) -> pd.DataFrame:
    """Summarize a staggered rollout as one row per cohort and period.

    Parameters
    ----------
    z
        Treatment indicator, shape ``(N, T)``. Must be absorbing.
    t
        Calendar time, length ``T``. Defaults to ``1..T``.
    labels
        Optional display names keyed by adoption time, e.g.
        ``{5: "W1", 6: "W2"}``. Unlabelled cohorts fall back to their adoption
        time.
    never_treated_label
        Name for the cohort that is never treated.

    Returns
    -------
    pandas.DataFrame
        Columns:

        ``cohort``
            Display name of the adoption cohort.
        ``first_treated``
            Adoption time, ``NaN`` for the never-treated cohort.
        ``n_units``
            Units in the cohort.
        ``period`` / ``period_index``
            Calendar time and its 0-based column in ``z``.
        ``status``
            One of ``'Not yet treated'``, ``'Treated'``, ``'Never treated'``.
        ``exposure``
            The exposure index the model would use for this cohort and period,
            from :func:`longbet.derive_exposure`, so a plot of the rollout and
            the event-time axis of an ATT cannot disagree.

        Rows are ordered by adoption time, with the never-treated cohort last.

    Examples
    --------
    >>> import numpy as np
    >>> z = np.zeros((4, 3)); z[:2, 1:] = 1          # two units adopt at t = 2
    >>> rollout_summary(z).query("period == 2")[["cohort", "status", "n_units"]]
               cohort   status  n_units
    1             2.0  Treated        2
    4   Never treated  Never treated  2
    """
    from longbet._model import derive_exposure  # local: avoids a cycle

    z, t = _validate(z, t)
    first_treated, cohort_times = adoption_cohorts(z, t)
    exposure = derive_exposure(z, t)
    labels = labels or {}

    rows = []

    def _label(value: float) -> str:
        for key, name in labels.items():
            if np.isclose(float(key), value):
                return name
        return f"{value:g}"

    for adopt in cohort_times:
        members = np.isclose(first_treated, adopt)
        n_units = int(members.sum())
        # Every member of a cohort shares the same status and exposure path, so
        # one representative row per period describes the whole cohort.
        rep = int(np.argmax(members))
        for j, period in enumerate(t):
            rows.append(
                {
                    "cohort": _label(adopt),
                    "first_treated": float(adopt),
                    "n_units": n_units,
                    "period": float(period),
                    "period_index": j,
                    "status": "Treated" if z[rep, j] == 1 else "Not yet treated",
                    "exposure": int(exposure[rep, j]),
                }
            )

    n_never = int(np.isnan(first_treated).sum())
    if n_never:
        for j, period in enumerate(t):
            rows.append(
                {
                    "cohort": never_treated_label,
                    "first_treated": np.nan,
                    "n_units": n_never,
                    "period": float(period),
                    "period_index": j,
                    "status": "Never treated",
                    "exposure": 0,
                }
            )

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    order = [_label(a) for a in cohort_times] + (
        [never_treated_label] if n_never else []
    )
    frame["cohort"] = pd.Categorical(frame["cohort"], categories=order, ordered=True)
    frame["status"] = pd.Categorical(
        frame["status"], categories=list(STATUSES), ordered=True
    )
    return frame.sort_values(["cohort", "period_index"]).reset_index(drop=True)


#: Default fill colours, keyed by status.
ROLLOUT_COLORS = {
    "Not yet treated": "#d9d9d9",
    "Treated": "#3b6ea5",
    "Never treated": "#9e9e9e",
}


def plot_rollout(
    z: np.ndarray,
    t: np.ndarray | None = None,
    labels: dict[Any, str] | None = None,
    never_treated_label: str = "Never treated",
    ax: Any = None,
    colors: dict[str, str] | None = None,
    title: str | None = "Adoption cohorts over time",
    show_counts: bool = True,
) -> Any:
    """Draw a staggered rollout as a tile chart, one row per adoption cohort.

    ``matplotlib`` is imported here rather than at module load, so it remains an
    optional dependency: :func:`rollout_summary` needs nothing beyond pandas.

    Parameters
    ----------
    z, t, labels, never_treated_label
        As in :func:`rollout_summary`.
    ax
        Axes to draw on. A new figure is created when omitted.
    colors
        Fill colour per status; defaults to :data:`ROLLOUT_COLORS`.
    title
        Axes title, or ``None`` for none.
    show_counts
        Append the unit count to each cohort's tick label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "plot_rollout() needs matplotlib. Install it with "
            "'pip install matplotlib', or use rollout_summary() and plot the "
            "table with whatever you already have."
        ) from exc

    frame = rollout_summary(
        z, t, labels=labels, never_treated_label=never_treated_label
    )
    if frame.empty:
        raise ValueError("z contains no units")

    palette = {**ROLLOUT_COLORS, **(colors or {})}
    cohorts = list(frame["cohort"].cat.categories)
    periods = np.sort(frame["period"].unique())

    if ax is None:
        _, ax = plt.subplots(figsize=(max(6.0, 0.45 * len(periods)),
                                      max(2.0, 0.5 * len(cohorts))))

    # Cohorts read top-to-bottom in adoption order, so the staircase descends.
    y_of = {c: len(cohorts) - 1 - i for i, c in enumerate(cohorts)}
    width = float(np.min(np.diff(periods))) if periods.size > 1 else 1.0

    for row in frame.itertuples():
        ax.add_patch(
            plt.Rectangle(
                (row.period - width / 2, y_of[row.cohort] - 0.45),
                width, 0.9,
                facecolor=palette.get(str(row.status), "#cccccc"),
                edgecolor="white", linewidth=1.0,
            )
        )

    ticks, tick_labels = [], []
    for cohort in cohorts:
        ticks.append(y_of[cohort])
        n = int(frame.loc[frame["cohort"] == cohort, "n_units"].iloc[0])
        tick_labels.append(f"{cohort} (n={n})" if show_counts else str(cohort))

    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels)
    ax.set_ylim(-0.6, len(cohorts) - 0.4)
    ax.set_xlim(periods.min() - width / 2, periods.max() + width / 2)
    ax.set_xticks(periods)
    ax.set_xlabel("Period")
    if title:
        ax.set_title(title)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    present = [s for s in STATUSES if (frame["status"] == s).any()]
    ax.legend(
        handles=[Patch(facecolor=palette[s], edgecolor="white", label=s)
                 for s in present],
        loc="upper center", bbox_to_anchor=(0.5, -0.18),
        ncol=len(present), frameon=False,
    )
    return ax
