"""Design-based encouragement contrasts with explicit randomization and targets.

The independent observations are randomization units, which can be clusters.
All estimators are linear Horvitz--Thompson (HT) contrasts with fixed population
weights. In particular, subgroup estimates do not divide by random subgroup arm
sizes, and unequal-size clusters are not silently given equal participant weight.

For blocked complete assignment let ``F[b,k]`` be the vector of weighted outcome
totals in randomization unit k, multiplied by the number of randomization units
in block b. The estimator is the sum of block arm-mean differences in F; its
Neyman covariance is ``sum_b(S[b,1]/n[b,1] + S[b,0]/n[b,0])``. Relative to the
true randomization covariance, its expectation adds the PSD sum of finite
population treatment-effect covariance matrices divided by block sizes.

For independent Bernoulli/categorical assignment let ``X[k]`` be the observed
vector of signed inverse-probability-weighted totals. We use ``sum_k X[k]X[k]'``.
Independence gives ``Var(sum X) = sum E[XX'] - sum E[X]E[X]'``: the estimator is
a conservative PSD covariance bound in expectation. It retains all cross-time,
cross-group and reused-control terms. This deliberately simple bound can be
loose and is not invariant to outcome location shifts. Neither covariance bound
implies finite-sample normal coverage. See the HT variance discussion at
https://www150.statcan.gc.ca/n1/pub/12-001-x/2013001/article/11831/section2-eng.htm.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from ._encourage import _binary_panel, _critical, _wald_set
from ._model import _check_time_vector, derive_exposure


@dataclass(frozen=True)
class EncouragementDesign:
    """Declared assignment scheme and fixed-population analysis weights.

    ``assignment`` is ``'complete'`` (fixed arm counts within each block),
    ``'bernoulli'`` (independent assignments with known probabilities), or
    ``'staggered'`` (independent categorical encouragement-cohort assignment).
    The last does not cover fixed cohort quotas or sequential rerandomization.

    ``blocks`` and ``clusters`` are length-N baseline labels. Clusters must be
    wholly within blocks and share an encouragement path. ``target='unit'``
    gives equal participant weight; ``'cluster'`` gives equal weight to clusters,
    and then equal weight to the subgroup members within each included cluster.
    Without custom block weights these targets apply over the whole subgroup.
    ``block_weights`` can instead map every block label to a nonnegative fixed
    weight summing to one; every positive-weight block must contain the subgroup.

    For Bernoulli designs ``probabilities`` is a scalar, length-N vector constant
    within clusters, or a vector in first-appearance order of the clusters. For
    staggered designs it is a mapping from zero-based cohort start indices and
    optionally ``'never'`` to such vectors, summing to one for every independent
    unit. A category can have zero probability. Observed assignments must have
    positive probability. ``cohort_weights`` is required for staggered designs:
    map the finite cohorts of interest to fixed positive weights summing to one.
    ``control='not_yet'`` uses later cohorts and permanent holdouts;
    ``'never'`` uses only holdouts. No anticipation is assumed in either case.
    Cohort contrasts and fixed-weight horizon aggregates are returned. Aggregates
    are unavailable if any specified cohort lacks the horizon or design support;
    the weights are never renormalized after controls disappear.

    Label/probability arrays describe the original experiment. Supplying them
    cannot establish randomization, no interference, exclusion, or monotonicity.
    """

    assignment: str = "complete"
    blocks: Any = None
    clusters: Any = None
    target: str = "unit"
    probabilities: Any = None
    block_weights: Any = None
    cohort_weights: Any = None
    control: str = "not_yet"


@dataclass
class EncouragementDesignResult:
    """Portable contrast table, joint covariance, and serializable provenance.

    ``covariance`` has a MultiIndex on both axes: ``(contrast_id, quantity)``
    with quantity ``itt_y`` or ``itt_d``. Unavailable rows/columns are NaN.
    Off-diagonal entries retain reuse of the same participants and controls.
    Intervals in ``table`` are pointwise normal intervals/AR-style sets, not
    simultaneous bands or posterior intervals.
    """

    table: pd.DataFrame
    covariance: pd.DataFrame
    metadata: dict[str, Any]


def _labels(value: Any, n: int, name: str, default: Any) -> tuple[np.ndarray, list[str]]:
    a = np.asarray(default if value is None else value, dtype=object)
    if a.shape != (n,) or np.any(pd.isna(a)):
        raise ValueError(f"{name} must be a complete length-{n} vector of baseline labels.")
    if any(not isinstance(v, (str, int, float, bool, np.integer, np.floating,
                              np.bool_)) for v in a):
        raise ValueError(f"{name} labels must be strings or finite numeric scalars.")
    if any(isinstance(v, (float, np.floating)) and not np.isfinite(v) for v in a):
        raise ValueError(f"{name} labels must be finite.")
    # String labels make metadata and R conversion portable. Refuse collisions.
    originals: dict[str, tuple[type, Any]] = {}
    for v in a:
        key = str(v)
        signature = (type(v), v)
        if key in originals and originals[key] != signature:
            raise ValueError(f"{name} labels collide after string conversion: {key!r}.")
        originals[key] = signature
    strings = np.array([str(v) for v in a])
    return strings, list(dict.fromkeys(map(str, strings)))


def _fixed_weights(value: Any, labels: list[str], name: str) -> np.ndarray:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must map labels to fixed nonnegative weights summing to one.")
    mapping = {str(k): v for k, v in value.items()}
    if len(mapping) != len(value) or set(mapping) != set(labels):
        raise ValueError(f"{name} must contain exactly these labels: {labels}.")
    try:
        weights = np.array([mapping[k] for k in labels], dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain finite numeric weights.") from exc
    if (not np.all(np.isfinite(weights)) or np.any(weights < 0)
            or not np.isclose(weights.sum(), 1, rtol=0, atol=1e-10)):
        raise ValueError(f"{name} must be nonnegative and sum to one.")
    return weights


def _cluster_probability(value: Any, members: list[np.ndarray], n: int,
                         name: str) -> np.ndarray:
    try:
        a = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain finite probabilities.") from exc
    if a.ndim == 0:
        result = np.full(len(members), float(a))
    elif a.shape == (n,):
        if any(not np.all(a[idx] == a[idx[0]]) for idx in members):
            raise ValueError(f"{name} must be constant within each randomization cluster.")
        result = np.array([a[idx[0]] for idx in members])
    elif a.shape == (len(members),):
        result = a.copy()
    else:
        raise ValueError(f"{name} must be scalar, length N, or length number of clusters.")
    if not np.all(np.isfinite(result)) or np.any((result < 0) | (result > 1)):
        raise ValueError(f"{name} must contain finite probabilities between zero and one.")
    return result


def _cohort_key(key: Any, n_periods: int) -> int:
    if key == "never":
        return n_periods
    # Numeric strings permit lossless JSON / R named-list round trips.
    if isinstance(key, (bool, np.bool_)):
        raise ValueError("Cohort keys must be zero-based integer indices or 'never'.")
    try:
        integer = int(key)
        if str(integer) != str(key) or not 0 <= integer < n_periods:
            raise ValueError
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError("Cohort keys must be zero-based integer indices or 'never'.") from exc
    return integer


def _panel(y: Any, d: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    z, d = _binary_panel(z, "z"), _binary_panel(d, "d")
    y = np.asarray(y)
    if y.shape != z.shape or y.dtype.kind not in "biuf" or d.shape != z.shape:
        raise ValueError("y, d, and z must be complete numeric panels with the same (N, T) shape.")
    y = y.astype(float)
    if not np.all(np.isfinite(y)):
        raise ValueError("y must be complete and finite; missing cells are not dropped.")
    t = np.arange(1, z.shape[1] + 1, dtype=float) if t is None else np.asarray(t, dtype=float)
    if t.shape != (z.shape[1],):
        raise ValueError("t must be a one-dimensional vector with one entry per period.")
    _check_time_vector(t)
    if not np.allclose(np.diff(t), np.rint(np.diff(t)), rtol=0, atol=1e-8):
        raise ValueError("t must have whole-unit gaps for exposure indexing.")
    starts = np.where(z[:, -1] == 1, np.argmax(z == 1, axis=1), z.shape[1])
    return y, d, z, t, starts


def design_encouragement_effects(
    y: Any, d: Any, z: Any, t: Any = None, *,
    design: EncouragementDesign | None = None, groups: Any = None, alpha: float = 0.05,
) -> EncouragementDesignResult:
    """Estimate weighted encouragement ITTs and joint reference uncertainty.

    See :class:`EncouragementDesign` for assignment assumptions and weighting.
    ``groups`` is a complete length-N vector of prespecified baseline subgroup
    labels. The result includes ``'all'`` and every supplied group; ``'all'`` is
    reserved. Data-adaptive group selection is not accounted for. Small groups
    retain HT point estimates, with unavailable uncertainty when fewer than two
    independent units in either arm contribute to an active complete-assignment
    block (or to a Bernoulli/categorical contrast).

    Both ITTs always use identical fixed participant weights, periods, and arms.
    No observations are pooled as independent across time. The covariance uses
    whole randomization clusters. Known probability HT estimates, including
    binary ITTs, can fall outside their population parameter bounds. They are
    not clipped. Blocked complete overall unit assignment reproduces the usual
    block-weighted differences in means; subgroup HT estimates generally differ
    from differences of subgroup arm means.

    Staggered rows explicitly distinguish ``design_support_unavailable`` from
    ``insufficient_arm_size``. Neither situation triggers extrapolation. Only
    horizons observed for all positively weighted cohorts can be aggregated.
    Complete/randomized-cohort quotas and adaptive sequential assignment are
    not approximated as independent categorical assignment.

    Wald sets invert the joint normal contrast using the same procedure as
    :func:`longbet.encouragement_effects`. These are asymptotic, pointwise
    confidence sets. A Wald/CACE interpretation requires additional IV and
    adoption-history assumptions; no test in this function establishes them.
    """
    critical = _critical(alpha)
    design = EncouragementDesign() if design is None else design
    if not isinstance(design, EncouragementDesign):
        raise TypeError("design must be an EncouragementDesign instance.")
    if design.assignment not in {"complete", "bernoulli", "staggered"}:
        raise ValueError("assignment must be 'complete', 'bernoulli', or 'staggered'.")
    if design.target not in {"unit", "cluster"}:
        raise ValueError("target must be 'unit' or 'cluster'.")
    if design.control not in {"not_yet", "never"}:
        raise ValueError("control must be 'not_yet' or 'never'.")
    if design.assignment != "staggered" and design.cohort_weights is not None:
        raise ValueError("cohort_weights only applies to staggered assignment.")
    if design.assignment == "complete" and design.probabilities is not None:
        raise ValueError("Complete assignment uses fixed block arm counts, not probabilities.")
    y, d, z, t, starts = _panel(y, d, z, t)
    n, nt = z.shape
    cluster, cluster_labels = _labels(design.clusters, n, "clusters", np.arange(n))
    block, block_labels = _labels(design.blocks, n, "blocks", np.repeat("all", n))
    members = [np.flatnonzero(cluster == key) for key in cluster_labels]
    nr = len(members)
    if any(np.unique(block[idx]).size != 1 for idx in members):
        raise ValueError("Each randomization cluster must lie wholly within one block.")
    if any(np.unique(starts[idx]).size != 1 for idx in members):
        raise ValueError("All members of each randomization cluster must share the same z path.")
    rb = np.array([block[idx[0]] for idx in members])
    ra = np.array([starts[idx[0]] for idx in members])
    block_members = [np.flatnonzero(rb == key) for key in block_labels]
    group_labels, masks = ["all"], [np.ones(n, dtype=bool)]
    group_by_unit = None
    if groups is not None:
        gl, keys = _labels(groups, n, "groups", None)
        if "all" in keys:
            raise ValueError("The group label 'all' is reserved for the overall target.")
        group_by_unit = list(map(str, gl))
        group_labels += keys
        masks += [gl == key for key in keys]
    ng = len(masks)
    weights = np.zeros((ng, n))
    group_reasons = [""] * ng
    weights_by_block: dict[str, dict[str, float]] = {}
    fixed = (None if design.block_weights is None else
             _fixed_weights(design.block_weights, block_labels, "block_weights"))
    for g, mask in enumerate(masks):
        mass = np.array([
            np.count_nonzero(mask & (block == label)) if design.target == "unit" else
            sum(np.any(mask[members[k]]) for k in ids)
            for label, ids in zip(block_labels, block_members, strict=True)
        ], dtype=float)
        wb = mass / mass.sum() if fixed is None else fixed
        weights_by_block[group_labels[g]] = dict(zip(block_labels, map(float, wb), strict=True))
        if np.any((wb > 0) & (mass == 0)):
            group_reasons[g] = "empty_target_block"
            continue
        for b, (label, ids) in enumerate(zip(block_labels, block_members, strict=True)):
            if wb[b] == 0:
                continue
            if design.target == "unit":
                weights[g, mask & (block == label)] = wb[b] / mass[b]
            else:
                for k in ids:
                    selected = members[k][mask[members[k]]]
                    if selected.size:
                        weights[g, selected] = wb[b] / mass[b] / selected.size

    probabilities: dict[int, np.ndarray] = {}
    cohort_weights: dict[int, float] = {}
    if design.assignment == "staggered":
        if not isinstance(design.probabilities, Mapping) or not design.probabilities:
            raise ValueError("Staggered assignment requires a mapping of known cohort probabilities.")
        for key, value in design.probabilities.items():
            c = _cohort_key(key, nt)
            if c in probabilities:
                raise ValueError("Duplicate cohort probability keys after index conversion.")
            probabilities[c] = _cluster_probability(value, members, n, f"probabilities[{key!r}]")
        if not np.allclose(sum(probabilities.values()), 1, rtol=0, atol=1e-10):
            raise ValueError("Cohort probabilities must sum to one for every randomization unit.")
        if any(c not in probabilities or probabilities[c][k] <= 0 for k, c in enumerate(ra)):
            raise ValueError("Every observed cohort assignment must have positive known probability.")
        if not isinstance(design.cohort_weights, Mapping) or not design.cohort_weights:
            raise ValueError("Staggered assignment requires explicit fixed cohort_weights.")
        for key, value in design.cohort_weights.items():
            c = _cohort_key(key, nt)
            if c == nt or c not in probabilities or c in cohort_weights:
                raise ValueError("cohort_weights keys must name distinct finite cohorts in probabilities.")
            try:
                cohort_weights[c] = float(value)
            except (TypeError, ValueError) as exc:
                raise ValueError("cohort_weights must contain finite positive weights.") from exc
        if (not all(np.isfinite(v) and v > 0 for v in cohort_weights.values())
                or not np.isclose(sum(cohort_weights.values()), 1, rtol=0, atol=1e-10)):
            raise ValueError("cohort_weights must be positive and sum to one.")
        cohorts = sorted(cohort_weights)
    else:
        cohorts = sorted(set(starts[starts < nt]))
        if len(cohorts) != 1:
            raise ValueError("Single-wave designs need exactly one observed start; use staggered for multiple waves.")
        if design.assignment == "bernoulli":
            if design.probabilities is None or isinstance(design.probabilities, Mapping):
                raise ValueError("Bernoulli assignment requires known scalar or vector probabilities.")
            p = _cluster_probability(design.probabilities, members, n, "probabilities")
            if np.any((p <= 0) | (p >= 1)):
                raise ValueError("Bernoulli probabilities must be strictly between zero and one.")
            probabilities = {cohorts[0]: p, nt: 1 - p}

    # Fixed weighted outcome totals for each independent randomization unit.
    f = np.zeros((nr, ng, nt, 2))
    paired = np.stack((y, d), axis=-1)
    for k, idx in enumerate(members):
        f[k] = np.einsum("gi,itq->gtq", weights[:, idx], paired[idx])
    active = np.array([[np.any(weights[g, idx] > 0) for g in range(ng)] for idx in members])
    rows: list[dict[str, Any]] = []
    contributions: list[np.ndarray] = []
    values: list[np.ndarray] = []
    point_reasons: list[str] = []
    variance_reasons: list[str] = []
    complete_vectors: list[np.ndarray] = []
    row_lookup: dict[tuple[int, int, int], int] = {}
    horizons: dict[int, dict[int, int]] = {}
    for c in cohorts:
        counterfactual = np.zeros((1, nt))
        counterfactual[:, c:] = 1
        exposure = derive_exposure(counterfactual, t)[0]
        horizons[c] = {int(exposure[j]): j for j in range(c, nt)}
        for g in range(ng):
            for j in range(c, nt):
                treated = ra == c
                eligible = (ra > j) if design.control == "not_yet" else (ra == nt)
                if design.assignment != "staggered":
                    eligible = ra == nt
                reason = group_reasons[g]
                variance_reason = ""
                x = np.zeros((nr, 2))
                point = np.zeros(2)
                complete_vector = np.zeros((nr, 2))
                if design.assignment == "complete":
                    for ids in block_members:
                        if not np.any(active[ids, g]):
                            continue
                        a, b = ids[treated[ids]], ids[eligible[ids]]
                        if min(a.size, b.size) == 0:
                            reason = reason or "design_support_unavailable"
                            continue
                        complete_vector[ids] = len(ids) * f[ids, g, j]
                        x[a] = complete_vector[a] / len(a)
                        x[b] = -complete_vector[b] / len(b)
                        # Preserve exact cancellation when binary arm rates are
                        # identical; summing interleaved signed contributions
                        # can create a spurious 1e-17 first stage and huge ratio.
                        point += complete_vector[a].mean(axis=0) - complete_vector[b].mean(axis=0)
                        if min(np.count_nonzero(active[a, g]), np.count_nonzero(active[b, g])) < 2:
                            variance_reason = "insufficient_arm_size"
                else:
                    pc = probabilities[c]
                    if design.assignment == "bernoulli" or design.control == "never":
                        pe = probabilities.get(nt, np.zeros(nr))
                    else:
                        pe = sum((p for k, p in probabilities.items() if k > j), np.zeros(nr))
                    if np.any(active[:, g] & ((pc <= 0) | (pe <= 0))):
                        reason = reason or "design_support_unavailable"
                    a = treated & active[:, g] & (pc > 0)
                    b = eligible & active[:, g] & (pe > 0)
                    x[a] = f[a, g, j] / pc[a, None]
                    x[b] = -f[b, g, j] / pe[b, None]
                    if min(np.count_nonzero(a), np.count_nonzero(b)) < 2:
                        variance_reason = "insufficient_arm_size"
                target_members = weights[g] > 0
                row = {
                    "contrast_id": f"c{len(rows)}", "group": group_labels[g],
                    "contrast_type": "cohort" if design.assignment == "staggered" else "single_wave",
                    "cohort_index": c, "cohort_period": float(t[c]),
                    "period": float(t[j]), "period_index": j, "horizon": int(exposure[j]),
                    "post_encouragement": True,
                    "n_units_target": int(target_members.sum()),
                    "n_randomization_units": int(active[:, g].sum()),
                    "n_encouraged": int(np.count_nonzero(target_members & (starts == c))),
                    "n_control": int(sum(np.count_nonzero(target_members[members[k]])
                                         for k in np.flatnonzero(eligible))),
                    "n_randomized_encouraged": int(np.count_nonzero(active[:, g] & treated)),
                    "n_randomized_control": int(np.count_nonzero(active[:, g] & eligible)),
                    "arm_count_basis": "contrast",
                    "cohort_weight": float(cohort_weights.get(c, 1)),
                }
                row_lookup[g, c, int(exposure[j])] = len(rows)
                rows.append(row)
                contributions.append(x)
                complete_vectors.append(complete_vector)
                values.append(point if design.assignment == "complete" else x.sum(axis=0))
                point_reasons.append(reason)
                variance_reasons.append(reason or variance_reason)

    aggregate_components: dict[str, list[dict[str, Any]]] = {}
    if design.assignment == "staggered":
        for g in range(ng):
            for h in sorted(set().union(*(set(x) for x in horizons.values()))):
                available = [(row_lookup[g, c, h], cohort_weights[c])
                             for c in cohorts if h in horizons[c]]
                missing = len(available) != len(cohorts)
                reason = "missing_cohort_horizon" if missing else ""
                variance_reason = ""
                x = np.zeros((nr, 2))
                for index, weight in available:
                    x += weight * contributions[index]
                    reason = reason or point_reasons[index]
                    variance_reason = variance_reason or variance_reasons[index]
                template = rows[available[0][0]]
                row = dict(template, contrast_id=f"c{len(rows)}", contrast_type="horizon_average",
                           cohort_index=None, cohort_period=np.nan, period=np.nan,
                           period_index=None, cohort_weight=1.0,
                           arm_count_basis="minimum_across_components")
                for name in ("n_encouraged", "n_control", "n_randomized_encouraged", "n_randomized_control"):
                    row[name] = min(rows[index][name] for index, _ in available)
                aggregate_components[row["contrast_id"]] = [
                    {"contrast_id": rows[index]["contrast_id"], "weight": weight}
                    for index, weight in available
                ]
                rows.append(row)
                contributions.append(x)
                values.append(x.sum(axis=0))
                point_reasons.append(reason)
                variance_reasons.append(reason or variance_reason)

    q = len(rows) * 2
    if design.assignment == "complete":
        vector = np.stack(complete_vectors, axis=1).reshape(nr, q)
        covariance = np.zeros((q, q))
        for ids in block_members:
            for arm in (ra[ids] < nt, ra[ids] == nt):
                selected = ids[arm]
                if selected.size >= 2:
                    covariance += np.atleast_2d(np.cov(vector[selected], rowvar=False, ddof=1)) / selected.size
    else:
        vector = np.stack(contributions, axis=1).reshape(nr, q)
        covariance = vector.T @ vector
    if not np.all(np.isfinite(covariance)) or not np.all(np.isfinite(values)):
        raise ValueError("Inputs produce nonfinite moments; rescale outcomes or avoid extreme probabilities.")
    for r, row in enumerate(rows):
        dy, dd = (np.full(2, np.nan) if point_reasons[r] else values[r])
        if variance_reasons[r]:
            covariance[2*r:2*r+2, :] = np.nan
            covariance[:, 2*r:2*r+2] = np.nan
        vy, vd = np.diag(covariance)[2*r:2*r+2]
        cyd = covariance[2*r, 2*r+1]
        row.update(itt_y_d_cov=float(cyd), wald=float(dy / dd) if dd != 0 else np.nan,
                   inference="normal_ar", estimator="horvitz_thompson",
                   design_support=not bool(point_reasons[r]),
                   inference_available=not bool(variance_reasons[r]),
                   inference_reason=variance_reasons[r],
                   negative_first_stage=bool(dd < 0) if np.isfinite(dd) else None)
        for name, delta, variance in (("itt_y", dy, vy), ("itt_d", dd, vd)):
            se = np.sqrt(max(0.0, variance)) if np.isfinite(variance) else np.nan
            row.update({name: float(delta), f"{name}_se": float(se),
                        f"{name}_lower": float(delta - critical * se),
                        f"{name}_upper": float(delta + critical * se)})
        row["first_stage_includes_zero"] = (
            bool(row["itt_d_lower"] <= 0 <= row["itt_d_upper"])
            if np.isfinite(vd) else None
        )
        row.update(_wald_set(dy, dd, vy, vd, cyd, critical))
        if variance_reasons[r]:
            row["wald_reason"] = variance_reasons[r]
    index = pd.MultiIndex.from_product(
        [[row["contrast_id"] for row in rows], ["itt_y", "itt_d"]],
        names=["contrast_id", "quantity"],
    )
    metadata = {
        "version": 1, "assignment": design.assignment, "target": design.target,
        "estimator": "horvitz_thompson", "n_units": n,
        "n_randomization_units": nr, "cluster_labels": cluster_labels,
        "cluster_by_unit": list(map(str, cluster)),
        "block_by_unit": list(map(str, block)),
        "group_by_unit": group_by_unit,
        "start_indices": list(map(int, starts)),
        "never_start_index": nt,
        "time": list(map(float, t)),
        "block_labels": block_labels, "group_labels": group_labels,
        "block_arm_counts": {
            label: {"n_randomization_units": int(len(ids)),
                    "n_encouraged": int(np.count_nonzero(ra[ids] < nt)),
                    "n_never": int(np.count_nonzero(ra[ids] == nt)),
                    "cohort_counts": {str(c): int(np.count_nonzero(ra[ids] == c))
                                      for c in sorted(set(map(int, ra)))}}
            for label, ids in zip(block_labels, block_members, strict=True)
        },
        "block_weights": weights_by_block,
        "unit_weights": {label: list(map(float, weights[g])) for g, label in enumerate(group_labels)},
        "cohort_weights": {str(c): w for c, w in cohort_weights.items()},
        "cohort_key_scale": "zero_based_observed_period_index",
        "probabilities": {("never" if c == nt else str(c)): list(map(float, p))
                          for c, p in probabilities.items()},
        "aggregate_components": aggregate_components,
        "control": design.control if design.assignment == "staggered" else "never",
        "covariance": ("blockwise_neyman" if design.assignment == "complete" else
                       "independent_assignment_uncentered_psd_bound"),
        "confidence": "pointwise_asymptotic_normal_ar",
        "alpha": float(alpha),
        "assumptions": [
            "Assignment follows the declared scheme at the randomization-unit level.",
            "No interference between randomization units; consistency and complete panels.",
            "Groups, blocks, target weights, and probabilities are fixed before inspecting effects.",
            "Wald interpretation additionally requires exclusion, monotonicity, relevance, and an appropriate adoption-history target.",
        ] + (["Cohort assignments are independent across randomization units, not fixed quotas.",
              "No anticipation: eligible later cohorts share the untreated potential outcome at the comparison period."]
             if design.assignment == "staggered" else []),
        "limitations": [
            "Normal critical values need sufficient independent units and suitable moments; covariance conservativeness does not imply finite-sample coverage.",
            "Intervals are pointwise, not adjusted for subgroup or horizon multiplicity.",
            "Fixed-denominator subgroup and unequal-size cluster HT contrasts need not be invariant to outcome location shifts.",
        ] + (["The uncentered covariance bound can be loose and depends on outcome location."]
             if design.assignment != "complete" else []),
    }
    return EncouragementDesignResult(pd.DataFrame(rows), pd.DataFrame(covariance, index=index, columns=index), metadata)
