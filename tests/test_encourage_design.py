"""Independent checks of generalized randomization estimands and covariance."""

import itertools
import json

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from longbet._encourage import encouragement_effects
from longbet._encourage_design import EncouragementDesign, design_encouragement_effects


def panel(assignment, periods=3, start=1):
    a = np.asarray(assignment, dtype=bool)
    z = np.zeros((a.size, periods))
    z[a, start:] = 1
    return z


def test_default_reproduces_single_wave_and_cross_horizon_covariance():
    rng = np.random.default_rng(71)
    z = panel([1] * 6 + [0] * 6)
    y = rng.normal(size=z.shape)
    d = panel(rng.binomial(1, .5, len(z)))
    result = design_encouragement_effects(y, d, z)
    reference = encouragement_effects(y, d, z)
    pd.testing.assert_frame_equal(result.table[reference.columns], reference)
    values = np.stack((y[:, 1:], d[:, 1:]), axis=-1).reshape(len(z), -1)
    expected = np.cov(values[:6], rowvar=False) / 6 + np.cov(values[6:], rowvar=False) / 6
    assert_allclose(result.covariance, expected, atol=1e-15)
    assert result.covariance.index.names == ["contrast_id", "quantity"]
    json.dumps(result.metadata, allow_nan=False)


def test_block_weights_are_fixed_and_joint_covariance_is_blockwise():
    rng = np.random.default_rng(9)
    a = np.tile([1, 1, 1, 0, 0, 0], 2)
    z = panel(a)
    y = rng.normal(size=z.shape)
    d = panel([1, 1, 0, 1, 0, 0] * 2)
    blocks = ["east"] * 6 + ["west"] * 6
    result = design_encouragement_effects(
        y, d, z, design=EncouragementDesign(blocks=blocks, block_weights={"east": .8, "west": .2}),
    )
    refs = [encouragement_effects(y[s], d[s], z[s]) for s in (slice(0, 6), slice(6, 12))]
    assert_allclose(result.table.itt_y, .8 * refs[0].itt_y + .2 * refs[1].itt_y)
    assert_allclose(result.table.itt_d, .8 * refs[0].itt_d + .2 * refs[1].itt_d)
    assert_allclose(result.table.itt_y_se**2, .8**2 * refs[0].itt_y_se**2 + .2**2 * refs[1].itt_y_se**2)
    assert_allclose(result.table.itt_y_d_cov, .8**2 * refs[0].itt_y_d_cov + .2**2 * refs[1].itt_y_d_cov)


def test_subgroups_use_fixed_population_weights_and_retain_overlap_covariance():
    a = np.array([1] * 6 + [0] * 6)
    z = panel(a)
    y = np.arange(36).reshape(12, 3)
    groups = np.array(["left"] * 4 + ["right"] * 2 + ["left"] * 2 + ["right"] * 4)
    result = design_encouragement_effects(y, z, z, groups=groups)
    left = result.table.query("group == 'left'")
    target = groups == "left"
    # Allocation is 1/2 in the original randomization block. The subgroup's
    # realized 4:2 allocation is not substituted for that known fraction.
    expected = 2 * (y[target & (a == 1), 1:].sum(axis=0) - y[target & (a == 0), 1:].sum(axis=0)) / target.sum()
    assert_allclose(left.itt_y, expected)
    assert not np.allclose(expected, y[target & (a == 1), 1:].mean(axis=0) - y[target & (a == 0), 1:].mean(axis=0))
    assert_allclose(result.metadata["unit_weights"]["left"], target / target.sum())
    assert np.any(result.covariance.loc[("c0", "itt_y"), left.contrast_id].to_numpy() != 0)


def test_tiny_and_unrepresented_groups_preserve_points_without_inference():
    z = panel([1] * 4 + [0] * 4)
    groups = ["tiny", "rare", "main", "main", "rare", "main", "main", "main"]
    result = design_encouragement_effects(np.arange(24).reshape(8, 3), z, z, groups=groups)
    tiny = result.table.query("group == 'tiny'")
    assert tiny.itt_y.notna().all()
    assert (tiny.n_control == 0).all()
    assert (tiny.wald_reason == "insufficient_arm_size").all()
    assert tiny.itt_y_se.isna().all()
    assert result.covariance.loc[tiny.contrast_id.tolist()].isna().all().all()
    assert result.table.query("group == 'all'").inference_available.all()


def test_empty_positive_weight_block_is_an_undefined_subgroup_target():
    z = panel([1, 1, 0, 0] * 2)
    result = design_encouragement_effects(
        z, z, z, groups=["left"] * 4 + ["right"] * 4,
        design=EncouragementDesign(blocks=["A"] * 4 + ["B"] * 4, block_weights={"A": .5, "B": .5}),
    )
    subgroup = result.table.query("group != 'all'")
    assert (subgroup.inference_reason == "empty_target_block").all()
    assert not subgroup.design_support.any()
    assert subgroup.itt_y.isna().all()


@pytest.mark.parametrize("target", ["unit", "cluster"])
def test_cluster_targets_exact_randomization_unbiasedness_and_neyman_bound(target):
    sizes = np.array([1, 2, 3, 2, 1, 3])
    clusters = np.repeat(np.arange(6), sizes)
    base = np.arange(12, dtype=float) / 7
    effect = np.repeat([-.7, .1, 1.8, 2.4, .3, -.1], sizes)
    # Two outcomes/horizons with heterogeneous effects, not identical vectors.
    y0 = np.column_stack((base, 1.2 * base, -.3 * base))
    y1 = y0 + np.column_stack((np.zeros(12), effect, effect * .8))
    d0 = panel(np.repeat([0, 1, 0, 0, 1, 0], sizes))
    d1 = panel(np.repeat([1, 1, 1, 0, 1, 0], sizes))
    weights = np.full(12, 1 / 12) if target == "unit" else 1 / (6 * sizes[clusters])
    truth = np.stack((weights @ (y1 - y0)[:, 1:], weights @ (d1 - d0)[:, 1:]), axis=-1).ravel()
    estimates, covariance = [], []
    for treated in itertools.combinations(range(6), 3):
        a = np.isin(clusters, treated)
        z = panel(a)
        result = design_encouragement_effects(
            np.where(a[:, None], y1, y0), np.where(a[:, None], d1, d0), z,
            design=EncouragementDesign(clusters=clusters, target=target),
        )
        estimates.append(result.table[["itt_y", "itt_d"]].to_numpy().ravel())
        covariance.append(result.covariance.to_numpy())
        assert (result.table.n_randomization_units == 6).all()
        assert (result.table.n_randomized_encouraged == 3).all()
    estimates = np.array(estimates)
    assert_allclose(estimates.mean(axis=0), truth, atol=1e-15)
    exact = np.cov(estimates, rowvar=False, ddof=0)
    excess = np.mean(covariance, axis=0) - exact
    assert np.linalg.eigvalsh(excess).min() > -1e-12
    # The omitted finite-population effect covariance supplies the exact gap.
    cluster_effects = np.stack([
        np.stack((weights[clusters == k] @ (y1 - y0)[clusters == k, 1:],
                  weights[clusters == k] @ (d1 - d0)[clusters == k, 1:]), axis=-1).ravel() * 6
        for k in range(6)
    ])
    assert_allclose(excess, np.cov(cluster_effects, rowvar=False) / 6, atol=1e-14)


def test_clusters_are_not_treated_as_independent_participants():
    cluster_assignment = [1, 1, 1, 0, 0, 0]
    z = panel(np.repeat(cluster_assignment, 10))
    y = np.repeat(np.arange(6), 10)[:, None] * np.ones((1, 3))
    clustered = design_encouragement_effects(y, z, z, design=EncouragementDesign(clusters=np.repeat(np.arange(6), 10)))
    naive = encouragement_effects(y, z, z)
    assert (clustered.table.itt_y_se > naive.itt_y_se * 2).all()
    assert_allclose(clustered.table.itt_y, naive.itt_y)


def test_cluster_subgroups_can_split_clusters_with_declared_weights():
    clusters = np.repeat(np.arange(6), [2, 3, 4, 3, 2, 4])
    groups = np.where(np.arange(len(clusters)) % 2, "odd", "even")
    z = panel(clusters < 3)
    result = design_encouragement_effects(z, z, z, groups=groups,
                                        design=EncouragementDesign(clusters=clusters, target="cluster"))
    weights = np.asarray(result.metadata["unit_weights"]["odd"])
    for k in range(6):
        assert_allclose(weights[clusters == k].sum(), 1 / 6)
    assert_allclose(weights[groups == "even"], 0)


def test_known_unequal_probabilities_and_uncentered_joint_bound():
    a = np.array([1, 1, 1, 0, 0, 0])
    z = panel(a)
    p = np.array([.3, .5, .7, .2, .4, .6])
    y = np.arange(18).reshape(6, 3) / 3
    d = panel([1, 0, 1, 0, 1, 0])
    result = design_encouragement_effects(y, d, z, design=EncouragementDesign(assignment="bernoulli", probabilities=p))
    observed = np.stack((y[:, 1:], d[:, 1:]), axis=-1).reshape(6, 4)
    contribution = observed * np.where(a, 1 / p, -1 / (1 - p))[:, None] / 6
    assert_allclose(result.table[["itt_y", "itt_d"]].to_numpy().ravel(), contribution.sum(axis=0))
    assert_allclose(result.covariance, contribution.T @ contribution)
    assert result.metadata["covariance"] == "independent_assignment_uncentered_psd_bound"


def test_bernoulli_covariance_bound_matches_exact_finite_population_identity():
    # Algebra over all categories for each independent unit, without Monte Carlo.
    p = np.array([.3, .5, .7, .2, .4, .6])
    f0 = np.arange(24).reshape(6, 4) / 20
    f1 = f0 + np.array([-.3, .8, 1.2, -.7])
    expected_bound = sum(np.outer(f1[k], f1[k]) / p[k] + np.outer(f0[k], f0[k]) / (1 - p[k]) for k in range(6))
    variance = sum(p[k] * (1 - p[k]) * np.outer(f1[k] / p[k] + f0[k] / (1 - p[k]),
                                               f1[k] / p[k] + f0[k] / (1 - p[k])) for k in range(6))
    omitted = (f1 - f0).T @ (f1 - f0)
    assert_allclose(expected_bound - variance, omitted, atol=1e-14)
    assert np.linalg.eigvalsh(omitted).min() > -1e-12


def test_cluster_probabilities_accept_unit_and_cluster_representations():
    clusters = np.repeat(np.arange(6), 2)
    z = panel(clusters < 3)
    p = [.3, .4, .5, .6, .7, .8]
    first = design_encouragement_effects(z, z, z, design=EncouragementDesign(assignment="bernoulli", clusters=clusters, probabilities=p))
    second = design_encouragement_effects(z, z, z, design=EncouragementDesign(assignment="bernoulli", clusters=clusters, probabilities=np.repeat(p, 2)))
    pd.testing.assert_frame_equal(first.table, second.table)
    pd.testing.assert_frame_equal(first.covariance, second.covariance)


def staggered_panel():
    starts = np.repeat([0, 1, 2, 3], 3)
    z = (np.arange(3)[None, :] >= starts[:, None]).astype(float)
    y = np.arange(36).reshape(12, 3) / 10
    return y, z, starts


@pytest.mark.parametrize("control", ["never", "not_yet"])
def test_staggered_eligible_controls_reuse_covariance_and_fixed_aggregation(control):
    y, z, starts = staggered_panel()
    result = design_encouragement_effects(y, z, z, design=EncouragementDesign(
        assignment="staggered", probabilities={"0": .25, "1": .25, "2": .25, "never": .25},
        cohort_weights={"0": .4, "1": .6}, control=control,
    ))
    rows = result.table.query("contrast_type == 'cohort'")
    contributions = []
    for row in rows.itertuples():
        c, j = int(row.cohort_index), int(row.period_index)
        eligible = (starts > j) if control == "not_yet" else (starts == 3)
        pe = (3 - j) / 4 if control == "not_yet" else .25
        coefficient = np.where(starts == c, 4, np.where(eligible, -1 / pe, 0)) / len(z)
        contributions.append(np.column_stack((y[:, j], z[:, j])) * coefficient[:, None])
    vectors = np.stack(contributions, axis=1).reshape(len(z), -1)
    assert_allclose(rows[["itt_y", "itt_d"]].to_numpy().ravel(), vectors.sum(axis=0))
    assert_allclose(result.covariance.iloc[:len(rows)*2, :len(rows)*2], vectors.T @ vectors)
    # c0/t1 and c1/t1 share controls: their covariance is retained and positive.
    first = rows.query("cohort_index == 0 and period_index == 1").iloc[0].contrast_id
    second = rows.query("cohort_index == 1 and period_index == 1").iloc[0].contrast_id
    assert result.covariance.loc[(first, "itt_y"), (second, "itt_y")] > 0
    aggregates = result.table.query("contrast_type == 'horizon_average'")
    for row in aggregates.itertuples():
        components = result.metadata["aggregate_components"][row.contrast_id]
        if row.horizon == 3:
            assert row.inference_reason == "missing_cohort_horizon"
            assert np.isnan(row.itt_y)
            continue
        keys = [(x["contrast_id"], "itt_y") for x in components]
        w = np.array([x["weight"] for x in components])
        expected = sum(x["weight"] * result.table.set_index("contrast_id").loc[x["contrast_id"], "itt_y"] for x in components)
        assert_allclose(row.itt_y, expected)
        assert_allclose(row.itt_y_se**2, w @ result.covariance.loc[keys, keys].to_numpy() @ w)


def test_staggered_controls_disappear_without_extrapolation_or_weight_renormalization():
    starts = np.repeat([0, 1, 2], 3)
    z = (np.arange(3)[None, :] >= starts[:, None]).astype(float)
    result = design_encouragement_effects(z, z, z, design=EncouragementDesign(
        assignment="staggered", probabilities={0: 1/3, 1: 1/3, 2: 1/3},
        cohort_weights={0: .5, 1: .5},
    ))
    late = result.table.query("contrast_type == 'cohort' and period_index == 2")
    assert (late.inference_reason == "design_support_unavailable").all()
    assert late.itt_y.isna().all()
    assert result.table.query("contrast_type == 'horizon_average' and horizon == 2").itt_y.isna().all()
    assert result.metadata["cohort_weights"] == {"0": .5, "1": .5}


def test_staggered_heterogeneous_probability_lack_of_common_support_is_retained():
    y, z, _ = staggered_panel()
    p0 = np.full(12, .25); p0[3] = 0
    p1 = np.full(12, .25); p1[3] = .5
    result = design_encouragement_effects(y, z, z, design=EncouragementDesign(
        assignment="staggered", probabilities={0: p0, 1: p1, 2: .25, "never": .25}, cohort_weights={0: .5, 1: .5},
    ))
    assert not result.table.query("cohort_index == 0").design_support.any()
    assert result.table.query("cohort_index == 1").design_support.all()


def test_uneven_times_use_existing_exposure_index_and_no_cohort_horizon_filling():
    y, z, _ = staggered_panel()
    result = design_encouragement_effects(y, z, z, t=[1, 3, 6], design=EncouragementDesign(
        assignment="staggered", probabilities={0: .25, 1: .25, 2: .25, "never": .25}, cohort_weights={0: .5, 1: .5},
    ))
    # Uses derive_exposure's preceding-period clock, not observed-column count.
    assert result.table.query("cohort_index == 1").horizon.tolist() == [2, 5]
    unavailable = result.table.query("contrast_type == 'horizon_average'")
    assert (unavailable.inference_reason == "missing_cohort_horizon").any()


@pytest.mark.parametrize("kwargs,message", [
    ({"assignment": "adaptive"}, "assignment"),
    ({"target": "treated"}, "target"),
    ({"control": "all"}, "control"),
    ({"probabilities": .5}, "Complete assignment"),
    ({"assignment": "bernoulli"}, "known"),
    ({"assignment": "bernoulli", "probabilities": 1}, "strictly"),
    ({"assignment": "bernoulli", "probabilities": [1, 2]}, "length"),
    ({"cohort_weights": {1: 1}}, "only applies"),
    ({"block_weights": {"all": .3}}, "sum to one"),
    ({"block_weights": {"wrong": 1}}, "exactly"),
    ({"clusters": [0, 0, 0, 0, 0, 1, 1, 1]}, "same z path"),
    ({"blocks": [0, 1, 0, 0, 0, 0, 0, 0], "clusters": [0, 0, 1, 1, 2, 2, 3, 3]}, "wholly"),
    ({"assignment": "staggered", "probabilities": {1: .5, "never": .5}}, "cohort_weights"),
    ({"assignment": "staggered", "probabilities": {1: .5, "never": .2}, "cohort_weights": {1: 1}}, "sum to one"),
    ({"assignment": "staggered", "probabilities": {1: 0, "never": 1}, "cohort_weights": {1: 1}}, "observed cohort"),
    ({"assignment": "staggered", "probabilities": {1: .5, "never": .5}, "cohort_weights": {"never": 1}}, "finite cohorts"),
])
def test_design_validation(kwargs, message):
    z = panel([1] * 4 + [0] * 4)
    with pytest.raises(ValueError, match=message):
        design_encouragement_effects(z, z, z, design=EncouragementDesign(**kwargs))


def test_no_missing_labels_or_post_assignment_cluster_splitting():
    z = panel([1] * 4 + [0] * 4)
    with pytest.raises(ValueError, match="groups"):
        design_encouragement_effects(z, z, z, groups=["g"] * 7 + [None])
    with pytest.raises(ValueError, match="reserved"):
        design_encouragement_effects(z, z, z, groups=["all"] * 8)
    with pytest.raises(ValueError, match="constant within"):
        design_encouragement_effects(z, z, z, design=EncouragementDesign(
            assignment="bernoulli", clusters=np.repeat(np.arange(4), 2), probabilities=np.linspace(.1, .9, 8),
        ))


def test_block_with_no_randomized_arm_has_no_design_support():
    z = panel([1] * 4 + [0] * 4)
    result = design_encouragement_effects(z, z, z, design=EncouragementDesign(blocks=["a"] * 4 + ["b"] * 4))
    assert not result.table.design_support.any()
    assert result.table.itt_y.isna().all()
    assert (result.table.wald_reason == "design_support_unavailable").all()
