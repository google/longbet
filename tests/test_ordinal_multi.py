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

"""Heterogeneous ordinal counts, SUR conditionals, sharing, and joint archives."""
import dataclasses
import json

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from scipy import stats

from longbet import LongBet, LongBetConfig, LongBetMulti, effect_draws, joint_prob
from longbet._multi_input import normalize_multi_inputs
from longbet._multi_state import init_multi_longbet
from longbet._multi_step import multi_single_step
from longbet._step import longbet_single_step
from longbet._sur import conditional_residual
from test_full_sur import _constant_state
from test_ordinal_model import assert_invariants, ordinal_config, ordinal_panel


def mixed_panel():
    panel = ordinal_panel(5)
    rng = np.random.default_rng(837)
    a = np.arange(48).reshape(12, 4)
    ys = dict(continuous=rng.normal(size=(12, 4)), ordinal3=(a % 3).astype(float),
              binary=(a % 2).astype(float), ordinal5=(a % 5).astype(float))
    ys["ordinal3"][0, 0] = np.nan
    ys["ordinal5"][1, 1] = np.nan  # Non-nested discrete missingness is valid.
    ys["continuous"][0, 0] = ys["continuous"][1, 1] = np.nan
    panel["y"] = ys
    panel["outcome"] = ("continuous", "ordinal", "binary", "ordinal")
    panel["num_categories"] = (None, 3, None, 5)
    return panel


def mixed_config(**kwargs):
    cfg = dataclasses.replace(ordinal_config(), outcome="continuous", num_categories=None,
                             sigma_prior_a=2., sigma_prior_b=1., num_chains=2)
    return dataclasses.replace(cfg, **kwargs)


def normalize(panel=None, config=None, **kwargs):
    p = mixed_panel() if panel is None else panel
    options = dict(outcome=p["outcome"], outcome_names=None,
                   num_categories=p["num_categories"], config=config or mixed_config())
    options.update(kwargs)
    return normalize_multi_inputs(p["y"], **options)


def test_category_metadata_and_discrete_stable_partition():
    norm = normalize()
    assert norm.order == (1, 2, 3, 0) and norm.inverse_order == (3, 0, 1, 2)
    assert norm.user_num_categories == (None, 3, None, 5)
    assert norm.internal_num_categories == (3, None, 5, None)
    assert norm.internal_outcomes == ("ordinal", "binary", "ordinal", "continuous")
    for m in (0, 1, 2):
        assert norm.meany[m] == 0 and norm.sdy[m] == 1
    assert normalize(num_categories=5).user_num_categories == (None, 5, None, 5)
    mapping = dict(continuous=None, ordinal3=3, binary=None, ordinal5=5)
    assert normalize(num_categories=mapping).user_num_categories == norm.user_num_categories
    cfg = dataclasses.replace(mixed_config(), outcome="ordinal", num_categories=5)
    assert normalize(config=cfg, num_categories=None).user_num_categories == (None, 5, None, 5)


@pytest.mark.parametrize("counts", [None, True, 3., "3", [None, 3, None],
    [None, 3.1, None, 5], [2, 3, None, 5], [None, 3, 2, 5], [None, 3, None, None],
    dict(ordinal3=3, ordinal5=5), dict(continuous=None, ordinal3=3, binary=None, wrong=5)])
def test_category_count_rejections(counts):
    with pytest.raises(ValueError, match="num_categories"):
        normalize(num_categories=counts)


def test_nonordinal_counts_strings_and_missingness_rejections():
    panel = mixed_panel()
    with pytest.raises(ValueError, match="at least one ordinal"):
        normalize(outcome=["continuous"]*4, num_categories=5)
    panel["y"]["ordinal3"] = panel["y"]["ordinal3"].astype(str)
    with pytest.raises(ValueError, match="numeric"):
        normalize(panel)
    panel = mixed_panel()
    panel["y"]["continuous"][0, 0] = 2.
    with pytest.raises(ValueError, match="missing predecessor"):
        normalize(panel)
    normalize(panel, config=mixed_config(sur=False))
    panel = mixed_panel()
    panel["y"]["ordinal3"][2, 0] = np.inf
    with pytest.raises(ValueError, match="infinity"):
        normalize(panel)


def test_ordinal_latent_uses_sur_conditional_mean_and_variance():
    n, K = 20000, 4
    labels = np.arange(n) % K
    ordinal = _constant_state(labels, np.ones(n, bool), outcome="ordinal", num_categories=K)
    full = np.r_[-np.inf, 0, np.asarray(ordinal.cutpoints), np.inf]
    G = jnp.array([[0., 0.], [.8, 0.]])
    raw = jnp.stack([ordinal.resid, jnp.full(n, .7)])
    resid, precision = conditional_residual(raw, jnp.ones((2, n), bool), G, jnp.array([1., .4]), 0)
    view = eqx.tree_at(lambda s: s.resid, ordinal, resid)
    updated = longbet_single_step(jax.random.key(662), view, precision)
    latent = np.asarray(updated.z)
    mean = np.asarray(ordinal.z - resid)
    sd = np.asarray(jax.lax.rsqrt(precision))
    assert (sd < .7).all()  # The oracle fails decisively if sd=1 is used.
    for k in range(K):
        mask = labels == k
        m, s = mean[mask].astype(float), sd[mask].astype(float)
        a, b = (full[k]-m)/s, (full[k+1]-m)/s
        u = stats.truncnorm.cdf(latent[mask], a, b, loc=m, scale=s)
        assert stats.kstest(u, "uniform").statistic < .04
        expected_mean, expected_var = stats.truncnorm.stats(a[0], b[0], loc=m[0], scale=s[0], moments="mv")
        assert abs(latent[mask].mean()-expected_mean) < 6*np.sqrt(expected_var/mask.sum())
        assert latent[mask].var() == pytest.approx(expected_var, rel=.1)
    assert float(updated.sigma2) == 1


@pytest.fixture(scope="module")
def fitted_multi():
    panel = mixed_panel()
    model = LongBetMulti(mixed_config()).fit(**panel)
    return model, panel


def test_joint_sampling_identification_residuals_and_prediction(fitted_multi):
    model, panel = fitted_multi
    assert model.order == (1, 2, 3, 0)
    assert model.num_categories == (None, 3, None, 5)
    assert model.internal_num_categories == (3, None, 5, None)
    G = np.asarray(model.trace.gamma_loadings)
    assert (G[..., :3, :] == 0).all()
    assert np.any(G[..., 3, :3] != 0)
    assert model.state.continuous_mask == (False, False, False, True)
    for m in range(3):
        assert_invariants(model.state.states[m])
        assert (np.asarray(model.trace.traces[m].sigma2) == 1).all()
    assert model.trace.traces[0].cutpoints.shape == (2, 5, 1)
    assert model.trace.traces[1].cutpoints is None
    assert model.trace.traces[2].cutpoints.shape == (2, 5, 3)
    pred = model.predict(x=panel["x"], z=panel["z"], t=panel["t"], summary_only=True)
    for name, K in [("ordinal3", 3), ("ordinal5", 5)]:
        child = pred[name]
        assert child.prob_y_summary.mean.shape == (12, 4, K)
        assert child.att_prob_full.shape == (2, K, 10)
        assert model[name].config.num_categories == K


def test_sur_off_sharing_off_matches_scalar_steps():
    panel = mixed_panel()
    cfg = mixed_config(sur=False, num_chains=1)
    norm = normalize(panel, cfg)
    n = 48
    st = init_multi_longbet(X_unified=jnp.zeros((1,n), jnp.uint8),
        unit_idx=jnp.repeat(jnp.arange(12),4), time_idx=jnp.tile(jnp.arange(4),12),
        exposure_idx=jnp.tile(jnp.array([0,0,1,2]),12), z_vec=jnp.array(panel["z"].ravel()),
        max_split_mu=jnp.zeros(1,jnp.uint8), max_split_nu=jnp.zeros(1,jnp.uint8),
        norm_input=norm, config=cfg)
    keys = jax.random.split(jax.random.key(301), 4)
    out = jax.jit(multi_single_step)(keys, keys, st, jnp.int32(1))
    for m, child in enumerate(st.states):
        scalar = longbet_single_step(keys[m], child)
        for name in ("z", "cutpoints", "beta", "alpha", "b0", "b1", "resid", "mu_fit", "nu_fit"):
            if getattr(scalar, name) is not None:
                np.testing.assert_array_equal(getattr(out.states[m], name), getattr(scalar, name))


def test_joint_and_extracted_child_archives_roundtrip(fitted_multi, tmp_path):
    model, panel = fitted_multi
    path = tmp_path / "joint.npz"
    model.save(path)
    with np.load(path) as file:
        meta = json.loads(str(file["_meta_json"]))
        assert meta["format_version"] == 3
        assert meta["ordinal_schema_version"] == 1
        assert "outcome_0_cutpoints" in file and "outcome_2_cutpoints" in file
    loaded = LongBetMulti.load(path)
    assert loaded.num_categories == model.num_categories
    assert loaded.internal_num_categories == model.internal_num_categories
    args = dict(x=panel["x"], z=panel["z"], t=panel["t"], summary_only=True)
    a, b = model.predict(**args), loaded.predict(**args)
    for name in ("ordinal3", "ordinal5"):
        for field in ("prob_y_summary", "att_prob_full", "cutpoints_samples"):
            np.testing.assert_array_equal(getattr(a[name], field), getattr(b[name], field))
        child_path = tmp_path / f"{name}.npz"
        loaded[name].save(child_path)
        child = LongBet.load(child_path)
        np.testing.assert_array_equal(child.predict(**args).att_prob_full, a[name].att_prob_full)
        assert child.config.num_categories == loaded[name].config.num_categories


def test_scalar_effect_utilities_reject_only_ordinal_selections(fitted_multi):
    model, panel = fitted_multi
    pred = model.predict(x=panel["x"], z=panel["z"], t=panel["t"])
    for name in ("ordinal3", "ordinal5"):
        with pytest.raises(ValueError, match="category or score"):
            effect_draws(pred, name)
    for name in ("continuous", "binary"):
        assert effect_draws(pred, name).shape == (12, 4, 10)
    with pytest.raises(ValueError, match="category or score"):
        joint_prob(pred, {name: lambda x: x > 0 for name in model.outcome_names})


@pytest.mark.parametrize("corruption", ["missing", "shape", "negative", "order", "count", "child_meta",
                                      "variance", "loading", "format", "schema"])
def test_joint_archive_rejects_malformed_ordinal(fitted_multi, tmp_path, corruption):
    model, _ = fitted_multi
    path = tmp_path / "bad.npz"
    model.save(path)
    with np.load(path) as file:
        arrays = dict(file)
    meta = json.loads(str(arrays["_meta_json"]))
    if corruption == "missing": del arrays["outcome_0_cutpoints"]
    elif corruption == "shape": arrays["outcome_2_cutpoints"] = arrays["outcome_0_cutpoints"]
    elif corruption == "negative": arrays["outcome_2_cutpoints"][:] = -1
    elif corruption == "order": meta["order"] = [2, 1, 3, 0]
    elif corruption == "count": meta["internal_num_categories"] = [5, None, 3, None]
    elif corruption == "child_meta": meta["ordinal_metadata"][0]["num_categories"] = 5
    elif corruption == "variance": arrays["outcome_0_sigma2"] *= 2
    elif corruption == "loading": arrays["gamma_loadings"][..., 2, 0] = .5
    elif corruption == "format": meta["format_version"] = 2
    elif corruption == "schema": meta["ordinal_schema_version"] = 2
    arrays["_meta_json"] = json.dumps(meta)
    np.savez(path, **arrays)
    with pytest.raises(ValueError):
        LongBetMulti.load(path)


def test_two_category_joint_archives_keep_empty_threshold_axes(tmp_path):
    panel = ordinal_panel(2)
    data = dict(y={"rating": panel["y"], "binary": panel["y"].copy()},
                x=panel["x"], z=panel["z"], t=panel["t"],
                outcome=["ordinal", "binary"], num_categories=[2, None])
    model = LongBetMulti(mixed_config(num_chains=1, num_sweeps=1,
        num_burnin=1)).fit(**data)
    path = tmp_path / "binary-limit.npz"
    model.save(path)
    loaded = LongBetMulti.load(path)
    args = {name: panel[name] for name in ("x", "z", "t")}
    original, restored = model.predict(**args)["rating"], loaded.predict(**args)["rating"]
    assert restored.cutpoints_samples.shape == (1, 0)
    assert restored.prob_y.shape == (12, 4, 2, 1)
    np.testing.assert_array_equal(original.prob_y, restored.prob_y)
    with np.load(path) as file:
        arrays = dict(file)
    assert arrays["outcome_0_cutpoints"].shape == (1, 0)
    del arrays["outcome_0_cutpoints"]
    np.savez(path, **arrays)
    with pytest.raises(ValueError, match="missing outcome_0_cutpoints"):
        LongBetMulti.load(path)
