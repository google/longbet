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

"""Draw-paired ordinal estimands, bounded summaries, and archive round-trips."""
import dataclasses
import json

import numpy as np
import pytest
from scipy import stats

from longbet import LongBet, effect_draws, effect_draws_from_arrays
from longbet._ordinal import category_probabilities
from longbet._summary import choose_ordinal_block_size
from test_ordinal_model import ordinal_config, ordinal_panel


@pytest.fixture(scope="module", params=[2, 3, 5])
def fitted(request):
    K = request.param
    data = ordinal_panel(K)
    model = LongBet(ordinal_config(K, num_chains=2)).fit(**data)
    return model, data


def test_probability_shapes_algebra_and_category_att(fitted):
    model, data = fitted
    pred = model.predict(**data, block_size=7)
    K, D = model.config.num_categories, 10
    N, T = data["y"].shape
    assert pred.num_categories == K
    np.testing.assert_array_equal(pred.categories, np.arange(K))
    assert pred.cutpoints_samples.shape == (D, K-2)
    np.testing.assert_array_equal(pred.cutpoints_samples, np.asarray(model.trace.cutpoints).reshape(D, K-2))
    assert pred.att_prob_full.shape == (2, K, D)
    for name, mean in [("prob_y", pred.yhats), ("prob_mu0", pred.muhats0)]:
        p = getattr(pred, name)
        assert p.shape == (N, T, K, D)
        expected = category_probabilities(mean.reshape(N*T, D).T, pred.cutpoints_samples)
        np.testing.assert_allclose(p, expected.transpose(1, 2, 0).reshape(N, T, K, D), atol=2e-15)
        np.testing.assert_allclose(p.sum(axis=2), 1, atol=4e-16)
    treated = category_probabilities((pred.muhats0+pred.tauhats).reshape(N*T, D).T,
                                     pred.cutpoints_samples).transpose(1, 2, 0).reshape(N, T, K, D)
    np.testing.assert_allclose(pred.prob_tau, treated-pred.prob_mu0, atol=2e-15)
    np.testing.assert_allclose(pred.prob_tau.sum(axis=2), 0, atol=4e-16)
    for s in (1, 2):
        membership = (pred.z == 1) & (pred.s == s)
        assert pred.att_counts[s-1] == membership.sum()
        np.testing.assert_allclose(pred.att_prob_full[s-1], pred.prob_tau[membership].mean(0), atol=1e-15)
    # One treated training cell has a missing y, but it belongs in prediction ATT.
    assert ((pred.z == 1) & np.isnan(data["y"])).any()
    if K == 2:
        np.testing.assert_allclose(pred.prob_y[:, :, 1], stats.norm.cdf(pred.yhats), atol=1e-15)


def test_full_and_summary_modes_match_exact_reductions(fitted):
    model, data = fitted
    full = model.predict(**data, block_size=48, alpha=.2)
    for block in (1, 7, 48):
        pred = model.predict(**data, block_size=block, summary_only=True, alpha=.2)
        assert pred.prob_y is pred.prob_mu0 is pred.prob_tau is None
        for arm, name in [("factual", "prob_y"), ("control", "prob_mu0"), ("effect", "prob_tau")]:
            draws = getattr(full, name)
            summary = pred.predict_probabilities(arm)
            assert summary.mean.shape == (*data["y"].shape, model.config.num_categories)
            want = [draws.mean(-1), draws.std(-1), *np.percentile(draws, [10, 90], axis=-1)]
            for actual, expected in zip(summary, want):
                np.testing.assert_allclose(actual, expected, atol=2e-15)
            np.testing.assert_array_equal(full.predict_probabilities(arm, summary=False), draws)
            with pytest.raises(ValueError, match="discarded"):
                pred.predict_probabilities(arm, summary=False)
        np.testing.assert_allclose(full.att_prob_full, pred.att_prob_full, atol=1e-15)
        for name in ("tau_summary", "mu0_summary", "y_summary", "att_full"):
            np.testing.assert_allclose(getattr(full, name), getattr(pred, name), atol=2e-6)


def test_att_expected_score_and_public_errors(fitted):
    model, data = fitted
    pred = model.predict(**data, summary_only=True)
    K = pred.num_categories
    cat = pred.att_probabilities(alpha=.1)
    assert cat["att"].shape == (2, K) and cat["intervals"].shape == (2, 2, K)
    np.testing.assert_allclose(cat["att"], pred.att_prob_full.mean(-1), atol=1e-15)
    np.testing.assert_allclose(cat["intervals"], np.percentile(pred.att_prob_full, [5, 95], axis=-1))
    for weights in (None, np.eye(K)[-1], np.arange(K)[::-1], np.ones(K)*3):
        result = pred.att_expected_score(weights, alpha=.1)
        w = np.arange(K) if weights is None else weights
        expected = np.sum(pred.att_prob_full * w[None, :, None], axis=1)
        np.testing.assert_allclose(result["att_full"], expected, atol=1e-15)
        np.testing.assert_allclose(result["intervals"], np.percentile(expected, [5, 95], axis=-1), atol=1e-15)
    np.testing.assert_allclose(pred.att_expected_score(np.ones(K))["att_full"], 0, atol=5e-16)
    for weights in ([0], np.ones((K, 1)), [np.inf]*K, [np.nan]*K, ["1"]*K):
        with pytest.raises(ValueError, match="weights"):
            pred.att_expected_score(weights)
    for alpha in (0, 1, np.nan, True, "0.1"):
        with pytest.raises(ValueError, match="alpha"):
            pred.att_probabilities(alpha)
    with pytest.raises(ValueError, match="arm"):
        pred.predict_probabilities("treated")
    with pytest.raises(ValueError, match="category or score"):
        effect_draws_from_arrays(None, outcome="ordinal", summary_only=True)


def test_empty_exposure_groups_are_nan(fitted):
    model, data = fitted
    pred = model.predict(x=data["x"], z=np.zeros_like(data["z"]), t=data["t"], summary_only=True)
    assert pred.att_prob_full.shape == (1, model.config.num_categories, 10)
    assert pred.att_counts[0] == 0
    for method in (pred.att_probabilities, pred.att_expected_score):
        result = method()
        for name in ("att", "att_full", "intervals"):
            assert np.isnan(result[name]).all()


def test_summary_memory_is_bounded_and_accounts_for_float64_categories(fitted, monkeypatch):
    import longbet._model as module
    import longbet._summary as summary_module
    model, data = fitted
    D, M, K = 10, data["y"].size, model.config.num_categories
    monkeypatch.setattr(summary_module, "_TARGET_BLOCK_ELEMENTS", 16000)
    budget_block = choose_ordinal_block_size(D, M, K)
    assert budget_block < M
    assert choose_ordinal_block_size(D, M*1000, K) == budget_block
    assert choose_ordinal_block_size(D, M, K+2) <= budget_block
    calls, allocations = [], []
    original = module.category_probabilities
    def transform(mean, cp):
        calls.append(mean.shape)
        assert mean.shape[0] == D and mean.shape[1] <= budget_block
        return original(mean, cp)
    monkeypatch.setattr(module, "category_probabilities", transform)
    for method in ("empty", "zeros"):
        orig = getattr(np, method)
        def allocator(shape, *args, _orig=orig, **kwargs):
            if isinstance(shape, tuple):
                allocations.append(shape)
            return _orig(shape, *args, **kwargs)
        monkeypatch.setattr(np, method, allocator)
    pred = model.predict(**data, summary_only=True, cache_forest_evaluations=False)
    assert calls and max(b for _, b in calls) == budget_block
    for shape in [(D, M), (M, D), (D, M, K), (M, K, D)]:
        assert shape not in allocations
    assert pred.prob_y is pred.yhats is None
    assert not hasattr(model, "_eval_cache")


def test_scalar_ordinal_roundtrip_including_projection_and_custom_score(fitted, tmp_path):
    model, data = fitted
    path = tmp_path / "ordinal.npz"
    model.save(path)
    loaded = LongBet.load(path)
    assert loaded.config == model.config
    z = np.c_[data["z"], data["z"][:, -1:]]
    kw = dict(x=data["x"], z=z, t=np.arange(1, 6), summary_only=True, block_size=7)
    a, b = model.predict(**kw), loaded.predict(**kw)
    for name in ("prob_y_summary", "prob_mu0_summary", "prob_tau_summary", "att_prob_full", "cutpoints_samples"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
    np.testing.assert_array_equal(a.att_expected_score(np.eye(a.num_categories)[-1])["att_full"],
                                  b.att_expected_score(np.eye(a.num_categories)[-1])["att_full"])


def test_even_binary_ordinal_archives_require_the_cutpoint_array(fitted, tmp_path):
    model, _ = fitted
    path = tmp_path / "missing-cutpoints.npz"
    model.save(path)
    with np.load(path) as file:
        arrays = dict(file)
    del arrays["cutpoints"]
    np.savez(path, **arrays)
    with pytest.raises(ValueError, match="missing cutpoints"):
        LongBet.load(path)


def test_prediction_pairs_both_forest_inputs_and_scales_the_offset_once(fitted, monkeypatch):
    import longbet._model as module
    model, data = fitted
    D, M = 10, data["y"].size
    row = 0
    for block in model.design_.blocks:
        if block.name == "s":
            break
        row += block.n_cols
    calls = []
    mu = (np.float32(model.offset_) + np.arange(D, dtype=np.float32)*.07)[:, None]
    def forest_values(X, lo, hi, width, trace):
        if trace is model.trace.mu_trace:
            return np.broadcast_to(mu, (D,hi-lo))
        calls.append(X.copy())
        return (.4 + np.arange(D,dtype=np.float32)[:,None]*.02 + .3*X[None,row,lo:hi]).astype(np.float32)
    monkeypatch.setattr(module, "_eval_block", forest_values)
    pred = model.predict(**data, block_size=M)
    # The treatment forest never sees the exposure index, so one evaluation serves both arms.
    assert len(calls) == 1
    factual = control = calls[0]
    tr = model.trace
    alpha = np.asarray(tr.alpha).reshape(D,1)
    beta = np.asarray(tr.beta).reshape(D,-1)
    b0, b1 = np.asarray(tr.b0).reshape(D,1), np.asarray(tr.b1).reshape(D,1)
    gamma = np.asarray(tr.gamma).reshape(D,12)[:,np.repeat(np.arange(12),4)]
    f0 = (.4 + np.arange(D,dtype=np.float32)[:,None]*.02 + .3*control[None,row]).astype(np.float32)
    f1 = (.4 + np.arange(D,dtype=np.float32)[:,None]*.02 + .3*factual[None,row]).astype(np.float32)
    trend = 0.0
    if getattr(tr, "trend_coef", None) is not None:
        from longbet._trend_basis import apply_unit_features, evaluate_trend
        trend = evaluate_trend(
            np.asarray(tr.trend_coef).reshape(D, -1),
            apply_unit_features(np.asarray(data["x"], dtype=np.float32), model.trend_feature_spec_),
            np.asarray(model.trend_time_basis_, dtype=np.float32),
            np.repeat(np.arange(12, dtype=np.int32), 4),
            np.tile(np.arange(4, dtype=np.int32), 12),
            period=np.asarray(model.trend_period_basis_, dtype=np.float32),
        )
    eta0 = alpha*mu + trend + b0*beta[:,:1]*f0 + gamma          # b0 = 0: the untreated surface is mu + trend + gamma
    eta1 = alpha*mu + trend + b1*beta[:,pred.s.ravel()]*f1*(pred.s.ravel() > 0) + gamma
    np.testing.assert_allclose(pred.muhats0.reshape(M,D).T, eta0, atol=2e-6)
    # Summing mu0+tau may differ by float32 rounding from computing eta1 directly.
    actual1 = pred.muhats0 + pred.tauhats
    np.testing.assert_allclose(actual1.reshape(M,D).T,eta1,atol=2e-6)
    p0 = category_probabilities(pred.muhats0.reshape(M,D).T,pred.cutpoints_samples)
    p1 = category_probabilities(actual1.reshape(M,D).T,pred.cutpoints_samples)
    np.testing.assert_allclose(pred.prob_tau.reshape(M,pred.num_categories,D).transpose(2,0,1),p1-p0,atol=3e-15)


@pytest.mark.parametrize("corruption", ["missing", "shape", "nan", "negative", "reversed",
    "variance", "precision", "chain", "count", "scale", "anchor", "schema", "old_cache"])
def test_scalar_archive_rejects_corruption(tmp_path, corruption):
    # Reuse one fit across corruption cases without sharing mutated archives.
    model, _ = _archive_fixture()
    path = tmp_path / "bad.npz"
    model.save(path)
    with np.load(path) as file:
        data = dict(file)
    meta = json.loads(str(data["_meta_json"]))
    if corruption == "missing": del data["cutpoints"]
    elif corruption == "shape": data["cutpoints"] = data["cutpoints"][:, 0]
    elif corruption == "nan": data["cutpoints"][:] = np.nan
    elif corruption == "negative": data["cutpoints"][:] = -1
    elif corruption == "reversed": data["cutpoints"] = data["cutpoints"][..., ::-1]
    elif corruption == "variance": data["sigma2"] *= 2
    elif corruption == "precision": data["mu_error_cov_inv"] *= 2
    elif corruption == "chain": meta["has_chains"] = False
    elif corruption == "count": meta["num_categories"] = 3
    elif corruption == "scale": meta["cutpoint_prior_scale"] = 2.
    elif corruption == "anchor": meta["cutpoint_anchor"] = "other"
    elif corruption == "schema": meta["ordinal_schema_version"] = 2
    elif corruption == "old_cache": del meta["precision_cache_version"]
    data["_meta_json"] = json.dumps(meta)
    np.savez(path, **data)
    with pytest.raises(ValueError, match="archive|Archive"):
        LongBet.load(path)


from functools import lru_cache

@lru_cache(maxsize=1)
def _archive_fixture():
    data = ordinal_panel(5)
    return LongBet(ordinal_config(5, num_chains=2)).fit(**data), data


def test_legacy_nonordinal_archive_and_methods(tmp_path):
    cfg = dataclasses.replace(ordinal_config(2), outcome="binary", num_categories=None)
    data = ordinal_panel(2)
    model = LongBet(cfg).fit(**data)
    path = tmp_path / "old.npz"
    model.save(path)
    with np.load(path) as file:
        arrays = dict(file)
    config = json.loads(str(arrays["_config_json"]))
    del config["num_categories"], config["cutpoint_prior_scale"]
    arrays["_config_json"] = json.dumps(config)
    np.savez(path, **arrays)
    loaded = LongBet.load(path)
    assert loaded.trace.cutpoints is None
    pred = loaded.predict(**data)
    for name in ("num_categories", "categories", "cutpoints_samples", "prob_y",
                 "prob_mu0", "prob_tau", "prob_y_summary", "att_prob_full"):
        assert getattr(pred, name) is None
    for method in (pred.predict_probabilities, pred.att_probabilities, pred.att_expected_score):
        with pytest.raises(ValueError, match="ordinal"):
            method()
    np.testing.assert_array_equal(pred.yhats, model.predict(**data).yhats)
