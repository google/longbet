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

"""Tests for NPZ persistence of LongBetMulti models, roundtrip invariance, and cross-kind rejection."""

from pathlib import Path
import json
import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, LongBetMulti, effect_draws, joint_prob, load_multi_npz, load_npz, outcome_correlation


def test_multi_npz_roundtrip(tmp_path: Path):
    """Verify that a fitted LongBetMulti model roundtrips through save/load with identical outputs."""
    N, T, P = 8, 4, 2
    rng = np.random.default_rng(333)

    y_bin = rng.choice([0.0, 1.0], size=(N, T)).astype(np.float32)
    y_cont = rng.standard_normal((N, T)).astype(np.float32)
    x = rng.standard_normal((N, P)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[:, 1:] = 1.0

    cfg = LongBetConfig(
        sigma_prior_a=2, sigma_prior_b=1,
        num_sweeps=6,
        num_burnin=3,
        num_chains=2,
        num_trees_pr=3,
        num_trees_trt=3,
        random_seed=7,
        sur=True,
    )

    y_data = {"c": y_cont, "b": y_bin}
    types = {"c": "continuous", "b": "binary"}

    model = LongBetMulti(cfg).fit(y=y_data, x=x, z=z, outcome=types)
    pred_before = model.predict(x=x, z=z)
    eff_before = effect_draws(pred_before, outcome="c")
    corr_before = outcome_correlation(model)

    save_path = tmp_path / "model_multi.npz"
    model.save(save_path)

    loaded = LongBetMulti.load(save_path)

    # Validate metadata
    assert loaded.outcome_names == model.outcome_names
    assert loaded.outcome == model.outcome
    assert loaded.order == model.order
    assert loaded.inverse_order == model.inverse_order
    assert loaded.sur_active == model.sur_active
    assert loaded.sampler_semantics == model.sampler_semantics
    assert loaded.provenance == model.provenance
    assert loaded.sampler_semantics == "full_precision_sur_v1"

    # Validate loading draws
    np.testing.assert_allclose(loaded.Gamma_draws, model.Gamma_draws)

    # Validate predictions from loaded model
    pred_after = loaded.predict(x=x, z=z)
    np.testing.assert_allclose(pred_after["c"].tauhats, pred_before["c"].tauhats)
    np.testing.assert_allclose(pred_after["b"].tauhats, pred_before["b"].tauhats)

    # Validate correlation
    corr_after = outcome_correlation(loaded)
    np.testing.assert_allclose(corr_after, corr_before)

    # A corrupt permutation/draw axis must fail at load, before wrong outcomes
    # or chains can be silently paired by a utility.
    with np.load(save_path) as archive:
        arrays = {key: archive[key] for key in archive.files}
    meta = json.loads(str(arrays["_meta_json"]))
    improper = dict(arrays)
    old_config = json.loads(str(improper["_config_json"]))
    old_config.update(sigma_prior_a=0, sigma_prior_b=0)
    improper["_config_json"] = json.dumps(old_config)
    old_path = tmp_path / "improper_target.npz"
    np.savez_compressed(old_path, **improper)
    with pytest.raises(ValueError):
        LongBetMulti.load(old_path)
    for name in model.outcome_names:
        child_path = tmp_path / f"improper_child_{name}.npz"
        model[name].save(child_path)
        with np.load(child_path) as child_archive:
            child = {key: child_archive[key] for key in child_archive.files}
        cfg = json.loads(str(child["_config_json"]))
        cfg.update(sigma_prior_a=0, sigma_prior_b=0)
        child["_config_json"] = json.dumps(cfg)
        np.savez_compressed(child_path, **child)
        with pytest.raises(ValueError):
            LongBet.load(child_path)
    assert meta["precision_cache_version"] == 1
    legacy = dict(meta)
    legacy.pop("precision_cache_version")
    legacy["sampler_semantics"] = "recursive_sur_v1"
    legacy_path = tmp_path / "legacy_precision_multi.npz"
    np.savez_compressed(legacy_path, **{**arrays, "_meta_json": json.dumps(legacy)})
    with pytest.raises(ValueError, match="Refit from the original data"):
        LongBetMulti.load(legacy_path)

    child_path = tmp_path / "legacy_precision_child.npz"
    model["c"].save(child_path)
    with np.load(child_path) as archive:
        child_arrays = {key: archive[key] for key in archive.files}
    child_meta = json.loads(str(child_arrays["_meta_json"]))
    child_meta.pop("precision_cache_version")
    child_arrays["_meta_json"] = json.dumps(child_meta)
    np.savez_compressed(child_path, **child_arrays)
    with pytest.raises(ValueError, match="Refit from the original data"):
        LongBet.load(child_path)
    for key, value in [("inverse_order", [0, 1]), ("sdy", [0, 1]),
                       ("outcome_names", ["c", "c"]), ("has_chains", False)]:
        corrupted = dict(arrays)
        corrupted["_meta_json"] = json.dumps({**meta, key: value})
        bad = tmp_path / f"bad_{key}.npz"
        np.savez_compressed(bad, **corrupted)
        with pytest.raises(ValueError, match="Invalid multi archive"):
            LongBetMulti.load(bad)
    corrupted = dict(arrays)
    corrupted["outcome_0_sigma2"] = arrays["outcome_0_sigma2"][:, :-1]
    bad = tmp_path / "bad_draw_shape.npz"
    np.savez_compressed(bad, **corrupted)
    with pytest.raises(ValueError, match="Invalid multi archive.*shape"):
        LongBetMulti.load(bad)

    # Legacy multi fits silently omitted the calendar and exposure cutpoints.
    # Refuse them (and extracted scalar children), rather than reuse bad draws.
    from copy import deepcopy
    for axis in ("t", "s"):
        broken = deepcopy(model.design_.to_dict())
        block = next(b for b in broken["blocks"] if b["name"] == axis)
        block["splits"] = [[]]
        block["max_split"] = [0]
        corrupted = dict(arrays)
        corrupted["_meta_json"] = json.dumps({**meta, "design_spec": broken})
        bad = tmp_path / f"legacy_multi_{axis}.npz"
        np.savez_compressed(bad, **corrupted)
        with pytest.raises(ValueError, match="Refit from the original data"):
            LongBetMulti.load(bad)

        child_path = tmp_path / f"legacy_child_{axis}.npz"
        model["c"].save(child_path)
        with np.load(child_path) as archive:
            child_arrays = {key: archive[key] for key in archive.files}
        child_meta = json.loads(str(child_arrays["_meta_json"]))
        child_arrays["_meta_json"] = json.dumps({**child_meta, "design": broken})
        np.savez_compressed(child_path, **child_arrays)
        with pytest.raises(ValueError):
            LongBet.load(child_path)


def test_cross_kind_archive_rejection(tmp_path: Path):
    """Verify that load_npz rejects multi archives and load_multi_npz rejects scalar archives."""
    N, T = 6, 4
    rng = np.random.default_rng(444)
    x = rng.standard_normal((N, 2)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)

    # 1. Scalar fit and save
    y_scalar = rng.standard_normal((N, T)).astype(np.float32)
    cfg_scalar = LongBetConfig(num_sweeps=3, num_burnin=2, num_chains=1, num_trees_pr=2, num_trees_trt=2)
    scalar_model = LongBet(cfg_scalar).fit(y=y_scalar, x=x, z=z)
    scalar_path = tmp_path / "scalar_model.npz"
    scalar_model.save(scalar_path)

    # load_multi_npz must reject scalar archive with clear error
    with pytest.raises(ValueError, match="single-outcome model"):
        load_multi_npz(scalar_path)

    # 2. Multi fit and save
    y_multi = {"y1": y_scalar, "y2": y_scalar + 1.0}
    multi_model = LongBetMulti(cfg_scalar, sigma_prior_a=2, sigma_prior_b=1).fit(y=y_multi, x=x, z=z)
    multi_path = tmp_path / "multi_model.npz"
    multi_model.save(multi_path)

    # load_npz must reject multi archive with clear error
    with pytest.raises(ValueError, match="multi-outcome LongBet model"):
        load_npz(multi_path)
