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

"""Tests for LongBetMulti prediction, child extraction, counterfactuals, and option forwarding."""

import numpy as np
import pytest

from longbet import LongBetConfig, LongBetMulti, get_att, get_catt


def test_predict_child_apis_and_access():
    """Verify child prediction indexing, ATT/CATT getters, and summary_only mode."""
    N, T, P = 8, 4, 3
    rng = np.random.default_rng(55)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    X = rng.standard_normal((N, P)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[0, 2:] = 1.0

    cfg = LongBetConfig(
        sigma_prior_a=2, sigma_prior_b=1,
        num_sweeps=8,
        num_burnin=4,
        num_chains=1,
        num_trees_pr=4,
        num_trees_trt=4,
        random_seed=42,
    )

    model = LongBetMulti(cfg).fit(y={"rev": y1, "cost": y2}, x=X, z=z)

    # Fits child indexing
    assert len(model) == 2
    assert model["rev"] is model.fits["rev"]
    assert model[0] is model.fits[0]
    with pytest.raises(KeyError):
        _ = model["unknown"]

    # Full prediction
    pred = model.predict(x=X, z=z)
    assert len(pred) == 2
    assert pred["rev"] is pred.preds["rev"]
    assert pred[1] is pred.preds[1]
    with pytest.raises(KeyError):
        _ = pred["unknown"]

    # Child getters
    att_rev = get_att(pred["rev"])
    assert att_rev is not None
    catt_cost = get_catt(pred["cost"])
    assert catt_cost["catt"].shape == (N, T)

    # Summary only prediction
    pred_summary = model.predict(x=X, z=z, summary_only=True)
    assert pred_summary["rev"].summary_only
    assert pred_summary["rev"].tauhats is None
    assert pred_summary["rev"].tau_summary.mean.shape == (N, T)


def test_predict_with_all_optional_design_blocks():
    """Verify prediction when x_trt, x_tv, x_trt_tv, and ps are all present."""
    N, T, P = 6, 4, 2
    rng = np.random.default_rng(66)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    x = rng.standard_normal((N, P)).astype(np.float32)
    x_trt = rng.standard_normal((N, P + 1)).astype(np.float32)
    x_tv = rng.standard_normal((N, T, 1)).astype(np.float32)
    x_trt_tv = rng.standard_normal((N, T, 1)).astype(np.float32)
    ps = rng.uniform(0.1, 0.9, size=(N, T)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[0, 1:] = 1.0

    cfg = LongBetConfig(
        sigma_prior_a=2, sigma_prior_b=1,
        num_sweeps=4,
        num_burnin=2,
        num_chains=1,
        num_trees_pr=3,
        num_trees_trt=3,
        random_seed=12,
    )

    model = LongBetMulti(cfg).fit(
        y={"y1": y1, "y2": y2},
        x=x,
        z=z,
        x_trt=x_trt,
        x_tv=x_tv,
        x_trt_tv=x_trt_tv,
        ps=ps,
    )

    pred = model.predict(
        x=x,
        z=z,
        x_trt=x_trt,
        x_tv=x_tv,
        x_trt_tv=x_trt_tv,
        ps=ps,
        cache_forest_evaluations=True,
    )
    assert pred["y1"].tauhats.shape == (N, T, 4)
    assert pred["y2"].tauhats.shape == (N, T, 4)
