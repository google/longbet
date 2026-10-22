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

"""Tests for multi-outcome missingness policies and placeholder invariance."""

import numpy as np
import pytest

from longbet import LongBetConfig
from longbet._multi_input import normalize_multi_inputs


def test_complete_panel():
    N, T = 8, 4
    rng = np.random.default_rng(10)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    cfg = LongBetConfig()
    norm = normalize_multi_inputs({"y1": y1, "y2": y2}, outcome="continuous", outcome_names=None, config=cfg)
    assert np.all(norm.obs_masks[0])
    assert np.all(norm.obs_masks[1])


def test_common_missingness():
    N, T = 8, 4
    rng = np.random.default_rng(11)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    # Mark cell (0, 0) and (2, 3) missing in both
    y1[0, 0] = np.nan
    y2[0, 0] = np.nan
    y1[2, 3] = np.nan
    y2[2, 3] = np.nan

    cfg = LongBetConfig()
    norm = normalize_multi_inputs({"y1": y1, "y2": y2}, outcome="continuous", outcome_names=None, config=cfg)
    assert not norm.obs_masks[0][0, 0]
    assert not norm.obs_masks[1][0, 0]
    assert not norm.obs_masks[0][2, 3]
    assert not norm.obs_masks[1][2, 3]


def test_nested_missingness_supported():
    N, T = 8, 4
    rng = np.random.default_rng(12)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    # y1 is observed everywhere
    # y2 has missing cells (subset of y1's observed cells)
    y2[1, 1] = np.nan
    y2[3, 2] = np.nan

    cfg = LongBetConfig()
    # Order: y1 (0), y2 (1). Since continuous y2 is downstream, y2's observed cells
    # are contained within y1's observed cells (which is all cells). This must succeed.
    norm = normalize_multi_inputs({"y1": y1, "y2": y2}, outcome="continuous", outcome_names=None, config=cfg)
    assert np.all(norm.obs_masks[0])
    assert not norm.obs_masks[1][1, 1]


def test_unsupported_predecessor_missingness_rejected():
    N, T = 8, 4
    rng = np.random.default_rng(13)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    # y1 has missing cell at (0, 1), but downstream y2 is observed at (0, 1)
    y1[0, 1] = np.nan

    cfg = LongBetConfig(sur=True, sur_prior_var=1.0)
    with pytest.raises(ValueError, match="missing predecessor"):
        normalize_multi_inputs({"y1": y1, "y2": y2}, outcome="continuous", outcome_names=None, config=cfg)


def test_arbitrary_missingness_allowed_when_sur_disabled():
    N, T = 8, 4
    rng = np.random.default_rng(14)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    # y1 is missing where y2 is observed
    y1[0, 1] = np.nan
    y2[1, 0] = np.nan

    # With sur=False, this must succeed!
    cfg_off = LongBetConfig(sur=False)
    norm = normalize_multi_inputs({"y1": y1, "y2": y2}, outcome="continuous", outcome_names=None, config=cfg_off)
    assert not norm.obs_masks[0][0, 1]
    assert not norm.obs_masks[1][1, 0]

    # With sur_prior_var=0, SUR is inactive, so this must also succeed
    cfg_zero = LongBetConfig(sur=True, sur_prior_var=0.0)
    norm_z = normalize_multi_inputs({"y1": y1, "y2": y2}, outcome="continuous", outcome_names=None, config=cfg_zero)
    assert not norm_z.obs_masks[0][0, 1]


def test_placeholder_invariance():
    N, T = 6, 4
    rng = np.random.default_rng(15)
    y1_base = rng.standard_normal((N, T)).astype(np.float32)
    y2_base = rng.standard_normal((N, T)).astype(np.float32)
    y1_base[0, 0] = np.nan
    y2_base[0, 0] = np.nan

    cfg = LongBetConfig()
    norm1 = normalize_multi_inputs(
        {"y1": y1_base, "y2": y2_base}, outcome="continuous", outcome_names=None, config=cfg
    )

    # In y_prepared, missing cells are replaced with 0.0
    assert norm1.y_prepared[0][0, 0] == 0.0
    assert norm1.y_prepared[1][0, 0] == 0.0
    # Observed statistics (meany, sdy) ignore missing cell
    assert norm1.meany[0] == pytest.approx(np.mean(y1_base[norm1.obs_masks[0]]))
    assert norm1.sdy[0] == pytest.approx(np.std(y1_base[norm1.obs_masks[0]], ddof=0))
