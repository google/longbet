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

"""Tests for input normalization, validation, ordering, and container parsing for LongBetMulti."""

import numpy as np
import pytest

from longbet import LongBetConfig
from longbet._multi_input import normalize_multi_inputs


def test_dict_list_3d_equivalence():
    N, T = 10, 5
    rng = np.random.default_rng(42)
    y1 = rng.standard_normal((N, T)).astype(np.float32)
    y2 = rng.standard_normal((N, T)).astype(np.float32)
    y3d = np.stack([y1, y2], axis=-1)
    cfg = LongBetConfig()

    res_dict = normalize_multi_inputs(
        {"y1": y1, "y2": y2}, outcome="continuous", outcome_names=None, config=cfg
    )
    res_list = normalize_multi_inputs(
        [y1, y2], outcome="continuous", outcome_names=["y1", "y2"], config=cfg
    )
    res_3d = normalize_multi_inputs(
        y3d, outcome="continuous", outcome_names=["y1", "y2"], config=cfg
    )

    assert res_dict.outcome_names == res_list.outcome_names == res_3d.outcome_names
    assert res_dict.order == res_list.order == res_3d.order
    for a, b, c in zip(res_dict.y_prepared, res_list.y_prepared, res_3d.y_prepared):
        np.testing.assert_allclose(a, b)
        np.testing.assert_allclose(a, c)


def test_n1_t1_support():
    # N=1, T=1 panels must be accepted
    y1 = np.array([[2.5]], dtype=np.float32)
    y2 = np.array([[1.0]], dtype=np.float32)
    cfg = LongBetConfig(standardize=False)
    res = normalize_multi_inputs(
        [y1, y2], outcome="continuous", outcome_names=None, config=cfg
    )
    assert res.N == 1
    assert res.T == 1
    assert res.M == 2
    assert res.outcome_names == ("outcome_1", "outcome_2")


def test_stable_binary_first_partition():
    # User order: cont1, bin1, cont2, bin2
    N, T = 12, 4
    rng = np.random.default_rng(1)
    y_c1 = rng.standard_normal((N, T)).astype(np.float32)
    y_b1 = rng.choice([0.0, 1.0], size=(N, T)).astype(np.float32)
    y_c2 = rng.standard_normal((N, T)).astype(np.float32)
    y_b2 = rng.choice([0.0, 1.0], size=(N, T)).astype(np.float32)

    cfg = LongBetConfig()
    y_dict = {"c1": y_c1, "b1": y_b1, "c2": y_c2, "b2": y_b2}
    types = {"c1": "continuous", "b1": "binary", "c2": "continuous", "b2": "binary"}

    res = normalize_multi_inputs(y_dict, outcome=types, outcome_names=None, config=cfg)
    assert res.outcome_names == ("c1", "b1", "c2", "b2")
    # Binary outcomes (b1, b2) must precede continuous (c1, c2) stably
    assert res.internal_names == ("b1", "b2", "c1", "c2")
    # order maps internal m -> user index
    # b1 is user index 1, b2 is user index 3, c1 is user index 0, c2 is user index 2
    assert res.order == (1, 3, 0, 2)
    # inverse_order maps user index -> internal m
    # c1(0)->2, b1(1)->0, c2(2)->3, b2(3)->1
    assert res.inverse_order == (2, 0, 3, 1)


def test_invalid_ranks_and_shapes():
    cfg = LongBetConfig()
    # Bare 2-D array is invalid
    with pytest.raises(ValueError, match="Bare 2-D panels are invalid"):
        normalize_multi_inputs(np.zeros((10, 5)), outcome="continuous", outcome_names=None, config=cfg)

    # 4-D array is invalid
    with pytest.raises(ValueError, match="must be 3-D"):
        normalize_multi_inputs(np.zeros((10, 5, 2, 2)), outcome="continuous", outcome_names=None, config=cfg)

    # M < 2
    with pytest.raises(ValueError, match="requires at least 2 outcomes"):
        normalize_multi_inputs([np.zeros((10, 5))], outcome="continuous", outcome_names=None, config=cfg)

    # Mismatched panel shapes
    with pytest.raises(ValueError, match="matching panel dimensions"):
        normalize_multi_inputs(
            [np.zeros((10, 5)), np.zeros((10, 6))],
            outcome="continuous",
            outcome_names=None,
            config=cfg,
        )


def test_infinities_rejection():
    cfg = LongBetConfig()
    y1 = np.arange(50, dtype=np.float32).reshape(10, 5)
    y2 = np.arange(50, dtype=np.float32).reshape(10, 5)
    y2[0, 0] = np.inf
    with pytest.raises(ValueError, match="contains infinite values"):
        normalize_multi_inputs({"y1": y1, "y2": y2}, outcome="continuous", outcome_names=None, config=cfg)


def test_duplicate_names_rejection():
    cfg = LongBetConfig()
    y1 = np.arange(50, dtype=np.float32).reshape(10, 5)
    y2 = np.arange(50, dtype=np.float32).reshape(10, 5)
    with pytest.raises(ValueError, match="unique"):
        normalize_multi_inputs([y1, y2], outcome="continuous", outcome_names=["a", "a"], config=cfg)


def test_invalid_binary_labels():
    cfg = LongBetConfig()
    y1 = np.arange(4, dtype=np.float32).reshape(2, 2)
    y2 = np.array([[0.0, 1.0], [2.0, 0.0]], dtype=np.float32)  # contains label 2
    with pytest.raises(ValueError, match="invalid labels"):
        normalize_multi_inputs([y1, y2], outcome=["continuous", "binary"], outcome_names=None, config=cfg)


def test_all_missing_or_zero_sd():
    cfg = LongBetConfig(standardize=True)
    y1 = np.full((10, 5), np.nan, dtype=np.float32)
    y2 = np.arange(50, dtype=np.float32).reshape(10, 5)
    with pytest.raises(ValueError, match="has no observed values"):
        normalize_multi_inputs([y1, y2], outcome="continuous", outcome_names=None, config=cfg)

    # Constant outcome with standardize=True has zero SD
    y1_const = np.ones((10, 5), dtype=np.float32)
    y2_norm = np.random.default_rng(0).standard_normal((10, 5)).astype(np.float32)
    with pytest.raises(ValueError, match="zero or non-finite standard deviation"):
        normalize_multi_inputs([y1_const, y2_norm], outcome="continuous", outcome_names=None, config=cfg)


def test_config_prior_validation():
    # sur must be bool
    with pytest.raises((TypeError, ValueError), match="sur must be a bool"):
        LongBetConfig(sur="yes")  # type: ignore

    # sur_prior_var must be finite nonnegative float, and not a bool
    with pytest.raises((TypeError, ValueError), match="sur_prior_var must be a finite nonnegative scalar"):
        LongBetConfig(sur_prior_var=True)  # type: ignore

    with pytest.raises(ValueError, match="sur_prior_var must be a finite nonnegative scalar"):
        LongBetConfig(sur_prior_var=-0.5)

    with pytest.raises(ValueError, match="sur_prior_var must be a finite nonnegative scalar"):
        LongBetConfig(sur_prior_var=float("inf"))

    # sur_active property
    cfg_on = LongBetConfig(sur=True, sur_prior_var=1.0)
    assert cfg_on.sur_active is True

    cfg_off = LongBetConfig(sur=False, sur_prior_var=1.0)
    assert cfg_off.sur_active is False

    cfg_zero = LongBetConfig(sur=True, sur_prior_var=0.0)
    assert cfg_zero.sur_active is False


@pytest.mark.parametrize("names", ["ab", ["a", None], ["a", "  "]])
def test_invalid_explicit_names(names):
    with pytest.raises(ValueError, match="strings"):
        normalize_multi_inputs([np.arange(4).reshape(2, 2)] * 2,
                               outcome="continuous", outcome_names=names, config=LongBetConfig())


@pytest.mark.parametrize("container", ["list", "dict", "array"])
def test_binary_labels_are_validated_before_float32_rounding(container):
    y = np.array([[0.0, 1.0], [0.0, 1.0 + 1e-9]])
    data = [y, y] if container == "list" else {"a": y, "b": y} if container == "dict" else np.stack([y, y], axis=-1)
    with pytest.raises(ValueError, match="invalid labels"):
        normalize_multi_inputs(data, outcome="binary", outcome_names=None, config=LongBetConfig())


def test_container_scaling_keeps_float64_precision():
    y = (1e9 + np.arange(12, dtype=float)).reshape(3, 4)
    cfg = LongBetConfig()
    as_list = normalize_multi_inputs([y, y], outcome="continuous", outcome_names=None, config=cfg)
    as_array = normalize_multi_inputs(np.stack([y, y], axis=-1), outcome="continuous", outcome_names=None, config=cfg)
    np.testing.assert_array_equal(as_list.y_prepared, as_array.y_prepared)
    assert as_list.sdy == as_array.sdy


@pytest.mark.parametrize("bad", [np.nan, np.inf, 0.5, 2.0])
def test_multi_rejects_nonbinary_treatment(bad):
    from longbet import LongBetMulti
    y = np.arange(12).reshape(3, 4)
    z = np.zeros((3, 4))
    z[:, -1] = bad
    with pytest.raises(ValueError, match="binary treatment"):
        LongBetMulti().fit([y, y], np.ones((3, 1)), z)


@pytest.mark.parametrize("t", [np.array([1, 2, 3, np.nan]), np.array([[1, 2, 3, 4]])])
def test_multi_rejects_invalid_time_vector(t):
    from longbet import LongBetMulti
    y = np.arange(12).reshape(3, 4)
    with pytest.raises(ValueError, match="finite 1-D"):
        LongBetMulti().fit([y, y], np.ones((3, 1)), np.zeros((3, 4)), t=t)
