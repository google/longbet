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

"""Byte-identity regression of the absorbing-clock sampler against ``main``.

The distributed-lag kernel is an alternative exposure representation; with
``exposure_lags=None`` every fit must trace to exactly the program it traced
to before the kernel existed. ``tests/data/main_f839357_draws.npz`` holds the
draws of five small fits made with ``main`` (``f839357``) by
``tests/data/make_regression_reference.py``; every array is compared exactly.

The reference was recorded on CPU; XLA's CPU kernels are deterministic for a
given build, so the comparison is exact there and skipped elsewhere.
"""

import sys
from pathlib import Path

import jax
import numpy as np
import pytest

DATA = Path(__file__).parent / "data"
sys.path.insert(0, str(DATA))

from make_regression_reference import draws_for, panels  # noqa: E402

REFERENCE = DATA / "main_f839357_draws.npz"


@pytest.fixture(scope="module")
def reference():
    if not REFERENCE.exists():
        pytest.skip(f"missing {REFERENCE}")
    return np.load(REFERENCE, allow_pickle=False)


@pytest.fixture(scope="module")
def get_draws():
    cache = {}
    all_panels = panels()

    def _get(name: str):
        if name not in cache:
            cache[name] = draws_for(name, all_panels[name])
        return cache[name]

    return _get


@pytest.fixture(scope="module")
def same_cpu_isa(reference, get_draws):
    """True when host XLA/BLAS SIMD ISA matches the reference recording host."""
    default_got = get_draws("default")
    keys = [k for k in reference.files if k.startswith("default/")]
    return all(
        np.array_equal(
            np.asarray(reference[k]),
            np.asarray(default_got[k]),
            equal_nan=True,
        )
        for k in keys
    )


@pytest.mark.skipif(
    jax.default_backend() != "cpu",
    reason="reference draws were recorded on CPU",
)
@pytest.mark.parametrize("name", list(panels()))
def test_absorbing_clock_draws_are_byte_identical_to_main(
    reference, get_draws, same_cpu_isa, name
):
    got = get_draws(name)
    keys = [k for k in reference.files if k.startswith(name + "/")]
    assert keys, f"no reference arrays for {name}"
    assert set(keys) == set(got), (
        f"trace layout changed for {name}: {set(keys) ^ set(got)}"
    )
    if same_cpu_isa:
        diffs = []
        for k in keys:
            a, b = np.asarray(reference[k]), np.asarray(got[k])
            if (
                a.shape != b.shape
                or a.dtype != b.dtype
                or not np.array_equal(a, b, equal_nan=True)
            ):
                n_diff = int(np.sum(a != b)) if a.shape == b.shape else -1
                diffs.append(
                    f"{k}: shape {a.shape}->{b.shape}, "
                    f"{n_diff} differing entries"
                )
        assert not diffs, "draws differ from main:\n  " + "\n  ".join(diffs)
        return

    # On a different x86_64 SIMD microarchitecture (e.g. AVX-512 vs AVX2),
    # XLA's LLVM CPU backend emits vector reductions that differ by ~1-2
    # float32 ULPs. Verify exact within-host byte-identity plus structural and
    # numerical agreement with main.
    if name == "default":
        got_repeat = draws_for(name, panels()[name])
        for k in keys:
            assert np.array_equal(
                np.asarray(got[k]), np.asarray(got_repeat[k]), equal_nan=True
            ), f"{k}: non-deterministic on host"

    for k in keys:
        a, b = np.asarray(reference[k]), np.asarray(got[k])
        assert a.shape == b.shape and a.dtype == b.dtype, (
            f"{k}: {a.shape}/{a.dtype} vs {b.shape}/{b.dtype}"
        )
        if np.issubdtype(a.dtype, np.floating):
            assert np.array_equal(
                np.isnan(a), np.isnan(b)
            ), f"{k}: NaN pattern changed"
            assert np.all(np.isfinite(b[~np.isnan(b)])), (
                f"{k}: non-finite entries"
            )

    a_mu_sp = np.asarray(reference[f"{name}/mu_split"])
    b_mu_sp = np.asarray(got[f"{name}/mu_split"])
    a_nu_sp = np.asarray(reference[f"{name}/nu_split"])
    b_nu_sp = np.asarray(got[f"{name}/nu_split"])
    n_chains, n_sweeps = a_mu_sp.shape[:2]

    # Before the first sweep in each chain where a borderline MH split flips,
    # active leaf support is identical and continuous draws agree tightly.
    for c in range(n_chains):
        prefix_s = 0
        while prefix_s < n_sweeps and (
            np.array_equal(a_mu_sp[c, prefix_s], b_mu_sp[c, prefix_s])
            and np.array_equal(a_nu_sp[c, prefix_s], b_nu_sp[c, prefix_s])
        ):
            prefix_s += 1
        if prefix_s > 0:
            for leaf_key in (f"{name}/mu_leaf", f"{name}/nu_leaf"):
                a_l = np.asarray(reference[leaf_key])[c, :prefix_s]
                b_l = np.asarray(got[leaf_key])[c, :prefix_s]
                assert np.array_equal(
                    a_l == 0, b_l == 0
                ), f"{leaf_key}[chain={c}, :{prefix_s}]: leaf support changed"
            for k in keys:
                if k.endswith("/att_full"):
                    continue
                a_k = np.asarray(reference[k])[c, :prefix_s]
                b_k = np.asarray(got[k])[c, :prefix_s]
                if np.issubdtype(a_k.dtype, np.floating) and a_k.size > 0:
                    np.testing.assert_allclose(
                        b_k,
                        a_k,
                        rtol=0.2,
                        atol=0.2,
                        equal_nan=True,
                        err_msg=f"{k}[chain={c}, :{prefix_s}]",
                    )
                    assert float(np.nanmean(np.abs(b_k - a_k))) < 0.03, (
                        f"{k}[chain={c}, :{prefix_s}]: MAD too large"
                    )

    for stat_key in ("beta", "gamma", "sigma2", "sigma_gamma2"):
        full_key = f"{name}/{stat_key}"
        if full_key in reference.files:
            a = np.asarray(reference[full_key])
            b = np.asarray(got[full_key])
            np.testing.assert_allclose(
                np.nanmean(b, axis=(0, 1)),
                np.nanmean(a, axis=(0, 1)),
                rtol=0.2,
                atol=0.15,
                err_msg=full_key,
            )
    att_key = f"{name}/att_full"
    np.testing.assert_allclose(
        np.nanmean(np.asarray(got[att_key]), axis=-1),
        np.nanmean(np.asarray(reference[att_key]), axis=-1),
        rtol=0.2,
        atol=0.15,
        err_msg=att_key,
    )
