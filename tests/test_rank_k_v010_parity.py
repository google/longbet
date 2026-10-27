# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Bit-for-bit parity of the v0.10 opt-outs against a stored v0.9.0 fixture.

``tests/data/rank_k_v090_parity.npz`` was generated with longbet-jax 0.9.0 (commit
``2805156``) by running this module as a script before any v0.10 change. Every
scenario below is re-run with the v0.10 opt-out keyword arguments that restore the
0.9.0 behavior (``composite_mode="bottom_up"`` and friends) and must reproduce the
stored draws exactly.

Regenerate (only on a 0.9.0 checkout):
``.venv/bin/python tests/test_rank_k_v010_parity.py``.
"""

from __future__ import annotations

from collections.abc import Callable
import inspect
import pathlib
import sys
from typing import Any
import warnings

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import _rank_k_v010_mocks as mocks  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet  # noqa: E402  pylint: disable=g-import-not-at-top

FIXTURE_PATH = pathlib.Path(__file__).resolve().parent / "data" / "rank_k_v090_parity.npz"

_ARMS4 = ("c", "b1", "b2", "b3")
_ARMS5 = ("c", "b1", "b2", "b3", "b4")
_MU_DENSE4 = np.array([0.0, 0.25, 0.40, 0.48])
_MU_SPARSE4 = np.array([0.0, 0.10, 0.15, 0.17])
_V_DENSE4 = np.array([0.0020, 0.0022, 0.0024, 0.0025])
_V_SPARSE4 = np.array([0.0012, 0.0014, 0.0015, 0.0016])
_CORR = np.array([[1.0, -0.3], [-0.3, 1.0]])


def _supported(func: Callable[..., Any], kwargs: dict[str, Any]) -> dict[str, Any]:
    """Keeps only the keyword arguments that ``func`` accepts (0.9.0 lacks the opt-outs)."""
    params = inspect.signature(func).parameters
    return {k: v for k, v in kwargs.items() if k in params}


def _two_outcome_fits(
    seed: int, pairs: list[tuple[int, int]], arms: tuple[str, ...]
) -> dict[str, dict[tuple[str, str], Any]]:
    rng = np.random.default_rng(seed)
    dense, _ = mocks.simulate_shared_arm_fits(rng, _MU_DENSE4, _V_DENSE4, pairs, arms, n_draws=400)
    sparse, _ = mocks.simulate_shared_arm_fits(rng, _MU_SPARSE4, _V_SPARSE4, pairs, arms, n_draws=400)
    return {"dense": dense, "sparse": sparse}


def _pack_mamo(prefix: str, res: Any, out: dict[str, np.ndarray]) -> None:
    for name in res.outcomes:
        rec = res.by_outcome[name]
        out[f"{prefix}/{name}/tau_mu"] = np.asarray(rec.tau_mu_draws)
        key = res.key_contrast
        out[f"{prefix}/{name}/key_traj"] = np.asarray(rec.contrast_trajectory_draws[key])


def scenario_mamo_bottom_up(opt_out: dict[str, Any]) -> dict[str, np.ndarray]:
    """A = 4 chain, all six pairs, two outcomes plus a bottom-up composite."""
    fits = _two_outcome_fits(11, mocks.all_pairs(4), _ARMS4)
    kwargs = dict(
        arms=_ARMS4,
        prior_mode="bundle_halfnormal",
        prior_scales="in_sample",
        cross_outcome_corr=_CORR,
        composite_outcomes={"total": ("dense", "sparse")},
        composite_mode="bottom_up",
        n_proposals=40_000,
        seed=5,
    )
    kwargs.update(opt_out)
    func = longbet.reconcile_multi_arm_multi_outcome_draws
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = func(fits, **_supported(func, kwargs))
    out: dict[str, np.ndarray] = {}
    _pack_mamo("mamo_bu", res, out)
    return out


def scenario_mamo_coherent_total(opt_out: dict[str, Any]) -> dict[str, np.ndarray]:
    """A = 4 chain, coherent_total with the synthetic composite route (no direct fits)."""
    fits = _two_outcome_fits(12, mocks.all_pairs(4), _ARMS4)
    kwargs = dict(
        arms=_ARMS4,
        prior_mode="bundle_halfnormal",
        prior_scales="in_sample",
        cross_outcome_corr=_CORR,
        composite_outcomes={"total": ("dense", "sparse")},
        composite_mode="coherent_total",
        n_proposals=40_000,
        seed=6,
    )
    kwargs.update(opt_out)
    func = longbet.reconcile_multi_arm_multi_outcome_draws
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = func(fits, **_supported(func, kwargs))
    out: dict[str, np.ndarray] = {}
    _pack_mamo("mamo_ct", res, out)
    return out


def scenario_nested3_bottom_up(opt_out: dict[str, Any]) -> dict[str, np.ndarray]:
    """Legacy nested 3-arm multi-outcome entry point with a bottom-up composite."""
    pairs = [(1, 0), (2, 1), (2, 0)]
    arms3 = ("C", "B", "A")
    rng = np.random.default_rng(13)
    res_by_m: dict[str, dict[tuple[str, str], Any]] = {}
    for m, mu, v in (
        ("dense", _MU_DENSE4[:3], _V_DENSE4[:3]),
        ("sparse", _MU_SPARSE4[:3], _V_SPARSE4[:3]),
    ):
        res_by_m[m], _ = mocks.simulate_shared_arm_fits(rng, mu, v, pairs, arms3, n_draws=400)
    kwargs = dict(
        prior_scales="in_sample",
        prior_mode="bundle_halfnormal",
        cross_outcome_corr=_CORR,
        composite_outcomes={"total": ("dense", "sparse")},
        composite_mode="bottom_up",
        n_proposals=40_000,
        seed=7,
    )
    kwargs.update(opt_out)
    func = longbet.reconcile_nested_3arm_multi_draws
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = func(
            {m: res_by_m[m][("B", "C")] for m in res_by_m},
            {m: res_by_m[m][("A", "B")] for m in res_by_m},
            {m: res_by_m[m][("A", "C")] for m in res_by_m},
            **_supported(func, kwargs),
        )
    out: dict[str, np.ndarray] = {}
    for name in res.outcomes:
        out[f"nested3_bu/{name}/post_att_ac"] = np.asarray(res.post_att_ac[name])
        out[f"nested3_bu/{name}/trt_catt_ac"] = np.asarray(res.trt_catt_ac[name])
    return out


def scenario_slice_a5(opt_out: dict[str, Any]) -> dict[str, np.ndarray]:
    """A = 5 chain, all ten pairs: the ``d = 4`` hit-and-run slice sampler path."""
    rng = np.random.default_rng(14)
    mu = np.array([0.0, 0.22, 0.36, 0.45, 0.50])
    v = np.full(5, 0.0022)
    fits, _ = mocks.simulate_shared_arm_fits(rng, mu, v, mocks.all_pairs(5), _ARMS5, n_draws=400)
    kwargs = dict(arms=_ARMS5, prior_mode="bundle_halfnormal", prior_scale="in_sample", seed=8)
    kwargs.update(opt_out)
    func = longbet.reconcile_multi_arm_draws
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = func(fits, **_supported(func, kwargs))
    return {
        "slice_a5/tau_mu": np.asarray(rec.tau_mu_draws),
        "slice_a5/top_traj": np.asarray(rec.contrast_trajectory_draws[("b4", "c")]),
        "slice_a5/ess": np.asarray(rec.sampler_diagnostics.ess),
        "slice_a5/r_hat": np.asarray(rec.sampler_diagnostics.r_hat),
    }


def scenario_bipartite_star(opt_out: dict[str, Any]) -> dict[str, np.ndarray]:
    """A = 4 vs-control star fits on a chain order (bipartite contrast graph)."""
    rng = np.random.default_rng(15)
    star = [(1, 0), (2, 0), (3, 0)]
    fits, _ = mocks.simulate_shared_arm_fits(rng, _MU_DENSE4, _V_DENSE4, star, _ARMS4, n_draws=400)
    edges = [("b1", "c"), ("b2", "b1"), ("b3", "b2"), ("b2", "c"), ("b3", "c")]
    kwargs = dict(
        arms=_ARMS4,
        topology="poset",
        edges=edges,
        prior_mode="bundle_halfnormal",
        prior_scale="in_sample",
        n_proposals=40_000,
        seed=9,
    )
    kwargs.update(opt_out)
    func = longbet.reconcile_multi_arm_draws
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = func(fits, **_supported(func, kwargs))
    return {
        "star/tau_mu": np.asarray(rec.tau_mu_draws),
        "star/arm_variances": np.asarray(rec.arm_variances, dtype=np.float64),
        "star/sigma_lik": np.asarray(rec.sigma_lik),
    }


def scenario_single_a4(opt_out: dict[str, Any]) -> dict[str, np.ndarray]:
    """A = 4 chain, all six pairs, single outcome (SIR path, in-sample scale)."""
    rng = np.random.default_rng(16)
    fits, _ = mocks.simulate_shared_arm_fits(
        rng, _MU_DENSE4, _V_DENSE4, mocks.all_pairs(4), _ARMS4, n_draws=400
    )
    kwargs = dict(
        arms=_ARMS4,
        prior_mode="bundle_halfnormal",
        prior_scale="in_sample",
        n_proposals=40_000,
        seed=10,
    )
    kwargs.update(opt_out)
    func = longbet.reconcile_multi_arm_draws
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = func(fits, **_supported(func, kwargs))
    return {
        "single_a4/tau_mu": np.asarray(rec.tau_mu_draws),
        "single_a4/mu_hat": np.asarray(rec.mu_hat),
        "single_a4/sigma_lik": np.asarray(rec.sigma_lik),
    }


SCENARIOS: dict[str, Callable[[dict[str, Any]], dict[str, np.ndarray]]] = {
    "mamo_bottom_up": scenario_mamo_bottom_up,
    "mamo_coherent_total": scenario_mamo_coherent_total,
    "nested3_bottom_up": scenario_nested3_bottom_up,
    "slice_a5": scenario_slice_a5,
    "bipartite_star": scenario_bipartite_star,
    "single_a4": scenario_single_a4,
}

# Keyword arguments that restore the 0.9.0 behavior in 0.10 (ignored by 0.9.0 itself).
V090_OPT_OUTS: dict[str, dict[str, Any]] = {
    "mamo_bottom_up": {},
    "mamo_coherent_total": {},
    "nested3_bottom_up": {},
    "slice_a5": {},
    "bipartite_star": {},
    "single_a4": {},
}


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_v090_parity_bit_for_bit(name: str) -> None:
    stored = np.load(FIXTURE_PATH)
    got = SCENARIOS[name](V090_OPT_OUTS[name])
    assert got, name
    for key, arr in got.items():
        assert key in stored.files, key
        np.testing.assert_array_equal(arr, stored[key], err_msg=key)


def main() -> None:
    """Regenerates the stored fixture (run on a 0.9.0 checkout only)."""
    if longbet.__version__ != "0.9.0":
        raise SystemExit(f"Refusing to regenerate the 0.9.0 fixture on {longbet.__version__}")
    out: dict[str, np.ndarray] = {}
    for name, func in SCENARIOS.items():
        out.update(func({}))
    FIXTURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(FIXTURE_PATH, **out)
    print(f"wrote {len(out)} arrays to {FIXTURE_PATH}")


if __name__ == "__main__":
    main()
