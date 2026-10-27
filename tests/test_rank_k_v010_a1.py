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

"""v0.10 item A1: ``coherent_total`` default composite mode and the ``bottom_up`` warning."""

from __future__ import annotations

import math
import pathlib
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import _rank_k_v010_mocks as mocks  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet  # noqa: E402  pylint: disable=g-import-not-at-top
import longbet._rank_k_multi as rkm  # noqa: E402  pylint: disable=g-import-not-at-top

_ARMS4 = ("c", "b1", "b2", "b3")
_ARMS3 = ("C", "B", "A")
_CORR = np.array([[1.0, -0.3], [-0.3, 1.0]])


def _fits(
    seed: int, mu_d: np.ndarray, mu_s: np.ndarray, v: float, arms: tuple[str, ...]
) -> dict[str, dict[tuple[str, str], longbet.RankKLongBetResult]]:
    rng = np.random.default_rng(seed)
    pairs = mocks.all_pairs(len(arms))
    var = np.full(len(arms), v)
    dense, _ = mocks.simulate_shared_arm_fits(rng, mu_d, var, pairs, arms, n_draws=400)
    sparse, _ = mocks.simulate_shared_arm_fits(rng, mu_s, var, pairs, arms, n_draws=400)
    return {"dense": dense, "sparse": sparse}


def _bottom_up_warnings(caught: list[warnings.WarningMessage]) -> list[str]:
    return [
        str(w.message)
        for w in caught
        if issubclass(w.category, longbet.RankKConfigWarning)
        and "composite_mode='bottom_up'" in str(w.message)
    ]


def test_floor_inflation_formula_matches_sqrt_m_for_independent_equal_components() -> None:
    diag = rkm._composite_null_floor_diagnostic([0.0, 0.0, 0.0], [0.1, 0.1, 0.1], np.eye(3))
    assert diag["expected_floor_inflation"] == pytest.approx(math.sqrt(3.0))
    assert diag["p_nonpositive"] == pytest.approx(0.5)
    corr = np.array([[1.0, 0.5], [0.5, 1.0]])
    diag2 = rkm._composite_null_floor_diagnostic([0.2, 0.1], [0.1, 0.2], corr)
    sd = math.sqrt(0.01 + 0.04 + 2 * 0.5 * 0.02)
    assert diag2["sd"] == pytest.approx(sd)
    assert diag2["expected_floor_inflation"] == pytest.approx(0.3 / sd)


def test_multi_arm_default_composite_mode_is_coherent_total() -> None:
    fits = _fits(1, np.array([0.0, 0.2, 0.3, 0.4]), np.array([0.0, 0.1, 0.15, 0.2]), 0.002, _ARMS4)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = longbet.reconcile_multi_arm_multi_outcome_draws(
            fits,
            arms=_ARMS4,
            cross_outcome_corr=_CORR,
            composite_outcomes={"total": ("dense", "sparse")},
            n_proposals=20_000,
            seed=3,
        )
    assert res.composite_mode == "coherent_total"
    assert res.max_abs_composite_discrepancy["total"] < 1e-12
    assert "total" in res.composite_floor_diagnostics


def test_nested3_default_composite_mode_is_coherent_total() -> None:
    fits = _fits(2, np.array([0.0, 0.2, 0.35]), np.array([0.0, 0.1, 0.15]), 0.002, _ARMS3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = longbet.reconcile_nested_3arm_multi_draws(
            {m: fits[m][("B", "C")] for m in fits},
            {m: fits[m][("A", "B")] for m in fits},
            {m: fits[m][("A", "C")] for m in fits},
            cross_outcome_corr=_CORR,
            composite_outcomes={"total": ("dense", "sparse")},
            n_proposals=20_000,
            seed=4,
        )
    assert res.composite_mode == "coherent_total"
    assert res.max_abs_composite_discrepancy["total"] < 1e-12


def test_default_falls_back_to_bottom_up_for_overlapping_composites() -> None:
    rng = np.random.default_rng(5)
    pairs = mocks.all_pairs(4)
    var = np.full(4, 0.002)
    fits = {
        m: mocks.simulate_shared_arm_fits(rng, np.array([0.0, 0.2, 0.3, 0.4]), var, pairs, _ARMS4, n_draws=300)[0]
        for m in ("x", "y", "z")
    }
    with pytest.warns(longbet.RankKConfigWarning, match="disjoint component sets"):
        res = longbet.reconcile_multi_arm_multi_outcome_draws(
            fits,
            arms=_ARMS4,
            composite_outcomes={"xy": ("x", "y"), "yz": ("y", "z")},
            n_proposals=10_000,
            seed=6,
        )
    assert res.composite_mode == "bottom_up"


def test_default_is_bottom_up_without_truncation() -> None:
    assert rkm._resolve_composite_mode(None, {"t": ("a", "b")}, truncated=False) == "bottom_up"
    assert rkm._resolve_composite_mode(None, {"t": ("a", "b")}, truncated=True) == "coherent_total"
    with pytest.raises(ValueError, match="composite_mode must be one of"):
        rkm._resolve_composite_mode("sum", None, truncated=True)


@pytest.mark.parametrize("entry", ["multi_arm", "nested3"])
def test_bottom_up_warning_fires_on_null_and_not_on_active_dgp(entry: str) -> None:
    n_arms = 4 if entry == "multi_arm" else 3
    arms = _ARMS4 if entry == "multi_arm" else _ARMS3
    null_mu = np.zeros(n_arms)
    active_d = np.array([0.0, 0.25, 0.40, 0.48])[:n_arms]
    active_s = np.array([0.0, 0.10, 0.15, 0.17])[:n_arms]

    def _run(fits: dict[str, dict[tuple[str, str], longbet.RankKLongBetResult]]) -> tuple[list[str], object]:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            if entry == "multi_arm":
                res = longbet.reconcile_multi_arm_multi_outcome_draws(
                    fits,
                    arms=arms,
                    composite_outcomes={"total": ("dense", "sparse")},
                    composite_mode="bottom_up",
                    n_proposals=10_000,
                    seed=7,
                )
            else:
                res = longbet.reconcile_nested_3arm_multi_draws(
                    {m: fits[m][("B", "C")] for m in fits},
                    {m: fits[m][("A", "B")] for m in fits},
                    {m: fits[m][("A", "C")] for m in fits},
                    composite_outcomes={"total": ("dense", "sparse")},
                    composite_mode="bottom_up",
                    n_proposals=10_000,
                    seed=7,
                )
        return _bottom_up_warnings(caught), res

    msgs_null, res_null = _run(_fits(8, null_mu, null_mu, 0.002, arms))
    assert len(msgs_null) == 1, msgs_null
    infl = res_null.composite_floor_diagnostics["total"]["expected_floor_inflation"]
    assert f"{infl:.3f}x" in msgs_null[0]
    assert res_null.composite_floor_diagnostics["total"]["p_nonpositive"] > 0.10

    msgs_active, res_active = _run(_fits(9, active_d, active_s, 0.002, arms))
    assert msgs_active == []
    assert res_active.composite_floor_diagnostics["total"]["p_nonpositive"] < 0.10


def test_synthetic_composite_matches_direct_composite_fits() -> None:
    """Without direct composite fits, coherent_total samples the composite from component moments.

    On shared-arm panels whose direct composite fit has the moments implied by its
    components, the synthetic route and the direct route agree up to Monte Carlo error
    panel by panel (the R = 500 Tier A comparison is reported in the v0.10 notes).
    """
    pairs = mocks.all_pairs(4)
    v_d = np.array([0.0020, 0.0022, 0.0024, 0.0025])
    v_s = np.array([0.0012, 0.0014, 0.0015, 0.0016])
    mu_d = np.array([0.0, 0.25, 0.40, 0.48])
    mu_s = np.array([0.0, 0.10, 0.15, 0.17])
    for panel in range(4):
        rng = np.random.default_rng(100 + panel)
        z_d = mu_d + np.sqrt(v_d) * rng.standard_normal(4)
        z_s = mu_s + np.sqrt(v_s) * rng.standard_normal(4)
        dense, y_d = mocks.simulate_shared_arm_fits(rng, mu_d, v_d, pairs, _ARMS4, n_draws=1000, z_bar=z_d)
        sparse, y_s = mocks.simulate_shared_arm_fits(rng, mu_s, v_s, pairs, _ARMS4, n_draws=1000, z_bar=z_s)
        v_tot = v_d + v_s
        total = {
            (_ARMS4[j], _ARMS4[l]): mocks.make_mock_fit(
                att_mean=float(y_d[p] + y_s[p]),
                att_sd=math.sqrt(float(v_tot[j] + v_tot[l])),
                n_draws=1000,
                n_units=8,
                t_post=6,
                seed=1000 + p,
                exact_moments=True,
            )
            for p, (j, l) in enumerate(pairs)
        }
        fits = {"dense": dense, "sparse": sparse}
        kwargs = dict(
            arms=_ARMS4,
            composite_outcomes={"total": ("dense", "sparse")},
            composite_mode="coherent_total",
            n_proposals=60_000,
            seed=11 + panel,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            synth = longbet.reconcile_multi_arm_multi_outcome_draws(fits, **kwargs)
            direct = longbet.reconcile_multi_arm_multi_outcome_draws(
                fits, composite_results={"total": total}, **kwargs
            )
        for pair in [("b3", "c"), ("b1", "c"), ("b3", "b1")]:
            d_syn = synth.by_outcome["total"].get_contrast_result(*pair).mode1_draws
            d_dir = direct.by_outcome["total"].get_contrast_result(*pair).mode1_draws
            sd = float(np.std(d_dir, ddof=1))
            assert abs(float(np.mean(d_syn)) - float(np.mean(d_dir))) <= 0.25 * sd, (panel, pair)
            assert float(np.std(d_syn, ddof=1)) == pytest.approx(sd, rel=0.2), (panel, pair)
