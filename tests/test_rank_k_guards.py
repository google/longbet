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

"""Guards for the Rank-K options: variant defaults, input validation and warnings.

No model is fitted here. The tests pin the variant-dependent defaults of
``RankKLongBetConfig`` and its ``ValueError``s (including the type checks), the
input validation and the ``RankKConfigWarning`` policy of
``fit_rank_k_longbet`` (through ``_validate_rank_k_inputs``, which runs before
any sampling), the dispatch of the default fit and the placebo audit's feature
rebuild and handling of ``features`` (with the samplers and the feature builder
stubbed out), and the robust kriging operator of the JAX sampler.
"""

from __future__ import annotations

import dataclasses
import os
from pathlib import Path
import re
import types
from typing import Any
import warnings

import jax
import numpy as np
import pytest

import longbet
from longbet import (
    PlaceboAuditResult,
    RankKConfigWarning,
    RankKLongBetConfig,
    fit_rank_k_longbet,
    rank_k_placebo_audit,
)
from longbet import _rank_k as rank_k_module
from longbet._rank_k import _PLUG_IN_DEFAULTS, _PROPAGATED_DEFAULTS, _validate_rank_k_inputs
from longbet._rank_k_jax import (
    _kriging_operator,
    build_spatial_kernel_cache,
    great_circle_distance_km,
    matern32_corr,
)
from longbet._x64 import enable_x64

N_UNITS, N_PERIODS = 20, 8
CLUSTERS = np.repeat(np.arange(5), 4)

# Validated budgets (companion paper): plug-in runs and the propagated benchmark runs.
PLUG_IN = {
    "num_chains": 1,
    "num_burnin": 200,
    "num_trees_mu": 75,
    "num_trees_tau_mode1": 35,
    "num_trees_tau_higher": 20,
}
PROPAGATED = {
    "num_chains": 2,
    "num_burnin": 100,
    "num_trees_mu": 50,
    "num_trees_tau_mode1": 30,
    "num_trees_tau_higher": 15,
}

EXPECTED_FIELDS = (
    "rank_k",
    "use_forecast",
    "engine",
    "uncertainty_propagation",
    "n_forecast_draws",
    "max_forecast_draws",
    "between_var_threshold",
    "num_chains",
    "draws_per_component",
    "num_burnin",
    "warm_burnin",
    "num_gfr",
    "num_mcmc",
    "num_trees_mu",
    "num_trees_tau_mode1",
    "num_trees_tau_higher",
    "theta_f_clip",
    "theta_prior_mean",
    "theta_prior_sd",
    "tau_prior_mean",
    "tau_prior_sd",
    "tau_clip",
    "catt_floor",
    "use_factor_loadings",
    "n_factor_loadings",
    "use_spatial_gp",
    "lengthscale_grid_km",
    "baseline_cells",
    "ridge_eps",
    "num_threads",
    "random_seed",
    # Remediation options (methodology review, Tier 1)
    "tau_prior_family",
    "short_cohort_policy",
    "thin",
    "sigma_nu_total",
    "sigma_nu_dynamic_total",
    "alpha_nu",
    "beta_nu",
    "nu_working_zeta_weight",
    "learn_sigma_nu",
    "sigma_nu_prior_scale",
    "treated_extra_variance",
    "basis_spec",
    "keep_unprojected_draws",
)


def _panel(seed: int = 0) -> types.SimpleNamespace:
    """Staggered panel: 6 units adopt at t = 4, 4 at t = 6 and 10 are never treated (T = 8)."""
    rng = np.random.default_rng(seed)
    cohorts = np.full(N_UNITS, N_PERIODS, dtype=np.int64)
    cohorts[:6] = 4
    cohorts[6:10] = 6
    z = (np.arange(N_PERIODS)[None, :] >= cohorts[:, None]).astype(np.int32)
    return types.SimpleNamespace(
        y=rng.normal(size=(N_UNITS, N_PERIODS)),
        z=z,
        x=rng.normal(size=(N_UNITS, 2)),
        cohorts=cohorts,
    )


def _features(**overrides: Any) -> types.SimpleNamespace:
    """Stand-in for the output of ``build_causal_features`` on the ``_panel()`` geometry."""
    fields = {
        "anchor_t": 4,
        "point_forecast": np.full((N_UNITS, N_PERIODS), 0.5),
        "point_forecast_draws": None,
        "unit_features": np.zeros((N_UNITS, 2)),
        "ps": np.full(N_UNITS, 0.5),
        "pre_anchor_loadings": np.zeros((N_UNITS, 3)),
    }
    fields.update(overrides)
    return types.SimpleNamespace(**fields)


def _coords(seed: int = 1) -> np.ndarray:
    """(latitude, longitude) of ``N_UNITS`` points in the contiguous US."""
    rng = np.random.default_rng(seed)
    return np.column_stack([rng.uniform(25.0, 49.0, N_UNITS), rng.uniform(-124.0, -67.0, N_UNITS)])


def _poke(arr: Any, index: Any, value: float) -> np.ndarray:
    """Return a float copy of ``arr`` with ``arr[index] = value``."""
    out = np.array(arr, dtype=np.float64)
    out[index] = value
    return out


def _validate(cfg: RankKLongBetConfig | None = None, **overrides: Any) -> None:
    """Run the input validation of ``fit_rank_k_longbet`` on ``_panel()`` with stub features."""
    panel = _panel()
    inputs = {"y": panel.y, "z": panel.z, "x": panel.x, "features": _features()}
    inputs.update(overrides)
    _validate_rank_k_inputs(
        inputs.pop("y"),
        inputs.pop("z"),
        inputs.pop("x"),
        cfg=RankKLongBetConfig() if cfg is None else cfg,
        **inputs,
    )


def _caught_warnings(cfg: RankKLongBetConfig | None = None, **overrides: Any) -> list[Any]:
    """Return every warning that ``_validate`` emits."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _validate(cfg, **overrides)
    return list(caught)


# ---- Variants and their defaults ------------------------------------------------------------


def test_default_config_is_the_propagated_variant():
    cfg = RankKLongBetConfig()
    assert cfg.uncertainty_propagation is True
    assert cfg.engine == "jax"
    assert {name: getattr(cfg, name) for name in PROPAGATED} == PROPAGATED
    assert (cfg.n_forecast_draws, cfg.max_forecast_draws, cfg.draws_per_component) == (8, None, 25)
    assert (cfg.use_spatial_gp, cfg.use_factor_loadings, cfg.n_factor_loadings) == (True, True, 3)
    assert cfg.baseline_cells == "untreated"
    assert (cfg.theta_prior_mean, cfg.theta_prior_sd, cfg.theta_f_clip) == (1.0, 0.5, (0.0, 1.5))
    assert cfg.lengthscale_grid_km == (100.0, 200.0, 400.0, 800.0, 1600.0)
    assert cfg.between_var_threshold == 0.20
    assert (cfg.rank_k, cfg.use_forecast) == (2, True)
    assert cfg.warm_burnin is None
    # Posterior draws: num_chains x n_forecast_draws x draws_per_component.
    assert cfg.num_chains * cfg.n_forecast_draws * cfg.draws_per_component == 400


def test_plug_in_variant_defaults():
    cfg = RankKLongBetConfig(uncertainty_propagation=False)
    assert cfg.uncertainty_propagation is False
    assert cfg.engine == "stochtree"
    assert {name: getattr(cfg, name) for name in PLUG_IN} == PLUG_IN
    assert (cfg.num_gfr, cfg.num_mcmc) == (10, 600)
    assert cfg.theta_f_clip == (0.0, 1.5)
    assert (cfg.ridge_eps, cfg.num_threads, cfg.random_seed) == (1e-4, 1, 2026)
    assert (cfg.rank_k, cfg.use_forecast) == (2, True)
    assert cfg.warm_burnin is None
    assert cfg.max_forecast_draws is None


def test_variant_budget_tables():
    assert _PLUG_IN_DEFAULTS == PLUG_IN
    assert _PROPAGATED_DEFAULTS == PROPAGATED


@pytest.mark.parametrize(
    ("kwargs", "propagated"),
    [
        pytest.param({}, True, id="default"),
        pytest.param({"uncertainty_propagation": False}, False, id="propagation-false"),
        pytest.param({"uncertainty_propagation": True}, True, id="propagation-true"),
        pytest.param({"engine": "stochtree"}, False, id="engine-stochtree"),
        pytest.param({"engine": "jax"}, True, id="engine-jax"),
        pytest.param({"engine": "stochtree", "uncertainty_propagation": False}, False, id="consistent-plug-in"),
        pytest.param({"engine": "jax", "uncertainty_propagation": True}, True, id="consistent-propagated"),
        pytest.param({"uncertainty_propagation": np.bool_(True)}, True, id="numpy-bool"),
        pytest.param({"uncertainty_propagation": 0}, False, id="integer-zero"),
    ],
)
def test_engine_is_derived_from_uncertainty_propagation(kwargs, propagated):
    cfg = RankKLongBetConfig(**kwargs)
    assert cfg.uncertainty_propagation is propagated
    assert cfg.engine == ("jax" if propagated else "stochtree")
    budget = PROPAGATED if propagated else PLUG_IN
    assert {name: getattr(cfg, name) for name in budget} == budget
    reference = RankKLongBetConfig(uncertainty_propagation=propagated)
    assert cfg == reference
    assert hash(cfg) == hash(reference)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        pytest.param(
            {"engine": "jax", "uncertainty_propagation": False},
            "engine='jax' with uncertainty_propagation=False is not a supported variant",
            id="jax-without-propagation",
        ),
        pytest.param(
            {"engine": "stochtree", "uncertainty_propagation": True},
            "engine='stochtree' with uncertainty_propagation=True is not implemented",
            id="stochtree-with-propagation",
        ),
    ],
)
def test_contradictory_variant_settings_raise(kwargs, match):
    with pytest.raises(ValueError, match=re.escape(match)):
        RankKLongBetConfig(**kwargs)


def test_replace_keeps_the_variant_and_cannot_switch_it():
    cfg = dataclasses.replace(RankKLongBetConfig(), rank_k=3)
    assert (cfg.engine, cfg.rank_k) == ("jax", 3)
    assert {name: getattr(cfg, name) for name in PROPAGATED} == PROPAGATED
    # max_forecast_draws stays None, so changing n_forecast_draws never turns the expansion on.
    assert dataclasses.replace(cfg, n_forecast_draws=4).max_forecast_draws is None
    assert dataclasses.replace(cfg, n_forecast_draws=16).max_forecast_draws is None
    switch = re.escape("build a new RankKLongBetConfig instead of calling dataclasses.replace")
    with pytest.raises(ValueError, match=switch):
        dataclasses.replace(RankKLongBetConfig(), uncertainty_propagation=False)
    with pytest.raises(ValueError, match=switch):
        dataclasses.replace(RankKLongBetConfig(uncertainty_propagation=False), uncertainty_propagation=True)
    # Replacing engine=None as well switches the engine but keeps the stored budget (the fit warns).
    trap = dataclasses.replace(RankKLongBetConfig(), engine=None, uncertainty_propagation=False)
    assert trap.engine == "stochtree"
    assert {name: getattr(trap, name) for name in PROPAGATED} == PROPAGATED


def test_explicit_budget_values_are_kept():
    cfg = RankKLongBetConfig(num_chains=4, num_burnin=300, draws_per_component=100)
    assert (cfg.num_chains, cfg.num_burnin, cfg.draws_per_component) == (4, 300, 100)
    assert (cfg.num_trees_mu, cfg.num_trees_tau_mode1, cfg.num_trees_tau_higher) == (50, 30, 15)
    assert cfg.num_chains * cfg.n_forecast_draws * cfg.draws_per_component == 3200
    cfg = RankKLongBetConfig(uncertainty_propagation=False, num_burnin=50, num_trees_mu=10)
    assert (cfg.engine, cfg.num_burnin, cfg.num_trees_mu, cfg.num_trees_tau_mode1) == ("stochtree", 50, 10, 35)
    # 0 is a value, not the None sentinel.
    assert RankKLongBetConfig(num_burnin=0).num_burnin == 0


def test_max_forecast_draws_none_means_no_expansion():
    cfg = RankKLongBetConfig(uncertainty_propagation=True, n_forecast_draws=16)
    assert cfg.max_forecast_draws is None
    assert rank_k_module._max_forecast_draws(cfg) == 16
    cfg = RankKLongBetConfig(uncertainty_propagation=True, max_forecast_draws=32)
    assert (cfg.max_forecast_draws, rank_k_module._max_forecast_draws(cfg)) == (32, 32)


def test_config_fields_are_the_documented_options():
    names = tuple(f.name for f in dataclasses.fields(RankKLongBetConfig))
    assert names == EXPECTED_FIELDS
    for name in names:
        assert f"  {name}: " in RankKLongBetConfig.__doc__, name
    with pytest.raises(TypeError):
        RankKLongBetConfig(use_forests=False)  # removed in 0.5.0
    assert not [f.name for f in dataclasses.fields(PlaceboAuditResult) if "moran" in f.name]


def test_public_api_and_packaging():
    assert longbet.__version__ == "0.8.0"
    assert "RankKConfigWarning" in longbet.__all__
    assert issubclass(RankKConfigWarning, UserWarning)
    tomllib = pytest.importorskip("tomllib")
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    project = tomllib.loads(pyproject.read_text())["project"]
    assert project["version"] == longbet.__version__
    # The plug-in and hybrid variants fit stochtree BCFs, so stochtree is a core dependency.
    assert any(dep.replace(" ", "").startswith("stochtree>=") for dep in project["dependencies"])


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        pytest.param({"rank_k": 0}, "rank_k must be >= 1", id="rank-k"),
        pytest.param({"engine": "numpy"}, "engine must be None, 'jax' or 'stochtree'", id="engine"),
        pytest.param(
            {"uncertainty_propagation": "yes"},
            "uncertainty_propagation must be None, True or False",
            id="propagation-string",
        ),
        pytest.param(
            {"uncertainty_propagation": 2},
            "uncertainty_propagation must be None, True or False",
            id="propagation-int",
        ),
        pytest.param(
            {"baseline_cells": "treated"}, "baseline_cells must be 'untreated' or 'all'", id="baseline-cells"
        ),
        pytest.param({"n_forecast_draws": -1}, "n_forecast_draws must be >= 0", id="n-forecast-draws"),
        pytest.param(
            {"uncertainty_propagation": True, "max_forecast_draws": 4},
            "max_forecast_draws (4) must be >= n_forecast_draws (8)",
            id="max-forecast-draws",
        ),
        pytest.param({"num_chains": 0}, "num_chains must be >= 1", id="num-chains"),
        pytest.param({"num_burnin": -1}, "num_burnin must be >= 0", id="num-burnin"),
        pytest.param({"warm_burnin": -1}, "warm_burnin must be >= 0", id="warm-burnin"),
        pytest.param({"draws_per_component": 1}, "draws_per_component must be >= 2", id="draws-per-component"),
        pytest.param({"num_mcmc": 1}, "num_mcmc must be >= 2", id="num-mcmc"),
        pytest.param({"num_trees_mu": 0}, "num_trees_mu must be >= 1", id="trees-mu"),
        pytest.param({"num_trees_tau_mode1": 0}, "num_trees_tau_mode1 must be >= 1", id="trees-tau-mode1"),
        pytest.param({"num_trees_tau_higher": 0}, "num_trees_tau_higher must be >= 1", id="trees-tau-higher"),
        pytest.param({"num_gfr": -1}, "num_gfr must be >= 0", id="num-gfr"),
        pytest.param({"theta_f_clip": (1.0, 0.5)}, "theta_f_clip must be a (low, high) pair", id="clip-order"),
        pytest.param({"theta_f_clip": (0.0,)}, "theta_f_clip must be a (low, high) pair", id="clip-length"),
        pytest.param({"theta_prior_sd": 0.0}, "theta_prior_sd must be > 0", id="prior-sd"),
        pytest.param({"tau_prior_sd": 0.0}, "tau_prior_sd must be > 0 or None", id="tau-prior-sd"),
        pytest.param({"tau_clip": (1.0, 0.0)}, "tau_clip must be a (low, high) pair", id="tau-clip-order"),
        pytest.param({"catt_floor": float("nan")}, "catt_floor must be a finite number", id="catt-floor-nan"),
        pytest.param({"n_factor_loadings": -1}, "n_factor_loadings must be >= 0", id="n-factor-loadings"),
        pytest.param(
            {"lengthscale_grid_km": ()},
            "lengthscale_grid_km must be a non-empty sequence of positive, finite ranges",
            id="grid-empty",
        ),
        pytest.param(
            {"lengthscale_grid_km": (100.0, -5.0)},
            "lengthscale_grid_km must be a non-empty sequence of positive, finite ranges",
            id="grid-negative",
        ),
        pytest.param({"between_var_threshold": 1.5}, "between_var_threshold must lie in [0, 1]", id="threshold-high"),
        pytest.param({"between_var_threshold": -0.1}, "between_var_threshold must lie in [0, 1]", id="threshold-low"),
        pytest.param({"ridge_eps": -1e-4}, "ridge_eps must be >= 0", id="ridge-eps"),
        pytest.param(
            {"num_threads": 0}, "num_threads must be >= 1, or -1 for all available threads; got 0", id="threads-zero"
        ),
        pytest.param({"num_threads": -2}, "num_threads must be >= -1", id="threads-negative"),
        pytest.param({"random_seed": -1}, "random_seed must be >= 0", id="seed-negative"),
        # Types: counts must be whole numbers (int() would truncate them), flags booleans, reals finite.
        pytest.param({"rank_k": 2.5}, "rank_k must be an integer >= 1; got 2.5", id="rank-k-fraction"),
        pytest.param({"rank_k": None}, "rank_k must be an integer >= 1; got None", id="rank-k-none"),
        pytest.param({"num_trees_mu": True}, "num_trees_mu must be an integer >= 1; got True", id="trees-bool"),
        pytest.param({"num_mcmc": float("nan")}, "num_mcmc must be an integer >= 2; got nan", id="num-mcmc-nan"),
        pytest.param(
            {"num_burnin": float("inf")}, "num_burnin must be an integer >= 0; got inf", id="num-burnin-inf"
        ),
        pytest.param(
            {"draws_per_component": "25"}, "draws_per_component must be an integer >= 2; got '25'", id="draws-string"
        ),
        pytest.param(
            {"n_forecast_draws": 2.7}, "n_forecast_draws must be an integer >= 0; got 2.7", id="draws-fraction"
        ),
        pytest.param(
            {"uncertainty_propagation": True, "max_forecast_draws": 8.5},
            "max_forecast_draws must be an integer >= 0; got 8.5",
            id="max-draws-fraction",
        ),
        pytest.param({"random_seed": 1.5}, "random_seed must be an integer >= 0; got 1.5", id="seed-fraction"),
        pytest.param({"use_forecast": "False"}, "use_forecast must be True or False; got 'False'", id="flag-string"),
        pytest.param({"use_spatial_gp": None}, "use_spatial_gp must be True or False; got None", id="flag-none"),
        pytest.param({"use_factor_loadings": 2}, "use_factor_loadings must be True or False; got 2", id="flag-int"),
        pytest.param({"theta_f_clip": "01"}, "theta_f_clip must be a (low, high) pair", id="clip-string"),
        pytest.param({"theta_f_clip": (0.0, float("nan"))}, "theta_f_clip must be a (low, high) pair", id="clip-nan"),
        pytest.param(
            {"lengthscale_grid_km": "100"},
            "lengthscale_grid_km must be a non-empty sequence of positive, finite ranges",
            id="grid-string",
        ),
        pytest.param(
            {"lengthscale_grid_km": (100.0, float("inf"))},
            "lengthscale_grid_km must be a non-empty sequence of positive, finite ranges",
            id="grid-inf",
        ),
        pytest.param(
            {"theta_prior_mean": float("inf")},
            "theta_prior_mean must be a finite number; got inf",
            id="prior-mean-inf",
        ),
        pytest.param({"theta_prior_sd": "0.5"}, "theta_prior_sd must be a finite number", id="prior-sd-string"),
        pytest.param(
            {"between_var_threshold": float("nan")},
            "between_var_threshold must be a finite number; got nan",
            id="threshold-nan",
        ),
    ],
)
def test_invalid_config_values_raise(kwargs, match):
    with pytest.raises(ValueError, match=re.escape(match)):
        RankKLongBetConfig(**kwargs)


def test_numeric_types_are_normalised():
    cfg = RankKLongBetConfig(
        rank_k=np.int64(3),
        num_mcmc=600.0,
        num_burnin=np.float64(150.0),
        use_forecast=np.bool_(True),
        use_spatial_gp=1,
        theta_f_clip=[0, 2],
        lengthscale_grid_km=np.array([100, 200]),
        between_var_threshold=1,
        ridge_eps=0,
    )
    assert (cfg.rank_k, cfg.num_mcmc, cfg.num_burnin) == (3, 600, 150)
    assert all(type(v) is int for v in (cfg.rank_k, cfg.num_mcmc, cfg.num_burnin))
    assert cfg.use_forecast is True and cfg.use_spatial_gp is True
    assert (cfg.theta_f_clip, cfg.lengthscale_grid_km) == ((0.0, 2.0), (100.0, 200.0))
    assert type(cfg.between_var_threshold) is float and type(cfg.ridge_eps) is float
    reference = RankKLongBetConfig(
        rank_k=3,
        num_burnin=150,
        theta_f_clip=(0.0, 2.0),
        lengthscale_grid_km=(100.0, 200.0),
        between_var_threshold=1.0,
        ridge_eps=0.0,
    )
    assert cfg == reference
    assert hash(cfg) == hash(reference)  # list and array inputs are stored as tuples


@pytest.mark.parametrize(
    "kwargs",
    [
        pytest.param({"n_forecast_draws": 0}, id="no-forecast-draws"),
        pytest.param({"uncertainty_propagation": True, "n_forecast_draws": 1}, id="one-forecast-draw"),
        pytest.param({"between_var_threshold": 0.0}, id="threshold-zero"),
        pytest.param({"between_var_threshold": 1.0}, id="threshold-one"),
        pytest.param({"num_threads": -1}, id="all-threads"),
        pytest.param({"random_seed": 0}, id="seed-zero"),
        pytest.param({"ridge_eps": 0.0}, id="no-ridge"),
        pytest.param({"num_gfr": 0, "num_burnin": 0}, id="no-gfr-no-burn-in"),
        pytest.param({"n_factor_loadings": 0}, id="no-loadings"),
        pytest.param({"theta_f_clip": (1.0, 1.0)}, id="fixed-slope"),
    ],
)
def test_boundary_config_values_are_accepted(kwargs):
    cfg = RankKLongBetConfig(**kwargs)
    for name, value in kwargs.items():
        assert getattr(cfg, name) == value, name


# ---- Input validation (errors are raised in both variants, before any sampling) -------------

_P = _panel()
_INVALID_INPUTS = [
    pytest.param({"y": _P.y[:, 0]}, "y must be a 2-D (N, T) outcome panel with T >= 2", id="y-1d"),
    pytest.param({"y": _poke(_P.y, (0, 0), np.nan)}, "y contains 1 non-finite value(s)", id="y-nan"),
    pytest.param({"z": _P.z[:, :-1]}, "z shape (20, 7) must match y shape (20, 8)", id="z-shape"),
    pytest.param({"z": 2 * _P.z}, "z must contain only 0/1 treatment indicators", id="z-not-binary"),
    pytest.param({"z": _poke(_P.z, (0, 6), 0.0)}, "z must be absorbing", id="z-not-absorbing"),
    pytest.param({"x": _P.x[:-1]}, "x must have one row per unit (N=20)", id="x-rows"),
    pytest.param({"x": _poke(_P.x, (0, 0), np.inf)}, "x contains 1 non-finite value(s)", id="x-inf"),
    pytest.param(
        {"cohorts": _P.cohorts[:-1]},
        "cohorts must be a 1-D array with one first-adoption period per unit (N=20)",
        id="cohorts-shape",
    ),
    pytest.param(
        {"cohorts": _poke(_P.cohorts, 12, N_PERIODS + 1)},
        "cohorts must be integer first-adoption periods in [0, T] = [0, 8]",
        id="cohorts-range",
    ),
    pytest.param(
        {"cohorts": _poke(_P.cohorts, 0, 4.5)},
        "cohorts must be integer first-adoption periods",
        id="cohorts-fraction",
    ),
    pytest.param(
        {"cohorts": _poke(_P.cohorts, 0, 5)}, "cohorts is inconsistent with z for 1 unit(s)", id="cohorts-mismatch"
    ),
    pytest.param(
        {"z": np.zeros_like(_P.z)},
        "Panel must contain both ever-treated and never-treated units",
        id="no-treated-units",
    ),
    pytest.param(
        {"z": _poke(_P.z, (slice(10, None), N_PERIODS - 1), 1.0)},
        "Panel must contain both ever-treated and never-treated units",
        id="no-never-treated-units",
    ),
    pytest.param(
        {"z": _poke(_P.z, (0, slice(None)), 1.0)},
        "Every ever-treated unit needs at least one pre-adoption period",
        id="treated-from-period-0",
    ),
    pytest.param({"t_0_panel": 5}, "t_0_panel must be a single integer in [1, 4]", id="t0-after-adoption"),
    pytest.param({"t_0_panel": 0}, "t_0_panel must be a single integer in [1, 4]", id="t0-zero"),
    pytest.param({"t_0_panel": 2.5}, "t_0_panel must be a single integer in [1, 4]", id="t0-fraction"),
    pytest.param({"t_0_panel": np.array([2, 3])}, "t_0_panel must be a single integer in [1, 4]", id="t0-vector"),
    pytest.param(
        {"features": _features(anchor_t=5)},
        "features.anchor_t = 5 is after the first adoption period 4",
        id="features-anchor",
    ),
    pytest.param(
        {"features": _features(point_forecast=np.zeros((N_UNITS, N_PERIODS - 1)))},
        "features.point_forecast has shape (20, 7), but y has shape (20, 8)",
        id="features-forecast-shape",
    ),
    pytest.param(
        {"features": _features(point_forecast_draws=np.zeros((3, N_UNITS, N_PERIODS - 1)))},
        "features.point_forecast_draws has shape (3, 20, 7); expected (M, 20, 8)",
        id="features-draws-shape",
    ),
    pytest.param(
        {"features": _features(ps=np.full(N_UNITS, 1.2))},
        "features.ps must hold selection (propensity) scores in [0, 1]",
        id="features-ps-range",
    ),
    pytest.param(
        {"ps": np.full(N_UNITS - 1, 0.5)}, "ps must hold one selection score per unit (N=20)", id="ps-length"
    ),
    pytest.param(
        {"ps": np.full(N_UNITS, -0.1)}, "ps must hold selection (propensity) scores in [0, 1]", id="ps-range"
    ),
    pytest.param(
        {"unit_features": np.zeros((N_UNITS - 1, 2))},
        "unit_features must have one row per unit (N=20)",
        id="unit-features-rows",
    ),
    pytest.param(
        {"baseline_forecast": np.zeros((N_UNITS, N_PERIODS - 1))},
        "baseline_forecast must have shape (N, T) = (20, 8) or (N, T, 1)",
        id="forecast-shape",
    ),
    pytest.param(
        {"baseline_forecast": _poke(np.zeros((N_UNITS, N_PERIODS)), (1, 1), np.nan)},
        "baseline_forecast contains 1 non-finite value(s)",
        id="forecast-nan",
    ),
    pytest.param(
        {"baseline_forecast_draws": np.zeros((N_UNITS, N_PERIODS))},
        "baseline_forecast_draws must have shape (M, N, T) = (M, 20, 8)",
        id="forecast-draws-shape",
    ),
    pytest.param(
        {"pre_anchor_loadings": np.zeros((N_UNITS - 1, 3))},
        "pre_anchor_loadings must have one row per unit (N=20)",
        id="loadings-rows",
    ),
    pytest.param({"coords": np.zeros((N_UNITS, 3))}, "coords must have shape (N, 2) = (20, 2)", id="coords-shape"),
    pytest.param(
        {"coords": _coords()[:, ::-1]}, "coords[:, 0] must be latitudes in [-90, 90] degrees", id="coords-swapped"
    ),
    pytest.param(
        {"coords": _poke(_coords(), (0, 1), 400.0)},
        "coords[:, 1] must be longitudes in [-180, 360] degrees",
        id="coords-longitude",
    ),
    pytest.param(
        {"coords": _poke(_coords(), (0, 0), np.nan)}, "coords contains 1 non-finite value(s)", id="coords-nan"
    ),
    pytest.param(
        {"cluster_ids": np.arange(N_UNITS - 1)},
        "cluster_ids must be a 1-D array with one cluster label per unit (N=20)",
        id="clusters-length",
    ),
    pytest.param(
        {"cluster_ids": _poke(CLUSTERS, 3, np.nan)}, "cluster_ids contains NaN or infinite labels", id="clusters-nan"
    ),
]


@pytest.mark.parametrize("propagated", [False, True], ids=["plug-in", "propagated"])
@pytest.mark.parametrize(("overrides", "match"), _INVALID_INPUTS)
def test_invalid_inputs_raise(overrides, match, propagated):
    with pytest.raises(ValueError, match=re.escape(match)):
        _validate(RankKLongBetConfig(uncertainty_propagation=propagated), **overrides)


# ---- Warnings: only for ignored or contradictory explicit inputs ----------------------------

_WARNING_CASES = [
    pytest.param(
        {"warm_burnin": 60},
        {},
        "cfg.warm_burnin=60 is deprecated and will be removed; in this fit it has no effect",
        id="warm-burnin-no-effect",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "warm_burnin": 150},
        {},
        "raises the mode-1 burn-in from num_burnin=100 to 150 sweeps",
        id="warm-burnin-raises-burn-in",
    ),
    pytest.param(
        {"uncertainty_propagation": False},
        {"coords": _coords()},
        "With the plug-in variant (uncertainty_propagation=False, engine 'stochtree') coords is ignored",
        id="plug-in-coords",
    ),
    pytest.param(
        {"uncertainty_propagation": False},
        {
            "cluster_ids": CLUSTERS,
            "pre_anchor_loadings": np.zeros((N_UNITS, 3)),
            "baseline_forecast_draws": np.zeros((8, N_UNITS, N_PERIODS)),
        },
        "cluster_ids, pre_anchor_loadings, baseline_forecast_draws are ignored",
        id="plug-in-propagated-inputs",
    ),
    pytest.param(
        {
            "uncertainty_propagation": False,
            "n_forecast_draws": 16,
            "num_chains": 4,
            "draws_per_component": 50,
            "use_spatial_gp": False,
        },
        {},
        "cfg.n_forecast_draws=16, cfg.num_chains=4, cfg.draws_per_component=50, cfg.use_spatial_gp=False are ignored",
        id="plug-in-propagated-options",
    ),
    pytest.param(
        {"uncertainty_propagation": False, "max_forecast_draws": 16},
        {},
        "cfg.max_forecast_draws=16 is ignored",
        id="plug-in-expansion",
    ),
    pytest.param(
        {"uncertainty_propagation": False, "baseline_cells": "all"},
        {},
        "cfg.baseline_cells='all' is ignored",
        id="plug-in-baseline-cells",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "num_gfr": 5, "ridge_eps": 1e-3, "num_threads": 4},
        {},
        "With the propagated variant (uncertainty_propagation=True, engine 'jax') cfg.num_gfr=5, "
        "cfg.ridge_eps=0.001, cfg.num_threads=4 are ignored",
        id="propagated-plug-in-options",
    ),
    pytest.param(
        {"use_forecast": False, "theta_f_clip": (0.0, 1.0)},
        {},
        "No forecast is used (use_forecast=False), so cfg.theta_f_clip=(0.0, 1.0) is ignored",
        id="no-forecast-clip",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "use_forecast": False},
        {"baseline_forecast": np.zeros((N_UNITS, N_PERIODS))},
        "No forecast is used (use_forecast=False), so baseline_forecast is ignored",
        id="no-forecast-baseline",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "theta_prior_sd": 1.0},
        {"features": None},
        "No forecast is used (neither features nor baseline_forecast supplies a forecast), so "
        "cfg.theta_prior_sd=1.0 is ignored",
        id="no-forecast-prior",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "n_forecast_draws": 1},
        {},
        "uncertainty_propagation=True contradicts n_forecast_draws=1",
        id="single-forecast-draw",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "n_forecast_draws": 0},
        {},
        "uncertainty_propagation=True contradicts n_forecast_draws=0: with fewer than 2 forecast draws",
        id="zero-forecast-draws",
    ),
    pytest.param(
        {"uncertainty_propagation": True, **PLUG_IN},
        {},
        "cfg selects the propagated variant (uncertainty_propagation=True) but holds the default budget of the "
        "plug-in variant (num_chains=1, num_burnin=200, num_trees_mu=75, num_trees_tau_mode1=35, "
        "num_trees_tau_higher=20)",
        id="plug-in-budget-with-propagation",
    ),
    pytest.param(
        {"uncertainty_propagation": True},
        {"baseline_forecast_draws": np.zeros((4, N_UNITS, N_PERIODS))},
        "baseline_forecast_draws holds M=4 draws, fewer than max_forecast_draws=8",
        id="too-few-forecast-draws",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "num_mcmc": 1000},
        {},
        "cfg.num_mcmc=1000 is ignored because forecast draws are propagated: the posterior has "
        "num_chains x n_forecast_draws x draws_per_component = 2 x 8 x 25 = 400 draws",
        id="num-mcmc-ignored",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "between_var_threshold": 0.3},
        {},
        "cfg.between_var_threshold=0.3 is ignored: it only triggers the adaptive forecast-draw expansion",
        id="threshold-without-expansion",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "use_forecast": False, "draws_per_component": 50},
        {},
        "No forecast draws are propagated (use_forecast=False), so cfg.draws_per_component=50 is ignored",
        id="no-forecast-draws",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "use_factor_loadings": False, "n_factor_loadings": 5},
        {"pre_anchor_loadings": np.zeros((N_UNITS, 3))},
        "use_factor_loadings=False, so pre_anchor_loadings, cfg.n_factor_loadings=5 are ignored",
        id="loadings-off",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "n_factor_loadings": 5},
        {},
        "cfg.n_factor_loadings=5 is ignored because the factor loadings come from "
        "features.pre_anchor_loadings (3 column(s))",
        id="loadings-from-features",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "n_factor_loadings": 4},
        {"pre_anchor_loadings": np.zeros((N_UNITS, 5))},
        "cfg.n_factor_loadings=4 is ignored because the factor loadings come from pre_anchor_loadings "
        "(5 column(s))",
        id="loadings-explicit-count-differs",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "n_factor_loadings": 5},
        {"features": None, "t_0_panel": 1},
        "cfg.n_factor_loadings=5 is ignored: a pre-period of 1 period is too short",
        id="pre-period-too-short",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "use_spatial_gp": False},
        {"coords": _coords()},
        "use_spatial_gp=False, so coords is ignored and no spatial random field is fitted",
        id="spatial-off-coords",
    ),
    pytest.param(
        {"uncertainty_propagation": True},
        {"coords": _coords(), "cluster_ids": CLUSTERS},
        "Both coords and cluster_ids were supplied: coords take precedence",
        id="coords-and-clusters",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "lengthscale_grid_km": (50.0, 100.0)},
        {"cluster_ids": CLUSTERS},
        "cfg.lengthscale_grid_km=(50.0, 100.0) is ignored: it is the range grid of the Matérn-3/2 field on "
        "coords, and the cluster random effect on cluster_ids has no range",
        id="grid-with-clusters",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "lengthscale_grid_km": (50.0, 100.0)},
        {},
        "without coords or cluster_ids no spatial random field is fitted",
        id="grid-without-spatial-inputs",
    ),
    pytest.param(
        {"rank_k": 3},
        {},
        "rank_k=3 exceeds the post-adoption length L_g of 1 treated cohort(s) (first-adoption period(s) [6] "
        "with L_g = [2])",
        id="short-cohorts",
    ),
]


@pytest.mark.parametrize(("cfg_kwargs", "overrides", "fragment"), _WARNING_CASES)
def test_ignored_or_conflicting_inputs_warn_once(cfg_kwargs, overrides, fragment):
    caught = _caught_warnings(RankKLongBetConfig(**cfg_kwargs), **overrides)
    messages = [str(w.message) for w in caught]
    assert len(caught) == 1, messages
    assert issubclass(caught[0].category, RankKConfigWarning)
    assert fragment in messages[0]


def test_single_forecast_draw_with_expansion_limit_warns_twice():
    cfg = RankKLongBetConfig(uncertainty_propagation=True, n_forecast_draws=1, max_forecast_draws=6)
    messages = [str(w.message) for w in _caught_warnings(cfg)]
    assert len(messages) == 2, messages
    assert "uncertainty_propagation=True contradicts n_forecast_draws=1" in messages[0]
    expected = "No forecast draws are propagated (n_forecast_draws=1), so cfg.max_forecast_draws=6 is ignored"
    assert expected in messages[1]


def test_budget_of_the_other_variant_warns():
    # dataclasses.replace(cfg, engine=None, uncertainty_propagation=...) keeps the stored budget.
    to_propagated = dataclasses.replace(
        RankKLongBetConfig(uncertainty_propagation=False), engine=None, uncertainty_propagation=True
    )
    assert to_propagated == RankKLongBetConfig(uncertainty_propagation=True, **PLUG_IN)
    messages = [str(w.message) for w in _caught_warnings(to_propagated)]
    assert len(messages) == 1, messages
    assert "but holds the default budget of the plug-in variant" in messages[0]
    assert "Build a new RankKLongBetConfig(uncertainty_propagation=True, ...) instead" in messages[0]
    to_plug_in = dataclasses.replace(RankKLongBetConfig(), engine=None, uncertainty_propagation=False)
    assert to_plug_in == RankKLongBetConfig(uncertainty_propagation=False, **PROPAGATED)
    messages = [str(w.message) for w in _caught_warnings(to_plug_in)]
    assert len(messages) == 2, messages
    assert (
        "cfg selects the plug-in variant (uncertainty_propagation=False) but holds the default budget of the "
        "propagated variant (num_chains=2, num_burnin=100, num_trees_mu=50, num_trees_tau_mode1=30, "
        "num_trees_tau_higher=15)"
    ) in messages[0]
    assert "cfg.num_chains=2 is ignored" in messages[1]
    # A budget that differs from the other variant's default in any field is a deliberate choice.
    assert not _caught_warnings(RankKLongBetConfig(uncertainty_propagation=True, **{**PLUG_IN, "num_burnin": 300}))


def test_warnings_point_at_the_caller():
    with pytest.warns(RankKConfigWarning, match=re.escape("coords is ignored")) as record:
        _validate(RankKLongBetConfig(uncertainty_propagation=False), coords=_coords())
    assert len(record) == 1
    assert os.path.basename(record[0].filename) == os.path.basename(__file__)


_QUIET_CASES = [
    pytest.param({}, {}, id="propagated-default"),
    pytest.param({"uncertainty_propagation": False}, {}, id="plug-in-default"),
    pytest.param(
        {},
        {
            "features": None,
            "unit_features": np.zeros((N_UNITS, 2)),
            "ps": np.full(N_UNITS, 0.5),
            "baseline_forecast": np.zeros((N_UNITS, N_PERIODS)),
        },
        id="propagated-explicit-arrays",
    ),
    pytest.param({"engine": "stochtree"}, {}, id="engine-stochtree"),
    pytest.param(
        {"uncertainty_propagation": False, "num_chains": 1, "num_burnin": 200},
        {},
        id="plug-in-explicit-defaults",
    ),
    pytest.param(
        {
            "uncertainty_propagation": False,
            "num_burnin": 400,
            "num_mcmc": 1000,
            "num_gfr": 20,
            "num_trees_mu": 100,
            "num_trees_tau_mode1": 50,
            "num_trees_tau_higher": 25,
            "theta_f_clip": (0.0, 2.0),
            "ridge_eps": 1e-3,
            "num_threads": 4,
            "random_seed": 1,
        },
        {},
        id="plug-in-budget",
    ),
    pytest.param({"use_forecast": False}, {}, id="propagated-no-forecast-default"),
    pytest.param({"uncertainty_propagation": False, "use_forecast": False}, {}, id="plug-in-no-forecast"),
    pytest.param({"rank_k": 1}, {}, id="rank-1"),
    pytest.param(
        {"uncertainty_propagation": False}, {"cohorts": _P.cohorts, "t_0_panel": 3}, id="plug-in-explicit-timing"
    ),
    pytest.param({"engine": "jax"}, {"coords": _coords()}, id="engine-jax-coords"),
    pytest.param({"uncertainty_propagation": True}, {"coords": _coords()}, id="propagated-coords"),
    pytest.param({"uncertainty_propagation": True}, {"cluster_ids": CLUSTERS}, id="propagated-clusters"),
    pytest.param(
        {"uncertainty_propagation": True, "num_chains": 4, "num_burnin": 300, "draws_per_component": 100},
        {"coords": _coords()},
        id="application-budget",
    ),
    pytest.param({"uncertainty_propagation": True, "use_spatial_gp": False}, {}, id="ablation-no-spatial-field"),
    pytest.param({"uncertainty_propagation": True, "use_factor_loadings": False}, {}, id="ablation-no-loadings"),
    pytest.param(
        {"uncertainty_propagation": True, "baseline_cells": "all"}, {"coords": _coords()}, id="ablation-all-cells"
    ),
    pytest.param(
        {"uncertainty_propagation": True, "max_forecast_draws": 16, "between_var_threshold": 0.3},
        {},
        id="adaptive-expansion",
    ),
    pytest.param({"uncertainty_propagation": True, "n_forecast_draws": 16}, {}, id="more-forecast-draws"),
    pytest.param(
        {"uncertainty_propagation": True},
        {"baseline_forecast_draws": np.zeros((8, N_UNITS, N_PERIODS))},
        id="explicit-forecast-draws",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "use_forecast": False, "num_mcmc": 1000},
        {"coords": _coords()},
        id="propagated-no-forecast",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "n_factor_loadings": 5}, {"features": None}, id="loadings-from-pre-period"
    ),
    pytest.param(
        {"uncertainty_propagation": True, "n_factor_loadings": 5},
        {"pre_anchor_loadings": np.zeros((N_UNITS, 5))},
        id="explicit-loadings-count-matches",
    ),
    pytest.param(
        {"uncertainty_propagation": True, "n_factor_loadings": 5},
        {"features": _features(pre_anchor_loadings=np.zeros((N_UNITS, 5)))},
        id="feature-loadings-count-matches",
    ),
    pytest.param(
        {"uncertainty_propagation": True}, {"cohorts": _P.cohorts, "t_0_panel": 3}, id="propagated-explicit-timing"
    ),
]


@pytest.mark.parametrize(("cfg_kwargs", "overrides"), _QUIET_CASES)
def test_supported_configurations_do_not_warn(cfg_kwargs, overrides):
    caught = _caught_warnings(RankKLongBetConfig(**cfg_kwargs), **overrides)
    assert not caught, [str(w.message) for w in caught]


# ---- Dispatch of the default fit and the placebo audit (samplers stubbed out) ---------------


def test_fit_without_cfg_uses_the_propagated_sampler_and_plug_in_cfg_uses_stochtree(monkeypatch):
    seen_jax = []
    seen_stochtree = {}

    def fake_build_cache(**kwargs):
        return types.SimpleNamespace(enabled=False)

    def fake_sample_mode(**kwargs):
        seen_jax.append(kwargs)
        n_chains = int(kwargs["num_chains"])
        n_draws_per_chain = int(kwargs["num_draws_per_chain"])
        df_draws = kwargs["df_k_draws"]
        m_comp = int(df_draws.shape[0]) if df_draws is not None else 1
        s_total = n_chains * m_comp * n_draws_per_chain
        n_units = int(kwargs["dy_k"].shape[0])
        rng = np.random.default_rng(0)
        return {
            "tau_draws": rng.normal(size=(n_units, s_total)),
            "comp_tau_draws": rng.normal(size=(m_comp, n_chains * n_draws_per_chain, n_units)),
            "tau_chains": rng.normal(size=(n_chains, m_comp * n_draws_per_chain, n_units)),
            "theta_f_draws": np.ones(s_total),
            "sigma_draws": np.full(s_total, 0.1),
            "s_zeta_draws": np.zeros(s_total),
            "ell_draws": np.zeros(s_total),
        }

    def fake_stochtree(**kwargs):
        seen_stochtree.update(kwargs)
        return "plug-in result"

    monkeypatch.setattr(rank_k_module, "_get_rank_k_jax", lambda: (fake_build_cache, fake_sample_mode))
    monkeypatch.setattr(rank_k_module, "_fit_rank_k_stochtree", fake_stochtree)
    panel = _panel()
    feats = _features(
        point_forecast_draws=np.zeros((8, N_UNITS, N_PERIODS)),
        forecast_draws={"draw_method": "posterior"},
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res_prop = fit_rank_k_longbet(panel.y, panel.z, panel.x, features=feats)
        res_plug = fit_rank_k_longbet(
            panel.y,
            panel.z,
            panel.x,
            features=feats,
            cfg=RankKLongBetConfig(uncertainty_propagation=False),
        )
    assert not caught, [str(w.message) for w in caught]
    assert len(seen_jax) == 2
    assert res_prop.att_draws.shape == (400,)
    assert res_plug == "plug-in result"
    assert seen_stochtree["cfg"] == RankKLongBetConfig(uncertainty_propagation=False)
    assert seen_stochtree["cfg"].engine == "stochtree"
    assert seen_stochtree["has_forecast"] is True


def _stub_placebo(monkeypatch: pytest.MonkeyPatch) -> dict[str, dict[str, Any]]:
    """Replace the feature builder and the fit of the placebo audit; return their recorded calls."""
    calls: dict[str, dict[str, Any]] = {}

    def fake_build(y, z, **kwargs):
        calls["build"] = {"y": y, "z": z, **kwargs}
        return _features()  # its content is irrelevant: the fit is stubbed as well

    def fake_fit(y, z, x, **kwargs):
        calls["fit"] = {"y": y, "z": z, **kwargs}
        return types.SimpleNamespace(att_hat=0.1, att_se=0.2, att_ci_low=-0.3, att_ci_high=0.5)

    monkeypatch.setattr(rank_k_module, "build_causal_features", fake_build)
    monkeypatch.setattr(rank_k_module, "_fit_rank_k_longbet", fake_fit)
    return calls


def test_placebo_audit_warns_that_features_are_ignored(monkeypatch):
    calls = _stub_placebo(monkeypatch)
    panel = _panel()
    y_hist = np.random.default_rng(2).normal(size=(N_UNITS, 6))
    with pytest.warns(RankKConfigWarning, match=re.escape("rank_k_placebo_audit ignores `features`")) as record:
        audit = rank_k_placebo_audit(panel.y, panel.z, panel.x, y_pre_history=y_hist, features=_features())
    assert len(record) == 1
    assert audit.passes_audit and not audit.flag_parallel_trends_violation
    # t0 = 4 and 6 history periods: L_pre = 10, h = 4 and T_pl = 8, so the pseudo adoption is at 4.
    assert (audit.placebo_anchor_t, audit.placebo_horizon) == (4, 4)
    fit = calls["fit"]
    assert fit["t_0_panel"] == 4
    assert fit["warn_short_cohorts"] is False
    assert fit["cfg"] == RankKLongBetConfig(rank_k=1)  # the placebo ATT depends on mode 1 only
    np.testing.assert_array_equal(fit["cohorts"], np.r_[np.full(6, 4), np.full(4, 6), np.full(10, 8)])
    np.testing.assert_array_equal(fit["y"], np.concatenate([y_hist, panel.y[:, :4]], axis=1)[:, 2:])
    assert calls["build"]["seed"] == 2026
    assert calls["build"]["n_factor_loadings"] == 3
    assert calls["build"]["y_pre_history"].shape == (N_UNITS, 2)


def test_placebo_audit_rebuilds_features_with_the_configured_loadings(monkeypatch):
    calls = _stub_placebo(monkeypatch)
    panel = _panel()
    y_hist = np.random.default_rng(2).normal(size=(N_UNITS, 6))
    cfg = RankKLongBetConfig(uncertainty_propagation=True, n_factor_loadings=5, random_seed=7)
    rank_k_placebo_audit(panel.y, panel.z, panel.x, y_pre_history=y_hist, cluster_ids=CLUSTERS, cfg=cfg)
    assert (calls["build"]["n_factor_loadings"], calls["build"]["seed"]) == (5, 7)
    assert calls["fit"]["cfg"] == dataclasses.replace(cfg, rank_k=1)
    np.testing.assert_array_equal(calls["fit"]["cluster_ids"], CLUSTERS)


def test_placebo_audit_without_features_is_quiet(monkeypatch):
    _stub_placebo(monkeypatch)
    panel = _panel()
    y_hist = np.random.default_rng(2).normal(size=(N_UNITS, 6))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rank_k_placebo_audit(panel.y, panel.z, panel.x, y_pre_history=y_hist)
    assert not caught, [str(w.message) for w in caught]


def test_placebo_audit_validates_before_rebuilding_features(monkeypatch):
    calls = _stub_placebo(monkeypatch)
    panel = _panel()
    y_hist = np.random.default_rng(2).normal(size=(N_UNITS, 6))
    with pytest.raises(ValueError, match="z must be absorbing"):
        rank_k_placebo_audit(panel.y, _poke(panel.z, (0, 6), 0.0), panel.x, y_pre_history=y_hist)
    assert not calls


@pytest.mark.parametrize("offset", [2.5, "4", True, None, float("nan"), float("inf")])
def test_placebo_audit_rejects_an_offset_that_is_not_a_whole_number(monkeypatch, offset):
    calls = _stub_placebo(monkeypatch)
    panel = _panel()
    y_hist = np.random.default_rng(2).normal(size=(N_UNITS, 6))
    with pytest.raises(ValueError, match="pseudo_anchor_offset must be a whole number"):
        rank_k_placebo_audit(panel.y, panel.z, panel.x, y_pre_history=y_hist, pseudo_anchor_offset=offset)
    assert not calls


@pytest.mark.parametrize(("offset", "horizon"), [(4.0, 4), (np.int64(3), 3), (1, 2)])
def test_placebo_audit_accepts_whole_number_offsets(monkeypatch, offset, horizon):
    _stub_placebo(monkeypatch)
    panel = _panel()
    y_hist = np.random.default_rng(2).normal(size=(N_UNITS, 6))
    audit = rank_k_placebo_audit(panel.y, panel.z, panel.x, y_pre_history=y_hist, pseudo_anchor_offset=offset)
    # t0 = 4, so the pseudo adoption stays at period 4 and the placebo post-period has h periods.
    assert (audit.placebo_anchor_t, audit.placebo_horizon) == (4, horizon)


# ---- Robust kriging operator of the JAX sampler -------------------------------------------


def _fibonacci_sphere(n: int = 40) -> np.ndarray:
    """``n`` (latitude, longitude) points spread evenly over the globe."""
    i = np.arange(n) + 0.5
    lat = np.degrees(np.arcsin(1.0 - 2.0 * i / n))
    lon = (np.degrees(np.pi * (1.0 + np.sqrt(5.0)) * i) % 360.0) - 180.0
    return np.column_stack([lat, lon])


_GLOBE = _fibonacci_sphere()
_CTRL, _TRT = np.arange(0, 40, 2), np.arange(1, 40, 2)


def _kriging_blocks(ell_km: float) -> tuple[np.ndarray, np.ndarray]:
    """Jittered control block ``C00`` and treated-control block ``C10`` as in the kernel cache."""
    corr = matern32_corr(great_circle_distance_km(_GLOBE, _GLOBE), ell_km)
    return corr[np.ix_(_CTRL, _CTRL)] + 1e-6 * np.eye(_CTRL.size), corr[np.ix_(_TRT, _CTRL)]


def _cholesky_kriging(c00: np.ndarray, c10: np.ndarray) -> np.ndarray:
    """The kriging operator as computed before 0.5.0 (Cholesky solve only)."""
    with enable_x64():
        chol = np.linalg.cholesky(c00)
        return np.asarray(jax.scipy.linalg.cho_solve((chol, True), c10.T), dtype=np.float64).T


def test_kriging_operator_is_unchanged_for_a_positive_definite_control_block():
    c00, c10 = _kriging_blocks(1600.0)
    np.testing.assert_array_equal(_kriging_operator(c00, c10), _cholesky_kriging(c00, c10))


def test_kriging_operator_handles_an_indefinite_control_block():
    c00, c10 = _kriging_blocks(20000.0)
    # Matérn-3/2 on great-circle distances is not positive definite on the sphere at long ranges.
    with pytest.raises(np.linalg.LinAlgError):
        np.linalg.cholesky(c00)
    kriging = _kriging_operator(c00, c10)
    assert kriging.shape == (_TRT.size, _CTRL.size)
    assert np.all(np.isfinite(kriging))
    vals, vecs = np.linalg.eigh(0.5 * (c00 + c00.T))
    floored = np.maximum(vals, 1e-6 * vals.max())
    np.testing.assert_allclose(kriging, ((c10 @ vecs) / floored[None, :]) @ vecs.T, rtol=1e-10, atol=1e-12)


def test_spatial_kernel_cache_builds_for_global_coordinates():
    cache = build_spatial_kernel_cache(
        ctrl_idx=_CTRL, trt_idx=_TRT, w_unit=np.ones(_GLOBE.shape[0]), coords=_GLOBE, ell_grid_km=(1600.0, 20000.0)
    )
    kriging = np.asarray(cache.krig_mat_10)
    assert kriging.shape == (2, _TRT.size, _CTRL.size)
    for name in ("krig_mat_10", "schur_chol_11", "eig_vals_00", "eig_vecs_00"):
        assert np.all(np.isfinite(np.asarray(getattr(cache, name)))), name
    np.testing.assert_array_equal(kriging[0], _cholesky_kriging(*_kriging_blocks(1600.0)))
    clusters = build_spatial_kernel_cache(
        ctrl_idx=_CTRL, trt_idx=_TRT, w_unit=np.ones(_GLOBE.shape[0]), cluster_ids=np.repeat(np.arange(8), 5)
    )
    assert clusters.enabled
    assert np.all(np.isfinite(np.asarray(clusters.krig_mat_10)))
