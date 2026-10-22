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

"""Serialization for multi-outcome LongBet models to self-contained NPZ archives."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

import jax.numpy as jnp
import numpy as np

from longbet._config import LongBetConfig
from longbet._design import Design
from longbet._io import (
    _FOREST_KEYS, _PARAM_KEYS, _rebuild_trace,
    PRECISION_CACHE_VERSION, _require_current_precision_cache,
    ORDINAL_SCHEMA_VERSION, ordinal_metadata, validate_ordinal_arrays,
)
from longbet._loop import LongBetTrace
from longbet._multi_loop import MultiLongBetTrace
from longbet._multi_input import resolve_category_counts
from longbet._sur import SAMPLER_SEMANTICS

if TYPE_CHECKING:
    from longbet._multi_model import LongBetMulti, OutcomeList


def _validate_archive(data: Any, meta: dict[str, Any], config: LongBetConfig) -> None:
    """Reject inconsistent axes/order instead of silently pairing unrelated draws."""
    def require(ok: bool, message: str) -> None:
        if not ok:
            raise ValueError(f"Invalid multi archive: {message}")

    M = meta.get("M")
    require(type(M) is int and M >= 2, "M must be an integer >= 2.")
    names = meta.get("outcome_names", [])
    types = meta.get("user_outcomes", [])
    require(len(names) == M and all(isinstance(n, str) and n.strip() for n in names), "invalid outcome names.")
    require(len(set(names)) == M, "duplicate outcome names.")
    require(len(types) == M and all(t in ("continuous", "binary", "ordinal") for t in types), "invalid outcome types.")
    ordinal = "ordinal" in types
    order = meta.get("order", [])
    inverse = meta.get("inverse_order", [])
    require(len(order) == M and all(type(i) is int for i in order) and sorted(order) == list(range(M)), "invalid order permutation.")
    require(inverse == np.argsort(order).tolist(), "inverse_order does not invert order.")
    expected_order = [i for discrete in (True, False) for i in range(M) if (types[i] != "continuous") == discrete]
    require(order == expected_order, "order is not the stable discrete-first (binary-first for legacy fits) permutation.")
    require(meta.get("internal_names") == [names[i] for i in order], "internal names disagree with order.")
    require(meta.get("internal_outcomes") == [types[i] for i in order], "internal types disagree with order.")
    for key in ("meany", "sdy", "offset_"):
        values = np.asarray(meta.get(key, []), dtype=float)
        require(values.shape == (M,) and np.all(np.isfinite(values)), f"invalid {key} scales.")
        if key == "sdy":
            require(np.all(values > 0), "sdy must be positive.")
    for key in ("N", "T", "S_max"):
        value = meta.get(key)
        require(type(value) is int and value >= (0 if key == "S_max" else 1), f"invalid {key}.")
    fitted_t = np.asarray(meta.get("fitted_t"), dtype=float)
    require(fitted_t.shape == (meta["T"],) and np.all(np.isfinite(fitted_t)) and np.all(np.diff(fitted_t) > 0), "invalid fitted_t.")
    require(meta.get("sur_active") == config.sur_active, "SUR settings disagree with config.")
    require(meta.get("sur_prior_var") == config.sur_prior_var, "SUR prior disagrees with config.")
    require(meta.get("sampler_semantics") == SAMPLER_SEMANTICS,
            "unsupported sampler semantics.")
    require(meta.get("format_version", 1) == (3 if ordinal else 1), "format does not match ordinal outcomes.")
    if ordinal:
        require(type(meta.get("ordinal_schema_version")) is int and
                meta["ordinal_schema_version"] == ORDINAL_SCHEMA_VERSION, "unsupported ordinal schema.")
        require("num_categories" in meta and "internal_num_categories" in meta,
                "missing ordinal category counts.")
        counts = resolve_category_counts(meta["num_categories"], names, types)
        require(meta["num_categories"] == list(counts), "invalid category metadata.")
        require(meta["internal_num_categories"] == [counts[i] for i in order],
                "internal category counts disagree with order.")
        require(len(meta.get("ordinal_metadata", [])) == M, "missing per-outcome ordinal metadata.")
    require(meta.get("rng_scheme_version") == 1, "unsupported RNG scheme.")
    chains = config.num_chains > 1
    require(meta.get("has_chains") is chains, "chain metadata disagrees with config.")
    draw_shape = ((config.num_chains,) if chains else ()) + (config.num_sweeps,)

    def array(key: str, shape: tuple[int, ...]) -> np.ndarray:
        require(key in data, f"missing {key}.")
        a = np.asarray(data[key])
        require(a.shape == shape and np.all(np.isfinite(a)), f"invalid shape or non-finite values in {key}; expected {shape}.")
        return a

    gamma = array("gamma_loadings", (*draw_shape, M, M))
    require(np.all(np.triu(gamma) == 0), "Gamma must be strictly lower triangular.")
    for m, typ in enumerate(meta["internal_outcomes"]):
        if typ != "continuous" or not config.sur_active:
            require(np.all(gamma[..., m, :] == 0), "inactive/discrete Gamma row must be zero.")
        if ordinal:
            child_meta = meta["ordinal_metadata"][m]
            if typ == "ordinal":
                require(isinstance(child_meta, dict), "missing ordinal child metadata.")
                validate_ordinal_arrays(data, f"outcome_{m}_", counts[order[m]],
                    config.cutpoint_prior_scale, child_meta, draw_shape)
                require(meta["meany"][m] == 0 and meta["sdy"][m] == 1, "ordinal scale must be unstandardized.")
            else:
                require(child_meta is None, "nonordinal child has ordinal metadata.")
        for key in _PARAM_KEYS:
            extra = (meta["S_max"] + 1,) if key == "beta" else (meta["N"],) if key == "gamma" else ()
            a = array(f"outcome_{m}_{key}", (*draw_shape, *extra))
            if key == "sigma2":
                require(np.all(a > 0), "sigma2 must be positive.")
                if typ != "continuous":
                    require(np.all(a == 1), "discrete sigma2 must equal one.")
        for forest in ("mu", "nu"):
            for key in _FOREST_KEYS:
                name = f"outcome_{m}_{forest}_{key}"
                require(name in data, f"missing {name}.")
                a = np.asarray(data[name])
                # Offsets and leaf units are shared constants, not traces.
                expected_prefix = () if key in ("offset", "leaf_unit") else draw_shape
                require(a.shape[:len(expected_prefix)] == expected_prefix and np.all(np.isfinite(a)), f"invalid draw axes or non-finite values in {name}.")
                if ordinal and key in ("leaf_tree", "var_tree", "split_tree"):
                    J = config.num_trees_pr if forest == "mu" else config.num_trees_trt
                    depth = config.max_depth_pr if forest == "mu" else config.max_depth_trt
                    slots = 2**depth if key == "leaf_tree" else 2**(depth-1)
                    require(a.shape == (*draw_shape, J, slots), f"invalid combined forest shape in {name}.")


def save_multi_npz(model: LongBetMulti, path: str | Path) -> None:
    """Save a fitted LongBetMulti model to an NPZ archive.

    Parameters
    ----------
    model
        Fitted LongBetMulti instance.
    path
        Output file path.
    """
    if model.trace is None or model.design_ is None:
        raise RuntimeError("Cannot save an unfitted LongBetMulti model.")

    # Local, like the loader's: longbet._model imports this module's siblings,
    # so a module-level import would be circular.
    from longbet._model import LongBet

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    M = len(model.outcome_names)
    has_chains = bool(np.asarray(model.trace.gamma_loadings).ndim > 3)

    meta: dict[str, Any] = {
        "model_kind": "multi",
        "format_version": 3 if "ordinal" in model.outcome else 1,
        "precision_cache_version": PRECISION_CACHE_VERSION,
        "M": M,
        "outcome_names": list(model.outcome_names),
        "user_outcomes": list(model.outcome),
        "order": list(model.order),
        "inverse_order": list(model.inverse_order),
        "internal_names": [model.outcome_names[i] for i in model.order],
        "internal_outcomes": [model.outcome[i] for i in model.order],
        "meany": [float(x) for x in model.meany],
        "sdy": [float(x) for x in model.sdy],
        "offset_": [float(x) for x in model.offset_],
        "sur_active": bool(model.sur_active),
        "sur_prior_var": float(model.config.sur_prior_var),
        "sampler_semantics": str(model.sampler_semantics),
        "rng_scheme_version": int(model.rng_scheme_version),
        "provenance": str(model.provenance),
        "N": int(model.N_),
        "T": int(model.T_),
        "S_max": int(model.S_max_),
        "has_chains": has_chains,
        "design_spec": model.design_.to_dict(),
        "fitted_t": (
            np.asarray(model.fitted_t_).tolist()
            if model.fitted_t_ is not None
            else None
        ),
        # The shared calendar-time block. Serialized by the *scalar* model's
        # own method -- it reads only ``trend_time_basis_``,
        # ``trend_period_basis_`` and ``trend_feature_spec_``, which
        # LongBetMulti carries under the same names -- so the two archive
        # formats cannot drift apart. ``None`` when the fit had no block, and
        # older archives without the key load as ``None`` too.
        "trend_basis": (
            LongBet._trend_basis_to_dict(model)
            if getattr(model, "trend_time_basis_", None) is not None
            else None
        ),
    }

    if "ordinal" in model.outcome:
        meta.update(ordinal_schema_version=ORDINAL_SCHEMA_VERSION,
            num_categories=list(model.num_categories),
            internal_num_categories=list(model.internal_num_categories),
            ordinal_metadata=[ordinal_metadata(model.internal_num_categories[m], model.config.cutpoint_prior_scale)
                              if model.outcome[model.order[m]] == "ordinal" else None for m in range(M)])

    arrays: dict[str, Any] = {
        "_config_json": json.dumps(model.config.to_dict()),
        "_meta_json": json.dumps(meta),
        "gamma_loadings": np.asarray(model.trace.gamma_loadings),
    }

    # Store traces for each internal equation
    for m in range(M):
        tr = model.trace.traces[m]
        if model.outcome[model.order[m]] == "ordinal":
            if tr.cutpoints is None:
                raise ValueError("Ordinal multi archives require cutpoint traces")
            arrays[f"outcome_{m}_cutpoints"] = np.asarray(tr.cutpoints)
        for k in _PARAM_KEYS:
            arrays[f"outcome_{m}_{k}"] = np.asarray(getattr(tr, k))
        # As in ``longbet._io.save_npz``: written only when the fit had a trend
        # block, so archives of fits without one stay byte-compatible with
        # readers that predate it.
        trend_coef = getattr(tr, "trend_coef", None)
        if trend_coef is not None and np.asarray(trend_coef).shape[-1] > 0:
            arrays[f"outcome_{m}_trend_coef"] = np.asarray(trend_coef)
        for prefix, sub in (("mu", tr.mu_trace), ("nu", tr.nu_trace)):
            for key in _FOREST_KEYS:
                arrays[f"outcome_{m}_{prefix}_{key}"] = np.asarray(getattr(sub, key))

    if "ordinal" in model.outcome:
        _validate_archive(arrays, meta, model.config)
    np.savez_compressed(path, **arrays)


def load_multi_npz(path: str | Path) -> LongBetMulti:
    """Read a LongBetMulti archive saved by :func:`save_multi_npz`.

    Parameters
    ----------
    path
        Path to the archive.

    Returns
    -------
    Reconstituted LongBetMulti instance.
    """
    from longbet._model import LongBet
    from longbet._multi_model import LongBetMulti, OutcomeList

    data = np.load(Path(path), allow_pickle=False)

    if "_meta_json" not in data:
        raise ValueError(
            f"Archive at {path} is not a valid LongBet archive (missing _meta_json)."
        )

    meta = json.loads(str(data["_meta_json"]))
    model_kind = meta.get("model_kind", "scalar")
    if model_kind != "multi":
        raise ValueError(
            f"Archive at {path} is a single-outcome model (kind={model_kind!r}). "
            f"Use LongBet.load(path) rather than LongBetMulti.load(path)."
        )

    format_version = meta.get("format_version", 1)
    if format_version not in (1, 2, 3):
        raise ValueError(
            f"Unsupported multi format_version {format_version}; expected 1, 2, or 3."
        )

    config = LongBetConfig.from_dict(json.loads(str(data["_config_json"])))
    _require_current_precision_cache(meta)
    _validate_archive(data, meta, config)
    M = int(meta["M"])
    has_chains = bool(meta.get("has_chains", False))
    design = Design.from_dict(meta["design_spec"])
    design.validate_panel_grids(int(meta["T"]), int(meta["S_max"]))

    # Rebuild child traces in internal order
    child_traces = []
    for m in range(M):
        key = f"outcome_{m}_trend_coef"
        tr = LongBetTrace(
            mu_trace=_rebuild_trace(data, f"outcome_{m}_mu", has_chains),
            nu_trace=_rebuild_trace(data, f"outcome_{m}_nu", has_chains),
            **{k: jnp.asarray(data[f"outcome_{m}_{k}"]) for k in _PARAM_KEYS},
            trend_coef=jnp.asarray(data[key]) if key in data.files else None,
            cutpoints=(jnp.asarray(data[f"outcome_{m}_cutpoints"])
                       if meta["internal_outcomes"][m] == "ordinal" else None),
        )
        child_traces.append(tr)

    gamma_loadings = jnp.asarray(data["gamma_loadings"])
    multi_trace = MultiLongBetTrace(
        traces=tuple(child_traces),
        gamma_loadings=gamma_loadings,
    )

    model = LongBetMulti(config)
    model.trace = multi_trace
    model.design_ = design
    model.outcome_names = tuple(meta["outcome_names"])
    model.outcome = tuple(meta["user_outcomes"])
    model.num_categories = tuple(meta.get("num_categories", [None]*M))
    model.internal_num_categories = tuple(meta.get("internal_num_categories", [None]*M))
    model.order = tuple(meta["order"])
    model.inverse_order = tuple(meta["inverse_order"])
    model.meany = tuple(meta["meany"])
    model.sdy = tuple(meta["sdy"])
    model.offset_ = tuple(meta["offset_"])
    model.sur_active = bool(meta["sur_active"])
    model.N_ = int(meta["N"])
    model.T_ = int(meta["T"])
    model.S_max_ = int(meta["S_max"])
    model.fitted_t_ = (
        np.asarray(meta["fitted_t"], dtype=np.float32)
        if meta.get("fitted_t") is not None
        else None
    )
    model.provenance = str(meta.get("provenance", ""))
    model.sampler_semantics = str(meta["sampler_semantics"])
    model.rng_scheme_version = int(meta.get("rng_scheme_version", 1))
    # Symmetric with the saver: the scalar restorer writes exactly the three
    # attributes LongBetMulti shares with LongBet, and tolerates a missing key
    # by clearing them.
    LongBet._trend_basis_from_dict(model, meta.get("trend_basis"))

    # Reconstitute child LongBet fits in user order
    child_fits = []
    for u in range(M):
        internal_idx = model.inverse_order[u]
        child_otype = model.outcome[u]
        child_cfg = dataclasses.replace(config, outcome=child_otype,
                                         num_categories=model.num_categories[u])

        child_model = LongBet(child_cfg)
        child_model.multi_origin = model._child_origin(u)
        child_model.design_ = design
        child_model.trace = multi_trace.traces[internal_idx]
        child_model.meany = model.meany[internal_idx]
        child_model.sdy = model.sdy[internal_idx]
        child_model.offset_ = model.offset_[internal_idx]
        child_model.N_ = model.N_
        child_model.T_ = model.T_
        child_model.S_max_ = model.S_max_
        child_model.t_fit_ = model.fitted_t_
        child_model.fitted_t_ = model.fitted_t_
        # Same sharing as LongBetMulti.fit: a reloaded child must predict with
        # the basis its coefficients were sampled against.
        child_model.trend_time_basis_ = model.trend_time_basis_
        child_model.trend_period_basis_ = model.trend_period_basis_
        child_model.trend_feature_spec_ = model.trend_feature_spec_

        child_fits.append(child_model)

    model.t_fit_ = model.fitted_t_
    model.fits = OutcomeList(child_fits, model.outcome_names)
    return model
