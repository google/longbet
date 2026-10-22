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

"""Model serialization to a self-contained ``.npz`` archive.

The archive holds plain arrays plus two JSON blobs (the config and the
metadata, which carries the design spec).  It contains no pickled objects, so it
loads with ``allow_pickle=False`` and can be embedded in an R object as a raw
vector and round-tripped through ``saveRDS`` without a live Python handle.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jax.numpy as jnp
import numpy as np

from bartz.mcmcloop._trace import MainTrace

from longbet._config import LongBetConfig
from longbet._loop import LongBetTrace

_PARAM_KEYS = ("beta", "gamma", "b0", "b1", "alpha", "sigma2", "sigma_gamma2")
#: Traces of the optional model components (linear moderator term, learned
#: trajectory scale, AR(1) errors, distributed-lag kernel and its scale, fixed
#: extra columns); written only when the fit has them.
_OPTIONAL_KEYS = ("delta", "tau_delta2", "beta_scale2", "ar1_phi", "kappa", "lag_scale2",
                  "extra_coef")
_FOREST_KEYS = ("leaf_tree", "split_tree", "var_tree", "offset", "leaf_unit", "error_cov_inv")

# Both scalar and coupled fits used stale bartz leaf precision sums before this
# contract was introduced. Loading trees cannot repair their posterior draws.
PRECISION_CACHE_VERSION = 1
ORDINAL_SCHEMA_VERSION = 1


def ordinal_metadata(num_categories: int, prior_scale: float) -> dict[str, Any]:
    """Link and identification contract shared by scalar and child archives."""
    return {"num_categories": num_categories, "cutpoint_prior_scale": prior_scale,
            "cutpoint_anchor": "first_finite_zero",
            "ordinal_schema_version": ORDINAL_SCHEMA_VERSION}


def validate_ordinal_arrays(data, prefix, num_categories, prior_scale, metadata, draw_shape):
    """Validate one ordinal equation's metadata, thresholds, and saved variances."""
    def require(ok, message):
        if not ok:
            raise ValueError(f"Invalid ordinal archive: {message}")
    for key, expected in ordinal_metadata(num_categories, prior_scale).items():
        value = metadata.get(key)
        require(type(value) is type(expected) and value == expected,
                f"{prefix}{key} disagrees with the ordinal model schema.")
    key = f"{prefix}cutpoints"
    require(key in data, f"missing {key} (required even for K=2).")
    cuts = np.asarray(data[key])
    require(cuts.dtype.kind == "f" and cuts.shape == (*draw_shape, num_categories - 2),
            f"invalid cutpoint shape/dtype in {key}.")
    require(np.isfinite(cuts).all(), "nonfinite cutpoints.")
    require((np.diff(np.concatenate((np.zeros((*draw_shape, 1)), cuts), axis=-1), axis=-1) > 0).all(),
            "free cutpoints must be strictly ordered and positive.")
    for suffix in ("sigma2", "mu_error_cov_inv", "nu_error_cov_inv"):
        key = f"{prefix}{suffix}"
        require(key in data, f"missing {key}.")
        arr = np.asarray(data[key])
        require(arr.shape == draw_shape and (arr == 1).all(),
                f"{key} must have aligned draws and equal one.")
    return cuts


def _validate_scalar_ordinal(data, config, meta):
    for key in ("N", "T", "S_max"):
        value = meta.get(key)
        if type(value) is not int or value < (0 if key == "S_max" else 1):
            raise ValueError(f"Invalid ordinal archive: invalid {key}")
    chains = config.num_chains > 1
    if meta.get("has_chains") is not chains:
        raise ValueError("Invalid ordinal archive: chain metadata disagrees with config")
    draw_shape = ((config.num_chains,) if chains else ()) + (config.num_sweeps,)
    cuts = validate_ordinal_arrays(data, "", config.num_categories,
                                   config.cutpoint_prior_scale, meta, draw_shape)
    if meta.get("meany") != 0 or meta.get("sdy") != 1:
        raise ValueError("Invalid ordinal archive: latent scale must be unstandardized")
    for key in _PARAM_KEYS:
        if key not in data:
            raise ValueError(f"Invalid ordinal archive: missing {key}")
        arr = np.asarray(data[key])
        extra = (meta["S_max"]+1,) if key == "beta" else (meta["N"],) if key == "gamma" else ()
        if arr.dtype.kind != "f" or arr.shape != (*draw_shape, *extra) or not np.isfinite(arr).all():
            raise ValueError(f"Invalid ordinal archive: invalid shape or values in {key}")
    for forest in ("mu", "nu"):
        trees = config.num_trees_pr if forest == "mu" else config.num_trees_trt
        depth = config.max_depth_pr if forest == "mu" else config.max_depth_trt
        for key in _FOREST_KEYS:
            name = f"{forest}_{key}"
            if name not in data:
                raise ValueError(f"Invalid ordinal archive: missing {name}")
            arr = np.asarray(data[name])
            shape = (() if key in ("offset", "leaf_unit") else draw_shape if key == "error_cov_inv"
                     else (*draw_shape, trees, 2**(depth if key == "leaf_tree" else depth-1)))
            if arr.shape != shape or arr.dtype.kind not in "iuf" or not np.isfinite(arr).all():
                raise ValueError(f"Invalid ordinal archive: invalid shape or values in {name}")
    return cuts


def _require_current_precision_cache(meta: dict[str, Any]) -> None:
    if meta.get("precision_cache_version") != PRECISION_CACHE_VERSION:
        raise ValueError(
            "Archive predates the dynamic leaf-precision cache repair or uses "
            "an unsupported sampler version. Its posterior draws cannot be "
            "upgraded by loading them. Refit from the original data."
        )


def save_npz(
    path: str | Path,
    trace: LongBetTrace,
    config: LongBetConfig,
    meany: float,
    sdy: float,
    metadata: dict[str, Any] | None = None,
) -> None:
    """Write posterior traces, configuration and design spec to ``path``.

    Parameters
    ----------
    path
        Output path; the ``.npz`` suffix is added by NumPy if absent.
    trace
        Posterior trace from ``run_longbet_mcmc``.
    config
        The configuration the model was fitted with.
    meany, sdy
        Outcome centring and scaling.
    metadata
        Extra JSON-serializable metadata; the design spec travels here.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    meta = dict(metadata or {})
    meta["precision_cache_version"] = PRECISION_CACHE_VERSION
    meta["meany"] = float(meany)
    meta["sdy"] = float(sdy)
    meta["has_chains"] = bool(np.asarray(trace.beta).ndim > 2)

    arrays: dict[str, Any] = {k: np.asarray(getattr(trace, k)) for k in _PARAM_KEYS}
    for prefix, sub in (("mu", trace.mu_trace), ("nu", trace.nu_trace)):
        for key in _FOREST_KEYS:
            arrays[f"{prefix}_{key}"] = np.asarray(getattr(sub, key))

    # The smooth trend block. Written only when the fit had one, so archives of
    # models fitted without it stay byte-compatible with older readers.
    trend_coef = getattr(trace, "trend_coef", None)
    if trend_coef is not None and np.asarray(trend_coef).shape[-1] > 0:
        arrays["trend_coef"] = np.asarray(trend_coef)

    for key in _OPTIONAL_KEYS:
        value = getattr(trace, key, None)
        if value is not None:
            arrays[key] = np.asarray(value)

    if config.outcome == "ordinal":
        meta.update(ordinal_metadata(config.num_categories, config.cutpoint_prior_scale))
        if trace.cutpoints is None:
            raise ValueError("Ordinal archives require cutpoints, including the K=2 empty array")
        arrays["cutpoints"] = np.asarray(trace.cutpoints)
        _validate_scalar_ordinal(arrays, config, meta)

    arrays["_config_json"] = json.dumps(config.to_dict())
    arrays["_meta_json"] = json.dumps(meta)

    np.savez_compressed(path, **arrays)


def _rebuild_trace(data: Any, prefix: str, has_chains: bool) -> MainTrace:
    """Rebuild one forest's ``MainTrace`` from the archive."""
    leaf_tree = jnp.asarray(data[f"{prefix}_leaf_tree"])
    shape = (leaf_tree.shape[0], leaf_tree.shape[1]) if has_chains else (leaf_tree.shape[0],)
    zeros = jnp.zeros(shape, dtype=jnp.int32)
    return MainTrace(
        has_chains=has_chains,
        mesh=None,
        leaf_tree=leaf_tree,
        split_tree=jnp.asarray(data[f"{prefix}_split_tree"]),
        var_tree=jnp.asarray(data[f"{prefix}_var_tree"]),
        offset=jnp.asarray(data[f"{prefix}_offset"]),
        leaf_unit=jnp.asarray(data[f"{prefix}_leaf_unit"]),
        error_cov_inv=jnp.asarray(data[f"{prefix}_error_cov_inv"]),
        theta=None,
        grow_prop_count=zeros,
        grow_acc_count=zeros,
        prune_prop_count=zeros,
        prune_acc_count=zeros,
        log_likelihood=None,
        log_trans_prior=None,
        varprob=None,
    )


def load_npz(
    path: str | Path,
) -> tuple[LongBetTrace, LongBetConfig, float, float, dict[str, Any]]:
    """Read an archive written by :func:`save_npz`.

    Returns
    -------
    trace, config, meany, sdy, metadata
    """
    data = np.load(Path(path), allow_pickle=False)

    config = LongBetConfig.from_dict(json.loads(str(data["_config_json"])))
    meta = json.loads(str(data["_meta_json"]))
    if meta.get("model_kind") == "multi":
        raise ValueError(
            f"Archive at {path} contains a multi-outcome LongBet model (model_kind='multi'). "
            "Use LongBetMulti.load() or load_multi_npz() to load this model."
        )
    _require_current_precision_cache(meta)
    cutpoints = _validate_scalar_ordinal(data, config, meta) if config.outcome == "ordinal" else None
    meany = float(meta.pop("meany", 0.0))
    sdy = float(meta.pop("sdy", 1.0))
    has_chains = bool(meta.pop("has_chains", np.asarray(data["beta"]).ndim > 2))

    trend_coef = (
        jnp.asarray(data["trend_coef"]) if "trend_coef" in data.files else None
    )

    trace = LongBetTrace(
        mu_trace=_rebuild_trace(data, "mu", has_chains),
        nu_trace=_rebuild_trace(data, "nu", has_chains),
        **{k: jnp.asarray(data[k]) for k in _PARAM_KEYS},
        trend_coef=trend_coef,
        cutpoints=None if cutpoints is None else jnp.asarray(cutpoints),
        **{k: jnp.asarray(data[k]) for k in _OPTIONAL_KEYS if k in data.files},
    )
    if config.linear_trt and "delta" not in data.files:
        raise ValueError(f"Archive at {path} was fitted with linear_trt but has no delta draws")
    return trace, config, meany, sdy, meta
