"""Reference draws for the byte-identity test of the masked event-time Rank-K path.

``tests/test_rank_k_masked.py`` fits the same single-wave panels with the masked
implementation (``mask=None`` and an all-``True`` ``mask``) and compares every
array exactly against this file, which must be generated from the source tree
*before* the masked event-time changes (``google/longbet-jax`` @ 6f714d5), e.g.::

    JAX_PLATFORMS=cpu .venv/bin/python tests/data/make_masked_reference.py \\
        tests/data/masked_reference_6f714d5.npz

The script only uses the public API that exists in both versions. ``*.npz`` is
git-ignored, so add the reference with ``git add -f``.
"""

from __future__ import annotations

import os
import sys
import warnings

os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np

import longbet
from longbet import RankKLongBetConfig, fit_rank_k_longbet

PROP = dict(
    uncertainty_propagation=True, n_forecast_draws=3, max_forecast_draws=3, num_chains=2,
    draws_per_component=4, num_burnin=15, num_trees_mu=20, num_trees_tau_mode1=12,
    num_trees_tau_higher=8, random_seed=2026,
)
PLUG = dict(
    uncertainty_propagation=False, num_gfr=3, num_burnin=20, num_mcmc=30,
    num_trees_mu=20, num_trees_tau_mode1=12, num_trees_tau_higher=8, random_seed=2026,
)


def make_panel(kind: str, seed: int = 20261023):
    """Single-wave (``kind='single'``) or two-wave (``kind='staggered'``) panel with a forecast."""
    rng = np.random.default_rng(seed)
    N, T, P = 60, 12, 3
    x = rng.normal(size=(N, P))
    cohorts = np.full(N, T, dtype=np.int32)
    if kind == "single":
        cohorts[:32] = 4
    elif kind == "staggered":
        cohorts[:16] = 4
        cohorts[16:32] = 6
    else:
        raise ValueError(kind)
    z = (np.arange(T)[None, :] >= cohorts[:, None]).astype(np.int32)
    f = 0.3 * np.cumsum(rng.normal(size=(N, T)), axis=1) + 0.2 * x[:, :1]
    s = np.maximum(np.arange(T)[None, :] - cohorts[:, None] + 1, 0)
    tau = (0.8 + 0.4 * np.tanh(x[:, 0]))[:, None] * (0.85 ** np.maximum(s - 1, 0)) * (s > 0)
    y = 0.8 * f + z * tau + 0.3 * rng.normal(size=(N, T))
    f_draws = f[None, :, :] + 0.15 * rng.normal(size=(4, N, T))
    coords = np.column_stack([rng.uniform(30.0, 45.0, size=N), rng.uniform(-110.0, -80.0, size=N)])
    cluster_ids = rng.integers(0, 5, size=N)
    return dict(x=x, z=z, y=y, f=f, f_draws=f_draws, coords=coords, cluster_ids=cluster_ids, cohorts=cohorts)


def cases():
    """``(name, panel kind, config kwargs, fit kwargs builder)``."""
    yield "single_prop_k3_gp", "single", dict(PROP, rank_k=3), lambda p: dict(
        baseline_forecast=p["f"], baseline_forecast_draws=p["f_draws"], coords=p["coords"])
    yield "single_prop_k2_nofc", "single", dict(PROP, rank_k=2, use_forecast=False), lambda p: dict()
    yield "single_prop_k2_clip_floor", "single", dict(
        PROP, rank_k=2, tau_clip=(0.0, float("inf")), catt_floor=0.0, tau_prior_sd=0.5), lambda p: dict(
        baseline_forecast=p["f"], baseline_forecast_draws=p["f_draws"])
    yield "single_prop_k2_cluster", "single", dict(PROP, rank_k=2), lambda p: dict(
        baseline_forecast=p["f"], baseline_forecast_draws=p["f_draws"], cluster_ids=p["cluster_ids"])
    yield "single_prop_k2_t0panel", "single", dict(PROP, rank_k=2, use_spatial_gp=False), lambda p: dict(
        baseline_forecast=p["f"], baseline_forecast_draws=p["f_draws"], t_0_panel=3)
    yield "single_plugin_k2", "single", dict(PLUG, rank_k=2), lambda p: dict(baseline_forecast=p["f"])
    yield "single_plugin_k3_nofc", "single", dict(PLUG, rank_k=3, use_forecast=False), lambda p: dict()
    yield "stag_prop_k2_gp", "staggered", dict(PROP, rank_k=2), lambda p: dict(
        baseline_forecast=p["f"], baseline_forecast_draws=p["f_draws"], coords=p["coords"])


def result_arrays(res) -> dict[str, np.ndarray]:
    """The arrays compared bit-for-bit."""
    out = {
        "att_draws": np.asarray(res.att_draws),
        "catt_draws": np.asarray(res.catt_draws),
        "catt_hat": np.asarray(res.catt_hat),
        "theta_f_draws_0": np.asarray(res.theta_f_draws[0]),
        "sigma_modes": np.asarray(res.sigma_modes),
        "att_summary": np.asarray([res.att_hat, res.att_se, res.att_ci_low, res.att_ci_high]),
        "att_var": np.asarray([res.att_var_within, res.att_var_between, res.att_var_slope]),
        "rhat_ess": np.asarray([res.rhat_att, res.ess_att]),
    }
    for k, td in enumerate(res.mode_tau_draws):
        out[f"mode_tau_draws_{k}"] = np.asarray(td)
    for g, basis in sorted(res.bases.items()):
        out[f"basis_{g}"] = np.asarray(basis)
    raw = getattr(res, "catt_draws_raw", None)
    if raw is not None:
        out["catt_draws_raw"] = np.asarray(raw)
    return out


def all_arrays() -> dict[str, np.ndarray]:
    out = {}
    for name, kind, cfg_kw, kw_builder in cases():
        panel = make_panel(kind)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = fit_rank_k_longbet(
                panel["y"], panel["z"], panel["x"], cfg=RankKLongBetConfig(**cfg_kw), **kw_builder(panel)
            )
        for key, val in result_arrays(res).items():
            out[f"{name}/{key}"] = val
    return out


if __name__ == "__main__":
    path = sys.argv[1]
    arrays = all_arrays()
    arrays["_source"] = np.array(longbet.__file__)
    arrays["_version"] = np.array(longbet.__version__)
    np.savez_compressed(path, **arrays)
    print(f"wrote {len(arrays)} arrays from {longbet.__file__} ({longbet.__version__}) to {path}")
