"""Generate bit-identity golden fixtures for the Rank-K estimator (paper repo or longbet-jax).

Usage:
  python make_golden.py --repo paper --src-first /path/to/orig_dir --out /path/to/fixtures
  python make_golden.py --repo pkg --out /path/to/fixtures

The cases cover only paper-default behaviour (no tau prior / clip / catt_floor, every cohort with
L_g >= K), which Design Invariant 1 of the remediation plan promises to keep bit-for-bit identical.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
from pathlib import Path

os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np


def make_panel(kind: str, seed: int = 2026):
    rng = np.random.default_rng(seed)
    N, T, P = 60, 12, 3
    x = rng.normal(size=(N, P))
    cohorts = np.full(N, T, dtype=np.int32)
    if kind == "stag_unequal":
        cohorts[:10] = 4
        cohorts[10:20] = 5
        cohorts[20:30] = 6
    elif kind == "balanced":
        cohorts[:30] = 4
    else:
        raise ValueError(kind)
    z = (np.arange(T)[None, :] >= cohorts[:, None]).astype(np.int32)
    f = 0.3 * np.cumsum(rng.normal(size=(N, T)), axis=1) + 0.2 * x[:, :1]
    tau = (0.8 + 0.4 * np.tanh(x[:, 0]))[:, None] * (1.0 - np.exp(-np.maximum(np.arange(T)[None, :] - cohorts[:, None] + 1, 0) / 3.0))
    y = 0.8 * f + z * tau + 0.3 * rng.normal(size=(N, T))
    f_draws = f[None, :, :] + 0.15 * rng.normal(size=(4, N, T))
    coords = np.column_stack([rng.uniform(30.0, 45.0, size=N), rng.uniform(-110.0, -80.0, size=N)])
    return dict(x=x, z=z, y=y, f=f, f_draws=f_draws, coords=coords, cohorts=cohorts)


def cases():
    prop = dict(
        uncertainty_propagation=True, n_forecast_draws=3, max_forecast_draws=3, num_chains=2,
        draws_per_component=4, num_burnin=15, num_trees_mu=20, num_trees_tau_mode1=12,
        num_trees_tau_higher=8, random_seed=2026,
    )
    plug = dict(
        uncertainty_propagation=False, num_gfr=3, num_burnin=20, num_mcmc=30,
        num_trees_mu=20, num_trees_tau_mode1=12, num_trees_tau_higher=8, random_seed=2026,
    )
    yield "stag_unequal_prop_k2_gp", "stag_unequal", dict(prop, rank_k=2), True
    yield "stag_unequal_prop_k3_nogp", "stag_unequal", dict(prop, rank_k=3, use_spatial_gp=False), False
    yield "stag_unequal_prop_k1_gp", "stag_unequal", dict(prop, rank_k=1), True
    yield "balanced_prop_k2_gp", "balanced", dict(prop, rank_k=2), True
    yield "stag_unequal_plugin_k2", "stag_unequal", dict(plug, rank_k=2), False
    yield "balanced_plugin_k3", "balanced", dict(plug, rank_k=3), False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", choices=["paper", "pkg"], required=True)
    ap.add_argument("--src-first", default=None, help="directory inserted first on sys.path (original sources)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--check", action="store_true", help="compare against existing fixtures instead of writing")
    args = ap.parse_args()

    if args.src_first:
        sys.path.insert(0, args.src_first)
    if args.repo == "paper":
        from rank_k_longbet import RankKLongBetConfig, fit_rank_k_longbet  # type: ignore
        import rank_k_longbet as mod
    else:
        from longbet import RankKLongBetConfig, fit_rank_k_longbet  # type: ignore
        import longbet._rank_k as mod
    import jax

    print("module:", mod.__file__)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = {"jax": jax.__version__, "platform": platform.platform(), "machine": platform.machine(), "repo": args.repo}
    n_fail = 0
    for name, kind, cfg_kw, use_coords in cases():
        if cfg_kw.get("uncertainty_propagation") is False:
            try:
                import stochtree  # noqa: F401
            except ImportError:
                print("skip (no stochtree):", name)
                continue
        p = make_panel(kind)
        cfg = RankKLongBetConfig(**cfg_kw)
        kw = dict(baseline_forecast=p["f"], cohorts=p["cohorts"])
        if cfg_kw.get("uncertainty_propagation"):
            kw["baseline_forecast_draws"] = p["f_draws"]
        if use_coords:
            kw["coords"] = p["coords"]
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = fit_rank_k_longbet(p["y"], p["z"], p["x"], cfg=cfg, **kw)
        arrays = {
            "att_draws": np.asarray(res.att_draws),
            "catt_draws": np.asarray(res.catt_draws),
            "theta_f_draws_0": np.asarray(res.theta_f_draws[0]),
            "sigma_modes": np.asarray(res.sigma_modes),
            "att_var": np.asarray([res.att_var_within, res.att_var_between, res.att_var_slope]),
            "rhat_ess": np.asarray([res.rhat_att, res.ess_att]),
        }
        for k, td in enumerate(res.mode_tau_draws):
            arrays[f"mode_tau_draws_{k}"] = np.asarray(td)
        for k, g in enumerate(res.gp_hyper_draws):
            arrays[f"s_zeta_{k}"] = np.asarray(g["s_zeta"])
            arrays[f"ell_{k}"] = np.asarray(g["ell"])
        path = out / f"{name}.npz"
        if args.check:
            ref = np.load(path)
            bad = []
            for key, val in arrays.items():
                r = ref[key]
                if r.shape != val.shape or not np.array_equal(r, val, equal_nan=True):
                    diff = float(np.max(np.abs(r - val))) if r.shape == val.shape else float("nan")
                    bad.append((key, diff))
            status = "OK" if not bad else f"MISMATCH {bad}"
            if bad:
                n_fail += 1
            print(f"{name}: {status}")
        else:
            np.savez_compressed(path, **arrays)
            print("wrote", path, {k: v.shape for k, v in arrays.items() if k in ("att_draws", "catt_draws")})
    (out / "metadata.json").write_text(json.dumps(meta, indent=2))
    if args.check and n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
