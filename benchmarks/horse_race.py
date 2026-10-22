#!/usr/bin/env python3
"""LongBet Horse Race Benchmark Script.

Measures JAX compilation, MCMC sampling, and prediction times on identical
staggered-adoption panel workloads.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import socket
import time
from pathlib import Path

# Prevent JAX plugin initialization issues when CUDA driver/library mismatch exists
if "JAX_PLATFORMS" not in os.environ:
    os.environ["JAX_PLATFORMS"] = "cpu"

import jax
import numpy as np

from longbet import LongBet, LongBetConfig, available_devices, get_att


def get_cpu_model() -> str:
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if "model name" in line:
                    return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return platform.processor() or "Unknown CPU"


def get_ram_gb() -> float:
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    return round(kb / (1024 * 1024), 1)
    except Exception:
        pass
    return 0.0


def generate_panel(n: int, T: int = 20, p: int = 5, seed: int = 42):
    """Staggered-adoption panel with an exposure-dependent effect."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, p)).astype(np.float32)
    t = np.arange(1, T + 1, dtype=np.float32)

    adopt = rng.integers(3, T - 1, size=n)
    adopt = np.where(rng.binomial(1, 0.5, size=n).astype(bool), adopt, 10**6)
    z = np.zeros((n, T), dtype=np.float32)
    for i in range(n):
        if adopt[i] <= T:
            z[i, adopt[i] - 1:] = 1.0
    s = np.where(z == 1, np.maximum(t[None, :] - (adopt[:, None] - 1), 0), 0)

    gamma = rng.normal(0.0, 0.5, size=n).astype(np.float32)
    y = (
        (0.5 * x[:, 0] + gamma)[:, None]
        + 0.05 * t[None, :]
        + 1.2 * np.sqrt(s)
        + rng.normal(0.0, 0.3, size=(n, T))
    ).astype(np.float32)
    return x, y, z, t


def run_benchmark(n: int, T: int, sweeps: int, burnin: int, trees: int, chains: int, device: str = "cpu"):
    x, y, z, t = generate_panel(n, T)
    config = LongBetConfig(
        num_sweeps=sweeps,
        num_burnin=burnin,
        num_chains=chains,
        num_trees_pr=trees,
        num_trees_trt=trees,
        device=device,
        random_seed=42,
    )

    # 1. Warmup / Compilation measurement
    probe = LongBetConfig(**{**config.to_dict(), "num_sweeps": 1, "num_burnin": 0})
    t0 = time.perf_counter()
    warm = LongBet(probe)
    warm.fit(y=y, x=x, z=z, t=t)
    jax.block_until_ready(warm.trace.beta)
    compile_time = time.perf_counter() - t0

    # 2. Sampling measurement
    model = LongBet(config)
    t0 = time.perf_counter()
    model.fit(y=y, x=x, z=z, t=t)
    jax.block_until_ready(model.trace.beta)
    sample_time = time.perf_counter() - t0
    cold_fit_time = compile_time + sample_time

    # 3. Prediction measurement
    t0 = time.perf_counter()
    pred = model.predict(x=x, z=z, t=t, summary_only=True)
    att = get_att(pred)
    predict_time = time.perf_counter() - t0

    # 4. MCMC diagnostics
    stab = pred.stability(warn=False)
    ess = float(stab.summary["ess_median"])
    total_time = cold_fit_time + predict_time
    total_iters = config.total_iterations()

    return {
        "N": n,
        "T": T,
        "M": n * T,
        "chains": chains,
        "sweeps": sweeps,
        "burnin": burnin,
        "trees": trees,
        "total_iterations": total_iters,
        "compile_time_s": round(compile_time, 3),
        "sample_time_s": round(sample_time, 3),
        "cold_fit_time_s": round(cold_fit_time, 3),
        "predict_time_s": round(predict_time, 3),
        "total_time_s": round(total_time, 3),
        "ms_per_iteration": round(1000.0 * sample_time / total_iters, 2),
        "att_ess_median": round(ess, 1),
        "att_rhat_max": round(float(stab.summary["rhat_max"]), 4),
        "seconds_per_effective_draw": round(total_time / ess, 4) if ess > 0 else None,
    }


def main():
    parser = argparse.ArgumentParser(description="LongBet Benchmark Runner")
    parser.add_argument("--sizes", type=int, nargs="+", default=[250, 1000])
    parser.add_argument("--periods", type=int, default=20)
    parser.add_argument("--sweeps", type=int, default=200)
    parser.add_argument("--burnin", type=int, default=100)
    parser.add_argument("--trees", type=int, default=20)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    host_info = {
        "hostname": socket.gethostname(),
        "cpu_model": get_cpu_model(),
        "cpu_logical_cores": os.cpu_count(),
        "ram_gb": get_ram_gb(),
        "os": f"{platform.system()} {platform.release()}",
        "jax_version": jax.__version__,
        "jax_devices": [str(d) for d in jax.devices()],
    }

    print("=" * 70)
    print(f"Host: {host_info['hostname']} | CPU: {host_info['cpu_model']}")
    print(f"Cores: {host_info['cpu_logical_cores']} threads | RAM: {host_info['ram_gb']} GB")
    print(f"JAX: {host_info['jax_version']} | Devices: {host_info['jax_devices']}")
    print("=" * 70)
    print(f"{'N':>6} {'M':>8} {'Compile':>9} {'Sample':>9} {'ms/iter':>8} {'Predict':>9} {'Total':>9} {'ESS':>7} {'s/eff-draw':>11}")
    print("-" * 79)

    results = []
    for n in args.sizes:
        res = run_benchmark(
            n=n,
            T=args.periods,
            sweeps=args.sweeps,
            burnin=args.burnin,
            trees=args.trees,
            chains=args.chains,
        )
        results.append(res)
        print(f"{res['N']:>6} {res['M']:>8,} "
              f"{res['compile_time_s']:>8.2f}s "
              f"{res['sample_time_s']:>8.2f}s "
              f"{res['ms_per_iteration']:>7.2f}ms "
              f"{res['predict_time_s']:>8.2f}s "
              f"{res['total_time_s']:>8.2f}s "
              f"{res['att_ess_median']:>7.1f} "
              f"{res['seconds_per_effective_draw']:>10.4f}s")

    payload = {
        "host": host_info,
        "config": vars(args) | ({"out": str(args.out)} if args.out else {}),
        "results": results,
    }

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2))
        print(f"\nSaved results to {args.out}")
    else:
        print("\nJSON Output:")
        print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
