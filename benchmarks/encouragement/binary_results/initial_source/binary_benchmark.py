"""Reproducible binary joint-model calibration; failures are first-class records.

Example: python binary_benchmark.py --output binary_results/cells --replicates 100
         python binary_benchmark.py --output binary_results/forest --models forest_monotone forest_unrestricted --scenarios strong weak --replicates 2 --draws 4000 --burnin 1000
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import time
import traceback
import warnings

import numpy as np

try:
    from .binary_cells import effects, fit_cells, identification_bounds
    from .binary_dgp import SCENARIOS, make_binary_data, targets
    from .binary_forest import BinaryConfig, fit_binary
except ImportError:
    from binary_cells import effects, fit_cells, identification_bounds
    from binary_dgp import SCENARIOS, make_binary_data, targets
    from binary_forest import BinaryConfig, fit_binary

MODELS = ("cells_joint", "cells_independent", "cells_monotone", "cells_null", "forest_unrestricted", "forest_monotone")
BOUND_NAMES = ("complier_encouragement", "treatment_z0", "treatment_z1")


def portable(value):
    if isinstance(value, dict):
        return {k: portable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [portable(v) for v in value]
    if isinstance(value, np.generic):
        return portable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None if np.isnan(value) else ("inf" if value > 0 else "-inf")
    return value


def write_json(path, value):
    path = Path(path)
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(portable(value), indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def diagnostics(blocks):
    """Diagnose each stored effect/function/parameter; never drop invalid draws."""
    import arviz as az
    import xarray as xr
    names, columns = [], []
    for name, value in blocks.items():
        for index in np.ndindex(value.shape[2:]):
            names.append(name + ("[" + ",".join(map(str, index)) + "]" if index else ""))
            columns.append(value[(slice(None), slice(None), *index)])
    values = np.stack(columns, axis=-1)
    finite = np.isfinite(values).all(axis=(0, 1))
    varied = np.ptp(values, axis=(0, 1)) > 0
    eligible = finite & varied
    stats = {}
    if eligible.any():
        ds = xr.Dataset({"v": (("chain", "draw", "quantity"), values[..., eligible])})
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            stats = dict(rhat=np.asarray(az.rhat(ds).v),
                         ess_bulk=np.asarray(az.ess(ds, method="bulk").v),
                         ess_tail=np.asarray(az.ess(ds, method="tail", prob=.025).v),
                         mcse_mean=np.asarray(az.mcse(ds, method="mean").v))
            for label, quantile in (("lower", .025), ("median", .5), ("upper", .975)):
                stats["mcse_" + label] = np.asarray(az.mcse(ds, method="quantile", prob=quantile).v)
    result, j = {}, 0
    for i, name in enumerate(names):
        value = dict(undefined_fraction=float(np.mean(~np.isfinite(values[..., i]))))
        if eligible[i]:
            value.update({key: float(val[j]) for key, val in stats.items()})
            if name == "wald":
                del value["mcse_mean"]  # raw ratio means need not exist
            value["pass"] = bool(value["rhat"] <= 1.01 and value["ess_bulk"] >= 400 and value["ess_tail"] >= 400)
            value["status"] = "sampled"
            j += 1
        else:
            value.update(status="constant" if finite[i] else "undefined_draws", **{"pass": False})
        result[name] = value
    return result


def summarize(archive, truth, *, full_diagnostics=True):
    draw = effects(archive)
    summaries = {}
    for name, value in draw.items():
        available = np.isfinite(value).all()
        lo, median, hi = np.quantile(value, [.025, .5, .975]) if available else [np.nan] * 3
        summaries[name] = dict(lower=lo, median=median, upper=hi, width=hi-lo,
                               mean=value.mean() if available and name != "wald" else None,
                               coverage={target: bool(lo <= row[name] <= hi) if np.isfinite(row[name]) and available else None
                                         for target, row in truth.items()})
    summaries["posterior_itt_covariance"] = np.cov(np.stack([draw["itt_y"].ravel(), draw["itt_d"].ravel()]))
    summaries["negative_stage_probability"] = np.mean(draw["itt_d"] < 0)
    summaries["stage_below_005_probability"] = np.mean(np.abs(draw["itt_d"]) < .05)
    blocks = {**draw, **{name: archive[name] for name in ("p", "q", "r")}}
    for name in ("eta", "leaves", "rules"):
        if name in archive:
            blocks[name] = archive[name]
    bounds = {}
    if str(archive["first_stage"]) in ("monotone", "null"):
        for delta in (0., .05, .1, 1.):
            value = identification_bounds(archive["p"], archive["q"], archive["weights"], delta=delta)
            bound = dict(compatible_fraction=value["compatible"].mean(), defined_fraction=value["defined"].mean(),
                         available_fraction=value["available"].mean())
            # These are posterior quantiles of bound endpoints, not a posterior
            # for a point inside the set or a guaranteed 95% confidence region.
            for name in BOUND_NAMES:
                endpoints = value[name]
                blocks[f"bounds_delta_{delta}_{name}"] = endpoints
                if value["available"].all():
                    endpoint_quantiles = np.quantile(endpoints, [.025, .5, .975], axis=(0, 1))
                    envelope = np.array([endpoint_quantiles[0, 0], endpoint_quantiles[2, 1]])
                    bound[name] = dict(endpoint_quantiles=endpoint_quantiles, outer_envelope=envelope,
                        truth_in_envelope={target: bool(envelope[0] <= row[name] <= envelope[1])
                                           if np.isfinite(row[name]) else None for target, row in truth.items()})
                else:
                    bound[name] = dict(status="withheld: some draws have empty sets or zero complier mass")
            bounds[str(delta)] = bound
    summaries["bounds"] = bounds
    # IID cell posterior simulations do not need an MCMC convergence claim.
    # Full draw diagnostics are still available on every forest run and cell
    # audit replicate; other IID replicates use exact independence and record
    # the Monte Carlo mean/covariance errors directly.
    summaries["diagnostics"] = diagnostics(blocks) if full_diagnostics else {}
    summaries["diagnostic_mode"] = "all_quantities" if full_diagnostics else "iid_exact_draws"
    if "eta" not in archive:
        summaries["iid_mean_mcse"] = {k: np.std(v, ddof=1) / np.sqrt(v.size) for k, v in draw.items() if k != "wald"}
    return summaries


def reference(data, truth, *, seed, permutations):
    from longbet import encouragement_effects
    try:
        from .randomization_ar import randomization_ar
    except ImportError:
        from randomization_ar import randomization_ar
    y, d, z = [data[k][:, None] for k in ("y", "d", "z")]
    row = encouragement_effects(y, d, z).iloc[0].to_dict()
    record = dict(row=row, targets=truth, itt_covariance=np.array([[row["itt_y_se"]**2, row["itt_y_d_cov"]],
                                                                   [row["itt_y_d_cov"], row["itt_d_se"]**2]]))
    coverage, ar = {}, {}
    for label, target in truth.items():
        coverage[label] = {name: bool(row[name + "_lower"] <= target[name] <= row[name + "_upper"]) for name in ("itt_y", "itt_d")}
        beta = target["wald"]
        coverage[label]["wald"] = (any(row[f"wald_lower_{k}"] <= beta <= row[f"wald_upper_{k}"] for k in (1, 2)) if np.isfinite(beta) else None)
        if np.isfinite(beta) and permutations:
            result = randomization_ar(y, d, z, beta=[beta], method="monte_carlo", seed=seed, permutations=permutations)
            ar[label] = dict(row=result.table.iloc[0].to_dict(), metadata=result.metadata)
    record.update(coverage=coverage, randomization_ar_at_truth=ar)
    return record


def run(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    config = BinaryConfig(trees=args.trees)
    specification = dict(n=args.n, chains=4, draws=args.draws, burnin=args.burnin, prior=args.prior,
                         seed=args.seed, forest=asdict(config), max_candidates=args.max_candidates,
                         permutations=args.permutations)
    source_files = list(Path(__file__).parent.glob("binary_*.py"))
    source_files += [Path(__file__).with_name(name) for name in ("direct_smooth_candidate.py", "randomization_ar.py")]
    source_files += [Path(__file__).parents[2] / "src/longbet" / name for name in ("_encourage.py", "_diagnostics.py")]
    source = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    manifest = dict(specification=specification, sources=source, python=platform.python_version(),
                    platform=platform.platform(), packages={name: importlib.metadata.version(name)
                    for name in ("numpy", "scipy", "arviz", "jax", "bartz", "pandas")})
    manifest_path = output / "manifest.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise ValueError("Output belongs to another configuration/source version; use a fresh directory.")
    write_json(manifest_path, manifest)
    for scenario in args.scenarios:
        for replicate in range(args.start, args.start + args.replicates):
            # Shared data across model comparisons and separate deterministic
            # population, assignment, sampler and randomization-test streams.
            root = np.random.SeedSequence([args.seed, SCENARIOS.index(scenario), replicate])
            population_seed, assignment_seed, sampler_seed, ar_seed = [int(v) for v in root.generate_state(4)]
            directory = output / f"{scenario}_{replicate:04d}"
            directory.mkdir(exist_ok=True)
            data = make_binary_data(scenario, n=args.n, population_seed=population_seed, assignment_seed=assignment_seed)
            np.savez_compressed(directory / "data.npz", **data)
            truth = targets(data)
            if not (directory / "reference.json").exists():
                write_json(directory / "reference.json", reference(data, truth, seed=ar_seed, permutations=args.permutations))
            for model in args.models:
                path = directory / (model + ".json")
                if path.exists():
                    continue
                started = time.monotonic()
                record = dict(model=model, scenario=scenario, replicate=replicate, specification=specification,
                              sampler_seed=sampler_seed, targets=truth, status="running")
                write_json(path, record)
                try:
                    input_data = [data[k] for k in ("y", "d", "z", "x")]
                    if model.startswith("forest_"):
                        archive = fit_binary(*input_data, seed=sampler_seed, chains=4, draws=args.draws, burnin=args.burnin,
                                             first_stage=model.removeprefix("forest_"), config=config)
                    else:
                        stage = {"cells_joint": "unrestricted", "cells_independent": "unrestricted", "cells_monotone": "monotone", "cells_null": "null"}[model]
                        archive = fit_cells(*input_data, seed=sampler_seed, chains=4, draws=args.draws, first_stage=stage,
                                            prior=args.prior, max_candidates=args.max_candidates, independent_reduced_forms=model == "cells_independent")
                    np.savez_compressed(directory / (model + ".npz"), **archive)
                    record["summary"] = summarize(archive, truth, full_diagnostics=model.startswith("forest_") or replicate == 0)
                    record["status"] = "complete"
                except Exception:
                    record.update(status="failed", traceback=traceback.format_exc())
                record["seconds"] = time.monotonic() - started
                write_json(path, record)
                print(f"{scenario}/{replicate}/{model}: {record['status']} ({record['seconds']:.1f}s)", flush=True)
    aggregate(output)


def rate(values):
    valid = np.asarray([v for v in values if v is not None], dtype=float)
    if not len(valid):
        return dict(rate=None, mcse=None, n=0)
    estimate = valid.mean()
    # Wilson interval remains informative for all-success/all-failure pilots.
    denominator = 1 + 1.96**2 / len(valid)
    center = (estimate + 1.96**2 / (2 * len(valid))) / denominator
    half = 1.96 * np.sqrt(estimate * (1 - estimate) / len(valid) + 1.96**2 / (4 * len(valid)**2)) / denominator
    return dict(rate=estimate, mcse=np.sqrt(estimate * (1 - estimate) / len(valid)), n=len(valid), wilson_95=[center-half, center+half])


def aggregate(output):
    output = Path(output)
    groups = {}
    for path in sorted(output.glob("*/*.json")):
        if path.stem == "reference":
            continue
        value = json.loads(path.read_text())
        groups.setdefault((value["scenario"], value["model"]), []).append(value)
    results = []
    for (scenario, model), runs in sorted(groups.items()):
        complete = [r for r in runs if r["status"] == "complete"]
        record = dict(scenario=scenario, model=model, attempts=len(runs), complete=len(complete),
                      failures=sum(r["status"] == "failed" for r in runs), pending=sum(r["status"] == "running" for r in runs),
                      estimates_condition_on_completed_fits=True,
                      failure_reasons=[r["traceback"] for r in runs if r["status"] == "failed"])
        record["targets"] = {}
        for label in ("conditional", "finite"):
            if not complete:
                continue
            errors = np.array([[r["summary"][q]["mean"] - r["targets"][label][q] for q in ("itt_y", "itt_d")] for r in complete])
            covariance = np.array([r["summary"]["posterior_itt_covariance"] for r in complete])
            centered = errors - errors.mean(axis=0)
            products = centered[:, :, None] * centered[:, None, :]
            record["targets"][label] = dict(bias=errors.mean(axis=0), mean_posterior_covariance=covariance.mean(axis=0),
                empirical_error_covariance=np.cov(errors.T) if len(complete) > 1 else None,
                empirical_covariance_mcse=products.std(axis=0, ddof=1) / np.sqrt(len(complete)) if len(complete) > 1 else None,
                coverage={q: rate([r["summary"][q]["coverage"][label] for r in complete]) for q in ("itt_y", "itt_d", "wald")})
        if complete:
            record["mean_interval_width"] = {q: np.mean([r["summary"][q]["width"] for r in complete if r["summary"][q]["width"] is not None])
                                              if any(r["summary"][q]["width"] is not None for r in complete) else None for q in ("itt_y", "itt_d", "wald")}
            record["negative_stage_probability"] = np.mean([r["summary"]["negative_stage_probability"] for r in complete])
            record["stage_below_005_probability"] = np.mean([r["summary"]["stage_below_005_probability"] for r in complete])
            # Count all diagnosed failures; undefined bound endpoints are
            # identified separately, never disguised as sampler failures.
            diagnosed = [r for r in complete if r["summary"]["diagnostics"]]
            record["diagnosed_runs"] = len(diagnosed)
            record["diagnostic_failures"] = [{"replicate": r["replicate"], "quantities": {k: d for k, d in r["summary"]["diagnostics"].items()
                if not d["pass"] and d["status"] == "sampled"}} for r in diagnosed]
            record["bound_availability"] = {delta: dict(compatible_fraction=np.mean([r["summary"]["bounds"][delta]["compatible_fraction"] for r in complete]),
                fully_available_datasets=sum(r["summary"]["bounds"][delta]["available_fraction"] == 1 for r in complete))
                for delta in complete[0]["summary"]["bounds"]}
        results.append(record)
    reference_groups = {}
    for path in sorted(output.glob("*/reference.json")):
        scenario = path.parent.name.rsplit("_", 1)[0]
        reference_groups.setdefault(scenario, []).append(json.loads(path.read_text()))
    references = []
    for scenario, runs in sorted(reference_groups.items()):
        record = dict(scenario=scenario, datasets=len(runs), targets={})
        for label in ("conditional", "finite"):
            errors = np.array([[r["row"][q] - r["targets"][label][q] for q in ("itt_y", "itt_d")] for r in runs])
            record["targets"][label] = dict(bias=errors.mean(axis=0), mean_estimated_covariance=np.mean([r["itt_covariance"] for r in runs], axis=0),
                empirical_error_covariance=np.cov(errors.T) if len(runs) > 1 else None,
                coverage={q: rate([r["coverage"][label][q] for r in runs]) for q in ("itt_y", "itt_d", "wald")},
                randomization_ar_acceptance_at_truth=rate([r["randomization_ar_at_truth"][label]["row"]["accepted"]
                    if label in r["randomization_ar_at_truth"] else None for r in runs]))
        references.append(record)
    write_json(output / "aggregate.json", dict(models=results, reference=references))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS[:4]))
    parser.add_argument("--scenarios", nargs="+", choices=SCENARIOS, default=list(SCENARIOS))
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--n", type=int, default=240)
    parser.add_argument("--seed", type=int, default=19037)
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--burnin", type=int, default=1000)
    parser.add_argument("--trees", type=int, default=3)
    parser.add_argument("--prior", type=float, default=1.)
    parser.add_argument("--max-candidates", type=int, default=500000)
    parser.add_argument("--permutations", type=int, default=499)
    parser.add_argument("--cpus", default="0,1")
    args = parser.parse_args()
    if args.replicates < 1 or args.start < 0 or args.draws < 1 or args.burnin < 0 or args.permutations < 0:
        parser.error("Invalid iteration, replicate or permutation count.")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {int(v) for v in args.cpus.split(",")})
    run(args)


if __name__ == "__main__":
    main()
