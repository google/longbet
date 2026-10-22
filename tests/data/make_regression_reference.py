"""Reference draws for the byte-identity regression test (``tests/test_lag_kernel_regression.py``).

Run with ``main``'s source first on ``sys.path`` to record, e.g.::

    PYTHONPATH=~/projects/longbet-jax/src JAX_PLATFORMS=cpu .venv/bin/python \
        tests/data/make_regression_reference.py tests/data/main_f839357_draws.npz

The test fits the same panels with the current source and compares every
array exactly.
"""
from __future__ import annotations

import sys
import warnings

import numpy as np

import longbet
from longbet import LongBet, LongBetConfig


def panels():
    rng = np.random.default_rng(20261010)
    N, T = 24, 6
    x = rng.normal(size=(N, 3)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[:8, 3:] = 1.0
    z[8:14, 4:] = 1.0
    y = (0.4 * x[:, :1] + 0.3 * rng.normal(size=(N, 1))
         + 0.6 * rng.normal(size=(N, T)) + 0.5 * z).astype(np.float64)
    y_nan = y.copy()
    y_nan[rng.random((N, T)) < 0.1] = np.nan
    y_bin = (y > np.median(y)).astype(np.float64)
    base = dict(num_sweeps=12, num_burnin=10, num_chains=2, num_trees_pr=6, num_trees_trt=8,
                trend_num_features=8, random_seed=10007)
    return {
        "default": (dict(base), y, x, z),
        "v2_unanchored": (dict(base, linear_trt=True, learn_beta_scale=True, ar1_errors=True,
                               interweave_gamma=True, split_calendar_mu=False), y, x, z),
        "no_trend_block": (dict(base, use_trend_block=False, split_calendar_mu=True), y, x, z),
        "unbalanced": (dict(base), y_nan, x, z),
        "binary": (dict(base, outcome="binary"), y_bin, x, z),
    }


def draws_for(name, spec):
    cfg_kw, y, x, z = spec
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = LongBet(LongBetConfig(**cfg_kw)).fit(y, x, z)
    tr = model.trace
    out = {
        f"{name}/beta": np.asarray(tr.beta), f"{name}/gamma": np.asarray(tr.gamma),
        f"{name}/sigma2": np.asarray(tr.sigma2), f"{name}/sigma_gamma2": np.asarray(tr.sigma_gamma2),
        f"{name}/trend_coef": np.asarray(tr.trend_coef),
        f"{name}/mu_leaf": np.asarray(tr.mu_trace.leaf_tree), f"{name}/nu_leaf": np.asarray(tr.nu_trace.leaf_tree),
        f"{name}/mu_split": np.asarray(tr.mu_trace.split_tree), f"{name}/nu_split": np.asarray(tr.nu_trace.split_tree),
    }
    for opt in ("delta", "beta_scale2", "ar1_phi"):
        v = getattr(tr, opt, None)
        if v is not None:
            out[f"{name}/{opt}"] = np.asarray(v)
    pred = model.predict(x, z, summary_only=True)
    out[f"{name}/att_full"] = np.asarray(pred.att_full)
    return out


def all_draws():
    out = {}
    for name, spec in panels().items():
        out.update(draws_for(name, spec))
    return out


if __name__ == "__main__":
    path = sys.argv[1]
    arrays = all_draws()
    arrays["_source"] = np.array(longbet.__file__)
    np.savez_compressed(path, **arrays)
    print(f"wrote {len(arrays)} arrays from {longbet.__file__} to {path}")
