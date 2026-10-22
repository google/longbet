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

"""MCMC diagnostics for LongBet, on the identified quantity.

This interface reports diagnostics for the **ATT series**. The treatment term
enters as a product ``beta_S * nu(X, S, t)``, and the likelihood can weakly
identify its factors separately. Proper priors still define a posterior scale;
that ridge does not excuse persistent chain disagreement. Diagnose the actual
reported effects, and inspect parameter traces when investigating poor mixing.

On the ``reliable`` verdict
---------------------------
The reference R implementation shipped an ``att_stability()`` whose ``reliable``
column was calibrated against *measured coverage of the XBART sweep sampler*.
That calibration does not transfer to a Metropolis sampler with real chains, and
recalibrating it requires a fresh coverage study.  Until that study exists this
function reports the numbers and **withholds the verdict**: ``reliable`` is
``None``.  ``ess_ok`` and ``rhat_ok`` are threshold checks on the numbers, and
are documented as exactly that -- not as evidence about coverage.

The split-half interval-width ratio is deliberately **not** reported.  It does
not detect the failure it appears to: on the reference simulation it read 0.93
for a fit that under-covered and 0.96 for one that did not.  Do not re-add it as
a verdict.

R-hat needs dispersed starts
----------------------------
R-hat compares within-chain to between-chain variance, and it only detects a
chain stuck in an unvisited region if the chains had a chance to start in
*different* regions.  Chains that begin at identical values and are separated
only by their random streams can agree with one another while all of them are
wrong.  LongBet therefore overdisperses the chains at initialisation, drawing
each one's ``beta``, ``gamma``, ``b0``, ``b1`` and variances from their priors
(see ``longbet._state.broadcast_to_chains``). Dispersed starts help reveal chain
disagreement, but neither agreement nor a small R-hat proves convergence.

Reading the output
------------------
Diagnostics are reported **per exposure time**, and the summary's ``rhat_max``
and ``ess_min`` are extremes over those.  Pass ``n_treated`` so the support
behind each one is reported alongside: in a staggered rollout only the earliest
adopters ever reach the longest exposures, so the thin end of the event-time
axis is estimated from a fraction of the panel. Support alone does not explain
poor mixing or excuse a failed diagnostic. The worst R-hat's location and
support land in ``summary`` as ``rhat_max_at_exposure`` and
``rhat_max_n_treated``.

Support is reported, never used to soften a threshold. ``ess_ok`` requires
finite bulk and tail ESS above the threshold at every exposure. MCSE describes
the posterior mean and uses mean ESS, not rank-normalized bulk ESS.

A useful cross-check is to compare per-chain ATT means and trace plots.
Agreement of means alone does not establish exploration of the full posterior,
and neither a likelihood ridge nor agreement on one estimand validates another.
"""

from __future__ import annotations

import warnings
from typing import Any, NamedTuple

import arviz as az
import numpy as np


class StabilityResult(NamedTuple):
    """Stability diagnostics for the ATT series.

    Attributes
    ----------
    summary
        Scalar diagnostics across exposure times.
    by_exposure
        Per-exposure-time diagnostics.
    reliable
        Always ``None``: see the module docstring. Present so callers have a
        stable field to read once a calibrated verdict exists.
    """

    summary: dict[str, Any]
    by_exposure: dict[str, np.ndarray]
    reliable: None


def _as_chains_draws(arr: np.ndarray) -> np.ndarray:
    """Coerce to ``(chains, draws)``."""
    a = np.asarray(arr)
    return a[None, :] if a.ndim == 1 else a


def compute_ess(
    draws: np.ndarray,
    method: str = "bulk",
    prob: float | tuple[float, float] | None = None,
) -> float | np.ndarray:
    """Effective sample size for the requested diagnostic.

    Parameters
    ----------
    draws
        Array of shape ``(draws,)`` or ``(chains, draws)``.
    method
        ``'bulk'`` (rank-normalized), ``'tail'``, or ``'mean'``.
    prob
        Tail probability or pair of quantile probabilities when ``method='tail'``
        (default 0.05, meaning 0.05 and 0.95). Use ``(.025, .975)`` to diagnose
        the endpoints of a 95% equal-tail interval.
    """
    x = _as_chains_draws(draws)
    if method == "tail":
        return az.ess(x, method="tail", prob=prob if prob is not None else 0.05)
    if method not in ("bulk", "mean"):
        raise ValueError(f"Unknown ESS method {method!r}; use 'bulk', 'tail', or 'mean'.")
    return az.ess(x, method=method)


def compute_rhat(draws: np.ndarray) -> float:
    """Rank-normalized split R-hat across chains; NaN for a single chain."""
    x = _as_chains_draws(draws)
    if x.shape[0] < 2:
        return float("nan")
    return float(az.rhat(x))


def att_stability(
    att_draws: np.ndarray,
    min_ess: float = 400.0,
    max_rhat: float = 1.01,
    alpha: float = 0.05,
    warn: bool = True,
    n_treated: np.ndarray | None = None,
) -> StabilityResult:
    """Diagnose sampling adequacy for the ATT.

    Parameters
    ----------
    att_draws
        ATT posterior draws, shaped ``(chains, draws, S)``, ``(draws, S)`` or
        ``(S, draws)``. Prefer the explicit 3-D form: the 2-D forms have to be
        disambiguated by comparing extents, which is ambiguous when the number
        of draws is close to the number of exposure times.
    min_ess
        Threshold that both bulk and tail ESS must meet at every exposure for
        ``ess_ok``. Undefined diagnostics fail the check. This is a rule of
        thumb, not a coverage guarantee.
    max_rhat
        R-hat threshold for ``rhat_ok`` (multi-chain fits only).
    alpha
        Tail probability used for the tail-ESS.
    warn
        Emit a ``UserWarning`` when a threshold is not met.
    n_treated
        Number of treated cells supporting each exposure time, if known.
        Reported alongside diagnostics, never used to explain away poor mixing
        or soften a threshold.

    Returns
    -------
    StabilityResult
        ``reliable`` is always ``None``; see the module docstring.
    """
    if hasattr(att_draws, "att_full"):
        att_draws = att_draws.att_full
    elif isinstance(att_draws, dict) and "att_full" in att_draws:
        att_draws = att_draws["att_full"]

    a = np.asarray(att_draws, dtype=np.float64)

    if a.ndim == 3:
        arr = a
    elif a.ndim == 2:
        # (S, draws) if the leading extent is the smaller one, else (draws, S).
        arr = (a.T if a.shape[0] < a.shape[1] else a)[None, ...]
    else:
        raise ValueError(f"Expected a 2-D or 3-D array of ATT draws, got shape {a.shape}")

    num_chains, num_draws, n_exposure = arr.shape
    if min(arr.shape) == 0:
        raise ValueError("ATT draws must have nonempty chain, draw and exposure axes.")
    total_draws = num_chains * num_draws

    ess_bulk = np.full(n_exposure, np.nan)
    ess_tail = np.full(n_exposure, np.nan)
    ess_mean = np.full(n_exposure, np.nan)
    rhat_vals = np.full(n_exposure, np.nan)
    mcse_vals = np.full(n_exposure, np.nan)

    for s in range(n_exposure):
        series = arr[..., s]
        if not np.all(np.isfinite(series)) or np.all(series == series.flat[0]):
            # A constant or non-finite series has no meaningful ESS.
            continue
        try:
            ess_bulk[s] = float(compute_ess(series, method="bulk"))
        except Exception:
            pass
        try:
            ess_tail[s] = float(compute_ess(series, method="tail", prob=alpha))
        except Exception:
            pass
        try:
            rhat_vals[s] = compute_rhat(series)
        except Exception:
            pass
        try:
            ess_mean[s] = float(compute_ess(series, method="mean"))
            mcse_vals[s] = float(az.mcse(series, method="mean"))
        except Exception:
            # Undefined diagnostics must stay undefined, rather than falling
            # back to the nominal number of retained draws.
            pass

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med_ess = float(np.nanmedian(ess_bulk))
        min_ess_obs = float(np.nanmin(ess_bulk))
        med_tail_ess = float(np.nanmedian(ess_tail))
        min_tail_ess = float(np.nanmin(ess_tail))
        max_rhat_obs = float(np.nanmax(rhat_vals)) if num_chains > 1 else float("nan")
        med_mcse = float(np.nanmedian(mcse_vals))

    # Locate the worst R-hat so callers can inspect that particular estimand.
    worst_idx = (
        int(np.nanargmax(rhat_vals)) if np.any(~np.isnan(rhat_vals)) else None
    )
    support = None
    if n_treated is not None:
        support_arr = np.asarray(n_treated, dtype=np.float64)
        if support_arr.size == n_exposure:
            support = support_arr

    ess_ok = bool(
        np.all(np.isfinite(ess_bulk)) and np.all(ess_bulk >= min_ess)
        and np.all(np.isfinite(ess_tail)) and np.all(ess_tail >= min_ess)
    )
    rhat_ok = (
        bool(np.all(np.isfinite(rhat_vals)) and np.all(rhat_vals <= max_rhat))
        if num_chains > 1
        else None
    )

    if warn:
        if not ess_ok:
            warnings.warn(
                f"att_stability(): ESS check failed. Every exposure needs finite "
                f"bulk and tail ESS of at least {min_ess:g}; minima among available "
                f"diagnostics are {min_ess_obs:.1f} (bulk) and {min_tail_ess:.1f} "
                "(tail). Inspect the affected traces and retain more effective "
                "draws before interpreting the reported uncertainty.",
                UserWarning,
                stacklevel=2,
            )
        if num_chains == 1:
            warnings.warn(
                "att_stability(): R-hat needs at least two chains. Set "
                "num_chains >= 2 to diagnose convergence rather than only "
                "autocorrelation.",
                UserWarning,
                stacklevel=2,
            )
        elif rhat_ok is False:
            where = ""
            if worst_idx is not None:
                where = f" at exposure time {worst_idx + 1}"
                if support is not None:
                    where += f", which {support[worst_idx]:.0f} treated cell(s) support"
            warnings.warn(
                f"att_stability(): maximum R-hat on the ATT is {max_rhat_obs:.3f}{where}, "
                f"and the check requires a finite R-hat <= {max_rhat:g} at every "
                "exposure. Inspect chain disagreement and both warmup and retained "
                "sampling. More burn-in alone does not guarantee adequate mixing; "
                "limited treated-cell support does not excuse a failed diagnostic.",
                UserWarning,
                stacklevel=2,
            )

    summary = {
        "num_chains": num_chains,
        "num_draws": num_draws,
        "total_draws": total_draws,
        "ess_median": med_ess,
        "ess_min": min_ess_obs,
        "ess_tail_median": med_tail_ess,
        "ess_tail_min": min_tail_ess,
        "rhat_max": max_rhat_obs,
        "mcse_median": med_mcse,
        "ess_ok": ess_ok,
        "rhat_ok": rhat_ok,
        "reliable": None,
        "verdict": "withheld",
        "verdict_note": (
            "No reliability verdict is issued. The reference implementation's "
            "threshold was calibrated against measured coverage of the XBART sweep "
            "sampler and does not transfer to this one. ess_ok and rhat_ok are "
            "threshold checks on the reported numbers, not evidence about coverage."
        ),
    }
    by_exposure = {
        "exposure": np.arange(1, n_exposure + 1),
        "ess_bulk": ess_bulk,
        "ess_tail": ess_tail,
        "ess_mean": ess_mean,
        "rhat": rhat_vals,
        "mcse": mcse_vals,
    }
    if worst_idx is not None:
        summary["rhat_max_at_exposure"] = worst_idx + 1
    if support is not None:
        by_exposure["n_treated"] = support
        summary["min_treated_cells"] = float(np.min(support))
        if worst_idx is not None:
            summary["rhat_max_n_treated"] = float(support[worst_idx])
    return StabilityResult(summary=summary, by_exposure=by_exposure, reliable=None)
