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

"""General A-arm coherent Mode-1 and trajectory reconciliation for Rank-K BCF posteriors.

Generalizes :func:`longbet._rank_k.reconcile_nested_3arm_draws` from a 3-arm nested chain
to an arbitrary number of arms ``A >= 2`` and any connected set of ``P >= A - 1`` fitted
pairwise contrasts under chain or poset treatment topologies:

1. **Graph validation and spanning tree**: verifies that the fitted contrasts form a
   connected graph on ``arms`` (with ``arms[0]`` as control), constructs the spanning tree
   ``T`` (adjacent chain edges for ``topology="chain"``, Hasse diagram edges for
   ``topology="poset"``), requires every tree edge to have a fitted contrast, and
   identifies redundant contrasts.
2. **Mode-1 likelihood (Proposition 1 & Proposition 2)**: forms the signed incidence
   matrix ``E`` and reduced cut matrix ``B = E[:, 1:]``, computes the Euclidean cut-space
   projection ``mu_hat = (B' B)^{-1} B' y_hat`` and cycle discrepancy diagnostic
   ``r_cycle = ||y_hat - B mu_hat||_2 / max(||y_hat||_2, 1e-12)``, checks whether the
   contrast graph is bipartite via 2-coloring, and recovers the arm sampling variances
   ``v in R_+^A`` via bounded NNLS (or sample-size-informed ridge-NNLS with a
   ``UserWarning`` when the contrast graph is bipartite).
3. **Prior families**: supports chain bundle priors with ``Dirichlet(alpha)`` increments,
   unique-top poset order-polytope priors via uniform linear extensions, independent
   truncated Gaussian/half-normal priors with optional poset monotonicity, and flat
   priors.
4. **Dimension-adaptive posterior sampling**: uses generalized defensive-mixture SIR for
   ``d = A - 1 <= 3`` with an effective sample size gate ``ess >= max(100, 0.02 * n_prop)``
   and automatic fallback (or for ``d >= 4`` / ``sampler="slice"``) to a vectorized
   multi-chain analytical-polyhedron hit-and-run slice sampler reporting split R-hat and
   bulk ESS.
5. **Exact trajectory cycle coherence (Proposition 3)**: projects each spanning-tree edge
   increment draw through that edge's fit via ``_project_cohort_damped_trajectory`` and
   synthesizes every contrast trajectory as the signed path sum on ``T`` so cycle
   identities hold to machine precision across all triples, draws, and periods.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, replace
import math
from typing import Any
import warnings

import numpy as np
from scipy import optimize as sp_optimize
from scipy import special as sp_special
from scipy import stats as sp_stats

from longbet._rank_k import (
    _BUNDLE_PRIOR_MODES,
    _bundle_prior_logpdf,
    _draw_bundle_prior,
    _get_rhat_ess,
    _infer_panel_masks_from_result,
    _project_cohort_damped_trajectory,
    _warn_rank_k,
    Nested3ArmReconciliationResult,
    RankKLongBetResult,
    reconcile_nested_3arm_draws,
)

__all__ = [
    "ContrastReconciliationResult",
    "MultiArmReconciliationResult",
    "RankKDesignWarning",
    "SamplerDiagnostics",
    "reconcile_multi_arm_draws",
]

_INDEPENDENT_PRIOR_MODES = (
    "independent",
    "independent_halfnormal",
    "independent_normal",
)
# Spike-and-slab bundle priors: tau_top ~ pi delta_0 + (1 - pi) slab, increments split by the
# slab's Dirichlet only on the slab component (v0.10 item B1).
_SPIKE_PRIOR_MODES = ("bundle_spike_halfnormal",)
_SPIKE_SLAB_MODE = {"bundle_spike_halfnormal": "bundle_halfnormal"}
_ALL_PRIOR_MODES = (*_BUNDLE_PRIOR_MODES, *_SPIKE_PRIOR_MODES, *_INDEPENDENT_PRIOR_MODES, "flat")
_TOPOLOGIES = ("chain", "poset")
_SAMPLERS = ("auto", "sir", "slice")
_BIPARTITE_POLICIES = ("allocation", "conservative", "marginalize")
_PROJECTIONS = ("ols", "huber")
_HUBER_K = 1.345
# Prior family used when ``prior_mode=None`` on a non-bipartite / bipartite contrast graph.
# The bipartite default only moves away from the bundle prior if the Tier A v2 C13' gate
# (reconciled RMSE <= naive, Cov95 >= 0.93) is met by the independent family.
_DEFAULT_PRIOR_MODE = "bundle_halfnormal"
_BIPARTITE_DEFAULT_PRIOR_MODE = "bundle_halfnormal"
_N_NULL_GRID = 9


_NORM_PDF_LOGC = np.log(np.sqrt(2 * np.pi))


def _fast_bundle_prior_logpdf(prior_mode: str, sd_ac: float, t: np.ndarray) -> np.ndarray:
    """Bit-identical, scipy-free fast path of ``_bundle_prior_logpdf`` for ``bundle_halfnormal``.

    Evaluates ``log 2 + (-(t / sd)^2 / 2 - log sqrt(2 pi)) - log sd`` in exactly the
    floating-point order of ``scipy.stats.norm.logpdf`` (about 7x faster on the small arrays
    of the slice sampler); other families defer to :func:`_bundle_prior_logpdf`.
    """
    if prior_mode != "bundle_halfnormal":
        return _bundle_prior_logpdf(prior_mode, sd_ac, t)
    t_arr = np.asarray(t, dtype=np.float64)
    x_arr = t_arr / sd_ac
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        val = np.log(2.0) + (-x_arr**2 / 2.0 - _NORM_PDF_LOGC) - np.log(sd_ac)
    return np.where(t_arr > 0.0, val, -np.inf)


def _spike_slab_mix(
    rng: np.random.Generator,
    post_mu_slab: np.ndarray,
    log_z_slab: float,
    mu_hat: np.ndarray,
    sigma_lik: np.ndarray,
    spike_prob: float,
    lik_logpdf: Any = None,
) -> tuple[np.ndarray, float]:
    """Mixes slab posterior draws with the point mass at ``mu = 0`` (spike-and-slab prior).

    The posterior spike weight is exact given the slab evidence:
    ``Pr(spike | y) = pi L(0) / (pi L(0) + (1 - pi) Z_slab)`` with ``L(0)`` the Gaussian
    cut-space likelihood at ``mu = 0`` (closed form) and ``Z_slab = int p_slab(mu) L(mu) dmu``.

    Args:
      rng: Random generator for the draw-wise spike indicators.
      post_mu_slab: ``(S, d)`` draws from the slab posterior.
      log_z_slab: Log evidence of the slab component.
      mu_hat: ``(d,)`` likelihood mean.
      sigma_lik: ``(d, d)`` likelihood covariance.
      spike_prob: Prior spike probability ``pi`` in ``(0, 1)``.
      lik_logpdf: Optional non-Gaussian likelihood (evaluated at ``mu = 0``).

    Returns:
      ``(post_mu, spike_weight)`` where rows of ``post_mu`` drawn from the spike are zero.
    """
    d_dim = int(mu_hat.shape[0])
    zero = np.zeros((1, d_dim), dtype=np.float64)
    if lik_logpdf is None:
        log_l0 = float(
            sp_stats.multivariate_normal(mean=mu_hat, cov=sigma_lik, allow_singular=True).logpdf(zero[0])
        )
    else:
        log_l0 = float(np.asarray(lik_logpdf(zero), dtype=np.float64).reshape(-1)[0])
    log_a = math.log(float(spike_prob)) + log_l0
    log_b = math.log1p(-float(spike_prob)) + float(log_z_slab)
    if not math.isfinite(log_b):
        weight = 1.0
    elif not math.isfinite(log_a):
        weight = 0.0
    else:
        weight = float(math.exp(log_a - np.logaddexp(log_a, log_b)))
    is_spike = rng.random(post_mu_slab.shape[0]) < weight
    post_mu = np.array(post_mu_slab, dtype=np.float64, copy=True)
    post_mu[is_spike] = 0.0
    return post_mu, weight


def _log_marginal_1d(
    sigma: float,
    t_hat: float,
    s_lik: float,
    family: str,
    *,
    low: float = 0.0,
    high: float = math.inf,
    mean: float = 0.0,
    spike_prob: float | None = None,
) -> float:
    """Log marginal likelihood ``log int p(tau; sigma) N(t_hat; tau, s_lik^2) dtau`` of one coordinate.

    Closed form for the half-normal (skew-normal convolution) and the untruncated normal;
    a fine-grid quadrature otherwise (half-Student-t / gamma bundle families, truncated or
    shifted normals). A spike at zero with probability ``spike_prob`` is mixed in exactly.

    Args:
      sigma: Prior scale.
      t_hat: Likelihood mean of the coordinate.
      s_lik: Likelihood standard deviation of the coordinate.
      family: ``"bundle_halfnormal"``, ``"bundle_halfstudentt"``, ``"bundle_gamma_spike"``,
        ``"halfnormal"`` or ``"normal"`` (truncated to ``[low, high]`` with ``mean``).
      low: Lower truncation of the normal families.
      high: Upper truncation of the normal families.
      mean: Prior mean of the normal families.
      spike_prob: Optional prior probability of a point mass at zero.

    Returns:
      The log marginal likelihood.
    """
    sig = max(float(sigma), 1e-300)
    s_l = max(float(s_lik), 1e-300)
    t_v = float(t_hat)
    halfnormal_like = family in ("bundle_halfnormal", "halfnormal") or (
        family == "normal" and low == 0.0 and math.isinf(high) and mean == 0.0
    )
    if halfnormal_like and low == 0.0 and math.isinf(high) and mean == 0.0:
        tot_sd = math.sqrt(sig**2 + s_l**2)
        log_m = (
            math.log(2.0)
            + float(sp_stats.norm.logpdf(t_v, loc=0.0, scale=tot_sd))
            + float(sp_special.log_ndtr(t_v * sig / (s_l * tot_sd)))
        )
    elif family == "normal" and math.isinf(low) and math.isinf(high):
        log_m = float(sp_stats.norm.logpdf(t_v, loc=float(mean), scale=math.sqrt(sig**2 + s_l**2)))
    else:
        if family.startswith("bundle_"):
            g_lo, g_hi = 0.0, max(t_v + 12.0 * s_l, 40.0 * sig, 1e-12)
        else:
            g_lo = max(low, min(t_v - 12.0 * s_l, mean - 12.0 * sig))
            g_hi = min(high, max(t_v + 12.0 * s_l, mean + 12.0 * sig))
        if not g_hi > g_lo:
            return -math.inf
        grid = np.linspace(g_lo, g_hi, 4001)[1:] if family.startswith("bundle_") else np.linspace(g_lo, g_hi, 4001)
        if family.startswith("bundle_"):
            log_p = _fast_bundle_prior_logpdf(family, sig, grid)
        else:
            a_std = (low - mean) / sig if math.isfinite(low) else -math.inf
            b_std = (high - mean) / sig if math.isfinite(high) else math.inf
            mass = max(float(sp_special.ndtr(b_std) - sp_special.ndtr(a_std)), 1e-300)
            log_p = sp_stats.norm.logpdf(grid, loc=mean, scale=sig) - math.log(mass)
        log_f = log_p + sp_stats.norm.logpdf(t_v, loc=grid, scale=s_l)
        log_m = float(sp_special.logsumexp(log_f) + math.log(float(grid[1] - grid[0])))
    if spike_prob is not None:
        log_spike = math.log(float(spike_prob)) + float(sp_stats.norm.logpdf(t_v, loc=0.0, scale=s_l))
        log_m = float(np.logaddexp(log_spike, math.log1p(-float(spike_prob)) + log_m))
    return log_m


def _empirical_bayes_scale(
    t_hat: float,
    s_lik: float,
    family: str,
    floor: float,
    **marginal_kwargs: Any,
) -> float:
    """Type-II maximum likelihood prior scale ``argmax_{sigma >= floor} log m(t_hat | sigma)``.

    Args:
      t_hat: Likelihood mean of the scale-defining coordinate (the top contrast for bundle
        priors).
      s_lik: Its likelihood standard deviation.
      family: Prior family passed to :func:`_log_marginal_1d`.
      floor: Lower bound on the scale (the in-sample noise scale by default).
      **marginal_kwargs: Extra arguments of :func:`_log_marginal_1d`.

    Returns:
      The empirical-Bayes scale (at least ``floor``).
    """
    flo = max(float(floor), 1e-12)
    upper = max(50.0 * flo, 20.0 * abs(float(t_hat)) + 20.0 * float(s_lik))
    grid = np.geomspace(flo, upper, 161)
    vals = np.array([_log_marginal_1d(g, t_hat, s_lik, family, **marginal_kwargs) for g in grid])
    k_best = int(np.argmax(vals))
    if 0 < k_best < grid.size - 1:
        res = sp_optimize.minimize_scalar(
            lambda lg: -_log_marginal_1d(math.exp(lg), t_hat, s_lik, family, **marginal_kwargs),
            bounds=(math.log(grid[k_best - 1]), math.log(grid[k_best + 1])),
            method="bounded",
            options={"xatol": 1e-6},
        )
        best = math.exp(float(res.x)) if res.success and -res.fun >= vals[k_best] else float(grid[k_best])
    else:
        best = float(grid[k_best])
    return max(best, flo)


class RankKDesignWarning(UserWarning):
    """The fitted contrast design leaves part of the reconciliation unidentified.

    Emitted by :func:`reconcile_multi_arm_draws` when the contrast graph is bipartite
    (has no odd cycle): the arm sampling variances are then identified only up to the
    alternating 2-coloring null direction, so the off-diagonal likelihood covariance
    ``v_0`` rests on an assumption (``bipartite_policy``) instead of on the data. The
    message names the one extra pairwise fit that would close an odd cycle.
    """


class _ArmVarianceArray(np.ndarray):
    """1-D float64 array of arm variances supporting both integer and arm-name indexing."""

    _arm_to_idx: dict[str, int]
    _arms: tuple[str, ...]

    def __new__(
        cls, values: Sequence[float] | np.ndarray, arms: Sequence[str]
    ) -> _ArmVarianceArray:
        obj = np.asarray(values, dtype=np.float64).view(cls)
        obj._arms = tuple(str(a) for a in arms)
        obj._arm_to_idx = {arm: idx for idx, arm in enumerate(obj._arms)}
        return obj

    def __array_finalize__(self, obj: Any) -> None:
        if obj is None:
            return
        self._arms = getattr(obj, "_arms", ())
        self._arm_to_idx = getattr(obj, "_arm_to_idx", {})

    def __getitem__(self, key: Any) -> Any:
        if isinstance(key, str):
            if key not in self._arm_to_idx:
                raise KeyError(f"Unknown arm {key!r}; available arms are {self._arms}")
            return float(super().__getitem__(self._arm_to_idx[key]))
        return super().__getitem__(key)

    def __contains__(self, item: Any) -> bool:
        if isinstance(item, str):
            return item in self._arm_to_idx
        return super().__contains__(item)

    def keys(self) -> tuple[str, ...]:
        """Return the tuple of arm names."""
        return self._arms

    def values(self) -> tuple[float, ...]:
        """Return the tuple of arm variances in arm order."""
        return tuple(float(super().__getitem__(i)) for i in range(len(self._arms)))

    def items(self) -> tuple[tuple[str, float], ...]:
        """Return ``(arm_name, variance)`` pairs in arm order."""
        return tuple(
            (arm, float(super().__getitem__(i))) for i, arm in enumerate(self._arms)
        )

    def get(self, arm: str, default: Any = None) -> Any:
        """Return the variance for ``arm`` or ``default`` if not present."""
        if arm in self._arm_to_idx:
            return float(super().__getitem__(self._arm_to_idx[arm]))
        return default

    def to_dict(self) -> dict[str, float]:
        """Convert to a standard ``{arm: variance}`` dictionary."""
        return dict(self.items())


@dataclass(frozen=True)
class SamplerDiagnostics:
    """Convergence and efficiency diagnostics for the Mode-1 reconciliation sampler.

    Supports both attribute access (``diag.method``, ``diag.ess``, ``diag.r_hat``) and
    dictionary-style lookup (``diag["method"]``, ``diag["ess"]``, ``diag["r_hat"]``).

    Attributes:
      method: Sampler used (``"sir"`` or ``"slice"``).
      ess: Scalar effective sample size (importance-weight ESS for SIR, minimum
        coordinate bulk ESS for slice sampling).
      r_hat: Maximum split R-hat across the ``A - 1`` coordinates.
      ess_per_param: Per-coordinate bulk ESS tuple of length ``A - 1``.
      r_hat_per_param: Per-coordinate split R-hat tuple of length ``A - 1``.
      n_proposals: Number of importance proposals or total MCMC evaluations.
      fallback_triggered: ``True`` if ``sampler="auto"`` fell back from SIR to slice.
      n_iterations: Slice sampler sweeps per chain actually run (burn-in included; ``0`` for
        SIR). Larger than the default when ``auto_continue`` appended blocks.
      converged: Whether the slice sampler met the ``auto_continue`` targets (split R-hat and
        bulk ESS); ``None`` when ``auto_continue`` was off or SIR was used.
    """

    method: str
    ess: float
    r_hat: float
    ess_per_param: tuple[float, ...] = ()
    r_hat_per_param: tuple[float, ...] = ()
    n_proposals: int = 0
    fallback_triggered: bool = False
    n_iterations: int = 0
    converged: bool | None = None

    @property
    def rhat(self) -> float:
        """Alias for ``r_hat``."""
        return self.r_hat

    def _as_dict(self) -> dict[str, Any]:
        return {
            "method": self.method,
            "ess": self.ess,
            "r_hat": self.r_hat,
            "rhat": self.r_hat,
            "ess_per_param": self.ess_per_param,
            "r_hat_per_param": self.r_hat_per_param,
            "n_proposals": self.n_proposals,
            "fallback_triggered": self.fallback_triggered,
            "n_iterations": self.n_iterations,
            "converged": self.converged,
        }

    def __getitem__(self, key: str) -> Any:
        data = self._as_dict()
        if key not in data:
            raise KeyError(key)
        return data[key]

    def __contains__(self, key: object) -> bool:
        return key in self._as_dict()

    def get(self, key: str, default: Any = None) -> Any:
        """Return diagnostic ``key`` or ``default``."""
        return self._as_dict().get(key, default)

    def keys(self) -> tuple[str, ...]:
        """Return canonical diagnostic keys."""
        return (
            "method",
            "ess",
            "r_hat",
            "ess_per_param",
            "r_hat_per_param",
            "n_proposals",
            "fallback_triggered",
            "n_iterations",
            "converged",
        )

    def items(self) -> tuple[tuple[str, Any], ...]:
        """Return ``(key, value)`` pairs for canonical diagnostic keys."""
        data = self._as_dict()
        return tuple((k, data[k]) for k in self.keys())


@dataclass(frozen=True)
class ContrastReconciliationResult:
    """Reconciled posterior draws and summaries for a single contrast ``(treated, ref)``.

    Attributes:
      treated_arm: Treated arm name.
      reference_arm: Reference arm name.
      mode1_draws: ``(S,)`` posterior Mode-1 ATT draws for ``treated - ref``.
      trajectory_draws: ``(S, T_post)`` posterior post-adoption trajectory draws.
      cumulative_draws: ``(S,)`` cumulative post-adoption effect draws over ``horizon``.
      prior_mode1_draws: Optional ``(S,)`` prior Mode-1 draws.
      prior_trajectory_draws: Optional ``(S, T_post)`` prior trajectory draws.
    """

    treated_arm: str
    reference_arm: str
    mode1_draws: np.ndarray
    trajectory_draws: np.ndarray
    cumulative_draws: np.ndarray
    prior_mode1_draws: np.ndarray | None = None
    prior_trajectory_draws: np.ndarray | None = None

    @property
    def att_draws(self) -> np.ndarray:
        """Alias for ``mode1_draws``."""
        return self.mode1_draws

    @property
    def post_att(self) -> np.ndarray:
        """Alias for ``mode1_draws``."""
        return self.mode1_draws

    @property
    def post_att_draws(self) -> np.ndarray:
        """Alias for ``mode1_draws``."""
        return self.mode1_draws

    @property
    def catt_draws(self) -> np.ndarray:
        """Alias for ``trajectory_draws``."""
        return self.trajectory_draws

    @property
    def trt_catt(self) -> np.ndarray:
        """Alias for ``trajectory_draws``."""
        return self.trajectory_draws

    @property
    def trt_catt_draws(self) -> np.ndarray:
        """Alias for ``trajectory_draws``."""
        return self.trajectory_draws

    @property
    def cum_att(self) -> np.ndarray:
        """Alias for ``cumulative_draws``."""
        return self.cumulative_draws

    @property
    def prior_att(self) -> np.ndarray | None:
        """Alias for ``prior_mode1_draws``."""
        return self.prior_mode1_draws

    @property
    def prior_catt(self) -> np.ndarray | None:
        """Alias for ``prior_trajectory_draws``."""
        return self.prior_trajectory_draws

    @property
    def mean(self) -> float:
        """Posterior mean of ``mode1_draws``."""
        return float(np.mean(self.mode1_draws))

    @property
    def sd(self) -> float:
        """Posterior standard deviation of ``mode1_draws``."""
        return float(np.std(self.mode1_draws, ddof=1))

    @property
    def ci_low(self) -> float:
        """2.5% posterior percentile of ``mode1_draws``."""
        return float(np.percentile(self.mode1_draws, 2.5))

    @property
    def ci_high(self) -> float:
        """97.5% posterior percentile of ``mode1_draws``."""
        return float(np.percentile(self.mode1_draws, 97.5))

    def __getitem__(self, key: str) -> Any:
        mapping = {
            "treated_arm": self.treated_arm,
            "reference_arm": self.reference_arm,
            "mode1_draws": self.mode1_draws,
            "att_draws": self.mode1_draws,
            "post_att_draws": self.mode1_draws,
            "trajectory_draws": self.trajectory_draws,
            "catt_draws": self.trajectory_draws,
            "trt_catt_draws": self.trajectory_draws,
            "cumulative_draws": self.cumulative_draws,
            "prior_mode1_draws": self.prior_mode1_draws,
            "prior_trajectory_draws": self.prior_trajectory_draws,
            "mean": self.mean,
            "sd": self.sd,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
        }
        if key not in mapping:
            raise KeyError(key)
        return mapping[key]


@dataclass
class MultiArmReconciliationResult:
    """Output of :func:`reconcile_multi_arm_draws` for ``A >= 2`` arms and ``P >= A - 1`` contrasts.

    Attributes:
      arms: Tuple of ``A`` arm labels, with ``arms[0]`` as the control arm.
      topology: ``"chain"`` or ``"poset"``.
      tree_edges: Tuple of ``A - 1`` spanning-tree edges ``(treated_arm, reference_arm)``.
      fitted_pairs: Tuple of ``P`` fitted contrast pairs ``(treated_arm, reference_arm)``.
      redundant_pairs: Tuple of ``P - (A - 1)`` non-tree fitted contrast pairs.
      tau_mu_draws: ``(S, A - 1)`` posterior Mode-1 draws of ``mu_j = tau_{(arms[j], arms[0])}``.
      tau_increment_draws: ``(S, A - 1)`` posterior Mode-1 draws along ``tree_edges``.
      contrast_mode1_draws: Mapping ``(treated, ref) -> (S,)`` Mode-1 draws across all pairs.
      contrast_trajectory_draws: Mapping ``(treated, ref) -> (S, T_post)`` trajectory draws.
      contrast_catt_draws: Mapping ``(treated, ref) -> (S, T_post)`` trajectory draws (alias).
      contrast_cumulative_draws: Mapping ``(treated, ref) -> (S,)`` cumulative trajectory sums.
      mu_hat: ``(A - 1,)`` Euclidean cut-space projection of the Mode-1 posterior means.
      sigma_lik: ``(A - 1, A - 1)`` Mode-1 Gaussian likelihood covariance on ``mu``.
      arm_variances: Length-``A`` array/mapping of recovered arm variances ``(v_0, ..., v_{A-1})``.
      nnls_residual: Euclidean norm ``||A_var v_hat - s^2||_2`` of the variance fit.
      is_bipartite: ``True`` iff the contrast graph is bipartite (contains no odd cycle).
      cycle_discrepancy_ratio: ``||y_hat - B mu_hat||_2 / max(||y_hat||_2, 1e-12)``.
      trajectory_additivity_max_gap: Maximum absolute difference between a redundant contrast's
        own projected trajectory and its spanning-tree-implied trajectory.
      sampler_diagnostics: :class:`SamplerDiagnostics` with ``method``, ``ess``, ``r_hat``.
      trajectory_additivity_gaps: Per-redundant-pair max trajectory discrepancy.
      prior_mode: Prior family name.
      prior_scale: Resolved prior scale(s).
      prior_mu_draws: Optional ``(S, A - 1)`` prior draws of ``mu``.
      prior_increment_draws: Optional ``(S, A - 1)`` prior draws of tree increments.
      prior_contrast_mode1_draws: Optional mapping ``(treated, ref) -> (S,)`` of prior Mode-1 draws.
      prior_contrast_trajectory_draws: Optional mapping ``(treated, ref) -> (S, T_post)`` of prior
        trajectory draws.
      legacy_3arm_result: Underlying :class:`Nested3ArmReconciliationResult` when ``A = 3``
        chain delegates to the legacy 3-arm path.
      bipartite_policy: Policy that resolved the unidentified arm-variance direction on a
        bipartite contrast graph (``None`` when the graph has an odd cycle).
      null_shift_interval: Feasible interval ``[xi_lo, xi_hi]`` of the shift
        ``v(xi) = v_alloc + xi z`` along the alternating null vector ``z`` (``z_a = +1`` on the
        control's color class, ``-1`` on the other) that keeps every arm variance above its
        floor (``None`` when the graph has an odd cycle).
      null_shift: Shift ``xi`` used by the policy (``0`` for ``"allocation"``; ``None`` for
        ``"marginalize"``, which integrates ``xi`` uniformly over ``null_shift_interval``).
      spike_posterior_prob: Posterior probability of the spike ``mu = 0`` under a
        spike-and-slab bundle prior (``None`` for other priors).
      projection: Cut-space projection used for ``mu_hat`` (``"ols"`` or ``"huber"``).
      projection_weights: ``{pair: weight}`` final IRLS Huber weights (``None`` for ``"ols"``).
      flagged_pairs: Fitted pairs down-weighted by the robust projection (weight ``< 1``),
        most down-weighted first.
    """

    arms: tuple[str, ...]
    topology: str
    tree_edges: tuple[tuple[str, str], ...]
    fitted_pairs: tuple[tuple[str, str], ...]
    redundant_pairs: tuple[tuple[str, str], ...]
    tau_mu_draws: np.ndarray
    tau_increment_draws: np.ndarray
    contrast_mode1_draws: dict[tuple[str, str], np.ndarray]
    contrast_trajectory_draws: dict[tuple[str, str], np.ndarray]
    contrast_catt_draws: dict[tuple[str, str], np.ndarray]
    contrast_cumulative_draws: dict[tuple[str, str], np.ndarray]
    mu_hat: np.ndarray
    sigma_lik: np.ndarray
    arm_variances: _ArmVarianceArray
    nnls_residual: float
    is_bipartite: bool
    cycle_discrepancy_ratio: float
    trajectory_additivity_max_gap: float
    sampler_diagnostics: SamplerDiagnostics
    trajectory_additivity_gaps: dict[tuple[str, str], float]
    prior_mode: str = "bundle_halfnormal"
    prior_scale: float | tuple[float, ...] | None = None
    horizon: int = 0
    own_vs_implied_trajectory_max_gap: float = 0.0
    prior_mu_draws: np.ndarray | None = None
    prior_increment_draws: np.ndarray | None = None
    prior_contrast_mode1_draws: dict[tuple[str, str], np.ndarray] | None = None
    prior_contrast_trajectory_draws: dict[tuple[str, str], np.ndarray] | None = None
    legacy_3arm_result: Nested3ArmReconciliationResult | None = None
    bipartite_policy: str | None = None
    null_shift_interval: tuple[float, float] | None = None
    null_shift: float | None = None
    spike_posterior_prob: float | None = None
    projection: str = "ols"
    projection_weights: dict[tuple[str, str], float] | None = None
    flagged_pairs: tuple[tuple[str, str], ...] = ()

    @property
    def all_pairs(self) -> tuple[tuple[str, str], ...]:
        """Tuple of all canonical ordered contrast pairs stored in ``contrast_mode1_draws``."""
        return tuple(self.contrast_mode1_draws.keys())

    @property
    def own_vs_implied_max_gap(self) -> dict[tuple[str, str], float]:
        """Alias for ``trajectory_additivity_gaps`` (own vs tree-implied trajectory gap per redundant pair)."""
        return self.trajectory_additivity_gaps

    @property
    def prior_tau_mu_draws(self) -> np.ndarray | None:
        """Alias for ``prior_mu_draws``."""
        return self.prior_mu_draws

    @property
    def prior_tau_increment_draws(self) -> np.ndarray | None:
        """Alias for ``prior_increment_draws``."""
        return self.prior_increment_draws

    @property
    def prior_att_bc(self) -> np.ndarray | None:
        """For 3-arm chains, prior Mode-1 draws of ``(arms[1], arms[0])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.prior_att_bc
        if self.prior_contrast_mode1_draws is None:
            return None
        return self.prior_contrast_mode1_draws.get((self.arms[1], self.arms[0]))

    @property
    def prior_att_ab(self) -> np.ndarray | None:
        """For 3-arm chains, prior Mode-1 draws of ``(arms[2], arms[1])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.prior_att_ab
        if self.prior_contrast_mode1_draws is None:
            return None
        return self.prior_contrast_mode1_draws.get((self.arms[2], self.arms[1]))

    @property
    def prior_att_ac(self) -> np.ndarray | None:
        """For 3-arm chains, prior Mode-1 draws of ``(arms[2], arms[0])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.prior_att_ac
        if self.prior_contrast_mode1_draws is None:
            return None
        return self.prior_contrast_mode1_draws.get((self.arms[2], self.arms[0]))

    @property
    def prior_catt_bc(self) -> np.ndarray | None:
        """For 3-arm chains, prior trajectory draws of ``(arms[1], arms[0])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.prior_catt_bc
        if self.prior_contrast_trajectory_draws is None:
            return None
        return self.prior_contrast_trajectory_draws.get((self.arms[1], self.arms[0]))

    @property
    def prior_catt_ab(self) -> np.ndarray | None:
        """For 3-arm chains, prior trajectory draws of ``(arms[2], arms[1])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.prior_catt_ab
        if self.prior_contrast_trajectory_draws is None:
            return None
        return self.prior_contrast_trajectory_draws.get((self.arms[2], self.arms[1]))

    @property
    def prior_catt_ac(self) -> np.ndarray | None:
        """For 3-arm chains, prior trajectory draws of ``(arms[2], arms[0])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.prior_catt_ac
        if self.prior_contrast_trajectory_draws is None:
            return None
        return self.prior_contrast_trajectory_draws.get((self.arms[2], self.arms[0]))

    @property
    def sir_ess(self) -> float:
        """Effective sample size from ``sampler_diagnostics.ess``."""
        return self.sampler_diagnostics.ess

    @property
    def ess(self) -> float:
        """Effective sample size from ``sampler_diagnostics.ess``."""
        return self.sampler_diagnostics.ess

    @property
    def r_hat(self) -> float:
        """Split R-hat from ``sampler_diagnostics.r_hat``."""
        return self.sampler_diagnostics.r_hat

    @property
    def lik_mean(self) -> np.ndarray:
        """Likelihood mean ``mu_hat`` (or legacy 3-arm increment mean if delegated)."""
        if self.legacy_3arm_result is not None and self.legacy_3arm_result.lik_mean is not None:
            return self.legacy_3arm_result.lik_mean
        return self.mu_hat

    @property
    def lik_cov(self) -> np.ndarray:
        """Likelihood covariance ``sigma_lik`` (or legacy 3-arm increment covariance if delegated)."""
        if self.legacy_3arm_result is not None and self.legacy_3arm_result.lik_cov is not None:
            return self.legacy_3arm_result.lik_cov
        return self.sigma_lik

    @property
    def trt_catt_bc(self) -> np.ndarray:
        """For 3-arm chains, posterior trajectory draws of ``(arms[1], arms[0])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.trt_catt_bc
        return self.contrast_trajectory_draws[(self.arms[1], self.arms[0])]

    @property
    def trt_catt_ab(self) -> np.ndarray:
        """For 3-arm chains, posterior trajectory draws of ``(arms[2], arms[1])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.trt_catt_ab
        return self.contrast_trajectory_draws[(self.arms[2], self.arms[1])]

    @property
    def trt_catt_ac(self) -> np.ndarray:
        """For 3-arm chains, posterior trajectory draws of ``(arms[2], arms[0])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.trt_catt_ac
        return self.contrast_trajectory_draws[(self.arms[2], self.arms[0])]

    @property
    def post_att_bc(self) -> np.ndarray:
        """For 3-arm chains, posterior Mode-1 draws of ``(arms[1], arms[0])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.post_att_bc
        return self.contrast_mode1_draws[(self.arms[1], self.arms[0])]

    @property
    def post_att_ab(self) -> np.ndarray:
        """For 3-arm chains, posterior Mode-1 draws of ``(arms[2], arms[1])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.post_att_ab
        return self.contrast_mode1_draws[(self.arms[2], self.arms[1])]

    @property
    def post_att_ac(self) -> np.ndarray:
        """For 3-arm chains, posterior Mode-1 draws of ``(arms[2], arms[0])``."""
        if self.legacy_3arm_result is not None:
            return self.legacy_3arm_result.post_att_ac
        return self.contrast_mode1_draws[(self.arms[2], self.arms[0])]

    def get_contrast_result(
        self, treated: str, ref: str
    ) -> ContrastReconciliationResult:
        """Return the :class:`ContrastReconciliationResult` for ``(treated, ref)``.

        Supports any ordered pair of distinct arms in ``self.arms``. When ``(ref, treated)``
        is stored as the canonical orientation, the returned draws are negated so that
        ``get_contrast_result(ref, treated) == -get_contrast_result(treated, ref)`` holds
        draw by draw.

        Args:
          treated: Treated arm label in ``self.arms``.
          ref: Reference arm label in ``self.arms`` (must differ from ``treated``).

        Returns:
          :class:`ContrastReconciliationResult`.
        """
        t_arm, r_arm = str(treated), str(ref)
        if t_arm not in self.arms or r_arm not in self.arms:
            raise ValueError(
                f"Both treated ({t_arm!r}) and ref ({r_arm!r}) must be in arms {self.arms}"
            )
        if t_arm == r_arm:
            raise ValueError(f"treated and ref must be distinct arms; got {t_arm!r}")
        pair = (t_arm, r_arm)
        rev_pair = (r_arm, t_arm)
        if pair in self.contrast_mode1_draws:
            mode1 = self.contrast_mode1_draws[pair]
            traj = self.contrast_trajectory_draws[pair]
            cum = self.contrast_cumulative_draws[pair]
            pri_m1 = (
                None
                if self.prior_contrast_mode1_draws is None
                else self.prior_contrast_mode1_draws.get(pair)
            )
            pri_tr = (
                None
                if self.prior_contrast_trajectory_draws is None
                else self.prior_contrast_trajectory_draws.get(pair)
            )
        else:
            mode1 = -self.contrast_mode1_draws[rev_pair]
            traj = -self.contrast_trajectory_draws[rev_pair]
            cum = -self.contrast_cumulative_draws[rev_pair]
            pri_m1 = (
                None
                if self.prior_contrast_mode1_draws is None
                or rev_pair not in self.prior_contrast_mode1_draws
                else -self.prior_contrast_mode1_draws[rev_pair]
            )
            pri_tr = (
                None
                if self.prior_contrast_trajectory_draws is None
                or rev_pair not in self.prior_contrast_trajectory_draws
                else -self.prior_contrast_trajectory_draws[rev_pair]
            )
        return ContrastReconciliationResult(
            treated_arm=t_arm,
            reference_arm=r_arm,
            mode1_draws=mode1,
            trajectory_draws=traj,
            cumulative_draws=cum,
            prior_mode1_draws=pri_m1,
            prior_trajectory_draws=pri_tr,
        )

    def __iter__(self) -> Iterator[np.ndarray]:
        if len(self.arms) == 3 and self.topology == "chain":
            yield self.trt_catt_bc
            yield self.trt_catt_ab
            yield self.trt_catt_ac
        else:
            for edge in self.tree_edges:
                yield self.contrast_trajectory_draws[edge]

    def __getitem__(self, item: int | tuple[str, str]) -> Any:
        if isinstance(item, tuple):
            return self.get_contrast_result(item[0], item[1])
        if len(self.arms) == 3 and self.topology == "chain":
            return (self.trt_catt_bc, self.trt_catt_ab, self.trt_catt_ac)[item]
        return self.contrast_trajectory_draws[self.tree_edges[item]]

    def __len__(self) -> int:
        if len(self.arms) == 3 and self.topology == "chain":
            return 3
        return len(self.tree_edges)


def _normalize_results_mapping(
    results: Mapping[Any, RankKLongBetResult],
    arms: tuple[str, ...],
) -> dict[tuple[str, str], RankKLongBetResult]:
    """Validate and normalize the ``results`` mapping keys to ``(treated_arm, ref_arm)``."""
    if not isinstance(results, Mapping) or not results:
        raise ValueError("results must be a non-empty mapping of {(treated, ref): RankKLongBetResult}")
    arm_set = set(arms)
    normalized: dict[tuple[str, str], RankKLongBetResult] = {}
    seen_undirected: set[frozenset[str]] = set()
    for raw_key, fit in results.items():
        if isinstance(raw_key, tuple) and len(raw_key) == 2:
            t_arm, r_arm = str(raw_key[0]), str(raw_key[1])
        elif isinstance(raw_key, str) and len(raw_key) == 2 and raw_key[0] in arm_set and raw_key[1] in arm_set:
            t_arm, r_arm = raw_key[0], raw_key[1]
        else:
            raise ValueError(
                f"results key {raw_key!r} must be a 2-tuple (treated_arm, reference_arm) of arms in {arms}"
            )
        if t_arm not in arm_set or r_arm not in arm_set:
            raise ValueError(
                f"Contrast pair {(t_arm, r_arm)} references an arm not in arms={arms}"
            )
        if t_arm == r_arm:
            raise ValueError(f"Contrast pair {(t_arm, r_arm)} compares an arm to itself")
        undir = frozenset((t_arm, r_arm))
        if undir in seen_undirected:
            raise ValueError(
                f"Duplicate or reverse-oriented contrast for arms {sorted(undir)} in results"
            )
        seen_undirected.add(undir)
        normalized[(t_arm, r_arm)] = fit
    return normalized


def _chain_order_as_poset(
    arms: tuple[str, ...],
    topology: str,
    edges: Sequence[tuple[str, str]] | None,
    fitted_pairs: Sequence[tuple[str, str]],
    split_alpha: Sequence[float] | float | None = None,
) -> tuple[str, Sequence[tuple[str, str]] | None]:
    """Maps a chain order whose adjacent edges are not all fitted to an equivalent poset.

    A chain ``arms[0] < ... < arms[A-1]`` fitted only through, e.g., a vs-control star has no
    fitted chain spanning tree. It is reconciled as the poset whose cover relations are the
    chain edges plus every fitted pair (oriented high -> low), so the spanning tree comes from
    the fitted pairs while the chain order stays the monotonicity constraint. With the
    default uniform split the poset's order-polytope bundle prior equals the chain's
    ``Dirichlet(1, ..., 1)`` bundle prior.

    Args:
      arms: Arm labels in chain order.
      topology: Requested topology.
      edges: Requested poset edges (must be ``None`` for chains).
      fitted_pairs: Fitted contrasts ``(treated, ref)``.
      split_alpha: Requested Dirichlet split (only uniform is supported on the mapped poset).

    Returns:
      ``(topology, edges)`` to use internally (unchanged unless the mapping applies).
    """
    if topology != "chain" or edges is not None:
        return topology, edges
    arm_to_idx = {arm: idx for idx, arm in enumerate(arms)}
    fitted_undir = {frozenset(p) for p in fitted_pairs}
    chain_edges = [(arms[a], arms[a - 1]) for a in range(1, len(arms))]
    if all(frozenset(e) in fitted_undir for e in chain_edges):
        return topology, edges
    if split_alpha is not None and not np.allclose(np.asarray(split_alpha, dtype=np.float64), 1.0):
        raise ValueError(
            "A chain whose adjacent edges are not all fitted is reconciled as a poset with the "
            "chain as cover relations; only the uniform split (split_alpha=None) is supported"
        )
    extra = []
    for t_arm, r_arm in fitted_pairs:
        hi, lo = (t_arm, r_arm) if arm_to_idx[t_arm] > arm_to_idx[r_arm] else (r_arm, t_arm)
        if frozenset((hi, lo)) not in {frozenset(e) for e in chain_edges}:
            extra.append((hi, lo))
    return "poset", [*chain_edges, *extra]


def _check_connected_and_bipartite(
    n_arms: int,
    pair_indices: Sequence[tuple[int, int]],
) -> tuple[bool, bool, np.ndarray]:
    """Check graph connectivity and 2-colorability (bipartiteness) via BFS from arm 0.

    Returns:
      ``(is_connected, is_bipartite, colors)`` where ``colors[a] in {0, 1}`` is the
      2-coloring when bipartite.
    """
    adj: list[list[int]] = [[] for _ in range(n_arms)]
    for j_idx, l_idx in pair_indices:
        adj[j_idx].append(l_idx)
        adj[l_idx].append(j_idx)

    colors = np.full(n_arms, -1, dtype=np.int32)
    colors[0] = 0
    queue: deque[int] = deque([0])
    visited = 1
    is_bipartite = True

    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if colors[v] == -1:
                colors[v] = 1 - colors[u]
                visited += 1
                queue.append(v)
            elif colors[v] == colors[u]:
                is_bipartite = False

    is_connected = visited == n_arms
    return is_connected, is_bipartite, colors


def _validate_and_orient_poset(
    arms: tuple[str, ...],
    topology: str,
    edges: Sequence[tuple[str, str]] | None,
    results_map: Mapping[tuple[str, str], RankKLongBetResult],
) -> tuple[
    tuple[tuple[str, str], ...],
    tuple[tuple[str, str], ...],
    list[tuple[int, int]],
    np.ndarray,
    tuple[int, ...],
]:
    """Validate topology, build spanning tree edges, poset reachability, and maximal arms.

    Returns:
      ``(tree_edges, redundant_pairs, poset_cover_pairs, reachability, maximal_arms)``
      where:
      - ``tree_edges`` has length ``A - 1`` with each edge ``(treated, ref)`` in ``arms``;
      - ``redundant_pairs`` are the fitted keys in ``results_map`` not in ``tree_edges``;
      - ``poset_cover_pairs`` is a list of ``(lower_idx, higher_idx)`` cover relations;
      - ``reachability`` is an ``(A, A)`` boolean matrix with ``reachability[l, j] == True``
        iff ``l <= j`` in the poset;
      - ``maximal_arms`` are the arm indices ``j >= 1`` that have no strict successor.
    """
    n_arms = len(arms)
    arm_to_idx = {arm: idx for idx, arm in enumerate(arms)}
    undir_to_fitted: dict[frozenset[str], tuple[str, str]] = {
        frozenset(pair): pair for pair in results_map.keys()
    }

    if topology == "chain":
        if edges is not None:
            raise ValueError("edges must be None when topology='chain'")
        tree_edges_list: list[tuple[str, str]] = []
        cover_pairs: list[tuple[int, int]] = []
        for a_idx in range(1, n_arms):
            t_arm, r_arm = arms[a_idx], arms[a_idx - 1]
            if frozenset((t_arm, r_arm)) not in undir_to_fitted:
                raise ValueError(
                    f"Missing tree-edge fit for chain edge {(t_arm, r_arm)!r}; every spanning-tree "
                    "edge must have a fitted contrast in results"
                )
            tree_edges_list.append((t_arm, r_arm))
            cover_pairs.append((a_idx - 1, a_idx))
    elif topology == "poset":
        if edges is None or len(edges) < n_arms - 1:
            raise ValueError(
                f"topology='poset' requires edges with at least A - 1 = {n_arms - 1} cover relations; "
                f"got {edges!r}"
            )
        raw_pairs: list[tuple[int, int]] = []
        for raw_u, raw_v in edges:
            u_str, v_str = str(raw_u), str(raw_v)
            if u_str not in arm_to_idx or v_str not in arm_to_idx or u_str == v_str:
                raise ValueError(
                    f"Invalid poset edge {(raw_u, raw_v)!r} for arms {arms}"
                )
            raw_pairs.append((arm_to_idx[u_str], arm_to_idx[v_str]))

        # Check connectivity of the supplied poset edges
        poset_conn, _, _ = _check_connected_and_bipartite(n_arms, raw_pairs)
        if not poset_conn:
            raise ValueError(
                f"Poset edges {list(edges)} form a disconnected graph on arms {arms}"
            )

        # Determine whether user passed (lower, higher) or (higher, lower) by checking arm 0
        zero_first = any(u == 0 for u, _ in raw_pairs)
        zero_second = any(v == 0 for _, v in raw_pairs)
        if zero_first and not zero_second:
            cover_pairs = [(u, v) for u, v in raw_pairs]
        elif zero_second and not zero_first:
            cover_pairs = [(v, u) for u, v in raw_pairs]
        else:
            # Orient from smaller arm index (closer to control 0) to larger arm index
            cover_pairs = [(min(u, v), max(u, v)) for u, v in raw_pairs]

        # Check acyclicity via transitive closure
        reach = np.eye(n_arms, dtype=bool)
        for lo, hi in cover_pairs:
            reach[lo, hi] = True
        for k in range(n_arms):
            reach |= np.outer(reach[:, k], reach[k, :])
        for i in range(n_arms):
            for j in range(i + 1, n_arms):
                if reach[i, j] and reach[j, i]:
                    raise ValueError(
                        f"Poset edges {list(edges)} contain a directed cycle involving "
                        f"{arms[i]!r} and {arms[j]!r}"
                    )
        if not np.all(reach[0, :]):
            raise ValueError(
                f"Control arm {arms[0]!r} must precede all other arms in the poset; "
                f"unreachable arms: {[arms[i] for i in range(n_arms) if not reach[0, i]]}"
            )

        # Every cover edge when len(edges) == A - 1 is a mandatory tree edge; when len(edges) > A - 1
        # (e.g. 2x2 diamond), pick a spanning tree of fitted cover edges (requiring every arm j >= 1
        # to have an incoming fitted cover edge, and requiring all edges if len(edges) == n_arms - 1).
        if len(cover_pairs) == n_arms - 1:
            tree_edges_list = []
            for lo, hi in cover_pairs:
                t_arm, r_arm = arms[hi], arms[lo]
                if frozenset((t_arm, r_arm)) not in undir_to_fitted:
                    raise ValueError(
                        f"Missing tree-edge fit for poset edge {(t_arm, r_arm)!r}; every "
                        "spanning-tree edge must have a fitted contrast in results"
                    )
                tree_edges_list.append((t_arm, r_arm))
        else:
            # Build BFS spanning tree from control 0 using cover edges present in results_map
            parent_of: dict[int, int] = {}
            visited_tree = {0}
            bfs_q: deque[int] = deque([0])
            while bfs_q:
                u = bfs_q.popleft()
                for lo, hi in cover_pairs:
                    if lo == u and hi not in visited_tree:
                        if frozenset((arms[hi], arms[lo])) in undir_to_fitted:
                            visited_tree.add(hi)
                            parent_of[hi] = lo
                            bfs_q.append(hi)
            if len(visited_tree) != n_arms:
                missing_arms = [arms[a] for a in range(1, n_arms) if a not in visited_tree]
                raise ValueError(
                    f"Missing tree-edge fit in poset: no fitted cover-edge path from control "
                    f"{arms[0]!r} to arm(s) {missing_arms}"
                )
            tree_edges_list = [(arms[a], arms[parent_of[a]]) for a in range(1, n_arms)]
    else:
        raise ValueError(f"topology must be one of {_TOPOLOGIES}; got {topology!r}")

    reachability = np.eye(n_arms, dtype=bool)
    for lo, hi in cover_pairs:
        reachability[lo, hi] = True
    for k in range(n_arms):
        reachability |= np.outer(reachability[:, k], reachability[k, :])

    maximal_arms = tuple(
        j for j in range(1, n_arms) if not any(reachability[j, k] for k in range(n_arms) if k != j)
    )

    tree_undir = {frozenset(edge) for edge in tree_edges_list}
    redundant_list = [
        pair for pair in results_map.keys() if frozenset(pair) not in tree_undir
    ]
    return (
        tuple(tree_edges_list),
        tuple(redundant_list),
        cover_pairs,
        reachability,
        maximal_arms,
    )


def _enumerate_linear_extensions(
    reachability_or_non_top: np.ndarray | Sequence[int],
    poset_edges_idx: Sequence[tuple[int, int]] | None = None,
) -> np.ndarray:
    """Enumerate all topological sort permutations of an ``m``-element interior poset.

    Accepts either an ``(m, m)`` boolean reachability matrix ``reachability_or_non_top``
    or ``(non_top_arms, poset_edges_idx)`` where ``poset_edges_idx`` are directed edges
    (either ``(child, parent)`` or ``(lower, higher)``).

    Returns:
      ``(n_ext, m)`` int32 array of linear extensions (each row is a permutation of ``0..m-1``).
    """
    if poset_edges_idx is not None:
        non_top = [int(a) for a in reachability_or_non_top]
        m = len(non_top)
        if m == 0:
            return np.zeros((1, 0), dtype=np.int32)
        n_total = max([0, *non_top, *(max(u, v) for u, v in poset_edges_idx)]) + 1
        reach_all = np.eye(n_total, dtype=bool)
        zero_first = any(u == 0 for u, _ in poset_edges_idx)
        zero_second = any(v == 0 for _, v in poset_edges_idx)
        for u, v in poset_edges_idx:
            if zero_second and not zero_first:
                lo, hi = int(v), int(u)
            elif zero_first and not zero_second:
                lo, hi = int(u), int(v)
            else:
                lo, hi = min(int(u), int(v)), max(int(u), int(v))
            reach_all[lo, hi] = True
        for k in range(n_total):
            reach_all |= np.outer(reach_all[:, k], reach_all[k, :])
        reachability_interior = reach_all[np.ix_(non_top, non_top)]
    else:
        reachability_interior = np.asarray(reachability_or_non_top, dtype=bool)
        m = int(reachability_interior.shape[0])
        if m == 0:
            return np.zeros((1, 0), dtype=np.int32)

    strict_pred_mask = [
        int(sum((1 << i) for i in range(m) if i != j and reachability_interior[i, j]))
        for j in range(m)
    ]
    extensions: list[list[int]] = []
    current: list[int] = []

    def _backtrack(placed_mask: int) -> None:
        if len(current) == m:
            extensions.append(current.copy())
            return
        for j in range(m):
            if (
                not (placed_mask & (1 << j))
                and (strict_pred_mask[j] & placed_mask) == strict_pred_mask[j]
            ):
                current.append(j)
                _backtrack(placed_mask | (1 << j))
                current.pop()

    _backtrack(0)
    return np.asarray(extensions, dtype=np.int32)


def _draw_prior_mu(
    rng: np.random.Generator,
    *,
    n: int,
    dim: int,
    topology: str,
    prior_mode: str,
    scale_val: float,
    alpha_vec: np.ndarray,
    monotone: bool,
    tau_clip: tuple[float, float],
    poset_edges_idx: Sequence[tuple[int, int]],
    non_top_arms: Sequence[int],
    linear_extensions: np.ndarray | None = None,
) -> np.ndarray:
    """Draw ``n`` samples of ``mu in R^{dim}`` from the specified multi-arm prior family."""
    if prior_mode in _BUNDLE_PRIOR_MODES:
        tot = _draw_bundle_prior(rng, prior_mode, float(scale_val), n)
        if topology == "chain" and monotone:
            u_mat = rng.dirichlet(np.asarray(alpha_vec, dtype=np.float64), size=n)
            theta_mat = u_mat * tot[:, None]
            return np.cumsum(theta_mat, axis=1)
        m_int = len(non_top_arms)
        pts = np.zeros((n, dim), dtype=np.float64)
        pts[:, -1] = tot
        if m_int > 0:
            interior_cols = [a - 1 for a in non_top_arms]
            if monotone:
                lin_exts = (
                    linear_extensions
                    if linear_extensions is not None
                    else _enumerate_linear_extensions(non_top_arms, poset_edges_idx)
                )
                n_exts = int(lin_exts.shape[0])
                ext_idx = rng.integers(0, n_exts, size=n)
                chosen_exts = lin_exts[ext_idx]
                u_sorted = np.sort(rng.uniform(0.0, 1.0, size=(n, m_int)), axis=1)
                w_mat = np.empty((n, m_int), dtype=np.float64)
                w_mat[np.arange(n)[:, None], chosen_exts] = u_sorted
            else:
                w_mat = rng.uniform(0.0, 1.0, size=(n, m_int))
            pts[:, interior_cols] = tot[:, None] * w_mat
        return pts

    low_c = max(float(tau_clip[0]), 0.0) if prior_mode == "independent_halfnormal" else float(tau_clip[0])
    high_c = float(tau_clip[1])
    sd = float(scale_val)
    a_std = low_c / sd if math.isfinite(low_c) else -np.inf
    b_std = high_c / sd if math.isfinite(high_c) else np.inf
    pts = np.zeros((n, dim), dtype=np.float64)
    for j in range(dim):
        if math.isfinite(a_std) or math.isfinite(b_std):
            pts[:, j] = sp_stats.truncnorm.rvs(a_std, b_std, loc=0.0, scale=sd, size=n, random_state=rng)
        else:
            pts[:, j] = sd * rng.standard_normal(n)
    return pts


def _prior_logpdf_mu(
    pts: np.ndarray,
    *,
    dim: int,
    topology: str,
    prior_mode: str,
    scale_val: float,
    alpha_vec: np.ndarray,
    monotone: bool,
    tau_clip: tuple[float, float],
    poset_edges_idx: Sequence[tuple[int, int]],
    non_top_arms: Sequence[int],
    linear_extensions: np.ndarray | None = None,
) -> np.ndarray:
    """Evaluate log prior density ``log p(mu)`` at rows of ``pts`` of shape ``(N, dim)``."""
    pts_arr = np.atleast_2d(np.asarray(pts, dtype=np.float64))
    n_pts = pts_arr.shape[0]
    low_c = (
        max(float(tau_clip[0]), 0.0)
        if prior_mode == "independent_halfnormal"
        else float(tau_clip[0])
    )
    high_c = float(tau_clip[1])
    valid = np.all((pts_arr >= low_c) & (pts_arr <= high_c), axis=1)

    if prior_mode in _BUNDLE_PRIOR_MODES:
        if topology == "chain" and monotone:
            theta_mat = np.empty_like(pts_arr)
            theta_mat[:, 0] = pts_arr[:, 0]
            if dim > 1:
                theta_mat[:, 1:] = np.diff(pts_arr, axis=1)
            valid &= np.all(theta_mat > 0.0, axis=1)
            out = np.full(n_pts, -np.inf, dtype=np.float64)
            if not np.any(valid):
                return out
            th_v = theta_mat[valid]
            tot_v = np.sum(th_v, axis=1)
            a_vec = np.asarray(alpha_vec, dtype=np.float64)
            u_v = th_v / np.maximum(tot_v[:, None], 1e-300)
            log_dir_norm = float(sp_special.gammaln(np.sum(a_vec)) - np.sum(sp_special.gammaln(a_vec)))
            dir_log = log_dir_norm + np.sum((a_vec[None, :] - 1.0) * np.log(np.maximum(u_v, 1e-300)), axis=1)
            out[valid] = (
                _fast_bundle_prior_logpdf(prior_mode, float(scale_val), tot_v)
                - float(dim - 1) * np.log(np.maximum(tot_v, 1e-300))
                + dir_log
            )
            return out

        m_int = len(non_top_arms)
        interior_cols = [a - 1 for a in non_top_arms]
        tot = pts_arr[:, -1]
        valid &= tot > 0.0
        if m_int > 0:
            w_int = pts_arr[:, interior_cols]
            valid &= np.all((w_int >= 0.0) & (w_int <= tot[:, None]), axis=1)
            if monotone:
                zero_first = any(u == 0 for u, _ in poset_edges_idx)
                zero_second = any(v == 0 for _, v in poset_edges_idx)
                for u, v in poset_edges_idx:
                    if zero_second and not zero_first:
                        lo, hi = int(v), int(u)
                    elif zero_first and not zero_second:
                        lo, hi = int(u), int(v)
                    else:
                        lo, hi = min(int(u), int(v)), max(int(u), int(v))
                    if lo >= 1 and hi >= 1:
                        valid &= pts_arr[:, hi - 1] >= pts_arr[:, lo - 1]
                    elif lo == 0 and hi >= 1:
                        valid &= pts_arr[:, hi - 1] >= 0.0
        out = np.full(n_pts, -np.inf, dtype=np.float64)
        if not np.any(valid):
            return out
        lin_exts = (
            linear_extensions
            if linear_extensions is not None
            else (_enumerate_linear_extensions(non_top_arms, poset_edges_idx) if (monotone and m_int > 0) else np.zeros((1, m_int), dtype=np.int32))
        )
        n_exts = max(int(lin_exts.shape[0]), 1)
        log_inv_vol = float(sp_special.gammaln(m_int + 1) - math.log(n_exts)) if (monotone and m_int > 0) else 0.0
        tot_v = tot[valid]
        out[valid] = (
            _fast_bundle_prior_logpdf(prior_mode, float(scale_val), tot_v)
            - float(m_int) * np.log(np.maximum(tot_v, 1e-300))
            + log_inv_vol
        )
        return out

    if monotone:
        if topology == "chain":
            theta_mat = np.empty_like(pts_arr)
            theta_mat[:, 0] = pts_arr[:, 0]
            if dim > 1:
                theta_mat[:, 1:] = np.diff(pts_arr, axis=1)
            valid &= np.all(theta_mat >= 0.0, axis=1)
        else:
            zero_first = any(u == 0 for u, _ in poset_edges_idx)
            zero_second = any(v == 0 for _, v in poset_edges_idx)
            for u, v in poset_edges_idx:
                if zero_second and not zero_first:
                    lo, hi = int(v), int(u)
                elif zero_first and not zero_second:
                    lo, hi = int(u), int(v)
                else:
                    lo, hi = min(int(u), int(v)), max(int(u), int(v))
                if lo >= 1 and hi >= 1:
                    valid &= pts_arr[:, hi - 1] >= pts_arr[:, lo - 1]
                elif lo == 0 and hi >= 1:
                    valid &= pts_arr[:, hi - 1] >= 0.0

    if prior_mode == "flat":
        return np.where(valid, 0.0, -np.inf)

    sd = float(scale_val)
    out = np.full(n_pts, -np.inf, dtype=np.float64)
    if not np.any(valid):
        return out
    pts_v = pts_arr[valid]
    a_std = low_c / sd if math.isfinite(low_c) else -np.inf
    b_std = high_c / sd if math.isfinite(high_c) else np.inf
    trunc_mass = max(float(sp_special.ndtr(b_std) - sp_special.ndtr(a_std)), 1e-300)
    log_norm_per_dim = -0.5 * math.log(2.0 * math.pi) - math.log(sd) - math.log(trunc_mass)
    out[valid] = -0.5 * np.sum((pts_v / sd) ** 2, axis=1) + float(dim) * log_norm_per_dim
    return out


def _infer_arm_sample_sizes(
    arms: tuple[str, ...],
    results_map: Mapping[tuple[str, str], RankKLongBetResult],
    explicit_sizes: Mapping[str, float] | Sequence[float] | None,
) -> np.ndarray:
    """Resolve per-arm sample sizes ``(n_0, ..., n_{A-1})`` for bipartite variance regularization."""
    n_arms = len(arms)
    if explicit_sizes is not None:
        if isinstance(explicit_sizes, Mapping):
            sizes = np.array([float(explicit_sizes[a]) for a in arms], dtype=np.float64)
        else:
            sizes = np.asarray(explicit_sizes, dtype=np.float64)
        if sizes.shape != (n_arms,) or np.any(~np.isfinite(sizes)) or np.any(sizes <= 0.0):
            raise ValueError(
                f"arm_sample_sizes must provide {n_arms} finite positive numbers; got {explicit_sizes!r}"
            )
        return sizes

    arm_to_idx = {arm: idx for idx, arm in enumerate(arms)}
    counts: list[list[float]] = [[] for _ in range(n_arms)]
    for (t_arm, r_arm), res in results_map.items():
        try:
            catt = np.asarray(res.catt_draws)
            trt_mask = np.asarray(res.treated_units, dtype=bool)
        except (AttributeError, TypeError, ValueError):
            continue
        n_total = int(catt.shape[1]) if catt.ndim == 3 else 0
        n_trt = int(np.sum(trt_mask)) if trt_mask.ndim == 1 else 0
        n_ref = n_total - n_trt
        if n_trt >= 2 and n_ref >= 2:
            counts[arm_to_idx[t_arm]].append(float(n_trt))
            counts[arm_to_idx[r_arm]].append(float(n_ref))

    if all(len(c) > 0 for c in counts):
        return np.array([float(np.median(c)) for c in counts], dtype=np.float64)
    return np.ones(n_arms, dtype=np.float64)


def _recover_arm_variances(
    arms: tuple[str, ...],
    pair_indices: Sequence[tuple[int, int]],
    s_sq: np.ndarray,
    is_bipartite: bool,
    colors: np.ndarray,
    arm_sample_sizes: np.ndarray,
    min_group_var_share: float,
    bipartite_ridge_weight: float,
    bipartite_policy: str = "allocation",
) -> tuple[_ArmVarianceArray, float]:
    """Recover arm sampling variances ``v in R_+^A`` from pairwise posterior variances ``s_p^2``."""
    n_arms = len(arms)
    n_pairs = len(pair_indices)
    share = float(min_group_var_share)
    if not (0.0 <= share < 0.5):
        raise ValueError(f"min_group_var_share must lie in [0, 0.5); got {min_group_var_share!r}")

    a_var = np.zeros((n_pairs, n_arms), dtype=np.float64)
    incident_vars: list[list[float]] = [[] for _ in range(n_arms)]
    for p_idx, (j_idx, l_idx) in enumerate(pair_indices):
        a_var[p_idx, j_idx] = 1.0
        a_var[p_idx, l_idx] = 1.0
        var_p = max(float(s_sq[p_idx]), 1e-12)
        incident_vars[j_idx].append(var_p)
        incident_vars[l_idx].append(var_p)

    v_floor = np.array(
        [share * min(inc) if inc else 1e-12 for inc in incident_vars],
        dtype=np.float64,
    )

    if not is_bipartite:
        # Exact reproduction of _backsolve_group_variances for A = 3 triangle
        if n_arms == 3 and n_pairs == 3:
            pair_var_map: dict[frozenset[int], float] = {
                frozenset((j_idx, l_idx)): float(s_sq[p_idx])
                for p_idx, (j_idx, l_idx) in enumerate(pair_indices)
            }
            v_bc = pair_var_map[frozenset((1, 0))]
            v_ab = pair_var_map[frozenset((2, 1))]
            v_ac = pair_var_map[frozenset((2, 0))]
            v_b = max(0.5 * (v_bc + v_ab - v_ac), share * min(v_bc, v_ab))
            v_c = max(v_bc - v_b, share * v_bc)
            v_a = max(v_ab - v_b, share * v_ab)
            v_vec = np.array([v_c, v_b, v_a], dtype=np.float64)
            residual = float(np.linalg.norm(a_var @ v_vec - s_sq))
            return _ArmVarianceArray(v_vec, arms), residual

        sol = sp_optimize.lsq_linear(
            a_var,
            s_sq,
            bounds=(v_floor, np.full(n_arms, np.inf)),
        )
        v_vec = np.maximum(np.asarray(sol.x, dtype=np.float64), v_floor)
        residual = float(np.linalg.norm(a_var @ v_vec - s_sq))
        return _ArmVarianceArray(v_vec, arms), residual

    # Bipartite graph: rank(A_var) = A - 1 with null vector z_a = (-1)^{colors[a]}.
    suggestion = _suggest_odd_cycle_fit(arms, pair_indices, colors)
    if suggestion is None:
        hint = "No single extra pairwise fit closes an odd cycle on these arms."
    else:
        (fit_j, fit_l), cycle = suggestion
        hint = (
            f"The missing odd cycle is {' -> '.join(arms[a] for a in cycle)}: fitting the "
            f"contrast ({arms[fit_j]!r}, {arms[fit_l]!r}) would close it and identify every arm "
            "variance."
        )
    _warn_rank_k(
        f"reconcile_multi_arm_draws: contrast graph on arms {arms} is bipartite (contains no odd "
        "cycle), so the unsigned variance incidence matrix A_var has rank A - 1 and arm variances "
        "are under-identified along the alternating 2-coloring null vector; the shared-arm "
        f"covariance of the likelihood then rests on bipartite_policy={bipartite_policy!r} "
        f"rather than on the data. {hint}",
        category=RankKDesignWarning,
    )
    inv_n = 1.0 / np.maximum(arm_sample_sizes, 1e-12)
    pair_weights = a_var @ inv_n
    sigma_sq_hat = float(
        np.dot(s_sq, pair_weights) / max(float(np.dot(pair_weights, pair_weights)), 1e-12)
    )
    v_target = np.maximum(sigma_sq_hat * inv_n, v_floor)

    null_vec = np.where(colors == 0, 1.0, -1.0).astype(np.float64)
    null_vec /= max(float(np.linalg.norm(null_vec)), 1e-12)
    scale_ref = max(float(np.mean(s_sq)), 1e-6)
    lam_null = max(float(bipartite_ridge_weight), 1e-4)
    # Regularize strongly along the 1-D null direction and lightly on the full vector for strict SPD
    a_aug = np.vstack([
        a_var,
        math.sqrt(scale_ref) * null_vec[None, :],
        math.sqrt(lam_null * scale_ref) * np.eye(n_arms, dtype=np.float64),
    ])
    b_aug = np.concatenate([
        s_sq,
        np.array([math.sqrt(scale_ref) * float(np.dot(null_vec, v_target))]),
        math.sqrt(lam_null * scale_ref) * v_target,
    ])
    sol = sp_optimize.lsq_linear(
        a_aug,
        b_aug,
        bounds=(v_floor, np.full(n_arms, np.inf)),
    )
    v_vec = np.maximum(np.asarray(sol.x, dtype=np.float64), v_floor)
    residual = float(np.linalg.norm(a_var @ v_vec - s_sq))
    return _ArmVarianceArray(v_vec, arms), residual


def _suggest_odd_cycle_fit(
    arms: tuple[str, ...],
    pair_indices: Sequence[tuple[int, int]],
    colors: np.ndarray,
) -> tuple[tuple[int, int], list[int]] | None:
    """Return the extra fit ``(j, l)`` that closes the shortest odd cycle, and that cycle.

    On a connected bipartite graph any unfitted pair of same-colored arms closes an odd
    cycle (the tree path between them has even length). Pairs adjacent in ``arms`` order
    are preferred, then pairs away from the control arm.

    Args:
      arms: Arm labels.
      pair_indices: Fitted contrasts as ``(treated_idx, ref_idx)``.
      colors: 2-coloring of the arms.

    Returns:
      ``((j, l), cycle)`` with ``j > l`` and ``cycle`` the arm-index path ``l -> ... -> j -> l``,
      or ``None`` when every same-colored pair is already fitted (cannot happen on a
      bipartite graph with ``A >= 3``).
    """
    n_arms = len(arms)
    fitted = {frozenset(p) for p in pair_indices}
    adj: list[list[int]] = [[] for _ in range(n_arms)]
    for j_idx, l_idx in pair_indices:
        adj[j_idx].append(l_idx)
        adj[l_idx].append(j_idx)
    candidates = [
        (j, l)
        for l in range(n_arms)
        for j in range(l + 1, n_arms)
        if colors[j] == colors[l] and frozenset((j, l)) not in fitted
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda p: (p[0] - p[1], p[1] == 0, p[1], p[0]))
    j_best, l_best = candidates[0]
    parent = {l_best: -1}
    queue: deque[int] = deque([l_best])
    while queue:
        u = queue.popleft()
        for v in sorted(adj[u]):
            if v not in parent:
                parent[v] = u
                queue.append(v)
    path = [j_best]
    while parent[path[-1]] != -1:
        path.append(parent[path[-1]])
    path.reverse()  # l_best -> ... -> j_best
    return (j_best, l_best), [*path, l_best]


def _bipartite_null_interval(
    arm_variances: np.ndarray,
    pair_indices: Sequence[tuple[int, int]],
    s_sq: np.ndarray,
    colors: np.ndarray,
    min_group_var_share: float,
) -> tuple[np.ndarray, float, float]:
    """Null vector ``z`` and the feasible shift interval keeping ``v + xi z >= floor``.

    Args:
      arm_variances: ``(A,)`` allocation solution ``v_alloc``.
      pair_indices: Fitted contrasts as ``(treated_idx, ref_idx)``.
      s_sq: ``(P,)`` pairwise posterior variances.
      colors: 2-coloring of the arms (control has color ``0``).
      min_group_var_share: Floor share ``eta`` of the minimum incident pairwise variance.

    Returns:
      ``(z, xi_lo, xi_hi)`` with ``z_a = +1`` on the control's color class and ``-1`` else.
    """
    n_arms = int(arm_variances.shape[0])
    incident: list[list[float]] = [[] for _ in range(n_arms)]
    for p_idx, (j_idx, l_idx) in enumerate(pair_indices):
        var_p = max(float(s_sq[p_idx]), 1e-12)
        incident[j_idx].append(var_p)
        incident[l_idx].append(var_p)
    floor = np.array(
        [float(min_group_var_share) * min(inc) if inc else 1e-12 for inc in incident],
        dtype=np.float64,
    )
    z_vec = np.where(colors == colors[0], 1.0, -1.0).astype(np.float64)
    v_arr = np.asarray(arm_variances, dtype=np.float64)
    lo_terms = [floor[a] - v_arr[a] for a in range(n_arms) if z_vec[a] > 0]
    hi_terms = [v_arr[a] - floor[a] for a in range(n_arms) if z_vec[a] < 0]
    xi_lo = max(lo_terms) if lo_terms else -np.inf
    xi_hi = min(hi_terms) if hi_terms else np.inf
    if xi_lo > xi_hi:
        xi_lo = xi_hi = 0.0
    return z_vec, float(min(xi_lo, 0.0)), float(max(xi_hi, 0.0))


def _sigma_lik_from_arm_variances(v_vec: np.ndarray) -> np.ndarray:
    """Cut-space likelihood covariance ``diag(v_1, ..., v_{A-1}) + v_0 1 1'``."""
    v_arr = np.asarray(v_vec, dtype=np.float64)
    d_dim = v_arr.shape[0] - 1
    return np.diag(v_arr[1:]) + float(v_arr[0]) * np.ones((d_dim, d_dim), dtype=np.float64)


def _huber_cut_space_projection(
    b_cut: np.ndarray,
    y_hat: np.ndarray,
    s_pair: np.ndarray,
    *,
    huber_k: float = _HUBER_K,
    max_iter: int = 200,
    tol: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Robust cut-space projection by IRLS with Huber weights on the cycle residual.

    Starts from the Euclidean projection ``mu = (B'B)^{-1} B' y`` and iterates
    ``mu <- (B' W B)^{-1} B' W y`` with ``w_p = min(1, k / |r_p / scale|)``,
    ``r = y - B mu``, ``k = 1.345`` and ``scale = 1.4826 * MAD(r)`` (floored at ``1e-6`` times
    the median pairwise posterior SD). Under the exact shared-arm model ``y`` lies in the
    column space of ``B``, every residual is zero, all weights are one and the Euclidean
    projection (Proposition 1) is returned unchanged.

    Args:
      b_cut: ``(P, A - 1)`` reduced cut matrix.
      y_hat: ``(P,)`` pairwise Mode-1 posterior means.
      s_pair: ``(P,)`` pairwise posterior SDs (only for the scale floor).
      huber_k: Huber tuning constant.
      max_iter: Maximum IRLS iterations.
      tol: Convergence tolerance on ``mu``.

    Returns:
      ``(mu_hat, weights, scale)``.
    """
    btb = b_cut.T @ b_cut
    mu = np.linalg.solve(btb, b_cut.T @ y_hat)
    weights = np.ones(y_hat.shape[0], dtype=np.float64)
    floor = 1e-6 * max(float(np.median(s_pair)), 1e-12)
    scale = floor
    for _ in range(int(max_iter)):
        resid = y_hat - b_cut @ mu
        scale = max(1.4826 * float(np.median(np.abs(resid - np.median(resid)))), floor)
        u_abs = np.abs(resid) / scale
        weights = np.where(u_abs <= huber_k, 1.0, huber_k / np.maximum(u_abs, 1e-300))
        bw = b_cut * weights[:, None]
        mu_new = np.linalg.solve(b_cut.T @ bw, bw.T @ y_hat)
        if float(np.max(np.abs(mu_new - mu))) <= tol * max(1.0, float(np.max(np.abs(mu)))):
            mu = mu_new
            break
        mu = mu_new
    return mu, weights, scale


def _build_support_inequalities(
    n_arms: int,
    topology: str,
    prior_mode: str,
    monotone: bool,
    tau_clip: tuple[float, float],
    t_tree: np.ndarray,
    cover_pairs: Sequence[tuple[int, int]],
    top_idx: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Construct linear inequality matrices ``(A_ineq, b_ineq)`` for ``A_ineq @ mu >= b_ineq``."""
    d_dim = n_arms - 1
    low_clip, high_clip = float(tau_clip[0]), float(tau_clip[1])
    rows_a: list[np.ndarray] = []
    rows_b: list[float] = []

    def _add_ineq(a_row: np.ndarray, b_val: float) -> None:
        if math.isfinite(b_val):
            rows_a.append(np.asarray(a_row, dtype=np.float64))
            rows_b.append(float(b_val))

    if topology == "chain" and (monotone or prior_mode == "independent_halfnormal"):
        # Chain increments theta = t_tree @ mu lie in [low_clip, high_clip]
        eff_low = max(low_clip, 0.0) if (monotone or prior_mode in _BUNDLE_PRIOR_MODES) else low_clip
        for e_idx in range(d_dim):
            row = t_tree[e_idx]
            if math.isfinite(eff_low):
                _add_ineq(row, eff_low)
            if math.isfinite(high_clip):
                _add_ineq(-row, -high_clip)
    else:
        eff_low = max(low_clip, 0.0) if prior_mode in _BUNDLE_PRIOR_MODES else low_clip
        for j_idx in range(d_dim):
            e_j = np.zeros(d_dim, dtype=np.float64)
            e_j[j_idx] = 1.0
            if math.isfinite(eff_low):
                _add_ineq(e_j, eff_low)
            if math.isfinite(high_clip):
                _add_ineq(-e_j, -high_clip)

        if prior_mode in _BUNDLE_PRIOR_MODES and top_idx is not None:
            top_col = top_idx - 1
            for j_idx in range(d_dim):
                if j_idx == top_col:
                    continue
                # 0 <= mu_j <= mu_top
                diff = np.zeros(d_dim, dtype=np.float64)
                diff[top_col] = 1.0
                diff[j_idx] = -1.0
                _add_ineq(diff, 0.0)

        if monotone:
            for lo_arm, hi_arm in cover_pairs:
                diff = np.zeros(d_dim, dtype=np.float64)
                if hi_arm >= 1:
                    diff[hi_arm - 1] = 1.0
                if lo_arm >= 1:
                    diff[lo_arm - 1] = -1.0
                _add_ineq(diff, max(low_clip, 0.0) if lo_arm == 0 else 0.0)

    if not rows_a:
        return np.zeros((0, d_dim), dtype=np.float64), np.zeros(0, dtype=np.float64)
    return np.vstack(rows_a), np.asarray(rows_b, dtype=np.float64)


def _inside_polyhedron(
    pts: np.ndarray, a_ineq: np.ndarray, b_ineq: np.ndarray, tol: float = 0.0
) -> np.ndarray:
    """Return boolean mask of rows in ``pts`` satisfying ``A_ineq @ x >= b_ineq + tol``."""
    if a_ineq.shape[0] == 0:
        return np.ones(pts.shape[0], dtype=bool)
    slacks = pts @ a_ineq.T - b_ineq[None, :]
    return np.all(slacks >= tol, axis=1)


def _sir_posterior_nd(
    rng: np.random.Generator,
    *,
    prior_draw: Any,
    prior_logpdf: Any,
    is_flat_prior: bool = False,
    mu_lik: np.ndarray | None = None,
    mu_hat: np.ndarray | None = None,
    sigma_lik: np.ndarray,
    a_ineq: np.ndarray | None = None,
    b_ineq: np.ndarray | None = None,
    n_draws: int,
    n_prop: int,
    lik_logpdf: Any = None,
    return_log_evidence: bool = False,
) -> Any:
    """Defensive-mixture SIR in ``d = A - 1`` dimensions on the linear polyhedron support.

    ``lik_logpdf`` optionally replaces the Gaussian ``N(mu_lik, sigma_lik)`` target likelihood
    (for example by a scale mixture); ``N(mu_lik, sigma_lik)`` remains the proposal. With
    ``return_log_evidence=True`` a fourth element is returned: the deterministic-mixture
    importance-sampling estimate of ``log int p(mu) L(mu) dmu`` (both densities normalized),
    ``-inf`` when no proposal has positive weight.
    """
    if mu_lik is None:
        if mu_hat is None:
            raise ValueError("Provide mu_lik or mu_hat")
        mu_lik = mu_hat
    d_dim = int(mu_lik.shape[0])
    if a_ineq is None or b_ineq is None:
        a_ineq = np.zeros((0, d_dim), dtype=np.float64)
        b_ineq = np.zeros(0, dtype=np.float64)
    n_half = max(n_prop // 2, 1000)
    lik_dist = sp_stats.multivariate_normal(mean=mu_lik, cov=sigma_lik, allow_singular=True)

    def _valid_support(cand_pts: np.ndarray) -> np.ndarray:
        mask = _inside_polyhedron(cand_pts, a_ineq, b_ineq)
        if a_ineq.shape[0] == 0:
            mask &= np.isfinite(prior_logpdf(cand_pts))
        return mask

    if is_flat_prior:
        wide_cov = 4.0 * sigma_lik
        wide_dist = sp_stats.multivariate_normal(mean=mu_lik, cov=wide_cov, allow_singular=True)
        pri_pts_list = []
        n_tried_a = 0
        n_kept_a = 0
        for _ in range(30):
            cand = rng.multivariate_normal(mu_lik, wide_cov, size=n_half)
            if cand.ndim == 1:
                cand = cand[:, None]
            keep = _valid_support(cand)
            n_tried_a += cand.shape[0]
            n_kept_a += int(np.sum(keep))
            pri_pts_list.append(cand[keep])
            if n_kept_a >= n_half:
                break
        pri = (
            np.concatenate(pri_pts_list, axis=0)[:n_half]
            if pri_pts_list
            else np.zeros((0, d_dim), dtype=np.float64)
        )
        z_a = max(n_kept_a / max(n_tried_a, 1), 1e-12)
    else:
        pri_all = prior_draw(n_half)
        keep_a = _valid_support(pri_all)
        z_a = max(float(np.mean(keep_a)), 1e-12)
        pri = pri_all[keep_a]

    lik_pts_list = []
    n_tried_b = 0
    n_kept_b = 0
    for _ in range(40):
        cand = rng.multivariate_normal(mu_lik, sigma_lik, size=n_half)
        if cand.ndim == 1:
            cand = cand[:, None]
        keep = _valid_support(cand)
        n_tried_b += cand.shape[0]
        n_kept_b += int(np.sum(keep))
        lik_pts_list.append(cand[keep])
        if n_kept_b >= n_half:
            break
    lik = (
        np.concatenate(lik_pts_list, axis=0)[:n_half]
        if lik_pts_list
        else np.zeros((0, d_dim), dtype=np.float64)
    )
    z_b = max(n_kept_b / max(n_tried_b, 1), 1e-12)

    pts = np.concatenate([pri, lik], axis=0)
    if pts.shape[0] < 10:
        if return_log_evidence:
            return np.zeros((0, d_dim), dtype=np.float64), pri, 0.0, -np.inf
        return np.zeros((0, d_dim), dtype=np.float64), pri, 0.0

    log_prior = prior_logpdf(pts)
    log_lik = np.asarray(lik_dist.logpdf(pts), dtype=np.float64).reshape(-1)
    log_target_lik = (
        log_lik
        if lik_logpdf is None
        else np.asarray(lik_logpdf(pts), dtype=np.float64).reshape(-1)
    )
    w_a = pri.shape[0] / float(pts.shape[0])
    w_b = 1.0 - w_a

    if is_flat_prior:
        log_q_a = np.asarray(wide_dist.logpdf(pts), dtype=np.float64).reshape(-1)
    else:
        log_q_a = log_prior

    if w_b <= 1e-12:
        log_q = log_q_a - math.log(z_a)
    elif w_a <= 1e-12:
        log_q = log_lik - math.log(z_b)
    else:
        log_q = np.logaddexp(
            math.log(w_a) + log_q_a - math.log(z_a),
            math.log(w_b) + log_lik - math.log(z_b),
        )
    log_w = log_prior + log_target_lik - log_q
    log_w = np.where(np.isfinite(log_w), log_w, -np.inf)
    max_log_w = float(np.max(log_w))
    if not math.isfinite(max_log_w):
        if return_log_evidence:
            return np.zeros((0, d_dim), dtype=np.float64), pri, 0.0, -np.inf
        return np.zeros((0, d_dim), dtype=np.float64), pri, 0.0
    log_w -= max_log_w
    weights = np.exp(log_w)
    sum_w = float(np.sum(weights))
    weights /= sum_w
    ess = float(1.0 / np.sum(weights**2))
    idx = rng.choice(pts.shape[0], size=n_draws, replace=True, p=weights)
    if return_log_evidence:
        log_z = max_log_w + math.log(sum_w) - math.log(float(pts.shape[0]))
        return pts[idx], pri, ess, log_z
    return pts[idx], pri, ess


def _find_interior_starting_points(
    rng: np.random.Generator,
    prior_draw: Any,
    prior_logpdf: Any,
    mu_lik: np.ndarray,
    sigma_lik: np.ndarray,
    a_ineq: np.ndarray,
    b_ineq: np.ndarray,
    num_chains: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Find ``num_chains`` strictly interior starting points and a preconditioning Cholesky factor."""
    d_dim = int(mu_lik.shape[0])
    lik_dist = sp_stats.multivariate_normal(mean=mu_lik, cov=sigma_lik, allow_singular=True)
    candidates: list[np.ndarray] = []

    for _ in range(10):
        c_pri = prior_draw(2000)
        ok_pri = _inside_polyhedron(c_pri, a_ineq, b_ineq, tol=1e-7)
        if np.any(ok_pri):
            candidates.append(c_pri[ok_pri])
        c_lik = rng.multivariate_normal(mu_lik, sigma_lik, size=2000)
        if c_lik.ndim == 1:
            c_lik = c_lik[:, None]
        ok_lik = _inside_polyhedron(c_lik, a_ineq, b_ineq, tol=1e-7)
        if np.any(ok_lik):
            candidates.append(c_lik[ok_lik])
        if sum(arr.shape[0] for arr in candidates) >= 500:
            break

    if not candidates:
        # Solve a linear program to find a strictly interior point maximizing minimum slack
        c_lp = np.zeros(d_dim + 1, dtype=np.float64)
        c_lp[-1] = -1.0
        a_ub = np.column_stack([-a_ineq, np.ones(a_ineq.shape[0], dtype=np.float64)])
        b_ub = -b_ineq
        bounds = [(None, None)] * d_dim + [(1e-6, 1.0)]
        res_lp = sp_optimize.linprog(c_lp, A_ub=a_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if not res_lp.success or not _inside_polyhedron(res_lp.x[:d_dim][None, :], a_ineq, b_ineq, tol=1e-9)[0]:
            raise ValueError(
                "The prior and tau_clip support contain no feasible interior points; widen tau_clip"
            )
        center = res_lp.x[:d_dim]
        chol = np.linalg.cholesky(sigma_lik + 1e-6 * np.eye(d_dim))
        starts = np.repeat(center[None, :], num_chains, axis=0)
        return starts, chol

    pool = np.concatenate(candidates, axis=0)
    log_p = prior_logpdf(pool) + np.asarray(lik_dist.logpdf(pool), dtype=np.float64).reshape(-1)
    finite_mask = np.isfinite(log_p)
    if not np.any(finite_mask):
        raise ValueError(
            "The prior and the likelihood put no finite density inside tau_clip; widen tau_clip or prior_scale"
        )
    pool = pool[finite_mask]
    log_p = log_p[finite_mask]
    log_p -= float(np.max(log_p))
    w = np.exp(log_p)
    w /= float(np.sum(w))

    chosen_idx = rng.choice(pool.shape[0], size=num_chains, replace=pool.shape[0] < num_chains, p=w)
    starts = pool[chosen_idx].copy()

    if pool.shape[0] >= d_dim + 2 and float(1.0 / np.sum(w**2)) >= d_dim + 2:
        emp_cov = np.cov(pool, rowvar=False, aweights=w)
        emp_cov = np.atleast_2d(emp_cov) + 1e-5 * sigma_lik + 1e-8 * np.eye(d_dim)
    else:
        emp_cov = sigma_lik + 1e-6 * np.eye(d_dim)
    chol = np.linalg.cholesky(0.5 * (emp_cov + emp_cov.T))
    return starts, chol


def _slice_posterior_nd(
    rng: np.random.Generator,
    *,
    prior_draw: Any,
    prior_logpdf: Any,
    mu_lik: np.ndarray,
    sigma_lik: np.ndarray,
    a_ineq: np.ndarray,
    b_ineq: np.ndarray,
    n_draws: int,
    num_chains: int = 4,
    num_burnin: int = 200,
    thin: int = 2,
    lik_logpdf: Any = None,
    auto_continue: bool = False,
    target_rhat: float = 1.01,
    target_ess: float = 400.0,
    max_iter_multiplier: float = 8.0,
    return_info: bool = False,
) -> Any:
    """Vectorized multi-chain analytical-polyhedron hit-and-run slice sampler in ``R^d``.

    With ``auto_continue=False`` (the 0.9.0 behavior) the chains run ``num_burnin`` sweeps
    plus ``ceil(n_draws / num_chains) * thin`` recorded sweeps. With ``auto_continue=True``
    the same first block is run (so a run that is already converged is bit-identical to
    0.9.0) and further blocks of the same length are appended until the returned draws
    reach split R-hat ``<= target_rhat`` and bulk ESS ``>= min(target_ess, 0.75 * n_draws)``
    on every coordinate, or the total number of sweeps reaches ``max_iter_multiplier`` times
    the default, in which case a :class:`RankKConfigWarning` is emitted. After ``k >= 2``
    blocks the first ``floor(k / 2)`` blocks are discarded as extra warm-up and the
    remaining draws are thinned evenly back to ``ceil(n_draws / num_chains)`` per chain.

    Args:
      rng: Random generator.
      prior_draw: Callable ``n -> (n, d)`` prior draws.
      prior_logpdf: Callable ``(N, d) -> (N,)`` prior log density.
      mu_lik: ``(d,)`` Gaussian likelihood mean.
      sigma_lik: ``(d, d)`` Gaussian likelihood covariance.
      a_ineq: ``(R, d)`` support inequalities ``a_ineq @ mu >= b_ineq``.
      b_ineq: ``(R,)`` inequality offsets.
      n_draws: Number of pooled posterior draws to return.
      num_chains: Number of parallel chains.
      num_burnin: Burn-in sweeps per chain.
      thin: Thinning interval of the first block.
      lik_logpdf: Optional replacement of the Gaussian ``N(mu_lik, sigma_lik)`` log likelihood
        in the target; the Gaussian is still used to initialize and precondition the chains.
      auto_continue: Continue in blocks until the convergence targets are met.
      target_rhat: Maximum split R-hat target for ``auto_continue``.
      target_ess: Minimum bulk ESS target for ``auto_continue``.
      max_iter_multiplier: Cap on total sweeps as a multiple of the default sweeps.
      return_info: Also return a dict with ``n_iterations`` (sweeps per chain),
        ``n_blocks`` and ``converged``.

    Returns:
      ``(pooled, prior_inside, ess_per_param, rhat_per_param, n_evals)`` and, when
      ``return_info`` is set, the info dict as a sixth element.
    """
    d_dim = int(mu_lik.shape[0])
    n_chains = max(int(num_chains), 2)
    draws_per_chain = max(int(math.ceil(n_draws / float(n_chains))), 4)
    total_sweeps = int(num_burnin) + draws_per_chain * max(int(thin), 1)
    prec_lik = np.linalg.inv(sigma_lik)

    # Probe prior draws to check if the prior lives on a radial order-polytope cone
    # {mu : 0 <= mu_j <= mu_{top_col}} (where Cartesian density has a tot^{-(d-1)} singularity at 0).
    probe_pri = prior_draw(2048)
    if a_ineq.shape[0] > 0:
        probe_ok = _inside_polyhedron(probe_pri, a_ineq, b_ineq, tol=0.0)
        probe_valid = probe_pri[probe_ok] if np.any(probe_ok) else probe_pri
    else:
        probe_ok = np.isfinite(prior_logpdf(probe_pri))
        probe_valid = probe_pri[probe_ok] if np.any(probe_ok) else probe_pri

    use_radial_cone = False
    top_col = d_dim - 1
    interior_cols: list[int] = []
    if d_dim >= 2 and probe_valid.shape[0] >= 16 and np.all(probe_valid > 0.0):
        cand_top = int(np.argmax(np.mean(probe_valid, axis=0)))
        cand_int = [c for c in range(d_dim) if c != cand_top]
        if np.all(probe_valid[:, cand_int] <= probe_valid[:, [cand_top]] * (1.0 + 1e-11)):
            test_pt = np.mean(probe_valid, axis=0, keepdims=True).copy()
            test_pt[0, cand_int[0]] = 2.0 * test_pt[0, cand_top] + 1.0
            outside_ineq = a_ineq.shape[0] > 0 and not _inside_polyhedron(test_pt, a_ineq, b_ineq)[0]
            outside_logp = not np.isfinite(float(prior_logpdf(test_pt)[0]))
            if outside_ineq or outside_logp:
                use_radial_cone = True
                top_col = cand_top
                interior_cols = cand_int

    def _mu_to_z(mu_pts: np.ndarray) -> np.ndarray:
        if not use_radial_cone:
            return mu_pts
        tot = np.maximum(mu_pts[:, top_col], 1e-15)
        w = np.clip(mu_pts[:, interior_cols] / tot[:, None], 1e-12, 1.0 - 1e-12)
        return np.column_stack([w, tot])

    def _z_to_mu(z_pts: np.ndarray) -> np.ndarray:
        if not use_radial_cone:
            return z_pts
        tot = z_pts[:, -1]
        mu_out = np.empty((z_pts.shape[0], d_dim), dtype=np.float64)
        mu_out[:, top_col] = tot
        mu_out[:, interior_cols] = z_pts[:, :-1] * tot[:, None]
        return mu_out

    if use_radial_cone:
        m_int = d_dim - 1
        rows_az: list[np.ndarray] = []
        rows_bz: list[float] = []
        for k in range(m_int):
            ek = np.zeros(d_dim, dtype=np.float64)
            ek[k] = 1.0
            rows_az.append(ek)
            rows_bz.append(0.0)
            rows_az.append(-ek)
            rows_bz.append(-1.0)
        etot = np.zeros(d_dim, dtype=np.float64)
        etot[-1] = 1.0
        rows_az.append(etot)
        rows_bz.append(0.0)

        # Include any pairwise ordering w_k1 <= w_k2 enforced by the prior
        mean_z = _mu_to_z(np.mean(probe_valid, axis=0, keepdims=True))[0]
        for k1 in range(m_int):
            for k2 in range(m_int):
                if k1 == k2:
                    continue
                if np.all(probe_valid[:, interior_cols[k1]] <= probe_valid[:, interior_cols[k2]] + 1e-12):
                    z_test = mean_z.copy()
                    z_test[k1] = min(0.8, max(mean_z[k2] + 0.1, 0.6))
                    z_test[k2] = max(0.2, min(mean_z[k1] - 0.1, 0.4))
                    mu_test = _z_to_mu(z_test[None, :])
                    viol_ineq = a_ineq.shape[0] > 0 and not _inside_polyhedron(mu_test, a_ineq, b_ineq)[0]
                    viol_logp = not np.isfinite(float(prior_logpdf(mu_test)[0]))
                    if viol_ineq or viol_logp:
                        diff_row = np.zeros(d_dim, dtype=np.float64)
                        diff_row[k2] = 1.0
                        diff_row[k1] = -1.0
                        rows_az.append(diff_row)
                        rows_bz.append(0.0)

        for r_idx in range(a_ineq.shape[0]):
            a_r = a_ineq[r_idx]
            b_r = float(b_ineq[r_idx])
            if abs(b_r) <= 1e-12:
                z_row = np.zeros(d_dim, dtype=np.float64)
                for k, col_idx in enumerate(interior_cols):
                    z_row[k] = a_r[col_idx]
                rows_az.append(z_row)
                rows_bz.append(-float(a_r[top_col]))
            elif np.all(np.abs(a_r[interior_cols]) <= 1e-12) and abs(a_r[top_col]) > 1e-12:
                z_row = np.zeros(d_dim, dtype=np.float64)
                z_row[-1] = a_r[top_col]
                rows_az.append(z_row)
                rows_bz.append(b_r)

        a_work = np.vstack(rows_az)
        b_work = np.asarray(rows_bz, dtype=np.float64)
    else:
        a_work = a_ineq
        b_work = b_ineq

    def _log_target(z_pts: np.ndarray) -> np.ndarray:
        inside = _inside_polyhedron(z_pts, a_work, b_work, tol=0.0)
        if use_radial_cone:
            inside &= z_pts[:, -1] > 0.0
        out = np.full(z_pts.shape[0], -np.inf, dtype=np.float64)
        if not np.any(inside):
            return out
        z_valid = z_pts[inside]
        mu_valid = _z_to_mu(z_valid)
        if a_ineq.shape[0] > 0 and not use_radial_cone:
            pass
        elif a_ineq.shape[0] > 0:
            in_mu = _inside_polyhedron(mu_valid, a_ineq, b_ineq, tol=0.0)
            if not np.all(in_mu):
                idx_in = np.nonzero(inside)[0]
                inside[idx_in[~in_mu]] = False
                if not np.any(inside):
                    return out
                z_valid = z_pts[inside]
                mu_valid = _z_to_mu(z_valid)
        lp = prior_logpdf(mu_valid)
        if use_radial_cone:
            lp = lp + float(d_dim - 1) * np.log(np.maximum(z_valid[:, -1], 1e-300))
        if lik_logpdf is not None:
            ll = np.asarray(lik_logpdf(mu_valid), dtype=np.float64).reshape(-1)
        else:
            diff = mu_valid - mu_lik[None, :]
            ll = -0.5 * np.sum((diff @ prec_lik) * diff, axis=1)
        out[inside] = lp + ll
        return out

    # Initialize chains via importance resampling of prior draws weighted by the likelihood
    if probe_valid.shape[0] >= n_chains:
        z_cand = _mu_to_z(probe_valid)
        ok_z = _inside_polyhedron(z_cand, a_work, b_work, tol=1e-9)
        if np.any(ok_z):
            z_pool = z_cand[ok_z]
            mu_pool = probe_valid[ok_z]
            diff_pool = mu_pool - mu_lik[None, :]
            ll_pool = -0.5 * np.sum((diff_pool @ prec_lik) * diff_pool, axis=1)
            lp_pool = _log_target(z_pool)
            fin = np.isfinite(ll_pool) & np.isfinite(lp_pool)
            if np.any(fin):
                z_pool = z_pool[fin]
                ll_pool = ll_pool[fin]
                ll_pool -= float(np.max(ll_pool))
                w_pool = np.exp(ll_pool)
                w_pool /= float(np.sum(w_pool))
                chosen = rng.choice(
                    z_pool.shape[0], size=n_chains, replace=z_pool.shape[0] < n_chains, p=w_pool
                )
                x_cur = z_pool[chosen].copy()
                if z_pool.shape[0] >= d_dim + 2 and float(1.0 / np.sum(w_pool**2)) >= d_dim + 2:
                    emp_cov = np.atleast_2d(np.cov(z_pool, rowvar=False, aweights=w_pool))
                else:
                    emp_cov = np.atleast_2d(np.cov(z_pool, rowvar=False))
                emp_cov = 0.5 * (emp_cov + emp_cov.T) + 1e-4 * np.eye(d_dim)
                chol = np.linalg.cholesky(emp_cov)
            else:
                x_cur, chol = _find_interior_starting_points(
                    rng, prior_draw, prior_logpdf, mu_lik, sigma_lik, a_ineq, b_ineq, n_chains
                )
                x_cur = _mu_to_z(x_cur)
        else:
            x_cur, chol = _find_interior_starting_points(
                rng, prior_draw, prior_logpdf, mu_lik, sigma_lik, a_ineq, b_ineq, n_chains
            )
            x_cur = _mu_to_z(x_cur)
    else:
        x_cur, chol = _find_interior_starting_points(
            rng, prior_draw, prior_logpdf, mu_lik, sigma_lik, a_ineq, b_ineq, n_chains
        )
        x_cur = _mu_to_z(x_cur)

    log_cur = _log_target(x_cur)
    base_scale = float(np.sqrt(np.mean(np.diag(chol @ chol.T))))
    n_evals = int(n_chains)

    chain_traces = np.zeros((n_chains, draws_per_chain, d_dim), dtype=np.float64)
    kept_idx = 0

    # Auto-continuation state: extra blocks of ``draws_per_chain * thin`` sweeps are appended
    # lazily by ``_sweep_schedule`` until the selected draws meet the convergence targets.
    n_thin = max(int(thin), 1)
    block_sweeps = draws_per_chain * n_thin
    max_sweeps = max(int(math.floor(float(max_iter_multiplier) * total_sweeps)), total_sweeps)
    ess_goal = min(float(target_ess), 0.75 * float(n_draws))
    extra_traces: list[np.ndarray] = []
    auto_state: dict[str, Any] = {"converged": None, "n_sweeps": total_sweeps}
    compute_rhat, compute_ess = _get_rhat_ess()

    def _select_traces() -> np.ndarray:
        blocks = [chain_traces, *extra_traces]
        if len(blocks) == 1:
            return chain_traces
        pool = np.concatenate(blocks[len(blocks) // 2 :], axis=1)
        sel = np.round(np.linspace(0.0, pool.shape[1] - 1.0, draws_per_chain)).astype(np.int64)
        return pool[:, sel, :]

    def _targets_met() -> bool:
        sel_traces = _select_traces()
        rh = [float(compute_rhat(sel_traces[:, :, j])) for j in range(d_dim)]
        es = [float(compute_ess(sel_traces[:, :, j])) for j in range(d_dim)]
        ok = all(math.isfinite(r) and r <= float(target_rhat) for r in rh) and all(
            math.isfinite(e) and e >= ess_goal for e in es
        )
        auto_state["converged"] = bool(ok)
        return bool(ok)

    def _sweep_schedule() -> Iterator[int]:
        yield from range(total_sweeps)
        if not auto_continue:
            return
        s_idx = total_sweeps
        while not _targets_met() and s_idx + block_sweeps <= max_sweeps:
            extra_traces.append(np.zeros((n_chains, draws_per_chain, d_dim), dtype=np.float64))
            for _ in range(block_sweeps):
                yield s_idx
                s_idx += 1
            auto_state["n_sweeps"] = s_idx

    for sweep in _sweep_schedule():
        for dir_step in range(d_dim + 1):
            if dir_step < d_dim:
                raw_z = rng.standard_normal((n_chains, d_dim))
                dirs = raw_z @ chol.T
            else:
                axis = sweep % d_dim
                dirs = np.zeros((n_chains, d_dim), dtype=np.float64)
                dirs[:, axis] = 1.0
            norms = np.linalg.norm(dirs, axis=1, keepdims=True)
            dirs = dirs / np.maximum(norms, 1e-12)

            if a_work.shape[0] > 0:
                slacks = x_cur @ a_work.T - b_work[None, :]
                rates = dirs @ a_work.T
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratios = -slacks / rates
                upper_candidates = np.where(rates < -1e-14, ratios, np.inf)
                lower_candidates = np.where(rates > 1e-14, ratios, -np.inf)
                t_max = np.min(upper_candidates, axis=1)
                t_min = np.max(lower_candidates, axis=1)
                chord = np.maximum(t_max - t_min, 0.0)
                pad = np.minimum(1e-11, 0.001 * chord)
                t_min = t_min + pad
                t_max = t_max - pad
            else:
                t_min = np.full(n_chains, -np.inf, dtype=np.float64)
                t_max = np.full(n_chains, np.inf, dtype=np.float64)

            log_y = log_cur - rng.exponential(1.0, size=n_chains)
            finite_chord = np.isfinite(t_min) & np.isfinite(t_max)
            w_step = np.where(
                finite_chord,
                np.maximum(0.5 * (t_max - t_min), 1e-8),
                2.0 * base_scale,
            )
            u_init = rng.uniform(0.0, 1.0, size=n_chains)
            left = np.maximum(-u_init * w_step, t_min)
            right = np.minimum(left + w_step, t_max)

            stepping_l = left > t_min + 1e-12
            for _ in range(12):
                if not np.any(stepping_l):
                    break
                act_l = np.nonzero(stepping_l)[0]
                pts_l = x_cur[act_l] + left[act_l, None] * dirs[act_l]
                log_l = _log_target(pts_l)
                n_evals += int(pts_l.shape[0])
                expand_l = log_l > log_y[act_l]
                still_l = act_l[expand_l]
                stepping_l[:] = False
                if still_l.size == 0:
                    break
                left[still_l] = np.maximum(left[still_l] - w_step[still_l], t_min[still_l])
                stepping_l[still_l] = left[still_l] > t_min[still_l] + 1e-12

            stepping_r = right < t_max - 1e-12
            for _ in range(12):
                if not np.any(stepping_r):
                    break
                act_r = np.nonzero(stepping_r)[0]
                pts_r = x_cur[act_r] + right[act_r, None] * dirs[act_r]
                log_r = _log_target(pts_r)
                n_evals += int(pts_r.shape[0])
                expand_r = log_r > log_y[act_r]
                still_r = act_r[expand_r]
                stepping_r[:] = False
                if still_r.size == 0:
                    break
                right[still_r] = np.minimum(right[still_r] + w_step[still_r], t_max[still_r])
                stepping_r[still_r] = right[still_r] < t_max[still_r] - 1e-12

            active = np.ones(n_chains, dtype=bool)
            for _ in range(60):
                if not np.any(active):
                    break
                act_idx = np.nonzero(active)[0]
                t_prop = rng.uniform(left[act_idx], right[act_idx])
                x_prop = x_cur[act_idx] + t_prop[:, None] * dirs[act_idx]
                log_prop = _log_target(x_prop)
                n_evals += int(act_idx.size)
                accepted = log_prop >= log_y[act_idx]
                if np.any(accepted):
                    acc_chains = act_idx[accepted]
                    x_cur[acc_chains] = x_prop[accepted]
                    log_cur[acc_chains] = log_prop[accepted]
                    active[acc_chains] = False
                if np.any(~accepted):
                    rej_chains = act_idx[~accepted]
                    t_rej = t_prop[~accepted]
                    neg_mask = t_rej < 0.0
                    left[rej_chains[neg_mask]] = t_rej[neg_mask]
                    right[rej_chains[~neg_mask]] = t_rej[~neg_mask]
                    collapsed = (right[rej_chains] - left[rej_chains]) < 1e-13
                    if np.any(collapsed):
                        active[rej_chains[collapsed]] = False

        if sweep >= total_sweeps:
            off_blk, off_pos = divmod(sweep - total_sweeps, block_sweeps)
            if off_pos % n_thin == 0:
                extra_traces[off_blk][:, off_pos // n_thin, :] = _z_to_mu(x_cur)
        elif sweep >= int(num_burnin) and (sweep - int(num_burnin)) % max(int(thin), 1) == 0:
            if kept_idx < draws_per_chain:
                chain_traces[:, kept_idx, :] = _z_to_mu(x_cur)
                kept_idx += 1

    info = {
        "n_iterations": int(auto_state["n_sweeps"]),
        "n_blocks": 1 + len(extra_traces),
        "converged": auto_state["converged"],
    }
    if auto_continue:
        chain_traces = _select_traces()
        if not auto_state["converged"]:
            _warn_rank_k(
                f"hit-and-run slice sampler did not reach split R-hat <= {float(target_rhat):g} and "
                f"bulk ESS >= {ess_goal:.0f} on every coordinate within {info['n_iterations']} sweeps "
                f"per chain ({float(max_iter_multiplier):g}x the default); increase num_burnin, thin "
                "or max_iter_multiplier and check sampler_diagnostics."
            )

    rhats: list[float] = []
    esses: list[float] = []
    for j_idx in range(d_dim):
        coord_chains = chain_traces[:, :, j_idx]
        rhats.append(float(compute_rhat(coord_chains)))
        esses.append(float(compute_ess(coord_chains)))

    pooled = chain_traces.reshape(-1, d_dim)[:n_draws]
    pri_samples = prior_draw(max(n_draws, 1000))
    pri_inside = pri_samples[_inside_polyhedron(pri_samples, a_ineq, b_ineq)]
    if return_info:
        return pooled, pri_inside, tuple(esses), tuple(rhats), n_evals, info
    return pooled, pri_inside, tuple(esses), tuple(rhats), n_evals


def _hit_and_run_slice_sampler(
    rng: np.random.Generator,
    *,
    prior_draw: Any,
    prior_logpdf: Any,
    mu_hat: np.ndarray,
    sigma_lik: np.ndarray,
    n_draws: int,
    num_chains: int = 4,
    num_warmup: int = 200,
    thin: int = 2,
    a_ineq: np.ndarray | None = None,
    b_ineq: np.ndarray | None = None,
    auto_continue: bool = False,
    target_rhat: float = 1.01,
    target_ess: float = 400.0,
    max_iter_multiplier: float = 8.0,
    return_info: bool = False,
) -> Any:
    """Convenience wrapper around :func:`_slice_posterior_nd` returning ``(draws, min_ess, max_rhat)``.

    With ``auto_continue=True`` the chains run in blocks until split R-hat ``<= target_rhat``
    and bulk ESS ``>= target_ess`` (capped at ``0.75 * n_draws``), up to
    ``max_iter_multiplier`` times the default sweeps; ``return_info=True`` appends the
    ``{"n_iterations", "n_blocks", "converged"}`` dict to the returned tuple.
    """
    d_dim = int(mu_hat.shape[0])
    if a_ineq is None or b_ineq is None:
        a_ineq = np.zeros((0, d_dim), dtype=np.float64)
        b_ineq = np.zeros(0, dtype=np.float64)
    pooled, _, esses, rhats, _, info = _slice_posterior_nd(
        rng,
        prior_draw=prior_draw,
        prior_logpdf=prior_logpdf,
        mu_lik=mu_hat,
        sigma_lik=sigma_lik,
        a_ineq=a_ineq,
        b_ineq=b_ineq,
        n_draws=n_draws,
        num_chains=num_chains,
        num_burnin=num_warmup,
        thin=thin,
        auto_continue=auto_continue,
        target_rhat=target_rhat,
        target_ess=target_ess,
        max_iter_multiplier=max_iter_multiplier,
        return_info=True,
    )
    if return_info:
        return pooled, float(min(esses)), float(max(rhats)), info
    return pooled, float(min(esses)), float(max(rhats))


def reconcile_multi_arm_draws(
    results: Mapping[tuple[str, str], RankKLongBetResult],
    arms: Sequence[str] = (),
    *,
    topology: str = "chain",
    edges: Sequence[tuple[str, str]] | None = None,
    prior_mode: str | None = None,
    prior_scale: float | str | Sequence[float] | Mapping[str, float] | None = "in_sample",
    prior_mean: float | Sequence[float] | Mapping[str, float] = 0.0,
    split_alpha: Sequence[float] | float | None = None,
    monotone: bool = True,
    tau_clip: tuple[float, float] | None = (0.0, float("inf")),
    catt_floor: float | None = 0.0,
    cohorts: Mapping[tuple[str, str], np.ndarray] | np.ndarray | None = None,
    ever: Mapping[tuple[str, str], np.ndarray] | np.ndarray | None = None,
    att_w: Mapping[tuple[str, str], np.ndarray] | np.ndarray | None = None,
    arm_sample_sizes: Mapping[str, float] | Sequence[float] | None = None,
    sample_sizes: Mapping[str, float] | Sequence[float] | None = None,
    min_group_var_share: float = 0.15,
    min_arm_var_share: float | None = None,
    bipartite_ridge_weight: float = 1e-2,
    bipartite_ridge: float | None = None,
    bipartite_policy: str = "allocation",
    key_contrast: tuple[str, str] | None = None,
    horizon: int | None = None,
    n_proposals: int | None = None,
    sampler: str = "auto",
    num_chains: int = 4,
    num_burnin: int = 200,
    num_warmup: int | None = None,
    thin: int = 2,
    auto_continue: bool = True,
    target_rhat: float = 1.01,
    target_ess: float = 400.0,
    max_iter_multiplier: float = 8.0,
    seed: int = 2026,
    return_prior_draws: bool = False,
    spike_prob: float = 0.5,
    prior_scale_floor: float | None = None,
    projection: str = "ols",
    _force_general: bool = False,
) -> MultiArmReconciliationResult:
    """Reconcile ``P >= A - 1`` pairwise Rank-K BCF posteriors across ``A >= 2`` arms.

    Implements the general multi-arm Mode-1 and trajectory reconciliation of Section 3.1:
    1. Validates that the contrast graph on ``arms`` is connected and that every spanning-tree
       edge has a fitted contrast in ``results``.
    2. Projects the ``P`` pairwise Mode-1 posterior means onto the ``(A - 1)``-dimensional cut
       space ``mu_hat = (B' B)^{-1} B' y_hat``, recovers arm variances ``v in R_+^A`` via
       bounded NNLS (or sample-size-informed ridge-NNLS with a ``UserWarning`` when the contrast
       graph is bipartite), and forms ``Sigma_lik = diag(v_1, ..., v_{A-1}) + v_0 1 1'``.
    3. Combines the Gaussian cut-space likelihood with the specified prior family (chain
       ``Dirichlet(alpha)`` bundle prior, unique-top poset order-polytope prior via uniform
       linear extensions, independent truncated Gaussian/half-normal prior, or flat prior).
    4. Samples the posterior via defensive-mixture SIR for ``d = A - 1 <= 3`` (with automatic
       fallback to hit-and-run slice sampling if ``ess < max(100, 0.02 * n_proposals)``) or
       multi-chain analytical-polyhedron hit-and-run slice sampling for ``d >= 4``.
    5. Projects each spanning-tree edge increment draw through that edge's fit via
       :func:`_project_cohort_damped_trajectory` and computes every contrast trajectory as the
       exact signed path sum on the spanning tree so cycle identities hold to ``< 1e-12``.

    When ``len(arms) == 3``, ``topology == "chain"``, ``split_alpha`` is ``None`` (or all ``1.0``),
    ``monotone=True``, ``sampler in ("auto", "sir")``, and ``_force_general=False`` with all 3
    pairs present, delegates directly to :func:`reconcile_nested_3arm_draws` for bit-for-bit
    backward compatibility.

    Args:
      results: Mapping ``{(treated_arm, reference_arm): RankKLongBetResult}`` of ``P >= A - 1``
        pairwise fits.
      arms: Sequence of ``A >= 2`` arm names, with ``arms[0]`` as the control arm.
      topology: ``"chain"`` (nested sequence ``arms[0] < arms[1] < ... < arms[A-1]``) or
        ``"poset"`` (partial order specified by ``edges``).
      edges: Cover edges of the Hasse diagram when ``topology="poset"``.
      prior_mode: ``None`` (default: ``"bundle_halfnormal"``; on a bipartite contrast graph the
        bipartite default family), ``"bundle_halfnormal"``, ``"bundle_halfstudentt"``,
        ``"bundle_gamma_spike"``, ``"bundle_spike_halfnormal"`` (spike-and-slab: point mass at
        ``mu = 0`` with probability ``spike_prob``, else the half-normal bundle prior),
        ``"independent"``, ``"independent_halfnormal"``, ``"independent_normal"``, or ``"flat"``.
      prior_scale: Positive float in outcome units, ``"in_sample"`` (the in-sample posterior SD
        of the scale-defining contrast), ``"empirical_bayes"`` (Type-II maximum likelihood of the
        scale under the Gaussian cut-space likelihood of the top contrast, or of each coordinate
        for the ``independent`` family: closed form for half-normal / normal priors, quadrature
        otherwise, maximized over ``sigma >= prior_scale_floor``), or per-arm sequence/mapping for
        the ``independent`` prior family.
      prior_mean: Prior mean(s) for the ``independent`` prior family (default ``0.0``).
      split_alpha: Optional ``(A - 1,)`` concentration parameters for the chain Dirichlet split
        (default ``None``, which uses ``(1.0, ..., 1.0)``).
      monotone: Enforce poset monotonicity ``mu_j >= mu_l`` for ``l <= j`` (default ``True``).
      tau_clip: ``(low, high)`` Mode-1 truncation bounds (default ``(0.0, inf)``; ``None`` for
        ``(-inf, inf)``).
      catt_floor: Optional per-period floor for zero-sum dynamic damping (default ``0.0``).
      cohorts: Optional adoption cohorts array or ``{pair: array}`` mapping.
      ever: Optional ever-treated mask array or ``{pair: array}`` mapping.
      att_w: Optional unit ATT weights array or ``{pair: array}`` mapping.
      arm_sample_sizes: Optional ``{arm: n_a}`` or length-``A`` sequence of arm unit counts used
        in bipartite ridge regularization.
      min_group_var_share: Floor share ``eta`` of the minimum incident pairwise variance per arm
        (default ``0.15``).
      bipartite_ridge_weight: Ridge penalty weight for bipartite variance recovery (default ``1e-2``).
      bipartite_policy: How a bipartite contrast graph (no odd cycle, so the arm variances are
        identified only up to the alternating null direction ``v + xi z``) resolves ``xi``:
        ``"allocation"`` (default; 0.9.0 behavior) ridge-shrinks toward the sample-size
        allocation ``v_a ~ sigma^2 / n_a``; ``"conservative"`` picks the ``xi`` in the feasible
        interval (every ``v_a`` above its floor) that maximizes the Gaussian-approximate
        posterior variance of the top contrast; ``"marginalize"`` integrates ``xi`` under a
        flat prior on the feasible interval (a scale-mixture likelihood on a 9-point grid).
        A :class:`RankKDesignWarning` names the extra fit that would close an odd cycle.
      key_contrast: Optional primary contrast ``(treated, ref)`` (defaults to ``(arms[-1], arms[0])``).
      horizon: Optional number of post-adoption periods summed into ``contrast_cumulative_draws``.
      n_proposals: Number of SIR importance proposals (default ``max(n_draws * 500, 100_000)``).
      sampler: ``"auto"`` (SIR for ``d <= 3`` with slice fallback; slice for ``d >= 4``),
        ``"sir"``, or ``"slice"``.
      num_chains: Number of parallel chains for the slice sampler (default ``4``).
      num_burnin: Burn-in sweeps per chain for the slice sampler (default ``200``).
      auto_continue: Slice sampler only (``d >= 4`` or SIR fallback): continue the chains in
        blocks until split R-hat ``<= target_rhat`` and bulk ESS ``>= target_ess`` (capped at
        ``0.75 * S``) on every coordinate, up to ``max_iter_multiplier`` times the default
        sweeps, then warn (default ``True`` since 0.10.0; ``False`` is bit-identical to 0.9.0).
        Sweeps used are reported in ``sampler_diagnostics.n_iterations``.
      target_rhat: Split R-hat target for ``auto_continue`` (default ``1.01``).
      target_ess: Bulk ESS target for ``auto_continue`` (default ``400``).
      max_iter_multiplier: Cap on slice sweeps relative to the default (default ``8``).
      seed: Random seed.
      return_prior_draws: Attach prior Mode-1 and trajectory draws to the result.
      projection: ``"ols"`` (default; Euclidean cut-space projection, 0.9.0) or ``"huber"``
        (IRLS with Huber weights, ``k = 1.345``, scale ``1.4826 * MAD`` of the cycle residual,
        so one biased pairwise fit cannot drag the whole reconciliation; the down-weighted
        fits are reported in ``flagged_pairs`` / ``projection_weights``). Identical to
        ``"ols"`` under the exact shared-arm model.
      prior_scale_floor: Lower bound of the ``"empirical_bayes"`` scale (default: the in-sample
        posterior SD of the scale-defining contrast, i.e. its noise scale).
      spike_prob: Prior spike probability ``pi`` of ``prior_mode="bundle_spike_halfnormal"``
        (``tau_top ~ pi delta_0 + (1 - pi) HalfNormal(prior_scale)``; default ``0.5``, which is
        also what a ``Beta(1, 1)`` prior on ``pi`` integrates to for one reconciliation). The
        posterior spike weight is exact given the slab evidence and is reported as
        ``spike_posterior_prob``; spike draws are exactly zero on every contrast and period.
      _force_general: Bypass the ``A = 3`` legacy delegation path (for testing equivalence).

    Returns:
      :class:`MultiArmReconciliationResult`.
    """
    _ = key_contrast
    if sample_sizes is not None and arm_sample_sizes is None:
        arm_sample_sizes = sample_sizes
    if min_arm_var_share is not None:
        min_group_var_share = float(min_arm_var_share)
    if bipartite_ridge is not None:
        bipartite_ridge_weight = float(bipartite_ridge)
    if num_warmup is not None:
        num_burnin = int(num_warmup)

    arms_tuple = tuple(str(a) for a in arms)
    n_arms = len(arms_tuple)
    if n_arms < 2 or len(set(arms_tuple)) != n_arms:
        raise ValueError(f"arms must contain at least 2 distinct arm names; got {arms!r}")
    if topology not in _TOPOLOGIES:
        raise ValueError(f"topology must be one of {_TOPOLOGIES}; got {topology!r}")
    if prior_mode is not None and prior_mode not in _ALL_PRIOR_MODES:
        raise ValueError(f"prior_mode must be one of {_ALL_PRIOR_MODES}; got {prior_mode!r}")
    if sampler not in _SAMPLERS:
        raise ValueError(f"sampler must be one of {_SAMPLERS}; got {sampler!r}")
    if projection not in _PROJECTIONS:
        raise ValueError(f"projection must be one of {_PROJECTIONS}; got {projection!r}")
    if bipartite_policy not in _BIPARTITE_POLICIES:
        raise ValueError(
            f"bipartite_policy must be one of {_BIPARTITE_POLICIES}; got {bipartite_policy!r}"
        )

    resolved_tau_clip = (
        (-float("inf"), float("inf"))
        if tau_clip is None
        else (float(tau_clip[0]), float(tau_clip[1]))
    )
    if not resolved_tau_clip[0] < resolved_tau_clip[1]:
        raise ValueError(f"tau_clip must satisfy low < high; got {tau_clip!r}")

    results_map = _normalize_results_mapping(results, arms_tuple)
    arm_to_idx = {arm: idx for idx, arm in enumerate(arms_tuple)}
    fitted_pairs = tuple(results_map.keys())
    n_pairs = len(fitted_pairs)
    topology_requested = topology
    topology, edges = _chain_order_as_poset(arms_tuple, topology, edges, fitted_pairs, split_alpha)

    # Validate draw counts across all fitted contrasts
    draw_counts: set[int] = set()
    for pair, fit in results_map.items():
        att_arr = np.asarray(fit.att_draws, dtype=np.float64)
        if att_arr.ndim != 1 or att_arr.size < 2 or not np.all(np.isfinite(att_arr)):
            raise ValueError(f"results[{pair!r}].att_draws must be a finite 1-D array with >= 2 draws")
        draw_counts.add(int(att_arr.size))
    if len(draw_counts) != 1:
        raise ValueError(
            f"Inconsistent draw counts across fitted contrasts: {sorted(draw_counts)}"
        )
    (n_draws,) = draw_counts

    pair_indices = [(arm_to_idx[t_arm], arm_to_idx[r_arm]) for t_arm, r_arm in fitted_pairs]
    is_connected, is_bipartite, colors = _check_connected_and_bipartite(n_arms, pair_indices)
    if not is_connected:
        raise ValueError(
            f"Contrast graph on arms {arms_tuple} is disconnected (not connected); fitted pairs are {fitted_pairs}"
        )
    if prior_mode is None:
        prior_mode = _BIPARTITE_DEFAULT_PRIOR_MODE if is_bipartite else _DEFAULT_PRIOR_MODE
    prior_mode_requested = prior_mode
    spike_mode = prior_mode in _SPIKE_PRIOR_MODES
    if spike_mode:
        if not 0.0 < float(spike_prob) < 1.0:
            raise ValueError(f"spike_prob must lie in (0, 1); got {spike_prob!r}")
        if resolved_tau_clip[0] > 0.0:
            raise ValueError(
                f"prior_mode={prior_mode!r} puts a spike at zero, which needs tau_clip[0] <= 0; "
                f"got tau_clip={tau_clip!r}"
            )
        prior_mode = _SPIKE_SLAB_MODE[prior_mode]

    tree_edges, redundant_pairs, cover_pairs, reachability, maximal_arms = _validate_and_orient_poset(
        arms_tuple, topology, edges, results_map
    )

    # Bundle priors require a unique top arm in the poset
    if prior_mode in _BUNDLE_PRIOR_MODES and len(maximal_arms) != 1:
        max_names = [arms_tuple[i] for i in maximal_arms]
        raise ValueError(
            f"Non-unique top arm in poset (maximal arms: {max_names}): bundle prior mode "
            f"{prior_mode!r} requires a unique top arm; use prior_mode='independent' instead"
        )
    top_idx = maximal_arms[0] if len(maximal_arms) == 1 else None

    d_dim = n_arms - 1
    if split_alpha is not None:
        if topology_requested != "chain":
            raise ValueError("split_alpha is only supported when topology='chain'")
        if np.ndim(split_alpha) == 0:
            alpha_vec = np.full(d_dim, float(split_alpha), dtype=np.float64)
        else:
            alpha_vec = np.asarray(split_alpha, dtype=np.float64)
        if alpha_vec.shape != (d_dim,) or np.any(~np.isfinite(alpha_vec)) or np.any(alpha_vec <= 0.0):
            raise ValueError(
                f"split_alpha must be a sequence of {d_dim} positive finite floats; got {split_alpha!r}"
            )
    else:
        alpha_vec = np.ones(d_dim, dtype=np.float64)

    # Step 2: Mode-1 cut-space projection and arm variance recovery
    e_inc = np.zeros((n_pairs, n_arms), dtype=np.float64)
    y_hat = np.zeros(n_pairs, dtype=np.float64)
    s_sq = np.zeros(n_pairs, dtype=np.float64)
    for p_idx, ((t_arm, r_arm), (j_idx, l_idx)) in enumerate(zip(fitted_pairs, pair_indices, strict=True)):
        fit = results_map[(t_arm, r_arm)]
        d_arr = np.asarray(fit.att_draws, dtype=np.float64)
        e_inc[p_idx, j_idx] = 1.0
        e_inc[p_idx, l_idx] = -1.0
        y_hat[p_idx] = float(np.mean(d_arr))
        s_sq[p_idx] = max(float(np.var(d_arr, ddof=1)), 1e-12)

    b_cut = e_inc[:, 1:]  # (P, A - 1)
    btb = b_cut.T @ b_cut
    mu_hat = np.linalg.solve(btb, b_cut.T @ y_hat)
    y_proj = b_cut @ mu_hat
    cycle_discrepancy_ratio = float(
        np.linalg.norm(y_hat - y_proj) / max(float(np.linalg.norm(y_hat)), 1e-12)
    )
    projection_weights: dict[tuple[str, str], float] | None = None
    flagged_pairs: tuple[tuple[str, str], ...] = ()
    if projection == "huber":
        mu_hat, w_proj, _ = _huber_cut_space_projection(b_cut, y_hat, np.sqrt(s_sq))
        projection_weights = {pair: float(w) for pair, w in zip(fitted_pairs, w_proj, strict=True)}
        order = np.argsort(w_proj, kind="stable")
        flagged_pairs = tuple(fitted_pairs[i] for i in order if w_proj[i] < 1.0 - 1e-12)

    sample_sizes_vec = _infer_arm_sample_sizes(arms_tuple, results_map, arm_sample_sizes)
    arm_variances, nnls_residual = _recover_arm_variances(
        arms_tuple,
        pair_indices,
        s_sq,
        is_bipartite=is_bipartite,
        colors=colors,
        arm_sample_sizes=sample_sizes_vec,
        min_group_var_share=min_group_var_share,
        bipartite_ridge_weight=bipartite_ridge_weight,
        bipartite_policy=bipartite_policy,
    )
    sigma_lik = np.diag(np.asarray(arm_variances[1:], dtype=np.float64)) + float(
        arm_variances[0]
    ) * np.ones((d_dim, d_dim), dtype=np.float64)

    # Helper to look up oriented fit and sign for any pair (t_arm, r_arm)
    undir_to_oriented: dict[frozenset[str], tuple[tuple[str, str], RankKLongBetResult, float]] = {}
    for (t_arm, r_arm), fit in results_map.items():
        undir_to_oriented[frozenset((t_arm, r_arm))] = ((t_arm, r_arm), fit, 1.0)

    def _get_pair_arg(arg_map: Any, pair: tuple[str, str]) -> Any:
        if isinstance(arg_map, Mapping):
            if pair in arg_map:
                return arg_map[pair]
            rev = (pair[1], pair[0])
            if rev in arg_map:
                return arg_map[rev]
            return None
        return arg_map

    # Check exact delegation to legacy reconcile_nested_3arm_draws at A = 3 chain
    can_delegate_3arm = (
        not _force_general
        and n_arms == 3
        and n_pairs == 3
        and topology == "chain"
        and bool(monotone)
        and np.allclose(alpha_vec, 1.0)
        and sampler in ("auto", "sir")
        and prior_mode in (*_BUNDLE_PRIOR_MODES, "independent_halfnormal")
        and not spike_mode
        and (arms_tuple[1], arms_tuple[0]) in results_map
        and (arms_tuple[2], arms_tuple[1]) in results_map
        and (arms_tuple[2], arms_tuple[0]) in results_map
        and isinstance(prior_scale, (int, float, str, type(None)))
        and prior_scale != "empirical_bayes"
        and projection == "ols"
        and isinstance(prior_mean, (int, float))
    )
    if can_delegate_3arm:
        c_arm, b_arm, a_arm = arms_tuple
        pair_bc = (b_arm, c_arm)
        pair_ab = (a_arm, b_arm)
        pair_ac = (a_arm, c_arm)
        leg_scale: float | str | None = (
            prior_scale if isinstance(prior_scale, (int, float, str)) else None
        )
        leg = reconcile_nested_3arm_draws(
            results_map[pair_bc],
            results_map[pair_ab],
            results_map[pair_ac],
            prior_scale=leg_scale,
            prior_mode=prior_mode,
            cohorts_bc=_get_pair_arg(cohorts, pair_bc),
            ever_bc=_get_pair_arg(ever, pair_bc),
            att_w_bc=_get_pair_arg(att_w, pair_bc),
            cohorts_ab=_get_pair_arg(cohorts, pair_ab),
            ever_ab=_get_pair_arg(ever, pair_ab),
            att_w_ab=_get_pair_arg(att_w, pair_ab),
            prior_mean_bc=float(prior_mean),
            prior_mean_ab=float(prior_mean),
            tau_clip=resolved_tau_clip,
            catt_floor=catt_floor,
            min_group_var_share=min_group_var_share,
            n_proposals=n_proposals,
            seed=seed,
            return_prior_draws=return_prior_draws,
        )
        t_post = int(leg.trt_catt_ac.shape[1])
        h_eff = t_post if horizon is None else int(horizon)
        if not (1 <= h_eff <= t_post):
            raise ValueError(f"horizon must lie in [1, {t_post}]; got {horizon!r}")

        cohorts_ac, ever_ac, _ = _infer_panel_masks_from_result(
            results_map[pair_ac],
            _get_pair_arg(cohorts, pair_ac),
            _get_pair_arg(ever, pair_ac),
            _get_pair_arg(att_w, pair_ac),
        )
        own_traj_ac = _project_cohort_damped_trajectory(
            leg.post_att_ac,
            results_map[pair_ac],
            ever_ac,
            cohorts_ac,
            resolved_tau_clip,
            catt_floor,
        )
        gap_ac = float(np.max(np.abs(own_traj_ac - leg.trt_catt_ac)))
        add_gaps = {pair_ac: gap_ac}

        tau_inc = np.column_stack([leg.post_att_bc, leg.post_att_ab])
        tau_mu = np.column_stack([leg.post_att_bc, leg.post_att_ac])
        c_m1 = {
            pair_bc: leg.post_att_bc,
            pair_ab: leg.post_att_ab,
            pair_ac: leg.post_att_ac,
        }
        c_traj = {
            pair_bc: leg.trt_catt_bc,
            pair_ab: leg.trt_catt_ab,
            pair_ac: leg.trt_catt_ac,
        }
        c_cum = {
            pair: traj[:, :h_eff].sum(axis=1) for pair, traj in c_traj.items()
        }
        compute_rhat, _ = _get_rhat_ess()
        rhat_bc = float(compute_rhat(leg.post_att_bc.reshape(2, -1))) if n_draws >= 8 and n_draws % 2 == 0 else 1.0
        rhat_ac = float(compute_rhat(leg.post_att_ac.reshape(2, -1))) if n_draws >= 8 and n_draws % 2 == 0 else 1.0
        diag = SamplerDiagnostics(
            method="sir",
            ess=float(leg.sir_ess),
            r_hat=max(rhat_bc, rhat_ac),
            ess_per_param=(float(leg.sir_ess), float(leg.sir_ess)),
            r_hat_per_param=(rhat_bc, rhat_ac),
            n_proposals=int(n_proposals) if n_proposals is not None else max(n_draws * 500, 250_000),
            fallback_triggered=False,
        )
        pri_mu = (
            np.column_stack([leg.prior_att_bc, leg.prior_att_ac])
            if leg.prior_att_bc is not None and leg.prior_att_ac is not None
            else None
        )
        pri_inc = (
            np.column_stack([leg.prior_att_bc, leg.prior_att_ab])
            if leg.prior_att_bc is not None and leg.prior_att_ab is not None
            else None
        )
        pri_m1_map = (
            {
                pair_bc: leg.prior_att_bc,
                pair_ab: leg.prior_att_ab,
                pair_ac: leg.prior_att_ac,
            }
            if leg.prior_att_bc is not None
            and leg.prior_att_ab is not None
            and leg.prior_att_ac is not None
            else None
        )
        pri_traj_map = (
            {
                pair_bc: leg.prior_catt_bc,
                pair_ab: leg.prior_catt_ab,
                pair_ac: leg.prior_catt_ac,
            }
            if leg.prior_catt_bc is not None
            and leg.prior_catt_ab is not None
            and leg.prior_catt_ac is not None
            else None
        )
        cycle_traj_gap = float(
            np.max(np.abs(leg.trt_catt_ac - (leg.trt_catt_bc + leg.trt_catt_ab)))
        )
        return MultiArmReconciliationResult(
            arms=arms_tuple,
            topology=topology,
            tree_edges=tree_edges,
            fitted_pairs=fitted_pairs,
            redundant_pairs=redundant_pairs,
            tau_mu_draws=tau_mu,
            tau_increment_draws=tau_inc,
            contrast_mode1_draws=c_m1,
            contrast_trajectory_draws=c_traj,
            contrast_catt_draws=dict(c_traj),
            contrast_cumulative_draws=c_cum,
            mu_hat=mu_hat,
            sigma_lik=sigma_lik,
            arm_variances=arm_variances,
            nnls_residual=nnls_residual,
            is_bipartite=is_bipartite,
            cycle_discrepancy_ratio=cycle_discrepancy_ratio,
            trajectory_additivity_max_gap=cycle_traj_gap,
            sampler_diagnostics=diag,
            trajectory_additivity_gaps=add_gaps,
            prior_mode=prior_mode,
            prior_scale=leg.prior_scale,
            horizon=h_eff,
            own_vs_implied_trajectory_max_gap=gap_ac,
            prior_mu_draws=pri_mu,
            prior_increment_draws=pri_inc,
            prior_contrast_mode1_draws=pri_m1_map,
            prior_contrast_trajectory_draws=pri_traj_map,
            legacy_3arm_result=leg,
        )

    # Build spanning-tree incidence T_tree: theta = T_tree @ mu and mu = M_tree @ theta
    t_tree = np.zeros((d_dim, d_dim), dtype=np.float64)
    for e_idx, (t_arm, r_arm) in enumerate(tree_edges):
        j_idx = arm_to_idx[t_arm]
        l_idx = arm_to_idx[r_arm]
        if j_idx >= 1:
            t_tree[e_idx, j_idx - 1] = 1.0
        if l_idx >= 1:
            t_tree[e_idx, l_idx - 1] = -1.0
    m_tree = np.linalg.inv(t_tree)

    # Step 3: Prior construction on mu in R^{A - 1}
    rng = np.random.default_rng(int(seed))
    a_ineq, b_ineq = _build_support_inequalities(
        n_arms=n_arms,
        topology=topology,
        prior_mode=prior_mode,
        monotone=bool(monotone),
        tau_clip=resolved_tau_clip,
        t_tree=t_tree,
        cover_pairs=cover_pairs,
        top_idx=top_idx,
    )

    def _resolve_scalar_scale(
        scale_arg: Any, fallback_sd: float, label: str, eb: Mapping[str, Any] | None = None
    ) -> float:
        if isinstance(scale_arg, str):
            if scale_arg == "empirical_bayes":
                if eb is None:
                    raise ValueError(f"prior_scale='empirical_bayes' is not supported for {prior_mode!r}")
                floor = float(prior_scale_floor) if prior_scale_floor is not None else float(fallback_sd)
                return _empirical_bayes_scale(floor=floor, **eb)
            if scale_arg != "in_sample":
                raise ValueError(
                    "prior_scale must be a positive float, 'in_sample' or 'empirical_bayes'; "
                    f"got {scale_arg!r}"
                )
            _warn_rank_k(
                f"prior_scale='in_sample' sets the {label} to the in-sample posterior SD ({fallback_sd:.4g}), "
                "which shrinks effects that are several posterior SDs from zero toward zero. Report it "
                "next to an out-of-sample calibrated scale, or use prior_scale='empirical_bayes'."
            )
            return float(fallback_sd)
        if scale_arg is None:
            raise ValueError(
                f"reconcile_multi_arm_draws requires a positive prior_scale or 'in_sample' for {prior_mode!r}"
            )
        val = float(scale_arg)
        if not (math.isfinite(val) and val > 0.0):
            raise ValueError(f"prior_scale must be > 0; got {scale_arg!r}")
        return val

    is_flat_prior = prior_mode == "flat"
    scale_used: float | tuple[float, ...] | None = None

    if prior_mode in _BUNDLE_PRIOR_MODES:
        assert top_idx is not None
        top_arm = arms_tuple[top_idx]
        ctrl_arm = arms_tuple[0]
        top_key = frozenset((top_arm, ctrl_arm))
        if top_key in undir_to_oriented:
            _, fit_top, _ = undir_to_oriented[top_key]
            in_sample_sd_top = float(np.std(np.asarray(fit_top.att_draws, dtype=np.float64), ddof=1))
        else:
            in_sample_sd_top = float(
                math.sqrt(max(float(sigma_lik[top_idx - 1, top_idx - 1]), 1e-12))
            )
        sd_top = _resolve_scalar_scale(
            prior_scale,
            max(in_sample_sd_top, 1e-12),
            "bundle prior scale",
            eb={
                "t_hat": float(mu_hat[top_idx - 1]),
                "s_lik": math.sqrt(max(float(sigma_lik[top_idx - 1, top_idx - 1]), 1e-24)),
                "family": prior_mode,
                "spike_prob": float(spike_prob) if spike_mode else None,
            },
        )
        scale_used = sd_top

        if topology == "chain" and monotone:
            log_dir_norm = float(
                sp_special.gammaln(np.sum(alpha_vec)) - np.sum(sp_special.gammaln(alpha_vec))
            )

            def _prior_draw(n: int) -> np.ndarray:
                tot = _draw_bundle_prior(rng, prior_mode, sd_top, n)
                u_mat = rng.dirichlet(alpha_vec, size=n)
                theta_mat = u_mat * tot[:, None]
                return theta_mat @ m_tree.T

            def _prior_logpdf(pts: np.ndarray) -> np.ndarray:
                theta_mat = pts @ t_tree.T
                valid = np.all(theta_mat > 0.0, axis=1)
                out = np.full(pts.shape[0], -np.inf, dtype=np.float64)
                if not np.any(valid):
                    return out
                th_v = theta_mat[valid]
                tot_v = np.sum(th_v, axis=1)
                u_v = th_v / np.maximum(tot_v[:, None], 1e-300)
                dir_log = log_dir_norm + np.sum(
                    (alpha_vec[None, :] - 1.0) * np.log(np.maximum(u_v, 1e-300)),
                    axis=1,
                )
                out[valid] = (
                    _fast_bundle_prior_logpdf(prior_mode, sd_top, tot_v)
                    - float(n_arms - 2) * np.log(np.maximum(tot_v, 1e-300))
                    + dir_log
                )
                return out

        else:
            # Unique-top poset order-polytope prior (or chain with monotone=False)
            top_col = top_idx - 1
            interior_cols = [c for c in range(d_dim) if c != top_col]
            interior_arms = [c + 1 for c in interior_cols]
            m_int = len(interior_cols)
            if monotone and m_int > 0:
                reach_int = reachability[np.ix_(interior_arms, interior_arms)]
                lin_exts = _enumerate_linear_extensions(reach_int)
                n_exts = int(lin_exts.shape[0])
                log_inv_vol = float(sp_special.gammaln(m_int + 1) - math.log(max(n_exts, 1)))
            else:
                lin_exts = np.zeros((1, m_int), dtype=np.int32)
                n_exts = 1
                log_inv_vol = 0.0

            def _prior_draw(n: int) -> np.ndarray:
                tot = _draw_bundle_prior(rng, prior_mode, sd_top, n)
                pts = np.zeros((n, d_dim), dtype=np.float64)
                pts[:, top_col] = tot
                if m_int > 0:
                    if monotone:
                        ext_idx = rng.integers(0, n_exts, size=n)
                        chosen_exts = lin_exts[ext_idx]  # (n, m_int)
                        u_sorted = np.sort(rng.uniform(0.0, 1.0, size=(n, m_int)), axis=1)
                        w_mat = np.empty((n, m_int), dtype=np.float64)
                        rows = np.arange(n)[:, None]
                        w_mat[rows, chosen_exts] = u_sorted
                    else:
                        w_mat = rng.uniform(0.0, 1.0, size=(n, m_int))
                    pts[:, interior_cols] = tot[:, None] * w_mat
                return pts

            def _prior_logpdf(pts: np.ndarray) -> np.ndarray:
                tot = pts[:, top_col]
                valid = tot > 0.0
                if m_int > 0:
                    w_int = pts[:, interior_cols]
                    valid &= np.all((w_int >= 0.0) & (w_int <= tot[:, None]), axis=1)
                    if monotone:
                        for lo_arm, hi_arm in cover_pairs:
                            if lo_arm >= 1 and hi_arm >= 1:
                                valid &= pts[:, hi_arm - 1] >= pts[:, lo_arm - 1]
                out = np.full(pts.shape[0], -np.inf, dtype=np.float64)
                if not np.any(valid):
                    return out
                tot_v = tot[valid]
                out[valid] = (
                    _fast_bundle_prior_logpdf(prior_mode, sd_top, tot_v)
                    - float(m_int) * np.log(np.maximum(tot_v, 1e-300))
                    + log_inv_vol
                )
                return out

    elif prior_mode in _INDEPENDENT_PRIOR_MODES:
        use_increment_coords = prior_mode == "independent_halfnormal" and topology == "chain"
        if isinstance(prior_scale, Mapping):
            sds = []
            for a_idx in range(1, n_arms):
                arm_name = arms_tuple[a_idx]
                sds.append(float(prior_scale[arm_name]))
            s_pri = np.asarray(sds, dtype=np.float64)
        elif isinstance(prior_scale, Sequence) and not isinstance(prior_scale, str):
            s_pri = np.asarray(prior_scale, dtype=np.float64)
            if s_pri.shape != (d_dim,):
                raise ValueError(f"prior_scale sequence must have length {d_dim}; got {s_pri.shape}")
        else:
            if use_increment_coords:
                fallback_sds = []
                for e_edge in tree_edges:
                    _, fit_e, _ = undir_to_oriented[frozenset(e_edge)]
                    fallback_sds.append(
                        max(float(np.std(np.asarray(fit_e.att_draws, dtype=np.float64), ddof=1)), 1e-12)
                    )
            else:
                fallback_sds = [
                    max(math.sqrt(float(sigma_lik[j, j])), 1e-12) for j in range(d_dim)
                ]
            coord_mat = t_tree if use_increment_coords else np.eye(d_dim, dtype=np.float64)
            if isinstance(prior_mean, Mapping):
                mean_eb = [float(prior_mean.get(arms_tuple[a], 0.0)) for a in range(1, n_arms)]
            elif isinstance(prior_mean, Sequence) and not isinstance(prior_mean, str):
                mean_eb = [float(x) for x in prior_mean]
            else:
                mean_eb = [float(prior_mean)] * d_dim
            low_eb = (
                max(resolved_tau_clip[0], 0.0)
                if prior_mode == "independent_halfnormal"
                else resolved_tau_clip[0]
            )
            s_pri = np.array(
                [
                    _resolve_scalar_scale(
                        prior_scale,
                        fb,
                        f"prior SD of coord {j}",
                        eb={
                            "t_hat": float(coord_mat[j] @ mu_hat),
                            "s_lik": math.sqrt(max(float(coord_mat[j] @ sigma_lik @ coord_mat[j]), 1e-24)),
                            "family": "normal",
                            "low": float(low_eb),
                            "high": float(resolved_tau_clip[1]),
                            "mean": mean_eb[j],
                        },
                    )
                    for j, fb in enumerate(fallback_sds)
                ],
                dtype=np.float64,
            )
        if np.any(~np.isfinite(s_pri)) or np.any(s_pri <= 0.0):
            raise ValueError(f"prior_scale entries must be positive and finite; got {s_pri}")
        scale_used = tuple(float(v) for v in s_pri)

        if isinstance(prior_mean, Mapping):
            m_pri = np.array([float(prior_mean.get(arms_tuple[a], 0.0)) for a in range(1, n_arms)])
        elif isinstance(prior_mean, Sequence) and not isinstance(prior_mean, str):
            m_pri = np.asarray(prior_mean, dtype=np.float64)
        else:
            m_pri = np.full(d_dim, float(prior_mean), dtype=np.float64)

        low_c = max(resolved_tau_clip[0], 0.0) if prior_mode == "independent_halfnormal" else resolved_tau_clip[0]
        high_c = resolved_tau_clip[1]
        a_std_vec = np.array(
            [(low_c - m_pri[j]) / s_pri[j] if math.isfinite(low_c) else -np.inf for j in range(d_dim)],
            dtype=np.float64,
        )
        b_std_vec = np.array(
            [(high_c - m_pri[j]) / s_pri[j] if math.isfinite(high_c) else np.inf for j in range(d_dim)],
            dtype=np.float64,
        )
        trunc_mass_vec = np.maximum(
            sp_special.ndtr(b_std_vec) - sp_special.ndtr(a_std_vec), 1e-300
        )
        log_norm_total = float(
            -0.5 * d_dim * math.log(2.0 * math.pi)
            - np.sum(np.log(s_pri))
            - np.sum(np.log(trunc_mass_vec))
        )

        def _prior_draw(n: int) -> np.ndarray:
            coords = np.zeros((n, d_dim), dtype=np.float64)
            for j in range(d_dim):
                a_std = float(a_std_vec[j])
                b_std = float(b_std_vec[j])
                if math.isfinite(a_std) or math.isfinite(b_std):
                    coords[:, j] = sp_stats.truncnorm.rvs(
                        a_std, b_std, loc=m_pri[j], scale=s_pri[j], size=n, random_state=rng
                    )
                else:
                    coords[:, j] = m_pri[j] + s_pri[j] * rng.standard_normal(n)
            return coords @ m_tree.T if use_increment_coords else coords

        def _prior_logpdf(pts: np.ndarray) -> np.ndarray:
            coords = pts @ t_tree.T if use_increment_coords else pts
            valid = np.all((coords >= low_c) & (coords <= high_c), axis=1)
            valid &= _inside_polyhedron(pts, a_ineq, b_ineq)
            out = np.full(pts.shape[0], -np.inf, dtype=np.float64)
            if not np.any(valid):
                return out
            zsc = (coords[valid] - m_pri[None, :]) / s_pri[None, :]
            out[valid] = -0.5 * np.sum(zsc**2, axis=1) + log_norm_total
            return out

    else:
        # prior_mode == "flat"
        scale_used = None
        chol_lik = np.linalg.cholesky(sigma_lik + 1e-9 * np.eye(d_dim))

        def _prior_draw(n: int) -> np.ndarray:
            return mu_hat[None, :] + 2.0 * (rng.standard_normal((n, d_dim)) @ chol_lik.T)

        def _prior_logpdf(pts: np.ndarray) -> np.ndarray:
            valid = _inside_polyhedron(pts, a_ineq, b_ineq)
            return np.where(valid, 0.0, -np.inf)

    # Step 3b: resolve the unidentified arm-variance direction on a bipartite contrast graph
    lik_logpdf_override: Any = None
    null_interval: tuple[float, float] | None = None
    null_shift: float | None = None
    if is_bipartite:
        z_null, xi_lo, xi_hi = _bipartite_null_interval(
            np.asarray(arm_variances, dtype=np.float64),
            pair_indices,
            s_sq,
            colors,
            min_group_var_share,
        )
        null_interval = (xi_lo, xi_hi)
        v_alloc = np.asarray(arm_variances, dtype=np.float64).copy()
        null_shift = 0.0
        if bipartite_policy != "allocation" and xi_hi > xi_lo:
            xi_grid = np.linspace(xi_lo, xi_hi, 41 if bipartite_policy == "conservative" else _N_NULL_GRID)
            if bipartite_policy == "conservative":
                top_vec = np.zeros(d_dim, dtype=np.float64)
                top_vec[(top_idx if top_idx is not None else d_dim) - 1] = 1.0
                if is_flat_prior:
                    prior_prec = np.zeros((d_dim, d_dim), dtype=np.float64)
                else:
                    pri_probe = _prior_draw(4000)
                    pri_probe = pri_probe[np.isfinite(_prior_logpdf(pri_probe))]
                    prior_cov = np.atleast_2d(np.cov(pri_probe, rowvar=False))
                    prior_prec = np.linalg.pinv(prior_cov + 1e-12 * np.eye(d_dim))
                post_vars = []
                for xi_val in xi_grid:
                    sig_xi = _sigma_lik_from_arm_variances(v_alloc + xi_val * z_null)
                    post_cov = np.linalg.inv(np.linalg.inv(sig_xi) + prior_prec)
                    post_vars.append(float(top_vec @ post_cov @ top_vec))
                null_shift = float(xi_grid[int(np.argmax(post_vars))])
                arm_variances = _ArmVarianceArray(v_alloc + null_shift * z_null, arms_tuple)
                sigma_lik = _sigma_lik_from_arm_variances(np.asarray(arm_variances))
            else:
                sig_grid = [_sigma_lik_from_arm_variances(v_alloc + xi_val * z_null) for xi_val in xi_grid]
                comp_dists = [
                    sp_stats.multivariate_normal(mean=mu_hat, cov=sig_k, allow_singular=True)
                    for sig_k in sig_grid
                ]
                log_k = math.log(float(len(comp_dists)))

                def lik_logpdf_override(pts: np.ndarray) -> np.ndarray:
                    stack = np.vstack(
                        [np.asarray(dist.logpdf(pts), dtype=np.float64).reshape(-1) for dist in comp_dists]
                    )
                    return sp_special.logsumexp(stack, axis=0) - log_k

                null_shift = None
                sigma_lik = np.mean(np.stack(sig_grid), axis=0)
                arm_variances = _ArmVarianceArray(v_alloc + 0.5 * (xi_lo + xi_hi) * z_null, arms_tuple)

    # Step 4: Dimension-adaptive posterior sampling
    n_prop = int(n_proposals) if n_proposals is not None else max(n_draws * 400, 100_000)
    ess_gate = max(100.0, 0.02 * float(n_prop))
    use_sir = sampler == "sir" or (sampler == "auto" and d_dim <= 3)
    fallback_triggered = False

    log_z_slab: float | None = None
    if use_sir:
        post_mu, pri_mu_pool, sir_ess, log_z_slab = _sir_posterior_nd(
            rng,
            prior_draw=_prior_draw,
            prior_logpdf=_prior_logpdf,
            is_flat_prior=is_flat_prior,
            mu_lik=mu_hat,
            sigma_lik=sigma_lik,
            a_ineq=a_ineq,
            b_ineq=b_ineq,
            n_draws=n_draws,
            n_prop=n_prop,
            lik_logpdf=lik_logpdf_override,
            return_log_evidence=True,
        )
        if sampler == "auto" and (post_mu.shape[0] < n_draws or sir_ess < ess_gate):
            fallback_triggered = True
            use_sir = False
        elif sampler == "sir" and post_mu.shape[0] < n_draws:
            raise ValueError(
                "The prior and the likelihood put no mass inside tau_clip; widen tau_clip or prior_scale"
            )
        else:
            compute_rhat, _ = _get_rhat_ess()
            rhats_sir = []
            for j_idx in range(d_dim):
                col = post_mu[:, j_idx]
                if n_draws >= 8 and n_draws % 2 == 0:
                    rhats_sir.append(float(compute_rhat(col.reshape(2, -1))))
                else:
                    rhats_sir.append(1.0)
            diag = SamplerDiagnostics(
                method="sir",
                ess=float(sir_ess),
                r_hat=float(max(rhats_sir)),
                ess_per_param=tuple(float(sir_ess) for _ in range(d_dim)),
                r_hat_per_param=tuple(rhats_sir),
                n_proposals=n_prop,
                fallback_triggered=False,
            )

    if not use_sir:
        post_mu, pri_mu_pool, ess_per_param, rhat_per_param, n_evals, slice_info = _slice_posterior_nd(
            rng,
            prior_draw=_prior_draw,
            prior_logpdf=_prior_logpdf,
            mu_lik=mu_hat,
            sigma_lik=sigma_lik,
            a_ineq=a_ineq,
            b_ineq=b_ineq,
            n_draws=n_draws,
            num_chains=num_chains,
            num_burnin=num_burnin,
            thin=thin,
            lik_logpdf=lik_logpdf_override,
            auto_continue=bool(auto_continue),
            target_rhat=target_rhat,
            target_ess=target_ess,
            max_iter_multiplier=max_iter_multiplier,
            return_info=True,
        )
        diag = SamplerDiagnostics(
            method="slice",
            ess=float(min(ess_per_param)),
            r_hat=float(max(rhat_per_param)),
            ess_per_param=ess_per_param,
            r_hat_per_param=rhat_per_param,
            n_proposals=n_evals,
            fallback_triggered=fallback_triggered,
            n_iterations=int(slice_info["n_iterations"]),
            converged=slice_info["converged"],
        )

    spike_weight: float | None = None
    if spike_mode:
        if log_z_slab is None:
            _, _, _, log_z_slab = _sir_posterior_nd(
                rng,
                prior_draw=_prior_draw,
                prior_logpdf=_prior_logpdf,
                mu_lik=mu_hat,
                sigma_lik=sigma_lik,
                a_ineq=a_ineq,
                b_ineq=b_ineq,
                n_draws=1,
                n_prop=n_prop,
                lik_logpdf=lik_logpdf_override,
                return_log_evidence=True,
            )
        post_mu, spike_weight = _spike_slab_mix(
            rng, post_mu, float(log_z_slab), mu_hat, sigma_lik, float(spike_prob), lik_logpdf_override
        )

    tau_increment_draws = post_mu @ t_tree.T

    # Step 5: Project spanning-tree edge increments and synthesize all contrast trajectories
    tree_edge_trajs: list[np.ndarray] = []
    t_post_set: set[int] = set()
    for e_idx, edge_pair in enumerate(tree_edges):
        orig_pair, fit_e, _ = undir_to_oriented[frozenset(edge_pair)]
        sign_e = 1.0 if orig_pair == edge_pair else -1.0
        cohorts_e, ever_e, _ = _infer_panel_masks_from_result(
            fit_e,
            _get_pair_arg(cohorts, orig_pair),
            _get_pair_arg(ever, orig_pair),
            _get_pair_arg(att_w, orig_pair),
        )
        inc_draws = tau_increment_draws[:, e_idx]
        fit_e_oriented = (
            replace(
                fit_e,
                mode_tau_draws=tuple(
                    -np.asarray(m, dtype=np.float64) for m in fit_e.mode_tau_draws
                ),
            )
            if sign_e < 0.0
            else fit_e
        )
        traj_e = _project_cohort_damped_trajectory(
            inc_draws, fit_e_oriented, ever_e, cohorts_e, resolved_tau_clip, catt_floor
        )
        tree_edge_trajs.append(traj_e)
        t_post_set.add(int(traj_e.shape[1]))

    if len(t_post_set) != 1:
        raise ValueError(f"Spanning-tree fits have inconsistent T_post lengths: {sorted(t_post_set)}")
    (t_post,) = t_post_set
    h_eff = t_post if horizon is None else int(horizon)
    if not (1 <= h_eff <= t_post):
        raise ValueError(f"horizon must lie in [1, {t_post}]; got {horizon!r}")

    # Arm trajectories vs control 0: mu_traj[j] = sum_e M_tree[j - 1, e] * tree_edge_trajs[e]
    arm_trajs: list[np.ndarray] = [np.zeros((n_draws, t_post), dtype=np.float64)]
    for j_idx in range(1, n_arms):
        traj_j = np.zeros((n_draws, t_post), dtype=np.float64)
        for e_idx in range(d_dim):
            coeff = float(m_tree[j_idx - 1, e_idx])
            if coeff != 0.0:
                traj_j = traj_j + coeff * tree_edge_trajs[e_idx]
        arm_trajs.append(traj_j)

    arm_mode1: list[np.ndarray] = [np.zeros(n_draws, dtype=np.float64)] + [
        post_mu[:, j_idx] for j_idx in range(d_dim)
    ]

    # Populate all canonical pairs (arms[j], arms[l]) for j > l, plus any reverse fitted pairs
    canonical_pairs: list[tuple[str, str]] = []
    for l_idx in range(n_arms):
        for j_idx in range(l_idx + 1, n_arms):
            canonical_pairs.append((arms_tuple[j_idx], arms_tuple[l_idx]))
    for pair in fitted_pairs:
        if pair not in canonical_pairs:
            canonical_pairs.append(pair)

    contrast_mode1_draws: dict[tuple[str, str], np.ndarray] = {}
    contrast_trajectory_draws: dict[tuple[str, str], np.ndarray] = {}
    contrast_cumulative_draws: dict[tuple[str, str], np.ndarray] = {}
    for t_arm, r_arm in canonical_pairs:
        j_idx = arm_to_idx[t_arm]
        l_idx = arm_to_idx[r_arm]
        m1_draws = arm_mode1[j_idx] - arm_mode1[l_idx]
        tr_draws = arm_trajs[j_idx] - arm_trajs[l_idx]
        contrast_mode1_draws[(t_arm, r_arm)] = m1_draws
        contrast_trajectory_draws[(t_arm, r_arm)] = tr_draws
        contrast_cumulative_draws[(t_arm, r_arm)] = tr_draws[:, :h_eff].sum(axis=1)

    # Own vs implied trajectory additivity diagnostic on redundant fitted contrasts
    trajectory_additivity_gaps: dict[tuple[str, str], float] = {}
    for red_pair in redundant_pairs:
        fit_r = results_map[red_pair]
        cohorts_r, ever_r, _ = _infer_panel_masks_from_result(
            fit_r,
            _get_pair_arg(cohorts, red_pair),
            _get_pair_arg(ever, red_pair),
            _get_pair_arg(att_w, red_pair),
        )
        implied_m1 = contrast_mode1_draws[red_pair]
        own_traj = _project_cohort_damped_trajectory(
            implied_m1, fit_r, ever_r, cohorts_r, resolved_tau_clip, catt_floor
        )
        implied_traj = contrast_trajectory_draws[red_pair]
        min_t = min(own_traj.shape[1], implied_traj.shape[1])
        gap = float(np.max(np.abs(own_traj[:, :min_t] - implied_traj[:, :min_t])))
        trajectory_additivity_gaps[red_pair] = gap

    own_vs_implied_max_gap = float(
        max(trajectory_additivity_gaps.values(), default=0.0)
    )
    # Cycle additivity max gap across all triples (l < k < j) of reconciled trajectories
    cycle_traj_max_gap = 0.0
    for l_idx in range(n_arms):
        for k_idx in range(l_idx + 1, n_arms):
            for j_idx in range(k_idx + 1, n_arms):
                tr_jl = contrast_trajectory_draws[(arms_tuple[j_idx], arms_tuple[l_idx])]
                tr_jk = contrast_trajectory_draws[(arms_tuple[j_idx], arms_tuple[k_idx])]
                tr_kl = contrast_trajectory_draws[(arms_tuple[k_idx], arms_tuple[l_idx])]
                cycle_traj_max_gap = max(
                    cycle_traj_max_gap,
                    float(np.max(np.abs(tr_jl - (tr_jk + tr_kl)))),
                )

    pri_mu_out: np.ndarray | None = None
    pri_inc_out: np.ndarray | None = None
    pri_m1_out: dict[tuple[str, str], np.ndarray] | None = None
    pri_traj_out: dict[tuple[str, str], np.ndarray] | None = None
    if return_prior_draws and pri_mu_pool.shape[0] > 0:
        sel = rng.choice(
            pri_mu_pool.shape[0],
            size=n_draws,
            replace=pri_mu_pool.shape[0] < n_draws,
        )
        pri_mu_out = pri_mu_pool[sel]
        if spike_mode:
            pri_mu_out = np.array(pri_mu_out, dtype=np.float64, copy=True)
            pri_mu_out[rng.random(n_draws) < float(spike_prob)] = 0.0
        pri_inc_out = pri_mu_out @ t_tree.T
        pri_edge_trajs: list[np.ndarray] = []
        for e_idx, edge_pair in enumerate(tree_edges):
            orig_pair, fit_e, _ = undir_to_oriented[frozenset(edge_pair)]
            sign_e = 1.0 if orig_pair == edge_pair else -1.0
            cohorts_e, ever_e, _ = _infer_panel_masks_from_result(
                fit_e,
                _get_pair_arg(cohorts, orig_pair),
                _get_pair_arg(ever, orig_pair),
                _get_pair_arg(att_w, orig_pair),
            )
            fit_e_oriented = (
                replace(
                    fit_e,
                    mode_tau_draws=tuple(
                        -np.asarray(m, dtype=np.float64) for m in fit_e.mode_tau_draws
                    ),
                )
                if sign_e < 0.0
                else fit_e
            )
            p_traj_e = _project_cohort_damped_trajectory(
                pri_inc_out[:, e_idx],
                fit_e_oriented,
                ever_e,
                cohorts_e,
                resolved_tau_clip,
                catt_floor,
            )
            pri_edge_trajs.append(p_traj_e)
        pri_arm_trajs: list[np.ndarray] = [np.zeros((n_draws, t_post), dtype=np.float64)]
        for j_idx in range(1, n_arms):
            p_tr_j = np.zeros((n_draws, t_post), dtype=np.float64)
            for e_idx in range(d_dim):
                coeff = float(m_tree[j_idx - 1, e_idx])
                if coeff != 0.0:
                    p_tr_j = p_tr_j + coeff * pri_edge_trajs[e_idx]
            pri_arm_trajs.append(p_tr_j)
        pri_arm_m1 = [np.zeros(n_draws, dtype=np.float64)] + [
            pri_mu_out[:, j_idx] for j_idx in range(d_dim)
        ]
        pri_m1_out = {}
        pri_traj_out = {}
        for t_arm, r_arm in canonical_pairs:
            j_idx = arm_to_idx[t_arm]
            l_idx = arm_to_idx[r_arm]
            pri_m1_out[(t_arm, r_arm)] = pri_arm_m1[j_idx] - pri_arm_m1[l_idx]
            pri_traj_out[(t_arm, r_arm)] = pri_arm_trajs[j_idx] - pri_arm_trajs[l_idx]

    return MultiArmReconciliationResult(
        arms=arms_tuple,
        topology=topology_requested,
        tree_edges=tree_edges,
        fitted_pairs=fitted_pairs,
        redundant_pairs=redundant_pairs,
        tau_mu_draws=post_mu,
        tau_increment_draws=tau_increment_draws,
        contrast_mode1_draws=contrast_mode1_draws,
        contrast_trajectory_draws=contrast_trajectory_draws,
        contrast_catt_draws=dict(contrast_trajectory_draws),
        contrast_cumulative_draws=contrast_cumulative_draws,
        mu_hat=mu_hat,
        sigma_lik=sigma_lik,
        arm_variances=arm_variances,
        nnls_residual=nnls_residual,
        is_bipartite=is_bipartite,
        cycle_discrepancy_ratio=cycle_discrepancy_ratio,
        trajectory_additivity_max_gap=cycle_traj_max_gap,
        sampler_diagnostics=diag,
        trajectory_additivity_gaps=trajectory_additivity_gaps,
        prior_mode=prior_mode_requested,
        prior_scale=scale_used,
        horizon=h_eff,
        own_vs_implied_trajectory_max_gap=own_vs_implied_max_gap,
        prior_mu_draws=pri_mu_out,
        prior_increment_draws=pri_inc_out,
        prior_contrast_mode1_draws=pri_m1_out,
        prior_contrast_trajectory_draws=pri_traj_out,
        legacy_3arm_result=None,
        bipartite_policy=bipartite_policy if is_bipartite else None,
        null_shift_interval=null_interval,
        null_shift=null_shift,
        spike_posterior_prob=spike_weight,
        projection=projection,
        projection_weights=projection_weights,
        flagged_pairs=flagged_pairs,
    )
