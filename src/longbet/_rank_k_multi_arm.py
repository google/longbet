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

"""Joint hierarchical Bayesian reconciliation across multi-arm posets.

Implements the unified v0.11 Bayesian reconciliation engine shared by both
single-outcome (:func:`reconcile_multi_arm_draws`) and multi-outcome
(:func:`longbet._rank_k_multi.reconcile_multi_arm_multi_outcome_draws`) studies:

1. **Window & exposure basis**: Extracts treated-average trajectory draws on
   ``t = 1..H`` from Stage-1 ``catt_draws`` (never 15-month ``att_draws``) and
   projects onto the Euclidean-orthonormal basis ``Psi`` of dimension ``K'``
   whose first column is ``1_H / sqrt(H)``.
2. **Collapsed SUR + shared-arm + composite likelihood**: Marginalizes shared
   arm-level sampling noise ``xi_{m,a} ~ N(0, V_{m,a})`` (with cross-outcome
   correlation ``R``) analytically when updating ``Theta``, while accounting
   for fit-specific forest errors ``eta_{q} | omega_{q} ~ N(0, F_{q} / omega_{q})``
   and exact linear composite measurements.
3. **Poset-cone prior & hierarchical scales**: Enforces cover-relation
   monotonicity and per-period positivity in prior support (a polyhedral cone)
   with exact cone-normalized scale updates: an independent per-outcome
   generalized-horseshoe prior ``(sigma_m / s_{0,m})^2 ~ BetaPrime(1/4, 1/4)`` on
   the level scale and a Huang--Wand Half-t_3 prior on the dynamic-to-level
   ratio ``kappa_m / sigma_m`` (``P1``, shipped), or local per-cover Half-t_3
   scales (``P2``: not shipped; it failed simulation-based calibration, see the
   companion pre-registration section 5, and is reachable only through the
   private ``_model_options``).
4. **Exact Pakman--Paninski (2014) truncated-Gaussian HMC**: Wall-reflecting
   harmonic HMC in posterior-whitened coordinates across 4 chains with
   unconditional auto-continuation until ``split-Rhat <= 1.01`` and
   ``ESS >= 400``.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
import math
from typing import Any
import warnings

import numpy as np
from scipy import interpolate as sp_interpolate
from scipy import linalg as sp_linalg
from scipy import optimize as sp_optimize
from scipy import special as sp_special
from scipy import stats as sp_stats
from scipy.stats import qmc as sp_qmc

from longbet._rank_k import (
    RankKLongBetResult,
    _infer_panel_masks_from_result,
    orthonormal_exposure_basis,
)

__all__ = [
    "ContrastReconciliationResult",
    "MultiArmReconciliationResult",
    "RankKDesignWarning",
    "SamplerDiagnostics",
    "reconcile_multi_arm_draws",
]

_ALLOWED_SCALE_PRIORS = ("global_halft", "local_halft")
_ALLOWED_ROBUST_ERRORS = ("t4", "mixture", "gaussian")
_ALLOWED_MODEL_OPTION_KEYS = frozenset({
    "scale_prior",
    "robust_error",
    "period_positivity",
    "use_sur",
    "dynamic_shrinkage",
    "n_chains",
    "n_warmup",
    "n_draws_per_chain",
    "max_continuation_rounds",
})

# Precomputed Gauss-Hermite & Gauss-Legendre rules for deterministic integrals.
_LOG_2 = math.log(2.0)
_GH16_X, _GH16_W = np.polynomial.hermite.hermgauss(16)
_GH16_Z = _GH16_X * math.sqrt(2.0)
_GH16_PROB = _GH16_W / math.sqrt(math.pi)
_GL16_X, _GL16_W = np.polynomial.legendre.leggauss(16)
_GL8_X, _GL8_W = np.polynomial.legendre.leggauss(3)
_GL8_P0_Z = sp_special.ndtri(np.clip(0.75 + 0.25 * _GL8_X, 0.5, 1.0 - 1e-15))
_GL64_RAW_X, _GL64_RAW_W = np.polynomial.legendre.leggauss(64)
_GL64_PSI_W = 0.5 * math.pi * _GL64_RAW_W
_GL64_SIN = np.sin(0.25 * math.pi * (_GL64_RAW_X + 1.0))
_GL64_LOG_SIN = np.log(_GL64_SIN)

# Contamination-mixture log-normal inflation hyperparameters (R2: median 16x var = 4x SD).
_R2_BETA_ALPHA = 1.0
_R2_BETA_BETA = 9.0
_R2_LOG_INFL_MU = math.log(16.0)
_R2_LOG_INFL_SD = 1.0

# Hierarchical Stage-1/Stage-2 variance-decomposition Beta prior hyperparameters:
# rho_m ~ Beta(5.0, 1.5) (shared-arm variance scale), phi_m ~ Beta(2.0, 10.0) (forest variance share floor).
_RHO_BETA_ALPHA = 5.0
_RHO_BETA_BETA = 1.5
_PHI_BETA_ALPHA = 2.0
_PHI_BETA_BETA = 10.0

# Independent per-outcome standardized poset-cone scale prior hyperparameters (P1 / "global_halft"):
#   tau_m = sigma_m / s_{0,m},  (tau_m / tau_0)^2 ~ F(1/2, 1/2) = BetaPrime(1/4, 1/4),  tau_0 = 1
#   gamma_m = kappa_m / sigma_m ~ Half-t_3(0.50)
# i.e. the shrinkage factor 1 / (1 + tau_m^2) ~ Beta(1/4, 1/4): a generalized horseshoe
# (Armagan, Dunson & Clyde 2011) with a single smooth family, an integrable pole at tau_m = 0
# (null outcomes shrink to zero) and a heavy tau_m^{-3/2} tail (active outcomes stay unshrunk).
# s_{0,m} = sqrt(median_a(sigma^2_{pre,m,a})) is the pre-treatment outcome SD (falling back to
# 2 * SE_{top,m} / tau_0 on synthetic mock inputs without pre-period metadata). Outcomes are
# independent in prior scale; cross-outcome borrowing arises purely from SUR correlation R and
# composite measurements.
_SCALE_PRIOR_ALPHA_TAU = 0.5
_SCALE_PRIOR_NU_TAU = 0.5
_SCALE_PRIOR_TAU0 = 1.0
_SCALE_PRIOR_NU_GAMMA = 3.0
_SCALE_PRIOR_GAMMA0 = 0.50
_SCALE_PRIOR_C0 = 1.0
_SCALE_PRIOR_W_TAU_BOUND = 11.0
_SCALE_PRIOR_U_GAMMA_BOUND = 9.5
_SCALE_PRIOR_TREE_LOCAL = False
_SCALE_PRIOR_NU_LOCAL = 3.0

# Process-level cache of 2,049-point P1 diamond log-normalizer lookup tables.
_P1_LUT_CACHE: dict[tuple[Any, ...], list[float]] = {}


class RankKDesignWarning(UserWarning):
    """Emitted when a contrast graph is bipartite and admits no odd-cycle check."""


@dataclass(frozen=True)
class SamplerDiagnostics:
    """Convergence and mixing diagnostics of the 4-chain Bayesian reconciliation sampler."""

    sampler: str
    r_hat: float
    ess: float
    ess_tail: float
    n_chains: int
    n_draws_per_chain: int
    n_warmup: int
    n_continuation_rounds: int
    converged: bool
    max_rhat_by_block: dict[str, float] = field(default_factory=dict)
    min_ess_by_block: dict[str, float] = field(default_factory=dict)

    @property
    def max_rhat(self) -> float:
        """Alias for :attr:`r_hat` (worst-case split-Rhat across all parameter blocks)."""
        return self.r_hat

    @property
    def min_ess_bulk(self) -> float:
        """Alias for :attr:`ess` (minimum bulk ESS across all parameter blocks)."""
        return self.ess


@dataclass(frozen=True)
class ContrastReconciliationResult:
    """Posterior summary and draw arrays for a single arm contrast ``(arm_high, arm_low)``."""

    contrast: tuple[str, str]
    mode1_draws: np.ndarray
    trajectory_draws: np.ndarray
    cumulative_draws: np.ndarray
    mean: float
    sd: float
    ci_95: tuple[float, float]


@dataclass
class MultiArmReconciliationResult:
    """Posterior draws and diagnostics for a reconciled multi-arm study."""

    arms: tuple[str, ...]
    order: tuple[str, ...] | None
    cover_edges: tuple[tuple[str, str], ...]
    topology: str
    tree_edges: tuple[tuple[str, str], ...]
    fitted_pairs: tuple[tuple[str, str], ...]
    redundant_pairs: tuple[tuple[str, str], ...]
    all_pairs: tuple[tuple[str, str], ...]
    horizon: int
    n_modes: int
    basis: np.ndarray
    monotone: bool
    prior_scale: float | None
    s0_scale: float
    theta_draws: np.ndarray
    sigma_draws: np.ndarray
    kappa_draws: np.ndarray
    local_scale_draws: np.ndarray
    omega_draws: dict[tuple[str, str], np.ndarray]
    tau_mu_draws: np.ndarray
    tau_increment_draws: np.ndarray
    arm_trajectory_draws: dict[str, np.ndarray]
    contrast_mode1_draws: dict[tuple[str, str], np.ndarray]
    contrast_trajectory_draws: dict[tuple[str, str], np.ndarray]
    contrast_cumulative_draws: dict[tuple[str, str], np.ndarray]
    contrast_catt_draws: dict[tuple[str, str], np.ndarray]
    contrast_scale_inflation: dict[tuple[str, str], float]
    posterior_corruption_prob: dict[tuple[str, str], float]
    ppc_cycle_pvalue: float | None
    ppc_composite_pvalue: float | None
    ppc_basis_pvalue: float | None
    max_cycle_discrepancy: float
    sampler_diagnostics: SamplerDiagnostics
    rho_draws: np.ndarray | None = None
    phi_draws: np.ndarray | None = None
    tau_draws: np.ndarray | None = None

    @property
    def diagnostics(self) -> SamplerDiagnostics:
        """Alias for :attr:`sampler_diagnostics`."""
        return self.sampler_diagnostics

    @property
    def trt_catt_bc(self) -> np.ndarray | None:
        """3-arm convenience accessor for ``("B", "C")`` treated trajectory draws."""
        return self.contrast_trajectory_draws.get(("B", "C"))

    @property
    def trt_catt_ab(self) -> np.ndarray | None:
        """3-arm convenience accessor for ``("A", "B")`` treated trajectory draws."""
        return self.contrast_trajectory_draws.get(("A", "B"))

    @property
    def trt_catt_ac(self) -> np.ndarray | None:
        """3-arm convenience accessor for ``("A", "C")`` treated trajectory draws."""
        return self.contrast_trajectory_draws.get(("A", "C"))

    def get_contrast_result(self, arm_high: str, arm_low: str) -> ContrastReconciliationResult:
        """Return posterior draws and summary for any distinct ordered pair ``(arm_high, arm_low)``."""
        u, v = str(arm_high), str(arm_low)
        if u not in self.arms or v not in self.arms:
            raise ValueError(f"Contrast ({u!r}, {v!r}) references arm outside {self.arms}")
        if u == v:
            raise ValueError(f"Contrast requires two distinct arms; got ({u!r}, {v!r})")
        if (u, v) in self.contrast_trajectory_draws:
            m1 = self.contrast_mode1_draws[(u, v)]
            traj = self.contrast_trajectory_draws[(u, v)]
            cum = self.contrast_cumulative_draws[(u, v)]
        else:
            m1 = -self.contrast_mode1_draws[(v, u)]
            traj = -self.contrast_trajectory_draws[(v, u)]
            cum = -self.contrast_cumulative_draws[(v, u)]
        q_lo, q_hi = np.quantile(m1, [0.025, 0.975])
        return ContrastReconciliationResult(
            contrast=(u, v),
            mode1_draws=m1,
            trajectory_draws=traj,
            cumulative_draws=cum,
            mean=float(np.mean(m1)),
            sd=float(np.std(m1, ddof=1)),
            ci_95=(float(q_lo), float(q_hi)),
        )


# ---- Linear algebra, diagnostics, and validation helpers -----------------------------------------


def _psd_clip(mat: np.ndarray, floor: float = 1e-12) -> np.ndarray:
    """Symmetrize ``mat`` and clip eigenvalues below ``floor``."""
    arr = np.atleast_2d(np.asarray(mat, dtype=np.float64))
    sym = 0.5 * (arr + arr.T)
    eigvals, eigvecs = np.linalg.eigh(sym)
    eigvals = np.maximum(eigvals, float(floor))
    out = (eigvecs * eigvals[None, :]) @ eigvecs.T
    return 0.5 * (out + out.T)


def _psd_sqrt(mat: np.ndarray) -> np.ndarray:
    """Symmetric positive-semidefinite square root ``V^{1/2}``."""
    sym = _psd_clip(mat, floor=0.0)
    eigvals, eigvecs = np.linalg.eigh(sym)
    sqrt_vals = np.sqrt(np.maximum(eigvals, 0.0))
    out = (eigvecs * sqrt_vals[None, :]) @ eigvecs.T
    return 0.5 * (out + out.T)


def _euclidean_exposure_basis(horizon: int, n_modes: int) -> np.ndarray:
    """Euclidean-orthonormal exposure basis ``Psi`` of shape ``(H, K')`` with ``Psi[:, 0] = 1/sqrt(H)``."""
    h_int = int(horizon)
    k_eff = min(int(n_modes), h_int)
    raw = orthonormal_exposure_basis(h_int, k_eff)
    psi = np.asarray(raw[:, :k_eff], dtype=np.float64) / math.sqrt(float(h_int))
    return psi


def _validate_prior_scale_scalar(value: Any, name: str = "prior_scale") -> float | None:
    """Validate ``prior_scale``: only ``None`` or a positive finite float is permitted."""
    if value is None:
        return None
    if isinstance(value, (str, bytes, bool, np.bool_)):
        raise ValueError(
            f"{name} must be None (hierarchical generalized-horseshoe default) or a positive float; "
            f"string/boolean modes such as {value!r} were removed in v0.11"
        )
    try:
        val = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be None or a positive float; got {value!r}") from exc
    if not (math.isfinite(val) and val > 0.0):
        raise ValueError(f"{name} must be finite and strictly positive; got {value!r}")
    return val


def _parse_model_options(options: Mapping[str, Any] | None) -> dict[str, Any]:
    """Validate private ``_model_options`` dictionary for candidate comparisons and ablations."""
    cfg: dict[str, Any] = {
        "scale_prior": "global_halft",
        "robust_error": "t4",
        "period_positivity": True,
        "use_sur": True,
        "dynamic_shrinkage": True,
        "n_chains": 4,
        "n_warmup": 200,
        "n_draws_per_chain": 400,
        "max_continuation_rounds": 4,
    }
    if options is None:
        return cfg
    if not isinstance(options, Mapping):
        raise ValueError(f"_model_options must be a mapping or None; got {type(options).__name__}")
    unknown = set(options.keys()) - _ALLOWED_MODEL_OPTION_KEYS
    if unknown:
        raise ValueError(f"Unsupported _model_options keys: {sorted(unknown)}")
    for k, v in options.items():
        cfg[k] = v
    if cfg["scale_prior"] not in _ALLOWED_SCALE_PRIORS:
        raise ValueError(f"scale_prior must be one of {_ALLOWED_SCALE_PRIORS}; got {cfg['scale_prior']!r}")
    if cfg["robust_error"] not in _ALLOWED_ROBUST_ERRORS:
        raise ValueError(f"robust_error must be one of {_ALLOWED_ROBUST_ERRORS}; got {cfg['robust_error']!r}")
    cfg["period_positivity"] = bool(cfg["period_positivity"])
    cfg["use_sur"] = bool(cfg["use_sur"])
    cfg["dynamic_shrinkage"] = bool(cfg["dynamic_shrinkage"])
    cfg["n_chains"] = max(int(cfg["n_chains"]), 2)
    cfg["n_warmup"] = max(int(cfg["n_warmup"]), 50)
    cfg["n_draws_per_chain"] = max(int(cfg["n_draws_per_chain"]), 50)
    cfg["max_continuation_rounds"] = max(int(cfg["max_continuation_rounds"]), 0)
    return cfg


def _split_rhat_and_ess_1d(chains_2d: np.ndarray) -> tuple[float, float, float]:
    """Compute Vehtari et al. (2021) split-Rhat, bulk ESS, and tail ESS for ``(C, N)`` chains."""
    arr = np.asarray(chains_2d, dtype=np.float64)
    n_c, n_s = arr.shape
    half = n_s // 2
    if half < 4:
        return 1.0, float(n_c * n_s), float(n_c * n_s)
    split = np.concatenate([arr[:, :half], arr[:, 2 * half - half : 2 * half]], axis=0)
    if float(np.std(split)) < 1e-14:
        return 1.0, float(n_c * n_s), float(n_c * n_s)

    def _rhat_ess_raw(x: np.ndarray) -> tuple[float, float]:
        m_chains, n_len = x.shape
        chain_means = np.mean(x, axis=1)
        chain_vars = np.var(x, axis=1, ddof=1)
        w_var = float(np.mean(chain_vars))
        b_var = float(n_len * np.var(chain_means, ddof=1)) if m_chains > 1 else 0.0
        if w_var <= 1e-30:
            return 1.0, float(m_chains * n_len)
        var_hat = ((n_len - 1.0) / n_len) * w_var + (1.0 / n_len) * b_var
        r_hat_val = math.sqrt(max(var_hat / w_var, 1.0))

        # FFT autocorrelation pooled across split chains (Geyer initial positive sequence)
        centered = x - chain_means[:, None]
        n_fft = 1 << (2 * n_len - 1).bit_length()
        spec = np.fft.rfft(centered, n=n_fft, axis=1)
        acov = np.fft.irfft(spec * np.conj(spec), n=n_fft, axis=1)[:, :n_len]
        acov /= np.arange(n_len, 0, -1, dtype=np.float64)[None, :]
        mean_acov = np.mean(acov, axis=0)
        rho_t = 1.0 - (w_var - mean_acov) / max(var_hat, 1e-30)
        rho_t[0] = 1.0
        tau = -1.0
        max_pair = n_len // 2
        for t_idx in range(max_pair):
            pair_sum = float(rho_t[2 * t_idx] + (rho_t[2 * t_idx + 1] if 2 * t_idx + 1 < n_len else 0.0))
            if pair_sum <= 0.0:
                break
            tau += 2.0 * pair_sum
        tau = max(tau, 1.0 / float(m_chains * n_len))
        ess_val = float(m_chains * n_len) / tau
        return r_hat_val, min(ess_val, float(m_chains * n_len))

    # Rank-normalize for bulk Rhat/ESS
    flat = split.reshape(-1)
    ranks = sp_stats.rankdata(flat, method="average")
    z_bulk = sp_special.ndtri((ranks - 0.375) / ( float(flat.size) + 0.25 )).reshape(split.shape)
    rhat_bulk, ess_bulk = _rhat_ess_raw(z_bulk)

    # Folded rank-normalize for tail Rhat/ESS
    folded = np.abs(split - np.median(split))
    f_ranks = sp_stats.rankdata(folded.reshape(-1), method="average")
    z_tail = sp_special.ndtri((f_ranks - 0.375) / (float(flat.size) + 0.25)).reshape(split.shape)
    rhat_tail, ess_tail = _rhat_ess_raw(z_tail)

    return max(rhat_bulk, rhat_tail), ess_bulk, ess_tail


# ---- Topology, poset transitive reduction, and bipartite graph check -----------------------------


def _transitive_reduction_dag(
    nodes: Sequence[str],
    directed_edges: Sequence[tuple[str, str]],
) -> tuple[tuple[str, str], ...]:
    """Validate acyclicity on ``(high, low)`` relations and return the unique cover edges."""
    node_set = set(nodes)
    adj: dict[str, set[str]] = {u: set() for u in nodes}
    dedup_edges: list[tuple[str, str]] = []
    seen_e: set[tuple[str, str]] = set()
    for raw_hi, raw_lo in directed_edges:
        hi, lo = str(raw_hi), str(raw_lo)
        if hi == lo:
            raise ValueError(f"Poset relation cannot be reflexive: ({hi!r}, {lo!r})")
        if hi not in node_set or lo not in node_set:
            raise ValueError(f"Relation ({hi!r}, {lo!r}) references arm outside {tuple(nodes)}")
        if (hi, lo) not in seen_e:
            seen_e.add((hi, lo))
            dedup_edges.append((hi, lo))
            adj[hi].add(lo)

    # Cycle check via 3-color DFS on high -> low
    state: dict[str, int] = {u: 0 for u in nodes}

    def _dfs_cycle(u: str) -> None:
        state[u] = 1
        for v in adj[u]:
            if state[v] == 1:
                raise ValueError(f"Poset relations contain a directed cycle involving {u!r} and {v!r}")
            if state[v] == 0:
                _dfs_cycle(v)
        state[u] = 2

    for u in nodes:
        if state[u] == 0:
            _dfs_cycle(u)

    # Reachability in >= 1 steps from each node
    reach: dict[str, set[str]] = {}
    for u in nodes:
        vis: set[str] = set()
        q: deque[str] = deque(adj[u])
        while q:
            cur = q.popleft()
            if cur not in vis:
                vis.add(cur)
                q.extend(adj[cur] - vis)
        reach[u] = vis

    # (hi, lo) is a cover edge iff lo is not reachable from any intermediate w in adj[hi]
    covers: list[tuple[str, str]] = []
    for hi, lo in dedup_edges:
        if not any(lo in reach[w] for w in adj[hi] if w != lo):
            covers.append((hi, lo))
    return tuple(covers)


def _topological_sort_low_to_high(
    nodes: Sequence[str],
    cover_edges: Sequence[tuple[str, str]],
    control_arm: str,
) -> tuple[str, ...]:
    """Topological order ``(control, ..., top)`` where ``low`` precedes ``high`` for every ``(high, low)``."""
    order_hint = {a: idx for idx, a in enumerate(nodes)}
    in_deg = {a: 0 for a in nodes}
    up_adj: dict[str, list[str]] = {a: [] for a in nodes}
    for hi, lo in cover_edges:
        up_adj[lo].append(hi)
        in_deg[hi] += 1
    if in_deg[control_arm] > 0:
        raise ValueError(f"Control arm {control_arm!r} cannot be strictly above another arm in poset")

    ready = sorted([a for a in nodes if in_deg[a] == 0], key=lambda a: (0 if a == control_arm else 1, order_hint[a]))
    result: list[str] = []
    while ready:
        cur = ready.pop(0)
        result.append(cur)
        for nxt in sorted(up_adj[cur], key=lambda a: order_hint[a]):
            in_deg[nxt] -= 1
            if in_deg[nxt] == 0:
                ready.append(nxt)
                ready.sort(key=lambda a: (0 if a == control_arm else 1, order_hint[a]))
    if len(result) != len(nodes):
        raise ValueError("Poset relations contain a cycle")
    if result[0] != control_arm:
        result.remove(control_arm)
        result.insert(0, control_arm)
    return tuple(result)


def _resolve_arms_and_poset(
    pair_keys: Sequence[tuple[str, str]],
    arms: Sequence[str] | None,
    order: Sequence[str] | None,
    poset: Sequence[tuple[str, str]] | None,
) -> dict[str, Any]:
    """Resolve arm list, cover edges, spanning tree, cut matrix, cycle basis, and bipartiteness."""
    if not pair_keys:
        raise ValueError("results must contain at least one fitted contrast pair")
    if order is not None and poset is not None:
        raise ValueError("Specify at most one of `order` or `poset`, not both")

    clean_pairs: list[tuple[str, str]] = []
    seen_pairs: set[tuple[str, str]] = set()
    discovered_arms: list[str] = []
    for item in pair_keys:
        if not (isinstance(item, tuple) and len(item) == 2):
            raise ValueError(f"Contrast keys must be 2-tuples (arm_high, arm_low); got {item!r}")
        u, v = str(item[0]), str(item[1])
        if not u or not v or u == v:
            raise ValueError(f"Contrast pair must name two distinct non-empty arms; got {item!r}")
        if (u, v) in seen_pairs or (v, u) in seen_pairs:
            raise ValueError(f"Duplicate or reversed contrast pair ({u!r}, {v!r}) in results")
        seen_pairs.add((u, v))
        clean_pairs.append((u, v))
        for a in (u, v):
            if a not in discovered_arms:
                discovered_arms.append(a)

    if arms is not None:
        arms_list = [str(a) for a in arms]
        if len(arms_list) < 2 or len(set(arms_list)) != len(arms_list):
            raise ValueError(f"`arms` must contain >= 2 unique arm names; got {arms!r}")
        unknown = set(discovered_arms) - set(arms_list)
        if unknown:
            raise ValueError(f"Fitted pairs reference arms {sorted(unknown)} not present in `arms`")
    else:
        arms_list = None

    resolved_order: tuple[str, ...] | None = None
    if order is not None:
        ord_list = [str(a) for a in order]
        if len(ord_list) < 2 or len(set(ord_list)) != len(ord_list):
            raise ValueError(f"`order` must contain >= 2 unique arm names; got {order!r}")
        if arms_list is not None and set(arms_list) != set(ord_list):
            raise ValueError(f"`arms` {arms_list} and `order` {ord_list} must contain the same arms")
        unknown = set(discovered_arms) - set(ord_list)
        if unknown:
            raise ValueError(f"Fitted pairs reference arms {sorted(unknown)} not in `order`")
        arms_tuple = tuple(ord_list)
        resolved_order = arms_tuple
        cover_edges = tuple((arms_tuple[i], arms_tuple[i - 1]) for i in range(1, len(arms_tuple)))
    elif poset is not None:
        poset_list = [(str(hi), str(lo)) for hi, lo in poset]
        if not poset_list:
            raise ValueError("`poset` must contain at least one (high, low) cover relation")
        poset_nodes = list(arms_list) if arms_list is not None else list(discovered_arms)
        for hi, lo in poset_list:
            for a in (hi, lo):
                if a not in poset_nodes:
                    if arms_list is not None:
                        raise ValueError(f"`poset` references arm {a!r} not in `arms`")
                    poset_nodes.append(a)
        cover_edges = _transitive_reduction_dag(poset_nodes, poset_list)
        if arms_list is not None:
            control_arm = arms_list[0]
            # Preserve caller's `arms` order if already topologically valid; otherwise sort cleanly.
            pos_of = {a: idx for idx, a in enumerate(arms_list)}
            if any(hi == control_arm for hi, _ in cover_edges):
                raise ValueError(f"Control arm {control_arm!r} (arms[0]) cannot dominate another arm in `poset`")
            if all(pos_of[lo] < pos_of[hi] for hi, lo in cover_edges):
                arms_tuple = tuple(arms_list)
            else:
                arms_tuple = _topological_sort_low_to_high(arms_list, cover_edges, control_arm)
        else:
            highs = {hi for hi, _ in cover_edges}
            minimals = [a for a in poset_nodes if a not in highs]
            if len(minimals) != 1:
                raise ValueError(
                    f"Cannot infer a unique control arm from `poset` (minimal arms: {minimals}); pass `arms=` explicitly"
                )
            arms_tuple = _topological_sort_low_to_high(poset_nodes, cover_edges, minimals[0])
    else:
        # Infer DAG directly from directed keys (arm_high, arm_low) in results
        base_nodes = list(arms_list) if arms_list is not None else list(discovered_arms)
        cover_edges = _transitive_reduction_dag(base_nodes, clean_pairs)
        if arms_list is not None:
            control_arm = arms_list[0]
            pos_of = {a: idx for idx, a in enumerate(arms_list)}
            if any(hi == control_arm for hi, _ in cover_edges):
                raise ValueError(
                    f"Inferred directed graph has {control_arm!r} (arms[0]) as treated; pass `order=` or `poset=`"
                )
            if all(pos_of[lo] < pos_of[hi] for hi, lo in cover_edges):
                arms_tuple = tuple(arms_list)
            else:
                arms_tuple = _topological_sort_low_to_high(arms_list, cover_edges, control_arm)
        else:
            highs = {hi for hi, _ in cover_edges}
            minimals = [a for a in base_nodes if a not in highs]
            if len(minimals) != 1:
                raise ValueError(
                    f"Cannot infer a unique control arm from pair keys (minimal arms: {minimals}); "
                    "pass `order=` or `poset=` explicitly"
                )
            arms_tuple = _topological_sort_low_to_high(base_nodes, cover_edges, minimals[0])

    n_arms = len(arms_tuple)
    arm_to_idx = {a: idx for idx, a in enumerate(arms_tuple)}

    # Cut matrix B_cut of shape (P, A - 1) on non-control arms 1..A-1
    n_pairs = len(clean_pairs)
    b_cut = np.zeros((n_pairs, n_arms - 1), dtype=np.float64)
    for p_idx, (u, v) in enumerate(clean_pairs):
        iu, iv = arm_to_idx[u], arm_to_idx[v]
        if iu > 0:
            b_cut[p_idx, iu - 1] = 1.0
        if iv > 0:
            b_cut[p_idx, iv - 1] = -1.0
    if np.linalg.matrix_rank(b_cut) < n_arms - 1:
        raise ValueError(
            f"Fitted contrast pairs {clean_pairs} do not form a connected graph across arms {arms_tuple}"
        )

    chain_covers = tuple((arms_tuple[i], arms_tuple[i - 1]) for i in range(1, n_arms))
    topology = "chain" if set(cover_edges) == set(chain_covers) and len(cover_edges) == n_arms - 1 else "poset"
    if topology == "chain" and resolved_order is None:
        resolved_order = arms_tuple

    # Spanning tree rooted at arms_tuple[0] (preferring cover edges, then fitted pairs)
    undir_adj: dict[str, list[str]] = {a: [] for a in arms_tuple}
    for hi, lo in (*cover_edges, *clean_pairs):
        if hi not in undir_adj[lo]:
            undir_adj[lo].append(hi)
        if lo not in undir_adj[hi]:
            undir_adj[hi].append(lo)
    parent_of: dict[str, str] = {}
    bfs_q: deque[str] = deque([arms_tuple[0]])
    visited_tree = {arms_tuple[0]}
    while bfs_q:
        cur = bfs_q.popleft()
        for nb in undir_adj[cur]:
            if nb not in visited_tree:
                visited_tree.add(nb)
                parent_of[nb] = cur
                bfs_q.append(nb)
    tree_edges = tuple((a, parent_of[a]) for a in arms_tuple[1:])

    # Partition fitted_pairs into a spanning forest vs redundant cycle-closing pairs
    uf_parent = list(range(n_arms))

    def _find(x: int) -> int:
        while uf_parent[x] != x:
            uf_parent[x] = uf_parent[uf_parent[x]]
            x = uf_parent[x]
        return x

    redundant_pairs: list[tuple[str, str]] = []
    for u, v in clean_pairs:
        ru, rv = _find(arm_to_idx[u]), _find(arm_to_idx[v])
        if ru == rv:
            redundant_pairs.append((u, v))
        else:
            uf_parent[ru] = rv

    # Left nullspace N_cycle of B_cut (N_cycle^T B_cut == 0)
    u_svd, s_svd, _ = np.linalg.svd(b_cut, full_matrices=True)
    rank_b = int(np.sum(s_svd > 1e-10))
    n_cycle = u_svd[:, rank_b:]  # (P, P - (A - 1))

    # Bipartiteness check on undirected fitted contrast graph
    f_adj: dict[int, list[int]] = {i: [] for i in range(n_arms)}
    fitted_undir: set[tuple[int, int]] = set()
    for u, v in clean_pairs:
        iu, iv = arm_to_idx[u], arm_to_idx[v]
        f_adj[iu].append(iv)
        f_adj[iv].append(iu)
        fitted_undir.add((max(iu, iv), min(iu, iv)))
    color = [-1] * n_arms
    color[0] = 0
    is_bipartite = True
    q_bip: deque[int] = deque([0])
    while q_bip and is_bipartite:
        cur = q_bip.popleft()
        for nb in f_adj[cur]:
            if color[nb] == -1:
                color[nb] = 1 - color[cur]
                q_bip.append(nb)
            elif color[nb] == color[cur]:
                is_bipartite = False
                break

    suggested_odd_pair: tuple[str, str] | None = None
    if is_bipartite and n_arms >= 3:
        # Pick the highest-indexed missing pair within the same color class to close an odd cycle
        candidates: list[tuple[int, int]] = []
        for j in range(n_arms - 1, -1, -1):
            for l in range(j - 1, -1, -1):
                if color[j] == color[l] and (j, l) not in fitted_undir:
                    # Prefer closing a 3-cycle with two already-fitted edges
                    closes_triangle = any(
                        (max(j, w), min(j, w)) in fitted_undir and (max(l, w), min(l, w)) in fitted_undir
                        for w in range(n_arms)
                        if w not in (j, l)
                    )
                    candidates.append((1 if closes_triangle else 0, j, l))
        if candidates:
            candidates.sort(reverse=True)
            _, hi_i, lo_i = candidates[0]
            suggested_odd_pair = (arms_tuple[hi_i], arms_tuple[lo_i])

    all_pairs = tuple(
        (arms_tuple[j], arms_tuple[l])
        for j in range(1, n_arms)
        for l in range(j)
    )

    # Cover matrix C_cover of shape (E_c, A - 1)
    c_cover = np.zeros((len(cover_edges), n_arms - 1), dtype=np.float64)
    for e_idx, (hi, lo) in enumerate(cover_edges):
        ihi, ilo = arm_to_idx[hi], arm_to_idx[lo]
        if ihi > 0:
            c_cover[e_idx, ihi - 1] = 1.0
        if ilo > 0:
            c_cover[e_idx, ilo - 1] = -1.0

    # Strict interior heights h_a >= 1 for cone feasibility
    heights = np.zeros(n_arms, dtype=np.float64)
    for a in arms_tuple[1:]:
        ia = arm_to_idx[a]
        lows_of_a = [arm_to_idx[lo] for hi, lo in cover_edges if hi == a]
        heights[ia] = (max(heights[lo_i] for lo_i in lows_of_a) + 1.0) if lows_of_a else 1.0

    return {
        "arms": arms_tuple,
        "arm_to_idx": arm_to_idx,
        "order": resolved_order,
        "cover_edges": cover_edges,
        "c_cover": c_cover,
        "topology": topology,
        "tree_edges": tree_edges,
        "fitted_pairs": tuple(clean_pairs),
        "redundant_pairs": tuple(redundant_pairs),
        "all_pairs": all_pairs,
        "b_cut": b_cut,
        "n_cycle": n_cycle,
        "is_bipartite": is_bipartite,
        "suggested_odd_pair": suggested_odd_pair,
        "interior_heights": heights[1:],  # (A - 1,)
    }


# ---- Stage-1 trajectory & arm residual covariance extraction ------------------------------------


def _extract_raw_trt_trajectory(res: RankKLongBetResult) -> np.ndarray:
    """Extract ``(S, T_post)`` treated-average post-adoption trajectory draws from ``res``."""
    meta = getattr(res, "arm_metadata", None)
    if isinstance(meta, Mapping):
        for k_traj in (
            "trt_trajectory_draws",
            "treated_trajectory_draws",
            "catt_trt_draws",
            "trt_post_traj_draws",
        ):
            if k_traj in meta and meta[k_traj] is not None:
                traj = np.asarray(meta[k_traj], dtype=np.float64)
                if (
                    traj.ndim == 2
                    and traj.shape[0] >= 2
                    and traj.shape[1] >= 1
                    and np.all(np.isfinite(traj))
                ):
                    return traj

    catt = np.asarray(res.catt_draws, dtype=np.float64)
    if catt.ndim != 3 or catt.shape[0] < 2:
        raise ValueError(f"RankKLongBetResult.catt_draws must be 3-D (S, N, T) with S >= 2; got {catt.shape}")
    cohorts_arr, ever_arr, _ = _infer_panel_masks_from_result(res, None, None, None)
    if not np.any(ever_arr):
        raise ValueError("RankKLongBetResult has no ever-treated units to extract post-adoption trajectory")
    n_draws, _, n_periods = catt.shape
    waves = sorted(int(g) for g in np.unique(cohorts_arr[ever_arr]))
    t0 = waves[0]
    t_post = max(n_periods - t0, 1)
    if len(waves) <= 1:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            traj = np.nanmean(catt[:, ever_arr, t0 : t0 + t_post], axis=1)
        return np.asarray(traj, dtype=np.float64)

    traj = np.zeros((n_draws, t_post), dtype=np.float64)
    trt_idx = np.nonzero(ever_arr)[0]
    trt_c = cohorts_arr[trt_idx]
    for s_step in range(t_post):
        cal_t = trt_c + s_step
        valid = cal_t < n_periods
        if np.any(valid):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                traj[:, s_step] = np.nanmean(catt[:, trt_idx[valid], cal_t[valid]], axis=1)
        elif s_step > 0:
            traj[:, s_step] = traj[:, s_step - 1]
    return traj


def _lookup_meta_cov(meta: Mapping[str, Any] | None, role: str) -> np.ndarray | None:
    """Read treated or control post-window residual covariance from ``arm_metadata``."""
    if not isinstance(meta, Mapping):
        return None
    if role == "treated":
        keys = ("post_residual_cov_treated", "treated_post_cov", "post_cov_treated")
    else:
        keys = ("post_residual_cov_control", "control_post_cov", "post_cov_control")
    for k in keys:
        if k in meta and meta[k] is not None:
            arr = np.asarray(meta[k], dtype=np.float64)
            if arr.ndim == 2 and arr.shape[0] == arr.shape[1] and arr.shape[0] >= 1 and np.all(np.isfinite(arr)):
                return arr
    nested = meta.get(role)
    if isinstance(nested, Mapping) and nested.get("post_cov") is not None:
        arr = np.asarray(nested["post_cov"], dtype=np.float64)
        if arr.ndim == 2 and arr.shape[0] == arr.shape[1] and arr.shape[0] >= 1 and np.all(np.isfinite(arr)):
            return arr
    return None


def _lookup_meta_pre_var(meta: Mapping[str, Any] | None, role: str) -> float | None:
    """Read treated or control pre-period residual variance from ``arm_metadata``."""
    if not isinstance(meta, Mapping):
        return None
    if role == "treated":
        keys = ("pre_residual_var_treated", "treated_pre_var", "pre_var_treated")
    else:
        keys = ("pre_residual_var_control", "control_pre_var", "pre_var_control")
    for k in keys:
        if k in meta and meta[k] is not None:
            try:
                val = float(meta[k])
            except (TypeError, ValueError):
                val = float("nan")
            if math.isfinite(val) and val > 0.0:
                return val
    nested = meta.get(role)
    if isinstance(nested, Mapping) and nested.get("pre_var") is not None:
        try:
            val = float(nested["pre_var"])
        except (TypeError, ValueError):
            val = float("nan")
        if math.isfinite(val) and val > 0.0:
            return val
    return None


def _lookup_meta_n(res: RankKLongBetResult, role: str) -> float:
    """Read arm unit count from ``arm_metadata`` or infer from ``catt_draws``."""
    meta = getattr(res, "arm_metadata", None)
    if isinstance(meta, Mapping):
        direct_key = "n_treated" if role == "treated" else "n_control"
        if direct_key in meta and meta[direct_key] is not None:
            return max(float(meta[direct_key]), 1.0)
        nested = meta.get(role)
        if isinstance(nested, Mapping):
            for nk in ("n", "n_units"):
                if nk in nested and nested[nk] is not None:
                    return max(float(nested[nk]), 1.0)
    cohorts_arr, ever_arr, _ = _infer_panel_masks_from_result(res, None, None, None)
    n_trt = float(np.sum(ever_arr))
    n_ctrl = float(max(cohorts_arr.size - int(n_trt), 1))
    return max(n_trt if role == "treated" else n_ctrl, 1.0)


def _extract_outcome_summaries(
    results_map: Mapping[tuple[str, str], RankKLongBetResult],
    arms: tuple[str, ...],
    arm_to_idx: Mapping[str, int],
    fitted_pairs: tuple[tuple[str, str], ...],
    psi: np.ndarray,
    horizon: int,
) -> dict[str, Any]:
    """Extract basis projections ``(y_p, S_p, V_a, F_p)`` and out-of-basis trajectory stats."""
    n_arms = len(arms)
    k_dim = int(psi.shape[1])
    y_by_pair: dict[tuple[str, str], np.ndarray] = {}
    s_by_pair: dict[tuple[str, str], np.ndarray] = {}
    raw_s_full_by_pair: dict[tuple[str, str], np.ndarray] = {}
    active_modes_by_pair: dict[tuple[str, str], tuple[int, ...]] = {}
    traj_mean_by_pair: dict[tuple[str, str], np.ndarray] = {}
    traj_cov_by_pair: dict[tuple[str, str], np.ndarray] = {}
    v_candidates: list[list[np.ndarray]] = [[] for _ in range(n_arms)]
    pre_var_candidates: list[list[float]] = [[] for _ in range(n_arms)]
    n_candidates: list[list[float]] = [[] for _ in range(n_arms)]
    fallback_ssq: list[list[float]] = [[] for _ in range(n_arms)]

    for pair in fitted_pairs:
        if pair not in results_map:
            rev = (pair[1], pair[0])
            if rev in results_map:
                raise ValueError(
                    f"Pair {pair!r} missing, found reversed {rev!r}; all outcomes must use identical pair orientation"
                )
            raise ValueError(f"Missing RankKLongBetResult for contrast pair {pair!r}")
        res = results_map[pair]
        raw_traj = _extract_raw_trt_trajectory(res)
        if raw_traj.shape[1] < horizon:
            raise ValueError(
                f"Contrast {pair!r} has only {raw_traj.shape[1]} post-adoption period(s), fewer than horizon={horizon}"
            )
        traj_h = np.asarray(raw_traj[:, :horizon], dtype=np.float64)
        if not np.all(np.isfinite(traj_h)):
            raise ValueError(f"Non-finite post-adoption trajectory draws in contrast {pair!r}")
        traj_mean_by_pair[pair] = np.mean(traj_h, axis=0)
        raw_cov_h = np.atleast_2d(np.cov(traj_h, rowvar=False, ddof=1))
        traj_cov_by_pair[pair] = 0.5 * (raw_cov_h + raw_cov_h.T)

        proj_draws = traj_h @ psi  # (S, K')
        y_full = np.mean(proj_draws, axis=0)
        raw_s = np.atleast_2d(np.cov(proj_draws, rowvar=False, ddof=1))
        raw_s = 0.5 * (raw_s + raw_s.T)
        raw_s_full_by_pair[pair] = raw_s

        rk_attr = getattr(res, "rank_k", None)
        rk_cap = (
            int(rk_attr)
            if (isinstance(rk_attr, (int, np.integer)) and not isinstance(rk_attr, (bool, np.bool_)) and int(rk_attr) > 0)
            else k_dim
        )
        var_thresh = 1e-10 * max(float(np.trace(raw_s)), 1e-12)
        act_list = [k for k in range(min(k_dim, rk_cap)) if float(raw_s[k, k]) > var_thresh]
        if not act_list:
            act_list = [0]
        act_modes = tuple(act_list)
        active_modes_by_pair[pair] = act_modes
        y_by_pair[pair] = y_full[list(act_modes)]
        s_by_pair[pair] = _psd_clip(raw_s[np.ix_(act_modes, act_modes)], floor=1e-12)

        ia, ib = arm_to_idx[pair[0]], arm_to_idx[pair[1]]
        meta = getattr(res, "arm_metadata", None)
        cov_trt = _lookup_meta_cov(meta, "treated")
        cov_ctrl = _lookup_meta_cov(meta, "control")
        if cov_trt is not None and cov_trt.shape[0] >= horizon:
            v_candidates[ia].append(_psd_clip(psi.T @ cov_trt[:horizon, :horizon] @ psi, floor=1e-12))
        if cov_ctrl is not None and cov_ctrl.shape[0] >= horizon:
            v_candidates[ib].append(_psd_clip(psi.T @ cov_ctrl[:horizon, :horizon] @ psi, floor=1e-12))

        pre_trt = _lookup_meta_pre_var(meta, "treated")
        pre_ctrl = _lookup_meta_pre_var(meta, "control")
        if pre_trt is not None:
            pre_var_candidates[ia].append(pre_trt)
        if pre_ctrl is not None:
            pre_var_candidates[ib].append(pre_ctrl)

        n_trt = _lookup_meta_n(res, "treated")
        n_ctrl = _lookup_meta_n(res, "control")
        n_candidates[ia].append(n_trt)
        n_candidates[ib].append(n_ctrl)

        # Stage-1 residual s_sq fallback strictly for synthetic/mock RankKLongBetResult without arm_metadata
        gp_list = getattr(res, "gp_hyper_draws", None) or []
        if gp_list and isinstance(gp_list[0], Mapping) and "s_sq" in gp_list[0]:
            ssq_arr = np.asarray(gp_list[0]["s_sq"], dtype=np.float64)
            if ssq_arr.size > 0 and np.all(np.isfinite(ssq_arr)) and float(np.mean(ssq_arr)) > 0.0:
                ssq_mean = float(np.mean(ssq_arr))
                fallback_ssq[ia].append(ssq_mean / n_trt)
                fallback_ssq[ib].append(ssq_mean / n_ctrl)

    # Resolve V_{m,a} for every arm a = 0..A-1
    v_arms: list[np.ndarray] = []
    missing_arms: list[int] = []
    for a in range(n_arms):
        if v_candidates[a]:
            v_arms.append(_psd_clip(np.mean(np.stack(v_candidates[a], axis=0), axis=0), floor=1e-12))
        elif fallback_ssq[a]:
            v_iso = float(np.mean(fallback_ssq[a]))
            v_arms.append(v_iso * np.eye(k_dim, dtype=np.float64))
        else:
            v_arms.append(np.zeros((k_dim, k_dim), dtype=np.float64))
            missing_arms.append(a)

    if missing_arms:
        # Fallback NNLS variance allocation from S_{m,p} inversely proportional to n_a (mock unit tests only)
        inv_n = np.array(
            [1.0 / max(float(np.mean(n_candidates[a])) if n_candidates[a] else 1.0, 1.0) for a in range(n_arms)],
            dtype=np.float64,
        )
        p_len = len(fitted_pairs)
        a_inc = np.zeros((p_len, n_arms), dtype=np.float64)
        rhs = np.zeros(p_len, dtype=np.float64)
        for p_i, pair in enumerate(fitted_pairs):
            ia, ib = arm_to_idx[pair[0]], arm_to_idx[pair[1]]
            a_inc[p_i, ia] = 1.0
            a_inc[p_i, ib] = 1.0
            k_p = max(len(active_modes_by_pair[pair]), 1)
            rhs[p_i] = 0.80 * float(np.trace(s_by_pair[pair])) / float(k_p)
        # Regularized NNLS toward sample-size-proportional baseline
        scale_guess = float(np.median(rhs / np.maximum(a_inc @ inv_n, 1e-12)))
        v_prior = scale_guess * inv_n
        reg = 0.15
        a_aug = np.vstack([a_inc, math.sqrt(reg) * np.eye(n_arms)])
        b_aug = np.concatenate([rhs, math.sqrt(reg) * v_prior])
        sol, _ = sp_optimize.nnls(a_aug, b_aug)
        for a in missing_arms:
            v_arms[a] = max(float(sol[a]), 1e-12) * np.eye(k_dim, dtype=np.float64)

    # Align mode-wise Stage-1 shrinkage strictly across pairs that actually estimate mode k
    gamma_modes = np.ones(k_dim, dtype=np.float64)
    for k in range(k_dim):
        ratios_k = [
            float(raw_s_full_by_pair[pair][k, k])
            / max(float(v_arms[arm_to_idx[pair[0]]][k, k] + v_arms[arm_to_idx[pair[1]]][k, k]), 1e-18)
            for pair in fitted_pairs
            if k in active_modes_by_pair[pair]
        ]
        if ratios_k:
            gamma_modes[k] = min(1.0, max(float(np.median(ratios_k)), 1e-6))
        else:
            gamma_modes[k] = 1.0
    sqrt_gamma = np.sqrt(gamma_modes)
    for a in range(n_arms):
        v_arms[a] = _psd_clip((sqrt_gamma[:, None] * v_arms[a]) * sqrt_gamma[None, :], floor=1e-12)

    # Precompute simultaneous diagonalization of (S_{m,p}, V^{(0)}_{m,a} + V^{(0)}_{m,b})
    # on each pair's active mode subspace so F_{m,p}(rho_m, phi_m), L_F, and F^{-1} update in O(K_p^2)
    rho0 = _RHO_BETA_ALPHA / (_RHO_BETA_ALPHA + _RHO_BETA_BETA)
    phi0 = _PHI_BETA_ALPHA / (_PHI_BETA_ALPHA + _PHI_BETA_BETA)
    f_by_pair: dict[tuple[str, str], np.ndarray] = {}
    w_f_by_pair: dict[tuple[str, str], np.ndarray] = {}
    winv_f_by_pair: dict[tuple[str, str], np.ndarray] = {}
    d_ratio_by_pair: dict[tuple[str, str], np.ndarray] = {}
    for pair in fitted_pairs:
        ia, ib = arm_to_idx[pair[0]], arm_to_idx[pair[1]]
        act = list(active_modes_by_pair[pair])
        s_mp = _psd_clip(s_by_pair[pair], floor=1e-14)
        v_sum = (v_arms[ia] + v_arms[ib])[np.ix_(act, act)]
        evals_s, evecs_s = np.linalg.eigh(s_mp)
        evals_s = np.maximum(evals_s, 1e-14)
        s_half = (evecs_s * np.sqrt(evals_s)[None, :]) @ evecs_s.T
        s_inv_half = (evecs_s * (1.0 / np.sqrt(evals_s))[None, :]) @ evecs_s.T
        sym_h = s_inv_half @ v_sum @ s_inv_half
        sym_h = 0.5 * (sym_h + sym_h.T)
        d_vals, u_h = np.linalg.eigh(sym_h)
        d_vals = np.maximum(d_vals, 0.0)
        w_q = s_half @ u_h
        winv_q = u_h.T @ s_inv_half
        w_f_by_pair[pair] = w_q
        winv_f_by_pair[pair] = winv_q
        d_ratio_by_pair[pair] = d_vals
        e0 = np.maximum(1.0 - rho0 * d_vals, phi0)
        f_by_pair[pair] = (w_q * e0[None, :]) @ w_q.T

    # Realization-independent baseline scale s_{0,m}: strictly from pre-period outcome SD when
    # available in Stage-1 arm_metadata (falling back to 2x top-contrast Mode-1 SE only on synthetic
    # mock results that lack pre-period variance metadata).
    top_arm = arms[-1]
    ctrl_arm = arms[0]
    if (top_arm, ctrl_arm) in raw_s_full_by_pair:
        s0_var = max(float(raw_s_full_by_pair[(top_arm, ctrl_arm)][0, 0]), 1e-16) / float(horizon)
    else:
        b_pinv = np.linalg.pinv(
            np.array(
                [
                    [(1.0 if arm_to_idx[u] == j else (-1.0 if arm_to_idx[v] == j else 0.0)) for j in range(1, n_arms)]
                    for u, v in fitted_pairs
                ],
                dtype=np.float64,
            )
        )
        top_row = b_pinv[-1]
        s0_var = (
            float(sum((top_row[i] ** 2) * max(float(raw_s_full_by_pair[p][0, 0]), 1e-16) for i, p in enumerate(fitted_pairs)))
            / float(horizon)
        )
    se_top_scale = max(2.0 * math.sqrt(max(s0_var, 1e-16)), 1e-6)

    any_pre_vars = [float(np.mean(pv)) for pv in pre_var_candidates if pv and float(np.mean(pv)) > 0.0]
    if any_pre_vars:
        pre_sd = math.sqrt(max(float(np.median(any_pre_vars)), 1e-16))
        s0_scale = max(_SCALE_PRIOR_C0 * pre_sd, 1e-6)
    else:
        s0_scale = max(_SCALE_PRIOR_C0 * se_top_scale / max(_SCALE_PRIOR_TAU0, 1e-6), 1e-6)
    sk_scale = s0_scale

    return {
        "y_by_pair": y_by_pair,
        "s_by_pair": s_by_pair,
        "active_modes_by_pair": active_modes_by_pair,
        "v_arms": v_arms,
        "f_by_pair": f_by_pair,
        "w_f_by_pair": w_f_by_pair,
        "winv_f_by_pair": winv_f_by_pair,
        "d_ratio_by_pair": d_ratio_by_pair,
        "gamma_modes": gamma_modes,
        "traj_mean_by_pair": traj_mean_by_pair,
        "traj_cov_by_pair": traj_cov_by_pair,
        "s0_scale": s0_scale,
        "sk_scale": sk_scale,
        "se_top_scale": se_top_scale,
    }


# ---- Polyhedral cone constraints & exact cone normalizer Z_m(r, lambda) -------------------------


def _build_cone_matrix_outcome(
    c_cover: np.ndarray,
    psi: np.ndarray,
    horizon: int,
    monotone: bool,
    period_positivity: bool,
) -> np.ndarray:
    """Build ``A_m`` such that ``A_m @ vec_colmajor(Theta_m) >= 0`` defines ``K_m``.

    Parameter ordering for outcome ``m`` is mode-major ``(Theta_{m,:,0}, ..., Theta_{m,:,K'-1})``
    of length ``(A - 1) * K'``. When ``K' == 2`` and ``period_positivity=True``, the ``H + 1``
    per-edge inequalities reduce algebraically to the 2 extreme envelope walls
    ``u_{e,0} + m_min * u_{e,1} >= 0`` and ``u_{e,0} + m_max * u_{e,1} >= 0``.
    """
    n_nonctrl = int(c_cover.shape[1])
    k_dim = int(psi.shape[1])
    d_m = n_nonctrl * k_dim
    if not monotone or c_cover.shape[0] == 0:
        return np.zeros((0, d_m), dtype=np.float64)

    if not period_positivity or k_dim == 1:
        a_mat = np.zeros((c_cover.shape[0], d_m), dtype=np.float64)
        a_mat[:, :n_nonctrl] = c_cover
        return a_mat

    psi_scaled = math.sqrt(float(horizon)) * psi  # column 0 == 1.0
    if k_dim == 2:
        p1 = psi_scaled[:, 1]
        m_min = float(np.min(p1))
        m_max = float(np.max(p1))
        if m_max - m_min <= 1e-12:
            a_mat = np.zeros((c_cover.shape[0], d_m), dtype=np.float64)
            a_mat[:, :n_nonctrl] = c_cover
            return a_mat
        rows: list[np.ndarray] = []
        for e_idx in range(c_cover.shape[0]):
            ce = c_cover[e_idx]
            rows.append(np.concatenate([ce, m_min * ce]))
            rows.append(np.concatenate([ce, m_max * ce]))
        return np.vstack(rows)

    # General K' >= 3: Mode-1 constraint + H period constraints
    rows = []
    for e_idx in range(c_cover.shape[0]):
        ce = c_cover[e_idx]
        r0 = np.zeros(d_m, dtype=np.float64)
        r0[:n_nonctrl] = ce
        rows.append(r0)
        for t_idx in range(horizon):
            rt = np.concatenate([psi_scaled[t_idx, k] * ce for k in range(k_dim)])
            rows.append(rt)
    return np.vstack(rows)


class _OutcomeConePriorSpec:
    """Precomputed prior precision structure and exact cone normalizer ``log Z_m(r, lambda)``."""

    def __init__(
        self,
        *,
        c_cover: np.ndarray,
        psi: np.ndarray,
        horizon: int,
        monotone: bool,
        period_positivity: bool,
        scale_prior: str,
        use_local_scales: bool | None = None,
    ) -> None:
        self.c_cover = np.asarray(c_cover, dtype=np.float64)
        self.n_edges, self.n_nonctrl = self.c_cover.shape
        self.psi = np.asarray(psi, dtype=np.float64)
        self.horizon = int(horizon)
        self.k_dim = int(psi.shape[1])
        self.monotone = bool(monotone)
        self.period_positivity = bool(period_positivity)
        self.scale_prior = str(scale_prior)

        base_prec = self.c_cover.T @ self.c_cover
        if np.linalg.matrix_rank(base_prec) < self.n_nonctrl:
            base_prec = base_prec + 1e-6 * np.eye(self.n_nonctrl)
        base_cov = np.linalg.inv(base_prec)
        self.c_scale = max(float(base_cov[-1, -1]), 1e-12)
        self.lambda0_unit = self.c_scale * base_prec
        _, logdet_u = np.linalg.slogdet(self.lambda0_unit)
        self.logdet_unit = float(logdet_u)

        self.is_tree_poset = (
            self.n_edges == self.n_nonctrl
            and np.linalg.matrix_rank(self.c_cover) == self.n_nonctrl
        )
        self.c_cover_inv: np.ndarray | None = (
            np.linalg.inv(self.c_cover) if self.is_tree_poset else None
        )
        if use_local_scales is None:
            self.use_local_scales = bool(
                self.n_edges > 0
                and (
                    self.scale_prior == "local_halft"
                    or (self.scale_prior == "global_halft" and _SCALE_PRIOR_TREE_LOCAL and self.is_tree_poset)
                )
            )
        else:
            self.use_local_scales = bool(use_local_scales and self.n_edges > 0)

        # Detect standard 2x2 diamond (A-1 == 3, E_c == 4, arms 10 and 01 conditionally independent given 11)
        self.is_2x2_diamond = (
            self.n_nonctrl == 3
            and self.n_edges == 4
            and abs(float(base_prec[0, 1])) < 1e-12
            and np.allclose(np.sort(np.abs(self.c_cover), axis=0), np.sort(np.abs(self.c_cover), axis=0))
        )
        if self.is_2x2_diamond:
            # Map which edges are (10, 00), (01, 00), (11, 10), (11, 01)
            self._diamond_map = self._identify_diamond_edges()
        else:
            self._diamond_map = None

        psi_scaled = math.sqrt(float(self.horizon)) * self.psi
        if self.k_dim >= 2:
            p1 = psi_scaled[:, 1]
            self.m_plus = max(float(np.max(p1)), 1e-12)
            self.m_minus = max(float(-np.min(p1)), 1e-12)
        else:
            self.m_plus = 0.0
            self.m_minus = 0.0

        self.a_cone = _build_cone_matrix_outcome(
            self.c_cover, self.psi, self.horizon, self.monotone, self.period_positivity
        )
        self.a_cone_0 = self.a_cone[:, : self.n_nonctrl]
        self.a_cone_dyn = self.a_cone[:, self.n_nonctrl :]

        # Deterministic QMC normal bank for general non-tree K' >= 3 fallback
        self._qmc_z: np.ndarray | None = None
        self._p1_lut: list[float] | None = None
        if self.monotone and not self.is_tree_poset:
            if not self.use_local_scales and self.period_positivity and self.k_dim >= 2:
                self._p1_lut = self._build_p1_logz_lut()

    def _identify_diamond_edges(self) -> tuple[int, int, int, int] | None:
        """Identify edge indices ``(e_10_0, e_01_0, e_11_10, e_11_01)`` if standard 2x2 diamond."""
        e_10_0 = e_01_0 = e_11_10 = e_11_01 = -1
        for idx, row in enumerate(self.c_cover):
            tup = tuple(int(round(x)) for x in row)
            if tup == (1, 0, 0):
                e_10_0 = idx
            elif tup == (0, 1, 0):
                e_01_0 = idx
            elif tup == (-1, 0, 1):
                e_11_10 = idx
            elif tup == (0, -1, 1):
                e_11_01 = idx
        if min(e_10_0, e_01_0, e_11_10, e_11_01) >= 0:
            return (e_10_0, e_01_0, e_11_10, e_11_01)
        return None

    def precision_and_logdet(self, local_scales: np.ndarray) -> tuple[np.ndarray, float]:
        """Return ``(Lambda_0(lambda), log det Lambda_0(lambda))``."""
        if not self.use_local_scales:
            return self.lambda0_unit, self.logdet_unit
        lam = np.maximum(np.asarray(local_scales, dtype=np.float64), 1e-12)
        inv_sq = 1.0 / (lam * lam)
        prec = self.c_scale * ((self.c_cover.T * inv_sq[None, :]) @ self.c_cover)
        if self.is_tree_poset:
            return prec, self.logdet_unit - 2.0 * float(np.sum(np.log(lam)))
        sign, logdet = np.linalg.slogdet(prec)
        if sign <= 0:
            prec = _psd_clip(prec, floor=1e-12)
            _, logdet = np.linalg.slogdet(prec)
        return prec, float(logdet)

    def _ensure_qmc_bank(self) -> np.ndarray:
        if self._qmc_z is None:
            d_tot = self.n_nonctrl * self.k_dim
            sampler = sp_qmc.Sobol(d=d_tot, scramble=True, seed=20260311)
            u = np.clip(sampler.random_base2(m=12), 1e-9, 1.0 - 1e-9)  # 4096 points
            z_half = sp_special.ndtri(u)
            self._qmc_z = np.concatenate([z_half, -z_half], axis=0)  # 8192 antithetic points
        return self._qmc_z

    def _diamond_logz_exact(
        self, ratio: float, local_scales: np.ndarray, *, fast_grid: bool = False
    ) -> float:
        """Deterministic Gauss-Legendre exact cone probability on 2x2 diamond."""
        assert self._diamond_map is not None
        e10, e01, e31, e32 = self._diamond_map
        lam = np.asarray(local_scales, dtype=np.float64)
        tau1 = 1.0 / max(float(lam[e10]) ** 2, 1e-20)
        tau2 = 1.0 / max(float(lam[e01]) ** 2, 1e-20)
        tau3 = 1.0 / max(float(lam[e31]) ** 2, 1e-20)
        tau4 = 1.0 / max(float(lam[e32]) ** 2, 1e-20)

        w1 = tau3 / (tau1 + tau3)
        w2 = tau4 / (tau2 + tau4)
        v1 = 1.0 / (self.c_scale * (tau1 + tau3))
        v2 = 1.0 / (self.c_scale * (tau2 + tau4))
        prec_33 = self.c_scale * (
            (tau1 * tau3) / (tau1 + tau3) + (tau2 * tau4) / (tau2 + tau4)
        )
        var_3 = 1.0 / max(prec_33, 1e-20)
        sd_3 = math.sqrt(var_3)
        sd_v1 = math.sqrt(v1)
        sd_v2 = math.sqrt(v2)

        if not self.period_positivity or self.k_dim == 1 or ratio <= 1e-12:
            # 1D integral over x_{11,0} >= 0: p in [0.5, 1.0)
            p_nodes = 0.75 + 0.25 * _GL16_X
            p_weights = 0.25 * _GL16_W
            x3 = sd_3 * sp_special.ndtri(np.clip(p_nodes, 0.5, 1.0 - 1e-15))
            p10 = sp_special.ndtr(((1.0 - w1) * x3) / sd_v1) - sp_special.ndtr((-w1 * x3) / sd_v1)
            p01 = sp_special.ndtr(((1.0 - w2) * x3) / sd_v2) - sp_special.ndtr((-w2 * x3) / sd_v2)
            prob = float(np.sum(p_weights * np.maximum(p10, 0.0) * np.maximum(p01, 0.0)))
            return math.log(max(prob, 1e-300))

        if self.k_dim == 2:
            gl_x = _GL8_X if fast_grid else _GL16_X
            gl_w = _GL8_W if fast_grid else _GL16_W
            w0 = 0.25 * gl_w
            if fast_grid:
                x3_0 = sd_3 * _GL8_P0_Z
            else:
                p0 = 0.75 + 0.25 * gl_x
                x3_0 = sd_3 * sp_special.ndtri(np.clip(p0, 0.5, 1.0 - 1e-15))

            mp, mm = self.m_plus, self.m_minus
            m_sum = mp + mm
            r_sd3 = ratio * sd_3
            q_lo = sp_special.ndtr((-x3_0 / mp) / r_sd3)
            q_hi = sp_special.ndtr((x3_0 / mm) / r_sd3)
            q_span = np.maximum(q_hi - q_lo, 0.0)
            q_grid = q_lo[:, None] + 0.5 * q_span[:, None] * (1.0 + gl_x[None, :])
            w_q = 0.5 * q_span[:, None] * gl_w[None, :]
            x3_1 = r_sd3 * sp_special.ndtri(np.clip(q_grid, 1e-15, 1.0 - 1e-15))

            # Exact piecewise-linear support [lb, ub] of x_{b,1} where low_0(x_{b,1}) < high_0(x_{b,1})
            lb = np.maximum(
                -x3_0[:, None] / mp,
                np.maximum(x3_1 - x3_0[:, None] / mm, (mm * x3_1 - x3_0[:, None]) / m_sum),
            )
            ub = np.minimum(
                x3_0[:, None] / mm,
                np.minimum(x3_1 + x3_0[:, None] / mp, (mp * x3_1 + x3_0[:, None]) / m_sum),
            )

            def _branch_prob(w_b: float, sd_b: float) -> np.ndarray:
                r_sdb = ratio * sd_b
                mu_1 = w_b * x3_1
                c_lo = sp_special.ndtr((lb - mu_1) / r_sdb)
                c_hi = sp_special.ndtr((ub - mu_1) / r_sdb)
                c_span = np.where(ub > lb, np.maximum(c_hi - c_lo, 0.0), 0.0)
                c_grid = c_lo[:, :, None] + 0.5 * c_span[:, :, None] * (1.0 + gl_x[None, None, :])
                w_c = 0.5 * c_span[:, :, None] * gl_w[None, None, :]
                xb_1 = mu_1[:, :, None] + r_sdb * sp_special.ndtri(np.clip(c_grid, 1e-15, 1.0 - 1e-15))
                diff_1 = x3_1[:, :, None] - xb_1
                low_0 = np.maximum(0.0, np.maximum(-mp * xb_1, mm * xb_1))
                high_0 = x3_0[:, None, None] - np.maximum(0.0, np.maximum(-mp * diff_1, mm * diff_1))
                mu_0 = w_b * x3_0[:, None, None]
                p_int = np.where(
                    high_0 > low_0,
                    np.maximum(
                        sp_special.ndtr((high_0 - mu_0) / sd_b) - sp_special.ndtr((low_0 - mu_0) / sd_b),
                        0.0,
                    ),
                    0.0,
                )
                return np.sum(w_c * p_int, axis=2)

            prob_10 = _branch_prob(w1, sd_v1)
            prob_01 = _branch_prob(w2, sd_v2)
            total = float(np.sum(w0[:, None] * w_q * prob_10 * prob_01))
            return math.log(max(total, 1e-300))

        return self._qmc_logz(ratio, local_scales)

    def _qmc_logz(self, ratio: float, local_scales: np.ndarray) -> float:
        """Smooth QMC fallback for general non-tree posets with K' >= 3."""
        z_bank = self._ensure_qmc_bank()
        prec, _ = self.precision_and_logdet(local_scales)
        cov = np.linalg.inv(prec)
        chol = np.linalg.cholesky(_psd_clip(cov, floor=1e-14))
        n_pts = z_bank.shape[0]
        # Map z_bank (N, (A-1)*K') -> Theta samples in mode-major order
        blocks = []
        for k in range(self.k_dim):
            zk = z_bank[:, k * self.n_nonctrl : (k + 1) * self.n_nonctrl]
            scale_k = 1.0 if k == 0 else ratio
            blocks.append(scale_k * (zk @ chol.T))
        pts = np.concatenate(blocks, axis=1)
        inside = np.all((pts @ self.a_cone.T) >= 0.0, axis=1)
        prob = max(float(np.mean(inside)), 0.5 / float(n_pts))
        return math.log(prob)

    def _build_p1_logz_lut(self) -> list[float]:
        """Precompute dense 2,049-point lookup table of ``log Z_m(exp(u), 1)`` on ``u in [-12, 12]``."""
        cache_key = (
            self.horizon,
            self.k_dim,
            round(self.m_plus, 12),
            round(self.m_minus, 12),
            round(self.c_scale, 12),
            self._diamond_map,
            self.c_cover.tobytes(),
        )
        cached = _P1_LUT_CACHE.get(cache_key)
        if cached is not None:
            return cached
        u_grid = np.linspace(-12.0, 12.0, 97)
        ones = np.ones(self.n_edges, dtype=np.float64)
        vals = np.array(
            [
                self._diamond_logz_exact(float(math.exp(u)), ones)
                if self._diamond_map is not None
                else self._qmc_logz(float(math.exp(u)), ones)
                for u in u_grid
            ],
            dtype=np.float64,
        )
        spline = sp_interpolate.CubicSpline(u_grid, vals, extrapolate=True)
        u_dense = np.linspace(-12.0, 12.0, 2049)
        lut = [float(x) for x in spline(u_dense)]
        _P1_LUT_CACHE[cache_key] = lut
        return lut

    def log_normalizer(self, ratio: float, local_scales: np.ndarray) -> float:
        """Return ``log Pr((X_0, ratio * X_{1..K'-1}) in K_m | Lambda_0(local_scales))``."""
        if not self.monotone or self.n_edges == 0:
            return 0.0
        r_f = float(ratio)
        r_clamped = 1e-8 if r_f < 1e-8 else (1e8 if r_f > 1e8 else r_f)
        if self.is_tree_poset:
            if not self.period_positivity or self.k_dim == 1:
                return -float(self.n_nonctrl) * math.log(2.0)
            if self.k_dim == 2:
                p_edge = (
                    math.atan(1.0 / (r_clamped * self.m_plus))
                    + math.atan(1.0 / (r_clamped * self.m_minus))
                ) / (2.0 * math.pi)
                return float(self.n_nonctrl) * math.log(max(p_edge, 1e-300))
            return self._qmc_logz(r_clamped, np.ones(self.n_edges, dtype=np.float64))

        if not self.use_local_scales:
            if not self.period_positivity or self.k_dim == 1:
                if self._diamond_map is not None:
                    return self._diamond_logz_exact(0.0, np.ones(self.n_edges))
                return self._qmc_logz(0.0, np.ones(self.n_edges))
            lut = self._p1_lut
            assert lut is not None
            log_r = math.log(r_clamped)
            if log_r >= 12.0:
                # Asymptotic tail: Z(r) ~ r^{-(A-1)(K'-1)} for large r
                return lut[2048] - float(self.n_nonctrl * (self.k_dim - 1)) * (log_r - 12.0)
            if log_r <= -12.0:
                return lut[0]
            pos = (log_r + 12.0) * (2048.0 / 24.0)
            idx = int(pos)
            if idx >= 2048:
                return lut[2048]
            y0 = lut[idx]
            return y0 + (pos - float(idx)) * (lut[idx + 1] - y0)

        if self._diamond_map is not None:
            return self._diamond_logz_exact(r_clamped, local_scales, fast_grid=True)
        return self._qmc_logz(r_clamped, local_scales)


# ---- Exact Pakman--Paninski (2014) truncated-Gaussian HMC & Gibbs samplers ----------------------


def _sample_univariate_truncnorm(
    mean: float,
    sd: float,
    low: float,
    high: float,
    rng: np.random.Generator,
) -> float:
    """Exact 1-D truncated normal draw on ``[low, high]`` stable deep into the tails."""
    if high <= low + 1e-15:
        return 0.5 * (low + high) if math.isfinite(low) and math.isfinite(high) else (low if math.isfinite(low) else high)
    a = (low - mean) / sd if math.isfinite(low) else -math.inf
    b = (high - mean) / sd if math.isfinite(high) else math.inf
    if a <= -6.0 and b >= 6.0:
        z = float(rng.standard_normal())
        return mean + sd * float(np.clip(z, a, b))
    if a >= 4.5:
        # Robert (1995) shifted-exponential tail rejection for far right tail [a, b]
        lam = 0.5 * (a + math.sqrt(a * a + 4.0))
        for _ in range(128):
            u_exp = float(rng.uniform(1e-15, 1.0))
            if math.isfinite(b):
                max_u = 1.0 - math.exp(-lam * (b - a))
                z = a - math.log(max(1.0 - u_exp * max_u, 1e-300)) / lam
            else:
                z = a - math.log(u_exp) / lam
            if math.log(float(rng.uniform(1e-15, 1.0))) <= -0.5 * (z - lam) ** 2:
                return mean + sd * float(np.clip(z, a, b))
        return mean + sd * a
    if b <= -4.5:
        return -_sample_univariate_truncnorm(-mean, sd, -high, -low, rng)

    cdf_a = float(sp_special.ndtr(a)) if math.isfinite(a) else 0.0
    cdf_b = float(sp_special.ndtr(b)) if math.isfinite(b) else 1.0
    diff = cdf_b - cdf_a
    if diff > 1e-12:
        u = float(rng.uniform(cdf_a, cdf_b))
        z = float(sp_special.ndtri(np.clip(u, 1e-15, 1.0 - 1e-15)))
        return mean + sd * float(np.clip(z, a, b))
    # Narrow interval fallback
    return float(rng.uniform(low, high)) if (math.isfinite(low) and math.isfinite(high)) else (mean + sd * float(np.clip(0.0, a, b)))


def _univariate_tn_gibbs_sweep(
    z: np.ndarray,
    f_walls: np.ndarray,
    g_walls: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """One exact coordinate-wise truncated-normal Gibbs sweep on ``z ~ N(0, I_D)`` s.t. ``F z + g >= 0``."""
    z_cur = z.copy()
    slack = f_walls @ z_cur + g_walls
    d_dim = z_cur.shape[0]
    for i in rng.permutation(d_dim):
        col = f_walls[:, i]
        slack_wo = slack - col * z_cur[i]
        pos_mask = col > 1e-13
        neg_mask = col < -1e-13
        low = float(np.max(-slack_wo[pos_mask] / col[pos_mask])) if np.any(pos_mask) else -math.inf
        high = float(np.min(-slack_wo[neg_mask] / col[neg_mask])) if np.any(neg_mask) else math.inf
        z_new_i = _sample_univariate_truncnorm(0.0, 1.0, low, high, rng)
        slack = slack_wo + col * z_new_i
        z_cur[i] = z_new_i
    return z_cur


_LAPACK_DPOTRF = sp_linalg.lapack.dpotrf
_LAPACK_DTRTRS = sp_linalg.lapack.dtrtrs
_LAPACK_DTRTRI = sp_linalg.lapack.dtrtri

# Fixed (state-independent) wall-bounce budget per Pakman--Paninski trajectory; trajectories
# exceeding it are rejected (the chain stays put), which preserves exact reversibility.
_PP_HMC_MAX_BOUNCES = 48


def _fast_cholesky(mat: np.ndarray, *, floor: float = 1e-12) -> np.ndarray:
    """Lower Cholesky factor via direct LAPACK ``dpotrf``, eigenvalue-clipping only on failure.

    ``mat`` must be symmetric, so its transpose (an F-contiguous view of a C-contiguous
    array) is passed to LAPACK to avoid a layout copy.
    """
    chol, info = _LAPACK_DPOTRF(mat.T, lower=1, clean=1)
    if info == 0:
        return chol
    return np.linalg.cholesky(_psd_clip(mat, floor=floor))


def _solve_lower(chol: np.ndarray, rhs: np.ndarray, *, transpose: bool = False) -> np.ndarray:
    """Solve ``L x = rhs`` (or ``L^T x = rhs`` if ``transpose``) for lower-triangular ``L``."""
    sol, _ = _LAPACK_DTRTRS(chol, rhs, lower=1, trans=1 if transpose else 0)
    return sol


def _inv_lower(chol: np.ndarray) -> np.ndarray:
    """Inverse of a lower-triangular Cholesky factor via LAPACK ``dtrtri``."""
    inv, _ = _LAPACK_DTRTRI(chol, lower=1)
    return inv


def _pakman_paninski_hmc_whitened(
    z_init: np.ndarray,
    f_walls: np.ndarray,
    g_walls: np.ndarray,
    rng: np.random.Generator,
    *,
    n_steps: int = 2,
    max_events: int | None = None,
) -> np.ndarray:
    """Exact Pakman & Paninski (2014) wall-reflecting HMC for ``z ~ N(0, I_D)`` s.t. ``F z + g >= 0``.

    Every row of ``f_walls`` must have unit Euclidean norm ``||f_j||_2 == 1``.
    Travel time is randomized ``T ~ Uniform(0.35 * pi, 0.65 * pi)`` (mean ``pi / 2``)
    to eliminate periodic harmonic orbits while preserving exact detailed balance.
    """
    m_walls, d_dim = f_walls.shape
    if m_walls == 0:
        return rng.standard_normal(d_dim)
    eff_max_events = _PP_HMC_MAX_BOUNCES if max_events is None else int(max_events)

    z = z_init.copy()
    # Ensure strict interior feasibility up to machine epsilon
    if float(np.min(f_walls @ z + g_walls)) < 1e-12:
        z = _univariate_tn_gibbs_sweep(z, f_walls, g_walls, rng)

    two_pi = 2.0 * math.pi
    abs_g_tol = np.abs(g_walls) + 1e-14
    for _ in range(n_steps):
        z_saved = z.copy()
        v = rng.standard_normal(d_dim)
        t_rem = float(rng.uniform(0.35 * math.pi, 0.65 * math.pi))
        last_wall = -1
        bounced_cleanly = True

        for _event in range(eff_max_events):
            a_vec = f_walls @ z
            b_vec = f_walls @ v
            u_vec = np.hypot(a_vec, b_vec)

            # Candidate outward wall hit: u_j * cos(t + phi_j) + g_j == 0 with sin(t + phi_j) > 0
            can_hit = u_vec > abs_g_tol
            if not np.any(can_hit):
                z = z * math.cos(t_rem) + v * math.sin(t_rem)
                break

            idx_cand = np.nonzero(can_hit)[0]
            u_c = u_vec[idx_cand]
            ratio_c = np.clip(-g_walls[idx_cand] / u_c, -1.0, 1.0)
            alpha_c = np.arccos(ratio_c)  # in [0, pi], sin(alpha_c) >= 0
            phi_c = np.arctan2(-b_vec[idx_cand], a_vec[idx_cand])
            t_c = np.mod(alpha_c - phi_c, two_pi)
            t_c = np.where((t_c < 1e-13) & (b_vec[idx_cand] >= 0.0), t_c + two_pi, t_c)
            if last_wall >= 0:
                lw_mask = idx_cand == last_wall
                if np.any(lw_mask):
                    t_c = np.where(
                        lw_mask,
                        np.where(alpha_c > 1e-11, 2.0 * alpha_c, two_pi),
                        t_c,
                    )

            valid_c = (u_c * np.sin(alpha_c) > 1e-13) & (t_c > 1e-14)
            if not np.any(valid_c):
                z = z * math.cos(t_rem) + v * math.sin(t_rem)
                break

            t_valid = np.where(valid_c, t_c, math.inf)
            k_best = int(np.argmin(t_valid))
            t_hit = float(t_valid[k_best])
            if t_hit >= t_rem:
                z = z * math.cos(t_rem) + v * math.sin(t_rem)
                break
            j_hit = int(idx_cand[k_best])

            # Advance to wall j_hit and reflect momentum elastically
            cos_h = math.cos(t_hit)
            sin_h = math.sin(t_hit)
            z_next = z * cos_h + v * sin_h
            v_next = -z * sin_h + v * cos_h
            fj = f_walls[j_hit]
            vn = float(np.dot(fj, v_next))
            if vn < 0.0:
                v_next = v_next - 2.0 * vn * fj
            # Snap exact wall coordinate to zero to prevent cumulative drift
            wall_err = float(np.dot(fj, z_next) + g_walls[j_hit])
            z_next = z_next - wall_err * fj

            z = z_next
            v = v_next
            t_rem -= t_hit
            last_wall = j_hit
            if t_rem <= 1e-14:
                break
        else:
            bounced_cleanly = False

        if not bounced_cleanly or float(np.min(f_walls @ z + g_walls)) < -1e-9:
            # Reject trajectories exceeding the fixed bounce budget (or numerically leaving K) and
            # stay at z_saved. The bounce count is invariant under the time-reversal involution
            # (z, v) -> (z', -v'), so the accept set is symmetric and the move is exactly
            # pi-reversible. (A state-dependent fallback move on rejection would not be.)
            z = z_saved

    return z


def _sample_truncated_mvn(
    mean: np.ndarray,
    precision: np.ndarray,
    a_cone: np.ndarray,
    theta_current: np.ndarray,
    interior_ray: np.ndarray,
    rng: np.random.Generator,
    *,
    n_hmc_steps: int = 2,
    chol_prec: np.ndarray | None = None,
) -> np.ndarray:
    """Sample ``theta ~ N(mean, precision^{-1})`` truncated to ``{theta : a_cone @ theta >= 0}``."""
    d_dim = mean.shape[0]
    if chol_prec is None:
        chol_prec = _fast_cholesky(precision, floor=1e-12)
    # precision = L L^T => cov = L^{-T} L^{-1}; set M_cov = L^{-T}
    m_cov = _inv_lower(chol_prec).T
    if a_cone.shape[0] == 0:
        return mean + m_cov @ rng.standard_normal(d_dim)

    f_raw = a_cone @ m_cov
    g_raw = a_cone @ mean
    row_norms = np.linalg.norm(f_raw, axis=1)
    keep = row_norms > 1e-300
    if not np.any(keep):
        return mean + m_cov @ rng.standard_normal(d_dim)
    f_walls = f_raw[keep] / row_norms[keep, None]
    g_walls = g_raw[keep] / row_norms[keep]

    z_init = chol_prec.T @ (theta_current - mean)
    min_w_slack = float(np.min(f_walls @ z_init + g_walls))
    if min_w_slack <= 0.0:
        ray_z = chol_prec.T @ interior_ray
        ray_w_slack = float(np.min(f_walls @ ray_z))
        z_init = z_init + ((-min_w_slack + 1e-12) / max(ray_w_slack, 1e-300)) * ray_z

    z_draw = _pakman_paninski_hmc_whitened(z_init, f_walls, g_walls, rng, n_steps=n_hmc_steps)
    theta_draw = mean + m_cov @ z_draw

    # Nudge tiny machine-rounding boundary violations (< 0) onto cone interior proportionally
    post_slack = float(np.min(a_cone @ theta_draw))
    if post_slack < 0.0:
        ray_slack = float(np.min(a_cone @ interior_ray))
        theta_draw = theta_draw + ((-post_slack * (1.0 + 1e-12) + 1e-300) / max(ray_slack, 1e-300)) * interior_ray
    return theta_draw


def _cone_qp_mode(
    mean: np.ndarray,
    precision: np.ndarray,
    a_cone: np.ndarray,
    interior_ray: np.ndarray,
) -> np.ndarray:
    """Compute exact cone-constrained posterior mode ``argmin 0.5 (x - mu)^T P (x - mu)`` s.t. ``A x >= 0``."""
    if a_cone.shape[0] == 0 or float(np.min(a_cone @ mean)) >= 1e-5:
        return mean.copy()
    chol_prec = _fast_cholesky(precision, floor=1e-12)
    m_cov = _inv_lower(chol_prec).T
    f_raw = a_cone @ m_cov
    g_raw = a_cone @ mean
    norms = np.linalg.norm(f_raw, axis=1)
    keep = norms > 1e-14
    f_w = f_raw[keep] / norms[keep, None]
    g_w = g_raw[keep] / norms[keep] - 1e-5  # small strict interior margin

    # Feasible starting point along interior_ray
    z0 = chol_prec.T @ (interior_ray - mean)
    s0 = float(np.min(f_w @ z0 + g_w))
    if s0 < 0.0:
        ray_z = chol_prec.T @ interior_ray
        rs = float(np.min(f_w @ ray_z))
        z0 = z0 + ((-s0 + 1e-3) / max(rs, 1e-12)) * ray_z

    res = sp_optimize.minimize(
        fun=lambda z: 0.5 * float(np.dot(z, z)),
        x0=z0,
        jac=lambda z: z,
        constraints={"type": "ineq", "fun": lambda z: f_w @ z + g_w, "jac": lambda _z: f_w},
        method="SLSQP",
        options={"ftol": 1e-10, "maxiter": 200},
    )
    z_opt = np.asarray(res.x, dtype=np.float64) if res.success else z0
    theta_opt = mean + m_cov @ z_opt
    min_sl = float(np.min(a_cone @ theta_opt))
    if min_sl <= 1e-6:
        rs = float(np.min(a_cone @ interior_ray))
        theta_opt = theta_opt + ((1e-5 - min_sl) / max(rs, 1e-12)) * interior_ray
    return theta_opt


# ---- 1D slice sampler & exact Huang--Wand Half-t_3 scale updates --------------------------------


def _slice_sample_1d(
    log_prob_fn: Any,
    x_init: float,
    rng: np.random.Generator,
    *,
    width: float = 1.5,
    max_steps: int = 12,
    bounds: tuple[float, float] = (-14.0, 14.0),
) -> float:
    """Exact Neal (2003) stepping-out + shrinkage 1-D slice sampler on ``bounds``."""
    x_new, _ = _slice_sample_1d_with_logp(
        log_prob_fn, x_init, rng, width=width, max_steps=max_steps, bounds=bounds
    )
    return x_new


def _slice_sample_1d_with_logp(
    log_prob_fn: Any,
    x_init: float,
    rng: np.random.Generator,
    *,
    width: float = 1.5,
    max_steps: int = 12,
    bounds: tuple[float, float] = (-14.0, 14.0),
    logp_init: float | None = None,
) -> tuple[float, float]:
    """Neal (2003) slice sampler returning ``(x_new, log_prob_fn(x_new))``.

    ``logp_init``, if given, must equal ``log_prob_fn(x_init)``; it is reused instead of
    re-evaluating the (possibly expensive) log density at the current state. The step-out
    width and step budget are fixed constants, so the transition is exactly invariant.
    """
    lo_b, hi_b = bounds
    x_f = float(x_init)
    x0 = lo_b + 1e-9 if x_f < lo_b + 1e-9 else (hi_b - 1e-9 if x_f > hi_b - 1e-9 else x_f)
    rand = rng.random
    fx0 = float(logp_init) if (logp_init is not None and x0 == x_f) else float(log_prob_fn(x0))
    log_y = fx0 - float(rng.standard_exponential())

    left_raw = x0 - float(rand()) * width
    right_raw = left_raw + width
    left = lo_b if left_raw < lo_b else left_raw
    right = hi_b if right_raw > hi_b else right_raw

    j_steps = int(float(rand()) * max_steps)
    k_steps = (max_steps - 1) - j_steps
    while j_steps > 0 and left > lo_b and float(log_prob_fn(left)) > log_y:
        left = max(lo_b, left - width)
        j_steps -= 1
    while k_steps > 0 and right < hi_b and float(log_prob_fn(right)) > log_y:
        right = min(hi_b, right + width)
        k_steps -= 1

    for _ in range(128):
        x_prop = left + float(rand()) * (right - left)
        f_prop = float(log_prob_fn(x_prop))
        if f_prop >= log_y:
            return x_prop, f_prop
        if x_prop < x0:
            left = x_prop
        else:
            right = x_prop
        if right - left < 1e-13:
            break
    return x0, fx0


def _sample_invgamma(shape: float, scale: float, rng: np.random.Generator) -> float:
    """Sample from ``InvGamma(shape, scale)`` with density ``propto x^{-shape - 1} exp(-scale / x)``."""
    g = float(rng.gamma(shape=max(shape, 1e-12), scale=1.0))
    return float(max(scale, 1e-30) / max(g, 1e-300))


def _log_halft_unit(diff: float, nu: float) -> float:
    """Exact log-density (up to constant in ``nu``) of ``diff = log(X / scale)`` for ``X ~ Half-t_nu(scale)``."""
    z = 2.0 * diff - math.log(nu)
    if z > 30.0:
        return diff - 0.5 * (nu + 1.0) * z
    return diff - 0.5 * (nu + 1.0) * math.log1p(math.exp(z))


def _log_scaled_f_logscale(w_val: float, alpha: float, nu: float, log_scale: float) -> float:
    """Exact normalized log-density of ``w = log(tau)`` where ``(tau / scale)^2 ~ F(alpha, nu)``."""
    u = w_val - log_scale
    log_ratio = math.log(alpha / nu)
    log_norm = (
        _LOG_2
        + 0.5 * alpha * log_ratio
        + math.lgamma(0.5 * (alpha + nu))
        - math.lgamma(0.5 * alpha)
        - math.lgamma(0.5 * nu)
    )
    z = 2.0 * u + log_ratio
    if z > 30.0:
        return log_norm + alpha * u - 0.5 * (alpha + nu) * z
    return log_norm + alpha * u - 0.5 * (alpha + nu) * math.log1p(math.exp(z))


def _log_outcome_tau_prior(w_tau: float) -> float:
    """Normalized log-density of ``w_tau = log(sigma_m / s_{0,m})`` under the per-outcome scale prior."""
    return _log_scaled_f_logscale(
        w_tau, _SCALE_PRIOR_ALPHA_TAU, _SCALE_PRIOR_NU_TAU, math.log(_SCALE_PRIOR_TAU0)
    )


def _sample_outcome_tau_prior(rng: np.random.Generator) -> float:
    """Exact forward sample of ``tau_m = sigma_m / s_{0,m}`` on ``exp((-w_bound, w_bound))``."""
    wt_b = _SCALE_PRIOR_W_TAU_BOUND
    while True:
        f_draw = float(rng.f(_SCALE_PRIOR_ALPHA_TAU, _SCALE_PRIOR_NU_TAU))
        tau = _SCALE_PRIOR_TAU0 * math.sqrt(max(f_draw, 1e-300))
        w_val = math.log(max(tau, 1e-300))
        if -wt_b < w_val < wt_b:
            return tau


def _update_scales_outcome(
    theta_m: np.ndarray,
    spec: _OutcomeConePriorSpec,
    sigma_sq: float,
    kappa_sq: float,
    a_sigma: float,
    a_kappa: float,
    local_scales: np.ndarray,
    b_local: np.ndarray,
    log_s0_eff: float,
    log_sk_eff: float,
    nu_sig: float,
    nu_kap: float,
    nu_loc: float,
    fixed_sigma: float | None,
    dynamic_shrinkage: bool,
    rng: np.random.Generator,
    coupled_p1: bool = False,
) -> tuple[float, float, float, float, np.ndarray, np.ndarray]:
    """Exact collapsed marginal scale updates for ``(sigma_m^2, kappa_m^2, lambda_m)``.

    Integrates out auxiliary variables analytically via Liu (1994) Rao--Blackwellization
    so centered scale updates have zero auxiliary-variable Gibbs lag while targeting the
    exact marginal posterior.
    """
    h_float = float(spec.horizon)
    n_nonctrl = spec.n_nonctrl
    k_dim = spec.k_dim
    theta_modes = theta_m.reshape(k_dim, n_nonctrl)

    # 1. Per-increment local scales lambda_{m,e} ~ Half-t_{nu_loc}(1) when active
    lam_cur = local_scales.copy()
    if spec.use_local_scales and spec.n_edges > 0:
        u_modes = theta_modes @ spec.c_cover.T  # (K', E_c)
        if spec.is_tree_poset:
            for e_idx in range(spec.n_edges):
                quad_e = float(u_modes[0, e_idx] ** 2) / max(sigma_sq, 1e-80)
                if k_dim > 1:
                    quad_e += float(np.sum(u_modes[1:, e_idx] ** 2)) / max(kappa_sq, 1e-80)
                c_quad = 0.5 * spec.c_scale * quad_e / h_float

                def _logp_we_tree(w_val: float, _cq: float = c_quad) -> float:
                    return (
                        _log_halft_unit(w_val, nu_loc)
                        - float(k_dim) * w_val
                        - _cq * math.exp(-2.0 * w_val)
                    )

                w_init = math.log(max(float(lam_cur[e_idx]), 1e-12))
                w_new = _slice_sample_1d(_logp_we_tree, w_init, rng, bounds=(-9.0, 9.0))
                lam_cur[e_idx] = math.exp(w_new)
        else:
            ratio_cur = math.sqrt(max(kappa_sq, 1e-80) / max(sigma_sq, 1e-80))
            for e_idx in rng.permutation(spec.n_edges)[:2]:
                quad_e = float(u_modes[0, e_idx] ** 2) / max(sigma_sq, 1e-80)
                if k_dim > 1:
                    quad_e += float(np.sum(u_modes[1:, e_idx] ** 2)) / max(kappa_sq, 1e-80)
                c_quad = 0.5 * spec.c_scale * quad_e / h_float

                def _logp_we_dag(w_val: float, _idx: int = int(e_idx), _cq: float = c_quad) -> float:
                    lam_try = lam_cur.copy()
                    lam_try[_idx] = math.exp(w_val)
                    _, logdet_l0 = spec.precision_and_logdet(lam_try)
                    log_z = spec.log_normalizer(ratio_cur, lam_try)
                    return (
                        _log_halft_unit(w_val, nu_loc)
                        - _cq * math.exp(-2.0 * w_val)
                        + 0.5 * float(k_dim) * logdet_l0
                        - log_z
                    )

                w_init = math.log(max(float(lam_cur[e_idx]), 1e-12))
                w_new = _slice_sample_1d(_logp_we_dag, w_init, rng, bounds=(-9.0, 9.0))
                lam_cur[e_idx] = math.exp(w_new)

        # Exact scale-ridge move (w_s + gamma, w_k + gamma, w_local - gamma) holding Sigma_0 invariant
        if fixed_sigma is None and (dynamic_shrinkage or k_dim == 1):
            ws0 = 0.5 * math.log(max(sigma_sq, 1e-80))
            wk0 = 0.5 * math.log(max(kappa_sq, 1e-80))
            w_loc0 = np.log(np.maximum(lam_cur, 1e-12))
            both_gk = bool(k_dim >= 2 and dynamic_shrinkage)
            if coupled_p1:
                wt_b_r = _SCALE_PRIOR_W_TAU_BOUND
                g_lo = max(log_s0_eff - wt_b_r - ws0, float(np.max(w_loc0 - 9.0)))
                g_hi = min(log_s0_eff + wt_b_r - ws0, float(np.min(w_loc0 + 9.0)))
                if g_hi - g_lo > 1e-9 and (g_lo < 0.0 < g_hi):
                    def _logp_gamma_p1(g_val: float) -> float:
                        lp = _log_outcome_tau_prior(ws0 + g_val - log_s0_eff)
                        for e_i in range(spec.n_edges):
                            lp += _log_halft_unit(float(w_loc0[e_i]) - g_val, nu_loc)
                        return lp

                    g_new = _slice_sample_1d(_logp_gamma_p1, 0.0, rng, width=2.5, bounds=(g_lo, g_hi))
                    sigma_sq = math.exp(2.0 * (ws0 + g_new))
                    if both_gk:
                        kappa_sq = math.exp(2.0 * (wk0 + g_new))
                    lam_cur = np.exp(w_loc0 - g_new)
            else:
                g_lo = max(-13.0 - ws0, float(np.max(w_loc0 - 9.0)))
                g_hi = min(13.0 - ws0, float(np.min(w_loc0 + 9.0)))
                if both_gk:
                    g_lo = max(g_lo, -13.0 - wk0)
                    g_hi = min(g_hi, 13.0 - wk0)
                if g_hi - g_lo > 1e-9 and (g_lo < 0.0 < g_hi):
                    def _logp_gamma(g_val: float) -> float:
                        lp = _log_halft_unit(ws0 + g_val - log_s0_eff, nu_sig)
                        if both_gk:
                            lp += _log_halft_unit(wk0 + g_val - log_sk_eff, nu_kap)
                        for e_i in range(spec.n_edges):
                            lp += _log_halft_unit(float(w_loc0[e_i]) - g_val, nu_loc)
                        return lp

                    g_new = _slice_sample_1d(_logp_gamma, 0.0, rng, bounds=(g_lo, g_hi))
                    sigma_sq = math.exp(2.0 * (ws0 + g_new))
                    if both_gk:
                        kappa_sq = math.exp(2.0 * (wk0 + g_new))
                    lam_cur = np.exp(w_loc0 - g_new)

    prec_l0, _ = spec.precision_and_logdet(lam_cur)
    q0 = float(theta_modes[0] @ prec_l0 @ theta_modes[0])
    q_dyn = float(sum(theta_modes[k] @ prec_l0 @ theta_modes[k] for k in range(1, k_dim))) if k_dim > 1 else 0.0
    c0 = 0.5 * q0 / h_float
    c_dyn = 0.5 * q_dyn / h_float
    d0 = float(n_nonctrl)
    d_dyn = float(n_nonctrl * (k_dim - 1))

    # 2. Global scale sigma_m^2 & dynamic scale kappa_m^2 via exact collapsed conditionals
    logz_depends_on_ratio = bool(spec.monotone and spec.period_positivity and k_dim >= 2 and spec.n_edges > 0)

    ws_cur = math.log(float(fixed_sigma)) if fixed_sigma is not None else 0.5 * math.log(max(sigma_sq, 1e-80))
    wk_cur = (
        (math.log(50.0) + log_sk_eff)
        if not dynamic_shrinkage
        else 0.5 * math.log(max(kappa_sq, 1e-80))
    )

    if coupled_p1:
        log_gamma0 = math.log(_SCALE_PRIOR_GAMMA0)
        nu_gamma = _SCALE_PRIOR_NU_GAMMA
        wt_b = _SCALE_PRIOR_W_TAU_BOUND
        ug_b = _SCALE_PRIOR_U_GAMMA_BOUND

        if fixed_sigma is None:
            ws_lo = log_s0_eff - wt_b
            ws_hi = log_s0_eff + wt_b
            if dynamic_shrinkage:
                ws_lo = max(ws_lo, wk_cur - ug_b)
                ws_hi = min(ws_hi, wk_cur + ug_b)
            if ws_hi - ws_lo > 1e-9 and (ws_lo < ws_cur < ws_hi):
                def _logp_ws_p1(ws: float) -> float:
                    lp = _log_outcome_tau_prior(ws - log_s0_eff) - d0 * ws - c0 * math.exp(-2.0 * ws)
                    if dynamic_shrinkage:
                        lp += _log_halft_unit(wk_cur - ws - log_gamma0, nu_gamma)
                    if logz_depends_on_ratio:
                        lp -= spec.log_normalizer(math.exp(wk_cur - ws), lam_cur)
                    return lp

                ws_cur = _slice_sample_1d(_logp_ws_p1, ws_cur, rng, width=2.5, bounds=(ws_lo, ws_hi))

        if dynamic_shrinkage:
            wk_lo = ws_cur - ug_b
            wk_hi = ws_cur + ug_b
            if k_dim == 1:
                def _logp_wk_prior_p1(wk: float) -> float:
                    return _log_halft_unit(wk - ws_cur - log_gamma0, nu_gamma)

                wk_cur = _slice_sample_1d(_logp_wk_prior_p1, wk_cur, rng, bounds=(wk_lo, wk_hi))
            else:
                def _logp_wk_p1(wk: float) -> float:
                    lp = (
                        _log_halft_unit(wk - ws_cur - log_gamma0, nu_gamma)
                        - d_dyn * wk
                        - c_dyn * math.exp(-2.0 * wk)
                    )
                    if logz_depends_on_ratio:
                        lp -= spec.log_normalizer(math.exp(wk - ws_cur), lam_cur)
                    return lp

                wk_cur = _slice_sample_1d(_logp_wk_p1, wk_cur, rng, bounds=(wk_lo, wk_hi))

        # Constant-ratio joint radial scale step (ws, wk) = (ws, ws + u) where gamma = exp(u) & log Z_m(gamma) stay fixed
        if fixed_sigma is None and dynamic_shrinkage:
            u_fixed = wk_cur - ws_cur
            ws_lo_u = log_s0_eff - wt_b
            ws_hi_u = log_s0_eff + wt_b
            if ws_hi_u - ws_lo_u > 1e-9 and (ws_lo_u < ws_cur < ws_hi_u):
                c_comb = c0 + c_dyn * math.exp(-2.0 * u_fixed)

                def _logp_ws_const_u_p1(ws: float) -> float:
                    return (
                        _log_outcome_tau_prior(ws - log_s0_eff)
                        - (d0 + d_dyn) * ws
                        - c_comb * math.exp(-2.0 * ws)
                    )

                ws_cur = _slice_sample_1d(_logp_ws_const_u_p1, ws_cur, rng, width=2.5, bounds=(ws_lo_u, ws_hi_u))
                wk_cur = ws_cur + u_fixed
    else:
        if fixed_sigma is None:
            def _logp_ws(ws: float) -> float:
                lp = _log_halft_unit(ws - log_s0_eff, nu_sig) - d0 * ws - c0 * math.exp(-2.0 * ws)
                if logz_depends_on_ratio:
                    lp -= spec.log_normalizer(math.exp(wk_cur - ws), lam_cur)
                return lp

            ws_cur = _slice_sample_1d(_logp_ws, ws_cur, rng, bounds=(-13.0, 13.0))

        if dynamic_shrinkage:
            if k_dim == 1:
                def _logp_wk_prior(wk: float) -> float:
                    return _log_halft_unit(wk - log_sk_eff, nu_kap)

                wk_cur = _slice_sample_1d(_logp_wk_prior, wk_cur, rng, bounds=(-13.0, 13.0))
            else:
                def _logp_wk(wk: float) -> float:
                    lp = _log_halft_unit(wk - log_sk_eff, nu_kap) - d_dyn * wk - c_dyn * math.exp(-2.0 * wk)
                    if logz_depends_on_ratio:
                        lp -= spec.log_normalizer(math.exp(wk - ws_cur), lam_cur)
                    return lp

                wk_cur = _slice_sample_1d(_logp_wk, wk_cur, rng, bounds=(-13.0, 13.0))

        # Constant-ratio joint radial scale step (ws, wk) = (ws, ws + u) where log Z_m(exp(u)) is constant
        if fixed_sigma is None and dynamic_shrinkage and k_dim >= 2:
            u_fixed = wk_cur - ws_cur
            ws_lo_u = max(-13.0, -13.0 - u_fixed)
            ws_hi_u = min(13.0, 13.0 - u_fixed)
            if ws_hi_u - ws_lo_u > 1e-9 and (ws_lo_u < ws_cur < ws_hi_u):
                c_comb = c0 + c_dyn * math.exp(-2.0 * u_fixed)

                def _logp_ws_const_u(ws: float) -> float:
                    return (
                        _log_halft_unit(ws - log_s0_eff, nu_sig)
                        + _log_halft_unit(ws + u_fixed - log_sk_eff, nu_kap)
                        - (d0 + d_dyn) * ws
                        - c_comb * math.exp(-2.0 * ws)
                    )

                ws_cur = _slice_sample_1d(_logp_ws_const_u, ws_cur, rng, bounds=(ws_lo_u, ws_hi_u))
                wk_cur = ws_cur + u_fixed

    new_sigma_sq = math.exp(2.0 * ws_cur)
    new_kappa_sq = math.exp(2.0 * wk_cur)
    return new_sigma_sq, new_kappa_sq, a_sigma, a_kappa, lam_cur, b_local


# ---- Collapsed SUR covariance, shared-arm noise draw, omega & R updates -------------------------


def _build_sigma_xi(
    sqrt_v_arms: list[list[np.ndarray]],
    corr_r: np.ndarray,
    n_arms: int,
    n_outcomes: int,
    k_dim: int,
    *,
    d_arms_stack: np.ndarray | None = None,
    rho_vec: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build block-diagonal ``Sigma_xi(R, rho)`` and its generative factor ``L_xi`` of shape ``(A*M*K', A*M*K')``."""
    chol_r = _fast_cholesky(corr_r, floor=1e-10)
    if rho_vec is not None:
        chol_r = np.sqrt(np.maximum(rho_vec, 1e-12))[:, None] * chol_r
    rk_factor = np.kron(chol_r, np.eye(k_dim, dtype=np.float64))  # (M*K', M*K')
    amk = n_arms * n_outcomes * k_dim
    sigma_xi = np.zeros((amk, amk), dtype=np.float64)
    l_xi = np.zeros((amk, amk), dtype=np.float64)
    mk = n_outcomes * k_dim
    if d_arms_stack is not None:
        la_stack = d_arms_stack @ rk_factor
        s_stack = la_stack @ np.swapaxes(la_stack, 1, 2)
        a_idx = np.arange(n_arms)
        l_xi.reshape(n_arms, mk, n_arms, mk)[a_idx, :, a_idx, :] = la_stack
        sigma_xi.reshape(n_arms, mk, n_arms, mk)[a_idx, :, a_idx, :] = s_stack
        return sigma_xi, l_xi
    for a in range(n_arms):
        d_a = np.zeros((mk, mk), dtype=np.float64)
        for m in range(n_outcomes):
            d_a[m * k_dim : (m + 1) * k_dim, m * k_dim : (m + 1) * k_dim] = sqrt_v_arms[m][a]
        la = d_a @ rk_factor
        sl = slice(a * mk, (a + 1) * mk)
        l_xi[sl, sl] = la
        sigma_xi[sl, sl] = la @ la.T
    return sigma_xi, l_xi


def _build_omega_cov(
    g_xi: np.ndarray,
    sigma_xi: np.ndarray,
    f_blocks: Sequence[np.ndarray] | np.ndarray,
    omega_vec: np.ndarray,
    k_dim: int | None = None,
) -> np.ndarray:
    """Collapsed measurement covariance ``Omega(omega, R) = G Sigma_xi(R) G^T + blkdiag(F_q / omega_q)``."""
    omega_mat = g_xi @ sigma_xi @ g_xi.T
    n_q = len(f_blocks)
    if (
        isinstance(f_blocks, np.ndarray)
        and k_dim is not None
        and omega_mat.shape[0] == n_q * k_dim
    ):
        q_idx = np.arange(n_q)
        omega_4d = omega_mat.reshape(n_q, k_dim, n_q, k_dim)
        omega_4d[q_idx, :, q_idx, :] += f_blocks / np.maximum(omega_vec, 1e-12)[:, None, None]
        return omega_mat
    off = 0
    for q in range(n_q):
        fq = f_blocks[q]
        kq = int(fq.shape[0])
        omega_mat[off : off + kq, off : off + kq] += fq / max(float(omega_vec[q]), 1e-12)
        off += kq
    return omega_mat


def _sample_omega_weights(
    eta_blocks: np.ndarray,
    f_inv_blocks: Sequence[np.ndarray] | np.ndarray,
    robust_error: str,
    pi_mix: float,
    rng: np.random.Generator,
    *,
    k_dims: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Update contrast precision weights ``omega_q`` and corruption indicators ``z_q``."""
    n_q, k_dim = eta_blocks.shape
    if robust_error == "gaussian":
        return np.ones(n_q, dtype=np.float64), np.zeros(n_q, dtype=np.float64), 0.0

    k_dims_arr = (
        np.asarray(k_dims, dtype=np.float64)
        if k_dims is not None
        else np.full(n_q, float(k_dim), dtype=np.float64)
    )
    f_inv_arr = f_inv_blocks if isinstance(f_inv_blocks, np.ndarray) else np.stack(f_inv_blocks, axis=0)
    d_sq = np.einsum("qi,qij,qj->q", eta_blocks, f_inv_arr, eta_blocks)
    if robust_error == "t4":
        # Student-t_4: omega_q ~ Gamma(2, 2) => omega_q | eta_q ~ Gamma(2 + K_q/2, 2 + d_q^2 / 2)
        shape_post = 2.0 + 0.5 * k_dims_arr
        rates = 2.0 + 0.5 * np.maximum(d_sq, 0.0)
        omega_out = np.maximum(rng.gamma(shape=shape_post, scale=1.0 / rates), 1e-10)
        flag_out = (omega_out < 0.40).astype(np.float64)
        return omega_out, flag_out, pi_mix

    # R2 ("mixture"): z_q ~ Bern(pi), log(1/omega_q) | z_q=1 ~ N(mu_infl, sd_infl^2), exact GH16 collapse
    s_nodes = _R2_LOG_INFL_MU + _R2_LOG_INFL_SD * _GH16_Z  # (16,)
    log_prior_odds = math.log(max(pi_mix, 1e-12) / max(1.0 - pi_mix, 1e-12))
    omega_out = np.ones(n_q, dtype=np.float64)
    z_out = np.zeros(n_q, dtype=np.float64)
    for q in range(n_q):
        dq2 = max(float(d_sq[q]), 0.0)
        kq_f = float(k_dims_arr[q])
        log_p0 = -0.5 * dq2
        log_terms = np.log(_GH16_PROB) - 0.5 * kq_f * s_nodes - 0.5 * dq2 * np.exp(-s_nodes)
        log_p1 = float(sp_special.logsumexp(log_terms))
        logit_q = float(np.clip(log_prior_odds + log_p1 - log_p0, -60.0, 60.0))
        prob_z1 = 1.0 / (1.0 + math.exp(-logit_q))
        if float(rng.uniform()) < prob_z1:
            z_out[q] = 1.0

            def _logp_s(s_val: float, _dq2: float = dq2, _kqf: float = kq_f) -> float:
                return (
                    -0.5 * _kqf * s_val
                    - 0.5 * _dq2 * math.exp(-s_val)
                    - 0.5 * ((s_val - _R2_LOG_INFL_MU) / _R2_LOG_INFL_SD) ** 2
                )

            s_draw = _slice_sample_1d(_logp_s, _R2_LOG_INFL_MU, rng, bounds=(-4.0, 14.0))
            omega_out[q] = math.exp(-s_draw)
        else:
            z_out[q] = 0.0
            omega_out[q] = 1.0

    n_z1 = float(np.sum(z_out))
    pi_new = float(rng.beta(_R2_BETA_ALPHA + n_z1, _R2_BETA_BETA + float(n_q) - n_z1))
    return omega_out, z_out, pi_new


def _sample_iw_correlation(
    df_iw: float,
    scale_iw: np.ndarray,
    rng: np.random.Generator,
    *,
    chol_inv: np.ndarray | None = None,
) -> np.ndarray:
    """Draw ``Sigma ~ InvWishart(df_iw, scale_iw)`` and rescale to correlation matrix ``R``."""
    m_dim = scale_iw.shape[0]
    if m_dim == 1:
        return np.ones((1, 1), dtype=np.float64)
    # Sample W ~ Wishart(df_iw, scale_iw^{-1}) via Bartlett decomposition
    if chol_inv is None:
        inv_scale = np.linalg.inv(_psd_clip(scale_iw, floor=1e-12))
        chol_inv = _fast_cholesky(inv_scale, floor=1e-12)
    a_bart = np.zeros((m_dim, m_dim), dtype=np.float64)
    for i in range(m_dim):
        a_bart[i, i] = math.sqrt(max(float(rng.chisquare(max(df_iw - float(i), 1.0))), 1e-15))
        for j in range(i):
            a_bart[i, j] = float(rng.standard_normal())
    la = chol_inv @ a_bart
    w_mat = la @ la.T
    try:
        sig_mat = np.linalg.inv(w_mat)
    except np.linalg.LinAlgError:
        sig_mat = np.linalg.inv(_psd_clip(w_mat, floor=1e-14))
    d_inv = 1.0 / np.sqrt(np.maximum(np.diag(sig_mat), 1e-24))
    r_mat = sig_mat * (d_inv[:, None] * d_inv[None, :])
    r_mat = 0.5 * (r_mat + r_mat.T)
    np.fill_diagonal(r_mat, 1.0)
    return r_mat


def _precompute_iw2_quadrature(df_iw: float) -> tuple[np.ndarray, np.ndarray]:
    """Precompute 64-point Gauss-Legendre log-sine terms and weights for ``IW_2(df_iw, Psi)``."""
    nu_eff = max(float(df_iw), 1.0 + 1e-9)
    return (nu_eff - 1.0) * _GL64_LOG_SIN, _GL64_PSI_W


def _log_iw2_marginal_z(
    z_val: float,
    df_iw: float,
    r_iw: float,
    iw2_log_sin_pow: np.ndarray,
    iw2_psi_w: np.ndarray,
) -> float:
    """Exact log-density (up to additive constant) of ``z = atanh(r)`` under ``IW_2(df_iw, Psi)``.

    Integrates out ``(sigma_1, sigma_2)`` analytically to obtain
    ``p(z) propto (1 - r^2)^{(nu - 1)/2} int_0^{pi/2} (sin psi)^{nu - 1} (1 - a sin psi)^{-nu} d psi``
    with ``a = r_iw * tanh(z)``. Factoring out the exact peak ``-nu * log(1 - a)`` at ``psi = pi/2``
    bounds every Gauss-Legendre summand in ``(0, 1]``, preventing overflow/underflow for all ``nu``.
    """
    r_val = math.tanh(z_val)
    a_val = r_val * r_iw
    az = abs(z_val)
    log_1mr2 = -2.0 * (az + math.log1p(math.exp(-2.0 * az)) - _LOG_2)
    log_1ma = math.log1p(-a_val)
    log_terms = iw2_log_sin_pow - df_iw * (np.log1p(-a_val * _GL64_SIN) - log_1ma)
    quad = float(np.dot(iw2_psi_w, np.exp(log_terms)))
    return 0.5 * (df_iw - 1.0) * log_1mr2 - df_iw * log_1ma + math.log(max(quad, 1e-300))


# ---- Posterior predictive checks (cycle, composite, basis) ---------------------------------------


def _compute_basis_ppc_pvalue(
    traj_mean_by_pair: Mapping[tuple[str, str], np.ndarray],
    traj_cov_by_pair: Mapping[tuple[str, str], np.ndarray],
    psi: np.ndarray,
    horizon: int,
) -> float | None:
    """Chi-square posterior predictive p-value for trajectory energy orthogonal to ``col(Psi)``.

    Returns ``None`` when ``col(Psi)`` spans ``R^H`` or when all Stage-1 trajectory
    draws already lie strictly inside ``col(Psi)`` (``total_df == 0``).
    """
    k_dim = int(psi.shape[1])
    if k_dim >= horizon:
        return None
    p_orth = np.eye(horizon, dtype=np.float64) - psi @ psi.T
    u_svd, s_svd, _ = np.linalg.svd(p_orth, full_matrices=False)
    q_perp = u_svd[:, s_svd > 0.5]  # (H, H - K')
    if int(q_perp.shape[1]) == 0:
        return None

    stat = 0.0
    total_df = 0
    for pair, mean_h in traj_mean_by_pair.items():
        cov_h = traj_cov_by_pair[pair]
        r_perp = q_perp.T @ mean_h
        w_perp_raw = q_perp.T @ cov_h @ q_perp
        w_perp_raw = 0.5 * (w_perp_raw + w_perp_raw.T)
        evals, evecs = np.linalg.eigh(w_perp_raw)
        eig_thresh = 1e-9 * max(float(np.trace(cov_h)), 1e-12)
        keep = evals > eig_thresh
        n_keep = int(np.count_nonzero(keep))
        if n_keep == 0:
            continue
        proj_r = evecs[:, keep].T @ r_perp
        stat += float(np.sum((proj_r * proj_r) / evals[keep]))
        total_df += n_keep
    if total_df == 0:
        return None
    return float(sp_stats.chi2.sf(max(stat, 0.0), df=total_df))


def _compute_cycle_ppc_pvalue(
    y_by_pair: Mapping[tuple[str, str], np.ndarray],
    f_by_pair: Mapping[tuple[str, str], np.ndarray],
    omega_draws_pair: Mapping[tuple[str, str], np.ndarray],
    fitted_pairs: tuple[tuple[str, str], ...],
    n_cycle: np.ndarray,
    k_dim: int,
    rng: np.random.Generator,
    *,
    active_modes_by_pair: Mapping[tuple[str, str], tuple[int, ...]] | None = None,
) -> float | None:
    """Posterior predictive p-value on cycle discrepancy ``(N_cycle^T otimes I_{K_c}) y_hat``."""
    n_cyc = int(n_cycle.shape[1])
    if n_cyc == 0:
        return None
    p_len = len(fitted_pairs)
    if active_modes_by_pair is None:
        common_modes = tuple(range(k_dim))
    else:
        cycle_pairs = [
            fitted_pairs[p_i]
            for p_i in range(p_len)
            if np.any(np.abs(n_cycle[p_i, :]) > 1e-12)
        ]
        if not cycle_pairs:
            return None
        common_modes = tuple(sorted(set.intersection(*(set(active_modes_by_pair[p]) for p in cycle_pairs))))
        if not common_modes:
            return None

    k_c = len(common_modes)
    y_parts: list[np.ndarray] = []
    f_parts: list[np.ndarray] = []
    for p_i, pair in enumerate(fitted_pairs):
        if active_modes_by_pair is None:
            y_parts.append(y_by_pair[pair])
            f_parts.append(f_by_pair[pair])
        else:
            act_p = active_modes_by_pair[pair]
            if set(common_modes).issubset(act_p):
                loc_idx = [act_p.index(k) for k in common_modes]
                y_parts.append(y_by_pair[pair][loc_idx])
                f_parts.append(f_by_pair[pair][np.ix_(loc_idx, loc_idx)])
            else:
                # Pair does not participate in any cycle (n_cycle[p_i, :] == 0)
                y_parts.append(np.zeros(k_c, dtype=np.float64))
                f_parts.append(1e-12 * np.eye(k_c, dtype=np.float64))

    y_stack = np.concatenate(y_parts, axis=0)  # (P * K_c,)
    nk_op = np.kron(n_cycle.T, np.eye(k_c, dtype=np.float64))  # (C * K_c, P * K_c)
    d_obs = nk_op @ y_stack

    m0_cov = np.zeros((n_cyc * k_c, n_cyc * k_c), dtype=np.float64)
    chol_f_list = []
    for p_i in range(p_len):
        fp = f_parts[p_i]
        chol_f_list.append(np.linalg.cholesky(_psd_clip(fp, floor=1e-14)))
        np_row = n_cycle[p_i : p_i + 1, :]  # (1, C)
        m0_cov += np.kron(np_row.T @ np_row, fp)
    m0_inv = np.linalg.inv(_psd_clip(m0_cov, floor=1e-13))
    t_obs = float(d_obs @ m0_inv @ d_obs)

    n_draws = next(iter(omega_draws_pair.values())).shape[0]
    eta_rep = np.zeros((n_draws, p_len * k_c), dtype=np.float64)
    for p_i, pair in enumerate(fitted_pairs):
        w_s = np.maximum(omega_draws_pair[pair], 1e-10)
        z_rep = rng.standard_normal((n_draws, k_c)) @ chol_f_list[p_i].T
        eta_rep[:, p_i * k_c : (p_i + 1) * k_c] = z_rep / np.sqrt(w_s)[:, None]
    d_rep = eta_rep @ nk_op.T  # (S, C * K_c)
    t_rep = np.einsum("si,ij,sj->s", d_rep, m0_inv, d_rep)
    return float(np.mean(t_rep >= t_obs))


def _compute_composite_ppc_pvalue(
    comp_y_by_pair: Mapping[tuple[str, str], np.ndarray],
    comp_f_by_pair: Mapping[tuple[str, str], np.ndarray],
    comp_omega_draws: Mapping[tuple[str, str], np.ndarray],
    struct_y_list: Sequence[Mapping[tuple[str, str], np.ndarray]],
    struct_f_list: Sequence[Mapping[tuple[str, str], np.ndarray]],
    struct_omega_list: Sequence[Mapping[tuple[str, str], np.ndarray]],
    comp_pairs: Sequence[tuple[str, str]],
    k_dim: int,
    rng: np.random.Generator,
    *,
    comp_active_modes: Mapping[tuple[str, str], tuple[int, ...]] | None = None,
    struct_active_modes_list: Sequence[Mapping[tuple[str, str], tuple[int, ...]]] | None = None,
) -> float | None:
    """Posterior predictive p-value for composite-vs-sum measurement discrepancy."""
    if not comp_pairs:
        return None
    n_draws = next(iter(comp_omega_draws.values())).shape[0]
    t_obs = 0.0
    t_rep = np.zeros(n_draws, dtype=np.float64)
    any_evaluated = False
    for pair in comp_pairs:
        if comp_active_modes is None or struct_active_modes_list is None:
            yc = comp_y_by_pair[pair]
            fc = comp_f_by_pair[pair]
            ys_list = [sy[pair] for sy in struct_y_list]
            fs_list = [sf[pair] for sf in struct_f_list]
            k_p = k_dim
        else:
            common_k = tuple(
                sorted(
                    set(comp_active_modes[pair]).intersection(
                        *(set(sm[pair]) for sm in struct_active_modes_list)
                    )
                )
            )
            if not common_k:
                continue
            c_idx = [comp_active_modes[pair].index(k) for k in common_k]
            yc = comp_y_by_pair[pair][c_idx]
            fc = comp_f_by_pair[pair][np.ix_(c_idx, c_idx)]
            ys_list = []
            fs_list = []
            for sy, sf, sm in zip(struct_y_list, struct_f_list, struct_active_modes_list, strict=True):
                s_idx = [sm[pair].index(k) for k in common_k]
                ys_list.append(sy[pair][s_idx])
                fs_list.append(sf[pair][np.ix_(s_idx, s_idx)])
            k_p = len(common_k)

        any_evaluated = True
        delta_obs = yc - sum(ys_list)
        m0_p = fc.copy()
        for sf_mat in fs_list:
            m0_p = m0_p + sf_mat
        m0_p_inv = np.linalg.inv(_psd_clip(m0_p, floor=1e-13))
        t_obs += float(delta_obs @ m0_p_inv @ delta_obs)

        chol_c = np.linalg.cholesky(_psd_clip(fc, floor=1e-14))
        wc = np.maximum(comp_omega_draws[pair], 1e-10)
        delta_rep = (rng.standard_normal((n_draws, k_p)) @ chol_c.T) / np.sqrt(wc)[:, None]
        for sf_mat, sw in zip(fs_list, struct_omega_list, strict=True):
            chol_m = np.linalg.cholesky(_psd_clip(sf_mat, floor=1e-14))
            wm = np.maximum(sw[pair], 1e-10)
            delta_rep = delta_rep - (rng.standard_normal((n_draws, k_p)) @ chol_m.T) / np.sqrt(wm)[:, None]
        t_rep += np.einsum("si,ij,sj->s", delta_rep, m0_p_inv, delta_rep)
    if not any_evaluated:
        return None
    return float(np.mean(t_rep >= t_obs))


# ---- Joint Bayesian multi-chain reconciliation core ---------------------------------------------


def _asis_interweave_outcome_scales(
    *,
    m_idx: int,
    theta: np.ndarray,
    w_res: np.ndarray,
    wx: np.ndarray,
    spec: _OutcomeConePriorSpec,
    sigma_sq_m: float,
    kappa_sq_m: float,
    a_sigma_m: float,
    a_kappa_m: float,
    local_scales_m: np.ndarray,
    log_s0_eff_m: float,
    log_sk_eff_m: float,
    nu_sig: float,
    nu_kap: float,
    nu_loc: float,
    fixed_sigma_m: float | None,
    dynamic_shrinkage: bool,
    rng: np.random.Generator,
    coupled_p1: bool = False,
) -> tuple[float, float, float, float, np.ndarray]:
    """Yu & Meng (2011) Ancillary Augmentation interweaving step for ``(sigma_m, kappa_m, lambda_m)``.

    Eliminates centered-Gibbs funnel slow-mixing when dynamic modes ``theta_{m,1..K'-1}``,
    the global level ``theta_{m,0}``, or individual tree increments concentrate near zero,
    while preserving exact detailed balance of the marginal cone prior + collapsed likelihood.
    """
    n_nonctrl = spec.n_nonctrl
    k_dim = spec.k_dim
    d_per_m = n_nonctrl * k_dim
    base_m = m_idx * d_per_m
    sl_m = slice(base_m, base_m + d_per_m)
    logz_dep = bool(spec.monotone and spec.period_positivity and k_dim >= 2 and spec.n_edges > 0)

    skip_expensive_ratio_aa = bool(spec.use_local_scales and not spec.is_tree_poset)
    log_gamma0 = math.log(_SCALE_PRIOR_GAMMA0)
    nu_gamma = _SCALE_PRIOR_NU_GAMMA
    wt_b = _SCALE_PRIOR_W_TAU_BOUND
    ug_b = _SCALE_PRIOR_U_GAMMA_BOUND

    # AA-0: Non-centered Mode-1 scale interweaving holding eta_0 = theta_0 / sigma_m fixed
    if fixed_sigma_m is None and k_dim >= 2 and dynamic_shrinkage and not skip_expensive_ratio_aa:
        sl_0 = slice(base_m, base_m + n_nonctrl)
        sl_dyn = slice(base_m + n_nonctrl, base_m + d_per_m)
        sig_0 = math.sqrt(max(sigma_sq_m, 1e-80))
        kap_0 = math.sqrt(max(kappa_sq_m, 1e-80))
        wk_0_log = math.log(kap_0)
        eta_0 = theta[sl_0] / sig_0
        if float(np.dot(eta_0, eta_0)) > 1e-50:
            if coupled_p1:
                ws_lo = max(log_s0_eff_m - wt_b, wk_0_log - ug_b)
                ws_hi = min(log_s0_eff_m + wt_b, wk_0_log + ug_b)
            else:
                ws_lo, ws_hi = -13.0, 13.0
            if spec.a_cone.shape[0] > 0:
                a0 = spec.a_cone_dyn @ theta[sl_dyn]
                a1 = spec.a_cone_0 @ eta_0
                neg_mask = a1 < -1e-15
                if np.any(neg_mask):
                    sig_max = float(np.min(-a0[neg_mask] / a1[neg_mask]))
                    ws_hi = min(ws_hi, math.log(max(sig_max, 1e-80)) - 1e-12)
                pos_neg = (a1 > 1e-15) & (a0 < 0.0)
                if np.any(pos_neg):
                    sig_min = float(np.max(-a0[pos_neg] / a1[pos_neg]))
                    ws_lo = max(ws_lo, math.log(max(sig_min, 1e-80)) + 1e-12)

            ws_cur = math.log(sig_0)
            if ws_hi - ws_lo > 1e-9 and (ws_lo < ws_cur < ws_hi):
                v_wx0 = wx[:, sl_0] @ eta_0
                r0_wy0 = w_res + sig_0 * v_wx0
                c1_0 = float(np.dot(r0_wy0, v_wx0))
                c2_0 = float(np.dot(v_wx0, v_wx0))

                def _logp_aa_ws(ws: float) -> float:
                    sig = math.exp(ws)
                    if coupled_p1:
                        lp = (
                            _log_outcome_tau_prior(ws - log_s0_eff_m)
                            + _log_halft_unit(wk_0_log - ws - log_gamma0, nu_gamma)
                            + sig * c1_0
                            - 0.5 * (sig * sig) * c2_0
                        )
                    else:
                        lp = _log_halft_unit(ws - log_s0_eff_m, nu_sig) + sig * c1_0 - 0.5 * (sig * sig) * c2_0
                    if logz_dep:
                        lp -= spec.log_normalizer(kap_0 / sig, local_scales_m)
                    return lp

                ws_new = _slice_sample_1d(_logp_aa_ws, ws_cur, rng, width=2.5, bounds=(ws_lo, ws_hi))
                sig_new = math.exp(ws_new)
                theta[sl_0] = sig_new * eta_0
                w_res = r0_wy0 - sig_new * v_wx0
                sigma_sq_m = sig_new * sig_new

    # AA-1: Non-centered dynamic-scale interweaving holding eta_dyn = theta_dyn / kappa_m fixed
    if dynamic_shrinkage and k_dim >= 2 and not skip_expensive_ratio_aa:
        sl_0 = slice(base_m, base_m + n_nonctrl)
        sl_dyn = slice(base_m + n_nonctrl, base_m + d_per_m)
        kap_0 = math.sqrt(max(kappa_sq_m, 1e-80))
        sig_0 = math.sqrt(max(sigma_sq_m, 1e-80))
        ws_0_log = math.log(sig_0)
        eta_dyn = theta[sl_dyn] / kap_0
        if float(np.dot(eta_dyn, eta_dyn)) > 1e-50:
            if coupled_p1:
                w_lo, w_hi = ws_0_log - ug_b, ws_0_log + ug_b
            else:
                w_lo, w_hi = -13.0, 13.0
            if spec.a_cone.shape[0] > 0:
                a0 = spec.a_cone_0 @ theta[sl_0]
                a1 = spec.a_cone_dyn @ eta_dyn
                neg_mask = a1 < -1e-15
                if np.any(neg_mask):
                    kap_max = float(np.min(-a0[neg_mask] / a1[neg_mask]))
                    w_hi = min(w_hi, math.log(max(kap_max, 1e-80)) - 1e-12)
                pos_neg = (a1 > 1e-15) & (a0 < 0.0)
                if np.any(pos_neg):
                    kap_min = float(np.max(-a0[pos_neg] / a1[pos_neg]))
                    w_lo = max(w_lo, math.log(max(kap_min, 1e-80)) + 1e-12)

            wk_cur = math.log(kap_0)
            if w_hi - w_lo > 1e-9 and (w_lo < wk_cur < w_hi):
                v_wx = wx[:, sl_dyn] @ eta_dyn
                r0_wy = w_res + kap_0 * v_wx
                c1 = float(np.dot(r0_wy, v_wx))
                c2 = float(np.dot(v_wx, v_wx))

                def _logp_aa_wk(wk: float) -> float:
                    kap = math.exp(wk)
                    if coupled_p1:
                        lp = _log_halft_unit(wk - ws_0_log - log_gamma0, nu_gamma) + kap * c1 - 0.5 * (kap * kap) * c2
                    else:
                        lp = _log_halft_unit(wk - log_sk_eff_m, nu_kap) + kap * c1 - 0.5 * (kap * kap) * c2
                    if logz_dep:
                        lp -= spec.log_normalizer(kap / sig_0, local_scales_m)
                    return lp

                wk_new = _slice_sample_1d(_logp_aa_wk, wk_cur, rng, bounds=(w_lo, w_hi))
                kap_new = math.exp(wk_new)
                theta[sl_dyn] = kap_new * eta_dyn
                w_res = r0_wy - kap_new * v_wx
                kappa_sq_m = kap_new * kap_new

    # AA-2: Joint radial cone-scale interweaving holding theta_m / sigma_m and kappa_m / sigma_m fixed
    if fixed_sigma_m is None and (dynamic_shrinkage or k_dim == 1 or not logz_dep):
        th_m_cur = theta[sl_m].copy()
        if float(np.dot(th_m_cur, th_m_cur)) > 1e-50:
            ws0 = 0.5 * math.log(max(sigma_sq_m, 1e-80))
            if k_dim >= 2 and not dynamic_shrinkage:
                # Only scale Mode-1 (theta_{m,0}, sigma_m) when log Z_m does not depend on ratio
                sl_0 = slice(base_m, base_m + n_nonctrl)
                th0_cur = theta[sl_0].copy()
                if coupled_p1:
                    d_lo = log_s0_eff_m - wt_b - ws0
                    d_hi = log_s0_eff_m + wt_b - ws0
                else:
                    d_lo = -13.0 - ws0
                    d_hi = 13.0 - ws0
                if float(np.dot(th0_cur, th0_cur)) > 1e-50 and d_hi - d_lo > 1e-9 and (d_lo < 0.0 < d_hi):
                    u_wx = wx[:, sl_0] @ th0_cur
                    r_base = w_res + u_wx
                    d1 = float(np.dot(r_base, u_wx))
                    d2 = float(np.dot(u_wx, u_wx))

                    def _logp_aa_d0(delta: float) -> float:
                        s = math.exp(delta)
                        if coupled_p1:
                            return _log_outcome_tau_prior(ws0 + delta - log_s0_eff_m) + s * d1 - 0.5 * (s * s) * d2
                        return _log_halft_unit(ws0 + delta - log_s0_eff_m, nu_sig) + s * d1 - 0.5 * (s * s) * d2

                    d_new = _slice_sample_1d(_logp_aa_d0, 0.0, rng, width=2.5, bounds=(d_lo, d_hi))
                    s_new = math.exp(d_new)
                    theta[sl_0] = s_new * th0_cur
                    w_res = r_base - s_new * u_wx
                    sigma_sq_m = math.exp(2.0 * (ws0 + d_new))
            else:
                # Scale entire outcome vector theta_m -> s * theta_m (cone-invariant & ratio-invariant!)
                wk0 = 0.5 * math.log(max(kappa_sq_m, 1e-80))
                if coupled_p1:
                    shift_wk = bool(dynamic_shrinkage)
                    d_lo = log_s0_eff_m - wt_b - ws0
                    d_hi = log_s0_eff_m + wt_b - ws0
                else:
                    shift_wk = bool(k_dim >= 2 and dynamic_shrinkage)
                    d_lo = max(-13.0 - ws0, -13.0 - wk0) if shift_wk else (-13.0 - ws0)
                    d_hi = min(13.0 - ws0, 13.0 - wk0) if shift_wk else (13.0 - ws0)
                if d_hi - d_lo > 1e-9 and (d_lo < 0.0 < d_hi):
                    u_wx = wx[:, sl_m] @ th_m_cur
                    r_base = w_res + u_wx
                    d1 = float(np.dot(r_base, u_wx))
                    d2 = float(np.dot(u_wx, u_wx))

                    def _logp_aa_joint(delta: float) -> float:
                        s = math.exp(delta)
                        if coupled_p1:
                            return _log_outcome_tau_prior(ws0 + delta - log_s0_eff_m) + s * d1 - 0.5 * (s * s) * d2
                        lp = _log_halft_unit(ws0 + delta - log_s0_eff_m, nu_sig) + s * d1 - 0.5 * (s * s) * d2
                        if shift_wk:
                            lp += _log_halft_unit(wk0 + delta - log_sk_eff_m, nu_kap)
                        return lp

                    d_new = _slice_sample_1d(_logp_aa_joint, 0.0, rng, width=2.5, bounds=(d_lo, d_hi))
                    s_new = math.exp(d_new)
                    theta[sl_m] = s_new * th_m_cur
                    w_res = r_base - s_new * u_wx
                    sigma_sq_m = math.exp(2.0 * (ws0 + d_new))
                    if shift_wk:
                        kappa_sq_m = math.exp(2.0 * (wk0 + d_new))

    # AA-3: Per-edge non-centered scale interweaving on tree posets when local scales are active
    if spec.use_local_scales and spec.is_tree_poset and spec.c_cover_inv is not None:
        c_inv = spec.c_cover_inv
        for e_idx in range(spec.n_edges):
            th_modes = theta[sl_m].reshape(k_dim, n_nonctrl)
            u_e = th_modes @ spec.c_cover[e_idx]  # (K',)
            if float(np.dot(u_e, u_e)) <= 1e-50:
                continue
            w_e0 = math.log(max(float(local_scales_m[e_idx]), 1e-12))
            de_lo = -9.0 - w_e0
            de_hi = 9.0 - w_e0
            if de_hi - de_lo <= 1e-9 or not (de_lo < 0.0 < de_hi):
                continue
            col_e = c_inv[:, e_idx]  # (A - 1,)
            d_th_e = np.concatenate([u_e[k] * col_e for k in range(k_dim)])
            u_wx_e = wx[:, sl_m] @ d_th_e
            r_base_e = w_res + u_wx_e
            d1_e = float(np.dot(r_base_e, u_wx_e))
            d2_e = float(np.dot(u_wx_e, u_wx_e))

            def _logp_aa_edge(
                delta: float, _we0: float = w_e0, _d1: float = d1_e, _d2: float = d2_e
            ) -> float:
                s = math.exp(delta)
                return _log_halft_unit(_we0 + delta, nu_loc) + s * _d1 - 0.5 * (s * s) * _d2

            de_new = _slice_sample_1d(_logp_aa_edge, 0.0, rng, bounds=(de_lo, de_hi))
            se_new = math.exp(de_new)
            theta[sl_m] = theta[sl_m] + (se_new - 1.0) * d_th_e
            w_res = r_base_e - se_new * u_wx_e
            local_scales_m[e_idx] = math.exp(w_e0 + de_new)

    return sigma_sq_m, kappa_sq_m, a_sigma_m, a_kappa_m, w_res


def _run_single_gibbs_chain(
    *,
    n_steps: int,
    thin: int = 1,
    state: dict[str, Any],
    record: bool,
    y_vec: np.ndarray,
    x_des: np.ndarray,
    g_xi: np.ndarray,
    f_blocks: list[np.ndarray],
    chol_f_blocks: list[np.ndarray],
    f_inv_blocks: list[np.ndarray],
    w_f_stack: np.ndarray,
    winv_f_stack: np.ndarray,
    d_ratio_stack: np.ndarray,
    q_member_weights: np.ndarray,
    sqrt_v_arms: list[list[np.ndarray]],
    prior_specs: list[_OutcomeConePriorSpec],
    a_cone_joint: np.ndarray,
    interior_ray_joint: np.ndarray,
    s0_scales: list[float],
    sk_scales: list[float],
    fixed_sigmas: list[float | None],
    n_arms: int,
    n_outcomes: int,
    k_dim: int,
    horizon: int,
    robust_error: str,
    dynamic_shrinkage: bool,
    sample_corr: bool,
    iw_df: float,
    iw_scale: np.ndarray,
    rng: np.random.Generator,
) -> tuple[dict[str, Any], dict[str, np.ndarray] | None]:
    """Execute ``n_steps * thin`` iterations of the collapsed blocked Gibbs + Pakman--Paninski HMC sampler."""
    n_nonctrl = n_arms - 1
    d_per_m = n_nonctrl * k_dim
    d_total = n_outcomes * d_per_m
    n_q = len(f_blocks)
    h_float = float(horizon)
    thin_eff = max(int(thin), 1)
    total_iters = int(n_steps) * thin_eff

    theta = state["theta"].copy()
    sigma_sq = state["sigma_sq"].copy()
    kappa_sq = state["kappa_sq"].copy()
    a_sigma = state["a_sigma"].copy()
    a_kappa = state["a_kappa"].copy()
    local_scales = [ls.copy() for ls in state["local_scales"]]
    b_local = [bl.copy() for bl in state["b_local"]]
    rho_vec = state["rho_vec"].copy()
    phi_vec = state["phi_vec"].copy()
    omega_vec = state["omega_vec"].copy()
    flag_vec = state["flag_vec"].copy()
    pi_mix = float(state["pi_mix"])
    corr_r = state["corr_r"].copy()

    log_s0_scales = [math.log(max(float(s0_scales[m]), 1e-16)) for m in range(n_outcomes)]
    log_sk_scales = [
        math.log(max(float(fixed_sigmas[m]) if fixed_sigmas[m] is not None else float(sk_scales[m]), 1e-16))
        for m in range(n_outcomes)
    ]

    tr_theta: np.ndarray | None = np.empty((n_steps, d_total), dtype=np.float64) if record else None
    tr_sigma: np.ndarray | None = np.empty((n_steps, n_outcomes), dtype=np.float64) if record else None
    tr_kappa: np.ndarray | None = np.empty((n_steps, n_outcomes), dtype=np.float64) if record else None
    tr_local: list[np.ndarray] = (
        [np.empty((n_steps, spec.n_edges), dtype=np.float64) for spec in prior_specs] if record else []
    )
    tr_rho: np.ndarray | None = np.empty((n_steps, n_outcomes), dtype=np.float64) if record else None
    tr_phi: np.ndarray | None = np.empty((n_steps, n_outcomes), dtype=np.float64) if record else None
    tr_omega: np.ndarray | None = np.empty((n_steps, n_q), dtype=np.float64) if record else None
    tr_flag: np.ndarray | None = np.empty((n_steps, n_q), dtype=np.float64) if record else None
    tr_corr: np.ndarray | None = (
        np.empty((n_steps, n_outcomes, n_outcomes), dtype=np.float64) if record else None
    )

    mk = n_outcomes * k_dim
    d_arms_stack = np.zeros((n_arms, mk, mk), dtype=np.float64)
    for a in range(n_arms):
        for m in range(n_outcomes):
            d_arms_stack[a, m * k_dim : (m + 1) * k_dim, m * k_dim : (m + 1) * k_dim] = sqrt_v_arms[m][a]

    k_q_list = [int(fq.shape[0]) for fq in f_blocks]
    k_q_vec = np.asarray(k_q_list, dtype=np.float64)
    all_full_rank = all(kq == k_dim for kq in k_q_list)
    n_meas = int(y_vec.shape[0])
    q_arange = np.arange(n_q)
    if all_full_rank:
        k_arange = np.arange(k_dim)
        blockdiag_flat_idx = (
            (q_arange[:, None, None] * k_dim + k_arange[None, :, None]) * n_meas
            + q_arange[:, None, None] * k_dim
            + k_arange[None, None, :]
        ).reshape(-1)
        active_f_flat: slice | np.ndarray = slice(None)
        active_vec_flat: slice | np.ndarray = slice(None)
    else:
        bd_parts: list[np.ndarray] = []
        act_f_mask = np.zeros((n_q, k_dim, k_dim), dtype=bool)
        act_v_mask = np.zeros((n_q, k_dim), dtype=bool)
        off_q = 0
        for q_i, kq in enumerate(k_q_list):
            rq = np.arange(kq)
            bd_parts.append(((off_q + rq[:, None]) * n_meas + (off_q + rq[None, :])).reshape(-1))
            act_f_mask[q_i, :kq, :kq] = True
            act_v_mask[q_i, :kq] = True
            off_q += kq
        blockdiag_flat_idx = np.concatenate(bd_parts, axis=0)
        active_f_flat = act_f_mask.reshape(-1)
        active_vec_flat = act_v_mask.reshape(-1)

    xy_des = np.column_stack([x_des, y_vec])
    g_xi_t = np.ascontiguousarray(g_xi.T)

    # Precompute constant bilinear basis B_{m1, m2} so G Sigma_xi(R, rho) G^T = sum_{m1, m2} sqrt(rho_{m1} rho_{m2}) R_{m1, m2} B_{m1, m2}
    xi_cov_basis = np.zeros((n_outcomes, n_outcomes, n_meas, n_meas), dtype=np.float64)
    for m1 in range(n_outcomes):
        for m2 in range(m1, n_outcomes):
            blk_mm = np.zeros((n_arms * mk, n_arms * mk), dtype=np.float64)
            for a in range(n_arms):
                r_sl = slice(a * mk + m1 * k_dim, a * mk + (m1 + 1) * k_dim)
                c_sl = slice(a * mk + m2 * k_dim, a * mk + (m2 + 1) * k_dim)
                blk_mm[r_sl, c_sl] = sqrt_v_arms[m1][a] @ sqrt_v_arms[m2][a].T
            b_12 = g_xi @ blk_mm @ g_xi_t
            xi_cov_basis[m1, m2] = b_12
            if m1 != m2:
                xi_cov_basis[m2, m1] = b_12.T

    xi_cov_basis_flat = xi_cov_basis.reshape(n_outcomes * n_outcomes, n_meas * n_meas)

    # Linear map from per-mode forest variance shares e_{q,k} / omega_q to the active block-diagonal
    # entries of blkdiag(F_q / omega_q), with F_q = W_q diag(e_q) W_q^T (rows ordered like
    # blockdiag_flat_idx), so Omega(rho, phi, omega, R) is assembled with two small GEMVs.
    f_design_5d = np.zeros((n_q, k_dim, k_dim, n_q, k_dim), dtype=np.float64)
    f_design_5d[q_arange, :, :, q_arange, :] = np.einsum("qik,qjk->qijk", w_f_stack, w_f_stack)
    f_design = np.ascontiguousarray(f_design_5d.reshape(n_q * k_dim * k_dim, n_q * k_dim)[active_f_flat, :])

    chol_iw_inv: np.ndarray | None = None
    r_iw_2 = 0.0
    iw2_t2: np.ndarray | None = None
    iw2_w_vt: np.ndarray | None = None
    if sample_corr and n_outcomes >= 2:
        inv_iw_scale = np.linalg.inv(_psd_clip(iw_scale, floor=1e-12))
        chol_iw_inv = _fast_cholesky(inv_iw_scale, floor=1e-12)
        if n_outcomes == 2:
            r_iw_2 = float(
                np.clip(
                    iw_scale[0, 1] / math.sqrt(max(iw_scale[0, 0] * iw_scale[1, 1], 1e-24)),
                    -0.99999999,
                    0.99999999,
                )
            )
            iw2_t2, iw2_w_vt = _precompute_iw2_quadrature(iw_df)

    def _compute_base_xi_cov(rv: np.ndarray, cr: np.ndarray) -> np.ndarray:
        if n_outcomes == 1:
            return float(rv[0]) * xi_cov_basis[0, 0]
        sr = np.sqrt(np.maximum(rv, 1e-12))
        rs = (sr[:, None] * cr) * sr[None, :]
        return (rs.reshape(-1) @ xi_cov_basis_flat).reshape(n_meas, n_meas)

    def _compute_f_stacks(rv: np.ndarray, pv: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        rho_q = q_member_weights @ rv
        phi_q = q_member_weights @ pv
        e_q = np.maximum(1.0 - rho_q[:, None] * d_ratio_stack, phi_q[:, None])
        sqrt_e_q = np.sqrt(e_q)
        cf_st = w_f_stack * sqrt_e_q[:, None, :]
        f_st = cf_st @ np.swapaxes(cf_st, 1, 2)
        winv_s = winv_f_stack / sqrt_e_q[:, :, None]
        finv_st = np.swapaxes(winv_s, 1, 2) @ winv_s
        return f_st, cf_st, finv_st

    def _make_omega_chol(base_cov: np.ndarray, f_st: np.ndarray, om_vec: np.ndarray) -> np.ndarray:
        om_mat = base_cov.copy()
        om_mat.reshape(-1)[blockdiag_flat_idx] += (
            f_st / np.maximum(om_vec, 1e-12)[:, None, None]
        ).reshape(-1)[active_f_flat]
        return _fast_cholesky(om_mat, floor=1e-13)

    sigma_xi_unit, l_xi_unit = _build_sigma_xi(
        sqrt_v_arms, corr_r, n_arms, n_outcomes, k_dim, d_arms_stack=d_arms_stack
    )
    if n_outcomes == 1:
        r0_f = max(float(rho_vec[0]), 1e-12)
        sigma_xi = r0_f * sigma_xi_unit
        l_xi = math.sqrt(r0_f) * l_xi_unit
    else:
        sigma_xi, l_xi = _build_sigma_xi(
            sqrt_v_arms, corr_r, n_arms, n_outcomes, k_dim, d_arms_stack=d_arms_stack, rho_vec=rho_vec
        )
    base_xi_cov = _compute_base_xi_cov(rho_vec, corr_r)
    f_stack, chol_f_stack, f_inv_stack = _compute_f_stacks(rho_vec, phi_vec)
    chol_omega = _make_omega_chol(base_xi_cov, f_stack, omega_vec)

    for it_idx in range(total_iters):
        # 1. Collapsed likelihood precision & mean marginalizing shared arm noise xi
        wxy = _solve_lower(chol_omega, xy_des)
        wx = wxy[:, :-1]
        wy = wxy[:, -1]
        prec_lik = wx.T @ wx
        h_lik = wx.T @ wy

        # 2. Prior precision across outcomes m = 0..M-1
        prec_prior = np.zeros((d_total, d_total), dtype=np.float64)
        for m in range(n_outcomes):
            l0_m, _ = prior_specs[m].precision_and_logdet(local_scales[m])
            base_m = m * d_per_m
            sl0 = slice(base_m, base_m + n_nonctrl)
            prec_prior[sl0, sl0] = l0_m / max(h_float * float(sigma_sq[m]), 1e-80)
            for k in range(1, k_dim):
                slk = slice(base_m + k * n_nonctrl, base_m + (k + 1) * n_nonctrl)
                prec_prior[slk, slk] = l0_m / max(h_float * float(kappa_sq[m]), 1e-80)

        prec_post = prec_lik + prec_prior
        chol_post = _fast_cholesky(prec_post, floor=1e-12)
        mean_post = _solve_lower(chol_post, _solve_lower(chol_post, h_lik), transpose=True)

        # 3. Exact Pakman--Paninski wall-reflecting HMC on cone K
        theta = _sample_truncated_mvn(
            mean_post,
            prec_post,
            a_cone_joint,
            theta,
            interior_ray_joint,
            rng,
            n_hmc_steps=1,
            chol_prec=chol_post,
        )

        # 4. Exact collapsed SA + Yu & Meng (2011) ASIS scale updates (independent per outcome)
        w_res = wy - wx @ theta
        for m in range(n_outcomes):
            spec_m = prior_specs[m]
            is_p1 = spec_m.scale_prior == "global_halft"
            log_s0_eff_m = log_s0_scales[m]
            log_sk_eff_m = log_sk_scales[m]
            nu_sig_m = 3.0
            nu_kap_m = _SCALE_PRIOR_NU_GAMMA if is_p1 else 3.0
            nu_loc_m = _SCALE_PRIOR_NU_LOCAL if is_p1 else 3.0

            theta_m = theta[m * d_per_m : (m + 1) * d_per_m]
            s_sq_m, k_sq_m, as_m, ak_m, lam_m, bl_m = _update_scales_outcome(
                theta_m=theta_m,
                spec=spec_m,
                sigma_sq=float(sigma_sq[m]),
                kappa_sq=float(kappa_sq[m]),
                a_sigma=float(a_sigma[m]),
                a_kappa=float(a_kappa[m]),
                local_scales=local_scales[m],
                b_local=b_local[m],
                log_s0_eff=log_s0_eff_m,
                log_sk_eff=log_sk_eff_m,
                nu_sig=nu_sig_m,
                nu_kap=nu_kap_m,
                nu_loc=nu_loc_m,
                fixed_sigma=fixed_sigmas[m],
                dynamic_shrinkage=dynamic_shrinkage,
                rng=rng,
                coupled_p1=is_p1,
            )
            s_sq_m, k_sq_m, as_m, ak_m, w_res = _asis_interweave_outcome_scales(
                m_idx=m,
                theta=theta,
                w_res=w_res,
                wx=wx,
                spec=spec_m,
                sigma_sq_m=s_sq_m,
                kappa_sq_m=k_sq_m,
                a_sigma_m=as_m,
                a_kappa_m=ak_m,
                local_scales_m=lam_m,
                log_s0_eff_m=log_s0_eff_m,
                log_sk_eff_m=log_sk_eff_m,
                nu_sig=nu_sig_m,
                nu_kap=nu_kap_m,
                nu_loc=nu_loc_m,
                fixed_sigma_m=fixed_sigmas[m],
                dynamic_shrinkage=dynamic_shrinkage,
                rng=rng,
                coupled_p1=is_p1,
            )
            sigma_sq[m] = s_sq_m
            kappa_sq[m] = k_sq_m
            a_sigma[m] = as_m
            a_kappa[m] = ak_m
            local_scales[m] = lam_m
            b_local[m] = bl_m

        res_vec = y_vec - x_des @ theta

        # 4b. Hierarchical Bayesian variance share updates (rho_m, phi_m) against collapsed marginal likelihood
        inv_om_vec = 1.0 / np.maximum(omega_vec, 1e-12)
        last_eval: dict[str, Any] = {"key": None, "chol": None}

        def _collapsed_loglik(rv: np.ndarray, eq_mat: np.ndarray, key: tuple[Any, ...]) -> float:
            """log N(res_vec; 0, Omega(rv, eq_mat, omega, R)) up to a constant; caches the factor."""
            sr = np.sqrt(np.maximum(rv, 1e-12))
            om_flat = ((sr[:, None] * corr_r) * sr[None, :]).reshape(-1) @ xi_cov_basis_flat
            om_flat[blockdiag_flat_idx] += f_design @ (eq_mat * inv_om_vec[:, None]).reshape(-1)
            ch = _fast_cholesky(om_flat.reshape(n_meas, n_meas), floor=1e-13)
            v_sol = _solve_lower(ch, res_vec)
            last_eval["key"] = key
            last_eval["chol"] = ch
            return -float(np.sum(np.log(np.diag(ch)))) - 0.5 * float(np.dot(v_sol, v_sol))

        # Collapsed log-likelihood at the current (rho, phi, omega, R): chol_omega is current here,
        # so it is reused instead of refactorizing Omega at the slice sampler's starting point.
        v_cur_4b = _solve_lower(chol_omega, res_vec)
        ll_cur = -float(np.sum(np.log(np.diag(chol_omega)))) - 0.5 * float(np.dot(v_cur_4b, v_cur_4b))

        def _log_beta_kernel(x_val: float, alpha: float, beta: float) -> float:
            return (alpha - 1.0) * math.log(max(x_val, 1e-300)) + (beta - 1.0) * math.log(
                max(1.0 - x_val, 1e-300)
            )

        for m in range(n_outcomes):
            phi_q_cur = q_member_weights @ phi_vec

            def _logp_rho(rho_try: float, _m: int = m) -> float:
                rv_try = rho_vec.copy()
                rv_try[_m] = rho_try
                rho_q_try = q_member_weights @ rv_try
                eq_try = np.maximum(1.0 - rho_q_try[:, None] * d_ratio_stack, phi_q_cur[:, None])
                lp_prior = _log_beta_kernel(rho_try, _RHO_BETA_ALPHA, _RHO_BETA_BETA)
                return lp_prior + _collapsed_loglik(rv_try, eq_try, ("rho", _m, rho_try))

            rho_old = float(rho_vec[m])
            rho_new, lp_rho_new = _slice_sample_1d_with_logp(
                _logp_rho,
                rho_old,
                rng,
                width=0.35,
                max_steps=3,
                bounds=(1e-5, 1.0 - 1e-5),
                logp_init=_log_beta_kernel(rho_old, _RHO_BETA_ALPHA, _RHO_BETA_BETA) + ll_cur,
            )
            rho_vec[m] = rho_new
            ll_cur = lp_rho_new - _log_beta_kernel(rho_new, _RHO_BETA_ALPHA, _RHO_BETA_BETA)

            rho_q_after = q_member_weights @ rho_vec
            min_unclipped = 1.0 - rho_q_after[:, None] * d_ratio_stack
            eq_base = np.maximum(min_unclipped, (q_member_weights @ phi_vec)[:, None])
            ll_base = ll_cur
            rho_vec_after = rho_vec.copy()

            def _logp_phi(phi_try: float, _m: int = m) -> float:
                pv_try = phi_vec.copy()
                pv_try[_m] = phi_try
                phi_q_try = q_member_weights @ pv_try
                eq_try = np.maximum(min_unclipped, phi_q_try[:, None])
                lp_prior = _log_beta_kernel(phi_try, _PHI_BETA_ALPHA, _PHI_BETA_BETA)
                if np.array_equal(eq_try, eq_base):
                    return lp_prior + ll_base
                return lp_prior + _collapsed_loglik(rho_vec_after, eq_try, ("phi", _m, phi_try))

            phi_old = float(phi_vec[m])
            phi_new, lp_phi_new = _slice_sample_1d_with_logp(
                _logp_phi,
                phi_old,
                rng,
                width=0.25,
                max_steps=3,
                bounds=(1e-5, 1.0 - 1e-5),
                logp_init=_log_beta_kernel(phi_old, _PHI_BETA_ALPHA, _PHI_BETA_BETA) + ll_cur,
            )
            phi_vec[m] = phi_new
            ll_cur = lp_phi_new - _log_beta_kernel(phi_new, _PHI_BETA_ALPHA, _PHI_BETA_BETA)

        f_stack, chol_f_stack, f_inv_stack = _compute_f_stacks(rho_vec, phi_vec)
        base_xi_cov = _compute_base_xi_cov(rho_vec, corr_r)
        # (sigma_xi, l_xi) are consumed only by Step 5 below, so they are rebuilt here once per
        # sweep from the current (rho, R) rather than after every R update in Step 6.
        if n_outcomes == 1:
            r0_f = max(float(rho_vec[0]), 1e-12)
            sigma_xi = r0_f * sigma_xi_unit
            l_xi = math.sqrt(r0_f) * l_xi_unit
        else:
            sigma_xi, l_xi = _build_sigma_xi(
                sqrt_v_arms, corr_r, n_arms, n_outcomes, k_dim, d_arms_stack=d_arms_stack, rho_vec=rho_vec
            )
        # The slice sampler's accepted point is the last density it evaluated, so the factor of
        # Omega at the final (rho, phi) is usually already available.
        if last_eval["key"] == ("phi", n_outcomes - 1, float(phi_vec[n_outcomes - 1])):
            chol_omega = last_eval["chol"]
        else:
            chol_omega = _make_omega_chol(base_xi_cov, f_stack, omega_vec)

        # 5. Exact conditional draw of shared arm noise xi | theta, omega, R, y via Matheron's identity
        if robust_error != "gaussian":
            for sub_o in range(2):
                if sub_o > 0:
                    chol_omega = _make_omega_chol(base_xi_cov, f_stack, omega_vec)
                xi_0 = l_xi @ rng.standard_normal(l_xi.shape[0])
                z_eta = rng.standard_normal((n_q, k_dim))
                eta_0_pad = (
                    np.einsum("qij,qj->qi", chol_f_stack, z_eta)
                    / np.sqrt(np.maximum(omega_vec, 1e-12))[:, None]
                )
                eta_0 = eta_0_pad.reshape(-1)[active_vec_flat]
                r_0 = g_xi @ xi_0 + eta_0
                delta_r = res_vec - r_0
                alpha_sol = _solve_lower(chol_omega, _solve_lower(chol_omega, delta_r), transpose=True)
                xi_draw = xi_0 + sigma_xi @ (g_xi_t @ alpha_sol)
                eta_vec = res_vec - g_xi @ xi_draw
                if all_full_rank:
                    eta_blocks = eta_vec.reshape(n_q, k_dim)
                else:
                    eta_blocks = np.zeros((n_q, k_dim), dtype=np.float64)
                    eta_blocks.reshape(-1)[active_vec_flat] = eta_vec
                omega_vec, flag_vec, pi_mix = _sample_omega_weights(
                    eta_blocks, f_inv_stack, robust_error, pi_mix, rng, k_dims=k_q_vec
                )
            chol_omega = _make_omega_chol(base_xi_cov, f_stack, omega_vec)

        # 6. Pre-period IW independence-MH + exact marginal IW_2 z-slice update for correlation R
        if sample_corr and n_outcomes >= 2:
            r_prop = _sample_iw_correlation(iw_df, iw_scale, rng, chol_inv=chol_iw_inv)
            base_xi_cov_prop = _compute_base_xi_cov(rho_vec, r_prop)
            chol_cur = chol_omega
            chol_prop = _make_omega_chol(base_xi_cov_prop, f_stack, omega_vec)

            v_cur = _solve_lower(chol_cur, res_vec)
            v_prop = _solve_lower(chol_prop, res_vec)
            loglik_cur = -float(np.sum(np.log(np.diag(chol_cur)))) - 0.5 * float(np.dot(v_cur, v_cur))
            loglik_prop = -float(np.sum(np.log(np.diag(chol_prop)))) - 0.5 * float(np.dot(v_prop, v_prop))
            if math.log(float(rng.uniform(1e-15, 1.0))) < (loglik_prop - loglik_cur):
                corr_r = r_prop
                base_xi_cov = base_xi_cov_prop
                chol_omega = chol_prop
                loglik_cur = loglik_prop

            if n_outcomes == 2 and iw2_t2 is not None and iw2_w_vt is not None:
                iw2_t2_arr = iw2_t2
                iw2_w_arr = iw2_w_vt
                r_cur = float(np.clip(corr_r[0, 1], -0.999995, 0.999995))
                z_cur = math.atanh(r_cur)
                chol_by_z: dict[float, np.ndarray] = {}
                # For M = 2, Omega(r) = Omega_fixed + r * B_off is affine in r = tanh(z).
                sr_z = np.sqrt(np.maximum(rho_vec, 1e-12))
                om_fixed_z = (sr_z[0] * sr_z[0]) * xi_cov_basis[0, 0] + (sr_z[1] * sr_z[1]) * xi_cov_basis[1, 1]
                om_fixed_z.reshape(-1)[blockdiag_flat_idx] += (
                    f_stack / np.maximum(omega_vec, 1e-12)[:, None, None]
                ).reshape(-1)[active_f_flat]
                b_off_z = (sr_z[0] * sr_z[1]) * (xi_cov_basis[0, 1] + xi_cov_basis[1, 0])

                def _logp_z_corr(z_val: float) -> float:
                    r_try = math.tanh(z_val)
                    ch_try = _fast_cholesky(om_fixed_z + r_try * b_off_z, floor=1e-13)
                    chol_by_z[z_val] = ch_try
                    v_try = _solve_lower(ch_try, res_vec)
                    ll_try = -float(np.sum(np.log(np.diag(ch_try)))) - 0.5 * float(np.dot(v_try, v_try))
                    return _log_iw2_marginal_z(z_val, iw_df, r_iw_2, iw2_t2_arr, iw2_w_arr) + ll_try

                z_new, _ = _slice_sample_1d_with_logp(
                    _logp_z_corr,
                    z_cur,
                    rng,
                    width=0.5,
                    max_steps=6,
                    bounds=(-6.5, 6.5),
                    logp_init=_log_iw2_marginal_z(z_cur, iw_df, r_iw_2, iw2_t2_arr, iw2_w_arr) + loglik_cur,
                )
                r_new = math.tanh(z_new)
                corr_r = np.array([[1.0, r_new], [r_new, 1.0]], dtype=np.float64)
                base_xi_cov = _compute_base_xi_cov(rho_vec, corr_r)
                cached_chol = chol_by_z.get(z_new)
                chol_omega = (
                    cached_chol
                    if cached_chol is not None
                    else _make_omega_chol(base_xi_cov, f_stack, omega_vec)
                )

        if (
            record
            and ((it_idx + 1) % thin_eff == 0)
            and tr_theta is not None
            and tr_sigma is not None
            and tr_kappa is not None
            and tr_rho is not None
            and tr_phi is not None
            and tr_omega is not None
            and tr_flag is not None
            and tr_corr is not None
        ):
            rec_idx = it_idx // thin_eff
            tr_theta[rec_idx] = theta
            tr_sigma[rec_idx] = np.sqrt(np.maximum(sigma_sq, 1e-80))
            tr_kappa[rec_idx] = np.sqrt(np.maximum(kappa_sq, 1e-80))
            for m in range(n_outcomes):
                tr_local[m][rec_idx] = local_scales[m]
            tr_rho[rec_idx] = rho_vec
            tr_phi[rec_idx] = phi_vec
            tr_omega[rec_idx] = omega_vec
            tr_flag[rec_idx] = flag_vec
            tr_corr[rec_idx] = corr_r

    new_state = {
        "theta": theta,
        "sigma_sq": sigma_sq,
        "kappa_sq": kappa_sq,
        "a_sigma": a_sigma,
        "a_kappa": a_kappa,
        "local_scales": local_scales,
        "b_local": b_local,
        "rho_vec": rho_vec,
        "phi_vec": phi_vec,
        "omega_vec": omega_vec,
        "flag_vec": flag_vec,
        "pi_mix": pi_mix,
        "corr_r": corr_r,
    }
    if not record:
        return new_state, None
    traces: dict[str, Any] = {
        "theta": tr_theta,
        "sigma": tr_sigma,
        "kappa": tr_kappa,
        "local": tr_local,
        "rho": tr_rho,
        "phi": tr_phi,
        "omega": tr_omega,
        "flag": tr_flag,
        "corr": tr_corr,
    }
    return new_state, traces


def _reconcile_joint_bayesian_core(
    *,
    structural_outcomes: tuple[str, ...],
    results_by_outcome: Mapping[str, Mapping[tuple[str, str], RankKLongBetResult]],
    composite_outcomes: Mapping[str, tuple[str, ...]],
    composite_results: Mapping[str, Mapping[tuple[str, str], RankKLongBetResult]],
    arms: Sequence[str] | None,
    order: Sequence[str] | None,
    poset: Sequence[tuple[str, str]] | None,
    monotone_by_outcome: Mapping[str, bool],
    prior_scale_by_outcome: Mapping[str, float | None],
    iw_spec: dict[str, Any] | None,
    horizon: int | None,
    n_modes: int,
    seed: int,
    model_options: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Unified Bayesian reconciliation engine for ``M >= 1`` outcomes + optional composites."""
    opts = _parse_model_options(model_options)
    k_req = int(n_modes)
    if k_req < 1:
        raise ValueError(f"n_modes must be >= 1; got {n_modes!r}")

    first_m = structural_outcomes[0]
    base_pairs = list(results_by_outcome[first_m].keys())
    for m_name in structural_outcomes[1:]:
        if set(results_by_outcome[m_name].keys()) != set(base_pairs):
            raise ValueError(
                f"All structural outcomes must share identical fitted contrast pairs; "
                f"{first_m!r} has {base_pairs} vs {m_name!r} has {list(results_by_outcome[m_name].keys())}"
            )

    topo = _resolve_arms_and_poset(base_pairs, arms=arms, order=order, poset=poset)
    arms_tuple: tuple[str, ...] = topo["arms"]
    arm_to_idx: dict[str, int] = topo["arm_to_idx"]
    fitted_pairs: tuple[tuple[str, str], ...] = topo["fitted_pairs"]
    n_arms = len(arms_tuple)
    n_nonctrl = n_arms - 1
    n_outcomes = len(structural_outcomes)

    # Infer or validate common window H across all fitted pairs & outcomes
    min_avail = min(
        _extract_raw_trt_trajectory(results_by_outcome[m][p]).shape[1]
        for m in structural_outcomes
        for p in fitted_pairs
    )
    for c_name, c_res_map in composite_results.items():
        for p, r_obj in c_res_map.items():
            min_avail = min(min_avail, _extract_raw_trt_trajectory(r_obj).shape[1])

    if horizon is None:
        h_eff = int(min_avail)
    else:
        h_eff = int(horizon)
        if h_eff < 1 or h_eff > min_avail:
            raise ValueError(
                f"horizon={horizon!r} must lie in [1, {min_avail}] (available post-adoption periods)"
            )

    def _has_stage1_post_cov(res_obj: Any) -> bool:
        meta_obj = getattr(res_obj, "arm_metadata", None)
        cov_trt_m = _lookup_meta_cov(meta_obj, "treated")
        cov_ctrl_m = _lookup_meta_cov(meta_obj, "control")
        return (
            cov_trt_m is not None
            and cov_ctrl_m is not None
            and cov_trt_m.shape[0] >= h_eff
            and cov_ctrl_m.shape[0] >= h_eff
        )

    missing_post_cov = any(
        not _has_stage1_post_cov(results_by_outcome[m][p])
        for m in structural_outcomes
        for p in fitted_pairs
    ) or any(
        not _has_stage1_post_cov(r_obj)
        for c_map in composite_results.values()
        for p, r_obj in c_map.items()
        if p in fitted_pairs
    )
    if topo["is_bipartite"] and n_arms >= 3 and topo["suggested_odd_pair"] is not None and missing_post_cov:
        sug_u, sug_v = topo["suggested_odd_pair"]
        warnings.warn(
            f"Contrast graph on arms {arms_tuple} is bipartite (no odd cycles) and one or more fitted pairs "
            f"lack Stage-1 arm_metadata post-period residual covariances. Pass Stage-1 arm_metadata or fit "
            f"({sug_u!r}, {sug_v!r}) to close an odd cycle.",
            RankKDesignWarning,
            stacklevel=3,
        )

    psi = _euclidean_exposure_basis(h_eff, k_req)
    k_dim = int(psi.shape[1])
    d_per_m = n_nonctrl * k_dim
    d_total = n_outcomes * d_per_m

    # Extract summaries for every structural outcome
    struct_summaries = {
        m: _extract_outcome_summaries(
            results_by_outcome[m], arms_tuple, arm_to_idx, fitted_pairs, psi, h_eff
        )
        for m in structural_outcomes
    }
    # Extract summaries for fitted composite outcomes (if composite_results supplied)
    comp_summaries: dict[str, dict[str, Any]] = {}
    for c_name, c_map in composite_results.items():
        c_pairs = tuple(p for p in fitted_pairs if p in c_map)
        if not c_pairs:
            continue
        comp_summaries[c_name] = _extract_outcome_summaries(
            c_map, arms_tuple, arm_to_idx, c_pairs, psi, h_eff
        )
        comp_summaries[c_name]["pairs"] = c_pairs

    # Assemble stacked measurement system (y_vec, X_des, G_xi, F_blocks) with partial-mode support
    y_blocks: list[np.ndarray] = []
    x_rows: list[np.ndarray] = []
    g_rows: list[np.ndarray] = []
    f_blocks: list[np.ndarray] = []
    w_f_list: list[np.ndarray] = []
    winv_f_list: list[np.ndarray] = []
    d_ratio_list: list[np.ndarray] = []
    member_weights_list: list[np.ndarray] = []
    meas_index: dict[tuple[str, str, tuple[str, str]], int] = {}

    amk = n_arms * n_outcomes * k_dim
    mk = n_outcomes * k_dim

    def _add_measurement(
        kind: str,
        out_name: str,
        pair: tuple[str, str],
        member_indices: Sequence[int],
        act_modes_q: tuple[int, ...],
        y_q: np.ndarray,
        f_q: np.ndarray,
        w_f_q: np.ndarray,
        winv_f_q: np.ndarray,
        d_ratio_q: np.ndarray,
    ) -> None:
        q_idx = len(f_blocks)
        meas_index[(kind, out_name, pair)] = q_idx
        kq = len(act_modes_q)
        y_blocks.append(y_q)
        f_blocks.append(f_q)
        w_pad = np.zeros((k_dim, k_dim), dtype=np.float64)
        winv_pad = np.zeros((k_dim, k_dim), dtype=np.float64)
        d_pad = np.zeros(k_dim, dtype=np.float64)
        w_pad[:kq, :kq] = w_f_q
        winv_pad[:kq, :kq] = winv_f_q
        d_pad[:kq] = d_ratio_q
        w_f_list.append(w_pad)
        winv_f_list.append(winv_pad)
        d_ratio_list.append(d_pad)
        mw = np.zeros(n_outcomes, dtype=np.float64)
        for m_i in member_indices:
            mw[m_i] = 1.0 / float(len(member_indices))
        member_weights_list.append(mw)

        ia, ib = arm_to_idx[pair[0]], arm_to_idx[pair[1]]
        x_q = np.zeros((kq, d_total), dtype=np.float64)
        g_q = np.zeros((kq, amk), dtype=np.float64)
        for m_i in member_indices:
            base_m = m_i * d_per_m
            for row_r, k in enumerate(act_modes_q):
                if ia > 0:
                    x_q[row_r, base_m + k * n_nonctrl + (ia - 1)] += 1.0
                if ib > 0:
                    x_q[row_r, base_m + k * n_nonctrl + (ib - 1)] -= 1.0
                g_q[row_r, ia * mk + m_i * k_dim + k] += 1.0
                g_q[row_r, ib * mk + m_i * k_dim + k] -= 1.0
        x_rows.append(x_q)
        g_rows.append(g_q)

    for m_i, m_name in enumerate(structural_outcomes):
        summ_m = struct_summaries[m_name]
        for pair in fitted_pairs:
            _add_measurement(
                "struct",
                m_name,
                pair,
                (m_i,),
                summ_m["active_modes_by_pair"][pair],
                summ_m["y_by_pair"][pair],
                summ_m["f_by_pair"][pair],
                summ_m["w_f_by_pair"][pair],
                summ_m["winv_f_by_pair"][pair],
                summ_m["d_ratio_by_pair"][pair],
            )

    m_name_to_idx = {m: idx for idx, m in enumerate(structural_outcomes)}
    for c_name, c_summ in comp_summaries.items():
        members = tuple(m_name_to_idx[m] for m in composite_outcomes[c_name])
        for pair in c_summ["pairs"]:
            _add_measurement(
                "comp",
                c_name,
                pair,
                members,
                c_summ["active_modes_by_pair"][pair],
                c_summ["y_by_pair"][pair],
                c_summ["f_by_pair"][pair],
                c_summ["w_f_by_pair"][pair],
                c_summ["winv_f_by_pair"][pair],
                c_summ["d_ratio_by_pair"][pair],
            )

    y_vec = np.concatenate(y_blocks, axis=0)
    x_des = np.vstack(x_rows)
    g_xi = np.vstack(g_rows)
    chol_f_blocks = [np.linalg.cholesky(_psd_clip(fq, floor=1e-14)) for fq in f_blocks]
    f_inv_blocks = [np.linalg.inv(_psd_clip(fq, floor=1e-13)) for fq in f_blocks]
    w_f_stack = np.stack(w_f_list, axis=0)
    winv_f_stack = np.stack(winv_f_list, axis=0)
    d_ratio_stack = np.stack(d_ratio_list, axis=0)
    q_member_weights = np.stack(member_weights_list, axis=0)
    sqrt_v_arms = [[_psd_sqrt(va) for va in struct_summaries[m]["v_arms"]] for m in structural_outcomes]

    # Build per-outcome cone prior specs and joint cone matrix
    prior_specs: list[_OutcomeConePriorSpec] = []
    s0_scales: list[float] = []
    sk_scales: list[float] = []
    fixed_sigmas: list[float | None] = []
    cone_blocks: list[np.ndarray] = []
    ray_parts: list[np.ndarray] = []

    for m_name in structural_outcomes:
        fixed_m = prior_scale_by_outcome[m_name]
        spec_m = _OutcomeConePriorSpec(
            c_cover=topo["c_cover"],
            psi=psi,
            horizon=h_eff,
            monotone=bool(monotone_by_outcome[m_name]),
            period_positivity=bool(opts["period_positivity"]),
            scale_prior=str(opts["scale_prior"]),
            use_local_scales=False if (opts["scale_prior"] == "global_halft" and fixed_m is not None) else None,
        )
        prior_specs.append(spec_m)
        s0_m = float(struct_summaries[m_name]["s0_scale"])
        sk_m = float(struct_summaries[m_name]["sk_scale"])
        s0_scales.append(s0_m)
        sk_scales.append(sk_m)
        fixed_sigmas.append(fixed_m)
        cone_blocks.append(spec_m.a_cone)

        ray_scale_m = fixed_m if fixed_m is not None else (_SCALE_PRIOR_TAU0 * s0_m)
        ray_m = np.zeros(d_per_m, dtype=np.float64)
        ray_m[:n_nonctrl] = math.sqrt(float(h_eff)) * ray_scale_m * topo["interior_heights"]
        ray_parts.append(ray_m)

    is_p1_global = bool(opts["scale_prior"] == "global_halft")

    total_walls = sum(cb.shape[0] for cb in cone_blocks)
    a_cone_joint = np.zeros((total_walls, d_total), dtype=np.float64)
    w_off = 0
    for m_i, cb in enumerate(cone_blocks):
        if cb.shape[0] > 0:
            a_cone_joint[w_off : w_off + cb.shape[0], m_i * d_per_m : (m_i + 1) * d_per_m] = cb
            w_off += cb.shape[0]
    interior_ray_joint = np.concatenate(ray_parts, axis=0)

    # Resolve cross-outcome correlation Inverse-Wishart parameters
    sample_corr = bool(opts["use_sur"] and n_outcomes >= 2 and iw_spec is not None)
    if sample_corr and iw_spec is not None:
        iw_df = float(iw_spec["df_post"])
        iw_scale = np.asarray(iw_spec["scale_post"], dtype=np.float64)
        r_init_mean = np.asarray(iw_spec["r_init"], dtype=np.float64)
    else:
        iw_df = float(n_outcomes + 1)
        iw_scale = np.eye(n_outcomes, dtype=np.float64)
        r_init_mean = np.eye(n_outcomes, dtype=np.float64)

    # Compute cone-constrained posterior mode at nominal hyperparameters for chain initialization
    rho_init_val = _RHO_BETA_ALPHA / (_RHO_BETA_ALPHA + _RHO_BETA_BETA)
    phi_init_val = _PHI_BETA_ALPHA / (_PHI_BETA_ALPHA + _PHI_BETA_BETA)
    sigma_xi_0, _ = _build_sigma_xi(
        sqrt_v_arms,
        r_init_mean,
        n_arms,
        n_outcomes,
        k_dim,
        rho_vec=np.full(n_outcomes, rho_init_val, dtype=np.float64),
    )
    omega_cov_0 = _build_omega_cov(g_xi, sigma_xi_0, f_blocks, np.ones(len(f_blocks)), k_dim)
    chol_w0 = np.linalg.cholesky(omega_cov_0)
    wx0 = np.linalg.solve(chol_w0, x_des)
    wy0 = np.linalg.solve(chol_w0, y_vec)
    prec_init = wx0.T @ wx0
    tau_init_nom = _SCALE_PRIOR_TAU0 if is_p1_global else 1.0
    for m_i in range(n_outcomes):
        l0_m, _ = prior_specs[m_i].precision_and_logdet(np.ones(prior_specs[m_i].n_edges))
        s_init = fixed_sigmas[m_i] if fixed_sigmas[m_i] is not None else (tau_init_nom * s0_scales[m_i])
        if not opts["dynamic_shrinkage"]:
            k_init = 50.0 * (fixed_sigmas[m_i] if fixed_sigmas[m_i] is not None else sk_scales[m_i])
        elif is_p1_global:
            k_init = _SCALE_PRIOR_GAMMA0 * s_init
        else:
            k_init = fixed_sigmas[m_i] if fixed_sigmas[m_i] is not None else (tau_init_nom * sk_scales[m_i])
        base_m = m_i * d_per_m
        prec_init[base_m : base_m + n_nonctrl, base_m : base_m + n_nonctrl] += l0_m / (float(h_eff) * s_init**2)
        for k in range(1, k_dim):
            slk = slice(base_m + k * n_nonctrl, base_m + (k + 1) * n_nonctrl)
            prec_init[slk, slk] += l0_m / (float(h_eff) * k_init**2)
    prec_init = _psd_clip(prec_init, floor=1e-12)
    mean_init = np.linalg.solve(prec_init, wx0.T @ wy0)
    theta_mode = _cone_qp_mode(mean_init, prec_init, a_cone_joint, interior_ray_joint)

    # Spawn 4 deterministic independent RNG streams from seed
    n_chains = int(opts["n_chains"])
    n_warmup = int(opts["n_warmup"])
    n_draws_chain = int(opts["n_draws_per_chain"])
    max_rounds = int(opts["max_continuation_rounds"])

    ss = np.random.SeedSequence(int(seed))
    chain_rngs = [np.random.default_rng(child) for child in ss.spawn(n_chains)]

    p1_tau_spread = (0.04, 0.25, 2.0, 12.0)
    chain_states: list[dict[str, Any]] = []
    for c_i in range(n_chains):
        rng_c = chain_rngs[c_i]
        jitter = float(rng_c.uniform(0.8, 1.25))
        tau_c_base = p1_tau_spread[c_i % len(p1_tau_spread)] if is_p1_global else tau_init_nom
        ray_mult_c = min(tau_c_base, 1.0)
        theta_c = theta_mode + 0.05 * ray_mult_c * jitter * interior_ray_joint
        if a_cone_joint.shape[0] > 0 and float(np.min(a_cone_joint @ theta_c)) <= 0.0:
            theta_c = theta_mode + 0.05 * ray_mult_c * interior_ray_joint
        sig_c_list = [
            float(fixed_sigmas[m])
            if fixed_sigmas[m] is not None
            else (
                tau_c_base * s0_scales[m] * float(rng_c.uniform(0.75, 1.35))
                if is_p1_global
                else (tau_init_nom * s0_scales[m] * float(rng_c.uniform(0.5, 1.8)))
            )
            for m in range(n_outcomes)
        ]
        sig_sq_c = np.array([sc * sc for sc in sig_c_list], dtype=np.float64)
        kap_sq_c = np.array(
            [
                ((50.0 * (fixed_sigmas[m] if fixed_sigmas[m] is not None else sk_scales[m])) ** 2)
                if not opts["dynamic_shrinkage"]
                else (
                    (sig_c_list[m] * _SCALE_PRIOR_GAMMA0 * float(rng_c.uniform(0.6, 1.6))) ** 2
                    if is_p1_global
                    else (
                        (
                            (fixed_sigmas[m] if fixed_sigmas[m] is not None else sk_scales[m])
                            * float(rng_c.uniform(0.6, 1.6))
                        )
                        ** 2
                    )
                )
                for m in range(n_outcomes)
            ],
            dtype=np.float64,
        )
        rho_c = np.clip(
            rng_c.beta(_RHO_BETA_ALPHA, _RHO_BETA_BETA, size=n_outcomes),
            0.05,
            0.95,
        )
        phi_c = np.clip(
            rng_c.beta(_PHI_BETA_ALPHA, _PHI_BETA_BETA, size=n_outcomes),
            0.02,
            0.85,
        )
        st_c = {
            "theta": theta_c,
            "sigma_sq": sig_sq_c,
            "kappa_sq": kap_sq_c,
            "a_sigma": np.array([s0_scales[m] ** 2 for m in range(n_outcomes)], dtype=np.float64),
            "a_kappa": np.array([sk_scales[m] ** 2 for m in range(n_outcomes)], dtype=np.float64),
            "local_scales": [np.ones(spec.n_edges, dtype=np.float64) for spec in prior_specs],
            "b_local": [np.ones(spec.n_edges, dtype=np.float64) for spec in prior_specs],
            "rho_vec": rho_c,
            "phi_vec": phi_c,
            "omega_vec": np.ones(len(f_blocks), dtype=np.float64),
            "flag_vec": np.zeros(len(f_blocks), dtype=np.float64),
            "pi_mix": 0.10,
            "corr_r": r_init_mean.copy(),
        }
        st_c, _ = _run_single_gibbs_chain(
            n_steps=n_warmup,
            state=st_c,
            record=False,
            y_vec=y_vec,
            x_des=x_des,
            g_xi=g_xi,
            f_blocks=f_blocks,
            chol_f_blocks=chol_f_blocks,
            f_inv_blocks=f_inv_blocks,
            w_f_stack=w_f_stack,
            winv_f_stack=winv_f_stack,
            d_ratio_stack=d_ratio_stack,
            q_member_weights=q_member_weights,
            sqrt_v_arms=sqrt_v_arms,
            prior_specs=prior_specs,
            a_cone_joint=a_cone_joint,
            interior_ray_joint=interior_ray_joint,
            s0_scales=s0_scales,
            sk_scales=sk_scales,
            fixed_sigmas=fixed_sigmas,
            n_arms=n_arms,
            n_outcomes=n_outcomes,
            k_dim=k_dim,
            horizon=h_eff,
            robust_error=str(opts["robust_error"]),
            dynamic_shrinkage=bool(opts["dynamic_shrinkage"]),
            sample_corr=sample_corr,
            iw_df=iw_df,
            iw_scale=iw_scale,
            rng=rng_c,
        )
        chain_states.append(st_c)

    # Post-warmup sampling + unconditional auto-continuation until Rhat <= 1.01 and ESS >= 400
    chain_traces: list[dict[str, Any]] = []
    for c_i in range(n_chains):
        chain_states[c_i], tr_c = _run_single_gibbs_chain(
            n_steps=n_draws_chain,
            state=chain_states[c_i],
            record=True,
            y_vec=y_vec,
            x_des=x_des,
            g_xi=g_xi,
            f_blocks=f_blocks,
            chol_f_blocks=chol_f_blocks,
            f_inv_blocks=f_inv_blocks,
            w_f_stack=w_f_stack,
            winv_f_stack=winv_f_stack,
            d_ratio_stack=d_ratio_stack,
            q_member_weights=q_member_weights,
            sqrt_v_arms=sqrt_v_arms,
            prior_specs=prior_specs,
            a_cone_joint=a_cone_joint,
            interior_ray_joint=interior_ray_joint,
            s0_scales=s0_scales,
            sk_scales=sk_scales,
            fixed_sigmas=fixed_sigmas,
            n_arms=n_arms,
            n_outcomes=n_outcomes,
            k_dim=k_dim,
            horizon=h_eff,
            robust_error=str(opts["robust_error"]),
            dynamic_shrinkage=bool(opts["dynamic_shrinkage"]),
            sample_corr=sample_corr,
            iw_df=iw_df,
            iw_scale=iw_scale,
            rng=chain_rngs[c_i],
        )
        assert tr_c is not None
        chain_traces.append(tr_c)

    def _evaluate_convergence(traces_list: list[dict[str, Any]]) -> tuple[bool, float, float, float, dict[str, float], dict[str, float]]:
        rhat_blocks: dict[str, float] = {}
        ess_blocks: dict[str, float] = {}

        theta_stack = np.stack([t["theta"] for t in traces_list], axis=0)  # (C, N, D)
        th_stats = [_split_rhat_and_ess_1d(theta_stack[:, :, d]) for d in range(d_total)]
        rhat_blocks["theta"] = max(s[0] for s in th_stats)
        ess_blocks["theta"] = min(s[1] for s in th_stats)
        ess_tail_min = min(s[2] for s in th_stats)

        if any(fs is None for fs in fixed_sigmas):
            sig_stack = np.stack([t["sigma"] for t in traces_list], axis=0)
            s_stats = [
                _split_rhat_and_ess_1d(sig_stack[:, :, m])
                for m in range(n_outcomes)
                if fixed_sigmas[m] is None
            ]
            rhat_blocks["sigma"] = max(s[0] for s in s_stats)
            ess_blocks["sigma"] = min(s[1] for s in s_stats)
            if is_p1_global:
                rhat_blocks["tau"] = rhat_blocks["sigma"]
                ess_blocks["tau"] = ess_blocks["sigma"]

        if opts["dynamic_shrinkage"]:
            kap_stack = np.stack([t["kappa"] for t in traces_list], axis=0)
            k_stats = [_split_rhat_and_ess_1d(kap_stack[:, :, m]) for m in range(n_outcomes)]
            rhat_blocks["kappa"] = max(s[0] for s in k_stats)
            ess_blocks["kappa"] = min(s[1] for s in k_stats)

        if any(spec.use_local_scales and spec.n_edges > 0 for spec in prior_specs):
            l_stats = []
            for m in range(n_outcomes):
                if not (prior_specs[m].use_local_scales and prior_specs[m].n_edges > 0):
                    continue
                loc_m = np.stack([t["local"][m] for t in traces_list], axis=0)
                for e_i in range(prior_specs[m].n_edges):
                    l_stats.append(_split_rhat_and_ess_1d(loc_m[:, :, e_i]))
            if l_stats:
                rhat_blocks["lambda"] = max(s[0] for s in l_stats)
                ess_blocks["lambda"] = min(s[1] for s in l_stats)

        rho_stack = np.stack([t["rho"] for t in traces_list], axis=0)
        rho_stats = [_split_rhat_and_ess_1d(rho_stack[:, :, m]) for m in range(n_outcomes)]
        rhat_blocks["rho"] = max(s[0] for s in rho_stats)
        ess_blocks["rho"] = min(s[1] for s in rho_stats)

        phi_stack = np.stack([t["phi"] for t in traces_list], axis=0)
        phi_stats = [_split_rhat_and_ess_1d(phi_stack[:, :, m]) for m in range(n_outcomes)]
        rhat_blocks["phi"] = max(s[0] for s in phi_stats)
        ess_blocks["phi"] = min(s[1] for s in phi_stats)

        if opts["robust_error"] == "t4":
            om_stack = np.stack([t["omega"] for t in traces_list], axis=0)
            o_stats = [_split_rhat_and_ess_1d(om_stack[:, :, q]) for q in range(om_stack.shape[2])]
            rhat_blocks["omega"] = max(s[0] for s in o_stats)
            ess_blocks["omega"] = min(s[1] for s in o_stats)

        if sample_corr and n_outcomes >= 2:
            c_stack = np.stack([t["corr"] for t in traces_list], axis=0)
            c_stats = [
                _split_rhat_and_ess_1d(c_stack[:, :, i, j])
                for i in range(n_outcomes)
                for j in range(i)
            ]
            rhat_blocks["corr"] = max(s[0] for s in c_stats)
            ess_blocks["corr"] = min(s[1] for s in c_stats)

        overall_rhat = max(rhat_blocks.values())
        overall_ess = min(ess_blocks.values())
        ok = (overall_rhat <= 1.01) and (overall_ess >= 400.0)
        return ok, overall_rhat, overall_ess, ess_tail_min, rhat_blocks, ess_blocks

    converged, max_rhat, min_ess, min_ess_tail, rhat_by_block, ess_by_block = _evaluate_convergence(chain_traces)
    cont_round = 0
    thin_stride = 1
    eff_max_rounds = (
        0
        if (opts["scale_prior"] == "local_halft" and not prior_specs[0].is_tree_poset)
        else max_rounds
    )
    while not converged and cont_round < eff_max_rounds:
        cont_round += 1
        for c_i in range(n_chains):
            chain_states[c_i], extra_tr = _run_single_gibbs_chain(
                n_steps=n_draws_chain,
                thin=thin_stride,
                state=chain_states[c_i],
                record=True,
                y_vec=y_vec,
                x_des=x_des,
                g_xi=g_xi,
                f_blocks=f_blocks,
                chol_f_blocks=chol_f_blocks,
                f_inv_blocks=f_inv_blocks,
                w_f_stack=w_f_stack,
                winv_f_stack=winv_f_stack,
                d_ratio_stack=d_ratio_stack,
                q_member_weights=q_member_weights,
                sqrt_v_arms=sqrt_v_arms,
                prior_specs=prior_specs,
                a_cone_joint=a_cone_joint,
                interior_ray_joint=interior_ray_joint,
                s0_scales=s0_scales,
                sk_scales=sk_scales,
                fixed_sigmas=fixed_sigmas,
                n_arms=n_arms,
                n_outcomes=n_outcomes,
                k_dim=k_dim,
                horizon=h_eff,
                robust_error=str(opts["robust_error"]),
                dynamic_shrinkage=bool(opts["dynamic_shrinkage"]),
                sample_corr=sample_corr,
                iw_df=iw_df,
                iw_scale=iw_scale,
                rng=chain_rngs[c_i],
            )
            assert extra_tr is not None
            # Concatenate n_draws_chain + n_draws_chain at thin_stride and thin by 2 -> spacing 2 * thin_stride
            merged_c: dict[str, Any] = {}
            for key in ("theta", "sigma", "kappa", "rho", "phi", "omega", "flag", "corr"):
                cat = np.concatenate([chain_traces[c_i][key], extra_tr[key]], axis=0)
                merged_c[key] = cat[1::2][:n_draws_chain]
            merged_local = []
            for m in range(n_outcomes):
                cat_l = np.concatenate([chain_traces[c_i]["local"][m], extra_tr["local"][m]], axis=0)
                merged_local.append(cat_l[1::2][:n_draws_chain])
            merged_c["local"] = merged_local
            chain_traces[c_i] = merged_c
        thin_stride *= 2
        converged, max_rhat, min_ess, min_ess_tail, rhat_by_block, ess_by_block = _evaluate_convergence(chain_traces)

    diagnostics = SamplerDiagnostics(
        sampler="pakman_paninski_hmc_gibbs",
        r_hat=float(max_rhat),
        ess=float(min_ess),
        ess_tail=float(min_ess_tail),
        n_chains=n_chains,
        n_draws_per_chain=int(chain_traces[0]["theta"].shape[0]),
        n_warmup=n_warmup,
        n_continuation_rounds=cont_round,
        converged=bool(converged),
        max_rhat_by_block=rhat_by_block,
        min_ess_by_block=ess_by_block,
    )

    # Pool chains into (S, ...) posterior arrays
    theta_all = np.concatenate([t["theta"] for t in chain_traces], axis=0)  # (S, M * (A-1) * K')
    sigma_all = np.concatenate([t["sigma"] for t in chain_traces], axis=0)  # (S, M)
    kappa_all = np.concatenate([t["kappa"] for t in chain_traces], axis=0)  # (S, M)
    rho_all = np.concatenate([t["rho"] for t in chain_traces], axis=0)  # (S, M)
    phi_all = np.concatenate([t["phi"] for t in chain_traces], axis=0)  # (S, M)
    omega_all = np.concatenate([t["omega"] for t in chain_traces], axis=0)  # (S, Q)
    flag_all = np.concatenate([t["flag"] for t in chain_traces], axis=0)  # (S, Q)
    corr_all = np.concatenate([t["corr"] for t in chain_traces], axis=0)  # (S, M, M)
    local_all = [
        np.concatenate([t["local"][m] for t in chain_traces], axis=0)
        for m in range(n_outcomes)
    ]
    n_total_draws = int(theta_all.shape[0])

    # Build per-outcome MultiArmReconciliationResult objects (structural + composite)
    theta_by_m: dict[str, np.ndarray] = {}
    for m_i, m_name in enumerate(structural_outcomes):
        raw_m = theta_all[:, m_i * d_per_m : (m_i + 1) * d_per_m]  # (S, K' * (A - 1))
        # Reshape (S, K', A - 1) -> transpose to (S, A - 1, K')
        theta_by_m[m_name] = np.swapaxes(raw_m.reshape(n_total_draws, k_dim, n_nonctrl), 1, 2).copy()

    for c_name, members_names in composite_outcomes.items():
        theta_by_m[c_name] = sum(theta_by_m[m] for m in members_names)  # exact linear sum!

    def _eval_f_by_pair_posterior_mean(
        summ_d: dict[str, Any], rho_bar: float, phi_bar: float, pairs_seq: Sequence[tuple[str, str]]
    ) -> dict[tuple[str, str], np.ndarray]:
        out_f: dict[tuple[str, str], np.ndarray] = {}
        for p in pairs_seq:
            w_q = summ_d["w_f_by_pair"][p]
            d_q = summ_d["d_ratio_by_pair"][p]
            e_q = np.maximum(1.0 - rho_bar * d_q, phi_bar)
            out_f[p] = (w_q * e_q[None, :]) @ w_q.T
        return out_f

    def _ppc_child_rng(kind_tag: int, name: str) -> np.random.Generator:
        """Deterministic outcome-order-independent PRNG stream for posterior predictive checks."""
        entropy = [int(seed) & 0xFFFFFFFF, 7919, kind_tag, *(ord(ch) for ch in name)]
        return np.random.default_rng(np.random.SeedSequence(entropy))

    # Compute composite PPC p-values where composite_results were provided
    comp_ppc_by_name: dict[str, float | None] = {}
    for c_name in composite_outcomes:
        if c_name in comp_summaries:
            c_pairs = comp_summaries[c_name]["pairs"]
            c_om = {p: omega_all[:, meas_index[("comp", c_name, p)]] for p in c_pairs}
            mem_names = composite_outcomes[c_name]
            mem_idxs = [m_name_to_idx[m] for m in mem_names]
            rho_bar_c = float(np.mean(rho_all[:, mem_idxs]))
            phi_bar_c = float(np.mean(phi_all[:, mem_idxs]))
            comp_f_post = _eval_f_by_pair_posterior_mean(comp_summaries[c_name], rho_bar_c, phi_bar_c, c_pairs)
            struct_f_post_list = [
                _eval_f_by_pair_posterior_mean(
                    struct_summaries[m],
                    float(np.mean(rho_all[:, m_name_to_idx[m]])),
                    float(np.mean(phi_all[:, m_name_to_idx[m]])),
                    c_pairs,
                )
                for m in mem_names
            ]
            comp_ppc_by_name[c_name] = _compute_composite_ppc_pvalue(
                comp_y_by_pair=comp_summaries[c_name]["y_by_pair"],
                comp_f_by_pair=comp_f_post,
                comp_omega_draws=c_om,
                struct_y_list=[struct_summaries[m]["y_by_pair"] for m in mem_names],
                struct_f_list=struct_f_post_list,
                struct_omega_list=[
                    {p: omega_all[:, meas_index[("struct", m, p)]] for p in c_pairs}
                    for m in mem_names
                ],
                comp_pairs=c_pairs,
                k_dim=k_dim,
                rng=_ppc_child_rng(1, c_name),
                comp_active_modes=comp_summaries[c_name]["active_modes_by_pair"],
                struct_active_modes_list=[struct_summaries[m]["active_modes_by_pair"] for m in mem_names],
            )
        else:
            comp_ppc_by_name[c_name] = None

    by_outcome_results: dict[str, MultiArmReconciliationResult] = {}
    all_outcome_names = (*structural_outcomes, *tuple(composite_outcomes.keys()))
    for out_name in all_outcome_names:
        is_struct = out_name in m_name_to_idx
        th_out = theta_by_m[out_name]  # (S, A - 1, K')

        # Full arm trajectories (S, A, H) with arm 0 identically zero
        arm_trajs_arr = np.zeros((n_total_draws, n_arms, h_eff), dtype=np.float64)
        arm_trajs_arr[:, 1:, :] = np.einsum("sak,hk->sah", th_out, psi)
        arm_traj_dict = {arms_tuple[a]: arm_trajs_arr[:, a, :] for a in range(n_arms)}
        tau_mu_arr = np.zeros((n_total_draws, n_arms), dtype=np.float64)
        tau_mu_arr[:, 1:] = th_out[:, :, 0] / math.sqrt(float(h_eff))

        # Contrasts for all_pairs AND any fitted_pairs with reversed orientation
        pair_set_to_store = list(topo["all_pairs"])
        for fp in fitted_pairs:
            if fp not in pair_set_to_store:
                pair_set_to_store.append(fp)

        c_traj: dict[tuple[str, str], np.ndarray] = {}
        c_mode1: dict[tuple[str, str], np.ndarray] = {}
        c_cum: dict[tuple[str, str], np.ndarray] = {}
        for u, v in pair_set_to_store:
            iu, iv = arm_to_idx[u], arm_to_idx[v]
            tr_uv = arm_trajs_arr[:, iu, :] - arm_trajs_arr[:, iv, :]
            c_traj[(u, v)] = tr_uv
            c_mode1[(u, v)] = tau_mu_arr[:, iu] - tau_mu_arr[:, iv]
            c_cum[(u, v)] = np.sum(tr_uv, axis=1)

        tau_inc = np.column_stack(
            [tau_mu_arr[:, arm_to_idx[ch]] - tau_mu_arr[:, arm_to_idx[pa]] for ch, pa in topo["tree_edges"]]
        )

        # Max cycle discrepancy across all 3-cycles (machine zero ~ 1e-16)
        max_cyc = 0.0
        for c_i in range(n_arms):
            for b_i in range(c_i + 1, n_arms):
                for a_i in range(b_i + 1, n_arms):
                    pa, pb, pc = arms_tuple[a_i], arms_tuple[b_i], arms_tuple[c_i]
                    gap = float(np.max(np.abs(c_traj[(pa, pc)] - c_traj[(pa, pb)] - c_traj[(pb, pc)])))
                    max_cyc = max(max_cyc, gap)

        if is_struct:
            m_i = m_name_to_idx[out_name]
            summ_o = struct_summaries[out_name]
            om_dict = {p: omega_all[:, meas_index[("struct", out_name, p)]] for p in fitted_pairs}
            fl_dict = {p: flag_all[:, meas_index[("struct", out_name, p)]] for p in fitted_pairs}
            infl_dict = {p: float(np.mean(1.0 / np.maximum(om_dict[p], 1e-12))) for p in fitted_pairs}
            prob_dict = {p: float(np.mean(fl_dict[p])) for p in fitted_pairs}
            f_post_o = _eval_f_by_pair_posterior_mean(
                summ_o, float(np.mean(rho_all[:, m_i])), float(np.mean(phi_all[:, m_i])), fitted_pairs
            )
            ppc_cyc = _compute_cycle_ppc_pvalue(
                summ_o["y_by_pair"],
                f_post_o,
                om_dict,
                fitted_pairs,
                topo["n_cycle"],
                k_dim,
                _ppc_child_rng(2, out_name),
                active_modes_by_pair=summ_o["active_modes_by_pair"],
            )
            ppc_bas = _compute_basis_ppc_pvalue(
                summ_o["traj_mean_by_pair"], summ_o["traj_cov_by_pair"], psi, h_eff
            )
            ppc_comp = None
            sig_draws_o = sigma_all[:, m_i]
            kap_draws_o = kappa_all[:, m_i]
            loc_draws_o = local_all[m_i]
            rho_draws_o = rho_all[:, m_i]
            phi_draws_o = phi_all[:, m_i]
            mono_o = bool(monotone_by_outcome[out_name])
            ps_o = prior_scale_by_outcome[out_name]
            s0_o = s0_scales[m_i]
        else:
            mem_names = composite_outcomes[out_name]
            mem_idxs = [m_name_to_idx[m] for m in mem_names]
            rho_draws_o = np.mean(rho_all[:, mem_idxs], axis=1)
            phi_draws_o = np.mean(phi_all[:, mem_idxs], axis=1)
            if out_name in comp_summaries:
                c_summ = comp_summaries[out_name]
                c_pairs = c_summ["pairs"]
                om_dict = {p: omega_all[:, meas_index[("comp", out_name, p)]] for p in c_pairs}
                fl_dict = {p: flag_all[:, meas_index[("comp", out_name, p)]] for p in c_pairs}
                infl_dict = {p: float(np.mean(1.0 / np.maximum(om_dict[p], 1e-12))) for p in c_pairs}
                prob_dict = {p: float(np.mean(fl_dict[p])) for p in c_pairs}
                f_post_c = _eval_f_by_pair_posterior_mean(
                    c_summ, float(np.mean(rho_draws_o)), float(np.mean(phi_draws_o)), c_pairs
                )
                ppc_cyc = (
                    _compute_cycle_ppc_pvalue(
                        c_summ["y_by_pair"],
                        f_post_c,
                        om_dict,
                        fitted_pairs,
                        topo["n_cycle"],
                        k_dim,
                        _ppc_child_rng(3, out_name),
                        active_modes_by_pair=c_summ["active_modes_by_pair"],
                    )
                    if set(c_pairs) == set(fitted_pairs)
                    else None
                )
                ppc_bas = _compute_basis_ppc_pvalue(
                    c_summ["traj_mean_by_pair"], c_summ["traj_cov_by_pair"], psi, h_eff
                )
                s0_o = float(c_summ["s0_scale"])
            else:
                om_dict = {p: np.ones(n_total_draws, dtype=np.float64) for p in fitted_pairs}
                infl_dict = {p: 1.0 for p in fitted_pairs}
                prob_dict = {p: 0.0 for p in fitted_pairs}
                ppc_cyc = None
                ppc_bas = None
                s0_o = math.sqrt(sum(s0_scales[m_name_to_idx[m]] ** 2 for m in mem_names))
            ppc_comp = comp_ppc_by_name[out_name]
            sig_draws_o = np.sqrt(sum(sigma_all[:, m_name_to_idx[m]] ** 2 for m in mem_names))
            kap_draws_o = np.sqrt(sum(kappa_all[:, m_name_to_idx[m]] ** 2 for m in mem_names))
            loc_draws_o = np.ones((n_total_draws, len(topo["cover_edges"])), dtype=np.float64)
            mono_o = all(bool(monotone_by_outcome[m]) for m in mem_names)
            ps_o = None

        by_outcome_results[out_name] = MultiArmReconciliationResult(
            arms=arms_tuple,
            order=topo["order"],
            cover_edges=topo["cover_edges"],
            topology=topo["topology"],
            tree_edges=topo["tree_edges"],
            fitted_pairs=fitted_pairs,
            redundant_pairs=topo["redundant_pairs"],
            all_pairs=topo["all_pairs"],
            horizon=h_eff,
            n_modes=k_dim,
            basis=psi,
            monotone=mono_o,
            prior_scale=ps_o,
            s0_scale=s0_o,
            theta_draws=th_out,
            sigma_draws=sig_draws_o,
            kappa_draws=kap_draws_o,
            local_scale_draws=loc_draws_o,
            rho_draws=rho_draws_o,
            phi_draws=phi_draws_o,
            tau_draws=sig_draws_o / max(s0_o, 1e-16),
            omega_draws=om_dict,
            tau_mu_draws=tau_mu_arr,
            tau_increment_draws=tau_inc,
            arm_trajectory_draws=arm_traj_dict,
            contrast_mode1_draws=c_mode1,
            contrast_trajectory_draws=c_traj,
            contrast_cumulative_draws=c_cum,
            contrast_catt_draws=c_traj,
            contrast_scale_inflation=infl_dict,
            posterior_corruption_prob=prob_dict,
            ppc_cycle_pvalue=ppc_cyc,
            ppc_composite_pvalue=ppc_comp,
            ppc_basis_pvalue=ppc_bas,
            max_cycle_discrepancy=max_cyc,
            sampler_diagnostics=diagnostics,
        )

    return {
        "by_outcome": by_outcome_results,
        "corr_draws": corr_all,
        "diagnostics": diagnostics,
        "topology": topo,
    }


def reconcile_multi_arm_draws(
    results: Mapping[tuple[str, str], RankKLongBetResult],
    *,
    arms: Sequence[str] | None = None,
    order: Sequence[str] | None = None,
    poset: Sequence[tuple[str, str]] | None = None,
    monotone: bool = True,
    prior_scale: float | None = None,
    horizon: int | None = None,
    n_modes: int = 2,
    seed: int = 2026,
    _model_options: Mapping[str, Any] | None = None,
) -> MultiArmReconciliationResult:
    """Reconcile pairwise Stage-1 Rank-K posterior fits across a multi-arm poset.

    Args:
      results: Mapping ``{(arm_high, arm_low): RankKLongBetResult}`` forming a
        connected contrast graph over the arms.
      arms: Optional explicit sequence of arm labels (``arms[0]`` is control).
      order: Optional total order ``(control, ..., top)`` for chain designs
        (mutually exclusive with ``poset``).
      poset: Optional sequence of ``(arm_high, arm_low)`` cover relations
        (transitive reduction computed automatically; mutually exclusive with ``order``).
      monotone: Enforce poset-cone monotonicity and per-period positivity in
        prior support (default ``True``).
      prior_scale: ``None`` (default) places the hierarchical per-outcome
        generalized-horseshoe prior ``(sigma / s_0)^2 ~ BetaPrime(1/4, 1/4)`` on
        the top-arm level scale (``s_0`` the pre-treatment outcome SD); a
        positive ``float`` fixes an externally calibrated prior scale in outcome
        units. Strings raise ``ValueError``.
      horizon: Estimand window ``1..H`` (default: full post-adoption window).
      n_modes: Exposure basis dimension ``K'`` (default ``2``: level + decay).
      seed: Deterministic PRNG seed for the 4-chain sampler.
      _model_options: Private internal options mapping for candidate comparisons
        and ablations; ``None`` is the shipped configuration.

    Returns:
      :class:`MultiArmReconciliationResult`.
    """
    if not isinstance(results, Mapping) or not results:
        raise ValueError("`results` must be a non-empty mapping of {(arm_high, arm_low): RankKLongBetResult}")
    if not isinstance(monotone, (bool, np.bool_)):
        raise ValueError(f"`monotone` must be a bool; got {monotone!r}")
    ps_val = _validate_prior_scale_scalar(prior_scale, "prior_scale")

    core = _reconcile_joint_bayesian_core(
        structural_outcomes=("_single",),
        results_by_outcome={"_single": results},
        composite_outcomes={},
        composite_results={},
        arms=arms,
        order=order,
        poset=poset,
        monotone_by_outcome={"_single": bool(monotone)},
        prior_scale_by_outcome={"_single": ps_val},
        iw_spec=None,
        horizon=horizon,
        n_modes=n_modes,
        seed=seed,
        model_options=_model_options,
    )
    return core["by_outcome"]["_single"]
