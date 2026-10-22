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

"""Coupled multi-outcome LongBet estimator and prediction API."""

from __future__ import annotations

import dataclasses
import uuid
import warnings
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
import scipy.linalg
import scipy.stats as stats
from jaxtyping import Array, Key

from longbet._config import LongBetConfig, resolve_model_defaults
from longbet._model import (
    LongBet,
    LongBetPrediction,
    _check_absorbing,
    _check_time_vector,
    derive_exposure,
    resolve_device,
)
from longbet._multi_input import NormalizedMultiInput, normalize_multi_inputs
from longbet._multi_io import load_multi_npz, save_multi_npz
from longbet._multi_loop import MultiLongBetTrace, run_multi_longbet_mcmc
from longbet._multi_state import MultiLongBetState, init_multi_longbet
from longbet._sur import SAMPLER_SEMANTICS
from longbet._tempering import swap_rates, warn_if_not_exchanging
from longbet._trend_basis import (
    build_period_basis,
    build_time_basis,
    build_trend_design,
    build_unit_features,
)


class OutcomeList(list):
    """Sequence of per-outcome objects that can be indexed by integer or name."""

    def __init__(self, items: Sequence[Any], names: Sequence[str]) -> None:
        super().__init__(items)
        self._names = tuple(names)
        self._by_name = {name: item for name, item in zip(names, items)}

    def __getitem__(self, key: Any) -> Any:
        if isinstance(key, str):
            if key not in self._by_name:
                raise KeyError(
                    f"Unknown outcome name {key!r}. Available: {list(self._names)}"
                )
            return self._by_name[key]
        return super().__getitem__(key)

    def __contains__(self, key: Any) -> bool:
        if isinstance(key, str):
            return key in self._by_name
        return super().__contains__(key)


class LongBetMultiPrediction:
    """Predictions across multiple outcomes.

    Attributes
    ----------
    preds
        Sequence of child ``LongBetPrediction`` objects in USER order, indexable
        by name (e.g. ``pred['y1']``) or position (``pred[0]``).
    outcome_names
        Names of the outcomes in user order.
    outcome
        Outcome types in user order ('continuous' or 'binary').
    num_chains
        Number of MCMC chains.
    num_sweeps
        Number of retained draws per chain.
    summary_only
        Whether full draws were suppressed.
    sampler_semantics
        Sampler semantic tag; joint-event inference rejects legacy semantics.
    provenance
        Unique fit/predict run identifier.
    """

    def __init__(
        self,
        *,
        preds: Sequence[LongBetPrediction],
        outcome_names: Sequence[str],
        outcome: Sequence[str],
        num_chains: int,
        num_sweeps: int,
        summary_only: bool,
        provenance: str,
        sampler_semantics: str = SAMPLER_SEMANTICS,
    ) -> None:
        self.outcome_names = tuple(outcome_names)
        self.outcome = tuple(outcome)
        self.num_chains = int(num_chains)
        self.num_sweeps = int(num_sweeps)
        self.summary_only = bool(summary_only)
        self.provenance = str(provenance)
        self.sampler_semantics = str(sampler_semantics)
        self.preds = OutcomeList(preds, self.outcome_names)

    def __getitem__(self, key: Any) -> LongBetPrediction:
        return self.preds[key]

    def __len__(self) -> int:
        return len(self.preds)

    def stability(
        self,
        min_ess: float = 400.0,
        max_rhat: float = 1.01,
        alpha: float = 0.05,
        warn: bool = True,
    ) -> dict[str, Any]:
        """Compute MCMC stability diagnostics per outcome without pooling outcomes.

        Returns
        -------
        Dictionary mapping each outcome name to its StabilityResult.
        """
        return {
            name: self.preds[name].stability(
                min_ess=min_ess, max_rhat=max_rhat, alpha=alpha, warn=warn
            )
            for name in self.outcome_names
        }


class LongBetMulti:
    """Coupled trees for multiple continuous and binary outcomes.

    Fits separate prognostic and treatment forests per outcome on a shared
    panel, using the full precision of an identified triangular SUR likelihood.
    Downstream observations inform earlier means and binary latent responses.
    Only the residual likelihood couples the outcome-specific forests.
    Convergence and calibration must not be inferred from aligned draws or a
    joint-event counting utility.

    Coupled fits containing continuous outcomes require explicitly positive
    ``sigma_prior_a`` and ``sigma_prior_b``. For standardized responses,
    ``LongBetMulti(sigma_prior_a=2, sigma_prior_b=1)`` selects a proper IG(2,1)
    innovation prior. This is a model choice, not a convergence guarantee.
    """

    def __init__(
        self,
        config: LongBetConfig | None = None,
        **kwargs: Any,
    ) -> None:
        if config is None:
            self.config = LongBetConfig(**kwargs)
        elif kwargs:
            self.config = dataclasses.replace(config, **kwargs)
        else:
            self.config = config
        self.config = resolve_model_defaults(self.config)

        self.swap_accepts_: np.ndarray | None = None
        self.swap_rates_: np.ndarray | None = None
        self.fits: OutcomeList | None = None
        self.outcome_names: tuple[str, ...] | None = None
        self.outcome: tuple[str, ...] | None = None
        self.num_categories: tuple[int | None, ...] | None = None
        self.internal_num_categories: tuple[int | None, ...] | None = None
        self.order: tuple[int, ...] | None = None
        self.inverse_order: tuple[int, ...] | None = None
        self.sur_active: bool = False
        self.trace: MultiLongBetTrace | None = None
        self.state: MultiLongBetState | None = None
        self.design_: Any = None
        self.N_: int | None = None
        self.T_: int | None = None
        self.S_max_: int | None = None
        self.fitted_t_: np.ndarray | None = None
        self.meany: tuple[float, ...] | None = None
        self.sdy: tuple[float, ...] | None = None
        self.offset_: tuple[float, ...] | None = None
        # Calendar-time trend block, shared by every outcome because it depends
        # only on (x, t). Mirrors the three attributes ``LongBet`` carries, and
        # is copied verbatim onto each child fit so ``child.predict`` can
        # rebuild the block.
        self.trend_time_basis_: np.ndarray | None = None
        self.trend_period_basis_: np.ndarray | None = None
        self.trend_feature_spec_: Any = None
        self.trend_forecast_spec_: Any = None
        self.sampler_semantics: str = SAMPLER_SEMANTICS
        self.rng_scheme_version: int = 1
        self.provenance: str = ""

    def __getitem__(self, key: Any) -> LongBet:
        if self.fits is None:
            raise RuntimeError("Model has not been fitted yet.")
        return self.fits[key]

    def __len__(self) -> int:
        if self.outcome_names is None:
            return 0
        return len(self.outcome_names)

    @property
    def Gamma_draws(self) -> np.ndarray:
        """Loading draws shaped (M*M, D) with row-major matrix entries per draw.

        Outcomes follow internal sampling order, with strict loadings and zero
        diagonal (unlike the reference package's unit-diagonal storage).
        Columns are chain-major draws (c*K + k), excluding burn-in.
        """
        if self.trace is None:
            raise RuntimeError("Model is not fitted.")
        gl = np.asarray(self.trace.gamma_loadings)
        M = len(self.outcome_names)
        if gl.ndim == 4:  # (C, K, M, M)
            C, K, _, _ = gl.shape
            flat = gl.reshape(C * K, M * M)
        else:  # (K, M, M)
            K, _, _ = gl.shape
            flat = gl.reshape(K, M * M)
        return flat.T

    def stability(
        self,
        min_ess: float = 400.0,
        max_rhat: float = 1.01,
        alpha: float = 0.05,
        warn: bool = True,
    ) -> dict[str, Any]:
        """Compute MCMC stability diagnostics per outcome without pooling outcomes.

        Returns
        -------
        Dictionary mapping each outcome name to its StabilityResult.
        """
        if self.fits is None:
            raise RuntimeError("Model is not fitted.")
        return {
            name: self.fits[name].stability(
                min_ess=min_ess, max_rhat=max_rhat, alpha=alpha, warn=warn
            )
            for name in self.outcome_names
        }

    def fit(
        self,
        y: Any,
        x: np.ndarray | Array,
        z: np.ndarray | Array,
        t: np.ndarray | Array | None = None,
        x_trt: np.ndarray | Array | None = None,
        x_tv: np.ndarray | Array | None = None,
        x_trt_tv: np.ndarray | Array | None = None,
        ps: np.ndarray | Array | None = None,
        key: Key[Array, ''] | None = None,
        *,
        outcome: str | Sequence[str] | Mapping[str, str] | None = None,
        outcome_names: Sequence[str] | None = None,
        num_categories: int | Sequence[int | None] | Mapping[str, int | None] | None = None,
    ) -> LongBetMulti:
        """Fit the coupled multi-outcome LongBet model."""
        if (self.config.linear_trt or self.config.learn_beta_scale or self.config.ar1_errors
                or self.config.interweave_gamma):
            raise NotImplementedError(
                "linear_trt, learn_beta_scale and ar1_errors are scalar-model options; "
                "LongBetMulti does not support them")
        cfg = self.config

        # 1. Validate and normalize multi-outcome responses
        norm_in = normalize_multi_inputs(
            y,
            outcome=outcome,
            outcome_names=outcome_names,
            config=cfg,
            num_categories=num_categories,
        )

        N, T, M = norm_in.N, norm_in.T, norm_in.M
        self.N_ = N
        self.T_ = T
        self.outcome_names = norm_in.outcome_names
        self.outcome = norm_in.user_outcomes
        self.num_categories = norm_in.user_num_categories
        self.internal_num_categories = norm_in.internal_num_categories
        self.order = norm_in.order
        self.inverse_order = norm_in.inverse_order
        self.sur_active = cfg.sur_active
        self.meany = norm_in.meany
        self.sdy = norm_in.sdy
        self.offset_ = norm_in.offset_
        self.provenance = str(uuid.uuid4())
        # A fresh fit uses the current sampler, regardless of metadata left by
        # an earlier fit on this object.
        self.sampler_semantics = SAMPLER_SEMANTICS

        # 2. Key schedule (Section 6.1, rng_scheme_version=1)
        if key is None:
            master = jax.random.key(cfg.random_seed)
        else:
            master = key

        k_design, k_init, k_mcmc = jax.random.split(master, 3)

        # 3. Covariate validation, exposure, and shared Design construction
        z_np = np.asarray(z)
        if z_np.ndim != 2:
            raise ValueError(f"z must be an (N, T) array, got shape {z_np.shape}")
        if z_np.shape != (N, T):
            raise ValueError(f"z has shape {z_np.shape}, expected {(N, T)}")
        _check_absorbing(z_np)

        x_np = np.asarray(x, dtype=np.float32)
        if x_np.ndim != 2 or x_np.shape[0] != N:
            raise ValueError(f"x must be (N, P) with N={N}, got {x_np.shape}")

        t_vec = (
            np.arange(1, T + 1, dtype=np.float32)
            if t is None
            else np.asarray(t, dtype=np.float32)
        )
        if t_vec.size != T:
            raise ValueError(f"t has length {t_vec.size}, expected {T}")
        _check_time_vector(t_vec)

        if cfg.random_intercept:
            n_always = int(np.sum(np.all(z_np == 1, axis=1)))
            if n_always:
                warnings.warn(
                    f"{n_always} unit(s) are treated in every period, so their random "
                    "intercept is not separable from their treatment effect. Consider "
                    "dropping them or setting random_intercept=False.",
                    UserWarning,
                    stacklevel=3,
                )

        S_mat = derive_exposure(z_np, t_vec)
        S_max = int(S_mat.max())
        self.S_max_ = S_max
        self.fitted_t_ = t_vec

        unit_idx = np.repeat(np.arange(N, dtype=np.int32), T)
        time_idx = np.tile(np.arange(T, dtype=np.int32), N)
        exposure_idx = S_mat.ravel().astype(np.int32)
        z_vec = z_np.ravel().astype(np.float32)

        temp_scalar = LongBet(cfg)
        # _build_design uses these dimensions to construct calendar/exposure
        # cutpoints. Leaving the temporary estimator's zero-valued defaults
        # silently bins BOTH time axes to zero, unlike the scalar fit.
        temp_scalar.N_ = N
        temp_scalar.T_ = T
        temp_scalar.S_max_ = S_max
        raw_covars = temp_scalar._raw_inputs(
            x=x_np,
            x_trt=x_trt,
            x_tv=x_tv,
            x_trt_tv=x_trt_tv,
            ps=ps,
            time_idx=time_idx,
            exposure_idx=exposure_idx,
            N=N,
            T=T,
        )

        shared_design = temp_scalar._build_design(
            raw_covars,
            x_trt_given=x_trt is not None,
            key=k_design,
        )
        self.design_ = shared_design

        # Assemble unified binned predictors (P, N*T)
        X_unified_np = shared_design.build(raw_covars)
        X_unified = jnp.asarray(X_unified_np, dtype=jnp.uint8)

        # Conjugate calendar-time block, built exactly as ``LongBet._fit_impl``
        # builds it and shared by every outcome: the bases depend only on
        # ``(T, x)``, both of which are common to the M equations, and no chain
        # key touches them. Without this the shipped default
        # ``split_calendar_mu=False`` would leave the coupled model with no way
        # to represent calendar time at all.
        trend_design = None
        trend_num_period = 0
        if cfg.use_trend_block:
            self.trend_time_basis_ = build_time_basis(T, cfg.trend_degree)
            self.trend_period_basis_ = (
                build_period_basis(T) if cfg.trend_period_effects
                else np.zeros((T, 0), dtype=np.float32)
            )
            features, self.trend_feature_spec_ = build_unit_features(
                x_np, cfg.trend_num_features
            )
            n_inter = self.trend_time_basis_.shape[1] * features.shape[1]
            trend_num_period = int(self.trend_period_basis_.shape[1])
            if n_inter + trend_num_period > 0:
                trend_design = build_trend_design(
                    features, self.trend_time_basis_, unit_idx, time_idx,
                    period=self.trend_period_basis_,
                )
            else:
                self.trend_time_basis_ = None
                self.trend_period_basis_ = None
                self.trend_feature_spec_ = None
                trend_num_period = 0
        else:
            self.trend_time_basis_ = None
            self.trend_period_basis_ = None
            self.trend_feature_spec_ = None

        unit_idx = jnp.asarray(unit_idx, dtype=jnp.int32)
        time_idx = jnp.asarray(time_idx, dtype=jnp.int32)
        exposure_idx = jnp.asarray(exposure_idx, dtype=jnp.int32)
        z_vec = jnp.asarray(z_vec, dtype=jnp.float32)

        max_split_mu = jnp.asarray(shared_design.max_split_mu, dtype=jnp.uint8)
        max_split_nu = jnp.asarray(shared_design.max_split_nu, dtype=jnp.uint8)

        # 4. Device placement and MCMC execution
        device = resolve_device(cfg.device)
        with jax.default_device(device):
            multi_state = init_multi_longbet(
                X_unified=X_unified,
                unit_idx=unit_idx,
                time_idx=time_idx,
                exposure_idx=exposure_idx,
                z_vec=z_vec,
                max_split_mu=max_split_mu,
                max_split_nu=max_split_nu,
                norm_input=norm_in,
                config=cfg,
                trend_design=trend_design,
                trend_num_period=trend_num_period,
                num_chains=cfg.num_chains * cfg.tempering_levels,
                init_key=k_init,
            )

            result = run_multi_longbet_mcmc(
                key=k_mcmc,
                state=multi_state,
                n_burn=cfg.num_burnin,
                n_save=cfg.num_sweeps,
                n_skip=cfg.n_skip,
                inner_loop_length=cfg.inner_loop_length,
            )

        self.state = result.final_state
        self.trace = result.main_trace
        self.swap_accepts_ = (None if result.swap_accepts is None
                              else np.asarray(result.swap_accepts))
        self.swap_rates_ = (None if self.swap_accepts_ is None else swap_rates(
            self.swap_accepts_, cfg.num_chains, cfg.tempering_levels,
            cfg.num_burnin + cfg.num_sweeps * cfg.n_skip))
        if self.swap_rates_ is not None:
            warn_if_not_exchanging(self.swap_rates_)
        if "ordinal" in self.outcome:
            jax.block_until_ready(self.trace)

        # 5. Build child LongBet fits in USER order
        child_fits = []
        for u in range(M):
            internal_idx = self.inverse_order[u]
            child_otype = self.outcome[u]
            child_cfg = dataclasses.replace(cfg, outcome=child_otype,
                                             num_categories=self.num_categories[u])

            child_model = LongBet(child_cfg)
            child_model.design_ = shared_design
            child_model.trace = self.trace.traces[internal_idx]
            child_model.state = self.state.states[internal_idx]
            child_model.multi_origin = self._child_origin(u)
            child_model.meany = self.meany[internal_idx]
            child_model.sdy = self.sdy[internal_idx]
            child_model.offset_ = self.offset_[internal_idx]
            child_model.N_ = N
            child_model.T_ = T
            child_model.S_max_ = S_max
            child_model.t_fit_ = t_vec
            child_model.fitted_t_ = t_vec
            # The trend block is shared, so the children share the objects
            # rather than copies: ``LongBet.predict`` reads all three (the time
            # basis, the period basis and the feature spec) and reproduces the
            # in-sample trend only if they are the ones the sampler was fitted
            # against. The per-outcome coefficients live on
            # ``child_model.trace.trend_coef``, which _multi_loop fills in.
            child_model.trend_time_basis_ = self.trend_time_basis_
            child_model.trend_period_basis_ = self.trend_period_basis_
            child_model.trend_feature_spec_ = self.trend_feature_spec_

            child_fits.append(child_model)

        self.t_fit_ = t_vec
        self.fits = OutcomeList(child_fits, self.outcome_names)
        return self

    def _child_origin(self, u: int) -> dict[str, Any]:
        """Joint provenance survives extracting/saving an outcome's marginal fit."""
        return {"provenance": self.provenance, "sampler_semantics": self.sampler_semantics,
                "outcome_name": self.outcome_names[u],
                "outcomes": list(self.outcome)}

    def predict(
        self,
        x: np.ndarray | Array,
        z: np.ndarray | Array,
        t: np.ndarray | Array | None = None,
        x_trt: np.ndarray | Array | None = None,
        x_tv: np.ndarray | Array | None = None,
        x_trt_tv: np.ndarray | Array | None = None,
        ps: np.ndarray | Array | None = None,
        y: Any = None,
        summary_only: bool = False,
        alpha: float = 0.05,
        key: Key[Array, ''] | None = None,
        block_size: int | None = None,
        sig_knl: float | None = None,
        lambda_knl: float | None = None,
        cache_forest_evaluations: bool = False,
    ) -> LongBetMultiPrediction:
        """Predict treatment effects and outcomes for all outcomes."""
        del y  # Preserves signature compatibility with call forwarding
        if self.fits is None or self.trace is None:
            raise RuntimeError("Model has not been fitted yet; call fit() first.")

        cfg = self.config
        M = len(self.outcome_names)

        # PRNG root for prediction
        if key is None:
            k_pred = jax.random.key(cfg.random_seed + 1)
        else:
            k_pred = key

        child_preds = []
        for u in range(M):
            internal_idx = self.inverse_order[u]
            # Forecast GP innovations need different keys across outcomes (Section 6.1)
            k_child = jax.random.fold_in(k_pred, internal_idx)

            pred_u = self.fits[u].predict(
                x=x,
                z=z,
                t=t,
                x_trt=x_trt,
                x_tv=x_tv,
                x_trt_tv=x_trt_tv,
                ps=ps,
                summary_only=summary_only,
                alpha=alpha,
                key=k_child,
                block_size=block_size,
                sig_knl=sig_knl,
                lambda_knl=lambda_knl,
                cache_forest_evaluations=cache_forest_evaluations,
            )
            child_preds.append(pred_u)

        return LongBetMultiPrediction(
            preds=child_preds,
            outcome_names=self.outcome_names,
            outcome=self.outcome,
            num_chains=cfg.num_chains,
            num_sweeps=cfg.num_sweeps,
            summary_only=summary_only,
            provenance=self.provenance,
            sampler_semantics=self.sampler_semantics,
        )

    def save(self, path: str | Path) -> None:
        """Save the fitted multi-outcome model to an NPZ archive."""
        save_multi_npz(self, path)

    @classmethod
    def load(cls, path: str | Path) -> LongBetMulti:
        """Load a multi-outcome model from an NPZ archive."""
        return load_multi_npz(path)


def effect_draws(
    pred: LongBetPrediction | LongBetMultiPrediction,
    outcome: str | int | None = None,
) -> np.ndarray:
    """Return draws of treatment effects on the natural outcome scale.

    Parameters
    ----------
    pred
        A ``LongBetPrediction`` or ``LongBetMultiPrediction``.
    outcome
        For multi predictions, the outcome name or 0-based index. For scalar
        predictions, an optional validation override matching the fitted outcome.

    Returns
    -------
    Array of shape ``(N, T, D)`` where D = C * K is the total number of draws.
    For continuous outcomes this is ``tauhats``; for binary outcomes this is
    ``Phi(mu0 + tau) - Phi(mu0)``.
    """
    if isinstance(pred, LongBetMultiPrediction):
        if outcome is None:
            raise ValueError(
                "For LongBetMultiPrediction, an outcome name or index must be specified."
            )
        target_pred: LongBetPrediction = pred[outcome]
    elif isinstance(pred, LongBetPrediction):
        if outcome is not None and outcome != pred.outcome:
            raise ValueError(
                f"Explicit outcome {outcome!r} conflicts with prediction outcome "
                f"{pred.outcome!r}."
            )
        target_pred = pred
    else:
        raise TypeError(
            f"Expected LongBetPrediction or LongBetMultiPrediction, got {type(pred)}"
        )

    return effect_draws_from_arrays(
        target_pred.tauhats, target_pred.muhats0,
        outcome=target_pred.outcome, summary_only=target_pred.summary_only,
    )


def effect_draws_from_arrays(
    tauhats: Any,
    muhats0: Any = None,
    *,
    outcome: str,
    summary_only: bool = False,
) -> np.ndarray:
    """Shared effect transformation for Python predictions and serialized R arrays.

    Continuous effects stay on the supplied response scale (including log scales).
    Binary effects are probability differences in [-1, 1], so 0.01 is one
    percentage point. No live prediction handle is required.
    """
    if outcome == "ordinal":
        raise ValueError("Ordinal effects require selecting a category or score with "
                         "att_probabilities() or att_expected_score().")
    if summary_only or tauhats is None:
        raise ValueError(
            "effect_draws requires full posterior draws. Re-run predict with "
            "summary_only=False."
        )

    if outcome not in ("continuous", "binary"):
        raise ValueError("outcome must be 'continuous' or 'binary'.")
    tau = np.asarray(tauhats)
    if tau.ndim != 3 or any(size == 0 for size in tau.shape):
        raise ValueError("tauhats must have nonempty shape (N, T, D).")
    if not np.all(np.isfinite(tau)):
        raise ValueError("tauhats contains non-finite values.")
    if outcome == "continuous":
        return tau

    if muhats0 is None:
        raise ValueError(
            "effect_draws for binary outcomes requires full muhats0 draws. "
            "Re-run predict with summary_only=False."
        )

    mu0 = np.asarray(muhats0)
    if mu0.shape != tau.shape or not np.all(np.isfinite(mu0)):
        raise ValueError("muhats0 must be finite and have the same (N, T, D) shape as tauhats.")
    return stats.norm.cdf(mu0 + tau) - stats.norm.cdf(mu0)


def joint_prob(
    pred: LongBetMultiPrediction,
    conditions: Sequence[Callable[[np.ndarray], np.ndarray]]
    | Mapping[str, Callable[[np.ndarray], np.ndarray]],
    cells: np.ndarray | None = None,
) -> np.ndarray:
    """Empirical joint probability of multiple treatment effect events.

    Parameters
    ----------
    pred
        A ``LongBetMultiPrediction`` with full draws.
    conditions
        A sequence of callables in user order, or a mapping from outcome name
        to callable. Each callable receives an ``(N, T, D)`` effect array and
        must return an ``(N, T, D)`` boolean array.
    cells
        Optional boolean panel of shape ``(N, T)``. When None, returns an
        ``(N, T)`` array of cellwise joint probabilities. When provided, returns
        a length-N vector where each unit's probability is averaged over its
        selected periods (or NaN if no periods are selected for that unit).

    Returns
    -------
    Array of shape ``(N, T)`` when cells is None, or shape ``(N,)`` when cells is given.
    """
    if not isinstance(pred, LongBetMultiPrediction):
        raise TypeError(f"pred must be a LongBetMultiPrediction, got {type(pred)}")

    if getattr(pred, "sampler_semantics", None) != SAMPLER_SEMANTICS:
        raise ValueError(
            "Joint prediction uses legacy or unsupported sampler semantics. "
            "Refit from the original data and regenerate predictions."
        )

    if pred.summary_only:
        raise ValueError(
            "joint_prob requires full posterior draws; predict with summary_only=False."
        )

    M = len(pred.outcome_names)
    cond_list: list[Callable[[np.ndarray], np.ndarray]]

    if isinstance(conditions, Mapping):
        expected = set(pred.outcome_names)
        actual = set(conditions.keys())
        if expected != actual:
            missing = expected - actual
            extra = actual - expected
            raise ValueError(
                f"Conditions mapping keys must exactly match fitted outcome names: "
                f"missing {sorted(missing)}, extra {sorted(extra)}"
            )
        cond_list = [conditions[name] for name in pred.outcome_names]
    elif isinstance(conditions, Sequence):
        if len(conditions) != M:
            raise ValueError(
                f"Number of conditions ({len(conditions)}) must match number of "
                f"outcomes ({M})."
            )
        cond_list = list(conditions)
    else:
        raise TypeError(
            f"conditions must be a sequence or mapping of callables, got {type(conditions)}"
        )

    for idx, fn in enumerate(cond_list):
        if not callable(fn):
            raise TypeError(
                f"Condition for outcome {pred.outcome_names[idx]!r} is not callable."
            )

    # Incrementally evaluate conditions and reduce across outcomes
    masks = []
    ref_shape = None

    for idx, (name, cond_fn) in enumerate(zip(pred.outcome_names, cond_list)):
        eff = effect_draws(pred[name])
        if ref_shape is None:
            ref_shape = eff.shape
        if eff.shape != ref_shape:
            raise ValueError(f"Effect draws for outcome {name!r} have shape {eff.shape}, expected {ref_shape}.")
        if (
            pred[name].num_chains != pred.num_chains
            or eff.shape[-1] != pred.num_chains * pred.num_sweeps
        ):
            raise ValueError(f"Draw/chain metadata for outcome {name!r} is inconsistent with the joint prediction.")

        if not np.all(np.isfinite(eff)):
            raise ValueError(f"effect_draws for outcome {name!r} contains non-finite values.")

        mask = cond_fn(eff)
        if not isinstance(mask, np.ndarray) or mask.dtype != np.bool_:
            raise ValueError(
                f"Condition for outcome {name!r} must return a boolean numpy array, "
                f"got {type(mask)} with dtype {getattr(mask, 'dtype', None)}"
            )
        if mask.shape != ref_shape:
            raise ValueError(
                f"Condition for outcome {name!r} returned array of shape {mask.shape}, "
                f"expected {ref_shape}."
            )
        masks.append(mask)

    return reduce_joint_masks(masks, cells=cells)


def reduce_joint_masks(
    masks: Sequence[np.ndarray],
    cells: np.ndarray | None = None,
) -> np.ndarray:
    """Validate boolean masks from multiple conditions and reduce to joint probabilities.

    Parameters
    ----------
    masks
        Sequence of boolean arrays of shape (N, T, D).
    cells
        Optional boolean panel of shape (N, T).

    Returns
    -------
    Array of shape (N, T) when cells is None, or shape (N,) when cells is given.
    """
    if len(masks) == 0:
        raise ValueError("At least one mask array is required.")

    ref_shape = None
    joint_ok = None

    for idx, mask in enumerate(masks):
        mask_arr = np.asarray(mask)
        if mask_arr.dtype != np.bool_:
            raise ValueError(
                f"Mask at index {idx} must be a boolean array, got dtype {mask_arr.dtype}."
            )
        if ref_shape is None:
            ref_shape = mask_arr.shape
            if mask_arr.ndim != 3 or any(size == 0 for size in mask_arr.shape):
                raise ValueError(
                    f"Masks must be nonempty 3-D arrays (N, T, D), got shape {mask_arr.shape}."
                )
        elif mask_arr.shape != ref_shape:
            raise ValueError(
                f"Mask at index {idx} has shape {mask_arr.shape}, expected {ref_shape}."
            )

        if joint_ok is None:
            joint_ok = mask_arr.copy()
        else:
            joint_ok &= mask_arr

    assert joint_ok is not None
    p_cell = np.mean(joint_ok, axis=-1)  # shape (N, T)

    if cells is None:
        return p_cell

    cells_arr = np.asarray(cells)
    if cells_arr.dtype != np.bool_ or cells_arr.shape != p_cell.shape:
        raise ValueError(
            f"cells must be a boolean array with shape (N, T)={p_cell.shape}, "
            f"got {cells_arr.shape} with dtype {cells_arr.dtype}."
        )
    if np.isnan(cells_arr.astype(float)).any():
        raise ValueError("cells contains NA/NaN values.")

    N = p_cell.shape[0]
    result = np.full(N, np.nan, dtype=np.float64)
    for i in range(N):
        sel = cells_arr[i]
        if np.any(sel):
            result[i] = float(np.mean(p_cell[i, sel]))

    return result


def outcome_correlation(model: LongBetMulti) -> np.ndarray:
    """Mean of per-draw implied innovation correlation matrices in user order.

    Parameters
    ----------
    model
        A fitted ``LongBetMulti`` model.

    Returns
    -------
    Float64 matrix of shape ``(M, M)`` in user outcome order with unit diagonal.
    """
    if not isinstance(model, LongBetMulti):
        raise TypeError(f"model must be a LongBetMulti instance, got {type(model)}")
    if model.trace is None:
        raise RuntimeError("Model has not been fitted yet.")

    M = len(model.outcome_names)
    if not model.sur_active:
        return np.eye(M, dtype=np.float64)

    gamma_loadings = np.asarray(model.trace.gamma_loadings, dtype=np.float64)
    # Shape is (*chains, K, M, M) -> flatten chains to (D, M, M)
    D = int(np.prod(gamma_loadings.shape[:-2]))
    gamma_flat = gamma_loadings.reshape(D, M, M)

    # Extract child innovation variances in internal order: each shape (D,)
    sigma2_list = []
    for tr in model.trace.traces:
        s2 = np.asarray(tr.sigma2, dtype=np.float64).reshape(D)
        if np.any(s2 <= 0.0) or not np.all(np.isfinite(s2)):
            raise ValueError("Child sigma2 trace contains non-finite or non-positive values.")
        sigma2_list.append(s2)

    sigma2_arr = np.column_stack(sigma2_list)  # shape (D, M)

    I_mat = np.eye(M, dtype=np.float64)
    corr_sum = np.zeros((M, M), dtype=np.float64)

    for d in range(D):
        G = gamma_flat[d]
        v = sigma2_arr[d]
        A = I_mat - G
        # solve_triangular for lower unit-diagonal A
        B = scipy.linalg.solve_triangular(A, I_mat, lower=True, unit_diagonal=True)
        Sigma = (B * v[None, :]) @ B.T
        sd = np.sqrt(np.diag(Sigma))
        if np.any(sd <= 0.0) or not np.all(np.isfinite(sd)):
            raise ValueError("Encountered non-finite standard deviations in Sigma.")
        R_draw = Sigma / (sd[:, None] * sd[None, :])
        corr_sum += R_draw

    R_internal = corr_sum / D
    # Symmetrize and enforce unit diagonal
    R_internal = 0.5 * (R_internal + R_internal.T)
    np.fill_diagonal(R_internal, 1.0)

    # Permute from internal order to user order: result[inverse_order][:, inverse_order]
    inv = list(model.inverse_order)
    R_user = R_internal[inv][:, inv]
    np.fill_diagonal(R_user, 1.0)
    return R_user
