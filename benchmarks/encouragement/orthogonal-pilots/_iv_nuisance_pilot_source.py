"""Experimental LongBet predictions for externally corrected IV estimation.

This learner reuses the direct LongBet Gaussian-process leaf forest and Gibbs
updates.  It deliberately exposes predictions, not causal posterior intervals.
An account's full post-assignment trajectory is one training observation;
``x`` must contain only covariates measured before assignment.  New-account
random intercepts are integrated out rather than estimated from test outcomes.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from typing import Any

import numpy as np

from longbet._direct_smooth import (
    DirectSmoothConfig,
    ForestDesign,
    SingleEquationState,
    TreeSpace,
    sample_tree,
    time_covariance,
)


@dataclass(frozen=True)
class LongBetIVNuisanceConfig:
    """Small prediction-only MCMC budget, selected before evaluation.

    Length scales are in units of the median spacing between post periods.
    ``0`` means independent time leaves (no temporal pooling).  Selection uses
    one assignment-stratified holdout entirely inside each training fold.
    """

    baseline_trees: int = 4
    effect_trees: int = 4
    cutpoints: int = 4
    burnin: int = 100
    draws: int = 200
    chains: int = 2
    length_scales: tuple[float, ...] = (0.0, 2.0, 8.0)
    validation_fraction: float = 0.2
    nugget: float = 0.1
    interaction_partitions: bool = True
    max_interaction_rules: int = 64
    seed: int = 0

    def __post_init__(self) -> None:
        for name, lower in (("baseline_trees", 1), ("effect_trees", 1),
                            ("cutpoints", 1), ("burnin", 0), ("draws", 1), ("chains", 1),
                            ("max_interaction_rules", 0),
                            ("seed", 0)):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < lower:
                raise ValueError(f"{name} must be an integer >= {lower}.")
        if not self.length_scales or any(not np.isfinite(v) or v < 0 for v in self.length_scales):
            raise ValueError("length_scales must contain finite nonnegative values.")
        if len(set(self.length_scales)) != len(self.length_scales):
            raise ValueError("length_scales must not contain duplicates.")
        if not np.isfinite(self.validation_fraction) or not 0 < self.validation_fraction < 0.5:
            raise ValueError("validation_fraction must be strictly between 0 and 0.5.")
        if not np.isfinite(self.nugget) or self.nugget <= 0:
            raise ValueError("nugget must be finite and positive.")
        if not isinstance(self.interaction_partitions, bool):
            raise ValueError("interaction_partitions must be boolean.")


def inner_validation_split(assignment: Any, seed: int,
                           validation_fraction: float = 0.2) -> tuple[np.ndarray, np.ndarray]:
    """Shared account-level tuning split; no response or covariate is inspected."""
    z = np.asarray(assignment)
    if z.ndim != 1 or not np.isin(z, [0, 1]).all():
        raise ValueError("assignment must be a binary vector.")
    if not np.isfinite(validation_fraction) or not 0 < validation_fraction < 0.5:
        raise ValueError("validation_fraction must be strictly between 0 and 0.5.")
    rng = np.random.default_rng(seed)
    validation = []
    for arm in (0, 1):
        indices = np.flatnonzero(z == arm)
        if len(indices) < 6:
            raise ValueError("Adaptive tuning requires at least six training accounts per assignment arm.")
        n_valid = min(len(indices) - 4, max(2, int(np.floor(len(indices) * validation_fraction))))
        validation.extend(rng.permutation(indices)[:n_valid])
    validation = np.sort(validation).astype(int)
    train = np.ones(len(z), dtype=bool)
    train[validation] = False
    return np.flatnonzero(train), validation


def nuisance_validation_score(prediction: Any, assignment: Any, responses: Any,
                              training_responses: Any) -> float:
    """Observed-arm mean squared loss, scaled using inner-training data only."""
    prediction = np.asarray(prediction, dtype=float)
    responses = np.asarray(responses, dtype=float)
    z = np.asarray(assignment, dtype=int)
    scale = np.asarray(training_responses, dtype=float).std(axis=(0, 1))
    scale = np.maximum(scale, np.array([1e-8, 0.05]))
    observed = prediction[np.arange(len(z)), :, z, :]
    return float(np.mean(((observed - responses) / scale)**2))


def _training_arrays(x: Any, assignment: Any, responses: Any, times: Any):
    x = np.asarray(x, dtype=float)
    z = np.asarray(assignment, dtype=float)
    responses = np.asarray(responses, dtype=float)
    times = np.asarray(times, dtype=float)
    if x.ndim != 2 or len(x) < 4 or not np.isfinite(x).all():
        raise ValueError("x_train must be a finite matrix with at least four accounts.")
    if z.shape != (len(x),) or not np.isin(z, [0, 1]).all():
        raise ValueError("assignment_train must be a binary vector matching x_train.")
    if min(np.bincount(z.astype(int), minlength=2)) < 2:
        raise ValueError("Each assignment arm must contain at least two training accounts.")
    if (responses.ndim != 3 or responses.shape[0] != len(x) or responses.shape[2] != 2
            or responses.shape[1] < 1 or not np.isfinite(responses).all()):
        raise ValueError("responses_train must be finite with shape (accounts, horizons, 2).")
    if not np.isin(responses[..., 1], [0, 1]).all():
        raise ValueError("The second response must be binary adoption.")
    if (times.shape != (responses.shape[1],) or not np.isfinite(times).all()
            or (len(times) > 1 and np.any(np.diff(times) <= 0))):
        raise ValueError("times_post must be finite, strictly increasing, and match the horizons.")
    return x, z, responses, times


def _rule_masks(x: np.ndarray, feature: np.ndarray, threshold: np.ndarray) -> np.ndarray:
    """Return the frozen training rules evaluated on new baseline covariates."""
    masks = np.zeros((len(feature), 4, len(x)))
    for rule, (columns, cuts) in enumerate(zip(feature, threshold)):
        if columns[0] < 0:
            masks[rule, 0] = 1
            continue
        left = x[:, columns[0]] <= cuts[0]
        if columns[1] < 0:
            masks[rule, 0], masks[rule, 1] = left, ~left
            continue
        bottom = x[:, columns[1]] <= cuts[1]
        masks[rule] = np.stack([left & bottom, left & ~bottom, ~left & bottom, ~left & ~bottom])
    return masks


def _prediction_tree_space(x: np.ndarray, config: LongBetIVNuisanceConfig) -> TreeSpace:
    """Finite prior over a root, stumps, and crossed-median depth-two trees.

    All cutpoints and admissible interactions depend on training covariates
    only.  When more than the rule cap are admissible, a seeded sample selects
    the finite support; Gibbs sampling subsequently updates the topology.
    """
    features = [(-1, -1)]
    thresholds = [(np.nan, np.nan)]
    for column in range(x.shape[1]):
        cuts = np.unique(np.quantile(x[:, column], np.arange(1, config.cutpoints + 1) / (config.cutpoints + 1)))
        for cut in cuts:
            left = x[:, column] <= cut
            if left.any() and (~left).any():
                features.append((column, -1))
                thresholds.append((float(cut), np.nan))
    n_stumps = len(features) - 1
    interactions = []
    if config.interaction_partitions:
        medians = np.median(x, axis=0)
        for first in range(x.shape[1]):
            for second in range(first + 1, x.shape[1]):
                candidate = _rule_masks(x, np.array([[first, second]]), np.array([[medians[first], medians[second]]]))
                if np.min(candidate.sum(axis=-1)) >= 2:
                    interactions.append(((first, second), (medians[first], medians[second])))
    if len(interactions) > config.max_interaction_rules:
        choose = np.sort(np.random.default_rng(config.seed + 63443).choice(
            len(interactions), config.max_interaction_rules, replace=False))
        interactions = [interactions[index] for index in choose]
    for columns, cuts in interactions:
        features.append(columns)
        thresholds.append(cuts)
    n_interactions = len(interactions)
    # Reserve half the prior for the constant topology and divide the other
    # half equally among whichever nonconstant topology families are present.
    probability = np.ones(len(features))
    families = int(n_stumps > 0) + int(n_interactions > 0)
    if families:
        probability[0] = 0.5
        if n_stumps:
            probability[1:1 + n_stumps] = 0.5 / families / n_stumps
        if n_interactions:
            probability[1 + n_stumps:] = 0.5 / families / n_interactions
    features = np.asarray(features)
    thresholds = np.asarray(thresholds)
    return TreeSpace(_rule_masks(x, features, thresholds), np.log(probability), features, thresholds)


def _initial_state(rng: np.random.Generator, baseline: ForestDesign,
                   effect: ForestDesign, config: DirectSmoothConfig,
                   y: np.ndarray) -> SingleEquationState:
    """DirectSmooth initial state generalized to four-leaf finite trees."""
    def forest(design: ForestDesign, trees: int):
        rules = rng.choice(len(design.space.log_prior), size=trees, p=np.exp(design.space.log_prior))
        leaves = ((rng.normal(size=(trees, design.space.mask.shape[1], len(design.eigenvalues)))
                   * np.sqrt(design.eigenvalues)) @ design.vectors.T)
        fits = np.einsum("jln,jlt->jnt", design.space.mask[rules], leaves) * design.weight[None, :, None]
        return rules, leaves, fits
    br, bl, bf = forest(baseline, config.baseline_trees)
    er, el, ef = forest(effect, config.effect_trees)
    sigma2 = config.sigma_scale / rng.gamma(config.sigma_shape)
    gamma2 = config.gamma_scale / rng.gamma(config.gamma_shape)
    gamma = rng.normal(scale=np.sqrt(gamma2), size=len(y))
    residual = y - bf.sum(axis=0) - ef.sum(axis=0) - gamma[:, None]
    return SingleEquationState(br, bl, er, el, gamma, sigma2, gamma2, bf, ef, residual)


def _intercept_whitening(sigma2: float, gamma2: float, periods: int) -> tuple[np.ndarray, np.ndarray]:
    """Symmetric square root and inverse for sigma² I + gamma² 11'."""
    within_sd = np.sqrt(sigma2)
    mean_sd = np.sqrt(sigma2 + periods * gamma2)
    inverse = np.eye(periods) / within_sd + np.ones((periods, periods)) * ((1 / mean_sd - 1 / within_sd) / periods)
    root = np.eye(periods) * within_sd + np.ones((periods, periods)) * ((mean_sd - within_sd) / periods)
    return inverse, root


def _marginal_forest_sweep(rng: np.random.Generator, state: SingleEquationState,
                         baseline: ForestDesign, effect: ForestDesign,
                         baseline_covariance: np.ndarray, effect_covariance: np.ndarray,
                         config: DirectSmoothConfig, y: np.ndarray) -> None:
    """Exact collapsed random-intercept update of the LongBet leaf forests.

    Tree topologies and leaves condition on the variances with account
    intercepts integrated out.  Intercepts are then redrawn before either
    variance update, so no step conditions on stale marginalized variables.
    """
    whitening, unwhitening = _intercept_whitening(state.sigma2, state.gamma2, y.shape[1])
    baseline_white = ForestDesign.build(baseline.space, baseline.weight,
                                       whitening @ baseline_covariance @ whitening)
    effect_white = ForestDesign.build(effect.space, effect.weight,
                                     whitening @ effect_covariance @ whitening)
    baseline_fits = state.baseline_fits @ whitening
    effect_fits = state.effect_fits @ whitening
    residual = y @ whitening - baseline_fits.sum(axis=0) - effect_fits.sum(axis=0)
    for index in range(config.baseline_trees):
        residual += baseline_fits[index]
        rule, leaves, fit = sample_tree(rng, baseline_white, residual, 1.0)
        state.baseline_rules[index] = rule
        state.baseline_leaves[index] = leaves @ unwhitening
        baseline_fits[index] = fit
        residual -= fit
    for index in range(config.effect_trees):
        residual += effect_fits[index]
        rule, leaves, fit = sample_tree(rng, effect_white, residual, 1.0)
        state.effect_rules[index] = rule
        state.effect_leaves[index] = leaves @ unwhitening
        effect_fits[index] = fit
        residual -= fit
    state.baseline_fits = baseline_fits @ unwhitening
    state.effect_fits = effect_fits @ unwhitening
    residual_raw = residual @ unwhitening
    variance = 1 / (1 / state.gamma2 + y.shape[1] / state.sigma2)
    mean = variance * residual_raw.sum(axis=1) / state.sigma2
    state.gamma = mean + np.sqrt(variance) * rng.normal(size=len(y))
    state.residual = residual_raw - state.gamma[:, None]
    state.gamma2 = ((config.gamma_scale + 0.5 * np.sum(state.gamma**2))
                    / rng.gamma(config.gamma_shape + len(y) / 2))
    state.sigma2 = ((config.sigma_scale + 0.5 * np.sum(state.residual**2))
                    / rng.gamma(config.sigma_shape + y.size / 2))


def _evaluate_forest(masks: np.ndarray, rule_leaves: np.ndarray) -> np.ndarray:
    """Posterior-mean forest prediction from accumulated rule/leaf coefficients."""
    return np.einsum("rln,rlhk->nhk", masks, rule_leaves)


class LongBetIVNuisance:
    """LongBet nuisance predictions for both randomized-assignment arms.

    ``fit(x, assignment, responses, times_post)`` consumes responses with shape
    ``(accounts, horizons, 2)`` in outcome/adoption order.  ``predict(x_new)``
    returns ``(accounts, horizons, 2 assignment arms, 2 responses)``.  The uptake
    equation remains a Gaussian working model, with posterior-mean predictions
    clipped to [0, 1].  Predictions do not assert adoption-hazard coherence.

    An external cross-fitting estimator must exclude every held-out account's
    trajectory from training.  This class cannot determine whether caller-
    supplied covariates were measured before assignment.
    """

    def __init__(self, config: LongBetIVNuisanceConfig | None = None, **kwargs: Any):
        if config is not None and (not isinstance(config, LongBetIVNuisanceConfig) or kwargs):
            raise TypeError("Pass either LongBetIVNuisanceConfig or configuration keyword arguments.")
        self.config = config if config is not None else LongBetIVNuisanceConfig(**kwargs)
        self.fitted_ = False
        self.metadata_: dict[str, Any] = {}

    def fit(self, x_train: Any, assignment_train: Any, responses_train: Any,
            times_post: Any) -> LongBetIVNuisance:
        x, z, responses, times = _training_arrays(x_train, assignment_train, responses_train, times_post)
        self.fitted_ = False
        scores: list[dict[str, float]] = []
        validation_indices = np.empty(0, dtype=int)
        selected = float(self.config.length_scales[0])
        if len(self.config.length_scales) > 1:
            train, validation_indices = inner_validation_split(z, self.config.seed, self.config.validation_fraction)
            for length_scale in self.config.length_scales:
                candidate = self._fit_fixed(x[train], z[train], responses[train], times,
                                            float(length_scale), self.config.seed + 12553)
                prediction = self._predict_fitted(x[validation_indices], candidate)
                score = nuisance_validation_score(prediction, z[validation_indices], responses[validation_indices], responses[train])
                scores.append({"length_scale": float(length_scale), "validation_mse": score})
            selected = scores[int(np.argmin([item["validation_mse"] for item in scores]))]["length_scale"]
        self._fitted = self._fit_fixed(x, z, responses, times, selected, self.config.seed)
        self.n_features_in_ = x.shape[1]
        self.selected_length_scale_ = selected
        self.metadata_ = {
            "learner": "longbet_direct_gp_finite_depth_two_forest",
            "role": "experimental_nuisance_prediction_only",
            "config": asdict(self.config),
            "selected_length_scale": selected,
            "length_scale_units": "median_post_period_spacing",
            "validation_scores": scores,
            "inner_validation_indices": validation_indices.tolist(),
            "inner_holdout_sha256": hashlib.sha256(np.asarray(validation_indices, dtype="<i8").tobytes()).hexdigest(),
            "inner_training_count": len(x) - len(validation_indices),
            "inner_validation_count": len(validation_indices),
            "new_unit_intercept": "integrated_prior_mean_zero",
            "adoption_prediction": "clipped_gaussian_working_mean",
            "posterior_interval_inference": False,
            "n_training_accounts": len(x),
            "n_horizons": len(times),
            "n_finite_tree_rules": len(self._fitted["feature"]),
            "prediction_chain_rms_discrepancy": self._fitted["chain_rms_discrepancy"].tolist(),
            "prediction_chain_rms_scale": "training_factual_predictions_per_response_training_sd",
            "random_intercept_sampling": "integrated_during_tree_updates_refreshed_before_variances",
        }
        self.fitted_ = True
        return self

    def _fit_fixed(self, x: np.ndarray, z: np.ndarray, raw: np.ndarray, times: np.ndarray,
                   length_scale: float, seed: int) -> dict[str, Any]:
        config = DirectSmoothConfig(
            baseline_trees=self.config.baseline_trees, effect_trees=self.config.effect_trees,
            cutpoints=self.config.cutpoints,
            joint_baseline_intercept=False,
            joint_effect_leaves=False,
            nugget=self.config.nugget,
        )
        center = raw.mean(axis=(0, 1))
        scale = raw.std(axis=(0, 1))
        scale = np.where(scale > 1e-10, scale, 1.0)
        standardized = (raw - center) / scale
        space = _prediction_tree_space(x, self.config)
        p = float(z.mean())
        normalized_times = (times - times[0]) / (np.median(np.diff(times)) if len(times) > 1 else 1.0)
        def covariance(trees: int) -> np.ndarray:
            if length_scale == 0:
                return np.eye(len(times)) / trees
            return time_covariance(normalized_times, sd=1 / np.sqrt(trees),
                                   length_scale=length_scale, nugget=config.nugget)
        baseline_covariance = covariance(config.baseline_trees)
        effect_covariance = covariance(config.effect_trees)
        baseline = ForestDesign.build(space, np.ones(len(x)), baseline_covariance)
        effect = ForestDesign.build(space, z - p, effect_covariance)
        chain_baseline = np.zeros((self.config.chains, len(space.feature), 4, len(times), 2))
        chain_effect = np.zeros_like(chain_baseline)
        for chain in range(self.config.chains):
            rng = np.random.default_rng(seed + 10007 * chain)
            states = [_initial_state(rng, baseline, effect, config, standardized[..., k]) for k in range(2)]
            for iteration in range(self.config.burnin + self.config.draws):
                for equation, state in enumerate(states):
                    _marginal_forest_sweep(rng, state, baseline, effect, baseline_covariance,
                                           effect_covariance, config, standardized[..., equation])
                    if iteration >= self.config.burnin:
                        np.add.at(chain_baseline[chain, ..., equation], state.baseline_rules, state.baseline_leaves)
                        np.add.at(chain_effect[chain, ..., equation], state.effect_rules, state.effect_leaves)
        chain_baseline /= self.config.draws
        chain_effect /= self.config.draws
        # Prediction disagreement is an algorithmic stability indicator, not
        # an MCMC convergence certificate or an inferential standard error.
        chain_fits = np.stack([_evaluate_forest(space.mask, b) +
                               _evaluate_forest(space.mask, e) * (z - p)[:, None, None]
                               for b, e in zip(chain_baseline, chain_effect)])
        chain_rms = np.sqrt(np.mean((chain_fits - chain_fits.mean(axis=0))**2, axis=(0, 1, 2)))
        return {
            "feature": space.feature.copy(), "threshold": space.threshold.copy(),
            "baseline_leaves": chain_baseline.mean(axis=0),
            "effect_leaves": chain_effect.mean(axis=0),
            "assignment_fraction": p, "center": center, "scale": scale,
            "chain_rms_discrepancy": chain_rms,
        }

    @staticmethod
    def _predict_fitted(x: np.ndarray, fitted: dict[str, Any]) -> np.ndarray:
        masks = _rule_masks(x, fitted["feature"], fitted["threshold"])
        baseline = _evaluate_forest(masks, fitted["baseline_leaves"])
        effect = _evaluate_forest(masks, fitted["effect_leaves"])
        centered_assignment = np.arange(2) - fitted["assignment_fraction"]
        prediction = ((baseline[:, :, None, :] + effect[:, :, None, :] * centered_assignment[None, None, :, None])
                      * fitted["scale"] + fitted["center"])
        prediction[..., 1] = np.clip(prediction[..., 1], 0.0, 1.0)
        return prediction

    def predict(self, x_test: Any) -> np.ndarray:
        if not self.fitted_:
            raise RuntimeError("Fit the nuisance learner before prediction.")
        x = np.asarray(x_test, dtype=float)
        if x.ndim != 2 or x.shape[1] != self.n_features_in_ or not np.isfinite(x).all():
            raise ValueError("x_test must be finite and match the training feature count.")
        return self._predict_fitted(x, self._fitted)


__all__ = ["LongBetIVNuisance", "LongBetIVNuisanceConfig", "inner_validation_split", "nuisance_validation_score"]
