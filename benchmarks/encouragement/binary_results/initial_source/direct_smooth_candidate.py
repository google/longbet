"""Research Gaussian encouragement forests with direct smooth effect leaves.

This is a deliberately different model from LongBet, not a sampler-only fix:

    standardized_y[i,t] = m(x[i],t) + gamma[i]
                         + (assignment[i] - p) * post[t] * tau(x[i],h[t]) + e[i,t].

Both m and tau are sums of genuinely sampled trees whose leaves contain Gaussian
time vectors. The finite topology space consists of a root-only tree and all
prespecified baseline-covariate stumps. A collapsed *exact* categorical Gibbs
update enumerates that space, then draws the Gaussian leaf conditional. No
adaptive coding, learned GP multiplier, tree-depth approximation to a larger
prior, or sign folding is present. The shallow topology prior is the model.

The two equations are independent conditional on their own parameters; binary
adoption is a Gaussian LPM working response. Neither a successful Geweke check
nor satisfactory MCMC diagnostics establishes coverage under a RED DGP.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json

import numpy as np
from scipy.special import logsumexp
from scipy.linalg import cho_solve, solve_triangular


@dataclass(frozen=True)
class ForestConfig:
    joint_baseline_intercept: bool = False
    joint_effect_leaves: bool = False
    baseline_trees: int = 6
    effect_trees: int = 6
    cutpoints: int = 4
    split_probability: float = .5
    length_scale: float = 2.
    nugget: float = .05
    baseline_sd: float = 1.
    effect_sd: float = 1.
    sigma_shape: float = 3.
    sigma_scale: float = 2.
    gamma_shape: float = 3.
    gamma_scale: float = 2.

    def __post_init__(self):
        if not isinstance(self.joint_baseline_intercept, bool):
            raise ValueError("joint_baseline_intercept must be boolean.")
        if not isinstance(self.joint_effect_leaves, bool):
            raise ValueError("joint_effect_leaves must be boolean.")
        if self.joint_effect_leaves and not self.joint_baseline_intercept:
            raise ValueError("joint_effect_leaves requires joint_baseline_intercept.")
        for name in ("baseline_trees", "effect_trees", "cutpoints"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer.")
        if not 0 < self.split_probability < 1:
            raise ValueError("split_probability must be strictly between zero and one.")
        for name in ("length_scale", "nugget", "baseline_sd", "effect_sd",
                     "sigma_shape", "sigma_scale", "gamma_shape", "gamma_scale"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive.")


@dataclass
class TreeSpace:
    """Finite prior support; mask[rule,leaf,unit], including empty root leaf 1."""

    mask: np.ndarray
    log_prior: np.ndarray
    feature: np.ndarray
    threshold: np.ndarray


def make_tree_space(x, config):
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or not np.isfinite(x).all() or len(x) < 2:
        raise ValueError("x must be a finite baseline covariate matrix.")
    n = len(x)
    left = [np.ones(n, dtype=bool)]
    features, thresholds = [-1], [np.nan]
    for column in range(x.shape[1]):
        cuts = np.unique(np.quantile(x[:, column],
                                    np.arange(1, config.cutpoints + 1)
                                    / (config.cutpoints + 1)))
        for cut in cuts:
            mask = x[:, column] <= cut
            if mask.any() and (~mask).any():
                left.append(mask)
                features.append(column)
                thresholds.append(cut)
    left = np.asarray(left)
    right = ~left
    prior = np.ones(len(left))
    if len(left) > 1:
        prior[0] = 1 - config.split_probability
        prior[1:] = config.split_probability / (len(left) - 1)
    return TreeSpace(np.stack([left, right], axis=1).astype(float), np.log(prior),
                     np.asarray(features), np.asarray(thresholds))


def time_covariance(times, *, sd, length_scale, nugget):
    times = np.asarray(times, dtype=float)
    delta = (times[:, None] - times[None, :]) / length_scale
    correlation = (np.exp(-.5 * delta**2) + nugget * np.eye(len(times))) / (1 + nugget)
    return sd**2 * correlation


def leaf_posterior(total, count, sigma2, covariance):
    """Gaussian leaf posterior for weighted rows with scalar residual variance.

    ``total[h]`` is sum_i w_i r_ih, and ``count`` is sum_i w_i**2.
    The covariance eigenbasis avoids forming a poorly conditioned GP precision.
    """
    eigenvalues, vectors = np.linalg.eigh(covariance)
    variances = eigenvalues * sigma2 / (sigma2 + count * eigenvalues)
    mean = vectors @ (variances * (vectors.T @ total) / sigma2)
    variance = (vectors * variances) @ vectors.T
    return mean, variance


@dataclass
class ForestDesign:
    space: TreeSpace
    weight: np.ndarray
    eigenvalues: np.ndarray
    vectors: np.ndarray
    counts: np.ndarray
    weighted_masks: np.ndarray

    @classmethod
    def build(cls, space, weight, covariance):
        values, vectors = np.linalg.eigh(covariance)
        if np.min(values) <= 0:
            raise ValueError("Leaf covariance must be positive definite.")
        masks = space.mask * np.asarray(weight)[None, None, :]
        return cls(space, np.asarray(weight), values, vectors,
                   np.sum(masks**2, axis=2), masks)

    def collapsed(self, residual, sigma2):
        total = self.weighted_masks @ residual
        projected = total @ self.vectors
        variance = (self.eigenvalues[None, None, :] * sigma2
                    / (sigma2 + self.counts[..., None] * self.eigenvalues))
        log_integral = .5 * np.sum(variance * (projected / sigma2)**2
                                  - np.log1p(self.counts[..., None]
                                             * self.eigenvalues / sigma2), axis=-1)
        # The unused second root leaf integrates to one automatically (count=0).
        score = self.space.log_prior + log_integral.sum(axis=1)
        return score, variance, projected


def sample_tree(rng, design, residual, sigma2):
    """Joint topology/leaf Gibbs draw, with leaves integrated in topology scores."""
    score, variances, projected = design.collapsed(residual, sigma2)
    probability = np.exp(score - logsumexp(score))
    rule = int(rng.choice(len(probability), p=probability))
    eigen_mean = variances[rule] * projected[rule] / sigma2
    leaves = ((eigen_mean + np.sqrt(variances[rule]) * rng.normal(size=eigen_mean.shape))
              @ design.vectors.T)
    # Store the otherwise unused prior draw too: archives represent fixed-size
    # augmentation, including an independent inactive leaf for root-only trees.
    fitted = (design.space.mask[rule].T @ leaves) * design.weight[:, None]
    return rule, leaves, fitted


@dataclass
class GaussianForestState:
    baseline_rules: np.ndarray
    baseline_leaves: np.ndarray
    effect_rules: np.ndarray
    effect_leaves: np.ndarray
    gamma: np.ndarray
    sigma2: float
    gamma2: float
    baseline_fits: np.ndarray
    effect_fits: np.ndarray
    residual: np.ndarray


def fitted(state, start):
    value = state.baseline_fits.sum(axis=0) + state.gamma[:, None]
    value[:, start:] += state.effect_fits.sum(axis=0)
    return value


def prior_state(rng, baseline, effect, config, y, start):
    def forest(design, trees):
        rules = rng.choice(len(design.space.log_prior), size=trees,
                           p=np.exp(design.space.log_prior))
        leaves = ((rng.normal(size=(trees, 2, len(design.eigenvalues)))
                   * np.sqrt(design.eigenvalues)) @ design.vectors.T)
        fits = np.stack([(design.space.mask[rule].T @ leaf) * design.weight[:, None]
                         for rule, leaf in zip(rules, leaves)])
        return rules, leaves, fits

    br, bl, bf = forest(baseline, config.baseline_trees)
    er, el, ef = forest(effect, config.effect_trees)
    sigma2 = config.sigma_scale / rng.gamma(config.sigma_shape)
    gamma2 = config.gamma_scale / rng.gamma(config.gamma_shape)
    gamma = rng.normal(scale=np.sqrt(gamma2), size=len(y))
    state = GaussianForestState(br, bl, er, el, gamma, sigma2, gamma2, bf, ef,
                                np.empty_like(y))
    state.residual = y - fitted(state, start)
    return state


def joint_baseline_intercept_draw(rng, state, baseline, y, start, *, effect=None):
    """Exact joint baseline-leaf/unit-intercept draw conditional on topologies.

    Each leaf is whitened with its fixed Gaussian time covariance. Integrating
    gamma leaves independent within-unit contrasts (variance sigma2) and unit
    mean directions (variance sigma2 + T gamma2). This block removes the
    baseline/intercept Gibbs ridge while leaving the posterior unchanged.
    Inactive leaves have zero design columns and retain their independent prior.
    Supplying ``effect`` adds all encouragement leaves to the same exact block,
    eliminating their shared unit-mean direction with the baseline and gamma.
    """
    n, periods = y.shape
    masks = baseline.space.mask[state.baseline_rules]  # tree, leaf, unit
    factor = baseline.vectors * np.sqrt(baseline.eigenvalues)
    matrix = np.einsum("jln,tq->ntjlq", masks, factor).reshape(n, periods, -1)
    baseline_columns = matrix.shape[-1]
    residual_without_baseline = y.copy()
    if effect is None:
        residual_without_baseline[:, start:] -= state.effect_fits.sum(axis=0)
    else:
        effect_masks = effect.space.mask[state.effect_rules]
        effect_factor = effect.vectors * np.sqrt(effect.eigenvalues)
        effect_matrix = np.zeros((n, periods, effect_masks.shape[0] * 2 * (periods - start)))
        effect_matrix[:, start:] = (np.einsum("jln,tq->ntjlq", effect_masks, effect_factor)
                                     .reshape(n, periods - start, -1) * effect.weight[:, None, None])
        matrix = np.concatenate([matrix, effect_matrix], axis=-1)
    x_mean = matrix.mean(axis=1)
    y_mean = residual_without_baseline.mean(axis=1)
    within = matrix - x_mean[:, None, :]
    within_flat = within.reshape(n * periods, -1)
    between_variance = state.sigma2 + periods * state.gamma2
    precision = (np.eye(matrix.shape[-1]) + within_flat.T @ within_flat / state.sigma2
                 + periods * x_mean.T @ x_mean / between_variance)
    rhs = (np.einsum("ntq,nt->q", within, residual_without_baseline - y_mean[:, None]) / state.sigma2
           + periods * x_mean.T @ y_mean / between_variance)
    chol = np.linalg.cholesky(precision)
    whitened = (cho_solve((chol, True), rhs, check_finite=False)
                + solve_triangular(chol.T, rng.normal(size=len(rhs)), lower=False, check_finite=False))
    state.baseline_leaves = whitened[:baseline_columns].reshape(len(masks), 2, periods) @ factor.T
    state.baseline_fits = np.einsum("jln,jlt->jnt", masks, state.baseline_leaves)
    if effect is not None:
        state.effect_leaves = (whitened[baseline_columns:].reshape(len(effect_masks), 2, periods - start)
                                @ effect_factor.T)
        state.effect_fits = (np.einsum("jln,jlt->jnt", effect_masks, state.effect_leaves)
                            * effect.weight[None, :, None])
        residual_without_baseline[:, start:] -= state.effect_fits.sum(axis=0)
    residual_without_gamma = residual_without_baseline - state.baseline_fits.sum(axis=0)
    variance = state.sigma2 * state.gamma2 / between_variance
    mean = state.gamma2 * residual_without_gamma.sum(axis=1) / between_variance
    state.gamma = mean + np.sqrt(variance) * rng.normal(size=n)
    state.residual = residual_without_gamma - state.gamma[:, None]


def sweep(rng, state, baseline, effect, config, y, start):
    """Systematic partially collapsed Gibbs sweep preserving the full residual."""
    for j in range(config.baseline_trees):
        state.residual += state.baseline_fits[j]
        rule, leaves, fit = sample_tree(rng, baseline, state.residual, state.sigma2)
        state.baseline_rules[j], state.baseline_leaves[j] = rule, leaves
        state.baseline_fits[j] = fit
        state.residual -= fit
    for j in range(config.effect_trees):
        state.residual[:, start:] += state.effect_fits[j]
        rule, leaves, fit = sample_tree(rng, effect, state.residual[:, start:], state.sigma2)
        state.effect_rules[j], state.effect_leaves[j] = rule, leaves
        state.effect_fits[j] = fit
        state.residual[:, start:] -= fit
    if config.joint_baseline_intercept:
        if config.joint_effect_leaves:
            joint_baseline_intercept_draw(rng, state, baseline, y, start, effect=effect)
        else:
            joint_baseline_intercept_draw(rng, state, baseline, y, start)
    else:
        state.residual += state.gamma[:, None]
        variance = 1 / (1 / state.gamma2 + y.shape[1] / state.sigma2)
        mean = variance * state.residual.sum(axis=1) / state.sigma2
        state.gamma = mean + np.sqrt(variance) * rng.normal(size=len(y))
        state.residual -= state.gamma[:, None]
    state.gamma2 = ((config.gamma_scale + .5 * np.sum(state.gamma**2))
                    / rng.gamma(config.gamma_shape + len(y) / 2))
    state.sigma2 = ((config.sigma_scale + .5 * np.sum(state.residual**2))
                    / rng.gamma(config.sigma_shape + y.size / 2))
    return state


def make_designs(x, assignment, times, start, config):
    space = make_tree_space(x, config)
    base_cov = time_covariance(times, sd=config.baseline_sd / np.sqrt(config.baseline_trees),
                               length_scale=config.length_scale, nugget=config.nugget)
    effect_cov = time_covariance(times[start:], sd=config.effect_sd / np.sqrt(config.effect_trees),
                                 length_scale=config.length_scale, nugget=config.nugget)
    centered = np.asarray(assignment, dtype=float) - np.mean(assignment)
    return (ForestDesign.build(space, np.ones(len(x)), base_cov),
            ForestDesign.build(space, centered, effect_cov))


def _validate_data(data, start):
    y, d, x = (np.asarray(data[name], dtype=float) for name in ("y", "d", "x"))
    assignment = np.asarray(data["assignment"], dtype=float)
    if y.ndim != 2 or d.shape != y.shape or x.ndim != 2 or len(x) != len(y):
        raise ValueError("Require aligned complete outcome/adoption panels and baseline x.")
    if not all(np.isfinite(value).all() for value in (y, d, x, assignment)):
        raise ValueError("Candidate inputs must be finite.")
    if assignment.shape != (len(y),) or not np.isin(assignment, (0, 1)).all():
        raise ValueError("assignment must be a binary unit vector.")
    if assignment.sum() in (0, len(y)):
        raise ValueError("Both randomized assignment arms must be present.")
    if isinstance(start, bool) or not isinstance(start, (int, np.integer)) or not 0 <= start < y.shape[1]:
        raise ValueError("start must index the first post-encouragement period.")
    if not np.isin(d, (0, 1)).all() or np.any(np.diff(d, axis=1) < 0):
        raise ValueError("Adoption must be a complete absorbing binary panel.")
    if "z" in data:
        expected = np.zeros_like(y)
        expected[:, start:] = assignment[:, None]
        if not np.array_equal(np.asarray(data["z"]), expected):
            raise ValueError("Candidate supports a single randomized encouragement wave.")
    times = np.asarray(data.get("t", np.arange(y.shape[1])), dtype=float)
    if times.shape != (y.shape[1],) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("t must contain strictly increasing finite observation times.")
    return np.stack([y, d], axis=-1), x, assignment, times


def _fingerprint(*arrays):
    digest = hashlib.sha256()
    for value in arrays:
        value = np.ascontiguousarray(value, dtype=float)
        digest.update(str(value.shape).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def fit_candidate(data, start, *, seed, chains=4, burnin=300, draws=300,
                  config=None, resume=None):
    """Fit independent direct forests; return natural-unit ITTs and audit arrays.

    Effect arrays have shape ``(1, horizon, chain, draw)`` and average over all
    original units. ``resume`` may be the previous audit dictionary/loaded NPZ;
    use burnin=0 to append draws from its exact saved states and RNG streams.
    No object arrays or pickle are used. Diagnosing the returned parameters and
    effects, and combining old/new draws on resume, belongs to the study runner.
    """
    config = config or ForestConfig()
    for name, value, lower in (("chains", chains, 1), ("burnin", burnin, 0), ("draws", draws, 1)):
        if isinstance(value, bool) or not isinstance(value, int) or value < lower:
            raise ValueError(f"{name} must be an integer at least {lower}.")
    raw, x, assignment, times = _validate_data(data, start)
    n, periods, outcomes = raw.shape
    center = raw.mean(axis=(0, 1))
    scale = raw.std(axis=(0, 1))
    scale = np.where(scale > 1e-10, scale, 1.)
    standardized = (raw - center) / scale
    baseline, effect = make_designs(x, assignment, times, start, config)
    h = periods - start
    settings = json.dumps(asdict(config), sort_keys=True)
    fingerprint = _fingerprint(raw, x, assignment, times, np.asarray(start))
    probe_indices = np.unique([0, n // 2, n - 1])
    archive = dict(
        baseline_rules=np.empty((chains, draws, outcomes, config.baseline_trees), dtype=np.int32),
        baseline_leaves=np.empty((chains, draws, outcomes, config.baseline_trees, 2, periods)),
        effect_rules=np.empty((chains, draws, outcomes, config.effect_trees), dtype=np.int32),
        effect_leaves=np.empty((chains, draws, outcomes, config.effect_trees, 2, h)),
        gamma=np.empty((chains, draws, n, outcomes)),
        sigma2=np.empty((chains, draws, outcomes)), gamma2=np.empty((chains, draws, outcomes)),
        itt=np.empty((chains, draws, h, outcomes)),
        baseline_probe=np.empty((chains, draws, len(probe_indices), periods, outcomes)),
        effect_probe=np.empty((chains, draws, len(probe_indices), h, outcomes)),
        probe_indices=probe_indices,
        center=center, scale=scale, config_json=np.asarray(settings),
        fingerprint=np.asarray(fingerprint), start=np.asarray(start),
        rule_feature=baseline.space.feature, rule_threshold=baseline.space.threshold,
        format_version=np.asarray(1), seed=np.asarray(seed),
        assignment_probability=np.asarray(assignment.mean()),
    )
    if resume is not None:
        saved_settings = json.dumps(asdict(ForestConfig(**json.loads(str(resume["config_json"])))),
                                    sort_keys=True)
        if burnin != 0 or saved_settings != settings or str(resume["fingerprint"]) != fingerprint:
            raise ValueError("Resume requires burnin=0 and identical data/configuration.")
        if np.asarray(resume["sigma2"]).shape[0] != chains:
            raise ValueError("Resume chain count must match the saved chains.")
        if int(resume["format_version"]) != 1:
            raise ValueError("Unsupported candidate archive version.")
    rng_states = []
    for c, child in enumerate(np.random.SeedSequence(seed).spawn(chains)):
        rng = np.random.default_rng(child)
        states = []
        for outcome in range(outcomes):
            state = prior_state(rng, baseline, effect, config, standardized[..., outcome], start)
            if resume is not None:
                state.baseline_rules = np.array(resume["baseline_rules"][c, -1, outcome], copy=True)
                state.baseline_leaves = np.array(resume["baseline_leaves"][c, -1, outcome], copy=True)
                state.effect_rules = np.array(resume["effect_rules"][c, -1, outcome], copy=True)
                state.effect_leaves = np.array(resume["effect_leaves"][c, -1, outcome], copy=True)
                state.gamma = np.array(resume["gamma"][c, -1, :, outcome], copy=True)
                state.sigma2 = float(resume["sigma2"][c, -1, outcome])
                state.gamma2 = float(resume["gamma2"][c, -1, outcome])
                for design, rules, leaves, fits in (
                        (baseline, state.baseline_rules, state.baseline_leaves, state.baseline_fits),
                        (effect, state.effect_rules, state.effect_leaves, state.effect_fits)):
                    for j, (rule, leaf) in enumerate(zip(rules, leaves)):
                        fits[j] = (design.space.mask[rule].T @ leaf) * design.weight[:, None]
                state.residual = standardized[..., outcome] - fitted(state, start)
            states.append(state)
        if resume is not None:
            rng.bit_generator.state = json.loads(str(resume["rng_state_json"][c]))
        for iteration in range(burnin + draws):
            for outcome, state in enumerate(states):
                sweep(rng, state, baseline, effect, config, standardized[..., outcome], start)
                if iteration < burnin:
                    continue
                k = iteration - burnin
                for name in ("baseline_rules", "baseline_leaves", "effect_rules", "effect_leaves"):
                    archive[name][c, k, outcome] = getattr(state, name)
                archive["gamma"][c, k, :, outcome] = state.gamma
                archive["sigma2"][c, k, outcome] = state.sigma2
                archive["gamma2"][c, k, outcome] = state.gamma2
                # Contrast A=1 versus A=0 is the unweighted tau forest. Gamma,
                # baseline, and centering probability cancel for each unit.
                tau = sum(effect.space.mask[rule].T @ leaf
                          for rule, leaf in zip(state.effect_rules, state.effect_leaves))
                archive["itt"][c, k, :, outcome] = tau.mean(axis=0) * scale[outcome]
                archive["baseline_probe"][c, k, :, :, outcome] = (
                    state.baseline_fits.sum(axis=0)[probe_indices] * scale[outcome] + center[outcome])
                archive["effect_probe"][c, k, :, :, outcome] = tau[probe_indices] * scale[outcome]
        rng_states.append(json.dumps(rng.bit_generator.state, sort_keys=True))
    archive["rng_state_json"] = np.asarray(rng_states)
    values = {name: archive["itt"][..., j].transpose(2, 0, 1)[None]
              for j, name in enumerate(("outcome", "takeup"))}
    return values, archive


def geweke(replications=3000, seed=93118, *, joint_baseline_intercept=False,
           joint_effect_leaves=False):
    """Independent prior/data/one-Gibbs check, with analytic prior moments."""
    rng = np.random.default_rng(seed)
    config = ForestConfig(baseline_trees=2, effect_trees=2, cutpoints=2,
                          joint_baseline_intercept=joint_baseline_intercept,
                          joint_effect_leaves=joint_effect_leaves)
    n, periods, start = 8, 4, 1
    x = np.arange(n)[:, None]
    assignment = np.arange(n) % 2
    baseline, effect = make_designs(x, assignment, np.arange(periods), start, config)
    values = []
    for _ in range(replications):
        state = prior_state(rng, baseline, effect, config, np.zeros((n, periods)), start)
        y = fitted(state, start) + rng.normal(size=(n, periods)) * np.sqrt(state.sigma2)
        state.residual = y - fitted(state, start)
        sweep(rng, state, baseline, effect, config, y, start)
        tau = sum(effect.space.mask[r].T @ leaf for r, leaf in
                  zip(state.effect_rules, state.effect_leaves))
        mu = state.baseline_fits.sum(axis=0)
        values.append([mu.mean(), np.mean(mu**2), tau.mean(), np.mean(tau**2),
                       np.mean(state.gamma**2), state.gamma2, state.sigma2,
                       np.mean(state.baseline_rules == 0), np.mean(state.effect_rules == 0)])
    values = np.asarray(values)
    expectation = np.array([0., config.baseline_sd**2, 0., config.effect_sd**2,
                            config.gamma_scale / (config.gamma_shape - 1),
                            config.gamma_scale / (config.gamma_shape - 1),
                            config.sigma_scale / (config.sigma_shape - 1),
                            1 - config.split_probability, 1 - config.split_probability])
    error = values.mean(axis=0) - expectation
    mcse = values.std(axis=0, ddof=1) / np.sqrt(replications)
    z = np.divide(error, mcse, out=np.zeros_like(error), where=mcse > 0)
    z[(mcse == 0) & (error != 0)] = np.copysign(np.inf, error[(mcse == 0) & (error != 0)])
    return dict(method="independent prior/data/one-Gibbs transitions", replications=replications,
                seed=seed, functions=["baseline_mean", "baseline_second_moment", "effect_mean",
                                     "effect_second_moment", "gamma_second_moment", "gamma_variance",
                                     "innovation_variance", "baseline_root", "effect_root"],
                expectation=expectation, estimate=values.mean(axis=0), mcse=mcse, z=z,
                passed=bool(np.all(np.abs(z) < 5)))
