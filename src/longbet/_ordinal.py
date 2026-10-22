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

"""Ordered-probit preparation, interval Gibbs draws, and probability transforms.

Host validation and prediction use float64 NumPy/SciPy. Sampling is entirely
JAX float32. The rejection kernel avoids CDF interpolation, which degenerates
in float32 normal tails (including [8, 9] on JAX 0.11.1).
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import equinox as eqx
import jax
import jax.numpy as jnp
import numpy as np
from jax import lax, random
from jax.scipy.special import log_ndtr as jax_log_ndtr
from scipy.special import log_ndtr, ndtri


@dataclass(frozen=True)
class OrdinalPreparation:
    labels: np.ndarray
    obs_mask: np.ndarray
    offset: float
    cutpoints: np.ndarray


def prepare_ordinal(y, num_categories) -> OrdinalPreparation:
    """Validate declared labels without recoding and prepare starting thresholds."""
    if (isinstance(num_categories, (bool, np.bool_))
            or not isinstance(num_categories, Integral) or num_categories < 2):
        raise ValueError("num_categories must be an integer >=2 for ordinal outcomes")
    K = int(num_categories)
    values = np.asarray(y)
    if values.dtype.kind not in "biuf":
        raise ValueError("ordinal y must contain numeric integer category labels")
    values = values.astype(np.float64)
    if np.isinf(values).any():
        raise ValueError("ordinal y cannot contain infinity; only NaN means missing")
    observed = ~np.isnan(values)
    if not observed.any():
        raise ValueError("y contains no observed cells")
    obs_y = values[observed]
    if np.any((obs_y != np.floor(obs_y)) | (obs_y < 0) | (obs_y >= K)):
        raise ValueError(f"ordinal y requires integer labels in [0, {K - 1}]")
    labels = np.where(observed, values, 0).astype(np.float32)
    if K == 2:
        # Same rate calculation and clipping as the legacy binary fit.
        offset = float(ndtri(np.clip(float(obs_y.mean()), 1e-4, 1 - 1e-4)))
        cutpoints = np.empty(0, np.float32)
    else:
        counts = np.bincount(obs_y.astype(np.int64), minlength=K)
        q = ndtri(np.clip(np.cumsum(counts)[:-1] / obs_y.size, 1e-4, 1 - 1e-4))
        offset = float(-q[0])
        cutpoints = np.empty(K - 2, np.float64)
        prev = 0.0
        for i in range(K - 2):
            prev = max(float(q[i + 1] - q[0]), prev + 0.1)
            cutpoints[i] = prev
        cutpoints = cutpoints.astype(np.float32)
        if (not np.isfinite(cutpoints).all()
                or not (np.diff(np.r_[np.float32(0), cutpoints]) > 0).all()):
            raise ValueError("initial ordinal cutpoints are not finite and strictly ordered")
    return OrdinalPreparation(labels, observed, offset, cutpoints)


def _interior(value, lower, upper):
    """Correct endpoint rounding only; reject unrepresentable intervals."""
    lo = jnp.nextafter(lower, jnp.float32(jnp.inf))
    hi = jnp.nextafter(upper, jnp.float32(-jnp.inf))
    valid = (lower < upper) & (lo <= hi) & jnp.isfinite(lo) & jnp.isfinite(hi)
    value = eqx.error_if(value, ~jnp.all(valid),
                         "ordinal numerical failure: interval has no representable interior")
    value = jnp.where(value <= lower, lo, jnp.where(value >= upper, hi, value))
    return eqx.error_if(value, ~jnp.all(jnp.isfinite(value)),
                        "ordinal numerical failure: nonfinite interval draw")


def sample_truncated_normal(key, lower, upper):
    """Standard normal on an open interval, stable in both float32 tails.

    Reflect negative intervals, use truncated exponential rejection on the
    positive side, uniform rejection for narrow intervals across zero, and
    ordinary normal rejection for wide intervals across zero. Each element
    retains its first accepted proposal. Inactive branches use benign inputs.
    Errors propagate through jit/vmap and surface when the host synchronizes.
    """
    lower, upper = jnp.broadcast_arrays(jnp.asarray(lower, jnp.float32),
                                        jnp.asarray(upper, jnp.float32))
    # Validate before entering a loop that could otherwise never terminate.
    lower = eqx.error_if(lower, ~jnp.all(
        (lower < upper) & (jnp.nextafter(lower, jnp.float32(jnp.inf)) < upper)),
        "ordinal numerical failure: invalid interval or no representable interior")
    reflect = upper <= 0
    a = jnp.where(reflect, -upper, lower)
    b = jnp.where(reflect, -lower, upper)
    positive = a >= 0
    narrow = (~positive) & ((b - a) <= 2)
    ap = jnp.where(positive, a, 0.0)
    bp = jnp.where(positive, b, 1.0)
    rate = ap / 2 + jnp.hypot(ap, jnp.float32(2)) / 2
    mass = -jnp.expm1(-rate * (bp - ap))
    au = jnp.where(narrow, a, -1.0)
    bu = jnp.where(narrow, b, 1.0)

    def body(carry):
        key, pending, draws = carry
        key, ku, kv, kn = random.split(key, 4)
        # Exclude both endpoints with ordinary representable float32 uniforms.
        u = random.uniform(ku, a.shape, jnp.float32, minval=jnp.finfo(jnp.float32).eps)
        logv = jnp.log(random.uniform(kv, a.shape, jnp.float32,
                                     minval=jnp.finfo(jnp.float32).eps))
        xp = ap - jnp.log1p(-u * mass) / rate
        xu = au + (bu - au) * u
        xn = random.normal(kn, a.shape, jnp.float32)
        proposal = jnp.where(positive, xp, jnp.where(narrow, xu, xn))
        accepted = jnp.where(positive, logv <= -0.5 * (xp - rate)**2,
                            jnp.where(narrow, logv <= -0.5 * xu**2,
                                      (xn > a) & (xn < b)))
        draws = jnp.where(pending & accepted, proposal, draws)
        return key, pending & ~accepted, draws

    _, _, draws = lax.while_loop(lambda c: jnp.any(c[1]), body,
                                (key, jnp.ones(a.shape, bool), jnp.zeros_like(a)))
    return _interior(jnp.where(reflect, -draws, draws), lower, upper)


def full_cutpoints(cutpoints):
    return jnp.concatenate((jnp.array([-jnp.inf, 0], jnp.float32),
                            cutpoints, jnp.array([jnp.inf], jnp.float32)))


def sample_ordinal_latents(key, labels, mean, sd, cutpoints, obs_mask, old_z):
    """Refresh latents with scalar or SUR conditional means and scales.

    Missing placeholders (including old_z) never enter interval arithmetic.
    old_z is accepted for the common augmentation interface; missing draws
    are always zero, irrespective of the previous latent.
    """
    full = full_cutpoints(cutpoints)
    codes = jnp.where(obs_mask, labels, 0).astype(jnp.int32)
    mean = jnp.where(obs_mask, mean, 0.0)
    sd = jnp.where(obs_mask, sd, 1.0)
    sd = eqx.error_if(sd, ~jnp.all(jnp.isfinite(sd) & (sd > 0)),
                      "ordinal numerical failure: invalid conditional scale")
    mean = eqx.error_if(mean, ~jnp.all(jnp.isfinite(mean)),
                        "ordinal numerical failure: nonfinite conditional mean")
    lower = jnp.where(obs_mask, full[codes], -jnp.inf)
    upper = jnp.where(obs_mask, full[codes + 1], jnp.inf)
    draw = sample_truncated_normal(key, (lower - mean) / sd, (upper - mean) / sd)
    latent = _interior(mean + sd * draw, lower, upper)
    return jnp.where(obs_mask, latent, 0.0)


def _log_cutpoint_target(u, labels, mean, sd, obs_mask, prior_scale):
    gaps = jnp.exp(u)
    cutpoints = jnp.cumsum(gaps)
    all_cuts = jnp.concatenate([jnp.array([0.0], jnp.float32), cutpoints])
    codes = jnp.where(obs_mask, labels, 0).astype(jnp.int32)
    m = jnp.where(obs_mask, mean, 0.0)
    s = jnp.where(obs_mask, sd, 1.0)
    y_safe = jnp.clip(codes, 1, len(all_cuts) - 1)

    hi_int = (all_cuts[y_safe] - m) / s
    lo_int = (all_cuts[y_safe - 1] - m) / s

    reflected = lo_int >= 0
    hi = jnp.where(reflected, -lo_int, hi_int)
    lo = jnp.where(reflected, -hi_int, lo_int)
    loghi = jax_log_ndtr(hi)
    loglo = jax_log_ndtr(lo)
    diff = jnp.minimum(loglo - loghi, -1e-15)
    log_p_int = loghi + jnp.log(-jnp.expm1(diff))

    log_p_0 = jax_log_ndtr(-m / s)
    log_p_last = jax_log_ndtr(-(all_cuts[-1] - m) / s)

    log_p = jnp.where(codes == 0, log_p_0, jnp.where(codes == len(all_cuts), log_p_last, log_p_int))
    cell_log_p = jnp.where(obs_mask, log_p, 0.0)

    log_lik = jnp.sum(cell_log_p)
    log_prior = -0.5 * jnp.sum(cutpoints ** 2) / (prior_scale ** 2)
    log_jac = jnp.sum(u)
    return log_lik + log_prior + log_jac


def sample_cutpoints_marginalized(key, cutpoints, labels, mean, sd, obs_mask, prior_scale):
    """Partially collapsed threshold update with marginalized latents.

    Proposes in log-gap space u_j = log(theta_j - theta_{j-1}), evaluating
    the multinomial ordered probit marginal log-likelihood with the prior and
    Jacobian. Bypasses the Albert-Chib sequential latent sandwich bottleneck.
    """
    if cutpoints.shape[0] == 0:
        return cutpoints
    prev = jnp.concatenate([jnp.array([0.0], jnp.float32), cutpoints[:-1]])
    gaps = jnp.maximum(cutpoints - prev, jnp.finfo(jnp.float32).tiny)
    u = jnp.log(gaps)
    lp = _log_cutpoint_target(u, labels, mean, sd, obs_mask, prior_scale)

    n_obs = jnp.sum(obs_mask)
    scales = jnp.array([0.5, 1.2, 2.5], jnp.float32) / jnp.sqrt(jnp.maximum(n_obs.astype(jnp.float32), 1.0))

    def mh_step(i, carry):
        k, u_curr, lp_curr = carry
        k, ks, kp, ku = random.split(k, 4)
        scale = random.choice(ks, scales)
        u_prop = u_curr + scale * random.normal(kp, u_curr.shape)
        lp_prop = _log_cutpoint_target(u_prop, labels, mean, sd, obs_mask, prior_scale)
        valid = jnp.isfinite(lp_prop)
        log_alpha = jnp.where(valid, lp_prop - lp_curr, -jnp.inf)
        accept = jnp.log(random.uniform(ku)) < log_alpha
        u_next = jnp.where(accept, u_prop, u_curr)
        lp_next = jnp.where(accept, lp_prop, lp_curr)
        return (k, u_next, lp_next)

    _, u_final, _ = lax.fori_loop(0, 3, mh_step, (key, u, lp))
    new_cuts = jnp.cumsum(jnp.exp(u_final))
    return jnp.where(jnp.isfinite(new_cuts), new_cuts, cutpoints)


def category_probabilities(mean_draws, cutpoint_draws):
    """Float64 host probabilities: means (D,B), free thresholds (D,K-2).

    Evaluate CDF differences in log space, reflecting positive tail intervals
    to survival differences. Thresholds stay paired with their own draw.
    """
    means = np.asarray(mean_draws, np.float64)
    cuts = np.asarray(cutpoint_draws, np.float64)
    if means.ndim != 2 or cuts.ndim != 2 or means.shape[0] != cuts.shape[0]:
        raise ValueError("means and cutpoints must have shapes (D,B) and (D,K-2)")
    if (not np.isfinite(means).all() or not np.isfinite(cuts).all()
            or not (np.diff(np.c_[np.zeros(len(cuts)), cuts], axis=1) > 0).all()):
        raise ValueError("ordinal probabilities require finite means and ordered positive cutpoints")
    full = np.c_[np.full(len(cuts), -np.inf), np.zeros(len(cuts)),
                 cuts, np.full(len(cuts), np.inf)]
    lower = full[:, None, :-1] - means[:, :, None]
    upper = full[:, None, 1:] - means[:, :, None]
    reflected = lower >= 0
    loghi = log_ndtr(np.where(reflected, -lower, upper))
    loglo = log_ndtr(np.where(reflected, -upper, lower))
    return np.exp(loghi) * -np.expm1(loglo - loghi)
