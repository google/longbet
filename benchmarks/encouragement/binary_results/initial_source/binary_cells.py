"""Exact single-horizon binary observed-data controls and identification bounds.

The factorization p(D|Z,X) p(Y|D,Z,X) is observational. It induces posterior
covariance between both reduced forms without treating adoption as randomized.
Sensitivity bounds concern stratum-conditional mean risk differences, not sharp
individual effects. No prior density is imposed inside the identified sets.
"""
from __future__ import annotations

import numpy as np


def validate(y, d, z, x):
    y, d, z = [np.asarray(v, dtype=float) for v in (y, d, z)]
    x = np.asarray(x, dtype=float)
    if y.ndim != 1 or d.shape != y.shape or z.shape != y.shape or len(y) < 4:
        raise ValueError("Require aligned single-horizon outcome, adoption and assignment vectors.")
    if not all(np.isin(v, [0, 1]).all() for v in (y, d, z)):
        raise ValueError("y, d and z must contain finite binary values.")
    if not 0 < z.sum() < len(z):
        raise ValueError("Both randomized arms are required.")
    if x.ndim != 2 or len(x) != len(y) or not np.isfinite(x).all():
        raise ValueError("x must be a finite baseline covariate matrix.")
    grid, inverse, counts = np.unique(x, axis=0, return_inverse=True, return_counts=True)
    return y, d, z, x, grid, inverse, counts / len(y)


def _counts(y, d, z, inverse, groups):
    uptake = np.zeros((groups, 2, 2))  # group,Z,D
    outcomes = np.zeros((groups, 2, 2, 2))  # group,D,Z,Y
    for g, yy, dd, zz in zip(inverse, y.astype(int), d.astype(int), z.astype(int)):
        uptake[g, zz, dd] += 1
        outcomes[g, dd, zz, yy] += 1
    return uptake, outcomes


def ordered_beta(rng, a, b, size, *, max_candidates=2000000):
    """Independent exact posterior pairs conditional on p1>=p0, with a hard cap."""
    kept, total, attempted = [], 0, 0
    while total < size and attempted < max_candidates:
        batch = min(max(2 * (size - total), 256), 65536, max_candidates - attempted)
        proposed = rng.beta(a, b, size=(batch, 2))
        accepted = proposed[proposed[:, 1] >= proposed[:, 0]]
        kept.append(accepted)
        total += len(accepted)
        attempted += batch
    if total < size:
        raise RuntimeError(f"Ordered Beta rejection budget exhausted: {total}/{size} accepted in {attempted} candidates.")
    return np.concatenate(kept, axis=0)[:size], attempted


def fit_cells(y, d, z, x, *, seed, chains=4, draws=1000, first_stage="unrestricted",
              prior=1., max_candidates=2000000, independent_reduced_forms=False):
    """IID posterior draws; p(C,K,G,Z), q(C,K,G,D,Z), empirical-X weights.

    `monotone` truncates the independent Beta uptake prior to p1>=p0 separately
    at each prespecified baseline cell. `null` pools uptake across assignment.
    `independent_reduced_forms` supplies an explicitly different comparison model
    for Y|Z and D|Z; it does not define a joint observed Y,D likelihood.
    """
    y, d, z, x, grid, inverse, weights = validate(y, d, z, x)
    if first_stage not in ("unrestricted", "monotone", "null"):
        raise ValueError("Unknown first_stage.")
    if not np.isfinite(prior) or prior <= 0:
        raise ValueError("prior must be positive.")
    for value in (chains, draws, max_candidates):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError("Counts must be positive integers.")
    uptake, outcomes = _counts(y, d, z, inverse, len(grid))
    rng = np.random.default_rng(seed)
    p = np.empty((chains, draws, len(grid), 2))
    proposals = []
    for g, count in enumerate(uptake):
        if first_stage == "monotone":
            value, attempted = ordered_beta(rng, count[:, 1] + prior, count[:, 0] + prior,
                                            chains * draws, max_candidates=max_candidates)
            p[..., g, :] = value.reshape(chains, draws, 2)
            proposals.append(attempted)
        elif first_stage == "null":
            value = rng.beta(count[:, 1].sum() + prior, count[:, 0].sum() + prior, size=(chains, draws))
            p[..., g, :] = value[..., None]
        else:
            p[..., g, :] = rng.beta(count[:, 1] + prior, count[:, 0] + prior, size=(chains, draws, 2))
    q = rng.beta(outcomes[..., 1] + prior, outcomes[..., 0] + prior,
                 size=(chains, draws, len(grid), 2, 2))
    r = p * q[..., 1, :] + (1 - p) * q[..., 0, :]
    if independent_reduced_forms:
        count_y = outcomes.sum(axis=1)  # group,Z,Y
        r = rng.beta(count_y[..., 1] + prior, count_y[..., 0] + prior, size=r.shape)
    return dict(p=p, q=q, r=r, weights=weights, grid=grid, uptake_counts=uptake,
                outcome_counts=outcomes, first_stage=np.asarray(first_stage),
                prior=np.asarray(prior), proposals=np.asarray(proposals),
                independent_reduced_forms=np.asarray(independent_reduced_forms), seed=np.asarray(seed))


def beta_joint_moments(uptake, outcomes, prior=1.):
    """Analytic unrestricted-cell moments of (r_z,p_z), independently of draws."""
    a, b = uptake[..., 1] + prior, uptake[..., 0] + prior
    mp = a / (a + b)
    vp = a * b / ((a + b)**2 * (a + b + 1))
    aq, bq = outcomes[..., 1] + prior, outcomes[..., 0] + prior
    mq = aq / (aq + bq)
    vq = aq * bq / ((aq + bq)**2 * (aq + bq + 1))
    mr = mp * mq[:, 1] + (1 - mp) * mq[:, 0]
    vr = ((vp + mp**2) * vq[:, 1] + (vp + (1 - mp)**2) * vq[:, 0]
          + vp * (mq[:, 1] - mq[:, 0])**2)
    return dict(p_mean=mp, r_mean=mr, p_variance=vp, r_variance=vr,
                covariance=vp * (mq[:, 1] - mq[:, 0]))


def effects(archive):
    p, r, w = archive["p"], archive["r"], archive["weights"]
    dy = np.sum(w * (r[..., 1] - r[..., 0]), axis=-1)
    dd = np.sum(w * (p[..., 1] - p[..., 0]), axis=-1)
    return dict(itt_y=dy, itt_d=dd, wald=np.divide(dy, dd, out=np.full_like(dy, np.nan), where=dd != 0))


def identification_bounds(p, q, weights, *, delta=1.):
    """Sharp mean-probability bounds under no defiers and direct-risk limits.

    delta is scalar or (always, never, complier), each in [0,1]. Return the
    joint encouragement/adoption complier contrast and treatment effects at
    fixed Z=0 and Z=1 separately. Bounds are weighted by complier masses.
    Empty sets or zero total complier mass yield NaN endpoints and explicit flags.
    This transforms every observed-data draw; it never drops incompatible draws.
    """
    p, q, w = np.asarray(p, dtype=float), np.asarray(q, dtype=float), np.asarray(weights, dtype=float)
    if p.shape[-1:] != (2,) or q.shape != (*p.shape[:-1], 2, 2):
        raise ValueError("Require p(...,G,Z) and q(...,G,D,Z).")
    if (not np.isfinite(p).all() or not np.isfinite(q).all()
            or np.any((p < 0) | (p > 1)) or np.any((q < 0) | (q > 1))):
        raise ValueError("Probabilities must be finite and in [0,1].")
    if w.shape != (p.shape[-2],) or np.any(w < 0) or not np.isfinite(w).all() or not np.isclose(w.sum(), 1):
        raise ValueError("weights must be nonnegative and sum to one.")
    if np.any(p[..., 1] < p[..., 0]):
        raise ValueError("Principal-stratum bounds require the no-defiers probability restriction.")
    delta = np.broadcast_to(np.asarray(delta, dtype=float), (3,))
    if not np.isfinite(delta).all() or np.any((delta < 0) | (delta > 1)):
        raise ValueError("Direct-risk bounds must lie in [0,1].")
    da, dn, dc = delta
    pa, pc, pn = p[..., 0], p[..., 1] - p[..., 0], 1 - p[..., 1]
    t11, t00 = p[..., 1] * q[..., 1, 1], (1 - p[..., 0]) * q[..., 0, 0]
    la, ua = np.maximum(0, q[..., 1, 0] - da), np.minimum(1, q[..., 1, 0] + da)
    ln, un = np.maximum(0, q[..., 0, 1] - dn), np.minimum(1, q[..., 0, 1] + dn)
    l11, u11 = np.maximum(0, t11 - pa * ua), np.minimum(pc, t11 - pa * la)
    l00, u00 = np.maximum(0, t00 - pn * un), np.minimum(pc, t00 - pn * ln)
    compatible = np.all(((l11 <= u11 + 1e-14) & (l00 <= u00 + 1e-14)) | (w == 0), axis=-1)
    share = np.sum(w * pc, axis=-1)
    defined = share > 0
    available = compatible & defined
    l10, u10 = np.maximum(0, l11 - pc * dc), np.minimum(pc, u11 + pc * dc)
    l01, u01 = np.maximum(0, l00 - pc * dc), np.minimum(pc, u00 + pc * dc)
    result = dict(complier_share=share, compatible=compatible, defined=defined, available=available)
    for name, low, high in (("complier_encouragement", l11 - u00, u11 - l00),
                             ("treatment_z0", l10 - u00, u10 - l00),
                             ("treatment_z1", l11 - u01, u11 - l01)):
        result[name] = np.stack([np.divide(np.sum(w * v, axis=-1), share,
            out=np.full_like(share, np.nan), where=available) for v in (low, high)], axis=-1)
    return result
