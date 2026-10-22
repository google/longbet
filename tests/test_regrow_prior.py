"""REGROW target probabilities must use the selected variable's rule prior."""
import jax.numpy as jnp
import numpy as np
from scipy.special import logsumexp
from longbet._regrow_move import _node_candidates


def candidates():
    X = jnp.array([[0, 1, 2, 3], [0, 0, 1, 1]], dtype=jnp.uint8)
    return _node_candidates(
        jnp.int32(1), jnp.ones(4), X, jnp.ones(4), jnp.zeros(4),
        jnp.zeros(2, jnp.uint8), jnp.zeros(2, jnp.uint8), jnp.array([3, 1], jnp.uint8),
        None, jnp.array([0., .25]), jnp.float32(0), jnp.float32(1), 1, 256, 2
    )


def test_every_variable_cut_gets_its_own_split_rule_prior():
    stop, splits, leaf, prior, _ = candidates()
    expected = [.75, .25 / 2 / 3, .25 / 2 / 3, .25 / 2 / 3, .25 / 2]
    indices = [1, 2, 3, 257]
    actual = np.r_[float(leaf), [float(prior.ravel()[i]) for i in indices]]
    np.testing.assert_allclose(actual, np.log(expected), atol=3e-7, rtol=0)
    assert prior.shape == splits.shape == (2, 256)


def test_enumerated_regrow_transition_preserves_the_intended_prior():
    stop, splits, leaf, prior, _ = candidates()
    indices = [1, 2, 3, 257]
    target = np.r_[float(leaf), [float(prior.ravel()[i]) for i in indices]]
    logits = np.r_[float(stop), [float(splits.ravel()[i]) for i in indices]]
    logq = logits - logsumexp(logits)
    q = np.exp(logq)
    intended = np.array([.75, .25 / 6, .25 / 6, .25 / 6, .25 / 2])
    P = np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            acc = min(1., np.exp(target[j] - logq[j] - target[i] + logq[i]))
            P[i, j] += q[j] * acc
            P[i, i] += q[j] * (1 - acc)
    np.testing.assert_allclose(intended @ P, intended, atol=3e-8, rtol=0)
    np.testing.assert_allclose(intended[:, None] * P, intended[None, :] * P.T, atol=3e-8, rtol=0)


def test_weighted_leaf_integral_matches_observation_space_gaussian():
    X = np.array([[0, 1, 2, 3], [0, 0, 1, 1]], dtype=jnp.uint8)
    w = np.array([1., 2., .5, 0.])
    y = np.array([.1, .2, -.5, .8])
    e, lam = 3., 2.
    stop, splits, leaf, prior, _ = _node_candidates(
        jnp.int32(1), jnp.ones(4), jnp.asarray(X),
        jnp.asarray(w, dtype=jnp.float32), jnp.asarray(y, dtype=jnp.float32),
        jnp.zeros(2, jnp.uint8), jnp.zeros(2, jnp.uint8), jnp.array([3, 1], jnp.uint8),
        None, jnp.array([0., .25]), jnp.float32(e), jnp.float32(lam), 1, 256, 2
    )
    choices = [(None, None, .75), (0, 1, .25 / 6), (0, 2, .25 / 6), (0, 3, .25 / 6), (1, 1, .25 / 2)]
    actual = [float(stop)] + [float(splits[v, c]) for v, c, _ in choices[1:]]
    reference = []
    mask = w > 0
    for v, c, p in choices:
        A = np.ones((4, 1)) if v is None else np.column_stack((X[v] < c, X[v] >= c)).astype(float)
        B = A[mask]
        Sigma = np.diag(1 / (e * w[mask])) + B @ B.T / lam
        reference.append(
            np.log(p) - .5 * np.linalg.slogdet(Sigma)[1]
            - .5 * y[mask] @ np.linalg.solve(Sigma, y[mask])
        )
    difference = np.asarray(actual) - reference
    np.testing.assert_allclose(difference, np.full(5, difference[0]), atol=3e-6, rtol=0)
