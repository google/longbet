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

"""Independent input, probability, and distributional ordered-probit oracles."""
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from scipy import integrate, stats

from longbet import LongBetConfig
from longbet._ordinal import (category_probabilities, prepare_ordinal,
    sample_cutpoints_marginalized, sample_ordinal_latents,
    sample_truncated_normal)


@pytest.mark.parametrize("K", [2, 3, 5])
@pytest.mark.parametrize("observed", [[0], [1], [0, 1]])
def test_preparation_empty_categories(K, observed):
    y = np.array([*observed, np.nan])
    prepared = prepare_ordinal(y, K)
    np.testing.assert_array_equal(prepared.labels[:-1], observed)
    assert prepared.labels[-1] == 0
    assert not prepared.obs_mask[-1]
    assert prepared.cutpoints.shape == (K - 2,)
    assert np.isfinite(prepared.cutpoints).all()
    assert (np.diff(np.r_[0, prepared.cutpoints]) > 0).all()


def test_preparation_sign_and_no_recoding():
    prepared = prepare_ordinal([0]*8 + [2, 4], 5)
    assert prepared.offset == pytest.approx(-0.8416212335729143)
    assert stats.norm.cdf(-prepared.offset) == pytest.approx(.8)
    np.testing.assert_array_equal(prepared.labels, [0]*8 + [2, 4])
    # The binary compatibility branch must use precisely the legacy calculation.
    for y in ([0]*10, [1]*10, [0]*8 + [1]*2):
        assert prepare_ordinal(y, 2).offset == float(stats.norm.ppf(
            np.clip(float(np.mean(y)), 1e-4, 1-1e-4)))


@pytest.mark.parametrize("K", [None, True, np.bool_(True), 1, 0, -3, 2., "3", np.nan])
def test_bad_category_count(K):
    with pytest.raises(ValueError, match="num_categories"):
        LongBetConfig(outcome="ordinal", num_categories=K)
    with pytest.raises(ValueError, match="num_categories"):
        prepare_ordinal([0, 1], K)


@pytest.mark.parametrize("scale", [True, np.bool_(True), 0, -1, np.inf, np.nan, "5"])
def test_bad_prior_scale(scale):
    with pytest.raises(ValueError, match="cutpoint_prior_scale"):
        LongBetConfig(cutpoint_prior_scale=scale)


def test_configuration_contract():
    for outcome in ("continuous", "binary"):
        with pytest.raises(ValueError, match="num_categories"):
            LongBetConfig(outcome=outcome, num_categories=3)
    with pytest.raises(ValueError, match="outcome"):
        LongBetConfig(outcome="nominal")
    cfg = LongBetConfig(outcome="ordinal", num_categories=np.int64(3))
    assert type(cfg.num_categories) is int


@pytest.mark.parametrize("y", [[0, 1.1], [-1, 2], [0, 3], [np.inf, 1],
                               [-np.inf, 0], [np.nan], ["0", "1"], [1j, 0]])
def test_bad_labels(y):
    with pytest.raises(ValueError):
        prepare_ordinal(y, 3)


@pytest.mark.parametrize("a,b", [(-np.inf, np.inf), (-1., 1.), (-.1, .2),
    (-2., 3.), (0., np.inf), (8., 9.), (12., np.inf), (-9., -8.),
    (-np.inf, -12.), (1., 1.0001), (-.0001, .0002)])
def test_interval_distribution(a, b):
    n = 30000
    draws = np.asarray(jax.jit(sample_truncated_normal)(
        jax.random.key(147), jnp.full(n, a), jnp.full(n, b)))
    assert np.isfinite(draws).all()
    assert (draws > a).all() and (draws < b).all()
    # Probability integral transform is independently evaluated by SciPy.
    # D=0.018 has a DKW tail bound below 8e-9 at this sample size.
    uniforms = stats.truncnorm.cdf(draws.astype(np.float64), a, b)
    assert stats.kstest(uniforms, "uniform").statistic < .018
    mean, var = stats.truncnorm.stats(a, b, moments="mv")
    assert abs(draws.astype(np.float64).mean() - mean) < 6*np.sqrt(var/n) + 2e-6


def test_jit_vmap_varying_latent_means_scales_and_missing_placeholders():
    n = 12000
    labels = jnp.arange(n) % 4
    mask = (jnp.arange(n) % 11) != 0
    mean = jnp.where(mask, jnp.sin(jnp.arange(n)), jnp.nan)
    sd = jnp.where(mask, .3 + (jnp.arange(n) % 3), jnp.nan)
    keys = jax.random.split(jax.random.key(843), 2)
    cuts = jnp.array([[1., 2.5], [.2, 1.3]])
    def sample(key, cp):
        return sample_ordinal_latents(key, labels, mean, sd, cp, mask, jnp.full(n, jnp.nan))
    draws = np.asarray(jax.jit(jax.vmap(sample))(keys, cuts))
    for c in range(2):
        full = np.r_[-np.inf, 0, np.array(cuts[c]), np.inf]
        sel = np.asarray(mask)
        lo = full[np.asarray(labels)[sel]]
        hi = full[np.asarray(labels)[sel] + 1]
        assert (draws[c, sel] > lo).all() and (draws[c, sel] < hi).all()
        assert (draws[c, ~sel] == 0).all()
        m, s = np.asarray(mean)[sel], np.asarray(sd)[sel]
        u = stats.truncnorm.cdf(draws[c, sel], (lo-m)/s, (hi-m)/s, loc=m, scale=s)
        assert stats.kstest(u, "uniform").statistic < .025


@pytest.mark.parametrize("lower,upper", [(1., 1.), (2., 1.), (np.nan, 2.),
    (1., np.nextafter(np.float32(1), np.float32(2)))])
def test_invalid_intervals_fail(lower, upper):
    with pytest.raises(Exception, match="ordinal numerical failure"):
        jax.jit(sample_truncated_normal)(jax.random.key(1), lower, upper).block_until_ready()


def test_empty_adjacent_categories_recover_ordered_prior():
    @jax.jit
    def run(key):
        def step(cp, k):
            cp = sample_cutpoints_marginalized(
                k, cp, jnp.array([0]), jnp.array([-.5]), 1.0, jnp.array([True]), 2.0
            )
            return cp, cp
        return jax.lax.scan(step, jnp.array([.2, .5, .9]), jax.random.split(key, 15000))[1]
    draws = np.asarray(run(jax.random.key(921)))[1000:]
    assert (np.diff(np.c_[np.zeros(len(draws)), draws], axis=1) > 0).all()
    # With only category zero observed, the posterior is exactly the ordered
    # half-normal prior. Sum-of-squares / scale^2 has expectation K-2.
    assert np.mean(np.sum((draws/2)**2, axis=1)) == pytest.approx(3., abs=.15)


def test_probabilities_pair_draws_and_tail_accuracy():
    means = np.array([[0, 1, -12, 12], [1, -1, -8, 8.]])
    cuts = np.array([[1.], [2.]])
    p = category_probabilities(means, cuts)
    assert p.shape == (2, 4, 3)
    np.testing.assert_allclose(p[0, 0], [.5, stats.norm.cdf(1)-.5, stats.norm.sf(1)])
    np.testing.assert_allclose(p.sum(-1), 1., atol=3e-16)
    assert (p >= 0).all() and np.isfinite(p).all()
    assert p[0, 2, 1] > 0  # Direct CDF subtraction would return zero here.
    np.testing.assert_allclose(p[0, 2, 1], stats.norm.sf(12)-stats.norm.sf(13), rtol=3e-14)
    p1 = category_probabilities(means+.2, cuts)
    assert ((p1-p) @ np.arange(3) >= -1e-15).all()
    np.testing.assert_allclose((p1-p).sum(-1), 0., atol=3e-16)
    binary = category_probabilities(means, np.empty((2, 0)))
    np.testing.assert_allclose(binary[:, :, 1], stats.norm.cdf(means), rtol=3e-14)


def test_sample_cutpoints_marginalized_contract():
    key = jax.random.key(123)
    # K=2 returns empty array
    k2 = sample_cutpoints_marginalized(key, jnp.empty(0, jnp.float32), jnp.array([0, 1]),
                                       jnp.zeros(2), 1.0, jnp.ones(2, bool), 5.0)
    assert k2.shape == (0,)

    # K=4 produces positive strictly ordered thresholds
    cp_init = jnp.array([1.0, 2.0], jnp.float32)
    labels = jnp.array([0, 1, 2, 3, 0, 1], jnp.float32)
    mean = jnp.array([-0.5, 0.2, 1.1, 2.3, -0.1, 0.4], jnp.float32)
    mask = jnp.array([True, True, True, True, True, False])
    cp_new = jax.jit(sample_cutpoints_marginalized)(key, cp_init, labels, mean, 1.0, mask, 5.0)
    assert cp_new.shape == (2,)
    assert (cp_new > 0).all()
    assert cp_new[0] < cp_new[1]
    assert np.isfinite(cp_new).all()

    # Masked values cannot affect the draw
    labels_alt = labels.at[5].set(999.0)
    cp_alt = jax.jit(sample_cutpoints_marginalized)(key, cp_init, labels_alt, mean, 1.0, mask, 5.0)
    np.testing.assert_array_equal(cp_new, cp_alt)


def _marginalized_chain_draws(labels, mean, sd, mask, initial, prior_scale):
    """Independent chains for fixed-surface threshold distribution checks."""
    @jax.jit
    def run(key, start):
        def step(cuts, subkey):
            cuts = sample_cutpoints_marginalized(
                subkey, cuts, labels, mean, sd, mask, prior_scale)
            return cuts, cuts
        return jax.lax.scan(step, start, jax.random.split(key, 14000))[1]
    keys = jax.random.split(jax.random.key(20260912), 4)
    starts = jnp.asarray(initial)[None, :] * jnp.array([.4, .8, 1.5, 2.5])[:, None]
    samples = np.asarray(jax.vmap(run)(keys, starts))[:, 2000:]
    assert np.isfinite(samples).all()
    return samples.reshape(-1, len(initial)).astype(np.float64)


@pytest.mark.parametrize("conditional_scale", [False, True])
def test_marginalized_threshold_posterior_matches_quadrature(conditional_scale):
    # K=3 has one free threshold, so its posterior can be integrated directly
    # in threshold coordinates, independently of the sampler's log-gap target.
    labels = np.array([0, 1, 1, 2, 2, 2, 999])
    mean = np.array([-.3, .1, .6, 1., 1.5, 2., np.nan])
    sd = np.array([.4, .7, .5, .9, .6, .8, np.nan]) if conditional_scale else np.ones(7)
    mask = np.arange(7) < 6
    prior_scale = 1.7

    def density(theta):
        bounds = np.array([-np.inf, 0., theta, np.inf])
        y, m, s = labels[mask], mean[mask], sd[mask]
        probabilities = stats.norm.cdf((bounds[y + 1] - m) / s) - stats.norm.cdf((bounds[y] - m) / s)
        return np.prod(probabilities) * stats.halfnorm.pdf(theta, scale=prior_scale)

    mass = integrate.quad(density, 0, np.inf, epsabs=1e-12)[0]
    draws = _marginalized_chain_draws(
        jnp.asarray(labels), jnp.asarray(mean), jnp.asarray(sd), jnp.asarray(mask),
        [1.], prior_scale)[:, 0]
    assert (draws > 0).all()
    # Conservative fixed tolerances allow serial dependence; IID KS critical
    # values are inappropriate for MCMC samples. Both moments and CDFs are
    # checked, so an omitted Jacobian or wrong prior cannot pass on shape alone.
    expected_mean = integrate.quad(lambda t: t * density(t), 0, np.inf, epsabs=1e-12)[0] / mass
    assert draws.mean() == pytest.approx(expected_mean, abs=.035)
    for theta in [.5, 1., 1.5, 2.]:
        expected_cdf = integrate.quad(density, 0, theta, epsabs=1e-12)[0] / mass
        assert np.mean(draws <= theta) == pytest.approx(expected_cdf, abs=.025)


def test_marginalized_empty_categories_recover_ordered_normal_prior():
    # Only category zero is observed, whose probability does not depend on the
    # free thresholds. K=4 must therefore recover two sorted half-normal draws.
    scale = 2.
    draws = _marginalized_chain_draws(
        jnp.array([0, 999]), jnp.array([-.2, jnp.nan]),
        jnp.array([1., jnp.nan]), jnp.array([True, False]), [1., 2.], scale)
    assert (draws[:, 0] > 0).all()
    assert (draws[:, 1] > draws[:, 0]).all()
    median = stats.halfnorm.ppf(.5, scale=scale)
    assert np.mean(draws[:, 0] <= median) == pytest.approx(.75, abs=.03)
    assert np.mean(draws[:, 1] <= median) == pytest.approx(.25, abs=.03)
    assert np.mean(np.sum((draws / scale)**2, axis=1)) == pytest.approx(2., abs=.12)
