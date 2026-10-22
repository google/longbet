"""Variance hyperparameters are proper by construction."""

import numpy as np
import pytest

from longbet import LongBet, LongBetConfig, LongBetMulti


@pytest.mark.parametrize("name", ["sigma_prior_a", "sigma_prior_b", "gamma_prior_a", "gamma_prior_b"])
@pytest.mark.parametrize("value", [-1, 0, np.nan, np.inf, True, "1"])
def test_invalid_variance_hyperparameters(name, value):
    with pytest.raises(ValueError, match=f"{name} must be a finite positive scalar"):
        LongBetConfig(**{name: value})


def test_binary_only_archives_and_children_keep_fixed_unit_variances(tmp_path):
    rng = np.random.default_rng(13)
    y = {name: rng.integers(0, 2, size=(6, 4)).astype(float) for name in ("a", "b")}
    x = rng.normal(size=(6, 2))
    z = np.zeros((6, 4)); z[:3, 2:] = 1
    model = LongBetMulti(num_burnin=2, num_sweeps=2, num_chains=1,
                        num_trees_pr=2, num_trees_trt=2,
                        max_depth_pr=2, max_depth_trt=2).fit(y, x, z, outcome="binary")
    path = tmp_path / "binary.npz"
    model.save(path)
    loaded = LongBetMulti.load(path)
    assert loaded.config == model.config
    for name in y:
        child = tmp_path / f"{name}.npz"
        loaded[name].save(child)
        scalar = LongBet.load(child)
        np.testing.assert_array_equal(scalar.trace.sigma2, 1)
        assert scalar.multi_origin["outcomes"] == ["binary", "binary"]
