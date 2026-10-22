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

"""The unified design matrix and its persisted spec.

These tests pin the two things that make a column silently disappear: the
private ``Binner._splits`` accessor, and the requirement that predict rebuild
exactly the columns fit used.
"""

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from bartz.prepcovars import UniqueQuantileBinner

from longbet import LongBet, LongBetConfig
from longbet._design import Block, Design, integer_grid_block, quantile_block


def test_quantile_block_matches_bartz_binner_on_unseen_data():
    """The persisted cutpoints must reproduce bartz's own binning exactly."""
    rng = np.random.default_rng(0)
    fit_raw = rng.normal(size=(3, 300)).astype(np.float32)
    new_raw = rng.normal(size=(3, 40)).astype(np.float32)

    block = quantile_block("x", "x", fit_raw, 16, mu_visible=True, nu_visible=True)
    reference = UniqueQuantileBinner(jnp.asarray(fit_raw), max_bins=16)

    assert np.array_equal(block.bin(new_raw), np.asarray(reference.bin(jnp.asarray(new_raw))))
    assert np.array_equal(block.max_split, np.asarray(reference.max_split))


def test_block_survives_serialization():
    rng = np.random.default_rng(1)
    fit_raw = rng.normal(size=(2, 200)).astype(np.float32)
    new_raw = rng.normal(size=(2, 25)).astype(np.float32)
    block = quantile_block("x", "x", fit_raw, 8, mu_visible=True, nu_visible=False)
    restored = Block.from_dict(block.to_dict())
    assert np.array_equal(restored.bin(new_raw), block.bin(new_raw))
    assert restored.mu_visible and not restored.nu_visible


def test_integer_grid_block_puts_a_cut_between_every_level():
    block = integer_grid_block("t", "t", 6, mu_visible=True, nu_visible=True)
    binned = block.bin(np.arange(6, dtype=np.float32).reshape(1, -1))
    assert np.array_equal(binned.ravel(), np.arange(6))


def test_empty_time_grid_is_rejected_not_silently_binned_to_zero():
    with pytest.raises(ValueError, match="at least one level"):
        integer_grid_block("t", "t", 0, mu_visible=True, nu_visible=True)
    # A genuinely single-period panel remains valid.
    block = integer_grid_block("t", "t", 1, mu_visible=True, nu_visible=True)
    assert block.splits.shape == (1, 0)


def test_saved_panel_grid_dimensions_are_validated():
    t = integer_grid_block("t", "t", 6, mu_visible=True, nu_visible=True)
    s = integer_grid_block("s", "s", 5, mu_visible=False, nu_visible=True)
    design = Design((t, s))
    design.validate_panel_grids(6, 4)
    with pytest.raises(ValueError, match="Refit from the original data"):
        design.validate_panel_grids(7, 4)
    with pytest.raises(ValueError, match="Refit from the original data"):
        design.validate_panel_grids(6, 5)


def test_max_split_blocks_the_right_columns():
    rng = np.random.default_rng(2)
    x = quantile_block("x", "x", rng.normal(size=(2, 50)).astype(np.float32), 8,
                       mu_visible=True, nu_visible=False)
    s = integer_grid_block("s", "s", 4, mu_visible=False, nu_visible=True)
    design = Design((x, s))
    assert np.all(design.max_split_mu[2:] == 0)   # mu cannot split on S
    assert np.all(design.max_split_nu[:2] == 0)   # nu cannot split on x


def test_missing_source_is_refused_not_silently_dropped():
    """A design fitted with a column must not build without it."""
    rng = np.random.default_rng(3)
    x = quantile_block("x", "x", rng.normal(size=(1, 20)).astype(np.float32), 4,
                       mu_visible=True, nu_visible=True)
    ps = quantile_block("ps", "ps", rng.uniform(size=(1, 20)).astype(np.float32), 4,
                        mu_visible=True, nu_visible=False)
    with pytest.raises(ValueError, match="needs input 'ps'"):
        Design((x, ps)).build({"x": rng.normal(size=(1, 5)).astype(np.float32)})


def _panel(seed=0, N=40, T=5):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, 3)).astype(np.float32)
    z = np.zeros((N, T), dtype=np.float32)
    z[: N // 2, 2:] = 1.0
    t = np.arange(1, T + 1, dtype=np.float32)
    y = (x[:, [0]] + 1.5 * z + rng.normal(0, 0.2, (N, T))).astype(np.float32)
    return x, y, z, t


def test_propensity_score_survives_the_round_trip_to_predict():
    """Regression: the ps column used to be dropped when rebuilding the design."""
    x, y, z, t = _panel()
    ps = (1.0 / (1.0 + np.exp(-2.0 * x[:, 0]))).astype(np.float32)
    model = LongBet(LongBetConfig(num_sweeps=6, num_burnin=3, num_trees_pr=5,
                                  num_trees_trt=5, num_chains=1, random_seed=1))
    model.fit(y=y, x=x, z=z, t=t, ps=ps)

    assert "ps" in [b.name for b in model.design_.blocks]
    # The fitted design has a ps column, so predicting without one must fail
    # loudly rather than quietly evaluating trees against shifted columns.
    with pytest.raises(ValueError, match="needs input 'ps'"):
        model.predict(x=x, z=z, t=t)

    pred = model.predict(x=x, z=z, t=t, ps=ps)
    assert pred.tauhats.shape == (x.shape[0], z.shape[1], 6)


def test_x_trt_gives_the_treatment_forest_its_own_columns():
    """Regression: x_trt used to be accepted and then ignored entirely."""
    x, y, z, t = _panel()
    rng = np.random.default_rng(7)
    x_trt = rng.normal(size=(x.shape[0], 2)).astype(np.float32)

    model = LongBet(LongBetConfig(num_sweeps=6, num_burnin=3, num_trees_pr=5,
                                  num_trees_trt=5, num_chains=1, random_seed=1))
    model.fit(y=y, x=x, z=z, t=t, x_trt=x_trt)

    names = [b.name for b in model.design_.blocks]
    assert "x_trt" in names
    by_name = {b.name: b for b in model.design_.blocks}
    assert by_name["x"].mu_visible and not by_name["x"].nu_visible
    assert by_name["x_trt"].nu_visible and not by_name["x_trt"].mu_visible

    # Neither forest may split on a column its max_split blocks.
    for trace, bound in ((model.trace.mu_trace, model.design_.max_split_mu),
                         (model.trace.nu_trace, model.design_.max_split_nu)):
        used = np.unique(np.asarray(trace.var_tree))
        used = used[used > 0]
        assert all(bound[v] > 0 for v in used if v < len(bound))


def test_split_time_flags_apply_to_their_own_forest():
    """split_time_ps governs mu; split_calendar_trt governs nu; s is never visible."""
    x, y, z, t = _panel()
    model = LongBet(LongBetConfig(num_sweeps=4, num_burnin=2, num_trees_pr=4,
                                  num_trees_trt=4, split_time_ps=False,
                                  split_calendar_trt=True,
                                  num_chains=1, random_seed=1))
    model.fit(y=y, x=x, z=z, t=t)
    by_name = {b.name: b for b in model.design_.blocks}
    assert not by_name["t"].mu_visible, "split_time_ps=False must blind mu to time"
    assert by_name["t"].nu_visible, "split_calendar_trt=True enables calendar time for nu"
    assert not by_name["s"].nu_visible, "nu never sees the exposure index"
    assert not by_name["s"].mu_visible, "mu is never shown the exposure index"


def test_split_calendar_trt_blinds_only_the_treatment_forest():
    """split_calendar_trt=False removes calendar time from nu and leaves mu alone.

    The block is switched off here on purpose. With it on, the prognostic
    forest is blinded to calendar time whatever ``split_calendar_trt`` says, so
    the assertion below would be testing the block rather than the flag.
    ``test_default_moves_calendar_time_out_of_the_forest`` covers that case.
    """
    x, y, z, t = _panel()
    model = LongBet(LongBetConfig(num_sweeps=4, num_burnin=2, num_trees_pr=4,
                                  num_trees_trt=4, use_trend_block=False,
                                  split_calendar_trt=False, num_chains=1, random_seed=1))
    model.fit(y=y, x=x, z=z, t=t)
    by_name = {b.name: b for b in model.design_.blocks}
    assert by_name["t"].mu_visible, "the prognostic forest still sees calendar time"
    assert not by_name["t"].nu_visible, "split_calendar_trt=False must blind nu to t"
    assert not by_name["s"].nu_visible
    max_split = model.design_.max_split_nu
    assert (max_split[-2:] == 0).all(), "both time columns are blocked for nu"
    assert (max_split[:-2] > 0).all(), "covariate columns stay available to nu"
    pred = model.predict(x, z, t=t, summary_only=True)
    assert np.all(np.isfinite(pred.tau_summary.mean))
    with pytest.raises(ValueError, match="split_calendar_trt"):
        LongBetConfig(split_calendar_trt=1)


def test_default_moves_calendar_time_out_of_the_forest():
    """Under the shipped defaults the block carries calendar time, not the forests.

    Shipping these two apart is the bug this test exists to catch. The forest
    blinded with no block is the configuration that scored coverage 0.00; the
    forest unblinded with the block on is the one that produced NaN draws.
    """
    x, y, z, t = _panel()
    cfg = LongBetConfig(num_sweeps=4, num_burnin=2, num_trees_pr=4,
                        num_trees_trt=4, num_chains=1, random_seed=1)
    assert cfg.use_trend_block and cfg.trend_period_effects
    assert cfg.split_calendar_mu is False
    assert cfg.split_calendar_trt is False
    model = LongBet(cfg)
    model.fit(y=y, x=x, z=z, t=t)
    by_name = {b.name: b for b in model.design_.blocks}
    assert not by_name["t"].mu_visible, "the block is on, so the forest must be blind to t"
    assert not by_name["t"].nu_visible, "default split_calendar_trt=False blinds nu to t"
    # The block has to actually exist, or the forest was blinded for nothing.
    assert model.state.use_trend_block
    assert model.state.trend_design.shape[1] > 0
    assert model.state.trend_num_period > 0, "the common-time factor must be present"



def test_time_varying_covariates_are_routed_to_their_forest():
    x, y, z, t = _panel()
    rng = np.random.default_rng(9)
    N, T = z.shape
    x_tv = rng.normal(size=(N, T, 2)).astype(np.float32)
    x_trt_tv = rng.normal(size=(N, T, 1)).astype(np.float32)
    model = LongBet(LongBetConfig(num_sweeps=4, num_burnin=2, num_trees_pr=4,
                                  num_trees_trt=4, num_chains=1, random_seed=1))
    model.fit(y=y, x=x, z=z, t=t, x_tv=x_tv, x_trt_tv=x_trt_tv)
    by_name = {b.name: b for b in model.design_.blocks}
    assert by_name["x_tv"].mu_visible and not by_name["x_tv"].nu_visible
    assert by_name["x_trt_tv"].nu_visible and not by_name["x_trt_tv"].mu_visible
    pred = model.predict(x=x, z=z, t=t, x_tv=x_tv, x_trt_tv=x_trt_tv)
    assert np.all(np.isfinite(pred.tauhats))


def test_large_panel_binning_needs_and_gets_a_key():
    """Regression: above 100,000 cells the quantile binner demands a PRNG key.

    That is the normal case for this package, not an edge case, so the failure
    was a hard error on exactly the panels it targets.
    """
    rng = np.random.default_rng(4)
    raw = rng.normal(size=(2, 120_000)).astype(np.float32)
    block = quantile_block("x", "x", raw, 32, mu_visible=True, nu_visible=True,
                           key=jax.random.key(3))
    assert block.n_cols == 2
    assert block.bin(raw[:, :10]).shape == (2, 10)


def test_binning_is_reproducible_from_the_seed():
    """Subsampled quantiles must still be deterministic given the model seed."""
    rng = np.random.default_rng(5)
    raw = rng.normal(size=(1, 150_000)).astype(np.float32)
    a = quantile_block("x", "x", raw, 16, mu_visible=True, nu_visible=True,
                       key=jax.random.key(1))
    b = quantile_block("x", "x", raw, 16, mu_visible=True, nu_visible=True,
                       key=jax.random.key(1))
    c = quantile_block("x", "x", raw, 16, mu_visible=True, nu_visible=True,
                       key=jax.random.key(2))
    assert np.array_equal(a.splits, b.splits)
    assert not np.array_equal(a.splits, c.splits)


def test_fit_is_reproducible_from_random_seed():
    x, y, z, t = _panel(seed=6, N=50, T=6)
    cfg = dict(num_sweeps=8, num_burnin=4, num_trees_pr=5, num_trees_trt=5,
               num_chains=1)
    a = LongBet(LongBetConfig(random_seed=11, **cfg)).fit(y=y, x=x, z=z, t=t)
    b = LongBet(LongBetConfig(random_seed=11, **cfg)).fit(y=y, x=x, z=z, t=t)
    c = LongBet(LongBetConfig(random_seed=12, **cfg)).fit(y=y, x=x, z=z, t=t)
    np.testing.assert_array_equal(np.asarray(a.trace.beta), np.asarray(b.trace.beta))
    assert not np.array_equal(np.asarray(a.trace.beta), np.asarray(c.trace.beta))
