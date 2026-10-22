# How to contribute

We'd love to accept your patches and contributions to this project.

## Before you begin

### Sign our Contributor License Agreement

Contributions to this project must be accompanied by a
[Contributor License Agreement](https://cla.developers.google.com/about) (CLA).
You (or your employer) retain the copyright to your contribution; this simply
gives us permission to use and redistribute your contributions as part of the
project.

If you or your current employer have already signed the Google CLA (even if it
was for a different project), you probably don't need to do it again.

Visit <https://cla.developers.google.com/> to see your current agreements or to
sign a new one.

### Review our community guidelines

This project follows
[Google's Open Source Community Guidelines](https://opensource.google/conduct/).

## Contribution process

### Code reviews

All submissions, including submissions by project members, require review. We
use GitHub pull requests for this purpose. Consult
[GitHub Help](https://help.github.com/articles/about-pull-requests/) for more
information on using pull requests.

## Working on LongBet

One engine, two front doors. The sampler lives in `src/longbet/` (Python + JAX);
the R package in `R/` marshals data in and results out. Nothing statistical is
implemented twice — if you find yourself writing an estimand in R, delegate to
the engine instead.

## Layout

```
src/longbet/
  _config.py       LongBetConfig: every user-facing option, with validation
  _design.py       the unified design matrix and the spec predict replays
  _state.py        LongBetState, initialization, the chain-axis partition
  _step.py         one Gibbs sweep
  _gp.py           kernels, the whitened conditional draw, the forecast
  _ridge.py        the Metropolis move along the beta-nu scale ridge
  _scales.py       per-sweep prec_scale rebuild, mirroring bartz
  _loop.py         the MCMC driver and trace buffers
  _model.py        the LongBet estimator: fit, predict, diagnostics
  _rollout.py      adoption cohorts: rollout_summary / plot_rollout
  _summary.py      blocked posterior summaries
  _diagnostics.py  ESS and R-hat on the ATT
  _io.py           .npz serialization
R/                 the R front door
contract/          the shared API contract, vendored into inst/contract/
```

## Running the tests

```bash
pytest -m "not slow"    # fast: ~1 minute
pytest -m slow          # Geweke, coverage, DGP recovery: ~5 minutes
pytest                  # everything

# The R side needs the engine installed; Docker is the reproducible way.
docker build -f Dockerfile.rtest -t longbet-rtest . && docker run --rm longbet-rtest
docker build -f Dockerfile.rtest-nopython -t longbet-rtest-nopy . \
  && docker run --rm -e RETICULATE_PYTHON=/nonexistent/python longbet-rtest-nopy
```

Set `LONGBET_BARTZ_BCF_PATH` to a checkout of
[bartz#189](https://github.com/bartz-org/bartz/pull/189)'s `src/bartz` to enable
the BCF cross-check; it skips otherwise.

## Rules that exist for a reason

**Adding an option** means editing `LongBetConfig`, `contract/longbet-api.yaml`,
`R/longbet.R`, and copying the contract to `inst/contract/`. `tests/test_contract.py`
fails if you skip the contract, and `tests/testthat/test-contract.R` fails if you
skip R. That is the mechanism that keeps the two front doors from drifting.

**The full-model residual is carried in data units.** `bartz` stores its residual
scaled by `resid_unit`, which differs between the two forests. Convert only at
the boundary of a `bartz` call, with `_load` / `_read`.
`tests/test_units.py::test_full_model_residual_invariant` checks that
`resid == y - fitted` after every sweep; it will catch a missed update.

**Recover forest fits from cached leaf memberships**, without traversing trees
again. Residual differencing after division by tiny coding/GP weights loses
small leaf updates to float32 cancellation. Use physical fit deltas when
updating the full-model residual, including cells with negligible weight.

**Never form the GP precision matrix.** The squared-exponential Gram matrix is
severely ill-conditioned and its inverse does not survive float32. Work through
the Cholesky factor; see the module docstring in `_gp.py` for the measurements.

**Chain membership is read from `bartz`'s field metadata**, not hardcoded.
`Forest.leaf_tree` carries a chain axis and `Forest.leaf_unit` does not; getting
that wrong silently breaks `evaluate_trace`. See `_state.chain_filter_spec`.

**Sweeps must be jitted at the boundary.** `bartz.mcmcstep.step` donates its
argument's buffers, so an un-jitted sweep deletes arrays the caller still holds.

**`predict` rebuilds the design from the persisted spec.** Adding an input means
adding a `Block` with a `source`, so a design fitted with a column refuses to
build without it rather than silently shifting every column index.

**When touching a conditional, run the Geweke test.** It is the only test that
catches a wrong constant, and `test_geweke_detects_a_broken_conditional` exists
to prove it still can.

**Diagnose every reported effect**, including subgroup and joint-event summaries.
The treatment factors have sign/scale symmetries, but their proper priors still
define a parameter posterior. Do not hide chain disagreement by clipping initial
GP draws or deterministically folding the sampled trajectory to one orientation.
Parameter traces help investigate poor mixing; agreement of ATT means alone
does not establish posterior exploration.

**Coupled continuous outcomes require explicit proper innovation priors.** The
reference `(0,0)` variance prior can produce an improper joint posterior and is
rejected at fit/load. Positive inverse-gamma hyperparameters address that defect;
convergence and calibration remain separate requirements.

**Chains start overdispersed**, from draws of the priors. R-hat from identical
starts measures how far two random streams drifted, which is a different and
easier question. If you change initialisation, keep the residual invariant: the
fit at initialisation is exactly `gamma_i`, so dispersing `gamma` means
correcting `resid`.

**Plotting stays optional.** `rollout_summary()` returns a table and needs only
pandas; `plot_rollout()` imports matplotlib when called. Do not move a plotting
import to module scope.
