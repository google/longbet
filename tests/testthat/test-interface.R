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

test_that("fit returns the documented object and prints", {
  skip_without_engine()
  d <- make_panel()
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 20, num_burnin = 10, num_chains = 1,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 42)

  expect_s3_class(fit, "longbet")
  expect_equal(nrow(fit$gamma_draws), nrow(d$x))
  expect_equal(ncol(fit$gamma_draws), 20L)
  expect_equal(ncol(fit$beta_values), 20L)
  expect_equal(dim(fit$sigma0_draws), c(1L, 20L))
  expect_gt(length(fit$raw_model), 0)
  expect_true(is.environment(fit$handle))
  expect_output(print(fit), "LongBet fit")

  # sigma0_draws is in standardized units: multiplying by sdy recovers the
  # residual standard deviation on the data scale, which is about 0.25 here.
  expect_lt(abs(mean(fit$sigma0_draws) * fit$sdy - 0.25), 0.15)
})

test_that("predict returns full draws and exact summaries", {
  skip_without_engine()
  d <- make_panel()
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 20, num_burnin = 10, num_chains = 1,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 42)

  full <- predict(fit, x = d$x, z = d$z, t = d$t, summary_only = FALSE)
  expect_s3_class(full, "longbet.pred")
  expect_equal(dim(full$tauhats), c(nrow(d$y), ncol(d$y), 20L))
  expect_equal(dim(full$preds), c(nrow(d$y), ncol(d$y), 20L))

  lean <- predict(fit, x = d$x, z = d$z, t = d$t, summary_only = TRUE)
  expect_null(lean$tauhats)
  expect_equal(lean$tauhats.mean, full$tauhats.mean, tolerance = 1e-5)
  expect_equal(lean$tauhats.lower, full$tauhats.lower, tolerance = 1e-5)
  # The ATT survives summary_only, so get_att() works in both modes.
  expect_equal(get_att(lean)$att, get_att(full)$att, tolerance = 1e-5)
})

test_that("estimand getters agree with the engine and recover the effect", {
  skip_without_engine()
  d <- make_panel(n = 120, effect = 1.2)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 40, num_burnin = 20, num_chains = 1,
                 num_trees_pr = 10, num_trees_trt = 10, random_seed = 1)
  pred <- predict(fit, x = d$x, z = d$z, t = d$t)

  att <- get_att(pred)
  expect_length(att$att, max(d$s))
  expect_equal(att$exposure, seq_len(max(d$s)))
  expect_equal(dim(att$intervals), c(2L, max(d$s)))
  # True ATT at exposure s is 1.2 * sqrt(s).
  truth <- d$effect * sqrt(seq_len(max(d$s)))
  expect_gt(cor(att$att, truth), 0.9)

  catt <- get_catt(pred)
  expect_equal(dim(catt$catt), dim(d$y))
  expect_true(all(catt$lower <= catt$upper))
  expect_equal(getTaus(pred), catt$catt)
  expect_equal(dim(getMus(pred)), dim(d$y))
})

test_that("R and the Python engine compute the same ATT", {
  skip_without_engine()
  d <- make_panel(seed = 3, n = 60)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 15, num_burnin = 8, num_chains = 1,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 11)
  pred <- predict(fit, x = d$x, z = d$z, t = d$t)

  # The R getters delegate to the engine rather than reimplementing the
  # alignment, so this must be exact, not merely close.
  np <- reticulate::import("numpy", convert = TRUE)
  py_att <- as.numeric(np$asarray(pred$py_pred$att()[["att"]]))
  expect_equal(get_att(pred)$att, py_att, tolerance = 1e-12)
})

test_that("saveRDS round-trips without a live Python handle", {
  skip_without_engine()
  d <- make_panel(seed = 5, n = 40)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 15, num_burnin = 8, num_chains = 1,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 9)
  before <- predict(fit, x = d$x, z = d$z, t = d$t)

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  reloaded <- readRDS(path)

  # No manual surgery: the cached handle lives in an environment that does not
  # survive serialization, and rehydration reads the raw vector.
  after <- predict(reloaded, x = d$x, z = d$z, t = d$t)
  expect_equal(after$tauhats, before$tauhats, tolerance = 1e-5)
  expect_equal(after$muhats0, before$muhats0, tolerance = 1e-5)
})

test_that("stability reports numbers and withholds the verdict", {
  skip_without_engine()
  d <- make_panel(seed = 7, n = 60)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 20, num_burnin = 10, num_chains = 2,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 4)
  pred <- predict(fit, x = d$x, z = d$z, t = d$t)

  stab <- att_stability(pred, warn = FALSE)
  expect_s3_class(stab$summary, "data.frame")
  expect_true(stab$summary$ess_median > 0)
  expect_equal(stab$summary$num_chains, 2)
  expect_true(is.na(stab$is_reliable))
  expect_match(stab$verdict_note, "coverage")
  expect_true(all(c("ess_bulk", "rhat", "mcse") %in% names(stab$by_exposure)))
})

test_that("multi-chain fits expose R-hat", {
  skip_without_engine()
  d <- make_panel(seed = 11, n = 50)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 15, num_burnin = 8, num_chains = 3,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 2)
  expect_equal(ncol(fit$beta_values), 45L)   # 3 chains x 15 draws
  pred <- predict(fit, x = d$x, z = d$z, t = d$t)
  expect_equal(dim(pred$tauhats)[3], 45L)
  stab <- att_stability(pred, warn = FALSE)
  expect_false(is.na(stab$summary$rhat_max))
})

test_that("att_stability reproduces the live numbers from a serialized prediction", {
  skip_without_engine()
  d <- make_panel(seed = 23, n = 60)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 20, num_burnin = 10, num_chains = 3,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 6)
  pred <- predict(fit, x = d$x, z = d$z, t = d$t)
  live <- att_stability(pred, warn = FALSE)

  # The Python handle does not survive serialization, so this exercises the
  # fallback that rebuilds the (chain, draw, exposure) array from att_full.
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(pred, path)
  dead <- readRDS(path)
  expect_true(reticulate::py_is_null_xptr(dead$py_pred))

  revived <- att_stability(dead, warn = FALSE)

  # Exact, not approximate: the fallback must recover the same chain layout.
  # A draw-major unflattening would mix the chains and deflate R-hat instead.
  expect_equal(revived$summary$num_chains, live$summary$num_chains)
  expect_equal(revived$summary$rhat_max, live$summary$rhat_max, tolerance = 1e-10)
  expect_equal(revived$summary$ess_min, live$summary$ess_min, tolerance = 1e-10)
  expect_equal(revived$by_exposure$rhat, live$by_exposure$rhat, tolerance = 1e-10)
  expect_equal(revived$by_exposure$n_treated, live$by_exposure$n_treated)
})

test_that("att_stability survives a summary_only prediction", {
  skip_without_engine()
  d <- make_panel(seed = 31, n = 50)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 15, num_burnin = 8, num_chains = 2,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 3)
  # summary_only drops the per-cell draws but keeps att_full, which is what the
  # diagnostic reads -- so the chapter's memory-saving advice stays compatible
  # with running the diagnostic.
  pred <- predict(fit, x = d$x, z = d$z, t = d$t, summary_only = TRUE)
  expect_null(pred$tauhats)
  stab <- att_stability(pred, warn = FALSE)
  expect_equal(stab$summary$num_chains, 2)
  expect_false(is.na(stab$summary$rhat_max))
})
