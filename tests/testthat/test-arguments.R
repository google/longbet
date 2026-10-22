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

test_that("pcat is refused rather than silently ignored", {
  skip_without_engine()
  d <- make_panel(n = 20, Tn = 5)
  expect_error(
    longbet(y = d$y, x = d$x, z = d$z, t = d$t, pcat = 2,
            num_sweeps = 4, num_burnin = 2, num_chains = 1),
    "pcat"
  )
  # pcat = 0 is what the reference interface's default meant, and is fine.
  expect_no_error(
    longbet(y = d$y, x = d$x, z = d$z, t = d$t, pcat = 0,
            num_sweeps = 4, num_burnin = 2, num_chains = 1, num_trees_pr = 3, num_trees_trt = 3)
  )
})

test_that("unknown arguments are refused", {
  skip_without_engine()
  d <- make_panel(n = 20, Tn = 5)
  expect_error(
    longbet(y = d$y, x = d$x, z = d$z, t = d$t, nonesuch = 1,
            num_sweeps = 4, num_burnin = 2, num_chains = 1),
    "Unsupported argument"
  )
})

test_that("x_trt reaches the treatment forest", {
  skip_without_engine()
  d <- make_panel(seed = 13, n = 40)
  x_trt <- matrix(rnorm(nrow(d$x) * 2), nrow(d$x), 2)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t, x_trt = x_trt,
                 num_sweeps = 10, num_burnin = 5, num_chains = 1,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 1)
  py <- rehydrate_jax_model(fit)
  names_ <- vapply(reticulate::py_to_r(py$design_$blocks), function(b) b$name, "")
  expect_true("x_trt" %in% names_)

  # And it is required at predict time, rather than silently dropped.
  expect_error(predict(fit, x = d$x, z = d$z, t = d$t), "x_trt")
  expect_no_error(predict(fit, x = d$x, z = d$z, t = d$t, x_trt = x_trt))
})

test_that("a propensity score must be supplied again at predict time", {
  skip_without_engine()
  d <- make_panel(seed = 17, n = 40)
  ps <- plogis(d$x[, 1])
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t, ps = ps,
                 num_sweeps = 10, num_burnin = 5, num_chains = 1,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 1)
  expect_error(predict(fit, x = d$x, z = d$z, t = d$t), "ps")
  expect_no_error(predict(fit, x = d$x, z = d$z, t = d$t, ps = ps))
})

test_that("non-absorbing treatment is refused", {
  skip_without_engine()
  d <- make_panel(n = 20, Tn = 5)
  z <- d$z
  z[1, ] <- c(0, 1, 1, 0, 0)
  expect_error(
    longbet(y = d$y, x = d$x, z = z, t = d$t, num_sweeps = 4, num_burnin = 2, num_chains = 1),
    "absorbing"
  )
})

test_that("the projection seed is honoured", {
  skip_without_engine()
  d <- make_panel(seed = 19, n = 30)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 10, num_burnin = 5, num_chains = 1,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 1)

  Tn <- ncol(d$z)
  z2 <- cbind(d$z, matrix(d$z[, Tn], nrow(d$z), 3))
  t2 <- seq_len(Tn + 3)

  a <- predict(fit, x = d$x, z = z2, t = t2, random_seed = 5)
  b <- predict(fit, x = d$x, z = z2, t = t2, random_seed = 5)
  c_ <- predict(fit, x = d$x, z = z2, t = t2, random_seed = 6)
  expect_equal(a$beta_values, b$beta_values)
  expect_false(isTRUE(all.equal(a$beta_values, c_$beta_values)))
})

test_that("device selection works and reports honestly", {
  skip_without_engine()
  devs <- longbet_devices()
  expect_gt(length(devs), 0)

  d <- make_panel(n = 20, Tn = 5)
  expect_no_error(
    longbet(y = d$y, x = d$x, z = d$z, t = d$t, device = "cpu",
            num_sweeps = 4, num_burnin = 2, num_chains = 1, num_trees_pr = 3, num_trees_trt = 3)
  )
  if (!any(grepl("gpu|cuda|rocm", devs))) {
    expect_error(
      longbet(y = d$y, x = d$x, z = d$z, t = d$t, device = "gpu",
              num_sweeps = 4, num_burnin = 2, num_chains = 1, num_trees_pr = 3, num_trees_trt = 3),
      "no GPU"
    )
  }
})

test_that("unbalanced panels are accepted", {
  skip_without_engine()
  d <- make_panel(seed = 23, n = 40)
  y <- d$y
  y[sample(length(y), 20)] <- NA
  fit <- longbet(y = y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 10, num_burnin = 5, num_chains = 1,
                 num_trees_pr = 5, num_trees_trt = 5, random_seed = 1)
  pred <- predict(fit, x = d$x, z = d$z, t = d$t)
  expect_true(all(is.finite(pred$tauhats.mean)))
})


test_that("internal ablation options are not part of the user-facing API", {
  # These stayed available in the engine for ablation work, but they are not
  # exposed here, so they must land in `...` and be refused by name rather than
  # be silently ignored. No engine is needed: the check is the first statement.
  d <- make_panel(n = 20, Tn = 5)
  internal <- c("split_time_ps", "split_calendar_trt", "split_exposure_trt",
                "tempering_levels", "tempering_beta_min")
  base <- list(y = d$y, x = d$x, z = d$z, t = d$t, num_chains = 1,
               num_burnin = 2, num_sweeps = 2, num_trees_pr = 2, num_trees_trt = 2)
  for (option in internal) {
    args <- base
    args[[option]] <- FALSE
    expect_error(do.call(longbet, args), option, fixed = TRUE)
    args$y <- list(d$y, d$y)
    expect_error(do.call(longbet_multi, args), option, fixed = TRUE)
  }
})
