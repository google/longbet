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

intensity_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      set.seed(77)
      n <- 24L; periods <- 6L; t0 <- 3L
      x <- matrix(rnorm(n * 2), n, 2)
      z <- matrix(0, n, periods); z[1:12, t0:periods] <- 1
      u <- rnorm(n)
      pulse <- c(1, 0.6, 0.3, 0.1)
      d <- 0.5 + 0.3 * matrix(rnorm(n * periods), n, periods) + 0.5 * u
      d[1:12, t0:periods] <- d[1:12, t0:periods] + rep(pulse, each = 12)
      d <- pmax(d, 0)
      lags <- function(l) cbind(matrix(0, n, l), d)[, seq_len(periods), drop = FALSE]
      y <- -0.3 * d + 0.5 * lags(1) + 0.4 * lags(2) + 0.2 * lags(3) + u + 0.5 * x[, 1] +
        0.3 * matrix(rnorm(n * periods), n, periods)
      strata <- rep(c("north", "south"), length.out = n)
      cached <<- list(x = x, y = y, d = d, z = z, strata = strata, t = seq_len(periods))
    }
    cached
  }
})

contains_python_handle <- function(x) {
  if (inherits(x, "python.builtin.object") || is.environment(x) || typeof(x) == "externalptr") return(TRUE)
  if (is.list(x) && any(vapply(x, contains_python_handle, logical(1)))) return(TRUE)
  attrs <- attributes(x)
  if (length(attrs) && any(vapply(attrs, contains_python_handle, logical(1)))) return(TRUE)
  FALSE
}

test_that("intensity_reference returns plain-R design-based tables", {
  skip_without_engine()
  p <- intensity_fixture()
  meta <- validate_intensity(p$z, p$d)
  expect_identical(meta$treatment, "intensity")
  expect_equal(meta$encouragement_index, 2L)
  expect_true(meta$dose_nonabsorbing)
  ref <- intensity_reference(p$y, p$d, p$z, strata = p$strata, x = p$x)
  expect_s3_class(ref, "intensity_reference")
  expect_false(contains_python_handle(ref))
  expect_equal(nrow(ref$itt), 6L)
  expect_equal(nrow(ref$kernel), 4L)
  expect_true(all(c("kappa", "kappa_se", "kappa_cum", "kappa_cum_se") %in% names(ref$kernel)))
  expect_true(all(c("omega_cum", "omega_cum_set_type", "omega_cum_lower_1") %in% names(ref$itt)))
  expect_equal(dim(ref$covariance), c(8L, 8L))
  expect_equal(names(ref$covariance_index), c("quantity", "period_index"))
  expect_true(is.finite(ref$condition_number))
  expect_true(is.finite(ref$dip$z_stat))
  expect_true(ref$fosd$pvalue >= 0 && ref$fosd$pvalue <= 1)
  expect_equal(ref$metadata$n_strata, 2L)
  expect_equal(ref$metadata$lag_horizon, 3L)
  with_draws <- intensity_reference(p$y, p$d, p$z, kappa_draws = matrix(rnorm(40), 10, 4))
  expect_equal(nrow(with_draws$discrepancy), 4L)
  expect_output(print(ref), "Design-based reference")
  expect_error(intensity_reference(p$y, p$d, p$z, lag_horizon = -1), "lag_horizon")
})

test_that("longbet_encourage(treatment = 'intensity') fits, predicts and round-trips", {
  skip_without_engine()
  p <- intensity_fixture()
  cfg <- list(num_chains = 2, num_sweeps = 6, num_burnin = 4, num_trees_pr = 2, num_trees_trt = 2,
              random_seed = 11)
  expect_error(longbet_encourage(p$y, p$d, p$z, p$x, identification = "unconfounded", config = cfg),
               "treatment = 'intensity'")
  fit <- longbet_encourage(p$y, p$d, p$z, p$x, treatment = "intensity",
                           identification = "control_function", stage2_backend = "linear",
                           strata = p$strata, config = cfg)
  expect_s3_class(fit, "longbet_encourage")
  expect_false(contains_python_handle(fit))
  expect_identical(fit$metadata$engine, "distributed_lag")
  expect_identical(fit$metadata$tier, "T1")
  expect_equal(fit$metadata$lag_horizon, 3L)
  expect_equal(fit$metadata$n_strata, 2L)
  expect_identical(fit$metadata$component_streams, "independent")
  expect_length(fit$metadata$extra_columns, 3L * (1L + ncol(p$x)))
  expect_output(print(fit), "continuous intensity")

  pred <- predict(fit, min_ess = 1, max_rhat = 100)
  expect_s3_class(pred, "longbet_encourage.pred")
  expect_false(contains_python_handle(pred))
  expect_setequal(names(pred$draws), c("itt_y", "itt_d", "wald", "omega_cum"))
  expect_equal(dim(pred$draws$itt_y), c(1L, 4L, 2L, 6L))
  expect_true(all(c("kappa", "kappa_acr") %in% pred$kernel$quantity))
  expect_equal(nrow(pred$kernel), 8L)
  expect_true(all(c("reference_kappa", "posterior_mean") %in% names(pred$kernel)))
  expect_equal(nrow(pred$discrepancy), 4L)
  expect_identical(unique(pred$discrepancy$label), "T1")
  expect_equal(nrow(pred$design_diagnostics), 1L)
  expect_true(all(c("condition_number", "z_dip", "fosd_pvalue", "z_discrepancy") %in%
                    names(pred$design_diagnostics)))
  expect_true(all(c("omega_cum_set_type", "posterior_median") %in% names(pred$cumulative_return)))
  expect_true(all(c("itt_y_se", "omega_cum", "condition_number") %in% names(pred$reference)))

  cond <- predict_conditional(fit, cost = 0.1)
  expect_s3_class(cond, "longbet_conditional_pred")
  expect_false(contains_python_handle(cond))
  expect_equal(dim(cond$citt_y$draws), c(24L, 4L, 12L))
  expect_equal(dim(cond$kernel$draws), c(24L, 4L, 12L))
  expect_equal(cond$kernel_cum$mean, t(apply(cond$kernel$mean, 1, cumsum)))
  expect_equal(nrow(cond$kernel_summary), 4L)
  expect_equal(cond$cace$mean, cond$omega_cum$mean)
  expect_null(cond$strata)
  expect_type(cond$decision, "logical")

  # Round trip through saveRDS in a fresh object and through the engine's archive.
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  again <- readRDS(path)
  pred2 <- predict(again, min_ess = 1, max_rhat = 100)
  expect_equal(pred2$draws$itt_y, pred$draws$itt_y, tolerance = 0)
  expect_equal(pred2$kernel, pred$kernel, tolerance = 0)
  engine <- longbet:::.encouragement_load(fit)
  py <- engine$predict()
  expect_equal(pred$draws, longbet:::.encouragement_plain(py$draws), tolerance = 0)
  expect_equal(pred$reference, longbet:::.encouragement_plain(py$reference), tolerance = 0)
})
