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

encourage_model_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      set.seed(56)
      n <- 8L; periods <- 4L
      x <- matrix(rnorm(n * 2), n, 2)
      z <- matrix(0, n, periods); z[1:4, 2:4] <- 1
      d <- z; d[1, ] <- 1; d[4, ] <- 0; d[6, 3:4] <- 1
      y <- 2 * d + x[, 1] + matrix(rnorm(n * periods), n, periods)
      cfg <- list(num_chains = 2, num_sweeps = 6, num_burnin = 4,
                  num_trees_pr = 2, num_trees_trt = 2, random_seed = 111)
      fit <- longbet_encourage(y, d, z, x, config = cfg)
      cached <<- list(fit = fit, x = x, y = y, d = d, z = z,
                      groups = rep(c("a", "b"), 4))
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

test_that("encouragement R model preserves Python posterior values and axes", {
  skip_without_engine()
  p <- encourage_model_fixture()
  fit <- p$fit
  expect_s3_class(fit, "longbet_encourage")
  expect_type(fit$raw_archive, "raw")
  expect_false(contains_python_handle(fit))
  expect_equal(fit$config$sigma_prior_a, 2)
  expect_equal(fit$config$sigma_prior_b, 1)
  expect_identical(fit$metadata$engine, "adoption_clock")
  pred <- predict(fit, groups = p$groups, block_size = 5,
                   practical_threshold = 0.1, min_ess = 1, max_rhat = 100)
  expect_s3_class(pred, "longbet_encourage.pred")
  expect_false(contains_python_handle(pred))
  expect_identical(pred$group_labels, c("all", "a", "b"))
  expect_equal(dim(pred$draws$itt_y), c(3L, 3L, 2L, 6L))
  expect_equal(dim(pred$itt_covariance), c(18L, 18L))
  expect_equal(names(pred$itt_covariance_index), c("group", "horizon", "quantity"))
  expect_equal(encouragement_wald(pred), pred$wald)
  expect_equal(encouragement_wald(pred, FALSE), pred$wald_unchecked)
  expect_equal(encouragement_first_stage(pred), pred$first_stage)
  expect_equal(encouragement_comparison(pred), pred$comparison)
  expect_equal(pred$first_stage$practical_threshold, rep(.1, 9))
  expect_true(all(is.na(pred$comparison$disagreement)))
  engine <- longbet:::.encouragement_load(fit)
  py <- engine$predict(groups = p$groups, block_size = 5L)
  expect_equal(pred$draws, longbet:::.encouragement_plain(py$draws), tolerance = 0)
  expect_equal(pred$reference, longbet:::.encouragement_plain(py$reference), tolerance = 0)
  expect_output(print(fit), "adoption clock")
  expect_output(print(pred), "3 groups, 3 horizons")
  weighted <- predict(fit, weights = rep(c(1, 2), 4), min_ess = 1, max_rhat = 100)
  expect_equal(dim(weighted$draws$itt_y), c(1L, 3L, 2L, 6L))
  expect_false(isTRUE(all.equal(weighted$draws$itt_y, pred$draws$itt_y[1, , , , drop = FALSE])))
})

test_that("saveRDS predictions replay without live Python model handles", {
  skip_without_engine()
  p <- encourage_model_fixture()
  before <- predict(p$fit, groups = p$groups)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(fit = p$fit, pred = before), path)
  reloaded <- readRDS(path)
  expect_identical(reloaded$pred, before)
  after <- predict(reloaded$fit, groups = p$groups, block_size = 5)
  expect_equal(after$draws, before$draws, tolerance = 1e-5)
  expect_equal(after$reference, before$reference, tolerance = 0)
  expect_identical(after$metadata, before$metadata)
  expect_false(contains_python_handle(reloaded))

  output <- tempfile(fileext = ".rds")
  on.exit(unlink(output), add = TRUE)
  loader <- if (requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("longbet")) {
    sprintf("pkgload::load_all(%s, quiet=TRUE)",
            deparse(getNamespaceInfo(asNamespace("longbet"), "path")))
  } else "library(longbet)"
  script <- sprintf(paste0(
    "%s; a <- readRDS(%s); p <- predict(a$fit, groups=rep(c('a','b'),4)); ",
    "saveRDS(list(draws=p$draws, reference=p$reference, wald=encouragement_wald(a$pred)), %s)"
  ), loader, deparse(path), deparse(output))
  logs <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(script)),
                  stdout = TRUE, stderr = TRUE)
  expect_null(attr(logs, "status"), info = paste(logs, collapse = "\n"))
  fresh <- readRDS(output)
  expect_equal(fresh$draws, before$draws, tolerance = 0)
  expect_equal(fresh$reference, before$reference, tolerance = 0)
  expect_equal(fresh$wald, before$wald, tolerance = 0)
})

test_that("config options reject accidental coercion and retired options", {
  skip_without_engine()
  p <- encourage_model_fixture()
  expect_error(longbet_encourage(p$y, p$d, p$z, p$x, config = list(num_sweeps = 1.5)), "finite integer")
  expect_error(longbet_encourage(p$y, p$d, p$z, p$x, config = list(sigma_prior_a = 0)), "positive")
  expect_error(longbet_encourage(p$y, p$d, p$z, p$x, config = list(unknown = 1)), "LongBetConfig options")
  expect_error(longbet_encourage(p$y, p$d, p$z, p$x, config = list(first_stage = "probit")), "LongBetConfig options")
  expect_error(longbet_encourage(p$y, p$d, p$z, p$x, outcome = "ordinal"), "arg")
  expect_error(predict(p$fit, block_size = 1.5), "positive integer")
  expect_error(predict(p$fit, unsupported = TRUE), "unsupported|Unsupported")
  expect_error(predict_conditional(p$fit, cace_shrinkage = "adaptive"), "unsupported|Unsupported")
  expect_error(encouragement_wald(list()), "encouragement prediction")
  expect_error(encouragement_wald(structure(list(), class = "longbet_encourage.pred"), NA), "TRUE or FALSE")
})

test_that("design-based R inference preserves covariance labels and serialization", {
  skip_without_engine()
  p <- encourage_model_fixture()
  specification <- encouragement_design(blocks = rep(c("a", "b"), 4),
                                         block_weights = c(a = .4, b = .6))
  result <- design_encouragement_effects(p$y, p$d, p$z, design = specification)
  expect_false(contains_python_handle(result))
  expect_identical(unserialize(serialize(result, NULL)), result)
  expect_equal(names(result$covariance_index), c("contrast_id", "quantity"))
  expect_equal(dim(result$covariance), c(6L, 6L))
  expect_equal(result$covariance, t(result$covariance))
  lb <- reticulate::import("longbet", convert = FALSE)
  py_design <- do.call(lb$EncouragementDesign, unclass(specification))
  py <- lb$design_encouragement_effects(p$y, p$d, p$z, design = py_design)
  expect_equal(result$table, longbet:::.encouragement_plain(py$table), tolerance = 0)
  expect_equal(result$covariance, longbet:::.encouragement_plain(py$covariance$to_numpy()), tolerance = 0)
  bernoulli <- encouragement_design(assignment = "bernoulli", probabilities = .5)
  expect_false(contains_python_handle(design_encouragement_effects(
    p$y, p$d, p$z, design = bernoulli, groups = p$groups)))
})

test_that("staggered design probabilities retain their named cohort mapping", {
  skip_without_engine()
  z <- matrix(0, 9, 4)
  z[1:3, 2:4] <- 1
  z[4:6, 3:4] <- 1
  declared <- encouragement_design(assignment = "staggered",
    probabilities = c("1" = .3, "2" = .3, never = .4),
    cohort_weights = c("1" = .5, "2" = .5))
  result <- design_encouragement_effects(2 * z + row(z), z, z, design = declared)
  expect_true(nrow(result$table) > 2)
  expect_false(contains_python_handle(result))
  expect_error(encouragement_design(assignment = "staggered",
                                   probabilities = c(.5, .5)), "nonempty names")
})

test_that("paired bootstrap R forwarding retains unsuccessful runs", {
  skip_without_engine()
  p <- encourage_model_fixture()
  result <- encouragement_bootstrap_comparison(p$fit, replicates = 2, random_seed = 53,
                                                min_ess = 1, max_rhat = 100)
  expect_false(contains_python_handle(result))
  expect_equal(dim(result$replicates), c(2L, 1L, 3L, 2L, 2L))
  expect_true(all(is.na(result$table$disagreement)))
  expect_identical(unserialize(serialize(result, NULL)), result)
})

test_that("predict_conditional supports out-of-sample prediction and decisions", {
  skip_without_engine()
  p <- encourage_model_fixture()
  x_new <- matrix(rnorm(6 * ncol(p$x)), 6, ncol(p$x))
  cond <- predict_conditional(p$fit, new_x = x_new, cost = 2.0, hurdle = 0.5)
  expect_s3_class(cond, "longbet_conditional_pred")
  expect_equal(dim(cond$citt_y$mean), c(6L, 3L))
  expect_equal(dim(cond$citt_d$mean), c(6L, 3L))
  expect_equal(dim(cond$cace$median), c(6L, 3L))
  expect_equal(length(cond$breakeven_prob), 6L)
  expect_equal(length(cond$decision), 6L)
  expect_true(all(cond$breakeven_prob >= 0 & cond$breakeven_prob <= 1))
  expect_false(contains_python_handle(cond))
  study <- predict_conditional(p$fit)
  expect_equal(dim(study$citt_y$mean), c(8L, 3L))
  copy <- predict_conditional(p$fit, new_x = p$x[1:2, , drop = FALSE])
  expect_equal(copy$citt_y$draws, study$citt_y$draws[1:2, , , drop = FALSE], tolerance = 1e-6)
})

test_that("predict_conditional enforces monotonicity and the CACE floor", {
  skip_without_engine()
  p <- encourage_model_fixture()
  cond_mono <- predict_conditional(p$fit, monotonic_first_stage = TRUE)
  expect_true(all(cond_mono$citt_d$draws >= 0))
  expect_true(all(cond_mono$citt_d$mean >= 0))
  expect_true(all(is.finite(cond_mono$cace$median)))
  expect_equal(cond_mono$cace$draws, cond_mono$citt_y$draws / (cond_mono$citt_d$draws + 0.02),
               tolerance = 1e-6)
  cond_floor <- predict_conditional(p$fit, cace_stabilization = 0.05)
  expect_true(all(is.finite(cond_floor$cace$median)))
})

test_that("knapsack_policy respects budget and capacity constraints", {
  skip_without_engine()
  p <- encourage_model_fixture()
  cond <- predict_conditional(p$fit, cost = 0.5)
  pol_cap <- knapsack_policy(cond, cost = 0.5, capacity = 2)
  expect_type(pol_cap, "logical")
  expect_lte(sum(pol_cap), 2)
  pol_bud <- knapsack_policy(cond, cost = 1.0, budget = 2.5)
  expect_lte(sum(pol_bud), 2)
  cond_bud <- predict_conditional(p$fit, cost = 1.0, budget = 2.0, capacity = 1)
  expect_lte(sum(cond_bud$decision), 1)
  expect_true(is.numeric(cond_bud$policy_value))
  for (m in c("expected_net_value", "certainty_adjusted", "breakeven_probability")) {
    pol_m <- knapsack_policy(cond, cost = 0.5, capacity = 2, ranking_metric = m)
    expect_lte(sum(pol_m), 2)
  }
})

test_that("binary outcomes give probability-scale offer effects", {
  skip_without_engine()
  p <- encourage_model_fixture()
  y <- (p$y > median(p$y)) * 1
  cfg <- list(num_chains = 1, num_sweeps = 6, num_burnin = 4,
              num_trees_pr = 2, num_trees_trt = 2, random_seed = 7)
  fit <- longbet_encourage(y, p$d, p$z, p$x, outcome = "binary", config = cfg)
  expect_identical(fit$metadata$outcome, "binary")
  pred <- predict(fit, min_ess = 1, max_rhat = 100)
  expect_true(all(abs(pred$draws$itt_y) <= 1))
  cond <- predict_conditional(fit, new_x = p$x[1:3, , drop = FALSE])
  expect_true(all(abs(cond$citt_y$draws) <= 1))
})

test_that("principal strata estimation computes compliers, never-takers, and always-takers", {
  skip_without_engine()
  p <- encourage_model_fixture()
  cond <- predict_conditional(p$fit)
  df <- principal_strata(cond)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("period", "horizon", "stratum", "prob_mean", "count_mean", "n_total") %in% names(df)))
  expect_setequal(unique(df$stratum), c("complier", "always_taker", "never_taker"))
  for (h in unique(df$horizon)) {
    sub <- df[df$horizon == h, ]
    expect_equal(sum(sub$prob_mean), 1.0, tolerance = 0.05)
    expect_equal(sum(sub$count_mean), sub$n_total[1], tolerance = 0.5)
  }
  first_h <- unique(df$horizon)[1]
  df_h0 <- principal_strata(cond, horizon = first_h)
  expect_equal(nrow(df_h0), 3L)
  expect_output(strata_summary(cond, horizon = first_h), "Principal Strata Estimates")
})
