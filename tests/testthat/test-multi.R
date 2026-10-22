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

test_that("variance priors must be proper", {
  skip_without_engine()
  d <- make_panel(n = 12)
  expect_error(longbet_multi(list(a = d$y, b = -d$y), d$x, d$z, sigma_prior_a = 0),
               "positive")
})

test_that("longbet_multi fits multiple continuous and binary outcomes", {
  skip_without_engine()
  d <- make_panel(n = 30)
  # The chapter needs two continuous and one binary outcome, in caller order.
  y1 <- d$y
  y2 <- matrix(as.numeric(d$y > mean(d$y)), nrow = nrow(d$y), ncol = ncol(d$y))

  fit <- longbet_multi(
    y = list(rev = y1, hours = -y1 + 0.1 * d$x[, 1], churn = y2),
    x = d$x,
    z = d$z,
    t = d$t,
    num_sweeps = 10,
    num_burnin = 5,
    num_chains = 2,
    num_trees_pr = 3,
    num_trees_trt = 3,
    sigma_prior_a = 2, sigma_prior_b = 1,
    outcome = c(churn = "binary", hours = "continuous", rev = "continuous"),
    random_seed = 42
  )

  expect_s3_class(fit, "longbet_multi")
  expect_identical(fit$sampler_semantics, "full_precision_sur_v1")
  expect_equal(fit$M, 3L)
  expect_equal(fit$outcome_names, c("rev", "hours", "churn"))
  expect_equal(fit$order, c(3L, 1L, 2L))
  expect_true(is.list(fit$fits))
  expect_s3_class(fit$fits$rev, "longbet")
  expect_s3_class(fit$fits$churn, "longbet")
  expect_equal(dim(fit$Gamma_draws), c(9L, 20L))
  # Scalar traces must flatten chain-major, just like vector and loading traces.
  np <- reticulate::import("numpy", convert = TRUE)
  py_model <- rehydrate_jax_model(fit)
  for (name in fit$outcome_names) {
    py_child <- py_model$fits$`__getitem__`(name)
    sigma2 <- as.array(np$asarray(py_child$trace$sigma2))
    expect_equal(as.numeric(fit$fits[[name]]$sigma0_draws)^2,
                 as.numeric(t(sigma2)), tolerance = 1e-6)
  }

  # Indexing
  expect_identical(fit["rev"]$sdy, fit$fits$rev$sdy)
  expect_identical(fit[1]$sdy, fit$fits$rev$sdy)

  # Prediction
  pred <- predict(fit, x = d$x, z = d$z, t = d$t)
  expect_s3_class(pred, "longbet_multi.pred")
  legacy_pred <- pred
  legacy_pred$sampler_semantics <- "recursive_sur_v1"
  expect_error(joint_prob(legacy_pred, list()), "Refit from the original data")
  expect_s3_class(pred["rev"], "longbet.pred")
  expect_s3_class(pred["churn"], "longbet.pred")

  # Effect draws
  eff_rev <- effect_draws(pred, outcome = "rev")
  expect_equal(dim(eff_rev), c(nrow(d$y), ncol(d$y), 20L))
  eff_churn <- effect_draws(pred, outcome = "churn")
  expect_equal(dim(eff_churn), c(nrow(d$y), ncol(d$y), 20L))
  expect_equal(eff_churn, pnorm(pred["churn"]$muhats0 + pred["churn"]$tauhats) -
                 pnorm(pred["churn"]$muhats0), tolerance = 1e-7)
  # Binary effect scale is in [-1, 1]
  expect_true(all(eff_churn >= -1.0 & eff_churn <= 1.0))

  # Joint probability
  conds <- list(
    rev = function(e) e > 0,
    hours = function(e) e < 0,
    churn = function(e) e < 0.1
  )
  p_joint <- joint_prob(pred, conds)
  expect_equal(dim(p_joint), c(nrow(d$y), ncol(d$y)))
  expect_true(all(p_joint >= 0 & p_joint <= 1))
  manual <- apply((eff_rev > 0) & (effect_draws(pred, "hours") < 0) &
                    (eff_churn < 0.1), c(1, 2), mean)
  expect_equal(p_joint, manual)
  expect_equal(joint_prob(pred, rev(conds)), manual)
  expect_error(joint_prob(pred, c(conds, conds[1])), "exactly match")
  bad_cells <- matrix(FALSE, nrow(d$y), ncol(d$y))
  bad_cells[1, 1] <- NA
  expect_error(joint_prob(pred, conds, cells = bad_cells), "without NA")
  expect_error(effect_draws(pred, 1.5), "finite integer")
  expect_error(effect_draws(predict(fit, d$x, d$z, summary_only = TRUE), "rev"), "full posterior")

  # Correlation
  R <- outcome_correlation(fit)
  expect_equal(dim(R), c(3L, 3L))
  expect_equal(unname(diag(R)), c(1.0, 1.0, 1.0))
  expect_equal(dimnames(R), list(fit$outcome_names, fit$outcome_names))

  # Rehydration / saveRDS roundtrip
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  loaded <- readRDS(path)
  pred_loaded <- predict(loaded, x = d$x, z = d$z, t = d$t)
  expect_equal(pred_loaded["rev"]$tauhats.mean, pred["rev"]$tauhats.mean, tolerance = 1e-5)

  # A genuinely fresh R process must use serialized bytes/arrays, not live handles.
  payload <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  on.exit(unlink(c(payload, output)), add = TRUE)
  saveRDS(list(fit = fit, child = fit["rev"], pred = pred, d = d, conds = conds), payload)
  loader <- if (requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("longbet")) {
    sprintf("pkgload::load_all(%s, quiet=TRUE)",
            deparse(getNamespaceInfo(asNamespace("longbet"), "path")))
  } else {
    "library(longbet)"
  }
  script <- sprintf(paste0(
    "%s; a <- readRDS(%s); ",
    "p <- predict(a$fit, a$d$x, a$d$z); ",
    "child <- predict(a$child, a$d$x, a$d$z); ",
    "saveRDS(list(joint=joint_prob(a$pred,a$conds), effect=effect_draws(a$pred,'churn'), ",
    "cor=outcome_correlation(a$fit), tau=p['rev']$tauhats, child=child$tauhats), %s)"
  ), loader, deparse(payload), deparse(output))
  logs <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(script)),
                  stdout = TRUE, stderr = TRUE)
  expect_null(attr(logs, "status"), info = paste(logs, collapse = "\n"))
  fresh <- readRDS(output)
  expect_equal(fresh$joint, manual)
  expect_equal(fresh$effect, eff_churn)
  expect_equal(fresh$cor, R)
  expect_equal(fresh$tau, pred["rev"]$tauhats)
  expect_equal(fresh$child, pred["rev"]$tauhats)
})
