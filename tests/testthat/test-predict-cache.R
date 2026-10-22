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

test_that("cached GP projections preserve R predictions for each kernel", {
  skip_without_engine()
  d <- make_panel(n = 12, Tn = 4, t0 = 2)
  fit <- longbet(y = d$y, x = d$x, z = d$z, t = d$t,
                 num_sweeps = 6, num_burnin = 3, num_chains = 1,
                 num_trees_pr = 2, num_trees_trt = 2,
                 max_depth_pr = 2, max_depth_trt = 2,
                 min_points_per_leaf_pr = 1, min_points_per_leaf_trt = 1,
                 sigma_prior_a = 2, sigma_prior_b = 1, random_seed = 42)
  z_future <- cbind(d$z, d$z[, ncol(d$z)], d$z[, ncol(d$z)])
  t_future <- c(d$t, max(d$t) + 1:2)
  project <- function(lengthscale, ...) {
    predict(fit, x = d$x, z = z_future, t = t_future,
            lambda_knl = lengthscale, random_seed = 7,
            summary_only = TRUE, ...)
  }

  first <- project(1)
  expect_false(reticulate::py_has_attr(fit$handle$py_model, "_eval_cache"))
  cached <- project(1, cache_forest_evaluations = TRUE)
  expect_true(reticulate::py_has_attr(fit$handle$py_model, "_eval_cache"))
  fields <- c("att_full", "beta_values", "tauhats.mean", "tauhats.lower",
              "tauhats.upper", "muhats0.mean", "preds.mean")
  expect_identical(cached[fields], first[fields])

  # Changing the GP projection must reuse only forests, never the beta draws.
  second <- project(3)
  cached_second <- project(3, cache_forest_evaluations = TRUE)
  expect_identical(cached_second[fields], second[fields])
  expect_false(identical(cached$beta_values, cached_second$beta_values))
  expect_null(cached_second$tauhats)
})
