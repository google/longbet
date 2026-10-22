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

test_that("rollout_summary derives cohorts from z", {
  skip_without_engine()
  set.seed(1); N <- 200; Tn <- 14
  LAUNCH <- c(W1 = 4, W2 = 6, W3 = 8, W4 = 10)
  wave <- sample(c(names(LAUNCH), "Holdout"), N, TRUE, prob = c(.15,.15,.15,.15,.4))
  z <- matrix(0, N, Tn)
  for (w in names(LAUNCH)) z[wave == w, LAUNCH[[w]]:Tn] <- 1

  df <- rollout_summary(z, 1:Tn,
                        labels = c(`4` = "W1", `6` = "W2", `8` = "W3", `10` = "W4"),
                        never_treated_label = "Holdout")

  expect_s3_class(df, "data.frame")
  expect_equal(levels(df$cohort), c("W1", "W2", "W3", "W4", "Holdout"))
  expect_equal(nrow(df), 5L * Tn)
  expect_setequal(names(df), c("cohort", "first_treated", "n_units", "period",
                               "period_index", "status", "exposure"))

  # Cohort sizes recover the assignment.
  counts <- tapply(df$n_units, df$cohort, function(v) v[1])
  expect_equal(sum(counts), N)
  expect_equal(unname(counts[["Holdout"]]), sum(wave == "Holdout"))

  # Each wave is untreated up to its launch week and treated from then on.
  for (w in names(LAUNCH)) {
    rows <- df[df$cohort == w, ]
    rows <- rows[order(rows$period), ]
    expect_equal(rows$status == "Treated", rows$period >= LAUNCH[[w]],
                 info = w)
  }
  hold <- df[df$cohort == "Holdout", ]
  expect_true(all(hold$status == "Never treated"))
  expect_true(all(hold$exposure == 0))
})

test_that("rollout_summary works without labels and without a holdout", {
  skip_without_engine()
  z <- matrix(0, 20, 6); z[1:10, 3:6] <- 1; z[11:20, 5:6] <- 1
  df <- rollout_summary(z)
  expect_equal(levels(df$cohort), c("3", "5"))
  expect_false("Never treated" %in% levels(droplevels(df$status)))
})

test_that("exposure matches the sampler's own definition", {
  skip_without_engine()
  # Uneven calendar time: exposure is elapsed time, not a count of periods.
  z <- matrix(0, 4, 3); z[1:2, 2:3] <- 1
  df <- rollout_summary(z, t = c(2, 5, 9))
  cohort <- df[df$cohort == "5", ]
  cohort <- cohort[order(cohort$period), ]
  expect_equal(cohort$exposure, c(0, 3, 7))
})

test_that("plot_rollout returns a ggplot", {
  skip_without_engine()
  skip_if_not_installed("ggplot2")
  z <- matrix(0, 40, 8); z[1:20, 4:8] <- 1
  p <- plot_rollout(z, 1:8, never_treated_label = "Holdout")
  expect_s3_class(p, "ggplot")
  # It builds, which is where a bad factor level or missing colour would show.
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("bad input is refused", {
  skip_without_engine()
  expect_error(rollout_summary(matrix(0, 3, 4), t = 1:3), "length")
  expect_error(rollout_summary(matrix(0, 3, 4), labels = c(nope = "X")),
               "adoption times")
})

test_that("a single-period panel works", {
  skip_without_engine()
  # as.matrix() turns a vector into one column, so this is a T = 1 panel --
  # and a length-1 t crosses to Python as a scalar, which used to crash.
  df <- rollout_summary(c(1, 0, 1))
  expect_equal(nrow(df), 2L)
  expect_equal(levels(df$cohort), c("1", "Never treated"))
  expect_equal(df$exposure[df$cohort == "1"], 1L)
})
