encourage_panel <- function() {
  z <- matrix(0, 6, 4); z[1:3, 2:4] <- 1
  d <- rbind(c(1,1,1,1), c(0,1,1,1), c(0,0,0,1),
             c(1,1,1,1), c(0,0,0,0), c(0,0,1,1))
  list(z = z, d = d, y = d * 2 + row(z))
}

test_that("encouragement metadata accepts full compliance and no pre-periods", {
  skip_without_engine()
  z <- matrix(c(1, 1, 0, 0), ncol = 1)
  metadata <- validate_encouragement(z, z, t = 7)
  expect_true(metadata$perfect_compliance)
  expect_false(metadata$has_pre_periods)
  expect_equal(metadata$encouragement_period, 7)
  expect_equal(metadata$n_encouraged, 2)
})

test_that("encouragement summaries and estimates equal the Python results", {
  skip_without_engine()
  p <- encourage_panel()
  lb <- longbet:::longbet_py()
  summary <- encouragement_summary(p$z, p$d, alpha = 0.1)
  py_summary <- lb$encouragement_summary(z = p$z, d = p$d, alpha = 0.1)
  py_summary <- as.data.frame(longbet:::.encouragement_as_r(py_summary)); rownames(py_summary) <- NULL
  attr(py_summary, "pandas.index") <- NULL
  expect_equal(summary, py_summary, tolerance = 0)
  expect_null(attr(summary, "pandas.index"))
  expect_identical(unserialize(serialize(summary, NULL)), summary)
  expect_equal(summary$median_observed_lag[2], -0.5)
  expect_equal(summary$n_encouraged_adopters[2], 2)

  result <- encouragement_effects(p$y, p$d, p$z, alpha = 0.1)
  py_result <- lb$encouragement_effects(y = p$y, d = p$d, z = p$z, alpha = 0.1)
  py_result <- as.data.frame(longbet:::.encouragement_as_r(py_result)); rownames(py_result) <- NULL
  attr(py_result, "pandas.index") <- NULL
  expect_equal(result, py_result, tolerance = 0)
  expect_null(attr(result, "pandas.index"))
  expect_identical(unserialize(serialize(result, NULL)), result)
  expect_equal(result$period_index, 1:3)
  expect_true(all(result$inference == "normal_ar"))
})

test_that("infinite and missing confidence-set components survive conversion", {
  skip_without_engine()
  z <- matrix(c(1, 1, 0, 0), ncol = 1)
  d <- matrix(0, 4, 1)
  result <- encouragement_effects(matrix(c(1, 2, 1, 2), ncol = 1), d, z)
  expect_equal(result$wald_set_type, "all_real")
  expect_equal(result$wald_lower_1, -Inf)
  expect_equal(result$wald_upper_1, Inf)
  expect_true(is.na(result$wald))
  expect_true(is.na(result$wald_lower_2))

  z <- matrix(c(1, 0, 0), ncol = 1)
  result <- encouragement_effects(3 * z, z, z)
  expect_equal(result$wald, 3)
  expect_equal(result$wald_set_type, "unavailable")
  expect_equal(result$wald_reason, "insufficient_arm_size")
  expect_true(is.na(result$itt_y_se))
})

test_that("negative stages are retained and incomplete outcomes rejected", {
  skip_without_engine()
  p <- encourage_panel()
  result <- encouragement_summary(p$z, p$d[c(4,5,6,1,2,3), ])
  expect_lt(result$first_stage[2], 0)
  p$y[1,2] <- NA_real_
  expect_error(encouragement_effects(p$y, p$d, p$z), "complete and finite")
  p$d[1, ] <- c(0,1,0,1)
  expect_error(validate_encouragement(p$z, p$d), "d:.*switches back")
})

test_that("encouragement utility formals and output columns follow the contract", {
  skip_if_not_installed("yaml")
  path <- system.file("contract", "longbet-api.yaml", package = "longbet")
  if (!nzchar(path)) path <- file.path("..", "..", "contract", "longbet-api.yaml")
  api <- yaml::read_yaml(path)$encouragement_api
  for (name in names(api$functions)) {
    spec <- api$functions[[name]]
    args <- formals(get(name, mode = "function"))
    expect_setequal(names(args), c(unlist(spec$required), names(spec$defaults)))
    for (key in names(spec$defaults)) {
      expect_equal(eval(args[[key]]), spec$defaults[[key]])
    }
  }
  skip_without_engine()
  p <- encourage_panel()
  expect_equal(names(validate_encouragement(p$z, p$d)), unlist(api$metadata_fields))
  expect_equal(names(encouragement_summary(p$z, p$d)), unlist(api$summary_columns))
  expect_equal(names(encouragement_effects(p$y, p$d, p$z)), unlist(api$effects_columns))
})

test_that("encouragement plot builds from assigned-arm summary", {
  skip_without_engine()
  skip_if_not_installed("ggplot2")
  p <- encourage_panel()
  plot <- plot_encouragement(p$z, p$d)
  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
})
