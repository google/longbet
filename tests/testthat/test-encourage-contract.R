test_that("encouragement model and design formals match the shared contract", {
  skip_if_not_installed("yaml")
  path <- system.file("contract", "longbet-api.yaml", package = "longbet")
  if (!nzchar(path)) path <- file.path("..", "..", "contract", "longbet-api.yaml")
  expect_true(file.exists(path))
  contract <- yaml::read_yaml(path)
  model <- contract$encouragement_model_api
  design <- contract$encouragement_design_api
  check_defaults <- function(fn, defaults) {
    actual <- formals(fn)
    for (name in names(defaults)) {
      expect_true(name %in% names(actual), info = name)
      value <- eval(actual[[name]])
      expected <- defaults[[name]]
      # An empty YAML map is a named empty list; R list() has no names.
      if (is.list(expected) && length(expected) == 0L) expected <- list()
      expect_equal(value, expected, info = paste("default", name))
    }
  }
  check_defaults(longbet_encourage, model$r_fit$defaults)
  expect_setequal(names(formals(longbet_encourage)),
                  c(model$r_fit$required, names(model$r_fit$defaults)))
  check_defaults(getS3method("predict", "longbet_encourage"), model$r_predict$defaults)
  expect_setequal(names(formals(getS3method("predict", "longbet_encourage"))),
                  c(model$r_predict$required, names(model$r_predict$defaults), "..."))
  check_defaults(encouragement_design, design$defaults)
  expect_setequal(names(formals(encouragement_design)), names(design$defaults))
  check_defaults(design_encouragement_effects, design$effects_defaults)
  expect_setequal(names(formals(design_encouragement_effects)),
                  c(design$required, names(design$effects_defaults)))
  defaults <- model$bootstrap_comparison$defaults
  defaults$callback <- NULL
  check_defaults(encouragement_bootstrap_comparison, defaults)
  expect_setequal(names(formals(encouragement_bootstrap_comparison)), c("object", names(defaults)))
  expect_identical(model$calibration_status, "not_established")
})
