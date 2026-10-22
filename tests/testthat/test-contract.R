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

test_that("longbet() formals match the shared API contract", {
  skip_if_not_installed("yaml")
  path <- system.file("contract", "longbet-api.yaml", package = "longbet")
  if (!nzchar(path)) path <- file.path("..", "..", "contract", "longbet-api.yaml")
  skip_if_not(file.exists(path), "contract file not found")

  contract <- yaml::read_yaml(path)
  formals_lb <- formals(longbet)

  for (arg in contract$arguments) {
    if (!isTRUE(arg$expose_in_r)) next
    r_name <- arg$r_name
    expect_true(r_name %in% names(formals_lb),
                info = sprintf("longbet() is missing the contract argument '%s'", r_name))

    default <- formals_lb[[r_name]]
    if (is.call(default) && identical(as.character(default[[1]]), "c")) {
      # match.arg style: the first choice is the default
      default <- eval(default)[[1]]
    } else if (is.name(default)) {
      next
    } else {
      default <- eval(default)
    }
    expected <- arg$default
    if (is.null(expected)) {
      expect_null(default, info = sprintf("default for '%s'", r_name))
    } else if (is.numeric(expected) && is.numeric(default)) {
      expect_equal(as.numeric(default), as.numeric(expected), tolerance = 1e-12,
                   info = sprintf("default for '%s'", r_name))
    } else {
      expect_equal(default, expected, info = sprintf("default for '%s'", r_name))
    }
  }
})

test_that("every exposed contract argument really reaches the engine", {
  skip_without_engine()
  skip_if_not_installed("yaml")
  path <- system.file("contract", "longbet-api.yaml", package = "longbet")
  if (!nzchar(path)) path <- file.path("..", "..", "contract", "longbet-api.yaml")
  skip_if_not(file.exists(path), "contract file not found")

  contract <- yaml::read_yaml(path)
  lb <- reticulate::import("longbet")
  cfg <- lb$LongBetConfig()
  py_defaults <- reticulate::py_to_r(cfg$to_dict())

  for (arg in contract$arguments) {
    py_name <- arg$python_name
    expect_true(py_name %in% names(py_defaults),
                info = sprintf("LongBetConfig has no field '%s'", py_name))
  }
})
