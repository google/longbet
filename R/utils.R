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

# Internal helpers shared by the R front door.

.as_np_vector <- function(x) {
  # reticulate converts length-one R vectors to Python scalars by default.
  reticulate::import("numpy", convert = FALSE)$atleast_1d(reticulate::r_to_py(x))
}

#' @keywords internal
.as_np_matrix <- function(x, name) {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric; got %s.", name, class(x)[1]), call. = FALSE)
  }
  reticulate::r_to_py(x)
}

#' @keywords internal
.as_np_array3 <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (length(dim(x)) != 3L) {
    stop(sprintf("`%s` must be a 3-D array with dimensions [N, T, P].", name), call. = FALSE)
  }
  reticulate::r_to_py(x)
}

#' @keywords internal
.as_np_forecast <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (is.matrix(x) || length(dim(x)) == 2L) {
    return(.as_np_matrix(x, name))
  }
  if (length(dim(x)) == 3L) {
    return(.as_np_array3(x, name))
  }
  stop(sprintf("`%s` must be a 2-D matrix [N, T] or 3-D array [N, T, K].", name), call. = FALSE)
}

#' @keywords internal
.reject_unsupported <- function(...) {
  extras <- list(...)
  if (length(extras) == 0L) return(invisible(NULL))
  nms <- names(extras)
  nms <- nms[nzchar(nms)]

  # pcat is the one argument of the reference C++ interface with no equivalent
  # here. Silently ignoring it would fit unordered categoricals as ordered
  # numerics, which is a different model, so refuse instead.
  if ("pcat" %in% nms && !identical(extras$pcat, 0)) {
    stop(
      "`pcat` is not supported. This engine has no notion of an unordered ",
      "categorical: every column is treated as ordered numeric.\n",
      "One-hot encode categorical variables with more than two levels before ",
      "calling longbet(), and note that this is a real modelling difference ",
      "from the reference C++ implementation.",
      call. = FALSE
    )
  }
  unknown <- setdiff(nms, "pcat")
  if (length(unknown)) {
    stop("Unsupported argument(s): ", paste(sQuote(unknown), collapse = ", "),
         ".\nSee ?longbet for the supported arguments.", call. = FALSE)
  }
  invisible(NULL)
}
