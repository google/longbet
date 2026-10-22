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

#' Rehydrate a saved LongBet model
#'
#' A fitted `longbet` object stores the engine state as a raw vector, never as a
#' live Python handle, so that `saveRDS()` and `readRDS()` round-trip. This
#' function returns a live Python model, deserialising it on first use and
#' caching the handle so repeated `predict()` calls do not re-read the state.
#'
#' R serializes environments, so a cached handle *does* survive `saveRDS()` --
#' as a dangling external pointer that would fail on first use. The cache is
#' therefore validated before being returned: a handle from another R session,
#' or one whose pointer has been nulled, is discarded and rebuilt.
#'
#' @param object A fitted `longbet` object.
#' @return A live Python `LongBet` instance.
#' @export
rehydrate_jax_model <- function(object) {
  if (!inherits(object, "longbet") && !inherits(object, "longbet_multi")) {
    stop("rehydrate_jax_model() requires a fitted longbet or longbet_multi object.", call. = FALSE)
  }
  cache <- object$handle
  if (is.environment(cache) && .handle_is_live(cache)) {
    return(cache$py_model)
  }
  if (is.null(object$raw_model)) {
    stop("Cannot rehydrate: the model object carries no serialized engine state.",
         call. = FALSE)
  }

  lb <- longbet_py()
  tmp <- tempfile(fileext = ".npz")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(object$raw_model, tmp)
  py_model <- if (inherits(object, "longbet_multi")) {
    lb$LongBetMulti$load(tmp)
  } else {
    lb$LongBet$load(tmp)
  }

  if (is.environment(cache)) {
    cache$py_model <- py_model
    cache$pid <- Sys.getpid()
  }
  py_model
}

#' @keywords internal
.handle_is_live <- function(cache) {
  if (is.null(cache$py_model)) return(FALSE)
  # A handle created in a different R session cannot be valid here.
  if (!identical(cache$pid, Sys.getpid())) return(FALSE)
  ok <- tryCatch(!reticulate::py_is_null_xptr(cache$py_model),
                 error = function(e) FALSE)
  isTRUE(ok)
}
