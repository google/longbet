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

#' Summarise a staggered rollout
#'
#' Describes who was treated and when, as one row per adoption cohort and
#' period. A cohort is the set of units sharing a first-treated period, derived
#' from `z` rather than declared, so this works whether or not the design has
#' named launch waves.
#'
#' This is a description of the *design*, so it needs no fitted model and can be
#' run before `longbet()`. The `exposure` column is computed by the same code
#' the sampler uses, so a plot of the rollout and the event-time axis of an ATT
#' cannot disagree.
#'
#' @param z Binary treatment indicator, `[N x T]`. Must be absorbing.
#' @param t Calendar time vector, length `T`. Defaults to `1:ncol(z)`.
#' @param labels Optional named vector of display names keyed by adoption time,
#'   e.g. `c("4" = "W1", "6" = "W2")`. Unlabelled cohorts fall back to their
#'   adoption time.
#' @param never_treated_label Name for the cohort that is never treated.
#' @return A data frame with columns `cohort`, `first_treated`, `n_units`,
#'   `period`, `period_index`, `status` and `exposure`. `cohort` and `status`
#'   are factors, ordered by adoption time and by legend order respectively.
#' @seealso [plot_rollout()]
#' @examples
#' \dontrun{
#' z <- matrix(0, 100, 10); z[1:50, 5:10] <- 1
#' rollout_summary(z)
#' }
#' @export
rollout_summary <- function(z, t = NULL,
                            labels = NULL,
                            never_treated_label = "Never treated") {
  lb <- longbet_py()
  z <- as.matrix(z)
  if (is.null(t)) t <- seq_len(ncol(z))
  t <- as.numeric(t)

  py_labels <- NULL
  if (!is.null(labels)) {
    keys <- suppressWarnings(as.numeric(names(labels)))
    if (anyNA(keys)) {
      stop("`labels` must be a named vector whose names are adoption times, ",
           "e.g. c(`4` = \"W1\").", call. = FALSE)
    }
    py_labels <- stats::setNames(as.list(as.character(labels)), as.character(keys))
    names(py_labels) <- as.character(keys)
  }

  df <- reticulate::py_to_r(
    lb$rollout_summary(
      z = reticulate::r_to_py(z),
      t = reticulate::r_to_py(t),
      labels = if (is.null(py_labels)) NULL else reticulate::r_to_py(py_labels),
      never_treated_label = never_treated_label
    )
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # reticulate hands back pandas Categoricals as plain character; restore the
  # ordering, which is what makes a plot come out in adoption order.
  cohort_order <- unique(df$cohort[order(df$first_treated, na.last = TRUE)])
  df$cohort <- factor(df$cohort, levels = cohort_order)
  df$status <- factor(df$status,
                      levels = c("Not yet treated", "Treated", "Never treated"))
  rownames(df) <- NULL
  df
}

#' Default fill colours for [plot_rollout()], by status
#' @export
ROLLOUT_COLORS <- c(
  "Not yet treated" = "#d9d9d9",
  "Treated"         = "#3b6ea5",
  "Never treated"   = "#9e9e9e"
)

#' Plot a staggered rollout
#'
#' Draws the adoption pattern as a tile chart, one row per adoption cohort. The
#' returned object is an ordinary `ggplot`, so titles, palettes and themes can
#' be changed in the usual way.
#'
#' @param z,t,labels,never_treated_label As in [rollout_summary()].
#' @param colors Named character vector of fill colours by status. Defaults to
#'   [ROLLOUT_COLORS].
#' @param title Plot title, or `NULL` for none.
#' @param show_counts Append each cohort's unit count to its axis label.
#' @return A `ggplot` object.
#' @seealso [rollout_summary()], which returns the underlying table and needs no
#'   plotting package.
#' @examples
#' \dontrun{
#' z <- matrix(0, 100, 10); z[1:50, 5:10] <- 1
#' plot_rollout(z, labels = c(`5` = "Wave 1"), never_treated_label = "Holdout")
#' }
#' @export
plot_rollout <- function(z, t = NULL,
                         labels = NULL,
                         never_treated_label = "Never treated",
                         colors = ROLLOUT_COLORS,
                         title = "Adoption cohorts over time",
                         show_counts = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_rollout() needs ggplot2.\n",
         "Install it with install.packages('ggplot2'), or use ",
         "rollout_summary() and plot the table with whatever you already have.",
         call. = FALSE)
  }
  df <- rollout_summary(z, t, labels = labels,
                        never_treated_label = never_treated_label)

  if (isTRUE(show_counts)) {
    counts <- tapply(df$n_units, df$cohort, function(v) v[1])
    levels(df$cohort) <- sprintf("%s (n=%d)", levels(df$cohort),
                                 as.integer(counts[levels(df$cohort)]))
  }
  # Cohorts read top-to-bottom in adoption order, so the staircase descends.
  df$cohort <- factor(df$cohort, levels = rev(levels(df$cohort)))

  present <- levels(droplevels(df$status))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$period, y = .data$cohort,
                                   fill = .data$status)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = colors[present], drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = unique(df$period)) +
    ggplot2::labs(x = "Period", y = NULL, fill = NULL, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid = ggplot2::element_blank())
}
