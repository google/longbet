% Generated from R/rollout.R -- regenerate with roxygen2::roxygenise().
\name{plot_rollout}
\alias{plot_rollout}
\title{Plot a staggered rollout}
\usage{plot_rollout(z, t = NULL, labels = NULL,
  never_treated_label = "Never treated", colors = ROLLOUT_COLORS,
  title = "Adoption cohorts over time", show_counts = TRUE)}
\arguments{
\item{z, t, labels, never_treated_label}{As in \code{\link{rollout_summary}}.}
\item{colors}{Named character vector of fill colours by status.}
\item{title}{Plot title, or \code{NULL} for none.}
\item{show_counts}{Append each cohort's unit count to its axis label.}
}
\value{A \code{ggplot} object.}
\description{
Draws the adoption pattern as a tile chart, one row per adoption cohort. The
result is an ordinary \code{ggplot}, so titles, palettes and themes can be
changed in the usual way.
}
\seealso{\code{\link{rollout_summary}}, which returns the underlying table and
needs no plotting package.}
