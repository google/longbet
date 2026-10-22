\name{plot_encouragement}
\alias{plot_encouragement}
\title{Plot Adoption by Assigned Encouragement Arm}
\usage{plot_encouragement(z, d, t = NULL)}
\arguments{\item{z,d,t}{As in \code{\link{validate_encouragement}}.}}
\description{Plots assigned-arm adoption curves, including baseline periods, with the gap shaded and encouragement time marked. Requires ggplot2.}
\value{A ggplot object built from the Python-produced summary table.}
\seealso{\code{\link{encouragement_summary}}}
