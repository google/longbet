% Generated from R/rehydrate.R -- regenerate with roxygen2::roxygenise().
\name{rehydrate_jax_model}
\alias{rehydrate_jax_model}
\title{Rehydrate a saved LongBet model}
\usage{rehydrate_jax_model(object)}
\arguments{\item{object}{A fitted \code{longbet} object.}}
\value{A live Python \code{LongBet} instance.}
\description{
A fitted \code{longbet} object stores the engine state as a raw vector, never as
a live Python handle, so that \code{saveRDS()} and \code{readRDS()} round-trip.
This function returns a live Python model, deserialising on first use and
caching the handle in the object's own environment.
}
