% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cea_forest.R
\name{predict.CEAforests}
\alias{predict.CEAforests}
\title{Predict with a CEA forest.}
\usage{
\method{predict}{CEAforests}(object, ...)
}
\arguments{
\item{object}{The trained CEA forest.}

\item{...}{Other options to be passed to grf::predict.causal_forest() or grf::predict.instrumental_forest(). See grf documentation for additional information.}
}
\value{
A matrix of predictions of conditional average treatment effects for the outcome and costs, along with variance estimates. Also returns debiased errors (estimates of the error of a forest with infinite size) and excess error due to Monte Carlo variability (estimated via jackknife). The latter provides an estimates of how unstable the estimates are if we grow forests of the same size on the same dataset. Increase the number of trees until the excess error becomes negligible. See grf::predict.causal_forest documentation for further details.
}
\description{
\code{predict.CEAforests} Gets estimates of conditional incremental outcomes and costs given x using a cea_forest object (a wrapper for grf::predict.causal_forest).
}
\examples{
\dontrun{
To be added...
}
}
