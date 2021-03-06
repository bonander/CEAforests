% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explain_forest.R
\name{cate.prepare}
\alias{cate.prepare}
\title{Internal function.}
\usage{
cate.prepare(forest, subset = NULL, compliance.score = NULL)
}
\arguments{
\item{forest}{A trained causal forest or instrumental forest object.}
}
\value{
Returns a vector of transformed outcomes that can be used in parametric or semi-parametric regressions to explain the results from the forest.
}
\description{
\code{cate.prepare} Provides de-biased tau estimates based on AIPW (for use in custom explainers, e.g., regression models or GAMs).
}
\examples{
\dontrun{
To be added...
}
}
\keyword{internal}
