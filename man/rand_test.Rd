% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_clean.R
\name{rand_test}
\alias{rand_test}
\title{Randomization tests}
\usage{
rand_test(x, y, n = 1000)
}
\arguments{
\item{x}{A numeric vector.}

\item{y}{A numeric vector.}

\item{n}{The number of randomization, default is 1000.}
}
\value{
a data frame with mean, rank, p.value, etc.
}
\description{
Perform randomization test between two vectors.
}
