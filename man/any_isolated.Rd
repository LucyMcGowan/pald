% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pald_functions.R
\name{any_isolated}
\alias{any_isolated}
\title{Any isolated}
\usage{
any_isolated(c)
}
\arguments{
\item{c}{A \code{cohesion_matrix} object, a matrix of cohesion values
(see \code{\link{cohesion_matrix}}).}
}
\value{
Logical, indicating whether any points are isolated.
}
\description{
Checks for isolated points.
}
\examples{
d <- data.frame(
  x1 = c(1, 2, 3, 6),
  x2 = c(2, 1, 3, 10)
  )
D <- dist(d)
C <- cohesion_matrix(D)
any_isolated(C)
}
