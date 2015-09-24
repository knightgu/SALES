\name{print.cpernet}
\alias{print.cpernet}
\title{Print an "cpernet"" object}

\description{
Print a summary of the cpernet path at each step along the path.
}

\usage{
\method{print}{cpernet}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{fitted \code{\link{cpernet}} object.}
  \item{digits}{significant digits in printout.}
  \item{\dots}{additional print arguments.}
}

\details{
The call that produced the \code{\link{cpernet}} object is printed, followed by a three-column matrix with columns \code{Df1}, \code{Df2} and \code{Lambda}. The \code{Df1} and \code{Df2} columns are the number of nonzero mean and scale coefficients respectively.
}

\value{
a three-column matrix, the first two columns are the number of nonzero mean and scale coefficients respectively and the third column is \code{Lambda}.}

\author{
Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu <guxxx192@umn.edu>
}

\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}
}

\keyword{models}
\keyword{regression}

