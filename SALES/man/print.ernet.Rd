\name{print.ernet}
\alias{print.ernet}
\title{Print an "ernet"" object}

\description{
Print a summary of the ernet path at each step along the path.
}

\usage{
\method{print}{ernet}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{fitted \code{\link{ernet}} object.}
  \item{digits}{significant digits in printout.}
  \item{\dots}{additional print arguments.}
}

\details{
The call that produced the \code{\link{ernet}} object is printed, followed by a two-column matrix with columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of nonzero coefficients.
}

\value{
a two-column matrix, the first columns is the number of nonzero coefficients and the second column is \code{Lambda}.
}

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

