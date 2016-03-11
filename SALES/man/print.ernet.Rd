\name{print.ernet}
\alias{print.ernet}
\title{Print an ernet object}

\description{
Print a summary of the ernet path at each step along the path.
}

\usage{
\method{print}{ernet}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{fitted \code{\link{ernet}} object.}
  \item{digits}{significant digits in the output.}
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

\examples{
set.seed(1)
n <- 100
p <- 400
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.90
pf <- abs(rnorm(p))
pf2 <- abs(rnorm(p))
lambda2 <- 1
m1 <- ernet(y = y, x = x, tau = tau, eps = 1e-8, pf = pf,
            pf2 = pf2, standardize = FALSE, intercept = FALSE, 
            lambda2 = lambda2)
print(m1)
}

\keyword{models}
\keyword{regression}

