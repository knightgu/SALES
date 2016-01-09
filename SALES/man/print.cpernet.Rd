\name{print.cpernet}
\alias{print.cpernet}
\title{Print a \code{cpernet} object}

\description{
Print a summary of the \code{cpernet} path at each step along the path.
}

\usage{
\method{print}{cpernet}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{fitted \code{\link{cpernet}} object.}
  \item{digits}{significant digits in the output.}
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

\examples{
set.seed(1)
n <- 100
p <- 400
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.30
pf <- abs(rnorm(p))
pf2 <- abs(rnorm(p))
w <- 2.0
lambda2 <- 1
m2 <- cpernet(y = y, x = x, w = w, tau = tau, eps = 1e-8, 
              pf.mean = pf, pf.scale = pf2,
              standardize = FALSE, lambda2 = lambda2)
print(m2)
}

\keyword{models}
\keyword{regression}

