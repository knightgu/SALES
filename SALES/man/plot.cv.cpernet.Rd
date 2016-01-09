\name{plot.cv.cpernet}
\alias{plot.cv.cpernet}
\title{Plot the cross-validated curve produced by \code{cv.cpernet}}
\description{
Plots the cross-validated curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv.glmnet} function from the \code{glmnet} package.}

\usage{
\method{plot}{cv.cpernet}(x, sign.lambda = 1, ...)
}

\arguments{
	\item{x}{fitted \code{\link{cv.cpernet}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}

\details{A plot is produced.}

\author{Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>}

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
m2.cv <- cv.cpernet(y = y, x = x, w = w, tau = tau, eps = 1e-8, 
              pf.mean = pf, pf.scale = pf2,
              standardize = FALSE, lambda2 = lambda2)
plot(m2.cv)
}

\seealso{\code{\link{cv.cpernet}}}

\keyword{models}
\keyword{regression}
