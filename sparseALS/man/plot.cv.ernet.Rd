\name{plot.cv.ernet}
\alias{plot.cv.ernet}
\title{Plot the cross-validated curve produced by cv.ernet}
\description{
Plots the cross-validated curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv.glmnet} function from the \code{glmnet} package.}

\usage{
\method{plot}{cv.ernet}(x, sign.lambda, ...)
}

\arguments{
	\item{x}{fitted \code{\link{cv.ernet}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}

\details{A plot is produced.}

\author{Yi Yang, Yuwen Gu and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.ernet}}.}

\keyword{models}
\keyword{regression}
