\name{plot.ernet}
\alias{plot.ernet}
\title{Plot coefficients from an "ernet" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{ernet}} object. This function is modified based on the \code{plot} method in the \code{glmnet} package.
}
\usage{
\method{plot}{ernet}(x, xvar = c("norm", "lambda"), color = FALSE, label = FALSE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{ernet}} model}
		\item{xvar}{what is on the X-axis. \code{"norm"} plots against the
		L1-norm of the coefficients, \code{"lambda"} against the log-lambda
		sequence.}
		\item{color}{if \code{TRUE}, plot the curves with rainbow colors. \code{FALSE} is gray colors. Default is \code{FALSE}}
		\item{label}{if \code{TRUE}, label the curves with variable sequence numbers. Default is \code{FALSE}}
		\item{\dots}{other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced.
}
\author{Yuwen Gu, Yi Yang, and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>}
\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{plot.cv.ernet}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
