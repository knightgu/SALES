\name{cv.ernet}
\alias{cv.ernet}
\title{Cross-validation for ernet}
\description{Does k-fold cross-validation for ernet, produces a plot, and returns a value for \code{lambda}. This function is based on the \code{cv} function from the \code{glmnet} package.}

\usage{
cv.ernet(x, y, lambda = NULL, pred.loss = "loss", 
         nfolds = 5, foldid, tau = 0.5, ...)
}

\arguments{
		\item{x}{\code{x} matrix as in \code{\link{ernet}}.}
		\item{y}{response variable \code{y} as in \code{\link{ernet}}.}
		\item{lambda}{optional user-supplied lambda sequence; default is \code{NULL}, and \code{\link{ernet}} chooses its own sequence.}
		\item{nfolds}{number of folds. Default value is 5. Although 
    \code{nfolds} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowed is 3.}
		\item{foldid}{an optional vector of values between 1 and \code{nfolds}, identifying what fold each observation is in. If supplied, \code{nfolds} will be supressed.}
		\item{pred.loss}{loss function used to calculate cross-validation error. The only option now is \code{"loss"}, which is the asymmetric squared error loss (ASEL).
		}
		\item{tau}{the asymmetry coefficient \eqn{\tau} used in the asymmetric squared error loss.}
		\item{\dots}{other arguments that can be passed to ernet.}
}

\details{The function runs \code{\link{ernet}} \code{nfolds}+1 times; the first to get the \code{lambda} sequence, and the remainder to compute the fit with each of the folds removed. The average error and standard deviation over the folds are computed. 
}

\value{an object of class \code{\link{cv.ernet}} is returned, which is a list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length \code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlower}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{number of non-zero coefficients at each \code{lambda}.}
		\item{name}{a text string indicating type of measure (for plotting purposes).}
		\item{ernet.fit}{a fitted \code{\link{gcdnet}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is within 1 standard error of the minimum.}
}

\author{Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>}
\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.} \url{http://www.jstatsoft.org/v33/i01/}}


\seealso{\code{\link{ernet}}, \code{\link{plot.cv.ernet}}, \code{\link{predict.cv.ernet}}, and \code{\link{coef.cv.ernet}} methods.}
\examples{
set.seed(2015)
}
\keyword{models}
\keyword{regression}

