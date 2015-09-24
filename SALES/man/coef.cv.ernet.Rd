\name{coef.cv.ernet}
\alias{coef.cv.ernet}
\title{Get coefficients from a "cv.ernet" object.}
\description{
This function gets coefficients from a cross-validated ernet model,
using the stored \code{"ernet.fit"} object, and the optimal value
chosen for \code{lambda}.
}

\usage{
\method{coef}{cv.ernet}(object, s = c("lambda.1se", "lambda.min"),...)
}

\arguments{
	\item{object}{fitted \code{\link{cv.ernet}} object.}
  
	\item{s}{value(s) of the penalty parameter \code{lambda} at which predictions are required. Default is the value \code{s="lambda.1se"} stored on the CV \code{object}, it is the largest value of \code{lambda} such that error is within 1 standard error of the minimum. Alternatively \code{s="lambda.min"} can be used, it is the optimal value of \code{lambda} that gives minimum cross validation error \code{cvm}. If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be used.}
  
	\item{\dots}{not used. Other arguments to predict. } 
}

\details{This function makes it easier to use the results of
cross-validation to get coefficients or make coefficient predictions.}
\value{The object returned depends the \dots argument which is passed on to the \code{\link{predict}} method for \code{\link{ernet}} objects.}

\author{ Yuwen Gu, Yi Yang, and Hui Zou\cr
Maintainer:  Yuwen Gu <guxxx192@umn.edu>}

\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.ernet}}, and \code{\link{predict.cv.ernet}} methods.} 

\keyword{models}
\keyword{regression}