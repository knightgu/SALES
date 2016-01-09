\name{predict.cv.ernet}
\alias{predict.cv.ernet}
\title{make predictions from a \code{cv.ernet} object.}
\description{
This function makes predictions from a cross-validated ernet model,
using the fitted \code{cv.ernet} object, and the optimal value
chosen for \code{lambda}.
}

\usage{
\method{predict}{cv.ernet}(object, newx, s=c("lambda.1se", "lambda.min"), ...)}

\arguments{
		\item{object}{fitted \code{\link{cv.ernet}} object.}
    
		\item{newx}{matrix of new values for \code{x} at which predictions are to be made. Must be a matrix. See documentation for \code{predict.ernet}.}
    
		\item{s}{value(s) of the penalty parameter \code{lambda} at which predictions are to be made. Default is the value \code{s = "lambda.1se"} stored on the CV object. Alternatively \code{s = "lambda.min"} can be used. If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be used.}
    
		\item{\dots}{not used. Other arguments to predict.}} \details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on to the \code{\link{predict}} method for \code{\link{ernet}} objects.}

\author{Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>}

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
m1.cv <- cv.ernet(y = y, x = x, tau = tau, eps = 1e-8, pf = pf,
            pf2 = pf2, standardize = FALSE, intercept = FALSE, 
            lambda2 = lambda2)
as.vector(predict(m1.cv, newx = x, s = "lambda.min"))
}

\seealso{\code{\link{cv.ernet}}, \code{\link{coef.cv.ernet}}}

\keyword{models}
\keyword{regression}
