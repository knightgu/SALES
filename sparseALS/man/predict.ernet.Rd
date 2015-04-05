\name{predict.ernet}
\alias{predict.ernet}
\alias{predict.alspath}
\title{Make predictions from an "ernet" object.}

\description{
Similar to other predict methods, this functions predicts fitted values from a fitted \code{\link{ernet}} object.
}

\usage{
\method{predict}{ernet}(object, newx, s = NULL, type = "response", ...)
}

\arguments{
		\item{object}{fitted \code{\link{ernet}} model object.}
    
		\item{newx}{matrix of new values for \code{x} at which predictions are to be made. NOTE: \code{newx} must be a matrix, \code{predict} function does not accept a vector or other formats of \code{newx}.}
    
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are to be made. Default is the entire sequence used to
		create the model.}
    
		\item{type}{
		type of prediction required. Only \code{response} is available. Gives predicted response for regression problems.}
    
		\item{\dots}{Not used. Other arguments to predict. } }
\details{\code{s} is the new vector at which predictions are to be made. If \code{s} is not in the lambda sequence used for fitting the model, the \code{predict} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \code{lambda} indices. }
\value{The object returned depends on type.}


\author{Yuwen Gu, Yi Yang, and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>}
\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{coef}} method}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
