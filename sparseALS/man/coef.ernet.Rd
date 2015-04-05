\name{coef.ernet}
\alias{coef.ernet}
\alias{coef.alspath}
\title{Get coefficients from an "ernet" object.}
\description{
Computes the coefficients or returns a list of the indices of the nonzero coefficients at the requested values for \code{lambda} from a fitted \code{\link{ernet}} object.
}
\usage{
\method{coef}{ernet}(object, s = NULL, type = c("coefficients", "nonzero"), ...)
}
\arguments{
		\item{object}{fitted \code{\link{ernet}} model object.}
    
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are to be made. Default is the entire sequence used to
		create the model.}
    
		\item{type}{type \code{"coefficients"} computes coefficients at the requested values for \code{s}. Type \code{"nonzero"} returns a list of the indices of nonzero coefficients for each value of \code{s}. Default is \code{"coefficients"}.}
    
		\item{\dots}{not used. Other arguments to predict. } 
}

\details{
\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{coef} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of coefficients from both left and right \code{lambda} indices. 
}
\value{The object returned depends on type.}


\author{Yuwen Gu, Yi Yang, and Hui Zou\cr
Maintainer: Yuwen Gu <guxxx192@umn.edu>}
\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{predict.ernet}} method}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
