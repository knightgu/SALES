#' Make predictions from a cv.cpernet object
#'
#' This function makes predictions from a cross-validated cpernet model, using
#' the fitted \code{cv.cpernet} object, and the optimal value chosen for
#' \code{lambda}.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction.
#'
#' @param object fitted \code{\link{cv.cpernet}} object.
#' @param newx matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix. See documentation for \code{predict.cpernet}.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#' predictions are to be made. Default is the value \code{s = "lambda.1se"}
#' stored on the CV object. Alternatively \code{s = "lambda.min"} can be used.
#' If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
#' used.
#' @param \dots not used. Other arguments to predict.
#' @return The object returned depends the \dots{} argument which is passed on
#' to the \code{\link{predict}} method for \code{\link{cpernet}} objects.
#' @author Yuwen Gu and Hui Zou\cr Maintainer: Yuwen Gu <guxxx192@@umn.edu>
#' @seealso \code{\link{cv.cpernet}}, \code{\link{coef.cv.cpernet}}
#' @keywords models regression
#' @examples
#'
#' set.seed(1)
#' n <- 100
#' p <- 400
#' x <- matrix(rnorm(n*p), n, p)
#' y <- rnorm(n)
#' tau <- 0.30
#' pf <- abs(rnorm(p))
#' pf2 <- abs(rnorm(p))
#' w <- 2.0
#' lambda2 <- 1
#' m2.cv <- cv.cpernet(y = y, x = x, w = w, tau = tau, eps = 1e-8,
#'               pf.mean = pf, pf.scale = pf2,
#'               standardize = FALSE, lambda2 = lambda2)
#' as.vector(predict(m2.cv, newx = x, s = "lambda.min"))
#'
#' @export
predict.cv.cpernet <- function(object, newx, s = c("lambda.1se",
    "lambda.min"), ...) {
    if (is.numeric(s))
        lambda <- s else if (is.character(s)) {
        s <- match.arg(s)
        lambda <- object[[s]]
    } else stop("Invalid form for s")
    predict(object$cpernet.fit, newx, s = lambda, ...)
}
