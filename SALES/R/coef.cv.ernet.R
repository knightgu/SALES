#' Get coefficients from a cv.ernet object
#'
#' This function gets coefficients from a cross-validated ernet model, using
#' the fitted \code{cv.ernet} object, and the optimal value chosen for
#' \code{lambda}.
#'
#' @param object fitted \code{\link{cv.ernet}} object.
#' @param s value(s) of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is the value \code{s="lambda.1se"} stored on the CV
#'   \code{object}, it is the largest value of \code{lambda} such that error is
#'   within 1 standard error of the minimum. Alternatively \code{s="lambda.min"}
#'   can be used, it is the optimal value of \code{lambda} that gives minimum
#'   cross validation error \code{cvm}. If \code{s} is numeric, it is taken as
#'   the value(s) of \code{lambda} to be used.
#' @param \dots not used. Other arguments to predict.
#'
#' @details This function makes it easier to use the results of cross-validation
#'   to get coefficients or make coefficient predictions.
#'
#' @return The object returned depends the \dots{} argument which is passed on
#'   to the \code{\link{predict}} method for \code{\link{ernet}} objects.
#'
#' @author Yuwen Gu and Hui Zou\cr
#'
#'   Maintainer: Yuwen Gu <yuwen.gu@uconn.edu>
#'
#' @seealso \code{\link{cv.ernet}}, \code{\link{predict.cv.ernet}}
#'
#' @keywords models regression
#'
#' @examples
#'
#' set.seed(1)
#' n <- 100
#' p <- 400
#' x <- matrix(rnorm(n * p), n, p)
#' y <- rnorm(n)
#' tau <- 0.90
#' pf <- abs(rnorm(p))
#' pf2 <- abs(rnorm(p))
#' lambda2 <- 1
#' m1.cv <- cv.ernet(y = y, x = x, tau = tau, eps = 1e-8, pf = pf,
#'                   pf2 = pf2, standardize = FALSE, intercept = FALSE,
#'                   lambda2 = lambda2)
#' as.vector(coef(m1.cv, s = "lambda.min"))
#'
#' @export
coef.cv.ernet <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s))
    lambda <- s else if (is.character(s)) {
                  s <- match.arg(s)
                  lambda <- object[[s]]
                } else stop("Invalid form for s")
  coef(object$ernet.fit, s = lambda, ...)
}
