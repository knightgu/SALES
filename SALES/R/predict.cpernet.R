#' Make predictions from a cpernet object
#'
#' Similar to other predict methods, this function predicts fitted values from
#' a cpernet object.
#'
#' @aliases predict.cpernet predict.cpalspath
#'
#' @param object fitted \code{\link{cpernet}} model object.
#' @param newx matrix of new values for \code{x} at which predictions are to be
#'   made. NOTE: \code{newx} must be a matrix, \code{predict} function does not
#'   accept a vector or other formats of \code{newx}.
#' @param s value(s) of the penalty parameter \code{lambda} at which predictions
#'   are to be made. Default is the entire sequence used to create the model.
#' @param type type of prediction required. Only \code{response} is available.
#'   Gives predicted response for regression problems.
#' @param \dots Not used. Other arguments to predict.
#'
#' @details \code{s} is the new vector at which predictions are to be made. If
#'   \code{s} is not in the lambda sequence used for fitting the model, the
#'   \code{predict} function will use linear interpolation to make predictions.
#'   The new values are interpolated using a fraction of predicted values from
#'   both left and right \code{lambda} indices.
#'
#' @return The object returned depends on type.
#'
#' @author Yuwen Gu and Hui Zou\cr
#'
#'   Maintainer: Yuwen Gu <yuwen.gu@uconn.edu>
#'
#' @seealso \code{\link{cpernet}}, \code{\link{coef.cpernet}},
#'   \code{\link{plot.cpernet}}, \code{\link{print.cpernet}}
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
#' tau <- 0.30
#' pf <- abs(rnorm(p))
#' pf2 <- abs(rnorm(p))
#' w <- 2.0
#' lambda2 <- 1
#' m2 <- cpernet(y = y, x = x, w = w, tau = tau, eps = 1e-8,
#'               pf.mean = pf, pf.scale = pf2,
#'               standardize = FALSE, lambda2 = lambda2)
#' predict(m2, newx = x, s = m2$lambda[50])
#'
#' @export
predict.cpernet <- function(object, newx, s = NULL, type = "response",
                            ...) NextMethod("predict")
