#' Make predictions from an ernet object
#'
#' Similar to other predict methods, this functions predicts fitted values from
#' a fitted ernet object.
#'
#' \code{s} is the new vector at which predictions are to be made. If \code{s}
#' is not in the lambda sequence used for fitting the model, the \code{predict}
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of predicted values from both left and
#' right \code{lambda} indices.
#'
#' @aliases predict.ernet predict.alspath
#' @param object fitted \code{\link{ernet}} model object.
#' @param newx matrix of new values for \code{x} at which predictions are to be
#' made. NOTE: \code{newx} must be a matrix, \code{predict} function does not
#' accept a vector or other formats of \code{newx}.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#' predictions are to be made. Default is the entire sequence used to create
#' the model.
#' @param type type of prediction required. Only \code{response} is available.
#' Gives predicted response for regression problems.
#' @param \dots Not used. Other arguments to predict.
#' @return The object returned depends on type.
#' @author Yuwen Gu and Hui Zou\cr Maintainer: Yuwen Gu <guxxx192@@umn.edu>
#' @seealso \code{\link{ernet}}, \code{\link{coef.ernet}}
#' @keywords models regression
#'
#' @export
predict.ernet <- function(object, newx, s = NULL,
    type = "response", ...) NextMethod("predict")
