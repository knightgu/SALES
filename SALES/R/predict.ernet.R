#' Make predictions from an ernet object
#'
#' Similar to other predict methods, this functions predicts fitted values from
#' a fitted ernet object.
#'
#' @aliases predict.ernet predict.alspath
#'
#' @param object fitted \code{\link{ernet}} model object.
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
#' @seealso \code{\link{ernet}}, \code{\link{coef.ernet}},
#'   \code{\link{plot.ernet}}, \code{\link{print.ernet}}
#'
#' @keywords models regression
#'
#' @export
predict.ernet <- function(object, newx, s = NULL, type = "response",
                          ...) NextMethod("predict")


##' Model predictions
##'
##' \code{predict} is a generic function for predictions from the results of
##' various model fitting functions. The function invokes particular
##' \emph{methods} which depend on the \code{\link{class}} of the first
##' argument.
##'
##' @param object a model object for which prediction is desired.
##' @param \dots additional arguments affecting the predictions produced.
##'
##' @return The form of the value returned by \code{predict} depends on the
##'   class of its argument. See the documentation of the particular methods for
##'   details of what is produced by that method.
##'
##' @seealso \code{\link{predict.ernet}}, \code{\link{predict.cpernet}}.
##'
##' @export
predict <- function(object,...) UseMethod("predict")
