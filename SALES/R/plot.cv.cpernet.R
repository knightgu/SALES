##############################################################################
## This function is adapted/modified based on the plot function
## from the `glmnet` package:
##
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate Descent.
## Journal of Statistical Software, 33(1), 1-22.
## URL: http://www.jstatsoft.org/v33/i01/.
##############################################################################

#' Plot the cross-validated curve produced by cv.cpernet
#'
#' Plots the cross-validated curve, and upper and lower standard deviation
#' curves, as a function of the \code{lambda} values used. This function is
#' modified based on the \code{plot.cv.glmnet} function from the \code{glmnet}
#' package.
#'
#' A plot is produced.
#'
#' @param x fitted \code{\link{cv.cpernet}} object
#' @param sign.lambda either plot against \code{log(lambda)} (default) or its
#'   negative if \code{sign.lambda=-1}.
#' @param \dots other graphical parameters to plot
#'
#' @author Yuwen Gu and Hui Zou\cr
#'
#'   Maintainer: Yuwen Gu <yuwen.gu@uconn.edu>
#'
#' @seealso \code{\link{plot.cpernet}}
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
#' m2.cv <- cv.cpernet(y = y, x = x, w = w, tau = tau, eps = 1e-8,
#'                     pf.mean = pf, pf.scale = pf2,
#'                     standardize = FALSE, lambda2 = lambda2)
#' plot(m2.cv)
#'
#' @export
plot.cv.cpernet <- function(x, sign.lambda = 1, ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  if (sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  plot.args <- list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
                    ylim = range(cvobj$cvupper, cvobj$cvlower), xlab = xlab,
                    ylab = cvobj$name, type = "n")
  new.args <- list(...)
  if (length(new.args))
    plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvupper,
             cvobj$cvlower, width = 0.01, col = "darkgrey")
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20, col = "red")
  ## axis(side = 3, at = sign.lambda * log(cvobj$lambda),
  ##      labels = paste(cvobj$nz),
  ## tick = FALSE, line = 0)
  abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
  abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
  invisible()
}
