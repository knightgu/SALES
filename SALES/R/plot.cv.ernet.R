###############################################################
## This function is adapted/modified based on the plot.cv
#    function from the glmnet package:
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via
#    Coordinate Descent.
## Journal of Statistical Software, 33(1), 1-22.
## URL http://www.jstatsoft.org/v33/i01/.
###############################################################



#' Plot the cross-validated curve produced by cv.ernet
#'
#' Plots the cross-validated curve, and upper and lower standard deviation
#' curves, as a function of the \code{lambda} values used. This function is
#' modified based on the \code{plot.cv.glmnet} function from the \code{glmnet}
#' package.
#'
#' A plot is produced.
#'
#' @param x fitted \code{\link{cv.ernet}} object
#' @param sign.lambda either plot against \code{log(lambda)} (default) or its
#' negative if \code{sign.lambda=-1}.
#' @param \dots other graphical parameters to plot
#' @author Yuwen Gu and Hui Zou\cr Maintainer: Yuwen Gu <guxxx192@@umn.edu>
#' @seealso \code{\link{cv.ernet}}
#' @keywords models regression
#' @examples
#'
#' set.seed(1)
#' n <- 100
#' p <- 400
#' x <- matrix(rnorm(n*p), n, p)
#' y <- rnorm(n)
#' tau <- 0.90
#' pf <- abs(rnorm(p))
#' pf2 <- abs(rnorm(p))
#' lambda2 <- 1
#' m1.cv <- cv.ernet(y = y, x = x, tau = tau, eps = 1e-8, pf = pf,
#'             pf2 = pf2, standardize = FALSE, intercept = FALSE,
#'             lambda2 = lambda2)
#' plot(m1.cv)
#'
#' @importFrom graphics points abline
#' @export
plot.cv.ernet <- function(x, sign.lambda = 1, ...) {
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
    axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz),
        tick = FALSE, line = 0)
    abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
    abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
    invisible()
}
