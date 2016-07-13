###############################################################
## This function is adapted/modified based on the plot
#    function from the glmnet package:
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via
#    Coordinate Descent.
## Journal of Statistical Software, 33(1), 1-22.
## URL http://www.jstatsoft.org/v33/i01/.
###############################################################



#' Plot coefficients from an ernet object
#'
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' ernet object. This function is modified based on the \code{plot} method in
#' the \code{glmnet} package.
#'
#' A coefficient profile plot is produced.
#'
#' @param x fitted \code{\link{ernet}} model
#' @param xvar what is on the x-axis. \code{"norm"} plots against the L1-norm
#' of the coefficients, \code{"lambda"} against the log-lambda sequence.
#' @param color if \code{TRUE}, plot the curves with rainbow colors. Otherwise,
#' plot the curves with gray colors. Default is \code{FALSE}.
#' @param label if \code{TRUE}, label the curves with variable sequence
#' numbers. Otherwise, do not put labels. Default is \code{FALSE}.
#' @param \dots other graphical parameters to plot.
#' @author Yuwen Gu and Hui Zou\cr Maintainer: Yuwen Gu <guxxx192@@umn.edu>
#' @seealso \code{\link{plot.cv.ernet}}
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
#' m1 <- ernet(y = y, x = x, tau = tau, eps = 1e-8, pf = pf,
#'             pf2 = pf2, standardize = FALSE, intercept = FALSE,
#'             lambda2 = lambda2)
#' plot(m1)
#'
#' @importFrom graphics matplot axis text
#' @importFrom grDevices gray.colors
#' @export
plot.ernet <- function(x, xvar = c("norm", "lambda"),
  color = FALSE, label = FALSE, ...) {
  beta <- x$beta
  lambda <- x$lambda
  df <- x$df
  if (all(df == 0)) stop("All coefficients are zero")
  xvar <- match.arg(xvar)
  ##beta should be in 'dgCMatrix' format
  which <- nonzero(beta)
  beta <- as.matrix(beta[which, , drop = FALSE])
  switch(xvar, norm = {
      index <- apply(abs(beta), 2, sum)
      iname <- "L1 Norm"
    }, lambda = {
         index <- log(lambda)
         iname <- "Log Lambda"})
  xlab <- iname
  ylab <- "Coefficients"
  dotlist <- list(...)
  type <- dotlist$type
  if (is.null(type)) {
    if (color == FALSE)
      matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
        type = "l", pch = 500, col = gray.colors(12, start = 0.05,
        end = 0.7, gamma = 2.2), ...) else matplot(index, t(beta),
        lty = 1, xlab = xlab, ylab = ylab, type = "l", pch = 500, ...)
  } else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, ...)
  atdf <- pretty(index)
  prettydf <- trunc(approx(x = index, y = df, xout = atdf, rule = 2)$y)
  axis(3, at = atdf, labels = prettydf, cex.axis = 1, tcl = NA)
  if (label) {
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    if (xvar == "lambda") {
      xpos <- min(index)
      pos <- 2
    }
    xpos <- rep(xpos, nnz)
    ypos <- beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
}
