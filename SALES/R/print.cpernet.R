#' Print a cpernet object
#'
#' Print a summary of the \code{cpernet} path at each step along the path.
#'
#' The call that produced the \code{\link{cpernet}} object is printed, followed
#' by a three-column matrix with columns \code{Df1}, \code{Df2} and
#' \code{Lambda}. The \code{Df1} and \code{Df2} columns are the number of
#' nonzero mean and scale coefficients respectively.
#'
#' @param x fitted \code{\link{cpernet}} object.
#' @param digits significant digits in the output.
#' @param \dots additional print arguments.
#' @return a three-column matrix, the first two columns are the number of
#' nonzero mean and scale coefficients respectively and the third column is
#' \code{Lambda}.
#' @author Yuwen Gu and Hui Zou\cr Maintainer: Yuwen Gu <guxxx192@@umn.edu>
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
#' m2 <- cpernet(y = y, x = x, w = w, tau = tau, eps = 1e-8,
#'               pf.mean = pf, pf.scale = pf2,
#'               standardize = FALSE, lambda2 = lambda2)
#' print(m2)
#'
#' @method print cpernet
#' @export
print.cpernet <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Df1 = x$df.beta, Df2 = x$df.theta,
                Lambda = signif(x$lambda, digits)))
}
