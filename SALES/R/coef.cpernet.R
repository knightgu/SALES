#' Get coefficients from a cpernet object
#'
#' Computes the coefficients or returns a list of the indices of the nonzero
#' coefficients at the requested values for \code{lambda} from a fitted cpernet
#' object.
#'
#' \code{s} is the new vector at which predictions are requested. If \code{s}
#' is not in the lambda sequence used for fitting the model, the \code{coef}
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' \code{lambda} indices.
#'
#' @aliases coef.cpernet coef.cpalspath
#' @param object fitted \code{\link{cpernet}} model object.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#' predictions are to be made. Default is the entire sequence used to create
#' the model.
#' @param type type \code{"coefficients"} computes coefficients at the
#' requested values for \code{s}. Type \code{"nonzero"} returns a list of the
#' indices of nonzero coefficients for each value of \code{s}. Default is
#' \code{"coefficients"}.
#' @param \dots not used. Other arguments to predict.
#' @return The object returned depends on type.
#' @author Yuwen Gu and Hui Zou\cr Maintainer: Yuwen Gu <guxxx192@@umn.edu>
#' @seealso \code{\link{cpernet}}, \code{\link{predict.cpernet}}
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
#'               pf.mean = pf, pf.scale = pf2, intercept = TRUE,
#'               standardize = FALSE, lambda2 = lambda2)
#'
#' mean.coef <- as.vector(coef(m2, s = m2$lambda[50])[[1]])
#' scale.coef <- as.vector(coef(m2, s = m2$lambda[50])[[2]])
#'
#' @export
coef.cpernet <- function(object, s = NULL, type = c("coefficients",
    "nonzero"), ...) NextMethod("coef")
