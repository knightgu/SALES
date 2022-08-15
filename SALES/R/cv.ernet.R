#' Cross-validation for ernet
#'
#' Does k-fold cross-validation for ernet, produces a plot, and returns a value
#' for \code{lambda}. This function is based on the \code{cv} function from the
#' \code{glmnet} package.
#'
#' @param x \code{x} matrix as in \code{\link{ernet}}.
#' @param y response variable \code{y} as in \code{\link{ernet}}.
#' @param lambda optional user-supplied lambda sequence; default is \code{NULL},
#'   and \code{\link{ernet}} chooses its own sequence.
#' @param nfolds number of folds. Default value is 5. Although \code{nfolds} can
#'   be as large as the sample size (leave-one-out CV), it is not recommended
#'   for large datasets. Smallest value allowed is 3.
#' @param foldid an optional vector of values between 1 and \code{nfolds},
#'   identifying what fold each observation is in. If supplied, \code{nfolds}
#'   will be supressed.
#' @param pred.loss loss function used to calculate cross-validation error. The
#'   only option now is \code{"loss"}, which is the asymmetric squared error
#'   loss (ASEL).
#' @param tau the asymmetry coefficient \eqn{\tau} used in the asymmetric
#'   squared error loss.
#' @param \dots other arguments that can be passed to ernet.
#'
#' @details The function runs \code{\link{ernet}} \code{nfolds}+1 times; the
#'   first to get the \code{lambda} sequence, and the remainder to compute the
#'   fit with each of the folds removed. The average error and standard
#'   deviation over the folds are computed.
#'
#' @return an object of class \code{\link{cv.ernet}} is returned, which is a
#'   list with the ingredients of the cross-validation fit.
#'
#'   \item{lambda}{the values of \code{lambda} used in the fits.}
#'
#'   \item{cvm}{the mean cross-validated error - a vector of length
#'   \code{length(lambda)}.}
#'
#'   \item{cvsd}{estimate of standard error of \code{cvm}.}
#'
#'   \item{cvupper}{upper curve = \code{cvm+cvsd}.}
#'
#'   \item{cvlower}{lower curve = \code{cvm-cvsd}.}
#'
#'   \item{nzero}{number of non-zero coefficients at each \code{lambda}.}
#'
#'   \item{name}{a text string indicating type of measure (for plotting
#'   purposes).}
#'
#'   \item{ernet.fit}{a fitted \code{\link{ernet}} object for the full data.}
#'
#'   \item{lambda.min}{The optimal value of \code{lambda} that gives minimum
#'   cross validation error \code{cvm}.}
#'
#'   \item{lambda.1se}{The largest value of \code{lambda} such that error is
#'   within 1 standard error of the minimum.}
#'
#' @author Yuwen Gu and Hui Zou\cr
#'
#'   Maintainer: Yuwen Gu <yuwen.gu@uconn.edu>
#'
#' @seealso \code{\link{ernet}}
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
#'
#' @export
cv.ernet <- function(x, y, lambda = NULL, pred.loss = "loss", nfolds = 5,
                     foldid, tau = 0.5, ...) {
  pred.loss <- match.arg(pred.loss)
  N <- nrow(x)
  ## Fit the model once to get dimensions etc of output
  y <- drop(y)
  ernet.object <- ernet(x, y, lambda = lambda, tau = tau, ...)
  lambda <- ernet.object$lambda
  ## predict -> coef
  nz <- sapply(coef(ernet.object, type = "nonzero"), length)
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N))
  } else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be at least 3; nfolds=10 recommended")
  outlist <- vector("list", length = nfolds)
  ## Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    whichfold <- foldid == i
    y_sub <- y[!whichfold]
    outlist[[i]] <- ernet(x = x[!whichfold, , drop = FALSE], y = y_sub,
                          lambda = lambda, tau = tau, ...)
  }
  ## What to do depends on the pred.loss and the model fit
  fun <- paste("cv", class(ernet.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, pred.loss, tau))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd,
              cvlower = cvm - cvsd, nzero = nz, name = cvname,
              ernet.fit = ernet.object)
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.ernet"
  obj
}
