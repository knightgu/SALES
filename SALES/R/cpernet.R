#' Regularization paths for the coupled sparse asymmetric least squares
#' (COSALES) regression (or the coupled sparse expectile regression)
#'
#' Fits regularization paths for coupled sparse asymmetric least squares
#' regression at a sequence of regularization parameters.
#'
#' @param x matrix of predictors, of dimension (nobs * nvars); each row is an
#'   observation.
#' @param y response variable.
#' @param w weight applied to the asymmetric squared error loss of the mean
#'   part. See details. Default is 1.0.
#' @param nlambda the number of \code{lambda} values (default is 100).
#' @param method a character string specifying the loss function to use. Only
#'   \code{cper} is available now.
#' @param lambda.factor The factor for getting the minimal lambda in the
#'   \code{lambda} sequence, where we set \code{min(lambda)} =
#'   \code{lambda.factor} * \code{max(lambda)} with \code{max(lambda)} being the
#'   smallest value of \code{lambda} that penalizes all coefficients to zero.
#'   The default value depends on the relationship between \eqn{N} (the number
#'   of observations) and \eqn{p} (the number of predictors). If \eqn{N < p},
#'   the default is \code{0.01}. If \eqn{N > p}, the default is \code{0.0001},
#'   closer to zero. A very small value of \code{lambda.factor} will lead to a
#'   saturated fit. The argument takes no effect if there is a user-supplied
#'   \code{lambda} sequence.
#' @param lambda a user-supplied \code{lambda} sequence. Typically, by leaving
#'   this option unspecified users can have the program compute its own
#'   \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor}. It
#'   is better to supply, if necessary, a decreasing sequence of \code{lambda}
#'   values than a single (small) value. The program will ensure that the
#'   user-supplied \code{lambda} sequence is sorted in decreasing order.
#' @param lambda2 regularization parameter \code{lambda2} for the quadratic
#'   penalty of the coefficients. Default is 0, meaning no L2 penalization.
#' @param pf.mean,pf.scale L1 penalty factor of length \eqn{p} used for adaptive
#'   LASSO or adaptive elastic net. Separate L1 penalty weights can be applied
#'   to each mean or scale coefficient to allow different L1 shrinkage. Can be 0
#'   for some variables, which imposes no shrinkage and results in that variable
#'   being always included in the model. Default is 1 for all variables (and
#'   implicitly infinity for variables listed in \code{exclude}).
#' @param pf2.mean,pf2.scale L2 penalty factor of length \eqn{p}{p} used for
#'   adaptive elastic net. Separate L2 penalty weights can be applied to each
#'   mean or scale coefficient to allow different L2 shrinkage. Can be 0 for
#'   some variables, which imposes no shrinkage. Default is 1 for all variables.
#' @param exclude indices of variables to be excluded from the model. Default is
#'   none. Equivalent to an infinite penalty factor.
#' @param dfmax limit the maximum number of variables in the model. Useful for
#'   very large \eqn{p}, if a partial path is desired. Default is \eqn{p+1}.
#' @param pmax limit the maximum number of variables ever to be nonzero. For
#'   example once \eqn{\beta} enters the model, no matter how many times it
#'   exits or re-enters the model through the path, it will be counted only
#'   once. Default is \code{min(dfmax*1.2, p)}.
#' @param standardize logical flag for variable standardization, prior to
#'   fitting the model sequence. The coefficients are always returned to the
#'   original scale. Default is \code{TRUE}.
#' @param intercept Should intercept(s) be fitted (default=TRUE) or set to zero
#'   (FALSE).
#' @param eps convergence threshold for coordinate descent. Each inner
#'   coordinate descent loop continues until the maximum change in any
#'   coefficient is less than \code{eps}. Defaults value is \code{1e-8}.
#' @param maxit maximum number of outer-loop iterations allowed at fixed lambda
#'   values. Default is 1e7. If the algorithm does not converge, consider
#'   increasing \code{maxit}.
#' @param tau the parameter \code{tau} in the coupled ALS regression model. The
#'   value must be in (0,1) and cannot be 0.5. Default is 0.8.
#'
#' @details Note that the objective function in \code{cpernet} is
#'   \deqn{w*1'\Psi(y-X\beta,0.5)/N + 1'\Psi(y-X\beta-X\theta,\tau)/N +
#'   \lambda_1*\Vert\beta\Vert_1 + 0.5\lambda_2\Vert\beta\Vert_2^2 +
#'   \mu_1*\Vert\theta\Vert +
#'   0.5\mu_2\Vert\theta\Vert_2^2,}{w*1'\Psi(y-X\beta,0.5)/N +
#'   1'\Psi(y-X\beta-X\theta,\tau)/N + \lambda1*|\beta| +
#'   0.5*\lambda2*||\beta||^2 + \mu1*|\theta| + 0.5*\mu2*||\theta||^2,} where
#'   \eqn{\Psi(u,\tau)=|\tau-I(u<0)|*u^2} denotes the asymmetric squared error
#'   loss and the penalty is a combination of L1 and L2 terms for both the mean
#'   and scale coefficients.
#'
#'   For faster computation, if the algorithm is not converging or running slow,
#'   consider increasing \code{eps}, decreasing \code{nlambda}, or increasing
#'   \code{lambda.factor} before increasing \code{maxit}.
#'
#' @return An object with S3 class \code{\link{cpernet}}. \item{call}{the call
#'   that produced this object.} \item{b0, t0}{intercept sequences both of
#'   length \code{length(lambda)} for the mean and scale respectively.}
#'   \item{beta, theta}{\code{p*length(lambda)} matrices of coefficients for the
#'   mean and scale respectively, stored as sparse matrices (\code{dgCMatrix}
#'   class, the standard class for sparse numeric matrices in the \code{Matrix}
#'   package). To convert them into normal R matrices, use \code{as.matrix()}.}
#'   \item{lambda}{the actual sequence of \code{lambda} values used}
#'   \item{df.beta, df.theta}{the number of nonzero mean and scale coefficients
#'   respectively for each value of \code{lambda}.} \item{dim}{dimensions of
#'   coefficient matrices.} \item{npasses}{total number of iterations summed
#'   over all lambda values.} \item{jerr}{error flag, for warnings and errors, 0
#'   if no error.}
#'
#' @author Yuwen Gu and Hui Zou\cr
#'
#'   Maintainer: Yuwen Gu <yuwen.gu@uconn.edu>
#'
#' @seealso \code{\link{plot.cpernet}}, \code{\link{coef.cpernet}},
#'   \code{\link{predict.cpernet}}, \code{\link{print.cpernet}}
#'
#' @references Gu, Y., and Zou, H. (2016).
#'   "High-dimensional generalizations of asymmetric least squares regression and their applications."
#'   \emph{The Annals of Statistics}, 44(6), 2661–2694.\cr
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
#'
#' @export
cpernet <- function(x, y, w = 1.0, nlambda = 100L, method = "cper",
                    lambda.factor = ifelse(2 * nobs < nvars, 1e-02, 1e-04),
                    lambda = NULL, lambda2 = 0, pf.mean = rep(1, nvars),
                    pf2.mean = rep(1, nvars), pf.scale = rep(1, nvars),
                    pf2.scale = rep(1, nvars), exclude, dfmax = nvars + 1,
                    pmax = min(dfmax * 1.2, nvars), standardize = TRUE,
                    intercept = TRUE, eps = 1e-08, maxit = 1000000L, tau = 0.80) {
  ##=======================================================================##
  ## DATA SETUP
  method <- match.arg(method)
  this.call <- match.call()
  y <- drop(y)
  x <- as.matrix(x)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")
  if (NROW(y) != nobs) stop("x and y have different number of observations")
  if (NCOL(y) > 1L) stop("Multivariate response is not supported now")
  ##=======================================================================##
  ## PARAMETER SETUP
  if (length(pf.mean) != nvars)
    stop("Size of L1 penalty factors for the mean does not match the number of input variables")
  if (length(pf2.mean) != nvars)
    stop("Size of L2 penalty factors for the mean does not match the number of input variables")
  if (length(pf.scale) != nvars)
    stop("Size of L1 penalty factors for the scale does not match the number of input variables")
  if (length(pf2.scale) != nvars)
    stop("Size of L2 penalty factors for the scale does not match the number of input variables")
  if (lambda2 < 0) {
    warning("lambda2 < 0; set to zero...")
    lambda2 <- 0
  }
  maxit <- as.integer(maxit)
  lam2 <- as.double(lambda2)
  pfmean <- as.double(pf.mean)
  pf2mean <- as.double(pf2.mean)
  pfscale <- as.double(pf.scale)
  pf2scale <- as.double(pf2.scale)
  isd <- as.integer(standardize)
  intr <- as.integer(intercept)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) stop("Some excluded variables out of range")
    jd <- as.integer(c(length(jd), jd))
  } else jd <- as.integer(0)
  ##=======================================================================##
  ## LAMBDA SETUP
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    flmin <- as.double(1) # flmin = 1 if user defines lambda
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  ##=======================================================================##
  fit <- cpalspath(x, y, w, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd,
                   pfmean, pf2mean, pfscale, pf2scale, maxit, lam2, tau, nobs,
                   nvars, vnames)
  if (is.null(lambda)) fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  ##=======================================================================##
  class(fit) <- c("cpernet", class(fit))
  fit
}
