#' Optimization Control Parameters Passed to \code{optim}
#'
#' Optimization parameters passed to \code{\link[stats]{optim}} for the fit of regression model to
#'     count data via \code{\link{mvreg}}. This function acts in the same spirit as
#'     \code{\link[betareg]{betareg.control}} from the \code{betareg} package. Its primary purpose
#'     is to gather all the optimization control arguments in a single function.
#'
#' @param method the method to be used. See 'Details' in \code{\link[stats]{optim}}. The default
#'     method (\code{"BFGS"}) is a quasi-Newton method (also known as a variable metric algorithm),
#'     specifically that published simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno.
#' @param maxit the maximum number of iterations of the algorithm. Defaults to \code{2000}.
#' @param hessian logical. Should a numerically differentiated Hessian matrix be returned?
#' @param start an optional vector with starting values for all parameters. It must be passed in the
#'     order: \code{(beta, gamma)}, where
#'     \code{beta} and \code{gamma} are the coefficients with the mean and variance submodels,
#'     respectively.
#' @param ... further arguments to be passed to \code{\link[stats]{optim}}.
#'
#' @references Cribari-Neto, F., and Zeileis, A. (2010). Beta regression in R.
#'     \emph{Journal of statistical software}, 34, 1-24.
#'
#'
#' @references Kokonendji, C. C., Medeiros, R. M. R., and Bourguignon, M. (2024). Mean and variance
#'     count regression models based on reparameterized distributions. Submitted.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @return A list with the arguments specified.
#' @export
control_fit <- function(method = "BFGS", maxit = 2000, hessian = FALSE, start = NULL, ...)
{
  rval <- list(method = method, maxit = maxit, hessian = hessian, start = start)

  rval <- c(rval, list(...))

  if (!is.null(rval$fnscale))
    warning("fnscale must not be modified\n")

  rval$fnscale <- -1


  rval
}

# Estimation
mle <- function(y, X = NULL, Z = NULL, count = ripo,
                control = control_fit(...), ...){

  ### Control fit list
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian
  trace <- control$trace
  start  <- control$start

  control$start <- control$no_grad <- NULL

  # Count distribution
  family <- as.count(count)

  if (is.null(start))
    start <- family$start(y, X, Z)

  ### Data setting
  n <- length(y)

  if (is.null(X)) X <- matrix(1, nrow = n)
  if (family$npar > 1 & is.null(Z)) Z <- matrix(1, nrow = n)

  p <- ncol(X); k <- ncol(Z)

  # Log-likelihood
  llaux <- function(par) ll(par, y, X, Z, count)

  # Score function
  if (family$has_dl) {
    Uaux <- function(par) Uscore(par, y, X, Z, count)
  }  else {
    Uaux <- NULL
  }

  opt <- suppressWarnings(stats::optim(par = start,
                                       fn = llaux,
                                       gr = Uaux,
                                       method = method,
                                       control = control,
                                       hessian = hessian))

  # Convergence status
  if (opt$convergence > 0)
    warning(cat("optimization failed to converge\n"))

  opt$logLik <- opt$value
  opt$start <- start

  opt

}
