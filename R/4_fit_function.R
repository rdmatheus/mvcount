#' Mean and Variance Regression Models for Count Data
#'
#' Fit a regression model to count data via maximum likelihood for a distributions indexed by
#'      mean and variance parameters
#'
#' @param formula a simbolic description of the model, of type
#'     \code{y ~ x} for covariates in the mean submodel only or \code{y ~ x | z}
#'     to enter covariates in the variance submodel. See details below.
#' @param data an optional data frame containing the variables in the formula.
#'     By default the variables are taken from environment(formula).
#' @param count a \code{"\link{count}"} object, which is used to define the distribution and the link
#'     functions of the mean and variance parameters. The current distribution families available
#'     in the \code{mvcount} package can be found in \code{\link{count}}.
#' @param y,x logical; if \code{TRUE} the corresponding components of the fit, response and model
#'     matrices, are returned.
#' @param control a list of control arguments specified via \code{\link{control_fit}}.
#'     (under development).
#' @param ... further arguments passed to \code{\link{control_fit}}.
#'
#' @return  The \code{mvreg} function returns an object of class \code{"mvreg"},
#'  which consists of a list with the following components:
#' \describe{
#'   \item{coefficients}{a list containing the elements \code{"mean"} and \code{"variance"}.
#'       which consists of the estimates of the coefficients associated with the
#'       mean and the variance, respectively.}
#'   \item{fitted.values}{a vector with the fitted means.}
#'   \item{sigma2}{a vector with the fitted variances.}
#'   \item{residuals}{a vector of raw residuals \code{(observed - fitted)}.}
#'   \item{mu.link, sigma2.link}{the specified link functions for the respective submodels of the
#'       mean and variance.}
#'   \item{count}{the count distribution.}
#'   \item{vcov}{asymptotic covariance matrix of the maximum likelihood
#'       estimator of the model parameters vector. Specifically, the inverse of
#'       the Fisher information matrix.}
#'   \item{logLik}{log-likelihood of the fitted model.}
#'   \item{freq}{expected frequencies after fitting the count regression.}
#'   \item{nobs}{the number of observations in the sample.}
#'   \item{df.null}{residual degrees of freedom in the null model
#'         (constant mean and variance), that is, \code{n - 2}.}
#'   \item{df.residual}{residual degrees of freedom in the fitted model.}
#'   \item{feasible}{logical; if \code{TRUE}, the estimates satisfy the constraints of the
#'       mean and variance parameters.}
#'   \item{upsilon}{a goodness-of-fit statistics, where the closer to zero the better fitted.}
#'   \item{optim.pars}{control optimization parameters used by \code{\link{control_fit}}.}
#'   \item{call}{the function call.}
#'   \item{formula}{the formula used to specify the model in \code{bergreg}.}
#'   \item{terms}{ a list with elements "mean", "variance" and "full"
#'          containing the terms objects for the respective models.}
#'   \item{y}{ the response vector (if \code{y = TRUE}).}
#'   \item{x}{ a list with elements "mean" and "variance" containing the
#'       model matrices from the respective models (if \code{X = TRUE}).}
#'  }
#'
#' @details The basic formula is of type \code{y ~ x1 + x2 + ... + xp} which
#'     specifies the model for the mean response only with \code{p} explanatory
#'     variables. Following the syntax of the \code{betareg} package
#'     (Cribari-Neto and Zeileis, 2010), the model for the variance, say in terms of
#'     z1, z2, ..., zk, is specified as \code{y ~ x1 + x2 + ... + xp | z1 + z2 +
#'     ... +zk} using functionalities inherited from package \code{Formula}
#'     (Zeileis and Croissant, 2010).
#'
#' @references Bourguignon, M., and de Medeiros, R. M. (2022). A simple and useful regression model
#'     for fitting count data. \emph{TEST}, \bold{31}, 790-827.
#'
#' @references Cribari-Neto F, Zeileis A (2010). Beta regression in R. \emph{Journal
#'  of Statistical Software}, \bold{34}, 1–24
#'
#' @references Kokonendji, C. C., Medeiros, R. M. R., and Bourguignon, M. (2024). Mean and variance
#'     count regression models based on reparameterized distributions. Submitted.
#'
#' @references Zeileis A., and Croissant, Y. (2010). Extended model formulas in
#'  R: multiple parts and multiple responses. \emph{Journal of Statistical Software},
#'  \bold{34}, 1–13.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
mvreg <- function(formula, data,
                  count = ripo,
                  y = FALSE, x = FALSE,
                  control = control_fit(...), ...)
{

  responseLogical <- y
  xLogical <- x

  # Count distribution
  family <- as.count(count)

  ## Call
  cl <- match.call()
  if (missing(data))  data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  }else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract terms, model matrix, response
  mt <- stats::terms(formula, data = data)

  y <- stats::model.response(mf, "numeric")
  n <- length(y)

  mtX <- stats::terms(formula, data = data, rhs = 1L)
  X <- stats::model.matrix(mtX, mf)
  p <- NCOL(X)

  Z <- NULL
  k <- 0
  if(family$npar > 1){
    mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2L))
    Z <- stats::model.matrix(mtZ, mf)
    k <- NCOL(Z)
  }

  if (length(y) < 1)
    stop("empty model")

  mu.link <- family$mu.link
  sigma2.link <- family$sigma2.link

  opt <- mle(y, X, Z, count = count, control = control)

  nu <- NULL
  if(family$npar > 1){

    ## Coefficients
    beta <- opt$par[1:p]
    gamma <- opt$par[1:k + p]

    nextra_par <- family$npar - 2
    if(family$npar > 2){
      nu <- opt$par[p + k + nextra_par]
    }

    ## Covariance matrix
    if(control$hessian){
      vcov <- solve(-opt$hessian)
    }else{
      vcov <- solve(-Jhessian(opt$par, y, X, Z, family))
    }

    ## Fitted values
    mu <- as.numeric(make.link.ldrm(mu.link)$inv(X%*%beta))
    sigma2 <- as.numeric(make.link.ldrm(sigma2.link)$inv(Z%*%gamma))

    ## Expected frequencies
    ddist <- get(paste0("d", family$abb))

    expected <- function(y){
      x <- sort(unique(y))
      n <- length(y)
      s <- rep(0, length(x))

      if (length(sigma2) == 1) sigma2 <- as.matrix(rep(sigma2, n))

      for(i in 1:n){
        s <- s + ddist(x, mu = mu[i], sigma2 = sigma2[i], nu = nu)
      }
      return(s)
    }

    freq <- expected(y)

    out <- list(coefficients = list(mean = beta, variance = gamma),
                fitted.values = structure(mu, .Names = names(y)),
                sigma2 = sigma2,
                nu = nu,
                residuals = as.numeric(y - mu),
                mu.link = mu.link,
                sigma2.link = sigma2.link,
                count = family$abb,
                vcov = vcov,
                logLik = opt$logLik,
                freq = freq,
                convergence = opt$convergence,
                message = opt$message,
                start = opt$start,
                control = control,
                feasible = all(family$constraint(mu, sigma2 = sigma2, nu = nu)),
                nobs = n,
                df.null = n - family$npar,
                df.residual = n - p - k - length(nu),
                optim.pars = opt,
                terms = list(mean = mtX, variance = mtZ, full = mt))

    if(responseLogical) out$y <- y
    if(xLogical) out$x <- list(mean = X, variance = Z)

  }else{

    ## Coefficients
    beta <- opt$par[1:p]

    ## Covariance matrix
    vcov <- solve(-Jhessian(opt$par, y, X, Z,  family))

    ## Fitted values
    mu <- as.numeric(make.link.ldrm(mu.link)$inv(X%*%beta))

    ## Expected frequencies
    ddist <- get(paste0("d", family$abb))

    expected <- function(y){
      x <- sort(unique(y))
      n <- length(y)
      s <- rep(0, length(x))

      for(i in 1:n){
        s <- s + ddist(x, mu[i])
      }
      return(s)
    }

    freq <- expected(y)

    out <- list(coefficients = list(mean = beta),
                fitted.values = structure(mu, .Names = names(y)),
                sigma2 = NULL,
                residuals = as.numeric(y - mu),
                mu.link = mu.link,
                sigma2.link = NULL,
                count = family$abb,
                vcov = vcov,
                logLik = opt$logLik,
                freq = freq,
                convergence = opt$convergence,
                message = opt$message,
                start = opt$start,
                control = control,
                feasible = all(family$constraint(mu)),
                nobs = n,
                df.null = n - family$npar,
                df.residual = n - p - k - length(nu),
                optim.pars = opt,
                terms = list(mean = mtX, variance = NULL, full = mtX))

    if(responseLogical) out$y <- y
    if(xLogical) out$x <- list(mean = X, variance = NULL)
  }

  ## Further model information
  out$call <- cl
  out$formula <- formula

  class(out) <- "mvreg"

  # Upsilon statistic
  upsilon <- vector()
  for(i in 1:100){
    aux_res <- stats::residuals(out, type = "quantile")
    upsilon[i] <- mean(abs(sort(aux_res) - EnvStats::evNormOrdStats(n = n)))
  }

  out$upsilon <- mean(upsilon)

  out
}
