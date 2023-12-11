# The Reparameterized Generalized Poisson Distribution -------------------------------------

#' @name gpo
#'
#' @title The Reparameterized Generalized Poisson Distribution
#'
#' @description Probability mass function, distribution function, quantile function, and random generation
#' for the reparameterized generalized Poisson distribution with mean \code{mu} and
#' variance \code{sigma2}.
#'
#' @param mu.link,sigma2.link defines the link functions for the mean and variance regression models,
#'     respectively, in the \code{\link{mvreg}} function, with \code{"log"} link as the default.
#'     Other links are \code{"identity"} and \code{"sqrt"}.
#' @param x vector of non-negative integer quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu numeric; non-negative mean.
#' @param sigma2 numeric; variance (greather than \code{abs(mu * (mu - 1))}),
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P(X <= x)}, otherwise, \code{P(X > x)}.
#' @param log.p 	logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... possible extra parameters, in addition to the mean and variance, in a distribution
#'    added to the \code{rigpsreg} package. Not used in generalized Poisson distribution.
#'
#' @return \code{dgpo} returns the probability mass function, \code{pgpo} gives the distribution
#'     function, \code{qgpo} gives the quantile function, and \code{rgpo} generates random
#'     observations.
#'
#'     \code{dgo} returns a \code{"count"} object mainly used to fit the generalized Poisson distribution
#'         via \code{\link{mvreg}} function. It returns a list with the following components:
#'         \describe{
#'           \item{abb}{abbreviation of the distribution in the \code{mvcount} package.}
#'           \item{name}{capitalized name of the distribution.}
#'           \item{npar}{number of parameters of the distribution.}
#'           \item{mu.link, sigma2.link}{specified link functions for the mean (\code{mu}) and the
#'               variance (\code{sigma2}) parameters, respectively.}
#'           \item{constraint}{a function \code{function(mu, sigma2)} that returns \code{TRUE} if
#'               the values of \code{mu} and \code{sigma} satisfy the constraint between the mean
#'               and variance of the distribution.}
#'           \item{start}{a function \code{function(x, X, Z)} which provides initial values for the
#'               regression coefficients of a double regression fit of the distribution.}
#'           \item{has_dl}{logical; if \code{TRUE}, it informs that the model has first and second
#'               order derivatives implemented; otherwise, they are obtained numerically via \code{\link{optim}}.}
#'           \item{dl.dmu,dl.dsigma}{functions \code{f(x, mu, sigma)} that provide the partial
#'               derivatives of the regression coefficients associated with the mean and variance,
#'               respectively. Its use is optional in the mvcount package. If they are not
#'               implemented, they must return \code{NULL}.}
#'          \item{dl2.dmu2, dl2,dsigma2, dl2.dmudsigma}{functions \code{f(x, mu, sigma)} that provide
#'               the second order partial derivatives of the regression coefficients associated with
#'               the mean and variance. Its use is optional in the mvcount package. If they are not
#'               implemented, they must return \code{NULL}.}
#'         }
#'
#' @details This set of functions represents the probability mass function, the cumulative
#'    distribution function, quantile function, and a random number generator for the generalized Poisson
#'    distribution parameterized in terms of the mean (\code{mu}) and the variance (\code{sigma2}).
#'    This distribution was applied by Kokonendji, Medeiros, and Bourguignon (2024).
#'
#' @references
#'
#' Consul, P. C., and Jain, G. C. (1973). A generalization of the Poisson distribution.
#'     \emph{Technometrics}, \bold{15}, 791-799.
#'
#' Kokonendji, C. C., Medeiros, R. M. R., and Bourguignon, M. (2024). Mean and variance
#'     count regression models based on reparameterized distributions. Submitted.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
#' @examples
#' n <- 10000
#'
#' mu <- 2
#' sigma2 <- 5
#'
#' y <- rgpo(n, mu, sigma2)
#'
#' # mean and variance theoretical vs empirical comparison
#' cat("\n mu:", mu, "--- sample mean:", round(mean(y), 3),
#'     "\n sigma2:", sigma2, "--- sample variance:", round(var(y), 3), "\n")
#'
#' # probability mass function theoretical vs empirical comparison
#' yvals <- sort(unique(y))
#' xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf", col = "white")
#' points(xcoords, dgpo(yvals, mu, sigma2), pch = 16, col = 2, type = "b")
#'
#' # distribution function theoretical vs empirical comparison
#' plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
#' curve(pgpo(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
#'
#' # quantile function theoretical vs empirical comparison
#' plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#'      pch = 16, xlab = "p", ylab = "Quantile")
#' curve(qgpo(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)

## Probability mass function
#' @rdname gpo
#' @export
dgpo <- function(x, mu, sigma2, log.p = FALSE, ...){

  phi <- (sqrt(sigma2) - sqrt(mu)) / (mu^(3/2))

  if (is.vector(x))
    x <- matrix(x, nrow = length(x))

  n <- dim(x)[1]
  mu <- matrix(mu, nrow = n)
  phi <- matrix(phi, nrow = n)

  pmf <- matrix(-Inf, nrow = n)

  # NaN index
  pmf[which(mu <= 0 | phi <= 0, arr.ind = TRUE)] <- NaN

  # Positive pmf index
  id <- which(!is.nan(pmf), arr.ind = TRUE)

  if (length(id) != 0)
    pmf[id] <- gamlss.dist::dGPO(x = x[id], mu = mu[id], sigma = phi[id], log = TRUE)

  if(!log.p) pmf <- exp(pmf)

  as.vector(pmf)


}

## Cumulative distribution function (inerited from gamlss)
#' @rdname gpo
#' @export
pgpo <- function(q, mu, sigma2, lower.tail = TRUE, log.p = FALSE, ...)
{

  phi <- (sqrt(sigma2) - sqrt(mu)) / (mu^(3/2))

  if (is.vector(q))
    q <- matrix(q, nrow = length(q))

  n <- dim(q)[1]

  mu <- matrix(mu, nrow = n)
  phi <- matrix(phi, nrow = n)

  cdf <- matrix(0, nrow = n)

  # NaN index
  cdf[which(mu <= 0 | phi <= 0, arr.ind = TRUE)] <- NaN

  # Positive pmf index
  id <- which(q >= 0 & !is.nan(cdf), arr.ind = TRUE)

  if (length(id) != 0)
    cdf[id] <- gamlss.dist::pGPO(q = q[id], mu = mu[id], sigma = phi[id])

  if(!lower.tail) cdf <- 1 - cdf
  if(log.p) cdf <- log(cdf)

  as.numeric(cdf)

}

## Quantile function (inerited from gamlss)
#' @rdname gpo
#' @export
qgpo <- function(p, mu, sigma2, lower.tail = TRUE, ...)
{

  phi <- (sqrt(sigma2) - sqrt(mu)) / (mu^(3/2))

  if(!lower.tail)
    p <- 1 - p

  if (is.vector(p))
    p <- matrix(p, nrow = length(p))

  n <- dim(p)[1]

  mu <- matrix(mu, nrow = n)
  phi <- matrix(phi, nrow = n)

  qtl <- matrix(0, nrow = n)

  # NaN index
  qtl[which(mu <= 0 | sigma2 <= 0, arr.ind = TRUE)] <- NaN

  id <- which(!is.nan(qtl), arr.ind = TRUE)
  if (length(id) != 0)
    qtl[id] <- gamlss.dist::qGPO(p[id], mu = mu[id], sigma = phi[id])

  as.numeric(qtl)
}

## Random generation
#' @rdname gpo
#' @export
rgpo <- function(n, mu, sigma2, ...){

  phi <- (sqrt(sigma2) - sqrt(mu)) / (mu^(3/2))

  gamlss.dist::rGPO(n, mu = mu, sigma = phi)

}

#' @rdname gpo
#' @export
gpo <- function(mu.link = "log", sigma2.link = "log"){

  out <- list()

  # Abbreviation
  out$abb <- "gpo"

  # Name
  out$name <- "Generalized Poisson"

  # Number of parameters
  out$npar <- 2

  # Link functions
  out$mu.link <- mu.link
  out$sigma2.link <- sigma2.link

  # Constraint -----------------------------------------------------------------
  out$constraint <- function(mu, sigma2, ...) (mu > 0) & (sigma2 > mu)

  # Initial values -------------------------------------------------------------
  out$start <- function(x, X, Z){

    g1 <- make.link.ldrm(mu.link)$fun
    g1.inv <- make.link.ldrm(mu.link)$inv
    g2 <- make.link.ldrm(sigma2.link)$fun
    g2.inv <- make.link.ldrm(sigma2.link)$inv

    k <- NCOL(Z)

    # Initial values
    beta0 <- solve(t(X)%*%X)%*%t(X)%*%g1(x + 0.1)
    mu0 <- g1.inv(X%*%beta0)

    sigma0 <- stats::var(x) + max(mu0)
    gamma0 <- c(g2(sigma0), rep(0, k - 1))
    c(beta0, gamma0)

  }

  # Derivatives ----------------------------------------------------------------
  out$has_dl <- FALSE

  out$dl.dmu <- function(x, mu, sigma, ...){NULL}
  out$dl2.dmu2 <- function(x, mu, sigma, ...){NULL}
  out$dl.dsigma <- function(x, mu, sigma, ...){NULL}
  out$dl2.dsigma2 <- function(x, mu, sigma, ...){NULL}
  out$dl2.dmudsigma <- function(x, mu, sigma, ...){NULL}

  structure(out, class = "count")

}


