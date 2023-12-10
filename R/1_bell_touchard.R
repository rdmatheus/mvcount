# The Reparametrized Bell-Touchard Distribution --------------------------------

#' @name beto
#'
#' @title The Reparameterized Bell-Touchard Distribution
#'
#' @description Probability mass function, distribution function, quantile function, and random generation
#' for the reparameterized Bell-Touchard distribution with mean \code{mu} and
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
#'    added to the \code{rigpsreg} package. Not used in Bell-Touchard distribution.
#'
#' @return \code{dbeto} returns the probability mass function, \code{pbeto} gives the distribution
#'     function, \code{qbeto} gives the quantile function, and \code{rbeto} generates random
#'     observations.
#'
#'     \code{beto} returns a \code{"count"} object mainly used to fit the Bell-Touchard distribution
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
#'    distribution function, quantile function, and a random number generator for the BerG
#'    distribution parameterized in terms of the mean (\code{mu}) and the variance (\code{sigma2}).
#'    This distribution was introduced by Bourguignon and Medeiros (2022) and
#'    also applied (in its reparametrized version) by Kokonendji, Medeiros, and Bourguignon (2024).
#'
#' @references Castellares, F., Lemonte, A. J., and Moreno–Arenas, G. (2020). On the two-parameter
#'     Bell–Touchard discrete distribution. \emph{Communications in Statistics-Theory and Methods},
#'     \bold{49}, 4834-4852.
#'
#' @references Kokonendji, C. C., Medeiros, R. M. R., and Bourguignon, M. (2024). Mean and variance
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
#' sigma2 <- 3
#'
#' y <- rbeto(n, mu, sigma2)
#'
#' # mean and variance theoretical vs empirical comparison
#' cat("\n mu:", mu, "--- sample mean:", round(mean(y), 3),
#'     "\n sigma2:", sigma2, "--- sample variance:", round(var(y), 3), "\n")
#'
#' # probability mass function theoretical vs empirical comparison
#' yvals <- sort(unique(y))
#' xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf", col = "white")
#' points(xcoords, dbeto(yvals, mu, sigma2), pch = 16, col = 2, type = "b")
#'
#' # distribution function theoretical vs empirical comparison
#' plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
#' curve(pbeto(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
#'
#' # quantile function theoretical vs empirical comparison
#' plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#'      pch = 16, xlab = "p", ylab = "Quantile")
#' curve(qbeto(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)


## Probability mass function
#' @rdname beto
#' @export
dbeto <- function(x, mu, sigma2, log.p = FALSE, ...){

  if (is.vector(x))
    x <- matrix(x, nrow = length(x))

  n <- dim(x)[1]
  mu <- matrix(mu, nrow = n)
  sigma2 <- matrix(sigma2, nrow = n)

  DI <- sigma2 / mu # transform from variance to dispersion index

  alpha <- DI - 1
  phi <- mu / (alpha * exp(alpha))

  pmf <- matrix(-Inf, nrow = n)

  # NaN index
  pmf[which(mu <= 0 | DI <= 1, arr.ind = TRUE)] <- NaN

  # Positive pmf index
  id1 <- which(x == 0 & !is.nan(pmf), arr.ind = TRUE)
  id2 <- which(x > 0 & !is.nan(pmf), arr.ind = TRUE)

  pmf[id1] <- phi[id1] * (1 - exp(alpha[id1])) + x[id1] * log(alpha[id1])
  pmf[id2] <- phi[id2] * (1 - exp(alpha[id2])) + x[id2] * log(alpha[id2]) +
    log(as.numeric(mapply(aux_fun_beto, x[id2], phi = phi[id2]))) -
    log(factorial(x[id2]))

  if(!log.p) pmf <- exp(pmf)
  as.vector(pmf)

}

## Cumulative distribution function (inerited from gamlss)
#' @rdname beto
#' @export
pbeto <- function(q, mu, sigma2, lower.tail = TRUE, log.p = FALSE, ...)
{

  if (is.vector(q))
    q <- matrix(q, nrow = length(q))

  n <- dim(q)[1]

  mu <- matrix(mu, nrow = n)
  sigma2 <- matrix(sigma2, nrow = n)

  cdf <- matrix(0, nrow = n)

  # NaN index
  cdf[which(mu <= 0 | sigma2 <= mu, arr.ind = TRUE)] <- NaN

  # Positive pmf index
  id <- which(q >= 0 & !is.nan(cdf), arr.ind = TRUE)

  fn <- function(q, mu, sigma2) sum(dbeto(0:q, mu, sigma2))

  Vcdf <- Vectorize(fn)
  cdf[id] <- Vcdf(q = q[id], mu = mu[id], sigma2 = sigma2[id])

  if(!lower.tail) cdf <- 1 - cdf
  if(log.p) cdf <- exp(cdf)

  as.numeric(cdf)

}

## Quantile function (inerited from gamlss)
#' @rdname beto
#' @export
qbeto <- function(p, mu, sigma2, lower.tail = TRUE, ...)
{

  max.value = 10000

  if(!lower.tail)
    p <- 1 - p

  if (is.vector(p))
    p <- matrix(p, nrow = length(p))

  n <- dim(p)[1]

  mu <- matrix(mu, nrow = n)
  sigma2 <- matrix(sigma2, nrow = n)

  qtl <- matrix(0, nrow = n)

  for (i in seq(along = p)) {
    cumpro <- 0
    if (p[i] + 1e-09 >= 1)
      qtl[i] <- Inf
    else {
      for (j in seq(from = 0, to = max.value)) {
        cumpro <- pbeto(j, mu, sigma2)
        qtl[i] <- j
        if (p[i] <= cumpro)
          break
      }
    }
  }
  qtl
}

## Random generation
#' @rdname beto
#' @export
rbeto <- function(n, mu, sigma2, ...){

  y <- rep(0, n)

  mu <- matrix(mu, n)
  sigma2 <- matrix(sigma2, n)

  DI <- sigma2 / mu # transform from variance to dispersion index

  alpha <- DI - 1
  phi <- mu / (alpha * exp(alpha))

  N <- stats::rpois(n, phi * (exp(alpha) - 1))

  id1 <- which(N > 0)
  id2 <- which(N == 0)

  X <- mapply(extraDistr::rtpois, n = N[id1], lambda = alpha[id1], a = 0)

  y[id1] <- mapply(sum, X)

  y
}

#' @rdname beto
#' @export
beto <- function(mu.link = "log", sigma2.link = "log"){

  out <- list()

  # Abbreviation
  out$abb <- "beto"

  # Name
  out$name <- "Reparametrized Bell-Touchard"

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

    sigma0 <- sqrt(stats::var(x) + max(mu0))
    gamma0 <- c(g2(sigma0), rep(0, k - 1))
    c(beta0, gamma0)

  }

  # Derivatives ----------------------------------------------------------------
  out$has_dl <- FALSE

  out$dl.dmu <- function(x, mu, sigma2, ...){NULL}
  out$dl2.dmu2 <- function(x, mu, sigma2, ...){NULL}
  out$dl.dsigma <- function(x, mu, sigma2, ...){NULL}
  out$dl2.dsigma2 <- function(x, mu, sigma2, ...){NULL}
  out$dl2.dmudsigma <- function(x, mu, sigma2, ...){NULL}

  structure(out, class = "count")

}


aux_fun_beto <- function(y, phi){
  sum(c(0, copula::Stirling2.all(y)) * phi^(0:y))
}



# ## Verification ----------------------------------------------------------------
# n <- 10000
# mu <- runif(1, 1, 6)
# sigma <- runif(1, mu, 8)
# y <- rbeto(n, mu, sigma)
# print(paste("mu:", round(mu, 3), "sample mean:", round(mean(y), 3)))
# print(paste("Standard dev.:", round(sigma, 3), "sample sd:", round(sd(y), 3)))
#
# #par(mfrow = c(1, 3))
# yvals <- sort(unique(y))
# xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf")
# points(xcoords, dbeto(yvals, mu, sigma), pch = 16, col = 2, type = "b")
#
# plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
# curve(pbeto(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
#
# plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#      pch = 16, xlab = "p", ylab = "Quantile")
# curve(qbeto(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
# par(mfrow = c(1, 1))


