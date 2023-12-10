# The Reparameterized Negative Binomial Distribution --------------------------------

#' @name nb
#'
#' @title The Reparameterized Negative Binomial Distribution
#'
#' @description Probability mass function, distribution function, quantile function, and random generation
#' for the reparameterized negative binomial (NB) distribution with mean \code{mu} and
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
#' @param sigma2 numeric; variance (greather than \code{mu}),
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P(X <= x)}, otherwise, \code{P(X > x)}.
#' @param log.p 	logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... possible extra parameters, in addition to the mean and variance, in a distribution
#'    added to the \code{rigpsreg} package. Not used in NB distribution.
#'
#' @return \code{dnb} returns the probability mass function, \code{pnb} gives the distribution
#'     function, \code{qnb} gives the quantile function, and \code{rnb} generates random
#'     observations.
#'
#'     \code{nb} returns a \code{"count"} object mainly used to fit the NB distribution
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
#'    distribution function, quantile function, and a random number generator for the RNB
#'    distribution parameterized in terms of the mean (\code{mu}) and the variance (\code{sigma2}).
#'    This parameterization was introduced by Joe and Zhu (2005) and
#'    also applied by Kokonendji, Medeiros, and Bourguignon (2024).
#'
#' @references Joe, H., and Zhu, R. (2005). Generalized Poisson distribution: the property of
#'     mixture of Poisson and comparison with negative binomial distribution.
#'     \emph{Biometrical Journal: Journal of Mathematical Methods in Biosciences}, \bold{47}, 219-229.
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
#' mu <- 2.5
#' sigma2 <- 4.3
#'
#' y <- rnb(n, mu, sigma2)
#'
#' # mean and variance theoretical vs empirical comparison
#' cat("\n mu:", mu, "--- sample mean:", round(mean(y), 3),
#'     "\n sigma2:", sigma2, "--- sample variance:", round(var(y), 3), "\n")
#'
#' # probability mass function theoretical vs empirical comparison
#' yvals <- sort(unique(y))
#' xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf", col = "white")
#' points(xcoords, dnb(yvals, mu, sigma2), pch = 16, col = 2, type = "b")
#'
#' # distribution function theoretical vs empirical comparison
#' plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
#' curve(pnb(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
#'
#' # quantile function theoretical vs empirical comparison
#' plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#'      pch = 16, xlab = "p", ylab = "Quantile")
#' curve(qnb(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)


## Probability mass function
#' @rdname nb
#' @export
dnb <- function(x, mu, sigma2, log.p = FALSE, ...){
  sigma <- sqrt(sigma2)
  p <- mu / (sigma^2)
  size <- mu^2 / (sigma^2 - mu)
  suppressWarnings(stats::dnbinom(x, size = size, prob = p, log = log.p))
}

## Cumulative distribution function (inerited from gamlss)
#' @rdname nb
#' @export
pnb <- function(q, mu, sigma2, lower.tail = TRUE, log.p = FALSE, ...)
{
  sigma <- sqrt(sigma2)
  p <- mu / (sigma^2)
  size <- mu^2 / (sigma^2 - mu)
  stats::pnbinom(q, size = size, prob = p, lower.tail = lower.tail, log.p = log.p)
}

## Quantile function (inerited from gamlss)
#' @rdname nb
#' @export
qnb <- function(p, mu, sigma2, lower.tail = TRUE, ...)
{
  sigma <- sqrt(sigma2)
  prob <- mu / (sigma^2)
  size <- mu^2 / (sigma^2 - mu)

  stats::qnbinom(p, size = size, prob = prob, lower.tail = lower.tail)
}

## Random generation
#' @rdname nb
#' @export
rnb <- function(n, mu, sigma2, ...){
  sigma <- sqrt(sigma2)
  p <- mu / (sigma^2)
  size <- mu^2 / (sigma^2 - mu)
  stats::rnbinom(n, size = size, prob = p)
}

#' @rdname nb
#' @export
nb <- function(mu.link = "log", sigma2.link = "log"){

  out <- list()

  # Abbreviation
  out$abb <- "nb"

  # Name
  out$name <- "Negative binomial"

  # Number of parameters
  out$npar <- 2

  # Link functions
  out$mu.link <- mu.link
  out$sigma2.link <- sigma2.link

  # Constraint -----------------------------------------------------------------
  out$constraint <- function(mu, sigma2, ...) (mu > 0) & (sigma2 > 0)

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

  out$dl.dmu <- function(x, mu, sigma2, ...){NULL}
  out$dl2.dmu2 <- function(x, mu, sigma2, ...){NULL}
  out$dl.dsigma <- function(x, mu, sigma2, ...){NULL}
  out$dl2.dsigma2 <- function(x, mu, sigma2, ...){NULL}
  out$dl2.dmudsigma <- function(x, mu, sigma2, ...){NULL}

  structure(out, class = "count")

}

# ## Verification ----------------------------------------------------------------
# n <- 10000
# mu <- runif(1, 1, 6)
# sigma <- runif(1, mu, 8)
# y <- nb(n, mu, sigma)
#
# print(paste("mu:", round(mu, 3), "sample mean:", round(mean(y), 3)))
# print(paste("Standard dev.:", round(sigma, 3), "sample sd:", round(sd(y), 3)))
#
# par(mfrow = c(1, 3))
# yvals <- sort(unique(y))
# xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf")
# points(xcoords, dnb(yvals, mu, sigma), pch = 16, col = 2, type = "b")
#
# plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
# curve(pnb(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
# plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#      pch = 16, xlab = "p", ylab = "Quantile")
# curve(qnb(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
# par(mfrow = c(1, 1))

# ## Verification ----------------------------------------------------------------
# n <- 10000
# mu <- runif(1, 1, 6)
# sigma <- runif(1, mu * (mu + 1), mu * (mu + 1) + 1)
# U <- SE <- matrix(NA, 1000, 2)
# COV <- c()
# for(i in 1:1000){
#  y <- rribe(n, mu, sigma)
#  dl <- dlribe(y, mu, sigma)
#  U[i, ] <- apply(matrix(c(dl$dl.dmu, dl$dl.dsigma), ncol = 2), 2, sum)
#  SE[i, ] <- sqrt(-apply(matrix(c(dl$dl2.dmu2, dl$dl2.dsigma2), ncol = 2), 2, sum))
#  COV[i] <- sum(-dl$dl2.dmudsigma)
#  print(i)
# }
# hist(U[,1], main = " ", xlab = expression(U[mu]), prob = TRUE)
# curve(dnorm(x, 0, mean(SE[,1])), add = TRUE, col = 2, lwd = 2)
# hist(U[,2], main = " ", xlab = expression(U[sigma]), prob = TRUE)
# curve(dnorm(x, 0, mean(SE[,2])), add = TRUE, col = 2, lwd = 2)
# cov(U)
# matrix(c(mean(SE[,1]^2), mean(COV), mean(COV), mean(SE[,2]^2)), 2, 2)
