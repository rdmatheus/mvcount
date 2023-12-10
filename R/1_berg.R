# The Reparameterized BerG Distribution ----------------------------

#' @name berg
#'
#' @title The Reparameterized BerG Distribution
#'
#' @description Probability mass function, distribution function, quantile function, and random generation
#' for the reparameterized BerG distribution with mean \code{mu} and
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
#'    added to the \code{rigpsreg} package. Not used in RIPO distribution.
#'
#' @return \code{dberg} returns the probability mass function, \code{pberg} gives the distribution
#'     function, \code{qberg} gives the quantile function, and \code{rberg} generates random
#'     observations.
#'
#'     \code{berg} returns a \code{"count"} object mainly used to fit the BerG distribution
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
#' @references Bourguignon, M., and de Medeiros, R. M. (2022). A simple and useful regression model
#'     for fitting count data. \emph{TEST}, \bold{31}, 790-827.
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
#' y <- rberg(n, mu, sigma2)
#'
#' # mean and variance theoretical vs empirical comparison
#' cat("\n mu:", mu, "--- sample mean:", round(mean(y), 3),
#'     "\n sigma2:", sigma2, "--- sample variance:", round(var(y), 3), "\n")
#'
#' # probability mass function theoretical vs empirical comparison
#' yvals <- sort(unique(y))
#' xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf", col = "white")
#' points(xcoords, dberg(yvals, mu, sigma2), pch = 16, col = 2, type = "b")
#'
#' # distribution function theoretical vs empirical comparison
#' plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
#' curve(pberg(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
#'
#' # quantile function theoretical vs empirical comparison
#' plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#'      pch = 16, xlab = "p", ylab = "Quantile")
#' curve(qberg(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)


## Probability mass function
#' @rdname berg
#' @export
dberg <- function(x, mu, sigma2, log.p = FALSE, ...)
{

  phi <- sigma2 / mu

  if (is.vector(x))
    x <- matrix(x, ncol = length(x))

  n <- dim(x)[2]
  mu <- matrix(mu, ncol = n)
  phi <- matrix(phi, ncol = n)

  pmf <- matrix(-Inf, ncol = n)

  # NaN index
  NaNid <- which(mu <= 0 | phi <= abs(mu - 1), arr.ind = TRUE)
  pmf[NaNid] <- NaN

  # Positive pmf index
  id1 <- which(x == 0 & mu > 0 & phi > abs(mu - 1) & !is.nan(pmf), arr.ind = TRUE)
  id2 <- which(x > 0 & mu > 0 & phi > abs(mu - 1) & !is.nan(pmf), arr.ind = TRUE)

  pmf[id1] <- log(1 - mu[id1] + phi[id1]) - log(1 + mu[id1] + phi[id1])
  pmf[id2]  <- log(4) + log(mu[id2]) + (x[id2] - 1) * log(mu[id2] + phi[id2] - 1) -
    (x[id2] + 1) * log(mu[id2] + phi[id2] + 1)

  if(nrow(pmf) == 1L)
    pmf <- as.vector(pmf)

  if(log.p) pmf else exp(pmf)
}

## Cumulative distribution function (inerited from gamlss)
#' @rdname berg
#' @export
pberg <- function(q, mu, sigma2, lower.tail = TRUE, log.p = FALSE, ...)
{

  phi <- sigma2 / mu

  q <- floor(q)

  if (is.vector(q))
    q <- matrix(q, ncol = length(q))

  n <- dim(q)[2]
  mu <- matrix(mu, ncol = n)
  phi <- matrix(phi, ncol = n)

  cdf <- matrix(0, ncol = n)

  # NaN index
  NaNid <- which(mu <= 0 | phi <= abs(mu - 1), arr.ind = TRUE)
  cdf[NaNid] <- NaN

  # Positive cdf index
  id <- which(q >= 0 & mu > 0 & phi > abs(mu - 1) & !is.nan(cdf), arr.ind = TRUE)
  cdf[id]  <- exp(log(1 - mu[id] + phi[id]) - log(1 + mu[id] + phi[id])) +
    exp(log(2 * mu[id]) - log(1 + mu[id] + phi[id]) +
          log(1 - ((mu[id] + phi[id] - 1) / (mu[id] + phi[id] + 1))^q[id]))

  if(nrow(cdf) == 1L)
    cdf <- as.vector(cdf)

  if(lower.tail) cdf else 1 - cdf
  if(log.p) cdf <- exp(cdf)

  as.numeric(cdf)

}

## Quantile function (inerited from gamlss)
#' @rdname berg
#' @export
qberg <- function(p, mu, sigma2, lower.tail = TRUE, ...)
{

  phi <- sigma2 / mu

  if(!lower.tail)
    p <- 1 - p

  if (is.vector(p))
    p <- matrix(p, ncol = length(p))

  n <- dim(p)[2]
  mu <- matrix(mu, ncol = n)
  phi <- matrix(phi, ncol = n)

  q <- matrix(0, ncol = n)

  # NaN index
  NaNid <- which(mu <= 0 | phi <= abs(mu - 1), arr.ind = TRUE)
  q[NaNid] <- NaN

  p0 <- exp(log(1 - mu + phi) - log(1 + mu + phi))

  id <- which(p >= p0 & mu > 0 & phi > abs(mu - 1) & !is.nan(q), arr.ind = TRUE)
  q[id] <- ceiling(
    (log(1 - p[id]) + log(1 + mu[id] + phi[id]) - log(2 * mu[id])) /
      (log(mu[id] + phi[id] - 1) - log(1 + mu[id] + phi[id]))
  )

  if(nrow(q) == 1L) as.vector(q) else q

}


## Random generation
#' @rdname berg
#' @export
rberg <- function(n, mu, sigma2, ...){

  u <- stats::runif(n)
  qberg(u, mu, sigma2)

}

#' @rdname berg
#' @export
berg <- function(mu.link = "log", sigma2.link = "log"){

  out <- list()

  # Abbreviation
  out$abb <- "berg"

  # Name
  out$name <- "Reparameterized BerG"

  # Number of parameters
  out$npar <- 2

  # Link functions
  out$mu.link <- mu.link
  out$sigma2.link <- sigma2.link

  # Constraint -----------------------------------------------------------------
  out$constraint <- function(mu, sigma2, ...) (mu > 0) & (sigma2 >= abs(mu * (mu - 1)))

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

# ## Distributional verifications
#
# ### Without regressors
# n <- 10000
#
# mu <- runif(1, 1, 6)
# sigma2 <- runif(1, mu, mu + 1)
#
# y <- rripo(n, mu, sigma2)
#
# print(paste("mu:", round(mu, 3), "sample mean:", round(mean(y), 3)))
# print(paste("sigma2:", sigma2, "sample var:", round(var(y), 3)))
#
# par(mfrow = c(1, 3))
# yvals <- sort(unique(y))
# xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf")
# points(xcoords, dripo(yvals, mu, sigma2), pch = 16, col = 2, type = "b")
#
# plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
# curve(pripo(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
# plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#     pch = 16, xlab = "p", ylab = "Quantile")
# curve(qripo(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
# par(mfrow = c(1, 1))
#
# ### With regressors
# n <- 100
#
# set.seed(1); X <- cbind(1, runif(n), rbinom(1, 1, 0.6))
# set.seed(1); Z <- cbind(1, runif(n))
#
# beta <- c(2, -0.6, 0.7)
# gamma <- c(2, 0.5)
#
# mu <- exp(X%*%beta)
# sigma2 <- exp(Z%*%gamma)
#
# rbind(`mu` = summary(as.numeric(mu)), `sigma2` = summary(as.numeric(sigma2)))
#
# all(ripo()$constraint(mu, sigma2))
#
# y <- matrix(NA, 1000, n)
# for(i in 1:1000){
#   y[i, ] <- rripo(n, mu, sigma2)
#   if(i %% 100 == 0) print(i)
# }
#
# i <- sample(1:n, 1)
#
# par(mfrow = c(1, 3))
# yvals <- sort(unique(y[, i]))
# xcoords <- barplot(prop.table(table(y[, i])), xlab = "y", ylab = "Pmf")
# points(xcoords, dripo(yvals, mu[i], sigma2[i]), pch = 16, col = 2, type = "b")
#
# plot(ecdf(y[, i]), xlab = "y", ylab = "Cdf", main = paste0("Replicate: ", i))
# curve(pripo(x, mu[i], sigma2[i]), add = TRUE, col = 2, lwd = 2)
#
# plot(seq(0.01, 0.99, 0.001), quantile(y[, i], seq(0.01, 0.99, 0.001)),
#      pch = 16, xlab = "p", ylab = "Quantile")
# curve(qripo(x, mu[i], sigma2[i]), add = TRUE, col = 2, lwd = 2)
# par(mfrow = c(1, 1))
#
# ## Derivative verifications ----------------------------------------------------
#
# n <- 500
# mu <- runif(1, 1, 6)
# sigma2 <- runif(1, mu, mu + 1)
#
# U <- SE <- matrix(NA, 1000, 2)
# COV <- c()
# for(i in 1:1000){
#   y <- rripo(n, mu, sigma2)
#   dl <- ripo()
#
#   U[i, ] <- apply(matrix(c(dl$dl.dmu(y, mu, sigma2),
#                            dl$dl.dsigma(y, mu, sigma2)), ncol = 2), 2, sum)
#   SE[i, ] <- sqrt(-apply(matrix(c(dl$dl2.dmu2(y, mu, sigma2),
#                                   dl$dl2.dsigma2(y, mu, sigma2)), ncol = 2), 2, sum))
#   COV[i] <- sum(-dl$dl2.dmudsigma(y, mu, sigma2))
#   print(i)
# }
#
# plot(density(U[,1]), main = " ", xlab = expression(U[mu]))
# curve(dnorm(x, 0, mean(SE[,1])), add = TRUE, col = 2, lwd = 2)
#
# plot(density(U[,2]), main = " ", xlab = expression(U[sigma]))
# curve(dnorm(x, 0, mean(SE[,2])), add = TRUE, col = 2, lwd = 2)
#
# cov(U)
# matrix(c(mean(SE[,1]^2), mean(COV), mean(COV), mean(SE[,2]^2)), 2, 2)
