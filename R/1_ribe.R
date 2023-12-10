# The inflated Bernoulli distribution ------------------------------------------

#' @name ribe
#'
#' @title The Reparameterized Inflated Bernoulli Distribution
#'
#' @description Probability mass function, distribution function, quantile function, and random generation
#' for the reparameterized inflated Bernoulli (RIBE) distribution with mean \code{mu} and
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
#'    added to the \code{rigpsreg} package. Not used in RIBE distribution.
#'
#' @return \code{dribe} returns the probability mass function, \code{pribe} gives the distribution
#'     function, \code{qribe} gives the quantile function, and \code{rribe} generates random
#'     observations.
#'
#'     \code{ribe} returns a \code{"count"} object mainly used to fit the RIBE distribution
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
#'    distribution function, quantile function, and a random number generator for the RIBE
#'    distribution parameterized in terms of the mean (\code{mu}) and the variance (\code{sigma2}).
#'    This parameterization was introduced by Kokonendji, Medeiros, and Bourguignon (2024).
#'
#' @references Kolev, N., Minkova, L. and Neytchev, P. (2000) Inflated-parameter family of generalized
#'     power series distributions and their application in analysis of overdispersed insurance
#'     data. \emph{ARCH Research Clearing House}, \bold{2}, 295-320.
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
#' mu <- 1.5
#' sigma2 <- 2
#'
#' y <- rribe(n, mu, sigma2)
#'
#' # mean and variance theoretical vs empirical comparison
#' cat("\n mu:", mu, "--- sample mean:", round(mean(y), 3),
#'     "\n sigma2:", sigma2, "--- sample variance:", round(var(y), 3), "\n")
#'
#' # probability mass function theoretical vs empirical comparison
#' yvals <- sort(unique(y))
#' xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf", col = "white")
#' points(xcoords, dribe(yvals, mu, sigma2), pch = 16, col = 2, type = "b")
#'
#' # distribution function theoretical vs empirical comparison
#' plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
#' curve(pribe(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
#'
#' # quantile function theoretical vs empirical comparison
#' plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#'      pch = 16, xlab = "p", ylab = "Quantile")
#' curve(qribe(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)


## Probability mass function
#' @rdname ribe
#' @export
dribe <- function(x, mu, sigma2, log.p = FALSE, ...){

  sigma <- sqrt(sigma2)

  if (is.vector(x))
    x <- matrix(x, nrow = length(x))

  n <- dim(x)[1]
  mu <- matrix(mu, nrow = n)
  sigma <- matrix(sigma, nrow = n)

  pmf <- matrix(-Inf, nrow = n)

  # NaN index
  pmf[which(mu <= 0 | sigma^2 <= abs(mu * (mu - 1)), arr.ind = TRUE)] <- NaN

  # Positive pmf index
  id1 <- which(x == 0  & !is.nan(pmf), arr.ind = TRUE)
  id2 <- which(x > 0 & !is.nan(pmf), arr.ind = TRUE)

  pmf[id1] <- log(sigma[id1]^2 + mu[id1] - mu[id1]^2) - log(sigma[id1]^2 + mu[id1] + mu[id1]^2)
  pmf[id2] <- log(4) + 3 * log(mu[id2]) - 2 * log(sigma[id2]^2 + mu[id2]^2 + mu[id2]) +
    (x[id2] - 1) * (log(sigma[id2]^2 + mu[id2]^2 - mu[id2]) -
                      log(sigma[id2]^2 + mu[id2]^2 + mu[id2]))


  if(!log.p) pmf <- exp(pmf)

  as.vector(pmf)

}

## Cumulative distribution function
#' @rdname ribe
#' @export
pribe <- function(q, mu, sigma2, lower.tail = TRUE, log.p = FALSE, ...){

  sigma <- sqrt(sigma2)

  if (is.vector(q))
    q <- matrix(q, nrow = length(q))

  n <- dim(q)[1]

  mu <- matrix(mu, nrow = n)
  sigma <- matrix(sigma, nrow = n)

  cdf <- matrix(0, nrow = n)

  # NaN index
  cdf[which(mu <= 0 | sigma^2 <= abs(mu * (mu - 1)), arr.ind = TRUE)] <- NaN

  # Positive pmf index
  id <- which(q >= 0 & !is.nan(cdf), arr.ind = TRUE)

  cdf[id] <- 1 - ((2 * mu[id]^2) / (sigma[id]^2 + mu[id] + mu[id]^2)) *
    ((sigma[id]^2 + mu[id] * (mu[id] - 1)) /
       (sigma[id]^2 + mu[id] * (mu[id] + 1)))^floor(q[id])

  if(!lower.tail) cdf <- 1 - cdf
  if(log.p) cdf <- log(cdf)

  as.numeric(cdf)

}

## Quantile function
#' @rdname ribe
#' @export
qribe <- function(p, mu, sigma2, lower.tail = TRUE, ...){

  sigma <- sqrt(sigma2)

  if(!lower.tail)
    p <- 1 - p

  if (is.vector(p))
    p <- matrix(p, nrow = length(p))

  n <- dim(p)[1]

  mu <- matrix(mu, nrow = n)
  sigma <- matrix(sigma, nrow = n)

  qtl <- matrix(0, nrow = n)

  # NaN index
  qtl[which(mu <= 0 | sigma^2 <= abs(mu * (mu - 1)), arr.ind = TRUE)] <- NaN

  p0 <- (sigma^2 + mu - mu^2) / (sigma^2 + mu + mu^2)

  qtl[which(p < p0 & !is.nan(qtl), arr.ind = TRUE)] <- 0
  qtl[which(p == 1 & !is.nan(qtl), arr.ind = TRUE)] <- Inf

  id <- which(p >= p0 & p < 1 & !is.nan(qtl), arr.ind = TRUE)
  qtl[id] <- ceiling(
    (log(1 - p[id]) - log(2 * mu[id]^2) + log(sigma[id]^2 + mu[id] * (mu[id] + 1))) /
      (log(sigma[id]^2 + mu[id] * (mu[id] - 1)) -
         log(sigma[id]^2 + mu[id] * (mu[id] + 1)))
  )

  as.numeric(qtl)

}

## Random generation
#' @rdname ribe
#' @export
rribe <- function(n, mu, sigma2, ...){

  sigma <- sqrt(sigma2)

  y <- rep(0, n)

  mu <- matrix(mu, n)
  sigma <- matrix(sigma, n)

  theta <- 2 * (mu^2) / (sigma^2 - mu^2 + mu)
  rho <- (sigma^2 + mu^2 - mu) / (sigma^2 + mu^2 + mu)

  # Quantile function of the zero-truncated geometric
  qgeom1 <- function(u, rho) ceiling(log(1 - u) / log(rho))

  N <- stats::rbinom(n, 1, theta / (1 + theta))

  id1 <- which(N > 0)
  id2 <- which(N == 0)

  u <- mapply(stats::runif, N[id1])
  X <- qgeom1(u, rho[id1])

  y[id1] <- mapply(sum, X)

  y
}

#' @rdname ribe
#' @export
ribe <- function(mu.link = "log", sigma2.link = "log"){

  out <- list()

  # Abbreviation
  out$abb <- "ribe"

  # Name
  out$name <- "Reparameterized inflated Bernoulli"

  # Number of parameters
  out$npar <- 2

  # Link functions
  out$mu.link <- mu.link
  out$sigma2.link <- sigma2.link

  # Constraint -----------------------------------------------------------------
  out$constraint <- function(mu, sigma2, ...){
    (mu > 0) & (sigma2 >= abs(mu * (mu - 1)))
  }

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

    sigma0 <- stats::var(x) + max(max(mu0) * (max(mu0) - 1), 0)
    gamma0 <- c(g2(sigma0), rep(0, k - 1))
    c(beta0, gamma0)
  }

  # Derivatives ----------------------------------------------------------------
  out$has_dl <- FALSE

  out$dl.dmu <- function(x, mu, sigma2, ...){

    sigma <- sqrt(sigma2)

    if (is.vector(x))
      x <- matrix(x, nrow = length(x))

    n <- dim(x)[1]

    mu <- matrix(mu, nrow = n)
    sigma <- matrix(sigma, nrow = n)

    dl <- matrix(NA, nrow = n)
    dl[which(mu <= 0 | sigma^2 <= abs(mu * (mu - 1)), arr.ind = TRUE)] <- NaN

    id1 <- which(x == 0  & !is.nan(dl), arr.ind = TRUE)
    id2 <- which(x > 0 & !is.nan(dl), arr.ind = TRUE)

    # Special constant and its derivatives
    cMS <- (sigma^2 + mu^2 - mu) / (sigma^2 + mu^2 + mu)

    # First
    dcMS.dmu <- 2 * (mu^2 - sigma^2) / ((mu^2 + mu + sigma^2)^2)

    dl[id1] <- (1 - 2 * mu[id1]) / (sigma[id1]^2 + mu[id1] - mu[id1]^2) -
      (1 + 2 * mu[id1]) / (sigma[id1]^2 + mu[id1] + mu[id1]^2)
    dl[id2] <- 3 / mu[id2] - 2 * (1 + 2 * mu[id2]) / (sigma[id2]^2 + mu[id2]^2 + mu[id2]) +
      (x[id2] - 1) * dcMS.dmu[id2] / cMS[id2]

    dl

  }

  out$dl2.dmu2 <- function(x, mu, sigma2, ...){

    sigma <- sqrt(sigma2)

    if (is.vector(x))
      x <- matrix(x, nrow = length(x))

    n <- dim(x)[1]

    mu <- matrix(mu, nrow = n)
    sigma <- matrix(sigma, nrow = n)

    dl <- matrix(NA, nrow = n)
    dl[which(mu <= 0 | sigma^2 <= abs(mu * (mu - 1)), arr.ind = TRUE)] <- NaN

    id1 <- which(x == 0  & !is.nan(dl), arr.ind = TRUE)
    id2 <- which(x > 0 & !is.nan(dl), arr.ind = TRUE)

    # Special constant and its derivatives
    cMS <- (sigma^2 + mu^2 - mu) / (sigma^2 + mu^2 + mu)

    # First
    dcMS.dmu <- 2 * (mu^2 - sigma^2) / ((mu^2 + mu + sigma^2)^2)

    # Second
    d2cMS.dmu2 <- 4 * (-mu^3 + 3 * mu * sigma^2 + sigma^2) / ((mu^2 + mu + sigma^2)^3)

    dl[id1] <- - 4 * (mu[id1]^5 + 3 * mu[id1]^4 * sigma[id1]^2 +
                        mu[id1] * sigma[id1]^4 + sigma[id1]^6) /
      ((- mu[id1]^4 + mu[id1]^2 + 2 * mu[id1] * sigma[id1]^2 + sigma[id1]^4)^2)
    dl[id2] <- -3/(mu[id2]^2) +
      2 * (2 * mu[id2]^2 + 2 * mu[id2] - 2 * sigma[id2]^2 + 1) /
      ((mu[id2]^2 + mu[id2] + sigma[id2]^2)^2) +
      (x[id2] - 1) *
      (d2cMS.dmu2[id2] * cMS[id2] - dcMS.dmu[id2]^2) / (cMS[id2]^2)

    dl

  }

  out$dl.dsigma <- function(x, mu, sigma2, ...){

    sigma <- sqrt(sigma2)

    if (is.vector(x))
      x <- matrix(x, nrow = length(x))

    n <- dim(x)[1]

    mu <- matrix(mu, nrow = n)
    sigma <- matrix(sigma, nrow = n)

    dl <- matrix(NA, nrow = n)
    dl[which(mu <= 0 | sigma^2 <= abs(mu * (mu - 1)), arr.ind = TRUE)] <- NaN

    id1 <- which(x == 0  & !is.nan(dl), arr.ind = TRUE)
    id2 <- which(x > 0 & !is.nan(dl), arr.ind = TRUE)

    # Special constant and its derivatives
    cMS <- (sigma^2 + mu^2 - mu) / (sigma^2 + mu^2 + mu)

    # First
    dcMS.dsigma <- 4 * mu * sigma / ((mu^2 + mu + sigma^2)^2)

    dl[id1] <- 2 * sigma[id1] / (sigma[id1]^2 + mu[id1] - mu[id1]^2) -
      2 * sigma[id1] / (sigma[id1]^2 + mu[id1] + mu[id1]^2)
    dl[id2] <- - 4 * sigma[id2] / (sigma[id2]^2 + mu[id2]^2 + mu[id2]) +
      (x[id2] - 1) * dcMS.dsigma[id2] / cMS[id2]

    dl

  }

  out$dl2.dsigma2 <- function(x, mu, sigma2, ...){

    sigma <- sqrt(sigma2)

    if (is.vector(x))
      x <- matrix(x, nrow = length(x))

    n <- dim(x)[1]

    mu <- matrix(mu, nrow = n)
    sigma <- matrix(sigma, nrow = n)

    dl <- matrix(NA, nrow = n)
    dl[which(mu <= 0 | sigma^2 <= abs(mu * (mu - 1)), arr.ind = TRUE)] <- NaN

    id1 <- which(x == 0  & !is.nan(dl), arr.ind = TRUE)
    id2 <- which(x > 0 & !is.nan(dl), arr.ind = TRUE)

    # Special constant and its derivatives
    cMS <- (sigma^2 + mu^2 - mu) / (sigma^2 + mu^2 + mu)

    # First
    dcMS.dsigma <- 4 * mu * sigma / ((mu^2 + mu + sigma^2)^2)

    # Second
    d2cMS.dsigma2 <- 4 * mu * (mu^2 + mu - 3 * sigma^2) / ((mu^2 + mu + sigma^2)^3)

    dl[id1] <- - 4 * mu[id1]^2 * (mu[id1]^4 - mu[id1]^2 +
                                    2 * mu[id1] * sigma[id1]^2 + 3 * sigma[id1]^4) /
      ((- mu[id1]^4 + mu[id1]^2 + 2 * mu[id1] * sigma[id1]^2 + sigma[id1]^4)^2)
    dl[id2] <-
      - 4 * (mu[id2]^2 + mu[id2] - sigma[id2]^2) /
      ((mu[id2]^2 + mu[id2] + sigma[id2]^2)^2) +
      (x[id2] - 1) *
      (d2cMS.dsigma2[id2] * cMS[id2] - dcMS.dsigma[id2]^2) / (cMS[id2]^2)

    dl

  }

  out$dl2.dmudsigma <- function(x, mu, sigma2, ...){

    sigma <- sqrt(sigma2)

    if (is.vector(x))
      x <- matrix(x, nrow = length(x))

    n <- dim(x)[1]

    mu <- matrix(mu, nrow = n)
    sigma <- matrix(sigma, nrow = n)

    dl <- matrix(NA, nrow = n)
    dl[which(mu <= 0 | sigma^2 <= abs(mu * (mu - 1)), arr.ind = TRUE)] <- NaN

    id1 <- which(x == 0  & !is.nan(dl), arr.ind = TRUE)
    id2 <- which(x > 0 & !is.nan(dl), arr.ind = TRUE)

    # Special constant and its derivatives
    cMS <- (sigma^2 + mu^2 - mu) / (sigma^2 + mu^2 + mu)

    # First
    dcMS.dmu <- 2 * (mu^2 - sigma^2) / ((mu^2 + mu + sigma^2)^2)
    dcMS.dsigma <- 4 * mu * sigma / ((mu^2 + mu + sigma^2)^2)

    # Second
    d2cMS.dmudsigma <- 4 * sigma * (-3 * mu^2 - mu + sigma^2) / ((mu^2 + mu + sigma^2)^3)

    dl[id1] <- 8 * mu[id1] * sigma[id1] * (mu[id1]^4 + mu[id1] * sigma[id1]^2 + sigma[id1]^4) /
      ((- mu[id1]^4 + mu[id1]^2 + 2 * mu[id1] * sigma[id1]^2 + sigma[id1]^4)^2)
    dl[id2] <- 4 * (2 * mu[id2] + 1) * sigma[id2] /
      ((mu[id2]^2 + mu[id2] + sigma[id2]^2)^2) +
      (x[id2] - 1) *
      (d2cMS.dmudsigma[id2] * cMS[id2] - dcMS.dmu[id2] * dcMS.dsigma[id2]) / (cMS[id2]^2)

    dl

  }

  structure(out, class = "count")

}

# # Distributional verificarions -------------------------------------------------
#
# ### Without regressors
# n <- 10000
#
# mu <- runif(1, 1, 3)
# sigma2 <- runif(1, abs(mu * (mu - 1)), abs(mu * (mu - 1)) + 1)
#
# y <- rribe(n, mu, sigma2)
#
# print(paste("mu:", round(mu, 3), "sample mean:", round(mean(y), 3)))
# print(paste("standard dev.:", round(sqrt(sigma2), 3), "sample sd:", round(sd(y), 3)))
#
# par(mfrow = c(1, 3))
# yvals <- sort(unique(y))
# xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf")
# points(xcoords, dribe(yvals, mu, sigma2), pch = 16, col = 2, type = "b")
#
# plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
# curve(pribe(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
# plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
#      pch = 16, xlab = "p", ylab = "Quantile")
# curve(qribe(x, mu, sigma2), add = TRUE, col = 2, lwd = 2)
# par(mfrow = c(1, 1))
#
# ### With regressors
# n <- 100
#
# set.seed(1); X <- cbind(1, runif(n), rbinom(1, 1, 0.6))
# set.seed(1); Z <- cbind(1, runif(n))
#
# beta <- c(1, -0.6, 0.7)
# gamma <- c(2, 0.6)
#
# mu <- exp(X%*%beta)
# sigma2 <- exp(Z%*%gamma)
#
# rbind(`mu` = summary(as.numeric(mu)), `sigma2` = summary(as.numeric(sigma2)))
#
# all(ribe()$constraint(mu, sigma2))
#
# y <- matrix(NA, 1000, n)
# for(i in 1:1000){
#   y[i, ] <- rribe(n, mu, sigma2)
#   if(i %% 100 == 0) print(i)
# }
#
# i <- sample(1:n, 1)

# par(mfrow = c(1, 3))
# yvals <- sort(unique(y[, i]))
# xcoords <- barplot(prop.table(table(y[, i])), xlab = "y", ylab = "Pmf")
# points(xcoords, dribe(yvals, mu[i], sigma2[i]), pch = 16, col = 2, type = "b")

# plot(ecdf(y[, i]), xlab = "y", ylab = "Cdf", main = paste0("Replicate: ", i))
# curve(pribe(x, mu[i], sigma2[i]), add = TRUE, col = 2, lwd = 2)

# plot(seq(0.01, 0.99, 0.001), quantile(y[, i], seq(0.01, 0.99, 0.001)),
#      pch = 16, xlab = "p", ylab = "Quantile")
# curve(qribe(x, mu[i], sigma2[i]), add = TRUE, col = 2, lwd = 2)
# par(mfrow = c(1, 1))


# ## Derivative verifications ----------------------------------------------------
#
# n <- 500
# mu <- runif(1, 1, 3)
# sigma2 <- runif(1, abs(mu * (mu - 1)), abs(mu * (mu - 1)) + 1)
#
# U <- SE <- matrix(NA, 2000, 2)
# COV <- c()
# for(i in 1:2000){
#   y <- rribe(n, mu, sigma2)
#   dl <- ribe()
#
#   U[i, ] <- apply(matrix(c(dl$dl.dmu(y, mu, sigma2),
#                            dl$dl.dsigma(y, mu, sigma2)), ncol = 2), 2, sum)
#   SE[i, ] <- sqrt(-apply(matrix(c(dl$dl2.dmu2(y, mu, sigma2), dl$dl2.dsigma2(y, mu, sigma2)), ncol = 2), 2, sum))
#   COV[i] <- sum(-dl$dl2.dmudsigma(y, mu, sigma2))
#   print(i)
# }
#
# plot(density(U[,1]), main = " ", xlab = expression(U[mu]))
# curve(dnorm(x, 0, mean(SE[, 1])), add = TRUE, col = 2, lwd = 2)
#
# plot(density(U[,2]), main = " ", xlab = expression(U[sigma]))
# curve(dnorm(x, 0, mean(SE[, 2])), add = TRUE, col = 2, lwd = 2)
#
# cov(U)
# matrix(c(mean(SE[, 1]^2), mean(COV), mean(COV), mean(SE[, 2]^2)), 2, 2)
#
#
#
