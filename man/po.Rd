% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_poisson.R
\name{po}
\alias{po}
\alias{dpo}
\alias{ppo}
\alias{qpo}
\alias{rpo}
\title{The Poisson Distribution}
\usage{
dpo(x, mu, log.p = FALSE, ...)

ppo(q, mu, lower.tail = TRUE, log.p = FALSE, ...)

qpo(p, mu, lower.tail = TRUE, ...)

rpo(n, mu, ...)

po(mu.link = "log", ...)
}
\arguments{
\item{x}{vector of non-negative integer quantiles.}

\item{mu}{numeric; non-negative mean.}

\item{log.p}{logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.}

\item{...}{possible extra parameters, in addition to the mean and variance, in a distribution
added to the \code{rigpsreg} package. Not used in Poisson distribution.}

\item{q}{vector of quantiles.}

\item{lower.tail}{logical; if TRUE (default), probabilities are \code{P(X <= x)}, otherwise, \code{P(X > x)}.}

\item{p}{vector of probabilities.}

\item{n}{number of random values to return.}

\item{mu.link}{defines the link function for the mean regression model in the \code{\link{mvreg}}
function, with \code{"log"} link as the default. Other links are \code{"identity"} and \code{"sqrt"}.}
}
\value{
\code{dpo} returns the probability mass function, \code{ppo} gives the distribution
    function, \code{qpo} gives the quantile function, and \code{rpo} generates random
    observations.

    \code{po} returns a \code{"count"} object mainly used to fit the Poisson distribution
        via \code{\link{mvreg}} function. It returns a list with the following components:
        \describe{
          \item{abb}{abbreviation of the distribution in the \code{mvcount} package.}
          \item{name}{capitalized name of the distribution.}
          \item{npar}{number of parameters of the distribution.}
          \item{mu.link, sigma2.link}{specified link functions for the mean (\code{mu}) and the
              variance (\code{sigma2}) parameters, respectively.}
          \item{constraint}{a function \code{function(mu, sigma2)} that returns \code{TRUE} if
              the values of \code{mu} and \code{sigma} satisfy the constraint between the mean
              and variance of the distribution.}
          \item{start}{a function \code{function(x, X, Z)} which provides initial values for the
              regression coefficients of a double regression fit of the distribution.}
          \item{has_dl}{logical; if \code{TRUE}, it informs that the model has first and second
              order derivatives implemented; otherwise, they are obtained numerically via \code{\link{optim}}.}
          \item{dl.dmu,dl.dsigma}{functions \code{f(x, mu, sigma)} that provide the partial
              derivatives of the regression coefficients associated with the mean and variance,
              respectively. Its use is optional in the mvcount package. If they are not
              implemented, they must return \code{NULL}.}
         \item{dl2.dmu2, dl2,dsigma2, dl2.dmudsigma}{functions \code{f(x, mu, sigma)} that provide
              the second order partial derivatives of the regression coefficients associated with
              the mean and variance. Its use is optional in the mvcount package. If they are not
              implemented, they must return \code{NULL}.}
        }
}
\description{
Probability mass function, distribution function, quantile function, and random generation
for the Poisson distribution with mean \code{mu}.
}
\details{
This set of functions are the same as the \code{dpois}, \code{ppois}, \code{qpois},
    and \code{rpois} \code{R} base functions. It provides the Poisson distribution for use
    in the \code{rigpsreg} package.
}
\examples{
n <- 10000

mu <- 10

y <- rpo(n, mu)

# mean and variance theoretical vs empirical comparison
cat("\n mu:", mu, "--- sample mean:", round(mean(y), 3),
    "\n sigma2:", mu, "--- sample variance:", round(var(y), 3), "\n")

# probability mass function theoretical vs empirical comparison
yvals <- sort(unique(y))
xcoords <- barplot(prop.table(table(y)), xlab = "y", ylab = "Pmf", col = "white")
points(xcoords, dpo(yvals, mu), pch = 16, col = 2, type = "b")

# distribution function theoretical vs empirical comparison
plot(ecdf(y), xlab = "y", ylab = "Cdf", main = " ")
curve(ppo(x, mu), add = TRUE, col = 2, lwd = 2)

# quantile function theoretical vs empirical comparison
plot(seq(0.01, 0.99, 0.001), quantile(y, seq(0.01, 0.99, 0.001)),
     pch = 16, xlab = "p", ylab = "Quantile")
curve(qpo(x, mu), add = TRUE, col = 2, lwd = 2)
}
