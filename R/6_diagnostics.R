#' @name envelope
#'
#' @title Envelope Plot
#'
#' @description Provides the simulated envelope plot of the residuals resulting from the
#'     \code{\link{mvreg}} regression fit.
#'
#' @param object a \code{"\link{mvreg}"} object.
#' @param x an \code{"envel_mvreg"} object.
#' @param type character; specifies which residual should be produced in the
#'     envelope plot. The available options are \code{"quantile"} (default),
#'     \code{"pearson"} ((y - mean) / sd), and \code{"response"} (raw residuals, y - mean).
#' @param nsim the number of replicates. The default is \code{nsim = 99}.
#' @param level level of the sample quantile of the residual used in the
#'     construction of confidence bands.
#' @param plot logical; if \code{TRUE}, the envelope plot of the residuals is displayed.
#' @param progressBar logical; if \code{TRUE}, a progress bar is displayed giving
#'     the progress of making the graph. It can slow down the function
#'     considerably in applications with a large sample size.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{mvreg_envel} returns an \code{"mvreg_envel"} object which consists of
#'     a list with the following components:
#' \describe{
#'   \item{residuals}{a list with the quantile, pearson, and raw residuals resulting from the
#'       fit of the regression model.}
#'   \item{simulation}{a list whose components are matrices containing the
#'       ordered quantile, pearson, and raw residuals of the simulation for the
#'       plot envelope.}
#'
#'  The method \code{plot} makes the envelope plot.
#'  }
#'
NULL

#' @rdname envelope
#' @export
envelope <- function(object, nsim = 99, progressBar = TRUE, plot = TRUE, ...)
{

  mu.link <- object$mu.link
  sigma2.link <- object$sigma2.link

  # Family
  family <- get(object$count)(mu.link, sigma2.link)
  rdist <- get(paste0("r", family$abb))

  ## Model specifications
  y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y
  n <- length(y)
  X <- if(is.null(object$x)) stats::model.matrix(object$terms$mean, stats::model.frame(object)) else object$y
  p <- NCOL(X)

  mu <- object$fitted.values
  nu <- object$nu

  Z <- NULL
  k <- 0
  sigma2 <- NULL
  if(family$npar > 1){
    Z <- if(is.null(object$x)) stats::model.matrix(object$terms$variance, stats::model.frame(object)) else object$y
    k <- NCOL(Z)
    sigma2 <- object$sigma2
  }

  ## Residuals
  rq <- stats::residuals(object, type = "quantile")
  rp <- stats::residuals(object, type = "pearson")
  rr <- stats::residuals(object, type = "response")

  ## Simulation
  rq_s <- matrix(0, n, nsim)
  rp_s <- matrix(0, n, nsim)
  rr_s <- matrix(0, n, nsim)

  pb <- utils::txtProgressBar(1, nsim, style = 3)
  for(i in 1:nsim){

    # Estimatives
    ctrl <- control_fit()
    ctrl$hessian <- FALSE
    ctrl$start <- c(object$coefficients$mean,
                    object$coefficients$variance,
                    object$nu)

    error <- TRUE
    count <- 0
    while (error == TRUE & count < 50){
      # Simulated sample
      y.tilde <- rdist(n, mu = mu, sigma2 = sigma2, nu = nu)

      out <- suppressWarnings(try(mle(y.tilde, X, Z, family, control = ctrl), silent = TRUE))

      error <- unique(grepl("Error", out))
    }

    ## Coefficients
    beta <- out$par[1:p]

    ## Fitted values
    out$fitted.values <- c(make.link.ldrm(mu.link)$inv(X%*%beta))
    out$sigma2 <- NULL

    if(family$npar > 1){
      gamma <- out$par[1:k + p]
      out$sigma2 <- c(make.link.ldrm(sigma2.link)$inv(Z%*%gamma))
    }

    out$residuals <- y.tilde - out$fitted.values
    out$nu <- nu
    out$count <- family$abb
    out$y <- y.tilde
    class(out) <- "mvreg"

    # Empirical residuals
    rq_s[, i] <- sort(stats::residuals(out, type = "quantile"))
    rp_s[, i] <- sort(stats::residuals(out, type = "pearson"))
    rr_s[, i] <- sort(stats::residuals(out, type = "response"))

    if(progressBar){ utils::setTxtProgressBar(pb, i)}
  }

  out <- list(residuals = list(quantile = rq,
                               pearson = rp,
                               response = rr),
              simulation = list(quantile = rq_s,
                                pearson = rp_s,
                                response = rr_s))

  class(out) <- "envelope"

  if (plot)
    plot(out)

  out

}

#' @rdname envelope
#' @export
print.envelope <- function(x, ...){
  cat("\nA 'mvreg_envel' object\n\n")
  utils::str(x)
}

#' @rdname envelope
#' @export
plot.envelope <- function(x, type = c("quantile", "pearson", "response"),
                            level = 0.95, ...){


  alpha <- (1 - level)/2
  type <- match.arg(type)

  if(type == "quantile"){
    res <- x$residuals$quantile
    n <- length(res)
    rs <- x$simulation$quantile

    # alpha and 1 - alpha quantiles
    lower <- apply(rs, 1, stats::quantile, probs = alpha)
    upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

    # Median
    md <- apply(rs, 1, stats::quantile, probs = 0.5)

    # Theoretical quantiles for the normal distribution
    qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

    gp <- cbind(qq, sort(res), md, lower, upper)
    graphics::plot(gp[,1], gp[,2],
                   xlab = "Normal quantiles", ylab = "Randomized quantile residuals", type = "n", ...)
    graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                       c(gp[,4], sort(gp[,5], decreasing = T)),
                       col = "lightgray", border=NA)
    graphics::lines(gp[,1], gp[,3], lty = 2)
    graphics::points(gp[,1], gp[,2], pch = "+")
    graphics::box()
  }

  if(type == "pearson"){

    res <- x$residuals$pearson
    n <- length(res)
    rs <- x$simulation$pearson

    # alpha and 1 - alpha quantiles
    lower <- apply(rs, 1, stats::quantile, probs = alpha)
    upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

    # Median
    md <- apply(rs, 1, stats::quantile, probs = 0.5)

    # Theoretical quantiles for the normal distribution
    qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

    gp <- cbind(qq, sort(res), md, lower, upper)
    graphics::plot(gp[,1], gp[,2],
                   xlab = "Normal quantiles", ylab = "Pearson residuals", type = "n", ...)
    graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                       c(gp[,4], sort(gp[,5], decreasing = T)),
                       col = "lightgray", border=NA)
    graphics::lines(gp[,1], gp[,3], lty = 2)
    graphics::points(gp[,1], gp[,2], pch = "+")
    graphics::box()
  }

  if(type == "response"){
    res <- x$residuals$response
    n <- length(res)
    rs <- x$simulation$response

    # alpha and 1 - alpha quantiles
    lower <- apply(rs, 1, stats::quantile, probs = alpha)
    upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

    # Median
    md <- apply(rs, 1, stats::quantile, probs = 0.5)

    # Theoretical quantiles for the normal distribution
    qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

    gp <- cbind(qq, sort(res), md, lower, upper)
    graphics::plot(gp[,1], gp[,2],
                   xlab = "Normal quantiles", ylab = "Raw residuals", type = "n", ...)
    graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                       c(gp[,4], sort(gp[,5], decreasing = T)),
                       col = "lightgray", border=NA)
    graphics::lines(gp[,1], gp[,3], lty = 2)
    graphics::points(gp[,1], gp[,2], pch = "+")
    graphics::box()
  }

  invisible(x)

}

