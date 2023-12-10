#' @name mvreg-methods
#' @title Methods for \code{"mvreg"} objects
#' @param x,object an object of class \code{"mvreg"}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
NULL

# Print
#' @rdname mvreg-methods
#' @export
print.mvreg <- function(x, ...)
{

  p <- length(x$coefficients$mean)
  n <- x$nobs

  mu.link <- x$mu.link
  sigma2.link <- x$sigma2.link
  family <- get(x$count)(mu.link, sigma2.link)

  cat("\n---------------------------------------------------------------\n")
  cat(family$name, "distribution")
  cat("\n---------------------------------------------------------------")
  cat("\nCall:\n")
  print(x$call)
  if(x$convergence > 0) {
    cat("\nmodel did not converge\n")
  }else{
    cat("\nMean coefficients with", mu.link, "link:\n")
    print((x$coefficients)$mean)
    k <- 0
    if(family$npar > 1){
      k <- length(x$coefficients$variance)
      cat("\nVariance coefficients with", sigma2.link, "link:\n")
      print((x$coefficients)$variance)
    }

    if(family$npar > 2){
      cat("\nNu estimate:", x$nu)
    }

    cat("\nLog-lik value:",x$logLik,
        "\nAIC:", 2 * (p + k - x$logLik),
        "and BIC:", log(n) * (p + k) - 2 * x$logLik)
  }

  invisible(x)
}

# Summary
#' @rdname mvreg-methods
#' @export
summary.mvreg <- function(object, ...)
{

  X <- stats::model.matrix(object$terms$mean, stats::model.frame(object))
  p <- ncol(X)
  n <- object$nobs

  mu.link <- object$mu.link
  sigma2.link <- object$sigma2.link
  family <- get(object$count)(mu.link, sigma2.link)

  ## Names
  if (p > 1) {
    mean_names <- colnames(X)
  }else{
    if (mu.link == "identity"){
      mean_names <- "mu"
    }else{
      mean_names <- "g1(mu)"
    }
  }

  ## Summary for quantile residuals
  res <- stats::residuals(object, type = "quantile")
  skewness <- mean((res - mean(res))^3) / (stats::sd(res)^3)
  kurtosis <- mean((res - mean(res))^4) / (stats::sd(res)^4)
  TAB.residuals <- round(cbind(mean(res), stats::sd(res),
                               skewness, kurtosis), 6)
  colnames(TAB.residuals) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(TAB.residuals) <- " "

  # Summary for mu
  est.beta <- object$coef$mean
  se.beta <- sqrt(diag(object$vcov))[1:p]
  zval.beta <- est.beta/se.beta
  pval.beta <- 2 * stats::pnorm(abs(zval.beta), lower.tail = FALSE)

  TAB.mu <- cbind(Estimate = est.beta,
                  `Std. Error` = se.beta,
                  `z value` = zval.beta,
                  `Pr(>|z|)` = pval.beta)
  rownames(TAB.mu) <- mean_names

  k <- 0
  TAB.sigma2 <- NULL
  if(family$npar > 1){

    Z <- stats::model.matrix(object$terms$variance, stats::model.frame(object))
    k <- ncol(Z)

    if (k > 1) {
      disp_names <- colnames(Z)
    }else{
      if (sigma2.link == "identity"){
        disp_names <- "sigma2"
      }else{
        disp_names <- "g2(sigma2)"
      }
    }

    # Summary for sigma2
    est.gamma <- object$coef$var
    se.gamma <- sqrt(diag(object$vcov)[1:k + p])
    zval.gamma <- est.gamma/se.gamma
    pval.gamma <- 2 * stats::pnorm(abs(zval.gamma), lower.tail = FALSE)

    TAB.sigma2 <- cbind(Estimate = est.gamma,
                       `Std. Error` = se.gamma,
                       `z value` = zval.gamma,
                       `Pr(>|z|)` = pval.gamma)
    rownames(TAB.sigma2) <- disp_names

  }

  nu <- nu.se <- NULL
  if(family$npar > 2){
    nu <- object$nu
    nu.se <- sqrt(diag(object$vcov)[k + p + family$npar - 2])
  }

  # Pseudo R2
  y <- if(is.null(object$y)){
    stats::model.response(stats::model.frame(object))
  } else{
    object$y
  }

  pseudoR2 <- stats::cor(table(y), object$freq)^2

  # logLik, AIC, and BIC
  logLik <- object$logLik
  AIC <- 2 * (p + k - logLik)
  BIC <- log(n) * (p + k) - 2 * logLik

  out <- list(name = family$name,
              call = object$call,
              residuals = TAB.residuals,
              mu.link = object$mu.link,
              sigma2.link = object$sigma2.link,
              mean = TAB.mu,
              variance = TAB.sigma2,
              nu = nu,
              nu.se = nu.se,
              logLik = logLik,
              AIC = AIC,
              BIC = BIC,
              pseudoR2 = pseudoR2,
              upsilon = object$upsilon)

  class(out) <- "summary.mvreg"
  out
}

# Print summary
#' @rdname mvreg-methods
#' @export
print.summary.mvreg <- function(x, ...)
{
  cat("\n---------------------------------------------------------------\n")
  cat(x$name, "distribution")
  cat("\n---------------------------------------------------------------\n")
  cat("Call:\n")
  print(x$call)
  cat("\nSummary for quantile residuals:\n")
  print(x$residuals)
  cat("\n---------------------------------------------------------------\n")
  cat("Mean model with", x$mu.link, "link:")
  cat("\n---------------------------------------------------------------\n")
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$mean)
  if(!is.null(x$variance)){
    cat("\n---------------------------------------------------------------\n")
    cat("Variance model with", x$sigma2.link, "link:")
    cat("\n---------------------------------------------------------------\n")
    cat("\nCoefficients:\n")
    stats::printCoefmat(x$variance)
  }

  cat("\n---------------------------------------------------------------")
  if(!is.null(x$nu)){
    cat("\nNu estimate:", x$nu, "Std. Error:", x$nu.se)
  }
  cat("\nLog-lik value:",x$logLik,
      "\nAIC:", x$AIC,
      "and BIC:", x$BIC,
      "\nUpsilon statistic:", x$upsilon,
      "\nPseudo-R2:", x$pseudoR2)

  invisible(x)
}

# Parameter estimates
#' @rdname mvreg-methods
#' @export
#' @param model a character indicating which parameter coefficients are
#'   required, parameters for the \code{"mean"} or for the
#'   \code{"variance"} submodel. If \code{"full"} (default), a list with
#'   coefficients for the \code{mean} and for the \code{variance}
#'   model is returned.
coef.mvreg <- function(object,
                      model = c("full", "mean", "variance"), ...) {
  model <- match.arg(model)
  if(!is.null(object$sigma2)){
    switch(model,
           "full"        = list(
             mean       = (object$coef)$mean,
             variance = (object$coef)$var),
           "mean"       = (object$coef)$mean,
           "variance" = (object$coef)$var)
  }else{
    (object$coef)$mean
  }

}

#  Variance-covariance matrix
#' @rdname mvreg-methods
#' @export
vcov.mvreg <- function(object,
                      model = c("full", "mean", "variance"), ...) {

  model <- match.arg(model)
  covm <- object$vcov

  mu.link <- object$mu.link
  sigma2.link <- object$sigma2.link
  family <- get(object$count)(mu.link, sigma2.link)

  X <- stats::model.matrix(object$terms$mean, stats::model.frame(object))
  p <- ncol(X)

  if (p > 1) {
    mean_names <- colnames(X)
  }else{
    if (mu.link == "identity"){
      mean_names <- "mu"
    }else{
      mean_names <- "g1(mu)"
    }
  }

  if(family$npar > 1){
    Z <- stats::model.matrix(object$terms$variance, stats::model.frame(object))
    k <- ncol(Z)

    if (k > 1) {
      disp_names <- colnames(Z)
    }else{
      if (sigma2.link == "identity"){
        disp_names <- "sigma2"
      }else{
        disp_names <- "g2(sigma2)"
      }
    }

    nu_names <- NULL
    if(family$npar == 3){
      nu_names <- "nu"
    }else if(family$npar > 3){
      nu_names <- paste0("nu", 1:(family$npar - 2))
    }

    switch(model,
           "mean" = {
             covm <- covm[seq.int(length.out = p),
                          seq.int(length.out = p), drop = FALSE]
             colnames(covm) <- rownames(covm) <- mean_names
             covm
           },
           "variance" = {
             covm <- covm[seq.int(length.out = k) + p,
                          seq.int(length.out = k) + p, drop = FALSE]
             colnames(covm) <- rownames(covm) <- disp_names
             covm
           },
           "full" = {
             colnames(covm) <- rownames(covm) <- c(mean_names,
                                                   disp_names,
                                                   nu_names)
             covm
           })
  }else{
    covm <- covm[seq.int(length.out = p),
                 seq.int(length.out = p), drop = FALSE]
    colnames(covm) <- rownames(covm) <- mean_names
    covm
  }
}


#' Predict Method for LDRM Fits
#'
#' Obtains predictions from a fitted SDL regression object.
#'
#' @param object an \code{'mvreg'} object.
#' @param newdata optionally, a data frame in which to look for variables
#'     with which to predict. If omitted, the fitted linear predictors are
#'     used.
#' @param type the type of prediction required. The default is on the scale of
#'     the response variable \code{("response")}, that is, the fitted values
#'     (fitted means). The alternative \code{"variance"}
#'     provides the fitted variances. Finally,
#'     the option \code{"quantile"} gives the fitted quantiles in the order
#'     specified via \code{'at'}.
#' @param at the order of the quantile to be predicted if
#'     \code{type = "quantile"}. The default is to predict the median,
#'     that is, \code{at = 0.5}.
#' @param na.action function determining what should be done with missing
#'     values in \code{newdata}. The default is to predict \code{NA}.
#' @param ...  arguments passed to or from other methods.
#'
#' @return A vector of predictions.
#' @export
#'
predict.mvreg <- function(object, newdata = NULL,
                         type = c("response", "link", "variance", "quantile"),
                         at = 0.5,
                         na.action = stats::na.pass, ...)
{
  type <- match.arg(type)

  mu.link <- object$mu.link
  sigma2.link <- object$sigma2.link
  family <- get(object$count)(mu.link, sigma2.link)

  qdist <- match.fun(paste0("q", family$abb))

  qfun <- function(at, mu, sigma2, nu) {

    aux <- function(arg){
      p <- arg[1]
      mu <- arg[2]
      sigma2 <- arg[3]
      qdist(p, mu = mu, sigma2 = sigma2, nu = nu)
    }

    c(apply(cbind(at, mu, sigma2, nu), 1, aux))

  }


  if(missing(newdata)) {

    rval <- switch(type,
                   "response" = {
                     as.numeric(object$fitted.values)
                   },

                   "link" = {
                     as.numeric(make.link.ldrm(object$mu.link)$fun(object$fitted.values))
                   },

                   "variance" = {
                     if(family$npar > 1){
                       object$sigma2
                     }else{
                       NULL
                     }
                   },
                   "quantile" = {
                     mu <- as.numeric(object$fitted.values)
                     sigma2 <- as.numeric(object$sigma2)
                     nu <- object$nu

                     qfun(at, mu, sigma2, nu)
                   })

    return(rval)

  } else {

    tnam <- switch(type,
                   "response" = "mean",
                   "variance" = "full",
                   "quantile" = "full")

    mf <- stats::model.frame(stats::delete.response(object$terms[[tnam]]),
                             newdata, na.action = na.action)
    newdata <- newdata[rownames(mf), , drop = FALSE]

    if(type %in% c("response", "variance", "quantile"))
      X <- stats::model.matrix(stats::delete.response(object$terms$mean), mf)

    if(type %in% c("variance", "quantile")){
      if(family$npar > 1){
        Z <- stats::model.matrix(object$terms$variance, mf)
      }else{
        Z <- NULL
      }
    }


    rval <- switch(type,
                   "response" = {
                     as.numeric(make.link.ldrm(object$mu.link)$inv(
                       drop(X %*% object$coefficients$mean)))
                   },

                   "link" = {
                     mu <- as.numeric(make.link.ldrm(object$mu.link)$inv(
                       drop(X %*% object$coefficients$mean)))

                     make.link.ldrm(object$mu.link)$fun(mu)
                   },

                   "variance" = {
                     mu <- make.link.ldrm(object$mu.link)$inv(drop(X %*% object$coefficients$mean))
                     nu <- object$nu
                     if(family$npar > 1){
                       make.link.ldrm(object$sigma2.link)$inv(drop(Z %*% object$coefficients$variance))
                     }else{
                       mu
                     }

                   },
                   "quantile" = {
                     mu <- make.link.ldrm(object$mu.link)$inv(drop(X %*% object$coefficients$mean))
                     nu <- object$nu
                     if(family$npar > 1){
                       sigma2 <- make.link.ldrm(object$sigma2.link)$inv(drop(Z %*% object$coefficients$variance))
                       qfun(at, mu, sigma2, nu)
                     }else{
                       qfun(at, mu, NULL, NULL)
                     }


                   }
    )
    return(rval)

  }
}


# Log-likelihood
#' @rdname mvreg-methods
#' @export
logLik.mvreg <- function(object, ...) {
  structure(object$logLik,
            df = object$nobs - object$df.residual,
            class = "logLik")
}


# AIC
#' @export
#' @rdname mvreg-methods
AIC.mvreg <- function(object, ..., k = 2) {

  AIC <- - 2 * object$logLik + k *
    (length(object$coefficients$mean) +
       length(object$coefficients$variance) +
       length(object$nu))

  class(AIC) <- "AIC"
  return(AIC)
}

# Residuals
#' @name residuals.mvreg
#' @title Extract Model Residuals for a BerG Regression
#'
#' @param object an \code{'mvreg'} object.
#' @param type character; specifies which residual should be extracted.
#'     The available arguments are "quantile" (default), "pearson",
#'     and "response" (raw residuals, y - mu).
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
residuals.mvreg <- function(object,
                           type = c("quantile", "pearson", "response"), ...)
{

  ## raw response residuals and desired type
  res <- object$residuals
  type <- match.arg(type, c("quantile", "pearson", "response"))
  if(type == "response") return(res)

  ## extract fitted information
  y <- if(is.null(object$y)){
    stats::model.response(stats::model.frame(object))
  } else{
    object$y
  }

  mu <- stats::fitted(object)
  sigma2 <- object$sigma2
  nu <- object$nu

  mu.link <- object$mu.link
  sigma2.link <- object$sigma2.link
  family <- get(object$count)(mu.link, sigma2.link)

  pdist <- match.fun(paste0("p", family$abb))

  rqr <- function(y){

      n <- length(y)
      a <- vector()
      b <- vector()
      u <- vector()
      for(i in 1:n){
        a[i] <- pdist(y[i] - 1, mu[i], sigma2 = sigma2[i], nu = nu)
        b[i] <- pdist(y[i], mu[i], sigma2 = sigma2[i], nu = nu)
        u[i] <- stats::runif(1, a[i], b[i])
      }

    stats::qnorm(u)

  }

  res <- switch(type,

                "pearson" = {

                  if (family$npar < 2) {
                    as.numeric(res / sqrt(mu))
                  } else {
                    as.numeric(res / sqrt(sigma2))
                  }


                },

                "quantile" = {
                  rqr(y)
                })

  res
}

## Model frame
model.frame.mvreg <- function(formula, ...) {

  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()

}

## Model matrix
#' @export
#' @rdname mvreg-methods
model.matrix.mvreg <- function(object,
                              model = c("mean", "variance"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
  else stats::model.matrix(object$terms[[model]], stats::model.frame(object))
  return(rval)
}

# Plot
#' Diagnostic plots for the SDL Regression
#'
#' Five plots (selectable by \code{which}) are currently available:
#' a plot of residuals against fitted values, a plot of residuals against
#' the indexes, a Normal Q-Q plot, a barplot with comparisons of the
#' observed and fitted frequencies, and a plot of the sample autocorelations
#' of the residuals.
#'
#' @param x an object of class \code{mvreg}.
#' @param type character; specifies which residual should be produced in the
#'     envelope plot. The available options are \code{"quantile"} (default),
#'     \code{"pearson"}, and \code{"response"} (raw residuals, y - mu).
#' @param which numeric; if a subset of the plots is required, specify a subset
#'     of the numbers 1:5.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
plot.mvreg <- function(x, type = c("quantile", "pearson", "response"),
                      which = 1:2,
                      ask = prod(graphics::par("mfcol")) < length(which) &&
                        grDevices::dev.interactive(),
                      ...)
{
  if(!is.numeric(which) || any(which < 1) || any(which > 5))
    stop("'which' must be in 1:5")

  ## Reading
  type <- match.arg(type)
  res <- stats::residuals(x, type = type)

  ## Legends
  types <- c("quantile", "pearson", "response")
  Types <- c("Randomized quantile residuals", "Pearson residuals", "Raw response residuals")
  Type <- Types[type == types]

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  ## Plots to shown
  show <- rep(FALSE, 5)
  show[which] <- TRUE

  ## Residuals versus Fitted values
  if (show[1]){
    #y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
    graphics::plot(stats::fitted(x), res,
                   xlab = "Fitted values", ylab = Type, pch = "+", ...)
    graphics::abline(h = 0, col = "gray50", lty = 3)
  }

  ## Residuals versus index observation
  if (show[2]){
    n <- x$nobs
    graphics::plot(1:n, res, xlab = "Index", ylab = Type, pch = "+", ...)
    graphics::abline(h = 0, col= "gray50", lty = 3)
  }

  ## Normal probability plot
  if(show[3]) {
    stats::qqnorm(res,
                  xlab = "Theoretical Quantiles",
                  ylab = "Sample Quantiles",
                  plot.it = TRUE,
                  frame.plot = TRUE, pch =  "+", ...)
    graphics::abline(0, 1, col = "gray50", lty = 2)
  }

  ## Expected frequencies
  if(show[4]) {
    y <- if(is.null(x$y)) stats::model.response(stats::model.frame(x)) else x$y
    obs <- as.numeric(table(y))
    xcoord <- graphics::barplot(table(y), xlab = "Count", ylab = "Frequency",
                                ylim = c(0, max(max(obs), max(x$freq)) + 5), ...)
    graphics::points(xcoord, x$freq, col = 2, type = "b", pch = 16)
    graphics::legend("topright", "Fitted frequencies", col = 2, lty = 1, pch = 16,
                     bty = "n")
  }

  ## ACF of residuals
  if(show[5]) {
    stats::acf(res, main = " ", xlab = "Lags",
               ylab = paste("Sample ACF of", type, "residuals"), ...)
  }

  invisible(x)

}
