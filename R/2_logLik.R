# Link function ----------------------------------------------------------------
make.link.ldrm <- function(link){

  switch(link,

         identity = {
           fun <- function(theta) theta
           inv <- function(eta) eta
           deriv. <- function(theta) rep.int(1, length(theta))
           deriv.. <- function(theta) rep.int(0, length(theta))
           valideta <- function(eta) TRUE
         },

         log = {
           fun <- function(theta) log(theta)
           inv <- function(eta) pmax(exp(eta), .Machine$double.eps)
           deriv. <- function(theta) 1 / theta
           deriv.. <- function(theta) -1 / (theta^2)
           valideta <- function(eta) TRUE
         },

         sqrt = {
           fun <- function(theta) sqrt(theta)
           inv <- function(eta) eta^2
           deriv. <- function(theta) 1 / (2 * sqrt(theta))
           deriv.. <- function(theta) -1 / (4 * (theta ^ (3 / 2)))
           valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
         },

         stop(gettextf("link %s not available", sQuote(link)), domain = NA))

  environment(fun) <- environment(inv) <- environment(deriv.) <-
    environment(deriv..) <- environment(valideta) <- asNamespace("stats")

  structure(list(fun = fun, inv = inv, deriv. = deriv.,
                 deriv.. = deriv.., valideta = valideta,
                 name = link))
}


# Log-likelihood  --------------------------------------------------------------
ll <- function(par, y, X, Z = NULL, count = ripo){

  # Count distribution
  family <- as.count(count)
  ddist <- get(paste0("d", family$abb))

  # Mean submodel
  X <- as.matrix(X); p <- NCOL(X)
  mu.link <- family$mu.link
  mu <- c(make.link.ldrm(mu.link)$inv(X%*%par[1:p]))

  # Variance submodel
  sigma2 <- NULL
  if(family$npar > 1){
    Z <- as.matrix(Z); k <- NCOL(Z)
    sigma2.link <- family$sigma2.link
    sigma2 <- c(make.link.ldrm(sigma2.link)$inv(Z%*%par[1:k + p]))
  }

  # Extra parameters (if any)
  nu <- NULL
  if(family$npar > 2){
    nu <- par[p + k + family$npar - 2]
  }

  if(any(!is.finite(mu)) | any(!is.finite(sigma2))){
    NaN
  }else {
    sum(ddist(y, mu = mu, sigma2 = sigma2, nu = nu, log.p = TRUE))
  }
}

# Score function ---------------------------------------------------------------
Uscore <- function(par, y, X, Z = NULL, count = ripo){

  # Count distribution
  family <- as.count(count)

  if(!family$has_dl){
    U <- try(numDeriv::grad(ll, par, y = y, X = X, Z = Z, family = family), silent = TRUE)
    error <- unique(grepl("Error", U))
    if (error) return(NULL) else return(U)
  }

  # Extra parameters
  nu <- NULL
  Un <- NULL
  if(family$npar > 2){
    nu <- par[p + k + family$npar - 2]
    Un <- family$dl.dnu(mu = mu, sigma2 = sigma2, nu = nu)
  }

  # Mean submodel
  X <- as.matrix(X); p <- NCOL(X)
  mu.link <- family$mu.link
  mu <- c(make.link.ldrm(mu.link)$inv(X%*%par[1:p]))

  g1. <- make.link.ldrm(mu.link)$deriv.
  D1 <- diag(as.numeric(1/g1.(mu)))

  # Variance submodel
  sigma2 <- NULL
  Ug <- NULL
  if(family$npar > 1){
    Z <- as.matrix(Z); k <- NCOL(Z)
    sigma2.link <- family$sigma2.link
    sigma2 <- c(make.link.ldrm(sigma2.link)$inv(Z%*%par[1:k + p]))

    g2. <- make.link.ldrm(sigma2.link)$deriv.
    D2 <- diag(as.numeric(1/g2.(sigma2)))
    u2 <- family$dl.dsigma(y, mu = mu, sigma2 = sigma2, nu = nu)

    Ug <- t(Z)%*%D2%*%u2
  }

  u1 <- family$dl.dmu(y, mu = mu, sigma2 = sigma2, nu = nu)

  Ub <- t(X)%*%D1%*%u1;

  c(Ub, Ug, Un)
}


# Hessian matrix  --------------------------------------------------------------
Jhessian <- function(par, y, X, Z = NULL, count = ripo){

  # Count distribution
  family <- as.count(count)

  if(!family$has_dl){
    J <- try(numDeriv::hessian(ll, par, y = y, X = X, Z = Z, family = family), silent = TRUE)
    error <- unique(grepl("Error", J))
    if (error) return(NA) else return(J)
  }

  # Extra parameters
  nu <- NULL
  if(family$npar > 2){
    nextra_par <- family$npar - 2
    nu <- par[p + k + nextra_par]
  }

  # Mean submodel
  X <- as.matrix(X)
  p <- NCOL(X)
  mu.link <- family$mu.link
  mu <- c(make.link.ldrm(mu.link)$inv(X%*%par[1:p]))

  g1. <- make.link.ldrm(mu.link)$deriv.
  g1.. <- make.link.ldrm(mu.link)$deriv..
  D1 <- diag(as.numeric(1/g1.(mu)))

  # Variance submodel
  sigma2 <- NULL
  Ug <- NULL
  k <- 0

  J <- matrix(NA, p + k, p + k)

  if(family$npar > 1){
    Z <- as.matrix(Z)
    k <- NCOL(Z)

    J <- matrix(NA, p + k, p + k)

    sigma2.link <- family$sigma2.link
    sigma2 <- c(make.link.ldrm(sigma2.link)$inv(Z%*%par[1:k + p]))

    g2. <- make.link.ldrm(sigma2.link)$deriv.
    g2.. <- make.link.ldrm(sigma2.link)$deriv..
    D2 <- diag(as.numeric(1/g2.(sigma2)))

    W2 <- diag(c(family$dl2.dmudsigma(y, mu = mu, sigma2 = sigma2, nu = nu)))

    W3 <- diag(c(family$dl2.dsigma2(y, mu = mu, sigma2 = sigma2, nu = nu) / g2.(sigma2) +
                   family$dl.dsigma(y, mu = mu, sigma2 = sigma2, nu = nu) * (g2..(sigma2) / (g2.(sigma2)^2))))

    J[1:p,(p + 1):(p + k)] <- t(X)%*%D1%*%W2%*%D2%*%Z
    J[(p + 1):(p + k),(1:p)] <- t(J[1:p,(p + 1):(p + k)])
    J[(p + 1):(p + k),(p + 1):(p + k)] <- t(Z)%*%W3%*%D2%*%Z
  }


  W1 <- diag(c(family$dl2.dmu2(y, mu, sigma2) / g1.(mu) +
                 family$dl.dmu(y, mu, sigma2) * (g1..(mu) / (g1.(mu)^2))))

  J[1:p,1:p] <- t(X)%*%W1%*%D1%*%X

  if(family$npar > 2){
    Jbnu <- family$dl2.dmudnu(mu, sigma2, nu)
    Jgnu <- family$dl2.dsigmadnu(mu, sigma2, nu)
    Jnu  <- family$dl2.dnu2(mu, sigma2, nu)

    J <- cbind(rbind(NA, J), NA)
    J[p + k + nextra_par, ] <- J[, p + k + nextra_par] <- c(Jbnu, Jgnu, Jnu)
  }

  J

}

# --------------------------------------------------------------------------------------------------
Jaux <- function(par, y, X, Z = NULL, count = ripo){

  # Count distribution
  family <- as.count(count)

  if(!family$has_dl){
    return(matrix(NA, length(par), length(par)))
  }

  if(family$npar == 1){

    X <- as.matrix(X)
    p <- NCOL(X)
    mu.link <- family$mu.link
    mu <- c(make.link.ldrm(mu.link)$inv(X%*%par[1:p]))

    g1. <- make.link.ldrm(mu.link)$deriv.

    u1 <- family$dl.dmu(y, mu = mu)
    D11 <- diag(as.numeric(u1 / g1.(mu))^2)

    return(t(X)%*%D11%*%X)

  }else{

    # Extra parameters
    nu <- NULL
    if(family$npar > 2){
      nextra_par <- family$npar - 2
      nu <- par[p + k + nextra_par]
    }

    X <- as.matrix(X); Z <- as.matrix(Z)
    p <- NCOL(X); k <- NCOL(Z)

    mu.link <- family$mu.link
    sigma2.link <- family$sigma2.link

    mu <- c(make.link.ldrm(mu.link)$inv(X%*%par[1:p]))
    sigma2 <- c(make.link.ldrm(sigma2.link)$inv(Z%*%par[1:k + p]))

    g1. <- make.link.ldrm(mu.link)$deriv.
    g2. <- make.link.ldrm(sigma2.link)$deriv.

    u1 <- family$dl.dmu(y, mu = mu, sigma2 = sigma2, nu = nu)
    u2 <- family$dl.dsigma(y, mu = mu, sigma2 = sigma2, nu = nu)

    D11 <- diag(as.numeric(u1 / g1.(mu))^2)
    D22 <- diag(as.numeric(u2 / g2.(sigma2))^2)
    D12 <- diag(as.numeric(u1 * u2 / (g1.(mu) * g2.(sigma2))))

    J <- matrix(NA, p + k, p + k)

    J[1:p,1:p] <- t(X)%*%D11%*%X
    J[1:p,(p + 1):(p + k)] <- t(X)%*%D12%*%Z
    J[(p + 1):(p + k),(1:p)] <- t(Z)%*%D12%*%X
    J[(p + 1):(p + k),(p + 1):(p + k)] <- t(Z)%*%D22%*%Z

    return(J)
  }



}

