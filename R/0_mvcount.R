#' @name count
#'
#' @title Methods for \code{"count"} objects
#'
#' @param x,object an object of class \code{"count"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The current available count distributions can be seen below.
#'   \tabular{llc}{
#'  \bold{Distribution}  \tab \bold{Abbreviation} \tab \bold{N. of parameters}\cr
#'  Bell-Touchard  \tab \code{"\link{beto}"}      \tab  2  \cr
#'  BerG  \tab \code{"\link{berg}"}      \tab  2  \cr
#'  Double Poisson  \tab \code{"\link{dbpo}"}      \tab  2  \cr
#'  Generalized Poisson \tab \code{"\link{gpo}"}      \tab  2  \cr
#'  Negative binomial  \tab \code{"\link{nb}"}      \tab  2  \cr
#'  Poisson  \tab \code{"\link{po}"}      \tab  1  \cr
#'  Inflated Bernoulli \tab \code{"\link{ribe}"}      \tab  2  \cr
#'  Inflated geometric \tab \code{"\link{rige}"}      \tab  2  \cr
#'  Inflated Poisson \tab \code{"\link{ripo}"}      \tab  2  \cr
#'  Zero inflated geometric \tab \code{"\link{zig}"}      \tab  2  \cr
#'  Zero inflated Poisson \tab \code{"\link{zip}"}      \tab  2  \cr
#'  }
#'
#'
#' @references
#'
#'  Bourguignon, M., and de Medeiros, R. M. (2022). A simple and useful regression model
#'     for fitting count data. \emph{TEST}, \bold{31}, 790-827.
#'
#'  Castellares, F., Lemonte, A. J., and Moreno–Arenas, G. (2020). On the two-parameter
#'     Bell–Touchard discrete distribution. \emph{Communications in Statistics-Theory and Methods},
#'     \bold{49}, 4834-4852.
#'
#' Consul, P. C., and Jain, G. C. (1973). A generalization of the Poisson distribution.
#'     \emph{Technometrics}, \bold{15}, 791-799.
#'
#'  Joe, H., and Zhu, R. (2005). Generalized Poisson distribution: the property of
#'     mixture of Poisson and comparison with negative binomial distribution.
#'     \emph{Biometrical Journal: Journal of Mathematical Methods in Biosciences}, \bold{47}, 219-229.
#'
#'  Kokonendji, C. C., Medeiros, R. M. R., and Bourguignon, M. (2024). Mean and variance
#'     count regression models based on reparameterized distributions. Submitted.
#'
#'  Rodríguez-Avi, J., and Olmo-Jiménez, M. J. (2017).
#'    A regression model for overdispersed data without too many zeros.
#'    \emph{Statistical Papers}, \bold{58}, 749-773.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' count("beto")
#' count("berg")
#' count("rige")
#' count("ripo")
#' count("zig")
#' count("zip")
NULL

# Class definition
#' @rdname count
#' @export
count <- function(object, ...) UseMethod("count")

# Default
#' @rdname count
#' @export
count.default <- function(object, ...) {
  cl <- data.class(object)[1]

  return(switch(cl,
                count = object,
                "function" = count(object()),
                character = count(get(object)),
                name = count(eval(object)),
                call = count(eval(object)),
                stop("The object argument is invalid")
  ))
}

# Transformation to bcs class
#' @rdname count
#' @export
as.count <- function(object) {
  if (inherits(object, "count")) object else count(object)
}

# Print method
#' @rdname count
#' @export
print.count <- function(x, ...) {
  cat(
    "Mean and Variance Parametrized Count Distribution\n",
    "\nName:", x$name,
    "\nAbbreviation:", x$abb,
    "\nNumber of parameters:", x$npar, "\n"
  )
}







