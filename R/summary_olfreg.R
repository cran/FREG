#' Summary of OLFREG model
#'
#' \code{summary.olfreg} produce summary of OLFREG model fitting function \code{olfreg}.
#'
#' @param object OLFREG model
#' @param ... additional arguments relevant for the generic method
#' @return call of the function
#' @return alpha regression coefficients and beta regression coefficient functions
#' @exportS3Method FREG::summary
#'

summary.olfreg = function(object, ...){
  name_object <- deparse(substitute(object))
  cat("----------------------------------------------------\n")
  cat("               Summary of :", name_object, "\n")
  cat("----------------------------------------------------\n")

  cat("Call:\n")
  print(object$call)
  cat("\n\n")

  cat("Coefficients :\n")
  cat("--------------  \n")
  print(object$coefficients)
  cat("\n\n")

  cat("Intercepts :\n")
  cat("--------------  \n")
  print(object$alpha)
  cat("\n\n")

}
