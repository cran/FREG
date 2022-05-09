#' Summary of LFREG model
#'
#' \code{summary.lfreg} produce summary of LFREG model fitting function \code{lfreg}.
#'
#' @param object LFREG model
#' @param ... additional arguments relevant for the generic method
#' @return call of the function
#' @return beta regression coefficient functions
#' @exportS3Method FREG::summary
#'

summary.lfreg = function(object, ...){

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

}
