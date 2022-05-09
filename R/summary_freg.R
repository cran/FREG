#' Summary of FREG model
#'
#' \code{summary.freg} produce summary of FREG model fitting function \code{freg}.
#'
#' @param object FREG model
#' @param ... additional arguments relevant for the generic method
#' @return call of the function
#' @return beta regression coefficient functions
#' @exportS3Method FREG::summary
#'
#'
#'
summary.freg = function(object, ...){

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


