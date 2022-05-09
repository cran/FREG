#' Cycling dataset
#'
#' A dataset containing information about cycling sessions
#'
#' @format A list of 9 elements. The first element of the list is time whereas the other 8 elements are functional variables.
#' The performance of 216 cyclists is recorded during one hour. Thus, all functional variables are expressed in a form of 216 x 3600 matrix
#' \describe{
#'   \item{SECS}{ 3600 seconds (1 hour)}
#'   \item{KM}{ kilometers}
#'   \item{WATTS}{ power}
#'   \item{CAD}{ cadence}
#'   \item{KPH}{ kilometers per hour}
#'   \item{HR}{ heart rate}
#'   \item{ALT}{ altitude}
#'   \item{SLOPE}{ slope}
#'   \item{TEMP}{ temperature}
#' }
"cycling"
