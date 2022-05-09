#' Plot coefficients from FREG model
#'
#' @param object FREG model
#'
#' @return plot of the beta coefficient regression functions for each variable
#' @importFrom graphics par
#' @importFrom graphics plot
#' @export plot_freg
#'
#' @examples
#' \donttest{
#' library(fda)
#' y = log10(apply(CanadianWeather$dailyAv[1:335,,2],2,sum))
#' x = CanadianWeather$dailyAv[1:335,,1] # temperature
#' xbasis = create.fourier.basis(c(1,335),5)
#' xfd = smooth.basis(c(1:335),x,xbasis)$fd
#' bbasis = create.fourier.basis(c(1,335),5)
#' betalist = list(bbasis)
#' formula = y ~ xfd
#' freg.model = freg(formula = formula, betalist = betalist)
#' plot_freg(freg.model)
#'}
#'
plot_freg = function(object){

  if(inherits(object, "freg")){

  betacoeff = list()
  p = length(object$xfdlist)

  beta = object$coefficients

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(p,1))

  j = 1
  for (i in 1:p){
    betafd = fdPar(object$betalist[[i]],0,0)$fd # create and extract fd object of each variable and store it
    nbeta = betafd$basis$nbasis # extract the number of basis aka coeffs
    m = j
    j = j + nbeta-1

    betacoef = beta[m:j]
    betafd$coefs = as.matrix(betacoef)

    betacoeff[[i]] = betafd

    plt = plot(betacoeff[[i]],
               xlab = "Time", ylab = "Value")

  }

  return(plt)

  }else stop('Model has to be of the class freg')
}
