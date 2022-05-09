#' Plot coefficients from LFREG model
#'
#' @param object LFREG model
#'
#' @return plot of the beta coefficient regression functions for each variable
#' @importFrom graphics plot
#' @export plot_lfreg
#'
#' @examples
#' library(fda)
#' precipitation_data = CanadianWeather$daily[1:334,,"Precipitation.mm"]
#' annualprec = apply(precipitation_data,2,sum) # without December
#' y = ifelse(annualprec<mean(annualprec), 0, 1)
#' y = as.factor(y)
#' x = CanadianWeather$daily[1:334,,"Temperature.C"]
#' xbasis = create.fourier.basis(c(1,334),5) # 5 basis functions
#' xfd = smooth.basis(c(1:334),x,xbasis)$fd
#' bbasis = create.fourier.basis(c(0,334),5)
#' betalist = list(bbasis)
#' formula = y ~ xfd
#' lfreg.model = lfreg(formula, betalist = betalist)
#' plot_lfreg(lfreg.model)
#'

plot_lfreg = function(object){

  if(inherits(object, "lfreg")){

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

  }else stop('Model has to be of the class lfreg')
}
