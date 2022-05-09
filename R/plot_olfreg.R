#' Plot coefficients from OLFREG model
#'
#' @param object OLFREG model
#'
#' @return plot of the beta coefficient regression functions for each intercept and for each variable
#' @importFrom graphics par
#' @importFrom graphics plot
#' @export plot_olfreg
#'
#' @examples
#' \donttest{
#' # cycling dataset
#' library(fda)
#' # creation of ordinal variable from HR variable
#' zoneHR=rep(0,216)
#' zoneHR[which(rowMeans(cycling$HR[,1:1700])<107)]=1
#' zoneHR[which((rowMeans(cycling$HR[,1:1700])<125)&(rowMeans(cycling$HR[,1:1700])>107))]=2
#' zoneHR[which((rowMeans(cycling$HR[,1:1700])<142)&(rowMeans(cycling$HR[,1:1700])>125))]=3
#' zoneHR[which((rowMeans(cycling$HR[,1:1700])<160)&(rowMeans(cycling$HR[,1:1700])>142))]=4
#' zoneHR[which((rowMeans(cycling$HR[,1:1700])>160))]=5
#' # first functional variable - power (WATTS)
#' watts = t(cycling$WATTS[,1:1700])
#' # set up a fourier basis system due to its cycling pattern
#' xbasis = create.fourier.basis(c(1,1700),50) # 50 basis functions for example
#' watts.fd = smooth.basis(c(1:1700),watts,xbasis)$fd
#' zoneHR = as.factor(zoneHR)
#' # additional functional variable - cadence (CAD)
#' cad = t(cycling$CAD[,1:1700])
#' # set up a functional variable for cad
#' xbasis2 = create.bspline.basis(c(1,1700), nbasis = 25, norder = 4)
#' cad.fd = smooth.basis(c(1:1700),cad,xbasis2)$fd
#' formula = zoneHR ~ watts.fd + cad.fd
#' olfreg.model = olfreg(formula = formula)
#' plot_olfreg(olfreg.model)
#'}

plot_olfreg = function(object){

  if(inherits(object, "olfreg")){

  betacoeff = list()
  p = length(object$xfdlist)

  beta = object$coefficients
  alphas = object$alpha
  nalphas = length(object$alpha)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(p,nalphas))


  j = 1
  for (i in 1:p){
    betafd = fdPar(object$betalist[[i]],0,0)$fd # create and extract fd object of each variable and store it
    nbeta = betafd$basis$nbasis # extract the number of basis aka coeffs
    m = j
    j = j + nbeta-1

    betacoef = beta[m:j]
    betafd$coefs = as.matrix(betacoef)

    betacoeff[[i]] = betafd

    for (n in 1:nalphas){
      plt = plot(betacoeff[[i]] + object$alpha[[n]],
                 xlab = "Time", ylab = "Value")
    }

  }

  return(plt)

  }else stop('Model has to be of the class olfreg')
}
