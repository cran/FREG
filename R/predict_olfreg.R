#' Predict OLFREG model
#'
#' Prediction of OLFREG model
#'
#' @param object OLFREG model for which predictions are computed
#' @param ... additional arguments relevant for the generic method
#' @param newdata an optional argument. Newdata should be organized as a list. The elements of the list
#' are covariates from OLFREG model, respectively. No data transformation is needed. Thus, functional covariates are entered in the list \code{newdata}
#' in their raw form. The \code{predict.olfreg} function will take care of the transformation of such covariates into the functional form of their equivalents from OLFREG model.
#' @param type c("probabilities", "labels")
#'
#' @return predictions of dependent variable \code{y}
#' @exportS3Method FREG::predict
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
#'
#' # Predict with new data included
#' watts_new = t(cycling$WATTS[,101:1800])
#' cad_new = t(cycling$CAD[,101:1800])
#' newdata = list(watts_new, cad_new) # could also be fd var instead of raw data
#' yhat = predict(olfreg.model, newdata = newdata, type = "labels")
#'}

predict.olfreg = function(object, ..., newdata = NULL, type = c("probabilities", "labels")){

  # Check

  if(inherits(object, "olfreg")){

  if (inherits(newdata, c("numeric", "factor", "matrix", "array", "fd"), FALSE))
    stop('Newdata should be organized as list')

  x.count = object$no.var

  if(length(newdata) != x.count)
    stop('There is a mismatch in number of predicted variables and initial covariates')

  xfdlist = vector('list', length = x.count) # stock them in the list

  for(i in 1:x.count){

    if(inherits(newdata[[i]], what = 'fd')){
      range2 = newdata[[i]]$basis$rangeval[2]
    }else if(inherits(newdata[[i]], what = 'matrix')){
      range2 = nrow(newdata[[i]])
    }
  }

  for(i in 1:x.count){

    no_basis = list()
    type_basis = list()
    if(inherits(object$xfdlist[[i]], what = 'fd')){
      type_basis[[i]] = object$xfdlist[[i]]$basis$type
      no_basis[[i]] = object$xfdlist[[i]]$basis$nbasis
    }

    if(inherits(newdata[[i]], what = 'fd')){
      x = newdata[[i]]
    }else if(inherits(newdata[[i]], what = 'matrix')){
      if(type_basis[[i]]=="fourier"){
        xbasis[[i]] = create.fourier.basis(c(1,range2),no_basis[[i]])
      }else if(type_basis[[i]]=="bspline"){
        xbasis[[i]] = create.bspline.basis(c(1,range2), nbasis = no_basis[[i]], norder = 4)
      }
      x = smooth.basis(c(1:range2),newdata[[i]],xbasis[[i]])$fd
    }else if(inherits(newdata[[i]], what = 'numeric')){
      cbasis = create.constant.basis(c(1,range2))
      x.len = length(newdata[[i]])
      x = fd(matrix(newdata[[i]],1,x.len),cbasis)
    }else stop("The variable has to be of class fd, matrix or numeric")

    xfdlist[[i]] = x
  }

  names(xfdlist) = object$x.var

  # Extraction of alpha and beta coeff from "olfreg" object

  alpha = object$alpha
  betacoef = object$coefficients
  betalist = object$betalist

  # Calculation of yhat

  p = length(xfdlist) # constant and independent functional variables
  Z  = NULL
  # for any number of covariates
  for (i in 1:p) {
    xfdi       = xfdlist[[i]]
    xcoef      = xfdi$coefs
    xbasis     = xfdi$basis
    bbasis     = betalist[[i]]
    basis.prod  = romberg_alg(xbasis,bbasis) # same as in the estimation
    Z          = cbind(Z, crossprod(xcoef, basis.prod)) # xcoef is new

  }

  n = nrow(Z)
  ylev = object$ylev
  q = length(ylev)-1L

  eta = as.vector(Z %*% betacoef) # Z instead of X
  cumprob = plogis(matrix(alpha, n, q, byrow = TRUE) - eta)
  yhat = as.matrix(t(apply(cumprob, 1, function(x) diff(c(0, x, 1)))))
  colnames(yhat) = ylev

  if(missing(type) || type == "probabilities"){
    return(yhat)
  }else if(type == "labels"){
    yhat.labels = apply(yhat,1,which.max)
    return(yhat.labels)
  }

  }else stop('Object is not of class "olfreg"')

}
