#' Predict FREG model
#'
#' Prediction of FREG model
#'
#' @param object FREG model for which predictions are computed
#' @param ... additional arguments relevant for the generic method
#' @param newdata an optional argument. Newdata should be organized as a list. The elements of the list
#' are covariates from FREG model, respectively. No data transformation is needed. Thus, functional covariates are entered in the list \code{newdata}
#' in their raw form. The \code{predict.freg} function will take care of the transformation of such covariates into the functional form of their equivalents from FREG model.
#'
#' @return predictions of dependent variable \code{y}
#' @exportS3Method FREG::predict
#'
#' @examples
#' \donttest{
#' library(fda)
#' y = log10(apply(CanadianWeather$dailyAv[1:334,,2],2,sum))
#' x = CanadianWeather$dailyAv[1:334,,1] # temperature
#' xbasis = create.fourier.basis(c(1,334),5)
#' xfd = smooth.basis(c(1:335),x,xbasis)$fd
#' bbasis = create.fourier.basis(c(1,334),5)
#' latitude = CanadianWeather$coordinates[,1]
#' longitude = CanadianWeather$coordinates[,2]
#' xfdlist = list(xfd, latitude, longitude)
#' cbasis = create.constant.basis(c(1,334))
#' betalist = list(bbasis, cbasis, cbasis)
#' formula = y ~ xfd + latitude + longitude
#' freg.model = freg(formula = formula, betalist = betalist)
#' # Prediction with new data included
#' newdata = list(CanadianWeather$dailyAv[1:365,,1], latitude, longitude)
#' # newdata = list(xfd_1, latitude, longitude) #funct. and scalar variable(s)
#' yhat = predict(freg.model, newdata = newdata)
#' }

predict.freg = function(object, ..., newdata = NULL){

  if(inherits(object, "freg")){

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

    #  names(xfdlist) = object$x.var

  # Extraction of beta coeff from "freg" object
  }

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

  yhat = Z %*% betacoef # Z instead of X

  return(yhat)

  }else stop('Object is not of class "freg"')
}
