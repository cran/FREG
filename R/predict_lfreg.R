#' Predict LFREG model
#'
#' Prediction of LFREG model
#'
#' @param object LFREG model for which predictions are computed
#' @param ... additional arguments relevant for the generic method
#' @param newdata an optional argument. Newdata should be organized as a list. The elements of the list
#' are covariates from LFREG model, respectively. No data transformation is needed. Thus, functional covariates are entered in the list \code{newdata}
#' in their raw form. The \code{predict.lfreg} function will take care of the transformation of such covariates into the functional form of their equivalents from LFREG model.
#' @param type c("probabilities", "labels")
#'
#' @return predictions of dependent variable \code{y}
#' @exportS3Method FREG::predict
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
#' # Prediction on new data
#' newdata = list(CanadianWeather$dailyAv[1:365,,1])
#' # newdata = list(xfd_1, latitude, longitude)
#' yhat = predict(lfreg.model, newdata = newdata, type = "labels")
#'

predict.lfreg = function(object, ..., newdata = NULL, type = c("probabilities", "labels")){

  # Check

  if(inherits(object, "lfreg")){

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

  # Extraction of beta coeff from "lfreg" object

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

  eta = Z %*% betacoef  # Z instead of X
  mu = function(x) 1 /(1 + exp(-x)) # inverse link
  yhat = mu(eta)

  if(missing(type) || type == "probabilities"){
    return(yhat)
  }else if(type == "labels"){
    yhat.labels = ifelse(yhat > 0.5, 1, 0)
    return(yhat.labels)

  }
  }else stop('Object is not of class "lfreg"')

}
