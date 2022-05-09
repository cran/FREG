#' Functional logistic regression model
#'
#' Functional logistic regression model in which the response variable is a factor variable
#' whereas the independent variables are functional variables. Independent variables could also be scalar variables.
#'
#' @param formula  a formula expression of the form \code{response ~ predictors}. On the left side of the formula, \code{y} is a factor variable whereas on the right side,
#' \code{X} can be either functional data object of class \code{fd} or a scalar variable of class \code{numeric}. The length of a scalar variable must equal the length
#' of a response variable. Similarly, the number of observations of a functional covariate must equal the length of a response variable.
#' @param betalist  an optional argument. A list which contains beta regression coefficient functions for independent variables.
#' If betalist is not provided, the number of estimated beta regression coefficient functions for one functional covariate would equal the number of basis functions used to represent that functional covariate.
#' For a scalar variable, beta regression coefficient function is also a functional object whose basis is constant.
#' Needless to say, for a scalar variable, there will be one beta regression coefficient.
#'
#' @return \item{call}{ call of the lfreg function}
#' @return \item{x.count}{ number of predictors}
#' @return \item{xfdlist}{ a list of functional data objects. The length of the list is equal to the number of predictors}
#' @return \item{betalist}{ a list of beta regression coefficient functions}
#' @return \item{coefficients}{  estimated beta regression coefficient functions}
#' @return \item{fitted.values}{  predicted values of a response variable \code{y}}
#' @return \item{loglik}{ a value of log-likelihood function at optimum}
#' @return \item{df}{ degrees of freedom}
#' @return \item{AIC}{ Akaike information criterion}
#' @return \item{iteration}{  number of iterations needed for convergence criterion to be met}
#' @export lfreg
#' @import fda
#' @import stats
#'
#' @examples
#'
#' library(fda)
#' precipitation_data = CanadianWeather$daily[,,"Precipitation.mm"]
#' annualprec = apply(precipitation_data,2,sum)
#' y = ifelse(annualprec<mean(annualprec), 0, 1)
#' y = as.factor(y)
#' x = CanadianWeather$daily[,,"Temperature.C"]
#' xbasis = create.fourier.basis(c(1,365),5) # 5 basis functions
#' # smoothing of the data and extraction of functional data object
#' xfd = smooth.basis(c(1:365),x,xbasis)$fd
#' # bbasis and betalist are optional arguments
#' bbasis = create.fourier.basis(c(0,365),3) # 3 bf
#' betalist = list(bbasis)
#' formula = y ~ xfd
#' lfreg.model = lfreg(formula, betalist = betalist)
#' # add scalar variables
#' latitude = CanadianWeather$coordinates[,1]
#' longitude = CanadianWeather$coordinates[,2]
#' # cbasis and betalist are optional arguments
#' cbasis = create.constant.basis(c(1,365))
#' betalist = list(bbasis, cbasis, cbasis)
#' formula = y ~ xfd + latitude + longitude
#' lfreg.model = lfreg(formula, betalist = betalist)


lfreg = function(formula, betalist = NULL){

  call = match.call()

  # extract y from formula

  y.name = formula[[2]]
  y = get(as.character(y.name)) # search y by name

  y.len = length(y)

  if (inherits(y, c("numeric", "matrix", "array"), FALSE))
    stop("Y has to be factor")

  # extract independent variables from formula

  x.var = all.vars(formula)[-1]
  x.count = length(x.var)

  xfdlist = vector('list', length = x.count) # stock them in the list
  names(xfdlist) = x.var

  type = c()
  nbasis = c()

  df = lapply(x.var, get)
  no = which(lapply(df, class)=="fd")
  if(length(no)>1){
    range = get(x.var[no[1]])$basis$range # take range from the first fd
  }else range = get(x.var[no])$basis$range

  # if (inherits(get(x.var), what = "fd")){
  #   x.fun = get(x.var)
  #   range = x.fun$basis$rangeval
  # }else stop("Please enter a functional covariate")

  bbasis.names = vector('list', length = x.count) # to stock beta basis names

  for(i in 1:x.count){

    x = get(x.var[i])

    if(inherits(x, what = "fd")){
      type[i] = x$basis$type
      nbasis[i] = x$basis$nbasis
      x.len = dim(x$coefs)[2]
      #range = x$basis$rangeval
    }else if(inherits(x, what = "numeric")){
      cbasis = create.constant.basis(rangeval = range)
      x.len = length(x)
      x = fd(matrix(x,1,y.len),cbasis)
    }

    if(x.len != y.len)
      stop('The number of observations of ',x.var[i], ' is ', x.len,
           ' and is not equal to the number of observations of y ', y.len)

    if(!class(x) %in% c("fd", "numeric"))
      stop('Variable ', x.var[i], ' has to be either fd or numeric')

    xfdlist[[i]] = x

    # create betalist

    if(is.null(betalist)){
      betalist = vector('list', length = x.count)
    }
    if(is.null(betalist[[i]])){
      if(class(x) %in% "fd"){
        bbasis = with(x, fd(basis = basis, fdnames = fdnames))$basis
      }else if(class(x) %in% "numeric"){
        bbasis = create.constant.basis(rangeval = range)
      }

      betalist[[i]] = bbasis

    }else if(length(betalist) != length(xfdlist)){

      stop('length(betalist) is ', length(betalist),
           ' but it  must be equal to the number of independent variables ', length(xfdlist))

      betaclass = sapply(betalist, class)
      wrong = which(betaclass != 'basisfd')

      if(length(wrong) > 0)
        stop('All components of betalist must have class basisfd')
    }
    bbasis.names[[i]] = betalist[[i]]$names
  }

  # estimation

  p = length(xfdlist) # independent functional variables
  y = as.matrix(as.numeric(as.character(y)))
  N = dim(y)[1] # number of observations
  Z  = NULL
  # for any number of covariates
  for (i in 1:p) {
    xfdi       = xfdlist[[i]]
    xcoef      = xfdi$coefs
    xbasis     = xfdi$basis
    bbasis     = betalist[[i]]
    basis_prod  = romberg_alg(xbasis,bbasis) # call of a function romberg_alg
    Z          = cbind(Z, crossprod(xcoef, basis_prod))

  }

  change <- Inf
  p = ncol(Z)
  beta_old <- rep(0,p)
  #intercept = 1 #without intercept when estimating with scalar variable...
  iter = 25
  tolerance = 1e-8


  for(i in 1:iter){

    eta = Z %*% beta_old
    mu = function(x) 1 /(1 + exp(-x)) # inverse link
    y.hat = mu(eta)
    wghts = y.hat * (1 - y.hat)
    z = eta + (y - y.hat) / wghts

    beta_new = lm(z ~ Z-1, weights = wghts)$coef
    change = sqrt(sum((beta_new - beta_old)^2))
    beta_old = beta_new

    if(change < tolerance) break

  }

  if (i==max(iter)) warning("The algorithm did not converge.")

  coefficients = as.matrix(beta_new)
  rownames(coefficients) = paste("X", unlist(bbasis.names), sep = ".")

  bbasis = betalist


  loglik = t(y)%*%eta-sum(log(1+exp(eta)))
  df = length(coefficients)
  AIC = -2*loglik + 2*df

  instance = list()

  instance$call = call
  instance$no.var = x.count
  instance$xfdlist = xfdlist
  instance$betalist = betalist
  instance$coefficients = coefficients
  #instance$bbasis = bbasis
  instance$fitted.values = y.hat
  instance$loglik = loglik
  instance$df = df
  instance$AIC = AIC
  instance$iteration = i

  class(instance) = "lfreg"

  return(instance)

}
