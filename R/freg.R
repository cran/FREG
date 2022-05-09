#' Functional linear regression model
#'
#' Functional linear regression model in which the response variable is a scalar variable
#' whereas the independent variables are functional variables. Independent variables could also be scalar variables.
#'
#' @param formula  a formula expression of the form \code{response ~ predictors}. On the left side of the formula, \code{y} is a numeric variable whereas on the right side,
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
#' @export freg
#'
#' @examples
#' library(fda)
#' y = log10(apply(daily$precav,2,sum))
#' x = daily$tempav
#' xbasis = create.fourier.basis(c(1,365),5) # 5 basis functions
#' # smoothing of the data and extraction of functional data object
#' xfd=smooth.basis(c(1:365),x,xbasis)$fd
#' formula = y ~ xfd
#' # betalist is an optional argument
#' bbasis = create.fourier.basis(c(1,365),5) # 5 basis functions
#' betalist = list(bbasis)
#' freg.model = freg(formula = formula, betalist = betalist)
#'
#' # Functional variable and two scalar variables
#' latitude = CanadianWeather$coordinates[,1]
#' longitude = CanadianWeather$coordinates[,2]
#' xfdlist = list(xfd, latitude, longitude)
#' cbasis = create.constant.basis(c(1,365))
#' betalist = list(bbasis, cbasis, cbasis)
#' formula = y ~ xfd + latitude + longitude
#' freg.model = freg(formula = formula, betalist = betalist)
#' print(freg.model$coefficients)
#'

freg = function(formula, betalist = NULL){

  # extract y from formula

  call = match.call()

  y.name = formula[[2]]
  y = get(as.character(y.name)) # search y by name

  y.len = length(y)

  if (inherits(y, c("factor", "matrix", "array"), FALSE))
    stop("Y has to be numeric")

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

  bbasis.names = vector('list', length = x.count) # to stock beta basis

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
      #x.len = dim(x$coefs)[2]
    }

    if(x.len != y.len)
      stop('The number of observations of ',x.var[i], ' is ', x.len,
           ' and is not equal to the number of observations of y ', y.len)

    if(!class(x) %in% c("fd", "numeric"))
      stop('Variable ', x.var[i], ' has to be either fd or numeric')

    # for freg include intercept
    #cbasis = create.constant.basis(c(0, x$basis$rangeval[2])) # what if first var is scalar...
    #xfdlist[[1]] = fd(matrix(1,1,y.len),cbasis)
    #xfdlist[[i+1]] = x

    xfdlist[[i]] = x # continue without intercept..

    # create betalist

    if(is.null(betalist)){
      betalist = vector('list', length = x.count)
      #bbasis.names = vector('list', length = x.count)
    }

    if(is.null(betalist[[i]])){
      if(class(x) %in% "fd"){
        bbasis = with(x, fd(basis = basis, fdnames = fdnames))$basis
        #betalist[[i]] = bbasis
      }else if(class(x) %in% "numeric"){
        bbasis = create.constant.basis(rangeval = range)
        #betalist[[i]] = bbasis
      }

      betalist[[i]] = bbasis
      #bbasis.names[[i]] = betalist[[i]]$names

      # }else if(!is.null(betalist[[i]])){
      #     bbasis.names[[i]] = betalist[[i]]$names

    }else if(length(betalist) != length(xfdlist)){

      stop('length(betalist) is ', length(betalist),
           ' but it must be equal to the number of independent variables ',length(xfdlist))

      betaclass = sapply(betalist, class)
      wrong = which(betaclass != 'basisfd')

      if(length(wrong) > 0)
        stop('All components of betalist must have class basisfd')
    }

    bbasis.names[[i]] = betalist[[i]]$names
  }

  # estimation

  p = length(xfdlist) # constant and independent functional variables
  y = as.matrix(y)
  N = dim(y)[1] # number of observations
  Z  = NULL
  betacoef = NULL
  # for any number of covariates
  for (i in 1:p) {
    xfdi       = xfdlist[[i]]
    xcoef      = xfdi$coefs
    xbasis     = xfdi$basis
    bbasis     = betalist[[i]]
    basis.prod  = romberg_alg(xbasis,bbasis) # call of a function romberg_alg
    Z          = cbind(Z, crossprod(xcoef, basis.prod))

    #  Beta coefficients

    Z_prim = as.matrix(t(Z) %*% Z)
    Z_prim.inv  = solve(Z_prim)
    DZ = t(Z) %*% y
    betacoef = Z_prim.inv %*% DZ

  }

  rownames(betacoef) = paste("X", unlist(bbasis.names), sep = ".")

  bbasis = betalist

  instance = list()

  instance$call = call
  instance$no.var = x.count
  instance$xfdlist = xfdlist
  instance$betalist = betalist
  instance$coefficients = betacoef


  class(instance) = "freg"

  return(instance)

}
