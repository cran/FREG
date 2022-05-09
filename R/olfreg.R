#' Functional ordinal logistic regression model
#'
#' Functional ordinal logistic regression model in which the response variable is a factor variable
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
#' @return \item{call}{ call of the olfreg function}
#' @return \item{x.count}{ number of predictors}
#' @return \item{xfdlist}{ a list of functional data objects. The length of the list is equal to the number of predictors}
#' @return \item{betalist}{ a list of beta regression coefficient functions}
#' @return \item{coefficients}{  estimated beta regression coefficient functions}
#' @return \item{alpha}{  estimated intercepts which represent boudaries of categories of dependent factor variable \code{y}}
#' @return \item{ylev}{ a number of categories of a response variable}
#' @return \item{fitted.values}{  fitted probabilities of a dependent factor variable \code{y}}
#' @return \item{loglik}{  a value of log-likelihood function at optimum}
#' @return \item{grd}{  a vector of gradient values at optimum}
#' @return \item{Hess}{  Hessian matrix at optimum}
#' @return \item{df}{ degrees of freedom}
#' @return \item{AIC}{ Akaike information criterion}
#' @return \item{iteration}{  number of iterations of Fisher Scoring algorithm needed for convergence}
#' @export olfreg
#'
#' @examples
#' # cycling dataset
#' library(fda)
#' # creation of ordinal variable from HR variable
#' zoneHR=rep(0,216)
#' zoneHR[which(rowMeans(cycling$HR[,1:60])<107)]=1
#' zoneHR[which((rowMeans(cycling$HR[,1:60])<125)&(rowMeans(cycling$HR[,1:60])>107))]=2
#' zoneHR[which((rowMeans(cycling$HR[,1:60])<142)&(rowMeans(cycling$HR[,1:60])>125))]=3
#' zoneHR[which((rowMeans(cycling$HR[,1:60])<160)&(rowMeans(cycling$HR[,1:60])>142))]=4
#' zoneHR[which((rowMeans(cycling$HR[,1:60])>160))]=5
#' # first functional variable - power (WATTS)
#' watts = t(cycling$WATTS[,1:60])
#' # set up a fourier basis system due to its cycling pattern
#' xbasis = create.fourier.basis(c(1,60),5) # 5 basis functions for example
#' watts.fd = smooth.basis(c(1:60),watts,xbasis)$fd
#' zoneHR = as.factor(zoneHR)
#' formula = zoneHR ~ watts.fd
#' olfreg.model = olfreg(formula = formula)
#' # additional functional variable - cadence (CAD)
#' cad = t(cycling$CAD[,1:60])
#' # set up a functional variable for cad
#' xbasis2 = create.bspline.basis(c(1,60), nbasis = 5, norder = 4)
#' cad.fd = smooth.basis(c(1:60),cad,xbasis2)$fd
#' formula = zoneHR ~ watts.fd + cad.fd
#' olfreg.model = olfreg(formula = formula)



olfreg = function(formula, betalist = NULL){

  call = match.call()

  # extract y from formula

  y.name = formula[[2]]
  y = get(as.character(y.name)) # search y by name

  y.len = length(y)

  if(inherits(y, c("numeric", "matrix", "array"), FALSE))
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

  p = length(xfdlist) # constant and independent functional variables
  y = as.matrix(y)
  #N = dim(y)[1] # number of observations
  Z  = NULL
  # for any number of covariates
  for (i in 1:p) {
    xfdi       = xfdlist[[i]]
    xcoef      = xfdi$coefs
    xbasis     = xfdi$basis
    bbasis     = betalist[[i]]
    basis.prod  = romberg_alg(xbasis,bbasis)
    Z          = cbind(Z, crossprod(xcoef, basis.prod))

  }

  x = Z
  n = nrow(x)
  xc = ncol(x)
  wt = rep(1, n) # weights
  ind_xc = seq_len(xc)
  if(!is.factor(y)) y = ordered(y) # as.factor(y)

  ylev = levels(y)
  lylev = length(ylev)
  q = length(ylev)-1L
  ind_q = seq_len(lylev-1L)
  y = unclass(y)

  coefs = rep(0, xc)
  logit = function(p) log(p/(1 - p)) # qlogis for quantiles
  taby = tabulate(y)
  alphas = cumsum(taby)[-length(taby)]/n # space
  initial = logit(alphas)

  start = c(coefs, initial)

  res = optimization(x, y, start, loglik, gradient, Hessian)$beta
  beta = res[seq_len(xc)]
  alpha = res[xc + ind_q]
  names(alpha) = paste("y <=", ylev[1:length(ylev)-1])
  names(beta) = paste("X", unlist(bbasis.names), sep = ".")

  # fitted values
  eta = as.vector(x %*% beta)
  cumprob = plogis(matrix(alpha, n, q, byrow = TRUE) - eta)
  fitted.values = as.matrix(t(apply(cumprob, 1, function(x) diff(c(0, x, 1)))))
  colnames(fitted.values) = ylev
  # additional output
  loglik = optimization(x, y, start, loglik, gradient, Hessian)$ll
  grd = optimization(x, y, start, loglik, gradient, Hessian)$grd
  hessian = optimization(x, y, start, loglik, gradient, Hessian)$hessian
  iteration =  optimization(x, y, start, loglik, gradient, Hessian)$iter

  # calculate degrees of freedom and AIC
  loglik = -loglik # return the actual value of log-likelihood and not the negative one
  df = length(alpha) + length(beta)
  AIC = -2*loglik + 2*df

  instance = list()

  instance$call = call
  instance$no.var = x.count
  instance$xfdlist = xfdlist
  instance$betalist = betalist
  instance$coefficients = beta
  instance$alpha = alpha
  instance$ylev = ylev
  instance$fitted.values = fitted.values
  instance$loglik = loglik
  instance$grd = grd
  instance$hessian = hessian
  instance$df = df
  instance$AIC = AIC
  instance$iteration = iteration

  class(instance) = "olfreg"

  return(instance)

}
