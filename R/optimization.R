#' Fisher Scoring algorithm
#'
#' Optimization algorithm for the estimation of beta regression coefficient functions and intercepts
#'
#' @param x  a design matrix which is a product of inner product of basis functions and basis coefficients of functional covariate \code{X}
#' @param y  a response variable of class \code{factor}
#' @param beta  initial values for beta regression coefficients and intercepts
#' @param loglik  log-likelihood function
#' @param gradient  function for the estimation of first derivative of log-likelihood function - gradient
#' @param Hessian  function for the estimation of second derivative of log-likelihood function - Hessian
#'
#' @return \item{beta}{  a vector with estimated beta regression coefficients and intercepts}
#' @return \item{ll}{  a value of the log-likelihood function at the estimated optimum}
#' @return \item{grd}{  a vector of gradient values at the estimated optimum}
#' @return \item{hessian}{  Hessian matrix at the estimated optimum}
#' @export
#'
optimization = function(x, y, beta, loglik, gradient, Hessian){

  epsilon = 1e-7
  step = 0
  iter = 0
  condition = NULL
  max.iter = 25

  while(all(condition > epsilon)){

    ll = loglik(x,y,beta)
    grd = gradient(x,y,beta)$grd
    Hess = Hessian(x,y,beta)

    ch = chol(Hess)
    step = c(backsolve(ch, backsolve(ch, grd, transpose=TRUE))) # to ensure that Hess is not negative definite
    #step_check = solve(Hess, grd)

    stepFactor = 1
    beta_old = beta
    beta = beta_old - stepFactor*step

    llcheck = loglik(x,y,beta)

    stephalf = (llcheck > ll) # if thresholds are not incresing activate step-halving
    if(stephalf){
      stepFactor = stepFactor/2
      beta <- beta + stepFactor * step
    }else{
      beta = beta_old - stepFactor*step
    }

    iter = iter + 1
    condition = abs(beta - beta_old)

    if(iter == max.iter) warning("The algorithm did not converge.") & break
  }

  return(list(beta = beta, ll = ll, grd = grd, hessian = Hess, iter = iter))

}

