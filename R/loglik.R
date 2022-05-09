#' Log-likelihood function
#'
#' @param x  a design matrix which is a product of inner product of basis functions and basis coefficients of functional covariate \code{X}
#' @param y  a response variable of class \code{factor}
#' @param beta  initial values for beta regression coefficients and intercepts
#'
#' @return \item{ll}{  a value of the log-likelihood function at the estimated optimum}
#' @export
#'

loglik = function(x,y,beta){ # score function

  #beta = start
  n = nrow(x)
  wt = rep(1, n) #weights
  xc = ncol(x)
  ind_xc = seq_len(xc)
  lylev = length(levels(y))
  q = lylev-1L
  ind_q = seq_len(q)

  inv_link = function(x) 1/(1+exp(-x)) # logit link function

  alphas <- c(-1e2, beta[xc + ind_q], 1e2)
  eta = x %*% beta[ind_xc] # linear predictor
  z1 = alphas[y+1L] - eta
  z2 = alphas[y] - eta
  fitted <- inv_link(z1) - inv_link(z2)
  if(all(fitted > 0)) ll = -sum(wt*log(fitted)) else ll = Inf # to ensure that thresholds are increasing
  return(ll)
}
