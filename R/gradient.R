#' Estimation of gradient of log-likelihood function
#'
#' First derivative of log-likelihood function
#'
#' @param x  a design matrix which is a product of inner product of basis functions and basis coefficients of functional covariate \code{X}
#' @param y  a response variable of class \code{factor}
#' @param beta  initial values for beta regression coefficients and intercepts
#'
#' @return \item{grd}{  a vector of gradient values at the estimated optimum}
#' @export
#'
#'
gradient <- function(x,y,beta){

  n = nrow(x)
  wt <- rep(1, n) #weights
  xc <- ncol(x)
  ind_xc <- seq_len(xc)
  lylev = length(levels(y))
  q = lylev-1L
  ind_q = seq_len(q)

  inv_link = function(x) 1/(1+exp(-x)) # logit link function
  den_fun = function(x) exp(-x)/(1+exp(-x))^2 # probability density function

  alphas <- c(-1e2, beta[xc + ind_q], 1e2)
  eta <- x %*% beta[ind_xc]
  z1 <- alphas[y+1L] - eta
  z2 <- alphas[y] - eta
  P <- inv_link(z1) - inv_link(z2)
  p1 <- as.vector(den_fun(z1))
  p2 <- as.vector(den_fun(z2))
  wfit = wt/P
  grd = rep(0, xc+q)
  grd[1:xc] = crossprod(x, wfit*(p1 - p2))
  Y <- matrix(0, n, q)
  ymat = 1 * (col(matrix(0, nrow(x), lylev)) == y)
  kr_delta = ymat[, -(lylev), drop = FALSE] # (y) #minus the last category
  kr_delta1 = ymat[, -1, drop = FALSE] # (y-1) #minus the first category
  nominator = kr_delta*p1 - kr_delta1*p2
  grd[xc+1:length(ind_q)] = -crossprod(nominator, wfit)

  return(list(inv_link = inv_link, grd = grd, z1 = z1, z2 = z2, p1 = p1, p2 = p2, kr_delta = kr_delta,
              kr_delta1 = kr_delta1))
}
