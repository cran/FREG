#' Romberg integration
#'
#' Romberg integration is a process of numerical integration. Composite Trapezoidal Rule is used for the approximation of an integral. Then, Richardson
#' extrapolation is used in order to improve previously computed approximations. The range over which the integral is defined is the range in which functional data are defined.
#' When the relative error is infinitesimally small, convergence criterion is fulfilled.
#'
#' @param xbasis  basis functional data object used to represent functional covariate \code{X}. If covariate is a scalar variable,
#' \code{xbasis} is a constant basis functional data object
#'
#' @param bbasis  basis functional data object used to represent beta regression coefficient function for each independent variable
#'
#' @return \item{S matrix}{  a matrix of inner product of two basis objects. The dimensions of the matrix
#' are determined by the number of basis functions. The number of rows is equal to the number of basis functions for \code{X},
#' and the number of columns is equal to the number of basis functions for \code{beta} regression coefficient functions.}
#' @export
#' @import fda
#'
#'
romberg_alg = function(xbasis, bbasis){

  l = NULL

  nrep1 = xbasis$nbasis
  nrep2 = bbasis$nbasis

  N = 20
  eps_step  = 1e-5

  a = xbasis$rangeval[[1]]
  b = xbasis$rangeval[[2]]
  h = b - a # width of the interval

  s = array(0,c(N,nrep1,nrep2))

  # Trapezoidal rule or fist column of Romberg integration matrix

  # 1st iteration (divide by width/2 pi, range is from a to b)
  basismat1 = eval.basis(c(a,b), xbasis) # eval.basis evaluate fourier basis function on a given range
  basismat2 = eval.basis(c(a,b), bbasis)
  s[1,,] = h*as.numeric(crossprod(basismat1, basismat2))/2

  # 2nd to Nth iteration
  for (i in 2:N){
    #h = h/2 #or mean # 2^(i-1)
    x = h/2^(i-1) # sequence
    rng_i = seq(a,b,x) # vector of values in which to evaluate functions
    basismat1 = eval.basis(rng_i, xbasis)
    basismat2 = eval.basis(rng_i, bbasis)
    s[i,,] = h*as.numeric(crossprod(basismat1, basismat2))/2^(i-1)
  }

  # Richardson's extrapolation

  #creation of arrays to stock matrices
  for(i in 1:(N-1)) { #-- Create arrays 'condo.1', 'condo.2' # ideally every next condo should have one less (N-1) place to store matrices but we only ever use first matrix of condo so the approximation is correct
    nam <- paste("condo", i, sep = ".")
    assign(nam, array(0,c(N-1,nrep1,nrep2)))
  }
  # stock condo 3D arrays in a list
  pattern = grep("condo",names(environment()),value=TRUE) # find a match for the name "condo" created in the env
  l = do.call("list",mget(pattern)) # create a list that contains condo arrays

  # iteration 1
  for (i in 1:(N-1)){

    l[[1]][i,,] = (4*s[i+1,,]-s[i,,])/3
  }

  for (j in 2:(N-1)){
    for (i in 1:(N-2)){

      l[[j]][i,,] = (4^j*l[[j-1]][i+1,,]-l[[j-1]][i,,])/(4^j-1)

    }
    err = abs(max(l[[j]][1,,])-max(l[[j-1]][1,,]))
    if(err<eps_step) break

    #print(paste("iteration", j+1))
  }

  if (xbasis$nbasis != 1 & j == (N-1)) warning("The algorithm did not converge.")

  #in case of a constant with value of 1
  if(xbasis$nbasis == 1 & bbasis$nbasis == 1) {
    S = s[1,,]
  }else {
    S = l[[j]][1,,]
  }
  return(S)
}
