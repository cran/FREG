#' Estimation of Hessian matrix
#'
#' Second derivative of log-likelihood function
#'
#' @param x  a design matrix which is a product of inner product of basis functions and basis coefficients of functional covariate \code{X}
#' @param y  a response variable of class \code{factor}
#' @param beta  initial values for beta regression coefficients and intercepts
#'
#' @return \item{Hess}{  a Hessian matrix at the estimated optimum}
#' @export
#'

Hessian = function(x,y,beta){

  inv_link = gradient(x,y,beta)$inv_link
  z1 = gradient(x,y,beta)$z1
  z2 = gradient(x,y,beta)$z2
  p1 = gradient(x,y,beta)$p1
  p2 = gradient(x,y,beta)$p2
  kr_delta = gradient(x,y,beta)$kr_delta
  kr_delta1 = gradient(x,y,beta)$kr_delta1
  ylev = levels(y)
  lylev = length(levels(y))
  q = lylev-1L

  P = inv_link(z1) - inv_link(z2)
  # the first derivative of distribution function is a density function
  # the second is as follows:
  ppa = inv_link(z1) * (1 - 3*inv_link(z1) + 2*inv_link(z1)*inv_link(z1))
  ppb = inv_link(z2) * (1 - 3*inv_link(z2) + 2*inv_link(z2)*inv_link(z2))

  # beta

  llbb = matrix(0, ncol(x), ncol(x))

  for(i in 1:ncol(x)){
    for(k in 1:ncol(x)){
      llbb[i,k] = sum((x[,i]*x[,k]/P)*((p2-p1)^2/P+ppb-ppa))
    }
  } #upper left

  # beta-alpha

  llba = NULL
  fst = ((p1-p2)*(p1*kr_delta-p2*kr_delta1))
  scnd = as.vector(ppa)*kr_delta-as.vector(ppb)*kr_delta1
  llba = t(x/(as.vector(P)))%*%(-fst/as.vector(P)+scnd) # upper right
  t(llba) # lower left

  # alpha

  llaa = matrix(0, q, q)

  for(i in 1:q){
    for(m in 1:q){
      a  = ppb*kr_delta1[,i]*kr_delta1[,m]-ppa*kr_delta[,i]*kr_delta[,m] # hors de la diagonale 0
      b = (p1*kr_delta[,i] - p2*kr_delta1[,i])*(p1*kr_delta[,m] - p2*kr_delta1[,m])
      llaa[i,m] = sum(b/P^2+a/P)
    }
  }

  block12 = cbind(llbb, llba)
  block34 = cbind(t(llba), llaa)
  Hess = rbind(block12, block34)
  colnames(Hess)[(ncol(x)+1):dim(Hess)[2]] = paste(ylev[1:length(ylev)-1], ylev[2:length(ylev)], sep = "|")
  rownames(Hess)[(ncol(x)+1):dim(Hess)[2]] = colnames(Hess)[(ncol(x)+1):dim(Hess)[2]]
  colnames(Hess)[1:ncol(x)] = paste("Z", 1:ncol(x))
  rownames(Hess)[1:ncol(x)] = paste("Z", 1:ncol(x))

  return(Hess)

}
