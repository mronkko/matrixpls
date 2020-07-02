# =========== Reliability estimators ===========

#'@title Reliabilities as products of weights and loadings
#'
#'@description
#'Calculates reliabilities as a matrix product of loadings and weights.
#'
#'@param S the data covariance matrix
#'
#'@param loadings matrix of factor loading estimates
#'
#'@param W matrix of weights
#'
#'@param ... All other arguments are ignored.
#'
#'@return a named vector of estimated composite reliabilities.
#'
#'@name reliabilityEstim
#'

NULL

#'@describeIn reliabilityEstim Reliability estimation based on weights and loadings.
#'@export

reliabilityEstim.weightLoadingProduct <- function(S, loadings, W, ...){
  
  Q <- diag(W %*% loadings)^2
  
  # Any composite with no reflective indicators is fixed to be perfectly reliable
  Q[apply(loadings==0,2,all)] <- 1
  
  Q
}

#'@describeIn reliabilityEstim Reliability estimation with Cronbach's alpha
#'@export

reliabilityEstim.alpha <- function(S, loadings, W, ...){
  
  S <- stats::cov2cor(S)
  
  Q <- apply(W,1,function(w){
    i <- which(w!=0)
    if(length(i) == 1) return(1)
    # Call the matrixpls alpha function
    alpha(S[i, i])
  })
  Q
}
