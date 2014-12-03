#
# This R code is contributed to matrixpls by Huang. It is a part of her dissertation.
#
# Huang, W. (2013). PLSe: Efficient Estimators and Tests for Partial Least Squares
# (Doctoral dissertation). University of California, Los Angeles.
#

# =========== Parameter estimators ===========

#'@title Parameter estimation with an adaptation of PLSc algorithm
#'
#'@description
#'Estimates the model parameters with an adapted version of Dijkstra's (2011) PLSc. The 
#'parameters between latent variables are estimated using two-stages least squares
#'estimation using a composite covariance matrix that has been disattenuated with estimated
#'composite reliabilities.
#'
#'@details
#'\code{params.plsc} estimates the statistical model described by \code{model} with the
#'following steps. If \code{model} is not in the native format, it is converted to the native
#'format containing matrices \code{inner}, \code{reflective}, and \code{formative}. The
#'weights \code{W} and the data covariance matrix \code{S} are used to calculate the composite
#'covariance matrix \code{C} and the indicator-composite covariance matrix \code{IC}.
#'
#'\code{C} is then disattenuated. The reliability estimates used
#'in the dissattenuation process are calculated by using Dijkstra's (2011)
#'correction or with a separate factor analysis for each indicator block. Indicators that are
#'used in a composite but are not specified as reflective indicators for the
#'latent variable that the composite approximates are assumed to be perfectly
#'reliable.
#'
#'The disattenuated \code{C} is used to estimate the the inner part of the model
#'with separate OLS regressions in the same way as in \code{\link{params.regression}}.
#'This differs from Dijkstra's (2011) PLSc estimator that uses
#'2SLS. The reason for not using 2SLS is that PLS is commonly used with models that
#'do not contain a sufficient number of exogenous variables that could be used 
#'as instruments in the 2SLS estimation.
#'
#'The results from the disattenuation process are used as estimates of the \code{reflective}
#'part of the model and the \code{formative} part of the model estimated with separate
#'OLS regressions using the \code{S} and \code{IC} matrices.
#'
#'Those elements of teh indicator-composite covariance matrix \code{IC} that correspond
#'to factor loadings are replaced with the factor loading estimates.
#'
#'The dissattenuation code for Dijkstra's method is adapted from Huang (2013), which is based on
#'Dijkstra (2011).
#'
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param fm factoring method for estimating the corrected factor loadings. \code{"dijkstra"}
#'will use the correction presented by Dijkstra (2011), where the PLS estimates of the 
#'factor loadings for a latent variable are multiplied with a scalar \code{a}, which is 
#'calculated by a simple formula that approximately minimizes squared residual covariances
#'of the factor model. \code{minres}, \code{wls}, \code{gls}, \code{pa}, and \code{ml}
#'will use a factor analysis and passing this parameter through to \code{\link[psych]{fa}}.
#'\code{"cfa"} estimates a maximum likelihood confirmatory factor analysis with \code{\link[lavaan]{lavaan}}
#'
#'@inheritParams matrixpls 
#'
#'@param tsls Should estimation be done with two stage least squares instead of regression
#'
#'@return A named vector of parameter estimates.
#'
#'@family parameter estimators
#'
#'@author Mikko Rönkkö, Wenjing Huang, Theo Dijkstra
#'
#'@references
#'
#' Huang, W. (2013). PLSe: Efficient Estimators and Tests for Partial Least Squares (Doctoral dissertation). University of California, Los Angeles.
#' Dijkstra, T. K. (2011). Consistent Partial Least Squares estimators for linear and polynomial factor models. A report of a belated, serious and not even unsuccessful attempt. Comments are invited. Retrieved from http://www.rug.nl/staff/t.k.dijkstra/consistent-pls-estimators.pdf

#'  
#'@example example/matrixpls.plsc-example.R
#'
#'@export

estimator.PLScLoadings <- function(S, model, W,  ...){
  
  ab <- nrow(W) #number of blocks
  ai <- ncol(W) #total number of indicators
  p <- lapply(1:nrow(W),function(x){which(W[x,]!=0)}) # indicator indices
  
  
  # Calculation of the correlations between the PLS proxies, C:
  C <- W %*% S %*% t(W)	
  
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  
  # Dijkstra's correction
  
  # Determination of the correction factors, based on (11) of Dijkstra, April 7, 2011.
  c2 <- rep(1,ab)
  for (i in 1:ab) {
    idx <- p[[i]]
    if (length(idx) > 1) { # only for latent factors, no need to correct for the single indicator for the phantom LV
      c2[i] <- t(W[i,idx])%*%(S[idx,idx]-diag(diag(S[idx,idx])))%*%W[i,idx]
      c2[i] <- c2[i]/(t(W[i,idx])%*%(W[i,idx]%*%t(W[i,idx])-diag(diag(W[i,idx]%*%t(W[i,idx]))))%*%W[i,idx])
    }
  }
  
  # Dijkstra's formula seems to have problems with negative weights
  c2 <- abs(c2)
  
  c <- sqrt(c2)
  
  # Determination of consistent estimates of the loadings, see (13) of Dijkstra, April 7, 2011.
  L <- model
  
  for (i in 1:ab) {
    idx <- p[[i]]
    L[idx,i] <- c[i]*W[i,idx]
  }
  
  attr(L,"c") <- c

  return(L)
}

estimator.EFALoadings <- function(S, model, W,  ...){
  
  Lp <- model
  
  # Loop over factors and use EFA
  
  # indicator indices based on the reflective model. Coerced to list to avoid the problem that
  # apply can return a list or a matrix depending on whether the number of indicators is equal
  # between the LVs
  
  p_refl <- apply(Lp, 2, function(x){list(which(x!=0))})
  
  for (i in 1:ncol(Lp)) {
    idx <- p_refl[[i]][[1]]
    Lp[idx,i] <- fa(S[idx,idx], fm = fm)$loadings
  }
  
  return(Lp)
}

estimator.CFALoadings <- function(S, model, W,  ..., estimator = "default",
                                  WLS.V = NULL, slotSampleStats = NULL){
  
  Lp <- model # Loading pattern
  
  # Loadings
  parTable <- data.frame(lhs = colnames(Lp)[col(Lp)[Lp!=0]], op = "=~",  rhs = rownames(Lp)[row(Lp)[Lp!=0]], stringsAsFactors = F)
  
  # Errors
  parTable <- rbind(parTable,data.frame(lhs = rownames(Lp)[row(Lp)[Lp!=0]], op = "~~",  rhs = rownames(Lp)[row(Lp)[Lp!=0]], stringsAsFactors = F))
  
  # Factor covariances
  a <- matrix(0,ncol(Lp),ncol(Lp))
  parTable <- rbind(parTable,data.frame(lhs = colnames(Lp)[col(a)[lower.tri(a)]], op = "~~",  rhs = colnames(Lp)[row(a)[lower.tri(a)]], stringsAsFactors = F))
  
  # Factor variances
  parTable <- rbind(parTable,data.frame(lhs = colnames(Lp), op = "~~",  rhs = colnames(Lp), stringsAsFactors = F))
  
  parTable <- cbind(id = as.integer(1:nrow(parTable)), 
                    parTable,
                    user = as.integer(ifelse(1:nrow(parTable)<=sum(Lp!=0),1,0)),
                    group = as.integer(1),
                    free = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-ncol(Lp)),1:nrow(parTable),0)),
                    ustart = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-ncol(Lp)),NA,1)),
                    exo = as.integer(0),
                    label = "",
                    eq.id = as.integer(0),
                    unco = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-ncol(Lp)),1:nrow(parTable),0)),
                    stringsAsFactors = FALSE)
  
  cfa.res <- lavaan::lavaan(parTable, sample.cov = S,
                            sample.nobs = 100, # this does not matter, but is required by lavaan
                            se="none",
                            sample.cov.rescale = FALSE,
                            estimator = estimator,
                            meanstructure = FALSE,
                            WLS.V = WLS.V,
                            slotSampleStats = slotSampleStats)
  
  Lp[Lp==1] <- coef(cfa.res)[1:sum(Lp!=0)]
  
  return(Lp)
  
}

