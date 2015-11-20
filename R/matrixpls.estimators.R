#'@title Parameter estimation with separate regression analyses
#'
#'@description
#'
#'Estimates the parameters of a model matrix with OLS regression.
#'
#'@details
#'
#'Providing \code{C} or \code{IC} allows for using disattenuated or otherwise
#'adjusted correlation matrices. If not provided, these matrices are calculated using \code{S} and
#'\code{W}.
#'
#'@inheritParams matrixpls-common
#'
#'@param modelMatrix A model matrix with dependent variables on rows, independent variables on colums, and
#'non-zero elements indicating relationships. Can be either \code{inner}, \code{reflective},
#'or \code{formative} matrix.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A matrix with estimated parameters.
#'
#'@family estimators
#'
#'@export
#'

estimator.ols <- function(S, model, W, ..., C = NULL, IC = NULL){
  
  covIV <- covDV <- NULL
  
  # Calculate the composite covariance matrix
  if(is.null(C)) C <- W %*% S %*% t(W)
  
  # Calculate the covariance matrix between indicators and composites
  if(is.null(IC)) IC <- W %*% S
  
  # Covariances between the independent variables
  for(m in list(S, C)){
    if(all(colnames(model) %in% colnames(m))){
      covIV <- m[colnames(model),colnames(model)]
      break()
    }
  }
  
  # Covariances between IVs and DVS
  
  for(m in list(S, C, IC)){
    if(all(rownames(model) %in% rownames(m)) &&
       all(colnames(model) %in% colnames(m))){
      covDV <- m[rownames(model),colnames(model)]
      break()
    }
  }
  
  for(row in 1:nrow(model)){
    
    independents <- which(model[row,]!=0, useNames = FALSE)
    
    if(length(independents)==1){
      # Simple regresion is the covariance divided by the variance of the predictor
      model[row,independents] <- covDV[row,independents]/covIV[independents,independents]
    }
    if(length(independents)>1){
      coefs <- solve(covIV[independents,independents],covDV[row,independents])
      model[row,independents] <- coefs
    }
  }
  
  return(model)
}

#'@title Parameter estimation with two-stage least squares regressions
#
#'@description
#'
#'Estimates the parameters of a model matrix (\code{inner}) with two-stage least squares regressions.
#'The exogenous variables are used as instruments. 
#'
#'@details
#'
#'Providing \code{C} allows for using disattenuated or otherwise
#'adjusted correlation matrices. If not provided, this matrix is calculated using \code{S} and
#'\code{W}.
#'
#'@inheritParams matrixpls-common
#'
#'@param modelMatrix A model matrix with dependent variables on rows, independent variables on colums, and
#'non-zero elements indicating relationships. Onyl the \code{inner} matrix can be used.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A matrix with estimated parameters.
#'
#'@family estimators
#'
#'@export
#'

estimator.tsls <- function(S, modelMatrix, W, ..., C){
  
  # Calculate the composite covariance matrix
  if(is.null(C)) C <- W %*% S %*% t(W)
  
  exog <- apply(model == 0, 1, all)
  endog <- ! exog
  
  # Use all exogenous variables as instruments
  
  instruments <- which(exog) 
  for(row in 1:nrow(model)){
    
    independents <- which(model[row,]!=0, useNames = FALSE)
    
    
    # Check if this variable is endogenous
    
    if(length(independents) > 0 ){
      
      # Stage 1
      
      needInstruments <- which(model[row,]!=0 & endog, useNames = FALSE)
      
      # This is a matrix that will transform the instruments into independent
      # variables. By default, all variables are instrumented by themselves
      
      stage1Model <- diag(nrow(model))
      
      for(toBeInstrumented in needInstruments){
        # Regress the variable requiring instruments on its predictors excluding the current DV
        
        coefs1 <- solve(C[instruments, instruments],C[toBeInstrumented, instruments])
        
        stage1Model[toBeInstrumented, toBeInstrumented] <- 0
        stage1Model[toBeInstrumented, instruments] <- coefs1
      }
      
      
      # Stage 2
      
      C2 <- stage1Model %*% C %*% t(stage1Model)
      
      if(length(independents)==1){
        # Simple regresion is the covariance divided by the variance of the predictor
        model[row,independents] <- C2[row,independents]/C2[independents,independents]
      }
      if(length(independents)>1){
        coefs2 <- solve(C2[independents,independents],C2[row,independents])
        model[row,independents] <- coefs2
      }
      
    }
    # Continue to the next equation
  }
  
  return(model)
}

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
#'\code{"dijkstra"}
#'will use the correction presented by Dijkstra (2011), where the PLS estimates of the 
#'factor loadings for a latent variable are multiplied with a scalar \code{a}, which is 
#'calculated by a simple formula that approximately minimizes squared residual covariances
#'of the factor model. \code{minres}, \code{wls}, \code{gls}, \code{pa}, and \code{ml}
#'will use a factor analysis and passing this parameter through to \code{\link[psych]{fa}}.
#'\code{"cfa"} estimates a maximum likelihood confirmatory factor analysis with \code{\link[lavaan]{lavaan}}
#'
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@inheritParams params.regression 
#'
#'@param ... All other arguments are ignored.
#'
#'@return A matrix with estimated parameters.
#'
#'@family estimators
#'
#'@export
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


#
# Parts of this R code is contributed to matrixpls by Huang. It is a part of her dissertation.
#
# Huang, W. (2013). PLSe: Efficient Estimators and Tests for Partial Least Squares
# (Doctoral dissertation). University of California, Los Angeles.
#


estimator.PLScLoadings <- function(S, model, W,  ...){
  
  ab <- nrow(W) #number of blocks
  ai <- ncol(W) #total number of indicators
  
  L <- model
  
  # Indicator indices based on the reflective model. Coerced to list to avoid the problem that
  # apply can return a list or a matrix depending on whether the number of indicators is equal
  # between the LVs
  
  p_refl <- apply(L, 2, function(x){list(which(x!=0))})
  
  # Calculation of the correlations between the PLS proxies, C:
  C <- W %*% S %*% t(W)	
  
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  
  # Dijkstra's correction
  
  # Determination of the correction factors, based on (11) of Dijkstra, April 7, 2011.
  c2 <- rep(1,ab)
  for (i in 1:ab) {
    idx <- p_refl[[i]][[1]]
    if (length(idx) > 1) { # only for latent factors, no need to correct for the single indicator for the phantom LV
      c2[i] <- t(W[i,idx])%*%(S[idx,idx]-diag(diag(S[idx,idx])))%*%W[i,idx]
      c2[i] <- c2[i]/(t(W[i,idx])%*%(W[i,idx]%*%t(W[i,idx])-diag(diag(W[i,idx]%*%t(W[i,idx]))))%*%W[i,idx])
    }
  }
  
  # Dijkstra's formula seems to have problems with negative weights
  c2 <- abs(c2)
  
  c <- sqrt(c2)
  
  # Determination of consistent estimates of the loadings, see (13) of Dijkstra, April 7, 2011.
  
  for (i in 1:ab) {
    idx <- p_refl[[i]][[1]]
    
    if(length(idx) > 1){
      L[idx,i] <- c[i]*W[i,idx]
    }
  }
  
  attr(L,"c") <- c
  return(L)
}

#'@title Parameter estimation with per-block exploratory factor analysis
#'
#'@description Estimates the factor loading parameters one latent variable at at time with exploratory
#'factor analysis using the \code{\link[psych]{fa}} function from the \code{psych} package.
#'
#'@details
#'Estimating an unconstrained single factor model requires three indicators. The loadings of 
#'single indicator factors are estimated as 1 and two indicator factors as estimated by the
#'square root of the indicator correlation.
#'
#'@param fm factoring method for estimating the factor loadings. Passed through to \code{\link[psych]{fa}}.
#'
#'@inheritParams params.regression 
#'
#'@return A named vector of parameter estimates.
#'
#'@family parameter estimators
#'
#'@export

estimator.EFALoadings <- function(S, model, W,  ... , fm = "minres"){
  
  L <- model
  
  # Indicator indices based on the reflective model. Coerced to list to avoid the problem that
  # apply can return a list or a matrix depending on whether the number of indicators is equal
  # between the LVs
  
  p_refl <- apply(L, 2, function(x){list(which(x!=0))})
  
  # Loop over factors and use EFA
  
  for (i in 1:ncol(L)) {
    idx <- p_refl[[i]][[1]]
    
    if(length(idx) == 1){ # Single indicator
      L[idx,i] <- 1
    }
    else if(length(idx) == 2){ # Two indicators
      L[idx,i] <- sqrt(S[idx,idx][2])
    }
    else if(length(idx) >= 3){ # Three or more indicators
      L[idx,i] <- psych::fa(S[idx,idx], fm = fm)$loadings
    }
  }
  return(L)
}

#'@title Parameter estimation with confirmatory factor analysis of the full model
#'
#'@description Estimates a maximum likelihood confirmatory factor analysis with \code{\link[lavaan]{lavaan}}
#'
#'@inheritParams params.regression 
#'
#'@return A named vector of parameter estimates.
#'
#'@family parameter estimators
#'
#'@export

estimator.CFALoadings <- function(S, model, W, ...){
  
  L <- model # Loading pattern
  
  hasReflIndicators <- which(apply(L!=0,2,any))
  # Loadings
  parTable <- data.frame(lhs = colnames(L)[col(L)[L!=0]], op = "=~",  rhs = rownames(L)[row(L)[L!=0]], stringsAsFactors = F)
  
  # Errors
  parTable <- rbind(parTable,data.frame(lhs = rownames(L)[row(L)[L!=0]], op = "~~",  rhs = rownames(L)[row(L)[L!=0]], stringsAsFactors = F))
  
  # Factor covariances
  if(length(hasReflIndicators)>1){
    a <- matrix(0,length(hasReflIndicators),length(hasReflIndicators))
    parTable <- rbind(parTable,data.frame(lhs = colnames(L)[hasReflIndicators][col(a)[lower.tri(a)]],
                                          op = "~~",  rhs = colnames(L)[hasReflIndicators][row(a)[lower.tri(a)]], stringsAsFactors = F))
  }
  
  # Factor variances
  parTable <- rbind(parTable,data.frame(lhs = colnames(L)[hasReflIndicators], op = "~~",  rhs = colnames(L)[hasReflIndicators], stringsAsFactors = F))
  
  parTable <- cbind(id = as.integer(1:nrow(parTable)), 
                    parTable,
                    user = as.integer(ifelse(1:nrow(parTable)<=sum(L!=0),1,0)),
                    group = as.integer(1),
                    free = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-length(hasReflIndicators)),1:nrow(parTable),0)),
                    ustart = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-length(hasReflIndicators)),NA,1)),
                    exo = as.integer(0),
                    label = "",
                    eq.id = as.integer(0),
                    unco = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-length(hasReflIndicators)),1:nrow(parTable),0)),
                    stringsAsFactors = FALSE)
  
  args <- list(model= parTable, sample.cov = S,
               sample.nobs = 100, # this does not matter, but is required by lavaan
               se="none",
               sample.cov.rescale = FALSE,
               meanstructure = FALSE)
  
  e <- list(...)
  f <- formals(lavaan::lavaan)
  include <- intersect(names(f), names(e))
  args[include] <- e[include]
  
  cfa.res <- do.call(lavaan::lavaan, args)
  
  L[L==1] <- lavaan::coef(cfa.res)[1:sum(L!=0)]
  return(L)
  
}
