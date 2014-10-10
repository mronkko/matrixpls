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

params.plsc <- function(S,  model, W, fm = "dijkstra", tsls = FALSE, ...){
  
  nativeModel <- parseModelToNativeFormat(model)

  ab <- nrow(W) #number of blocks
  ai <- ncol(W) #total number of indicators
  p <- lapply(1:nrow(W),function(x){which(W[x,]!=0)}) # indicator indices
  
  
  # Calculation of the correlations between the PLS proxies, C:
  C <- W %*% S %*% t(W)	
  
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  
  # Dijkstra's correction
  
  if(fm == "dijkstra"){
    
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
    L <- matrix(0,ai,ab)
    
    for (i in 1:ab) {
      idx <- p[[i]]
      L[idx,i] <- c[i]*W[i,idx]
    }
    
    # Determination of the quality of the proxies, see (15) of Dijkstra, April 7, 2011.
    
    Q <- c2
    for (i in 1:ab) {
      idx <- p[[i]]
      Q[i] <- c2[i]*(t(W[i,idx])%*%W[i,idx])^2
      # Debug inadmissible Qs
      # if(abs(Q[i]) > 1) browser()
    }
  }
  
  
  # Else use factor analysis
  
  else{
    
    L <- matrix(0,ai,ab)
    Q <- rep(1,ab)
    
    # Use confirmatory factor analysis
    
    if(fm == "cfa"){
      
      Lp <- nativeModel$reflective # Loading pattern
      
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
                                sample.cov.rescale = FALSE)
      
      L[Lp==1] <- coef(cfa.res)[1:sum(Lp!=0)]
      Q <- rowSums(t(L)*W)^2
      
    }
    
    # Else loop over factors and use EFA
    else{
      
      # indicator indices based on the reflective model. Coerced to list to avoid the problem that
      # apply can return a list or a matrix depending on whether the number of indicators is equal
      # between the LVs
      
      p_refl <- apply(nativeModel$reflective, 2, function(x){list(which(x!=0))})
      
      for (i in 1:ab) {
        idx <- p_refl[[i]][[1]]
        L[idx,i] <- fa(S[idx,idx], fm = fm)$loadings
        # Non-factor indicators are assumed to be perfectly reliable and not corrected
        indicator_reliabilities <- L
        indicator_reliabilities[indicator_reliabilities == 0] <- 1
        Q[i] <- sum(indicator_reliabilities[idx,i] * W[i,idx])^2
        
      }
    }
  }
  # Determination of consistent estimates for the correlation between the
  # latent variables, see (15) and (16) of Dijkstra, April 7, 2011.
  R <- C / sqrt(Q) %*% t(sqrt(Q))
  diag(R) <- 1
  
  
  # Choose the specified values and add names
  
  innerRegressionIndices <- which(nativeModel$inner==1, useNames = FALSE)
  
  if(tsls){
    inner <- TwoStageLeastSquaresWithCovarianceMatrixAndModelPattern(R,nativeModel$inner)
  }
  else{
    inner <- regressionsWithCovarianceMatrixAndModelPattern(R,nativeModel$inner)
  }
  
  innerVect <- inner[innerRegressionIndices]
  names(innerVect) <- paste(rownames(inner)[row(inner)[innerRegressionIndices]],"~",
                            colnames(inner)[col(inner)[innerRegressionIndices]], sep="")		
  
  reflective <- nativeModel$reflective
  
  loadingIndices <- which(t(W)!=0)
  loadingVect <- L[loadingIndices]
  names(loadingVect) <- paste(colnames(reflective)[col(reflective)[loadingIndices]], "=~",
                              rownames(reflective)[row(reflective)[loadingIndices]], sep = "")
  
  
  results <- c(innerVect, loadingVect)
  results <- c(results,	params.internal_formative(S, IC, nativeModel))
  
  # Store these in the result object
  attr(results,"C") <- R
  
  if(fm == "dijkstra"){
    attr(results,"c") <- c
  } 
  
  # Fix the IC matrix. Start by replacing correlations with the corrected loadings
  tL <- t(L)
  IC[tL!=0] <- tL[tL!=0]
  
  # Disattenuate the remaining correlations
  IC[tL==0] <- (IC/sqrt(Q))[tL==0]
  
  attr(results,"IC") <- IC
  attr(results,"beta") <- inner
  
  names(Q) <- rownames(inner)
  attr(results,"Q") <- Q
  
  return(results)
  
}


#
# Runs 2SLS defined by model
# using covariance matrix S and places the results in model
#

TwoStageLeastSquaresWithCovarianceMatrixAndModelPattern <- function(S,model){

  # Ensure that S and model have right order
  S <- S[rownames(model),colnames(model)]
  
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
        
        coefs1 <- solve(S[instruments, instruments],S[toBeInstrumented, instruments])
        
        stage1Model[toBeInstrumented, toBeInstrumented] <- 0
        stage1Model[toBeInstrumented, instruments] <- coefs1
      }
      
      
      # Stage 2
      
      S2 <- stage1Model %*% S %*% t(stage1Model)
        
      if(length(independents)==1){
        # Simple regresion is the covariance divided by the variance of the predictor
        model[row,independents] <- S2[row,independents]/S2[independents,independents]
      }
      if(length(independents)>1){
        coefs2 <- solve(S2[independents,independents],S2[row,independents])
        model[row,independents] <- coefs2
      }
      
    }
    # Continue to the next equation
  }
  
  return(model)
}



# =========== Original code by Huang ===========

#This R code is based on Dijkstra's PLSc_UCLA.m with addition of situations when 
#there is only one observed variable loading on a single phantom latent variable,
#so that this one observed variable can serve as both predictor or outcome later on
#when we estimate the structural model regression coefficients.

#define an index function 
PLSc.idx <- function(p,i) {
  (1+p[i]):p[i+1]
}

# Determination of the sign-factors, procedure not described in the paper
# Wold's adjacency concept is not used: all latent variables are taken to be adjacent.
#ab is number of blocks/factors, 
PLSc.sign <- function(ab,p,W,S) {
  s <- matrix(0,ab,ab) #to set up a factor corr matrix
  for (i in 1:(ab-1)) {
    for (j in (i+1):ab) {
      idxi <- PLSc.idx(p,i)
      idxj <- PLSc.idx(p,j)
      s[i,j] <- sign(t(W[idxi])%*%S[idxi,idxj]%*%W[idxj])
    }
  }
  s <- t(s)+s
  return(s)
}

# Normalization of the weight vectors, equation 6 of Dijkstra, April 7, 2011.
PLSc.normalize <- function(ab,p,W,S) {
  for (i in 1:ab) {
    idx <- PLSc.idx(p,i)
    if (length(idx) > 1) { # only for latent factors
      W[idx] <- W[idx]/sqrt(t(W[idx])%*%S[idx,idx]%*%W[idx])
    }
  }
  return(W)
}

# Run through equation (4) of Dijkstra, April 7, 2011.
PLSc.weights <- function(ab,p,Wold,s,S) {
  Wnew <- Wold
  for (i in 1:ab) {
    # Example: idxi=1:3 when i=1, p=c(0 3 6 7), note p is not the p user defined
    idxi <- PLSc.idx(p,i)
    if (length(idxi) > 1) { # if more than one indicator per latent factor, then update the weight.
      # when there is only one indicator per factor, length(idxi)=1, don't update, weight is 1.
      Wnew[idxi] <- 0
      for (j in 1:ab) {
        idxj <- PLSc.idx(p,j)
        Wnew[idxi] <- Wnew[idxi]+s[i,j]*(S[idxi,idxj]%*%matrix(Wold[idxj]))
      }
    }
  }
  return(Wnew)
}

PLSc <- function(S,p,tol=1e-6,tmax=500) {
  ab <- length(p) #number of blocks
  ai <- sum(p) #total number of indicators
  p <- c(0,cumsum(p)) # column indices, e.g if p=c(3,3,1) then p=c(0,3,6,7)
  
  # Normalization of the first weight vectors
  Wold <- PLSc.normalize(ab,p,rep(1,ai),S) #The starting values for W are vectors proportional to(1,1,...,1)',
  # iteration to convergence or max iteration exceeded
  for (t in 1:tmax) {
    s <- PLSc.sign(ab,p,Wold,S) # Determination of the sign-factors
    Wnew <- PLSc.weights(ab,p,Wold,s,S) # weights
    Wnew <- PLSc.normalize(ab,p,Wnew,S) # normalize again
    Wnew <- abs(Wnew) # make loadings positive: see note 2 in Matlab code
    if (max(abs(Wnew-Wold)) < tol) { # test whether convergence has occurred
      break # if yes, break out of the iterations
    } else { # if not, keep on iterating
      Wold <- Wnew
    }
  }
  
  # End of iterations, determination of correlations, loadings and the quality of the proxies
  wdak <- Wnew # wdak contains the PLS mode A estimates for the weights.
  
  
  # Calculation of the correlations between the PLS mode A proxies, Rdak:
  Rdak <- matrix(0,ab,ab)
  for (i in 1:(ab-1)) {
    for (j in (i+1):ab) {
      idxi <- PLSc.idx(p,i)
      idxj <- PLSc.idx(p,j)
      Rdak[i,j] <- t(wdak[idxi])%*%S[idxi,idxj]%*%wdak[idxj]
    }
  }
  Rdak <- Rdak+t(Rdak)+diag(ab)
  ###may use Rdak as PLS estimates to compare with PLSc estimates R###
  
  
  # Determination of the correction factors, based on (11) of Dijkstra, April 7, 2011.
  c2 <- rep(1,ab)
  for (i in 1:ab) {
    idx <- PLSc.idx(p,i)
    if (length(idx) > 1) { # only for latent factors, no need to correct for the single indicator for the phantom LV
      c2[i] <- t(wdak[idx])%*%(S[idx,idx]-diag(diag(S[idx,idx])))%*%wdak[idx]
      c2[i] <- c2[i]/(t(wdak[idx])%*%(wdak[idx]%*%t(wdak[idx])-diag(diag(wdak[idx]%*%t(wdak[idx]))))%*%wdak[idx])
    }
  }
  c <- sqrt(c2)
  
  
  # Determination of consistent estimates of the loadings, see (13) of Dijkstra, April 7, 2011.
  L <- wdak
  for (i in 1:ab) {
    idx <- PLSc.idx(p,i)
    L[idx] <- c[i]*wdak[idx]
  }
  
  # Determination of the quality of the proxies, see (15) of Dijkstra, April 7, 2011.
  Q <- c2
  for (i in 1:ab) {
    idx <- PLSc.idx(p,i)
    Q[i] <- c2[i]*(t(wdak[idx])%*%wdak[idx])^2
  }
  
  # Determination of consistent estimates for the correlation between the
  # latent variables, see (15) and (16) of Dijkstra, April 7, 2011.
  R <- matrix(0,ab,ab)
  # R start with a matrix with all 0, then define the upper triangle of the corr matrix
  for (i in 1:(ab-1)) {
    for (j in (i+1):ab) {
      R[i,j] <- Rdak[i,j]/(sqrt(Q[i]*Q[j]))
    }
  }
  R <- R+t(R)+diag(ab) # to form the complete corr matrix of latent variables
  return(list(R0=Rdak,R=R,L=L,Q=Q,t=t))
} # R0 is PLS estimates for LV corr without Dijkstra's correction.



##Use Two Stage Least Squares to find the structural relationship between the laten factors
TSLS_findindex <- function(v) {
  I <- NULL
  for (j in 1:length(v)) {
    if (v[j] == 1) {
      I <- c(I,j)
    }
  }
  return(I)
}

# R is the covariance matrix of composites
# IB is the model 
#

TSLS_general <- function(R,IB,IC) {
  n <- nrow(IB) # number of etas
  x <- ncol(IC) # number of ksis
  B <- matrix(0,n,n)
  C <- matrix(0,n,x)
  H <- solve(R[1:x,1:x])
  
  for (i in 1:n) {
    I <- TSLS_findindex(IB[i,])
    J <- TSLS_findindex(IC[i,])
    HH <- rbind(cbind(R[x+I,1:x]%*%H%*%R[1:x,x+I],matrix(R[x+I,J],length(x+I),length(J))),
                cbind(matrix(R[J,x+I],length(J),length(x+I)),matrix(R[J,J],length(J),length(J))))
    HH <- ginv(HH)
    h <- rbind(R[x+I,1:x]%*%H%*%R[1:x,i+x],matrix(R[J,i+x]))
    theta <- HH%*%h
    B[i,I] <- theta[1:sum(IB[i,])]
    C[i,J] <- theta[(sum(IB[i,])+1):length(theta)]
  } 
  cov_zeta <- (diag(1,n)-B)%*%R[(x+1):(x+n),(x+1):(x+n)]%*%t(diag(1,n)-B)-C%*%R[(1:x),(1:x)]%*%t(C)
  return(list(B=B,C=C,cov=cov_zeta))
}



# takes covariance matrix return normal theory information matrix
FisherInfo <- function(V) {
  p <- nrow(V)
  q <- p*(p+1)/2
  invV <- solve(V) #solve returns inverse of V
  I <- matrix(0,q,q)
  
  row <- 0
  for (i in 1:p) {
    for (j in 1:i) {
      row <- row + 1
      col <- 0
      # derivative of V with respect to element (i,j)
      dV1 <- matrix(0,p,p)
      dV1[i,j] <- 1
      dV1[j,i] <- 1
      
      for (k in 1:p) {
        for (l in 1:k) {
          col <- col + 1
          # derivative of V wrt element (k,l)
          dV2 <- matrix(0,p,p)
          dV2[k,l] <- 1
          dV2[l,k] <- 1
          I[row,col] <- 0.5*sum(diag(invV%*%dV1%*%invV%*%dV2))
          #          print(c(row,col))
        }
      }
    }
  }
  return(I)
}

# utilities

vech <- function(V) {
  d <- nrow(V)
  v <- rep(0,d*(d+1)/2)
  count <- 1
  for (i in 1:d) {
    for (j in 1:i) {
      v[count] <- V[i,j]
      count <- count+1
    }
  }
  return(v)
}


vech2 <- function(V) {
  d <- nrow(V)
  v <- rep(0,d*(d+1)/2)
  count <- 1
  for (i in 1:d) {
    for (j in i:d) {
      v[count] <- V[j,i]
      count <- count+1
    }
  }
  return(v)
}


vec <- function(V) {
  d <- nrow(V)
  v <- rep(0,d*d)
  count <- 1
  for (i in 1:d) {
    for (j in 1:d) {
      v[count] <- V[i,j]
      count <- count+1
    }
  }
  return(v)
}

# generate file name base
outnamebase <- function(rep) {
  o <- format(rep)
  if (rep < 10) {
    o <- paste("000",o,sep="")
  } else {
    if (rep < 100) {
      o <- paste("00",o,sep="")
    } else {
      if (rep < 1000) {
        o <- paste("0",o,sep="")
      }
    }
  }
  return(o)
}
