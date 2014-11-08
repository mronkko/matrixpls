#'@title GSCA inner estimation
#'
#'@description
#'
#'First step of the GSCA algorithm
#'
#'@template GSCA-details
#'
#'@inheritParams inner.centroid
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@references
#'Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. Psychometrika, 69(1), 81–99. doi:10.1007/BF02295841
#'
#'@example example/gsca-example.R
#'
#'@family inner estimators  
#'@family GSCA functions
#'
#'@export

inner.GSCA <- function(S, W, inner.mod, ...){
  
  # Calculate the composite covariance matrix
  C <- W %*% S %*% t(W)
  
  E <- regressionsWithCovarianceMatrixAndModelPattern(C,inner.mod)
  
  return(E)
}


#'@title GSCA outer estimation
#'
#'@description
#'
#'This implements the second step of the GSCA algorithm.
#'
#'@template GSCA-details
#'
#'@inheritParams outer.modeA
#'
#'@param model A matrixpls model. See \code{\link{matrixpls} for details}
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@references
#'Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. \emph{Psychometrika}, 69(1), 81–99. doi:10.1007/BF02295841
#'
#'@example example/gsca-example.R
#'
#'@family outer estimators
#'@family GSCA functions
#'
#'
#'@export

outer.GSCA <- function(S, W, E, W.mod, model, ...){
  
  nativeModel <- parseModelToNativeFormat(model)
  
  # Update the weights one composite at a time
  
  for(row in 1:nrow(W.mod)){
    
    # We will start by creating a function that returns the sum of residual
    # variances for all regressions where the composite is a predictor
    
    ssFun <- function(Wvect){
      
      resid <- NULL
      
      W[row,W.mod[row,]!=0] <- Wvect
      
      # Start by standardizing the weights because optim does not know that the
      # weights must result in a composite with unit variance
      
      W <- scaleWeights(S, W)
      
      # Calculate the covariance matrix between indicators and composites
      IC <- W %*% S
      
      # Calculate the composite covariance matrix
      C <- IC %*% t(W)
      
      # Regressions from composites to indicators, one indicator at a time
      
      for(i in which(nativeModel$reflective[,row] == 1)){
        c <- which(nativeModel$reflective[i,] == 1)
        coefs <- solve(C[c,c],IC[c,i])
        resid <- c(resid,1 - sum(coefs*IC[c,i]))
      }
      
      # Regressions from composites to composites, one composite at a time
      # The current composite can be either an IV or a DV
      
      DVs <- which(nativeModel$inner[,row] == 1) # current is IV
      if(any(nativeModel$inner[row,] != 0)) DVs <- c(DVs, row) # current is DV
      
      for(i in DVs){
        c <- which(nativeModel$inner[i,] == 1)
        coefs <- solve(C[c,c],C[c,i])
        resid <- c(resid,1 - sum(coefs*C[c,i]))
      }
      
      return(sum(resid))
    }  
    
    # The previous weights are used as starting values
    
    start <- W[row,W.mod[row,]!=0]
    
    # Then we will find the minimum using optim
    
    Wvect <- optim(start, ssFun)$par
    
    # Update the weights based on the estimated parameters
    
    W[row,W.mod[row,]!=0] <- Wvect
    
    # Finally standardize the weights because optim does not know that the
    # weights must result in a composite with unit variance
    
    W <- scaleWeights(S, W)
    
    # Proceed to next composite
  }
  
  return(W)
}

#'@title GSCA optimization criterion
#'
#'@details Optimization criterion for minimizing the sum of all residual
#'variances in the model. 
#'
#'The matrixpls implementation of the GSCA criterion extends the criterion
#'presented by Huang and Takane by including also the minimization of the
#'residual variances of the formative part of the model. The formative
#'regressions in a model are typically specified to be identical to the 
#'weight pattern \code{W.mod} resulting zero residual variances by definition.
#'However, it is possible to specify a formative model that does not completely
#'overlap with the weight pattern leading to non-zero residuals that can be
#'optimizes.
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@return Sum of residual variances.
#'
#'@family Weight optimization criteria
#'@family GSCA functions
#'
#'@export
#'
#'@references
#'
#'Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. 
#'Psychometrika, 69(1), 81–99. doi:10.1007/BF02295841
#'

optim.GSCA <- function(matrixpls.res){
  
  C <- attr(matrixpls.res,"C")
  IC <- attr(matrixpls.res,"IC")
  nativeModel <- attr(matrixpls.res,"model")
  
  reflective <- nativeModel$reflective
  reflective[which(reflective==1)] <- matrixpls.res[grep("=~",names(matrixpls.res))]
  
  formative <- nativeModel$formative
  formative[which(formative==1)] <- matrixpls.res[grep("<~",names(matrixpls.res))]
  
  f <- apply(nativeModel$formative != 0,1,any)
  r <- apply(nativeModel$reflective != 0,1,any)
  endo <- apply(nativeModel$inner != 0,1,any)
  
  inner_resid <- (1 - R2(matrixpls.res)[endo])
  form_resid <- (1 - rowSums(IC[f,] * formative[f,]))
  refl_resid <- (1 - rowSums(t(IC[,r]) * reflective[r,]))
  
  sum(inner_resid, form_resid, refl_resid)
  
}
