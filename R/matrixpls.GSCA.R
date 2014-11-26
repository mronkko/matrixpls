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
  
  inner <- nativeModel$inner
  
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  
  # Estimate the reflective and formative parts of the model
  
  formative <- nativeModel$formative
  
  for(row in which(rowSums(formative!=0)>0)){
    independents <- which(formative[row,] != 0)
    formative[row,independents] <- solve(S[independents,independents],IC[row,independents])
  }
  
  reflective <- nativeModel$reflective
  
  for(row in which(rowSums(reflective!=0)>0)){
    independents <- which(reflective[row,] != 0)
    reflective[row,independents] <- solve(S[independents,independents],IC[independents, row])
  }  
  
  # E, formative, and reflective form the A matrix of GSCA. In this implementation,
  # the full A matrix is never calculated, but we use the elements of 
  # E, formative, and reflective directly.

  # Composite correlations implied by A
  C <- W %*% S %*% t(W)
  
  # Update the weights one composite at a time
  
  for(row in 1:nrow(W.mod)){
    
    # TODO: Compare with the old implementation
    
    # Indicator indices
    i <- which(W.mod[row,]!=0)
    
    # correlations between the indicators and dvs
    ic <- NULL

    # Calculate the covariance matrix between indicators and composites
    # This needs to be recalculated for each composite because W is updated
    # immediately
    
    IC <- E %*% W %*% S
    
    # The composite is dependent in inner model
    
    if(any(inner[row,]!=0)){
      ic <- cbind(ic,IC[row,i])
    }
    
    # Composites that this composite depends on
    
    for(j in which(inner[,row]!=0)){
      ic <- cbind(ic,IC[j,i] * C[row,j])
    }

    # The composite is dependent in the formative model
    # and the weight pattern and formative indicators do not
    # overlap completely
    
    if(any(formative[row,]!=0) && 
         ! identical((formative[row,] ==0), (W.mod[row,] == 0))){
      ic <- cbind(ic, W[row,i] %*% S[i,i])
    }
    
    # Indicators that this composite depends on
    
    for(j in which(reflective[row,]!=0)){
      ic <- cbind(ic, S[i,j])
    }

    if(is.null(ic)) stop(paste("Composite '",rownames(inner)[row],
                               "' is neither a dependent nor an independent variable in a GSCA outer model regressions. Estimation fails.",sep=""))
    
    Wvect <- solve(S[i,i],apply(ic,1,mean))
    
    # Update the weights based on the estimated parameters
    
    W[row,W.mod[row,]!=0] <- Wvect
    
    # Finally standardize the weights 
    
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
