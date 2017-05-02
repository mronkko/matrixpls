#'@title Predict method for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{predict} predict.
#'Predicts the reflective indicators of endogenous latent variables using
#'estimated model and data for the indicators of exogenous latent variables
#'
#'@template postestimationFunctions
#'
#'@param newData A data frame or a matrix containing data used for prediction.
#'
#'@param predictionType "exogenous" (default) predicts indicators from exogenous composites. "redundancy"
#'and "communality" are alternative strategies described by Chin (2010).
#'
#'@param means A vector of means of the original data used to calculate intercepts for the
#'linear prediciton equations. If not provided, calculated from the new data or assumed zero.
#'
#'@return a matrix of predicted values for reflective indicators of endogenous latent variables.
#'
#'@references
#'
#'Wold, H. (1974). Causal flows with latent variables: Partings of the ways in the light of NIPALS modelling. \emph{European Economic Review}, 5(1), 67–86. doi:10.1016/0014-2921(74)90008-7
#'
#'Chin, W. W. (2010). How to write up and report PLS analyses. In V. Esposito Vinzi, W. W. Chin, J.
#'Henseler, & H. Wang (Eds.), Handbook of partial least squares (pp. 655–690). Berlin Heidelberg: Springer.
#'
#'@export
#'
#'@method predict matrixpls
#'
#'@export


predict.matrixpls <- function(object, newData, 
                              predictionType = "exogenous",
                              means = NULL, ...){
  
  nativeModel <- attr(object,"model")
  exog <- rowSums(nativeModel$inner)==0
  W <- attr(object,"W")
  inner <- attr(object,"inner")
  
  
  reflective <- attr(object,"reflective")
  
  # Reorder the variables in newData
  data <- matrix(NA, nrow(newData), 0)
  
  for(name in colnames(W)){
    if(! name %in% colnames(newData)){
      data <- cbind(data,NA)
    }
    else{
      data <- cbind(data,newData[,name])
    }
  }
  
  colnames(data) <- colnames(W)
  
  if(is.null(means)){
    means <- apply(data,2,mean,na.rm = TRUE)
    means[is.nan(means)] <- 0
  }
  else{
    means <- means[colnames(data)]
  }
  
  # Calculate prediction matrix 
  
  
  if(predictionType == "exogenous"){
    
    beta <- diag(exog)
    
    # Predict endog LVs using reduced form equations
    Beta <- inner[! exog, ! exog]
    Gamma <- inner[! exog, exog]
    
    beta[! exog, exog] <- t(Gamma) %*% solve(diag(nrow(Beta))- t(Beta))
    
  } else if(predictionType == "redundancy"){
    
    beta <- inner
    diag(beta)[exog]<-1
    
  } else if(predictionType == "communality"){
    
    beta <- diag(nrow(inner))
    
  } else stop('predictionType must be "exogenous", "redundancy", or "communality"')
  
  predictionMatrix <- t(W) %*% t(beta) %*% t(reflective)
  
  intercepts <- means - means %*% predictionMatrix
  
  # Do the predictions
  pred <- data %*% predictionMatrix
  sweep(pred, 2, intercepts, "+")

}
