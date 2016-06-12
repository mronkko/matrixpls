#'@title Cross-validation of predictions from matrixpls results
#'
#'@description
#'\code{matrixpls.crossvalidate} Calculates cross-validation predictions using \code{matrixpls}.
#'
#'@details In cross-validation, some elements of the data matrix are set to missing and then
#'predicted based on the remaining observations. The process is repeated until all elements
#'of the data have been predicted. 
#'
#'Cross-validation is typically applied by dividing the data into two groups, the training
#'sample and the validation sample. The prediction model is calculated based on the 
#'training sample and used to calculate predictions for the validation sample.
#'
#'In blindfolding, the data are not omitted case wise, but elements of the data are omitted
#'diagonally. After this, imputation is applied to missing data and the prediction model
#'is calibrated with the dataset containing also the imputations. The imputed values are then 
#'predicted with the model.
#'
#'@param data Matrix or data frame containing the raw data.
#'
#'@param nGroup The number of groups to divide the data into.
#'
#'@param predictFun The function used to calculate the predictions.
#'
#'@param blindfold Whether blindfolding should be used instead of holdout sample cross-validation.
#'
#'@param imputationFun The function used to impute missing data before blindfold prediction. 
#'If \code{NULL}, simple mean substitution is used.
#'
#'@param ... All other arguments are passed through to \code{\link{matrixpls}} and \code{predictFun}.
#'
#'@inheritParams matrixpls-common
#'
#'@return A matrix containing predictions calculated with cross-validation.
#'
#'@export
#'
#'@reference
#'Chin, W. W. (2010). How to write up and report PLS analyses. In V. Esposito Vinzi, W. W. Chin, J.
#'Henseler, & H. Wang (Eds.), Handbook of partial least squares (pp. 655–690). Berlin Heidelberg: Springer.
#'
#'Chin, W. W. (1998). The partial least squares approach to structural equation modeling. 
#'In G. A. Marcoulides (Ed.), Modern methods for business research (pp. 295–336). 
#'Mahwah, NJ: Lawrence Erlbaum Associates Publishers.


matrixpls.crossvalidate <- function(data, model, ..., predictFun = stats::predict,nGroup = 4,
                                    blindfold = FALSE,
                                    imputationFun = NULL){
  
  nativeModel <- parseModelToNativeFormat(model)
  
  data <- as.matrix(data)
  data <- data[,rownames(nativeModel$reflective)]

  assertive::assert_all_are_in_closed_range(nGroup,2,nrow(data))
  
  
  # Construct the omission pattern
  if(blindfold){
    pattern <- data
    pattern[] <- NA

    counter <- 0
    
    for(i in 1:ncol(nativeModel$reflective)){
      # identify reflective blocks
      j <- which(nativeModel$reflective[,i] != 0)
      if(length(j) > 0){
        pattern[,j] <- rep_len(1:nGroup + counter, length(j)*nrow(data))
        counter <- counter + nGroup
      }
      # Set to the number of groups taking the number of blocks into consideration
      nGroup <- counter
    }
    
    if(is.null(imputationFun)){
      # Use mean substitution
      imputationFun <- function(d){apply(d,2,function(x){x[is.na(x)] <- mean(x,na.rm = TRUE);x})}
    }
  }
  else{
    pattern <- (row(data)-1)%%nGroup+1
  }
  
  crossvalidate.out <- matrix(NA, nrow(data), ncol(data),
                              dimnames = dimnames(data))
  
  for(group in 1:nGroup){
    bdata <- data
    bdata[pattern == group] <- NA
    
    if(blindfold){
      bdata <- imputationFun(bdata)
    }
    
    S <- stats::cov(bdata, use = "complete.obs")
    matrixpls.res <- matrixpls(S, model = model, ...)
    
    if(!blindfold){
      bdata <- data
    }
    
    pred <- predictFun(matrixpls.res, bdata, ...)
    crossvalidate.out[pattern == group] <- pred[pattern == group]
  }
  
  attr(crossvalidate.out, "groups") <- pattern
  crossvalidate.out
}

#'@title Q2	predictive relevance statistics
#'
#'@description Calculates Q2 predictive relevance statistics based on comparing predictions
#'and real data.
#'
#'@details The Q2 statistic is calculated as \code{1-sse/sso} where \code{sse} is the sum of 
#'squared prediction errors based on comparison of the \code{originalData} and 
#'\code{predictedData} and \code{sso} is based on prediction with mean. If the predicted
#'data contain the \code{groups} attribute, which indicates the groups used in blindfolding
#'or cross-validation, the means are calculated separately for each group excluding the 
#'predicted group from the calculation.
#'
#'@param originalData A matrix or a data.frame containing the original data.
#'
#'@param predictedData A matrix or a data.frame containing the predicted data that are compared
#'against the original data to calculate the predictive relevance statistic.
#'
#'@inheritParams matrixpls-common
#'
#'@return A list with \code{total}, \code{block}, and \code{indicator} elements containing
#'the Q2 predictive relevance statistics for the full dataset, for each indicator block, and
#'for each indicator
#'
#'@export
#'

q2 <- function(originalData, predictedData, model = NULL){
  
  pattern <- attr(predictedData,"groups")
  
  # If there is no pattern, assume that each observation is predicted separately from
  # all other observations
  if(is.null(pattern)) pattern <- matrix(1,nrow(predictedData), ncol(originalData))
  
  # Choose only variables that exists in both data frames / matrices and 
  # sort the variables into same order
  
  v <- intersect(colnames(originalData), colnames(predictedData))
  originalData <- originalData[,v]
  predictedData <- predictedData[,v]
  
  means <- matrix(NA,nrow(predictedData), ncol(predictedData))
  
  groups <- unique(as.vector(pattern))
  for(group in groups){
    for(i in 1:ncol(originalData)){
      if(length(groups)!= 1){
        j <- pattern[,i]==group
        if(any(j)){
          k <- which(j)
          l <- which(!j)
          means[k,i] <- mean(originalData[l,i])
          
        }
      }
      else{
        means[,i] <- mean(originalData[,i])
      }
    }
  }
  
  sse <- apply((originalData-predictedData)^2,2,sum)
  sso <- apply((originalData-means)^2,2,sum)
  
  if(! is.null(model)){
    
    model <- parseModelToNativeFormat(model)
    reflective <- model$reflective
    q2 <- list(total = 1-sum(sse)/sum(sso),
               block = apply(reflective[v,],2,function(inBlock){
                 1-sum(sse[inBlock==1])/sum(sso[inBlock==1])
               }),
               indicator = 1-sse/sso)
  }
  else{
    q2 <- list(total = 1-sum(sse)/sum(sso),
               indicator = 1-sse/sso)
  }
  
  class(q2) <- "matrixplsq2"
  q2
}

#'@S3method print matrixplsq2

print.matrixplsq2 <- function(x, ...){
  cat("\n Q2 predictive relevance statistics\n")
  cat("\n Overall Q2\n")
  cat(x$total, ...)
  if(! is.null(x$block)){
    cat("\n Block Q2\n")
    print(x$block, ...)
  }
  cat("\n Indicator Q2\n")
  print(x$indicator, ...)
}
