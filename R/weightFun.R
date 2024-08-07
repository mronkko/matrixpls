# =========== Indicator weighting algorithms  ===========

#'@title Indicator weight algoritms
#'
#'@description
#'Estimates a weight matrix using Partial Least Squares or a related algorithm.
#'
#'@template modelSpecification
#'
#'@template weightSpecification
#'
#'@inheritParams matrixpls-common
#'@inheritParams matrixpls-functions
#'
#'@param validateInput A boolean indicating whether the validity of the parameter values should be tested.
#'
#'@param ... All other arguments are passed through to other estimation functions.
#'
#'@param standardize A boolean indicating whether \code{S} the weights should be scaled to produce
#'standardized composites.
#'
#'@template weights-return
#'
#'@templateVar attributes weightFun.pls,S E iterations converged history
#'@template attributes
#'
#'@name weightFun
NULL


#'@describeIn weightFun partial Least Squares and other iterative two-stage weight algorithms.
#'
#'@details
#'
#'\code{weightFun.pls} calculates indicator weights by calling the 
#'\code{innerEstim} and \code{outerEstim} iteratively until either the convergence criterion or
#'maximum number of iterations is reached and provides the results in a matrix.

#'@param outerEstim A function or a list of functions used for outer estimation. If
#'the value of this parameter is a function, the same function is applied to all
#'composites. If the value is a list, the composite \code{n} is estimated
#'with the estimator in the \code{n}th position in the list. If this argument is
#'\code{NULL} the \code{\link{outerEstim.modeA}} is used for all composites that are linked to at least
#'one indicator in the \code{reflective} matrix.\code{\link{outerEstim.modeB}} is used for all other
#'composites. See \code{\link{outerEstim}}.
#'
#'
#'@param tol Decimal value indicating the tolerance criterion for convergence. 
#'
#'@param iter An integer indicating the maximum number of iterations.

#'@param variant Choose either Lohmöller's (\code{"lohmoller"}, default) or Wold's (\code{"wold"}) 
#'variant of PLS. In Wold's variant the inner and outer estimation steps are repeated for each
#'indicator block whereas in Lohmöller's variant the weights for all composites are calculated
#'simultaneously. 

#'@export

weightFun.pls <- function(S, model, W.model,
                       outerEstim = NULL, 
                       innerEstim = innerEstim.path, ..., 
                       convCheck = convCheck.absolute,
                       variant = "lohmoller",
                       tol = 1e-05, iter = 100, validateInput = TRUE) {
  
  if(validateInput){
    
    # All parameters must have values
    assertthat::assert_that(all(!is.na(formals())))
    
    
    # S must be symmetric and a valid covariance matrix
    assertthat::assert_that(is.matrix(S))
    assertthat::assert_that(is.symmetric.matrix(S))
    assertthat::assert_that(is.positive.semi.definite(S))
    
    # W.model must be a real matrix and each indicators must be
    # linked to at least one composite and each composite at least to one indicator

    assertthat::assert_that(is.matrix(W.model))
    assertthat::assert_that(all(is.numeric(W.model)),
                            all(! is.na(W.model)),
                            all(! is.nan(W.model)))
    
    if(! all(apply(W.model!=0,1,any))){
      print(W.model)	
      stop("All composites must have at least one indicator")	
    }
    if(! all(apply(W.model!=0,2,any))){
      print(W.model)	
      stop("All indicators must be linked to at least one composite")	
    }
    
    
    if(ncol(S)!=ncol(W.model)){
      print(list(S=S,W.model = W.model))
      stop("Data matrix column count does not match weight patter column count")
    }
    
    if(!variant %in% c("lohmoller","wold")){
      stop("Variant must be \"lohmoller\" or \"wold\"")
    }

    # outerEstim must be a list of same length as number of rows in inner.mod or
    # a function
    if(! is.null(outerEstim)){
      if(is.list(outerEstim)){
        assertthat::assert_that(length(outerEstim) == nrow(W.model))
        assertthat::assert_that(all(lapply(outerEstim, is.function)))
      }
      else{
        assertthat::assert_that(is.function(outerEstim))
      }
    }
    
    if(! is.null(innerEstim)){
      assertthat::assert_that(is.function(innerEstim))
    }
    # tol must be non negative
    assertthat::assert_that(tol >= 0)
    
    #iter must not be negative
    assertthat::assert_that(iter >= 0)
  }
  
  nativeModel <- parseModelToNativeFormat(model)
  inner.mod <- nativeModel$inner
  
  # If the outer estimators (tpyically Mode A and Mode B) are not defined, default to using
  # Mode A for reflective composites and Mode B for formative composites
  
  if(is.null(outerEstim)){
    hasFormativeIndicators <- any(nativeModel$formative == 1)
    hasReflectiveIndicators <- any(nativeModel$reflective == 1)
    
    if(! hasFormativeIndicators) outerEstim = outerEstim.modeA
    else if (! hasReflectiveIndicators) outerEstim = outerEstim.modeB
    else{
      # composites with at least one reflective indicator are ModeA and others are ModeB
      outerEstim <- list()
      for(composite in 1:ncol(nativeModel$reflective)){
        if(any(nativeModel$reflective[,composite] == 1)) outerEstim[[composite]] <- outerEstim.modeA
        else outerEstim[[composite]] <- outerEstim.modeB
      }
    }
  }
  
  
  ##################################################################################################
  #
  # Start of the estimation process
  #
  ##################################################################################################
  
  # The initial weight matrix
  
  weightPattern <- W.model!=0
  W <- scaleWeights(S, W.model)
  iteration <- 0
  
  weightHistory <- matrix(NA,iter+1,sum(weightPattern))
  weightHistory[1,] <- W[weightPattern]
  
  if(iter > 0) rownames(weightHistory) <- c("start",1:iter)
  else rownames(weightHistory) <- "start"
  
  # Set up outer estimators
  
  if(is.list(outerEstim)){
    uniqueOuterEstimators <- unique(outerEstim)
    outerEstimIndices <- lapply(uniqueOuterEstimators, function(x){
      sapply(outerEstim, function(y){
        identical(y,x)})
    })
  }
  
  E <- NULL
  
  # =========== Start of iterative procedure ===========
  
  # Lohmöller's variant updates all weights at the same time
  # Wold's variant updates one composite at a time
  
  if(variant =="wold"){
    compositeIndices <- 1:nrow(W)
  }
  else{
    compositeIndices <- NA
  }
  
  repeat {
    
    if(iteration == iter){
      converged <- FALSE;
      break;
    }
    
    W_old <- W
    
    # Loop over the composites or calculate all as one pass if NA
    
    for(k in compositeIndices){
      # Get new inner weights from inner estimation
      
      if(! is.null(innerEstim)){
        E <- innerEstim(S, W, inner.mod, model = model, ...)
      }
      
      # Get new weights from outer estimation
      
      # Lohmöller
      
      if(is.na(k)){
        if(is.list(outerEstim)){
          
          # Run each estimator separately
          
          for(i in 1:length(uniqueOuterEstimators)){
            W.modelForThisEstimator <- W.model
            W.modelForThisEstimator[!outerEstimIndices[[i]],] <- 0
            W[outerEstimIndices[[i]],] <- uniqueOuterEstimators[[i]](S, W_old, E, W.modelForThisEstimator,...)[outerEstimIndices[[i]],]
          }
        }
        else{
          W <- outerEstim(S, W_old, E, W.model, model = model, ...)
        }	
      }
      
      # Wold
      
      else{
        if(is.list(outerEstim)) outerEstimForThisComposite <- outerEstim[[k]]
        else outerEstimForThisComposite <- outerEstim
        
        W.modelForThisEstimator <- W.model
        W.modelForThisEstimator[-k,] <- 0
        W[W.modelForThisEstimator != 0] <- outerEstimForThisComposite(S, W_old, E, W.modelForThisEstimator,...)[W.modelForThisEstimator != 0]        
      }
      W <- scaleWeights(S, W)
    }
    
    iteration <- iteration +1 
    weightHistory[iteration+1,] <- W[weightPattern]
    
    # Check convergence. If we are not using inner estimator, converge to the first iteration
    
    if(is.null(innerEstim) || convCheck(W,W_old) < tol){
      converged <- TRUE
      break;
    }
    
  }
  
  if(!converged) warning(paste("Iterative weight algorithm did not converge."))
  
  attr(W,"S") <- S
  attr(W,"E") <- E
  attr(W,"iterations") <- iteration
  attr(W,"converged") <- converged
  attr(W,"history") <- weightHistory[1:(iteration+1),]
  class(W) <-("matrixplsweights")
  rownames(W) <- rownames(inner.mod)
  
  return(W)
}

#'@details
#'
#'\code{weightFun.optim} calculates indicator weights by optimizing the indicator
#'weights against the criterion function using \code{\link[stats]{optim}}. The
#'algorithm works by first estimating the model with the starting weights. The
#'resulting \code{matrixpls} object is passed to the \code{optimCrit}
#'function, which evaluates the optimization criterion for the weights. The
#'weights are adjusted and new estimates are calculated until the optimization
#'criterion converges.
#'
#'@inheritParams matrixpls
#'
#'@param method The minimization algorithm to be used. See \code{\link[stats]{optim}}
#' for details. Default is \code{"BFGS"}.
#'
#'@example example/matrixpls.optim-example.R
#'@describeIn weightFun calculates a set of weights to minimize an optimization criterion.
#'@export

weightFun.optim <- function(S, model, W.model,
                         parameterEstim = parameterEstim.separate, 
                         optimCrit = optimCrit.maximizeInnerR2, method = "BFGS",
                         ..., 
                         validateInput = TRUE,
                         standardize = TRUE) {
  
  W <- W.model
  
  optim.res <- stats::optim(W.model[W.model != 0], fn = function(par, ...){
    W[W.model != 0] <- par
    # Use fixed weights estimation
    matrixpls.res <- matrixpls(S, model, W, weightFun = weightFun.fixed,
                               parameterEstimator = parameterEstim.separate,
                               ..., validateInput = FALSE, standardize = standardize)
    
    optimCrit(matrixpls.res)
  }, method = method, ...)
  
  W[W.model != 0] <- optim.res$par
  W <- scaleWeights(S, W)
  
  if(optim.res$convergence) warning(paste("Weight optimization did not converge. Optim returned",optim.res$convergence))
  
  attr(W,"S") <- S
  attr(W,"iterations") <- optim.res$counts[1]
  attr(W,"converged") <- optim.res$convergence == 0
  class(W) <-("matrixplsweights")
  
  return(W)
  
}

#'@describeIn weightFun returns the starting weights.
#'@export


weightFun.fixed <- function(S, model, W.model = NULL,
                         ..., 
                         standardize = TRUE) {
  
  W <- W.model 
  
  if(standardize){
    W <- scaleWeights(S,W)
  }
  
  attr(W,"S") <- S
  class(W) <-("matrixplsweights")
  
  return(W)
  
}

#'@describeIn weightFun blockwise factor score weights.
#'@description
#'
#'\code{weightFun.factor} calculates weights by estimating a common factor analysis model with a single factor for each 
#'indicator block and using the resulting estimates to calculate factor score weights
#'
#'@param fm factoring method for estimating the common factor model. Possible values are
#'\code{minres}, \code{wls}, \code{gls}, \code{pa}, and \code{ml}. The parameter is passed through to
#' to \code{\link[psych]{fa}}.
#'
#'@export



weightFun.factor <- function(S, model, W.model = NULL, ..., fm ="minres",
                          standardize = TRUE) {
  
  # Set up a weight pattern
  W <- ifelse(W.model==0,0,1)
  
  # Do the factor analyses
  
  for(row in which(rowSums(W)>0, useNames = FALSE)){
    indicatorIndices <- W[row,]==1
    fa.res <- psych::fa(S[indicatorIndices,indicatorIndices], fm=fm)
    # Calculate the factor score weights based on the loadings and indicator covariance matrix
    W[row,indicatorIndices] <- solve(S[indicatorIndices,indicatorIndices])%*%fa.res$loading
    
  }
  
  if(standardize){
    W <- scaleWeights(S,W)
  }
  
  attr(W,"S") <- S
  class(W) <-("matrixplsweights")
  
  return(W)
  
}

#'@describeIn weightFun blockwise principal component weights.
#'
#'@description
#'
#'\code{weightFun.principal} calculates weights by calculating a principal component analysis for each 
#'indicator block and returning the weights for the first principal component.
#'
#'@export

weightFun.principal <- function(S, model, W.model = NULL, ..., 
                             standardize = TRUE) {
  
  # Set up a weight pattern
  W <- ifelse(W.model==0,0,1)
  
  # Do the factor analyses
  
  for(row in which(rowSums(W)>0, useNames = FALSE)){
    indicatorIndices <- W[row,]==1
    principal.res <- psych::principal(S[indicatorIndices,indicatorIndices], ...)
    W[row,indicatorIndices] <- principal.res$weights
    
  }
  
  if(standardize){
    W <- scaleWeights(S,W)
  }
  
  attr(W,"S") <- S
  class(W) <-("matrixplsweights")
  
  return(W)
}


#'@export

print.matrixplsweights <- function(x, ...){

  cat("\n matrixpls weights\n")

  t <- x
  class(t) <- "matrix"
  attributes(t) <- attributes(t)[c("dim", "dimnames")]
  print(t, ...)

  if(! is.null(attr(x,"converged")))
    cat("\nWeight algorithm",ifelse(attr(x,"converged"),"converged","did not converge"),"in",attr(x,"iterations"),"iterations.\n")
}
