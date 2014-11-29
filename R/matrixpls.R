# =========== Main functions ===========

#'Partial Least Squares and other composite variable models.
#'
#'Estimates a weight matrix using Partial Least Squares or a related algorithm and then
#'uses the weights to estimate the parameters of a statistical model.
#'
#'\code{matrixpls} is the main function of the matrixpls package. This function
#'parses a model object and then uses the results to call \code{weightFunction} to
#'to calculate indicator weight. After this the \code{parameterEstimator} function is
#'applied to the indicator weights, the data covariance matrix,
#'and the model object and the resulting parameter estimates are returned.
#'
#'Model can be specified in the lavaan format or the native matrixpls format.
#'The native model format is a list of three binary matrices, \code{inner}, \code{reflective},
#'and \code{formative} specifying the free parameters of a model: \code{inner} (\code{l x l}) specifies the 
#'regressions between composites, \code{reflective} (\code{k x l}) specifies the regressions of observed
#'data on composites, and \code{formative} (\code{l x k}) specifies the regressions of composites on the
#'observed data. Here \code{k} is the number of observed variables and \code{l} is the number of composites.
#'
#'If the model is specified in lavaan format, the native
#'format model is derived from this model by assigning all regressions between latent
#'variables to \code{inner}, all factor loadings to \code{reflective}, and all regressions
#'of latent variables on observed variables to \code{formative}. Regressions between
#'observed variables and all free covariances are ignored. All parameters that are
#'specified in the model will be treated as free parameters. If model is specified in
#'lavaan syntax, the model that is passed to the \code{parameterEstimator} will be that
#'model and not the native format model.
#'
#'The original papers about Partial Least Squares, as well as many of the current PLS
#'implementations, impose restrictions on the matrices \code{inner},
#'\code{reflective}, and \code{formative}: \code{inner} must be a lower triangular matrix,
#'\code{reflective} must have exactly one non-zero value on each row and must have at least
#'one non-zero value on each column, and \code{formative} must only contain zeros.
#'Some PLS implementations allow \code{formative} to contain non-zero values, but impose a
#'restriction that the sum of \code{reflective} and \code{t(formative)} must satisfy
#'the original restrictions of \code{reflective}. The only restrictions that matrixpls
#'imposes on \code{inner}, \code{reflective}, and \code{formative} is that these must be
#'binary matrices and that the diagonal of \code{inner} must be zeros.
#'
#'The argument \code{W.mod} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites. The original papers about
#'Partial Least Squares as well as all current PLS implementations define this as
#'\code{t(reflective) | formative}, which means that the weight patter must match the
#'model specified in \code{reflective} and \code{formative}. Matrixpls does not
#'require that \code{W.mod} needs to match \code{reflective} and \code{formative}, but
#'accepts any numeric matrix. If this argument is not specified, \code{W.mod} is
#'defined as \code{t(reflective) | formative}.
#'
#'@param S Covariance matrix of the data.
#
#'@param model There are two options for this argument: 1. lavaan script or lavaan parameter
#'table, or 2. a list containing three matrices
#'\code{inner}, \code{reflective}, and \code{formative} defining the free regression paths
#'in the model.
#'
#'@param W.mod An optional numeric matrix representing the weight patter and starting weights
#'(i.e. the how the indicators are combined to form the composite variables). If this argument is not specified,
#'the weight patter is defined based on the relationships in the \code{reflective} and  \code{formative}
#'elements of \code{model}.

#'@param weightFunction A function that takes three or more arguments, the data covariance matrix \code{S},
#'a model specification \code{model}, and a weight pattern \code{W.mod} and 
#' returns a named vector of parameter estimates. The default is \code{\link{weight.pls}}
#'
#'@param parameterEstimator A function that takes three  or more arguments, the data covariance matrix \code{S},
#'model specification \code{model}, and weights \code{W} and returns a named vector of parameter estimates. The default is \code{\link{params.regression}}
#'
#'@param standardize \code{TRUE} (default) or \code{FALSE} indicating whether S should be converted
#' to a correlation matrix.
#'
#'@param validateInput A boolean indicating whether the arguments should be validated.
#' 
#'@param ... All other arguments are passed through to \code{weightFunction} and \code{parameterEstimator}.
#'
#'@seealso 
#'Weight algorithms: \code{\link{weight.pls}}; \code{\link{weight.fixed}}; \code{\link{weight.optim}}
#'
#'Parameter estimators: \code{\link{params.regression}}; \code{\link{params.plsregression}}; \code{\link{params.plsc}}; \code{\link{params.tsls}}
#'
#'@return A named numeric vector of class \code{matrixpls} containing parameter estimates followed by weight.
#'
#'@references Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. \emph{Journal of Statistical Software}, 48(2), 1–36. Retrieved from http://www.jstatsoft.org/v48/i02
#'
#'Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. In K. G. Jöreskog & S. Wold (Eds.),\emph{Systems under indirect observation: causality, structure, prediction} (pp. 1–54). Amsterdam: North-Holland.
#'
#'@export


matrixpls <- function(S, model, W.mod = NULL, weightFunction = weight.pls, parameterEstimator = params.regression,
                      ..., validateInput = TRUE, standardize = TRUE) {
  
  
  # TODO: Validate input.
  
  nativeModel <- parseModelToNativeFormat(model)
  
  lvNames <- colnames(nativeModel$inner)
  
  if(lvNames != rownames(nativeModel$inner) ||
       lvNames != colnames(nativeModel$reflective) ||
       lvNames != rownames(nativeModel$formative)){
    print(nativeModel)
    stop("Names of composites are inconsistent between inner, reflective, and formative models.")
  }
  
  if(! identical(rownames(nativeModel$reflective), colnames(nativeModel$formative))){
    print(nativeModel)
    stop("Names of observed variables are inconsistent between inner, reflective, and formative models.")
  }
  
  if(! identical(colnames(S), rownames(S))) stop("S must have identical row and column names.")
  
  if(! identical(colnames(S), colnames(nativeModel$formative))){
    # Do all variables of the model exists in S
    d <- setdiff(colnames(nativeModel$formative), colnames(S))
    if(length(d)>0) stop(paste("Variable(s)",paste(d,collapse=", "),"are not included in S."))
    
    # Choose the subset that we uses. This also guarantees the correct order.
    S <- S[colnames(nativeModel$formative), colnames(nativeModel$formative)]
  }
  
  # If the weight patter is not defined, calculate it based on the model.
  if(is.null(W.mod)) W.mod <- defaultWeightModelWithModel(nativeModel)
  if(! identical(colnames(W.mod), colnames(S))) stop("Column names of W.mod do not match the variable names")
  if(! identical(rownames(W.mod), lvNames)) stop("Row names of W.mod do not match the latent variable names")
  
  ##################################################################################################
  #
  # Start of the estimation process
  #
  ##################################################################################################
  
  if(standardize){
    v <- diag(1/sqrt(diag(S)))
    S <- cov2cor(S)
    # cov2cor can result in non.symmetrix matrix because of computational inaccuracy
    S[lower.tri(S)] <- t(S)[lower.tri(S)]
  }
  
  # Calculate weight.
  
  W <- weightFunction(S, W.mod = W.mod,
                      model = nativeModel, parameterEstimator = params.regression,
                      ..., validateInput = validateInput, standardize = standardize)
  
  
  # Apply the parameter estimator and return the results
  
  estimates <- parameterEstimator(S, model, W, ...)
  
  ##################################################################################################
  #
  # Construct the return object
  #
  ##################################################################################################
  
  # Construct the return vector by combining estimates with weights in a named vector
  indices <- W.mod!=0
  WVect <- W[indices]
  
  names(WVect) <- paste(rownames(nativeModel$formative)[row(W)[indices]],"=+",colnames(nativeModel$formative)[col(W)[indices]], sep="")
  
  ret <- c(estimates,WVect)
  
  # Copy all non-standard attributes from the weight and estimation algorithms
  
  allAttributes <- c(attributes(W),attributes(estimates))
  
  for(a in setdiff(names(allAttributes), c("dim", "dimnames", "class", "names"))){
    attr(ret,a) <- allAttributes[[a]]
    attr(W,a) <- NULL
  }
  
  # Make W a plain numeric matrix instead of class matrixplweights
  class(W) <- "numeric"
  attr(ret,"W") <- W
  attr(ret,"model") <- nativeModel
  attr(ret,"call") <- match.call()
  class(ret) <-("matrixpls")
  
  return(ret)
}

#'@S3method print matrixpls

print.matrixpls <- function(x, ...){
  
  cat("\n matrixpls parameter estimates\n")
  
  indices <- ! grepl("=+", names(x), fixed=TRUE)
  toPrint <- x[indices]
  
  estimates <- data.frame(Estimate = toPrint)
  
  boot.out <- attr(x,"boot.out")
  if(! is.null(boot.out)){
    estimates$Std.Err. <- apply(boot.out$t[,indices],2,sd)
  }
  
  print(estimates, ...)
  
  if(! is.null(boot.out)){
    cat("\n Standard errors based on",boot.out$R,"bootstrap replications\n")
  }
  
  W <- attr(x,"W")
  
  attr(W,"iterations") <- attr(x,"iterations")
  attr(W,"converged") <- attr(x,"converged")
    
  class(W) <- "matrixplsweights"
  print(W, ...)
  
}

#'@S3method summary matrixpls

summary.matrixpls <- function(object, ...){
  
  ret <- list(estimates = object,
              effects = effects(object),
              R2 = R2(object),
              residuals = residuals(object),
              gof = GoF(object),
              CR = CR(object),
              AVE = AVE(object))
  
  class(ret) <-("matrixplssummary")
  
  return(ret)
}

#'@S3method print matrixplssummary

print.matrixplssummary <- function(x, ...){
  for(element in x){
    print(element, ...)
  }
}


# =========== Indicator weighting algorithms  ===========

#'@title Partial Least Squares and other iterative two-stage weight algorithms
#'
#'@description
#'Estimates a weight matrix using Partial Least Squares or a related algorithm.
#
#'@details
#'
#'\code{weight.pls} calculates indicator weights by calling the 
#'\code{innerEstimator} and \code{outerEstimators} iteratively until either the convergence criterion or
#'maximum number of iterations is reached and provides the results in a matrix.
#'
#'
#'Model can be specified in the lavaan format or the native matrixpls format.
#'The native model format is a list of three binary matrices, \code{inner}, \code{reflective},
#'and \code{formative} specifying the free parameters of a model: \code{inner} (\code{l x l}) specifies the 
#'regressions between composites, \code{reflective} (\code{k x l}) specifies the regressions of observed
#'data on composites, and \code{formative} (\code{l x k}) specifies the regressions of composites on the
#'observed data. Here \code{k} is the number of observed variables and \code{l} is the number of composites.
#'
#'If the model is specified in lavaan format, the native
#'format model is derived from this model by assigning all regressions between latent
#'variables to \code{inner}, all factor loadings to \code{reflective}, and all regressions
#'of latent variables on observed variables to \code{formative}. Regressions between
#'observed variables and all free covariances are ignored. All parameters that are
#'specified in the model will be treated as free parameters. If model is specified in
#'lavaan syntax, the model that is passed to the \code{parameterEstimator} will be that
#'model and not the native format model.
#'
#'The original papers about Partial Least Squares, as well as many of the current PLS
#'implementations, impose restrictions on the matrices \code{inner},
#'\code{reflective}, and \code{formative}: \code{inner} must be a lower triangular matrix,
#'\code{reflective} must have exactly one non-zero value on each row and must have at least
#'one non-zero value on each column, and \code{formative} must only contain zeros.
#'Some PLS implementations allow \code{formative} to contain non-zero values, but impose a
#'restriction that the sum of \code{reflective} and \code{t(formative)} must satisfy
#'the original restrictions of \code{reflective}. The only restrictions that matrixpls
#'imposes on \code{inner}, \code{reflective}, and \code{formative} is that these must be
#'binary matrices and that the diagonal of \code{inner} must be zeros.
#'
#'
#'The argument \code{W.mod} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites. The original papers about
#'Partial Least Squares as well as all current PLS implementations define this as
#'\code{t(reflective) | formative}, which means that the weight patter must match the
#'model specified in \code{reflective} and \code{formative}. Matrixpls does not
#'require that \code{W.mod} needs to match \code{reflective} and \code{formative}, but
#'accepts any numeric matrix. If this argument is not specified, \code{W.mod} is
#'defined as \code{t(reflective) | formative}.
#'
#'@inheritParams matrixpls
#'
#'@param outerEstimators A function or a list of functions used for outer estimation. If
#'the value of this parameter is a function, the same function is applied to all
#'composites. If the value of the parameter is a list, the composite \code{n} is estimated
#'with the estimator in the \code{n}th position in the list. If this argument is
#'\code{NULL} the starting weights specified in \code{W.mod} will be returned. The default is \code{\link{outer.modeA}} 
#'(PLS Mode A estimation).
#'
#'@param innerEstimator A function used for inner estimation. Setting this argument to
#'\code{null} will use identity matrix as inner estimates and causes the algorithm to converge
#'after the first iteration. This is useful when using 
#'\code{\link{outer.fixedWeights}} or some other outer estimation function that
#'does not use inner estimation results. The default is \code{\link{inner.path}}
#'(PLS path weighting scheme).

#'@param convCheck A function that takes the old and new weight matrices and
#'returns a scalar that is compared against \code{tol} to check for convergence. The default
#'function calculates the differences between each cell of old and new weights and returns the largest
#'absolute difference.
#'
#'@param tol Decimal value indicating the tolerance criterion for the
#'iterations. 
#'
#'@param iter An integer indicating the maximum number of iterations.
#'
#'@param validateInput A boolean indicating whether the validity of the parameter values should be tested.
#'
#'@param ... All other arguments are passed through to \code{outerEstimators} and \code{innerEstimator}.
#'
#'@template weights-return
#'
#'@seealso 
#'Inner estimators: \code{\link{inner.path}}; \code{\link{inner.centroid}}; \code{\link{inner.factor}}; \code{\link{inner.GSCA}}; \code{\link{inner.identity}}
#'
#'Outer estimators: \code{\link{outer.modeA}}; \code{\link{outer.modeB}}; \code{\link{outer.GSCA}}; \code{\link{outer.factor}}; \code{\link{outer.fixedWeights}}
#'
#'Convergence checks: \code{\link{convCheck.absolute}}, \code{\link{convCheck.square}}, and \code{\link{convCheck.relative}}.
#'
#'@export

weight.pls <- function(S, model, W.mod,
                       outerEstimators = NULL, 
                       innerEstimator = inner.path, ..., 
                       convCheck = convCheck.absolute,
                       tol = 1e-05, iter = 100, validateInput = TRUE) {
  
  if(validateInput){
    
    # All parameters must have values
    assertive::assert_all_are_not_na(formals())
    
    # S must be symmetric and a valid covariance matrix
    assertive::assert_is_matrix(S)
    assertive::assert_is_symmetric_matrix(S)
    assertive::assert_is_identical_to_true(is.positive.semi.definite(S))
    
    # W.mod must be a real matrix and each indicators must be
    # linked to at least one composite and each composite at least to one indicator
    assertive::assert_is_matrix(W.mod)
    assertive::assert_all_are_real(W.mod)
    
    if(! all(apply(W.mod!=0,1,any))){
      print(W.mod)	
      stop("All constructs must have at least one indicator")	
    }
    if(! all(apply(W.mod!=0,2,any))){
      print(W.mod)	
      stop("All indicators must be linked to at least one construct")	
    }
    
    
    if(ncol(S)!=ncol(W.mod)){
      print(list(S=S,W.mod = W.mod))
      stop("Data matrix column count does not match weight patter column count")
    }
    
    # outerEstimators must be a list of same length as number of rows in inner.mod or
    # a function
    if(! is.null(outerEstimators)){
      if(is.list(outerEstimators)){
        assertive::assert_is_identical_to_true(length(outerEstimators) == nrow(W.mod))
        for(outerEstimator in outerEstimators){
          assertive::assert_is_function(outerEstimator)
        }
      }
      else{
        assertive::assert_is_function(outerEstimators)
      }
    }
    
    if(! is.null(innerEstimator)){
      assertive::assert_is_function(innerEstimator)
    }
    # tol must be non negative
    assertive::assert_all_are_non_negative(tol)
    
    #iter must not be negative
    assertive::assert_all_are_non_negative(iter)
  }
  
  nativeModel <- parseModelToNativeFormat(model)
  inner.mod <- nativeModel$inner
  
  # If the outer estimators (tpyically Mode A and Mode B) are not defined, default to using
  # Mode A for reflective constructs and Mode B for formative constructs
  
  if(is.null(outerEstimators)){
    hasFormativeIndicators <- any(nativeModel$formative == 1)
    hasReflectiveIndicators <- any(nativeModel$reflective == 1)
    
    if(! hasFormativeIndicators) outerEstimators = outer.modeA
    else if (! hasReflectiveIndicators) outerEstimators = outer.modeB
    else{
      # Constructs with at least one reflective indicator are ModeA and others are ModeB
      outerEstimators <- list()
      for(construct in 1:ncol(nativeModel$reflective)){
        if(any(nativeModel$reflective[,construct] == 1)) outerEstimators[[construct]] <- outer.modeA
        else outerEstimators[[construct]] <- outer.modeB
      }
    }
  }
  
  
  ##################################################################################################
  #
  # Start of the estimation process
  #
  ##################################################################################################
  
  # The initial weight matrix
  
  weightPattern <- W.mod!=0
  W <- scaleWeights(S, W.mod)
  iteration <- 0
  
  weightHistory <- matrix(NA,iter+1,sum(weightPattern))
  weightHistory[1,] <- W[weightPattern]
  
  if(iter > 0) rownames(weightHistory) <- c("start",1:iter)
  else rownames(weightHistory) <- "start"
  
  
  # If we are not using an outer estimator, do not perform iterative estimation
  
  if(!is.null(outerEstimators)){
    
    # Set up outer estimators
    
    if(is.list(outerEstimators)){
      uniqueOuterEstimators <- unique(outerEstimators)
      outerEstimatorIndices <- lapply(uniqueOuterEstimators, function(x){
        sapply(outerEstimators, function(y){
          identical(y,x)})
      })
    }
    
    E <- NULL
    
    # =========== Start of iterative procedure ===========
    
    repeat {
      
      if(iteration == iter){
        converged <- FALSE;
        break;
      }
      
      # Get new inner weights from inner estimation
      
      if(! is.null(innerEstimator)){
        E <- innerEstimator(S, W, inner.mod, model = model, ...)
      }
      
      # Get new weights from outer estimation
      
      W_old <- W
      if(is.list(outerEstimators)){
        
        # Run each estimator separately
        
        for(i in 1:length(uniqueOuterEstimators)){
          W.modForThisEstimator <- W.mod
          W.modForThisEstimator[!outerEstimatorIndices[[i]],] <- 0
          W[outerEstimatorIndices[[i]],] <- uniqueOuterEstimators[[i]](S, W_old, E, W.modForThisEstimator,...)[outerEstimatorIndices[[i]]]
        }
      }
      else{
        W <- outerEstimators(S, W_old, E, W.mod, model = model, ...)
      }	
      
      W <- scaleWeights(S, W)
      
      iteration <- iteration +1 
      weightHistory[iteration+1,] <- W[weightPattern]
      
      # Check convergence. If we are not using inner estimator, converge to the first iteration
      
      if(is.null(innerEstimator) || convCheck(W,W_old) < tol){
        converged <- TRUE
        break;
      }
      
    }
  }
  
  # Mark the estimation as converged if not running iterations
  
  else converged <- TRUE	
  
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

#'@title Optimized weights
#'
#'@description
#'Calculates a set of weights to minimize an optimization criterion
#
#'@details
#'
#'\code{weight.optim} calculates indicator weights by optimizing the indicator
#'weights against the criterion function using \code{\link[stats]{optim}}. The
#'algoritmh works by first estimating the model with the starting weights. The
#'resulting \code{matrixpls} object is passed to the \code{optimCriterion}
#'function, which evaluates the optimization criterion for the weights. The
#'weights are adjusted and new estimates are calculated until the optimization
#'criterion converges.
#'
#'
#'Model can be specified in the lavaan format or the native matrixpls format.
#'The native model format is a list of three binary matrices, \code{inner}, \code{reflective},
#'and \code{formative} specifying the free parameters of a model: \code{inner} (\code{l x l}) specifies the 
#'regressions between composites, \code{reflective} (\code{k x l}) specifies the regressions of observed
#'data on composites, and \code{formative} (\code{l x k}) specifies the regressions of composites on the
#'observed data. Here \code{k} is the number of observed variables and \code{l} is the number of composites.
#'
#'If the model is specified in lavaan format, the native
#'format model is derived from this model by assigning all regressions between latent
#'variables to \code{inner}, all factor loadings to \code{reflective}, and all regressions
#'of latent variables on observed variables to \code{formative}. Regressions between
#'observed variables and all free covariances are ignored. All parameters that are
#'specified in the model will be treated as free parameters. If model is specified in
#'lavaan syntax, the model that is passed to the \code{parameterEstimator} will be that
#'model and not the native format model.
#'
#'The original papers about Partial Least Squares, as well as many of the current PLS
#'implementations, impose restrictions on the matrices \code{inner},
#'\code{reflective}, and \code{formative}: \code{inner} must be a lower triangular matrix,
#'\code{reflective} must have exactly one non-zero value on each row and must have at least
#'one non-zero value on each column, and \code{formative} must only contain zeros.
#'Some PLS implementations allow \code{formative} to contain non-zero values, but impose a
#'restriction that the sum of \code{reflective} and \code{t(formative)} must satisfy
#'the original restrictions of \code{reflective}. The only restrictions that matrixpls
#'imposes on \code{inner}, \code{reflective}, and \code{formative} is that these must be
#'binary matrices and that the diagonal of \code{inner} must be zeros.
#'
#'
#'The argument \code{W.mod} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites. The original papers about
#'Partial Least Squares as well as all current PLS implementations define this as
#'\code{t(reflective) | formative}, which means that the weight patter must match the
#'model specified in \code{reflective} and \code{formative}. Matrixpls does not
#'require that \code{W.mod} needs to match \code{reflective} and \code{formative}, but
#'accepts any numeric matrix. If this argument is not specified, \code{W.mod} is
#'defined as \code{t(reflective) | formative}.
#'
#'@inheritParams matrixpls
#'
#'@param method The minimization algorithm to be used. See \code{\link[stats]{optim}}
#' for details. Default is "BFGS".
#' 
#'@param optimCriterion A function that taking an object of class class 
#'\code{matrixpls} and returning a scalar. The default is \code{\link{optim.maximizeInnerR2}}
#'
#'@param ... All other arguments are passed through to \code{\link[stats]{optim}} and \code{parameterEstimator}.
#'
#'@template weights-return
#'
#'@seealso 
#'Optimization criteria: \code{\link{optim.maximizeInnerR2}}
#'
#'@example example/matrixpls.optim-example.R
#'
#'@export



weight.optim <- function(S, model, W.mod,
                         parameterEstimator = params.regression, 
                         optimCriterion = optim.maximizeInnerR2, method = "BFGS",
                         ..., 
                         validateInput = TRUE,
                         standardize = TRUE) {
  
  W <- W.mod
  
  optim.res <- optim(W.mod[W.mod != 0], fn = function(par, ...){
    W[W.mod != 0] <- par
    # Use fixed weights estimation
    matrixpls.res <- matrixpls(S, model, W, weightFunction = weight.fixed,
                               parameterEstimator = params.regression,
                               ..., validateInput = FALSE, standardize = standardize)
    
    optimCriterion(matrixpls.res)
  }, method = method, ...)
  
  W[W.mod != 0] <- optim.res$par
  W <- scaleWeights(S, W)
  
  if(optim.res$convergence) warning(paste("Weight optimization did not converge. Optim returned",optim.res$convergence))
  
  attr(W,"S") <- S
  attr(W,"iterations") <- optim.res$counts[1]
  attr(W,"converged") <- optim.res$convergence == 0
  attr(W,"history") <- NA
  class(W) <-("matrixplsweights")
  
  return(W)
  
}

#'@title Fixed weights
#'
#'@description
#'Returns fixed weights
#
#'@details
#'
#'Returns the starting weights specified
#'in \code{W.mod}. If \code{standardize} is \code{TRUE} the weights are
#'standardized so that the composites have unit variances.
#'
#'
#'Model can be specified in the lavaan format or the native matrixpls format.
#'The native model format is a list of three binary matrices, \code{inner}, \code{reflective},
#'and \code{formative} specifying the free parameters of a model: \code{inner} (\code{l x l}) specifies the 
#'regressions between composites, \code{reflective} (\code{k x l}) specifies the regressions of observed
#'data on composites, and \code{formative} (\code{l x k}) specifies the regressions of composites on the
#'observed data. Here \code{k} is the number of observed variables and \code{l} is the number of composites.
#'
#'If the model is specified in lavaan format, the native
#'format model is derived from this model by assigning all regressions between latent
#'variables to \code{inner}, all factor loadings to \code{reflective}, and all regressions
#'of latent variables on observed variables to \code{formative}. Regressions between
#'observed variables and all free covariances are ignored. All parameters that are
#'specified in the model will be treated as free parameters. If model is specified in
#'lavaan syntax, the model that is passed to the \code{parameterEstimator} will be that
#'model and not the native format model.
#'
#'The original papers about Partial Least Squares, as well as many of the current PLS
#'implementations, impose restrictions on the matrices \code{inner},
#'\code{reflective}, and \code{formative}: \code{inner} must be a lower triangular matrix,
#'\code{reflective} must have exactly one non-zero value on each row and must have at least
#'one non-zero value on each column, and \code{formative} must only contain zeros.
#'Some PLS implementations allow \code{formative} to contain non-zero values, but impose a
#'restriction that the sum of \code{reflective} and \code{t(formative)} must satisfy
#'the original restrictions of \code{reflective}. The only restrictions that matrixpls
#'imposes on \code{inner}, \code{reflective}, and \code{formative} is that these must be
#'binary matrices and that the diagonal of \code{inner} must be zeros.
#'
#'
#'The argument \code{W.mod} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites. The original papers about
#'Partial Least Squares as well as all current PLS implementations define this as
#'\code{t(reflective) | formative}, which means that the weight patter must match the
#'model specified in \code{reflective} and \code{formative}. Matrixpls does not
#'require that \code{W.mod} needs to match \code{reflective} and \code{formative}, but
#'accepts any numeric matrix. If this argument is not specified, \code{W.mod} is
#'defined as \code{t(reflective) | formative}.
#'
#'@inheritParams matrixpls
#'
#'@param ... All other arguments are ignored.
#'
#'@template weights-return
#'
#'@export


weight.fixed <- function(S, model, W.mod = NULL,
                         ..., 
                         standardize = TRUE) {
  
  W <- W.mod 
  
  if(standardize){
    W <- scaleWeights(S,W)
  }

  attr(W,"S") <- S
  attr(W,"iterations") <- 0
  attr(W,"converged") <- TRUE
  attr(W,"history") <- W
  class(W) <-("matrixplsweights")
  
  return(W)
  
}

#'@S3method print matrixplsweights

print.matrixplsweights <- function(x, ...){
  cat("\n matrixpls weights\n")
  print.table(x, ...)
  cat("\nWeight algorithm",ifelse(attr(x,"converged"),"converged","did not converge"),"in",attr(x,"iterations"),"iterations.\n")
}

# =========== Convergence checks for weight functions ===========

#'@title Relative difference convergence check criterion
#'
#'@description
#'Compares two weight matrices and returns the value of the relative 
#'difference convergence criterion.
#'
#'@details
#'Relative difference convergence criterion is the maximum of relative differences between
#'weights from two iterations
#'
#'@param Wnew Current weight matrix
#'@param Wold Weight matrix from the previous iteration
#'
#'@return Scalar value of the convergence criterion
#'
#'@family Convergence checks

convCheck.relative <- function(Wnew, Wold){
  max(abs((Wold[Wnew != 0]-Wnew[Wnew != 0])/Wnew[Wnew != 0]))
}

#'@title Squared difference convergence check criterion
#'
#'@description
#'Compares two weight matrices and returns the value of the square difference
#'convergence criterion.
#'
#'@details
#'Squared difference convergence criterion is the maximum of squared absolute 
#'differences between weights from two iterations
#'
#'@inheritParams convCheck.relative
#'
#'@return Scalar value of the convergence criterion
#'
#'@family Convergence checks
#'

convCheck.square <- function(Wnew, Wold){
  max((Wold-Wnew)^2)
}

#'@title Absolute difference convergence check criterion
#'
#'@description
#'Compares two weight matrices and returns the value of the absolute
#'difference convergence criterion.
#'
#'@details
#'Absolute difference convergence criterion is the maximum of absolute differences
#'between weights from two iterations
#'
#'@inheritParams convCheck.relative
#'
#'@return Scalar value of the convergence criterion
#'
#'@family Convergence checks
#'

convCheck.absolute <- function(Wnew, Wold){
  max(abs(Wold-Wnew))
}


# =========== Parameter estimators ===========

#'@title Parameter estimation with separate regression analyses
#'
#'@description
#'Estimates the model parameters with weighted composites using separate OLS regressions.
#'
#'@details
#'\code{params.regression} estimates the statistical model described by \code{model} with the
#'following steps. If \code{model} is not in the native format, it is converted to the native
#'format containing matrices \code{inner}, \code{reflective}, and \code{formative}. The
#'weights \code{W} and the data covariance matrix \code{S} are used to calculate the composite
#'covariance matrix \code{C} and the indicator-composite covariance matrix \code{IC}. These
#'are used to estimate multiple OLS regression models.
#'
#'The OLS regressions are estimated separately for each of the three model parts \code{inner},
#'\code{reflective}, and \code{formative}. These matrices are analyzed one row at a time
#'so that the row specifies the index of the dependent variable in the OLS regression and 
#'the non-zero elements on the row specify the indices of the independent variables.
#'
#'This approach of estimating the inner and outer models separately with separate 
#'OLS regression analyses is the standard way of estimation in the PLS literature.
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@inheritParams matrixpls
#'
#'@return A named vector of parameter estimates.
#'
#'@family parameter estimators
#'
#'@export

params.regression <- function(S, model, W, ...){
  
  return(params.internal_generic(S,model, W,
                                 regressionsWithCovarianceMatrixAndModelPattern))
}


#'@title Parameter estimation with two-stage least squares
#'
#'@description
#'Estimates the model parameters with weighted composites using two-stage least squares
#'
#'@details
#'\code{params.tsls} estimates the statistical model described by \code{model} with the
#'following steps. If \code{model} is not in the native format, it is converted to the native
#'format containing matrices \code{inner}, \code{reflective}, and \code{formative}. The
#'weights \code{W} and the data covariance matrix \code{S} are used to calculate the composite
#'covariance matrix \code{C} and the indicator-composite covariance matrix \code{IC}. These
#'are used to estimate multiple regression models with two-stage least squares for
#'inner model and OLS regressions for the outer model.
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@inheritParams matrixpls
#'
#'@return A named vector of parameter estimates.
#'
#'@family parameter estimators
#'
#'@export

params.tsls <- function(S, model, W, ...){
  
  return(params.internal_generic(S,model, W, 
                                 TwoStageLeastSquaresWithCovarianceMatrixAndModelPattern))
}

params.internal_generic <- function(S, model, W, pathEstimator){
  
  nativeModel <- parseModelToNativeFormat(model)
  
  results <- c()
  
  # Calculate the composite covariance matrix
  C <- W %*% S %*% t(W)
  
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  
  innerRegressionIndices <- which(nativeModel$inner==1, useNames = FALSE)
  
  # Choose the specified values and add names
  if(length(innerRegressionIndices)>0){
    inner <- pathEstimator(C,nativeModel$inner)
    innerVect <- inner[innerRegressionIndices]
    names(innerVect) <- paste(rownames(inner)[row(inner)[innerRegressionIndices]],"~",
                              colnames(inner)[col(inner)[innerRegressionIndices]], sep="")
    
    results <- c(results, innerVect)
  }
  else{
    inner <- matrix(0,nrow(nativeModel$inner), ncol(nativeModel$inner))
  }
  
  results <- c(results, params.internal_reflective(C, IC, nativeModel))
  results <- c(results,	params.internal_formative(S, IC, nativeModel))
  
  # Store these in the result object
  attr(results,"C") <- C
  attr(results,"IC") <- IC
  attr(results,"beta") <- inner
  
  return(results)
}
params.internal_formative <- function(S, IC, nativeModel){
  
  results <-c()
  
  for(row in 1:nrow(nativeModel$formative)){
    
    independents <- which(nativeModel$formative[row,]!=0, useNames = FALSE)
    
    if(length(independents)>0){
      if(length(independents)==1){
        # Simple regresion is the covariance divided by the variance of the predictor, which are standardized
        thisResults <- IC[row,independents]
      }
      else{
        thisResults <- solve(S[independents,independents],IC[row,independents])
      }
      
      names(thisResults) <- paste(rownames(nativeModel$formative)[row], "<~",
                                  names(thisResults), sep = "")
      
      results <- c(results, thisResults)
    }
  }
  return(results)
}

params.internal_reflective <- function(C, IC, nativeModel){	
  
  results <-c()
  
  for(row in 1:nrow(nativeModel$reflective)){
    
    independents <- which(nativeModel$reflective[row,]!=0, useNames = FALSE)
    
    if(length(independents)>0){
      if(length(independents)==1){
        # Simple regresion is the covariance divided by the variance of the predictor, which are standardized
        coefs <- IC[independents,row]
      }
      else{
        coefs <- solve(C[independents,independents],IC[independents,row])
      }
      
      names(coefs) <- paste(colnames(nativeModel$reflective)[independents], "=~",
                            rownames(nativeModel$reflective)[row], sep = "")
      
      results <- c(results,coefs)
      
    }
  }
  
  return(results)
}


# =========== Inner estimators ===========

#'@title PLS inner estimation with the centroid scheme
#'
#'@description
#'Calculates a set of inner weights based on the centroid scheme. 
#
#'@details
#'In the centroid scheme, inner weights are set to the signs (1 or -1) of correlations between
#'composites that are connected in the model specified in \code{inner.mod} and zero otherwise.
#'
#'Falls back to to identity scheme for composites that are not connected to any other composites.
#'
#'@inheritParams inner.factor
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@family inner estimators
#'
#'@references
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.

#'@export

inner.centroid <- function(S, W, inner.mod, ignoreInnerModel = FALSE, ...){
  
  # Centroid is just the sign of factor weighting
  
  E <- sign(inner.factor(S, W, inner.mod, ignoreInnerModel, ...))
  
  return(E)
}

#'@title PLS inner estimation with the path scheme
#'
#'@description
#'Calculates a set of inner weights based on the path scheme.
#
#'@details
#'In the path scheme, inner weights are based on regression estimates of the relationships between
#'composites that are connected in the model specified in \code{inner.mod}, and correlations for
#'the inverse relationships. If a relationship is reciprocal, regression is used for both directions.
#'
#'Falls back to to identity scheme for composites that are not connected to any other composites.
#'
#'@inheritParams inner.centroid
#'
#'@family inner estimators
#'
#'@references
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.

#'@export

inner.path <- function(S, W, inner.mod, ...){
  
  # Calculate the composite covariance matrix	
  C <- W %*% S %*% t(W)
  
  E <- regressionsWithCovarianceMatrixAndModelPattern(C, inner.mod)
  
  # Use correlations for inverse relationships for non-reciprocal paths
  inverseRelationships <- t(inner.mod) & ! inner.mod
  E[inverseRelationships] <- C[inverseRelationships]
  
  # If we have LVs that are not connected to any other LVs, use identity scheme as fallback
  diag(E)[rowSums(E) == 0] <- 1
  
  return(E)
}

#'@title PLS inner estimation with the factor scheme
#'
#'@description
#'Calculates a set of inner weights based on the factor scheme.
#
#'@details
#'In the factor scheme, inner weights are set to the correlations between
#'composites that are connected in the model specified in \code{inner.mod} and zero otherwise.
#'
#'Falls back to to identity scheme for composites that are not connected to any other composites.
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param inner.mod A square matrix specifying the relationships of the composites in the model.

#'@param ignoreInnerModel Should the inner model be ignored and all correlations be used.
#'
#'@param ... Other parameters are ignored
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@family inner estimators
#'
#'@references
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.

#'@export

inner.factor <- function(S, W, inner.mod, ignoreInnerModel = FALSE, ...){
  
  # Calculate the composite covariance matrix
  
  C <- W %*% S %*% t(W)
  
  # Calculate unscaled inner weights using the factor weighting scheme
  if(ignoreInnerModel){
    E <- C
    diag(E) <- 0
  } 
  else E <- C * (inner.mod | t(inner.mod))
  
  # If we have LVs that are not connected to any other LVs, use identity scheme as fallback
  diag(E)[rowSums(E) == 0] <- 1
  
  return(E)
} 

#'@title PLS inner estimation with the Horst scheme
#'
#'@description
#'Calculates a set of inner weights based on the Horst scheme.
#
#'@details
#'In the Hosrt scheme, inner weights are unit weights for
#'composites that are connected in the model specified in \code{inner.mod} and zero otherwise.
#'
#'Falls back to to identity scheme for composites that are not connected to any other composites.
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param inner.mod A square matrix specifying the relationships of the composites in the model.

#'@param ignoreInnerModel Should the inner model be ignored and all correlations be used.
#'
#'@param ... Other parameters are ignored
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@family inner estimators
#'
#'@references
#'Tenenhaus, A., & Tenenhaus, M. (2011). Regularized Generalized Canonical 
#'Correlation Analysis. Psychometrika, 76(2), 257–284. doi:10.1007/s11336-011-9206-8

#'@export

inner.Horst <- function(S, W, inner.mod, ignoreInnerModel = FALSE, ...){
  
  if(ignoreInnerModel){
    E <- lower.tri(inner.mod) * 1
  } 
  else E <- inner.mod 
  
  # If we have LVs that are not connected to any other LVs, use identity scheme as fallback
  diag(E)[rowSums(E) == 0] <- 1
  
  return(E)
} 


#'@title PLS inner estimation with the identity scheme
#'
#'@description
#'
#'Returns a set of inner weights using the identity scheme. 
#'
#'@details
#'
#'This scheme is not commonly discussed in the current PLS literature, but it is a special case 
#'that is analyzed in the early PLS literature and currently implemented in at least the WarpPLS software.
#'In the identity scheme identity matrix is used as the inner weight matrix \code{E}.
#'
#'The identity scheme with Mode A outer estimation converges toward the first principal component of each indicator block.
#'
#'@inheritParams inner.centroid
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@references
#'
#'Wold, H. (1966). Nonlinear estimation by iterative least squares procedures. \emph{Research Papers in Statistics: Festschrift for J. Neyman}, 411–444.
#'
#'Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. In K. G. Jöreskog & S. Wold (Eds.), \emph{Systems under indirect observation: causality, structure, prediction} (pp. 1–54). Amsterdam: North-Holland.
#'
#'@family inner estimators
#'
#'@export



inner.identity <- function(S, W, inner.mod, ...){
  return(diag(nrow(inner.mod)))
}

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



# =========== Outer estimators ===========

#'@title PLS outer estimation with Mode A
#'
#'@description
#'
#'Performs Mode A outer estimation by calculating correlation matrix between the indicators and composites.
#
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param E Inner weight matrix. A square matrix of inner estimates between the composites.
#'
#'@param W.mod A matrix specifying the weight relationships and their starting values.
#'
#'@inheritParams inner.centroid
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@family outer estimators  
#'
#'@export


outer.modeA <- function(S, W, E, W.mod, ...){
  
  # Calculate the covariance matrix between indicators and composites
  W_new <- E %*% W %*% S
  
  # Set the non-existing weight relations to zero
  W_new[W.mod == 0] <- 0
  
  return(W_new)
}

#'@title PLS outer estimation with Mode B
#'
#'@description
#'
#'Performs Mode B outer estimation by regressing the composites on the observed variables.
#'
#'@inheritParams outer.modeA
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@family outer estimators
#'
#'@export


outer.modeB <- function(S, W, E, W.mod, ...){
  
  # Calculate the covariance matrix between indicators and composites
  IC <- E %*% W %*% S
  
  # Set up a weight pattern
  W_new <- ifelse(W.mod==0,0,1)
  
  # Do the outer model regressions
  
  for(row in which(rowSums(W_new)>0, useNames = FALSE)){
    indicatorIndices <- W_new[row,]==1
    W_new[row,indicatorIndices] <- solve(S[indicatorIndices,indicatorIndices],IC[row,indicatorIndices])
  }
  
  return(W_new)
  
}

#'@title PLS outer estimation with fixed weights
#'
#'@description
#'
#'Returns the starting weights as outer estimates.
#'
#'@inheritParams outer.modeA
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@family outer estimators
#'
#'@export

outer.fixedWeights <- function(S, W, E, W.mod, ...){
  return(W.mod)
}


#'@title RGCCA outer estimation (Experimental)
#'
#'@description
#'
#'Implements the second step of the RGCCA estimation describe by Tenenhaus & Tenehaus (2011)
#'
#'@details
#'
#'The second step of RGCCA estimation method, as describe by Tenenhaus & Tenehaus (2011)
#'@inheritParams outer.modeA
#'
#'@param tau The shrinkage parameter (0-1) 
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@references
#'Tenenhaus, A., & Tenenhaus, M. (2011). Regularized Generalized Canonical 
#'Correlation Analysis. Psychometrika, 76(2), 257–284. 
#'doi:10.1007/s11336-011-9206-8
#'
#'@example example/matrixpls.RGCCA-example.R
#'
#'@family outer estimators
#'
#'@export


outer.RGCCA <- function(S, W, E, W.mod, ..., tau = 1){
  
  if(length(tau) == 1) tau <- rep(tau, nrow(W.mod))
  else if(length(tau) != nrow(W.mod)) stop("Parameter tau of outer.RGCCA must be of length one or the number of composites.")
  
  # Calculate the covariance matrix between indicators and composites
  IC <- E %*% W %*% S
  
  # Set up a weight pattern
  W_new <- ifelse(W.mod==0,0,1)
  
  # Do the outer model regressions
  
  for(row in which(rowSums(W_new)>0, useNames = FALSE)){
    
    indicatorIndices <- W_new[row,]==1
    # Tikhonov matrix
    Gamma <- tau[row] * diag(sum(indicatorIndices))
    
    # Ridge regression
    W_new[row,indicatorIndices] <- solve(S[indicatorIndices,indicatorIndices] + t(Gamma) %*% Gamma,IC[row,indicatorIndices])
    
  }
  
  return(W_new)
  
}

#'@title Blockwise factor score outer estimation
#'
#'@description
#'
#'Calculates outer weights by estimating a common factor analysis model with a single factor for each 
#'indicator block and using the resulting estimates to calculate factor score weights
#'
#'
#'@inheritParams outer.modeA
#'
#'@param fm factoring method for estimating the common factor model. Possible values are
#'\code{minres}, \code{wls}, \code{gls}, \code{pa}, and \code{ml}. The parameter is passed through to
#' to \code{\link[psych]{fa}}.

#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@family outer estimators
#'
#'@export

outer.factor <- function(S, W, E, W.mod, fm="minres", ...){
  
  # Set up a weight pattern
  W_new <- ifelse(W.mod==0,0,1)
  
  # Do the factor analyses
  
  for(row in which(rowSums(W_new)>0, useNames = FALSE)){
    indicatorIndices <- W_new[row,]==1
    fa.res <- fa(S[indicatorIndices,indicatorIndices], fm = fm)
    # Calculate the factor score weights based on the loadings and indicator covariance matrix
    W_new[row,indicatorIndices] <- solve(S[indicatorIndices,indicatorIndices])%*%fa.res$loading
    
  }
  
  return(W_new)
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
    
    #
    # In Hwang, H., & Takane, Y. (2004), the weights for one composite are
    # defined by specifying a series of regression analyses where the
    # indicators are independent variables. These regressions are then estimated
    # simulataneously by stacking the data so that the system of equations
    # can be estimated with OLS estimator in one go.
    #
    # This is equivalent to collecting all covariances between the IVs and
    # DVs into a matrix and then taking a mean over all DVs so that we 
    # have a vector of mean covariances for each IV. 
    #
    
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
    
    # Indicators that depend on this composite
    
    for(j in which(reflective[,row]!=0)){
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


# =========== Optimization criterion functions ===========

#'@title Optimization criterion to maximize the inner model mean R2
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@return Negative mean R2.
#'
#'@family Weight optimization criteria
#'
#'@export
#'

optim.maximizeInnerR2 <- function(matrixpls.res){
  -mean(R2(matrixpls.res))
}

#'@title Optimization criterion for maximal prediction
#'
#'@details Calculates the predicted variances of reflective indicators. The
#'prediction criterion is negative of the sum of predicted variances.
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@return Mean squared prediction error.
#'
#'@family Weight optimization criteria
#'
#'@export
#'

optim.maximizePrediction <- function(matrixpls.res){
  
  C <- attr(matrixpls.res,"C")
  nativeModel <- attr(matrixpls.res,"model")
  exog <- rowSums(nativeModel$inner)==0
  W <- attr(matrixpls.res,"W")
  beta <- attr(matrixpls.res,"beta")
  
  reflective <- nativeModel$reflective
  reflective[which(reflective==1)] <- matrixpls.res[grep("=~",names(matrixpls.res))]
  
  # Predict endog LVs using reduced from equations
  Beta <- beta[! exog, ! exog]
  Gamma <- beta[! exog, exog]
  
  invBeta <- solve(diag(nrow(Beta)) - Beta)
  endogC <- invBeta %*% Gamma %*% C[exog,exog] %*% t(Gamma) %*% t(invBeta)
  C[!exog, !exog] <- endogC
  
  -sum(diag(reflective %*% C %*% t(reflective)))
}


#'@title GSCA optimization criterion
#'
#'@details Optimization criterion for maximixing the weighted sum of composite
#'correlations
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@inheritParams weight.pls
#'
#'@return Sum of residual variances.
#'
#'@family Weight optimization criteria
#'
#'@export
#'
#'@references
#'
#'Tenenhaus, A., & Tenenhaus, M. (2011). Regularized Generalized Canonical 
#'Correlation Analysis. Psychometrika, 76(2), 257–284. doi:10.1007/s11336-011-9206-8
#'

optim.GCCA <- function(matrixpls.res, innerEstimator, ...){
  
  C <- attr(matrixpls.res,"C")
  S <- attr(matrixpls.res,"S")
  W <- attr(matrixpls.res,"W")
  inner <- attr(matrixpls.res,"model")$inner.mod
  
  E <- function(S, W, inner.mod, ...)
  
  -sum(C * E)
  
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


# =========== Utility functions ===========

#
# Scales the weight matrix so that the resulting composites have unit variance
#

scaleWeights <- function(S, W){
  
  # Calculate the variances of the unscaled composites
  
  var_C_unscaled <- diag(W %*% S %*% t(W))
  
  # Scaling the unscaled weights and return
  
  ret <- diag(x = 1/sqrt(var_C_unscaled)) %*% W
  
  colnames(ret) <- colnames(W)
  rownames(ret) <- rownames(W)
  
  return(ret)
}

lavaanParTableToNativeFormat <- function(partable){
  
  factorLoadings <- partable[partable$op == "=~",]
  regressions <- partable[partable$op == "~",]
  formativeLoadings <- partable[partable$op == "<~",]
  
  # Parse the variables
  latentVariableNames <- unique(c(factorLoadings$lhs, formativeLoadings$lhs))
  observedVariableNames <- setdiff(unique(c(partable$rhs,partable$lhs)), latentVariableNames)
  
  # Set up empty model tables
  
  inner <- matrix(0,length(latentVariableNames),length(latentVariableNames))
  colnames(inner)<-rownames(inner)<-latentVariableNames
  
  formative <- matrix(0,length(latentVariableNames),length(observedVariableNames))
  colnames(formative) <- observedVariableNames
  rownames(formative) <- latentVariableNames
  
  reflective <- t(formative)
  
  # Set the relationships in the tables
  
  rows <- match(factorLoadings$rhs, observedVariableNames)
  cols <- match(factorLoadings$lhs,latentVariableNames)
  indices <- rows + (cols-1)*nrow(reflective) 
  
  reflective[indices] <- 1
  
  latentRegressions <- regressions[regressions$rhs %in% latentVariableNames & 
                                     regressions$lhs %in% latentVariableNames,]
  
  rows <- match(regressions$lhs, latentVariableNames)
  cols <- match(regressions$rhs,latentVariableNames)
  indices <- rows + (cols-1)*nrow(inner) 
  
  inner[indices] <- 1
  
  formativeRegressions <- rbind(regressions[regressions$rhs %in% observedVariableNames & 
                                              regressions$lhs %in% latentVariableNames,],
                                formativeLoadings)
  
  
  rows <- match(formativeRegressions$lhs, latentVariableNames)
  cols <- match(formativeRegressions$rhs, observedVariableNames)
  indices <- rows + (cols-1)*nrow(formative) 
  
  formative[indices] <- 1
  
  return(list(inner=inner, reflective=reflective, formative=formative))
}

parseModelToNativeFormat <- function(model){
  
  if(is.matrixpls.model(model)){
    #Already in native format
    return(model)
  }
  else if(is.character(model)) {
    # Remove all multigroup specifications because we do not support multigroup analyses
    model <- gsub("c\\(.+?\\)","NA",model)
    
    return(lavaanParTableToNativeFormat(lavaanify(model)))
  } else if (is.partable(model)) {
    return(lavaanParTableToNativeFormat(model))
  } else if (is.lavaancall(model)) {
    browser()
  } else if (is(model, "lavaan")) {
    return(lavaanParTableToNativeFormat(model@ParTable))
    #	} else if (is(model, "MxModel")) {
    #		browser()
    #	} else if (is(model, "SimSem")) {
    #		browser()
  } else {
    stop("Please specify an appropriate object for the 'model' argument: simsem model template, lavaan script, lavaan parameter table, list of options for the 'lavaan' function, or a list containing three matrices in the matrixpls native model format.")
  }
}

#
# Returns a matrix with composites on rows and observed on columns
#

defaultWeightModelWithModel <- function(model){
  
  nativeModel <- parseModelToNativeFormat(model)
  W.mod <- matrix(0,nrow(nativeModel$formative), ncol(nativeModel$formative))
  
  colnames(W.mod) <- colnames(nativeModel$formative)
  rownames(W.mod) <- rownames(nativeModel$formative)
  
  W.mod[nativeModel$formative!=0] <- 1 
  W.mod[t(nativeModel$reflective)!=0] <- 1 
  
  return(W.mod)
  
}

is.matrixpls.model <- function(model) {
  return (is.list(model) &&
            setequal(names(model),c("inner","reflective","formative")) &&
            all(sapply(model,is.matrix)))
  
}

#
# Runs regressions defined by model
# using covariance matrix S and places the results in model
#

regressionsWithCovarianceMatrixAndModelPattern <- function(S,model){
  
  # Ensure that S and model have right order
  S <- S[rownames(model),colnames(model)]
  
  for(row in 1:nrow(model)){
    
    independents <- which(model[row,]!=0, useNames = FALSE)
    
    if(length(independents)==1){
      # Simple regresion is the covariance divided by the variance of the predictor
      model[row,independents] <- S[row,independents]/S[independents,independents]
    }
    if(length(independents)>1){
      coefs <- solve(S[independents,independents],S[row,independents])
      model[row,independents] <- coefs
    }
  }
  
  return(model)
}

#
# Functions adapted from simsem
#

is.partable <- function(object) {
  is.list(object) && all(names(object) %in% c("id", "lhs", "op", "rhs", "user", "group", "free", "ustart", "exo", "label", "eq.id", "unco"))
}

is.lavaancall <- function(object) {
  is.list(object) && ("model" %in% names(object))
}


