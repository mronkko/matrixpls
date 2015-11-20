#' Matrix-based Partial Least Squares estimation
#'
#'\pkg{matrixpls} calculates composite variable models using partial least squares (PLS)
#'algorithm and related methods. In contrast to most other PLS software which implement
#'the raw data version of the algorithm, \pkg{matrixpls} works with data covariance matrices. The 
#'algorithms are designed to be computationally efficient, modular in programming, and well documented.
#'\pkg{matrixpls} integrates with \pkg{simsem} to enable Monte Carlo simulations with as little custom
#'programming as possible.
#'
#'\pkg{matrixpls} calculates models where sets of indicator variables are combined as weighted composites. 
#'These composites are then used to estimate a statistical model describing the relationships between
#'the composites and composites and indicators. While a number of such methods exists, the partial
#'least squares (PLS) technique is perhaps the most widely used.
#'
#'The \pkg{matrixpls} package implements a collection of PLS techniques as well as the more recent 
#'GSCA and PLSc techniques and older methods based on analysis with composite variables,
#'such as regression with
#'unit weighted composites or factor scores. The package provides a unified
#'framework that enables the comparison and analysis of these algorithms. In contrast to previous R
#'packages for PLS, such as \pkg{plspm} and \pkg{semPLS} and all currently
#'available commercial PLS software, which work with
#'raw data, \pkg{matrixpls} calculates the indicator weights and model estimates from data covariance 
#'matrices. Working with covariance data allows for reanalyzing covariance matrices that are sometimes
#'published as appendices of articles, is computationally more efficient, and lends itself more easily
#'for formal analysis than implementations based on raw data. 
#'
#'@section Overview of package usage:
#'
#'\pkg{matrixpls} has modular design that is easily expanded and contains more calculation options
#'than the two other PLS packages for R. To allow validation of the algorithms by end users
#'and to help porting existing analysis files from the two other R packages to 
#'\pkg{matrixpls}, the package contains compatibility functions for both \pkg{plspm} and \pkg{semPLS}.
#'
#'The desing principles and functionality of the package is best explained by first explaining the main
#'function \code{matrixpls}. The function performs two tasks. It first calculates a set of indicator 
#'weights to form composites based on data covariance matrix and then estimates a statistical model
#'with the indicators and composites using the weights. The main function takes the following arguments:
#'
#'\preformatted{
#'matrixpls(S, model, W.model = NULL, 
#'           weightFunction = weight.pls, 
#'           parameterEstimator = params.regression,
#'           weightSignCorrection = NULL, ..., 
#'           validateInput = TRUE, standardize = TRUE)
#'}
#'
#'The first five arguments of \code{\link{matrixpls}} are most relevant for understanding how the package
#'works. \code{S}, is the data covariance or correlation matrix. \code{model} defines the model
#'which is estimated in the second stage and \code{W.model} defines how the indicators are to be
#'aggregated as composites. If \code{W.model} is left undefined, it will be constructed based on
#'\code{model} following rules that are explained elsewhere in the documentation.
#' \code{weightFunction} and 
#'\code{parameterEstimator} are functions that 
#'implement the first and second task of the function respectively. All other arguments are passed 
#'down to these two functions, which in turn can pass arguments to other functions that they call.
#'
#'@template modelSpecification
#'
#'@section Specifying estimation algorithms:
#'
#'Many of the commonly used arguments of \code{matrixpls} function are functions themselves. For 
#'example, executing a PLS analysis with Mode B outer estimation for all indicator blocks and centroid inner 
#'estimation could be specified as follows:
#'
#'\preformatted{
#'matrixpls(S, model, 
#'           outerEstimators = outer.modeB,
#'           innerEstimator = inner.centroid)
#'}
#'
#'The arguments \code{outerEstimators} and \code{innerEstimator} are not defined by the
#'\code{matrixpls} function, but are passed down to \code{weight.pls} which is used as the default
#'\code{weightFunction}. \code{outer.modeB} and \code{inner.centroid} are themselves functions provided
#'by the \pkg{matrixpls} package, which perform the actual inner and outer estimation stages of the
#'PLS algorithm. Essentially, all parts of the estimation algorithm can be provided as arguments for
#'the main function. This allows for adjusting the inner workings of the algorithm in a way that is
#'currently not possible with any other PLS software.
#'
#'It is also possible to define custom functions. For example, we could define a new Mode B outer
#'estimator that only produces positive weights by creating a custom function:
#'
#'\preformatted{
#'myModeB <- function(...)\{
#'   abs(outer.ModeB(...))
#'\}
#'
#'matrixpls(S, model,
#'           outerEstimators = myModeB,
#'           innerEstimator = inner.centroid)
#'}
"_PACKAGE"


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
#'The argument \code{W.model} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites. The original papers about
#'Partial Least Squares as well as all current PLS implementations define this as
#'\code{t(reflective) | formative}, which means that the weight patter must match the
#'model specified in \code{reflective} and \code{formative}. Matrixpls does not
#'require that \code{W.model} needs to match \code{reflective} and \code{formative}, but
#'accepts any numeric matrix. If this argument is not specified, \code{W.model} is
#'defined as \code{t(reflective) | formative}.
#'
#'@param S Covariance matrix of the data.
#
#'@param model There are two options for this argument: 1. lavaan script or lavaan parameter
#'table, or 2. a list containing three matrices
#'\code{inner}, \code{reflective}, and \code{formative} defining the free regression paths
#'in the model.
#'
#'@param W.model An optional numeric matrix representing the weight patter and starting weights
#'(i.e. the how the indicators are combined to form the composite variables). If this argument is not specified,
#'the weight patter is defined based on the relationships in the \code{reflective} and  \code{formative}
#'elements of \code{model}.

#'@param weightFunction A function that takes three or more arguments, the data covariance matrix \code{S},
#'a model specification \code{model}, and a weight pattern \code{W.model} and 
#' returns a named vector of parameter estimates. The default is \code{\link{weight.pls}}
#'
#'@param parameterEstimator A function that takes three or more arguments, the data covariance matrix \code{S},
#'model specification \code{model}, and weights \code{W} and returns a named vector of parameter estimates. The default is \code{\link{params.regression}}
#'
#'@param weightSignCorrection A function that is applied to weights to correct their signs before parameter estimation.
#'
#'@param validateInput A boolean indicating whether the arguments should be validated.
#'
#'@param standardize A boolean indicating whether \code{S} is standardized
#' 
#'@param ... All other arguments are passed through to \code{weightFunction} and \code{parameterEstimator}.
#'
#'@seealso 
#'Weight algorithms: \code{\link{weight.pls}}; \code{\link{weight.fixed}}; \code{\link{weight.optim}}
#'
#'Weight sign corrections:\code{\link{weightSign.Wold1985}}; \code{\link{weightSign.dominantIndicator}}
#'
#'Parameter estimators: \code{\link{params.regression}}
#'
#'@return A named numeric vector of class \code{matrixpls} containing parameter estimates followed by weight.
#'
#'@references Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. \emph{Journal of Statistical Software}, 48(2), 1–36. Retrieved from http://www.jstatsoft.org/v48/i02
#'
#'Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. In K. G. Jöreskog & S. Wold (Eds.),\emph{Systems under indirect observation: causality, structure, prediction} (pp. 1–54). Amsterdam: North-Holland.
#'
#'@export


matrixpls <- function(S, model, W.model = NULL, weightFunction = weight.pls,
                      parameterEstimator = params.regression, weightSignCorrection = NULL,
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
  if(is.null(W.model)) W.model <- defaultWeightModelWithModel(nativeModel)
  if(! identical(colnames(W.model), colnames(S))) stop("Column names of W.model do not match the variable names")
  if(! identical(rownames(W.model), lvNames)) stop("Row names of W.model do not match the latent variable names")
  
  ##################################################################################################
  #
  # Start of the estimation process
  #
  ##################################################################################################
  
  if(standardize){
    S <- cov2cor(S)
    # cov2cor can result in nonsymmetrix matrix because of computational inaccuracy
    S[lower.tri(S)] <- t(S)[lower.tri(S)]
  }
  
  # Calculate weight.
  
  W <- weightFunction(S, W.model = W.model,
                      model = nativeModel, parameterEstimator = params.regression,
                      ..., validateInput = validateInput, standardize = standardize)
  
  # Correct the signs of the weights
  
  if(!is.null(weightSignCorrection)) W <- weightSignCorrection(W)
  
  # Apply the parameter estimator and return the results
  
  estimates <- parameterEstimator(S, model, W, ...)
  
  ##################################################################################################
  #
  # Construct the return object
  #
  ##################################################################################################
  
  # Construct the return vector by combining estimates with weights in a named vector
  indices <- W.model!=0
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
              r2 = r2(object),
              residuals = residuals(object),
              gof = gof(object),
              cr = cr(object),
              ave = ave(object))
  
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
#'The argument \code{W.model} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites. The original papers about
#'Partial Least Squares as well as all current PLS implementations define this as
#'\code{t(reflective) | formative}, which means that the weight patter must match the
#'model specified in \code{reflective} and \code{formative}. Matrixpls does not
#'require that \code{W.model} needs to match \code{reflective} and \code{formative}, but
#'accepts any numeric matrix. If this argument is not specified, \code{W.model} is
#'defined as \code{t(reflective) | formative}.
#'
#'@inheritParams matrixpls
#'
#'@param outerEstimators A function or a list of functions used for outer estimation. If
#'the value of this parameter is a function, the same function is applied to all
#'composites. If the value of the parameter is a list, the composite \code{n} is estimated
#'with the estimator in the \code{n}th position in the list. If this argument is
#'\code{NULL} the starting weights specified in \code{W.model} will be returned. The default is \code{\link{outer.modeA}} 
#'(PLS Mode A estimation).
#'
#'@param innerEstimator A function used for inner estimation. The default is \code{\link{inner.path}}
#'(PLS path weighting scheme).

#'@param convCheck A function that takes the old and new weight matrices and
#'returns a scalar that is compared against \code{tol} to check for convergence. The default
#'function calculates the differences between each cell of old and new weights and returns the largest
#'absolute difference.
#'
#'
#'@param tol Decimal value indicating the tolerance criterion for the
#'iterations. 
#'
#'@param iter An integer indicating the maximum number of iterations.
#'
#'@param validateInput A boolean indicating whether the validity of the parameter values should be tested.
#'
#'@param variant Choose either Lohmöller's (\code{"lohmoller"}, default) or Wold's (\code{"wold"}) 
#'variant of PLS. In Wold's variant the inner and outer estimation steps are repeated for eac
#'indicator block whereas in Lohmöller's variant the weights for all composites are calculated
#'simultaneously. 
#'
#'@param ... All other arguments are passed through to \code{outerEstimators} and \code{innerEstimator}.
#'
#'@template weights-return
#'
#'@seealso 
#'Inner estimators: \code{\link{inner.path}}; \code{\link{inner.centroid}}; \code{\link{inner.factor}}; \code{\link{inner.gsca}}; \code{\link{inner.identity}};
#'
#'Outer estimators: \code{\link{outer.modeA}}; \code{\link{outer.modeB}}; \code{\link{outer.gsca}}
#'
#'Convergence checks: \code{\link{convCheck.absolute}}, \code{\link{convCheck.square}}, and \code{\link{convCheck.relative}}.
#'
#'@export

weight.pls <- function(S, model, W.model,
                       outerEstimators = NULL, 
                       innerEstimator = inner.path, ..., 
                       convCheck = convCheck.absolute,
                       variant = "lohmoller",
                       tol = 1e-05, iter = 100, validateInput = TRUE) {
  
  if(validateInput){
    
    # All parameters must have values
    assertive::assert_all_are_not_na(formals())
    
    # S must be symmetric and a valid covariance matrix
    assertive::assert_is_matrix(S)
    assertive::assert_is_symmetric_matrix(S)
    assertive::assert_is_identical_to_true(matrixcalc::is.positive.semi.definite(S))
    
    # W.model must be a real matrix and each indicators must be
    # linked to at least one composite and each composite at least to one indicator
    assertive::assert_is_matrix(W.model)
    assertive::assert_all_are_real(W.model)
    
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
    
    # outerEstimators must be a list of same length as number of rows in inner.mod or
    # a function
    if(! is.null(outerEstimators)){
      if(is.list(outerEstimators)){
        assertive::assert_is_identical_to_true(length(outerEstimators) == nrow(W.model))
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
  # Mode A for reflective composites and Mode B for formative composites
  
  if(is.null(outerEstimators)){
    hasFormativeIndicators <- any(nativeModel$formative == 1)
    hasReflectiveIndicators <- any(nativeModel$reflective == 1)
    
    if(! hasFormativeIndicators) outerEstimators = outer.modeA
    else if (! hasReflectiveIndicators) outerEstimators = outer.modeB
    else{
      # composites with at least one reflective indicator are ModeA and others are ModeB
      outerEstimators <- list()
      for(composite in 1:ncol(nativeModel$reflective)){
        if(any(nativeModel$reflective[,composite] == 1)) outerEstimators[[composite]] <- outer.modeA
        else outerEstimators[[composite]] <- outer.modeB
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
  
  if(is.list(outerEstimators)){
    uniqueOuterEstimators <- unique(outerEstimators)
    outerEstimatorIndices <- lapply(uniqueOuterEstimators, function(x){
      sapply(outerEstimators, function(y){
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
      
      if(! is.null(innerEstimator)){
        E <- innerEstimator(S, W, inner.mod, model = model, ...)
      }
      
      # Get new weights from outer estimation
      
      # Lohmöller
      
      if(is.na(k)){
        if(is.list(outerEstimators)){
          
          # Run each estimator separately
          
          for(i in 1:length(uniqueOuterEstimators)){
            W.modelForThisEstimator <- W.model
            W.modelForThisEstimator[!outerEstimatorIndices[[i]],] <- 0
            W[outerEstimatorIndices[[i]],] <- uniqueOuterEstimators[[i]](S, W_old, E, W.modelForThisEstimator,...)[outerEstimatorIndices[[i]],]
          }
        }
        else{
          W <- outerEstimators(S, W_old, E, W.model, model = model, ...)
        }	
      }
      
      # Wold
      
      else{

        if(is.list(outerEstimators)) outerEstimator <- outerEstimators[[k]]
        else outerEstimator <- outerEstimators[[k]]
        
        W.modelForThisEstimator <- W.model
        W.modelForThisEstimator[-k,] <- 0
        W[W.modelForThisEstimator != 0] <- outerEstimator(S, W_old, E, W.modelForThisEstimator,...)[W.modelForThisEstimator != 0]        
      }
      W <- scaleWeights(S, W)
    }
    
    iteration <- iteration +1 
    weightHistory[iteration+1,] <- W[weightPattern]
    
    # Check convergence. If we are not using inner estimator, converge to the first iteration
    
    if(is.null(innerEstimator) || convCheck(W,W_old) < tol){
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
#'The argument \code{W.model} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites. The original papers about
#'Partial Least Squares as well as all current PLS implementations define this as
#'\code{t(reflective) | formative}, which means that the weight patter must match the
#'model specified in \code{reflective} and \code{formative}. Matrixpls does not
#'require that \code{W.model} needs to match \code{reflective} and \code{formative}, but
#'accepts any numeric matrix. If this argument is not specified, \code{W.model} is
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



weight.optim <- function(S, model, W.model,
                         parameterEstimator = params.regression, 
                         optimCriterion = optim.maximizeInnerR2, method = "BFGS",
                         ..., 
                         validateInput = TRUE,
                         standardize = TRUE) {
  
  W <- W.model
  
  optim.res <- optim(W.model[W.model != 0], fn = function(par, ...){
    W[W.model != 0] <- par
    # Use fixed weights estimation
    matrixpls.res <- matrixpls(S, model, W, weightFunction = weight.fixed,
                               parameterEstimator = params.regression,
                               ..., validateInput = FALSE, standardize = standardize)
    
    optimCriterion(matrixpls.res)
  }, method = method, ...)
  
  W[W.model != 0] <- optim.res$par
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
#'in \code{W.model}. If \code{standardize} is \code{TRUE} the weights are
#'standardized so that the composites have unit variances.
#
#'@inheritParams matrixpls
#'
#'@param model Ignored by \code{weight.fixed}
#'
#'@param ... All other arguments are ignored.
#'
#'@template weights-return
#'
#'@export


weight.fixed <- function(S, model, W.model = NULL,
                         ..., 
                         standardize = TRUE) {
  
  W <- W.model 
  
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

#'@title Blockwise factor score weights
#'
#'@description
#'
#'Calculates weights by estimating a common factor analysis model with a single factor for each 
#'indicator block and using the resulting estimates to calculate factor score weights
#'
#'@inheritParams matrixpls
#'
#'@param model Ignored by \code{weight.factor}
#'
#'@param fm factoring method for estimating the common factor model. Possible values are
#'\code{minres}, \code{wls}, \code{gls}, \code{pa}, and \code{ml}. The parameter is passed through to
#' to \code{\link[psych]{fa}}.
#'
#'@param ... All other arguments are ignored.
#'
#'@template weights-return
#'
#'@export



weight.factor <- function(S, model, W.model = NULL, ..., fm ="minres",
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
  attr(W,"iterations") <- 0
  attr(W,"converged") <- TRUE
  attr(W,"history") <- W
  class(W) <-("matrixplsweights")
  
  return(W)
  
}

#'@title Blockwise principal component weights
#'
#'@description
#'
#'Calculates weights by calculating a principal component analysis for each 
#'indicator block and returning the weights for the first principal component.
#'
#'@inheritParams matrixpls
#'
#'@param model Ignored by \code{weight.principal}
#'
#'@param ... All other arguments are ignored.
#'
#'@template weights-return
#'
#'@export

weight.principal <- function(S, model, W.model = NULL, ..., 
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

# =========== Weight sign corrections ===========

#'@title Wold's sign ambiguity correction
#'
#'@description Adjust the signs of W so that the majority of the indicators are positively 
#'correlated with the composite as proposed by Wold (1985).
#'
#'@inheritParams inner.factor
#'
#'@return W after sign correction.
#'
#'@family Sign corrections
#'@seealso 
#'\code{\link{matrixpls}};
#'
#'
#'@references
#'Wold, H. (1985). Partial Least Squares. In S. Kotz & N. L. Johnson (Eds.), Encyclopedia  of 
#'statistical sciences (Vol. 6, pp. 581–591). New York: Wiley.
#'
#'@export

weightSign.Wold1985 <- function(W,S){
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  signSums <- rowSums(sign(IC * (W!=0)))
  sweep(W, 1,ifelse(signSums<0, -1, 1),"*")
}

#'@title Dominant indicator sign ambiguity correction
#'
#'@description Adjust the signs of W so that the first indicator of each composite has positive
#'weight.
#'
#'@inheritParams inner.factor
#'
#'@return W after sign correction.
#'
#'@family Sign corrections
#'@seealso 
#'\code{\link{matrixpls}};
#'
#'@export

weightSign.dominantIndicator <- function(W,S){
  # Calculate the covariance matrix between indicators and composites
  signs <- apply(W,1,function(x){
    i <- min(which(x!=0))
    sign(x[i])
  })
  sweep(W, 1,signs,"*")
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
#'@param parametersInner A function that takes three or more arguments, the data covariance matrix \code{S},
#'model specification \code{model}, and weights \code{W} and returns a named vector of parameter estimates. The default is \code{\link{params.regression}}
#'
#'@param parametersReflective A function that takes three or more arguments, the data covariance matrix \code{S},
#'model specification \code{model}, and weights \code{W} and returns a named vector of parameter estimates. The default is \code{\link{params.regression}}
#'
#'@param parametersFormative A function that takes three or more arguments, the data covariance matrix \code{S},
#'model specification \code{model}, and weights \code{W} and returns a named vector of parameter estimates. The default is \code{\link{params.regression}}
#'
#'@param disattenuate \code{TRUE} or \code{FALSE} (default) indicating whether C should be disattenuated before applying \code{parametersInner}.
#'
#'@param reliabilities A function that takes three or more arguments, the data covariance matrix \code{S},
#'matrix of estimated loadings \code{loadings}, and weights \code{W} and returns a named vector of 
#'estimated composite reliabilities. The default is \code{\link{reliability.weightLoadingProduct}}
#'@inheritParams matrixpls
#'
#'@return A named vector of parameter estimates.
#'
#'@family parameter estimators
#'
#'@export

params.regression <- function(S, model, W, ...,
                              parametersInner = estimator.ols,
                              parametersReflective = estimator.ols,
                              parametersFormative = estimator.ols,
                              disattenuate = FALSE,
                              reliabilities = reliability.weightLoadingProduct){
  
  nativeModel <- parseModelToNativeFormat(model)
  
  results <- c()
  
  # Calculate the composite covariance matrix
  C <- W %*% S %*% t(W)
  
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  
  reflectiveEstimates <- parametersReflective(S, nativeModel$reflective, W, ..., C = C, IC = t(IC))  
  
  if(disattenuate){
    
    Q <- reliabilities(S, reflectiveEstimates, W, ...)
    
    C <- C / sqrt(Q) %*% t(sqrt(Q))
    diag(C) <- 1
    
    # Fix the IC matrix. Start by replacing correlations with the corrected loadings
    tL <- t(reflectiveEstimates)
    IC[tL!=0] <- tL[tL!=0]
    
    # Disattenuate the remaining correlations
    IC[tL==0] <- (IC/sqrt(Q))[tL==0]
  }
  
  formativeEstimates <- parametersFormative(S, nativeModel$formative, W, ..., IC = IC)  
  innerEstimates <- parametersInner(S, nativeModel$inner, W, ..., C = C)  
  
  results <- c(estimatesMatrixToVector(innerEstimates, nativeModel$inner, "~"),
               estimatesMatrixToVector(reflectiveEstimates, nativeModel$reflective, "=~", reverse = TRUE),
               estimatesMatrixToVector(formativeEstimates, nativeModel$formative, "<~"))
  
  # Copy all non-standard attributes from the estimates objects
  
  for(object in list(innerEstimates, reflectiveEstimates, formativeEstimates)){
    for(a in setdiff(names(attributes(object)), c("dim", "dimnames", "class", "names"))){
      attr(results,a) <- attr(object,a)
      attr(W,a) <- NULL
    }
  }
  
  # Store these in the result object
  attr(results,"C") <- C
  attr(results,"IC") <- IC
  attr(results,"inner") <- innerEstimates
  attr(results,"reflective") <- reflectiveEstimates
  attr(results,"formative") <- formativeEstimates
  
  if(disattenuate){
    attr(results,"Q") <- Q
  }
  
  return(results)
}

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
#'@param S Covariance matrix of the data.
#'
#'@param model A model matrix with dependent variables on rows, independent variables on colums, and
#'non-zero elements indicating relationships. Can be either \code{inner}, \code{reflective},
#'or \code{formative} matrix.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param C Correlation matrix of the composites.
#'
#'@param IC Correlation matrix of the composites and indicators.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A matrix with estimated parameters.
#'
#'@family regression estimators
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
#'@param S Covariance matrix of the data.
#'
#'@param model A model matrix with dependent variables on rows, independent variables on colums, and
#'non-zero elements indicating relationships. Onyl the \code{inner} matrix can be used.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param C Correlation matrix of the composites.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A matrix with estimated parameters.
#'
#'@family regression estimators
#'
#'@export
#'

estimator.tsls <- function(S, model, W, ..., C){
  
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

# =========== Reliability estimators ===========

#'@title Reliabilities as products of weights and loadings
#'
#'@description
#'Calculates reliabilities as a matrix product of loadings and weidhts
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
#'@family reliability estimators
#'
#'@export

reliability.weightLoadingProduct <- function(S, loadings, W, ...){
  
  Q <- diag(W %*% loadings)^2
  
  # Any composite with no reflective indicators is fixed to be perfectly reliable
  Q[apply(loadings==0,2,all)] <- 1
  
  Q
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
  
  E <- estimator.ols(S, inner.mod,W, C = C)
  
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

inner.gsca <- function(S, W, inner.mod, ...){
  E <- estimator.ols(S, inner.mod,W)
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
#'@param W.model A matrix specifying the weight relationships and their starting values.
#'
#'@inheritParams inner.centroid
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.model}.
#'
#'@family outer estimators  
#'
#'@export


outer.modeA <- function(S, W, E, W.model, ...){
  
  # Calculate the covariance matrix between indicators and composites
  W_new <- E %*% W %*% S
  
  # Set the non-existing weight relations to zero
  W_new[W.model == 0] <- 0
  
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
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.model}.
#'
#'@family outer estimators
#'
#'@export


outer.modeB <- function(S, W, E, W.model, ...){
  
  # Calculate the covariance matrix between indicators and composites
  IC <- E %*% W %*% S
  
  # Set up a weight pattern
  W_new <- ifelse(W.model==0,0,1)
  
  # Do the outer model regressions
  
  for(row in which(rowSums(W_new)>0, useNames = FALSE)){
    indicatorIndices <- W_new[row,]==1
    W_new[row,indicatorIndices] <- solve(S[indicatorIndices,indicatorIndices],IC[row,indicatorIndices])
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
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.model}.
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

outer.gsca <- function(S, W, E, W.model, model, ...){
  
  nativeModel <- parseModelToNativeFormat(model)
  
  inner <- nativeModel$inner
  
  # Calculate the covariance matrix between indicators and composites and between composites
  IC <- W %*% S
  C <- W %*% S %*% t(W)
  
  # Estimate the reflective parts of the model
  
  reflective <- nativeModel$reflective
  
  for(row in which(rowSums(reflective!=0)>0)){
    independents <- which(reflective[row,] != 0)
    reflective[row,independents] <- solve(C[independents,independents],IC[independents, row])
  }  
  
  # Number of composites and indicators, and their sum
  P <- nrow(W.model)
  J <- ncol(W.model)
  JP <- J + P
  
  # E, and reflective form the A matrix of GSCA. 
  # Indicators first, then composites
  
  A <- rbind(reflective, E)
  V <- rbind(diag(J), W)
  
  # The following code is based on the ASGSCA package (licensed
  # under GPL-3). All matrices are transposed from the original
  # ASGSCA code
  
  # Step 2: Update W
  
  tr_w <- 0
  for(p in 1:P){
    t <- J + p
    windex_p <- which(W.model[p, ] != 0)
    m <- matrix(0, 1, JP)
    m[t] <- 1
    a <- A[, p]
    beta <- m - a
    H1 <- diag(P)
    H2 <- diag(JP)
    H1[p,p] <- 0
    H2[t,t] <- 0
    
    Delta <- A%*%H1%*%W - H2%*%V 
    Sp <- S[windex_p , windex_p]
    if (length(windex_p)!=0){        
      
      theta <- MASS::ginv(as.numeric(beta%*%t(beta))*S[windex_p,windex_p]) %*%
        t(beta %*% Delta %*% S[,windex_p])
      
      # Update the weights based on the estimated parameters and standardize
      W[p,windex_p] <- theta
      W <- scaleWeights(S, W)
      
    }
    
    # Proceed to next composite
  }
  
  
  return(W)
}


# =========== Optimization criterion functions ===========

#'@title Optimization criterion for maximizing inner R2
#'
#'@description Calculates the sum of R2s of \code{inner} matrix. The
#'optimization criterion is negative of this sum.
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@return Negative of sum of R2 values of the \code{inner} matrix.
#'
#'@family Weight optimization criteria
#'
#'@export
#'
#'

optim.maximizeInnerR2 <- function(matrixpls.res){
  -sum(R2(matrixpls.res))
}

# #'@title Optimization criterion for maximal prediction 
# #'
# #'@details Calculates the predicted variances of reflective indicators. The
# #'prediction criterion is negative of the sum of predicted variances.
# #'
# #'@param matrixpls.res An object of class \code{matrixpls} from which the
# #'criterion function is calculated
# #'
# #'@return Mean squared prediction error.
# #'
# #'@family Weight optimization criteria
# #'
# #'@export
# #'
# 
# optim.maximizePrediction <- function(matrixpls.res){
#   
#   # TODO: convert to using the predict function
#   
#   C <- attr(matrixpls.res,"C")
#   nativeModel <- attr(matrixpls.res,"model")
#   exog <- rowSums(nativeModel$inner)==0
#   W <- attr(matrixpls.res,"W")
#   inner <- attr(matrixpls.res,"inner")
#   
#   reflective <- nativeModel$reflective
#   reflective[which(reflective==1)] <- matrixpls.res[grep("=~",names(matrixpls.res))]
#   
#   # Predict endog LVs using reduced from equations
#   Beta <- inner[! exog, ! exog]
#   Gamma <- inner[! exog, exog]
#   
#   invBeta <- solve(diag(nrow(Beta)) - Beta)
#   endogC <- invBeta %*% Gamma %*% C[exog,exog] %*% t(Gamma) %*% t(invBeta)
#   C[!exog, !exog] <- endogC
#   
#   -sum(diag(reflective %*% C %*% t(reflective)))
# }

#'@title Optimization criterion for maximizing reflective R2
#'
#'@description Calculates the sum of R2s of \code{reflective} matrix. The
#'optimization criterion is negative of this sum.
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@return Negative of sum of R2 values of the \code{reflective} matrix.
#'
#'@family Weight optimization criteria
#'
#'@export
#'
#'

optim.maximizeIndicatorR2 <- function(matrixpls.res){
  lambda <- loadings(matrixpls.res)
  IC <- attr(matrixpls.res,"IC")
  -sum(diag(lambda %*% IC))
}

#'@title Optimization criterion for maximizing reflective and innner R2
#'
#'@description Calculates the sum of R2s of \code{inner} matrix and 
#'\code{reflective} matrix. The
#'optimization criterion is negative of this sum.
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@return Negative of sum of R2 values of the \code{inner} and
#'\code{reflective} matrices.
#'
#'@family Weight optimization criteria
#'
#'@export
#'

optim.maximizeFullR2 <- function(matrixpls.res){
  optim.maximizeIndicatorR2(matrixpls.res) + optim.maximizeInnerR2(matrixpls.res)
}

#'@title GSCA optimization criterion
#'
#'@description Optimization criterion for minimizing the sum of all residual
#'variances in the model. 
#'
##'The matrixpls implementation of the GSCA criterion extends the criterion
##'presented by Huang and Takane by including also the minimization of the
##'residual variances of the formative part of the model. The formative
##'regressions in a model are typically specified to be identical to the 
##'weight pattern \code{W.model} resulting zero residual variances by definition.
##'However, it is possible to specify a formative model that does not completely
##'overlap with the weight pattern leading to non-zero residuals that can be
##'optimizes.
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

optim.gsca <- function(matrixpls.res){
  
  C <- attr(matrixpls.res,"C")
  IC <- attr(matrixpls.res,"IC")
  nativeModel <- attr(matrixpls.res,"model")
  
  reflective <- nativeModel$reflective
  reflective[which(reflective==1)] <- matrixpls.res[grep("=~",names(matrixpls.res))]
  
  #  formative <- nativeModel$formative
  #  formative[which(formative==1)] <- matrixpls.res[grep("<~",names(matrixpls.res))]
  
  #  f <- apply(nativeModel$formative != 0,1,any)
  r <- apply(nativeModel$reflective != 0,1,any)
  endo <- apply(nativeModel$inner != 0,1,any)
  
  inner_resid <- (1 - R2(matrixpls.res)[endo])
  #  form_resid <- (1 - rowSums(IC[f,] * formative[f,]))
  refl_resid <- (1 - rowSums(t(IC[,r]) * reflective[r,]))
  
  #  sum(inner_resid, form_resid, refl_resid)
  sum(inner_resid, refl_resid)
}


# =========== Utility functions ===========

estimatesMatrixToVector <- function(est, m, sep, reverse = FALSE){
  
  ind <- which(m != 0)
  ret <- est[ind]
  
  if(reverse){
    names(ret) <- paste(colnames(m)[col(m)[ind]],
                        rownames(m)[row(m)[ind]], sep=sep)
  }
  else{
    names(ret) <- paste(rownames(m)[row(m)[ind]],
                        colnames(m)[col(m)[ind]], sep=sep)
  }
  
  ret
}


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
  observedVariableNames <- setdiff(unique(c(factorLoadings$rhs,formativeLoadings$rhs,
                                            regressions$lhs,regressions$rhs)), latentVariableNames)
  
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
    
    return(lavaanParTableToNativeFormat(lavaan::lavaanify(model)))
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
  W.model <- matrix(0,nrow(nativeModel$formative), ncol(nativeModel$formative))
  
  colnames(W.model) <- colnames(nativeModel$formative)
  rownames(W.model) <- rownames(nativeModel$formative)
  
  W.model[nativeModel$formative!=0] <- 1 
  W.model[t(nativeModel$reflective)!=0] <- 1 
  
  return(W.model)
  
}

is.matrixpls.model <- function(model) {
  return (is.list(model) &&
            setequal(names(model),c("inner","reflective","formative")) &&
            all(sapply(model,is.matrix)))
  
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


