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
#'packages for PLS, such as plspm and semPLS and all currently
#'available commercial PLS software, which work with
#'raw data, \pkg{matrixpls} calculates the indicator weights and model estimates from data covariance 
#'matrices. Working with covariance data allows for reanalyzing covariance matrices that are sometimes
#'published as appendices of articles, is computationally more efficient, and lends itself more easily
#'for formal analysis than implementations based on raw data. 
#'
#'\pkg{matrixpls} has modular design that is easily expanded and contains more calculation options
#'than the two other PLS packages for R. To allow validation of the algorithms by end users
#'and to help porting existing analysis files from the two other R packages to 
#'\pkg{matrixpls}, the package contains compatibility functions for both plspm and semPLS.
#'
#'The design principles and functionality of the package is best explained by first explaining the main
#'function \code{matrixpls}. The function performs two tasks. It first calculates a set of indicator 
#'weights to form composites based on data covariance matrix and then estimates a statistical model
#'with the indicators and composites using the weights. The main function takes the following arguments:
#'
#'\preformatted{
#'matrixpls(S, model, W.model = NULL, 
#'           weightFun = weightFun.pls, 
#'           parameterEstim = parameterEstim.separate,
#'           weightSign = NULL, ..., 
#'           validateInput = TRUE, standardize = TRUE)
#'}
#'
#'The first five arguments of \code{\link{matrixpls}} are most relevant for understanding how the package
#'works. \code{S}, is the data covariance or correlation matrix. \code{model} defines the model
#'which is estimated in the second stage and \code{W.model} defines how the indicators are to be
#'aggregated as composites. If \code{W.model} is left undefined, it will be constructed based on
#'\code{model} following rules that are explained elsewhere in the documentation.
#' \code{weightFun} and 
#'\code{parameterEstim} are functions that 
#'implement the first and second task of the function respectively. All other arguments are passed 
#'down to these two functions, which in turn can pass arguments to other functions that they call.
#'
#'@template modelSpecification
#'
#'@details
#'Many of the commonly used arguments of \code{matrixpls} function are functions themselves. For 
#'example, executing a PLS analysis with Mode B outer estimation for all indicator blocks and centroid inner 
#'estimation could be specified as follows:
#'
#'\preformatted{
#'matrixpls(S, model, 
#'           outerEstim = outerEstim.modeB,
#'           innerEstim = innerEstim.centroid)
#'}
#'
#'The arguments \code{outerEstim} and \code{innerEstim} are not defined by the
#'\code{matrixpls} function, but are passed down to \code{weightFun.pls} which is used as the default
#'\code{weightFun}. \code{outerEstim.modeB} and \code{innerEstim.centroid} are themselves functions provided
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
#'   abs(outerEstim.ModeB(...))
#'\}
#'
#'matrixpls(S, model,
#'           outerEstim = myModeB,
#'           innerEstim = innerEstim.centroid)
#'}
#'
#'@references
#'
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial least squares.} Heidelberg: Physica-Verlag.
#'
#'Rönkkö, M., McIntosh, C. N., & Antonakis, J. (2015). On the adoption of partial
#' least squares in psychological research: Caveat emptor. 
#' \emph{Personality and Individual Differences}, (87), 76–84. 
#' \doi{10.1016/j.paid.2015.07.019}
#'
#'Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. 
#'In K. G. Jöreskog & S. Wold (Eds.),\emph{Systems under indirect observation:
#'causality, structure, prediction} (pp. 1–54). Amsterdam: North-Holland.
#'
#'
#'@aliases matrixpls-package
#'
"_PACKAGE"

#'Commonly used arguments
#'
#'The following describes commonly used arguments in the \pkg{matrixpls} package.
#'
#'@param S Covariance matrix of the data.
#
#'@param model There are two options for this argument: 1. lavaan script or lavaan parameter
#'table, or 2. a list containing three matrices
#'\code{inner}, \code{reflective}, and \code{formative} defining the free regression paths
#'in the model.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param C Correlation matrix of the composites.
#'
#'@param IC Correlation matrix of the composites and indicators.
#'
#'@param E Inner weight matrix. A square matrix of inner estimates between the composites.
#'
#'@param W.model A matrix specifying the weight relationships and their starting values.
#'


#'@name matrixpls-common
NULL

#'All estimation function types
#'
#'The following describes all estimation function types used in \pkg{matrixpls} package.
#'
#'@param weightFun A function for calculating indicator weights using the data covariance matrix
#' \code{S}, a model specification \code{model}, and a weight pattern \code{W.model}. Returns
#' a weight matrix \code{W}. The default is \code{\link{weightFun.pls}}
#'
#'@param parameterEstim A function for estimating the model parameters using
#'the data covariance matrix \code{S}, model specification \code{model}, 
#'and weight matrix \code{W}. Returns a named vector of parameter estimates.
#'The default is \code{\link{parameterEstim.separate}}
#'
#'@param estimator A function for estimating the parameters of one model matrix using
#'the data covariance matrix \code{S}, a model matrix \code{modelMatrix}, and a weight matrix
#'\code{W}. Disattenuated composite correlation matrix \code{C} and indicator composite 
#'covariance matrix \code{IC} are optional. Returns matrix of parameter estimates.
#'The default is \code{\link{estimator.ols}}
#'
#'@param weightSign A function for resolving weight sign ambiguity based on the data covariance matrix
#' \code{S} and a weight matrix \code{W}.  Returns
#' a weight matrix \code{W}. See \code{\link{weightSign}}
#' for details. 
#' 
#'@param outerEstim A function for calculating outer weights using the data covariance matrix
#' \code{S}, a weight matrix \code{W}, an inner weight matrix \code{E},
#'and a weight pattern \code{W.model}. Returns
#'a weight matrix \code{W}. See \code{\link{outerEstim}}.
#'
#'@param innerEstim A function for calculating inner weights using  the data covariance matrix
#'\code{S}, a weight matrix \code{W}, and an inner model matrix \code{inner.mod}. Returns
#'an inner weight matrix \code{E}. The default is \code{\link{innerEstim.path}}.

#'@param convCheck A function that takes the old \code{Wold} and new weight \code{Wold} matrices and
#'returns a scalar that is compared against \code{tol} to check for convergence. The default
#'is \code{\link{convCheck.absolute}}.
#'
#'@param optimCrit A function for calculating value for an optimization criterion based on a
#'\code{matrixpls} result object. Returns a scalar. The default is \code{\link{optimCrit.maximizeInnerR2}}. 
#'
#'@param reliabilities A function for calculating reliability estimates based on the 
#'data covariance matrix \code{S}, factor loading matrix \code{loadings}, and a weight matrix \code{W}.
#'Returns a vector of reliability estimates. The default is
#' \code{\link{reliabilityEstim.weightLoadingProduct}}
#' 


#'@name matrixpls-functions
NULL

# =========== Main functions ===========

#'Partial Least Squares and other composite variable models.
#'
#'Estimates a weight matrix using Partial Least Squares or a related algorithm and then
#'uses the weights to estimate the parameters of a statistical model.
#'
#'\code{matrixpls} is the main function of the matrixpls package. This function
#'parses a model object and then uses the results to call \code{weightFun} to
#'to calculate indicator weight. After this the \code{parameterEstim} function is
#'applied to the indicator weights, the data covariance matrix,
#'and the model object and the resulting parameter estimates are returned.
#'
#'@template modelSpecification
#'
#'@template weightSpecification
#'
#'@param W.model An optional numeric matrix representing the weight pattern and starting weights
#'(i.e. the how the indicators are combined to form the composite variables). If this argument is not specified,
#'the weight patter is defined based on the relationships in the \code{reflective} and  \code{formative}
#'elements of \code{model}.
#'
#'@inheritParams matrixpls-common
#'@inheritParams matrixpls-functions
#'
#'@param validateInput If \code{TRUE}, the arguments are validated.
#'
#'@param standardize If \code{TRUE}, \code{S} is converted to a correlation matrix before analysis.
#' 
#'@param ... All other arguments are passed through to \code{weightFun} and \code{parameterEstim}.
#'
#'@seealso 
#'Weight algorithms: \code{\link{weightFun.pls}}; \code{\link{weightFun.fixed}}; \code{\link{weightFun.optim}}; \code{\link{weightFun.principal}}; \code{\link{weightFun.factor}}
#'
#'Weight sign corrections:\code{\link{weightSign.Wold1985}}; \code{\link{weightSign.dominantIndicator}}
#'
#'Parameter estimators: \code{\link{parameterEstim.separate}}
#'
#'@return A named numeric vector of class \code{matrixpls} containing the parameter estimates followed by weights.
#'
#'@templateVar attributes matrixpls,W model call,S E iterations converged history C IC inner reflective formative Q c
#'@template attributes
#'
#'@references 
#'
#'Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. In K. G. Jöreskog & S. Wold (Eds.),\emph{Systems under indirect observation: causality, structure, prediction} (pp. 1–54). Amsterdam: North-Holland.
#'
#'@export
#'


matrixpls <- function(S, model, W.model = NULL, weightFun = weightFun.pls,
                      parameterEstim = parameterEstim.separate, weightSign = NULL,
                      ..., validateInput = TRUE, standardize = TRUE) {
  
  
  # TODO: Validate input.
  
  nativeModel <- parseModelToNativeFormat(model)
  
  lvNames <- colnames(nativeModel$inner)
  if(any(duplicated(lvNames)))
    stop(paste("Each composite must have a unique name. The names contain duplicates:",
               lvNames))
  
  if(! identical(lvNames, colnames(nativeModel$reflective)) ||
     ! identical(lvNames,rownames(nativeModel$formative))){
    print(nativeModel)
    stop("Names of composites are inconsistent between inner, reflective, and formative models.")
  }
  
  if(! identical(rownames(nativeModel$reflective), colnames(nativeModel$formative))){
    print(nativeModel)
    stop("Names of observed variables are inconsistent between reflective and formative models.")
  }
  
  if(! identical(colnames(S), rownames(S))) stop("S must have identical row and column names.")
  if(! is.symmetric.matrix(S)) stop("S must be symmetric to be a valid covariance matrix.")
  if(! is.positive.semi.definite(S)) stop("S must be positive semi-definite to be a valid covariance matrix.")
  
  if(! identical(colnames(S), colnames(nativeModel$formative))){
    # Do all variables of the model exists in S
    d <- setdiff(colnames(nativeModel$formative), colnames(S))
    if(length(d)>0) stop(paste("Variable(s)",paste(d,collapse=", "),"are not included in S."))
    
    # Choose the subset that we uses. This also guarantees the correct order.
    S <- S[colnames(nativeModel$formative), colnames(nativeModel$formative)]
  }
  
  # If the weight pattern is not defined, calculate it based on the model.
  if(is.null(W.model)) W.model <- defaultWeightModelWithModel(nativeModel)
  if(! identical(colnames(W.model), colnames(S))) stop("Column names of W.model do not match the variable names")
  if(! identical(rownames(W.model), lvNames)) stop("Row names of W.model do not match the latent variable names")
  
  
  
  
  ##################################################################################################
  #
  # Start of the estimation process
  #
  ##################################################################################################
  
  if(standardize){
    S <- stats::cov2cor(S)
    # cov2cor can result in nonsymmetrix matrix because of computational inaccuracy
    S[lower.tri(S)] <- t(S)[lower.tri(S)]
  }
  
  # Calculate weight.
  
  W <- weightFun(S, W.model = W.model,
                 model = nativeModel, parameterEstim = parameterEstim.separate,
                 ..., validateInput = validateInput, standardize = standardize)
  
  # Correct the signs of the weights
  
  if(!is.null(weightSign)) W <- weightSign(W, S)
  
  # Apply the parameter estimator and return the results
  
  estimates <- parameterEstim(S, model, W, ...)
  
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

#'@export

print.matrixpls <- function(x, ...){
  
  cat("\n matrixpls parameter estimates\n")
  
  indices <- ! grepl("=+", names(x), fixed=TRUE)
  toPrint <- x[indices]
  
  estimates <- as.matrix(toPrint)
  colnames(estimates)[1] <- "Est."
  
  
  
  se <- attr(x,"se")
  boot.out <- attr(x,"boot.out")
  
  if(! is.null(boot.out)){
    estimates <- cbind(estimates, apply(boot.out$t[,indices],2,stats::sd))
    colnames(estimates)[2] <- "SE"
  }
  else if(! is.null(se)){
    estimates <- cbind(estimates,se)
    colnames(estimates)[2] <- "SE"
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

#'@export

summary.matrixpls <- function(object, ...){
  
  ret <- list(estimates = object,
              effects = effects(object),
              r2 = r2(object),
              residuals = stats::residuals(object),
              gof = gof(object),
              cr = cr(object),
              ave = ave(object),
              htmt = htmt(object),
              cei = cei(object))
  
  class(ret) <-("matrixplssummary")
  
  return(ret)
}

#'@export

print.matrixplssummary <- function(x, ...){
  for(element in x){
    print(element, ...)
  }
}




