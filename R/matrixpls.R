# =========== Main functions ===========

#'@title Partial Least Squares estimation
#'
#'@description
#'Estimates a weight matrix using Partial Least Squares or a related algorithm and then
#'uses the weights to estimate the parameters of a statistical model.
#'
#'@details
#'\code{matrixpls} is the main function of the matrixpls package. This function
#'parses a model object and then uses the results to call \code{\link{matrixpls.weights}} to
#'to calculate indicator weights. After this the \code{parameterEstimator} function is
#'applied to the indicator weights, the data covariance matrix,
#'and the model object and the resulting parameter estimates are returned.
#'
#'The model object specifies the model that is estimated with \code{parameterEstimator}
#'and is used to construct \code{inner.mod} argument for \code{\link{matrixpls.weights}}.
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
#'\code{t(reflective) | formative}, which means that the weight model must match the
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
#'@param W.mod An optional numeric matrix representing the weight model and starting weights
#'(i.e. the how the indicators are combined to form the composite variables). If this argument is not specified,
#'the weight model is defined based on the relationships in the \code{reflective} and  \code{formative}
#'elements of \code{model}.
#'
#'@param parameterEstimator A function that takes three arguments, the data covariance matrix \code{S},
#'weights \code{W}, and model specification \code{model} and returns a named vector of parameter estimates.
#'
#'@param ... All other arguments are passed through to \code{matrixpls.weights}.
#'
#'@inheritParams matrixpls.weights
#'
#'@seealso \code{\link{matrixpls.weights}}
#'
#'@return A named numeric vector of class \code{matrixpls} containing parameter estimates followed by weights.
#'
#'@references Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. \emph{Journal of Statistical Software}, 48(2), 1–36. Retrieved from http://www.jstatsoft.org/v48/i02
#'
#'Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. In K. G. Jöreskog & S. Wold (Eds.),\emph{Systems under indirect observation: causality, structure, prediction} (pp. 1–54). Amsterdam: North-Holland.
#'
#'@export


matrixpls <- function(S, model, W.mod = NULL, parameterEstimator = params.regression, ..., outerEstimators = NULL, validateInput = TRUE) {
	
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
	
	# If the weight model is not defined, calculate it based on the model.
	if(is.null(W.mod)) W.mod <- defaultWeightModelWithModel(nativeModel)
	
	if(! identical(colnames(W.mod), colnames(S))) stop("Column names of W.mod do not match the variable names")
	if(! identical(rownames(W.mod), lvNames)) stop("Row names of W.mod do not match the latent variable names")
	
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
	
	# Calculate weights
	W <- matrixpls.weights(S, nativeModel$inner, W.mod, outerEstimators, ... , validateInput = validateInput)
	
	# Apply the parameter estimator and return the results
	estimates <- parameterEstimator(S, W, model)
	
	# Construct the return vector by combining estimates with weights in a named vector
	indices <- W.mod!=0
	WVect <- W[indices]
	
	names(WVect) <- paste(rownames(nativeModel$formative)[row(W)[indices]],"=+",colnames(nativeModel$formative)[col(W)[indices]], sep="")
	
	ret <- c(estimates,WVect)
	
	# Copy all non-standard attributes from the weight and estimation algorithms
	
	allAttributes <- c(attributes(W),attributes(estimates))
	
	for(a in setdiff(names(allAttributes), c("dim", "dimnames", "class", "names"))){
		attr(ret,a) <- allAttributes[[a]]
	}
	
	attr(ret,"W") <- W
	attr(ret,"model") <- nativeModel
	
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
	
	print(attr(x,"W"), ...)
	
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


#'@title Partial Least Squares weight algorithm
#'
#'@description
#'Estimates a weight matrix using Partial Least Squares or a related algorithm.
#
#'@details
#'\code{matrixpls.weights} performs the PLS weighting algorithm by calling the 
#'inner and outer estimators iteratively until either the convergence criterion or
#'maximum number of iterations is reached and provides the results in a matrix.
#'
#'The argument \code{inner.mod} is a binary (\code{l x l}) matrix, where
#'\code{l} is the number of composites. This matrix indicates how the composites
#'should be linked in the inner estimation. The matrix \code{inner.mod} will
#'contain a 1 when a variable on the column \code{j}
#'has a link toward the variable on row \code{i} and 0 otherwise. 
#'
#'The argument \code{W.mod} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites, where \code{k} is the 
#'number of observed variables. \code{W.mod} will contain a non-zero value
#'when a variable on the
#'column \code{j} is a part of the composites variable on row \code{i}, zero otherwise. 
#'The values are used as starting values for the weighting algorithm.
#'
#'@param S Covariance matrix of the data.
#
#'@param inner.mod A square matrix of ones and zeros representing the
#'inner model (i.e. the path relationships between the composite variables).
#
#'@param W.mod numeric matrix representing the weight model and starting weights
#'(i.e. the how the indicators are combined to form the composite variables).
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
#'@return An object of class \code{"matrixplsweights"}, which is a matrix containing the weights with the following attributes: 
#'@return \item{iterations}{Number of iterations performed}.
#'@return \item{converged}{A boolean indicating if the algorithm converged}.
#'@return \item{history}{A data.frame containing the weights for each iteration}.
#'@export

matrixpls.weights <- function(S, inner.mod, W.mod,
															outerEstimators = outer.modeA, 
															innerEstimator = inner.path,
															convCheck = function(W,W_new){max(abs(W_new - W))},
															tol = 1e-05, iter = 100, validateInput = TRUE) {
	
	if(validateInput){
		
		library(assertive)
		library(matrixcalc)
		
		# All parameters must have values
		assert_all_are_not_na(formals())
		
		# S must be symmetric and a valid covariance matrix
		assert_is_matrix(S)
		assert_is_symmetric_matrix(S)
		assert_is_identical_to_true(is.positive.semi.definite(S))
		
		# inner.mod must be a square matrix consisting of ones and zeros and zero diagonal
		# and all variables must be linked to at least one other variable
		assert_is_matrix(inner.mod)
		assert_is_identical_to_true(nrow(inner.mod)==ncol(inner.mod))
		assert_all_are_true(inner.mod==0 | inner.mod == 1)
		assert_all_are_true(diag(inner.mod)==0)
		
		# W.mod must be a real matrix and each indicators must be
		# linked to at least one composite and each composite at least to one indicator
		assert_is_matrix(W.mod)
		assert_all_are_real(W.mod)
		
		if(! all(apply(W.mod!=0,1,any))){
			print(W.mod)	
			stop("All constructs must have at least one indicator")	
		}
		if(! all(apply(W.mod!=0,2,any))){
			print(W.mod)	
			stop("All indicators must be linked to at least one construct")	
		}
		
		# the number of rows in W.mod must match the dimensions of other matrices
		if(nrow(inner.mod)!=nrow(W.mod)){
			print(list(inner.mod=inner.mod,W.mod = W.mod))
			stop("Inner model row count does not match weight model row count")
		}
		
		if(ncol(S)!=ncol(W.mod)){
			print(list(S=S,W.mod = W.mod))
			stop("Data matrix column count does not match weight model column count")
		}
		
		# outerEstimators must be a list of same length as number of rows in inner.mod or
		# a function
		if(! is.null(outerEstimators)){
			if(is.list(outerEstimators)){
				assert_is_identical_to_true(length(outerEstimators) == nrow(W.mod))
				for(outerEstimator in outerEstimators){
					assert_is_function(outerEstimator)
				}
			}
		}
		else{
			assert_is_function(outerEstimators)
		}
		
		if(! is.null(innerEstimator)){
			assert_is_function(innerEstimator)
		}
		# tol must be positive
		assert_all_are_positive(tol)
		
		#iter must not be negative
		assert_all_are_positive(iter)
	}
	
	# Done checking arguments
	
	# The initial weight matrix
	
	weightPattern <- W.mod!=0
	W <- scaleWeights(S, W.mod)
	iteration <- 0
	
	weightHistory <- matrix(NA,iter+1,sum(weightPattern))
	weightHistory[1,] <- W[weightPattern]
	rownames(weightHistory) <- c("start",1:100)
	
	
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
			
			# Get new inner weights from inner estimation
			
			if(! is.null(innerEstimator)){
				E <- innerEstimator(S, W, inner.mod)
			}
			
			# Get new weights from outer estimation
			
			W_old <- W
			if(is.list(outerEstimators)){
				
				# Run each estimator separately
				
				for(i in 1:length(uniqueOuterEstimators)){
					W.modForThisEstimator <- W.mod
					W.modForThisEstimator[!outerEstimatorIndices[[i]],] <- 0
					W[outerEstimatorIndices[[i]],] <- uniqueOuterEstimators[[i]](S, W_old, E, W.modForThisEstimator)[outerEstimatorIndices[[i]]]
				}
			}
			else{
				W <- outerEstimators(S, W_old, E, W.mod)
			}	
			
			W <- scaleWeights(S, W)
			
			iteration <- iteration +1 
			weightHistory[iteration+1,] <- W[weightPattern]
			
			# Check convergence. If we are not using inner estimator, converge to the first iteration
			
			if(is.null(innerEstimator) || convCheck(W,W_old) < tol){
				converged <- TRUE
				break;
			}
			else if(iteration == iter){
				converged <- FALSE;
				break;
			}
			
		}
	}
	
	# Mark the estimation as converged if not running iterations
	
	else converged <- TRUE	
	
	attr(W,"S") <- S
	attr(W,"iterations") <- iteration
	attr(W,"converged") <- converged
	attr(W,"history") <- weightHistory[1:iteration+1,]
	class(W) <-("matrixplsweights")
	rownames(W) <- rownames(inner.mod)
	
	return(W)
}

#'@S3method print matrixplsweights

print.matrixplsweights <- function(x, ...){
	cat("\n matrixpls weights\n")
	print.table(x, ...)
	cat("\nWeight algorithm",ifelse(attr(x,"converged"),"converged","did not converge"),"in",attr(x,"iterations"),"iterations.\n")
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

params.regression <- function(S, W, model){
	
	return(params.internal_generic(S,W,model,
																 regressionsWithCovarianceMatrixAndModelPattern))
}

params.internal_generic <- function(S, W, model, pathEstimator){
	
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
	browser()
	return(results)
}


params.internal_generic <- function(S, W, model, pathEstimator){
	
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
				results <- c(results, IC[independents,row])
			}
			else{
				coefs <- solve(C[independents,independents],IC[independents,row])
				results <- c(results,coefs)
			}
			names(results)[length(results)] <- paste(colnames(nativeModel$reflective)[independents], "=~",
																							 rownames(nativeModel$reflective)[row], sep = "")
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
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param inner.mod A square matrix specifying the relationships of the composites in the model.
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@family inner estimators
#'
#'@references
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.

#'@export

inner.centroid <- function(S, W, inner.mod){
	
	# Centroid is just the sign of factor weighting
	
	E <- sign(inner.factor(S, W, inner.mod))
	
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
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param inner.mod A square matrix specifying the relationships of the composites in the model.
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@family inner estimators
#'
#'@references
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.

#'@export

inner.path <- function(S, W, inner.mod){
	
	# Calculate the composite covariance matrix	
	C <- W %*% S %*% t(W)
	
	E <- regressionsWithCovarianceMatrixAndModelPattern(C, inner.mod)
	
	# Use  correlations for inverse relationships for non-reciprocal paths
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
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@family inner estimators
#'
#'@references
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.

#'@export

inner.factor <- function(S, W, inner.mod){
	
	# Calculate the composite covariance matrix
	
	C <- W %*% S %*% t(W)
	
	# Calculate unscaled inner weights using the centroid weighting scheme
	
	E <- C * (inner.mod | t(inner.mod))
	
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
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param inner.mod A square matrix specifying the relationships of the composites in the model.
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



inner.identity <- function(S, W, inner.mod){
	return(diag(nrow(inner.mod)))
}

#'@title GSCA inner estimation
#'
#'@description
#'
#'This implements the first step of the GSCA estimation describe by Hwang & Takane (2004). GSCA inner estimation should
#'be used only with GSCA outer estimation. 
#
#'@details
#'
#'The first step of GSCA estimation method, as describe by Hwang & Takane (2004), involves estimating all model regressions,
#'including also the relationships from composites toward indicators in the first step. 
#'
#'The implementation of GSCA in matrixpls differs from the Hwang & Takane (2004) version in that during the
#'first step, only regressions between composites are estimated. The reason for this is that the
#'relationhips from the composites to indicators are simple regressions that are simply the covariances between
#'the indicators and compositess. Since these covariances need to be calculated in the second step, it is more
#'efficient to not calculate them during the first step.
#'
#'This algorithm is therefore identical to the PLS path weighting scheme with the exception that correlations
#'are not used for inverse relationships and there is no falling back to identity scheme for composites
#'that are not connected to other composites.
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param inner.mod A square matrix specifying the relationships of the composites in the model
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'@references
#'Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. Psychometrika, 69(1), 81–99. doi:10.1007/BF02295841
#'
#'@example example/gsca-example.R
#'
#'@family inner estimators  
#'      
#'@export

inner.GSCA <- function(S, W, inner.mod){
	
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
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@family outer estimators  
#'
#'@export


outer.modeA <- function(S, W, E, W.mod){
	
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
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param E Inner weight matrix. A square matrix of inner estimates between the composites.
#'
#'@param W.mod A matrix specifying the weight relationships and their starting values.
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@family outer estimators
#'
#'@export


outer.modeB <- function(S, W, E, W.mod){
	
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
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param E Inner weight matrix. A square matrix of inner estimates between the composites.
#'
#'@param W.mod A matrix specifying the weight relationships and their starting values.
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@family outer estimators
#'
#'@export

outer.fixedWeights <- function(S, W, E, W.mod){
	return(W.mod)
}

#'@title GSCA outer estimation
#'
#'@description
#'
#'This implements the second step of the GSCA estimation describe by Hwang & Takane (2004). GSCA outer estimation should
#'be used only with GSCA inner estimation. 
#'
#'@details
#'
#'The second step of GSCA estimation method, as describe by Hwang & Takane (2004), involves calculation of new
#'weights given the regression estimates form the first step. In the second step, the following function is
#'minimized (Hwang & Takane, 2004, eq. 7, first row):
#'
#'\deqn{SS(Z[V-\Lambda])}{SS(Z[V-Lambda])}
#'
#'Because \eqn{\Lambda}{Lambda} is defined as \eqn{WA}{WA}, the function to be minimized is
#'identical to the first step function (Hwang & Takane, 2004, eq. 4, first row):
#'
#'\deqn{SS(ZV-ZWA)}{SS(ZV-ZWA)}
#'
#'In the second step, this function is minimized in respect to weights \eqn{W}{W} and \eqn{V}{V}.
#'This involves estimating each regression ananalysis in the model including
#'regressions between the composites and from composites to indicators and to minimize the sum of all OLS
#'dicrepancy functions simultaneously. Because one weight can be included in many regressions, these equations
#'must be estimated simultaneously. The minimization algoritm is the Nelder-Mead algorithm implemented in the
#'\code{\link{optim}} function.
#'
#'The GSCA algoritm described by Hwang & Takane (2004) allows some indicators to be excluded from the second
#'step, but in this implementation all indicators are always used so that each weight relation described in
#'\code{W.mod} has always a corresponding regression relationship from a composite to a variable in the second
#'step of GSCA estimation.
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param E Inner weight matrix. A square matrix of inner estimates between the composites.
#'
#'@param W.mod A matrix specifying the weight relationships and their starting values.
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.mod}.
#'
#'@references
#'Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. \emph{Psychometrika}, 69(1), 81–99. doi:10.1007/BF02295841
#'
#'@example example/gsca-example.R
#'
#'@family outer estimators
#'
#'@export

outer.GSCA <- function(S, W, E, W.mod){
	
	
	Wpattern <- W.mod!=0
	
	# We will start by creating a function that returns the sum of squares of all the regressions
	
	ssFun <- function(Wvect){
		
		W[Wpattern] <- Wvect
		
		# Start by standardizing the weights because optim does not know that the
		# weights must result in a composite with unit variance
		
		W <- scaleWeights(S, W)
		
		# Calculate the covariance matrix between indicators and composites
		IC <- W %*% S
		
		# Calculate the composite covariance matrix
		C <- IC %*% t(W)
		
		# Sum of squares from regressions from composites to indicators
		
		# Sum of squares is the sum of residual variances. For indicators, residual variance is the
		# difference between covariance between the composite and the indicator variances. 
		# we will take a sum of these 
		
		indicatorIndices <- col(W.mod)[Wpattern]
		ss_indicators <-sum (diag(S)[indicatorIndices]-IC[Wpattern])
		
		# Sum of squares from regressions between composites
		
		# Sum of squares is the sum of residual variances. For composites, residual variance is 
		# 1 - R2 and we will take a sum of these 
		
		R2 <- rowSums(E * C)
		ss_composites <- sum(1-R2)
		
		return(ss_composites + ss_indicators)
	}	
	
	# The previous weights are used as starting values
	
	start <- W[Wpattern]
	
	# Then we will find the minimum using optim
	
	Wvect <- optim(start, ssFun)$par
	
	# Update the weights based on the estimated parameters
	
	W[Wpattern] <- Wvect
	
	# Finally standardize the weights because optim does not know that the
	# weights must result in a composite with unit variance
	
	W <- scaleWeights(S, W)
	
	return(W)
}

# =========== Utility functions ===========

#
# Scales the weight matrix so that the resulting composites have unit variance
#

scaleWeights <- function(S, W){
	
	# Calculate the variances of the unscaled composites
	
	var_C_unscaled <- diag(W %*% S %*% t(W))
	
	# Scaling the unscaled weights and return
	
	return(diag(x = 1/sqrt(var_C_unscaled)) %*% W)
}

lavaanParTableToNativeFormat <- function(partable){
	
	library(lavaan)
	
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


