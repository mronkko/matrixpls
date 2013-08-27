# =========== Main functions ===========

#'@title Basic results for Partial Least Squares Path Modeling as a vector
#'
#'@description
#'Estimates a weight matrix using the PLS algorithm
#
#'@details
#'\code{matrixpls.weights} performs the PLS weighting algorithm by calling the 
#'inner and outer estimators iteratively until either the convergence criterion or
#'maximum number of iterations is reached and provides the results in a matrix.
#'
#'The argument \code{weightRelations} is a matrix of zero and non-zeros that indicates
#'how the indicators are combined to form the composites. In Wold's original papers, 
#'each indicator was allowed to contribute to only one composite, but MatrixPLS does
#'not have this restriction.
#'\code{weightRelations} will contain a non-zero value when a indicator variable on the
#'column \code{j} is a part of the composites variable on row \code{i}, zero otherwise. 
#'The values are used as starting values for the weighting algorithm.\cr
#'
#'@param S Covariance matrix of the data.
#
#'@param Model There are four options for this argument: 1. SimSem object created by model
#'command from simsem package, 2. lavaan script, lavaan parameter table, or a list that
#'contains all argument that users use to run lavaan (including cfa, sem, lavaan), 3.
#'MxModel object from the OpenMx package, or 4. A list containing three matrices
#'\code{inner}, \code{reflective}, and \code{formative} defining the free regression paths
#'in the model.
#'
#'@param weightRelations real matrix representing the weight relations and starting weights
#'(i.e. the how the indicators are combined to form the composite variables).
#'
#'@param ... All other parameters are passed through to \code{matrixpls.weights}
#'
#'@param validateInput A boolean indicating whether the validity of the input matrix
#'and the parameter values should be tested
#
# 
#'@return A named vector of class \code{matrixpls} containing parameter estimates followed by weights
#'@export


matrixpls <- function(S, model, weightRelations = NULL, parameterEstimator = matrixpls.parameterEstimator.regression, ..., validateInput = TRUE) {
	
	# TODO: Validate input.
	
	nativeModel <- parseModelToNativeFormat(model)
	
	if(is.null(weightRelations)) weightRelations <- defaultWeightRelationsWithModel(nativeModel)
	
	# Calculate weights
	W <- matrixpls.weights(S, nativeModel$inner, weightRelations, ..., validateInput = validateInput)
	
	# Apply the parameter estimator and return the results
	estimates <- parameterEstimator(S, W, model)
	
	# Construct the return vector by combining estimates with weights in a named vector
	indices <- weightRelations!=0
	WVect <- W[indices]
	names(WVect) <- paste(rownames(nativeModel$formative)[row(W)[indices]],"=+",colnames(nativeModel$formative)[col(W)[indices]], sep="")
	
	ret <- c(estimates,WVect)
	
	attr(ret,"iterations") <- attr(W,"iterations")
	attr(ret,"converged") <- attr(W,"converged")
	class(ret) <-("matrixpls")
	
	return(ret)
}

#'@title Basic results for Partial Least Squares Path Modeling as a vector
#'
#'@description
#'Estimates a weight matrix using the PLS algorithm
#
#'@details
#'\code{matrixpls.weights} performs the PLS weighting algorithm by calling the 
#'inner and outer estimators iteratively until either the convergence criterion or
#'maximum number of iterations is reached and provides the results in a matrix.
#'
#'The argument \code{innerMatrix} is a matrix of zeros and ones that indicates
#'the regression relationships between composites. In Wold's original papers, 
#'the model was constrained to be recursive implying a lower
#'triangular matrix, but MatrixPLS does not have this restriction.
#'\code{innerMatrix} will contain a 1 when a variable on the column \code{j}
#'has a regression path toward the variable on row \code{i}, 0 otherwise. \cr
#'
#'The argument \code{weightRelations} is a matrix of zero and non-zeros that indicates
#'how the indicators are combined to form the composites. In Wold's original papers, 
#'each indicator was allowed to contribute to only one composite, but MatrixPLS does
#'not have this restriction.
#'\code{weightRelations} will contain a non-zero value when a indicator variable on the
#'column \code{j} is a part of the composites variable on row \code{i}, zero otherwise. 
#'The values are used as starting values for the weighting algorithm.\cr
#'
#'@param S Covariance matrix of the data.
#
#'@param innerMatrix A square matrix of ones and zeros representing the
#'inner model (i.e. the path relationships between the composite variables).
#
#'@param weightRelations real matrix representing the weight relations and starting weights
#'(i.e. the how the indicators are combined to form the composite variables).
#'
#'@param outerEstimators A function or a list of functions used for outer estimation. If
#'the value of this parameter is a function, the same function is applied to all
#'composites. If the value of the parameter is a list, the composite \code{n} is estimated
#'with the estimator in the \code{n}th position in the list.
##
#'@param tol Decimal value indicating the tolerance criterion for the
#'iterations (\code{tol=0.00001} by default). Must be a positive number.
#
#'@param iter An integer indicating the maximum number of iterations
#'(\code{iter=100} by default).
#'
#'@param validateInput A boolean indicating whether the validity of the input matrix
#'and the parameter values should be tested
#'
#'@return An object of class \code{"matrixpls.weights"}, which is a matrix consisting of weights and the following attributes: 
#'@return \item{iterations}{Number of iterations performed}
#'@return \item{converged}{A boolean indicating if the algoritm converged}
#'@export

matrixpls.weights <- function(S, innerMatrix, weightRelations,
															outerEstimators = matrixpls.outerEstimator.modeA, 
															innerEstimator = matrixpls.innerEstimator.path,
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
		
		# innerMatrix must be a square matrix consisting of ones and zeros and zero diagonal
		# and all variables must be linked to at least one other variable
		assert_is_matrix(innerMatrix)
		assert_is_identical_to_true(nrow(innerMatrix)==ncol(innerMatrix))
		assert_all_are_true(innerMatrix==0 | innerMatrix == 1)
		assert_all_are_true(diag(innerMatrix)==0)
		
		# weightRelations must be a real matrix and each indicators must be
		# linked to at least one composite and each composite at least to one indicator
		assert_is_matrix(weightRelations)
		assert_all_are_real(weightRelations)
		assert_all_are_true(apply(weightRelations!=0,1,any))
		assert_all_are_true(apply(weightRelations!=0,2,any))
		
		# the number of rows in weightRelations must match the dimensions of other matrices
		assert_is_identical_to_true(nrow(innerMatrix)==nrow(weightRelations))
		assert_is_identical_to_true(ncol(S)==ncol(weightRelations))
		
		# outerEstimators must be a list of same length as number of rows in innerMatrix or
		# a function
		
		if(is.list(outerEstimators)){
			assert_is_identical_to_true(length(outerEstimators) == nrow(weightRelations))
			for(outerEstimator in outerEstimators){
				assert_is_function(outerEstimator)
			}
		}
		else{
			assert_is_function(outerEstimators)
		}
		
		assert_is_function(innerEstimator)
		
		# tol must be positive
		assert_all_are_positive(tol)
		
		#iter must not be negative
		assert_all_are_positive(iter)
	}
	
	# =========== Start of iterative procedure ===========
	
	# The weight matrix
	
	W <- scaleWeights(S, weightRelations)
	
	for (iteration in 1:iter) {
		
		# Get new inner weights from inner estimation
		
		E <- innerEstimator(S, W, innerMatrix)
		
		# Get new weights from outer estimation
		
		if(is.list(outerEstimators)){
			stop("Using multiple outer estimators has not been implemented yet")
		}
		else{
			W_new<-outerEstimators(S, W, E, weightRelations)
			W_new <- scaleWeights(S, W_new)
		}	
		
		# Check convergence
		
		if(max(abs(W_new-W)) <= tol){
			attr(W_new,"iterations") <- iteration
			attr(W_new,"converged") <- TRUE
			class(W_new) <-("matrixpls.weights")
			return(W_new)			
		} 
		
		# Prepare for the next round
		
		W <- W_new
		
	}
	
	# Reaching the end of the loop without convergence
	
	attr(W,"iterations") <- iter
	attr(W,"converged") <- FALSE
	class(W) <-("matrixpls.weights")
	return(W)
}

# =========== Parameter estimators ===========

#'@title PLS parameter estimation with separate regression analyses
#'
#'@description
#'Estimates the model parameters with weighted composites using separate OLS regression and estimating
#'measurement model and latent variable model separately. This is the standard way of estimation in the PLS literature.
#'
#'
#'@details
#'~~~ Explain how PLS estimates are calculated after weights have been ~~~
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param Model There are four options for this argument: 1. SimSem object created by model
#'command from simsem package, 2. lavaan script, lavaan parameter table, or a list that
#'contains all argument that users use to run lavaan (including cfa, sem, lavaan), 3.
#'MxModel object from the OpenMx package, or 4. A list containing three matrices
#'\code{inner}, \code{reflective}, and \code{formative} defining the free regression paths
#'in the model.
#'
#'@return A named vector of parameter estimates
#'
#'@export

matrixpls.parameterEstimator.regression <- function(S, W, model){

	nativeModel <- parseModelToNativeFormat(model)
	
	results <- c()

	# Create the composite covariance matrix
	C <- W %*% S %*% t(W)

	# Calculate the covariance matrix between indicators and composites
	IC <- W %*% S
	
	# Join these to form the covariance matrix of indicators and composites
	
	IC_full <- rbind(cbind(C,IC),
									 cbind(t(IC),S))
	
	innerRegressionIndices <- which(nativeModel$inner==1)
	
	
	# Choose the specified values and add names
	if(length(innerRegressionIndices)>0){
				
		inner <- regressionsWithCovarianceMatrixAndModelPattern(C,nativeModel$inner)
		innerVect <- inner[innerRegressionIndices]
		names(innerVect) <- paste(rownames(inner)[row(inner)[innerRegressionIndices]],"~",
														colnames(inner)[col(inner)[innerRegressionIndices]], sep="")
		
		results <- c(results, innerVect)
	}
	
	
	
	if(any(nativeModel$reflective ==1)){
		
		reflectiveModel <- rbind(matrix(0,nrow(C),ncol(IC_full)),
														 cbind(nativeModel$reflective, matrix(0,nrow(S),ncol(S))))
		
		reflective <- regressionsWithCovarianceMatrixAndModelPattern(IC_full, reflectiveModel)

		reflectiveRegressionIndices <- which(reflectiveModel==1)
		
		# Choose the specified values and add names
		reflectiveVect <- reflective[reflectiveRegressionIndices]
		names(reflectiveVect) <- paste(colnames(reflective)[col(reflective)[reflectiveRegressionIndices]],"=~",
														rownames(reflective)[row(reflective)[reflectiveRegressionIndices]], sep="")
		
		results <- c(results,reflectiveVect)
	}
	
	formativeRegressionIndices <- which(nativeModel$formative==1)
	
	if(any(nativeModel$formative ==1)){

		formativeModel <- rbind(cbind(matrix(0,nrow(C),ncol(C)), nativeModel$formative),
														 cbind(matrix(0,nrow(S),ncol(IC_full))))
		
		formative <- regressionsWithCovarianceMatrixAndModelPattern(IC_full,formativeModel)
	
		formativeRegressionIndices <- which(formativeModel==1)
		
		# Choose the specified values and add names
		formativeVect <- formative[formativeRegressionIndices]
		names(formativeVect) <- paste(rownames(formative)[row(formative)[formativeRegressionIndices]],"~",
														colnames(formative)[col(formative)[formativeRegressionIndices]], sep="")
		
		results <- c(results,formativeVect)
	}
	
	return(results)
}


# =========== Inner estimators ===========

#'@title PLS inner estimation with the centroid scheme
#'
#'@description
#'Calculates a set of inner weights based on the centroid scheme
#
#'@details
#'~~~ Explain  the centroid scheme here ~~~
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param innerMatrix A square matrix specifying the relationships of the latent variables in the model
#'
#'@return A matrix of unscaled inner weights with the same dimesions as \code{innerMatrix}
#'
#'@export

matrixpls.innerEstimator.centroid <- function(S, W, innerMatrix){

	# Centroid is just the sign of factor weighting
	
	E <- sign(matrixpls.innerEstimator.factor(S, W, innerMatrix))
	
	return(E)
}

#'@title PLS inner estimation with the path scheme
#'
#'@description
#'Calculates a set of inner weights based on the path scheme
#
#'@details
#'~~~ Explain  the path scheme here ~~~
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param innerMatrix A square matrix specifying the relationships of the latent variables in the model
#'
#'@return A matrix of unscaled inner weights with the same dimesions as \code{innerMatrix}
#'
#'@export

matrixpls.innerEstimator.path <- function(S, W, innerMatrix){
	
	# Start with the factor weights
	E <- matrixpls.innerEstimator.factor(S, W, innerMatrix)
	
	E <- regressionsWithCovarianceMatrixAndModelPattern(S,E)
	
	return(E)
}

#'@title PLS inner estimation with the factor scheme
#'
#'@description
#'Calculates a set of inner weights based on the factor scheme
#
#'@details
#'~~~ Explain  the factor scheme here ~~~
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param innerMatrix A square matrix specifying the relationships of the latent variables in the model
#'
#'@return A matrix of unscaled inner weights with the same dimesions as \code{innerMatrix}
#'
#'@export

matrixpls.innerEstimator.factor <- function(S, W, innerMatrix){

	# Create the composite covariance matrix
		
	C <- W %*% S %*% t(W)

	# Calculate unscaled inner weights using the centroid weighting scheme
	
	E <- C * (innerMatrix | t(innerMatrix))
	
	# If we have LVs that are not connected to any other LVs, use identity scheme as fallback
	diag(E)[rowSums(E) == 0] <- 1
	
	return(E)
} 

#'@title PLS inner estimation with the identity scheme
#'
#'@description
#'
#'The identity scheme is an inner weighting scheme where the inner weight matrix is an identity matrix.
#
#'@details
#'
#'This scheme is not commonly discussed in the current PLS literature, but it is a special case 
#'that is analyzed in the early PLS literature and currently implmented in at least WarpPLS.
#'
#'The identity scheme with Mode A outer estimation converges toward the first principal component of each indicator block.
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param innerMatrix A square matrix specifying the relationships of the latent variables in the model
#'
#'@return A matrix of unscaled inner weights with the same dimesions as \code{innerMatrix}
#'
#'@references
#'Wold, H. (1966). Nonlinear estimation by iterative least squares procedures. Research Papers in Statistics: Festschrift for J. Neyman, 411–444.
#'
#'Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. In K. G. Jöreskog & S. Wold (Eds.), Systems under indirect observation: causality, structure, prediction (pp. 1–54). Amsterdam: North-Holland.
#'@export



matrixpls.innerEstimator.identity <- function(C, innerMatrix){
	return(diag(nrow(innerMatrix)))
}

# =========== Outer estimators ===========

#'@title PLS outer estimation with Mode A
#'
#'@description
#'
#'Performs Mode A outer estimation by calculating correlation matrix between the indicators and composites.
#
#'@details
#'
#'~~~ Explain details here ~~
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param E Inner weight matrix. A square matrix of inner estimates between the composites.
#'
#'@param weightRelations A matrix specifying the weight relationships and their starting values.
#'
#'@return A matrix of unscaled outer weights with the same dimesions as \code{weightRelations}
#'
#'@export


matrixpls.outerEstimator.modeA <- function(S, W, E, weightRelations){

	# Calculate the covariance matrix between indicators and composites
	W_new <- E %*% W %*% S

	# Set the non-existing weight relations to zero
	W_new[weightRelations == 0] <- 0
	
	return(W_new)
}

#'@title PLS outer estimation with Mode B
#'
#'@description
#'
#'Performs Mode B outer estimation by calculating correlation matrix between the indicators and composites.
#
#'@details
#'
#'~~~ Explain details here ~~
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param E Inner weight matrix. A square matrix of inner estimates between the composites.
#'
#'@param weightRelations A matrix specifying the weight relationships and their starting values.
#'
#'@return A matrix of unscaled outer weights with the same dimesions as \code{weightRelations}
#'
#'@export


matrixpls.outerEstimator.modeB <- function(S, W, E, weightRelations){

	# Calculate the covariance matrix between indicators and composites
	IC <- E %*% W %*% S
	
	# Set up a weight pattern
	W_new <- iflse(weightRelations==0,0,1)
	
	# Do the outer model regressions
	W_new <- regressionsWithCovarianceMatrixAndModelPattern(IC,W_new)
	
	return(W_new)
	
}

#'@title PLS outer estimation with fixed weights
#'
#'@description
#'
#'Returns the starting weights as outer estimates.
#
#'@details
#'
#'~~~ Explain details here ~~
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param E Inner weight matrix. A square matrix of inner estimates between the composites.
#'
#'@param weightRelations A matrix specifying the weight relationships and their starting values.
#'
#'@return A matrix of unscaled outer weights with the same dimesions as \code{weightRelations}
#'
#'@export

matrixpls.outerEstimator.fixedWeights <- function(S, W, E, weightRelations){
	return(weightRelations)
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
	
	# Parse the variables
	latentVariableNames <- unique(factorLoadings$lhs)
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
	
	formativeRegressions <- regressions[regressions$rhs %in% observedVariableNames & 
																	 	regressions$lhs %in% latentVariableNames,]


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
	} else if (is(model, "MxModel")) {
		browser()
	} else if (is(model, "SimSem")) {
		browser()
	} else {
		stop("Please specify an appropriate object for the 'model' argument: simsem model template, lavaan script, lavaan parameter table, list of options for the 'lavaan' function, or a list containing three matrices in the MatrixPLS native model format.")
	}
}

#
# Returns a matrix with composites on rows and observed on columns
#

defaultWeightRelationsWithModel <- function(model){
		
	nativeModel <- parseModelToNativeFormat(model)
	weightRelations <- matrix(0,nrow(nativeModel$formative), ncol(nativeModel$formative))
	weightRelations[nativeModel$formative!=0] <- 1 
	weightRelations[t(nativeModel$reflective)!=0] <- 1 
	
	return(weightRelations)
	
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

	assert_is_symmetric_matrix(S)
	
	for(row in 1:nrow(model)){
		
		independents <- which(model[row,]!=0)
		
		if(length(independents)>0){
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


