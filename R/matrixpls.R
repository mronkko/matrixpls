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
	
	# Copy all non-standard attributes from the weight and estimation algorithms
	
	allAttributes <- c(attributes(W),attributes(estimates))
	 
	for(a in setdiff(names(allAttributes), c("dim", "dimnames", "class", "names"))){
		attr(ret,a) <- allAttributes[[a]]
	}

	attr(ret,"W") <- W
	
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
#'@param convergenceCheckFunction A function that takes the old and new weight matrices and
#'returns a scalar that is compared against \code{tol} to check for convergence
#'
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
															convergenceCheckFunction = function(W,W_new){max(abs(W_new) - abs(W))},
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
		
		if(convergenceCheckFunction(W,W_new) < tol){
			attr(W_new,"iterations") <- iteration
			attr(W_new,"converged") <- TRUE
			class(W_new) <-("matrixpls.weights")
			rownames(W_new) <- rownames(innerMatrix)
			return(W_new)			
		} 
		
		# Prepare for the next round
		
		W <- W_new
		
	}
	
	# Reaching the end of the loop without convergence
	
	attr(W,"iterations") <- iter
	attr(W,"converged") <- FALSE
	class(W) <-("matrixpls.weights")
	rownames(W) <- rownames(innerMatrix)
	
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

	return(matrixpls.parameterEstimator.internal_generic(S,W,model,
																											 regressionsWithCovarianceMatrixAndModelPattern,
																											 regressionsWithCovarianceMatrixAndModelPattern,
																											 regressionsWithCovarianceMatrixAndModelPattern))
}

#'@title Generic parameter estimation that prepares the data and applies the estimation functions
#'
#'@description
#'Estimates the model parameters with weighted composites using separate estimations for reflective, formative, and latent path relationships
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
#'@param pathEstimator A function for estimating the path relationships
#'
#'@param formativeEstimator A function for estimating the formative relationships
#'
#'@param reflectiveEstimator A function for estimating the reflective relationships
#'
#'@return A named vector of parameter estimates
#'

matrixpls.parameterEstimator.internal_generic <- function(S, W, model, pathEstimator, formativeEstimator, reflectiveEstimator){
	
	nativeModel <- parseModelToNativeFormat(model)
	
	results <- c()
	
	# Calculate the composite covariance matrix
	C <- W %*% S %*% t(W)
	
	# Calculate the covariance matrix between indicators and composites
	IC <- W %*% S
	
	# Join these to form the covariance matrix of indicators and composites
	
	IC_full <- rbind(cbind(C,IC),
									 cbind(t(IC),S))
	
	innerRegressionIndices <- which(nativeModel$inner==1)
	
	
	# Choose the specified values and add names
	if(length(innerRegressionIndices)>0){
		
		inner <- pathEstimator(C,nativeModel$inner)
		innerVect <- inner[innerRegressionIndices]
		names(innerVect) <- paste(rownames(inner)[row(inner)[innerRegressionIndices]],"~",
															colnames(inner)[col(inner)[innerRegressionIndices]], sep="")
		
		results <- c(results, innerVect)
	}
	
	
	
	if(any(nativeModel$reflective ==1)){
		
		reflectiveModel <- rbind(matrix(0,nrow(C),ncol(IC_full)),
														 cbind(nativeModel$reflective, matrix(0,nrow(S),ncol(S))))
		
		reflective <- reflectiveEstimator(IC_full, reflectiveModel)
		
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
		
		formative <- formativeEstimator(IC_full,formativeModel)
		
		formativeRegressionIndices <- which(formativeModel==1)
		
		# Choose the specified values and add names
		formativeVect <- formative[formativeRegressionIndices]
		names(formativeVect) <- paste(rownames(formative)[row(formative)[formativeRegressionIndices]],"~",
																	colnames(formative)[col(formative)[formativeRegressionIndices]], sep="")
		
		results <- c(results,formativeVect)
	}
	
	# Store these in the result object
	attr(results,"C") <- C
	attr(results,"IC") <- IC
	attr(results,"beta") <- inner
	
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
	
	# Calculate the composite covariance matrix	
	C <- W %*% S %*% t(W)
	
	E <- regressionsWithCovarianceMatrixAndModelPattern(C, innerMatrix)
	
	# Use  correlations for inverse relationships for non-reciprocal paths
	inverseRelationships <- t(innerMatrix) & ! innerMatrix
	E[inverseRelationships] <- C[inverseRelationships]
	
	# If we have LVs that are not connected to any other LVs, use identity scheme as fallback
	diag(E)[rowSums(E) == 0] <- 1
	
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

	# Calculate the composite covariance matrix
		
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



matrixpls.innerEstimator.identity <- function(S, W, innerMatrix){
	return(diag(nrow(innerMatrix)))
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
#'including also the relationships from composites toward indicators in the first step. The term component
#'is used in a very generic way to refer to a variable that is a weighted composite of the raw data.
#'
#'The implementation of GSCA in MatrixPLS differens from the Hwang & Takane (2004) version in that during the
#'first step, only regressions between components are estimated. The reason for this is that the
#'relationhips from the components to indicators are simple regressions that are simply the covariances between
#'the indicators and components. Since these covariances need to be calculated in the second step, it is more
#'efficient to not calculate them during the first step.
#'
#'This algorithm is therefore identical to the PLS path weighting scheme with the exception that correlations
#'are not used for inverse relationships and there is no falling back to identity scheme for components
#'that are not connected to other components
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
#'Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. Psychometrika, 69(1), 81–99. doi:10.1007/BF02295841
#'
#'@examples
#'  
#'  # Run the example from plspm package using GSCA estimation
#'  
#'  library(plspm)
#'  
#'  # load dataset satisfaction
#'  data(satisfaction)
#'  
#'  # inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0)
#'  LOY = c(1,0,0,0,1,0)
#'  
#'  inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'  
#'  # reflective relationships
#'  reflective <- matrix(0, 27, 6)
#'  for(lv in 1:6){
#'  	reflective[list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)[[lv]],lv] <-1
#'  }
#'  # no formative relationships
#'  formative <- matrix(0, 6, 27)
#'  
#'  model <- list(inner=inner, reflective=reflective, formative=formative)
#'  
#'  # apply MatrixPLS with GSCA estimators
#'  
#'  matrixpls(cor(satisfaction[,1:27]), model, 
#'  					outerEstimators = matrixpls.outerEstimator.GSCA, 
#'  					innerEstimator = matrixpls.innerEstimator.GSCA)
#'
#'@export

matrixpls.innerEstimator.GSCA <- function(S, W, innerMatrix){
		
	# Calculate the composite covariance matrix
	C <- W %*% S %*% t(W)
	
	E <- regressionsWithCovarianceMatrixAndModelPattern(C,innerMatrix)
		
	return(E)
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
#'weights given the regression estimates form the first step.The term component
#'is used in a very generic way to refer to a variable that is a weighted composite of the raw data.
#'
#'In the first step (Hwang & Takane, 2004, eq. 7)
#'
#'\deqn{latex}{SS(\textbf{Z}[\textbf{V}-\textbf{\Lambda}])}
#'
#'Because \eqn{latex}{\textbf{\Lambda}} is defined as \eqn{latex}{\textbf{WA}}, the function to be minimized is
#'identical to the first step function (Hwang & Takane, 2004, eq. 4)
#'
#'\deqn{latex}{SS(\textbf{ZV}-\textbf{ZWA})}
#'
#'In the second step, this function was minimized in respect to W and V, which are returned and used as fixed
#'parameters the next estimation step. This involves estimating each regression ananalysis in the model including
#'regressions between the components and from components to indicators and to minimize the sum of all OLS
#'dicrepancy functions simultaneously. Because one weight can be included in many regressions, these equations
#'must be estimated simultaneously.
#'
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
#'@references
#'Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. Psychometrika, 69(1), 81–99. doi:10.1007/BF02295841
#'
#'@examples
#'
#'  # Run the example from plspm package using GSCA estimation
#'  
#'  library(plspm)
#'  
#'  # load dataset satisfaction
#'  data(satisfaction)
#'  
#'  # inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0)
#'  LOY = c(1,0,0,0,1,0)
#'  
#'  inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'  
#'  # reflective relationships
#'  reflective <- matrix(0, 27, 6)
#'  for(lv in 1:6){
#'  	reflective[list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)[[lv]],lv] <-1
#'  }
#'  # no formative relationships
#'  formative <- matrix(0, 6, 27)
#'  
#'  model <- list(inner=inner, reflective=reflective, formative=formative)
#'  
#'  # apply MatrixPLS with GSCA estimators
#'  
#'  matrixpls(cor(satisfaction[,1:27]), model, 
#'  					outerEstimators = matrixpls.outerEstimator.GSCA, 
#'  					innerEstimator = matrixpls.innerEstimator.GSCA)
#'
#'
#'@export

matrixpls.outerEstimator.GSCA <- function(S, W, E, weightRelations){
	
	
	Wpattern <- weightRelations!=0
	
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
		
		indicatorIndices <- col(weightRelations)[Wpattern]
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


