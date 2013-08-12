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
#'@return An object of class \code{"matrixpls.weights"}. 
#'@return \item{weights}{Weight matrix where composites are on rows and indicators are on columns}
#'@return \item{iterations}{Number of iterations performed}
#'@return \item{converged}{A boolean indicating if the algoritm converged}
#'@export

matrixpls.weights <- function(S, innerMatrix, weightRelations,
	outerEstimators = matrixpls.outerEstimator.modeA, 
	innerEstimator = matrixpls.innerEstimator.centroid,
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
		assert_all_are_true(rowSums(innerMatrix)+colSums(innerMatrix)>0)
		
		# weightRelations must be a real matrix and each indicators must be
		# linked to at least one composite and each composite at least to one indicator
		assert_is_matrix(weightRelations)
		assert_all_are_real(weightRelations)
		assert_all_are_true(apply(1,innerMatrix!=0,any))
		assert_all_are_true(apply(2,innerMatrix!=0,any))
		
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
			assert_is_function(outerEstimator)
		}
		
		assert_is_function(innerEstimator)

		# tol must be positive
		assert_all_are_positive(tol)
		
		#iter must not be negative
		assert_all_are_positive(iter)
	}
	
	# =========== Start of iterative procedure ===========

	# The weight matrix
	
	W <- scaleWeights(weightRelations)

	for (iteration in 1:iter) {
	
		# Get new inner weights from inner estimation

		E <- innerEstimator(S, W, innerMatrix)
		
		# Get new weights from outer estimation
		
		if(is.list(outerEstimator)){
			stop("Not implemented")
		}
		else{
			W_new<-outerEstimator()
			W_new <- scaleWeights(outerEstimator())
		}	
		
		# Check convergence
		
		if(max(abs(W_new-W)) <= tol) return(list(weightMatrix = W, iterations = iteration,
			converged = TRUE))
		
		# Prepare for the next round
		
		W <- W_new
	
	}
	
	# Reaching the end of the loop without convergence
	
	return(list(weightMatrix = W, iterations = iter, converged = FALSE))
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
#'MxModel object from the OpenMx package, or 4. A square matrix of ones and zeros 
#'representing the free regression paths in the model.
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


matrixpls <- function(S, model, weightRelations, ..., validateInput = TRUE) {

	stop("Not implemented")
	
#	# =========== inputs setting ===========
#	
#	# Names of the composites
#	lvs.names <- colnames(innerMatrix)
#	mv.names <- colnames(S)
#	
#	# A binary vector showing which composites are endogenous
#	endo <- rowSums(innerMatrix) > 0
#	
#	# The estimation uses the following matrices
#	# 
#	# innerMatrix: a p by p matrix specifying the inner model weightRelations: a n by p matrix specifying the outer model
#	# 
#	# S: n by n indicator covariance matrix
#	# 
#	# W: n by p matrix of outer weights C: is a p by p covariance matrix of composites E: is a p by p matrix of inner weights
#	# 
#	# where
#	# 
#	# n: number of indicators p: number of composites
#	
#	n <- nrow(S)
#	p <- nrow(innerMatrix)
#	
#	weightRelations <- sapply(outer_list, function(x) match(1:n, x, nomatch = 0) > 0)
#	
#	S <- cov2cor(S)
#	
#	# =========== Stage 2: Path coefficients ===========
#	
#	
#	if (unbiasedLoadings) {
#	load.est <- matrix(0, nrow(S), ncol(inner))
#	# Loop over the indicator blocks and perform EFAs
#	
#	for (col in 1:length(outer)) {
#				load.est[outer[[col]], col] <- fa(S[outer[[col]], outer[[col]]])$loadings[, 1]
#				
#			}
#	} else {
#	load.est = W_unscaled
#	}
#	
#	if (disattenuate) {
#	
#	if (unbiasedCompositeReliability) {
#	# Calculate composite reliabilities
#	# 
#	# Aguirre-Urreta, M. I., Marakas, G. M., & Ellis, M. E. (in press). Measurement of Composite Reliability in Research Using Partial Least Squares: Some Issues and an Alternative Approach. The DATA BASE
#	# for Advances in Information Systems.
#	
#	reliabilities <- colSums(W * load.est)^2
#	disattenuationMatrix <- sqrt(reliabilities %o% reliabilities)
#	diag(disattenuationMatrix) <- C <- C * disattenuationMatrix
#	} else {
#	stop("Disattennuation has been implemented only for unbiased composite reliability")
#	}
#	}
#	
#	# A p x p lower diagonal matrix of path estimates
#	Path <- innerMatrix
#	R2 <- rep(0, p)
#	
#	# Loop over the endogenous variables and do the regressions
#	
#	for (aux in which(endo)) {
#	indeps <- which(innerMatrix[aux, ] == 1)
#	coefs <- solve(C[indeps, indeps], C[indeps, aux])
#	Path[aux, indeps] <- coefs
#	}
#	
#	
#	# =========== Results ===========
#	
#	
#	return(c(W[which(weightRelations==1)], load.est[which(weightRelations==1)], Path[which(innerMatrix==1)], C[lower.tri(innerMatrix)]))
	
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
	
	# Loop over each row of innerMatrix
	
	for(row in 1:nrow(innerMatrix)){

		# Which LVs does this LV depend on
		
		independents <- which(innerMatrix[row,]!=0)
		
		# If there are at least (regression with one independent equals correlation, which
		# we already have from the factor weights)

		if(length(independents)>1){
			# Use the regression coefficients as weights
			coefs <- solve(C[independents,independents],C[independents,row])
			E[row,independents] <- coefs
		}
	}
	
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
		
	C <- t(W) %*% S %*% W

	# Calculate unscaled inner weights using the centroid weighting scheme
	
	E <- C * (innerMatrix | t(innerMatrix))
	
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
	W_new <- S %*% W %*% E

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
	IC <- S %*% W %*% E
	
	W_new = matrix(0,nrow(weightRelations),ncol(weightRelations))
	
	for(row in 1:nrow(weightRelations)){

		# Which indicators contribute to this LV
		
		independents <- which(weightRelations[row,]!=0)
		
		# If there are at least (regression with one independent equals correlation, which
		# we already have from the factor weights)

		if(length(independents)>1){
			# Use the regression coefficients as weights
			coefs <- solve(IC[independents,independents],IC[independents,row])
			W_new[row,independents] <- coefs
		}
	}
	
	return(E)
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
		
	var_C_unscaled <- diag(t(W) %*% S %*% W)
		
	# Scaling the unscaled weights and return
		
	return(W  %*% diag(x = 1/sqrt(var_C_unscaled)))
}

