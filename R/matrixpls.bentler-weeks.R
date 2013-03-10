#'@title Bentler-Weeks structural equation model
#'
#'@description
#'Calculates the matrices of Bentler-Weeks type structural equation model.

#
#'@details
#'
#'@param beta beta matrix of the Bentler-Weeks model (regression paths between dependent
#'variables)
#'@param gamma gamma matrix of the Bentler-Weeks model (regression paths between dependent
#'and independent variables)
#'@param Phi Phi matrix of the Bentler-Weeks model (independent variable covariance matrix)
#'@param observed.independents A vector containing indices of observed independent variables
#'@param observed.dependents A vector containing indices of observed dependent variables

#'@return An object of class \code{'matrixpls.bentler-weeks'} containing.
#'@return \item{estimates}{A matrix containing the PLS results as columns in the following order: weights, loading
#'estimates, path estimates, lower triangular of the composite correlation matrix}
#'@return \item{beta}{The beta matrix}
#'@return \item{gamma}{The gamma matrix}
#'@return \item{Phi}{The Phi matrix}
#'@return \item{Sigma}{The population covariance matrix for the data based on the Bentler-Weeks model}
#'@author Mikko R\303\266nkk\303\266
#'
#'@references # Bentler, P. M., & Weeks, D. G. (1980). \emph{Linear structural equations
#'with latent variables.} Psychometrika, 45(3), 289\342\200\223308. doi:10.1007/BetaF02293905
#'


matrixpls.bentler-weeks <- function(beta = NULL, gamma = NULL, Phi = NULL, observed.independents = c(), observed.dependents = c()) {
	
	# Calculate  Sigma based on the Bentler-Weeks model
	
	# Matrix orders
	
	p = length(observed.dependents)
	q = length(observed.independents)
	n = nrow(Phi)
	m = nrow(beta)
	r = p + q
	s = m + n
	
	# The matrices
	
	Beta <- matrix(0, s, s)
	Gamma <- matrix(0, s, n)
	G <- matrix(0, r, s)
	
	# Identity matrix
	
	I <- diag(nrow(Beta))
	
	# Populate the parameter matrices
	
	Gamma[1:m, ] <- gamma
	Gamma[(m + 1):s, ] <- diag(n)
	
	Beta[1:m, 1:m] <- beta
	
	# Populate the selection matrix
	
	
	if (p != 0) {
		Gy <- outer(observed.dependents, 1:m, function(x, y) x == y)
		G[1:p, 1:m] <- Gy
	}
	if (q != 0) {
		Gx <- outer(observed.independents, 1:n, function(x, y) x == y)
		G[(p + 1):r, (m + 1):s] <- Gx
	}
	
	# Calculate the model implied covariance matrix
	
	Sigma <- G %*% solve(I - Beta) %*% Gamma %*% Phi %*% t(Gamma) %*% t(solve(I - Beta)) %*% t(G)
	
	ret <- list(beta = beta, gamma = gamma, Phi = Phi, Sigma = Sigma)
	class(ret) <- "matrixpls.bentler-weeks"
	return(ret)
	
	}
	
}
	
