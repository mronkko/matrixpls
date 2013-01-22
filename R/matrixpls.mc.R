#'@title Monte Carlo simulations with MatrixPLS
#'
#'@description
#'Generates samples of data using Bentler-Weeks type structural equation model and
#'estimates these using MatrixPLS and optionally PLSPM.
#
#'@details
#'
#'@param rep Number of replications
#'@param beta beta matrix of the Bentler-Weeks model (regression paths between dependent
#'variables)
#'@param gamma gamma matrix of the Bentler-Weeks model (regression paths between dependent
#'and independent variables)
#'@param Phi Phi matrix of the Bentler-Weeks model (independent variable covariance matrix)
#'@param observed.independents A vector containing indices of observed independent variables
#'@param observed.dependents A vector containing indices of observed dependent variables
#'@param nobs Number of observations to generate
#'@param plspm A boolean determining if the results from PLSPM should be compared with
#'MatrixPLS
#'@param ... All other parameters are passed through to MatrixPLS (and PLSPM if plspm is
#'TRUE)
#'@return Nothing at the moment
#'@author Mikko R\303\266nkk\303\266
#'
#'@references # Bentler, P. M., & Weeks, D. G. (1980). \emph{Linear structural equations
#'with latent variables.} Psychometrika, 45(3), 289\342\200\223308. doi:10.1007/BetaF02293905
#'@export
#'@examples

library(MASS)

matrixpls.mc <- function(rep = 1, beta = NULL, gamma = NULL, Phi = NULL, observed.independents = c(), observed.dependents = c(), nobs = NULL, plspm = FALSE, ...) {
    
    # Calculate and print Sigma based on the Bentler-Weeks model
    
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
    # Calculate the model implied covariance matrix and print it
    
    Sigma <- G %*% solve(I - Beta) %*% Gamma %*% Phi %*% t(Gamma) %*% t(solve(I - Beta)) %*% t(G)
    
    times.matrixpls <- rep(0,5)
    times.plspm <- rep(0,5)
    
    for (replication in 1:rep) {
        data <- mvrnorm(n = nobs, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
        times.matrixpls <- times.matrixpls + system.time(results.matrixpls <- matrixpls(data, ...))[1:5]
        if (plspm) {
            times.plspm <- times.plspm + system.time(results.plspm <- plspm(data, ...))[1:5]
            
            # Compare the results and stop if they are too different
            for(index in 1:length(results.plspm$inner.mod)){
	            if (max(abs(results.plspm$inner.mod[[index]]$value-results.matrixpls$inner.mod[[index]]$value)) > 0.001) {
    	            print("Results from MatrixPLS and PLSPM do not match")
        	        print("PLSPM")
            	    print(results.plspm$inner.mod)
                	print("MatrixPLS")
	                print(results.matrixpls$inner.mod)
    	            stop(paste("Validation failed for replication", replication))
            }
			} 
        }
    }
    
    print("Time elapsed: MatrixPLS")
    print(times.matrixpls)
    if (plspm) {
        print("Time elapsed: PLSPM")
        print(times.plspm)
    }
}
 
