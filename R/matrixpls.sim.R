#'@title Monte Carlo simulations with MatrixPLS
#'
#'@description
#'Generates samples of data using a structural equation model and
#'estimates these using MatrixPLS.
#
#'@details
#
#'@param R Number of replications
#'@param inner_matrix
#'@param outer_list
#'@param modes
#'@param nobs Number of observations to generate
#'@param matrixpls.arguments
#'@param generate.arguments
#'@param boot.R = 0
#'@param boot.arguments
#'@param FUN
#'@param printProgress 

#'@return An object of class \code{'matrixpls.mc'} containing.
#'@return \item{estimates}{A matrix containing the PLS results as columns in the following order: weights, loading
#'estimates, path estimates, lower triangular of the composite correlation matrix}
#'@return \item{boot}{A list of \code{'boot'} objects if bootstrapping was enabled}
#'@return \item{beta}{The beta matrix}
#'@return \item{gamma}{The gamma matrix}
#'@return \item{Phi}{The Phi matrix}
#'@return \item{Sigma}{The population covariance matrix for the data based on the Bentler-Weeks model}
#'@author Mikko R\303\266nkk\303\266
#'
#'@references # Bentler, P. M., & Weeks, D. G. (1980). \emph{Linear structural equations
#'with latent variables.} Psychometrika, 45(3), 289\342\200\223308. doi:10.1007/BetaF02293905
#'
#'@examples
#'
#'#
#'# Creates Monte Carlo study with four latent variables, A-D each with 3 indicators and estimates
#'#
#'# Set up parallel bootstrapping
#'
#'options(boot.parallel='multicore')
#'options(boot.ncpus=detectCores())
#'
#'# Create a Monte Carlo study
#'
#'Phi<-diag(c(1,1, #Variances of A and B
#'\t\t\t1-(.5*(.5+.3*.1)+.1*(.1+.3*.5)), #Error variance of C
#'\t\t\t1-2*.3*(.5+.3*.1), # Error variance of D
#'\t\t\trep(1-.7^2,12))) # Indicator error variances
#'
#'Phi[1,2]<-Phi[2,1]<-.3 #Correlation between A and B
#'
#'colnames(Phi)<-c('A','B','eC','eD',paste('ex',1:12,sep=''))
#'rownames(Phi)<-colnames(Phi)
#'
#'gamma <- matrix(0,14,16)
#'gamma[1,1] <- .5 # C on A
#'gamma[1,2] <- .1 # C on B
#'gamma[2,1] <- .3 # D on A
#'gamma[3:5,1] <- .7 # x1-x3 on A
#'gamma[6:8,2] <- .7 # x4-x6 on B
#'
#'gamma[1:14,3:16] <- diag(14) # errors
#'
#'colnames(gamma)<-colnames(Phi)
#'rownames(gamma)<-c('C','D',paste('x',1:12,sep=''))
#'
#'beta <- matrix(0,14,14)
#'beta[2,1] <- .3 # D on C
#'beta[9:11,1] <- .7 # x7-x9 on C
#'beta[12:14,2] <- .7 # x10-x12 on D
#'
#'rownames(beta)<-rownames(gamma)
#'colnames(beta)<-rownames(beta)
#'
#'inner<-matrix(0,4,4)
#'inner[c(3,4),1]<-1
#'inner[c(3,4),2]<-1
#'inner[4,3]<-1
#'
#'matrixpls.mc(rep=10, beta=beta, gamma=gamma, Phi=Phi, observed.dependents=c(3:14), n=100, plspm = TRUE,
#'\tinner = inner, outer = split(1:12,rep(1:4, each=3)), modes = rep('A',4), boot.val=TRUE, br=500)
#'


library(simsem)
library()

A function that takes a data set and returns a list. The list must
contain at least three objects: a vector of parameter estimates (coef),
a vector of standard error (se), and the convergence status as TRUE or
FALSE (converged). There are five optional objects in the list: a vector
of fit indices (fit), a vector of standardized estimates (std), any
extra output (extra), fraction missing type I (FMI1), and fraction
missing type II (FMI2). Note that the coef, se, std, FMI1, and FMI2 must
be a vector with names. The name of those vectors across different
objects must be the same.

matrixpls.sim <- function(R = 1, inner_matrix, outer_list, modes, nobs = NULL,
	matrixpls.arguments = list(),
	generate.arguments = list(),
	boot.R = 0, boot.arguments=list(),
	FUN = NULL,
	printProgress = FALSE) {
    
    #
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
    
    if (plspm) {
        if (!require(plspm)) 
            stop("PLSPM is required for comparing the performance of PLSPM and MatrixPLS")
        
        times.matrixpls <- rep(0, 5)
        times.plspm <- rep(0, 5)
        
        # List to store the results
        
        resultList <- vector("list", rep)
        
        for (replication in 1:rep) {
            data <- mvrnorm(n = nobs, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
            times.matrixpls <- times.matrixpls + system.time(results.matrixpls <- matrixpls(data, ...))[1:5]
            if (plspm) {
                times.plspm <- times.plspm + system.time(results.plspm <- plspm(data, ...))[1:5]
                
                # Compare the results and stop if they are too different
                for (index in 1:length(results.plspm$inner.mod)) {
                  if (max(abs(results.plspm$inner.mod[[index]]$value - results.matrixpls$inner.mod[[index]]$value)) > 0.001) {
                    print("Results from MatrixPLS and PLSPM do not match")
                    print("PLSPM")
                    print(results.plspm$inner.mod)
                    print("MatrixPLS")
                    print(results.matrixpls$inner.mod)
                    stop(paste("Validation failed for replication", replication))
                  }
                }
            }
            resultList[[replication]] <- results.matrixpls
        }
        
        print("Time elapsed: MatrixPLS")
        print(times.matrixpls)
        if (plspm) {
            print("Time elapsed: PLSPM")
            print(times.plspm)
        }
    } else {
        
        results <- NA
        boot.results <- NA
        
        if (boot.val) {
            # Do bootstrapping
            
            resultList <- mclapply(1:rep, function(replication) {
            		if(printProgress && (replication %% 50 == 0)) print(paste("Starting replication",replication))
	                data <- mvrnorm(n = nobs, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
	                return(matrixpls.boot(data, br, ...))
            })
            
            # loop over the results and combine as a matrix
            
            for (i in 1:length(resultList)) {
            	tryCatch({
	                result <- resultList[[i]]
    	            if (i==1) {
        	          results <- matrix(NA, length(resultList), length(result$t0))
            	      boot.results <- list(rep(NA, length(resultList)))
                	}
	                results[i, ] <- result$t0
    	            boot.results[[i]] <- result
                }, error = function(e){return(NA)})

            }
        } else {
            
            # No bootstrapping
            
            resultList <- mclapply(1:rep, function(replication) {
            	if(printProgress) print(paste("Starting replication",replication))
                data <- mvrnorm(n = nobs, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
                return(matrixpls.fit(cov(data), ...))
            })
            
            # loop over the results and combine as a matrix
            
            for (i in 1:length(resultList)) {
                result <- resultList[[i]]
                if (i==1) {
                  results <- matrix(NA, length(resultList), length(result))
                }
                results[i, ] <- result
            }
        }
        
        
        
        ret <- list(estimates = results, boot = boot.results, beta = beta, gamma = gamma, Phi = Phi, Sigma = Sigma)
        class(ret) <- "matrixpls.mc"
        return(ret)
        
    }
    
}
 
