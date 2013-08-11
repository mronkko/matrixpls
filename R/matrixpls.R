#'@title Basic results for Partial Least Squares Path Modeling as a vector
#'
#'@description
#'Estimates a PLS model as described by Wold (1982) and provides a minimum set of estimates as a named vector.
#
#'@details
#'\code{matrixpls} performs the basic PLS algorithm and provides
#'the results in a vector without names. The function is designed to be as efficient
#'as possible. This function does not do any input validation and typically should not
#'be called directly\cr
#'
#'The argument \code{inner_matrix} is a matrix of zeros and ones that indicates
#'the structural relationships between composites. This must be a lower
#'triangular matrix. \code{inner_matrix} will contain a 1 when column \code{j}
#'affects row \code{i}, 0 otherwise. \cr
#'
#'@param S Covariance matrix that is used to estimate the PLS model.
#'@param inner_matrix A square (lower triangular) boolean matrix representing the
#'inner model (i.e. the path relationships between the composite variables).
#'@param outer_list List of vectors with column indices from \code{x} indicating
#'the sets of manifest variables associated to the composites
#'Length of \code{outer_list} must be equal to the number of rows of \code{inner_matrix}.
#'@param modes A character vector indicating the type of estimation for each
#'composite, Either \code{'A'} or \code{'B'}. The length of \code{modes}
#'must be equal to the length of \code{outer_list}).
#'@param scheme A string of characters indicating the type of inner weighting
#'scheme. Possible values are \code{'centroid'}, \code{'factor'}, or
#'\code{'path'}.
#'@param tol Decimal value indicating the tolerance criterion for the
#'iterations (\code{tol=0.00001}). Must be a positive number.
#'@param iter An integer indicating the maximum number of iterations
#'(\code{iter=100} by default).'
#'@param unbiasedLoadings PLS estimates factor loadings with composite loadings severely biased estimates of factor loadings. Setting this option to TRUE will use maximum likelihood common factor analysis instead of the component loadings (default: FALSE)
#'@param unbiasedCompositeReliability The Composite Reliability statistic that is typically presented with PLS results assumes equal weights for the indicators. This is never correct in a PLS analysis resulting in a biased estimate. Setting this option to TRUE will use a formulate that takes indicator weighting into account (default: FALSE)
#'@param disattenuate Setting this option to TRUE will dissattenuate the composite covariance matrix with the estimated composite reliabilities before estimating the inner model regressions. (default: FALSE)
#'@return A vector containing the PLS results in the following order: weights, loading
#'estimates, path estimates, lower triangular of the composite correlation matrix.
#'@author Mikko R\303\266nkk\303\266
#'
#'@references Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. In K. G. J\303\266reskog & S. Wold (Eds.), \emph{Systems under indirect observation\342\200\257: causality, structure, prediction} (pp. 1\342\200\22354). Amsterdam: North-Holland.
#'
#'Lohm\303\266ller, J. B. (1989). \emp{composite path modeling with partial least squares}. Physica-Verlag.
#'
#'Aguirre-Urreta, M. I., Marakas, G. M., & Ellis, M. E. (in press). Measurement of Composite Reliability in Research Using Partial Least Squares: Some Issues and an Alternative Approach. \emp{The DATA BASE for Advances in Information Systems.}
#'
#'Dijkstra, T. K., & Henseler, J. (2012). Consistent and Asymptotically Normal PLS Estimators for Linear Structural Equations. RU Groningen working paper. Retrieved from http://www.rug.nl/staff/t.k.dijkstra/Dijkstra-Henseler-PLSc-linear.pdf
#'

library(psych)

matrixpls <- function(S, inner_matrix, outer_list, modes, scheme = "centroid", tol = 1e-05, iter = 100, unbiasedLoadings = FALSE, unbiasedCompositeReliability = FALSE, disattenuate = FALSE) {
    
    # ==================================================================================== inputs setting ====================================================================================
    
    # Names of the composites
    lvs.names <- colnames(inner_matrix)
    mv.names <- colnames(S)
    
    # A binary vector showing which composites are endogenous
    endo <- rowSums(inner_matrix) > 0
    
    # The estimation uses the following matrices
    # 
    # inner_matrix: a p by p matrix specifying the inner model outer_matrix: a n by p matrix specifying the outer model
    # 
    # S: n by n indicator covariance matrix
    # 
    # W: n by p matrix of outer weights C: is a p by p covariance matrix of composites E: is a p by p matrix of inner weights
    # 
    # where
    # 
    # n: number of indicators p: number of composites
    
    n <- nrow(S)
    p <- nrow(inner_matrix)
    
    outer_matrix <- sapply(outer_list, function(x) match(1:n, x, nomatch = 0) > 0)
    
    S <- cov2cor(S)
    
    # ==================================================================================== Stage 1: Iterative procedure ====================================================================================
    
    for (iteration in 1:max(1, iter)) {
        
        # Outer estimation
        
        if (iteration == 1) {
            
            # If this is the first iteration, start with equal weights instead of doing outer estimation
            W_unscaled <- outer_matrix
        } else {
            
            # We have already performed inner estimation, do outer estimation of weights Element-wise multiplication by outer_matrix chooses the relevant correlations
            IC <- S %*% W %*% E
            W_unscaled <- IC * outer_matrix
        }
        
        
        # Calculate the variances of the unscaled composites
        
        var_C_unscaled <- diag(t(W_unscaled) %*% S %*% W_unscaled)
        
        # Create the weights that are used to form the standardized composites after the outer estimation by scaling the unscaled weights
        
        W_new <- W_unscaled %*% diag(x = 1/sqrt(var_C_unscaled))
        
        
        converged <- ifelse(iteration > 1, sum((abs(W) - abs(W_new))^2) < tol, 0)
        
        W <- W_new
        
        # Create the composite covariance matrix
        
        C <- t(W) %*% S %*% W
        
        # Are the weights converged?
        
        if (converged) 
            break
        
        # Inner estimation
        
        # Calculate unscaled inner weights using the centroid weighting scheme
        
        E_unscaled <- sign(C * (inner_matrix + t(inner_matrix)))
        
        # Calculate the variances of the unscaled composites.  TODO: Optimize this as a single matrix operation
        var_C_unscaled <- rep(1, ncol(C))
        
        for (x in 1:ncol(C)) {
            var_C_unscaled[x] = sum(E_unscaled[x, ] %o% E_unscaled[x, ] * C)
        }
        
        # Create the weights that are used to form the standardized composites after the inner estimation
        
        E <- E_unscaled %*% diag(x = 1/var_C_unscaled)
        
        # In rare cases, we have underidentified regressions that the algorithm cannot handle
        if (max(is.nan(E))) {
            for (matName in c("S", "C", "W", "W_unscaled", "E", "E_unscaled", "var_C_unscaled")) {
                print(matName)
                print(get(matName))
            }
            print(paste("Iteration:", iteration, max(is.nan(E))))
            stop("An inner model regression was empirically underidentified.")
        }
    }
    
    # Check if we got a convergence
    
    if (!converged) 
        return(NULL)
    
    # ==================================================================================== Stage 2: Path coefficients ====================================================================================
    
    
    if (unbiasedLoadings) {
        load.est <- matrix(0, nrow(S), ncol(inner))
        # Loop over the indicator blocks and perform EFAs

        for (col in 1:length(outer)) {
			load.est[outer[[col]], col] <- fa(S[outer[[col]], outer[[col]]])$loadings[, 1]
			
		}
    } else {
        load.est = W_unscaled
    }
    
    if (disattenuate) {
        
        if (unbiasedCompositeReliability) {
            # Calculate composite reliabilities
            # 
            # Aguirre-Urreta, M. I., Marakas, G. M., & Ellis, M. E. (in press). Measurement of Composite Reliability in Research Using Partial Least Squares: Some Issues and an Alternative Approach. The DATA BASE
            # for Advances in Information Systems.
            
            reliabilities <- colSums(W * load.est)^2
            disattenuationMatrix <- sqrt(reliabilities %o% reliabilities)
            diag(disattenuationMatrix) <- C <- C * disattenuationMatrix
        } else {
            stop("Disattennuation has been implemented only for unbiased composite reliability")
        }
    }
    
    # A p x p lower diagonal matrix of path estimates
    Path <- inner_matrix
    R2 <- rep(0, p)
    
    # Loop over the endogenous variables and do the regressions
    
    for (aux in which(endo)) {
        indeps <- which(inner_matrix[aux, ] == 1)
        coefs <- solve(C[indeps, indeps], C[indeps, aux])
        Path[aux, indeps] <- coefs
    }
    
    
    # ==================================================================================== Results ====================================================================================
    
    
    return(c(W[which(outer_matrix==1)], load.est[which(outer_matrix==1)], Path[which(inner_matrix==1)], C[lower.tri(inner_matrix)]))
    
}
 
