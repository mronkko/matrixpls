#'@title Basic results for Partial Least Squares Path Modeling
#'
#'@description
#'Estimate path models with latent variables by partial least squares approach
#'without providing the full list of results as \code{matrixpls}. This might be helpful
#'when doing simulations, intensive computations, or when you don't want
#'the whole enchilada.
#'
#'@details
#'\code{matrixpls.fit} performs the basic PLS algorithm and provides
#'limited results (e.g. outer weights, LVs scores, path coefficients, R2, and
#'loadings). \cr
#'
#'The argument \code{inner_matrix} is a matrix of zeros and ones that indicates
#'the structural relationships between latent variables. This must be a lower
#'triangular matrix. \code{inner_matrix} will contain a 1 when column \code{j}
#'affects row \code{i}, 0 otherwise. \cr
#'
#'@param Data A numeric matrix or data frame containing the manifest variables.
#'@param inner_matrix A square (lower triangular) boolean matrix representing the
#'inner model (i.e. the path relationships betwenn latent variables).
#'@param outer_list List of vectors with column indices from \code{x} indicating 
#'the sets of manifest variables asociated to the latent variables
#'(i.e. which manifest variables correspond to the latent variables).
#'Length of \code{outer_list} must be equal to the number of rows of \code{inner_matrix}.
#'@param modes A character vector indicating the type of measurement for each
#'latent variable. \code{"A"} for reflective measurement or \code{"B"} for
#'formative measurement (\code{NULL} by default). The length of \code{modes}
#'must be equal to the length of \code{outer_list}).
#'@param scheme A string of characters indicating the type of inner weighting
#'scheme. Possible values are \code{"centroid"}, \code{"factor"}, or
#'\code{"path"}.
#'@param scaled A logical value indicating whether scaling data is performed
#'When (\code{TRUE} data is scaled to standardized values (mean=0 and variance=1)
#'The variance is calculated dividing by \code{N} instead of \code{N-1}).
#'@param tol Decimal value indicating the tolerance criterion for the
#'iterations (\code{tol=0.00001}). Can be specified between 0 and 0.001.
#'@param iter An integer indicating the maximum number of iterations
#'(\code{iter=100} by default). The minimum value of \code{iter} is 100.
#'@return An object of class \code{"matrixpls"}. 
#'@return \item{outer.mod}{Results of the outer (measurement) model. Includes:
#'outer weights, standardized loadings, communalities, and redundancies}
#'@return \item{inner.mod}{Results of the inner (structural) model. Includes: path
#'coefficients and R-squared for each endogenous latent variable}
#'@return \item{latents}{Matrix of standardized latent variables (variance=1
#'calculated divided by \code{N}) obtained from centered data (mean=0)}
#'@return \item{scores}{Matrix of latent variables used to estimate the inner
#'model. If \code{scaled=FALSE} then \code{scores} are latent variables
#'calculated with the original data (non-stardardized). If \code{scaled=TRUE}
#'then \code{scores} and \code{latents} have the same values}
#'@return \item{out.weights}{Vector of outer weights}
#'@return \item{loadings}{Vector of standardized loadings (i.e. correlations with
#'LVs)}
#'@return \item{path.coefs}{Matrix of path coefficients (this matrix has a similar
#'form as \code{inner_matrix})}
#'@return \item{r.sqr}{Vector of R-squared coefficients}
#'@author Gaston Sanchez
#'
#'@references Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro C.
#'(2005) PLS path modeling. \emph{Computational Statistics & Data Analysis},
#'\bold{48}, pp. 159-205.
#'
#'Lohmoller J.-B. (1989) \emph{Latent variables path modelin with partial
#'least squares.} Heidelberg: Physica-Verlag.
#'
#'Wold H. (1985) Partial Least Squares. In: Kotz, S., Johnson, N.L. (Eds.),
#'\emph{Encyclopedia of Statistical Sciences}, Vol. 6. Wiley, New York, pp.
#'581-591.
#'
#'Wold H. (1982) Soft modeling: the basic design and some extensions. In: K.G.
#'Joreskog & H. Wold (Eds.), \emph{Systems under indirect observations:
#'Causality, structure, prediction}, Part 2, pp. 1-54. Amsterdam: Holland.
#'@export
#'@examples
#'
#'  \dontrun{
#'  ## typical example of PLS-PM in customer satisfaction analysis
#'  ## model with six LVs and reflective indicators
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
#'  sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'
#'  # outer model list
#'  sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'
#'  # vector of reflective modes
#'  sat_mod = rep("A", 6)
#'
#'  # apply matrixpls.fit
#'  satpls = matrixpls.fit(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE)
#'  
#'  # summary of results
#'  summary(satpls)
#'
#'  # default plot (inner model)
#'  plot(satpls)
#'  }
#'
matrixpls.fit <-
function(Data, inner_matrix, outer_list, modes = NULL, scheme="centroid", 
         scaled = TRUE, tol = 0.00001, iter = 100)
{
	
	# =======================================================
    # inputs setting
    # =======================================================  

	# Names of the composites
	lvs.names <- colnames(inner_matrix)
	mv.names  <- colnames(Data)
	
	# A binary vector showing which composites are endogenous
	endo <- rowSums(inner_matrix)>0
	
    # The estimation uses the following matrices
    #
    # inner_matrix: a p by p matrix specifying the inner model
    # outer_matrix: a n by p matrix specifying the outer model
	#
    # S: n by n indicator covariance matrix
	#
    # W: n by p matrix of outer weights
    # C: is a p by p covariance matrix of composites
    # E: is a p by p matrix of inner weights
    # 
    # where
    #
	# n: number of indicators
    # p: number of composites

	n<-nrow(Data)
	p<-nrow(inner_matrix)

	outer_matrix <- sapply(outer_list,function(x) match(1:n,x,nomatch=0)>0)
	
	# Standardize the data
	if(scaled) Data <-cov2cor(Data)

    # =======================================================
    # Stage 1: Iterative procedure
    # =======================================================  

	for(iteration in 1:max(1,iter)){
	
		#
		# Outer estimation
		#
		
		if(iteration == 1 ){

			# If this is the first iteration, start with equal weights instead of doing
			# outer estimation
			W_unscaled <- outer_matrix
		}
		else{
			
			# We have already performed inner estimation, do outer estimation of weights
			# Element-wise multiplication by outer_matrix chooses the relevant
			# correlations
			IC <- S %*% W %*% E
			W_unscaled <- IC * outer_matrix
		}

		
		# Calculate the variances of the unscaled composites
		
		var_C_unscaled <- diag(t(W_unscaled) %*% S %*% W_unscaled)
		
		# Create the weights that are used to form the standardized composites after
		# the outer estimation by scaling the unscaled weights
		
		W_new <- W_unscaled %*% diag(x = 1/sqrt(var_C_unscaled))
	
	
		converged <- ifelse(iteration > 1, max(abs(W_new-W)) < tol, 0)
		
		W <- W_new
		
		# Create the composite covariance matrix
	
	 	C <- t(W) %*% S %*% W
	 		
		# Are the weights converged?
		
		if(converged) break;
		
		#
		# Inner estimation
		#
		
		# Calculate unscaled inner weights using the centroid weighting scheme
		
		E_unscaled <- sign(C * (inner_matrix + t(inner_matrix)))
		
		# Calculate the variances of the unscaled composites. 
		
		var_C_unscaled <- diag(E_unscaled * C * E_unscaled)
		
		# Create the weights that are used to form the standardized composites after the
		# inner estimation
		
		E <- E_unscaled %*% diag(x = 1/var_C_unscaled)
		
		# In rare cases, we have underidentified regressions that the algorithm cannot handle
		if(max(is.nan(E))) stop("An inner model regression was empirically underidentified.")
	}	
	
	# Check if we got a convergence
	
	if (!converged) return(NULL)

    # =======================================================
    # Stage 2: Path coefficients
    # =======================================================  

    # A p x p lower diagonal matrix of path estimates
    Path <- inner_matrix
    R2 <- rep(0,p)
    
    # Loop over the endogenous variables and do the regressions

    for (aux in which(endo)){
    	indeps <- inner_matrix[aux,]
    	coefs <- solve(C[indeps,indeps],C[indeps,aux])
    	Path[indeps,aux] <- coefs
    	R2[aux] <- colSums(coefs * C[indeps,aux])
    }

    # Composite loadings and weights are available already from the iterative algorithm
    

    # =======================================================
    # Results
    # =======================================================  

    res = list(out.weights = W, 
               loadings = W_unscaled, 
               path.coefs = Path,
               r.sqr = R2,
               effects = C
               outer.cor = IC)
               
    class(res) = c("matrixpls.fit")
    return(res)
}

