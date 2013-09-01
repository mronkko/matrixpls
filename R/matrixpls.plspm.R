library(plspm)
library(psych)
library(boot)
library(parallel)

# =========== Main functions ===========

#'@title A PLS-PM compatibility wrapper for MatrixPLS
#'
#'@description
#'Estimate path models with latent variables by partial least squares approach
#'
#'@details
#'The function \code{matrixpls.plspm} estimates a path model by partial least squares
#'approach providing the full set of results. \cr
#'
#'The argument \code{inner_matrix} is a matrix of zeros and ones that indicates
#'the structural relationships between latent variables. This must be a lower
#'triangular matrix. \code{inner_matrix} will contain a 1 when column \code{j}
#'affects row \code{i}, 0 otherwise. \cr
#'
#'@param Data A numeric matrix or data frame containing the manifest variables.
#'@param inner_matrix A square (lower triangular) boolean matrix representing the
#'inner_matrix model (i.e. the path relationships betwenn latent variables).
#'@param outer_list List of vectors with column indices from \code{x} indicating 
#'the sets of manifest variables asociated to the latent variables
#'(i.e. which manifest variables correspond to the latent variables).
#'Length of \code{outer_list} must be equal to the number of rows of \code{inner_matrix}.
#'@param modes A character vector indicating the type of measurement for each
#'latent variable. \code{"A"} for reflective measurement or \code{"B"} for
#'formative measurement (\code{NULL} by default). The length of \code{modes}
#'must be equal to the length of \code{outer_list}).
#'@param scheme A string of characters indicating the type of inner_matrix weighting
#'scheme. Possible values are \code{"centroid"}, \code{"factor"}, or
#'\code{"path"}.
#'@param scaled A logical value indicating whether scaling data is performed
#'When (\code{TRUE} data is scaled to standardized values (mean=0 and variance=1)
#'The variance is calculated dividing by \code{N} instead of \code{N-1}).
#'@param tol Decimal value indicating the tolerance criterion for the
#'iterations (\code{tol=0.00001}). Can be specified between 0 and 0.001.
#'@param iter An integer indicating the maximum number of iterations
#'(\code{iter=100} by default). The minimum value of \code{iter} is 100.
#'@param boot.val A logical value indicating whether bootstrap validation is
#'performed (\code{FALSE} by default). 
#'@param br An integer indicating the number bootstrap resamples. Used only
#'when \code{boot.val=TRUE}. When \code{boot.val=TRUE}, the default number of 
#'re-samples is 100, but it can be specified in a range from 100 to 1000.
#'@param plsr A logical value indicating whether pls regression is applied
#'to calculate path coefficients (\code{FALSE} by default).
#'@param dataset A logical value indicating whether the data matrix should be
#'retrieved (\code{TRUE} by default).
#'@return An object of class \code{"plspm"}. 
#'@return \item{outer.mod}{Results of the outer (measurement) model. Includes:
#'outer weights, standardized loadings, communalities, and redundancies}
#'@return \item{inner_matrix.mod}{Results of the inner_matrix (structural) model. Includes: path
#'coefficients and R-squared for each endogenous latent variable}
#'@return \item{latents}{Matrix of standardized latent variables (variance=1
#'calculated divided by \code{N}) obtained from centered data (mean=0)}
#'@return \item{scores}{Matrix of latent variables used to estimate the inner_matrix
#'model. If \code{scaled=FALSE} then \code{scores} are latent variables
#'calculated with the original data (non-stardardized). If \code{scaled=TRUE}
#'then \code{scores} and \code{latents} have the same values}
#'@return \item{out.weights}{Vector of outer weights}
#'@return \item{loadings}{Vector of standardized loadings (i.e. correlations with
#'LVs)}
#'@return \item{path.coefs}{Matrix of path coefficients (this matrix has a similar
#'form as \code{inner_matrix})}
#'@return \item{r.sqr}{Vector of R-squared coefficients}
#'@return \item{outer.cor}{Correlations between the latent variables and the
#'manifest variables (also called crossloadings)}
#'@return \item{inner_matrix.sum}{Summarized results by latent variable of the inner_matrix
#'model. Includes: type of LV, type of measurement, number of indicators,
#'R-squared, average communality, average redundancy, and average variance
#'extracted}
#'@return \item{effects}{Path effects of the structural relationships. Includes:
#'direct, indirect, and total effects}
#'@return \item{unidim}{Results for checking the unidimensionality of blocks
#'(These results are only meaningful for reflective blocks)}
#'@return \item{gof}{Goodness-of-Fit index}
#'@return \item{data}{Data matrix containing the manifest variables used in the
#'model. Only when \code{dataset=TRUE}}
#'@return \item{boot}{List of bootstrapping results; only available when argument
#'\code{boot.val=TRUE}}
#'
#'@references Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro C.
#'(2005) PLS path modeling. \emph{Computational Statistics & Data Analysis},
#'\bold{48}, pp. 159-205.
#'
#'Tenenhaus M., Pages J. (2001) Multiple factor analysis combined with
#'PLS path modelling. Application to the analysis of relationships between
#'physicochemical variables, sensory profiles and hedonic judgements.
#'\emph{Chemometrics and Intelligent Laboratory Systems}, \bold{58}, pp.
#'261-273.
#'
#'Tenenhaus M., Hanafi M. (2010) A bridge between PLS path modeling and
#'multi-block data analysis. \emph{Handbook on Partial Least Squares (PLS):
#'Concepts, methods, and applications.} Springer.
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
#'@seealso \code{\link[plspm]{plspm}}, \code{\link[plspm]{plspm.fit}}, \code{\link[plspm]{plot.plspm}}
#'@export
#'@examples
#'
#'  \dontrun{
#'  library(plspm)
#'  
#'  # Run the example from plspm package
#'  
#'  # load dataset satisfaction
#'  data(satisfaction)
#'  # inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0)
#'  LOY = c(1,0,0,0,1,0)
#'  sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'  # outer model list
#'  sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'  # vector of modes (reflective indicators)
#'  sat_mod = rep("A", 6)
#'  
#'  # apply plspm
#'  plspm.res <- plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=FALSE)
#'  
#'  # apply MatrixPLS
#'  matrixpls.res <- matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=FALSE)
#'  
#'  # plspm scales latent variable scores based on estimated population sd whereas matrixpls scales
#'  # standardizes them in the sample. This creates small differences in some of the results.
#'  
#'  print(plspm.res)
#'  print(plspm.res)
#'  
#'  # If RUnit is installed check that the results are identical
#'  
#'  if(is.element("RUnit", installed.packages()[,1])){
#'  	library(RUnit)
#'  	checkEquals(plspm.res, plspm.res, tol = 0.0001)
#'  }
#'}

matrixpls.plspm <-
function(Data, inner_matrix, outer_list, modes = NULL, scheme = "centroid", 
         scaled = TRUE, tol = 0.00001, iter = 100, boot.val = FALSE, 
         br = NULL, plsr = FALSE, dataset = TRUE){
    
		# =======================================================
		# checking arguments
		# =======================================================
		params = get_params(x=Data, inner=inner_matrix, outer=outer_list, modes=modes, 
										 scheme=scheme, scaled=scaled, tol=tol, iter=iter,
										 boot.val=boot.val, br=br, plsr=plsr, dataset=dataset)

		
		#
		# Setup names etc
		#
		
		if(is.null(colnames(params$inner))) colnames(params$inner) <- rownames(params$inner)
		
		# 
		# Set up  matrixpls parameters
		#
		
		lvNames <- rownames(params$inner)
		reflective <- matrix(0, max(unlist(params$outer)),nrow(params$inner))
		rownames(reflective) <- colnames(params$x)[1:nrow(reflective)]
		colnames(reflective) <- lvNames
			
		col <- 1		
		
		for(outerSpec in params$outer){
			reflective[outerSpec,col] <- 1
			col <- col + 1 
		}
		
		# Only keep the variables that are used in the model

		varsToUse <- which(rowSums(reflective) > 0)
		reflective <- reflective[varsToUse,]
		dataToUse <- as.matrix(params$x[,varsToUse])
		
		# plspm does not support formative relationships, so we create a matrix of zero
		formative <- t(reflective)
		formative[] <- 0
		
		nativeModel <- list(inner=params$inner, 
												reflective=reflective, 
												formative=formative) 

		weightRelations <- t(reflective)

		if(params$plsr) parameterEstimator <- matrixpls.parameterEstimator.plsregression
		else parameterEstimator <- matrixpls.parameterEstimator.regression

		modeA <- params$modes == "A"
		
		if(max(modeA) == 0) outerEstimators <- matrixpls.outerEstimator.modeB	
		else if(min(modeA) == 1) outerEstimators <- matrixpls.outerEstimator.modeA
		else{
			outerEstimators <- list(rep(NA,length(modeA)))
			outerEstimators[!modeA] <- list(matrixpls.outerEstimator.modeB)
			outerEstimators[modeA] <- list(matrixpls.outerEstimator.modeA)
		}
		
		
		innerEstimator <- c(matrixpls.innerEstimator.centroid,
												matrixpls.innerEstimator.factor,
												matrixpls.innerEstimator.path)[[params$scheme]]
	
		# Convergence is checked with the following function in plspm
		
		convergenceCheckFunction <- function(W,W_new){sum((abs(W_new) - abs(W))^2)} 
		
		#
		# Perform estimation
		# 

		if(scaled){
			S <- cor(dataToUse)
			S_std <- S
		}
		else{
			S <- cov(dataToUse)
			S_std <- cov2cor(S)
		} 
				
		
		
		if(params$boot.val){
			
			#
			# The matrixpls.boot function is not used here because we need to support data scaling
			#
			
			boot.res <- boot(dataToUse, function(originalData,indices){
				
				# plspm reports always standardized loadings even if the data were not standardized
				# We mimic this by standardizing the data
				
				if(scaled){
					S_boot <- cor(originalData[indices,])
				}
				else{
					S_boot <- cov(originalData[indices,])
				}
				
				tryCatch(
					matrixpls(S_boot, model = nativeModel, weightRelations = weightRelations, parameterEstimator = parameterEstimator,
										outerEstimators = outerEstimators, innerEstimator = innerEstimator,
										tol = params$tol, iter = params$iter, convergenceCheckFunction = convergenceCheckFunction,
										validateInput = FALSE), 
					
					error = function(e){
						print(e)
						print(S)
						print(nativeModel)
						print(weightRelations)
					})
			}, params$br)
			
			matrixpls.res <- boot.res$t0
		}
		else{
			matrixpls.res <- matrixpls(S, model = nativeModel, weightRelations = weightRelations, parameterEstimator = parameterEstimator,
																 outerEstimators = outerEstimators, innerEstimator = innerEstimator,
																 tol = params$tol, iter = params$iter, convergenceCheckFunction = convergenceCheckFunction,
																 validateInput = FALSE)
		}
		
		#
		# Compile a result object in the plspm format
		#
		
		C <- attr(matrixpls.res,"C") 
		IC <- attr(matrixpls.res,"IC")
		W <- attr(matrixpls.res,"W") 
		beta <- attr(matrixpls.res,"beta")
		
		IC_std <- IC %*% (diag(1/sqrt(diag(S))))
		colnames(IC_std) <- colnames(IC)
		
		# Inner model R2s
		R2 <- rowSums(beta * C)
		names(R2) <- lvNames

		# PLSPM does LV score scaling differently, so we need a correction
		
		sdv = sqrt(nrow(dataToUse)/(nrow(dataToUse)-1))   
		
		# Composite values, fittes values, and intercepts
		
		lvScores <- dataToUse %*% t(W) * sdv
		lvScoreFittedValues <- lvScores %*% t(beta)
		intercepts <- apply(lvScores-lvScoreFittedValues,2,mean)
		rownames(lvScores) <- 1:nrow(lvScores)
		
		
		Modes <- ifelse(params$modes == "A", "Reflective","Formative")
		
		exogenousLVs <- rowSums(nativeModel$inner) == 0
				
		outer.mod <- sapply(lvNames, function(lvName){
			row <- which(lvNames == lvName)
			weights <- W[row,params$outer[[row]]]*sdv
			std.loads <- IC_std[row,params$outer[[row]]]
			communal <- std.loads^2
			redundan <- communal * R2[row]
			
			ret <- cbind(weights, std.loads, communal, redundan)
			rownames(ret) <- gsub(paste(lvName,"=+",sep=""),"",names(weights), fixed=TRUE)

			return(ret)
		}, simplify= FALSE)
		
		inner.mod <- sapply(lvNames[!exogenousLVs], function(lvName){
			row <- which(lvNames == lvName)
			regressors <- beta[row,]!=0
			ret <- data.frame(concept =c ("R2","Intercept",paste("path_",lvNames[regressors],sep="")),
												 value = c(R2[row],intercepts[row],beta[row,regressors]))
			return(ret)
		}, simplify= FALSE)
		
		pathIndices <- lower.tri(nativeModel$inner)
		
		pathLabels <- paste(colnames(nativeModel$inner)[col(nativeModel$inner)[pathIndices]],"->",
												rownames(nativeModel$inner)[row(nativeModel$inner)[pathIndices]], sep="")
		
		if(params$boot.val){
			
			pathCount <- sum(nativeModel$inner)
			weightCount <- sum(weightRelations)
			bootPathIndices <- 1:pathCount
			bootLoadingIndices <- 1:weightCount + pathCount
		  bootWeightIndices <- bootLoadingIndices + weightCount
			
			bootIndices <- boot.array(boot.res, indices = TRUE)
						
			boot <- list(weights = get_bootDataFrame(boot.res$t0[bootWeightIndices]*sdv, boot.res$t[,bootWeightIndices]*sdv, rownames(S)),
									 loadings = get_bootDataFrame(IC_std[weightRelations == 1], 
									 														 t(mcmapply(function(x){
									 														 	
									 														 	# Recover the data that were used in the bootrap replications and calculate indicator variances
									 														 	sdX <- apply(dataToUse[bootIndices[x,],],2,sd)
									 														 	boot.res$t[x,bootLoadingIndices] / sdX									 														 	
									 														 },1:params$br)),
									 														 rownames(S)),
									 paths = get_bootDataFrame(boot.res$t0[bootPathIndices], boot.res$t[,bootPathIndices], 
									 													pathLabels[nativeModel$inner[lower.tri(nativeModel$inner)]==1]),
									 
									 # Use parallel processing because of the large number of elements
									 
									 rsq = get_bootDataFrame(R2[!exogenousLVs],
									 												t(mcmapply(function(x){
									 	 
									 	 # Recover the data that were used in the bootrap replication
									 												if(scaled)
																					   S <- co(dataToUse[bootIndices[x,],])
									 												else
									 													S <- cov(dataToUse[bootIndices[x,],])
									   
									   # Reconstruct the matrices
									   W <- matrix(0,nrow(weightRelations),ncol(weightRelations))
									   W[weightRelations==1] <- boot.res$t[x,bootWeightIndices]
									   C <- W %*% S %*% t(W)
									   beta <- matrix(0,nrow(nativeModel$inner),ncol(nativeModel$inner))
										 beta[nativeModel$inner == 1] <- boot.res$t[x,bootPathIndices]
									   
									   # Calculate R2 values
									   rowSums(beta[!exogenousLVs,]*C[!exogenousLVs,])
									   	
									 },1:params$br)),
									 												lvNames[!exogenousLVs]),
									 
									 # Again, parellel proceesing
									 
									 total.efs = get_bootDataFrame(effects(matrixpls.res)$Total[pathIndices[-1,]],
									 																t(mcmapply(function(x){
									 	
									 	beta <- matrix(0,nrow(nativeModel$inner),ncol(nativeModel$inner))
									 	beta[nativeModel$inner == 1] <- boot.res$t[x,bootPathIndices]
									 	
									 	effects.matrixpls(beta = beta, innerModel = nativeModel$inner)$Total[pathIndices[-1,]]
									 	
									 },1:params$br)),
									 															pathLabels))
		}
		else{
			boot = FALSE
		}
		
		blocks <- unlist(lapply(params$outer, length))
		names(blocks) <- lvNames
		
		model <- list(IDM=params$inner,
									blocks=blocks,
									scheme=switch(scheme, "centroid"="centroid", "factor"="factor", "path"="path"),
									modes=modes, scaled=params$scaled, 
									boot.val=params$boot.val, plsr=params$plsr, obs=nrow(params$x), br=params$br, 
									tol=params$tol, iter=params$iter,
									n.iter=attr(matrixpls.res,"iterations"), outer=params$outer)
		
		if(params$dataset){
			data <- as.matrix(params$x[,unlist(params$outer)])
		}
		else{
			data <- FALSE
		}
		
		lvScores_std <- apply(lvScores, 2, scale) * sdv
		rownames(lvScores_std) <- 1:nrow(lvScores_std)
				
		Windices <-W!=0
		out.weights <- W[Windices] * sdv
		names(out.weights) <- colnames(W)[col(W)[Windices]]

		# Weight and loading pattern is always identical in plspm
		loadings <- IC_std[Windices]
		names(loadings) <- names(out.weights)
		
		outer.cor <- sapply(lvNames, function(lvName){
			row <- which(lvNames == lvName)
			vars <- params$outer[[row]]	
			return(t(IC_std[,vars]))
		}, simplify= FALSE)
		
		
		inner.sum <- data.frame(LV.Type = ifelse(exogenousLVs,"Exogen","Endogen"),
														Measure = abbreviate(Modes, 5),
														MVs = blocks,
														R.square = R2, 
														Av.Commu = sapply(outer.mod,function(x){mean(x[,"communal"])}), 
														Av.Redun = sapply(outer.mod,function(x){mean(x[,"redundan"])}), 
														AVE = ifelse(params$modes == "A",
																				 sapply(outer.mod,function(x){mean(x[,"std.loads"]^2)}),
																				 0))
		
		effs <- effects(matrixpls.res)
				
		effects <- data.frame(relationships = pathLabels,
														dir.effects = effs$Direct[pathIndices[-1,]],
													ind.effects = effs$Indirect[pathIndices[-1,]], 
													tot.effects = effs$Total[pathIndices[-1,]])

		unidim <- data.frame(Type.measure = Modes,
												 MVs = blocks,
												 C.alpha = ifelse(params$modes == "A",
												 								 sapply(params$outer,function(indices){
												 								 	alpha(S[indices,indices])$total[[2]]
												 								 	}, simplify = TRUE),
												 								 0),
												 DG.rho = ifelse(params$modes == "A",
												 								sapply(params$outer,function(indices){
												 									pc <- principal(S[indices,indices])
												 									std.loads <- pc$loadings
												 									numer.rho <- sum(std.loads)^2
												 									denom.rho <- numer.rho + (length(indices) - sum(std.loads^2))
												 									numer.rho / denom.rho
												 									}, simplify = TRUE),
												 								0),
												 
												 eig.1st = sapply(params$outer,function(indices){
												 	 eigen(S_std[indices,indices])$values[1]
												 }, simplify = TRUE),   
												 eig.2nd = sapply(params$outer,function(indices){
												 	eigen(S_std[indices,indices])$values[2]
												 }, simplify = TRUE))
		
		# Goodness of Fit is square root of product of mean communality and mean R2
		
		gof <- sqrt(mean(IC_std[t(nativeModel$reflective)==1]^2) * 
									mean(R2[!exogenousLVs]))
		
		res = list(outer.mod = outer.mod, 
							 inner.mod = inner.mod, 
							 latents = lvScores_std, 
							 scores = lvScores,
							 out.weights = out.weights, 
							 loadings = loadings, 
							 path.coefs = beta, 
							 r.sqr = R2,
							 outer.cor = outer.cor, 
							 inner.sum = inner.sum, 
							 effects = effects,
							 unidim = unidim, 
							 gof = gof, 
							 boot = boot, 
							 data = data, 
							 model = model)
		
		class(res) = "plspm"

		return(res)

}





# =========== Parameter estimators ===========

#'@title PLS parameter estimation with separatePLS regression analyses
#'
#'@description
#'Estimates the model parameters with weighted composites using separate PLS regression and estimating
#'measurement model and latent variable model separately.
#'
#'
#'@details
#'The implementation of PLS regression is ported from the raw data version included in plspm
#'
#'~~~ Explain how this works ~~~
#'
#'@param S Covariance matrix of the data.
#'
#'@param W Weight matrix, where the indicators are on colums and composites are on the rows.
#'
#'@param model There are four options for this argument: 1. SimSem object created by model
#'command from simsem package, 2. lavaan script, lavaan parameter table, or a list that
#'contains all argument that users use to run lavaan (including cfa, sem, lavaan), 3.
#'MxModel object from the OpenMx package, or 4. A list containing three matrices
#'\code{inner}, \code{reflective}, and \code{formative} defining the free regression paths
#'in the model.
#'
#'@param nc number of components
#'@param scaled logical indicating whether to scale the data
#'
#'@return A named vector of parameter estimates
#'
#'@export

matrixpls.parameterEstimator.plsregression <- function(S, W, model, nc=2, scaled = TRUE){
	
	return(matrixpls.parameterEstimator.internal_generic(S,W,model,
																											 plsregressionsWithCovarianceMatrixAndModelPattern))
	
}

# =========== Utility functions ===========

plsregressionsWithCovarianceMatrixAndModelPattern <- function(S,model){

	assert_is_symmetric_matrix(S)
	
	for(row in 1:nrow(model)){
		
		independents <- which(model[row,]!=0)
		
		if(length(independents)>0){
			vars <- c(row,independents)
			coefs <- get_plsr1(S[vars,vars])
			model[row,independents] <- coefs
		}
	}
	
	return(model)
	
}

#
# Run PLS regression. Ported from PLSPM
#

get_plsr1 <-function(C, nc=NULL, scaled=TRUE)
{
	# ============ checking arguments ============
	p <- ncol(X)
	if (is.null(nc))
		nc <- p
	# ============ setting inputs ==============
	if (scaled) C <- cov2cor(C)
	C.old <- C

	Ph <- matrix(NA, p, nc)# matrix of X-loadings
	Wh <- matrix(NA, p, nc)# matrix of raw-weights
	ch <- rep(NA, nc)# vector of y-loadings
	
	# ============ pls regression algorithm ==============
	
	for (h in 1:nc)
	{
		# Covariancese between the independent and dependent vars
		w.old <- C[2:nrow(C),1]
		w.new <- w.old / sqrt(sum(w.old^2)) # normalization
		
		# Covariances between the component and the variables
		cv.new <- C %*% c(0,w.new)
		
		# Covariances between the independents and the component
		p.new <- cv.new[1]
		
		# Covariance between the dependent and the component
		c.new <- cv.new[2:length(cv.new)]
		
		# Deflation
		C.old <- C.old - cv.new

		Ph[,h] <- p.new
		Wh[,h] <- w.new
		ch[h] <- c.new
		
		
	}
	Ws <- Wh %*% solve(t(Ph)%*%Wh)# modified weights
	Bs <- as.vector(Ws %*% ch) # std beta coeffs    
	Br <- Bs * C[1,1]/diag(C)[2:nrow(C)]   # beta coeffs

	return(Br)
}

#
# Validate parameters. Copied from PLSPM
#

get_params <-
	function(x, inner, outer, modes=NULL, scheme="centroid", scaled=TRUE,
					 boot.val=FALSE, br=NULL, plsr=FALSE, tol=0.00001, iter=100, dataset=TRUE)
	{
		# =======================================================
		# checking arguments
		# =======================================================
		if (!is.matrix(x) && !is.data.frame(x))
			stop("Invalid object 'x'. Must be a numeric matrix or data frame.")
		if (is.null(rownames(x))) 
			rownames(x) <- 1:nrow(x)
		if (is.null(colnames(x))) 
			colnames(x) <- paste("MV", 1:ncol(x), sep="")
		if (!is.matrix(inner))
			stop("Invalid argument 'inner'. Must be a matrix.")
		if (nrow(inner)!=ncol(inner))
			stop("Invalid argument 'inner'. Must be a square matrix.")
		for (j in 1:ncol(inner))
			for (i in 1:nrow(inner)) {
				if (i<=j)
					if (inner[i,j]!=0) 
						stop("argument 'inner' must be a lower triangular matrix")
				if (length(intersect(inner[i,j], c(1,0)))==0)
					stop("elements in 'inner' must be '1' or '0'")
			}
		if (is.null(dimnames(inner)))
			lvs.names <- paste("LV", 1:ncol(inner), sep="")
		if (!is.null(rownames(inner)))
			lvs.names <- rownames(inner)
		if (!is.null(colnames(inner)))
			lvs.names <- colnames(inner)
		if (!is.list(outer))
			stop("Invalid argument 'outer'. Must be a list.")
		if (length(outer) != nrow(inner))
			stop("Number of rows of 'inner' does not coincide with length of 'outer'.")
		if (is.null(modes)) {
			modes <- rep("A",length(outer))
			warning("Argument 'modes' missing. Default reflective 'modes' is used.")
		}
		if (length(outer) != length(modes)) {
			warning("Warning: Invalid length of 'modes'. Default reflective 'modes' is used.")
			modes <- rep("A", length(outer))
		}
		for (i in 1:length(modes))
			if (modes[i]!="A" && modes[i]!="B") modes[i]<-"A"
		if (!is.na(pmatch(scheme, "centroid"))) 
			scheme <- "centroid"
		SCHEMES <- c("centroid", "factor", "path")
		scheme <- pmatch(scheme, SCHEMES)
		if (is.na(scheme)) {
			warning("Warning: Invalid argument 'scheme'. Default 'scheme=centroid' is used.")   
			scheme <- "centroid"
		}
		if (!is.logical(scaled)) {
			warning("Warning: Invalid argument 'scaled'. Default 'scaled=TRUE' is used.")
			scaled <- TRUE
		}
		if (!is.logical(boot.val)) {
			warning("Warning: Invalid argument 'boot.val'. No bootstrap validation is done.")
			boot.val <- FALSE
		}   
		if (boot.val) {
			if (!is.null(br)) {        
				if (mode(br)!="numeric" || length(br)!=1 || (br%%1)!=0 ||
							br<100 || br>1000) {
					warning("Warning: Invalid argument 'br'. Default 'br=100' is used.")   
					br <- 100
				} 
			} else
				br <- 100
		}
		if (!is.logical(plsr)) plsr<-FALSE
		if (mode(tol)!="numeric" || length(tol)!=1 || tol<=0 || tol>0.001) {
			warning("Warning: Invalid argument 'tol'. Default 'tol=0.00001' is used.")   
			tol <- 0.00001
		} 
		if (mode(iter)!="numeric" || length(iter)!=1 || iter<100) {
			warning("Warning: Invalid argument 'iter'. Default 'iter=100' is used.")   
			iter <- 100
		} 
		if (!is.logical(dataset)) 
			dataset <- TRUE
		
		# =======================================================
		# list with verified arguments
		# =======================================================
		res = list(
			x = x,
			inner = inner,
			outer = outer,
			modes = modes,
			scheme = scheme,
			SCHEMES = SCHEMES,
			scaled = scaled,
			boot.val = boot.val,
			br = br,
			plsr = plsr,
			tol = tol,
			iter = iter,
			dataset = dataset)
		res
}

get_bootDataFrame <- function(t0, t, rownames){
	
	ret <- data.frame(Original = t0, 
										Mean.Boot = apply(t, 2, mean), 
										Std.Error = apply(t, 2, sd), 
										perc.025 = apply(t, 2, quantile, 0.025),
										perc.975 = apply(t, 2, quantile, 0.975))
	
	rownames(ret) <- rownames
	
	return(ret)
}
 
