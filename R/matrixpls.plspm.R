library(plspm)

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
#'@author Mikko Rönkkö
#'@author Gaston Sanchez
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
#'  ## typical example of PLS-PM in customer satisfaction analysis
#'  ## model with six LVs and reflective indicators
#'
#'  # load dataset satisfaction
#'  data(satisfaction)
#'
#'  # inner_matrix model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0) 
#'  LOY = c(1,0,0,0,1,0)
#'  sat_inner_matrix = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'
#'  # outer model list
#'  sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'
#'  # vector of modes (reflective indicators)
#'  sat_mod = rep("A", 6)
#'
#'  # apply plspm
#'  satpls = matrixpls.plspm(satisfaction, sat_inner_matrix, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)
#'  
#'  # summary of results
#'  summary(satpls)
#'
#'  # default plot (inner_matrix model)
#'  plot(satpls)
#'  }
#'

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
			outerEstimators <- list(rep(NA,length(outerEstimators)))
			outerEstimators[! modeA] <- matrixpls.outerEstimator.modeB 
			outerEstimators[modeA] <- matrixpls.outerEstimator.modeA
		}
		
		
		innerEstimator <- c(matrixpls.innerEstimator.centroid,
												matrixpls.innerEstimator.factor,
												matrixpls.innerEstimator.path)[[params$scheme]]
		
		# Convergence is checked with the following function in plspm
		
		convergenceCheckFunction <- function(W,W_new){sum((abs(W_new) - abs(W))^2)} 
		
		#
		# Perform estimation
		# 

		if(scaled) S <- cor(dataToUse)
		else  S <- cov(dataToUse)
				
		
		
		if(params$boot.val){
			
			boot.res <- boot(dataToUse, function(originalData,indices){
				S_boot <- cov(originalData[indices,])
				if(scaled) S_boot <- cov2cor(S) 
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
			
			matrixpls.res <- boot.out$t0
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
		
		# Composite values, fitte values, and intercepts
		
		lvScores <- dataToUse %*% t(W)
		lvScoreFittedValues <- lvScores %*% t(beta)
		intercepts <- apply(lvScores-lvScoreFittedValues,2,mean)
		rownames(lvScores) <- 1:nrow(lvScores)
		
		outer.mod <- sapply(lvNames, function(lvName){
			row <- which(lvNames == lvName)
			weights <- W[row,params$outer[[row]]]
			std.loads <- IC_std[row,params$outer[[row]]]
			communal <- std.loads^2
			redundan <- communal * R2[row]
			
			ret <- cbind(weights, std.loads, communal, redundan)
			rownames(ret) <- gsub(paste(lvName,"=+",sep=""),"",names(weights), fixed=TRUE)

			return(ret)
		}, simplify= FALSE)
		
		inner.mod <- sapply(lvNames[rowSums(nativeModel$inner) > 0], function(lvName){
			row <- which(lvNames == lvName)
			regressors <- beta[row,]!=0
			ret <- data.frame(concept =c ("R2","Intercept",paste("path_",lvNames[regressors],sep="")),
												 value = c(R2[row],intercepts[row],beta[row,regressors]))
			return(ret)
		}, simplify= FALSE)
		
		if(params$boot.val){
			stop("Bootstrapping not implemented")
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
		
		lvScores_std <- apply(lvScores,2,scale)
		rownames(lvScores_std) <- 1:nrow(lvScores_std)
		
		Windices <-W!=0
		out.weights <- W[Windices]
		names(out.weights) <- colnames(W)[col(W)[Windices]]

		# Weight and loading pattern is always identical in plspm
		loadings <- IC_std[Windices]
		names(loadings) <- names(out.weights)
		
		outer.cor <- sapply(lvNames, function(lvName){
			row <- which(lvNames == lvName)
			vars <- params$outer[[row]]	
			return(t(IC_std[,vars]))
		}, simplify= FALSE)
		
		> a$inner.sum
		LV.Type Measure MVs  R.square  Av.Commu  Av.Redun       AVE
		IMAG  Exogen   Rflct   5 0.0000000 0.5822708 0.0000000 0.5822708
		EXPE Endogen   Rflct   5 0.3351926 0.6164204 0.2066196 0.6164204
		QUAL Endogen   Rflct   5 0.7196882 0.6585720 0.4739665 0.6585720
		VAL  Endogen   Rflct   4 0.5900843 0.6644158 0.3920613 0.6644158
		SAT  Endogen   Rflct   4 0.7073207 0.7588907 0.5367791 0.7588907
		LOY  Endogen   Rflct   4 0.5099190 0.6390545 0.3258660 0.6390545
		
		res = list(outer.mod = outer.mod, 
							 inner.mod = inner.mod, 
							 latents = lvScores_std, 
							 scores = lvScores,
							 out.weights = out.weights, 
							 loadings = loadings, 
							 path.coefs = beta, 
							 r.sqr = R2,
							 outer.cor = outer.cor, 
#							 inner.sum = inner.sum, 
#							 effects = effects,
#							 unidim = unidim, 
#							 gof = gof, 
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
																											 plsregressionsWithCovarianceMatrixAndModelPattern,
																											 regressionsWithCovarianceMatrixAndModelPattern,
																											 regressionsWithCovarianceMatrixAndModelPattern))
	
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

#' @title PLS regression for \code{MatrixPLS}
#' 
#' @description
#' Internal function. \code{get_plsr1} is called by \code{MatrixPLS} to do PLS-R1. Ported from \code{plspm}
#' 
#' @param C covariance matrix of the pls composites. The first composite is regressed on the other composites
#' @param nc number of components
#' @param scaled logical indicating whether to scale the data
#' @return A list with pls regression results
#' @keywords internal
#' @export

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

#' @title Check parameters for \code{plspm} and \code{plspm.fit}
#' 
#' @description
#' Internal function. \code{get_params} is called by \code{plspm}.
#'
#' @param x numeric matrix or data frame containing the manifest variables.
#' @param inner square (lower triangular) boolean matrix for inner model.
#' @param outer List of vectors with column indices from \code{x} indicating 
#' the sets of manifest variables asociated to the latent variables.
#' @param modes character vector indicating the type of measurement.
#' @param scheme string indicating the type of inner weighting scheme.
#' @param scaled logical indicating whether scaling data is performed.
#' @param boot.val logical indicating whether bootstrap validation is performed.
#' @param br integer indicating the number bootstrap resamples.
#' @param plsr logical indicating whether to use pls regression for path coefs.
#' @param tol decimal value indicating the tolerance criterion for covergence.
#' @param iter integer indicating the maximum number of iterations.
#' @param dataset logical indicating whether the data matrix should be retrieved.
#' @return list of validated parameters for \code{plspm} and \code{plspm.fit}.
#' @keywords internal
#' @export
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

 
