#
# Various post-estimation functions
#

#'@title 	Total, Direct, and Indirect Effects for PLS model
#'
#'@description
#'
#'The \code{matrixpls} method for the standard generic function \code{effects} computes total, direct, 
#'and indirect effects for a PLS model according to the method described in Fox (1980).
#'
#'Adapted from the \code{\link[sem]{effects}} function of the \code{sem} package
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A list with \code{Total}, \code{Direct}, and \code{Indirect} elements.
#'
#'@references
#'
#'Fox, J. (1980) Effect analysis in structural equation models: Extensions and simplified methods of computation. \emph{Sociological Methods and Research}
#'9, 3--28.
#'
#'@family post-estimation functions
#'
#'@importFrom stats effects
#'
#'@method effects matrixpls
#'
#'@S3method effects matrixpls

effects.matrixpls <- function(object = NULL,  ...) {
	
	A <- attr(object,"beta")
	endog <- rowSums(attr(object,"model")$inner)!=0 
	
	I <- diag(endog)
	AA <- - A
	diag(AA) <- 1
	Total <- solve(AA) - I
	Indirect <-  Total - A
	result <- list(Total=Total[endog, ], Direct=A[endog, ], Indirect=Indirect[endog, ])
	class(result) <- "matrixplseffects"
	result
}

#'@S3method print matrixplseffects

print.matrixplseffects <- function(x, ...){
	cat("\n Total Effects (column on row)\n")
	Total <- x$Total
	Direct <- x$Direct
	Indirect <- x$Indirect
	select <- !(apply(Total, 2, function(x) all( x == 0)) & 
								apply(Direct, 2, function(x) all( x == 0)) & 
								apply(Indirect, 2, function(x) all( x == 0)))
	print(Total[, select], ...)
	cat("\n Direct Effects\n")
	print(Direct[, select], ...)
	cat("\n Indirect Effects\n")
	print(Indirect[, select], ...)
	invisible(x)
}

#'@title Residual diagnostics for PLS model
#'
#'@description
#'
#'The \code{matrixpls} method for generic function \code{residuals} computes the residual
#'covariance matrix and various fit indices presented by Lohmöller (1989, ch 2.4).
#''
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A list with three elements: \code{inner}, \code{outer}, and \code{indices} elements
#' containing the residual covariance matrix of regressions of composites on other composites,
#' the residual covariance matrix of indicators on composites, and various fit indices
#' calculated based on the residuals.
#'
#'@references
#'
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.
#'
#'@export
#'
#'@family post-estimation functions
#'
#'@method residuals matrixpls
#'
#'@S3method residuals matrixpls

residuals.matrixpls <- function(object, ...) {
	
	library(psych)
	library(MASS)
	
	# Equation numbers in parenthesis refer to equation number in Lohmoller 1989
	
	nativeModel <- attr(object,"model")
	
	W <- attr(object,"W")
	S <- attr(object,"S")

	# Number of reflective indicators
	reflectiveIndicators<- rowSums(nativeModel$reflective)>0
	k <- sum(reflectiveIndicators)
	
	# Number of endog LVs
	endog <- rowSums(nativeModel$inner)>0
	h <- sum(endog)
		
	
	# Factor loading matrix
	P <- nativeModel$reflective
	P[P==1] <- object[grepl("=~", names(object), fixed=TRUE)]
	
	
	B <- attr(object,"beta")
	
	# Lohmoller is not clear whether R should be based on the estimated betas or calculated scores
	# The scores are used here because this results in less complex code
	
	R <- W %*% S %*% t(W)
	R_star <- (B %*% R %*% t(B))[endog,endog] # e. 2.99
	
	# Model implied indicator covariances
	H <- P %*% R %*% t(P) # eq 2.96
	
	I <- diag(S)
	H2 <- (I * H)  %*% ginv(I * S) # eq 2.9	
	
	F <- P %*% B %*% R %*% t(B) %*% t(P) # eq 2.104
	
	F2 <- (I * F) %*% ginv(I * S) # eq 2.105
	
	R2 <- R2(object)
	
	# C in Lohmoller 1989 is different from the matrixpls C
	# residual covariance matrix of indicators
	
	C <- S-H
	C_std <- diag(1/diag(S)) %*% C

	Q <- (W %*% S %*% t(W))[endog,endog] - R_star
	

	indices <- list(Communality = tr(H2)/k, # eq 2.109
									Redundancy = tr(F2)/k,  # eq 2.110
									SMC = sum(R2)/h,         # eq 2.111
									"RMS outer residual covariance" = sqrt(sum(C[lower.tri(C)]^2)/length(C[lower.tri(C)])), # eq 2.118
									"RMS inner residual covariance" = sqrt(sum(Q[lower.tri(Q)]^2)/length(Q[lower.tri(Q)])), # eq 2.118
									SRMR = sqrt(sum(C_std[lower.tri(C_std)]^2)/length(C[lower.tri(C_std, diag=TRUE)])))
	
	result<- list(inner = Q, outer = C, indices = indices)
	
	class(result) <- "matrixplsresiduals"
	result
}

#'@S3method print matrixplsresiduals

print.matrixplsresiduals <- function(x, ...){
	cat("\n Inner model (composite) residual covariance matrix\n")
	print(x$inner, ...)
	cat("\n Outer model (indicator) residual covariance matrix\n")
	print(x$outer, ...)
	cat("\n Residual-based fit indices\n")
	print(data.frame(Value = unlist(x$indices)), ...)
}

#'@title R2	for PLS model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{R2} computes the squared multiple correlation (R2)
#'for composites predicted by other composites in the model.
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A named numberic vector containing the R2 values.
#'
#'@family post-estimation functions
#'
#'@export
#'

R2 <- function(object, ...){
	UseMethod("R2")
}

#'@S3method R2 matrixpls

R2.matrixpls <- function(object, ...){
	
	R2 <- rowSums(attr(object,"beta") * attr(object,"C"))
	names(R2) <- colnames(attr(object,"model")$inner)
	class(R2) <- "matrixplsr2"
	R2
}

#'@S3method print matrixplsr2

print.matrixplsr2 <- function(x, ...){
	cat("\n Inner model squared multiple correlations (R2)\n")
	print.table(x, ...)
}

#'@title Goodness of Fit indices for PLS model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{GoF} computes the Goodness of Fit index for PLS model.
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A list with \code{Total}, \code{Direct}, and \code{Indirect} elements.
#'
#'@references
#'
#'Henseler, J., & Sarstedt, M. (2013). Goodness-of-fit indices for partial least squares path modeling. \emph{Computational Statistics}, 28(2), 565–580. doi:10.1007/s00180-012-0317-1
#'
#'
#'@family post-estimation functions
#'
#'@export
#'

GoF <- function(object, ...){
	UseMethod("GoF")
}

#'@S3method GoF matrixpls

GoF.matrixpls <- function(object, ...) {
	
	nativeModel <- attr(object,"model")
	IC <- attr(object,"IC")
	S <- attr(object,"S")

	exogenousLVs <- rowSums(nativeModel$inner) == 0
	
	
	IC_std <- IC %*% (diag(1/sqrt(diag(S))))
	
	result <- sqrt(mean(IC_std[t(nativeModel$reflective)==1]^2) * 
			 	mean(R2(object)[!exogenousLVs]))
	class(result) <- "matrixplsgof"
	result
}

#'@S3method print matrixplsgof

print.matrixplsgof <- function(x, digits=getOption("digits"), ...){
	cat("\n Absolute goodness of fit:", round(x, digits = digits))
	cat("\n")
}

#'@title Factor loadings matrix from PLS model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{loadings} computes standardized factor loading matrix
#'from the model results.
#'
#'@param x An object of class \code{matrixpls}
#'
#'@param ... All other arguments are ignored.
#'
#'@return A matrix of factor loadings.
#'
#'@family post-estimation functions
#'
#'@importFrom stats loadings
#'
#'@export
#
#'@S3method loadings matrixpls

loadings.matrixpls <- function(x, ...) {
	nativeModel <- attr(x,"model")
	IC <- attr(x,"IC")
	
	#Standardize
	S <- attr(x,"S")
	IC_std <- IC %*% (diag(1/sqrt(diag(S))))
	
	res <- nativeModel$reflective
	res[res==1] <- t(IC_std)[res==1]
	res
}

#'@title Composite Reliability indices for PLS model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{CR} computes Composite Reliability 
#'indices for the model using the formula presented by Fornell and Larcker (1981).
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A named numeric vector containing the Composite Reliability indices.
#'
#'@references
#'
#'Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of marketing research}, 18(1), 39–50.
#'
#'@family post-estimation functions
#'
#'@export
#'

CR <- function(object, ...){
	UseMethod("CR")
}

#'@S3method CR matrixpls

CR.matrixpls <- function(object, ...) {

	loadings <- loadings.matrixpls(object)
	reflectiveModel <- attr(object, "model")$reflective
	
	result <- unlist(lapply(1:ncol(loadings), function(col){	
		useLoadings <- loadings[,col][reflectiveModel[,col]==1]
		
		crvalue <- sum(useLoadings)^2/
			(sum(useLoadings)^2 + sum(1-useLoadings^2))
		crvalue
	}))
	
	class(result) <- "matrixplscr"
	names(result) <- colnames(loadings)
	result
}

#'@S3method print matrixplscr

print.matrixplscr <- function(x, ...){
	cat("\n Composite Reliability indices\n")
	print.table(x, ...)
}


#'@title Average Variance Extracted indices for PLS model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{AVE} computes Average Variance Extracted 
#'indices for the model using the formula presented by Fornell and Larcker (1981).
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@param ... All other arguments are ignored.
#'
#'@return A list containing the Average Variance Extracted indices in the first position and the differences
#'between AVEs and largest squared correlations with other composites in the second position.
#'
#'@references
#'
#'Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of marketing research}, 18(1), 39–50.
#'
#'@family post-estimation functions
#'
#'@export
#'

AVE <- function(object, ...){
	UseMethod("AVE")
}


#'@S3method AVE matrixpls

AVE.matrixpls <- function(object, ...) {
	
	loadings <- loadings.matrixpls(object, standardized = TRUE)
	reflectiveModel <- attr(object, "model")$reflective
	
	aves <- unlist(lapply(1:ncol(loadings), function(col){
		useLoadings <- loadings[,col][reflectiveModel[,col]==1]
		
		sum(useLoadings^2)/
			(sum(useLoadings^2) + sum(1-useLoadings^2))
		
	}))

	names(aves) <- colnames(loadings)
	
	C <- attr(object,"C")
	C2 <- C^2
	diag(C2) <- 0
	maxC <- apply(C2,1,max)

	aves_correlation <- aves - maxC
	
	names(aves_correlation) <- colnames(loadings)
	
	result = list(AVE = aves,
								AVE_correlation = aves_correlation)
	
	class(result) <- "matrixplsave"

	
	C <- attr(object,"C")
	result
}

#'@S3method print matrixplsave

print.matrixplsave <- function(x, ...){
	cat("\n Average Variance Extracted indices\n")
	print.table(x$AVE, ...)
	cat("\n AVE - largest squared correlation\n")
	print.table(x$AVE_correlation, ...)
}
