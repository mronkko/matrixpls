#
# Various post-estimation functions
#

#'@title 	Total, Direct, and Indirect Effects for PLS latent variable model
#'
#'@description
#'
#'The \code{matrixpls} method for the standard generic function \code{effects} computes total, direct, 
#'and indirect effects for a PLS latent variable model according to the method described in Fox (1980).
#'
#'Adapted from the \code{\link[sem]{effects}} function of the \code{sem} package
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@return A list with \code{Total}, \code{Direct}, and \code{Indirect} elements.
#'
#'@references
#'
#'Fox, J. (1980) Effect analysis in structural equation models: Extensions and simplified methods of computation. \emph{Sociological Methods and Research}
#'9, 3--28.
#'
#'
#'@export

effects <- function(x, ...){
	UseMethod("effects")
}

#'@S3method effects matrixpls

effects.matrixpls <- function(object = NULL, beta = NULL, innerModel = NULL, ...) {
	
	if(!is.null(object)){
		A <- attr(object,"beta")
		endog <- rowSums(attr(object,"model")$inner)!=0 
	}
	else{
		A <- beta
		endog <- rowSums(innerModel)!=0 		
	}
	
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

print.matrixplseffects <- function(x, digits=getOption("digits"), ...){
	cat("\n Total Effects (column on row)\n")
	Total <- x$Total
	Direct <- x$Direct
	Indirect <- x$Indirect
	select <- !(apply(Total, 2, function(x) all( x == 0)) & 
								apply(Direct, 2, function(x) all( x == 0)) & 
								apply(Indirect, 2, function(x) all( x == 0)))
	print(Total[, select], digits=digits)
	cat("\n Direct Effects\n")
	print(Direct[, select], digits=digits)
	cat("\n Indirect Effects\n")
	print(Indirect[, select], digits=digits)
	invisible(x)
}

#'@title Residual diagnostics for PLS latent variable model
#'
#'@description
#'
#'The \code{matrixpls} method for generic function \code{residuals} computes the residual
#'covariance matrix and various fit indices presented by Lohm\\"oller (1989, ch 2.4).
#''
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@return A list with \code{Total}, \code{Direct}, and \code{Indirect} elements.
#'
#'@references
#'
#'Lohm\\"oller J.-B. (1989) \emph{Latent variables path modelin with partial
#'least squares.} Heidelberg: Physica-Verlag.
#'
#'@export
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
	
	R2 <- r2(object)
	
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

print.matrixplsresiduals <- function(x, digits=getOption("digits"), ...){
	cat("\n Inner model (latent variable) residual covariance matrix\n")
	print(x$inner, digits = digits)
	cat("\n Outer model (indicator) residual covariance matrix\n")
	print(x$outer, digits = digits)
	cat("\n Residual-based fit indices\n")
	print(data.frame(Value = unlist(x$indices)))
}

#'@title R2	for PLS latent variable model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{r2} computes the squared multiple correlation (R2)
#'for latent variables predicted by other latent variables in the model
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@return A named numberic vector
#'
#'
#'@export
#'

r2 <- function(x, ...){
	UseMethod("r2")
}

#'@S3method r2 matrixpls

r2.matrixpls <- function(object, ...){
	
	R2 <- rowSums(attr(object,"beta") * attr(object,"C"))
	names(R2) <- colnames(attr(object,"model")$inner)
	class(R2) <- "matrixplsr2"
	R2
}

#'@S3method print matrixplsr2

print.matrixplsr2 <- function(x, digits=getOption("digits"), ...){
	cat("\n Inner model squared multiple correlations (R2)\n")
	print.table(x, digits = digits)
}

#'@title Goodness of Fit indices for PLS latent variable model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{gof} computes Goodness of Fit 
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@return A list with \code{Total}, \code{Direct}, and \code{Indirect} elements.
#'
#'@references
#'
#'Henseler, J., & Sarstedt, M. (2013). Goodness-of-fit indices for partial least squares path modeling. \emph{Computational Statistics}, 28(2), 565–580. doi:10.1007/s00180-012-0317-1

#'@export
#'

gof <- function(x, ...){
	UseMethod("gof")
}

#'@S3method gof matrixpls

gof.matrixpls <- function(object, ...) {
	
	nativeModel <- attr(object,"model")
	IC <- attr(object,"IC")
	S <- attr(object,"S")

	exogenousLVs <- rowSums(nativeModel$inner) == 0
	
	
	IC_std <- IC %*% (diag(1/sqrt(diag(S))))
	
	result <- sqrt(mean(IC_std[t(nativeModel$reflective)==1]^2) * 
			 	mean(r2(object)[!exogenousLVs]))
	class(result) <- "matrixplsgof"
	result
}

#'@S3method print matrixplsgof

print.matrixplsgof <- function(x, digits=getOption("digits"), ...){
	cat("\n Absolute goodness of fit:", round(x, digits = digits))
	cat("\n")
}

#'@title Factor loading matrix from PLS latent variable model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{loading} computes standardized factor loading matrix
#'from the model results.
#'
#'@param standarized logical indicating whether the factor loadings should be standardized
#'
#'@return A matrix of factor loadings
#'
#'@export
#
loadings <- function(x, ...){
	UseMethod("loadings")
}

#'@S3method loadings matrixpls

loadings.matrixpls <- function(object, ...) {

	nativeModel <- attr(object,"model")
	IC <- attr(object,"IC")
	
	#Standardize
	S <- attr(object,"S")
	IC_std <- IC %*% (diag(1/sqrt(diag(S))))
	
	res <- nativeModel$reflective
	res[res==1] <- t(IC_std)[res==1]
	res
}

#'@title Composite Reliability indices for PLS latent variable model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{cr} computes Composite Reliability 
#'indices for the model
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@return A named numeric vector containing the Composite Reliability indices
#'
#'@references
#'
#'Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of marketing research}, 18(1), 39–50.
#'
#'@export
#'

cr <- function(x, ...){
	UseMethod("cr")
}

#'@S3method cr matrixpls

cr.matrixpls <- function(object, ...) {
	
	loadings <- loadings(object)
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

print.matrixplscr <- function(x, digits=getOption("digits"), ...){
	cat("\n Composite Reliability indices\n")
	print.table(x, digits = digits)
}


#'@title Average Variance Extracted indices for PLS latent variable model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{cr} computes Average Variance Extracted 
#'indices for the model
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@return A list containing the Average Variance Extracted indices in first position and the differences
#'between AVEs and largest correlations with other latent variables in the second position
#'
#'@references
#'
#'Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of marketing research}, 18(1), 39–50.
#'
#'@export
#'

ave <- function(x, ...){
	UseMethod("ave")
}


#'@S3method ave matrixpls

ave.matrixpls <- function(object, ...) {
	
	loadings <- loadings(object, standardized = TRUE)
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

print.matrixplsave <- function(x, digits=getOption("digits"), ...){
	cat("\n Average Variance Extracted indices\n")
	print.table(x$AVE, digits = digits)
	cat("\n AVE - largest squared correlation\n")
	print.table(x$AVE_correlation, digits = digits)
}

#'@title Summary of model fit of PLS latent variable model
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{fitSummary} computes a list of statistics
#'that can be used to access the overall fit of the model.
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@return A list containing a selection of fit indices
#'
#'@export
#'

fitSummary <- function(x, ...){
	UseMethod("fitSummary")
}

#'@S3method fitSummary matrixpls

fitSummary.matrixpls <- function(x, ...){

	
	
	ret <- list("Min CR" = min(cr(x)),
							"Min AVE" = min(ave(x)$AVE),
							"Min AVE - sq. cor" = min(ave(x)$AVE_correlation),
							"Goodness of Fit" = gof(x),
							SRMR = residuals(x)$indices$SRMR)
	
	class(ret) <- "matrixplfitsummary"
	ret
}

#'@S3method print matrixplfitsummary

print.matrixplfitsummary <- function(x, digits=getOption("digits"), ...){
	cat("\n Summary indices of model fit\n")
	print(data.frame(Value = unlist(x)))
}