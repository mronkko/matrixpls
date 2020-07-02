#'@title Parameter estimation with PLSe2
#'
#'@description
#'Estimates the model parameters with weighted least squares estimator using the weight matrix
#'\code{W} as the WLS weights.
#'
#'@details
#'\code{params.PLSe2} estimates the statistical model described by \code{model} using the PLSe2
#'algorithm. This algorithm is simply WLS estimation with the weights \code{W}. The model is 
#'estimated with \code{\link[lavaan]{lavaan}}.
#'
#'@param model Lavaan script or lavaan parameter table defining the model.
#'
#'@inheritParams params.separate
#'@inheritParams matrixpls
#'
#'@references
#' Bentler, P. M., & Huang, W. (submitted). On Components, Latent Variables, PLS and Simple Methods:
#' Reactions to Ridgon’s Rethinking of PLS. \emph{Long Range Planning}.
#'
#'@return A named vector of parameter estimates.
#'
#'@family parameter estimators
#'
#'@export
#'

params.plse2 <- function(S, W, model){
	
	stop("Not implemented")
	
	W.vect <- apply(W,2,function(x){
		mean(x[x!=0])
	})
	
	lavaan.out = lavaan(model, sample.cov = S, estimator = "GLS", WLS.V = diag(W.vect))	
	
}
