library(simsem)

#'@title Monte Carlo simulations with MatrixPLS
#'
#'@description
#' A wrapper function to use MatrixPLS with simsem
#
#'@details
#' This funtion calls the \code{\link[simsem]{sim}} function from the \code{simsem} package to perform Monte
#' Carlo simulations with MatrixPLS. The function parses the model parameter and replaces it with
#' a function call that estimates the same model with MatrixPLS.
#' 
#'@export


matrixpls.sim <- function(nRep = NULL, model = NULL, n = NULL, generate = NULL, ..., cilevel = 0.95, citype=c("norm","basic", "stud", "perc", "bca"), bootstrap = 500){

	# Basic verification of the arguments
	assert_all_are_positive(bootstrap)	
	assert_all_are_positive(cilevel)
	assert_all_are_true(cilevel<1)
	citype <- match.arg(citype)

	# Beacuse we are using a custom estimator function, the generate model must be set.
	
	if(is.null(generate)) generate <- model
	
	nativeModel <- parseModelToNativeFormat(model)	
	weightRelations <- defaultWeightRelationsWithModel(model)
		
  # A function that takes a data set and returns a list. The list must
	# contain at least three objects: a vector of parameter estimates (coef),
	# a vector of standard error (se), and the convergence status as TRUE or
	# FALSE (converged). There are five optional objects in the list: a vector
	# of fit indices (fit), a vector of standardized estimates (std), any
	# extra output (extra), fraction missing type I (FMI1), and fraction
	# missing type II (FMI2). Note that the coef, se, std, FMI1, and FMI2 must
	# be a vector with names. The name of those vectors across different
	# objects must be the same.

	model  <- function(data){
			
			
			boot.out <- boot(data, function(originalData,indices){
					S <- cor(originalData[indices,])
					matrixpls(S, nativeModel, weightRelations)
			}, bootstrap) 
			
			# Indices for parameters excluding weights
			parameterIndices <- 1:(sum(nativeModel$inner) + sum(nativeModel$reflective) + sum(nativeModel$formative))
			
			boot.ci.out <- boot.ci(boot.out, conf = cilevel, type = citype)
			cis <- boot.ci.out[[citype]]
			cis <- cis[,ncol(cis)-1:0]

			return(list(coef = boot.out$t0[parameterIndices],
									se = apply(boot.out$t,2,sd,na.rm=TRUE),
									converged = attr(boot.out$t0,"converged"),
									cilower = cis[,1],
									ciupper = cis[,2]))
	}
	

	
	sim(nRep = nRep, model = model, n = n, generate = generate)
}
 
