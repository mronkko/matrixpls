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
#'@usage sim(nRep, model, ..., boot.R = 500)
#'
#'@param nRep Number of replications. If any of the \code{n}, \code{pmMCAR}, or \code{pmMAR} arguments are specified as lists, the number of replications will default to the length of the list(s), and \code{nRep} need not be specified.
#'
#'@param model \code{lavaan} script or \code{lavaan} parameter table. If the \code{generate} argument is not specified, then the object in the \code{model} argument will be used for both data generation and analysis. If \code{generate} is specified, then the \code{model} argument will be used for data analysis only. 
#'
#'@param ... All other parameters are passed through to \code{\link[simsem]{sim}}
#'
#'@param boot.R Number of bootstrap replications to use to estimate standard errors
#'
#'@export


matrixpls.sim <- function(nRep = NULL, model = NULL, n = NULL, generate = NULL, ..., rawData = NULL, cilevel = 0.95, citype=c("norm","basic", "stud", "perc", "bca"), boot.R = 500, stopOnError = FALSE){

	library(simsem)
	library(assertive)
	
	# Basic verification of the arguments
	assert_all_are_positive(boot.R)	
	assert_all_are_positive(cilevel)
	assert_all_are_true(cilevel<1)
	citype <- match.arg(citype)

	# Beacuse we are using a custom estimator function, the generate model must be set.
	
	if(is.null(generate) && is.null(rawData)) generate <- model
	
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
			
		# Indices for parameters excluding weights
		
		parameterIndices <- 1:(sum(nativeModel$inner) + sum(nativeModel$reflective) + sum(nativeModel$formative))
		
		# Convert the data to matrix for efficiency
		
		boot.out <- matrixpls.boot(data, boot.R, model = nativeModel, weightRelations = weightRelations)
		
		cis <- sapply(parameterIndices, FUN = function(index) {
			boot.ci.out <- boot.ci(boot.out, conf = cilevel, type = citype, index=index)

				# The cis start from the fourth slot and we only have one type of ci. 
				# The list names do not match the type parameter exactly (e.g. "norm" vs. "normal")
				
			cis <- boot.ci.out[[4]]
			cis[,ncol(cis)-1:0]
		})
		
		ses <- apply(boot.out$t[,parameterIndices],2,sd,na.rm=TRUE)
		names(ses) <- names(boot.out$t0[parameterIndices])
		colnames(cis) <- names(boot.out$t0[parameterIndices])
		
		ret <- list(coef = boot.out$t0[parameterIndices],
								se = ses,
								converged = attr(boot.out$t0,"converged"),
								cilower = cis[1,],
								ciupper = cis[2,], 
								fit = unlist(fitSummary(boot.out$t0)))
		
		return(ret)
	}
	

	
	sim(nRep = nRep, model = model, n = n, generate = generate, ..., rawData = rawData, stopOnError = stopOnError)
}
 
