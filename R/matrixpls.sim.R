#'@title Monte Carlo simulations with matrixpls
#'
#'@description
#'Performs Monte Carlo simulations of \code{\link{matrixpls}} with the \code{\link[simsem]{sim}} function of the \code{simsem} package.
#'The standard errors and confidence intervals are estimated with the \code{\link[boot]{boot}} and \code{\link[boot]{boot.ci}} functions
#'of the \code{boot} package.
#
#'@details
#' This funtion calls the \code{\link[simsem]{sim}} function from the \code{simsem} package to perform Monte
#' Carlo simulations with matrixpls. The function parses the model parameters and replaces it with
#' a function call that estimates the model and bootstrapped standard errors and confidence
#' intervals with \link{matrixpls.boot}.
#' 
#' If the \code{generate} or \code{rawdata} arguments are not specified in the \code{\link[simsem]{sim}} arguments
#'  then the \code{model} argument will be used for data generation and must be specified in lavaan format.
#'
#'@param nRep Number of replications. If any of the \code{n}, \code{pmMCAR}, or \code{pmMAR} arguments are specified as lists, the number of replications will default to the length of the list(s), and \code{nRep} need not be specified.
#'
#'@inheritParams matrixpls
#'
#'@param ... All other arguments are forwared through to \code{\link[simsem]{sim}} or \code{\link{matrixpls.boot}}.
#'
#'@param cilevel Confidence level. This argument will be forwarded to the \code{\link[boot]{boot.ci}} when calculating the confidence intervals.
#'
#'@param citype Type of confidence interval. This argument will be forwarded to the \code{\link[boot]{boot.ci}} when calculating the confidence intervals.
#'
#'@param boot.R Number of bootstrap replications to use to estimate standard errors.
#'
#'@param fitIndices A function that returns a list of fit indices for the model. Setting this argument to \code{NULL} disables fit indices.
#'
#'@return An object of class \code{\link[simsem]{SimResult-class}}.
#'
#'@inheritParams simsem::sim
#'
#'@example example/matrixpls.sim-example.R
#'
#'@seealso
#'
#'\code{\link[simsem]{sim}}, \code{\link[simsem]{SimResult-class}}, \code{\link{matrixpls.boot}}, \code{\link[boot]{boot}}, \code{\link[boot]{boot.ci}}
#'
#'@include matrixpls.R
#'
#'@export


matrixpls.sim <- function(nRep = NULL, model = NULL, n = NULL, ..., cilevel = 0.95, citype=c("norm","basic", "stud", "perc", "bca"), boot.R = 500, fitIndices = fitSummary){

	library(simsem)
	library(assertive)
	
	# Basic verification of the arguments
	assert_all_are_positive(boot.R)	
	assert_all_are_positive(cilevel)
	assert_all_are_true(cilevel<1)
	citype <- match.arg(citype)

	# Decide which arguments to route to simsem and which to matrispls.boot
	
	allArgs <- list(...)	
	simsem.argNames <- names(formals(simsem::sim))
	
	whichArgsToSimsem <- names(allArgs) %in% simsem.argNames
	simsemArgs <- allArgs[whichArgsToSimsem]
	matrixplsArgs <- allArgs[! whichArgsToSimsem]
	
	# Beacuse we are using a custom estimator function, the generate model must be set.
	
	if(!"generate" %in% names(simsemArgs) &&
		 	!"rawData" %in% names(simsemArgs)) simsemArgs$generate <- model
	
	nativeModel <- parseModelToNativeFormat(model)
	
	if(!"W.mod" %in% names(matrixplsArgs)) matrixplsArgs$W.mod<- defaultWeightModelWithModel(model)
		
	matrixplsArgs <- c(list(R= boot.R, model = nativeModel),matrixplsArgs)
	
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

		boot.out <- do.call(matrixpls.boot, c(list(as.matrix(data)), matrixplsArgs))
		
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
								ciupper = cis[2,])

		if(! is.null(fitIndices)){
			assert_is_function(fitIndices)
			fitlist <- unlist(fitIndices(boot.out$t0))
			ret$fit <- fitlist
		}
		return(ret)
	}
	simsemArgs <- c(list(nRep = nRep, model = model, n = n), simsemArgs)
	do.call(simsem::sim, simsemArgs)

}

#'@title Summary of model fit of PLS model
#'
#'@description
#'
#'\code{fitSummary} computes a list of statistics
#'that are commonly used to access the overall fit of the PLS model.
#'
#'@param object PLS estimation result object produced by the \code{\link{matrixpls}} function.
#'
#'@return A list containing a selection of fit indices.
#'
#'@include matrixpls.postestimation.R
#'
#'@export
#'

fitSummary <- function(object){
	
	ret <- list("Min CR" = min(CR(object)),
							"Min AVE" = min(AVE(object)$AVE),
							"Min AVE - sq. cor" = min(AVE(object)$AVE_correlation),
							"Goodness of Fit" = GoF(object),
							SRMR = residuals(object)$indices$SRMR)
	
	class(ret) <- "matrixplfitsummary"
	ret
}

#'@S3method print matrixplfitsummary

print.matrixplfitsummary <- function(x, ...){
	cat("\n Summary indices of model fit\n")
	print(data.frame(Value = unlist(x)), ...)
}