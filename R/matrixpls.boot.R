#'@title Bootstrapping of matrixpls function
#'
#'@description
#' Runs bootstrap replications from rawdata using matrixplis
#
#'@details
#'
#' Uses the \code{\link[boot]{boot}} method from the 
#'
#'@param data Matrix or data frame containing the raw data
#'
#'@param R Number of bootstrap samples to draw
#'
#'@param ... All other parameters are passed through to \code{matrixpls}
#'
#'@return An object of class \code{\link[boot]{boot}}.
#'
#'#'  \item{t0}{
#'    The observed value of \code{statistic} applied to \code{data}. 
#'  }
#'  \item{t}{
#'    A matrix with \code{sum(R)} rows each of which is a bootstrap replicate
#'    of the result of calling \code{statistic}.
#'  }
#'  \item{R}{
#'    The value of \code{R} as passed to \code{boot}.
#'  }
#'  \item{data}{
#'    The \code{data} as passed to \code{boot}.
#'  }
#'  \item{seed}{
#'    The value of \code{.Random.seed} when \code{boot} was called.  
#'  }
#'  \item{statistic}{
#'    The function \code{statistic} as passed to \code{boot}.
#'  }
#'  \item{sim}{
#'    Simulation type used.
#'  }
#'  \item{stype}{
#'    Statistic type as passed to \code{boot}.
#'  }
#'  \item{call}{
#'    The original call to \code{boot}.
#'  }
#'  \item{strata}{
#'    The strata used.  This is the vector passed to \code{boot}, if it
#'    was supplied or a vector of ones if there were no strata.  It is not
#'    returned if \code{sim} is \code{"parametric"}.
#'  }
#'  \item{weights}{
#'    The importance sampling weights as passed to \code{boot} or the empirical 
#'    distribution function weights if no importance sampling weights were
#'    specified.  It is omitted if \code{sim} is not one of
#'    \code{"ordinary"} or \code{"balanced"}.
#'  }
#'  \item{pred.i}{
#'    If predictions are required (\code{m > 0}) this is the matrix of
#'    indices at which predictions were calculated as they were passed to
#'    statistic.  Omitted if \code{m} is \code{0} or \code{sim} is not
#'    \code{"ordinary"}. 
#'  }
#'  \item{L}{
#'    The influence values used when \code{sim} is \code{"antithetic"}.
#'    If no such values were specified and \code{stype} is not \code{"w"}
#'    then \code{L} is returned as consecutive integers corresponding to
#'    the assumption that data is ordered by influence values. This
#'    component is omitted when \code{sim} is not \code{"antithetic"}.
#'  }
#'  \item{ran.gen}{
#'    The random generator function used if \code{sim} is
#'    \code{"parametric"}. This component is omitted for any other value
#'    of \code{sim}.
#'  }
#'  \item{mle}{
#'    The parameter estimates passed to \code{boot} when \code{sim} is
#'    \code{"parametric"}.  It is omitted for all other values of
#'    \code{sim}.
#'  }
#'
#'
#'@export
#'

matrixpls.boot <- function(data, R = 500, ..., parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 1L), validateInput = FALSE, stopOnError = FALSE){
	
	library(boot)
	
	#
	# The precalculation of the covariance matrix products is not  efficient if implemented with R code
	# This is left here to save the logic if this is ever implemented in C code
	#
	
	PRECALCULATE <- FALSE
	
	data <- as.matrix(data)
	
	if(PRECALCULATE){
		# Precalculate the products used in forming the covariance matrix 
		
		data.centered <- apply(data, 2, scale, scale=FALSE, center = TRUE)
		
		zeroMatrix <- matrix(0,ncol(data),ncol(data))
		lt <- which(lower.tri(zeroMatrix, diag = TRUE), useNames = FALSE)
		ut <- matrix(1:ncol(data)^2,ncol(data),ncol(data), byrow = TRUE)[lower.tri(zeroMatrix, diag = TRUE)]
		
		rows <- row(zeroMatrix)[lt]
		cols <- col(zeroMatrix)[lt]
		
		# The rows are observations and columns are the elements in the covariance matrix 
		
		products <- mcmapply(function(row,col){
			data.centered[,row]*data.centered[,col]/(nrow(data)-1)
		}, rows, cols)
	}
	
	boot.out <- boot(data,
									 function(data, indices){
									 	
									 	if(PRECALCULATE){
									 		covs <- apply(products[indices,],2,sum)
									 		S <- matrix(0,ncol(data),ncol(data))
									 		S[lt] <- covs
									 		S[ut] <- covs
									 	}
									 	else{
									 		S <- cov(data[indices,])
									 	}
									 	
									 	if(stopOnError){
										 	tryCatch(
										 		matrixpls(S, ...,  validateInput = validateInput)
										 	)
									 	}
										 else{
										 	matrixpls(S, ..., , validateInput = validateInput)
										 }
									 },
									 R, parallel = parallel, ncpus = ncpus)
	
	class(boot.out) <- c("matrixplsboot", class(boot.out))
	boot.out
}

#'@S3method print matrixplsboot

print.matrixplsboot <- function(object, digits=getOption("digits"), ...){
	matrixpls.out <- object$t0
	attr(matrixpls.out,"boot.out") <- object
	print(matrixpls.out, digits = digits)
}

#'@S3method summary matrixplsboot

summary.matrixplsboot <- function(object, ..., digits = max(3, getOption("digits")-3)){
	matrixpls.out <- object$t0
	attr(matrixpls.out,"boot.out") <- object
	summary(matrixpls.out, digits = digits)
}
