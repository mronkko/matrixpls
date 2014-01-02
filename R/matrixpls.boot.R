#'@title Bootstrapping of matrixpls function
#'
#'@description
#'\code{matrixpls.boot} is a convenience method that implements bootstrapping of \code{matrixpls}
#'with \code{\link[boot]{boot}} method of the \code{boot} package.
#'
#'@param data Matrix or data frame containing the raw data.
#'
#'@param R Number of bootstrap samples to draw.
#'
#'@param ... All other parameters are passed through to \code{matrixpls}.
#'
#'@param stopOnError A logical indicating whether boostrapping should be continued when error occurs in a replication.
#'
#'@inheritParams boot::boot
#'@inheritParams simsem::sim
#'@inheritParams matrixpls.boot
#'
#'@return An object of class \code{\link[boot]{boot}}.
#'
#'@seealso
#'\code{\link[boot]{boot}}
#'
#'@export
#'

matrixpls.boot <- function(data, R = 500, ..., parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 1L), stopOnError = FALSE){
	
	library(boot)
	
	if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
	
	data <- as.matrix(data)
	
	boot.out <- boot(data,
									 function(data, indices){
									 	
									 	S <- cov(data[indices,])

									 	if(stopOnError){
									 		matrixpls(S, ...)
									 	}
									 	else{
									 		tryCatch(
									 			matrixpls(S, ...)
									 		)
									 	}
									 },
									 R, parallel = parallel, ncpus = ncpus)
	
	#	class(boot.out) <- c("matrixplsboot", class(boot.out))
	boot.out
}

# These are not in use 

#'#@S3method print matrixplsboot

#print.matrixplsboot <- function(x, ...){
#	matrixpls.out <- x$t0
#	attr(matrixpls.out,"boot.out") <- x
#	print(matrixpls.out, ...)
#}

#'#@S3method summary matrixplsboot

#summary.matrixplsboot <- function(object, ...){
#	matrixpls.out <- object$t0
#	attr(matrixpls.out,"boot.out") <- object
#	summary(matrixpls.out, ...)
#}

