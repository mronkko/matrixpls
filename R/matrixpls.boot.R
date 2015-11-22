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
#'@param ... All other arguments are passed through to \code{\link{matrixpls}}.
#'
#'@param signChangeCorrection Sign change correction function.
#'
#'@param extraFun A function that takes a \code{matrixpls} object and returns a numeric vector. The
#'vector is appended to bootstrap replication. Can be used for boostrapping additional
#'statistics calculated based on the estimation results.
#'
#'@param stopOnError A logical indicating whether boostrapping should be continued when error occurs
#' in a replication.
#'
#'@inheritParams boot::boot
#'
#'
#'@return An object of class \code{\link[boot]{boot}}.
#'
#'@seealso
#'\code{\link[boot]{boot}}
#'
#'Sign change corrections: \code{\link{signChange.individual}}; \code{\link{signChange.construct}}
#'
#'@export
#'
#'
matrixpls.boot <- function(data, ..., R = 500,
                           signChangeCorrection = NULL,
                           parallel = c("no", "multicore", "snow"),
                           ncpus = getOption("boot.ncpus", 1L),
                           stopOnError = FALSE,
                           extraFun = NULL){
  
  if(! requireNamespace("boot")) stop("matrixpls.boot requires the boot package")
  
  if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
  
  data <- as.matrix(data)

  arguments <- list(...)
  
  # Prepare sign change corrections
  
  if(! is.null(signChangeCorrection)){
    Worig <- attr(matrixpls(cov(data),...), "W")
    
    # Get the original weight function
    weightFunction <- arguments[["weightFunction"]]
    if(is.null(weightFunction)) weightFunction <- weight.pls
    
    # Wrap inside sign change correction
    arguments[["weightFunction"]] <- function(S, ...){
      Wrep <- weightFunction(S, ...)
      W <- signChangeCorrection(Worig, Wrep)
      W <- scaleWeights(S, W)
      attributes(W) <- attributes(Wrep)
      W
    }
    
  }
  
  # Bootstrap
  
  boot.out <- boot::boot(data,
                   function(data, indices){
                     
                     S <- cov(data[indices,])
                     arguments <- c(list(S),arguments)
                     
                     if(stopOnError){
                       boot.rep <- do.call(matrixpls, arguments)
                     }
                     else{
                       tryCatch(
                         boot.rep <- do.call(matrixpls, arguments)
                       )
                     }
                     
                     # Add additional statistics
                     
                     if(!is.null(extraFun)){
                       a <- attributes(boot.rep)
                       boot.rep <-c(boot.rep,extraFun(boot.rep))
                       a$names <- names(boot.rep)
                       attributes(boot.rep) <- a
                     }
                     
                     # If the indices are not sorted, then this is not the original sample
                     # and we can safely omit all attributes to save memory
                     
                     if(is.unsorted(indices)) attributes(boot.rep) <- NULL
                     
                     boot.rep
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

