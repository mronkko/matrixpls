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
#'@param ... All other parameters are passed through to \code{\link{matrixpls}}.
#'
#'@param signChangeCorrection Sign change correction function.
#'
#'@param stopOnError A logical indicating whether boostrapping should be continued when error occurs
#' in a replication.
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
#'Sign change corrections: \code{\link{signChange.individual}}; \code{\link{signChange.construct}}
#'@export
#'

matrixpls.boot <- function(data, R = 500, ..., 
                           signChangeCorrection = NULL,
                           parallel = c("no", "multicore", "snow"),
                           ncpus = getOption("boot.ncpus", 1L),
                           stopOnError = FALSE){
  
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
                       do.call(matrixpls, arguments)
                     }
                     else{
                       tryCatch(
                         do.call(matrixpls, arguments)
                       )
                     }
                   },
                   R, parallel = parallel, ncpus = ncpus)
  
  #	class(boot.out) <- c("matrixplsboot", class(boot.out))
  boot.out
}

#'@title Individual indicator sign change correction for boostrapping
#'
#'@description Changes the signs of \code{W} to match \code{Worig}
#'
#'@param Worig the original weight matrix
#'
#'@param W a weight matrix of a bootstrap sample
#'
#'@return A weight matrix with the same dimensions as \code{W} after
#'
#'@seealso
#'\code{\link{matrixpls.boot}}
#'
#'@family signChangeCorrection
#'
#'@references Tenenhaus, M., Esposito Vinzi, V., Chatelin, Y.-M., & Lauro, C. (2005). PLS Path Modeling.
#'\emph{Computational Statistics & Data Analysis}, 48(1), 159–205. doi:10.1016/j.csda.2004.03.005
#'
#'@export
#'

signChange.individual <- function(Worig,W){
  W * sign(Worig) * sign(W)
}

#'@title Construct level sign change correction for boostrapping
#'
#'@description Changes the signs of \code{W} on all rows where the sign of the
#'sum of the row differs between \code{Worig} and \code{W}
#'
#'@param Worig the original weight matrix
#'
#'@param W a weight matrix of a bootstrap sample
#'
#'@return A weight matrix with the same dimensions as \code{W} after
#'
#'@seealso
#'\code{\link{matrixpls.boot}}
#'
#'@family signChangeCorrection
#'
#'@references Tenenhaus, M., Esposito Vinzi, V., Chatelin, Y.-M., & Lauro, C. (2005). PLS Path Modeling.
#'\emph{Computational Statistics & Data Analysis}, 48(1), 159–205. doi:10.1016/j.csda.2004.03.005
#'
#'@export
#'
signChange.construct <- function(Worig,W){
  origSums <- apply(Worig,1,sum)
  curSums <- apply(W,1,sum)
  signs <- ifelse(abs(origSums-curSums) > abs(origSums+curSums), -1,1)
  sweep(W,MARGIN=1,signs,`*`)
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

