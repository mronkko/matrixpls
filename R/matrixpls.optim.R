# =========== Optimization criterion functions ===========

#'@title Optimization criteria functions
#'
#'@description Optimization criterion functions calculate various optimization criterion values
#'from \code{matrixpls} objects. 
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@return Value of the optimization criterion.
#'
#'@seealso \code{\link{weightFun.optim}}
#'@name optimCrit
NULL


#'@describeIn optimCrit maximizes the sum of R2 statistics of the \code{inner} matrix
#'@export

optimCrit.maximizeInnerR2 <- function(matrixpls.res){
  -sum(r2(matrixpls.res))
}


#'@describeIn optimCrit maximizes the sum of R2 statistics of the \code{reflective} matrix.
#'@export

optimCrit.maximizeIndicatorR2 <- function(matrixpls.res){
  lambda <- loadings(matrixpls.res)
  IC <- attr(matrixpls.res,"IC")
  -sum(diag(lambda %*% IC))
}

#'@describeIn optimCrit maximizes the sum of R2 statistics of the \code{inner} and \code{reflective} matrices.
#'@export


optimCrit.maximizeFullR2 <- function(matrixpls.res){
  optimCrit.maximizeIndicatorR2(matrixpls.res) + optimCrit.maximizeInnerR2(matrixpls.res)
}

#'@describeIn optimCrit minimies the generalized structured component analysis criterion. See \link{GSCA}
#'@export

optimCrit.gsca <- function(matrixpls.res){
  
  C <- attr(matrixpls.res,"C")
  IC <- attr(matrixpls.res,"IC")
  nativeModel <- attr(matrixpls.res,"model")
  
  reflective <- nativeModel$reflective
  reflective[which(reflective==1)] <- matrixpls.res[grep("=~",names(matrixpls.res))]
  
  #  formative <- nativeModel$formative
  #  formative[which(formative==1)] <- matrixpls.res[grep("<~",names(matrixpls.res))]
  
  #  f <- apply(nativeModel$formative != 0,1,any)
  r <- apply(nativeModel$reflective != 0,1,any)
  endo <- apply(nativeModel$inner != 0,1,any)
  
  inner_resid <- (1 - r2(matrixpls.res)[endo])
  #  form_resid <- (1 - rowSums(IC[f,] * formative[f,]))
  refl_resid <- (1 - rowSums(t(IC[,r]) * reflective[r,]))
  
  #  sum(inner_resid, form_resid, refl_resid)
  sum(inner_resid, refl_resid)
}
