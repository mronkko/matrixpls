#
# Various post-estimation functions
#

#'@title 	Total, Direct, and Indirect Effects for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the standard generic function \code{effects} computes total, direct, 
#'and indirect effects for a matrixpls results according to the method described in Fox (1980).
#'
#'Adapted from the \code{\link[sem]{effects}} function of the \code{sem} package
#'
#'@template postestimationFunctions
#'
#'@return A list with \code{Total}, \code{Direct}, and \code{Indirect} elements.
#'
#'@references
#'
#'Fox, J. (1980) Effect analysis in structural equation models: Extensions and simplified methods of computation. \emph{Sociological Methods and Research}
#'9, 3--28.
#'
#'@importFrom stats effects
#'
#'@method effects matrixpls
#'
#'@export

effects.matrixpls <- function(object,  ...) {
  
  A <- attr(object,"inner")
  endog <- rowSums(attr(object,"model")$inner)!=0 
  
  I <- diag(endog)
  AA <- - A
  diag(AA) <- 1
  Total <- solve(AA) - I
  Indirect <-  Total - A
  result <- list(Total=Total[endog, , drop = FALSE], Direct=A[endog, , drop = FALSE], Indirect=Indirect[endog, , drop = FALSE])
  class(result) <- "matrixplseffects"
  result
}

#'@export

print.matrixplseffects <- function(x, ...){
  cat("\n Total Effects (column on row)\n")
  Total <- x$Total
  Direct <- x$Direct
  Indirect <- x$Indirect
  select <- !(apply(Total, 2, function(x) all( x == 0)) & 
                apply(Direct, 2, function(x) all( x == 0)) & 
                apply(Indirect, 2, function(x) all( x == 0)))
  print(Total[, select, drop = FALSE], ...)
  cat("\n Direct Effects\n")
  print(Direct[, select, drop = FALSE], ...)
  cat("\n Indirect Effects\n")
  print(Indirect[, select, drop = FALSE], ...)
  invisible(x)
}

#'@title Residual diagnostics for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for generic function \code{residuals} computes the residual
#'covariance matrix and various fit indices. 
#'
#'@details The residuals can be
#'either observed residuals from the regressions of indicators on composites and composites 
#'on composites
#'(i.e. the \code{reflective} and \code{inner} models) as presented by Lohmöller (1989, ch 2.4) or 
#'model implied residuals calculated by subtracting model implied covariance matrix from the 
#'sample covariance matrix as done by Henseler et al. (2014). 
#'
#'The root mean squared residual indices (Lohmöller, 1989, eq 2.118) are calculated from the
#'off diagonal elements of the residual covariance matrix. The
#'standardized root mean squared residual (SRMR) is calculated based on the standardized residuals
#'of the \code{reflective} model matrix. 
#'
#'Following Hu and Bentler (1999, Table 1), the SRMR index is calculated by dividing with 
#'\eqn{p(p+1)/2}, where \eqn{p} is the number of indicator variables. In typical SEM applications,
#'the diagonal of residual covariance matrix consists of all zeros because error term variances
#'are freely estimated. To make the SRMR more comparable with the index produced by 
#'SEM software, the SRMR is calculated by summing only the squares of off-diagonal elements,
#'which is equivalent to including a diagonal of all zeros. 
#'
#'Two versions of the 
#'SRMR index are rovided, the traditional SRMR that includes all residual covariances, and the 
#'version proposed by Henseler et al. (2014) where the within-block residual covariances are 
#'ignored. 
#'
#'@template postestimationFunctions
#'
#'@param observed If \code{TRUE} (default) the observed residuals from the outerEstim.model regressions
#'(indicators regressed on composites) are returned. If \code{FALSE}, the residuals are calculated
#'by combining  \code{inner}, \code{reflective}, and \code{formative} as a simultaneous equations
#'system and subtracting the covariances implied by this system from the observed covariances.
#'The error terms are constrained to be uncorrelated and covariances between exogenous observed
#'values are fixed at their sample values.
#'
#'@return A list with three elements: \code{inner}, \code{outer}, and \code{indices} elements
#' containing the residual covariance matrix of regressions of composites on other composites,
#' the residual covariance matrix of indicators on composites, and various indices
#' calculated based on the residuals.
#'
#'@references
#'
#'Henseler, J., Dijkstra, T. K., Sarstedt, M., Ringle, C. M., Diamantopoulos, A., Straub, D. W., …
#'Calantone, R. J. (2014). Common Beliefs and Reality About PLS Comments on Rönkkö and Evermann 
#'(2013). \emph{Organizational Research Methods}, 17(2), 182–209. \doi{10.1177/1094428114526928}
#'
#'Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure 
#'analysis: Conventional criteria versus new alternatives. \emph{Structural Equation Modeling: 
#'A Multidisciplinary Journal}, 6(1), 1–55.
#'
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.

#'@export
#'
#'@method residuals matrixpls
#'
#'@export

residuals.matrixpls <- function(object, ..., observed = TRUE) {
  
  RMS <- function(num) sqrt(sum(num^2)/length(num))
  
  S <- attr(object,"S")
  
  # Lohmöller defines quite a few statistics based on correlations.
  # Because S is a covariance matrix, we need to calculate the 
  # corresponding correlation matrix as well.
  
  Scor <- stats::cov2cor(S)
  
  nativeModel <- attr(object,"model")
  
  # Equation numbers in parenthesis refer to equation number in Lohmoller 1989
  
  if(observed){
    
    W <- attr(object,"W")
    
    # Number of reflective indicators
    reflectiveIndicators<- rowSums(nativeModel$reflective)>0
    k <- sum(reflectiveIndicators)
    
    # Number of endog LVs
    endog <- rowSums(nativeModel$inner)>0
    h <- sum(endog)	
    
    
    
    # Factor loading matrix
    P <- nativeModel$reflective
    
    P[P==1] <- object[grepl("=~", names(object), fixed=TRUE)]
    
    # Standardized loadings
    Pstd <- sweep(P,MARGIN=1,sqrt(diag(S)),`/`)
    
    # This is always standardized, so no need to rescale 
    B <- attr(object,"inner")
    
    # Lohmoller is not clear whether R should be based on the estimated betas or calculated scores
    # The scores are used here because this results in less complex code
    
    R <- attr(object,"C")
    R_star <- (B %*% R %*% t(B))[endog,endog] # e. 2.99
    
    
    # Model implied indicator correlations
    H <- Pstd %*% R %*% t(Pstd) # eq 2.96
    
    I <- diag(ncol(Scor))
    H2 <- (I * H)  %*% solve(I * Scor) # eq 2.97
    
    F <- Pstd %*% B %*% R %*% t(B) %*% t(Pstd) # eq 2.104
    
    F2 <- (I * F) %*% solve(I * Scor) # eq 2.105
    
    r2 <- r2(object)
    
    
    # Lohmoller 1989 uses C for the residual covariance matrix of indicators
    # matrixpls uses C for the composite correlation matrix, but from here on
    # until the end of the function, C is used for the residual covariance matrix
    
    # Lohmoller does not define C in covariance form, so we need to do it ourselfs.
    # 
    # Start with the observed residuals:
    #
    # e = X-XW'P'
    # 
    # because residuals have a mean of zero cov(e) can be defined as
    # 
    # cov(e) = (X-XW'P')’(X-XW'P')
    # cov(e) = (X’-(XW'P')')(X-XW'P')
    # cov(e) = (X’-PWX')(X-XW'P')
    # cov(e) = X’X-X’XW'P'-PWX'X+PWX'XW'P'
    # cov(e) = S-SW'P'-PWS+PWSXW'P'
    # cov(e) = S-SW'P'-(SW'P')’+PRP'

    C <- S - S%*%t(W)%*%t(P) - t(S-S%*%t(W)%*%t(P)) + P%*%R%*%t(P)
    
    Q <- (W %*% S %*% t(W))[endog,endog] - R_star
    
    
    indices <- c(Communality = psych::tr(H2)/k, # eq 2.109
                    Redundancy = psych::tr(F2)/k,  # eq 2.110
                    SMC = sum(r2)/h,         # eq 2.111
                    "RMS outer residual covariance" = RMS(C[lower.tri(C)]), # eq 2.118
                    "RMS inner residual covariance" = RMS(C[lower.tri(Q)]) # eq 2.118
                    )
  }
  else{
    C <- S-stats::fitted(object)
    indices <- c()
  }
  
  C_std <- diag(1/diag(S)) %*% C
  
  indices <- c(indices,c(
                  
                  # SRMR as calculated in SEM. (Hu and Bentler, 1999, p. 3)
                  
                  SRMR = sqrt(sum(C_std[lower.tri(C_std)]^2)/length(C[lower.tri(C_std, diag=TRUE)])),
                  
                  # SRMR calculated ignoring within block residuals from Henserler et al 2014.
                  
                  "SRMR (Henseler)" = RMS(C_std[nativeModel$reflective %*% t(nativeModel$reflective)==0]))
  )
  
  if(observed){
    result<- list(inner = Q, outer = C, indices = indices)
  }
  else{
    result<- list(outer = C, indices = indices)
  }
  
  class(result) <- "matrixplsresiduals"
  result
}

#'@export

print.matrixplsresiduals <- function(x, ...){

    if("inner" %in% names(x)){
    cat("\n Inner model (composite) residual covariance matrix\n")
    print(x$inner, ...)
    cat("\n Outer model (indicator) residual covariance matrix\n")
    print(x$outer, ...)
  }
  else{
    cat("\n Model implied residual covariance matrix\n")
    print(x$outer, ...)
  }

  cat("\n Residual-based fit indices\n")
  i <- as.matrix(x$indices)
  colnames(i) <- "Value"
  print(i, ...)
}

#'@title Model implied covariance matrix based on matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for generic function \code{fitted} computes the model implied
#'covariance matrix by combining  \code{inner}, \code{reflective}, and \code{formative} as a 
#'simultaneous equations system. The error terms are constrained to be uncorrelated and 
#'covariances between exogenous variables are fixed at their sample values. Defining a
#'composite as dependent variable in both inner and formative creates an impossible model
#'and results in an error.
#'
#'@template postestimationFunctions
#'
#'@return a matrix containing the model implied covariances.
#'
#'@export
#'
#'@method fitted matrixpls
#'
#'@export

fitted.matrixpls <- function(object, ...) {
  
  # Contains: $inner (J x J), $reflective (K x J), $formative (J x K)
  
  nativeModel <- attr(object,"model")
  
  # Check that the matrices form a valid model
  
  if(any(rowSums(nativeModel$inner) > 0 & rowSums(nativeModel$formative) > 0))
    stop("Cannot calculate model implied covariance matrix. A composite is a dependent variable in both formative and inner matrices resulting in an impossible model. See http://urn.fi/URN:NBN:fi:aalto-201605031907, p. 823 for an explanation and pp. 70-76 for further technical details.")
  
  ### Matrices -------------------------------------------------------------------------------------
  ## S      := (K x K) Empirical indicator VCV matrix: Cov(x)
  ## IC     := (J x K) Empirical indicator-construct VCV matrix
  ## C      := (J x J) Empirical construct covariance/correlation matrix (V(eta) = WSW')
  ## B      := (J x J) Matrix of (estimated) path coefficients (zero if there is no path)
  ## F      := (J x K) Matrix containing the estimated parameters of the formative (composite) model
  ## P      := (K X J) Matrix of factor and/or composite loadings (usually: Lambda)
  ##
  ## ---- Bentler & Weeks notation
  ##
  ## y      := (p x 1) Vector of observed/measured dependent indicators (manifest variables)
  ## x      := (q x 1) Vector of observed independent indicators (usually 0 in a "normal" framework)
  ## p      := Number of observed dependent variables (usually all indicators)
  ## q      := Number of observed independent variables (0 in a complete latent-variable model)
  ## m      := Number of dependent variables ("endogenous" constructs + indicators)
  ## n      := Number of independent variables ("exogenous" constructs + zetas + deltas)
  ## r      := Number of observed variables, r = p + q
  ## s      := Total number of variables, s = m + n
  ## beta   := (m x m) Coefficient matrix of dep. variables on dep. variables.
  ## gamma  := (m x n) Coefficient matrix of indep. variables on dep. variables.
  ## Phi    := (n x n) VCV of independent variables .
  ## Gamma  := ()
  ## Beta   := ()
  ## G      := (r x s) 
  ### ----------------------------------------------------------------------------------------------
  
  ## Get relevant matrices
  S <- attr(object,"S")
  IC <- attr(object,"IC")
  C <- attr(object,"C")
  B <- attr(object,"inner")
  F <- attr(object,"formative")
  L <- attr(object,"reflective")
  
  
  #### Define full matrices to select from ---------------------------------------------------------
  
  # Matrices containing all regressions and covariances
  # indicators first, then composites
  
  fullB <- rbind(cbind(matrix(0,nrow(S),nrow(S)),L),
                 cbind(F,B))
  
  fullC <- rbind(cbind(S,t(IC)),
                 cbind(IC,C))
  
  ## Selector for all exogenous variables
  
  exog <- rowSums(fullB) == 0
  
  ## Compute outer and inner errors variance and put them in a vector (Var(delta) and Var(zeta)).
  
  e <- c(diag(S) - rowSums(L * t(IC)),
         1-r2(object))
  
  fullC <- rbind(cbind(fullC,matrix(0, nrow(fullC), sum(!exog))),
                 cbind(matrix(0, sum(!exog), nrow(fullC)), diag(e)[!exog,!exog]))
  
  fullB <- cbind(fullB,diag(length(e))[,! exog])
  fullB <- rbind(fullB, matrix(0,sum(!exog), ncol(fullB)))
  
  # Update exog because fullB is now larger
  exog <- rowSums(fullB) == 0
  
  #
  # Derive the implied covariance matrix using the Bentler-Weeks model
  #
  # Bentler, P. M., & Weeks, D. G. (1980). Linear structural equations with 
  # latent variables. Psychometrika, 45(3), 289–308. doi:10.1007/BetaF02293905
  #
  
  # beta (regression paths between dependent variables)
  
  beta <- fullB[!exog,!exog]
  
  
  # gamma (regression paths between dependent and independent  variables)
  
  gamma <- fullB[!exog,exog]
  
  # exogenous variable covariance matrix
  Phi <- fullC[exog,exog]
  
  #
  # Calculate sigma
  #
  
  # Matrix orders
  
  # number of observed dependent
  p = sum(!exog[1:nrow(S)])
  # number of observed independent
  q = nrow(S)-p
  
  # number of independent variables
  n = nrow(Phi)
  # number of dependent variables
  m = nrow(beta)
  
  # Number of observed variables
  r = p + q
  # Number of variables
  s = m + n
  
  # The matrices
  
  Beta <- matrix(0,s,s)
  Gamma <- matrix(0,s,n)
  G <- matrix(0,r,s)
  
  # Identity matrix
  
  I<-diag(nrow(Beta))
  
  # Populate the parameter matrices
  
  Gamma[1:m,] <- gamma
  Gamma[(m+1):s,] <- diag(n)
  
  Beta[1:m,1:m]<-beta
  
  # Populate the selection matrix (columns of all variables  selected to rows of observed variables)
  observed <- NULL
  if(p>0){
    observed <- 1:p
  }
  if(q>0){
    observed <- c(observed,(1:q)+n)
  }
  
  G <- diag(s)[observed,]
  
  # G has dependent variables followed by independent variables. Reorder to match the variables
  # in S
  
  G <- G[order(exog[1:nrow(S)]),]
  
  # Calculate the model implied covariance matrix and return
  
  Sigma = G %*% solve(I-Beta) %*% Gamma %*% Phi %*% t(Gamma) %*% t(solve(I-Beta))%*% t(G)
  rownames(Sigma) <- colnames(Sigma) <- rownames(S)
  
  Sigma
}

#'@title R2	for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{r2} computes the squared multiple correlation (R2)
#'for composites predicted by other composites in the model.
#'
#'@template postestimationFunctions
#'
#'@return A named numeric vector containing the R2 values.
#'
#'@export
#'

r2 <- function(object, ...){
  UseMethod("r2")
}

#'@export

r2.matrixpls <- function(object, ...){
  
  r2 <- rowSums(attr(object,"inner") * attr(object,"C"))
  names(r2) <- colnames(attr(object,"model")$inner)
  class(r2) <- "matrixplsr2"
  r2
}

#'@export

print.matrixplsr2 <- function(x, ...){
  cat("\n Inner model squared multiple correlations (R2)\n")
  print.table(x, ...)
}

#'@title Goodness of Fit indices for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{gof} computes the Goodness of Fit index for matrixpls results.
#'
#'@template postestimationFunctions
#'
#'@return The Goodness of Fit index.
#'
#'@references
#'
#'Henseler, J., & Sarstedt, M. (2013). Goodness-of-fit indices for partial least
#'squares path modeling. \emph{Computational Statistics}, 28(2), 565–580. 
#'\doi{10.1007/s00180-012-0317-1}
#'
#'
#'@export
#'

gof <- function(object, ...){
  UseMethod("gof")
}

#'@export

gof.matrixpls <- function(object, ...) {
  
  nativeModel <- attr(object,"model")
  IC <- attr(object,"IC")
  S <- attr(object,"S")
  
  exogenousLVs <- rowSums(nativeModel$inner) == 0
  
  
  IC_std <- IC %*% (diag(1/sqrt(diag(S))))
  
  result <- sqrt(mean(IC_std[t(nativeModel$reflective)==1]^2) * 
                   mean(r2(object)[!exogenousLVs]))
  
  # TODO: Relative gof
  # http://www.stat.wmich.edu/wang/561/codes/R/CanCor.R
  
  class(result) <- "matrixplsgof"
  result
}

#'@export

print.matrixplsgof <- function(x, digits=getOption("digits"), ...){
  cat("\n Absolute goodness of fit:", round(x, digits = digits))
  cat("\n")
}

#'@title Factor loadings matrix from matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{loadings} computes standardized factor loading matrix
#'from the model results.
#'
#'@template postestimationFunctions
#'
#'@return A matrix of estimated factor loadings.
#'
#'@importFrom stats loadings
#'
#'@export
#'
loadings <- function(object, ...){
  UseMethod("loadings")
}

#
#'@export

loadings.matrixpls <- function(object, ...) {

  attr(object,"reflective")

}

#'@title Composite Reliability indices for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{cr} computes Composite Reliability 
#'indices for the model using the formula presented by Fornell and Larcker (1981).
#'
#'@template postestimationFunctions
#'
#'@return A named numeric vector containing the Composite Reliability indices.
#'
#'@references
#'
#'Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of marketing research}, 18(1), 39–50.
#'
#'Aguirre-Urreta, M. I., Marakas, G. M., & Ellis, M. E. (2013). Measurement of 
#'composite reliability in research using partial least squares: Some issues and
#'an alternative approach. 
#' \emph{The DATA BASE for Advances in Information Systems}, 44(4), 11–43. 
#' \doi{10.1145/2544415.2544417}
#'
#'@export
#'

cr <- function(object, ...){
  UseMethod("cr")
}

#'@export

cr.matrixpls <- function(object, ...) {
  
  object <- standardize(object)
  
  loadings <- loadings(object)
  reflectiveModel <- attr(object, "model")$reflective
  
  result <- unlist(lapply(1:ncol(loadings), function(col){
    useLoadings <- abs(loadings[,col][reflectiveModel[,col]==1])
    
    crvalue <- sum(useLoadings)^2/
      (sum(useLoadings)^2 + sum(1-useLoadings^2))
    crvalue
  }))
  
  class(result) <- "matrixplscr"
  names(result) <- colnames(loadings)
  result
}

#'@export

print.matrixplscr <- function(x, ...){
  cat("\n Composite Reliability indices\n")
  print.table(x, ...)
}



#'@title Average Variance Extracted indices for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{ave} computes Average Variance Extracted 
#'indices for the model using the formula presented by Fornell and Larcker (1981).
#'
#'@template postestimationFunctions
#'
#'@return A list containing the Average Variance Extracted indices in the first position and the differences
#'between aves and largest squared correlations with other composites in the second position.
#'
#'
#'@references
#'
#'Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of marketing research}, 18(1), 39–50.
#'
#'@export
#'

ave <- function(object, ...){
  UseMethod("ave")
}


#'@export

ave.matrixpls <- function(object, ...) {
  
  object <- standardize(object)

  loadings <- loadings(object)
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
  
  result = list(ave = aves,
                ave_correlation = aves_correlation)
  
  class(result) <- "matrixplsave"
  result
}

#'@export

print.matrixplsave <- function(x, ...){
  cat("\n Average Variance Extracted indices\n")
  print.table(x$ave, ...)
  cat("\n AVE - largest squared correlation\n")
  print.table(x$ave_correlation, ...)
}

#'@title Heterotrait-monotrait ratio 
#'
#'@description 
#'
#'The \code{matrixpls} method for the generic function \code{htmt} computes Heterotrait-monotrait ratio 
#'for the model using the formula presented by Henseler et al (2014).
#'
#'@template postestimationFunctions
#'
#'@return Heterotrait-monotrait ratio as a scalar
#'
#'
#'@references
#'
#'Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for assessing discriminant validity in variance-based structural equation modeling. \emph{Journal of the Academy of Marketing Science}, 43(1), 115–135.
#'
#'
#'@export
#'

htmt <- function(object, ...){
  
  #HTMT is calculated from correlation matrix
  S <- stats::cov2cor(attr(object,"S"))
  
  reflective <- attr(object,"reflective")!=0
  
  # Mean within- and between-block correlations. Calculated as a sum of 
  # all correlations within block divided by the number of elements 
  
  excludeDiag <- 1-diag(nrow(S))
    
  meanBlockCor <- (t(reflective) %*% (S*excludeDiag) %*% reflective) /
    (t(reflective) %*% excludeDiag %*% reflective) 
  
  # Choose the blocks that have at least two reflective indicators
  i <- which(colSums(reflective) > 1)
  meanBlockCor <- meanBlockCor[i,i]

  if(length(i)<2){
    warning("HTMT can only be calculated if there are at least two composites with reflective parameters to at least two indicators.")
    return(NULL)
  }
  
  # Calculate the indices
  htmt <- meanBlockCor*lower.tri(meanBlockCor) /
    sqrt(diag(meanBlockCor) %o% diag(meanBlockCor))
  
  class(htmt) <- "matrixplshtmt"
  htmt
}

#'@export

print.matrixplshtmt <- function(x, ...){
  cat("\n Heterotrait-monotrait matrix\n")
  print.table(x, ...)
}

#'@title Composite Equivalence Indices 
#'
#'@description 
#'
#'The \code{matrixpls} method for the generic function \code{cei} computes 
#'composite equivalence indices (CEI) for the \code{matrixpls} object. By
#'default, the composites are compared against unit-weighted composites.
#'
#'@details
#'
#'Composite equivalence indices quantify if two sets of composites calculated
#'from the same data using different weight algorithms differ. Composites are 
#'matched by name and correlations for each pair are reported.
#'
#'@param object2 Another \code{matrixpls} object that \code{matrixpls} is
#'compared against.
#'
#'@template postestimationFunctions
#'
#'@return Composite equivalence indices as a vector
#'
#'@example example/matrixpls.cei-example.R
#'
#'@export
#'

cei <- function(object, object2 = NULL, ...){
  
  
  S <- attr(object,"S")
  W <- attr(object,"W")
  model <- attr(object,"model")
  
  if(is.null(object2)){
    object2 <- matrixpls(S, model, weightFun = weightFun.fixed)  
  }
  # Validate a user-provided object
  else{
    if(! identical(S, attr(object2,"S"))) stop("Both objects must use the same data covariance matrix S")
  }
  
  W2 <- attr(object2,"W")[rownames(W), colnames(W)]
  cei <- diag(W %*% S %*% t(W2))
  class(cei) <- "matrixplscei"
  cei
}

#'@export

print.matrixplscei <- function(x, digits=getOption("digits"), ...){
  cat("\n Composite equilevance indices\n\n CEI individual\n")
  print.table(x, digits = digits, ...)
  cat("\n CEI total:", round(min(x), digits = digits))
  cat("\n")
  
}


#'@title Summary of model fit of PLS model
#'
#'@description
#'
#'\code{fitSummary} computes a list of statistics
#'that are commonly used to access the overall fit of the PLS model.
#'
#'@template postestimationFunctions
#
#'@return A list containing a selection of fit indices.
#'
#'@export
#'

fitSummary <- function(object, ...){
  
  # Returns the minumum or NA if all are NA
  
  m <- function(x){
    x <- stats::na.exclude(x)
    if(length(x) == 0 ) NA
    else min(x)
  }
  
  ret <- c("Min CR" = m(cr(object)),
              "Min AVE" = m(ave(object)$ave),
              "Min AVE - sq. cor" = m(ave(object)$ave_correlation),
              "Goodness of Fit" = gof(object),
              stats::residuals(object)$indices[6:7],
              "CEI total" = cei(object)
           )
  
  class(ret) <- "matrixplfitsummary"
  ret
}

#'@export

print.matrixplfitsummary <- function(x, ...){
  cat("\n Summary indices of model fit\n")
  
  x <- as.matrix(x)
  colnames(x) <- "Value"
  print(x, ...)

}
