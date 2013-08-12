
#
#'@title A semPLS compatibility wrapper for MatrixPLS
#'
#'@description
#'  \code{matrixpls.sempls} fits structural equation models by the patial least
#'  squares (PLS) method. The estimation is based on the raw data and
#'  requires no distributional assumptions.
#'
#'@param model 
#'    An object inheriting from class \code[sempls]{plsm} as returned from
#'    \code{\link[sempls]{plsm}} or \code{\link[sempls]{read.splsm}}.
#' 
#'@param data 
#'    A \code{data.frame} containing the observed variables
#'    (MVs). The storage mode for all the MVs included in the model must
#'    be \code{numeric}.
#' 
#'@param maxit 
#'    A \code{numeric} value, which determines the maximum number of
#'    iterations performed by the PLS algorithm. The default is \eqn{20}
#'    iterations.
#' 
#'@param tol 
#'    A \code{numeric} value, specifying the tolerance for the maximum relative
#'    differences in the outer weights. The default value is
#'    \eqn{10^{-7}}.
#' 
#'@param scaled 
#'    A \code{logical} value indicating, whether the observed
#'    data shall be scaled to zero mean and unit variance. The default is
#'    \code{TRUE}.
#' 
#'@param sum1 
#'    A \code{logical} value indicating, whether the outer
#'    weights foreach latent variable (LV) shall be standardized to sum up
#'    to one. The default is \code{FALSE}. Since the factor scores are
#'    scaled in each step of the PLS algorithm, changing this value to
#'    \code{TRUE} does not affect the results.
#' 
#'
#'@param wscheme A \code{character} naming the weighting scheme to
#'    use. Possible values are:
#'    \itemize{
#'      \item \code{"A"} or \code{"centroid"} for the centroid scheme, the default,
#'      \item \code{"B"} or \code{"factorial"}for the factorial scheme and
#'      \item \code{"C"}, \code{"pw"} or \code{"pathWeighting"} for the path weighting scheme.
#'    }
#'  
#'
#'@param pairwise A \code{logical} value indicating, whether
#'    correlations shall be calculated pairwise. If the observed data
#'    does not contain missing values, the results are not affected.
#'    The default is \code{FALSE}. For more details the R help,
#'    \code{?cor}, can be consulted.}
#'    
#'@param method A \code{character} naming the method to calculate
#'    the correlations. Possible values are:
#'     \itemize{
#'       \item \code{"pearson"} , the default,
#'       \item \code{"kendall"},
#'       \item \code{"spearman"}.
#'     }
#'    For more details on the method, the R help, \code{?cor}, can be
#'    consulted. Note, that despite of the \code{method} argument, pearson
#'    correlations are always used for the inner approximation (step 2).
#' 
#'@param convCrit
#'    The convergence criteria to use:
#'    \itemize{
#'      \item \code{"relative"}, the default,
#'      \item \code{"square"}.
#'    }
#'  
#'@param verbose Logical: If \code{FALSE} no status messages are printed.
#'  
#'@param object
#'    An object of class \code[sempls]{sempls}.
#'    
#'@param x
#'    An object of the according class.
#'  
#'@param type
#'    If the argument \code{what="loadings"}, \code{type} describes the
#'    loadings to be extracted -- those are:
#'    \itemize{
#'      \item \code{"discriminant"}, the default, contrasts outer against cross
#'      loadings to check for discrimant validity of the measurement model,
#'      \item \code{"outer"} for the outer loadings and
#'      \item \code{"cross"} for the cross loadings.
#'    }
#'  
#'@param cutoff
#'    A numerical value at which to cutoff the loadings -- this means
#'    loadings smaller than the cutoff value will not be printed.
#'  
#'@param reldiff
#'    The argument is only effectiv when \code{type="discriminant"}. It is
#'    a \code{numeric} value, specifying the relative difference between
#'    outer and cross loadings at which cross loadings will still be
#'    printed.
#'  
#'@param na.print
#'    A \code{character} substituting values not to be printed.
#'  
#'@param digits
#'    minimal number of _significant_ digits, see \code{\link{print.default}}.
#' 
#'@param use
#'    The values for which the density plots are created. If
#'    \itemize{
#'      \item \code{"fscores"}: the factor scores are used,
#'      \item \code{"prediction"}: the estimated factor scores are used,
#'      \item \code{"residuals"}: the residuals are used.
#'    }
#'  
#'@param abbreviate
#'    A logical indicating whether dimnames should be abbreviated. For
#'    Details see \code{\link{abbreviate}}. The default is \code{FALSE}.
#'  
#'
#'
#'@return
#'  \code{matrixpls.sempls} returns an object of class \code[sempls]{sempls}, with the following elements:
#'
#'  \item{coefficients}{
#'    A \code{data.frame} containing the estimates for
#'    all the arcs in the path model, those are the outer loadings for
#'    mode \sQuote{A} type LVs and outer weights for mode \sQuote{B} type LVs and path
#'    coefficients for those belonging to the structural model.
#'  }
#'  \item{path_coefficient}{
#'    The \code{matrix} of path coefficients.
#'  }
#'  \item{outer_loadings}{
#'    The \code{matrix} of outer loadings.
#'  }
#'  \item{cross_loadings}{
#'    The \code{matrix} of cross loadings.
#'  }
#'  \item{total_effects}{
#'    The \code{matrix} of total effects.
#'  }
#'  \item{inner_weights}{
#'    The \code{matrix} of inner weights.
#'  }
#'  \item{outer_weights}{
#'    The \code{matrix} of outer weights.
#'  }
#'  \item{factor_scores}{
#'    A \code{data.frame} containing the estimated factor scores for the
#'    LVs.
#'  }
#'  \item{data}{
#'    A \code{data.frame} containing the preprocessed obseravtions of the
#'    MVs.
#'  }
#'  \item{incomplete}{
#'    The index of the incomplete observations.
#'  }
#'  \item{...}{
#'    All the other values are just storing information used in
#'    the \code{call}.
#'  }
#'}
#'
#'@author Mikko Rönkkö
#'@author Armin Monecke
#'
#'@seealso
#'   \code{\link[sempls]{plsm}}, \code{\link[sempls]{read.splsm}},
#'   \code{\link[sempls]{rSquared}}, \code{\link[sempls]{pathDiagram}},
#'   \code{\link[sempls]{bootsempls}}, \code{\link[sempls]{plsm2sem}},
#'   \code{\link[sem]{sem}}
#'
#'
#'@examples
#'\dontrun{
#'data(ECSImobi)
#'ecsi <- sempls(model=ECSImobi, data=mobi, wscheme="pathWeighting")
#'ecsi
#'
#'## create plots
#'densityplot(ecsi)
#'densityplot(ecsi, use="prediction")
#'densityplot(ecsi, use="residuals")
#'
#'## Values of 'sempls' objects
#'names(ecsi)
#'ecsi$outer_weights
#'ecsi$outer_loadings
#'ecsi$path_coefficients
#'ecsi$total_effects
#'
#'
#'### using convenience methods to sempls results
#'## path coefficients
#'pathCoeff(ecsi)
#'
#'## total effects
#'totalEffects(ecsi)
#'
#'## get loadings and check for discriminant validity
#'(l <- plsLoadings(ecsi))
#'# outer loadings
#'print(l, type="outer", digits=2)
#'# outer loadings greater than 0.5
#'print(l,type="outer", cutoff=0.5, digits=2)
#'# cross loadings greater than 0.5
#'print(l, type="cross", cutoff=0.5, digits=2)
#'
#'
#'### R-squared
#'rSquared(ecsi)
#'
#'
#'### Create .dot representation of the path diagram and
#'### create .pdf file if graphviz is available.
#'
#'pathDiagram(ecsi, file="ecsiPLS1", edge.labels="both",
#'            output.type="graphics", digits=3, graphics.fmt = "pdf")
#'
#'# include R-squared values
#'pathDiagram(ecsi, file="ecsiPLS2", edge.labels="both",
#'            output.type="graphics", digits=3, graphics.fmt = "pdf",
#'            rSquared=rSquared(ecsi))
#'
#'# only the structural model
#'pathDiagram(ecsi, file="ecsiPLS3", edge.labels="both",
#'            output.type="graphics", digits=3, graphics.fmt = "pdf",
#'            rSquared=rSquared(ecsi), full=FALSE)
#'
#'}
#'

matrixpls.sempls <-
function(model, data, maxit=20, tol=1e-7, scaled=TRUE, sum1=FALSE, wscheme="centroid", pairwise=FALSE,
         method=c("pearson", "kendall", "spearman"),
         convCrit=c("relative", "square"), verbose=TRUE){
         
         stop("Not implemented")
  method <- match.arg(method)
  convCrit <- match.arg(convCrit)
  result <- list(coefficients=NULL, path_coefficients=NULL,
                 outer_loadings=NULL ,cross_loadings=NULL,
                 total_effects=NULL,inner_weights=NULL, outer_weights=NULL,
                 blocks=NULL, factor_scores=NULL, data=NULL, scaled=scaled,
                 model=model, weighting_scheme=NULL, weights_evolution=NULL,
                 sum1=sum1, pairwise=pairwise, method=method, iterations=NULL,
                 convCrit=convCrit, verbose=verbose, tolerance=tol, maxit=maxit, N=NULL,
                 incomplete=NULL, Hanafi=NULL)
  class(result) <- "sempls"

  # checking the data
  data <- data[, model$manifest]
  N <- nrow(data)
  missings <- which(complete.cases(data)==FALSE)
  if(length(missings)==0 & verbose){
    cat("All", N ,"observations are valid.\n")
    if(pairwise){
      pairwise <- FALSE
      cat("Argument 'pairwise' is reset to FALSE.\n")
    }
  }
  else if(length(missings)!=0 & !pairwise & verbose){
    # Just keeping the observations, that are complete.
    data <- na.omit(data[, model$manifest])
    cat("Data rows:", paste(missings, collapse=", "),
        "\nare not taken into acount, due to missings in the manifest variables.\n",
        "Total number of complete cases:", N-length(missings), "\n")
  }
  else if(verbose){
     cat("Data rows", paste(missings, collapse=", "),
         " contain missing values.\n",
         "Total number of complete cases:", N-length(missings), "\n")
  }
  ## check the variances of the data
  if(!all(apply(data, 2, sd, na.rm=TRUE) != 0)){
     stop("The MVs: ",
          paste(colnames(data)[which(apply(data, 2, sd)==0)], collapse=", "),
          "\n  have standard deviation equal to 0.\n",
          "  Recheck model!\n")
  }

  ## scale data?
  # Note: scale() changes class(data) to 'matrix'
  if(scaled) data <- scale(data)

  ## compute PLS approximation of LV scores

  #############################################
  # step 1: Initialisation
  stp1 <- step1(model, data, sum1=sum1, pairwise, method)
  factor_scores <- stp1$latent
  if(!sum1){
      # to ensure: w'Sw=1
      sdYs <- rep(attr(factor_scores, "scaled:scale"),
                  each=length(model$manifest))
      Wold <- stp1$outerW / sdYs
  }
  else Wold <- stp1$outerW
  weights_evolution <- reshape(as.data.frame(Wold),
                               v.names="weights",
                               ids=rownames(Wold),
                               idvar="MVs",
                               times=colnames(Wold),
                               timevar="LVs",
                               varying=list(colnames(Wold)),
                               direction="long")
  ## fixme: find something more efficient than 'cbind'
  weights_evolution <- cbind(weights_evolution, iteration=0)
  Hanafi <- cbind(f=sum(abs(cor(factor_scores)) * model$D),
                  g=sum(cor(factor_scores)^2 * model$D),
                  iteration=0)


  #############################################
  # Select the function according to the weighting scheme
  if(wscheme %in% c("A", "centroid")) {
    innerWe <- centroid
    result$weighting_scheme <- "centroid"
  }
  else if(wscheme %in% c("B", "factorial")) {
    innerWe <- factorial
    result$weighting_scheme <- "factorial"
  }
  else if(wscheme %in% c("C", "pw", "pathWeighting")) {
    innerWe <- pathWeighting
    result$weighting_scheme <- "path weighting"
  }
  else {stop("The argument E can only take the values 'A', 'B' or 'C'.\n See ?sempls")}

  converged <- c()
  i <- c()
  Wnew <- c()
  innerWeights <- c()
  eval(plsLoop)

  ## print
  if(converged & verbose){
      cat(paste("Converged after ", (i-1), " iterations.\n",
                "Tolerance: ", tol ,"\n", sep=""))
      if (wscheme %in% c("A", "centroid")) cat("Scheme: centroid\n")
      if (wscheme %in% c("B", "factorial")) cat("Scheme: factorial\n")
      if (wscheme %in% c("C", "pw", "pathWeighting")) cat("Scheme: path weighting\n")
  }
  else if(!converged){
      stop("Result did not converge after ", result$maxit, " iterations.\n",
           "\nIncrease 'maxit' and rerun.", sep="")
  }

  weights_evolution <- weights_evolution[weights_evolution!=0,]
  weights_evolution$LVs <- factor(weights_evolution$LVs,  levels=model$latent)
  # create result list
  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  result$path_coefficients <- pathCoeff(model=model, factor_scores, method, pairwise)
  result$cross_loadings <- cor(data, factor_scores, use, method)
  result$outer_loadings <- result$cross_loadings
  result$outer_loadings[Wnew==0] <- 0
  result$total_effects <- totalEffects(result$path_coefficients)
  result$inner_weights <- innerWeights
  result$outer_weights <- Wnew
  result$weights_evolution <- weights_evolution
  result$Hanafi <- Hanafi
  result$factor_scores <- factor_scores
  result$data <- data
  result$N <- N
  result$incomplete <- missings
  result$iterations <- (i-1)
  result$coefficients <- coefficients(result)
  return(result)
}