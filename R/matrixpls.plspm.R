# =========== Main functions ===========

#'@title A plspm compatibility wrapper for matrixpls
#'
#'@description
#'\code{matrixpls.plspm} mimics \code{plspm} function of the \code{plspm} package.
#'The arguments and their default values and the output of the function are identical with \code{plspm} function,
#'but internally the function uses matrixpls estimation.
#'
#'@details
#'The function \code{matrixpls.plspm} calculates indicator weights and estimates a model
#'identically to the  \code{plspm} function. In contrast to the \code{\link{matrixpls}} function
#'that provides only weights and parameter estimates, this function also reports multiple post-estimation
#'statistics. Because of this \code{matrixpls.plspm} is substantially less efficient than the \code{\link{matrixpls}}
#'function.
#'
#'The argument \code{path_matrix} is a matrix of zeros and ones that indicates
#'the structural relationships between composites. This must be a lower
#'triangular matrix. \code{path_matrix} will contain a 1 when column \code{j}
#'affects row \code{i}, 0 otherwise. \cr
#'
#' @param Data matrix or data frame containing the manifest variables.
#' @param path_matrix A square (lower triangular) boolean matrix representing 
#' the inner model (i.e. the path relationships between latent variables).
#' @param blocks list of vectors with column indices or column names
#' from \code{Data} indicating the sets of manifest variables forming 
#' each block (i.e. which manifest variables correspond to each block).
#' @param modes character vector indicating the type of measurement for each
#' block. Possible values are: \code{"A", "B"}. 
#' The length of \code{modes} must be equal to the length of \code{blocks}.
#' @param scheme string indicating the type of inner weighting
#' scheme. Possible values are \code{"centroid"}, \code{"factorial"}, or
#' \code{"path"}.
#' @param scaled whether manifest variables should be standardized. 
#' When (\code{TRUE}, data is 
#' scaled to standardized values (mean=0 and variance=1). 
#' The variance is calculated dividing by \code{N} instead of \code{N-1}).
#' @param tol decimal value indicating the tolerance criterion for the
#' iterations (\code{tol=0.000001}). Can be specified between 0 and 0.001.
#' @param maxiter integer indicating the maximum number of iterations
#' (\code{maxiter=100} by default). The minimum value of \code{maxiter} is 100.
#' @param boot.val whether bootstrap validation should be performed. 
#' (\code{FALSE} by default). 
#' @param br number bootstrap resamples. Used only
#' when \code{boot.val=TRUE}. When \code{boot.val=TRUE}, the default number of 
#' re-samples is 100.
#' @param dataset whether the data matrix used in the computations should be
#' retrieved (\code{TRUE} by default).
#' @return An object of class \code{"plspm"}. 
#' @return \item{outer_model}{Results of the outer model. Includes:
#' outer weights, standardized loadings, communalities, and redundancies}
#' @return \item{inner_model}{Results of the inner (structural) model. 
#' Includes: path coeffs and R-squared for each endogenous latent variable}
#' @return \item{scores}{Matrix of latent variables used to estimate the inner
#' model. If \code{scaled=FALSE} then \code{scores} are latent variables
#' calculated with the original data (non-standardized).}
#' @return \item{path_coefs}{Matrix of path coefficients 
#' (this matrix has a similar form as \code{path_matrix})}
#' @return \item{crossloadings}{Correlations between the latent variables 
#' and the manifest variables (also called crossloadings)}
#' @return \item{inner_summary}{Summarized results of the inner model. 
#' Includes: type of LV, type of measurement, number of indicators, R-squared,
#' average communality, average redundancy, and average variance
#' extracted}
#' @return \item{effects}{Path effects of the structural relationships. 
#' Includes: direct, indirect, and total effects}
#' @return \item{unidim}{Results for checking the unidimensionality of blocks
#' (These results are only meaningful for reflective blocks)}
#' @return \item{gof}{Goodness-of-Fit index}
#' @return \item{data}{Data matrix containing the manifest variables used in the
#' model. Only available when \code{dataset=TRUE}}
#' @return \item{boot}{List of bootstrapping results; only available 
#' when argument \code{boot.val=TRUE}}
#'
#'@references Sanchez, G. (2013). \emph{PLS Path Modeling with R.} Retrieved from http://www.gastonsanchez.com/PLS Path Modeling with R.pdf
#'#'
#'@export
#'@example example/matrixpls.plspm-example.R
#'



matrixpls.plspm <-
  function(Data, path_matrix, blocks, modes = NULL, scheme = "centroid", 
           scaled = TRUE, tol = 0.000001, maxiter = 100, boot.val = FALSE, 
           br = NULL, dataset = TRUE){
    
    if(! requireNamespace("boot")) stop("matrixpls.plspm requires the boot package")
    if(! requireNamespace("parallel")) stop("matrixpls.plspm requires the parallel package")
    

    
    # =======================================================
    # checking arguments
    # =======================================================
    params = get_params(x=Data, path_matrix, outer=blocks, modes=modes, 
                        scheme=scheme, scaled=scaled, tol=tol, maxiter=maxiter,
                        boot.val=boot.val, br=br, plsr=FALSE, dataset=dataset)
    
    
    #
    # Setup names etc
    #
    
    if(is.null(colnames(params$inner))) colnames(params$inner) <- rownames(params$inner)
    
    # 
    # Set up  matrixpls parameters
    #
    
    lvNames <- rownames(params$inner)
    reflective <- matrix(0, max(unlist(params$outer)),nrow(params$inner))
    rownames(reflective) <- colnames(params$x)[1:nrow(reflective)]
    colnames(reflective) <- lvNames
    
    col <- 1		
    
    for(outerSpec in params$outer){
      reflective[outerSpec,col] <- 1
      col <- col + 1 
    }
    
    # Only keep the variables that are used in the model
    
    varsToUse <- which(rowSums(reflective) > 0)
    reflective <- reflective[varsToUse,]
    dataToUse <- as.matrix(params$x[,varsToUse])
    
    # plspm does not support formative relationships, so we create a matrix of zero
    formative <- t(reflective)
    formative[] <- 0
    
    nativeModel <- list(inner=params$inner, 
                        reflective=reflective, 
                        formative=formative) 
    
    W.model <- t(reflective)
    
    if(params$plsr) paramsInner <- estimator.plsreg
    else paramsInner <- estimator.plsreg
    
    modeA <- params$modes == "A"
    
    if(max(modeA) == 0) outerEstim <- outerEstim.modeB	
    else if(min(modeA) == 1) outerEstim <- outerEstim.modeA
    else{
      outerEstim <- list(rep(NA,length(modeA)))
      outerEstim[!modeA] <- list(outerEstim.modeB)
      outerEstim[modeA] <- list(outerEstim.modeA)
    }
    
    
    innerEstim <- c(innerEstim.centroid,
                        innerEstim.factor,
                        innerEstim.path)[[params$scheme]]
    
    # Convergence is checked with the following function in plspm
    
    convCheck <- function(W,W_new){sum((abs(W_new) - abs(W))^2)} 
    
    #
    # Perform estimation
    # 
    
    if(scaled){
      dataToUse <- scale(dataToUse)
      S <- stats::cor(dataToUse, use = "complete.obs")
      S_std <- S
    }
    else{
      S <- stats::cov(dataToUse, use = "complete.obs")
      S_std <- stats::cov2cor(S)
    } 
    
    
    
    if(params$boot.val){
      
      #
      # The matrixpls.boot function is not used here because we need to support data scaling
      #
      
      boot.res <- boot::boot(dataToUse, function(originalData,indices){
        
        # plspm reports always standardized loadings even if the data were not standardized
        # We mimic this by standardizing the data
        
        if(scaled){
          S_boot <- stats::cor(originalData[indices,])
        }
        else{
          S_boot <- stats::cov(originalData[indices,])
        }
        
        tryCatch(
          matrixpls(S_boot, model = nativeModel, W.model = W.model, paramsInner = paramsInner,
                    outerEstim = outerEstim, innerEstim = innerEstim,
                    tol = params$tol, iter = params$maxiter, convCheck = convCheck,
                    validateInput = FALSE, standardize = FALSE), 
          
          error = function(e){
            print(e)
            print(S)
            print(nativeModel)
            print(W.model)
          })
      }, params$br)
      
      matrixpls.res <- boot.res$t0
    }
    else{
      matrixpls.res <- matrixpls(S, model = nativeModel, W.model = W.model, paramsInner = paramsInner,
                                 outerEstim = outerEstim, innerEstim = innerEstim,
                                 tol = params$tol, iter = params$maxiter, convCheck = convCheck,
                                 validateInput = FALSE, standardize = FALSE)
    }

    #
    # Compile a result object in the plspm format
    #
    
    C <- attr(matrixpls.res,"C") 
    IC <- attr(matrixpls.res,"IC")
    W <- attr(matrixpls.res,"W") 
    inner <- attr(matrixpls.res,"inner")
    
    IC_std <- IC %*% (diag(1/sqrt(diag(S))))
    colnames(IC_std) <- colnames(IC)
    
    # Inner model R2s
    R2 <- r2(matrixpls.res)
    class(R2) <- "numeric"
    
    # PLSPM does LV score scaling differently, so we need a correction
    sdv = sqrt(nrow(dataToUse)/(nrow(dataToUse)-1))   

    # Composite values, fittes values, and intercepts
    
    lvScores <- dataToUse %*% t(W) * sdv
    rownames(lvScores) <- 1:nrow(lvScores)
    
    lvScores_std <- apply(lvScores, 2, scale) * sdv
    rownames(lvScores_std) <- 1:nrow(lvScores_std)
    
    
    Modes <- ifelse(params$modes == "A", "Reflective","Formative")
    
    exogenousLVs <- rowSums(inner) == 0
    endogenousLVs <- rowSums(inner) != 0
    
    outer.mod <- sapply(lvNames, function(lvName){
      row <- which(lvNames == lvName)
      colnames <-colnames(Data)[params$outer[[row]]]
      weights <- W[row,colnames] * sdv
      std.loads <- IC_std[row,colnames]
      communal <- std.loads^2
      redundan <- communal * R2[row]
      
      ret <- cbind(weights, std.loads, communal, redundan)
      rownames(ret) <- colnames(Data)[params$outer[[row]]]
      return(ret)
    }, simplify= FALSE)

    outer.mod.dataframe <- do.call(rbind,sapply(lvNames,function(lvName){
      temp <- outer.mod[[lvName]]
      # as.factor is required to maintain compatibility with the orphaned plspm package
      data.frame(name = as.factor(rownames(temp)),
                 block = as.factor(lvName),
                 weight = temp[,1],
                 loading = temp[,2],
                 communality = temp[,3],
                 redundancy = temp[,4])
    }, simplify= FALSE))
    
    rownames(outer.mod.dataframe) <- NULL
    
    lvScoreFittedValues <- lvScores_std %*% t(inner)
    intercepts <- apply(lvScores_std-lvScoreFittedValues,2,mean)
    
    inner.mod <- sapply(lvNames[!exogenousLVs], function(lvName){
      
      # Recalculate the regressions to get standard errors and intercepts
      row <- which(lvNames == lvName)
      regressors <- inner[row,]!=0
      
      inner <- c(intercepts[lvName],inner[row,regressors]) 
      
      # Calculating of SEs adapted from mat.regress from the psych package
      
      if(sum(regressors) == 1){
        uniq <- 1
      }
      else{
        uniq <- (1-psych::smc(C[regressors,regressors]))
      }
      df <- nrow(params$x)-sum(regressors)-1
      se <- sqrt((1-R2[row])/(df))*c(1,sqrt(1/uniq))
      tvalue <- inner/se
      prob <- 2*(1- stats::pt(abs(tvalue),df))
      
      ret <- cbind(inner, se, tvalue, prob)
      
      rownames(ret) <- c("Intercept",lvNames[regressors])
      colnames(ret) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      
      return(ret)
    }, simplify= FALSE)
    
    # Indices of existing paths from the square model matrix
    pathIndices <- which(inner!=0)
    

    pathLabels <- paste(colnames(inner)[col(inner)[pathIndices]]," -> ",
                        rownames(inner)[row(inner)[pathIndices]], sep="")

    effectLabels <- paste(colnames(inner)[col(inner)[lower.tri(inner)]]," -> ",
                        rownames(inner)[row(inner)[lower.tri(inner)]], sep="")
    
    # A function for converting effects
    
    getEffects <- function(x){
      e <- inner
      e[endogenousLVs,] <- x
      e[lower.tri(e)]
    }
    
    if(params$boot.val){
      
      pathCount <- sum(nativeModel$inner)
      weightCount <- sum(W.model)
      bootPathIndices <- 1:pathCount
      bootLoadingIndices <- 1:weightCount + pathCount
      bootWeightIndices <- bootLoadingIndices + weightCount
      
      bootIndices <- boot::boot.array(boot.res, indices = TRUE)

      boot <- list(weights = get_bootDataFrame(boot.res$t0[bootWeightIndices]*sdv, boot.res$t[,bootWeightIndices]*sdv, rownames(S)),
                   loadings = get_bootDataFrame(IC_std[W.model == 1], 
                                                t(parallel::mcmapply(function(x){
                                                  
                                                  # Recover the data that were used in the bootrap replications and calculate indicator variances
                                                  sdX <- apply(dataToUse[bootIndices[x,],],2,stats::sd)
                                                  boot.res$t[x,bootLoadingIndices] / sdX									 														 	
                                                },1:params$br)),
                                                rownames(S)),
                   paths = get_bootDataFrame(boot.res$t0[bootPathIndices], boot.res$t[,bootPathIndices], 
                                             pathLabels),
                   
                   # Use parallel processing because of the large number of elements
                   
                   rsq = get_bootDataFrame(R2[!exogenousLVs],
                                           t(parallel::mcmapply(function(x){
                                             
                                             # Recover the data that were used in the bootrap replication
                                             if(scaled)
                                               S <- stats::cor(dataToUse[bootIndices[x,],])
                                             else
                                               S <- stats::cov(dataToUse[bootIndices[x,],])
                                             
                                             # Reconstruct the matrices
                                             W <- matrix(0,nrow(W.model),ncol(W.model))
                                             W[W.model==1] <- boot.res$t[x,bootWeightIndices]
                                             C <- W %*% S %*% t(W)
                                             inner <- matrix(0,nrow(nativeModel$inner),ncol(nativeModel$inner))
                                             inner[nativeModel$inner == 1] <- boot.res$t[x,bootPathIndices]
                                             
                                             # Calculate R2 values
                                             rowSums(inner[!exogenousLVs,]*C[!exogenousLVs,])
                                             
                                           },1:params$br)),
                                           lvNames[!exogenousLVs]),
                   
                   # Again, parellel processing
                   
                   total.efs = get_bootDataFrame(getEffects(effects(matrixpls.res)$Total),
                                                 t(parallel::mcmapply(function(x){
                                                   
                                                   inner <- matrix(0,nrow(nativeModel$inner),ncol(nativeModel$inner))
                                                   inner[nativeModel$inner == 1] <- boot.res$t[x,bootPathIndices]

                                                   # Create a fake matrixpls object to calculate total effects
                                                   
                                                   obj <- numeric()
                                                   class(obj) <- "matrixpls"
                                                   attr(obj,"inner") <- inner
                                                   attr(obj,"model") <- nativeModel
                                                   getEffects(effects.matrixpls(obj)$Total)
                                                   
                                                 },1:params$br)),
                                                 effectLabels))
    }
    else{
      boot = FALSE
    }
    
    blocks <- unlist(lapply(params$outer, length))
    names(blocks) <- lvNames
    
    blocksList <- params$outer
    names(blocksList) <- lvNames
    
    mvsToUse <- unlist(params$outer)
    
    model <- list(IDM=params$inner,
                  blocks= blocksList,
                  specs = list(scaling = NULL,
                               modes = modes, 
                               scheme=switch(scheme, "centroid"="centroid", "factor"="factor", "path"="path"),
                               scaled=params$scaled, 
                               tol=params$tol, 
                               maxiter = params$maxiter,
                               plscomp = NULL),
                  iter = attr(matrixpls.res,"iterations") + 1,
                  boot.val=params$boot.val, 
                  br=params$br, 
                  gens = list(obs=nrow(params$x),
                              obs_names=rownames(params$x),
                              mvs=length(mvsToUse),
                              mvs_names = colnames(params$x)[mvsToUse],
                              lvs = length(lvNames),
                              lvs_names = lvNames))
    
    
    if(params$dataset){
      data <- as.matrix(params$x[,unlist(params$outer)])
    }
    else{
      data <- FALSE
    }
    
    Windices <-W!=0
    out.weights <- W[Windices] * sdv
    names(out.weights) <- colnames(W)[col(W)[Windices]]
    
    Wlist <- apply(Windices, 1, which)

    # Weight and loading pattern is always identical in plspm
    loadings <- IC_std[Windices]
    names(loadings) <- names(out.weights)
    
    outer.cor <- sapply(lvNames, function(lvName){
      row <- which(lvNames == lvName)
      vars <- Wlist[[row]]
      return(t(IC_std[,vars]))
    }, simplify= FALSE)
    
    
    inner.sum <- data.frame(Type = ifelse(exogenousLVs,"Exogenous","Endogenous"),
                            R2 = R2, 
                            Block_Communality = sapply(outer.mod,function(x){mean(x[,"communal"])}), 
                            Mean_Redundancy = sapply(outer.mod,function(x){mean(x[,"redundan"])}), 
                            AVE = ifelse(params$modes == "A",
                                         sapply(outer.mod,function(x){mean(x[,"std.loads"]^2)}),
                                         0))
    
    effs <- effects(matrixpls.res)
    
    effects <- data.frame(relationships = effectLabels,
                          direct = getEffects(effs$Direct),
                          indirect = getEffects(effs$Indirect), 
                          total = getEffects(effs$Total))

    unidim <- data.frame(Mode = params$modes,
                         MVs = blocks,
                         C.alpha = ifelse(params$modes == "A",
                                          sapply(Wlist,function(indices){
                                            if(length(indices) == 1) return(1)
                                            Ssub <- stats::cov2cor(S[indices,indices])
                                            alpha(Ssub)
                                          }, simplify = TRUE),
                                          0),
                         DG.rho = ifelse(params$modes == "A",
                                         sapply(Wlist,function(indices){
                                           if(length(indices) == 1) return(1)
                                           pc <- psych::principal(S[indices,indices])
                                           std.loads <- pc$loadings
                                           numer.rho <- sum(std.loads)^2
                                           denom.rho <- numer.rho + (length(indices) - sum(std.loads^2))
                                           numer.rho / denom.rho
                                         }, simplify = TRUE),
                                         0),
                         
                         eig.1st = sapply(Wlist,function(indices){
                           if(length(indices) == 1) return(1)
                           eigen(S_std[indices,indices])$values[1]
                         }, simplify = TRUE),   
                         eig.2nd = sapply(Wlist,function(indices){
                           if(length(indices) == 1) return(0)
                           eigen(S_std[indices,indices])$values[2]
                         }, simplify = TRUE))
    
    # Goodness of Fit is square root of product of mean communality and mean R2
    # plspm calulates this a bit differently

    gof <- gof(matrixpls.res)
    class(gof) <- "numeric"
    
    # Crossloadings are the IC matrix

    crossloadings <- cbind(outer.mod.dataframe[,1:2],
                           t(IC[,as.character(outer.mod.dataframe$name)])/
                             sqrt(diag(S))[as.character(outer.mod.dataframe$name)])

    rownames(crossloadings) <- NULL
    
    data <- as.data.frame(data)
    attr(data,"row.names") <- as.character(attr(data,"row.names"))
    
    manifests <- scale(data, scale=FALSE)

    res = list(outer_model = outer.mod.dataframe, 
               inner_model = inner.mod, 
               path_coefs = inner, 
               scores = lvScores_std,
               crossloadings = crossloadings,
               inner_summary = inner.sum, 
               effects = effects,
               unidim = unidim, 
               gof = gof, 
               boot = boot, 
               data = data, 
               manifests = manifests,
               model = model)
    
    #		out.weights = out.weights, 
    #		loadings = loadings, 
    #		r.sqr = R2,
    #		outer.cor = outer.cor, 
    #							 latents = lvScores_std, 
    
    class(res) = "plspm"
    
    return(res)
    
  }


#
# Validate parameters. Copied from PLSPM
#

get_params <-
  function(x, inner, outer, modes=NULL, scheme="centroid", scaled=TRUE,
           boot.val=FALSE, br=NULL, plsr=FALSE, tol=0.00001, maxiter=100, dataset=TRUE)
  {
    # =======================================================
    # checking arguments
    # =======================================================
    if (!is.matrix(x) && !is.data.frame(x))
      stop("Invalid object 'x'. Must be a numeric matrix or data frame.")
    if (is.null(rownames(x))) 
      rownames(x) <- 1:nrow(x)
    if (is.null(colnames(x))) 
      colnames(x) <- paste("MV", 1:ncol(x), sep="")
    if (!is.matrix(inner))
      stop("Invalid argument 'inner'. Must be a matrix.")
    if (nrow(inner)!=ncol(inner))
      stop("Invalid argument 'inner'. Must be a square matrix.")
    for (j in 1:ncol(inner))
      for (i in 1:nrow(inner)) {
        if (i<=j)
          if (inner[i,j]!=0) 
            stop("argument 'inner' must be a lower triangular matrix")
        if (length(intersect(inner[i,j], c(1,0)))==0)
          stop("elements in 'inner' must be '1' or '0'")
      }
    if (is.null(dimnames(inner)))
      lvs.names <- paste("LV", 1:ncol(inner), sep="")
    if (!is.null(rownames(inner)))
      lvs.names <- rownames(inner)
    if (!is.null(colnames(inner)))
      lvs.names <- colnames(inner)
    if (!is.list(outer))
      stop("Invalid argument 'outer'. Must be a list.")
    if (length(outer) != nrow(inner))
      stop("Number of rows of 'inner' does not coincide with length of 'outer'.")
    if (is.null(modes)) {
      modes <- rep("A",length(outer))
      warning("Argument 'modes' missing. Default reflective 'modes' is used.")
    }
    if (length(outer) != length(modes)) {
      warning("Warning: Invalid length of 'modes'. Default reflective 'modes' is used.")
      modes <- rep("A", length(outer))
    }
    for (i in 1:length(modes))
      if (modes[i]!="A" && modes[i]!="B") modes[i]<-"A"
    if (!is.na(pmatch(scheme, "centroid"))) 
      scheme <- "centroid"
    SCHEMES <- c("centroid", "factor", "path")
    scheme <- pmatch(scheme, SCHEMES)
    if (is.na(scheme)) {
      warning("Warning: Invalid argument 'scheme'. Default 'scheme=centroid' is used.")   
      scheme <- "centroid"
    }
    if (!is.logical(scaled)) {
      warning("Warning: Invalid argument 'scaled'. Default 'scaled=TRUE' is used.")
      scaled <- TRUE
    }
    if (!is.logical(boot.val)) {
      warning("Warning: Invalid argument 'boot.val'. No bootstrap validation is done.")
      boot.val <- FALSE
    }   
    if (boot.val) {
      if (!is.null(br)) {        
        if (mode(br)!="numeric" || length(br)!=1 || (br%%1)!=0 ||
              br<100 || br>1000) {
          warning("Warning: Invalid argument 'br'. Default 'br=100' is used.")   
          br <- 100
        } 
      } else
        br <- 100
    }
    if (!is.logical(plsr)) plsr<-FALSE
    if (mode(tol)!="numeric" || length(tol)!=1 || tol<=0 || tol>0.001) {
      warning("Warning: Invalid argument 'tol'. Default 'tol=0.00001' is used.")   
      tol <- 0.00001
    } 
    if (mode(maxiter)!="numeric" || length(maxiter)!=1 || maxiter<100) {
      warning("Warning: Invalid argument 'maxiter'. Default 'maxiter=100' is used.")   
      maxiter <- 100
    } 
    if (!is.logical(dataset)) 
      dataset <- TRUE
    
    # =======================================================
    # list with verified arguments
    # =======================================================
    res = list(
      x = x,
      inner = inner,
      outer = outer,
      modes = modes,
      scheme = scheme,
      SCHEMES = SCHEMES,
      scaled = scaled,
      boot.val = boot.val,
      br = br,
      plsr = plsr,
      tol = tol,
      maxiter = maxiter,
      dataset = dataset)
    res
  }

get_bootDataFrame <- function(t0, t, rownames){
  ret <- data.frame(Original = t0, 
                    Mean.Boot = apply(t, 2, mean), 
                    Std.Error = apply(t, 2, stats::sd), 
                    perc.025 = apply(t, 2, stats::quantile, 0.025),
                    perc.975 = apply(t, 2, stats::quantile, 0.975))
  
  rownames(ret) <- rownames
  
  return(ret)
}

