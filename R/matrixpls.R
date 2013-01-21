#'@title PLS-PM: Partial Least Squares Path Modeling
#'
#'@description
#'Estimate path models with latent variables by partial least squares approach
#'
#'@details
#'The function \code{matrixpls} estimates a path model by partial least squares
#'approach providing the full set of results. \cr
#'
#'The argument \code{inner_matrix} is a matrix of zeros and ones that indicates
#'the structural relationships between latent variables. This must be a lower
#'triangular matrix. \code{inner_matrix} will contain a 1 when column \code{j}
#'affects row \code{i}, 0 otherwise. \cr
#'
#'@param Data A numeric matrix or data frame containing the manifest variables.
#'@param inner_matrix A square (lower triangular) boolean matrix representing the
#'inner model (i.e. the path relationships betwenn latent variables).
#'@param outer_list List of vectors with column indices from \code{x} indicating 
#'the sets of manifest variables asociated to the latent variables
#'(i.e. which manifest variables correspond to the latent variables).
#'Length of \code{outer_list} must be equal to the number of rows of \code{inner_matrix}.
#'@param modes A character vector indicating the type of measurement for each
#'latent variable. \code{"A"} for reflective measurement or \code{"B"} for
#'formative measurement (\code{NULL} by default). The length of \code{modes}
#'must be equal to the length of \code{outer_list}).
#'@param scheme A string of characters indicating the type of inner weighting
#'scheme. Possible values are \code{"centroid"}, \code{"factor"}, or
#'\code{"path"}.
#'@param scaled A logical value indicating whether scaling data is performed
#'When (\code{TRUE} data is scaled to standardized values (mean=0 and variance=1)
#'The variance is calculated dividing by \code{N} instead of \code{N-1}).
#'@param tol Decimal value indicating the tolerance criterion for the
#'iterations (\code{tol=0.00001}). Can be specified between 0 and 0.001.
#'@param iter An integer indicating the maximum number of iterations
#'(\code{iter=100} by default). The minimum value of \code{iter} is 100.
#'@param boot.val A logical value indicating whether bootstrap validation is
#'performed (\code{FALSE} by default). 
#'@param br An integer indicating the number bootstrap resamples. Used only
#'when \code{boot.val=TRUE}. When \code{boot.val=TRUE}, the default number of 
#'re-samples is 100, but it can be specified in a range from 100 to 1000.
#'@param plsr A logical value indicating whether pls regression is applied
#'to calculate path coefficients (\code{FALSE} by default).
#'@param dataset A logical value indicating whether the data matrix should be
#'retrieved (\code{TRUE} by default).
#'@return An object of class \code{"matrixpls"}. 
#'@return \item{outer.mod}{Results of the outer (measurement) model. Includes:
#'outer weights, standardized loadings, communalities, and redundancies}
#'@return \item{inner.mod}{Results of the inner (structural) model. Includes: path
#'coefficients and R-squared for each endogenous latent variable}
#'@return \item{latents}{Matrix of standardized latent variables (variance=1
#'calculated divided by \code{N}) obtained from centered data (mean=0)}
#'@return \item{scores}{Matrix of latent variables used to estimate the inner
#'model. If \code{scaled=FALSE} then \code{scores} are latent variables
#'calculated with the original data (non-stardardized). If \code{scaled=TRUE}
#'then \code{scores} and \code{latents} have the same values}
#'@return \item{out.weights}{Vector of outer weights}
#'@return \item{loadings}{Vector of standardized loadings (i.e. correlations with
#'LVs)}
#'@return \item{path.coefs}{Matrix of path coefficients (this matrix has a similar
#'form as \code{inner_matrix})}
#'@return \item{r.sqr}{Vector of R-squared coefficients}
#'@return \item{outer.cor}{Correlations between the latent variables and the
#'manifest variables (also called crossloadings)}
#'@return \item{inner.sum}{Summarized results by latent variable of the inner
#'model. Includes: type of LV, type of measurement, number of indicators,
#'R-squared, average communality, average redundancy, and average variance
#'extracted}
#'@return \item{effects}{Path effects of the structural relationships. Includes:
#'direct, indirect, and total effects}
#'@return \item{unidim}{Results for checking the unidimensionality of blocks
#'(These results are only meaningful for reflective blocks)}
#'@return \item{gof}{Goodness-of-Fit index}
#'@return \item{data}{Data matrix containing the manifest variables used in the
#'model. Only when \code{dataset=TRUE}}
#'@return \item{boot}{List of bootstrapping results; only available when argument
#'\code{boot.val=TRUE}}
#'@author Gaston Sanchez
#'
#'@references Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro C.
#'(2005) PLS path modeling. \emph{Computational Statistics & Data Analysis},
#'\bold{48}, pp. 159-205.
#'
#'Tenenhaus M., Pages J. (2001) Multiple factor analysis combined with
#'PLS path modelling. Application to the analysis of relationships between
#'physicochemical variables, sensory profiles and hedonic judgements.
#'\emph{Chemometrics and Intelligent Laboratory Systems}, \bold{58}, pp.
#'261-273.
#'
#'Tenenhaus M., Hanafi M. (2010) A bridge between PLS path modeling and
#'multi-block data analysis. \emph{Handbook on Partial Least Squares (PLS):
#'Concepts, methods, and applications.} Springer.
#'
#'Lohmoller J.-B. (1989) \emph{Latent variables path modelin with partial
#'least squares.} Heidelberg: Physica-Verlag.
#'
#'Wold H. (1985) Partial Least Squares. In: Kotz, S., Johnson, N.L. (Eds.),
#'\emph{Encyclopedia of Statistical Sciences}, Vol. 6. Wiley, New York, pp.
#'581-591.
#'
#'Wold H. (1982) Soft modeling: the basic design and some extensions. In: K.G.
#'Joreskog & H. Wold (Eds.), \emph{Systems under indirect observations:
#'Causality, structure, prediction}, Part 2, pp. 1-54. Amsterdam: Holland.
#'@seealso \code{\link{matrixpls.fit}}, \code{\link{plot.matrixpls}}
#'@export
#'@examples
#'
#'  \dontrun{
#'  ## typical example of PLS-PM in customer satisfaction analysis
#'  ## model with six LVs and reflective indicators
#'
#'  # load dataset satisfaction
#'  data(satisfaction)
#'
#'  # inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0) 
#'  LOY = c(1,0,0,0,1,0)
#'  sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'
#'  # outer model list
#'  sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'
#'  # vector of modes (reflective indicators)
#'  sat_mod = rep("A", 6)
#'
#'  # apply matrixpls
#'  satpls = matrixpls(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)
#'  
#'  # summary of results
#'  summary(satpls)
#'
#'  # default plot (inner model)
#'  plot(satpls)
#'  }
#'
matrixpls <-
function(Data, inner_matrix, outer_list, modes = NULL, scheme = "centroid", 
         scaled = TRUE, tol = 0.00001, iter = 100, boot.val = FALSE, 
         br = NULL, plsr = FALSE, dataset = TRUE, matrixData = FALSE)
{
  # =======================================================
  # checking arguments
  # =======================================================
  
  valid = get_params(x=Data, inner=inner_matrix, outer=outer_list, modes=modes, 
                     scheme=scheme, scaled=scaled, tol=tol, iter=iter,
                     boot.val=boot.val, br=br, plsr=plsr, dataset=dataset)
  x = valid$x
  inner = valid$inner
  outer = valid$outer
  modes = valid$modes
  scheme = valid$scheme
  SCHEMES = valid$SCHEMES
  scaled = valid$scaled
  boot.val = valid$boot.val
  br = valid$br
  plsr = valid$plsr
  tol = valid$tol
  iter = valid$iter
  dataset = valid$dataset
  
  # Stop if there are unimplemented options
  if(sum(modes=="B")>0) stop("Mode B is not currently implemented in matrixpls")
  if(scheme!="centroid") stop("Only centroid scheme is currently implemented in matrixpls")
  if(plsr) stop("PLS regression is currently not implemented in matrixpls")
  

  if(matrixData){
  	indicatorCovariances <- x
  }
  else{
  	indicatorCovariances <- cov(x)
  }
  indicatorCorrelations <- cov2cor(indicatorCovariances)
  
  doBootstrap = boot.val
  if (doBootstrap && nrow(Data) <= 10) {
      warning("Bootstrapping stopped: very few cases.") 
      doBootstrap<-FALSE
  }
  
  #
  # Boot will provide also estimates for the full data
  #
  
  if(doBootstrap){
  	if(matrixData) stop("Bootstrapping requires raw data")
  	library(boot)
  	boot.results<-boot(x,function(d,p) {
  		
  		cov.matrix<-cov(d[p,],d[p,])
  		
  		pls<-matrixpls.fit(cov.matrix,inner,outer,modes,scheme,scaled,tol,iter)
  		
  		# Selection matrix
  		t<-lower.tri(inner)

  		if(is.null(pls)){
  			# if the method failed to converge, return a vector of NAs with appropriate lenght
  			rep(NA,sum(t)*2+nrow(cov.matrix)*2+nrow(inner))
  		}
  		else{
	  		c(pls$loading,pls$out.weights,pls$r.sqr,pls$path.coef[t]pls$path.effects[t])
	  	}
  	},br)
  	
  	# Unpack the results from t0
	boot.loadIndices <- 1:ncol(Data)
	boot.weightIndices <- boot.loadIndices+ncol(Data)
	boot.R2Indices <- (tail(boot.weightIndices,1)+1):(tail(boot.weightIndices,1)+ncol(inner))
	boot.pathIndices <- (tail(boot.R2Indices,1)+1):(tail(boot.R2Indices,1)+sum(lower.tri(inner)))
	boot.effectIndices <- boot.pathIndices + sum(lower.tri(inner))
	
	effects<-matrix(0, nrow(inner),nrow(inner))
	effects[lower.tri(effects)] <- boot$t0[boot.effectIndices]
	path.coef<-matrix(0, nrow(inner),nrow(inner))
	path.coef[lower.tri(path.coef)] <- boot$t0[boot.pathIndices]
	
  	pls<-list(path.coef = path.coef,
  		 loadings = boot$t0[boot.loadIndices],
  		 out.weights =  boot$t0[boot.weightIndices],
  		 r.sqr = boot$t0[boot.R2Indices],
  		 effects = effects)
  }
  
  #
  # If bootstrapping is not done, just call matrixpls.fit directly
  #
  
  else{
  	pls<-matrixpls.fit(indicatorCovariances,inner,outer,modes,scheme,scaled,tol,iter)
  }

  if (is.null(pls) || max(is.na(pls$loadings))) {
    print(paste("Iterative process is non-convergent with 'iter'=", 
                iter, " and 'tol'=", tol, sep=""))
    message("Algorithm stops") 
    stop("")
  }


  # We will be storing the results in the pls object
  class(pls)<-"matrixpls"
  

  # =======================================================
  # Define variables that are needed for reporting
  # =======================================================  
  IDM = inner
  lvs.names = rownames(IDM)
  dimnames(IDM) = list(lvs.names, lvs.names)
  lvs = nrow(IDM)
  blocks = unlist(lapply(outer, length))
  mvs = sum(blocks)
  names(blocks) = lvs.names
  blocklist = outer
  for (k in 1:length(blocks))
    blocklist[[k]] = rep(k,blocks[k])
  blocklist = unlist(blocklist)
  Mode = ifelse(modes=="A","Mode A","Mode B")


  # =======================================================
  # Stage 2: Path coefficients and total effects
  # =======================================================  

  endo = rowSums(IDM)
  endo[endo!=0] <- 1  # vector indicating endogenous LVs
  innmod <- as.list(1:sum(endo))

  for (aux in 1:sum(endo)) 
  {
    k1 <- which(endo==1)[aux]    # index for endo LV
    k2 <- which(IDM[k1,]==1)     # index for indep LVs
    inn.val <- c(pls$r.sqr[k1], 0, pls$Path(k1,k2))
    inn.lab <- c("R2", "Intercept", paste("path_",names(k2),sep=""))
    names(inn.val) <- NULL
    innmod[[aux]] <- data.frame(concept=inn.lab, value=round(inn.val,4))
  }
  
  names(innmod) <- lvs.names[endo!=0]  
  
  pls["inner.mod"] <- innmod
  
  # initialize
  efs.labs <- dir.efs <- ind.efs <- tot.efs <- NULL
  for (j in 1:lvs) 
  {
    for (i in j:lvs)
    {
      if (i != j) 
      {
        # labels
        efs.labs = c(efs.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
        # direct effects
        dir.efs = c(dir.efs, pls$Path[i,j])
        # indirect effects
        ind.efs = c(ind.efs, pls$effects[i,j]-pls$Path[i,j])
        # total effects
        tot.efs = c(tot.efs, pls$effects[i,j])
      }
    }
  }
  # results
  Effects = data.frame(relationships = efs.labs, 
                       dir.effects = dir.efs, 
                       ind.effects = ind.efs, 
                       tot.effects = tot.efs)

  pls["effects"] <- Effects
    
  # Composites scores
  
  if(matrixData){
    pls["latents"] <- NA
    pls["scores"] <- NA
  }
  else{
    outer_matrix <- sapply(outer_list,function(x) match(1:n,x,nomatch=0)>0)
    composites <- x %*% outer_matrix
    dimnames(composites) = list(rownames(x), lvs.names)
    pls["latents"] <- composites
    pls["scores"] <- composites

  }
  
  # =======================================================
  # Stage 3: Measurement loadings and communalities
  # =======================================================  
  
  # Communalities are needed later
  comu<-pls$loadings^2
  
  # =======================================================
  # Measurement model
  # =======================================================  
  outcor <- outmod <- as.list(1:lvs)
  for (j in 1:lvs)
  {
    aux <- which(blocklist==j)
    outmod[[j]] <- round(cbind(weights=out.weights[aux], std.loads=loads[aux], 
                               communal=comu[aux], redundan=redun[aux]), 4)
  }
  names(outmod) = lvs.names
  
  pls["outer.mod"] <- outmod

  # =======================================================
  # Unidimensionality
  # ======================================================= 

  Alpha = rep(1, lvs)   # Cronbach's Alpha for each block
  Rho = rep(1, lvs)     # D.G. Rho for each block
  eig.1st = rep(1, lvs) # first eigenvalue
  eig.2nd = rep(0, lvs) # second eigenvalue

  for (aux in 1:lvs) 
  {      
    if (blocks[aux] != 1) 
    { 
      if (modes[aux]=="A") 
      {
      	# Cronbach´s alpha
      	# PLSMP does not use information about weights  to determine if indicators are
      	# reverse codes, so we will not do it either
      	
        Alpha[aux] <- alpha(indicatorCovariances[outer_list[[aux]],outer_list[[aux]]])

        # dillon-goldstein rho
        
        # Equation 2.5 in Vinzi, V. E., Amato, S., & Trinchera, L. (2010). PLS path modeling: Recent developments and open issues for model assessment and improvement. V. Esposito Vinzi, W. Chin, J. Henseler & H. Wang, eds,“Handbook of Partial Least Squares-Concepts, Methods and Applications”, Springer, Berlin, Heidelberg, New York.

		loadings <- loadings[outer_list[[aux]]]
		
		# Standardize the loadings if not already standardized
		
		if(! scaled){
			loadings <- loadings * diag(indicatorCovariances)[outer_list[[aux]]]^2
		}
		squaredLoadings <- loadings^2
        Rho[aux] <- sum(squaredLoadings) /
        	 (sum(squaredLoadings) + sum(1-squaredLoadings))
      } else {  # modes[aux]=="B"
        Alpha[aux] = 0
        Rho[aux] = 0
      }
      eigenValues<-eigen(indicatorCorrelations[outer_list[[aux]],outer_list[[aux]]])
      eig.1st[aux] = eigenValues$values[1]
      eig.2nd[aux] = eigenValues$values[2]
    }
  }
  unidim = data.frame(Type.measure = Mode, 
                      MVs = blocks,
                      C.alpha = Alpha, 
                      DG.rho = Rho,
                      eig.1st, 
                      eig.2nd)
  rownames(unidim) = lvs.names

  pls["unidim"] <- unidim

  # =======================================================
  # Summary Inner model
  # =======================================================  
  exo.endo = ifelse(rowSums(IDM)==0
  exo.endo[rowSums(IDM)==0] = "Exogen"
  exo.endo[rowSums(IDM)!=0] = "Endogen"
  av.comu = rep(0, lvs)   # average communality
  av.redu = rep(0, lvs)   # average redundancy
  ave = rep(0, lvs)      # average variance extracted
  for (k in 1:lvs)
  {
    av.comu[k] = mean(comu[which(blocklist==k)])
    av.redu[k] = mean(redun[which(blocklist==k)])
    if (modes[k]=="A")
    {
      ave.num = sum(comu[which(blocklist==k)])
      ave.denom = sum(comu[which(blocklist==k)]) + sum(1-(comu[which(blocklist==k)]))
      ave[k] = ave.num / ave.denom
    }
  }
  names(ave) = lvs.names
  innsum = data.frame(LV.Type = exo.endo, 
                      Measure = Mode, 
                      MVs = blocks, 
                      R.square = pls["r.sqr"], 
                      Av.Commu = av.comu, 
                      Av.Redun = av.redu, 
                      AVE = ave)
  rownames(innsum) = lvs.names
  pls["inner.sum"] <- innsum

  
  # =======================================================
  # GoF Index
  # =======================================================  
  
  pls["gof"] <- get_gof(comu, pls["r.sqr"], blocks, IDM)

  # =======================================================
  # Results
  # =======================================================  
  
  if(matrixData) obs<-NA
  else obs <- nrow(matrixData)
  
  pls["model"] <- list(IDM=IDM, blocks=blocks, scheme=scheme, modes=modes, scaled=scaled, 
                boot.val=boot.val, plsr=plsr, obs=obs, br=br, 
                tol=tol, iter=iter, n.iter=NA, outer=outer)
  # deliver dataset?
  if (dataset && ! matrixData){
	pls["data"] <- as.matrix(Data)

  }
  else{
    pls["data"] <- NULL
  }
  # deliver bootstrap validation results? 
  if (doBootstrap) 
  {
	# TODO: Refactor this to use the boot package to calculate the descriptives
	
    WEIGS <- boots$t[,boot.weightIndices]
    LOADS <- boots$t[,boot.loadIndices]
    PATHS <- boots$t[,boot.pathIndices]
    TOEFS <- boots$t[,boot.effectIndices]
    RSQRS <- boots$t[,boot.R2Indices]

    colnames(WEIGS) <- mvs.names
    WB <- data.frame(Original = wgs.orig, 
                     Mean.Boot = apply(WEIGS, 2, mean), 
                     Std.Error = apply(WEIGS, 2, sd), 
                     perc.025 = apply(WEIGS, 2, function(x) quantile(x, 0.025)),
                     perc.975 = apply(WEIGS, 2, function(x) quantile(x, 0.975)))
    # Loadings
    colnames(LOADS) <- mvs.names
    LB <- data.frame(Original = loads.orig, 
                     Mean.Boot = apply(LOADS, 2, mean),
                     Std.Error = apply(LOADS, 2, sd), 
                     perc.025 = apply(LOADS, 2, function(x) quantile(x, 0.025)),
                     perc.975 = apply(LOADS, 2, function(x) quantile(x, 0.975)))
    # Path coefficients
    colnames(PATHS) <- path.labs
    PB <- data.frame(Original = path.orig, 
                     Mean.Boot = apply(PATHS, 2, mean),
                     Std.Error = apply(PATHS, 2, sd), 
                     perc.025 = apply(PATHS, 2, function(x) quantile(x, 0.025)),
                     perc.975 = apply(PATHS, 2, function(x) quantile(x, 0.975)))

    # Total effects
    colnames(TOEFS) <- Path.efs[, 1]
    TE <- data.frame(Original = Path.efs[, 4], 
                     Mean.Boot = apply(TOEFS, 2, mean), 
                     Std.Error = apply(TOEFS, 2, sd),
                     perc.025 = apply(TOEFS, 2, function(x) quantile(x, 0.025)), 
                     perc.975 = apply(TOEFS, 2, function(x) quantile(x, 0.975)))
    # R-squared
    colnames(RSQRS) <- lvs.names[endo == 1]
    RB <- data.frame(Original = r2.orig, 
                     Mean.Boot = apply(RSQRS, 2, mean),
                     Std.Error = apply(RSQRS, 2, sd), 
                     perc.025 = apply(RSQRS, 2, function(x) quantile(x, 0.025)),
                     perc.975 = apply(RSQRS, 2, function(x) quantile(x, 0.975)))
    # Bootstrap Results
    pls["boot"] <- list(weights = WB, 
                     loadings = LB, 
                     paths = PB, 
                     rsq = RB, 
                     total.efs = TE)
    }
    pls["boot"] <- NA

  } else {
      pls["boot"] <- FALSE
  }

  return(pls)
}





