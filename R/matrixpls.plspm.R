library(plspm)
library(psych)

#'@title A PLS-PM compatibility wrapper for MatrixPLS
#'
#'@description
#'Estimate path models with latent variables by partial least squares approach
#'
#'@details
#'The function \code{matrixpls.plspm} estimates a path model by partial least squares
#'approach providing the full set of results. \cr
#'
#'The argument \code{inner_matrix} is a matrix of zeros and ones that indicates
#'the structural relationships between latent variables. This must be a lower
#'triangular matrix. \code{inner_matrix} will contain a 1 when column \code{j}
#'affects row \code{i}, 0 otherwise. \cr
#'
#'@param Data A numeric matrix or data frame containing the manifest variables.
#'@param inner_matrix A square (lower triangular) boolean matrix representing the
#'inner_matrix model (i.e. the path relationships betwenn latent variables).
#'@param outer_list List of vectors with column indices from \code{x} indicating 
#'the sets of manifest variables asociated to the latent variables
#'(i.e. which manifest variables correspond to the latent variables).
#'Length of \code{outer_list} must be equal to the number of rows of \code{inner_matrix}.
#'@param modes A character vector indicating the type of measurement for each
#'latent variable. \code{"A"} for reflective measurement or \code{"B"} for
#'formative measurement (\code{NULL} by default). The length of \code{modes}
#'must be equal to the length of \code{outer_list}).
#'@param scheme A string of characters indicating the type of inner_matrix weighting
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
#'@return An object of class \code{"plspm"}. 
#'@return \item{outer.mod}{Results of the outer (measurement) model. Includes:
#'outer weights, standardized loadings, communalities, and redundancies}
#'@return \item{inner_matrix.mod}{Results of the inner_matrix (structural) model. Includes: path
#'coefficients and R-squared for each endogenous latent variable}
#'@return \item{latents}{Matrix of standardized latent variables (variance=1
#'calculated divided by \code{N}) obtained from centered data (mean=0)}
#'@return \item{scores}{Matrix of latent variables used to estimate the inner_matrix
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
#'@return \item{inner_matrix.sum}{Summarized results by latent variable of the inner_matrix
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
#'@author Mikko Rönkkö
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
#'@seealso \code{\link[plspm]{plspm}}, \code{\link[plspm]{plspm.fit}}, \code{\link[plspm]{plot.plspm}}
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
#'  # inner_matrix model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0) 
#'  LOY = c(1,0,0,0,1,0)
#'  sat_inner_matrix = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'
#'  # outer model list
#'  sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'
#'  # vector of modes (reflective indicators)
#'  sat_mod = rep("A", 6)
#'
#'  # apply plspm
#'  satpls = matrixpls.plspm(satisfaction, sat_inner_matrix, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)
#'  
#'  # summary of results
#'  summary(satpls)
#'
#'  # default plot (inner_matrix model)
#'  plot(satpls)
#'  }
#'

matrixpls.plspm <-
function(Data, inner_matrix, outer_list, modes = NULL, scheme = "centroid", 
         scaled = TRUE, tol = 0.00001, iter = 100, boot.val = FALSE, 
         br = NULL, plsr = FALSE, dataset = TRUE){
         
    # ======================================================= checking arguments =======================================================
    
    if (!is.matrix(x) && !is.data.frame(x)) 
        stop("Invalid object 'x'. Must be a numeric matrix or data frame.")
    if (is.null(rownames(x))) 
        rownames(x) <- 1:nrow(x)
    if (is.null(colnames(x))) 
        colnames(x) <- paste("MV", 1:ncol(x), sep = "")
    if (!is.matrix(inner_matrix)) 
        stop("Invalid argument 'inner_matrix'. Must be a matrix.")
    if (nrow(inner_matrix) != ncol(inner_matrix)) 
        stop("Invalid argument 'inner_matrix'. Must be a square matrix.")
    for (j in 1:ncol(inner_matrix)) for (i in 1:nrow(inner_matrix)) {
        if (i <= j) 
            if (inner_matrix[i, j] != 0) 
                stop("argument 'inner_matrix' must be a lower triangular matrix")
        if (length(intersect(inner_matrix[i, j], c(1, 0))) == 0) 
            stop("elements in 'inner_matrix' must be '1' or '0'")
    }
    if (is.null(dimnames(inner_matrix))) 
        lvs.names <- paste("LV", 1:ncol(inner_matrix), sep = "")
    if (!is.null(rownames(inner_matrix))) 
        lvs.names <- rownames(inner_matrix)
    if (!is.null(colnames(inner_matrix))) 
        lvs.names <- colnames(inner_matrix)
    if (!is.list(outer)) 
        stop("Invalid argument 'outer'. Must be a list.")
    if (length(outer) != nrow(inner_matrix)) 
        stop("Number of rows of 'inner_matrix' does not coincide with length of 'outer'.")
    if (is.null(modes)) {
        modes <- rep("A", length(outer))
        warning("Argument 'modes' missing. Default reflective 'modes' is used.")
    }
    if (length(outer) != length(modes)) {
        warning("Warning: Invalid length of 'modes'. Default reflective 'modes' is used.")
        modes <- rep("A", length(outer))
    }
    for (i in 1:length(modes)) if (modes[i] != "A" && modes[i] != "B") 
        modes[i] <- "A"
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
            if (mode(br) != "numeric" || length(br) != 1 || (br%%1) != 0 || br < 100 || br > 1000) {
                warning("Warning: Invalid argument 'br'. Default 'br=100' is used.")
                br <- 100
            }
        } else br <- 100
    }
    if (!is.logical(plsr)) 
        plsr <- FALSE
    if (mode(tol) != "numeric" || length(tol) != 1 || tol <= 0 || tol > 0.001) {
        warning("Warning: Invalid argument 'tol'. Default 'tol=0.00001' is used.")
        tol <- 1e-05
    }
    if (mode(iter) != "numeric" || length(iter) != 1 || iter < 100) {
        warning("Warning: Invalid argument 'iter'. Default 'iter=100' is used.")
        iter <- 100
    }
    if (!is.logical(dataset)) 
        dataset <- TRUE
    
    # Stop if there are unimplemented options
    if (plsr) 
        stop("PLS regression is currently not implemented in matrixpls")
    
    
	indicatorCovariances <- cov(x)
    
    indicatorCorrelations <- cov2cor(indicatorCovariances)
    
    doBootstrap = boot.val
    if (doBootstrap && nrow(Data) <= 10) {
        warning("Bootstrapping stopped: very few cases.")
        doBootstrap <- FALSE
    }
    
    # Boot will provide also estimates for the full data
    
    if (doBootstrap) {
        
        # Selection matrix
        selectWeightsFromMatrix <- matrix(0, nrow(indicatorCovariances), ncol(inner_matrix))
        for (i in 1:ncol(inner_matrix)) {
            selectWeightsFromMatrix[outer[[i]], i] <- 1
        }
        selectWeightsFromMatrix <- which(selectWeightsFromMatrix == 1)
                
        boot.results <- boot(x, function(d, p) {
            
            cov.matrix <- cov(d[p, ], d[p, ])
            
            # TODO: Fix this call
            pls <- matrixpls(cov.matrix, inner_matrix, outer, modes, scheme, scaled, tol, iter)
            
            if (is.null(pls)) {
                return(rep(NA, length(selectWeightsFromMatrix) * 2 + nrow(inner_matrix) + sum(inner_matrix) + sum(lower.tri(inner_matrix))))
            } else {
                return(c(pls$loadings[selectWeightsFromMatrix], pls$out.weights[selectWeightsFromMatrix], pls$r.sqr, pls$path.coef[which(inner_matrix == 1)], pls$effects[lower.tri(inner_matrix)]))
            }
        }, br)
        
        # Unpack the results from t0
        boot.loadIndices <- 1:ncol(Data)
        boot.weightIndices <- ncol(Data) + boot.loadIndices
        boot.R2Indices <- tail(boot.weightIndices, 1) + 1:ncol(inner_matrix)
        boot.pathIndices <- tail(boot.R2Indices, 1) + 1:sum(inner_matrix)
        boot.effectIndices <- tail(boot.pathIndices, 1) + 1:sum(lower.tri(inner_matrix))
        
        path.coefs <- matrix(0, ncol(inner_matrix), ncol(inner_matrix))
        path.coefs[which(inner_matrix == 1)] <- boot.results$t0[boot.pathIndices]
        
        loads <- matrix(0, nrow(indicatorCovariances), ncol(inner_matrix))
        loads[selectWeightsFromMatrix] <- boot.results$t0[boot.loadIndices]
        
        out.weights <- matrix(0, nrow(indicatorCovariances), ncol(inner_matrix))
        out.weights[selectWeightsFromMatrix] <- boot.results$t0[boot.weightIndices]
        
        effs <- matrix(0, ncol(inner_matrix), ncol(inner_matrix))
        effs[lower.tri(inner_matrix)] <- boot.results$t0[boot.effectIndices]
        
        pls <- list(path.coefs = path.coefs, loadings = loads, out.weights = out.weights, r.sqr = boot.results$t0[boot.R2Indices], effects = effs)
    } else {
        # If bootstrapping is not done, just call matrixpls directly
        # TODO: Fix this call
        pls <- matrixpls(indicatorCovariances, inner_matrix, outer, modes, scheme, scaled, tol, iter)
    }
    
    if (is.null(pls) || max(is.na(pls$loadings))) {
        print(paste("Iterative process is non-convergent with 'iter'=", iter, " and 'tol'=", tol, sep = ""))
        message("Algorithm stops")
        stop("")
    }
    
    # ======================================================= Define variables that are needed for reporting =======================================================
    IDM = inner_matrix
    if (is.null(dimnames(inner_matrix))) 
        lvs.names <- paste("LV", 1:ncol(inner_matrix), sep = "")
    if (!is.null(rownames(inner_matrix))) 
        lvs.names <- rownames(inner_matrix)
    if (!is.null(colnames(inner_matrix))) 
        lvs.names <- colnames(inner_matrix)
    
    dimnames(IDM) = list(lvs.names, lvs.names)
    lvs = nrow(IDM)
    blocks = unlist(lapply(outer, length))
    mvs = sum(blocks)
    mvs.names = 1:mvs
    names(blocks) = lvs.names
    blocklist = outer
    for (k in 1:length(blocks)) blocklist[[k]] = rep(k, blocks[k])
    blocklist = unlist(blocklist)
    Mode = ifelse(modes == "A", "Mode A", "Mode B")
    
    # ======================================================= Stage 2: Path coefficients and total effects =======================================================
    
    colnames(pls$path.coefs) <- lvs.names
    rownames(pls$path.coefs) <- lvs.names
    
    endo = rowSums(IDM)
    endo[endo != 0] <- 1  # vector indicating endogenous LVs
    innmod <- as.list(1:sum(endo))
    
    for (aux in 1:sum(endo)) {
        k1 <- which(endo == 1)[aux]  # index for endo LV
        k2 <- which(IDM[k1, ] == 1)  # index for indep LVs
        inn.val <- c(pls$r.sqr[k1], 0, pls$path.coefs[k1, k2])
        inn.lab <- c("R2", "Intercept", paste("path_", names(k2), sep = ""))
        names(inn.val) <- NULL
        innmod[[aux]] <- data.frame(concept = inn.lab, value = round(inn.val, 4))
    }
    
    
    names(innmod) <- lvs.names[endo != 0]
    
    pls[["inner_matrix.mod"]] <- innmod
    
    # initialize
    efs.labs <- dir.efs <- ind.efs <- tot.efs <- NULL
    for (j in 1:lvs) {
        for (i in j:lvs) {
            if (i != j) {
                # labels
                efs.labs = c(efs.labs, paste(lvs.names[j], "->", lvs.names[i], sep = ""))
                # direct effects
                dir.efs = c(dir.efs, pls$path.coefs[i, j])
                # indirect effects
                ind.efs = c(ind.efs, pls$effects[i, j] - pls$path.coefs[i, j])
                # total effects
                tot.efs = c(tot.efs, pls$effects[i, j])
            }
        }
    }
    
    Effects = data.frame(relationships = efs.labs, dir.effects = dir.efs, ind.effects = ind.efs, tot.effects = tot.efs)
    
    pls[["effects"]] <- Effects
    
    # Composites scores
    
        outer_matrix <- sapply(outer_list, function(x) match(1:mvs, x, nomatch = 0) > 0)
        composites <- x %*% outer_matrix
        dimnames(composites) = list(rownames(x), lvs.names)
        pls[["latents"]] <- composites
        pls[["scores"]] <- composites
        
    
    # ======================================================= Stage 3: Measurement loadings and communalities =======================================================
    
    comu <- rowSums(pls$loadings)^2
    redun = rep(0, mvs)
    for (j in 1:lvs) {
        if (endo[j] == 1) {
            redun[blocklist == j] <- comu[blocklist == j] * pls$r.sqr[j]
        }
    }
    # ======================================================= Measurement model =======================================================
    outcor <- outmod <- as.list(1:lvs)
    
    for (j in 1:lvs) {
        aux <- which(blocklist == j)
        outmod[[j]] <- round(cbind(weights = pls$out.weights[aux, j], std.loads = pls$loadings[aux, j], communal = comu[aux], redundan = redun[aux]), 4)
    }
    names(outmod) = lvs.names
    pls[["outer.mod"]] <- outmod
    
    # ======================================================= Unidimensionality =======================================================
    
    Alpha = rep(1, lvs)  # Cronbach's Alpha for each block
    Rho = rep(1, lvs)  # D.G. Rho for each block
    eig.1st = rep(1, lvs)  # first eigenvalue
    eig.2nd = rep(0, lvs)  # second eigenvalue
    
    for (aux in 1:lvs) {
        if (blocks[aux] != 1) {
            if (modes[aux] == "A") {
                # Cronbach\302\264s alpha PLSMP does not use information about weights to determine if indicators are reverse codes, so we will not do it either
                
                Alpha[aux] <- alpha(indicatorCovariances[outer_list[[aux]], outer_list[[aux]]])$total$std.alpha
                
                # dillon-goldstein rho
                
                # Equation 2.5 in Vinzi, V. E., Amato, S., & Trinchera, L. (2010). PLS path modeling: Recent developments and open issues for model assessment and improvement. V. Esposito Vinzi, W. Chin, J. Henseler
                # & H. Wang, eds,\342\200\234Handbook of Partial Least Squares-Concepts, Methods and Applications\342\200\235, Springer, Berlin, Heidelberg, New York.
                
                loads <- pls$loadings[outer_list[[aux]], aux]
                
                # Standardize the loads if not already standardized
                
                if (!scaled) {
                  loads <- loads * diag(indicatorCovariances)[outer_list[[aux]]]^2
                }
                squaredLoads <- loads^2
                Rho[aux] <- sum(squaredLoads)/(sum(squaredLoads) + sum(1 - squaredLoads))
            } else {
                Alpha[aux] = 0
                Rho[aux] = 0
            }
            eigenValues <- eigen(indicatorCorrelations[outer_list[[aux]], outer_list[[aux]]])
            eig.1st[aux] = eigenValues$values[1]
            eig.2nd[aux] = eigenValues$values[2]
        }
    }
    unidim = data.frame(Type.measure = Mode, MVs = blocks, C.alpha = Alpha, DG.rho = Rho, eig.1st, eig.2nd)
    rownames(unidim) = lvs.names
    
    pls[["unidim"]] <- unidim
    
    # ======================================================= Summary inner_matrix model =======================================================
    
    exo.endo = ifelse(rowSums(IDM) == 0, "Exogen", "Endogen")
    
    av.comu = rep(0, lvs)  # average communality
    av.redu = rep(0, lvs)  # average redundancy
    ave = rep(0, lvs)  # average variance extracted
    for (k in 1:lvs) {
        av.comu[k] = mean(comu[which(blocklist == k)])
        av.redu[k] = mean(redun[which(blocklist == k)])
        if (modes[k] == "A") {
            ave.num = sum(comu[which(blocklist == k)])
            ave.denom = sum(comu[which(blocklist == k)]) + sum(1 - (comu[which(blocklist == k)]))
            ave[k] = ave.num/ave.denom
        }
    }
    names(ave) = lvs.names
    innsum = data.frame(LV.Type = exo.endo, Measure = Mode, MVs = blocks, R.square = pls[["r.sqr"]], Av.Commu = av.comu, Av.Redun = av.redu, AVE = ave)
    rownames(innsum) = lvs.names
    pls[["inner_matrix.sum"]] <- innsum
    
    
    # ======================================================= GoF Index =======================================================

	# this is just to avoid compilation warnings. Needs to be fixed. 
	R2 <- NULL
    # average of communalities

    R2.aux <- R2[endo == 1]
    comu.aux <- n.comu <- 0
    for (j in 1:lvs) {
        if (length(which(blocklist == j)) > 1) {
            comu.aux = comu.aux + sum(comu[which(blocklist == j)])
            n.comu = n.comu + length(which(blocklist == j))
        }
    }
    pls[["gof"]] <- sqrt((comu.aux/n.comu) * mean(R2.aux))
    
    
    # ======================================================= Results =======================================================
    
    
    pls[["model"]] <- list(IDM = IDM, blocks = blocks, scheme = scheme, modes = modes, scaled = scaled, boot.val = boot.val, plsr = plsr, obs = nrow(Data), br = br, tol = tol, iter = iter, n.iter = NA, outer = outer)
    # deliver dataset?
    pls[["data"]] <- as.matrix(Data)
        
    # deliver bootstrap validation results?
    if (doBootstrap) {
        
        # TODO: Refactor this to use the boot package to calculate the descriptives
        
        WEIGS <- boot.results$t[, boot.weightIndices]
        LOADS <- boot.results$t[, boot.loadIndices]
        PATHS <- boot.results$t[, boot.pathIndices]
        TOEFS <- boot.results$t[, boot.effectIndices]
        RSQRS <- boot.results$t[, boot.R2Indices]
        
        colnames(WEIGS) <- mvs.names
        WB <- data.frame(Original = boot.results$t0[boot.weightIndices], Mean.Boot = apply(WEIGS, 2, mean, na.rm = TRUE), Std.Error = apply(WEIGS, 2, sd, na.rm = TRUE), perc.025 = apply(WEIGS, 2, function(x) quantile(x, 
            0.025, na.rm = TRUE)), perc.975 = apply(WEIGS, 2, function(x) quantile(x, 0.975, na.rm = TRUE)))
        # Loadings
        colnames(LOADS) <- mvs.names
        LB <- data.frame(Original = boot.results$t0[boot.loadIndices], Mean.Boot = apply(LOADS, 2, mean, na.rm = TRUE), Std.Error = apply(LOADS, 2, sd, na.rm = TRUE), perc.025 = apply(LOADS, 2, function(x) quantile(x, 
            0.025, na.rm = TRUE)), perc.975 = apply(LOADS, 2, function(x) quantile(x, 0.975, na.rm = TRUE)))
        # Path coefficients colnames(PATHS) <- path.labs
        PB <- data.frame(Original = boot.results$t0[boot.pathIndices], Mean.Boot = apply(PATHS, 2, mean, na.rm = TRUE), Std.Error = apply(PATHS, 2, sd, na.rm = TRUE), perc.025 = apply(PATHS, 2, function(x) quantile(x, 
            0.025, na.rm = TRUE)), perc.975 = apply(PATHS, 2, function(x) quantile(x, 0.975, na.rm = TRUE)))
        
        # Total effects colnames(TOEFS) <- Path.efs[, 1]
        TE <- data.frame(Original = boot.results$t0[boot.effectIndices], Mean.Boot = apply(TOEFS, 2, mean, na.rm = TRUE), Std.Error = apply(TOEFS, 2, sd, na.rm = TRUE), perc.025 = apply(TOEFS, 2, function(x) quantile(x, 
            0.025, na.rm = TRUE)), perc.975 = apply(TOEFS, 2, function(x) quantile(x, 0.975, na.rm = TRUE)))
        # R-squared colnames(RSQRS) <- lvs.names[endo == 1]
        RB <- data.frame(Original = boot.results$t0[boot.R2Indices], Mean.Boot = apply(RSQRS, 2, mean, na.rm = TRUE), Std.Error = apply(RSQRS, 2, sd, na.rm = TRUE), perc.025 = apply(RSQRS, 2, function(x) quantile(x, 
            0.025, na.rm = TRUE)), perc.975 = apply(RSQRS, 2, function(x) quantile(x, 0.975, na.rm = TRUE)))
        # Bootstrap Results
        pls[["boot"]] <- list(weights = WB, loadings = LB, paths = PB, rsq = RB, total.efs = TE)
        
        pls[["boot.raw"]] <- boot.results
        
    } else {
        pls[["boot"]] <- FALSE
    }
    
    # Set the class to be plspm so that the results can be printed with plspm functions
    class(pls) <- "plspm"
    return(pls)
}




 
