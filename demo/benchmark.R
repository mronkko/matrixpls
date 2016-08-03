library(matrixpls)
library(parallel)

if(require(plspm) & require(semPLS)){
	
	# Benchmark the performance of plspm, semPLS, and matrixpls with R bootstrap samples
  # using the customer satisfaction example from plspm
	
	# load dataset satisfaction
	data(satisfaction)
	
	# inner model matrix
	IMAG = c(0,0,0,0,0,0)
	EXPE = c(1,0,0,0,0,0)
	QUAL = c(0,1,0,0,0,0)
	VAL = c(0,1,1,0,0,0)
	SAT = c(1,1,1,1,0,0)
	LOY = c(1,0,0,0,1,0)
	sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
	colnames(sat_inner) <- rownames(sat_inner)
	
	# outer model list
	sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
	# vector of modes (reflective indicators)
	sat_mod = rep(c("A"), 6)
	
	# Set up the native matrixpls model
  
	reflective <- matrix(0, max(unlist(sat_outer)),nrow(sat_inner))
	rownames(reflective) <- colnames(satisfaction)[1:nrow(reflective)]
	colnames(reflective) <- colnames(sat_inner)
	
	col <- 1		
	
	for(outerSpec in sat_outer){
		reflective[outerSpec,col] <- 1
		col <- col + 1 
	}
	
	formative <- t(reflective)
	formative[] <- 0
	
	nativeModel <- list(inner=sat_inner, 
											reflective=reflective, 
											formative=formative) 
	
	W.model <- t(reflective)
	
  # Setup the semPLS model
  
  i <- which(sat_inner == 1)
  j <- which(reflective == 1)
  
  semPLSmodel <- plsm(data = satisfaction[,1:27], 
                      strucmod = cbind(source = colnames(sat_inner)[col(sat_inner)[i]],
                                            target = rownames(sat_inner)[row(sat_inner)[i]]),
                      measuremod = cbind(source = colnames(reflective)[col(reflective)[j]],
                                              target = rownames(reflective)[row(reflective)[j]]))
  
  # Number of replications and iterations
	maxiter <- 100
  R <- 100
	Rboot <- 1000
  
	#
	# Single thread comparisons
	#
	
	# Store the old option values
	
	opt_boot.ncpus <- getOption("boot.ncpus")
	opt_boot.parallel <- getOption("boot.parallel")
	opt_cores <- getOption("mc.cores")
	
	options(boot.ncpus = 1)
	options(boot.parallel = "no")
	options(mc.cores = 1)
	
	plspm_boot <- system.time(plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=Rboot))
	plspm_no_boot <- system.time(replicate(R,plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE)))
	
	matrixpls.plspm_boot <- system.time(matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=Rboot))
	matrixpls.plspm_no_boot <- system.time(replicate(R,matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE)))

	sempls_boot <- system.time(bootsempls(sempls(semPLSmodel,satisfaction[,1:27], maxit = maxiter, verbose = FALSE), nboot = Rboot, verbose = FALSE))
	sempls_no_boot <- system.time(replicate(R,system.time(sempls(semPLSmodel,satisfaction[,1:27], maxit = maxiter, verbose = FALSE))))
	
	matrixpls.sempls_no_boot <- system.time(replicate(R,matrixpls.sempls(semPLSmodel,satisfaction[,1:27], maxit = maxiter)))
  
	matrixpls_boot <- system.time(matrixpls.boot(satisfaction[,1:27], R = Rboot, model = nativeModel, W.model = W.model,
																							 innerEstim = innerEstim.centroid,
																							 tol = 0.00001, iter = maxiter, validateInput = FALSE))
	S <- cov(satisfaction[,1:27])
	
	matrixpls_no_boot <- system.time(replicate(R,matrixpls(S, model = nativeModel, W.model = W.model,
																													 innerEstim = innerEstim.centroid,
																													 tol = 0.00001, iter = maxiter, validateInput = FALSE)))
	#
	# Multithreaded bootstrapping. Only supported by matrixpls.
	#
	
	cores = detectCores()
	
	options(boot.ncpus = cores)
	options(mc.cores = cores)
  
	if(.Platform$OS.type == "windows") {
		options(boot.parallel = "snow")
	} else {
		options(boot.parallel = "multicore")
	}
	
	
	matrixpls.plspm_boot_multicore <- system.time(matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=Rboot))
	
	matrixpls_boot_multicore <- system.time(matrixpls.boot(satisfaction[,1:27], R = Rboot, model = nativeModel, W.model = W.model,
																												 innerEstim = innerEstim.centroid,
																												 tol = 0.00001, iter = maxiter, validateInput = FALSE))
	
	#
	# Set the options back to default values
	#
	
	
	options(boot.ncpus = opt_boot.ncpus)
	options(boot.parallel = opt_boot.parallel)
	options(mc.cores = opt_cores)
	
	# Print the times
	
	print(rbind(plspm_no_boot,matrixpls.plspm_no_boot, sempls_no_boot, matrixpls.sempls_no_boot, matrixpls_no_boot))
	
	print(rbind(plspm_boot,matrixpls.plspm_boot, matrixpls_boot, sempls_boot,
        matrixpls.plspm_boot_multicore, matrixpls_boot_multicore))
	
} else{
	print("This example requires the plspm and semPLS packages")
}
