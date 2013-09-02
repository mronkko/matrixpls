library(plspm)

# Benchmark the performance of plspm and MatrixPLS with 1000 bootstrap samples using the customer satisfaction example

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

#
# Test 1: Apply plspm
#

plspm <- system.time(a <- plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=1000))


# Store the old option values

opt_boot.ncpus <- getOption("boot.ncpus")
opt_boot.parallel <- getOption("boot.parallel")
opt_cores <- getOption("cores")

options(boot.ncpus = 1)
options(boot.parallel = "no")
options(cores = 1)


#
# Test 2: Apply MatrixPLS, single thread
#

matrixpls.plspm <- system.time(b <- matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=1000))

#
# Test 3: Apply MatrixPLS, multiple threads
#

cores = detectCores()

options(boot.ncpus = cores)

if(.Platform$OS.type == "windows") {
	options(boot.parallel = "snow")
} else {
	options(boot.parallel = "multicore")
}
options(cores = cores)

matrixpls.plspm_multithread <- system.time(c <- matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val = TRUE, br=1000))

#
# Test 4: Apply MatrixPLS, multiple threads, native call
#

# Set up the native model
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

weightRelations <- t(reflective)

# Do the estimation 

matrixpls <- system.time(d <- boot(as.matrix(satisfaction[,1:27]),
								 function(data, indices){
								 	S <- cov(data[indices,])
								 	matrixpls(S, model = nativeModel, weightRelations = weightRelations,
								 						innerEstimator = matrixpls.innerEstimator.centroid,
								 						tol = 0.00001, iter = 100, validateInput = FALSE)
								 },
								 1000))

#
# 
#

#
# Set the options back to default values
#


options(boot.ncpus = opt_boot.ncpus)
options(boot.parallel = opt_boot.parallel)
options(cores = opt_cores)

# Print the times

rbind(plspm,matrixpls.plspm,matrixpls.plspm_multithread,matrixpls)
