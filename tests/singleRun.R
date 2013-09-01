
library(plspm)

# Do a single run the performance of plspm and MatrixPLS with 1000 bootstrap samples using the customer satisfaction example

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

convergenceCheckFunction <- function(W,W_new){sum((abs(W_new) - abs(W))^2)} 

						S <- cov(satisfaction[,1:27])
						m <- matrixpls(S, model = nativeModel, weightRelations = weightRelations,
											innerEstimator = matrixpls.innerEstimator.centroid,
											tol = 0.00001, iter = 100, convergenceCheckFunction = convergenceCheckFunction,
											validateInput = FALSE)

