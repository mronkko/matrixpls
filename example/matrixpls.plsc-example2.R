# Simulates the direct effects model from 
#
# Dijkstra, T. K., & Schermelleh-Engel, K. (2013). Consistent Partial Least Squares for Nonlinear Structural Equation Models. Psychometrika, 1–20. doi:10.1007/s11336-013-9370-0
#
# using PLSc, disattenuated summed scales, and disattenuated factor scores
#
# This simulation requires features not included in matrixpls 0.2.0, which means that you need to 
# use the development version from github. See this page for installing instructions
# https://github.com/mronkko/matrixpls

####################################################################################################
#
# Simulation parameters
#
####################################################################################################

SAMPLE <- 100
REPLICATIONS <- 500

# Set this to FALSE allows using breakpoints and traceback. This is useful if you want to experiment
# with the analysis, and encounter problems

MULTICORE <- TRUE

MODEL <- "
Eta1 =~ .8*y1 + .8*y2 + .8*y3
Eta2 =~ .8*y4 + .8*y5 + .8*y6 
Eta3 =~ .8*y7 + .8*y8 + .8*y9 

y1~~.36*y1
y2~~.51*y2
y3~~.64*y3
y4~~.36*y4
y5~~.51*y5
y6~~.64*y6
y7~~.36*y7
y8~~.51*y8
y9~~.64*y9

Eta1~~1*Eta1
Eta2~~1*Eta2
Eta2~~-.3*Eta1

Eta3~.5*Eta1 + -.3*Eta2
Eta3~~.57*Eta3  # error tem is calculated as 1-(.5*(.5 + -.3*-.3) + -.3*(-.3 + -.3*.5))
"


####################################################################################################
#
# Run simulations
#
####################################################################################################

library(matrixpls)

sim.PLSf1 <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, #General setup
													 parameterEstimator = params.plsc, fm = "minres", # Setup disattenuation with minres factor analysis
													 outerEstimator = outer.factor, innerEstimator = NULL, # Use factor scores as proxies. There is no need for inner estimation
													 boot.R = FALSE, # We are not interested in bootstrap
													 multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE) # Control parameters

sim.PLSf2 <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, #General setup
													 parameterEstimator = params.plsc, fm = "minres", # Setup disattenuation with minres factor analysis
													 outerEstimator = outer.fixedWeights, innerEstimator = NULL, # Use equal weights.  There is no need for inner estimation
													 boot.R = FALSE, # We are not interested in bootstrap
													 multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE) # Control parameters
	
sim.PLSc <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, #General setup
													parameterEstimator = params.plsc, # Setup PLSc disattenuation 
													boot.R = FALSE, # We are not interested in bootstrap
													multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE) # Control parameters


####################################################################################################
#
# Report results
#
####################################################################################################


# SimSem summary methods have a bug that causes an error with inadmissible results. (therefore 'try')

try(summary(sim.PLSf1))
try(summary(sim.PLSf2))
try(summary(sim.PLSc))

# Plot the sampling distribution of the parameters

par(mfrow=c(2,3))
for(i in 1:6){
	estimatorIndex <- (i-1) %% 3 +1
	obj <- switch(estimatorIndex, sim.PLSf1,sim.PLSf2,sim.PLSc)
	var <- colnames(obj@coef)[ceiling(i/3)]
	title <- paste(switch(estimatorIndex, "Factor scores", "Summed scales", "PLSc"),var,sep=": ")
	plot(density(obj@coef[,var], na.rm = TRUE), main = title)	
}

