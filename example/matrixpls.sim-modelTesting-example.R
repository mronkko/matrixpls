# This example contains code cotributed by Florian Schuberth (florian.schuberth@uni-wuerzburg.de)

library(matrixpls)
library(boot) # packages for bootstrap
library(simsem)  # package for simulation of data sets
library(lavaan)  # package for latent variable analysis
library(expm) # matrix square roots
 
set.seed(2015)

### Data generation  #####
# The model of Henseler & Sarstedt (2013) is used page 573
MODEL <- "
# measurement model
XI1 =~ .6*x1 + .7*x2 + .8*x3
XI2 =~ .6*y1 + .7*y2 + .8*y3
XI3 =~ .6*y4 + .7*y5 + .8*y6

# Variances of the indicators
x1~~.64*x1 ; x2~~.51*x2; x3~~.36*x3

y1~~.64*y1; y2~~.51*y2; y3~~.36*y3

y4~~.64*y4; y5~~.51*y5; y6~~.36*y6

#structural model
XI2~.6*XI1
XI3~.6*XI2 +.0*XI1

# (Co-)Variance of the common factors
XI1~~1*XI1
XI2~~0.64*XI2
XI3~~0.64*XI3
"

ndraws <- 1000        # number of draws
nobservations <- 1200  # number of observations
multi <- TRUE       # parallel computing TRUE or FALSE

# Create data sets
datmat.factor <-  sim(nRep = ndraws,
                      model = MODEL,
                      n = nobservations,
                      dataOnly = TRUE,
                      multicore = multi)

# Estimate the model with PLSc (Dijkstra & Henseler (2015)) using the data sets created earlier

sim.res <- matrixpls.sim(model = MODEL,
                         disattenuate = TRUE,
                         reliabilities = reliabilityEstim.weightLoadingProduct,
                         parametersReflective = estimator.plscLoadings,
                         rawData = datmat.factor,
                         multicore = multi,
                         boot.R = FALSE)

summary(sim.res)

# Extract estimates
plsc.res <- sim.res@extraOut

# Function that calculates the squared Euclidean distance between the empirical and
# the model implied correlation matrix for every estimation
calculateTestStat <- function(matrixpls.res) {
  S <- attr(matrixpls.res, "S")
  Sigma <- fitted(matrixpls.res)
  sum((S - Sigma)[lower.tri(S, diag = FALSE)] ^ 2)
}

# calculate the distances for the initial data sets
teststat <- lapply(plsc.res, calculateTestStat)


# Function that tranforms the data sets in the way proposed by Yuan & Hayashi (2003), 
# see also Dijkstra & Henseler (2015)

# The bootstrap approach only works correctly if  the empirical variance-covariance matrix of the
# transformed data set (Z) equals the model-implied covariance matrix, cov(Z)=Sigma_hat 
# (Beran & Srivastava, 1985). Therefore, we need to check that that the data set X corresponds to
# S and Sigma_hat. If X is not scaled to means zero and unit variances, S and Sigma_hat need to be
# covariance matrices, otherwise, if X is scaled, S and Sigma_hat need to be correlation matrices.
#
# Beran, R. & Srivastava, M. S.
# Bootstrap tests and confidence regions for functions of a covariance matrix
# The Annals of Statistics, JSTOR, 1985, 13, 95-115  
#
# Because PLS weights are not scale free, that is, they depend on the scale of the raw data,
# the function will check that S matches the covariance matrix of the data and automatically
# rescale the dataset

scaleDataSet <- function(data, matrixpls.res) {
  
  S   <- attr(matrixpls.res, "S")
  Sig <- fitted(matrixpls.res)
  
  ## S^{-1/2}
  S_half   <- solve(expm::sqrtm(S))
  Sig_half <- expm::sqrtm(Sig)
  
  ## Scale factor Sig_hat^(1/2)%*%S^(-1/2)
  sf <- S_half %*% Sig_half  
  
  sdata <- as.matrix(data)
  
  if(any(cov(data) != S)) {sdata <- scale(sdata)}

  sdata <- sdata %*% sf
  colnames(sdata) <- colnames(data)
  sdata
}

# Create transformed data sets
scaleddata <- mapply(scaleDataSet, 
                     data=datmat.factor,
                     matrixpls.res=plsc.res,
                     SIMPLIFY = FALSE)


# Calculate the empirical reference distribution
sim.res2 <- matrixpls.sim(model = MODEL, disattenuate = TRUE,
                          reliabilities = reliabilityEstim.weightLoadingProduct,
                          parametersReflective = estimator.plscLoadings,
                          rawData = scaleddata,
                          multicore = F,
                          boot.R = 500,
                          # Calculate the test statistic for each bootstrap replication
                          extraFun = calculateTestStat)


rejection=mapply(function(x,y){
  # The last column contains the test statistic
  crit=quantile(x$t[,ncol(x$t)],0.9) # critical value for comparison
  return(y>crit)      
},x=sim.res2@extraOut,y=teststat)

# Present results
# FALSE: non-rejection
# TRUE: rejection

table(rejection)

