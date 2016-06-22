library(matrixpls) # package to estimate VBSEM
library(boot) # packages for bootstrap
library(simsem)  # package for simulation of data sets
library(lavaan)  # package for latent variable analysis

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

ndraws = 1000        # number of draws
nobservations = 1200  # number of observations
multi = FALSE       # parallel computing TRUE or FALSE

# Create data sets
datmat.factor <-  sim(nRep = ndraws,
                      model = MODEL,
                      n = nobservations,
                      dataOnly = TRUE,
                      multicore = F)

# Estimate the model with PLSc (Dijkstra & Henseler (2015)) using the data sets created earlier
sim.res <- matrixpls.sim(model = MODEL,
                         disattenuate = TRUE,
                         reliabilities = reliabilityEstim.weightLoadingProduct,
                         parametersReflective = estimator.plscLoadings,
                         rawData = datmat.factor,
                         multicore = multi,
                         boot.R = FALSE)

summary(sim.res)

# save estimates
plsc.res <- sim.res@extraOut

# Calculate the model implied correlation matrix for every estimation
Sighats <- lapply(plsc.res, function(x)
  fitted(x))

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
scaleDataSet <- function(data,Sig){    
  S <- cor(data)
  # singular value decomposition of S and Sigma
  Ssvd=svd(S); Sigsvd=svd(Sig)
  dS=Ssvd$d; dSig=Sigsvd$d
  uS=Ssvd$u; uSig=Sigsvd$u
  vS=Ssvd$v; vSig=Sigsvd$v
  #S^{-1/2}
  S_half=uS%*%(diag(dS^(-1/2)))%*%t(vS)
  # Sigma_hat^(1/2)
  Sig_half=uSig%*%(diag(dSig^(1/2)))%*%t(vSig)
  # scale factor Sig_hat^(1/2)%*%S^(-1/2)
  sf=S_half%*%Sig_half
  
  sdata <- as.matrix(data)%*%sf
  colnames(sdata) <- colnames(data)
  sdata
}

# Create transformed data sets
scaleddata <- mapply(scaleDataSet, 
                     data=datmat.factor,
                     Sig=Sighats,
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
#FALSE: non-rejection
# TRUE: rejection
table(rejection)

# Results for different sample sizes:
# n=300
# FALSE  TRUE 
#  914    86 

# n=600
# FALSE  TRUE 
# 921    79 

# n=1200
# FALSE  TRUE 
#  901    99

# The results support the findings by Dijkstra & Henseler (2015) (Table 7) 

# References:
#Yuan, Ke-Hai, and Kentaro Hayashi. "Bootstrap approach to inference and power analysis based on three test statistics for covariance structure models." British Journal of Mathematical and Statistical Psychology 56.1 (2003): 93-110.

#Dijkstra, Theo K., and Henseler, Joerg. "Consistent and asymptotically normal PLS estimators for linear structural equations." Computational statistics & data analysis 81 (2015): 10-23.

#Henseler, Joerg, and Sarstedt, Marko. "Goodness-of-fit indices for partial least squares path modeling." Computational Statistics 28.2 (2013): 565-580.

