# Run the example from ASGSCA package using GSCA estimation

library(ASGSCA) 
library(gdata)
library(matrixpls)

# Run the GSCA example from the ASGSCA package

data(GenPhen)
W0 <- matrix(c(rep(1,2),rep(0,8),rep(1,2),rep(0,8),rep(1,3),rep(0,7),rep(1,2)),nrow=8,ncol=4)
B0 <- matrix(c(rep(0,8),rep(1,2),rep(0,3),1,rep(0,2)),nrow=4,ncol=4)

# ASGSCA uses random numbers as starting values 
# set.seed(1)

GSCA.res <- GSCA(GenPhen,W0, B0,estim=TRUE,path.test=FALSE, 
                             latent.names=c("Gene1","Gene2","Clinical pathway 1","Clinical pathway 2"))

# Setup matrixpls to estimate the same model. Note GSCA places dependent variables
# on columns but matrixpls uses rows for dependent variables

inner <- t(B0)
# formative <- matrix(0,4,8)
# reflective <- W0
formative <- t(W0)
reflective <- matrix(0,8,4)

colnames(formative) <- rownames(reflective) <- names(GenPhen)
colnames(inner) <- rownames(inner) <- rownames(formative) <- colnames(reflective) <-
  c("Gene1","Gene2","Clinical pathway 1","Clinical pathway 2")

model <- list(inner = inner, 
              reflective = reflective,
              formative = formative)

# Estimate using alternating least squares

matrixpls.res1 <- matrixpls(cov(GenPhen),  model,
                            outerEstimators = outer.GSCA,
                            innerEstimator = inner.GSCA)

# Estimate using direct minimization of the estimation criterion

matrixpls.res2 <- matrixpls(cov(GenPhen),  model,
                            weightFunction = weight.optim,
                            optimCriterion = optim.GSCA,
                            # Set the convergence criterion to be slightly stricter than normally
                            control = list(reltol = 1e-12))

# Compare the weights
do.call(cbind,lapply(list(GSCA.res$Weight,t(attr(matrixpls.res1,"W")), t(attr(matrixpls.res2,"W"))),
       function(W){
         W <- unmatrix(W)
         W[W!=0]
       }))

# Compare the paths


# Check the criterion function values
print(optim.GSCA(matrixpls.res1))
print(optim.GSCA(matrixpls.res2))

