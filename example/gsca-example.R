# Run the example from ASGSCA package using GSCA estimation

# library(Matrix)
# library(MASS)
# 
# source("~/Downloads/ASGSCA/R/GSCAestim.R")
# source("~/Downloads/ASGSCA/R/GSCAtest.R")

 if(require(ASGSCA)) {
  
  # Run the GSCA example from the ASGSCA package
  
  #Scenario (g) in Romdhani et al. (2014): 4 SNPs mapped to 2 genes and 4 
  #traits involved in 2 clinical pathways 
  #In total: 8 observed variables and 4 latent variables.
  #One of the traits is involved in both clinical pathways.
  #One gene is connected to one of the clinical pathways and
  #the other to both of them.
  
  data(GenPhen)
  W0 <- matrix(c(rep(1,2),rep(0,8),rep(1,2),rep(0,8),rep(1,3),rep(0,7),rep(1,2)),nrow=8,ncol=4)
  B0 <- matrix(c(rep(0,8),rep(1,2),rep(0,3),1,rep(0,2)),nrow=4,ncol=4)
  
  #Estimation only
  GSCA.res <- GSCA(GenPhen,W0, B0,estim=TRUE,path.test=FALSE, latent.names=c("Gene1","Gene2","Clinical pathway 1","Clinical pathway 2"))
  
  # Setup matrixpls to estimate the same model. Note GSCA places dependent variables
  # on columns but matrixpls uses rows for dependent variables
  
  inner <- t(B0)
  reflective <- matrix(0,8,4)
  formative <- t(W0)
  
  colnames(formative) <- rownames(reflective) <- names(GenPhen)
  colnames(inner) <- rownames(inner) <- rownames(formative) <- colnames(reflective) <- c("Gene1","Gene2","Clinical pathway 1","Clinical pathway 2")
  
  model <- list(inner = inner, 
                reflective = reflective,
                formative = formative)
    
  # Estimate using alternating least squares
  
  matrixpls.res1 <- matrixpls(cov(GenPhen),  model,
                              outerEstimators = outer.GSCA,
                              innerEstimator = inner.GSCA,
                              tol = 0.000000000001)

  print("Comparing alternating least squares estimates to original estimats")
  # Compare thw weights
  print(GSCA.res$Weight)
  print(t(attr(matrixpls.res1,"W")))

  print(GSCA.res$Weight - t(attr(matrixpls.res1,"W")))
  
  
  # Compare the paths
  print(GSCA.res$Path)
  print(t(attr(matrixpls.res1,"beta")))
  
  print(GSCA.res$Path-t(attr(matrixpls.res1,"beta")))
  
  stop("done")
  
  # Estimate using direct minimization of the estimation criterion
  
  matrixpls.res2 <- matrixpls(cov(GenPhen),  model,
                              weightFunction = weight.optim,
                              optimCriterion = optim.GSCA)

  print("Comparing numeric minimization estimates to original estimats")
  
  # Compare thw weights
  print(GSCA.res$Weight)
  print(t(attr(matrixpls.res2,"W")))
  
  print(GSCA.res$Weight - t(attr(matrixpls.res2,"W")))
  
  
  # Compare the paths
  print(GSCA.res$Path)
  print(t(attr(matrixpls.res2,"beta")))
  
  print(GSCA.res$Path-t(attr(matrixpls.res2,"beta")))
  
  
  # Check the criterion function values
  print(optim.GSCA(matrixpls.res1))
  print(optim.GSCA(matrixpls.res2))
  
} else{
  print("This example requires the ASGSCA package")
}

