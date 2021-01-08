# Load the Tenenhaus et al 2005 model and data from semPLS
library(semPLS)
data(ECSImobi)
data(mobi)

# Reflective and empty formative model
reflective<- ECSImobi$M
formative <- t(reflective)
formative[] <- 0

# Estimation using covariance matrix
model <- list(inner =  t(ECSImobi$D),
              reflective = reflective,
              formative = formative)


S <- cor(mobi)

matrixpls.ModeA <- matrixpls(S, model, innerEstim = innerEstim.centroid)
matrixpls.Fixed <- matrixpls(S, model, weightFun = weightFun.fixed)

cei(matrixpls.ModeA)
cei(matrixpls.ModeA, matrixpls.Fixed)