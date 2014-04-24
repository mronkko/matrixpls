library(semPLS)

data(ECSImobi)
ecsi.sempls <- sempls(model=ECSImobi, data=mobi, wscheme="pathWeighting")

ecsi <- matrixpls.sempls(model=ECSImobi, data=mobi, wscheme="pathWeighting")

checkEquals(ecsi.sempls,ecsi)

ecsi

## create plots
densityplot(ecsi)
densityplot(ecsi, use="prediction")
densityplot(ecsi, use="residuals")

## Values of 'sempls' objects
names(ecsi)
ecsi$outer_weights
ecsi$outer_loadings
ecsi$path_coefficients
ecsi$total_effects


### using convenience methods to sempls results
## path coefficients
pathCoeff(ecsi)

## total effects
totalEffects(ecsi)

## get loadings and check for discriminant validity
(l <- plsLoadings(ecsi))
# outer loadings
print(l, type="outer", digits=2)
# outer loadings greater than 0.5
print(l,type="outer", cutoff=0.5, digits=2)
# cross loadings greater than 0.5
print(l, type="cross", cutoff=0.5, digits=2)


### R-squared
rSquared(ecsi)