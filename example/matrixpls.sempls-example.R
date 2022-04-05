data(ECSImobi)

ecsi <- matrixpls.sempls(model=ECSImobi, data=mobi, wscheme="pathWeighting")

## Values of 'sempls' objects
names(ecsi)
ecsi$outer_weights
ecsi$outer_loadings
ecsi$path_coefficients
ecsi$total_effects
