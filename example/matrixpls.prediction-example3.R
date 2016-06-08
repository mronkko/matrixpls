library(semPLS)
data(ECSImobi)

# Construct the model based on the ECSImobi example

model <- list(inner = t(ECSImobi$D),
     reflective = ECSImobi$M,
     formative = t(ECSImobi$M))
model$formative[] <- 0


# Estimation using covariance matrix

matrixpls.out <- matrixpls(cov(mobi),
                           model = model,
                           standardize = FALSE)

# Calculate within-sample predictions

predictions <- predict(matrixpls.out, mobi)

# Calculate root mean squared prediction errors

sqrt(apply((predictions-mobi[61:75,])**2,2,mean))

# Q2 predictive relevance

q2(mobi, predictions, model)

# Crossvalidation predictions using holdout samples and blindfolding

predictions.holdout <- matrixpls.crossvalidate(mobi,
                                               model = model,
                                               standardize = FALSE)

q2(mobi, predictions.holdout, model)

# Mimic the blindfolding procedure used in semPLS

predictions.blindfold <- matrixpls.crossvalidate(mobi,
                                                 model = model,
                                                 blindfold = TRUE,
                                                 predictionType = "redundancy",
                                                 groups = 4)


q2(mobi, predictions.blindfold, model)

# The results are similar to semPLS q2 values, but not exactly identical due to differences 
# in how the two packages apply standardization when calculating predictions.

ecsi <- sempls(model=ECSImobi, data=mobi, E="C")
qSquared(ecsi, d=4, dlines = FALSE)




