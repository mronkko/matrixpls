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


# Mimic the blindfolding procedure used in semPLS

predictions.blindfold <- matrixpls.crossvalidate(mobi,
                                                 model = model,
                                                 blindfold = TRUE,
                                                 predictionType = "redundancy",
                                                 groups = 4)

# Q2 predictive relevance

q2(mobi, predictions.blindfold, model)



