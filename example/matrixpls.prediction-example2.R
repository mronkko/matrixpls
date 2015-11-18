library(lavaan)

## The industrialization and Political Democracy example 
## Bollen (1989), page 332. (Adopted from the lavaan package.)

model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

matrixpls.out <- matrixpls(cov(PoliticalDemocracy[1:60,]),model)

# Calculate predictions using the hold-out sample

predictions <- predict(matrixpls.out, PoliticalDemocracy[61:75,])

summary(predictions)

# Calculate root mean squared prediction errors

sqrt(apply((predictions - PoliticalDemocracy[61:75,])^2,2,mean))

