# Runs the power analysis example from 

# Aguirre-Urreta, M. I., & Rönkkö, M. (2015). Sample size determination and statistical power 
# analysis in PLS using R: An annotated tutorial. Communications of the Association for 
# Information Systems, 36(1). Retrieved from http://aisel.aisnet.org/cais/vol36/iss1/3

library(matrixpls)
library(simsem)
library(lavaan)

model<-"! regressions
A=~0.7*x1
A=~0.7*x2
A=~0.7*x3
B=~0.7*x4
B=~0.7*x5
B=~0.8*x6
B=~0.8*x7
C=~0.6*x8
C=~0.6*x9
C=~0.6*x10
C=~0.8*x11
C=~0.8*x12
D=~0.8*x13
D=~0.8*x14
D=~0.8*x15
D ~ 0.3*A
C ~ 0.1*B
C ~ 0.5*A
D ~ 0.3*C

! error variances, variances and covariances
A ~~ 1.0*A
B ~~ 1.0*B
C ~~ 0.71*C
D ~~ 0.725*D
B ~~ 0.3*A
x1 ~~ 0.51*x1
x2 ~~ 0.51*x2
x3 ~~ 0.51*x3
x4 ~~ 0.51*x4
x5 ~~ 0.51*x5
x6 ~~ 0.36*x6
x7 ~~ 0.36*x7
x8 ~~ 0.64*x8
x9 ~~ 0.64*x9
x10 ~~ 0.64*x10
x11 ~~ 0.36*x11
x12 ~~ 0.36*x12
x13 ~~ 0.36*x13
x14 ~~ 0.36*x14
x15 ~~ 0.36*x15"


# Normal data

output <- matrixpls.sim(1000, model, n=100, # Basic parameters
                        multicore = TRUE, # Additional parameters
                        completeRep = TRUE)

summary(output)

# Non-normal data

distributions <- bindDist(skewness = c(0.50, 0.50, 0.50, 0.50, 0.50,
                                       0.50, 0.50, 0.50, 0.75, 0.75,
                                       0.75, 0.75, 0.75, 0.75, 0.75),
                          kurtosis = c(1, 1, 1, 1, 1,
                                       1, 1, 1, 1.50, 1.50,
                                       1.50, 1.50, 1.50, 1.50, 1.50))

nonNorm <- matrixpls.sim(1000, model, n=100, # Basic parameters 
                         indDist = distributions, # Indicator distr
                         multicore = FALSE, # Additional params
                         completeRep = TRUE)

summary(nonNorm)

# CBSEM, normal data

cbsem <- sim(nRep = 1000,
             model = gsub("(0-9.)+\\*","",model), # Free all params
             generate = model, n=100, # Basic parameters
             lavaanfun = "sem", # Do a SEM analysis
             multicore = TRUE, # Additional parameters
             completeRep = TRUE)

summary(cbsem)