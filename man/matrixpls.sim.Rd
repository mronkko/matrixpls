\name{matrixpls.sim}
\alias{matrixpls.sim}
\title{Monte Carlo simulations with matrixpls}
\usage{
  matrixpls.sim(nRep = NULL, model = NULL, n = NULL, ...,
    cilevel = 0.95,
    citype = c("norm", "basic", "stud", "perc", "bca"),
    boot.R = 500, fitIndices = fitSummary)
}
\arguments{
  \item{nRep}{Number of replications. If any of the
  \code{n}, \code{pmMCAR}, or \code{pmMAR} arguments are
  specified as lists, the number of replications will
  default to the length of the list(s), and \code{nRep}
  need not be specified.}

  \item{...}{All other arguments are forwared through to
  \code{\link[simsem]{sim}} or
  \code{\link{matrixpls.boot}}.}

  \item{cilevel}{Confidence level. This argument will be
  forwarded to the \code{\link[boot]{boot.ci}} when
  calculating the confidence intervals.}

  \item{citype}{Type of confidence interval. This argument
  will be forwarded to the \code{\link[boot]{boot.ci}} when
  calculating the confidence intervals.}

  \item{boot.R}{Number of bootstrap replications to use to
  estimate standard errors or \code{FALSE} to disable
  bootstrapping.}

  \item{fitIndices}{A function that returns a list of fit
  indices for the model. Setting this argument to
  \code{NULL} disables fit indices.}

  \item{model}{There are two options for this argument: 1.
  lavaan script or lavaan parameter table, or 2. a list
  containing three matrices \code{inner},
  \code{reflective}, and \code{formative} defining the free
  regression paths in the model.}

  \item{n}{ Sample size. Either a single value, or a list
  of values to vary sample size across replications. The
  \code{n} argument can also be specified as a random
  distribution object; if any resulting values are
  non-integers, the decimal will be rounded.  }
}
\value{
  An object of class \code{\link[simsem]{SimResult-class}}.
}
\description{
  Performs Monte Carlo simulations of
  \code{\link{matrixpls}} with the
  \code{\link[simsem]{sim}} function of the \code{simsem}
  package. The standard errors and confidence intervals are
  estimated with the \code{\link[boot]{boot}} and
  \code{\link[boot]{boot.ci}} functions of the \code{boot}
  package.
}
\details{
  This funtion calls the \code{\link[simsem]{sim}} function
  from the \code{simsem} package to perform Monte Carlo
  simulations with matrixpls. The function parses the model
  parameters and replaces it with a function call that
  estimates the model and bootstrapped standard errors and
  confidence intervals with \link{matrixpls.boot}.

  If the \code{generate} or \code{rawdata} arguments are
  not specified in the \code{\link[simsem]{sim}} arguments
  then the \code{model} argument will be used for data
  generation and must be specified in lavaan format.
}
\examples{
#
#  Runs an example from the SimSem website
# 
#  https://github.com/simsem/simsem/wiki/Example-2:-Getting-Started-\%28lavaan-input\%29
#

if(require(simsem)){
  
  popModel <- "
f1 =~ 1*y1 + 0.6*y2 + 0.7*y3
f2 =~ 1*y4 + 1.1*y5 + 0.9*y6
f3 =~ 1*y7 + 1.2*y8 + 1.1*y9
f1 ~~ 0.8*f1
f2 ~~ 0.9*f2
f3 ~~ 0.4*f3
f1 ~~ 0.4*f2
f1 ~~ 0.2*f3
f2 ~~ 0.3*f3
y1 ~~ 0.5*y1
y2 ~~ 1.1*y2
y3 ~~ 0.8*y3
y4 ~~ 0.4*y4
y5 ~~ 0.4*y5
y6 ~~ 0.8*y6
y7 ~~ 0.8*y7
y8 ~~ 0.5*y8
y9 ~~ 0.6*y9
"
  
  analyzeModel <- "
f1 =~ y1 + y2 + y3
f2 =~ y4 + y5 + y6
f3 =~ y7 + y8 + y9
"
  
  # Only run 10 replications with 50 bootstrap replications each so that the 
  #example runs faster
  
  Output <- matrixpls.sim(10, analyzeModel, n=200, generate=popModel, boot.R=50,
                          multicore = FALSE, stopOnError = TRUE)
  getCutoff(Output, 0.05)
  plotCutoff(Output, 0.05)
  summary(Output)
  
  
} else{
  print("This example requires the simsem package")
}

#
# Runs the second model from
#
# Aguirre-Urreta, M., & Marakas, G. (2013). Partial Least Squares and Models with Formatively 
# Specified Endogenous Constructs: A Cautionary Note. Information Systems Research. 
# doi:10.1287/isre.2013.0493

if(require(simsem)){
  
  library(MASS)
  
  X <- diag(15) 
  X[upper.tri(X, diag=TRUE)] <- c(
    1.000,
    0.640, 1.000,
    0.640, 0.640, 1.000,
    0.640, 0.640, 0.640, 1.000,
    0.096, 0.096, 0.096, 0.096, 1.000,
    0.096, 0.096, 0.096, 0.096, 0.640, 1.000,
    0.096, 0.096, 0.096, 0.096, 0.640, 0.640, 1.000,
    0.096, 0.096, 0.096, 0.096, 0.640, 0.640, 0.640, 1.000,
    0.115, 0.115, 0.115, 0.115, 0.192, 0.192, 0.192, 0.192, 1.000,
    0.115, 0.115, 0.115, 0.115, 0.192, 0.192, 0.192, 0.192, 0.640, 1.000,
    0.115, 0.115, 0.115, 0.115, 0.192, 0.192, 0.192, 0.192, 0.640, 0.640,
    1.000,
    0.115, 0.115, 0.115, 0.115, 0.192, 0.192, 0.192, 0.192, 0.640, 0.640,
    0.640, 1.000,
    0.000, 0.000, 0.000, 0.000, 0.271, 0.271, 0.271, 0.271, 0.325, 0.325,
    0.325, 0.325, 1.000,
    0.000, 0.000, 0.000, 0.000, 0.271, 0.271, 0.271, 0.271, 0.325, 0.325,
    0.325, 0.325, 0.300, 1.000,
    0.000, 0.000, 0.000, 0.000, 0.271, 0.271, 0.271, 0.271, 0.325, 0.325,
    0.325, 0.325, 0.300, 0.300, 1.000
  )
  X <- X + t(X) - diag(diag(X)) 
  
  colnames(X) <- rownames(X) <- c(paste("Y",1:12,sep=""),paste("X",1:3,sep=""))
  
  # Print the population covariance matrix X to see that it is correct
  
  X
  
  # The estimated model in Lavaan syntax
  
  analyzeModel <- "
ksi =~ Y1 + Y2 + Y3 + Y4
eta1 <~ X1 + X2 + X3
eta2 =~ Y5 + Y6 + Y7 + Y8
eta3 =~ Y9 + Y10 + Y11 + Y12

eta1 ~ ksi
eta2 ~ eta1
eta3 ~ eta1
"
  
  # Only run 100 replications with 2 bootstrap replications each so that the 
  # example runs faster. Generate the data outside simsem because simsem
  # does not support drawing samples from population matrix
  
  dataSets <- lapply(1:100, function(x){
    mvrnorm(300,                 # Sample size
            rep(0,15),           # Means
            X)                   # Population covarancematrix
  })
  
  Output <- matrixpls.sim(model = analyzeModel, rawData = dataSets, boot.R=2,
                          multicore = FALSE, stopOnError = TRUE)
  
  summary(Output)
  
  
} else{
  print("This example requires the simsem package")
}

# Simulates the direct effects model from 
#
# Dijkstra, T. K., & Schermelleh-Engel, K. (2013). Consistent Partial Least Squares for Nonlinear
# Structural Equation Models. Psychometrika, 1–20. doi:10.1007/s11336-013-9370-0
#
# using PLSc, disattenuated summed scales, and disattenuated factor scores
#

####################################################################################################
#
# Simulation parameters
#
####################################################################################################

SAMPLE <- 100
REPLICATIONS <- 50
MULTICORE <- FALSE

MODEL <- "
Eta1 =~ .8*y1 + .8*y2 + .8*y3
Eta2 =~ .8*y4 + .8*y5 + .8*y6 
Eta3 =~ .8*y7 + .8*y8 + .8*y9 

y1~~.36*y1
y2~~.51*y2
y3~~.64*y3
y4~~.36*y4
y5~~.51*y5
y6~~.64*y6
y7~~.36*y7
y8~~.51*y8
y9~~.64*y9

Eta1~~1*Eta1
Eta2~~1*Eta2
Eta2~~-.3*Eta1

Eta3~.5*Eta1 + -.3*Eta2
Eta3~~.57*Eta3  # error tem is calculated as 1-(.5*(.5 + -.3*-.3) + -.3*(-.3 + -.3*.5))
"


####################################################################################################
#
# Run simulations
#
####################################################################################################

library(matrixpls)

sim.PLSf1 <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, #General setup
                           parameterEstimator = params.plsc,
                           fm = "minres", # Setup disattenuation with minres factor analysis
                           outerEstimator = outer.factor, # Use factor scores as proxies.
                           innerEstimator = NULL, # Factor score proxies do not use inner estimation
                           boot.R = FALSE, # We are not interested in bootstrap
                           multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE)


sim.PLSf2 <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, #General setup
                           parameterEstimator = params.plsc, 
                           fm = "minres", # Setup disattenuation with minres factor analysis
                           outerEstimator = outer.fixedWeights, # Use equal weights
                           innerEstimator = NULL, # Equal weightproxies do not use inner estimation
                           boot.R = FALSE, # We are not interested in bootstrap
                           multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE)

sim.PLSc <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, #General setup
                          parameterEstimator = params.plsc, # Setup PLSc disattenuation 
                          boot.R = FALSE, # We are not interested in bootstrap
                          multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE)

####################################################################################################
#
# Report results
#
####################################################################################################


# SimSem summary methods have a bug that causes a crash with inadmissible results. (therefore 'try')

try(summary(sim.PLSf1))
try(summary(sim.PLSf2))
try(summary(sim.PLSc))

# Plot the sampling distribution of the parameters

par(mfrow=c(2,3))
for(i in 1:6){
  estimatorIndex <- (i-1) \%\% 3 +1
  obj <- switch(estimatorIndex, sim.PLSf1,sim.PLSf2,sim.PLSc)
  var <- colnames(obj@coef)[ceiling(i/3)]
  title <- paste(switch(estimatorIndex, "Factor scores", "Summed scales", "PLSc"),var,sep=": ")
  plot(density(obj@coef[,var], na.rm = TRUE), main = title)	
}

}
\seealso{
  \code{\link[simsem]{sim}},
  \code{\link[simsem]{SimResult-class}},
  \code{\link{matrixpls.boot}}, \code{\link[boot]{boot}},
  \code{\link[boot]{boot.ci}}
}

