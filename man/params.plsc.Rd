\name{params.plsc}
\alias{params.plsc}
\title{Parameter estimation with an adaptation of PLSc algorithm}
\usage{
  params.plsc(S, W, model, fm = "dijkstra", ...)
}
\arguments{
  \item{S}{Covariance matrix of the data.}

  \item{W}{Weight matrix, where the indicators are on
  colums and composites are on the rows.}

  \item{fm}{factoring method for estimating the corrected
  factor loadings. \code{"dijkstra"} will use the
  correction presented by Dijkstra (2011), where the PLS
  estimates of the factor loadings for a latent variable
  are multiplied with a scalar \code{a}, which is
  calculated by a simple formula that approximately
  minimizes squared residual covariances of the factor
  model. \code{minres}, \code{wls}, \code{gls}, \code{pa},
  and \code{ml} will use a factor analysis and passing this
  parameter through to \code{\link[psych]{fa}}.}

  \item{model}{There are two options for this argument: 1.
  lavaan script or lavaan parameter table, or 2. a list
  containing three matrices \code{inner},
  \code{reflective}, and \code{formative} defining the free
  regression paths in the model.}

  \item{...}{All other arguments are passed through to
  \code{matrixpls.weights} and \code{parameterEstimator}.}
}
\value{
  A named vector of parameter estimates.
}
\description{
  Estimates the model parameters with an adapted version of
  Dijkstra's (2011) PLSc. The parameters between latent
  variables are estimated using two-stages least squares
  estimation using a composite covariance matrix that has
  been disattenuated with estimated composite
  reliabilities.
}
\details{
  \code{params.plsc} estimates the statistical model
  described by \code{model} with the following steps. If
  \code{model} is not in the native format, it is converted
  to the native format containing matrices \code{inner},
  \code{reflective}, and \code{formative}. The weights
  \code{W} and the data covariance matrix \code{S} are used
  to calculate the composite covariance matrix \code{C} and
  the indicator-composite covariance matrix \code{IC}.

  \code{C} is then disattenuated. The reliability estimates
  used in the dissattenuation process are calculated by
  using Dijkstra's (2011) correction or with a separate
  factor analysis for each indicator block. Indicators that
  are used in a composite but are not specified as
  reflective indicators for the latent variable that the
  composite approximates are assumed to be perfectly
  reliable.

  The disattenuated \code{C} is used to estimate the the
  inner part of the model with separate OLS regressions in
  the same way as in \code{\link{params.regression}}. This
  differs from Dijkstra's (2011) PLSc estimator that uses
  2SLS. The reason for not using 2SLS is that PLS is
  commonly used with models that do not contain variables
  that could be used as instruments in the 2SLS estimation.

  The results from the disattenuation process are used as
  estimates of the \code{reflective} part of the model and
  the \code{formative} part of the model estimated with
  separate OLS regressions using the \code{S} and \code{IC}
  matrices.

  Those elements of teh indicator-composite covariance
  matrix \code{IC} that correspond to factor loadings are
  replaced with the factor loading estimates.

  The dissattenuation code for Dijkstra's method is adapted
  from Huang (2013), which is based on Dijkstra (2011).
}
\examples{
library(plspm)

# Run the customer satisfaction examle form plspm

# load dataset satisfaction
data(satisfaction)
# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
colnames(inner) <- rownames(inner)

# Reflective model
list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

reflective<- matrix(
  c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
  27,6, dimnames = list(colnames(satisfaction)[1:27],colnames(inner)))

# empty formative model

formative <- matrix(0, 6, 27, dimnames = list(colnames(inner), colnames(satisfaction)[1:27]))

# Estimation using covariance matrix

matrixpls.out <- matrixpls(cov(satisfaction[,1:27]), model = list(inner = inner,
                                                                  reflective = reflective,
                                                                  formative = formative),
                           parameterEstimator = params.plsc)

summary(matrixpls.out)


# Simulates the direct effects model from 
#
# Dijkstra, T. K., & Schermelleh-Engel, K. (2013). Consistent Partial Least Squares for Nonlinear Structural Equation Models. Psychometrika, 1–20. doi:10.1007/s11336-013-9370-0
#
# using PLSc, disattenuated summed scales, and disattenuated factor scores
#
# This simulation requires features not included in matrixpls 0.2.0, which means that you need to 
# use the development version from github. See this page for installing instructions
# https://github.com/mronkko/matrixpls

####################################################################################################
#
# Simulation parameters
#
####################################################################################################

SAMPLE <- 100
REPLICATIONS <- 500

# Set this to FALSE allows using breakpoints and traceback. This is useful if you want to experiment
# with the analysis, and encounter problems

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
													 parameterEstimator = params.plsc, fm = "minres", # Setup disattenuation with minres factor analysis
													 outerEstimator = outer.factor, innerEstimator = NULL, # Use factor scores as proxies. There is no need for inner estimation
													 boot.R = FALSE, # We are not interested in bootstrap
													 multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE) # Control parameters

sim.PLSf2 <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, #General setup
													 parameterEstimator = params.plsc, fm = "minres", # Setup disattenuation with minres factor analysis
													 outerEstimator = outer.fixedWeights, innerEstimator = NULL, # Use equal weights.  There is no need for inner estimation
													 boot.R = FALSE, # We are not interested in bootstrap
													 multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE) # Control parameters
	
sim.PLSc <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, #General setup
													parameterEstimator = params.plsc, # Setup PLSc disattenuation 
													boot.R = FALSE, # We are not interested in bootstrap
													multicore = MULTICORE, fitIndices = NULL, stopOnError = TRUE) # Control parameters


####################################################################################################
#
# Report results
#
####################################################################################################


# SimSem summary methods have a bug that causes an error with inadmissible results. (therefore 'try')

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
\author{
  Mikko Rönkkö, Wenjing Huang, Theo Dijkstra
}
\references{
  Huang, W. (2013). PLSe: Efficient Estimators and Tests
  for Partial Least Squares (Doctoral dissertation).
  University of California, Los Angeles. Dijkstra, T. K.
  (2011). Consistent Partial Least Squares estimators for
  linear and polynomial factor models. A report of a
  belated, serious and not even unsuccessful attempt.
  Comments are invited. Retrieved from
  http://www.rug.nl/staff/t.k.dijkstra/consistent-pls-estimators.pdf
}
\seealso{
  Other parameter estimators:
  \code{\link{params.plsregression}},
  \code{\link{params.regression}}
}

