#'@details
#'
#'The two step GSCA algorithm is designed to minimize.
#'
#'\code{SS(ZV-ZWA)}
#'
#'The parameter matrix \code{A} contains all model parameters including
#'\code{inner}, reflective \code{inner}, and \code{formative}. The weight
#'matrices \code{V} and \code{W}, which can contain duplicate elements,
#'contain the indicator weights.
#'
#'The first step of GSCA estimation method is calculation of regressions
#'coefficients \code{A} given weights \code{W} and \code{V}. The function
#'\code{\link{inner.GSCA}} update the part of \code{A} corresponding to 
#'regressions between the composites, corresponding to \code{E} matrix in 
#'matrixpls. The regressions between composites and indicators are estimated
#'in \code{\link{outer.GSCA}}.
#'
#'This algorithm for estimating the relationships between the composites
#'is therefore identical to the PLS path weighting scheme with
#'the exception that correlations are not used for inverse relationships and
#'there is no falling back to identity scheme for composites that are not
#'connected to other composites

#'The second step of GSCA is calculating a new set of weights conditional on
#'the regression coeffcients \code{A} to minimize the sum of error terms in
#'the regressions. This step is implemented in \code{\link{outer.GSCA}} after
#'updating the regresions between indicators and composites.

#'The implementation of GSCA in matrixpls differs from the Hwang & Takane (2004)
#'version in that during the first step, only regressions between composites are
#'estimated. The regressions between composites and indicators are estimated by
#'the second stage 
#'the indicators and compositess. Since these covariances need to be calculated in the second step, it is more
#'efficient to not calculate them during the first step.
#'
