#'@title Data generation with MatrixPLS
#'
#'@description
#'Generates a sample of data using a covariance matrix and distributions.
#
#'@details
#'
#'@param Sigma Population covariance matrix
#'@param nobs Number of observations to generate
#'@param dist A data distribution object created by the bindDist function of simsem package (optional, if omitted normal distribution is used)
#'@return A 
#'@author Mikko R\303\266nkk\303\266
#'
#'@references Sunthud Pornprasertmanit, Patrick Miller and Alexander Schoemann (2013). simsem: SIMulated Structural Equation Modeling.. R package version 0.5-0.
#'  http://CRAN.R-project.org/package=simsem
#'

library("simsem")

matrixpls.gendata <- function(Sigma = NULL, nobs = NULL, dist = NULL) {
	
}
 
