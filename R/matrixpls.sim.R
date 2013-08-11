#'@title Monte Carlo simulations with MatrixPLS
#'
#'@description
#'A wrapper function to use MatrixPLS with SimSem
#
#'@details
#
# This function
#
# All the parameters are identical to SimSem sim function

library(simsem)


matrixpls.sim <- function(nRep, model, ...)

    # A function that takes a data set and returns a list. The list must
	# contain at least three objects: a vector of parameter estimates (coef),
	# a vector of standard error (se), and the convergence status as TRUE or
	# FALSE (converged). There are five optional objects in the list: a vector
	# of fit indices (fit), a vector of standardized estimates (std), any
	# extra output (extra), fraction missing type I (FMI1), and fraction
	# missing type II (FMI2). Note that the coef, se, std, FMI1, and FMI2 must
	# be a vector with names. The name of those vectors across different
	# objects must be the same.

}
 
