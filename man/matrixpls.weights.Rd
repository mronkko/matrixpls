\name{matrixpls.weights}
\alias{matrixpls.weights}
\title{Partial Least Squares weight algorithm}
\usage{
  matrixpls.weights(S, inner.mod, W.mod,
    outerEstimators = outer.modeA,
    innerEstimator = inner.path, ...,
    convCheck = function(W, W_old) {
        max(abs(W_old - W)) }, tol = 1e-05, iter = 100,
    validateInput = TRUE)
}
\arguments{
  \item{S}{Covariance matrix of the data.}

  \item{inner.mod}{A square matrix of ones and zeros
  representing the inner model (i.e. the path relationships
  between the composite variables).}

  \item{W.mod}{numeric matrix representing the weight model
  and starting weights (i.e. the how the indicators are
  combined to form the composite variables).}

  \item{outerEstimators}{A function or a list of functions
  used for outer estimation. If the value of this parameter
  is a function, the same function is applied to all
  composites. If the value of the parameter is a list, the
  composite \code{n} is estimated with the estimator in the
  \code{n}th position in the list. If this argument is
  \code{NULL} the starting weights specified in
  \code{W.mod} will be returned. The default is
  \code{\link{outer.modeA}} (PLS Mode A estimation).}

  \item{innerEstimator}{A function used for inner
  estimation. Setting this argument to \code{null} will use
  identity matrix as inner estimates and causes the
  algorithm to converge after the first iteration. This is
  useful when using \code{\link{outer.fixedWeights}} or
  some other outer estimation function that does not use
  inner estimation results. The default is
  \code{\link{inner.path}} (PLS path weighting scheme).}

  \item{convCheck}{A function that takes the old and new
  weight matrices and returns a scalar that is compared
  against \code{tol} to check for convergence. The default
  function calculates the differences between each cell of
  old and new weights and returns the largest absolute
  difference.}

  \item{tol}{Decimal value indicating the tolerance
  criterion for the iterations.}

  \item{iter}{An integer indicating the maximum number of
  iterations.}

  \item{validateInput}{A boolean indicating whether the
  validity of the parameter values should be tested.}

  \item{...}{All other arguments are passed through to
  \code{outerEstimators} and \code{innerEstimator}.}
}
\value{
  An object of class \code{"matrixplsweights"}, which is a
  matrix containing the weights with the following
  attributes:

  \item{iterations}{Number of iterations performed}.

  \item{converged}{A boolean indicating if the algorithm
  converged}.

  \item{history}{A data.frame containing the weights for
  each iteration}.
}
\description{
  Estimates a weight matrix using Partial Least Squares or
  a related algorithm.
}
\details{
  \code{matrixpls.weights} performs the PLS weighting
  algorithm by calling the inner and outer estimators
  iteratively until either the convergence criterion or
  maximum number of iterations is reached and provides the
  results in a matrix.

  The argument \code{inner.mod} is a binary (\code{l x l})
  matrix, where \code{l} is the number of composites. This
  matrix indicates how the composites should be linked in
  the inner estimation. The matrix \code{inner.mod} will
  contain a 1 when a variable on the column \code{j} has a
  link toward the variable on row \code{i} and 0 otherwise.

  The argument \code{W.mod} is a (\code{l x k}) matrix that
  indicates how the indicators are combined to form the
  composites, where \code{k} is the number of observed
  variables. \code{W.mod} will contain a non-zero value
  when a variable on the column \code{j} is a part of the
  composites variable on row \code{i}, zero otherwise. The
  values are used as starting values for the weighting
  algorithm.
}

