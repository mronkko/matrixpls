\name{matrixpls}
\alias{matrixpls}
\title{Partial Least Squares estimation}
\usage{
  matrixpls(S, model, W.mod = NULL,
    parameterEstimator = params.regression,
    outerEstimators = NULL, validateInput = TRUE, ...)
}
\arguments{
  \item{S}{Covariance matrix of the data.}

  \item{model}{There are two options for this argument: 1.
  lavaan script or lavaan parameter table, or 2. a list
  containing three matrices \code{inner},
  \code{reflective}, and \code{formative} defining the free
  regression paths in the model.}

  \item{W.mod}{An optional numeric matrix representing the
  weight model and starting weights (i.e. the how the
  indicators are combined to form the composite variables).
  If this argument is not specified, the weight model is
  defined based on the relationships in the
  \code{reflective} and \code{formative} elements of
  \code{model}.}

  \item{parameterEstimator}{A function that takes three
  arguments, the data covariance matrix \code{S}, weights
  \code{W}, and model specification \code{model} and
  returns a named vector of parameter estimates.}

  \item{...}{All other arguments are passed through to
  \code{matrixpls.weights} and \code{parameterEstimator}.}

  \item{outerEstimators}{A function or a list of functions
  used for outer estimation. If the value of this parameter
  is a function, the same function is applied to all
  composites. If the value of the parameter is a list, the
  composite \code{n} is estimated with the estimator in the
  \code{n}th position in the list. If this argument is
  \code{NULL} the starting weights specified in
  \code{W.mod} will be returned. The default is
  \code{\link{outer.modeA}} (PLS Mode A estimation).}

  \item{validateInput}{A boolean indicating whether the
  validity of the parameter values should be tested.}
}
\value{
  A named numeric vector of class \code{matrixpls}
  containing parameter estimates followed by weights.
}
\description{
  Estimates a weight matrix using Partial Least Squares or
  a related algorithm and then uses the weights to estimate
  the parameters of a statistical model.
}
\details{
  \code{matrixpls} is the main function of the matrixpls
  package. This function parses a model object and then
  uses the results to call \code{\link{matrixpls.weights}}
  to to calculate indicator weights. After this the
  \code{parameterEstimator} function is applied to the
  indicator weights, the data covariance matrix, and the
  model object and the resulting parameter estimates are
  returned.

  The model object specifies the model that is estimated
  with \code{parameterEstimator} and is used to construct
  \code{inner.mod} argument for
  \code{\link{matrixpls.weights}}. Model can be specified
  in the lavaan format or the native matrixpls format. The
  native model format is a list of three binary matrices,
  \code{inner}, \code{reflective}, and \code{formative}
  specifying the free parameters of a model: \code{inner}
  (\code{l x l}) specifies the regressions between
  composites, \code{reflective} (\code{k x l}) specifies
  the regressions of observed data on composites, and
  \code{formative} (\code{l x k}) specifies the regressions
  of composites on the observed data. Here \code{k} is the
  number of observed variables and \code{l} is the number
  of composites.

  If the model is specified in lavaan format, the native
  format model is derived from this model by assigning all
  regressions between latent variables to \code{inner}, all
  factor loadings to \code{reflective}, and all regressions
  of latent variables on observed variables to
  \code{formative}. Regressions between observed variables
  and all free covariances are ignored. All parameters that
  are specified in the model will be treated as free
  parameters. If model is specified in lavaan syntax, the
  model that is passed to the \code{parameterEstimator}
  will be that model and not the native format model.

  The original papers about Partial Least Squares, as well
  as many of the current PLS implementations, impose
  restrictions on the matrices \code{inner},
  \code{reflective}, and \code{formative}: \code{inner}
  must be a lower triangular matrix, \code{reflective} must
  have exactly one non-zero value on each row and must have
  at least one non-zero value on each column, and
  \code{formative} must only contain zeros. Some PLS
  implementations allow \code{formative} to contain
  non-zero values, but impose a restriction that the sum of
  \code{reflective} and \code{t(formative)} must satisfy
  the original restrictions of \code{reflective}. The only
  restrictions that matrixpls imposes on \code{inner},
  \code{reflective}, and \code{formative} is that these
  must be binary matrices and that the diagonal of
  \code{inner} must be zeros.

  The argument \code{W.mod} is a (\code{l x k}) matrix that
  indicates how the indicators are combined to form the
  composites. The original papers about Partial Least
  Squares as well as all current PLS implementations define
  this as \code{t(reflective) | formative}, which means
  that the weight model must match the model specified in
  \code{reflective} and \code{formative}. Matrixpls does
  not require that \code{W.mod} needs to match
  \code{reflective} and \code{formative}, but accepts any
  numeric matrix. If this argument is not specified,
  \code{W.mod} is defined as \code{t(reflective) |
  formative}.
}
\references{
  Rosseel, Y. (2012). lavaan: An R Package for Structural
  Equation Modeling. \emph{Journal of Statistical
  Software}, 48(2), 1–36. Retrieved from
  http://www.jstatsoft.org/v48/i02

  Wold, H. (1982). Soft modeling - The Basic Design And
  Some Extensions. In K. G. Jöreskog & S. Wold
  (Eds.),\emph{Systems under indirect observation:
  causality, structure, prediction} (pp. 1–54). Amsterdam:
  North-Holland.
}
\seealso{
  \code{\link{matrixpls.weights}}
}

