#'@section Model specification:
#'
#'Model can be specified in the lavaan format or the native matrixpls format.
#'The native model format is a list of three binary matrices, \code{inner}, \code{reflective},
#'and \code{formative} specifying the free parameters of a model: \code{inner} (\code{l x l}) specifies the 
#'regressions between composites, \code{reflective} (\code{k x l}) specifies the regressions of observed
#'data on composites, and \code{formative} (\code{l x k}) specifies the regressions of composites on the
#'observed data. Here \code{k} is the number of observed variables and \code{l} is the number of composites.
#'
#'If the model is specified in lavaan format, the native
#'format model is derived from this model by assigning all regressions between latent
#'variables to \code{inner}, all factor loadings to \code{reflective}, and all regressions
#'of latent variables on observed variables to \code{formative}. Regressions between
#'observed variables and all free covariances are ignored. All parameters that are
#'specified in the model will be treated as free parameters. If model is specified in
#'lavaan syntax, the model that is passed to the \code{parameterEstimator} will be that
#'model and not the native format model.
#'
#'The original papers about Partial Least Squares, as well as many of the current PLS
#'implementations, impose restrictions on the matrices \code{inner},
#'\code{reflective}, and \code{formative}: \code{inner} must be a lower triangular matrix,
#'\code{reflective} must have exactly one non-zero value on each row and must have at least
#'one non-zero value on each column, and \code{formative} must only contain zeros.
#'Some PLS implementations allow \code{formative} to contain non-zero values, but impose a
#'restriction that the sum of \code{reflective} and \code{t(formative)} must satisfy
#'the original restrictions of \code{reflective}. The only restrictions that matrixpls
#'imposes on \code{inner}, \code{reflective}, and \code{formative} is that these must be
#'binary matrices and that the diagonal of \code{inner} must be zeros.
