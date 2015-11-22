#'@details
#'
#'The argument \code{W.model} is a (\code{l x k}) matrix that indicates
#'how the indicators are combined to form the composites. The original papers about
#'Partial Least Squares as well as all current PLS implementations define this as
#'\code{t(reflective) | formative}, which means that the weight patter must match the
#'model specified in \code{reflective} and \code{formative}. Matrixpls does not
#'require that \code{W.model} needs to match \code{reflective} and \code{formative}, but
#'accepts any numeric matrix. If this argument is not specified, all elements of \code{W.model} that
#'correspond to non-zero elements in the \code{reflective} or \code{formative} formative
#'matrices receive the value 1.
