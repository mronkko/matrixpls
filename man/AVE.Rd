\name{AVE}
\alias{AVE}
\title{Average Variance Extracted indices for PLS model}
\usage{
  AVE(object, ...)
}
\arguments{
  \item{object}{PLS estimation result object produced by
  the \code{\link{matrixpls}} function.}

  \item{...}{All other arguments are ignored.}
}
\value{
  A list containing the Average Variance Extracted indices
  in the first position and the differences between AVEs
  and largest squared correlations with other composites in
  the second position.
}
\description{
  The \code{matrixpls} method for the generic function
  \code{AVE} computes Average Variance Extracted indices
  for the model using the formula presented by Fornell and
  Larcker (1981).
}
\references{
  Fornell, C., & Larcker, D. F. (1981). Evaluating
  structural equation models with unobservable variables
  and measurement error. \emph{Journal of marketing
  research}, 18(1), 39–50.
}
\seealso{
  Other post-estimation functions: \code{\link{CR}},
  \code{\link{GoF}}, \code{\link{R2}},
  \code{\link{effects.matrixpls}},
  \code{\link{loadings.matrixpls}},
  \code{\link{residuals.matrixpls}}
}

