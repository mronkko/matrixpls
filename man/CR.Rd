\name{CR}
\alias{CR}
\title{Composite Reliability indices for PLS model}
\usage{
  CR(object, ...)
}
\arguments{
  \item{object}{PLS estimation result object produced by
  the \code{\link{matrixpls}} function.}

  \item{...}{All other arguments are ignored.}
}
\value{
  A named numeric vector containing the Composite
  Reliability indices.
}
\description{
  The \code{matrixpls} method for the generic function
  \code{CR} computes Composite Reliability indices for the
  model using the formula presented by Fornell and Larcker
  (1981).
}
\references{
  Fornell, C., & Larcker, D. F. (1981). Evaluating
  structural equation models with unobservable variables
  and measurement error. \emph{Journal of marketing
  research}, 18(1), 39–50.
}
\seealso{
  Other post-estimation functions: \code{\link{AVE}},
  \code{\link{GoF}}, \code{\link{R2}},
  \code{\link{effects.matrixpls}},
  \code{\link{loadings.matrixpls}},
  \code{\link{residuals.matrixpls}}
}

