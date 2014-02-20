\name{residuals.matrixpls}
\alias{residuals.matrixpls}
\title{Residual diagnostics for PLS model}
\usage{
  \method{residuals}{matrixpls} (object, ...)
}
\arguments{
  \item{object}{PLS estimation result object produced by
  the \code{\link{matrixpls}} function.}

  \item{...}{All other arguments are ignored.}
}
\value{
  A list with three elements: \code{inner}, \code{outer},
  and \code{indices} elements containing the residual
  covariance matrix of regressions of composites on other
  composites, the residual covariance matrix of indicators
  on composites, and various fit indices calculated based
  on the residuals.
}
\description{
  The \code{matrixpls} method for generic function
  \code{residuals} computes the residual covariance matrix
  and various fit indices presented by Lohmöller (1989, ch
  2.4). '
}
\references{
  Lohmöller J.-B. (1989) \emph{Latent variable path
  modeling with partial least squares.} Heidelberg:
  Physica-Verlag.
}
\seealso{
  Other post-estimation functions: \code{\link{AVE}},
  \code{\link{CR}}, \code{\link{GoF}}, \code{\link{R2}},
  \code{\link{effects.matrixpls}},
  \code{\link{loadings.matrixpls}}
}

