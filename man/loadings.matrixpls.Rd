\name{loadings.matrixpls}
\alias{loadings.matrixpls}
\title{Factor loadings matrix from PLS model}
\usage{
  loadings.matrixpls(x, ...)
}
\arguments{
  \item{x}{An object of class \code{matrixpls}}

  \item{...}{All other arguments are ignored.}
}
\value{
  A matrix of factor loadings.
}
\description{
  The \code{matrixpls} method for the generic function
  \code{loadings} computes standardized factor loading
  matrix from the model results.
}
\seealso{
  Other post-estimation functions: \code{\link{AVE}},
  \code{\link{CR}}, \code{\link{GoF}}, \code{\link{R2}},
  \code{\link{effects.matrixpls}},
  \code{\link{residuals.matrixpls}}
}

