\name{effects.matrixpls}
\alias{effects.matrixpls}
\title{Total, Direct, and Indirect Effects for PLS model}
\usage{
  \method{effects}{matrixpls} (object = NULL, ...)
}
\arguments{
  \item{object}{PLS estimation result object produced by
  the \code{\link{matrixpls}} function.}

  \item{...}{All other arguments are ignored.}
}
\value{
  A list with \code{Total}, \code{Direct}, and
  \code{Indirect} elements.
}
\description{
  The \code{matrixpls} method for the standard generic
  function \code{effects} computes total, direct, and
  indirect effects for a PLS model according to the method
  described in Fox (1980).

  Adapted from the \code{\link[sem]{effects}} function of
  the \code{sem} package
}
\references{
  Fox, J. (1980) Effect analysis in structural equation
  models: Extensions and simplified methods of computation.
  \emph{Sociological Methods and Research} 9, 3--28.
}
\seealso{
  Other post-estimation functions: \code{\link{AVE}},
  \code{\link{CR}}, \code{\link{GoF}}, \code{\link{R2}},
  \code{\link{loadings.matrixpls}},
  \code{\link{residuals.matrixpls}}
}

