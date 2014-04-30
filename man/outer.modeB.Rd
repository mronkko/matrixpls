\name{outer.modeB}
\alias{outer.modeB}
\title{PLS outer estimation with Mode B}
\usage{
  outer.modeB(S, W, E, W.mod, ...)
}
\arguments{
  \item{S}{Covariance matrix of the data.}

  \item{W}{Weight matrix, where the indicators are on
  colums and composites are on the rows.}

  \item{E}{Inner weight matrix. A square matrix of inner
  estimates between the composites.}

  \item{W.mod}{A matrix specifying the weight relationships
  and their starting values.}

  \item{...}{Other parameters are ignored}
}
\value{
  A matrix of unscaled outer weights \code{W} with the same
  dimesions as \code{W.mod}.
}
\description{
  Performs Mode B outer estimation by regressing the
  composites on the observed variables.
}
\seealso{
  Other outer estimators: \code{\link{outer.GSCA}},
  \code{\link{outer.factor}},
  \code{\link{outer.fixedWeights}},
  \code{\link{outer.modeA}}
}

