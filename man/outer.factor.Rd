\name{outer.factor}
\alias{outer.factor}
\title{Blockwise factor score outer estimation}
\usage{
  outer.factor(S, W, E, W.mod, fm = "minres", ...)
}
\arguments{
  \item{fm}{factoring method for estimating the common
  factor model. Possible values are \code{minres},
  \code{wls}, \code{gls}, \code{pa}, and \code{ml}. The
  parameter is passed through to to
  \code{\link[psych]{fa}}.}

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
  Calculates outer weights by estimating a common factor
  analysis model with a single factor for each indicator
  block and using the resulting estimates to calculate
  factor score weights
}
\seealso{
  Other outer estimators: \code{\link{outer.GSCA}},
  \code{\link{outer.fixedWeights}},
  \code{\link{outer.modeA}}, \code{\link{outer.modeB}}
}

