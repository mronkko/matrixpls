\name{inner.centroid}
\alias{inner.centroid}
\title{PLS inner estimation with the centroid scheme}
\usage{
  inner.centroid(S, W, inner.mod, ...)
}
\arguments{
  \item{S}{Covariance matrix of the data.}

  \item{W}{Weight matrix, where the indicators are on
  colums and composites are on the rows.}

  \item{inner.mod}{A square matrix specifying the
  relationships of the composites in the model.}

  \item{...}{Other parameters are ignored}
}
\value{
  A matrix of unscaled inner weights \code{E} with the same
  dimesions as \code{inner.mod}.
}
\description{
  Calculates a set of inner weights based on the centroid
  scheme.
}
\details{
  In the centroid scheme, inner weights are set to the
  signs (1 or -1) of correlations between composites that
  are connected in the model specified in \code{inner.mod}
  and zero otherwise.

  Falls back to to identity scheme for composites that are
  not connected to any other composites.
}
\references{
  Lohmöller J.-B. (1989) \emph{Latent variable path
  modeling with partial least squares.} Heidelberg:
  Physica-Verlag.
}
\seealso{
  Other inner estimators: \code{\link{inner.GSCA}},
  \code{\link{inner.factor}}, \code{\link{inner.identity}},
  \code{\link{inner.path}}
}

