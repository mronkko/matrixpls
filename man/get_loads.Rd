\name{get_loads}
\alias{get_loads}
\title{Calculate loadings and communality for \code{matrixpls}}
\usage{
  get_loads(X, Y.lvs, blocks)
}
\arguments{
  \item{X}{Data matrix}

  \item{Y.lvs}{Matrix of latent variables}

  \item{blocks}{vector with number of manifest variabls per
  block}
}
\value{
  oadinsg and communalities
}
\description{
  Internal function. \code{get_loads} is called by
  \code{matrixpls}.
}
\keyword{internal}

