\name{get_plsr1}
\alias{get_plsr1}
\title{PLS regression for \code{matrixpls}}
\usage{
  get_plsr1(predictors, response, nc = NULL, scaled = TRUE)
}
\arguments{
  \item{predictors}{matrix of predictors}

  \item{response}{response variable}

  \item{nc}{number of components}

  \item{scaled}{logical indicating whether to scale the
  data}
}
\value{
  A list with pls regression results
}
\description{
  Internal function. \code{get_plsr1} is called by
  \code{matrixpls} to do PLS-R1.
}
\keyword{internal}

