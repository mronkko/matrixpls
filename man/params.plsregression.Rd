\name{params.plsregression}
\alias{params.plsregression}
\title{Parameter estimation with PLS regression}
\usage{
  params.plsregression(S, W, model)
}
\arguments{
  \item{S}{Covariance matrix of the data.}

  \item{W}{Weight matrix, where the indicators are on
  colums and composites are on the rows.}

  \item{model}{There are two options for this argument: 1.
  lavaan script or lavaan parameter table, or 2. a list
  containing three matrices \code{inner},
  \code{reflective}, and \code{formative} defining the free
  regression paths in the model.}
}
\value{
  A named vector of parameter estimates.
}
\description{
  Estimates the model parameters with weighted composites
  using separate OLS regressions for outer model and
  separate PLS regressions for inner model.
}
\details{
  \code{params.plsregression} estimates the model
  parameters similarly to \code{\link{params.regression}}
  with the exception that instead of separate OLS
  regressions the \code{inner} part of the model is
  estimated with separate PLS regressions using the PLS1
  algorithm with two rounds of estimation.

  The implementation of PLS regression is ported from the
  raw data version implemented in
  \code{\link[plspm]{get_plsr1}} funtion of the
  \code{plspm} package.
}
\references{
  Sanchez, G. (2013). \emph{PLS Path Modeling with R.}
  Retrieved from http://www.gastonsanchez.com/PLS Path
  Modeling with R.pdf

  Bjørn-Helge Mevik, & Ron Wehrens. (2007). The pls
  Package: Principal Component and Partial Least Squares
  Regression in R. \emph{Journal of Statistical Software},
  18. Retrieved from http://www.jstatsoft.org/v18/i02/paper
}
\seealso{
  Other parameter estimators: \code{\link{params.plsc}},
  \code{\link{params.regression}}
}

