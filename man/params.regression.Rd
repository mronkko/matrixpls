\name{params.regression}
\alias{params.regression}
\title{Parameter estimation with separate regression analyses}
\usage{
  params.regression(S, W, model, ...)
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

  \item{...}{All other arguments are passed through to
  \code{matrixpls.weights} and \code{parameterEstimator}.}
}
\value{
  A named vector of parameter estimates.
}
\description{
  Estimates the model parameters with weighted composites
  using separate OLS regressions.
}
\details{
  \code{params.regression} estimates the statistical model
  described by \code{model} with the following steps. If
  \code{model} is not in the native format, it is converted
  to the native format containing matrices \code{inner},
  \code{reflective}, and \code{formative}. The weights
  \code{W} and the data covariance matrix \code{S} are used
  to calculate the composite covariance matrix \code{C} and
  the indicator-composite covariance matrix \code{IC}.
  These are used to estimate multiple OLS regression
  models.

  The OLS regressions are estimated separately for each of
  the three model parts \code{inner}, \code{reflective},
  and \code{formative}. These matrices are analyzed one row
  at a time so that the row specifies the index of the
  dependent variable in the OLS regression and the non-zero
  elements on the row specify the indices of the
  independent variables.

  This approach of estimating the inner and outer models
  separately with separate OLS regression analyses is the
  standard way of estimation in the PLS literature.
}
\seealso{
  Other parameter estimators: \code{\link{params.plsc}},
  \code{\link{params.plsregression}}
}

