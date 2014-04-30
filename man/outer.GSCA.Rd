\name{outer.GSCA}
\alias{outer.GSCA}
\title{GSCA outer estimation}
\usage{
  outer.GSCA(S, W, E, W.mod, ...)
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
  This implements the second step of the GSCA estimation
  describe by Hwang & Takane (2004). GSCA outer estimation
  should be used only with GSCA inner estimation.
}
\details{
  The second step of GSCA estimation method, as describe by
  Hwang & Takane (2004), involves calculation of new
  weights given the regression estimates form the first
  step. In the second step, the following function is
  minimized (Hwang & Takane, 2004, eq. 7, first row):

  \deqn{SS(Z[V-\Lambda])}{SS(Z[V-Lambda])}

  Because \eqn{\Lambda}{Lambda} is defined as \eqn{WA}{WA},
  the function to be minimized is identical to the first
  step function (Hwang & Takane, 2004, eq. 4, first row):

  \deqn{SS(ZV-ZWA)}{SS(ZV-ZWA)}

  In the second step, this function is minimized in respect
  to weights \eqn{W}{W} and \eqn{V}{V}. This involves
  estimating each regression ananalysis in the model
  including regressions between the composites and from
  composites to indicators and to minimize the sum of all
  OLS dicrepancy functions simultaneously. Because one
  weight can be included in many regressions, these
  equations must be estimated simultaneously. The
  minimization algoritm is the Nelder-Mead algorithm
  implemented in the \code{\link{optim}} function.

  The GSCA algoritm described by Hwang & Takane (2004)
  allows some indicators to be excluded from the second
  step, but in this implementation all indicators are
  always used so that each weight relation described in
  \code{W.mod} has always a corresponding regression
  relationship from a composite to a variable in the second
  step of GSCA estimation.
}
\examples{
# Run the example from plspm package using GSCA estimation

if(require(plspm)) {
  
  # Run the customer satisfaction examle form plspm
  
  # load dataset satisfaction
  data(satisfaction)
  # inner model matrix
  IMAG = c(0,0,0,0,0,0)
  EXPE = c(1,0,0,0,0,0)
  QUAL = c(0,1,0,0,0,0)
  VAL = c(0,1,1,0,0,0)
  SAT = c(1,1,1,1,0,0)
  LOY = c(1,0,0,0,1,0)
  inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
  colnames(inner) <- rownames(inner)
  
  # Reflective model
  list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
  
  reflective<- matrix(
    c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
    27,6, dimnames = list(colnames(satisfaction)[1:27],colnames(inner)))
  
  # empty formative model
  
  formative <- matrix(0, 6, 27, dimnames = list(colnames(inner), colnames(satisfaction)[1:27]))
  
  # Estimation using covariance matrix and the GSCA estimators
  
  matrixpls(cov(satisfaction[,1:27]),  model = list(inner = inner,
                                                    reflective = reflective,
                                                    formative = formative),
            outerEstimators = outer.GSCA,
            innerEstimator = inner.GSCA)
  
} else{
  print("This example requires the plspm package")
}

}
\references{
  Hwang, H., & Takane, Y. (2004). Generalized structured
  component analysis. \emph{Psychometrika}, 69(1), 81–99.
  doi:10.1007/BF02295841
}
\seealso{
  Other outer estimators: \code{\link{outer.factor}},
  \code{\link{outer.fixedWeights}},
  \code{\link{outer.modeA}}, \code{\link{outer.modeB}}
}

