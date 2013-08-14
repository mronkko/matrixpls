\name{matrixpls-package}
\alias{matrixpls-package}
\docType{package}
\title{
  Matrix-based Partial Least Squares estimation
}
\description{
MatrixPLS does Partial Least Squares estimation using covariance data. The algorithms are based on matrix multiplication
and equation group solving (using the \code{\link{solve}} function).
}
\details{

The matrices are constructed in such a way that they can be interpreted as a set of multiple regression equations where
the variables on rows are regressed on variables on columns. 

	\tabular{lll}{
	Matrix \tab Order \tab Description\cr\cr
	\code{S} \tab k x k \tab Covariance matrix of the observed variables.\cr\cr
	\code{W} \tab l x k \tab Outer weight matrix. Returned by outer estimator.\cr\cr
	\code{E} \tab l x l \tab Inner weight matrix. Returned by inner estimator.\cr\cr
	         \tab       \tab Covariance matrix of the composites formed by combining\cr
	\code{C} \tab l x l \tab the observed variables into weighted composites using the\cr
	         \tab       \tab weight matrix \code{W}. Calculated as \code{W \%*\% S \%*\% t(W)}.\cr\cr
	         \tab       \tab Covariance matrix between composites and indicators.\cr
	\code{IC}\tab l x k \tab Calculated as \code{W \%*\% S} after outer estimation or\cr
	         \tab       \tab as \code{IC <- E \%*\% W \%*\% S} after inner estimation.\cr
	}

Where \code{k} is the number of observed variables and \code{l} is the number of latent variables. 

\tabular{ll}{
Package: \tab MatrixPLS\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2013-08-12\cr
License: \tab GPL-3\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{

}
