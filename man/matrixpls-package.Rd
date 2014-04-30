\name{matrixpls-package}
\alias{matrixpls-package}
\docType{package}
\title{
  Matrix-based Partial Least Squares estimation
}
\description{
Matrixpls estimates composite variable models using Partial Least Squares estimation and related methods. The algorithms are based on matrix multiplication
and equation group solving (using the \code{\link{solve}} function).
}
\details{

\tabular{ll}{
Package: \tab matrixpls\cr
Type: \tab Package\cr
Version: \tab 0.3.0\cr
Date: \tab 2014-04-28\cr
License: \tab GPL-3\cr
}

Matrixpls estimates models where sets of indicator variables are combined as weighted composites. These composites are then used to estimate a statistical model describing the relationships between the composites and composites and indicators. The most typical use for these method is approximating structural equation models
with latent variables.

The estimation algoritm consists of two steps. In the first step the indicators are combined as weighted composites and in the second step the composites are used to estimate a statistical model. The indicator weights are determined with the following iterative algorithm.

\tabular{ll}{
Inner      \tab Inner estimation function is applied to the data covariance matrix \code{S},\cr
estimation \tab weight matrix \code{W}, and composite variable model matrix \code{inner.mod}.\cr
					 \tab The function returns an inner weight matrix \code{E}.\cr\cr
Outer      \tab Outer estimation function is applied to the data covariance matrix \code{S},\cr
estimation \tab weight matrix \code{W}, inner weight matrix \code{W}, and weight model matrix \code{W.mod}. \cr
					 \tab The function returns a weight matrix \code{W}.\cr\cr
Convergence\tab Convergence check function is applied to the weight matrix \code{W} before and\cr
check      \tab after outer estimation. This function returns a scalar that is compared \cr
					 \tab against the tolerance value. If the scalar is smaller than the tolerance\cr
					 \tab value, the algorithm converges. Otherwise, a new iteration is started.\cr\cr
}

After the weights have been calculated, the parameter estimator function is applied to the data covariance matrix \code{S}, the weight matrix \code{W}, and model object.

The parameter matrices are constructed in such a way that they can be interpreted either as a set of multiple regression equations where
the variables on rows are regressed on variables on columns or as covariance matrices. 

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

Where \code{k} is the number of observed variables and \code{l} is the number of composites. 

}

\references{

Wold, H. (1982). Soft modeling - The Basic Design And Some Extensions. In K. G. Jöreskog & S. Wold (Eds.),\emph{Systems under indirect observation: causality, structure, prediction} (pp. 1–54). Amsterdam: North-Holland.

Dijkstra, T. K. (1983). Some comments on maximum likelihood and partial least squares methods. \emph{Journal of Econometrics}, 22(1-2), 67–90.

Lohmöller, J. B. (1989). \emph{Latent variable path modeling with partial least squares}. Physica-Verlag.

Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. \emph{Psychometrika}, 69(1), 81–99. doi:10.1007/BF02295841

Rönkkö, M., & Evermann, J. (2013). A Critical Examination of Common Beliefs About Partial Least Squares Path Modeling. \emph{Organizational Research Methods}, 16(3), 425–448. doi:10.1177/1094428112474693

}
\keyword{ package }
\seealso{
	\code{\link[plspm]{plspm}}, \code{\link[semPLS]{sempls}},
	\code{\link[simsem]{sim}}
}