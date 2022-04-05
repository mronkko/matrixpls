\name{mobi}
\Rdversion{1.1}
\alias{mobi}
\alias{ECSImobi}
\alias{ECSImm}
\alias{ECSIsm}
\docType{data}
\title{
  Measurement Instrument for the Mobile Phone Industry
}
\description{
  The data set is used as measurement instrument for the european customer
  satisfaction index (ECSI) adapted to the mobile phone market, see
  Tenenhaus et al. (2005). Adapted from semPLS package.
}
\usage{
data(mobi)
data(ECSImobi)
data(ECSImm)
data(ECSIsm)
}

\details{
  The data frame \code{mobi} contains the observed data for the model
  specified by \code{ECSImobi}.
  The from-to-matrizes \code{ECSImm} and code{ECSIsm} represent directed edges connecting the
  variables in outer/measurement and inner/structural model. They are
  needed for model specification by mean of the \code{plsm} method.
}

\format{
  The data frame \code{mobi} has 250 observations on the following 24 items. All the
  items are scaled from \eqn{1} to \eqn{10}.
  \describe{
    \item{\code{CUEX1}}{Expectations for the overall quality of \dQuote{your mobile phone provider} at the moment you became customer of this provider}
    \item{\code{CUEX2}}{Expectations for \dQuote{your mobile phone provider} to provide products and services to meet your personal need}
    \item{\code{CUEX3}}{How often did you expect that things could go wrong at \dQuote{your mobile phone provider}}
    \item{\code{CUSA1}}{Overall satisfaction}
    \item{\code{CUSA2}}{Fulfillment of expectations}
    \item{\code{CUSA3}}{How well do you think \dQuote{your mobile phone provider} compares with your ideal mobile phone provider?}
    \item{\code{CUSCO}}{You complained about \dQuote{your mobile phone provider} last year. How well, or poorly, was your most recent complaint handled or You did not complain about \dQuote{your mobile phone provider} last year. Imagine you have to complain to \dQuote{your mobile phone rovider} because of a bad quality of service or product. To what extent do you think that \dQuote{your mobile phone provider} will care about your complaint?}
    \item{\code{CUSL1}}{If you would need to choose a new mobile phone provider how likely is it that you would choose \dQuote{your provider} again?}
    \item{\code{CUSL2}}{Let us now suppose that other mobile phone providers decide to lower their fees and prices, but \dQuote{your mobile phone provider} stays at the same level as today. At which level of difference (in \%) would you choose another mobile phone provider?}
    \item{\code{CUSL3}}{If a friend or colleague asks you for advice, how likely is it that you would recommend \dQuote{your mobile phone provider}?}
    \item{\code{IMAG1}}{It can be trusted what it says and does}
    \item{\code{IMAG2}}{It is stable and firmly established}
    \item{\code{IMAG3}}{It has a social contribution to society}
    \item{\code{IMAG4}}{It is concerned with customers}
    \item{\code{IMAG5}}{It is innovative and forward looking}
    \item{\code{PERQ1}}{Overall perceived quality}
    \item{\code{PERQ2}}{Technical quality of the network}
    \item{\code{PERQ3}}{Customer service and personal advice offered}
    \item{\code{PERQ4}}{Quality of the services you use}
    \item{\code{PERQ5}}{Range of services and products offered}
    \item{\code{PERQ6}}{Reliability and accuracy of the products and services provided}
    \item{\code{PERQ7}}{Clarity and transparency of information provided}
    \item{\code{PERV1}}{Given the quality of the products and services offered by \dQuote{your mobile phone provider} how would you rate the fees and prices that you pay for them?}
    \item{\code{PERV2}}{Given the fees and prices that you pay for \dQuote{your mobile phone provider} how would you rate the quality of the products and services offered by \dQuote{your mobile phone provider}?}
  }
}

\references{
  Tenenhaus, M., V. E. Vinzi, Y.-M. Chatelin, and C. Lauro (2005)
  PLS path modeling.
  Computational Statistics & Data Analysis 48, 159-205.
}

\seealso{
  \code{\link{ECSImobi}}, \code{\link{ECSIsm}}, \code{\link{ECSImm}}
}

\examples{
data(mobi)
data(ECSImobi)
ecsi <- matrixpls.sempls(model=ECSImobi, data=mobi, E="C")
ecsi
}

\keyword{datasets}
