\name{vas5}
\alias{vas5}
\docType{data}
\title{Visual analog scale (VAS) data}
\description{
In the original data 368 patients, measured at 18 times after
treatment with one of 7 drug treatments (including placebo), plus
a baseline measure (time=0) and one or more pre-baseline measures
(time=-1). Here for illustration we will ignore the repeated measure nature of the
data and we shall use data from time 5 only (364 observations).
The VAS scale response variable, Y, is assumed to be distributed
as  \code{BEINF(mu,sigma,nu,tau)} where any of the
distributional parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau} are
modelled as a constant or as a function of the treatment,
}
\usage{data(vas5)}
\format{
  A data frame with 364 observations on the following 3 variables.
  \describe{
    \item{\code{patient}}{a factor indicationg the patient}
    \item{\code{treat}}{the treatment  factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7}}
    \item{\code{vas}}{the response variable}
  }
}
\details{
The Visual analog scale is used to measure pain and quality of
life. For example patients are required to indicate in a scale
from 0 to 100 the amount of discomfort they have. This can be
easily translated to a value from 0 to 1 and consequently analyzed
using the beta distribution. Unfortunately if 0's or 100's are
recorded the beta distribution is not appropriate since the values
0 and 1 are not allowed in the definition of the beta
distribution. Note that the inflated beta distribution
allows values at 0 and 1. This is a mixed distribution
(continuous and discrete) having four parameters, \code{nu} for
modelling the probability at zero p(Y=0) relative to p(0<Y<1), \code{tau} for modelling
the probability at one p(Y=1) relative to p(0<Y<1), and \code{mu} and \code{sigma} for
modelling the between values, $0<Y<1$, using a beta distributed
variable \code{BE(mu,sigma)} with mean \code{mu} and variance
\code{sigma*mu*(1-mu)}.
}
\source{
The data were provided by Dr. Peter Lane
}

\examples{
data(vas5)
}
\keyword{datasets}
