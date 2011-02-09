\name{acidity}
\alias{acidity}
\docType{data}
\title{The Acidity Data files for GAMLSS}
\description{
The data shows the acidity index for 155 lakes in the Northeastern United States
(previously analysed as a mixture of gaussian distributions on the log scale by Crawford \emph{et al.}(1992, 1994)). 
These 155 observations are the log acidity indices for the lakes.
}
\usage{data(acidity)}
\format{
  A data frame with 155 observations on the following variable.
  \describe{
    \item{\code{y}}{a numeric vector showing  the acidity index for 155 lakes in the Northeastern United States}
  }
}

\references{

Crawford S.L.,  DeGroot M.H.,  Kadane  J.B., and  Small M.J. (1992), Modeling lake-chemistry 
distributions: Approximate Bayesian methods for estimating a finite-mixture model,   \emph{Technometrics}, 
34, pp 441-450.

Crawford S.L. (1994) An application of the Laplace method to finite mixture distributions, \emph{JASA}, 89. pp 269-278.  


McLachlan G. and  Peel D., \emph{Finite Mixture Models}, Wiley, New York.

}
\examples{
data(acidity)
with( acidity, hist(y))
}
\keyword{datasets}
