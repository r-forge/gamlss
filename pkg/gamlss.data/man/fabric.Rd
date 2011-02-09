\name{fabric}
\alias{fabric}
\docType{data}
\title{ The Fabric Data }
\description{The data are 32 observations on faults in rolls of fabric
}
\usage{data(fabric)}
\format{
  A data frame with 32 observations on the following 3 variables.
  \describe{
    \item{leng}{the length of the roll : a  numeric vector}
    \item{y}{the number of faults in the roll of fabric  : a  discrete vector}
    \item{x}{the log of the length of the roll : a  numeric vector}
  }
}
\details{
 The data are 32 observations on faults in rolls of fabric taken from Hinde (1982) 
 who used the EM algorithm to fit a Poisson-normal model. 
 The response variable is the number of faults in the roll of fabric and the explanatory variable is the log of the length of the roll. 
}
\source{ John Hinde
  
}
\references{
  Hinde,  J. (1982) Compound Poisson regression models: in \emph{GLIM} 82, \emph{Proceedings of the
    International Conference on Generalized Linear Models}, ed. Gilchrist, R., 109--121, Springer: New York.

}
\examples{
data(fabric)
attach(fabric)
plot(x,y)
detach(fabric)
}
\keyword{datasets}
