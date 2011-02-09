\name{usair}
\alias{usair}
\docType{data}
\title{ US air pollution data set }
\description{US air pollution data set taken from  Hand et al. (1994) data set 26, USAIR.DAT, originally from 
Sokal and Rohlf (1981). 
}
\usage{data(usair)}
\format{
  A data frame with 41 observations on the following 7 variables.
  \describe{
    \item{y}{a numeric vector: sulpher dioxide concentration in air mgs. per cubic metre in 41 cities in the USA}
    \item{x1}{a numeric vector: average annual temperature in degrees F}
    \item{x2}{a numeric vector: number of manufacturers employing >20 workers}
    \item{x3}{a numeric vector: population size in thousands}
    \item{x4}{a numeric vector: average annual wind speed in miles per hour}
    \item{x5}{a numeric vector: average annual rainfall in inches}
    \item{x6}{a numeric vector: average number of days rainfall per year}
  }
}
\details{
}
\source{
 Hand et al. (1994) data set 26, USAIR.DAT, originally from 
 Sokal and Rohlf (1981)
}
\references{
Hand, D. J., Daly, F., Lunn, A. D., McConway, K. J. and Ostrowski, E. (1994),
A handbook of small data sets, Chapman and Hall, London.
}
\examples{
data(usair)
str(usair) 
plot(usair)
# a possible gamlss model
# gamlss(library)
#ap<-gamlss(y~cs(x1,2)+x2+x3+cs(x4,2)+x5+cs(x6,3)+x4:x5, 
#               data=usair, family=GA(mu.link="inverse"))
# 
}
\keyword{datasets}
