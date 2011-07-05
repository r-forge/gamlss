\name{gamlss.ts-package}
\alias{gamlss.ts-package}
\alias{gamlss.ts}
\docType{package}
\title{
Time Series and GAMLSS
}
\description{
This package will contain several functions  fitting time series models. 
Currently only contains \code{garmaFit()} to fit a GARMA model, see Benjamin et al. (2003). [Other functions will appear soon. The idea is to be able to use 
techniques design for normally distributed time series data to all the \code{gamlss.family} distributions].   
}
\details{
\tabular{ll}{
Package: \tab gamlss.ts\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-04-28\cr
License: \\tab GPL (version 2 or later)\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Mikis Stasinopoulos <\email{d.stasinopoulos@londonmet.ac.uk}>, Bob Rigby <\email{r.rigby@londonmet.ac.uk}> 
Vlasis Voudouris

Maintainer: Mikis Stasinopoulos <d.stasinopoulos@londonmet.ac.uk>

}
\references{

Benjamin M. A., Rigby R. A. and Stasinopoulos D.M. (2003) Generalised Autoregressive Moving  Average Models.  \emph{J. Am. Statist. Ass.}, 98, 214-223.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

}

\keyword{ package }

\examples{
data(polio)
ti <- as.numeric(time(polio))
mo <- as.factor(cycle(polio))
x1 <- 0:167    #Index used in Tutz p197
x2 <- cos(2*pi*x1/12)
x3 <- sin(2*pi*x1/12)
x4 <- cos(2*pi*x1/6)
x5 <- sin(2*pi*x1/6)
# all the data here 
da <-data.frame(polio,x1,x2,x3,x4,x5, ti, mo)
m02 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(0,2), data=da, family=NBI, tail=3)
summary(m02)
}
