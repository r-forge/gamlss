\name{gamlss.tr-package}
\alias{gamlss.tr-package}
\alias{gamlss.tr}
\docType{package}
\title{The GAMLSS Add On Package for Truncated Distributions}
\description{
The main purpose of this package is to allow the user of the GAMLSS models to fit truncated distributions.
}
\details{
\tabular{ll}{
Package: \tab gamlss.tr\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2005-11-21\cr
License: \tab GPL (version 2 or later) \cr
}
The user can take any \code{gamlss.family} and create a truncated version of it. Left, 
right and both sides truncation is allowed. 
For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value).    
}

\author{
Mikis Stasinopoulos <\email{d.stasinopoulos@londonmet.ac.uk}> and 
Bob Rigby <\email{r.rigby@londonmet.ac.uk}> 

Maintainer: Mikis Stasinopoulos <d.stasinopoulos@londonmet.ac.uk>

}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).
}
\keyword{ package }
\keyword{ distribution }
\keyword{ regression }

\examples{
gen.trun(par=c(0,100),family="TF", name="0to100", type="both")
op <- par(mfrow = c(2, 2)) 
plot(function(x) dTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100)
plot(function(x) pTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100)
plot(function(x) qTF0to100(x, mu=80 ,sigma=20, nu=5), 0.01, .999)
hist(rTF0to100(1000, mu=80 ,sigma=20, nu=5))
op
}
