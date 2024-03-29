\name{EGB2}
\alias{EGB2}
\alias{dEGB2}
\alias{pEGB2}
\alias{qEGB2}
\alias{rEGB2}


\title{The exponential generalized Beta type 2 distribution for fitting a GAMLSS}
\description{
This function defines the generalized t distribution, a four parameter distribution. The response variable is 
in the range from minus infinity to plus infinity.  
The functions \code{dEGB2}, 
\code{pEGB2}, \code{qEGB2} and \code{rEGB2} define the density, 
distribution function, quantile function and random
generation for the generalized beta type 2 distribution.   
}
\usage{
EGB2(mu.link = "identity", sigma.link = "identity", nu.link = "log", 
      tau.link = "log")
dEGB2(x, mu = 0, sigma = 1, nu = 1, tau = 0.5, log = FALSE)
pEGB2(q, mu = 0, sigma = 1, nu = 1, tau = 0.5, lower.tail = TRUE, 
      log.p = FALSE)
qEGB2(p, mu = 0, sigma = 1, nu = 0, tau = 0.5, lower.tail = TRUE, 
      log.p = FALSE)
rEGB2(n, mu = 0, sigma = 1, nu = 0, tau = 0.5)
}
%
\arguments{
 \item{mu.link}{Defines the \code{mu.link}, with "identity" link as the default for the \code{mu} parameter. }
 \item{sigma.link}{  Defines the  \code{sigma.link}, with "log" link as the default for the \code{sigma} parameter. }
 \item{nu.link}{Defines the  \code{nu.link}, with "log" link as the default for the \code{nu} parameter. }
 \item{tau.link}{Defines the  \code{tau.link}, with "log" link as the default for the \code{tau} parameter.}
 \item{x,q}{vector of quantiles}
 \item{mu}{vector of location parameter values }
 \item{sigma}{vector of scale parameter values}
 \item{nu}{vector of skewness \code{nu} parameter values}
 \item{tau}{vector of kurtosis \code{tau} parameter values}
 \item{log, log.p}{ logical; if TRUE, probabilities p are given as log(p).}
 \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x],
          otherwise, P[X > x] }
 \item{p}{vector of probabilities. }
 \item{n}{ number of observations. If \code{length(n) > 1}, the length is
          taken to be the number required}      
}
\details{   
 The probability density function of the Generalized Beta type 2, (\code{GB2}), is defined as 
  
\deqn{f(y|\mu,\sigma\,\nu,\tau)= e^{\mbox{\hspace{0.01cm}}\nu \mbox{\hspace{0.01cm}}z } \{|\sigma|\mbox{\hspace{0.05cm}} B(\nu,\tau) \mbox{\hspace{0.05cm}} [1+e^z]^{\nu+\tau}\}^{-1} }{f(y|mu,sigma,nu,tau)=exp{nu*z}(abs(sigma)*Beta(nu.tau)*(1+exp(z))^(nu+tau) )^(-1)}


for \eqn{-\infty<y<\infty}{-Inf<y<Inf}, where \eqn{z=(y-\mu)/\sigma}{z=(y-mu)/sigma} and 
\eqn{-\infty<\mu<\infty}{-Inf<mu<Inf}, \eqn{-\infty<\sigma<\infty}{-Inf<sigma<Inf}, 
\eqn{\nu>0}{nu>0} and \eqn{\tau>0}{tau>0},  McDonald and Xu (1995).
}
\value{
\code{EGB2()}  returns a \code{gamlss.family} object which can be used to fit the EGB2 distribution in the 
      \code{gamlss()} function. 
\code{dEGB2()} gives the density, \code{pEGB2()} gives the distribution
     function, \code{qEGB2()} gives the quantile function, and \code{rEGB2()}
     generates random deviates. 
}
\references{


Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R. A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).  

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Bob Rigby \email{r.rigby@londonmet.ac.uk} and Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} }



\seealso{  \code{\link{gamlss.family}}, \code{\link{JSU}}, \code{\link{BCT}}}
\examples{
EGB2()   # 
y<- rEGB2(200, mu=5, sigma=2, nu=1, tau=4)
library(MASS)
truehist(y)
fx<-dEGB2(seq(min(y), 20, length=200), mu=5 ,sigma=2, nu=1, tau=4)
lines(seq(min(y),20,length=200),fx)
# something funny here
# library(gamlss)
# histDist(y, family=EGB2, n.cyc=60)
integrate(function(x) x*dEGB2(x=x, mu=5, sigma=2, nu=1, tau=4), -Inf, Inf)
curve(dEGB2(x, mu=5 ,sigma=2, nu=1, tau=4), -10, 10, main = "The EGB2  density mu=5, sigma=2, nu=1, tau=4")
}
\keyword{distribution}
\keyword{ regression}
