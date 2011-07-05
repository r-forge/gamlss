\name{dMX}
\alias{dMX}
\alias{pMX}
\title{Evaluate the d (pdf) and p (cdf) functions from GAMLSS mixtures}
\description{
The functions \code{dMX} and \code{pMX} can be used to evalaluate the 
pdf (p function) and the cdf (p function) repectively from a gamlss.family mixture.
}   
\usage{
dMX(y, mu = list(mu1 = 1, mu2 = 5), sigma = list(sigma1 = 1, sigma2 = 1), 
       nu = list(nu1 = 1, nu2 = 1), tau = list(tau1 = 1, tau2 = 1), 
       pi = list(pi1 = 0.2, pi2 = 0.8), family = list(fam1 = "NO", fam2 = "NO"), 
       log = FALSE, ...)
pMX(q, mu = list(mu1 = 1, mu2 = 5), sigma = list(sigma1 = 1, sigma2 = 1), 
       nu = list(nu1 = 1, nu2 = 1), tau = list(tau1 = 1, tau2 = 1), 
       pi = list(pi1 = 0.2, pi2 = 0.8), family = list(fam1 = "NO", fam2 = "NO"), 
       log = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y,q}{vector of quantiles }
  \item{mu}{a vector of \code{mu}'s}
  \item{sigma}{a vector of \code{sigma}'s}
  \item{nu}{a vector of \code{nu}'s }
  \item{tau}{a vector of \code{tau}'s }
  \item{pi}{ a vector of \code{pi}'s }
  \item{family}{a vector of GAMLSS \code{family}'s }
  \item{log}{whether the  \code{log} of the function or not}
  \item{\dots}{for extra arguments}
}

\value{
Returns values or pdf  or cdf. 
}

\author{Mikis Stasinopoulos}
\examples{
fyNO<-dMX(y=seq(0,3,.01), mu=list(1.253, 0.1876), sigma=list(exp(-0.6665 ), exp(-2.573 )),
                  pi=list(0.4079609, 0.5920391 ), family=list("NO","NO") )
plot(fyNO~seq(0,3,.01), type="l")                  
FyNO<-pMX(q=seq(0,3,.01), mu=list(1.253, 0.1876), sigma=list(exp(-0.6665 ), exp(-2.573 )),
                  pi=list(0.4079609, 0.5920391 ), family=list("NO","NO") )
plot(FyNO~seq(0,3,.01), type="l")  
}
\keyword{distribution}
\keyword{regression}
