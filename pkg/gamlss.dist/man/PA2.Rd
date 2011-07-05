\name{PA2}
\alias{PA2}
\alias{dPA2}
\alias{pPA2}
\alias{qPA2}
\alias{rPA2}
\title{
Pareto Type 2 distribution for fitting a GAMLSS
}
\description{
The function \code{PA2()} defines the Pareto Type 2 distribution, a two parameter distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting using the function \code{gamlss()}, with parameters \code{mu} and \code{sigma}. The functions \code{dPA2}, \code{pPA2}, \code{qPA2} and \code{rPA2} define the density, distribution function, quantile function and random generation for the \code{PA2} parameterization of the Pareto type 2 distribution. 
}
\usage{

PA2(mu.link = "log", sigma.link = "log")
dPA2(x, mu = 2, sigma = 0.5, log = FALSE)
pPA2(q, mu = 2, sigma = 0.5, lower.tail = TRUE, log.p = FALSE)
qPA2(p, mu = 2, sigma = 0.5, lower.tail = TRUE, log.p = FALSE)
rPA2(n, mu = 2, sigma = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{mu.link}{Defines the \code{mu.link}, with "log" link as the default for the mu parameter}
  \item{sigma.link}{Defines the \code{sigma.link}, with "log" as the default for the sigma parameter}
  \item{x, q}{vector of quantiles}
  \item{mu}{vector of location parameter values}
  \item{sigma}{vector of scale parameter values}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p)}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], otherwise P[X > x]}
  \item{p}{vector of probabilities}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required}
}
\details{
The parameterization of the Pareto Type 2 distribution in the function \code{PA2} is
\deqn{f(y|mu, sigma) = (1/sigma)*mu^{(1/sigma)}*(y+mu)^{(-(1/sigma+1))}}
for \eqn{y>=0}, \eqn{mu>0} and \eqn{sigma>0}.
}
\value{
returns a gamlss.family object which can be used to fit a Pareto type 2 distribution in the \code{gamlss()} function.}
\author{
Bob Rigby \email{r.rigby@londonmet.ac.uk}, Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Fiona McElduff \email{F.Mcelduff@londonmet.ac.uk} and Kalliope Akantziliotou
}
\references{
Rigby, R. A. and Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), \emph{Appl. Statist.}, 54, part 3, pp 507--554. 

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R. Accompanying documentation in the current GAMLSS help files, (see also \url{http://www.gamlss.com/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R. \emph{Journal of Statistical Software}, Vol. 23, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}. 

}

\seealso{
\code{\link{gamlss.family}}
}
\examples{
par(mfrow=c(2,2))
y<-seq(0.2,20,0.2)
plot(y, dPA2(y), type="l" , lwd=2)
q<-seq(0,20,0.2)
plot(q, pPA2(q), ylim=c(0,1), type="l", lwd=2) 
p<-seq(0.0001,0.999,0.05)
plot(p, qPA2(p), type="l", lwd=2)
dat <- rPA2(100)
hist(rPA2(100), nclass=30)
#summary(gamlss(dat~1, family="PA2"))
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
