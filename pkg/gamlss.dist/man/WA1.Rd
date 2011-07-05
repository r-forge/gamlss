\name{WA1}
\alias{WA1}
\alias{dWA1}
\alias{pWA1}
\alias{qWA1}
\alias{rWA1}
\title{
Waring distribution for fitting a GAMLSS
}
\description{
The function \code{WA1()} defines the Waring distribution, a two parameter distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting using the function \code{gamlss()}, with parameters \code{mu} and \code{sigma}. The functions \code{dWA1}, \code{pWA1}, \code{qWA1} and \code{rWA1} define the density, distribution function, quantile function and random generation for the \code{WA1} parameterization of the Waring distribution. 
}
\usage{

WA1(mu.link = "log", sigma.link = "log")
dWA1(x, mu = 2, sigma = 2, log = FALSE)
pWA1(q, mu = 2, sigma = 2, lower.tail = TRUE, log.p = FALSE)
qWA1(p, mu = 2, sigma = 2, lower.tail = TRUE, log.p = FALSE, max.value = 10000)
rWA1(n, mu = 2, sigma = 2)
}

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
  \item{max.value}{constant; generates a sequence of values for the cdf function}
}
\details{
The parameterization of the Pareto Type 2 distribution in the function \code{WA1} is
\deqn{f(y|mu, sigma) = Beta(sigma+y, mu+1)/Beta(sigma, mu)}
for \eqn{y>=0}, \eqn{mu>0} and \eqn{sigma>0}.
}
\value{
returns a gamlss.family object which can be used to fit a Waring distribution in the \code{gamlss()} function.}
\author{
Bob Rigby \email{r.rigby@londonmet.ac.uk}, Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Fiona McElduff \email{F.Mcelduff@londonmet.ac.uk} and Kalliope Akantziliotou
}
\references{
Rigby, R. A. and Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), \emph{Appl. Statist.}, 54, part 3, pp 507-554. 

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R. Accompanying documentation in the current GAMLSS help files, (see also \url{http://www.gamlss.com/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R. \emph{Journal of Statistical Software}, Vol. 23, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}. 

}

\seealso{
\code{\link{gamlss.family}}
}
\examples{
par(mfrow=c(2,2))
y<-seq(0,20,1)
plot(y, dWA1(y), type="h")
q <- seq(0, 20, 1)
plot(q, pWA1(q), type="h")
p<-seq(0.0001,0.999,0.05)
plot(p , qWA1(p), type="s")
dat <- rWA1(100)
hist(dat)
#summary(gamlss(dat~1, family="WA1"))
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
