\name{ING}
\alias{ING}
\alias{dING}
\alias{pING}
\alias{qING}
\alias{rING}
\title{
Inverse Gamma distribution for fitting a GAMLSS
}
\description{
The function \code{ING()} defines the Inverse Gamma distribution, a two parameter distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting using the function \code{gamlss()}, with parameters \code{mu} (the mode) and \code{sigma}. The functions \code{dING}, \code{pING}, \code{qING} and \code{rING} define the density, distribution function, quantile function and random generation for the \code{ING} parameterization of the Inverse Gamma distribution. 
}
\usage{

ING(mu.link = "log", sigma.link="log")
dING(x, mu = 2, sigma = 1, log = FALSE)
pING(q, mu = 2, sigma = 1, lower.tail = TRUE, log.p = FALSE)
qING(p, mu = 2, sigma = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1000)
rING(n, mu = 2, sigma = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{mu.link}{Defines the \code{mu.link}, with "log" link as the default for the mu parameter}
  \item{sigma.link}{Defines the \code{sigma.link}, with "log" as the default for the sigma parameter }
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
The parameterization of the Inverse Gamma distribution in the function \code{ING} is
\deqn{f(y|mu, sigma) = ([mu*(alpha+1)]^{alpha})/Gamma(alpha)*y^{(-(alpha+1))}*exp{-(mu*(alpha+1)}/y)}
where \eqn{alpha = 1/(sigma^2)}
for \eqn{y>0}, \eqn{mu>0} and \eqn{sigma>0}.
}
\value{
returns a gamlss.family object which can be used to fit an Inverse Gamma distribution in the \code{gamlss()} function.}
\author{
Bob Rigby \email{r.rigby@londonmet.ac.uk}, Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Fiona McElduff \email{F.Mcelduff@londonmet.ac.uk} and Kalliope Akantziliotou
}
\references{
Rigby, R. A. and Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), \emph{Appl. Statist.}, 54, part 3, pp 507--554. 

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R. Accompanying documentation in the current GAMLSS help files, (see also \url{http://www.gamlss.com/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R. \emph{Journal of Statistical Software}, Vol. 23, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}. 

}
\note{
For the function \code{ING()}, \emph{mu} is the mode of the Inverse Gamma distribution.
}

\seealso{
\code{\link{gamlss.family}}
}
\examples{
par(mfrow=c(2,2))
y<-seq(0.2,20,0.2)
plot(y, dING(y), type="l")
q <- seq(0.2, 20, 0.2)
plot(q, pING(q), type="l")
p<-seq(0.0001,0.999,0.05)
plot(p , qING(p), type="l")
dat <- rING(50)
hist(dat)
#summary(gamlss(dat~1, family="ING"))
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
