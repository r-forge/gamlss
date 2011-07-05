\name{gen.trun}
\alias{gen.trun}

\title{Generates a Truncate Distribution from a gamlss.family}
\description{
  The \code{gen.trun()} function allows the user to generate the \code{d}, the \code{p}, the \code{q}, the \code{r} and the 
  \code{gamlss} fitting functions of a truncated  distribution given that the input \code{family} is a current GAMLSS family distribution. 

For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value).    
}
\usage{
gen.trun(par = c(0), family = "NO", name = "tr", type = c("left", "right", "both"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{ a vector with one (for \code{left} or \code{right} truncation) or two elements for \code{both} }
  \item{family}{ a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   Functions such as BI() (binomial) produce a family object.}
  \item{name}{the characters you want to add to the name of new truncated distribution, by default adds \code{tr}}
  \item{type}{whether \code{left}, \code{right} or in \code{both} sides truncation is required  }
  \item{\dots}{for extra arguments}
}

\value{
 Returns  the \code{d}, the \code{p}, the \code{q} and the \code{r} functions of 
  a truncated  distribution 
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk}}

\seealso{ \code{\link{trun.d}}, \code{\link{trun.p}}, \code{\link{trun.q}}, \code{\link{trun.r}} }
\examples{
gen.trun(par=c(0,100),family="TF", name="0to100", type="both")
plot(function(x) dTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100)
plot(function(x) pTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100)
plot(function(x) qTF0to100(x, mu=80 ,sigma=20, nu=5), 0.01, .999)
hist(s1<-rTF0to100(1000, mu=80 ,sigma=20, nu=5))
m1<-gamlss(s1~1, family=TF0to100)
}
\keyword{distribution}
\keyword{regression}
