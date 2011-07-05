\name{trun.q}
\alias{trun.q}

\title{Truncated Inverse Cumulative Density Function of a gamlss.family Distribution}
\description{
  Creates a function to produce the inverse of a truncated cumulative density function generated from a current 
  GAMLSS family distribution
  
  For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value). 
}
\usage{
trun.q(par, family = "NO", type = c("left", "right", "both"), ...)
}

\arguments{
  \item{par}{ a vector with one (for \code{left} or \code{right} truncation) or two elements for \code{both}  }
  \item{family}{ a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   Functions such as BI() (binomial) produce a family object. }
  \item{type}{ whether \code{left}, \code{right} or in \code{both} sides truncation is required, (left is the default)  }
  \item{\dots}{for extra arguments }
}

\value{
Returns a q family function 
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).}
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk} }

\seealso{  \code{\link{trun.d}}, \code{\link{trun.q}}, \code{\link{trun.r}}, \code{\link{gen.trun}} }
\examples{
# continuous
#----------------------------------------------------------------------------------------
# left
test1<-trun.q(par=c(0), family="TF", type="left")
test1(.6)
qTF(pTF(0)+0.6*(1-pTF(0)))
#----------------------------------------------------------------------------------------
# right
test2 <- trun.q(par=c(10), family="BCT", type="right")
test2(.6)
qBCT(0.6*pBCT(10))
#----------------------------------------------------------------------------------------
# both
test3<-trun.q(par=c(-3,3), family="TF", type="both")
test3(.6)
qTF(0.6*(pTF(3)-pTF(-3))+pTF(-3))
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# FOR DISCRETE DISTRIBUTIONS
# trucated q function
# left
test4<-trun.q(par=c(0), family="PO", type="left")
test4(.6)
qPO(pPO(0)+0.6*(1-pPO(0)))
#----------------------------------------------------------------------------------------
# right
test5 <- trun.q(par=c(10), family="NBI", type="right")
test5(.6)
qNBI(0.6*pNBI(10))
test5(.6, mu=10, sigma=2)
qNBI(0.6*pNBI(10, mu=10, sigma=2), mu=10, sigma=2)
#----------------------------------------------------------------------------------------
# both
test6<-trun.q(par=c(0,10), family="NBI", type="both")
test6(.6)
qNBI(0.6*(pNBI(10)-pNBI(0))+pNBI(0))
#----------------------------------------------------------------------------------------
}
\keyword{distribution}
\keyword{regression}
