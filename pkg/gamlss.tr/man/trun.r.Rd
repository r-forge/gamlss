\name{trun.r}
\alias{trun.r}

\title{Generates Random Values from a Truncated Density Function of a gamlss.family Distribution}
\description{
 Creates a function to generate randon values from a truncated probability density function created from a current GAMLSS family distribution

For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value). 
}
\usage{
trun.r(par, family = "NO", type = c("left", "right", "both"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{par}{ a vector with one (for \code{left} or \code{right} truncation) or two elements for \code{both}  }
  \item{family}{ a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   Functions such as BI() (binomial) produce a family object. }
  \item{type}{ whether \code{left}, \code{right} or in \code{both} sides truncation is required, (left is the default)  }
  \item{\dots}{for extra arguments }
}

\value{
Returns a r family function  
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}). }
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk}  }

\seealso{ \code{\link{trun.p}}, \code{\link{trun.q}}, \code{\link{trun.d}}, \code{\link{gen.trun}} }
\examples{
# trucated r function
# continuous
#----------------------------------------------------------------------------------------
# left
test1<-trun.r(par=c(0), family="TF", type="left")
rr<-test1(1000)
hist(rr)
#----------------------------------------------------------------------------------------
# right
test2 <- trun.r(par=c(10), family="BCT", type="right")
rr<-test2(1000)
hist(rr)
#----------------------------------------------------------------------------------------
# both
test3<-trun.r(par=c(-3,3), family="TF", type="both")
rr<-test3(1000)
hist(rr)
#----------------------------------------------------------------------------------------
# discrete 
# trucated r function
# left
test1<-trun.r(par=c(0), family="PO", type="left")
tN <- table(Ni <- test1(1000))
 r <- barplot(tN, col='lightblue')
#----------------------------------------------------------------------------------------
# right
test2 <- trun.r(par=c(10), family="NBI", type="right")
tN <- table(Ni <- test2(1000))
 r <- barplot(tN, col='lightblue')
tN <- table(Ni <- test2(1000,mu=5))
  r <- barplot(tN, col='lightblue')
tN <- table(Ni <- test2(1000,mu=10, sigma=.1))
  r <- barplot(tN, col='lightblue')
#----------------------------------------------------------------------------------------
# both 
test3<-trun.r(par=c(0,10), family="NBI", type="both")
 tN <- table(Ni <- test3(1000,mu=5))
 r <- barplot(tN, col='lightblue')
}
\keyword{distribution}
\keyword{regression}
