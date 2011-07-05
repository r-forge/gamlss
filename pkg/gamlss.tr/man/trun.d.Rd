\name{trun.d}
\alias{trun.d}

\title{Truncated Probability Density Function of a gamlss.family Distribution}
\description{
  Creates a truncated probability density function version from a current GAMLSS family distribution
  
  For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value). 
}
\usage{
trun.d(par, family = "NO", type = c("left", "right", "both"), ...)
}

\arguments{
  \item{par}{a vector with one (for \code{left} or \code{right} truncation) or two elements for \code{both}   }
  \item{family}{a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   Functions such as BI() (binomial) produce a family object. }
  \item{type}{whether \code{left}, \code{right} or in \code{both} sides truncation is required, (left is the default).}
  \item{\dots}{ for extra arguments }
}

\value{
Returns a d family function 
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}). }
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk}  }

\seealso{  \code{\link{trun.p}}, \code{\link{trun.q}}, \code{\link{trun.r}}, \code{\link{gen.trun}} }
\examples{
# continuous distribution 
# left truncation 
test1<-trun.d(par=c(0), family="TF", type="left")
test1(1)
dTF(1)/(1-pTF(0))
if(abs(test1(1)-(dTF(1)/pTF(0)))>0.00001) stop("error in left trucation")
test1(1, log=TRUE)
log(dTF(1)/(1-pTF(0)))
if(abs(test1(1, log=TRUE)-log(dTF(1)/pTF(0)))>0.00001) stop("error in left trucation")
#if (!identical(test1(1, log=TRUE), log(dTF(1)/(1-pTF(0)))))   stop("error in left trucation")
 integrate(function(x) test1(x, mu=-2, sigma=1, nu=1),0,Inf) # the pdf is defined even with negative mu
 integrate(function(x) test1(x, mu=0, sigma=10, nu=1),0,Inf) 
 integrate(function(x) test1(x, mu=5, sigma=5, nu=10),0,Inf)
plot(function(x) test1(x, mu=-3, sigma=1, nu=1),0,10)
plot(function(x) test1(x, mu=3, sigma=5, nu=10),0,10)
#----------------------------------------------------------------------------------------
# right truncation
test2<-trun.d(par=c(10), family="BCT", type="right")
test2(1)
dBCT(1)/(pBCT(10))
#if(abs(test2(1)-(dBCT(1)/pBCT(10)))>0.00001) stop("error in right trucation")
test2(1, log=TRUE)
log(dBCT(1)/(pBCT(10)))
if(abs(test2(1, log=TRUE)-log(dBCT(1)/(pBCT(10))))>0.00001) stop("error in right trucation")
integrate(function(x) test2(x, mu=2, sigma=1, nu=1),0,10) 
integrate(function(x) test2(x, mu=2, sigma=.1, nu=1),0,10) 
integrate(function(x) test2(x, mu=2, sigma=.1, nu=10),0,10) 
plot(function(x) test2(x, mu=2, sigma=.1, nu=1),0,10)
plot(function(x) test2(x, mu=2, sigma=1, nu=1),0,10)
#----------------------------------------------------------------------------------------
# both left and right truncation
test3<-trun.d(par=c(-3,3), family="TF", type="both")
test3(0)
dTF(0)/(pTF(3)-pTF(-3))
if(abs(test3(0)-dTF(0)/(pTF(3)-pTF(-3)))>0.00001) stop("error in right trucation")
test3(0, log=TRUE)
log(dTF(0)/(pTF(3)-pTF(-3)))
if(abs(test3(0, log=TRUE)-log(dTF(0)/(pTF(3)-pTF(-3))))>0.00001) stop("error in both trucation")
plot(function(x) test3(x, mu=0, sigma=1, nu=1),-3,3)
integrate(function(x) test3(x, mu=2, sigma=1, nu=1),-3,3)
#----------------------------------------------------------------------------------------
# discrete distribution
# left 
# Poisson truncated at zero means zero is excluded
test4<-trun.d(par=c(0), family="PO", type="left")
test4(1)
dPO(1)/(1-pPO(0))
if(abs(test4(1)-dPO(1)/(1-pPO(0)))>0.00001) stop("error in left trucation")
test4(1, log=TRUE)
log(dPO(1)/(1-pPO(0)))
if(abs(test4(1, log=TRUE)-log(dPO(1)/(1-pPO(0))))>0.00001) stop("error in left trucation")
 sum(test4(x=1:20, mu=2)) # 
 sum(test4(x=1:200, mu=80)) #
plot(function(x) test4(x, mu=20), from=1, to=51, n=50+1, type="h") # pdf 
# right
# right truncated at 10 means 10 is excluded
test5<-trun.d(par=c(10), family="NBI", type="right")
test5(2)
dNBI(2)/(pNBI(9))
if(abs(test5(1)-dNBI(1)/(pNBI(9)))>0.00001) stop("error in right trucation")
test5(1, log=TRUE)
log(dNBI(1)/(pNBI(9)))
if(abs(test5(1, log=TRUE)-log(dNBI(1)/(pNBI(9))))>0.00001) stop("error in right trucation")
sum(test5(x=0:9, mu=2,   sigma=2)) # 
sum(test5(x=0:9, mu=300, sigma=5)) # can have mu > parameter
plot(function(x) test5(x, mu=20, sigma=3), from=0, to=9, n=10, type="h") # pdf
plot(function(x) test5(x, mu=300, sigma=5), from=0, to=9, n=10, type="h") # pdf
#----------------------------------------------------------------------------------------
# both
test6<-trun.d(par=c(0,10), family="NBI", type="both")
test6(2)
dNBI(2)/(pNBI(9)-pNBI(0))
if(abs(test6(2)-dNBI(2)/(pNBI(9)-pNBI(0)))>0.00001) stop("error in right trucation")
test6(1, log=TRUE)
log(dNBI(1)/(pNBI(9)-pNBI(0)))
if(abs(test6(1, log=TRUE)-log(dNBI(1)/(pNBI(9)-pNBI(0))))>0.00001) 
           stop("error in right trucation")
sum(test6(x=1:9, mu=2,   sigma=2)) # 
 sum(test6(x=1:9, mu=100, sigma=5)) # can have mu > parameter
plot(function(x) test6(x, mu=20, sigma=3), from=1, to=9, n=9, type="h") # pdf
plot(function(x) test6(x, mu=300, sigma=.4), from=1, to=9, n=9, type="h") # pdf
#----------------------------------------------------------------------------------------
}
\keyword{distribution}
\keyword{regression}
