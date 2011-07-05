\name{trun.p}
\alias{trun.p}

\title{Truncated Cumulative Density Function of a gamlss.family Distribution}
\description{
 Creates a truncated cumulative density function version from a current GAMLSS family distribution.
 
 For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value). 
}
\usage{
trun.p(par, family = "NO", type = c("left", "right", "both"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{ a vector with one (for \code{left} or \code{right} truncation) or two elements for \code{both}  }
  \item{family}{ a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   Functions such as BI() (binomial) produce a family object. }
  \item{type}{ whether \code{left}, \code{right} or in \code{both} sides truncation is required, (left is the default)}
  \item{\dots}{for extra arguments }
}

\value{
   Return a p family function 
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}). }
\author{  Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk}  }


\seealso{  \code{\link{trun.d}}, \code{\link{trun.q}}, \code{\link{trun.r}}, \code{\link{gen.trun}}}
\examples{
# trucated p  continuous function
# continuous
#----------------------------------------------------------------------------------------
# left
test1<-trun.p(par=c(0), family="TF", type="left")
test1(1)
(pTF(1)-pTF(0))/(1-pTF(0))
if(abs(test1(1)-(pTF(1)-pTF(0))/(1-pTF(0)))>0.00001) stop("error in left trucation of p")
plot(function(x) test1(x, mu=2, sigma=1, nu=2),0,10)
#----------------------------------------------------------------------------------------
# right
test5 <- trun.p(par=c(10), family="BCT", type="right")
test5(1)
pBCT(1)/pBCT(10)
if(abs(test5(1)-pBCT(1)/pBCT(10))>0.00001) stop("error in right trucation")
test5(1, lower.tail=FALSE)
1-pBCT(1)/pBCT(10)
if(abs(test5(1, lower.tail=FALSE)-(1-pBCT(1)/pBCT(10)))>0.00001) stop("error in right trucation")
test5(1, log.p=TRUE)
log(pBCT(1)/pBCT(10))
if(abs(test5(1, log.p=TRUE)-log(pBCT(1)/pBCT(10)))>0.00001) stop("error in right trucation")
plot(function(x) test5(x, mu=2, sigma=1, nu=2, tau=2),0,10)
plot(function(x) test5(x, mu=2, sigma=1, nu=2, tau=2, lower.tail=FALSE),0,10)
#----------------------------------------------------------------------------------------
# both 
test3<-trun.p(par=c(-3,3), family="TF", type="both")
test3(1)
(pTF(1)-pTF(-3))/(pTF(3)-pTF(-3))
if(abs(test3(1)-(pTF(1)-pTF(-3))/(pTF(3)-pTF(-3)))>0.00001) stop("error in right trucation")
test3(1, lower.tail=FALSE)
1-(pTF(1)-pTF(-3))/(pTF(3)-pTF(-3))
if(abs(test3(0,  lower.tail=FALSE)-(1-(pTF(0)-pTF(-3))/(pTF(3)-pTF(-3))))>0.00001) 
           stop("error in right trucation")
plot(function(x) test3(x, mu=2, sigma=1, nu=2, ),-3,3)
plot(function(x) test3(x, mu=2, sigma=1, nu=2, lower.tail=FALSE),-3,3)
#----------------------------------------------------------------------------------------
# Discrete
#----------------------------------------------------------------------------------------
# trucated p function
# left
test4<-trun.p(par=c(0), family="PO", type="left")
test4(1)
(pPO(1)-pPO(0))/(1-pPO(0))
if(abs(test4(1)-(pPO(1)-pPO(0))/(1-pPO(0)))>0.00001) stop("error in left trucation of p")
plot(function(x) test4(x, mu=2), from=1, to=10, n=10, type="h")
cdf <- stepfun(1:40, test4(1:41, mu=5), f = 0)
plot(cdf, main="cdf", ylab="cdf(x)", do.points=FALSE )
#----------------------------------------------------------------------------------------
# right
test2<-trun.p(par=c(10), family="NBI", type="right")
test2(2)
pNBI(2)/(pNBI(9))
if(abs(test2(2)-(pNBI(2)/(pNBI(9))))>0.00001) stop("error in right trucation of p")
plot(function(x) test2(x, mu=2), from=0, to=9, n=10, type="h")
cdf <- stepfun(0:8, test2(0:9, mu=5), f = 0)
plot(cdf, main="cdf", ylab="cdf(x)", do.points=FALSE )
#----------------------------------------------------------------------------------------
# both 
test6<-trun.p(par=c(0,10), family="NBI", type="both")
test6(2)
(pNBI(2)-pNBI(0))/(pNBI(9)-pNBI(0))
if(abs(test6(2)-(pNBI(2)-pNBI(0))/(pNBI(9)-pNBI(0)))>0.00001) stop("error in the both trucation")
test6(1, log=TRUE)
log((pNBI(1)-pNBI(0))/(pNBI(9)-pNBI(0)))
if(abs(test6(1, log=TRUE)-log((pNBI(1)-pNBI(0))/(pNBI(9)-pNBI(0))))>0.00001) stop("error in both trucation")
plot(function(y) test6(y, mu=20, sigma=3), from=1, to=9, n=9, type="h") # cdf
plot(function(y) test6(y, mu=300, sigma=.4), from=1, to=9, n=9, type="h") # cdf
cdf <- stepfun(1:8, test6(1:9, mu=5), f = 0)
plot(cdf, main="cdf", ylab="cdf(x)", do.points=FALSE )
#----------------------------------------------------------------------------------------  
}
\keyword{distribution}
\keyword{regression}
