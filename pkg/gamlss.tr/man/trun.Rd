\name{trun}
\alias{trun}

\title{Fits a Truncate Distribution from a gamlss.family}
\description{
This function can be used to fit truncated distributions. It takes as an argument an existing GAMLSS family distribution and 
a parameter vector,  of the type c(left.value, right.value), and  generates a \code{gamlss.family} object which then can be used to fit 
 a truncated distribution.    
}
\usage{
trun(par = c(0), family = "NO",  type = c("left", "right", "both"), name = "tr", 
        local = TRUE, delta=NULL, ...)
}

\arguments{
  \item{par}{a scalar for left and right truncation or a vector of the type c(left.value, right.value) for interval truncation}
  \item{family}{an existing \code{gamlss.family} distribution}
  \item{type}{what type of truncation is required, \code{left}, \code{right} or \code{both}. If \code{both} the \code{par} should be a vector of length two. (the default is left truncation)}
  \item{name}{a character string to be added to name of the created object i.e. with \code{family=TF} and \code{name=trZero} the \code{gamlss.family} object will be called \code{TFtrZero}}
  \item{local}{if TRUE the function will try to find the environment of \code{gamlss} to generate the \code{d} and \code{p} functions required for the fitting, if FALSE the functions will be generated in the global environment}
  \item{delta}{the delta increment used in the numerical derivatives} 
  \item{\dots}{for extra arguments }
}
\details{
This function is created to help the user to fit a truncated form of existing \code{gamlss} distribution.
It does this by taking an existing \code{gamlss.family} and changing some of the components of the distribution to help the fitting process. 
It particular it i) creates a pdf (\code{d}) and a cdf (\code{p}) function within \code{gamlss}, 
ii) changes  the global deviance function \code{G.dev.incr}, the first derivative functions (see note below) and the quantile residual function.   
}
\value{
It returns a \code{gamlss.family} object which has all the components needed for fitting a distribution in \code{gamlss}. 
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk}}

\note{This function is experimental and could be changed. The function \code{trun}  changes 
      the first derivatives of the original gamlss family \code{d} function to numerical derivatives 
      for the new truncated \code{d} function. The default increment \code{delta},
      for this numerical derivatives function, is \code{eps * pmax(abs(x), 1)} where 
      \code{eps<-sqrt(.Machine$double.eps)}.  The default \code{delta} could be inappropriate 
      for specific applications and can be overwritten by using the argument \code{delta}.}
\seealso{\code{\link{trun.d}}, \code{\link{trun.p}}, \code{\link{trun.q}}, \code{\link{trun.r}},  \code{\link{gen.trun}} }
\examples{
# generate a left truncated zero t family
gen.trun(0,family="TF")
# take a random sample of 1000 observations
sam<-rTFtr(1000,mu=10,sigma=5, nu=5 )
hist(sam)
# fit the distribution to the data
mod1<-gamlss(sam~1, family=trun(0,TF))
mod1
# now create a gamlss.family object before the fitting 
Ttruc.Zero<- trun(par=0,family=TF, local=FALSE)
mod2<-gamlss(sam~1, family=Ttruc.Zero)
# now check the sensitivity of delta 
Ttruc.Zero<- trun(par=0,family=TF, local=FALSE, delta=c(0.01,0.01, 0.01))
mod3<-gamlss(sam~1, family=Ttruc.Zero)
}
\keyword{distribution}
\keyword{regression}
