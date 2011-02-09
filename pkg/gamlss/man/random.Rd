\name{random}
\alias{random}

\title{Specify a simple random effect in a GAMLSS Formula}
\description{ Includes random effect terms in an GAMLSS model.}
\usage{random(xvar, df = NULL, lambda = 0)}
\arguments{
  \item{xvar}{a factor }
  \item{df}{the target degrees of freedom}
  \item{lambda}{the smoothing parameter lambda which can be viewed as a shrinkage parameter.}
  }
\details{
This is an experimental smoother for use with factors in gamlss(). 
It allows the fitted values for a factor predictor to be shrunk towards the overall mean, 
where the amount of shrinking depends either on lambda, or on the equivalent degrees of freedom. 
Similar in spirit to smoothing splines, this fitting method can be justified on Bayesian grounds or by a random effects model.

Since factors are coded by model.matrix() into a set of contrasts, care has been taken to add an appropriate "contrast" 
attribute to the output of random(). This zero contrast results in a column of zeros in the model matrix, 
which is aliased with any column and is hence ignored
}
\value{
 x is returned with class "smooth", with an attribute named "call" which is to be evaluated in the backfitting  \code{additive.fit()} 
   called by \code{gamlss()}
}
\references{Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 


Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
            
\author{ Trevor Hastie (amended by Mikis Stasinopoulos)}


\seealso{\code{\link{gamlss}}, \code{\link{gamlss.random}}}

\examples{
data(aids)
attach(aids)
# fitting a loess curve with span=0.4 plus the a quarterly  effect 
aids1<-gamlss(y~lo(x,span=0.4)+qrt,data=aids,family=PO) # 
# now we string the quarterly  effect using random 
aids2<-gamlss(y~lo(x,span=0.4)+random(qrt,df=2),data=aids,family=PO) # 
plot(x,y)
lines(x,fitted(aids1),col="red")
lines(x,fitted(aids2),col="purple")
rm(aids1,aids2)
detach(aids)
}
\keyword{regression}% 
