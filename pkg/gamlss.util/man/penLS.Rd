\name{penLS}
\alias{penLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit penalised least squares}
\description{
The function \code{penLS()} can be used to fit a penalised least square to a response variable y. 
There is no explanatory variable here. 
The underline model is a random walk.  
}
\usage{
penLS(y, w = rep(1, length(y)), df = NULL, lambda = NULL, 
      start = 10, order = 1, plot = FALSE, 
      type = c("level", "trend"), 
      method = c("ML", "GAIC", "GCV"), k = 2, ...)
}
\arguments{
  \item{y}{the response variable usually a time series}
  \item{w}{prior weights if needed otherwise 1}
  \item{df}{effective degrees of freedom}
  \item{lambda}{the smoothing parameter}
  \item{start}{the lambda starting value if the local methods are use}
  \item{order}{the required difference in the vector of coefficients, see below}
  \item{plot}{whether to plot the data and  the fitted function}
  \item{type}{the type of X matrix, if "level" X is a diagonel matrix of 1's if "trend" X is a 
   diagonal matrix 1:n}
  \item{method}{The method used in the estimation of the smoothing parameter}
 \item{k}{the penalty used in "GAIC" and "GCV"}
  \item{\dots}{for extra arguments}
}
\details{
The order refers to differences in the penalty matrix, (i)  order = 0 : white noise random effects
(ii) order = 1 : random walk (iii) order = 2 : random walk of order 2
(iv) order = 3 : random walk of order 3 
}
\value{
Returns a fitted object of class \code{penLS}. The object contains 
1) the fitted \code{coefficients} 
2) the  \code{fitted.values} 
3) the response variable \code{y}, 
4) the smoothing parameter  \code{lambda}, 8) the effective degrees of freedom \code{df},
5) the estimete for sigma \code{sigma},  
6) the residual sum of squares   \code{rss}, 
7) the Akaike information criterion \code{aic},
8) the Bayesian information criterion \code{sbc} and 
9) the \code{deviance} 
}
\references{
Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties (with comments and rejoinder). \emph{Statist. Sci},
\bold{11}, 89-121.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk } and Paul Eilers
}
\section{Warning }{The function can be slow for large response variable} 

\seealso{
\code{\link{penReg}}
}
\examples{
set.seed(1234)
 x<-seq(0,10,length=200); y<-(yt<-1+2*x+.6*x^2-.1*x^3)+rnorm(200,0, 16)
yts <- ts(y)
plot(yts)
#---------------------------------------------------------
#lambda fix
 m1<-penLS(yts,lambda=1) ; deviance(m1) 
#--------------------------------------------------------- 
# fixing df 
 m2<-penLS(yts, df=10) ; deviance(m2)
#---------------------------------------------------------  
# estimating lambda - ML
m3<-penLS(yts) ; deviance(m3)
#---------------------------------------------------------
# estimating lambda - GAIC
m4<-penLS(yts, method="GAIC", k=3) ; deviance(m4)
#---------------------------------------------------------
# different order
PPP <- par(mfrow=c(2,2))
penLS(yts, plot=TRUE, order=0, main="order=0")
penLS(yts, plot=TRUE, order=1, main="order=1")
penLS(yts, plot=TRUE, order=2, main="order=2")
penLS(yts, plot=TRUE, order=3, main="order=3")
par(PPP) 
}
\keyword{ regression }
