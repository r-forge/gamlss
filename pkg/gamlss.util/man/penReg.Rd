\name{penReg}
\alias{penReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to fit penalised regression}
\description{
The function \code{penReg()} can be used to fit a P-spline. 
It can be used as demonstration of how the penalised B-splines can be fitted to one explanatory variable.
For more that explanatory variable use the function \code{pb()} in \pkg{gamlss}.  
}
\usage{
penReg(y, x, w = rep(1, length(y)), df = NULL, lambda = NULL, start = 10, 
      inter = 20, order = 2, degree = 3,  plot = FALSE,
      method = c("ML", "ML-1", "GAIC", "GCV", "EM"), k = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the response variable}
  \item{x}{the unique explanatory variable}
  \item{w}{prior weights}
  \item{df}{effective degrees of freedom}
  \item{lambda}{the smoothing parameter}
  \item{start}{the lambda starting value if the local methods are used}
  \item{inter}{the no of break points (knots) in the x-axis}
  \item{order}{the required difference in the vector of coefficients}
  \item{degree}{the degree of the piecewise polynomial}
  \item{plot}{whether to plot the data and  the fitted function}
  \item{method}{The method used in the (local) performance iterations. Available methods are "ML", "ML-1", "EM", "GAIC" and "GCV"}
  \item{k}{the penalty used in "GAIC" and "GCV"}
  \item{\dots}{for extra arguments}
}

\value{Returns a fitted object of class \code{penReg}. The object contains 1) the fitted 
\code{coefficients} 2) the   \code{fitted.values} 3) the response variable \code{y}, 
4) the label of the response variable \code{ylabel}
5) the explanatory variable \code{x}, 6) the lebel of the explanatory variable
7) the smoothing parameter  \code{lambda}, 8) the effective degrees of freedom \code{df},
9) the estimete for sigma \code{sigma},  
10) the residual sum of squares   \code{rss}, 11) the Akaike information criterion \code{aic},
12) the Bayesian information criterion \code{sbc} and 13) the \code{deviance} 
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
\seealso{
\code{\link{penLS}}
}

\examples{
set.seed(1234)
 x <- seq(0,10,length=200); y<-(yt<-1+2*x+.6*x^2-.1*x^3)+rnorm(200, 4)
 library(gamlss)
#------------------ 
# df fixed
 g1<-gamlss(y~pb(x, df=4))
 m1<-penReg(y,x, df=4) 
 cbind(g1$mu.coefSmo[[1]]$lambda, m1$lambda)
 cbind(g1$mu.df, m1$df)
 cbind(g1$aic, m1$aic)
 cbind(fitted(g1), fitted(m1))[1:10,]
# identical
#------------------
# estimate lambda using ML
g2<-gamlss(y~pb(x))
m2<-penReg(y,x) 
cbind(g2$mu.df, m2$df)
cbind(g2$mu.lambda, m2$lambda) 
cbind(g2$aic, m2$aic) # different lambda
cbind(fitted(g2), fitted(m2))[1:10,]
# identical
#------------------
#  estimate lambda using GCV
g3 <- gamlss(y~pb(x, method="GCV"))
m3 <- penReg(y,x, method="GCV") 
cbind(g3$mu.df, m3$df)
cbind(g3$mu.lambda, m3$lambda)
cbind(g3$aic, m3$aic)
cbind(fitted(g3), fitted(m3))[1:10,]
# almost identical
#------------------
# estimate lambda using EM
g4<-gamlss(y~pb(x, method="EM"))
m4<-penReg(y,x, method="EM") 
cbind(g4$mu.df, m4$df )
cbind(g4$mu.lambda, m4$lambda)
cbind(g4$aic, m4$aic)
cbind(fitted(g4), fitted(m4))[1:10,]
# almost identical
#------------------
#  estimate lambda using  GAIC(#=3)
g5<-gamlss(y~pb(x, method="GAIC", k=3))
m5<-penReg(y,x, method="GAIC", k=3) 
cbind(g5$mu.df, m5$df )
cbind(g5$mu.lambda, m5$lambda)
cbind(g5$aic, m5$aic)
cbind(g5$mu.df, m5$df)
cbind(g5$mu.lambda, m5$lambda)
cbind(fitted(g5), fitted(m5))[1:10,]
#-------------------
plot(y~x)
lines(fitted(m1)~x, col="green")
lines(fitted(m2)~x, col="red")
lines(fitted(m3)~x, col="blue")
lines(fitted(m4)~x, col="yellow")
lines(fitted(m4)~x, col="grey")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
