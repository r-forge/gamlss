\name{disSmo}
\alias{disSmo}
\alias{disSmoA}
\alias{disSmoQ}
\alias{simRW}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Function to Fit Discrete Smoothing see for example Eilers( 2003)
}
\description{There are three function design to fit discrete Smoothers in time series data. 
 i)   \code{disSmoA()} for Alternate fitting when the smoothing parameters is estimated
 ii)  \code{disSmoQ()} for Q function fitting when the smoothing parameters is estimated
 iii) \code{disSmo} which combines the above two functions. Both of functions \code{disSmoA()} and \code{disSmoQ()}  should give identical results when lambda is estimated. The lambda can be fixed but only in function \code{sisSmoA()}. 
}
\usage{
disSmo(y, method = c("Alter", "Qfun"), ...)
disSmoA(y, weights = rep(1, length(y)), lambda = NULL, order = 1, 
        start = 10, plot = FALSE)
disSmoQ(y, weights = rep(1, length(y)), order = 1, start = 10, 
         plot = FALSE)
simRW(N = 1000, mu = 100, sige = 10, sigb = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the response variable usually a time series}
  \item{method}{method of finding the lambda, i) alternated, \code{Alter} ii) using the Q function \code{Qfun}. Note that both should give the   maximum likelihood estimator for code{lambda}}
  \item{weights}{prior weights}
  \item{lambda}{use this if you want to fix \code{lambda} (it works only with the \code{disSmoA()} function)}
  \item{order}{the required difference in the vector of coefficients see below}
  \item{start}{starting value for \code{lambda}}
  \item{plot}{whether to plot the data and the line after fit}
  \item{N}{the number of simulated observations using \code{simRW}}
  \item{mu}{The mean of the random walk simulated by \code{simRW}}
  \item{sige}{the standard deviation for the simulated data}
  \item{sigb}{the standard deviation for the random walk}
  \item{\dots}{for additional arguments}
}
\details{
The estimation of smoothing parameter \code{lambda} is done my maximum likelihood estimation.
The method in the function \code{disSmoA}  "alternates between estimates of \code{sigma}, (the variance of the data)  and \code{tau} 
(the variance of the random effects).
This method in  order to work   requires an estimate of the degrees of freedom of the model. The  degrees of freedom are the trace of the smoothing matrix   which is difficult to calculate for large number of observations. An approximation is used here within the \code{disSmoA)()} function.  For unweighted observations  This approximation is very accurate (compared with the Q function method which is exact).


}
\value{
Returns a fitted object of class \code{disSmo}. The object contains 
1)  the  \code{fitted.values} 
2) the effective degrees of freedom  \code{df} 
3) the smoothing parameter  \code{lambda},
4) the order of the difference \code{order},
5) \code{sige} the estimate of \code{sigma^2},
6) \code{sigma} the estimate of \code{sigma},
7) \code{tau2} the estimate of the random effect error,
8) the response variable \code{y}, 
9) \code{weights} the prior weights,
10) \code{N} the number of observations,
11) the call of the function \code{call},
12) the residual sum of squares   \code{rss}, 
13) the Akaike information criterion \code{aic},
14) the Bayesian information criterion \code{sbc} and 
15) the \code{deviance}
}

\references{   
Eilers, P. (2003). A perfect smoother. \emph{Analytical Chemistry}, 75: 3631-3636.
Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties (with comments and rejoinder). \emph{Statist. Sci},
\bold{11}, 89-121.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

}
\author{
Mikis Stasinopoulos, Paul Eilers, Bob Rigby and Majid  Djennad
}


\examples{
# test 500 obs -----------------------

y500 <- simRW(500)
bA500<- disSmo(y500, method="Alter")
bQ500<- disSmo(y500, method="Qfun")
a500 <- disSmoA(y500)
q500 <- disSmoQ(y500)

identical(fitted(a500), fitted(bA500))
identical(fitted(q500), fitted(bQ500))
identical(fitted(bA500), fitted(bQ500)) # false
plot(fitted(bA500), fitted(bQ500)) # but practivaly identical to 5 decimal places
nn<- abs(fitted(bA500)-fitted(bQ500))
 any(nn>0.00001) 
# very similar df's'
bA500$df
bQ500$df
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

