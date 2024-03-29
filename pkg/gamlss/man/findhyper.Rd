\name{find.hyper}
\alias{find.hyper}
\title{A function to select values of hyperparameters in a GAMLSS model }
\description{
 This function selects the values of hyper parameters and/or non-linear parameters in a GAMLSS model. It uses the \R{} function \code{optim} 
 which then minimised the generalized Akaike information criterion (GAIC) with a user defined penalty. 
}
\usage{
find.hyper(model = NULL, parameters = NULL, other = NULL, penalty = 2,
        steps = c(0.1), lower = -Inf, upper = Inf, method = "L-BFGS-B",
        ...)
}

\arguments{
  \item{model}{this is a GAMLSS model. e.g. \cr
                     \code{model=gamlss(y~cs(x,df=p[1]),sigma.fo=~cs(x,df=p[2]),data=abdom)} 
                     where p[1] and p[2] denore the parametets to be estimated}
  \item{parameters}{the starting values in the search of the optimum hyperparameters and/or non-linear parameters e.g. \code{parameters=c(3,3)}  }
  \item{other}{this is used to optimize other non-parameters, for example a transformation of the explanatory variable of the kind \eqn{x^{p[3]}}{x^p[3]},
      \code{others=quote(nx<-x^p[3])} where nx is now in the model formula }
  \item{penalty}{specifies the penalty in the GAIC, (the default is 2 e.g. \code{penalty=3}}
  \item{steps}{the steps taken in the optimization procedure [see the ndeps option in \code{optim()}], by default is set to 0.1 for all hyper parameters and non-linear parameters}
  \item{lower}{the lower permissible level of the parameters  i.e. \code{lower=c(1,1)} this does not apply if a method other than the default method "L-BFGS-B" is used }
  \item{upper}{the upper permissible level of the parameters  i.e. \code{upper=c(30,10)}, this is not apply if a method other than the default method "L-BFGS-B" is used }
  \item{method}{the method used in \code{optim()} to numerically minimize the GAIC over the hyperparameters and/or non-linear parameters.  By default this is "L-BFGS-B" to allow box-restriction on the parameters }
  \item{\dots}{for extra arguments to be passed to the \code{R} function \code{optim()} used in the optimization}
}
\details{
  This is an experimental function which appears to work well for the search of the optimum degrees of freedom and non-linear parameters 
  (e.g. power parameter \eqn{\lambda}{lambda} used to transform x to \eqn{x^\lambda}{x^lambda}). 
  Further investigation will check whether this function is reliable in general. 
        }
\value{
  The function turns the same output as the function \code{optim()}
  \item{par}{the optimum hyperparameter values}
  \item{value}{the minimized value of the GAIC}
  \item{counts}{A two-element integer vector giving the number of calls to `fn' and `gr' respectively }
  \item{convergence}{An integer code. `0' indicates successful convergence. see the function \code{optim()} for other errors}
  \item{message}{A character string giving any additional information returned by the optimizer, or `NULL'}
     }
\references{ 
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.


Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{ Mikis Stasinopoulos }


\section{Warning }{It may be slow to find the optimum} 
\seealso{ \code{\link{gamlss}}, \code{\link{plot.gamlss}}, \code{\link{optim}}  }
\examples{
data(abdom)
attach(abdom)
# declare the model
mod1<-quote(gamlss(y~cs(nx,df=p[1]),family=BCT,data=abdom,
                        control=gamlss.control(trace=FALSE)))
# we want also to check for a transformation in x 
# so we use the other option
op<-find.hyper(model=mod1, other=quote(nx<-x^p[2]), parameters=c(3,0.5), 
            lower=c(1,0.001), steps=c(0.1,0.001))
# the optimum parameters found are 
# p=(p[1],p[2]) = (2.944836 0.001000) = (df for mu, lambda)
# so it needs df = 3 on top of the constant and linear 
# in  the cubic spline model for mu since p[1] is approximately  3
# and log transformation for x since p[2] is approximately  0 
op
rm(op)
}
\keyword{regression}% 
