\name{gamlssMX}
\alias{gamlssMX}
\alias{gamlssMXfits}

\title{Function to fit finite mixture of gamlss family distributions}
\description{
The  function \code{gamlssMX} is design for fitting a K fold non parametric mixture of gamlss family distributions.
}
\usage{
gamlssMX(formula = formula(data), pi.formula = ~1, 
         family = "NO", weights, K = 2, prob = NULL, 
         data = sys.parent(), control = MX.control(), 
         g.control = gamlss.control(trace = FALSE), 
         zero.component = FALSE,   ...)
gamlssMXfits(n = 5, formula = formula(data), pi.formula = ~1, 
         family = "NO", weights, K = 2, prob = NULL, 
         data = sys.parent(), control = MX.control(), 
         g.control = gamlss.control(trace = FALSE),
         zero.component = FALSE, ... )
}
\arguments{
  \item{formula}{This argument it should be a formula (or a list of formulea of length
  K) for modelling the \code{mu} parameter of the model. Note that
  modelling the rest of the distributional parameters it can be done
  by using the usual \code{\dots} which passes the arguments to
  \code{gamlss()} }
  \item{pi.formula}{This should be a formula for modelling the prior probabilities as a
  function of explanatory variables. Note that no smoothing of other
  additive terms are allowed here only the usual linear terms. The
  modelling here is done using the \code{multinom()} function from
  package \code{nnet}}
  \item{family}{This should be a \code{gamlss.family} distribution (or a list of
 distributions). Note that if different distributions are used here
 their parameters should be comparable for ease of interpretation.}
  \item{weights}{prior weights if needed}
  \item{K}{the number of finite mixtures with default \code{K=2} }
  \item{prob}{prior probabilities if required for starting values}
  \item{data}{the data frame nedded for the fit. Note that this is compulsory if \code{pi.formula} is used.}
  \item{control}{ This argument sets the control parameters for the EM iterations algorithm.
               The default setting are given in the \code{MX.control} function  }
  \item{g.control}{This argument can be used to pass to \code{gamlss()} control parameters, as in
  \code{gamlss.control} }
  \item{n}{the number of fits required in \code{gamlssMXfits()}}
  \item{zero.component}{whether zero component models exist, default is \code{FALSE}}
  \item{\dots}{ for extra arguments} 
}

\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,
(with discussion), \emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).
 }
\author{Mikis Stasinopoulos and Bob Rigby}

\seealso{\code{\link[gamlss]{gamlss}}, \code{\link[gamlss.dist]{gamlss.family}}}
\examples{
library(MASS)
data(geyser)
# fitting 2 finite normal mixtures 
m1<-gamlssMX(waiting~1,data=geyser,family=NO, K=2)
#fitting 2 finite gamma mixtures 
m2<-gamlssMX(waiting~1,data=geyser,family=GA, K=2)
# fitting a model for pi
# first create a data frame
geyser1<-matrix(0,ncol=2, nrow=298)
geyser1[,1] <-geyser$waiting[-1]
geyser1[,2] <-geyser$duration[-299]
colnames(geyser1)<- c("waiting", "duration")
geyser1 <-data.frame(geyser1)
# get the best of 5 fits
m3<-gamlssMXfits(n=5, waiting~1, pi.formula=~duration, data=geyser1,family=NO, K=2)
m3

}
\keyword{regression}

