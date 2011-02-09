\name{gamlss.cs}
\alias{gamlss.cs}
\alias{gamlss.scs}
\alias{gamlss.vc}


\title{Support for Function cs(), s() and vc()}
\description{ This is support for the functions cs(), s() and vc().
It is not intended to be called directly by users. The function \code{gamlss.cs} is based on the R function \code{smooth.spline}
}
\usage{
gamlss.cs(x, y = NULL, w = NULL, df = 5, spar = NULL, cv = FALSE,  nknots = NULL,
         all.knots = TRUE, df.offset = 0, penalty = 1, 
         control.spar = list(low = -1.5, high = 2), xeval = NULL)
gamlss.scs(x, y = NULL, w = NULL, df = NULL, spar = NULL, cv = FALSE, all.knots = FALSE, 
          nknots = NULL, keep.data = TRUE, df.offset = 0, penalty = 1.4, 
          control.spar = list(low = -1.5, high = 2), xeval = NULL)
gamlss.vc(x, y = NULL, w = NULL, df = 5, spar = NULL, cv = FALSE, all.knots = TRUE, 
          df.offset = 0, penalty = 1, control.spar = list(low = -1.5, high = 2), 
          xeval = NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{the design matrix}
  \item{y}{the response variable}
  \item{w}{ prior weights}
  \item{df}{ effective degrees of freedom }
  \item{spar}{spar the smoothing parameter}
  \item{cv}{options for the \code{smooth.spline} function not to use here }
  \item{nknots}{options for the \code{smooth.spline} function not to use here }
  \item{all.knots}{options for the \code{smooth.spline} function not to use here }
  \item{df.offset}{options for the \code{smooth.spline} function not to use here }
  \item{keep.data}{whether to keep data}
  \item{penalty}{options for the \code{smooth.spline} function not to use here }
  \item{control.spar}{control for spar.  It can be changed through \code{\link{cs}}}
  \item{xeval}{used in prediction } 
}

\value{
  Returns a class "smooth.spline" object with 
  \item{residuals}{The residuals of the fit}
  \item{fitted.values}{The smoothing values }
  \item{var}{the variance for the fitted smoother}
  \item{lambda}{the final value for spar}
  \item{nl.df}{the smoothing degrees of freedom excluding the constant and linear terms, i.e. (df-2)}
  \item{coefSmo}{this is a list containing among others the knots and the coefficients} 
  ...
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} }

\seealso{ \code{\link{gamlss}}, \code{\link{cs}} }

\keyword{regression}% 
