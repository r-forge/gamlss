\name{par.plot}
\alias{par.plot}

\title{A function to plot parallel plot for repeated measurement data}
\description{
  This function can be used to plot parallel plots for each individual in a repeated measurement study. 
  It is based on the \code{coplot()} function of \code{R}. 
}
\usage{
par.plot(formula = NULL, data = NULL, subjects = NULL, 
          color = TRUE, show.given = TRUE, ...)
}

\arguments{
  \item{formula}{ a formula describing the form of conditioning plot.  A
          formula of the form \code{y ~ x | a} indicates that plots of \code{y}
          versus \code{x} should be produced conditional on the variable
          \code{a}.  A formula of the form \code{y ~ x| a * b} indicates that
          plots of \code{y} versus \code{x} should be produced conditional on the
          two variables \code{a} and \code{b}.}
  \item{data}{ a data frame containing values for any variables in the
          formula.  By default the environment where \code{par.plot} was
          called from is used.}
  \item{subjects}{a factor which distinguish between the individual participants }
  \item{color}{whether the parallel plot are shown in colour, \code{color=TRUE} (the default) or not \code{color=FALSE}  }
  \item{show.given}{ logical (possibly of length 2 for 2 conditioning variables): 
         should conditioning plots be shown for the  corresponding conditioning variables (default 'TRUE')}
 \item{\dots}{for extra arguments}
}

\value{
It returns a plot. 
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.
 
 Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} }

\note{Note that similar plot can be fount in the library \code{nlme} by Pinheiro and Bates}

\seealso{\code{\link{gamlss}}}
\examples{

library(nlme)
data(Orthodont)
par.plot(distance~age,data=Orthodont,sub=Subject)
par.plot(distance~age|Sex,data=Orthodont,sub=Subject)
par.plot(distance~age|Subject,data=Orthodont,sub=Subject,show.given=FALSE)

}
\keyword{regression}%
