\name{VictimsOfCrime}
\alias{VictimsOfCrime}
\docType{data}
\title{Reported victims of crime  data}
\description{
The data shows whether victims of crime were reported in the local media. 
}
\usage{data(VictimsOfCrime)}
\format{
  A data frame with 10590 observations on the following 2 variables.
  \describe{
    \item{\code{reported}}{ Whether the crime was reported in local media.}
    \item{\code{age}}{the age of the victim}
  }
}
\details{
 Whether the crime was reported in local media.
}
\source{
The data were given  by Prof Brian Francis of Lancaster University. 
They can be used to demonstrate the usefulness of smoothing techniques with a binary response variable. 
}

\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\examples{
data(VictimsOfCrime)
}
\keyword{datasets}
