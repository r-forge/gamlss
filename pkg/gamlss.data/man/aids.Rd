\name{aids}
\alias{aids}
\docType{data}
\title{ Aids Cases in England and Wales }
\description{The quarterly reported AIDS cases in the U.K. from
January 1983 to March 1994 obtained from the Public Health
Laboratory Service, Communicable Disease Surveillance Centre,
London.  }
\usage{data(aids)}
\format{
  A data frame with 45 observations on the following 3 variables.
  \describe{
    \item{y}{the number of quarterly aids cases in England and Wales: a numeric vector}
    \item{x}{time in months from January 1983, 1:45 : a numeric vector}
    \item{qrt}{the quarterly seasonal effect a factor with 4 levels, [1=Q1 (Jan-March), 2=Q2 (Apr-June), 3=Q3 (July-Sept), 4=Q4 (Oct-Dec)]}
  }
}
\details{
  The counts y can be modelled using  a (smooth) Poisson regression  model in time x with the quarterly effects i.e. 
  cs(x,df=7)+qrt. Overdispersion persists,  so use a Negative Binomial distribution of type I or II.  
  The data also can be used to find a break point in time, see Rigby and Stasinopoulos (1992). 
}
\source{Public Health
Laboratory Service, Communicable Disease Surveillance Centre, London. 
}
\references{
 Stasinopoulos, D.M. and Rigby, R. A. (1992). Detecting break points in generalized linear models.
 \emph{Computational Statistics and Data Analysis}, \bold{13}, 461--471.  
}
\examples{
data(aids)
attach(aids)
plot(x,y,pch=21,bg=c("red","green3","blue","yellow")[unclass(qrt)])
detach(aids)
}
\keyword{datasets}
