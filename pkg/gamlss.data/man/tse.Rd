\name{tse}
\alias{tse}
\docType{data}
\title{The Turkish stock exchange index}
\description{The Turkish stock exchange index, was recorded daily from 
1/1/1988 to 31/12/1998. 
The daily returns, \code{ret=log(I_(i+1)/I_(i))}, were obtained for i = 1,2,...,2868.}
\usage{data(tse)}
\format{
  A data frame with 2868 observations on the following 4 variables.
  \describe{
    \item{\code{year}}{the year}
    \item{\code{month}}{the month}
    \item{\code{day}}{the day}
    \item{\code{ret}}{day returns \code{ret[t]=ln(currency[t])-ln(currency[t-1])}}
    \item{\code{currency}}{the currency exchange rate}
    \item{\code{tl}}{day return \code{ret[t]=log10(currency[t])-log10(currency[t-1])}}
  }
}
\references{Ricard D. F. Harris and C. Coskun Kucukozen 
The Empirical Distribution of Stock returns: Evidence from a Emerging European Market,
Applied Economic Letters, 2001,8, pages 367-371.}

\examples{
data(tse)
plot(ts(tse$ret))
}
\keyword{datasets}
