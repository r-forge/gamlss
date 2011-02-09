\name{aep}
\alias{aep}
\docType{data}
\title{ The Hospital Stay Data}
\description{
The data, 1383 observations, are from a study at the Hospital del
Mar, Barcelona during the years 1988 and 1990, Gange \emph{et
al.} (1996).
}
\usage{data(aep)}
\format{
  A data frame with 1383 observations on the following 8 variables.
  \describe{
    \item{los}{the total number of days patients spent in hospital: a discrete vector}
    \item{noinap}{the number of inappropriate days spent in hospital: a discrete vector }
    \item{loglos}{the log(los/10): a numeric vector}
    \item{sex}{the gender of patient: a factor with levels \code{1}=male, \code{2}=female}
    \item{ward}{the type of ward in the hospital: a factor with levels \code{1}=medical \code{2}=surgical, \code{3}=others}
    \item{year}{the specific year 1988 or 1990: a factor with levels \code{88} and  \code{90}}
    \item{age}{the age of the patient subtracted from 55: a numeric vector}
    \item{y}{the response variable a matrix with 2 columns, the first is noinap the second is equal to (los-noinap) }
  }
}
\details{
Gange \emph{et al.} (1996) used a logistic regression model for the number of inappropriate
days (noinap) out of the total number of days spent in hospital (los), with binomial and beta binomial errors and found that the later
provided a better fit to the data. They modelled both the mean and
the dispersion of the beta binomial distribution (BB) as functions
of explanatory variables  }
\source{ Gange, S. J. Munoz, A.  Saez, M. and  Alonso, J. (1996)}
\references{
 Gange, S. J. Munoz, A.  Saez, M. and  Alonso, J. (1996) Use of the beta-binomial distribution to model
        the effect of policy changes on appropriateness of hospital stays. \emph{Appl. Statist}, \bold{45}, 371--382
}
\examples{
data(aep)
attach(aep)
pro<-noinap/los
plot(ward,pro)
rm(pro)
detach(aep)
}
\keyword{datasets}
