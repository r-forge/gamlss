\name{dbbmi}
\alias{dbbmi}
\docType{data}
\title{BMI of Dutch Boys}
\description{
The data are comming from the Fourth Dutch Growth Study, Fredriks et al. (2000a,
2000b),  which is a cross-sectional study that measures growth and development of the
Dutch population between the ages 0 and 21 years. The study
measured, among other variables, height, weight, head
circumference and age for 7482 males and 7018 females. Here we have the only the BMI
of Dutch boys. 
}
\usage{data(dbbmi)}
\format{
  A data frame with 7294 observations on the following 2 variables.
  \describe{
    \item{\code{age}}{a numeric vector}
    \item{\code{bmi}}{a numeric vector}
  }
}
\source{ The data were kindly given by professor Stef. van Buuren. }
\references{ 
Fredriks, A.M. van Buuren, S. Burgmeijer, R.J.F. Meulmeester, J.F. Beuker, R.J. Brugman, E.
                   Roede, M.J. Verloove-Vanhorick, S.P.  and Wit, J. M. (2000a), Continuing 
                   positive secular change in The Netherlands, 1955-1997, \emph{Pediatric Research}, \bold{47}, 316--323

Fredriks, A.M. van Buuren, S. Wit, J.M. and Verloove-Vanhorick, S. P. (2000b) 
Body index measurments in 1996-7 compared with 1980,
\emph{Archives of Childhood Diseases}, \bold{82}, 107--112

van Buuren and Fredriks M. (2001) Worm plot: simple diagnostic device for modelling growth reference curves. 
            \emph{Statistics in Medicine}, \bold{20}, 1259--1277
}
\examples{
data(dbbmi)
plot(bmi~age, data=dbbmi)
}
\keyword{datasets}
