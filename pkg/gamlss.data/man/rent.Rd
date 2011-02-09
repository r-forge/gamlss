\name{rent}
\alias{rent}
\docType{data}
\title{Rent data}
\description{
 A survey was conducted in April 1993 by Infratest Sozialforschung.
A random sample of accommodation with new tenancy agreements or
increases of rents within the last four years in Munich was
selected including: i) single rooms, ii) small apartments, iii)
flats, iv) two-family houses. Accommodation subject to price
control rents, one family houses and special houses, such as
penthouses, were excluded because they are rather different from
the rest and are considered a separate market. For the purpose of
this study, 1967 observations of the variables listed below were
used, i.e. the rent response variable R followed by the
explanatory variables found to be appropriate for a regression
analysis approach by Fahrmeir \emph{et al.}  (1994, 1995):

}
\usage{data(rent)}
\format{
  A data frame with 1969 observations on the following 9 variables.
  \describe{
    \item{R}{: rent response variable, the monthly net rent in DM, i.e.
                the monthly rent minus calculated or estimated utility cost}
    \item{Fl}{: floor space in square meters}
    \item{A}{: year of construction}
    \item{Sp}{: a variable indicating whether the location is above 
              average, \code{1}, (550 observations) or not, \code{0}, (1419 observations)}
    \item{Sm}{: a variable indicating whether the location is below, \code{1},
               average (172 obs.) or not, \code{0}, (1797 obs.)}
    \item{B}{: a factor with levels indicating whether there is a bathroom, \code{1}, (1925
                 obs.)  or not, \code{0}, (44 obs.) }
    \item{H}{: a factor with levels indicating whether there is central heating, \code{1},
              (1580 obs.)  or not, \code{0}, (389 obs.) }
    \item{L}{: a factor with levels indicating whether the kitchen equipment is
above average, \code{1},  (161 obs.) or not, \code{0}, (1808 obs.) }
    \item{loc}{: a factor (combination of Sp and Sm) indicating whether the location is below, \code{1}, average,  \code{2}, or above average \code{3}}
  }
}
\details{
  This set of data were used by Stasinopoulos \emph{et al.}  (2000) to fit a model where both the mean 
  and the dispersion parameter of a Gamma distribution were modelled using the explanatory variables.   
}
\source{Provide by Prof. L. Fahrmeir
  
}
\references{
 Fahrmeir L., Gieger C., Mathes H. and Schneeweiss H. (1994)
Gutachten zur Erstellung des Mietspiegels fur Munchen 1994, Teil
B: Statistiche Analyse der Nettomieten. Hrsg: Landeshaupttstadt
Munchen, Sozialreferat-Amt fur Wohnungswesen.

Fahrmeir L., Gieger C., and Klinger, A.  (1995) Additive,
dynamic and multiplicative regression. In \emph{ Applied Statistics:
Recent Developments}, Vandenhoeck and Ruprecht, Gottingen.

Stasinopoulos, D. M., Rigby, R. A. and Fahrmeir, L., (2000), 
Modelling rental guide data  using mean and dispersion additive models, \emph{Statistician}, \bold{49 },
479-493.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
 
\examples{
data(rent)
attach(rent)
plot(Fl,R)
}
\keyword{datasets}
