\name{hodges}
\alias{hodges}
\alias{hodges1}
\docType{data}


\title{Hodges data}
\description{
There two data sets contain data used in Hodges (1998).  In addition 
to the data used in that manuscript, it contains other data items.

The original data consists of two matrices of dimensions of 341x6 and a 
45x4 respectively.
  
The first matrix \code{hodges} describes plans.  The information for each plan is:  
the state, a two-character code that identifies plans within state, the total premium for 
an individual, the total premium for a family, the total enrollment of 
federal employees as individuals, and the total enrollment of federal 
employees as families.

The second matrix, \code{hodges}, describes states.  The information for each state is:  
its two-letter abbreviation, the state average expenses per admission 
(from American Medical Association 1991 Annual Survey of Hospitals), 
population (1990 Census), and the region (from the Marion Merrill Dow 
Managed Care Digest 1991).  

The Hodges manuscript used these variables:
Plan level:  individual premium, individual enrollment. 
State level:  expenses per admission, region.  

}
\usage{data(hodges)}

\format{
  Two data frames the first with 341 observations on the following 6 variables.
  \describe{
    \item{\code{state}}{a factor with 45 levels \code{AL} \code{AZ} \code{CA} \code{CO} 
                                             \code{CT} \code{DC} \code{DE} \code{FL} 
                                             \code{GA} \code{GU} \code{HI} \code{IA} 
                                             \code{ID} \code{IL} \code{IN} \code{KS} 
                                             \code{KY} \code{LA} \code{MA} \code{MD} 
                                             \code{ME} \code{MI} \code{MN} \code{MO} 
                                             \code{NC} \code{ND} \code{NE} \code{NH} 
                                             \code{NJ} \code{NM} \code{NV} \code{NY} 
                                             \code{OH} \code{OK} \code{OR} \code{PA} 
                                             \code{PR} \code{RI} \code{SC} \code{TN} 
                                             \code{TX} \code{UT} \code{VA} \code{WA} 
                                             \code{WI}}
    \item{\code{plan}}{a two-character code that identifies plans within state declared here as factor with 325 levals.}
    \item{\code{prind}}{a numeric vector showing the total premium for an individual}
    \item{\code{prfam}}{a numeric vector showing the total premium for a family}
    \item{\code{enind}}{a numeric vector showing the total enrollment of federal employees as individuals}
    \item{\code{enfam}}{a numeric vector showing the total enrollment of federal employees as families.}
   }
   and the second with 45 observations on the following 4 variables
    \describe{
   \item{\code{State}}{a factor with levels same as state above }
    \item{\code{expe}}{a numeric vector showing the state average expenses per admission (from American Medical Association 1991 Annual Survey of Hospitals)}
    \item{\code{pop}}{a numeric vector shoing the population (1990 Census)}
    \item{\code{region}}{ the region (from the Marion Merrill Dow Managed Care Digest 1991), 
    a factor with levels \code{MA} \code{MT} \code{NC} \code{NE} \code{PA} \code{SA} \code{SC}}
  }
}

\source{
   \url{http://www.biostat.umn.edu/~hodges/}}
\references{
Hodges, J. S. (1998). Some algebra and geometry for hierarchical
models, applied to diadnostics. \emph{J. R. Statist. Soc. B.}, \bold{60} pp 497:536.
}
\examples{
data(hodges)
attach(hodges)
plot(prind~state, cex=1, cex.lab=1.5, cex.axis=1, cex.main=1.2)
str(hodges)
data(hodges1)
str(hodges1)
}
\keyword{datasets}
