\name{lo}
\alias{lo}
\title{ Specify a loess fit in a GAMLSS formula}
\description{
Allows the user to specify a loess fit in a GAMLSS formula. This function is similar to the \code{lo} function in the \code{gam} implementation of S-plus  
}
\usage{
lo(..., span = 0.5, df = NULL, degree = 1)
}

\arguments{
  \item{\dots}{ the unspecified ...  can be a comma-separated list of numeric vectors, numeric matrix, 
  or expressions that evaluate to either of these. If it is a list of vectors, they must all have the same length.
    }
  \item{span}{ the number of observations in a neighborhood. This is the smoothing parameter for a loess fit.
     }
  \item{df}{the effective degrees of freedom can be specified instead of span, e.g. \code{df=5}  }
  \item{degree}{the degree of local polynomial to be fit; can be 1 or 2. }
}
\details{
 Note that \code{lo} itself does no smoothing; it simply sets things up for the function \code{gamlss.lo()} which is used by the backfitting
 function  \code{gamlss.add()}.
}
\value{
 a numeric matrix is returned.  The simplest case is when there is a single argument to lo and degree=1; 
 a one-column matrix is returned, consisting of a normalized version of the vector.  
 If degree=2 in this case, a two-column matrix is returned, consisting of a 2d-degree orthogonal-polynomial basis.  
 Similarly, if there are two arguments, or the single argument is a two-column matrix, either a two-column matrix is returned if degree=1, 
 or a five-column matrix consisting of powers and products up to degree 2.  
 Any dimensional argument is allowed, but typically one or two vectors are used in practice. 
 The matrix is endowed with a number of attributes; the matrix itself is used in the construction of the model matrix, 
 while the attributes are needed for the backfitting algorithms all.wam or lo.wam (weighted additive model). 
 Local-linear curve or surface fits reproduce linear responses, while local-quadratic fits reproduce quadratic curves or surfaces.
 These parts of the loess fit are computed exactly together with the other parametric linear parts of the model. 
}
\references{Chambers, J. M. and Hastie, T. J. (1991). \emph{Statistical Models in S}, Chapman and Hall, London. 


Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}, based on the Trevor Hastie S-plus \code{lo} function }
\note{Note that \code{lo} itself does no smoothing; it simply sets things up for \code{gamlss.lo()} to do the backfitting.}

\section{Warning }{For user wanted to compare the \code{gamlss()} results with the equivalent \code{gam()} 
results in S-plus: make sure that the convergence criteria epsilon and bf.epsilon in S-plus are decreased sufficiently to ensure proper convergence in S-plus} 

\seealso{ \code{\link{cs}}, \code{\link{random}},  }

\examples{
data(aids)
attach(aids)
# fitting a loess curve with span=0.4 plus the a quarterly  effect 
aids1<-gamlss(y~lo(x,span=0.4)+qrt,data=aids,family=PO) # 
plot(x,y)
lines(x,fitted(aids1))
rm(aids1)
detach(aids)
}
\keyword{regression}% 
