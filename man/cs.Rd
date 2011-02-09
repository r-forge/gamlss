\name{cs}

\alias{cs}
\alias{vc}
\alias{scs}

\title{Specify a Smoothing Cubic Spline Fit in a GAMLSS Formula}

\description{The functions \code{cs()} \code{scs()} and \code{vc} are using cubic smoothing splines. 
They take a vector and return it with several attributes. 
The vector is used in the construction of the model matrix.  The functions  do no smoothing, 
but assigns the attributes to the vector to aid gamlss in the smoothing.
The function doing the smoothing are \code{\link{gamlss.cs}()}, \code{gamlss.scs()} and \code{gamlss.vc()}.
The \code{\link{gamlss.cs}()}   is a modified version of the R function 
   \code{smooth.spline()} which is used by the backfitting function \code{\link{additive.fit}}.
The function  \code{gamlss.scs()}  just uses   \code{smooth.spline()}. It differs from the function \code{cs()} 
in that allows cross valdation of the smoothing parameters.   
The (experimental) function \code{vc} can be use to fit varying coefficient models, see Hastie and Tibshirani (1993).   
}
\usage{
cs(x, df = 3, spar = NULL, c.spar = NULL)
vc(r, x, df = 3, spar = NULL, c.spar = NULL)
scs(x, df = NULL, spar = NULL, control.spar = NULL, all.knots = TRUE, nknots = NULL, penalty =1.4)
}

\arguments{
  \item{x}{the univariate predictor, (or expression, that evaluates to a numeric vector). 
           For the function \code{vc} the x argument is the vector which has its (linear) coefficient change with \code{r}}
  \item{df}{the desired equivalent number of degrees of freedom (trace of the smoother matrix minus two for the constant and linear fit).
  The real smoothing parameter (spar below) is found such that df=tr(S)-2, 
  where S is the implicit smoother matrix. Values for df should be greater than 0, with 0 implying a linear fit.}
  \item{spar}{ smoothing parameter, typically (but not necessarily) in (0,1].   
               The coefficient lambda of the integral of the squared second derivative in the fit (penalized log likelihood)
               criterion is a monotone function of `spar', see the details in \code{smooth.spline}.}
  \item{c.spar}{This is an option to be used when the degrees of freedom of the fitted gamlss object are different from the ones given 
           as input in the option \code{df}.  The default values used are the ones given  the  option \code{control.spar}   
           in the R function \code{smooth.spine()} and they are \code{c.spar=c(-1.5, 2)}. 
   For  very large data sets e.g. 10000 observations, the upper limit may have to increase for example to  \code{c.spar=c(-1.5, 2.5)}.
   Use this option if you have received the warning 'The output df are different from the input, change the control.spar'.
   \code{c.spar} can take both vectors or lists of length 2, for example \code{c.spar=c(-1.5, 2.5)}
    or \code{c.spar=list(-1.5, 2.5)} would have the same effect.}
   \item{r}{for the function \code{vc}, \code{r} represent the vector of the explanatory variable which effects the coefficients of \code{x}
    i.e. beta(r)*x.  Both the \code{x} and \code{r} vectors should be adjusted by subtracting the their mean} 
    \item{control.spar}{see above \code{c.spar} or the equivalent argument in the function \code{smooth.spline}.}
    \item{all.knots}{if TRUE, all distinct points in x are used as knots. If FALSE (default), a subset of x is used, see the same argument for function 
      \code{smooth.spline}}
  \item{nknots}{integer giving the number of knots to use when all.knots=FALSE. Per default, this is less than n, the number of unique x values for n > 49}    
   \item{penalty}{the penalty applied to cross validation}  
}
\details{
Note that \code{cs} itself does no smoothing; it simply sets things up for the function \code{gamlss()} which in turn uses the function 
  \code{additive.fit()} for backfitting which in turn uses \code{gamlss.cs()}

Note that  \code{cs()} and  \code{s()} functions behave differently at their default values that is if df and lambda are not specified.
\code{cs(x)} by default will use 3 extra degrees of freedom for smoothing for \code{x}.
\code{ps(x)} by default will estimate lambda (and the degrees of freedom) automatically using generalised cross validation (GCV).
Note that if GCV is used the convergence of the gamlss model  can be  less stable
 compared to a model where the degrees of freedom are fixed. This can be true especially for small data sets.         
  }


\value{
  the vector x is returned, endowed with a number of attributes. The vector itself is used in the construction of the model matrix, 
  while the attributes are needed for the backfitting algorithms \code{additive.fit()}. 
  Since smoothing splines includes linear fits, the linear part will be efficiently computed with the other parametric linear parts of the model. 
}
\references{
Hastie, T. J. and Tibshirani, R. J. (1993), Varying coefficient models (with discussion),J. R. Statist. Soc. B., \bold{55},
    757-796.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{Mikis Stasinopoulos and Bob Rigby}

\note{The degrees of freedom df are defined differently from that of the gam() function in S-plus. Here df are the additional degrees of freedom 
excluding the constant and the linear part of x. For example \code{df=4} in \code{gamlss()} is equivalent to \code{df=5} in \code{gam()} in S-plus }

\section{Warning}{For a user who wishes to compare the \code{gamlss()} results with the equivalent \code{gam()} 
results in S-plus: make sure when using S-plus that the convergence criteria epsilon and bf.epsilon in \code{control.gam()} are decreased sufficiently
to ensure proper convergence in S-plus. 
Also note that the degrees of freedom are defined on top of the linear term in \code{gamlss}, but on top of the constant term in S-plus, 
(so use an extra degrees of freedom in S-plus in order to obtain comparable results to those in \code{galmss}).

Change the upper limit of spar if you received the warning 'The output df are different from the input, change the control.spar'.

For large data sets do not use expressions, e.g. \code{cs(x^0.5)} inside the \code{gamlss} function command but evaluate the expression, 
e.g. nx=\eqn{x^0.5}{x^0.5}, first and then use \code{cs(nx)}.  
} 

\seealso{ \code{\link{gamlss}}, \code{\link{gamlss.cs}}, \code{\link{lo}}  }

\examples{
# cubic splines example
data(aids)
# fitting a smoothing cubic spline with 7 degrees of freedom
# plus the a quarterly  effect  
aids1<-gamlss(y~cs(x,df=7)+qrt,data=aids,family=PO) # 
aids2<-gamlss(y~scs(x,df=5)+qrt,data=aids,family=PO) # 
aids3<-gamlss(y~scs(x)+qrt,data=aids,family=PO) # using GCV 
with(aids, plot(x,y))
lines(aids$x,fitted(aids1), col="red")
lines(aids$x,fitted(aids3), col="green")
rm(aids1, aids2, aids3)
# varying-coefficient example
data(rent)
attach(rent)
# adjusting the variables
Flbar<-Fl-mean(Fl)
Abar<-A-mean(A)
# additive model
 m1<-gamlss(R~cs(Flbar, df=3)+cs(Abar))
# varying-coefficient model
 m2<-gamlss(R~cs(Flbar, df=3)+cs(Abar)+vc(r=Abar,x=Flbar))
AIC(m1,m2)
detach(rent)
}

\keyword{regression}% 
