\name{garmaFit}
\alias{garmaFit}
\title{
A function to fit a GARMA model 
}
\description{
This function is fitting a GARMA model, see  Benjamin et al. (2003) for details. That is, a model in which the linear predictor of the mean (or location parameter) of a \code{gamlss.family} distribution  contains, apart from linear terms,  Autoregressive (AR) and/or Moving average (MA) terms. Note that the the model fitted here is more general than the one suggested in   Benjamin et al. (2003) in the sense that it is not restricted to the exponential family.  
}
\usage{
garmaFit(formula = formula(data), order = c(0, 0), 
         weights = NULL, data = sys.parent(), 
         family = NO(), alpha = 0.1, 
         phi.start = NULL, theta.start = NULL, 
         tail = max(order), control = list())
}

\arguments{
  \item{formula}{A formula for the linear terms i.e. like in linear model, \code{lm()}.
}
  \item{order}{
\code{order} specify the order of the generalised ARMA model i.e. \code{c(1,2)} is an AR 1 and MA 2 model. 
}
  \item{weights}{
prior weights, they are working like in \code{gamlss()} function.
}
  \item{data}{
the relevant \code{data.frame}.
}
  \item{family}{
A \code{gamlss.family} distribution.
}
  \item{alpha}{
This parameter is used in the definition of the link function of the response variable eg. for a log link use \eqn{log(y*)} will be \eqn{y*=max(y, \alpha)}.
}
  \item{phi.start}{
starting values for the AR parameters
}
  \item{theta.start}{
starting values for the MA part
}
  \item{tail}{
how many values, in the beginning of time, of the response variable should be weighted out from the fit (and used for initialising the process)   
}
  \item{control}{
control for \code{optim()}  or \code{nlminb()} function use for optimisation.
}
}
\details{
The model is described in  Benjamin et al. (2003). The implementation here is more general because it allows all the \code{gamlss.family} distributions to be fitted rather than only for the exponential family distributions which was described in the original formulation.
}
\value{
It returns a fitted \code{garma} model.
}
\references{
Benjamin M. A., Rigby R. A. and Stasinopoulos D.M. (2003) Generalised Autoregressive Moving  Average Models.  \emph{J. Am. Statist. Ass.}, 98, 214-223.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.com/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} and Vlasios Voudouris  
}
\note{
The is no check done whether the fitted model is stationary. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[gamlss.dist]{gamlss.family}}, \code{\link[gamlss]{gamlss}}
}
\examples{
data(polio)
ti <- as.numeric(time(polio))
mo <- as.factor(cycle(polio))
x1 <- 0:167    #Index used in Tutz p197
x2 <- cos(2*pi*x1/12)
x3 <- sin(2*pi*x1/12)
x4 <- cos(2*pi*x1/6)
x5 <- sin(2*pi*x1/6)
# all the data here 
da <-data.frame(polio,x1,x2,x3,x4,x5, ti, mo)
rm(ti,mo,x1,x2,x3,x4,x5)

#-------------------------------------------------------------------
# with linear trend 
m00 <- garmaFit(polio~x1+x2+x3+x4+x5, data=da, order=c(0,0), family=NBI, tail=3) # 
m10 <- garmaFit(polio~x1+x2+x3+x4+x5, data=da, order=c(1,0), family=NBI, tail=3) # OK
m01 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(0,1), data=da, family=NBI, tail=3)
m20 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(2,0), data=da, family=NBI, tail=3)
m11 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(1,1), data=da, family=NBI, tail=3)
m02 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(0,2), data=da, family=NBI, tail=3)
m30 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(3,0), data=da, family=NBI, tail=3)
m21 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(2,1), data=da, family=NBI, tail=3)
m12 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(1,2), data=da, family=NBI, tail=3)
m03 <- garmaFit(polio~x1+x2+x3+x4+x5, order=c(0,3), data=da, family=NBI, tail=3)
AIC(m00,m10,m01,m20,m11,m02,m30,m21,m12,m03 , k=0)
AIC(m00,m10,m01,m20,m11,m02,m30,m21,m12,m03 , k=log(168))

# without linear trend 
n00 <- garmaFit(polio~x2+x3+x4+x5, data=da, order=c(0,0), family=NBI, tail=3) # 
n10 <- garmaFit(polio~x2+x3+x4+x5, data=da, order=c(1,0), family=NBI, tail=3) # OK
n01 <- garmaFit(polio~x2+x3+x4+x5, order=c(0,1), data=da, family=NBI, tail=3)
n20 <- garmaFit(polio~x2+x3+x4+x5, order=c(2,0), data=da, family=NBI, tail=3)
n11 <- garmaFit(polio~x2+x3+x4+x5, order=c(1,1), data=da, family=NBI, tail=3)
n02 <- garmaFit(polio~x2+x3+x4+x5, order=c(0,2), data=da, family=NBI, tail=3)
n30 <- garmaFit(polio~x2+x3+x4+x5, order=c(3,0), data=da, family=NBI, tail=3)
n21 <- garmaFit(polio~x2+x3+x4+x5, order=c(2,1), data=da, family=NBI, tail=3)
n12 <- garmaFit(polio~x2+x3+x4+x5, order=c(1,2), data=da, family=NBI, tail=3)
n03 <- garmaFit(polio~x2+x3+x4+x5, order=c(0,3), data=da, family=NBI, tail=3)

AIC(m00,n10,n01,n20,n11,n02,n30,n21,n12,  k=0)
AIC(m00,n10,n01,n20,n11,n02,n30,n21,n12, k=log(168))}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
