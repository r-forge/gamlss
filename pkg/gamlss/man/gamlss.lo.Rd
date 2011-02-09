\name{gamlss.lo}
\alias{gamlss.lo}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Support for Function lo()}
\description{
 This is support for the \code{loess} function \code{lo()}.
It is not intended to be called directly by users. The function \code{gamlss.lo} is based on the R function \code{loess}
}
\usage{
gamlss.lo(x, y, w = NULL, span, df = NULL, degree = 1, ncols =
                 FALSE, wspan = TRUE, parametric = FALSE, drop.square
                 = FALSE, normalize = FALSE, family = "gaussian",
                 method = "loess", control = loess.control(...), xeval
                 = NULL, ...)
}
\arguments{
  \item{x}{the design matrix}
  \item{y}{the response variable}
  \item{w}{ prior weights}
  \item{span}{the smoothing parameter}
  \item{df}{effective degrees of freedom }
  \item{degree}{the order of the polynomial}
  \item{ncols}{the number of columns of the x matrix}
  \item{wspan}{argument for the \code{loess} function not to use here}
  \item{parametric}{argument for the \code{loess} function not to use here}
  \item{drop.square}{argument for the \code{loess} function not to use here}
  \item{normalize}{argument for the \code{loess} function not to use here}
  \item{family}{argument for the \code{loess} function not to use here}
  \item{method}{argument for the \code{loess} function not to use here}
  \item{control}{argument for the \code{loess} function not to use here}
   \item{xeval}{used in prediction } 
  \item{...}{further arguments passed to or from other methods.}
}

\value{

 Returns an object with 
 
  \item{fitted}{the smooth values}
  \item{residuals}{the residuals}
  \item{var}{the variance of the smoother}
  \item{nl.df}{the non-linear degrees of freedom}
  \item{coefSmo}{with value NULL}
  \item{lambda}{the value of span}
 
}

\author{Mikis Stasinopoulos based on Brian Ripley loess function in R}


\seealso{ \code{\link{gamlss}}, \code{\link{lo}} }


\keyword{regression}% 
