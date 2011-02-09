\name{Mums}
\alias{Mums}
\docType{data}
\title{Mothers encouragement data}
\description{Mothers encouragement for participation in Higher Education. The response variable is \code{mums} a three level factor 
which can be used in a multinomial Logistic model or mumsB a two level factor suitable for binary logistic model.  
}
\usage{data(Mums)}
\format{
  A data frame with 871 observations on the following 7 variables.
  \describe{
    \item{mums}{mothers encouragement: factor with levels \code{1} is for strong encouragement, 
                                     \code{2} is for some encouragement and 
                                     \code{3} for no encouragement/discouragement} 
    \item{class}{social class: a factor with levels \code{1}is C1, \code{2} is C2, \code{3} is D and   \code{4} is E}
    \item{age}{age of the participants: a factor with levels \code{1} is 16-18,  \code{2} is 19-20 and \code{3} is 20-30}
    \item{gender}{a factor with levels \code{1} is male  and \code{2} is female}
    \item{ethn}{ethnicity of the participants: a factor with levels \code{1} is white, \code{2} is black,  \code{3} is asian and  \code{4} is other}
    \item{qual}{qualifications of the participants: a factor with levels, \code{1} is greater or equal to 2 A levels, 
                                                                          \code{2} is HND or more than 5 GCSE's,
                                                                          \code{3} is less than 5 GSCSE's ar none above and
                                                                          \code{4} no formal qualification}
    \item{mumsb}{mothers encouragement: a factor with levels, \code{0} is no encouragement or some  encouragement \code{1} is for strong encouragement}
  }
}
\details{
 The data were collected as part of the Social Class and widening Participation in Higher Education Project based at the 
University of North London (now London Metropolitan University) and supported by the 
University's Development and Diversity Fund over the period 1998-2000.
}
\source{
 Professor Robert Gilchrist director of STORM at London Metropolitan
}
\references{
Collier T., Gilchrist R. and Phillips D. (2003), Who Plans to Go to University? Statistical Modelling of 
potential Working-Class Participants, Education Research and Evaluation, Vol 9, No 3, pp 239-263.  
}
\examples{
data(Mums)
MM<-xtabs(~mums+qual, data=Mums)
mosaicplot(MM, color=TRUE)
MM<-xtabs(~mums+ethn+gender, data=Mums)
mosaicplot(MM, color=TRUE)
}
\keyword{datasets}
