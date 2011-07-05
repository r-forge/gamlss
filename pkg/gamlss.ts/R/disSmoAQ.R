# AUTHORS 
# Mikis Stasinopoulos, Bob Rigby, Paul Eilers, Majid  Djennad
# this file contains  functions to fit discrete Smoother in time series data 
# I has following  main functions 
# i)   disSmoA for Alternate fitting when the smoothing parameters is estimated
# ii)  disSmoQ for Q function fiting when the smoothing parameters is estimated
# iii) disSmo which combines the above two functions 
# both functions disSmoA and disSmoQ  should give identical results
# Both functions are ussing  the new spam 0.24 package (Not in CRAN) 
# NOTE AT THE MOMENT CAN NOT RUN WITHOUT spam 0.24
# 
# The method in disSmoA "A"lternating between sigma and sigma_b 
# note that this method requires and esimate of the degres of freedom
# The  degrees of freedom are the trace of the smoothing matrix
# (W+lambdaG)^(-1) W 
# which is difficult to calculate for large N 
# therefore an approximation is used here.
# this approximation is very accurate (comapred with the Q method)
# for unweighted obsevations but can result to discprepancies if the data are weighted
# both methods create an disSmo object which has the following methods
# print(), fitted(), resid(), AIC(), deviance(), plot(), coef()
# Not that wp() and dtop() can also be used
#  TO DO
# predict()-forecast()
# missing values ??
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# METHOD ALTER : alternate ML estimation
disSmoA <- function(y, 
             weights = rep(1,length(y)),  
              lambda = NULL,
               order = 1,
               start = 10,
                plot = FALSE)                    
 {
 # local functions
 #------------------------------- 
      regpen <- function(y, weights, lambda) 
       {
               fv <- solve(W+lambda*G,  weights*y)
               if (order==0)
                 { # if order zero we can get them easy  
                    edf <- sum(weights/(weights+lambda))    
                 }
                else
                 { # here is more difficult to get dfs
                    if (N<=500) edf <- sum(diag(solve(W+lambda*G,  W)))
                    else  # otherwise get a sample of N1 
                     {    # and try to estimate it see   P. Eilers (2003) method (A Perfect Smoother)
                                     II <- seq(1,N, ceiling( N/(log(N)*50))) # get a sample 
                                     N1 <- length(II)                        # of length N1
                                     z2 <- rep(0,N)                          # create two dummy variable
                                     z1 <- rep(0,N1)
                    z2[as.integer(N/2)] <- 1                                 # put 1 in the middle
                   z1[as.integer(N1/2)] <- 1                                  
                                lambda1 <- lambda #((N1/N)^(2*order))*(lambda) # 
                                     #yy <- y[II]                             # get sample data
                                     ww <- weights[II] 
                                     WW <- diag.spam(x=ww)
                                     EE <- diag.spam(N1)     
                                     DD <- diff(EE, diff = order)   
                                     GG <- t(DD)%*%DD
                                   edf1 <- sum(H1<-diag(solve(WW+lambda1*GG, WW))) # get edf from the sample
                                    edf <- edf1 +(N-N1)* mean(H1[c(as.integer((N1/2)-100):as.integer((N1/2)+100))])
                  
                     } 
                  }      
             fit <- list(fv=fv, edf=edf)
            fit 
        }
 #-------------------------------	          
   require(spam); require(gamlss)
        scall <- deparse(sys.call())
            N <-  length(y)
            W <- diag.spam(x=weights)
            E <- diag.spam(N)     
            D <- diff(E, diff = order)   
            #D <- if (order==0) E
            #     else  -1*E[-N, ] + E[-1, ]
                 #as.spam(diff(E, differences=order))
            G <- t(D)%*%D 
            sig2 <- tau2 <- 0
            
         if (is.null(lambda)) 
           {
                    lambda <- start
                for (it in 1:200) 
                {
                          fit <- regpen(y, weights, lambda)
                       gamma. <- diff(fit$fv, differences=order)
                           fv <- fit$fv          
                         sig2 <- sum(weights * (y - fv) ^ 2) / (N - fit$edf)
                         tau2 <- sum(gamma.^2)/(fit$edf - order)
                   lambda.old <- lambda
                       lambda <- sig2/tau2
                  if (lambda < 1e-07) lambda <- 1e-07
                  if (lambda > 1e+07) lambda <- 1e+07
                  if (abs(lambda - lambda.old) < 1e-07 || lambda > 1e+10) break
                }
              }            
           fit <- regpen(y, weights, lambda) 
          sig2 <- sum(weights * (y - fit$fv) ^ 2) / (N - fit$edf)
          tau2 <- sum(diff(fit$fv, differences=order)^2)/(fit$edf - order)          
# saving staff 
           fit <- list(fitted = fit$fv, 
                           df = fit$edf, 
                       lambda = lambda, 
                        order = order, 
                         sig2 = sig2,
                        sigma = sqrt(sig2),
                         tau2 = tau2,
                            y = y,
                      weights = weights,
                            N = N,
                         call = scall,
                          rss = sum(weights*(y-fit$fv)^2),
                          aic = sum(weights*(-2*dNO(y, mu=fit$fv, sigma=sqrt(sig2), log=TRUE)))+2*(fit$edf+1) , 
                          sbc = sum(weights*(-2*dNO(y, mu=fit$fv, sigma=sqrt(sig2), log=TRUE)))+log(N)*(fit$edf+1),
                    deviance = sum(weights*(-2*dNO(y, mu=fit$fv, sigma=sqrt(sig2), log=TRUE))))
          class(fit) <- "disSmo"
           if (plot) 
           {
            plot(y, col = gray(0.7), type = 'l') 
           lines(fitted(fit), col = 'red', lwd = 2)
            }
           return(fit)     
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# disSmoQ
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------      

#-------------------------------------------------------------------------------------
# disctere smoothing using the Q function
# that is the evaluation of the smoothing parameters is done
# by maximising the Q-function of Pawitan 2002
# which is a profile likelihood for the sigma and tau parameters 
# it uses the spam package version 0.24 which is NOT in the CRAN yet
# this was based on the function disSmo3() of Mikis Stasinopoulos
#-------------------------------------------------------------------------------------
disSmoQ <- function(y, 
                weights = rep(1,length(y)), # for weighted observations 
                  order = 1,
                  start = 10,
                   plot = FALSE)                       
 {          
  require(spam); require(gamlss)
        scall <- deparse(sys.call())
      lambda  <- start 
            N <- length(y)
            W <- diag.spam(x=weights)
            E <- diag.spam(N)
            D <- diff(E, diff = order)   
            G <- t(D)%*%D  
            #       BD <- bdiag(W+lambda*G) 
           fv <- solve(W+lambda*G,  weights*y)        
       params <- c(sige <- sum(weights*(y-fv)^2)/N, sigb<-sum((D%*%fv)^2)/N)        
    Qf<- function(par)
     {
       lambda <- par[1]/par[2]
           fv <-solve(W+lambda*G,  weights*y)   
            b <-  diff(fv, differences=order)
           D3 <-  determinant((1/par[1])*W+(1/par[2])*G)$modulus
            f <- -(N/2)*log(2*pi*par[1]) -sum(weights*(y-fv)^2)/(2*par[1])-((N-order)/2)*log(2*pi*par[2])-sum(b^2)/(2*par[2]) -.5*D3+(N/2)*log(2*pi)
           -f
     }
    # out <- optim (params, Qf, method = c("L-BFGS-B"),lower = c(0,0), upper = c(1000000, 1000000), )
    out<-nlminb(start = params, objective = Qf, lower = c(0.00001,0.000001),  upper = c(Inf, Inf))
      sige <- out$par[1]
      sigb <- out$par[2]
    lambda <- sige/sigb
        fv <- solve(W+lambda*G,  weights*y)   
         b <-  diff(fv, differences=order)
       tr1 <- order + sum(b^2)/(sigb)
       tr2 <- N-(sum(weights*(y-fv)^2))/(sige)
    if (plot) {plot(y, col = gray(0.7), type = 'l'); lines(fv, col = 'red', lwd = 2)}
    # get the output 
           fit <- list(fitted = fv, 
                           df = tr1,
                           df2= tr2,
                       lambda = lambda, 
                        order = order, 
                         sig2 =  sige,
                        sigma = sqrt(sige),
                         tau2 = sigb,
                            y = y,
                      weights = weights,
                            N = N,
                         call = scall,
                          rss = sum(weights*(y-fv)^2),
                          aic = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE)))+2*(tr1) , 
                          sbc = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE)))+log(N)*(tr1),
                     deviance = sum(weights*(-2*dNO(y, mu=fv, sigma=sqrt(sige), log=TRUE))))
          class(fit) <- "disSmo"
  fit
  }



#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# disSmo 
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
disSmo<- function(y, 
              method = c("Alter", "Qfun"),
                   ...)
{
  method <- match.arg(method)
   scall <- deparse(sys.call())
     fit <- switch(method, "Alter"= disSmoA(y,...),
                    "Qfun" = disSmoQ(y,...)) 
fit$call <- scall
 fit
}

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# disSmo methods
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
fitted.disSmo<-function(object,...) 
{
object$fitted
}
#------------------------------------------
coef.disSmo<-function(object,...) 
{
list(sigma=object$sigma, tau=sqrt(object$tau2))
}
#------------------------------------------
residuals.disSmo<-function(object,...) 
{
qNO(pNO(object$y, mu=object$fitted, sigma=sqrt(object$sig2)))
}
#-------------------------------------------
AIC.disSmo <- function(object, ...,k=2)
{
 val <- if (is(object, "disSmo")) 
            object$deviance + (object$df+1) * k
        else stop(paste("this is not a disSmo object"))
val
}
#------------------------------------------
deviance.disSmo<-function(object,...) 
{
object$deviance
}
#------------------------------------------
plot.disSmo <- function (x, xvar=NULL, parameters=NULL, ts=FALSE, summaries=TRUE, ...) 
{
    residx <- resid(x) # get the residuals 
         w <- x$weights
       x$x <- if (is(x$y , "ts"))  as.numeric(time(x$y))
              else                 1:length(x$y)
    xlabel <- if(!missing(xvar)) deparse(substitute(xvar)) else deparse(substitute(index))
## plotting parameters
    if(is.null(parameters))
          op <- par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col.main="blue4", col.lab="blue4",  col="darkgreen", bg="beige" )
    else  op <- parameters
## now the two top  figures 
## if time series plot acf and pacf  
    if(identical(ts, TRUE))
     {  # get the acf and pacf
     acf.new<-acf(residx,plot=FALSE)
     plot(acf.new,xlim=c(2,length(acf.new$acf)),ylim=range(acf.new$acf[-1]))   # ms Tuesday, August 19, 2003 at 11:04
     pacf(residx)
     }
     else 
     {# otherwise 
     ## I am assuming that is x$noObs!=x$N then we have weights (with frequencies)
     if (length(residx)==x$N)
        {
         fittedvalues <- if(is.null(fitted(x))) fitted(x,"sigma") else fitted(x) # MS Wednesday, September 10, 2003 at 21:20
         ## whether index or x-variable
         if(is.null(xvar))     xvar <- seq(1,length(residx),1) # MS
        }
        else
        { # if weights
         fittedvalues <- rep( if(is.null(fitted(x))) fitted(x,"sigma") else fitted(x), w)
          xvar <- if(is.null(xvar))  seq(1,length(residx),1) else rep(xvar,w)
        } 
    # top left
    plot(fittedvalues , residx,
         xlab = "Fitted Values",  
         ylab = "Quantile Residuals", 
         main = "Against Fitted Values",
         frame.plot = TRUE) 
    # top right  
    plot(xvar, residx, 
         ylab = "Quantile Residuals",
         xlab = xlabel, 
         main = paste("Against ", xlabel), 
         frame.plot = TRUE) #  points(par(col="blue4"))
     }    
    plot(density(residx), 
         xlab = "Quantile. Residuals", 
         ylab = "Density", 
         main = "Density Estimate",
         frame.plot = TRUE, 
         col="black", 
         lwd=0.4 ) #col="deepskyblue4", col="darkgreen", 
         rug(residx, col="red")
 
    qqnorm(residx, main = "Normal Q-Q Plot",
            xlab = "Theoretical Quantiles",
            ylab = "Sample Quantiles", 
            plot.it = TRUE, 
            frame.plot = TRUE, 
            col="darkgreen")
       lines(as.numeric(residx), as.numeric(residx), col="red" , lwd=.4, cex=.4 )
 
     if ( identical(summaries, TRUE))
               { 
                     qq <- as.data.frame(qqnorm(residx, plot = FALSE))
               Filliben <- cor(qq$y,qq$x)
                    # mr <- as.matrix(residx)
                    m.1 <- mean(residx)
                    m.2 <- var(residx) # cov.wt(mr,w)$cov
                  n.obs <- sum(w) 
                    m.3 <- sum((residx-m.1)**3)/n.obs 
                    m.4 <- sum((residx-m.1)**4)/n.obs 
                    b.1 <- m.3^2/m.2^3
                sqrtb.1 <- sign(m.3)*sqrt(abs(b.1))
                    b.2 <- m.4/m.2^2 
                     cat("*******************************************************************")
                     cat("\n")
                     if (identical(x$type,"Continuous")) 
                         {cat("\t","     Summary of the Quantile Residuals")}
                     else{cat("\t","Summary of the Randomised Quantile Residuals")}    
                     cat("\n")
                     cat("                           mean   = ", m.1, "\n")
                     cat("                       variance   = ", m.2, "\n")
                     cat("               coef. of skewness  = ", sqrtb.1, "\n")
                     cat("               coef. of kurtosis  = ", b.2, "\n")
                     cat("Filliben correlation coefficient  = ", Filliben, "\n")
                     cat("*******************************************************************")
                     cat("\n")

               }
    par(op)
}
#----------------------------------------------------------------------------------------------

print.disSmo  <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{   
    cat("\nDiscrete Smoothing fit")
    cat("Fitting method:", deparse(x$method), "\n")
    cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
    cat("Coefficients")
    co <- coef(x)
    cat("\n  ", names(co), "\n")
    cc <-simplify2array(co, higher=TRUE)
    cat(cc, " \n")
    cat("\n Degrees of Freedom for the fit:", x$df, "Residual Deg. of Freedom  ", 
        x$N-x$df, "\n")
    cat("Global Deviance:    ", format(signif(x$deviance)), 
        "\n            AIC:    ", format(signif(x$aic)), "\n            SBC:    ", 
        format(signif(x$sbc)), "\n")
    invisible(x)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#-------------------------------------------------------------------
# function to generate y values from a random walk senario
simRW <- function(N=1000, mu=100, sige=10, sigb=2)
{
  wn1 <- ts(rnorm(N, 0, sigb))
    y <- ts(mu+cumsum(wn1))+rnorm(N, m=0,s=sige)
    y
}
#--------------------------------------------------------------------
