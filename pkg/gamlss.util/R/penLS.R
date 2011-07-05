#penalised Least squares
# the penLS needs checking against 
# the penReg function seems OK (it have been checked agaist gamlss+pb())
#---------------------------------------------------------------------------------------- 
penLS <- function(y, 
                  w = rep(1,length(y)), 
                 df = NULL, 
             lambda = NULL,  
              start = 10, 
              order = 1,
               plot = FALSE,
               type = c("level", "trend"), 
             method = c("ML","GAIC", "GCV"), 
                  k = 2,
                  ...)
{
library(Matrix)
#----------------------------------------------------
# lambda : the smoothing parameter
# df : the effective df's  
# if both lambda=NULL  and df=NULL then lambda is estimated using the Schall method
# if df is not NULL but lambda is NULL then df are used for smoothing
# if lambda is not NULL (whether df=NULL  or not) lambda is used for smoothing
# ---------------------------------------------------
# creates the basis for p-splines
# Paul Eilers' function
# -------------------------------------------------- 
 # a siple penalized regression     
regpen <- function(y, X, w, lambda, D)                
  {                                                   
             G <- lambda * t(D) %*% D                 
            XW <- w * X                               
           XWX <- t(XW) %*% X                      
          beta <- solve(XWX + G, t(XW) %*% y)         
            fv <- as.vector(X %*%beta)                           
             H <- solve(XWX + G, XWX)                 
           fit <- list(beta = beta, edf = sum(diag(H)))
            return(fit)                                         
  }                                                   
#---------------------------------------------------  
# a similar as obove but extra saving                 
regpenEM <- function(y, X, w, lambda, D)
  {
             G <- lambda * t(D) %*% D
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + G, t(XW) %*% y)
            fv <- as.vector(X %*%beta)
             H <- solve(XWX + G, XWX)
             V <- solve(XWX + G)
           fit <- list(beta = beta, edf = sum(diag(H)), V=V)
  return(fit)  
  }
#--------------------------------------------------
## function to find lambdas miimizing the local GAIC        
     fnGAIC <- function(lambda, k, D)
    {
    #cat("lambda, k", lambda, k, "\n")
       fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
        fv <- as.vector(X %*% fit$beta)
       sig2 <- sum(w*(y-fv)^2)/(length(y))
       NOd <- -2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)    
      GAIC <- sum(w*NOd)+k*fit$edf 
    #cat("GAIC", GAIC, "\n")
      GAIC   
    }
#--------------------------------------------------
## function to find the lambdas wich minimise the local GCV 
#      fnGCV <- function(lambda, k)
#           {
#    I.lambda.D <- (1+lambda*UDU$values)
#           edf <- sum(1/I.lambda.D)
#         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
#           GCV <- (n*y_Hy2)/(n-k*edf)
#           GCV
#           }  
#--------------------------------------------------
## local function to get df using eigen values
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }  
#--------------------------------------------------
        method <- match.arg(method)
         type <- match.arg(type) 
             y <- as.numeric(y)
             n <- length(y) 
             X <- if (type=="level") diag(n) else diag(1:n) # create the basis
             D <- if(order==0) diag(n) else diff(diag(n), diff=order) # the penalty matrix
             X <- as(X, "dgCMatrix")
             D <- as(D, "dgCMatrix")
            if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (df < 1)  1  else  df+1
              if (df < 1)  warning("the df are set to 1")    
             }
#---------------------------------------------------
# case 1: if lambda is known just fit
 if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, X, w, lambda,  D)
           fv <- as.vector(X %*% fit$beta)
         sig2 <- sum(w*(y-fv)^2)/(length(y))     
 } # case 2: if lambda is estimated ------------------------------------------- 
 else if (is.null(df) && is.null(lambda)) 
 { #   
  # cat("----------------------------","\n")
        lambda <- start# get(startLambdaName, envir=gamlss.env) ## geting the starting value
  # if ML ----------------------
  switch(method,
  "ML"={ 
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, w, lambda, D) # fit model
         gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
             fv <- as.vector(X %*% fit$beta)             # fitted values
           sig2 <- sum(w * (y - fv) ^ 2) / (n - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
     lambda.old <- lambda
         lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
     if (lambda>1.0e+7) lambda<-1.0e+7 # DS Wednesday, July 29, 2009 at 19:00
      #    cat("iter tau2 sig2",it,tau2, sig2, '\n')
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e10) break
   #   assign(startLambdaName, lambda, envir=gamlss.env)
     #cat("lambda",lambda, '\n')
         }
        sig2 <- sum(w*(y-fv)^2)/(length(y))  
       },
#  "EM"={
#        tau2 <- 0
#      for (it in 1:500) 
#         {
#             fit  <- regpenEM(y, X, w, lambda, D)
#           gamma. <- D %*% as.vector(fit$beta)
#           vgamma <- sum(diag(D%*%fit$V%*%t(D))) # this is crucial for estimating the variance of gamma Monday, March 23, 2009
#               fv <- as.vector(X %*% fit$beta)
#               sig2 <- sum(w * (y - fv) ^ 2) / (n - fit$edf) 
#             tau2.old <- tau2  
#             tau2 <- ((sum(gamma.^ 2))+vgamma)/length(gamma.) 
#       lambda.old <- lambda
#           lambda <- sig2 / tau2
#           # cat("it tau2 sig2 lambda",it, tau2, sig2, lambda, '\n')
#         if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009    
#         if (lambda>1.0e+7) lambda<-1.0e+7 # DS Wednesday, July 29, 2009 at 19:00
#       if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e20) break
#         }
#       },
  "GAIC"=
       {
        lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e20, k=k,  D=D)$par 
           fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
            fv <- as.vector(X %*% fit$beta) 
           sig2 <- sum(w * (y - fv) ^ 2) / (n - fit$edf)  
       },
 # "GCV"={
 #          QR <-suppressWarnings(qr(sqrt(w)*X))
 #          wy <- sqrt(w)*y
 #         y.y <- sum(wy^2)
 #        Rinv <- solve(qr.R(QR))
 #           S <- t(D)%*%D
 #         UDU <- eigen(t(Rinv)%*%S%*%Rinv)
 #        browser()
 #          yy <- t(UDU$vectors)%*%t(qr.Q(QR))%*%wy
 #      lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e20, k=k)$par
 #         fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
 #          fv <- as.vector(X %*% fit$beta)      
 #      }
       )
   }
  else # if df are required
  { 
 #           XW <- w * X
 #          XWX <- t(XW) %*% X
 #    edf_df <- function(lambda)
 #        {
 #         G <- lambda * t(D) %*% D
 #         H <- solve(XWX + G, XWX)
 #       edf <- sum(diag(H))
 #       #cat("edf", edf, "\n")
 #        (edf-df)
 #        }  
 #      lambda <- uniroot(edf_df, c(0,100000))$root   
 #       fit  <- regpen(y, X, w, lambda, D)
 #       fv <- as.vector(X %*% fit$beta)
        #
        QR <- suppressWarnings(qr(sqrt(w)*X))
         Rinv <- suppressWarnings(solve(qr.R(QR)))
          S   <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
       lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
                 else  uniroot(edf1_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, w, lambda, D)
            fv <- as.vector(X %*% fit$beta)
          sig2 <- sum(w*(y-fv)^2)/(length(y))       
  }             
#---------
       pfit <- list(coefficients = as.vector(fit$beta), 
                   fitted.values = fv, 
                               y = y, 
                          lambda = lambda,
                          sigma = sqrt(sig2), 
                              df = fit$edf, 
                             rss = sum(w*(y-fv)^2),
                             aic = sum(w*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)))+2*(fit$edf+1) , 
                             sbc = sum(w*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)))+log(n)*(fit$edf+1),
                        deviance = sum(w*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE))))
  class(pfit) <- "penLS"
   if (plot==TRUE)
  {
   plot(ts(y), ...)
   lines(1:length(y),fv, col="red", ...)
  }
  return(pfit)
}

#----------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
fitted.penLS<-function(object,...) 
{
object$fitted.values
}
#------------------------------------------
coef.penLS<-function(object,...) 
{
object$coefficients
}
#------------------------------------------
residuals.penLS<-function(object,...) 
{
object$y-object$fitted
}
#-------------------------------------------
AIC.penLS <- function(object, ...,k=2)
{
 val <- if (is(object, "penLS")) 
            object$deviance + (object$df+1) * k
        else stop(paste("this is not a penReg object"))
val
}
#------------------------------------------
deviance.penLS<-function(object,...) 
{
object$deviance
}
#------------------------------------------
plot.penLS <- function(x, ...)
{
    residx <- (x$y-fitted(x))/x$sigma
  x$x <- 1:length(x$y)
  x$xlabel <- "index"
 op <- par(mfrow=c(2,2), mar=par("mar")+c(0,1,0,0), col.axis="blue4", col.main="blue4", col.lab="blue4",  col="darkgreen", bg="beige" ) 
 plot(fitted(x) , residx,
         xlab = "Fitted Values",  
         ylab = "Residuals", 
         main = "Against Fitted Values",
         frame.plot = TRUE) 
    # top right  
    plot(x$x, residx, 
         ylab = "Residuals",
         xlab = x$xlabel, 
         main = paste("Against ", x$xlabel), 
         frame.plot = TRUE) #  points(par(col="blue4"))  
    plot(density(residx), 
         xlab = "Residuals", 
         ylab = "Density", 
         main = "Density Estimate",
         frame.plot = TRUE, 
         col="black", 
         lwd=0.4 ) #col="deepskyblue4", col="darkgreen", 
         rug(residx, col="red")
 browser()
    qqnorm(residx, main = "Normal Q-Q Plot",
            xlab = "Theoretical Quantiles",
            ylab = "Sample Quantiles", 
            plot.it = TRUE, 
            frame.plot = TRUE, 
            col="darkgreen"), 
            points(par(col="darkgreen")))
     lines(residx, residx, col="red" , lwd=.4, cex=.4 )
par(op)
}
#------------------------------------------
#---------------------------------------------------------------------------------------
