#---------------------------------------------------------------------------------------
# Mikis Stasinopoulos 2-11-09
#----------------------------------------------------------------------------------------
# the penalised regression function is ckecked against gamlss+pb()
# the results are identical or very similar 
# the only I am concern is EM (ML-1)
#----------------------------------------------------------------------------------------
# this is a simple smoother using P-splines
# Paul Eilers and Mikis Stasinopoulos 
penReg <- function(y, x, w = rep(1,length(y)), 
                 df = NULL, 
             lambda = NULL,  
              start = 10,
              inter = 20, 
              order = 2,
             degree = 3,
               plot = FALSE,
             method = c("ML","ML-1","GAIC", "GCV", "EM"), 
                  k = 2,
                  ...)
{
#----------------------------------------------------
# inter : is the number of equal space intervals in x
# degree: is the degree of the polynomial
# order refers to differences in the penalty matrix
# order = 0 : white noise
# order = 1 : random effect
# order = 2 : random walk
# order = 3 : random walk of order 2
# lambda : the smoothing parameter
# df : the effective df's  
# if both lambda=NULL  and df=NULL then lambda is estimated using the Schall method
# if df is not NULL but lambda is NULL then df are used for smoothing
# if lambda is not NULL (whether df=NULL  or not) lambda is used for smoothing
# ---------------------------------------------------
# creates the basis for p-splines
 bbase <- function(x, xl, xr, ndx, deg, quantiles=FALSE)
  {
 tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)
# DS xl= min, xr=max, ndx= number of points within 
# Construct B-spline basis
     dx <- (xr - xl) / ndx # DS increment 
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx) 
      P <- outer(x, knots, tpower, deg)# calculate the power in the knots
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) # 
      B <- (-1) ^ (deg + 1) * P %*% t(D) 
      B 
  }
#-------------------------------------------------  
# a siple penalized regression
regpen <- function(y, X, w, lambda, D)
  {
             G <- lambda * t(D) %*% D
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + G, t(XW) %*% y)
            fv <- X %*%beta
             H <- solve(XWX + G, XWX)
         #  edf <- sum(diag(H))
           fit <- list(beta = beta, edf = sum(diag(H)))
  return(fit)  
  }
#--------------------------------------------------
# a similar as obove but extra saving
regpenEM <- function(y, X, w, lambda, order, D)
  {
             G <- lambda * t(D) %*% D
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- solve(XWX + G, t(XW) %*% y)
            fv <- X %*%beta
             H <- solve(XWX + G, XWX)
             V <- solve(XWX + G)
           fit <- list(beta = beta, edf = sum(diag(H)), V=V)
  return(fit)  
  }
#--------------------------------------------------
## function to find lambdas miimizing the local GAIC        
     fnGAIC <- function(lambda, k)
    {
    #cat("lambda, k", lambda, k, "\n")
       fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
        fv <- X %*% fit$beta
       sig2 <- sum(w*(y-fv)^2)/(length(y))
       NOd <- -2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)    
      GAIC <- sum(w*NOd)+k*fit$edf 
    #cat("GAIC", GAIC, "\n")
      GAIC   
    }
#--------------------------------------------------
## function to find the lambdas wich minimise the local GCV 
      fnGCV <- function(lambda, k)
           {
    I.lambda.D <- (1+lambda*UDU$values)
           edf <- sum(1/I.lambda.D)
         y_Hy2 <- y.y-2*sum((yy^2)/I.lambda.D)+sum((yy^2)/((I.lambda.D)^2))
           GCV <- (n*y_Hy2)/(n-k*edf)
           GCV
           }  
#--------------------------------------------------
## local function to get edf from lambda 
#   edf_df <- function(lambda)
#         {
#             G <- lambda * t(D) %*% D
#             H <- solve(XWX + G, XWX)
#           edf <- sum(diag(H))
#          # cat("edf", edf, "\n")
#           (edf-df)
#          }
## local function to get df using eigen values
    edf1_df <- function(lambda)
           {
           edf <-  sum(1/(1+lambda*UDU$values))
           (edf-df)
           }  
#------------------------------------------------------------------
# the main function starts here
#------------------------------------------------------------------
        method <- match.arg(method)
            lx <- n <- length(x)
         inter <- if (lx<100) 10 else inter # this is to prevent singularities when length(x) is small
            xl <- min(x)
            xr <- max(x)
          xmax <- xr + 0.01 * (xr - xl)
          xmin <- xl - 0.01 * (xr - xl)   
             X <- bbase(x, xmin, xmax, inter, degree) # create the basis
             r <- ncol(X)
             D <- if(order==0) diag(r) else diff(diag(r), diff=order) # the penalty matrix
             if(!is.null(df)) # degrees of freedom
             {
             if (df>(dim(X)[2]-2)) 
              {df <- 3;  
              warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
              df <- if (df < 1)  1  else  df+2
              if (df < 1)  warning("the df are set to 1")    
             }
#---------------------------------------------------
# case 1: if lambda is known just fit
 if (is.null(df)&&!is.null(lambda)||!is.null(df)&&!is.null(lambda))
 {
          fit <- regpen(y, X, w, lambda,  D)
           fv <- X %*% fit$beta
         sig2 <- sum(w*(y-fv)^2)/(length(y))     
 } # case 2: if lambda is estimated ------------------------------------------- 
 else if (is.null(df)&&is.null(lambda)) 
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
             fv <- X %*% fit$beta             # fitted values
           sig2 <- sum(w * (y - fv) ^ 2) / (n - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
     lambda.old <- lambda
         lambda <- sig2 / tau2 # maybe only 1/tau2 will do since it gives exactly the EM results see LM-1
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009 at 14:18
     if (lambda>1.0e+7) lambda<-1.0e+7 # DS Wednesday, July 29, 2009 at 19:00
      #    cat("iter tau2 sig2",it,tau2, sig2, '\n')
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e20) break
   #   assign(startLambdaName, lambda, envir=gamlss.env)
     #cat("lambda",lambda, '\n')
         }
        sig2 <- sum(w*(y-fv)^2)/(length(y))  
       },
  "ML-1"={
       for (it in 1:50) 
         {
           fit  <- regpen(y, X, w, lambda, D) # fit model
         gamma. <- D %*% as.vector(fit$beta)  # get the gamma differences
             fv <- X %*% fit$beta             # fitted values
         #  sig2 <- 1 # sum(w * (y - fv) ^ 2) / (n - fit$edf)
           tau2 <- sum(gamma. ^ 2) / (fit$edf-order)# Monday, March 16, 2009 at 20:00 see LNP page 279
     lambda.old <- lambda
         lambda <- 1 / tau2 # 1/tau2 
     if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009
     if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e20) break
    #  assign(startLambdaName, lambda, envir=gamlss.env)
         }
       },
  "EM"={
      for (it in 1:500) 
         {
             fit  <- regpenEM(y, X, w, lambda, order, D)
           gamma. <- D %*% as.vector(fit$beta)
           vgamma <- sum(diag(D%*%fit$V%*%t(D))) # this is crucial for estimating the variance of gamma Monday, March 23, 2009
               fv <- X %*% fit$beta
             tau2 <- ((sum(gamma.^ 2))+vgamma)/length(gamma.) 
       lambda.old <- lambda
           lambda <- 1 / tau2
         if (lambda<1.0e-7) lambda<-1.0e-7 # DS Saturday, April 11, 2009    
       #    cat("iter sigma_t^2",it, tau2, "lambda",lambda, '\n')
       if (abs(lambda-lambda.old) < 1.0e-7||lambda>1.0e20) break
         }
          sig2 <- sum(w*(y-fv)^2)/(length(y))
    #cat("lambda",lambda, '\n')
     # assign(startLambdaName, lambda, envir=gamlss.env)
       },
  "GAIC"=
       {
        lambda <- nlminb(lambda, fnGAIC,  lower = 1.0e-7, upper = 1.0e20, k=k)$par 
           fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
            fv <- X %*% fit$beta 
          sig2 <- sum(w*(y-fv)^2)/(length(y))    
        #assign(startLambdaName, lambda, envir=gamlss.env)
       },
  "GCV"={
  # 
           QR <-qr(sqrt(w)*X)
           wy <- sqrt(w)*y
          y.y <- sum(wy^2)
         Rinv <- solve(qr.R(QR))
            S <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)
           yy <- t(UDU$vectors)%*%t(qr.Q(QR))%*%wy
       lambda <- nlminb(lambda, fnGCV,  lower = 1.0e-7, upper = 1.0e20, k=k)$par
          fit <- regpen(y=y, X=X, w=w, lambda=lambda, D)
           fv <- X %*% fit$beta
         sig2 <- sum(w*(y-fv)^2)/(length(y))     
      #  assign(startLambdaName, lambda, envir=gamlss.env) 
       })
  }
  else # case 3 : if df are required---------------------------------
  { 
      #method 1
      #      XW <- w * X
      #     XWX <- t(XW) %*% X
      #  lambda <- if (sign(edf_df(0))==sign(edf_df(100000))) 100000  # in case they have the some sign
      #            else  uniroot(edf_df, c(0,100000))$root
      #method 2 from Simon Wood (2006) pages 210-211, and 360 
           QR <- qr(sqrt(w)*X)
         Rinv <- solve(qr.R(QR))
          S   <- t(D)%*%D
          UDU <- eigen(t(Rinv)%*%S%*%Rinv)           
       lambda <- if (sign(edf1_df(0))==sign(edf1_df(100000))) 100000  # in case they have the some sign
                 else  uniroot(edf1_df, c(0,100000))$root
      # if (any(class(lambda)%in%"try-error")) {lambda<-100000}   
           fit <- regpen(y, X, w, lambda, D)
            fv <- X %*% fit$beta
          sig2 <- sum(w*(y-fv)^2)/(length(y))
  }#--------------------------------------------------------------------------end of case 3
#---------
         pfit <- list(coefficients = as.vector(fit$beta),
                     fitted.values = as.vector(fv), 
                                 y = y, 
                            ylabel =  substitute(y),
                                 x = x,
                            xlabel =   substitute(x),
                            lambda = lambda, 
                                df = fit$edf,
                             sigma = sqrt(sig2),
                               rss = sum(w*(y-fv)^2),
                               aic = sum(w*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)))+2*(fit$edf+1) , 
                               sbc = sum(w*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE)))+log(n)*(fit$edf+1),
                          deviance = sum(w*(-2*dNO(y, mu=fv, sigma=sqrt(sig2), log=TRUE))))
  class(pfit) <- "penReg"
  if (plot==TRUE)
  {
  plot(y~x, ...)
   lines(fv~x, ...)
  }
  return(pfit)
}
#----------------------------------------------------------------------------------------
fitted.penReg<-function(object,...) 
{
object$fitted.values
}
#------------------------------------------
coef.penReg<-function(object,...) 
{
object$coefficients
}
residuals.penReg<-function(object,...) 
{
object$y-object$fitted
}
AIC.penReg <- function(object, ...,k=2)
{
 val <- if (is(object, "penReg")) 
            object$deviance + (object$df+1) * k
        else stop(paste("this is not a penReg object"))
val
}
deviance.penReg<-function(object,...) 
{
object$deviance
}

plot.penReg <- function(x, ...)
{
  residx <- (x$y-fitted(x))/x$sigma
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
         rug(residx, col="red", points(par(col="blue4")))
 
    qqnorm(residx, main = "Normal Q-Q Plot",
            xlab = "Theoretical Quantiles",
            ylab = "Sample Quantiles", 
            plot.it = TRUE, 
            frame.plot = TRUE, 
            col="darkgreen", 
            points(par(col="darkgreen")))
     lines(residx, residx, col="red" , lwd=.4, cex=.4 )
par(op)
}
#------------------------------------------
#---------------------------------------------------------------------------------------
