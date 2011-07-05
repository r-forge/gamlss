#----------------------------------------------------------------------------------------
# for fixed knots it fits a piecewise beta spline
# created by MS Tuesday, July 7, 2009 
#----------------------------------------------------------------------------------------
fitFixBP <- function(x,y,w = NULL, knots=NULL, degree=3, fixed = NULL, ...)
{
#local functions ---
# A least square fit---
ls.wr <- function(y, X, w)
  {
            XW <- w * X
           XWX <- t(XW) %*% X
          beta <- as.vector(solve(XWX , crossprod(XW, y)))
            fv <- as.vector(X %*%beta)
           rss <- sum((w*(y-fv))^2)
           fit <- list(coefficients = beta, fv=fv, rss=rss, df=length(beta) ) #edf = sum(diag(H)))
  return(fit)
  }  
#--------------------------
# main function starts here
             lx <- length(x)
             xl <- min(x)
             xr <- max(x)
           xmax <- xr + 0.01 * (xr - xl)
           xmin <- xl - 0.01 * (xr - xl)
           if (is.null(knots)) warning("No knots (break points) are declared")   
             kn <-  sort(c(xmin, xmax, fixed, knots)) 
              # splineDesign(knots, x, ord = 4, derivs, outer.ok = FALSE)
             X2 <- splineDesign(knots=kn, x, ord= degree + 1, derivs=0 * x, outer.ok=TRUE)#
                   #bbase(x=x, xl=xmin, xr=xmax, knots=knots, deg=degree) # the problem here is the difference between break points
             form<- switch(as.character(degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
                           "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
             X1 <- model.matrix(form) # model.matrix(~poly(x,degree)) 
             X <-cbind(X1,X2) 
             w <- if (is.null(w)) rep(1, lx)
                   else rep(w, length = lx)
           #fit <- ls.wr(y,X,w)
            # or  
           fit <- lm.wfit(X,y,w) 
          out  <- list(fitted.values=fitted(fit), residuals=resid(fit), df=fit$rank, rss=sum(resid(fit)^2), knots=kn,
                        coef=fit$coefficients,  
                        fixed=fixed, breakPoints=knots,  degree=degree, 
                        y=y, x=x, w=w)
          #out <- list(fitted.values=fit$fv,       residuals=y-fit$fv,  df=fit$df,   rss=fit$rss,         knots=kn,
                       #, X=X 
   class(out) <- "FixBreakPointsReg"
         out
}
#------------------------------------------
fitted.FixBreakPointsReg<-function(object,...) 
{
object$fitted.values
}
#------------------------------------------
residuals.FixBreakPointsReg<-function(object,...) 
{
object$residuals
}
#------------------------------------------
coef.FixBreakPointsReg<-function(object,...) 
{
object$coef
}
#------------------------------------------
knots.FixBreakPointsReg<-function(Fn,...)
{
Fn$breakPoints
}
#------------------------------------------
predict.FixBreakPointsReg<-function(object, newdata=NULL, old.x.range=TRUE,...)
{
if (is.null(newdata))  #
    {
    predictor<- fitted(object)
    return(predictor)
    }
  lox <- length(object$x)
    x <- c(object$x, newdata)    
   lx <- length(x)
# here we have two choices
# the first is to create new end points including the newdata
# the second is to only include the old data end points
# If we take the fist choice old.x.range=TRUE
# the prediction could be possible better outside the x range
# but would not coincide with the original predictions i.e. fitted(model)
# here the default includes the old data point range with the consequence that predictions 
# outside the original x-range will linear, quadratic, cubic  etc, depending on the 
# value of the degree of the polynomial fitted
if (old.x.range)
 {
   X2 <- splineDesign(knots=object$knots, x, ord= object$degree + 1, derivs=0 * x, outer.ok=TRUE)#
 }
 else
 {
   xl <- min(x)
   xr <- max(x)
 xmax <- xr + 0.01 * (xr - xl)
 xmin <- xl - 0.01 * (xr - xl)
   kn <- sort(c(xmin, xmax, object$fixed, knots(object))) 
   X2 <- splineDesign(knots=kn, x, ord= object$degree + 1, derivs=0 * x, outer.ok=TRUE)#
}   
 form <- switch(as.character(object$degree), "0"=~1, "1"=~x, "2"=~x+I(x^2), "3"=~x+I(x^2)+I(x^3), 
                           "4"=~x+I(x^2)+I(x^3)+I(x^4), "5"=~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) 
   X1 <- model.matrix(form) # model.matrix(~poly(x,degree))
    X <- cbind(X1,X2) 
   fv <- X%*%coef(object) 
   if (any(abs(fv[1:lox]-fitted(object))>1e-005)) 
warning(paste("There is a discrepancy  between the original prediction and the re-fit"," \n used to achieve 'safe' predictions \n ", sep = "" ))  
 pred <-  fv[(lox+1):lx]
 pred
}

#----------------------------------------------------------------------------------------
fitFreeKnots <- function(x,y,w = NULL, knots=NULL, degree=3, fixed = NULL, trace = 0,...)
{
#------------------------------------------------------------------
    penalty.opt <- function(kn, x, y, w, k, fixed = NULL, degree, ...) 
        {
       
       # kn <- sort(c(kn, fixed))
        if (length(u <- unique(knots)) < length(knots)) 
            stop(sprintf("%d coincident knot(s) detected", length(kn) - length(u)))
        #       fitFixBP(x, y, knots = knots, degree=degree, fixed = fixed)
        #cat("knots", kn, "\n")
        sp <-  fitFixBP(x=x, y, w, knots = kn, degree=degree, ...)
        sp$rss
       }
#-------------------------------------------------------------------
       lx <- length(x)
     xmin <- min(x)
     xmax <- max(x)
    #    w <- if (is.null(w))  rep(1, lx)
    #        else rep(w, length = lx)
      eps <- 5e-04
    shift <- .Machine$double.eps^0.25
        g <- length(knots)
    # initial fit
      sp0 <- fitFixBP(x=x, y=y, w=w, knots = knots, degree=degree, fixed = fixed)
   sigma0 <- sp0$rss
        lambda <- if (length(knots) > 1) 
        {
            optim(knots, penalty.opt, #if (is.null(fixed)) penalty.gr, 
                  x = x, y = y, w = w, degree = degree, 
                lower = xmin + shift, upper = xmax - shift, method = "L-BFGS-B", 
                eps = eps, fixed = fixed )$par #, control = list(trace = trace - 1))$par
        }
        else
        {
            optimize(penalty.opt, c(xmin, xmax), x = x, y = y, w=w,  degree=degree, fixed = fixed)$minimum
        }
        
        fit <- fitFixBP(x=x, y=y, w=w, knots = lambda, degree=degree, fixed = fixed) 
        out <- list(fitted.values=fitted(fit), residuals=resid(fit), df=fit$df+length(lambda), rss=fit$rss, knots= fit$knots,
                    fixed=fit$fixed, breakPoints=lambda,  coef=coef(fit), degree=degree, y=y, x=x, w=w )
 class(out) <- c("FreeBreakPointsReg", "FixBreakPointsReg")
         out
 }
#----------------------------------------------------------------------------------------
fitted.FreeBreakPointsReg<-function(object,...) 
{
object$fitted.values
}
#------------------------------------------
residuals.FreeBreakPointsReg<-function(object,...) 
{
object$residuals
}
#------------------------------------------
coef.FreeBreakPointsReg<-function(object,...) 
{
object$coef
}
#----------------------------------------------------------------------------------------
