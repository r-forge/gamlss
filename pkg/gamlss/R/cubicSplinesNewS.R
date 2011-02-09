# Non linear 
# This is an experimental smoothing cubic spline function 
# created  Monday, April 28, 2008 at 08:42 
# Mikis Stasinopoulos
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
scs <-function(x, df = NULL, spar = NULL, 
        control.spar = NULL, 
           all.knots = TRUE,
              nknots = NULL,
              penalty =1.4) 
{   
  scall <- deparse(sys.call()) 
     if(ncol(as.matrix(x)) > 1)
    stop(paste("The default smoother is bivariate; you gave a matrix as an argument in ", scall, "\n"))
   # len <- if(is.null(dim(x))) length(x) else dim(x)[1]
    if(!is.null(levels(x))) 
      {
        if(inherits(x, "ordered"))
            x <- codes(x)
        else stop("unordered factors cannot be used as smoothing variables")
      }
      if(is.null(control.spar)) control.spar <- list(low=-1.5,high=2) 
    else
      { control.spar <- if (is.list(control.spar) & length(control.spar)==2)  
                      {  list(low=control.spar[[1]],high=control.spar[[2]])}
                  else if (is.vector(control.spar) & length(control.spar)==2)   
                      {  list(low=control.spar[1],high=control.spar[2])}
        else stop("c.spar is not defined properly") 
      }   
     #if(!is.null(df)||df < 0) warning("the df are set to 0")
    #df <- if (is.null(df))   NULL
    #      else df
    #attr(x, "spar") <- spar
    #  attr(x, "df") <- df
       real.call <- substitute(gamlss.scs(data[[scall]], z, w, spar = spar, df = df, 
                            control.spar=control.spar, all.knots= all.knots,
                            nknots = nknots, penalty=penalty))
    attr(x, "call") <- real.call
   attr(x, "class") <- "smooth"
                  a <- is.na(x)
    if(any(a))
        attr(x, "NAs") <- seq(along = x)[a]
    x
}
#----------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.scs <-function(x, y = NULL, w = NULL, df=NULL, spar = NULL, cv = FALSE, 
    all.knots = FALSE, nknots = NULL, keep.data = TRUE, df.offset = 0, 
    penalty = 1.4, control.spar = list(low=-1.5, high=2),  xeval = NULL)
{      
             x <- signif(x, 6)
           pox <- order(x)
          freq <- table(x)
          df <- if (!is.null(df))  df+2
    #       lambda <- as.vector(attr(x,"lambda")) # lambda
    #   df <- as.vector(attr(x,"df")) # degrees of freedom
    if (is.null(df)&&is.null(spar))
      {
       fit <- smooth.spline(y=y, x=x, w=w, cv=cv,
                             all.knots=all.knots, nknots=nknots, 
                             control.spar=control.spar ,  penalty = penalty)
      }
    else if (is.null(df))
      {  fit <- smooth.spline(y=y, x=x, w=w, spar=spar, 
                              all.knots=all.knots, nknots=nknots, 
                              control.spar=control.spar, penalty = penalty)
      }
     else 
     {  fit <- smooth.spline(y=y, x=x, w=w, df=df,  
                             all.knots=all.knots, nknots=nknots, 
                             control.spar=control.spar,  penalty = penalty)
      }
        longfv  <- rep(fit$y,freq)
   longfv[pox] <- longfv # OK
          llev <- rep(fit$lev,freq) # the leverage   calculations    
     llev[pox] <- llev           
          sumw <- rep(fit$w,freq)     
     sumw[pox] <- sumw
            wt <-  (w * sum(w > 0))/sum(w)           
       longlev <- llev*wt/sumw        
          levl <- (longlev-.hat.WX(wt,x))
           var <- levl/w # MS Tuesday, June 22, 2004 at 20:58
if (is.null(xeval)) # if no prediction  
    {
       obj.out <- list(residuals=y-longfv, fitted.values = longfv, 
                   var = var,  nl.df = fit$df-2, lambda = fit$lambda, 
                   coefSmo = list(knot = fit$fit[["knot"]], #
                                    nk = fit$fit[["nk"]], 
                                   min = fit$fit[["ux[1]"]], 
                                 range = fit$fit[["r.ux"]], 
                                  coef = fit$fit[["coef"]], 
                              pen.crit = fit$pen.crit, 
                                   lev = levl,
                   lambda1 = fit$lambda))
class(obj.out) <- "smooth.spline"
    obj.out
    }
else 
   { # if  prediction  
                pred <- predict(fit,x = xeval)
                pred$y  
     }    
}







 
