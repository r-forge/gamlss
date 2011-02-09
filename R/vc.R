# varying coefficient
#----------------------------------------------------------------------------------------
# MS 
vc <-function(r, x, df = 3, spar = NULL, c.spar=NULL) 
{
    scall <- deparse(sys.call())
    if(ncol(as.matrix(x)) > 1)
    stop(paste("The default smoother is bivariate; you gave a matrix as an argument in ",
      scall, "\n"))
    if(ncol(as.matrix(r)) > 1)
    stop(paste("The default smoother is bivariate; you gave a matrix as an argument in ",
      scall, "\n"))  
    if(!is.null(levels(r))) 
      {
        if(inherits(r, "ordered"))
       r <- codes(r)
        else stop("unordered factors cannot be used as smoothing variables")
      }
      
     if(!is.null(levels(x))) 
      {
       x <-  if(inherits(x, "ordered")) codes(x)
             else if  (nlevels(x) == 2 ) ifelse(x==levels(x)[2], 1,0)
             else  stop("only two level factors with 0 and 1's are allowed as covariates in vc()")
      } 
    if(is.null(c.spar)) c.spar <- list(low=-1.5,high=2) 
    else                c.spar <- list(low=c.spar[1],high=c.spar[2] )
    if(df < 0) warning("the df are set to 0")
    df <- if (df < 0)   0
          else          df
    #-----from the mean
    xvar <- if (range(x)[1]==0&range(x)[2]==1) x else x-mean(x)
    rvar <- r#-mean(r)
      xr <- model.matrix(~rvar*xvar-1,contrast="")
    attr(xr, "spar") <- spar
      attr(xr, "df") <- df
    attr(xr, "rvar") <- rvar
    attr(xr, "xvar") <- xvar
   #   attr(xr, "lambda") <- NULL
    real.call <- substitute(gamlss.vc(data[[scall]], z, w, spar = spar, df = df+2,# df+1
                            control.spar=c.spar))
    attr(xr, "call") <- real.call
   attr(xr, "class") <- "smooth"
                   a <- is.na(x)
    if(any(a))
     attr(xr, "NAs") <- seq(along = r)[a]
    xr
}

#------------------------------------------------------------------------------------
# this is a modified version of the smooth.spline in R stats
#  MS Thursday, Monday, November 10, 2003 at 10:50
gamlss.vc <- function (x, y = NULL, w = NULL, df = 5, spar = NULL, 
                       cv = FALSE,  all.knots = TRUE, df.offset = 0, 
                       penalty = 1, control.spar = list(low=-1.5,high=2),
                       xeval = NULL) 
{
    sknotl <- function(x) 
         { #function 1
         n.kn <- function(n) 
            { # function 2
             if (n < 50) 
                 n
             else trunc(
              {
                a1 <- log(50, 2)
                a2 <- log(100, 2)
                a3 <- log(140, 2)
                a4 <- log(200, 2)
                if (n < 200) 2^(a1 + (a2 - a1) * (n - 50)/150) else if (n < 
                  800) 2^(a2 + (a3 - a2) * (n - 200)/600) else if (n < 
                  3200) 2^(a3 + (a4 - a3) * (n - 800)/2400) else 200 + 
                  (n - 3200)^0.2
              })
            } # end of function 2
           nk <- n.kn(n <- length(x))
           c(rep(x[1], 3), x[seq(1, n, len = nk)], rep(x[n], 3))
          } # end of function 1
  
  
   #--- control for spar---
    contr.sp <- list(low = -1.5, high = 1.5, tol = 1e-04, eps = 2e-08, 
                     maxit = 500, trace = getOption("verbose"))
    contr.sp[(namc <- names(control.spar))] <- control.spar
    if (!all(sapply(contr.sp[1:4], is.double)) || contr.sp$tol < 
        0 || contr.sp$eps <= 0 || contr.sp$maxit <= 0) 
        stop("invalid `control.spar'")
# 
   #---get y x r w---  
     #browser()
       if (is.null(xeval))
        {
       rvar <- as.vector(attr(x,"rvar"))  
       xvar <- as.vector(attr(x,"xvar"))  
        }
       else
        {
       rvar <-  x[,"rvar"]
       xvar <-  x[,"xvar"]    
        }
         ry <- xy.coords(rvar, y)
       yvar <- ry$y #  -eval.parent(expression(s[,j])) 
         wt <- w #ms Tuesday, June 22, 2004 at 23:16  
          y <- (yvar/ifelse(xvar==0,0.0001,xvar)) # we should check fo zeros 
          r <- ry$x
          n <- length(r)
          w <- if (is.null(w)) 
                  xvar^2
               else 
               {
               if (n != length(w)) 
               stop("lengths of x and w must match")
               if (any(w < 0)) 
               stop("all weights should be non-negative")
               if (all(w == 0)) 
               stop("some weights should be positive")
               (w * xvar^2*sum(w > 0))/sum(w)
               }
      #   print(cbind(y,w))      
         W <- rep(1, n) # this for counting repeated x's
         r <- signif(r, 6)
       por <- order(r)
        ur <- unique(sort(r)) # unique values of x
        or <- match(r, ur)
        tmp <- matrix(unlist(tapply(
                              seq(along = y), or, 
                                   function(i, y, w, W) 
                                  c(
                                   sum(w[i]), 
                                   sum(w[i] * y[i]), 
                                   sum(w[i] * y[i]^2), 
                                   sum(W[i])
                                   ), 
                                   y = y, 
                                   w = w, 
                                   W = W
                                 )
                           ), 
                            ncol = 4, 
                            byrow = TRUE)
       wbar <- tmp[, 1] # weights for fitting 
       ybar <- tmp[, 2]/ifelse(wbar > 0, wbar, 1) 
       yssw <- sum(tmp[, 3] - wbar * ybar^2)
       Wbar <- tmp[,4]
         nr <- length(ur) # no unique x values
       if (nr <= 3) 
        stop("need at least for unique `r' values")
       if (cv && nr < n) 
        warning("crossvalidation with non-unique `x' seems doubtful")
       r.ur <- ur[nr] - ur[1]
       rbar <- (ur - ur[1])/r.ur # 0<=xbar<=1       
       if (all.knots) 
       {
        knot <- c(rep(rbar[1], 3), rbar, rep(rbar[nr], 3))
          nk <- nr + 2
       }
       else 
       {
        knot <- sknotl(rbar)
          nk <- length(knot) - 4
       }
      ispar <- if (is.null(spar) || missing(spar)) 
                  {
                   if (contr.sp$trace) 
                      -1
                  else 0
                  }
               else 1
      spar <- if (ispar == 1) 
                 as.double(spar)
             else double(1)
     icrit <- if (cv) 
                   2
             else 1
    dofoff <- df.offset
    if (!missing(df)) 
       {
        if (df > 1 && df <= nr) 
           {
             icrit <- 3
            dofoff <- df
           }
        else warning(paste("you must supply 1 < df <= n,  n = #{unique x} =", 
            nr))
       }
    iparms <- as.integer(c(icrit, ispar, contr.sp$maxit))
    names(iparms) <- c("icrit", "ispar", "iter")
    fit <- .Fortran("qsbart", as.double(penalty), as.double(dofoff), 
                    x = as.double(rbar), y = as.double(ybar), w = as.double(wbar), 
                    ssw = as.double(yssw), as.integer(nr), as.double(knot), 
                    as.integer(nk), coef = double(nk), ty = double(nr), 
                    lev = double(nr), crit = double(1), iparms = iparms, 
                    spar = spar, parms = unlist(contr.sp[1:4]), 
                    isetup = as.integer(0), scrtch = double((17 + nk) * nk), 
                    ld4 = as.integer(4), ldnk = as.integer(1), ier = integer(1), 
                    DUP = FALSE, PACKAGE = "stats")[c("coef", "ty", "lev", 
                    "spar", "parms", "crit", "iparms", "ier")]
    lev <- fit$lev
    df <- sum(lev)
    if (is.na(df))   
       stop("NA lev[]; probably smoothing parameter `spar' way too large!")
    if (abs(df-dofoff)>1) 
    warning("The output df ", df,"  are different for the input ", dofoff,
             "\n","change the control.spar" )
    if (fit$ier > 0) 
       {
         sml <- fit$spar < 0.5
        wtxt <- paste("smoothing parameter value too", 
                      if (sml) 
                        "small"
                      else "large")
        if (sml) 
         {
            stop(wtxt)
         }
        else 
         {
            fit$ty <- rep(mean(y), nr)
            df <- 1
            warning(paste(wtxt, "setting df = 1  __use with care!__", 
                sep = "\n"))
         }
       }
    cv.crit <- if (cv) 
                 {
                  ww <- wbar
                  ww[!(ww > 0)] <- 1
                  weighted.mean(((y - fit$ty[or])/(1 - (lev[or] * w)/ww[or]))^2, 
                  w)
                 }
             else weighted.mean((y - fit$ty[or])^2, w)/(1 - (df.offset + 
                     penalty * df)/n)^2
      pen.crit <- sum(wbar * (ybar - fit$ty) * ybar)
        longy  <- rep(fit$ty, Wbar)
    longy[por] <- longy 
          llev <- rep(lev,Wbar) # the leverage   calculations    
     llev[por] <- llev           
          sumw <- rep(wbar,Wbar)     
     sumw[por] <- sumw           
       longlev <- llev*w/sumw        
          levl <- (longlev-.hat.WX(w,x))
           var <- levl/wt # MS Tuesday, June 22, 2004 at 20:58              
          fval <-  longy*xvar 
if (is.null(xeval))
    {   
     object <- list(residuals=yvar-fval, fitted.values = fval, var = var, 
                   nl.df = df-2, lambda = fit$spar, #df-4
                   coefSmo = list(knot = knot, nk = nk, min = ur[1], range = r.ur, 
                   coef = fit$coef, pen.crit = pen.crit, lev = levl, 
                   lambda1 = unname(fit$parms["low"]) ))
    class(object) <- "smooth.spline"
    object
    }
else {
  #------------------------
           fit.object <- list(knot = knot, nk = nk, min = ur[1], range = r.ur, 
                    coef = fit$coef)
  class(fit.object) <- "smooth.spline.fit"
             object <- list(x = ur, y = fit$ty, w = wbar, yin = ybar, 
                           lev = lev, cv.crit = cv.crit, pen.crit = pen.crit, 
                           crit = fit$crit, df = df, spar = fit$spar, 
                           lambda = unname(fit$parms["low"]), 
                           iparms = fit$iparms, fit = fit.object, 
                           call = match.call())
       class(object) <- "smooth.spline"   
                pred <- predict(object,x=xeval[,"rvar"])
                pred$y*xeval[,"xvar"]  
  #------------------------
     }    
}
