##---------------------------------------------------------------------------------------
##  CUBIC SPLINES 
##  the cs() function is based on Hastie's gam cubic spline function s()
##  the gamlss.cs() is based on Brian Ripley function smooth.spline()
##  author Mikis Stasinopoulos
##  last change MS Tuesday, April 29, 2008  
cs <-function(x, df = 3, spar = NULL, c.spar=NULL) 
{
    scall <- deparse(sys.call())
    if(ncol(as.matrix(x)) > 1)
    stop(paste("The default smoother is bivariate; you gave a matrix as an argument in ",
     scall, "\n"))
    if(!is.null(levels(x))) 
      {
        if(inherits(x, "ordered"))
            x <- as.numeric(x) # MS Saturday, July 18, 2009 
        else stop("unordered factors cannot be used as smoothing variables")
      }
    if(is.null(c.spar)) c.spar <- list(low=-1.5,high=2) 
    else
      { c.spar <- if (is.list(c.spar) & length(c.spar)==2)  
                      {  list(low=c.spar[[1]],high=c.spar[[2]])}
                  else if (is.vector(c.spar) & length(c.spar)==2)   
                      {  list(low=c.spar[1],high=c.spar[2])}
        else stop("c.spar is not defined properly") 
      }   
    if(df < 0) warning("the df are set to 0")
    df <- if (df < 0)   0
          else          df
    attr(x, "spar") <- spar
      attr(x, "df") <- df
   #   attr(x, "lambda") <- NULL
    real.call <- substitute(gamlss.cs(data[[scall]], z, w, spar = spar, df = df+2, 
                            control.spar=c.spar))
    attr(x, "call") <- real.call
   attr(x, "class") <- "smooth"
                  a <- is.na(x)
    if(any(a))
        attr(x, "NAs") <- seq(along = x)[a]
    x
}
##---------------------------------------------------------------------------------------
## this is a modified version of the smooth.spline() in R stats
##  MS Thursday, Monday, November 10, 2003 at 10:50
gamlss.cs <- function (x, y = NULL, w = NULL, df = 5, spar = NULL, 
                       cv = FALSE,  nknots = NULL, all.knots = TRUE, df.offset = 0, 
                       penalty = 1, control.spar = list(low=-1.5,high=2),
                       xeval = NULL) 
{
 sknotl <- function(x, nk = NULL) {
        n.kn <- function(n) {
            if (n < 50) 
                n
            else trunc({
                a1 <- log(50, 2)
                a2 <- log(100, 2)
                a3 <- log(140, 2)
                a4 <- log(200, 2)
                if (n < 200) 2^(a1 + (a2 - a1) * (n - 50)/150) else if (n < 
                  800) 2^(a2 + (a3 - a2) * (n - 200)/600) else if (n < 
                  3200) 2^(a3 + (a4 - a3) * (n - 800)/2400) else 200 + 
                  (n - 3200)^0.2
            })
        }
        n <- length(x)
        if (is.null(nk)) 
            nk <- n.kn(n)
        else if (!is.numeric(nk)) 
            stop("'nknots' must be numeric <= n")
        else if (nk > n) 
            stop("cannot use more inner knots than unique 'x' values")
        c(rep(x[1], 3), x[seq.int(1, n, length.out = nk)], rep(x[n], 
            3))
    }
   #--- control for spar---
    contr.sp <- list(low = -1.5, high = 1.5, tol = 1e-04, eps = 2e-08, 
                     maxit = 500, trace = getOption("verbose"))
    contr.sp[(namc <- names(control.spar))] <- control.spar
    if (!all(sapply(contr.sp[1:4], is.double)) || contr.sp$tol < 
        0 || contr.sp$eps <= 0 || contr.sp$maxit <= 0) 
        stop("invalid `control.spar'")
   #---get y x w---   
         wt <- w #ms Tuesday, June 22, 2004 at 23:16  
         xy <- xy.coords(x, y)
          y <- xy$y
          x <- xy$x
          n <- length(x)
          w <- if (is.null(w)) 
                  rep(1, n)
               else 
               {
               if (n != length(w)) 
               stop("lengths of x and w must match")
               if (any(w < 0)) 
               stop("all weights should be non-negative")
               if (all(w == 0)) 
               stop("some weights should be positive")
               (w * sum(w > 0))/sum(w)
               }
         W <- rep(1, n) # this for counting repeated x's
         x <- signif(x, 6)
       pox <- order(x)
        ux <- unique(sort(x)) # unique values of x
        ox <- match(x, ux)
       tmp <- matrix(
                    unlist(tapply(
                              seq(along = y), ox, 
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
                     byrow = TRUE
                     )
       wbar <- tmp[, 1] # weights for fitting 
       ybar <- tmp[, 2]/ifelse(wbar > 0, wbar, 1) #y/weights the y-variable for smooting
       yssw <- sum(tmp[, 3] - wbar * ybar^2)
    freq.ux <- tmp[,4] # frequency of unique x's
         nx <- length(ux) # no unique x values
       if (nx <= 3) 
        stop("need at least for unique `x' values")
       if (cv && nx < n) 
        warning("crossvalidation with non-unique `x' seems doubtful")
       r.ux <- ux[nx] - ux[1]
       xbar <- (ux - ux[1])/r.ux # 0<=xbar<=1
    if (all.knots) 
    {
        knot <- c(rep(xbar[1], 3), xbar, rep(xbar[nx], 3))
        nk <- nx + 2
    }
    else 
    {
        knot <- sknotl(xbar, nknots)
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
        if (df > 1 && df <= nx) 
           {
             icrit <- 3
            dofoff <- df
           }
        else warning(paste("you must supply 1 < df <= n,  n = #{unique x} =", 
            nx))
       }
    iparms <- as.integer(c(icrit, ispar, contr.sp$maxit))
    names(iparms) <- c("icrit", "ispar", "iter")
    fit <- .Fortran("qsbart", as.double(penalty), as.double(dofoff), 
                    x = as.double(xbar), y = as.double(ybar), w = as.double(wbar), 
                    ssw = as.double(yssw), as.integer(nx), as.double(knot), 
                    as.integer(nk), coef = double(nk), ty = double(nx), 
                    lev = double(nx), crit = double(1), iparms = iparms, 
                    spar = spar, parms = unlist(contr.sp[1:4]), 
                    isetup = as.integer(0), scrtch = double((17 + nk) * nk), 
                    ld4 = as.integer(4), ldnk = as.integer(1), ier = integer(1), 
                    DUP = FALSE, PACKAGE = "stats")[c("coef", "ty", "lev", 
                    "spar", "parms", "crit", "iparms", "ier")]
    wbar <- tmp[, 1] # new added by MS Tuesday, April 29, 2008 
    lev <- fit$lev
    df <- sum(lev)
    if (is.na(df))   
       stop("NA lev[]; probably smoothing parameter `spar' way too large!")
    if (abs(df-dofoff)>.5) # ms Saturday, February 7, 2004 at 16:06 it used to be 1
    if (ispar != 1)
      {warning("The output df ", df,"  are different for the input ", dofoff,
             "\n","change the control.spar" )
      }       
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
            fit$ty <- rep(mean(y), nx)
            df <- 1
            warning(paste(wtxt, "setting df = 1  __use with care!__", 
                sep = "\n"))
         }
       }
    cv.crit <- if (cv) 
                 {
                  ww <- wbar
                  ww[!(ww > 0)] <- 1
                  weighted.mean(((y - fit$ty[ox])/(1 - (lev[ox] * w)/ww[ox]))^2, 
                  w)
                 }
               else weighted.mean((y - fit$ty[ox])^2, w)/(1 - (df.offset + 
                     penalty * df)/n)^2
      pen.crit <- sum(wbar * (ybar - fit$ty) * ybar)
        longy  <- rep(fit$ty,freq.ux)
    longy[pox] <- longy 
          llev <- rep(lev,freq.ux) # the leverage   calculations    
     llev[pox] <- llev           
          sumw <- rep(wbar,freq.ux)     
     sumw[pox] <- sumw           
       longlev <- llev*w/sumw   
          levl <- (longlev-.hat.WX(w,x))
           var <- levl/wt # MS Tuesday, June 22, 2004 at 20:58
if (is.null(xeval)) # if no prediction  
    {
       obj.out <- list(residuals=y-longy, fitted.values = longy, 
                   var = var,  nl.df = df-2, lambda = fit$spar, 
                   coefSmo = list(knot = knot, nk = nk, min = ux[1], range = r.ux, 
                   coef = fit$coef, pen.crit = pen.crit, lev = levl,
                   lambda1 = unname(fit$parms["low"]) ))
class(obj.out) <- "smooth.spline"
    obj.out
    }
else { # if  prediction 
         fit.object <- list(knot = knot, nk = nk, min = ux[1], range = r.ux, 
                    coef = fit$coef)
  class(fit.object) <- "smooth.spline.fit"
             object <- list(x = ux, y = fit$ty, w = wbar, yin = ybar, 
                           lev = lev, cv.crit = cv.crit, pen.crit = pen.crit, 
                           crit = fit$crit, df = df, spar = fit$spar, 
                           lambda = unname(fit$parms["low"]), 
                           iparms = fit$iparms, fit = fit.object, 
                           call = match.call())
       class(object) <- "smooth.spline" 
                pred <- predict(object,x = xeval)
                pred$y  
     }    
}




    
