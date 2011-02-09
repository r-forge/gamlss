##---------------------------------------------------------------------------------------
## the lo() and gamlss.lo() functions
## are based on R loess() and S-plus lo() and gam.lo() functions
## author Mikis Stasinopoulos
## note that for more that one variable say lo(x1,x2,x3)
## gamlss() and Splus gam() give identical results
## with one variable say lo(x1) the result are differ
## but usually this is reflected into the df's  
##---------------------------------------------------------------------------------------
lo <-function(..., span = 0.5, df = NULL, degree = 1) 
 {
    lodummy <- function(span = 0.5, df = NULL, degree = 1)
    list(span = span, df = df, degree = degree)
    vars <- list(...)
    locall <- sys.call()
    chcall <- deparse(locall)
    mcall <- match.call(expand = FALSE)
    wspan <- !is.null(mcall$span) # 
    mcall$... <- NULL
    nvars <- length(vars)
    if(nvars > 1) 
    {
## a bit of freedom in giving the span and degree
        scalars <- sapply(vars, length) == 1
        if(any(scalars)) 
        {
            nvars <- nvars - sum(scalars)
            mcall <- c(mcall, as.call(vars[scalars]))
            vars <- vars[!scalars]
        }
    }
    mcall[[1]] <- as.name("lodummy")
    m <- eval(mcall)
    degree <- m$degree
    span <- m$span
    df   <- m$df
    if(degree > 2)
        stop("degrees 1 or 2 are implemented")
    if(nvars == 1) 
    {
        xvar <- vars[[1]]
        xnames <- deparse(locall[[2]])
    }
    else 
    {
        nobs <- length(vars[[1]])
        xvar <- matrix(0, nobs, nvars)
        xnames <- character(nvars)
        for(i in seq(nvars)) 
        {
            tt <- vars[[i]]
            if(!is.null(dd <- dim(tt)) && dd[2] > 1)
                stop(
        "either call lo with a matrix argument, or else a comma separated list x1, x2"
                  )
            exptt <- locall[[i + 1]]
            xnames[i] <- deparse(exptt)
            xvar[, i] <- as.numeric(tt)
        }
        dimnames(xvar) <- list(NULL, xnames)
    }
    xvar <- poly.matrix(xvar, degree = degree)
    pd <- attr(xvar, "degree")
    opd <- order(pd)
    if(length(pd) > 1) 
    {
        xvar <- xvar[, opd]
        p <- sum(pd == 1)
    }
    else p <- 1
    nobs <- length(xvar)[1]
    if(is.matrix(xvar))
    nas <- is.na(xvar[, 1:p])
    else nas<-is.na(xvar)
    if(any(nas)) 
    {
        if(p > 1)
            nas <- nas %*% array(1, c(p, 1))
        attr(xvar, "NAs") <- seq(nobs)[nas > 0]
    }
    if(span * nobs < 1)
        stop(paste("span is too small; the minimum is 1/n =", format(round(1/
            nobs, 4))))
    real.call <- substitute(gamlss.lo(data[[chcall]], z, w, span = span, df=df, 
                      degree = degree, ncols = p,  wspan=wspan), 
                      list(span = span, df=df, degree = degree, chcall = chcall, 
                      p = p, wspan=wspan))
      atts <- c(attributes(xvar), list(span = span, degree = degree, df=df, 
              ncols = p, wspan=wspan, call = real.call, class = c("smooth", "matrix")))
    attributes(xvar) <- atts
    xvar
 }
##---------------------------------------------------------------------------------------
## created by MS : Thursday, September 12, 2002 at 14:46
## based on Brian Ripley loess functions
## in this implementation of loess x must be a lo() object and ncols must be given
## note that the available options for loess.control are
## surface="interpolate", statistics="approximate",  trace.hat ="exact"
## compare to the available  
## surface = c("interpolate", "direct"),
## statistics = c("approximate", "exact"),
## trace.hat = c("exact", "approximate"),
gamlss.lo <- function (x, y, w = NULL, span , df=NULL, 
                      degree = 1, ncols=FALSE, wspan=TRUE, 
                      parametric = FALSE, # w=rep(1, length(y))
                      drop.square = FALSE, normalize = FALSE, family="gaussian", 
                      method = "loess", control = loess.control(...), 
                      xeval = NULL, ...)
{
 if (!length(y)) stop("invalid `y'")
 if(ncols==FALSE) stop("ncols must be given")
## checking the structure os x
 if (any(sapply(x, is.factor))) 
        stop("predictors must all be numeric") 
# wt <- w # save for later 
## checking the dimensions of  and picking only the fist order polynomial 
 if(is.null(np <- dim(x))) 
    {
        N <- as.integer(length(x))
        p <- as.integer(1)
        x <- as.matrix(x)
        attr(x,"dimnames")[2]<-list("x")
    }
    else 
    {
       np <- as.integer(np)
        N <- np[1]
        p <- ncols
        x <- as.matrix(x[,1:p])
        if (ncols==1) {
                      attr(x,"dimnames")[2]<-list("x")
                      }
    }
        D <- ncol(x)
    if (!match(degree, 0:2, 0)) 
        stop("degree must be 0, 1 or 2")
    if (!is.null(df)) 
        if (wspan) 
            warning("both span and df specified: span will be used")
        else 
        {
            tau <- switch(degree + 1, 1, D + 1, (D + 1) * (D + 
                2)/2) - sum(drop.square)
            span <- 1.2 * tau/df
        }
    #---------------------------------------
              cell <- control$cell
           surface <- control$surface
        statistics <- control$statistics
         trace.hat <- control$trace.hat
        iterations <- 1 # always one in GAMLSS 
            max.kd <- max(N, 200)
            robust <- rep(1, N)
               nmx <- colnames(x)
        names(nmx) <- nmx
       drop.square <- match(nmx, nmx[drop.square], 0) > 0
        parametric <- match(nmx, nmx[parametric], 0) > 0 
      sum.drop.sqr <- sum(drop.square)
    sum.parametric <- sum(parametric)
     nonparametric <- sum(!parametric)
  order.parametric <- order(parametric)
    order.drop.sqr <- (2 - drop.square)[order.parametric]
    if (iterations) 
        for (j in 1:iterations) 
        {
            robust <-  robust # w *
            if (j > 1) 
                statistics <- "none"
            else if (surface == "interpolate" && statistics == 
                "approximate") 
                statistics <- if (trace.hat == "exact") 
                  "1.approx"
                else "2.approx"
            surf.stat <- paste(surface, statistics, sep = "/")
            fit <- .C("loess_raw", as.double(y), as.double(x), 
                as.double(w), as.double(robust), as.integer(D), 
                as.integer(N), as.double(span), as.integer(degree), 
                as.integer(nonparametric), as.integer(order.drop.sqr), 
                as.integer(sum.drop.sqr), as.double(span * cell), 
                as.character(surf.stat), fitted.values = double(N), 
                parameter = integer(7), a = integer(max.kd), 
                xi = double(max.kd), vert = double(2 * D), vval = double((D + 
                  1) * max.kd), diagonal = double(N), trL = double(1), 
                delta1 = double(1), delta2 = double(1), as.integer(surf.stat == 
                  "interpolate/exact"), PACKAGE = "stats")         
        }       
##  before we get out we want to map the fitted values to the space orthogonal to M(X)
## This is also what gam() function in S-plus is doing
## MS weights are new Friday, April 4, 2003 
#cat("the df",fit$trL,"\n")
#cat("fitted fit", fitted(fit),  " \n")
#print(x)
 linfit <-  if(dim(x)[2]==1)    
                  lm(fitted(fit)~poly.matrix(x[,1],degree=degree), weights=w)
                else 
                  lm(fitted(fit)~poly.matrix(x,degree=degree),weights=w)
# cat("coef ", coef(linfit),  " \n")
if (is.null(xeval)) # if no prediction
    {
         lev <- (fit$diagonal-hat(linfit$qr)) 
         var <- lev/w
         out <- list(fitted=resid(linfit), residuals=y-resid(linfit), var=var, 
                nl.df=fit$trL-linfit$rank, coefSmo=list(lev=lev), lambda=span) 
                # ms Tuesday, March 25, 2003 at 18:00
  out
    }
else 
   { 
    if (surface == "interpolate") 
      {
        pars <- fit$parameter
        names(pars) <- c("d", "n", "vc", "nc", "nv", "liv", "lv")
        enough <- (D + 1) * pars["nv"]
          a <- fit$a[1:pars[4]] 
         xi <- fit$xi[1:pars[4]]
       vert <- fit$vert
       vval <- fit$vval[1:enough]
        #kd <- list(parameter = pars, a = fit$a[1:pars[4]], 
        #    xi = fit$xi[1:pars[4]], vert = fit$vert, vval = fit$vval[1:enough])
       M <- NROW(xeval)
       newx <- as.matrix(xeval) 
       inside <- matrix(TRUE, M, ncol = D)
        #ranges <- apply(x, 2, range)
        #inside <- (xeval <= rep(ranges[2, ], rep(M, D))) & 
        #    (xeval >= rep(ranges[1, ], rep(M, D)))
        inside <- inside %*% rep(1, D) == D
        M1 <- sum(inside)
        nfit <- rep(as.numeric(NA), M)
        if (any(inside)) 
            nfit[inside] <- .C("loess_ifit", as.integer(pars), 
                as.integer(a), as.double(xi), as.double(vert), 
                as.double(vval), as.integer(M1), as.double(xeval[inside, 
                  ]), fit = double(M1), PACKAGE = "stats")$fit               
      }
  if (surface == "direct") 
      {
        M <- NROW(xeval)
            fit <- .C("loess_dfit", as.double(y), as.double(x), 
                as.double(xeval), as.double(w), 
                as.double(span), as.integer(degree), as.integer(nonparametric), 
                as.integer(order.drop.sqr), as.integer(sum.drop.sqr), 
                as.integer(D), as.integer(N), as.integer(M), 
                fit = double(M), PACKAGE = "stats")$fit
        }
  plin <- cbind(1,xeval[,, drop = FALSE]) %*% coef(linfit)
  fit  <- nfit-plin  
  fit
   }
}
