trun.r <- function(par, family = "NO", type = c("left", "right", "both"), ...)
  { 
   type <- match.arg(type)
if (type=="both" && length(par)!= 2)  stop(paste("the length of par should be 2 \n")) 
if (type!="both" && length(par)!= 1)  stop(paste("the length of par should be 1 \n")) 
  fname <- family
  if (mode(family) != "character" && mode(family) != "name")
  fname <- as.character(substitute(family))
   qfun <- paste("q",fname,sep="")
   pfun <- paste("p",fname,sep="")
 invcdf <- eval(parse(text=qfun))
    cdf <- eval(parse(text=pfun))
fun <- if (type=="left")  
       function(n,...)
    {
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
     n <- ceiling(n)
     p <- runif(n)
    pp <- cdf(par,...)
     r <- invcdf((pp+p*(1-pp)),...)
     r
    }
     else if (type=="right")
     function(n,...)
    {
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
     n <- ceiling(n)
     p <- runif(n)
    pp <- cdf(par,...)
     r <- invcdf(p*pp,...)
     r
    }
     else if (type=="both")    
      function(n,...)
    {
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))   
     n <- ceiling(n)
     p <- runif(n)
   pp1 <- cdf(par[1],...)
   pp2 <- cdf(par[2],...)
     r <- invcdf(p*(pp2-pp1)+pp1,...)
     r
    }
  }
