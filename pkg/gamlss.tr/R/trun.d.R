trun.d <-function(par, family = "NO", type = c("left", "right", "both"), ...)
  {
    type <- match.arg(type)
if (type=="both" && length(par)!= 2)  stop(paste("the length of par should be 2 \n"))
if (type!="both" && length(par)!= 1)  stop(paste("the length of par should be 1 \n"))
   fname <- family
   if (mode(family) != "character" && mode(family) != "name")
   fname <- as.character(substitute(family))
 distype <- eval(call(family))$type
    dfun <- paste("d",fname,sep="")
    pfun <- paste("p",fname,sep="")
     pdf <- eval(parse(text=dfun))
     cdf <- eval(parse(text=pfun))
fun <- if (type=="left")  
       function(x, log = FALSE, ...)
        {
        if (distype=="Discrete" &&  any(x <= par))  
          stop(paste("x must be greater than ", par, "\n", ""))
        if (distype!="Discrete" && any(x < par))
          stop(paste("x must be greater or equal than ", par, "\n", ""))
        dfun <- pdf(x,log = TRUE,...)-log(1-cdf(par,...))
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
       }
     else if (type=="right")
      function(x, log = FALSE, ...)
       {
        if (distype=="Discrete" &&  any(x >= par))  
          stop(paste("x must be less than ", par, "\n", ""))
        if (distype!="Discrete" && any(x > par))
          stop(paste("x must be less or equal than ", par, "\n", ""))   
        dfun <-  if (distype=="Discrete") pdf(x, log = TRUE,...)-log(cdf(par-1,...))
                 else pdf(x, log = TRUE,...)-log(cdf(par,...))
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
       } 
     else if (type=="both")    
      function(x, log = FALSE, ...)
       {
        if (distype=="Discrete" &&  (any(x <= par[1]) || any(x >= par[2])) )  
          stop(paste("x must be greater than", par[1], "and less than", par[2], 
                      "\n", ""))
        if (distype!="Discrete" && (any(x < par[1]) || any(x > par[2])) )
         stop(paste("x must be greater than", par[1], "and less or equal to", par[2], 
                      "\n", ""))  
        dfun <-   if (distype=="Discrete") pdf(x, log = TRUE,...) - log(cdf(par[2]-1,...)-cdf(par[1],...)) 
                 else pdf(x, log = TRUE,...) - log(cdf(par[2],...)-cdf(par[1],...))      
         #pdf(x, log = TRUE,...) - log(cdf(par[2],...)-cdf(par[1],...))
        dfun <- if (log == TRUE) dfun else exp(dfun)
        dfun
       }   
#    formals(fun) <- formals(dfun) 
#environment(fun) <- .GlobalEnv # it is not working since it can not find par
  fun
  }
