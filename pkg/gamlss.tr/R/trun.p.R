trun.p <- function(par, family = "NO", type = c("left", "right", "both"), ...)
  {
   type <- match.arg(type)
if (type=="both" && length(par)!= 2)  stop(paste("the length of par should be 2 \n"))
if (type!="both" && length(par)!= 1)  stop(paste("the length of par should be 1 \n"))  
  fname <- family
  if (mode(family) != "character" && mode(family) != "name")
  fname <- as.character(substitute(family))
distype <- eval(call(family))$type
   pfun <- paste("p",fname,sep="")
    cdf <- eval(parse(text=pfun))
fun <- if (type=="left")  
       function(q, lower.tail = TRUE, log.p = FALSE, ...)
      {
      if (distype=="Discrete" &&  any(q <= par))  
          stop(paste("q must be greater than ", par, "\n", ""))
      if (distype!="Discrete" && any(q < par))
          stop(paste("q must be greater or equal than ", par, "\n", ""))
      cof <- (cdf(q,...)-cdf(par,...))/(1-cdf(par,...))
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof
      }
     else if (type=="right")
       function(q, lower.tail = TRUE, log.p = FALSE, ...)
      {
      if (distype=="Discrete" &&  any(q >= par))  
          stop(paste("q must be less than ", par, "\n", ""))
        if (distype!="Discrete" && any(q > par))
          stop(paste("q must be less or equal than ", par, "\n", ""))
      cof <- if (distype=="Discrete") cdf(q,...)/cdf(par-1,...)
             else cdf(q,...)/cdf(par,...) # added Friday, February 26, 2010 
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof
      }
     else if (type=="both")    
      function(q, lower.tail = TRUE, log.p = FALSE, ...)
      {
     if (distype=="Discrete" &&  (any(q <= par[1]) || any(q >= par[2])) )  
          stop(paste("q must be greater than", par[1], "and less than", par[2], "\n", ""))
        if (distype!="Discrete" && (any(q < par[1]) || any(q > par[2])) )
         stop(paste("q must be greater or equal than", par[1], "and less or equal to", par[2], "\n", ""))  
      cof <- if (distype=="Discrete") (cdf(q,...)-cdf(par[1],...))/(cdf(par[2]-1,...)-cdf(par[1],...)) 
             else (cdf(q,...)-cdf(par[1],...))/(cdf(par[2],...)-cdf(par[1],...))   # added Friday, February 26, 2010     
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof
      }
  fun
  }
