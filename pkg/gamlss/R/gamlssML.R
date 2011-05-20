# this is an attempt to create a functon which 
# fits a GAMLSS distribition using non linear maximation like optim()
# in fact here we use MLE() which is a copy of the mle() function of stat4
# The reason fot doing this 
#     i) to impove the speed of fitting 
#     ii)  because it would be usefull to extent this to Hidden Markov Models
#          and other models as ARMA and GARCH type
# Author Mikis Stasinopoulos 
# TO DO
# i)   weights implementation OK 
# ii)  residuals  OK
# iii) do I need mu.fv and mu.lp?? probably not
# iv)  mu.coef etc OK
# v) if implemented may histDist() should use it (Not working at the moment)
# vi) check all the methods OK
# vii) do I need data argument OK
# viii) vcoc() ?
#######################################################################################
#names(m1)          
#[5]  "weights"                   
#        "iter"              "method"         
#[13] "converged"       "residuals"       "noObs"          
#[17] "mu.fv"           "mu.lp"            
#[21] "mu.link"              
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# this should be a GAMLSS based mle procedure for fitting distributions
# the main problem is to create the likelihood function automatcally 
# not to be given by the user as in MLE below
# also it would be better to use the link function parameters rather the original
# to avoid boundary problems in the parameters
# this probably can be done by 
#   i)    using the eta.par as parameters in the likelihood fun say LogLikelihood(tparameters)
#   ii)   within the function use the inverse link to go to the original parameters
#   iii)  evalute the likelihod 
#   vi)   after exit transfer back (also use this for the se's)?? 
#---------------------------------------------------------------------------------------- 
#require(gamlss)
gamlssML<-function(y, 
		     family = NO,  
		    weights = NULL, 
		   mu.start = NULL, 
        sigma.start = NULL, 
	       nu.start = NULL, 
		  tau.start = NULL,
	         mu.fix = FALSE,
		  sigma.fix = FALSE,
			 nu.fix = FALSE,
			tau.fix = FALSE,
			 data = NULL, 
			 ...) 
{
#------------------------------------------------------------------------------------------
# local functions
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# this is  copy of the stats::mle()
# but it exports a S3 object rather than an S4 object
# see  the examples in lme() for its use
# here we use it as convenient way of calling optim()
# and to cover the case that some parameters can be fixed
#-------------------------------------------------------------------------------------------
	MLE <- function (minuslogl, start = formals(minuslogl), method = "BFGS", 
			fixed = list(), ...) 
	{ 
		  call <- match.call()
		     n <- names(fixed)
	  fullcoef <- formals(minuslogl)
		if (any(!n %in% names(fullcoef))) 
			stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
   fullcoef[n] <- fixed
		if (!missing(start) && (!is.list(start) || is.null(names(start)))) 
			stop("'start' must be a named list")
	 start[n] <- NULL
		start <- sapply(start, eval.parent)
		   nm <- names(start)
		   oo <- match(nm, names(fullcoef))
		if (any(is.na(oo))) 
			stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
		start <- start[order(oo)]
		   nm <- names(start)
		    f <- function(p) 
		      {
			           l <- as.list(p)
			    names(l) <- nm
			        l[n] <- fixed
			    do.call("minuslogl", l)
		      }
  		 oout <- if (length(start)) 
					optim(start, f, method = method, hessian = TRUE, ...)
				else list(par = numeric(0L), value = f(start))
		 coef <- oout$par
		 vcov <- if (length(coef)) 
					solve(oout$hessian)
				else matrix(numeric(0L), 0L, 0L)
		 min <- oout$value
fullcoef[nm] <- coef
		#new("mle", call = call, coef = coef, fullcoef = unlist(fullcoef), 
		#    vcov = vcov, min = min, details = oout, minuslogl = minuslogl, 
		#    method = method)
		out <- list(call = call, coef = coef, fullcoef = unlist(fullcoef), 
				vcov = vcov, min = min, details = oout, minuslogl = minuslogl, 
				method = method)
class(out) <- "MLE"
		out
	}
# end of MLE
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# this is to replicate rqres within gamlss enviroment DS Friday, March 31, 2006 at 10:30
# it is used as in gamlss()
rqres <- function (pfun = "pNO", 
                   type = c("Continuous", "Discrete", "Mixed"),
               censored = NULL,  
                   ymin = NULL, 
                 mass.p = NULL, 
                prob.mp = NULL,
                      y = y,
                         ... )
{ }
body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# main function starts here
  mlFitcall <- match.call()  #   the function call    
       fam  <- as.gamlss.family(family)
       if (!is.null(data)) {attach(data); on.exit(detach(data))}
      fname <- fam$family[[1]] 
       dfun <- paste("d",fname,sep="")
	  # pfun <- paste("p",fname,sep="")
        PDF <- eval(parse(text=dfun))
	#	CDF <- eval(parse(text=pfun))
      nopar <- fam$nopar
          N <- length(y)
          w <- if(is.null(weights))    rep(1, N) else weights
               if(any(w < 0)) stop("negative weights not allowed") # 
# get the initial values if are not set by the user
           if ("mu"%in%names(fam$parameters))
           {
             mu <- if(is.null(mu.start))  mean(eval(fam$mu.initial, list(y=y)))
			       else mu.start 
               eta.mu <- fam$mu.linkfun(mu)
			  #eta.par <- c(eta.mu)          
           }
           if ("sigma"%in%names(fam$parameters))
           {
             sigma <- if(is.null(sigma.start)) mean(eval(fam$sigma.initial, list(y=y, mu=mu)))
			          else sigma.start
            eta.sigma <- fam$sigma.linkfun(sigma)
			  #eta.par <- c(eta.mu, eta.sigma)
           }
           if ("nu"%in%names(fam$parameters))
           {
			   nu <- if(is.null(nu.start))  mean(eval(fam$nu.initial, list(y=y, mu=mu, sigma=sigma)))
			         else nu.start
               eta.nu <- fam$nu.linkfun(nu)
			  #eta.par <- c(eta.mu, eta.sigma, eta.nu)
           }
           if ("tau"%in%names(fam$parameters))
           {  
			   tau <- if(is.null(tau.start)) mean(eval(fam$tau.initial, list(y=y, mu=mu, sigma=sigma, nu=nu)))
			   else tau.start
           	eta.tau <- fam$tau.linkfun(tau)
			#eta.par <- c(eta.mu, eta.sigma, eta.nu, eta.tau)
           }
# whether to fix parameters
   fixed <- list()
   if (mu.fix)      fixed <- c(fixed, eta.mu=mu)
   if (sigma.fix)   fixed <- c(fixed, eta.sigma=sigma)
   if (nu.fix)      fixed <- c(fixed, eta.nu=nu)
   if (tau.fix)     fixed <- c(fixed, eta.tau=tau)
                  noFixed <- sum(c(mu.fix, sigma.fix, nu.fix, tau.fix))
    switch(nopar, 
			{# one parameter 
				ll1 <- function(eta.mu)
				{
					mu <- fam$mu.linkinv(eta.mu)
					-sum(w*PDF(y, mu=mu, log=TRUE))
				}
			   fit <-	MLE(ll1, start=list(eta.mu=eta.mu), fixed=fixed, ...)
			},
			{# two paremeters
				ll2 <- function(eta.mu, eta.sigma)
				{
					mu <- fam$mu.linkinv(eta.mu)
				 sigma <- fam$sigma.linkinv(eta.sigma)
					-sum(w*PDF(y, mu=mu, sigma=sigma, log=TRUE))
				}
				fit <-	MLE(ll2, start=list(eta.mu=eta.mu, eta.sigma=eta.sigma), fixed=fixed, ...)
			},
			{# three parameters
				ll3 <- function(eta.mu, eta.sigma, eta.nu)
				{
					   mu <- fam$mu.linkinv(eta.mu)
					sigma <- fam$sigma.linkinv(eta.sigma)
					   nu <- fam$nu.linkinv(eta.nu)
					-sum(w*PDF(y, mu=mu, sigma=sigma, nu=nu, log=TRUE))
				}
				fit <-	MLE(ll3, start=list(eta.mu=eta.mu, eta.sigma=eta.sigma, eta.nu=eta.nu), fixed=fixed, ...)
				
			},
			{# four parameters
				ll4 <- function(eta.mu, eta.sigma, eta.nu, eta.tau)
				{
				    	mu <- fam$mu.linkinv(eta.mu)
				     sigma <- fam$sigma.linkinv(eta.sigma)
					    nu <- fam$nu.linkinv(eta.nu)
                       tau <- fam$tau.linkinv(eta.tau)
					-sum(w*PDF(y, mu=mu, sigma=sigma, nu=nu, tau=tau, log=TRUE))
				}
				fit <-	MLE(ll4, start=list(eta.mu=eta.mu, eta.sigma=eta.sigma, eta.nu=eta.nu, eta.tau=eta.tau), fixed=fixed, ...)
			}
	  )
    df.fit <- nopar-noFixed
     out <-list(family = fam$family,  
			parameters =  as.character(names(fam$par)), 
			      type = fam$type, 
				  call = mlFitcall, 
				     y = y,
			   weights = w,
	        G.deviance = 2*fit$min, 
			         N = N,  
				df.fit = df.fit, 
		   df.residual = N-df.fit, 
				   aic = 2*fit$min+2*nopar,
				   sbc = 2*fit$min+log(N)*nopar,
				method = "BFGS",
		   		  vcov = fit$vcov) 
		  
    if ("mu"%in%names(fam$parameters))
		   {
			    mu <- fam$mu.linkinv(fit$coef["eta.mu"])
				names(mu) <- "mu"
				mu.coefficients = fit$coef["eta.mu"]
				names(mu.coefficients) <- "mu.coefficients"
				
				out <- c(out, mu, mu.coefficients )
				out$mu.link <- fam$mu.link
		   }
		   ## Output for sigma model: ---------------------------------------------------------------
		   if ("sigma"%in%names(fam$parameters))
		   {
			    sigma <- fam$sigma.linkinv(fit$coef["eta.sigma"])
				names(sigma) <- "sigma"
				sigma.coefficients = fit$coef["eta.sigma"]
				names(sigma.coefficients) <- "sigma.coefficients"
				
				 out <- c(out, sigma, sigma.coefficients )
				 out$sigma.link <- fam$sigma.link
		   }
		   ##  output for nu ------------------------------------------------------------------------
		   if ("nu"%in%names(fam$parameters))
		   {
			     nu <- fam$nu.linkinv(fit$coef["eta.nu"])
				names(nu) <- "nu"
				nu.coefficients = fit$coef["eta.nu"]
				names(nu.coefficients) <- "nu.coefficients"			    
				 out <- c(out, nu, nu.coefficients )
				 out$nu.link <- fam$nu.link

		   }
		   ##  output for tau -----------------------------------------------------------------------
		   if ("tau"%in%names(fam$parameters))
		   {
			   tau <- fam$tau.linkinv(fit$coef["eta.tau"])
				names(tau) <- "tau"
				tau.coefficients = fit$coef["eta.tau"]
				names(tau.coefficients) <- "tau.coefficients"			    
				 out <- c(out, tau, tau.coefficients )
				 out$tau.link <- fam$tau.link
		   }
		    out$residuals <- eval(fam$rqres)
		        out$rqres <- fam$rqres
		  if (!is.null(data) ) out$call$data <- substitute(data)
     class(out) <- c("gamlssML", "gamlss")
	 out
 }                                        
    
######################################################################################
fitted.gamlssML<-function (object, what = c("mu", "sigma", "nu", "tau"), ... ) 
{
what <- match.arg(what)
if (! what%in%object$par) stop(paste(what,"is not a parameter in the gamlss object","\n"))
x <- rep(object[[what]], object$N)
x
}
######################################################################################
######################################################################################
vcov.gamlssML<-function (object, type = c("vcov", "cor", "se",  "all"), ... ) 
{
type <- match.arg(type)
switch(type,
       "se"={x <- sqrt(diag(object$vcov))},
       "cor"={x <- cov2cor(object$vcov)},
       "vcov"={x<- object$vcov},
       	"all"={x <-list(vcov=object$vcov, cor=cov2cor(object$vcov), se=sqrt(diag(object$vcov)))}
       	)
x
}
######################################################################################
  

   
 
 


