#----------------------------------------------------------------------------------------
# MS + BR last change Friday, February 2, 2007 
NOF <- function (mu.link="identity", sigma.link="log", nu.link ="identity")
{
    mstats <- checklink("mu.link", "normal Family", substitute(mu.link), c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "normal Family", substitute(sigma.link), c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "normal Family", substitute(nu.link), c("1/mu^2", "log", "identity"))
    
    structure(
          list(family = c("NOF", "normal Family"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                nopar = 3, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)), 
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv, 
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta, 
                 dldm = function(y,mu,sigma,nu) {   
                                 #     mu1 <- mu
                                 #  sigma1 <- sigma*mu^(nu/2)
                                 #  dldm1 <- (1/sigma1^2)*(y-mu1) 
                                 #  dldd1 <- ((y-mu1)^2-sigma1^2)/(sigma1^3)
                                 #   dldm2 <- dldm1+dldd1*sigma*(nu/2)*mu^((nu/2)-1) 
                                 #dldm
             dldm <- -(nu/(4*mu))+(y-mu)/((sigma^2)*(mu^nu))+(((y-mu)^2)*nu)/(2*(sigma^2)*mu^(nu+1))
                                   dldm
                                    },
               d2ldm2 = function(mu,sigma,nu) {
                                     # mu1 <- mu
                                   sigma1 <- sigma*mu^(nu/2)
                                  d2ldm21 <- -(1/sigma1^2)
                                  d2ldd21 <- -(2/(sigma1^2))
                                   d2ldm2 <- d2ldm21+d2ldd21*((sigma*nu/2)^2)*mu^(nu-2)
                                   d2ldm2
                                   },
                 dldd = function(y,mu,sigma,nu) {
                                 #    mu1 <- mu
                                 # sigma1 <- sigma*mu^(nu/2) 
                                 #  dldd1 <- ((y-mu1)^2-sigma1^2)/(sigma1^3)
                                 #  dldd2  <- dldd1*mu^(nu/2) 
                                   dldd <- -1/(sigma)+((y-mu)^2)/((sigma^3)*(mu^nu))
                                   dldd  
                                    },
               d2ldd2 = function(mu,sigma,nu) {
                                  #    mu1 <- mu
                                   sigma1 <- sigma*mu^(nu/2)
                                  d2ldd21 <- -(2/(sigma1^2))
                                   d2ldd2 <- d2ldd21*mu^(nu)
                                   d2ldd2
                                    },
                 dldv = function(y,mu,sigma,nu) {
                                   #dldv1 <- -0.5*log(mu) + (((y-mu)*log(mu)/sigma)^2)/(2*(mu^nu))
                                   dldv <- -0.5*log(mu)+(((y-mu)^2)*log(mu))/(2*(sigma^2)*mu^nu)
                                   dldv     
                                    },
               d2ldv2 = function(mu)  {
                                 d2ldv2 <- -0.5*(log(mu))^2
                                 d2ldv2  
                                    },
              d2ldmdd = function(y,mu,sigma,nu)  {
                                   #   mu1 <- mu
                                   sigma1 <- sigma*mu^(nu/2)
                                  d2ldd21 <- -(2/(sigma1^2))
                                  d2ldmdd <- d2ldd21*sigma*(nu/2)*mu^(nu-1)
                                  d2ldmdd
                                    },
              d2ldmdv = function(mu,nu)  {
                                  d2ldmdv <- -nu*log(mu)/(2*mu)
                                  d2ldmdv
                                    },
              d2ldddv = function(mu,sigma)  {
                                  d2ldddv <- -log(mu)/sigma
                                  d2ldddv
                                    },
          G.dev.incr  = function(y,mu,sigma,nu,...) -2*dNOF(y,mu,sigma,nu,log=TRUE),                           
                rqres = expression(
                rqres(pfun="pNOF", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu) 
                                   ),
            mu.initial =  expression(  mu <- (y+mean(y))/2), 
         sigma.initial =  expression({ sigma <- rep((0.2*sd(y))/sqrt(mean(y)),length(y)) }), 
            nu.initial =  expression(   nu <- rep(1, length(y))), 
              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) TRUE, 
               y.valid = function(y)  all(y>0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dNOF<-function(x, mu=0, sigma=1, nu=0, log=FALSE)
 { 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    #if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
               mu1 <- mu
            sigma1 <- sigma*mu^(nu/2)
               fy1 <- dnorm(x, mean=mu1, sd=sigma1, log=log)
     # fy <- -0.5*(log(2*pi))-log(sigma)-(nu/2)*log(mu)-((x-mu)^2)/(2*(sigma^2)*(mu^nu))
     # fy <- if (log==TRUE) fy else exp(fy) 
    fy1 
  }
#---------------------------------------------------------------------------------------- 
pNOF <- function(q, mu=0, sigma=1, nu=0, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    #if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
               mu1 <- mu
            sigma1 <- sigma*mu^(nu/2)
            cdf <- pnorm(q, mean=mu1, sd=sigma1, lower.tail = lower.tail, log.p = log.p)
      cdf
   }
#----------------------------------------------------------------------------------------
qNOF <- function(p, mu=0, sigma=1, nu=0, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", ""))  
    #if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
      mu1 <- mu
            sigma1 <- sigma*mu^(nu/2)
    q <- qnorm(p, mean=mu1, sd=sigma1, lower.tail = lower.tail )
    q
   }
#----------------------------------------------------------------------------------------
rNOF <- function(n, mu=0, sigma=1, nu=0)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
   # if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qNOF(p, mu=mu, sigma=sigma, nu=nu)
    r
  }
