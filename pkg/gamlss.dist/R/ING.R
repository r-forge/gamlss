######################################
##### Inverse Gamma Distribution #####
######################################

##ING (alpha, mu/(alpha+1))?NO
# where alpha = 1/sigma^2
# x>0, alpha > 0, mu > 0

#Probability Density function
dING <- function(x, mu = 2, sigma = 1, log = FALSE)
{
   if (any(mu < 0))
	stop(paste("mu must be greater than 0", "\n", ""))
   if (any(sigma <= 0))
	stop(paste("sigma must be greater than 0", "\n", ""))
   if (any(x < 0))
	stop(paste("x must be greater than 0", "\n", ""))
   alpha <- 1/(sigma^2)
   lfy <- alpha*log(mu) + alpha*log(alpha+1) - lgamma(alpha) -
	  (alpha + 1)*log(x) - ((mu*(alpha + 1))/x)
   if (log == FALSE) fy <- exp(lfy)
   else fy <-lfy
   fy
}

#Cumulative density function
pING <- function(q, mu = 2, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
   if (any(mu <= 0))
	stop(paste("mu must be greater than 0", "\n", ""))
   if (any(sigma <= 0))
	stop(paste("sigma must be greater than 0", "\n", ""))
   if (any(q < 0))
	stop(paste("q must be greater than 0", "\n", ""))
   alpha <- 1/(sigma^2)
   lcdf <- pgamma(((mu*(alpha + 1))/q), alpha, lower=FALSE, log.p = TRUE)  
   #cdf <- (pgamma(((mu*(alpha + 1))/q), alpha, lower=FALSE)*gamma(alpha))/ gamma(alpha)
   if (lower.tail == TRUE) lcdf <- lcdf
   else lcdf <- 1 - lcdf
   if (log.p == FALSE) cdf <- exp(lcdf)
   else cdf <- lcdf 
   cdf
} 

#Quantile 
qING <- function(p, mu = 2, sigma = 1, lower.tail = TRUE, log.p = FALSE, 
                 max.value=1000)
{
   if (any(mu <= 0))
	stop(paste("mu must be greater than 0", "\n", ""))
   if (any(sigma <= 0))
	stop(paste("sigma must be greater than 0", "\n", ""))
   if (any(p < 0)|any(p > 1))
	stop(paste("p must be between 0 and 1", "\n", ""))  
    QQQ <- rep(0, length = length(p))   
   for (i in seq(along = p)) {
        cumpro <- 0
        if (p[i] + 1e-09 >= 1) 
            QQQ[i] <- Inf
        else {
            for (j in seq(from = 0, to = max.value, by=0.001)) {
                cumpro <- pING(j, mu = mu, sigma = sigma, 
				log.p = FALSE)
                QQQ[i] <- j
                if (p[i] <= cumpro) 
                  break
            }
        }
    }
    QQQ
}

#Random Generating
rING <- function (n, mu = 2, sigma = 1)
{
   if (any(mu <= 0))
	stop(paste("mu must be greater than 0", "\n", ""))
   if (any(sigma <= 0))
	stop(paste("sigma must be greater than 0", "\n", ""))
   if (any(n <= 0)) 
        stop(paste("n must be a positive integer", "\n", ""))  
   n <- ceiling(n)
   p <- runif(n)
   r <- qING(p, mu = mu, sigma = sigma)
   r 
}

#browser()

#Gamlss Family Function
ING <- function (mu.link = "log", sigma.link = "log") 
{
    mstats <- checklink("mu.link", "Inverse Gamma", substitute(mu.link), 
        c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Inverse Gamma", substitute(sigma.link), 
        c("inverse", "log", "identity", "own"))
    structure(list(family = c("ING", "Inverse Gamma"), 
    parameters = list(mu = TRUE, sigma = TRUE), 
    nopar = 2, type = "Continuous", 
    mu.link = as.character(substitute(mu.link)), 
    sigma.link = as.character(substitute(sigma.link)), 
    mu.linkfun = mstats$linkfun, 
    sigma.linkfun = dstats$linkfun, 
    mu.linkinv = mstats$linkinv, 
    sigma.linkinv = dstats$linkinv, 
    mu.dr = mstats$mu.eta, 
    sigma.dr = dstats$mu.eta, 
    dldm = function(y, mu, sigma){ 
       alpha <- 1/(sigma^2)
       dldm <- (alpha/mu) - ((alpha + 1)/y)
       dldm
    }, 
    d2ldm2 = function(y, mu, sigma){
       alpha <- 1/(sigma^2)
       d2ldm2 <- -(alpha/(mu^2)) 
       d2ldm2
    }, 
    dldd = function(y, mu, sigma){
        alpha <- 1/(sigma^2)
        dldd <- (-2/(sigma^3))*(log(mu) + (alpha/(alpha+1)) + log(alpha+1) - 
                 digamma(alpha) - log(y) - (mu/y))
        dldd
    }, 
    d2ldd2 = function(y, mu, sigma){
       # alpha <- 1/(sigma^2)
        #d2ldd2<- ((4*sigma^6)*((1/((alpha+1)^2)) + (1/(alpha+1)) -
        #         trigamma(alpha)))
        #d2ldd2
        alpha <- 1/(sigma^2)
        dldd <- (-2/(sigma^3))*(log(mu) + (alpha/(alpha+1)) + log(alpha+1) - 
                 digamma(alpha) - log(y) - (mu/y))
        d2ldd2<- -dldd^2
        d2ldd2
    }, 
    d2ldmdd = function(y, mu, sigma){ 
       alpha <- 1/(sigma^2)
       d2ldmdd <- -2/(mu*(sigma^3)*(alpha+1))
       d2ldmdd 
    }, 
    G.dev.incr = function(y, mu, sigma, ...) -2 * dING(y, mu, sigma, log = TRUE), 
    rqres = expression(rqres(pfun = "pING", type = "Continuous", y = y, 
                             mu = mu, sigma = sigma)), 
    mu.initial = expression({ mu <- rep(mean(y), length(y))}), 
    sigma.initial = expression({ sigma <- rep(((mean(y)^2)/var(y))+2, length(y)) }), 
    mu.valid = function(mu) all(mu > 0), 
    sigma.valid = function(sigma) all(sigma > 0), 
    y.valid = function(y) TRUE), 
    class = c("gamlss.family", "family"))
}






