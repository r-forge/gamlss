######################################
##### PARETO TYPE 2 DISTRIBUTION #####
######################################

#Probability distribution
dPA2 <- function(x, mu = 2, sigma = 0.5, log = FALSE)
{
   if (any(mu <= 0))
        stop(paste("mu must be between greater than 0", "\n", ""))
   if (any(sigma <= 0))
        stop(paste("sigma must be between greater than 0", "\n", ""))
   if (any(x < 0)) 
        stop(paste("x must be greater than 0", "\n", ""))
   lfy <- -log(sigma) + (1/sigma)*log(mu) - (1 + (1/sigma))*log(x+mu)
   if (log == FALSE) fy <- exp(lfy)
   else fy <- lfy
   fy
}

#Density 
pPA2 <- function(q, mu = 2, sigma = 0.5, lower.tail = TRUE, log.p = FALSE)
{
   if (any(mu <= 0))
        stop(paste("mu must be between greater than 0", "\n", ""))
   if (any(sigma <= 0))
        stop(paste("sigma must be between greater than 0", "\n", ""))
   if (any(q < 0)) 
        stop(paste("q must be be greater than 0", "\n", ""))   
   cdf <- 1 - (mu/(mu+q))^(1/sigma)
   if (lower.tail == TRUE) cdf <- cdf  
   else cdf <- 1 - cdf
   if (log.p == FALSE) cdf <- cdf
   else cdf < - log(cdf)
   cdf
}   

#Quantile 
qPA2 <- function(p, mu = 2, sigma = 0.5, lower.tail = TRUE, log.p = FALSE, 
		 max.value=1000)
{
   if (any(mu <= 0))
        stop(paste("mu must be between greater than 0", "\n", ""))
   if (any(sigma <= 0))
        stop(paste("sigma must be between greater than 0", "\n", ""))
   if (any(p < 0) | any(p > 1)) 
        stop(paste("p must be between 0 and 1", "\n", ""))
   QQQ <- rep(0, length = length(p))   
   for (i in seq(along = p)) {
        cumpro <- 0
        if (p[i] + 1e-09 >= 1) 
            QQQ[i] <- Inf
        else {
            for (j in seq(from = 0, to = max.value, by=0.001)) {
                cumpro <- pPA2(j, mu = mu, sigma = sigma, 
				log.p = FALSE)
                QQQ[i] <- j
                if (p[i] <= cumpro) 
                  break
            }
        }
    }
    QQQ
}

#Quantile 
qPA2 <- function(p, mu = 2, sigma = 0.5, lower.tail = TRUE, log.p = FALSE)
{
   if (any(mu <= 0))
        stop(paste("mu must be between greater than 0", "\n", ""))
   if (any(sigma <= 0))
        stop(paste("sigma must be between greater than 0", "\n", ""))
   if (any(p < 0) | any(p > 1)) 
        stop(paste("p must be between 0 and 1", "\n", ""))
   if (log.p == TRUE) 
        p <- exp(p)
   else p <- p
   if (lower.tail == TRUE) 
        p <- p
   else p <- 1 - p
   q <- mu*((1-p)^(-sigma)-1) 
   q   
}

#Random generation 
rPA2 <- function(n, mu = 2, sigma = 0.5)
{
   if (any(mu <= 0))
        stop(paste("mu must be between greater than 0", "\n", ""))
   if (any(sigma <= 0))
        stop(paste("sigma must be between greater than 0", "\n", ""))
   if (any(n <= 0)) 
        stop(paste("n must be a positive integer", "\n", ""))  
   n <- ceiling(n)
   p <- runif(n)
   r <- qPA2(p, mu = mu, sigma = sigma)
   r 
}

#Gamlss Family Function
PA2 <- function (mu.link = "log", sigma.link = "log") 
{
    mstats <- checklink("mu.link", "Pareto Type 2", substitute(mu.link), 
        c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Pareto Type 2", substitute(sigma.link), 
        c("inverse", "log", "identity", "own"))
    structure(list(family = c("PA2", "Pareto Type 2"), 
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
        dldm = function(y, mu, sigma) (1/(mu*sigma))-(1+(1/sigma))/(y+mu), 
        d2ldm2 = function(y, mu, sigma) -((1/(mu*sigma))-(1+(1/sigma))/(y+mu))*
                                        ((1/(mu*sigma))-(1+(1/sigma))/(y+mu)), 
        dldd = function(y, mu, sigma) -(1/sigma)-(1/sigma^2)*log(mu) + 
					1/sigma^2*log(y+mu), 
        d2ldd2 = function(y, mu, sigma) -(-(1/sigma)-(1/sigma^2)*log(mu) + 
					1/sigma^2*log(y+mu))*(-(1/sigma)-(1/sigma^2)*log(mu) + 
					1/sigma^2*log(y+mu)), 
        d2ldmdd = function(y, mu, sigma) -((1/(mu*sigma))-(1+(1/sigma))/(y+mu))*(-(1/sigma)-(1/sigma^2)*log(mu) + 
					1/sigma^2*log(y+mu)), 
        G.dev.incr = function(y, mu, sigma, ...) -2 * 
            dPA2(y, mu, sigma, log = TRUE), 
        rqres = expression(rqres(pfun = "pPA2", 
            type = "Continuous", y = y, mu = mu, sigma = sigma)), 
        mu.initial = expression({mu <- (y + mean(y))/2}), 
        sigma.initial = expression({sigma <- rep(sd(y), length(y))}), 
        mu.valid = function(mu) all(mu > 0), 
        sigma.valid = function(sigma) all(sigma > 0), 
        y.valid = function(y) TRUE), 
        class = c("gamlss.family", "family"))
}
