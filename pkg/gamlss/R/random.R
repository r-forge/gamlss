#-------------------------------------------------------------------------------------------
# this function is take from Hastie's S-plus gam 
random <- function(xvar, df = NULL, lambda = 0) 
{
 scall <- deparse(sys.call())
 if(!inherits(xvar, "factor")) # | !is.category(xvar))
    stop("random() expects a factor or category as its first argument")
 xvar <- C(xvar, rep(0, length(levels(xvar))), 1) # puts zero in the X matrix 
 attr(xvar, "call") <- substitute( gamlss.random(data[[scall]], z, w, df = df, lambda))
 class(xvar) <- c("smooth", class(xvar))
 xvar
}
#-------------------------------------------------------------------------------------------
# this function is take from Hastie's gam 
# last change MS Thursday, April 25, 2002 at 13:31
gamlss.random <- function(x, y, w, df = NULL, lambda = 0) 
{
    df.inv <- function(n, df, lambda = sum(n)/df - mean(n), iterations = 10)
       { # given df find lambda
        if(df > length(n))
            return(0)
        current.df <- sum(n/(n + lambda))
        if(abs((df - current.df)/df) < 0.0001 | iterations == 1)
            lambda
        else {
            lambda <- exp(log(lambda) + (current.df - df)/(sum((n * lambda)/(n +lambda)^2)))
            Recall(n, df, lambda, iterations - 1)
             }
        }
    nw <- tapply(w, x, sum)
    non.zero <- !is.na(nw)
    if(is.null(df))     df <- sum(non.zero)
    if(lambda == 0) lambda <- df.inv(nw[non.zero], df)
    df <- sum(nw[non.zero]/(nw[non.zero] + lambda))
    fit <- tapply(w * y, x, sum)/(nw + lambda)
    var <- as.vector(w/(nw[x] + lambda))
    residuals <- as.vector(y - fit[x])
    list(x = seq(along = nw), y = fit, residuals = residuals, var = var, 
         nl.df = df-1, lambda=lambda, coefSmo=lambda) # MS 
}
