#-----------------------------------------------
# a function to create lags 
# author MS
# last change 2-2-11 
#-----------------------------------------------
createLags <- function(y,lag=1, omit.na=FALSE)
{
 ## local function 	
 lag1 <- function(y)
   {
    l <- c(NA, y[1:(length(y)-1)])
    l
   }
  ## main function starts here
     d <- matrix(0, nrow=length(y), ncol=lag)
 d[,1] <- lag1(y)
 yname <- deparse(substitute(y))
 cname <- paste("lag1", yname, sep="")
 if (lag==1) 
  {names(d[,1]) <- c("lag1")
  if (omit.na) d <-na.omit(d)
  return(d)
  }
for(i in 2:lag) 
 {
 d[,i] <- lag1(d[,i - 1])
 cname <- c(cname, paste("lag", paste(yname,i, sep=""), sep=""))
 }
colnames(d ) <- cname
if (omit.na) d <-na.omit(d)
d
 }
#####
#------------------------------------------------
