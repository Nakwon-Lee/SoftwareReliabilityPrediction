MSE2 <- function(n,k,observed,estimated){
  
  if(length(observed)!=length(estimated)){
    stop("length inconsistency: length(observed)!=length(estimated)")
  }
  
  sum <- 0
  for(i in (n+1):k){
    sum <- sum + ((observed[i] - estimated[i])^2)
  }
  return(sum/(k-n))
} 

FC.MSE <- function(k,p,observed,estimated){
  obsrvd <- observed[1:k]
  estmtd <- estimated[1:k]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(Inf)
  }
  
  tsum <- sum((obsrvd-estmtd)^2)
  
  ret <-NULL
  
  if((length(obsrvd)-p)<1){
    ret <- tsum/1
  }else{
    ret <- tsum/(length(obsrvd)-p)
  }
  
  stopifnot(ret>=0)
  
  return(ret)
}

FC.MAE <- function(k,p,observed,estimated){
  obsrvd <- observed[1:k]
  estmtd <- estimated[1:k]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(Inf)
  }
  
  tsum <- sum(abs(obsrvd-estmtd))
  ret <-NULL
  
  if((length(obsrvd)-p)<1){
    ret <- tsum/1
  }else{
    ret <- tsum/(length(obsrvd)-p)
  }
  
  stopifnot(ret>=0)
  
  return(ret)
}

FC.Rsquare <- function(k,observed,estimated){
  obsrvd <- observed[1:k]
  estmtd <- estimated[1:k]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(-Inf)
  }
  
  thismean <- sum(obsrvd)/length(obsrvd)
  SST <- sum((obsrvd - thismean)^2)
  SSR <- sum((obsrvd-estmtd)^2)
  
  return(1-(SSR/SST))
}

FC.Noise <- function(k,frestimated){
  frestmtd <- frestimated[1:k]
  
  if(any(is.na(frestmtd))){
    return(NA)
  }
  
  if(any(is.infinite(frestmtd))){
    return(Inf)
  }
  
  fleft <- frestmtd[2:length(frestmtd)]
  fright <- frestmtd[1:(length(frestmtd)-1)]
  tsum <- 0
  for(i in 1:length(fleft)){
    if(!is.infinite(fright[i])){
      if(fright[i]>0){
        tsum <- tsum + abs((fleft[i] - fright[i])/fright[i])
      }else{
        tsum <- tsum + 1
      }
    }
  }
  return(tsum)
}

FC.Bias <- function(n,k,observed,estimated){
  obsrvd <- observed[n:k]
  estmtd <- estimated[n:k]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(Inf)
  }
  
  tsum <- 0
  count <- 0
  for (i in 1:length(obsrvd)) {
    if(!is.infinite(estmtd[i])){
      tsum <- tsum + (estmtd[i] - obsrvd[i])
      count <- count + 1
    }
  }
  
  if(count>0){
    return(tsum/count)
  }else{
    return(Inf)
  }
}

FC.Variation <- function(k,observed,estimated){
  obsrvd <- observed[1:k]
  estmtd <- estimated[1:k]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(Inf)
  }
  
  pbias <- FC.Bias(1,k,observed,estimated)
  
  tsum <- 0
  count <- 0
  for (i in 1:length(obsrvd)) {
    if(!is.infinite(estmtd[i])){
      tsum <- tsum + (((obsrvd[i] - estmtd[i]) - pbias)^2)
      count <- count + 1
    }
  }
  
  return(sqrt(tsum/(count-1)))
}

FC.PRR <- function(k,observed,estimated){
  obsrvd <- observed[1:k]
  estmtd <- estimated[1:k]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(Inf)
  }
  
  tsum <- 0
  for(i in 1:length(estmtd)){
    if(estmtd[i]>0){
      tsum <- tsum + (((estmtd[i] - obsrvd[i])/estmtd[i])^2)
    }else{
      tsum <- tsum + 1
    }
  }
  return(tsum)
}

FC.WLSE <- function(k,observed,estimated,pc){
  obsrvd <- observed[1:k]
  estmtd <- estimated[1:k]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(Inf)
  }
  
  ivec <- c(1:length(obsrvd))
  si <- sum(ivec)
  e <- (length(obsrvd)*(1-pc))/si
  expr1 <- pc + (e * ivec)
  ret <- sum(((estmtd - obsrvd)^2)*expr1)
  return(ret)
}

FC.EP <- function(k,observed,estimated){
  obsrvd <- observed[k]
  estmtd <- estimated[k]
  
  if(is.na(estmtd)){
    return(NA)
  }
  
  if(is.infinite(estmtd)){
    return(Inf)
  }
  
  return(abs(obsrvd - estmtd))
}

FC.MEOP <- function(n,k,observed,estimated){
  obsrvd <- observed[n:k]
  estmtd <- estimated[n:k]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(Inf)
  }
  tsum <- sum(abs(obsrvd - estmtd))
  return(tsum/(length(obsrvd)+1))
}

FC.TOUP <- function(t,observed,estimated){
  obsrvd <- observed[1:t]
  estmtd <- estimated[1:t]
  
  if(any(is.na(estmtd))){
    return(NA)
  }
  
  if(any(is.infinite(estmtd))){
    return(Inf)
  }
  alpha <- 0.5 #[0,1]
  m <- length(obsrvd) - 2
  ivec <- c(0:m)
  errvec <- abs(estmtd - obsrvd)
  expr1 <- alpha * ((1 - alpha)^ivec)
  errleft <- errvec[2:length(obsrvd)]
  errright <- errvec[1:(length(obsrvd)-1)]
  stopifnot(length(expr1)==length(errleft))
  tsum <- 0
  for(i in 0:m){
    tsum <- tsum + (expr1[(i+1)] * (errleft[(length(errleft)-i)] - errright[(length(errright)-i)]))
  }
  return(tsum)
}

GROUP <- function(pbias){
  
  if(is.na(pbias)){
    return(NA)
  }
  
  if(is.infinite(pbias)){
    return("O")
  }
  
  ret <- NULL
  if(pbias >= 0 ){
    ret <- "O"
  }else{
    ret <- "U"
  }
  return(ret)
}
