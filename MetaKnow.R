#observed: number of observed failures for each time step
Variance <- function(observed){
  
  meanf <- mean(observed)
  nums <- length(observed)
  
  tsum <- sum(((observed-meanf)^2))
  
  return((tsum/(nums-1)))
}

#observed: number of observed failures for each time step
#pi: time variable
Inclination <- function(pi,observed){
  
  stopifnot(length(observed)==length(pi))
  
  obmin <- min(observed)
  obmax <- max(observed)
  
  normvec <- ((observed-obmin)/(obmax-obmin))
  
  stopifnot(length(observed)==length(normvec))
  
  nums <- length(normvec)
  meanf <- mean(normvec)
  meani <- mean(pi)
  
  sumtop <- sum(((normvec-meanf)*(pi-meani)))
  sumbot <- sum(((normvec-meanf)^2))
  
  return(sumtop/sumbot)
}

#observed: number of observed failures for each time step
#pi: time variable
AutoCorr <- function(pi,observed){
  
  stopifnot(length(observed)==length(pi))
  
  obmin <- min(observed)
  obmax <- max(observed)
  
  normvec <- ((observed-obmin)/(obmax-obmin))
  
  stopifnot(length(observed)==length(normvec))
  
  nums <- length(normvec)
  meanx <- mean(normvec[c(2:length(normvec))])
  meany <- mean(normvec[c(1:(length(normvec)-1))])
  
  sumx <- sum(((normvec[c(2:length(normvec))]-meanx)^2))
  varx <- sumx/(nums-1)
  
  sumy <- sum(((normvec[c(1:(length(normvec)-1))]-meany)^2))
  vary <- sumy/(nums-1)
  
  sumf <- 0
  for(i in 2:length(normvec)){
    sumf <- sumf + ((normvec[i]-meanx)*(normvec[i-1]-meany))
  }
  
  return((sumf/nums)/(varx*vary))
}

#observed: number of observed failures for each time step
MetaNO <- function(observed){
  
  nums <- length(observed)
  stopifnot(nums>=3)
  
  obmin <- min(observed)
  obmax <- max(observed)
  
  normvec <- ((observed-obmin)/(obmax-obmin))
  
  stopifnot(length(observed)==length(normvec))
  
  term <- 5
  
  if (nums<11){
    term <-((nums-1) %/% 2)
  }
  
  sumf <- 0
  
  for(i in (term+1):(length(normvec)-term)){
    sump <- 0
    for(j in (i-term):(i+term)){
      sump <- sump + (normvec[j]/(term*2))
    }
    sumf <- sumf + ((sump - normvec[i])^2)
  }
  
  ret <- sqrt((sumf/nums))
  
  return(ret)
}

PRatio <- function(pPredH, pTrainH){
  
  ret <- (pTrainH/pPredH)
  
  stopifnot(ret<1)
  
  return(ret)
}

NumP <- function(pTrainH){
  return(pTrainH)
}

SubAddi <- function(observed){
  nums <- length(observed)
  
  ivec <- c(1:nums)
  
  top1vec <- vector()
  
  for(i in 1:nums){
    top1vec[i] <- ((ivec[i]-1)*observed[i])
  }
  
  top1 <- sum(top1vec)
  
  top2 <- (((nums-1)/2)*sum(observed))
  
  top <- (top1-top2)
  
  return(top)
}

LapFact <- function(observed){
  nums <- length(observed)
  
  ivec <- c(1:nums)
  
  top1vec <- vector()

  for(i in 1:nums){
    top1vec[i] <- ((ivec[i]-1)*observed[i])
  }
  
  # top1vec <- ((ivec-1)*observed)
  
  top1 <- sum(top1vec)
  
  top2 <- (((nums-1)/2)*sum(observed))
  
  top <- (top1-top2)
  
  bot <- sqrt((((nums^2)-1)/12)*sum(observed))
  
  return((top/bot))
}
