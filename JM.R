funcJM <- function(ppar,nn){
  return(ppar$k*(ppar$a-(nn-1)))
}

residFuncJM <- function(ppar,observed,nn){
  return(observed - funcJM(ppar,nn))
}

parInitJM <- function(){
  return(list(a=100,k=0.01))
}

parLowerJM<- function(){
  return(c(0,0))
}

parUpperJM <- function(){
  return(c(Inf,Inf))
}

parGridJM <- function(){
  veca <- 2^seq(from=0,to=14,length = 30)
  veck <- 2^seq(from=-8,to=7,length = 30)
  df <- list(a=veca,k=veck)
  return(df)
}

parmEstJM <- function(pparInit,observed){
  ivec <- c(1:length(observed))
  pcA <- length(observed)*sum(observed)
  pcB <- sum(ivec*observed)
  pcC <- sum(observed)
  
  # mle.out <- nlm(likelihood.a,c(pparInit$a),ivec,pcA,pcB,pcC)
  # estA <- mle.out$estimate
  # estK <- sum(1/(estA-ivec+1))/sum(observed)
  # return(list(a=estA,k=estK))
  
  estA <- bisectionMethodJM(ivec,pcA,pcB,pcC)

  if(estA != -1){
    topsum <- sum(1/(estA-(ivec-1)))
    estK <- topsum/sum(observed)
    return(list(a=estA,k=estK))
  }else{
    return(NULL)
  }
}

likelihood.a <- function(ppar,pivec,conA,conB,conC){
  abs((((ppar[1]*conC)-conB+conC)*sum(1/(ppar[1]-pivec+1)))-conA)
}

bisectionMethodJM <- function(pivec,conA,conB,conC){
  eps <- 0.0000000001
  from <- 0.001
  to <- sum(conC)/length(pivec)
  b1 <- from
  b2 <- to
  for(i in 1:10000000){
    b3 <- b1 + ((b2-b1)/2)
    
    val <- abs(likelihood.a(c(b3),pivec,conA,conB,conC))
    
    if((abs(val)<eps) || (((b2-b1)/2)<eps)){
      return(b3)
    }
    
    if((likelihood.a(c(b1),pivec,conA,conB,conC)*likelihood.a(c(b3),pivec,conA,conB,conC))<0){
      b2 <- b3
    }else{
      b1 <- b3
    }
  }
  return(-1)
}

formulaJM <- Mu ~ k*(a-(t-1))
