funcGEO <- function(ppar,nn){
  return(ppar$b*(ppar$p^(nn-1)))
}

residFuncGEO <- function(ppar,observed,nn){
  return(observed - funcJM(ppar,nn))
}

parInitGEO <- function(){
  return(list(p=0.8,b=1))
}

parLowerGEO<- function(){
  return(c(0,0))
}

parUpperGEO <- function(){
  return(c(1,Inf))
}

parGridGEO <- function(){
  vecp <- 2^seq(from=-8,to=0,length = 30)
  vecb <- 2^seq(from=-5,to=5,length = 30)
  df <- list(p=vecp,b=vecb)
  return(df)
}

parmEstGEO <- function(pparInit,observed){
  counts <- c(1:length(observed))
  mle.out <- nlm(likelihood.p,c(pparInit$p),observed,counts,length(observed))
  estP <- mle.out$estimate
  estB <- (estP*length(observed))/sum((estP^counts)*observed)
  return(list(p=estP,b=estB))
  
  # estP <- bisectionMethodGEO(observed,counts,length(observed))
  # if(estP != -1){
  #   estB <- (estP*length(observed))/sum((estP^counts)*observed)
  #   return(list(p=estP,b=estB))
  # }else{
  #   return(NULL)
  # }
}

likelihood.p <- function(ppar,obsrvd,ivec,nlen){
  (((ppar[1]*nlen)/sum((ppar[1]^ivec)*obsrvd))*sum(ivec*(ppar[1]^ivec)*obsrvd))-sum(ivec/ppar[1])
}

bisectionMethodGEO <- function(ob,ivec,nlen){
  eps <- 0.0000000001
  from <- 0.001
  to <- sum(ob)/length(ivec)
  b1 <- from
  b2 <- to
  for(i in 1:10000000){
    b3 <- b1 + ((b2-b1)/2)
    
    val <- abs(likelihood.p(c(b3),ob,ivec,nlen))
    
    if((abs(val)<eps) || (((b2-b1)/2)<eps)){
      return(b3)
    }
    
    if((likelihood.p(c(b1),ob,ivec,nlen)*likelihood.p(c(b3),ob,ivec,nlen))<0){
      b2 <- b3
    }else{
      b1 <- b3
    }
  }
  return(-1)
}

formulaGEO <- Mu ~ b*(p^(t-1))