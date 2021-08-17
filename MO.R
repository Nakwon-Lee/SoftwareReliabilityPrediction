#p*q: initial failure intensity p^-1 dubbed failure intensity decay param 

funcMO <- function(ppar,tt){
  return(ppar$p*log(1+(ppar$q*tt)))
}

funcFRMO <- function(ppar,tt){
  expr2 <- 1 + (ppar$q * tt)
  return(ppar$p*(ppar$q/expr2))
}

residFuncMO <- function(ppar,observed,tt){
  observed - funcMO(ppar,tt)
}

parInitMO <- function(){
  return(list(p=100,q=0.0001))
}

parLowerMO <- function(){
  return(c(0,0))
}

parUpperMO <- function(){
  return(c(Inf,Inf))
}

parGridMO <- function(){
  vecp <- 2^seq(from=0,to=14,length = 30)
  vecq <- 2^seq(from=-3,to=3,length = 30)
  df <- list(p=vecp,q=vecq)
  return(df)
}

formulaMO <- Mu ~ p*log(1+(q*t))