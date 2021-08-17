funcMB <- function(ppar,nn){
  return(ppar$b*(1-(nn/ppar$a)))
}

residFuncMB <- function(ppar,observed,nn){
  return(observed - funcMB(ppar,nn))
}

parInitMB <- function(){
  return(list(a=150,b=0.01))
}

parLowerMB<- function(){
  return(c(0,0))
}

parUpperMB <- function(){
  return(c(Inf,Inf))
}

parGridMB <- function(){
  veca <- 2^seq(from=0,to=14,length = 30)
  vecb <- 2^seq(from=-8,to=7,length = 30)
  df <- list(a=veca,b=vecb)
  return(df)
}

formulaMB <- Mu ~ b * (1 - (t/a))