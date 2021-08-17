funcLV <- function(ppar,nn){
  return((ppar$a - 1)/((ppar$b * nn) + (ppar$c * (nn^2))))
}

residFuncLV <- function(ppar,observed,nn){
  return(observed - funcLV(ppar,nn))
}

parInitLV <- function(){
  return(list(a=10,b=10,c=5))
}

parLowerLV<- function(){
  return(c(0,0,0))
}

parUpperLV <- function(){
  return(c(Inf,Inf,Inf))
}

parGridLV <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=0,to=5,length = 10)
  vecc <- 2^seq(from=0,to=5,length = 10)
  df <- list(a=veca,b=vecb,c=vecc)
  return(df)
}

formulaLV <- Mu ~ (a - 1)/((b * t) + (c * (t^2)))