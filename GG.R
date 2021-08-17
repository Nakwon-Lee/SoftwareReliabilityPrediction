#a: total failures, b: faults per unit time, h: s-shaped over 1

funcGG<- function(ppar,tt){
  return(ppar$a*(1-exp(-ppar$b*(tt^ppar$h))))
}

residFuncGG <- function(ppar,observed,tt){
  observed - funcGG(ppar,tt)
}

funcFRGG<- function(ppar,tt){
  expr4 <- exp(-ppar$b * (tt^ppar$h))
  return(ppar$a * (expr4 * (ppar$b * (tt^(ppar$h-1) * ppar$h))))
}

parInitGG <- function(){
  return(list(a=50,b=0.01,h=1))
}

parLowerGG <- function(){
  return(c(0,0,0))
}

parUpperGG <- function(){
  return(c(Inf,Inf,Inf))
}

parGridGG <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=-8,to=7,length = 10)
  vech <- 2^seq(from=-8,to=7,length = 10)
  df <- list(a=veca,b=vecb,h=vech)
  return(df)
}

formulaGG <- Mu ~ a*(1-exp(-b*(t^h)))