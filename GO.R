#a: total failures, b: rate of change of failure rate (estimated per failure failure rate)

funcGO <- function(ppar,tt){
  return(ppar$a*(1-exp(-ppar$b*tt)))
}

residFuncGO <- function(ppar,observed,tt){
  return(observed - funcGO(ppar,tt))
}

funcFRGO <- function(ppar,tt){
  expr3 <- exp(-ppar$b*tt)
  return(ppar$a * (expr3 * ppar$b))
}

parInitGO <- function(){
  return(list(a=100,b=0.01))
}

parLowerGO <- function(){
  return(c(0,0))
}

parUpperGO <- function(){
  return(c(Inf,Inf))
}

parGridGO <- function(){
  veca <- 2^seq(from=0,to=14,length = 30)
  vecb <- 2^seq(from=-8,to=7,length = 30)
  df <- list(a=veca,b=vecb)
  return(df)
}

formulaGO <- Mu ~ a*(1-exp(-b*t))