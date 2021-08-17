#a: total failures, 
#g(t): fault detection rate with imperfection (b^2*t/1+bt)-(d/1+dt) 
#b: peak failure rate, d: degradation of failure rate 0 for no degradation

funcPZI <- function(ppar,tt){
  left <- ppar$a*exp(-ppar$b*tt)
  right <- 1+((ppar$b+ppar$d)*tt)+(ppar$b*ppar$d*(tt^2))
  return(ppar$a-left*right)
}

funcFRPZI <- function(ppar,tt){
  expr3 <- exp(-ppar$b * tt)
  expr4 <- ppar$a * expr3
  expr5 <- ppar$b + ppar$d
  expr8 <- ppar$b * ppar$d
  expr11 <- 1 + expr5 * tt + expr8 * tt^2
  return(-(expr4 * (expr5 + expr8 * (2 * tt)) - 
             ppar$a * (expr3 * ppar$b) * expr11))
}

residFuncPZI <- function(ppar,observed,tt){
  observed - funcPZI(ppar,tt)
}

parInitPZI <- function(){
  return(list(a=50,b=0.01,d=0.01))
}

parLowerPZI <- function(){
  return(c(0,0,0))
}

parUpperPZI <- function(){
  return(c(Inf,Inf,Inf))
}

parGridPZI <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=-8,to=7,length = 10)
  vecd <- 2^seq(from=-8,to=7,length = 10)
  df <- list(a=veca,b=vecb,d=vecd)
  return(df)
}

formulaPZI <- Mu ~ a-(a*exp(-b*t)*(1+((b+d)*t)+(b*d*(t^2))))

