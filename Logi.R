#a: total failures, b: failure rate, d: s-shaped when bigger than 1 (almost same a for valid estimation at t=0)

funcLogi <- function(ppar,tt){
  bot <- 1+(ppar$d*exp(-ppar$b*tt))
  return(ppar$a/bot)
}

funcFRLogi <- function(ppar,tt){
  expr3 <- exp(-ppar$b * tt)
  expr5 <- 1 + ppar$d * expr3
  return(ppar$a * (ppar$d * (expr3 * ppar$b))/expr5^2)
}

residFuncLogi <- function(ppar,observed,tt){
  observed - funcLogi(ppar,tt)
}

parInitLogi <- function(){
  return(list(a=50,b=0.01,d=45))
}

parLowerLogi <- function(){
  return(c(0,0,0))
}

parUpperLogi <- function(){
  return(c(Inf,Inf,Inf))
}

parGridLogi <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=-8,to=7,length = 10)
  vecd <- 2^seq(from=-5,to=5,length = 10)
  df <- list(a=veca,b=vecb,d=vecd)
  return(df)
}

formulaLogi <- Mu ~ a/(1+(d*exp(-b*t)))