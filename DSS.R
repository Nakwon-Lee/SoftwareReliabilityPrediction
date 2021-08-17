#a: total failures, b(t): ((b^2)*t)/(1+b*t) (b: peak failure rate)

funcDSS <- function(ppar,tt){
  return(ppar$a*(1-((1+(ppar$b*tt)*exp(-ppar$b*tt)))))
}

funcFRDSS <- function(ppar,tt){
  expr1 <- ppar$b * tt
  expr4 <- exp(-ppar$b * tt)
  return(-(ppar$a * (ppar$b * expr4 - expr1 * (expr4 * ppar$b))))
}

residFuncDSS <- function(ppar,observed,tt){
  observed - funcDSS(ppar,tt)
}

parInitDSS <- function(){
  return(list(a=100,b=0.1))
}

parLowerDSS <- function(){
  return(c(0,0))
}

parUpperDSS <- function(){
  return(c(Inf,Inf))
}

parGridDSS <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=-8,to=7,length = 10)
  df <- list(a=veca,b=vecb)
  return(df)
}

formulaDSS <- Mu ~ a*(1-((1+(b*t)*exp(-b*t))))