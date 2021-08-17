#a: total failures, f: time for 50% detection, g: weight (~1~)

funcMD <- function(ppar,tt){
  return(ppar$a*(1-((ppar$f/(ppar$f+tt))^ppar$g)))
}

funcFRMD <- function(ppar,tt){
  expr1 <- ppar$f+tt
  expr2 <- ppar$f/expr1
  return(ppar$a*((expr2^(ppar$g-1))*(ppar$g*(ppar$f/(expr1^2)))))
}

residFuncMD <- function(ppar,observed,tt){
  observed - funcMD(ppar,tt)
}

parInitMD <- function(){
  return(list(a=100,f=50,g=1))
}

parLowerMD <- function(){
  return(c(0,0,0))
}

parUpperMD <- function(){
  return(c(Inf,Inf,Inf))
}

parGridMD <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecf <- 2^seq(from=0,to=10,length = 10)
  vecg <- 2^seq(from=-5,to=5,length = 10)
  df <- list(a=veca,f=vecf,g=vecg)
  return(df)
}

formulaMD <- Mu ~ a*(1-((f/(f+t))^g))