#a(t): a(1+ct) (a:total failures, c:failure inducing rate), b: failure rate

funcYID2 <- function(ppar,tt){
  return((ppar$a*(1-exp(-ppar$b*tt))*(1-(ppar$c/ppar$b)))+(ppar$c*ppar$a*tt))
}

funcFRYID2 <- function(ppar,tt){
  expr3 <- exp(-ppar$b * tt)
  expr7 <- 1 - ppar$c/ppar$b
  expr9 <- ppar$c * ppar$a
  return(ppar$a * (expr3*ppar$b) *expr7 + expr9)
}

residFuncYID2 <- function(ppar,observed,tt){
  observed - funcYID2(ppar,tt)
}

parInitYID2 <- function(df){
  
  return(list(a=50,b=0.1,c=0))
}

parLowerYID2 <- function(){
  return(c(0,0,0))
}

parUpperYID2 <- function(){
  return(c(Inf,Inf,Inf))
}

parGridYID2 <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=-8,to=7,length = 10)
  vecc <- 2^seq(from=-5,to=5,length = 10)
  df <- list(a=veca,b=vecb,c=vecc)
  return(df)
}

formulaYID2 <- Mu ~ (a*(1-exp(-b*t))*(1-(c/b)))+(c*a*t)