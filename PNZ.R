#a(t) = a(1+ct), b(t) = b/(1+d*exp(-bt))
#a: total failures, b: failure rate
#c: rate of introduced faults per detected faults
#d: 0=fixed failure rate and 1=s-shaped rate start

funcPNZ <- function(ppar,tt){
  left <- ppar$a/(1+(ppar$d*exp(-ppar$b*tt)))
  right <- ((1-exp(-ppar$b*tt))*(1-(ppar$c/ppar$b)))+(ppar$c*ppar$a*tt)
  return(left*right)
}

funcFRPNZ <- function(ppar,tt){
  expr3 <- exp(-ppar$b * tt)
  expr5 <- 1 + ppar$d * expr3
  expr6 <- ppar$a/expr5
  expr9 <- 1 - ppar$c/ppar$b
  expr11 <- ppar$c * ppar$a
  expr13 <- (1 - expr3) * expr9 + expr11 * tt
  expr15 <- expr3 * ppar$b
  return(ppar$a * (ppar$d * expr15)/expr5^2 * expr13 + 
           expr6 * (expr15 * expr9 + expr11))
}

residFuncPNZ <- function(ppar,observed,tt){
  observed - funcPNZ(ppar,tt)
}

parInitPNZ <- function(){
  return(list(a=50,b=0.1,c=0.001,d=1))
}

parLowerPNZ <- function(){
  return(c(0,0,0,0))
}

parUpperPNZ <- function(){
  return(c(Inf,Inf,Inf,Inf))
}

parGridPNZ <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=-8,to=7,length = 10)
  vecc <- 2^seq(from=-5,to=5,length = 3)
  vecd <- 2^seq(from=-5,to=5,length = 3)
  df <- list(a=veca,b=vecb,c=vecc,d=vecd)
  return(df)
}

formulaPNZ <- Mu ~ (a/(1+(d*exp(-b*t))))*(((1-exp(-b*t))*(1-(c/b)))+(c*a*t))


