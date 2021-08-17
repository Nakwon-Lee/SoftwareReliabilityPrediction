#a(t) = c+a(1-exp(-e*t))  b(t) = b/(1+(d*exp(-b*t)))
#h: decreasing rate of the error introduction (might be significantly smaller than b)
#d: 0: fixed failure rate, 1: s-shaped rate start

funcPZ <- function(ppar,tt){
  left <- 1/(1+(ppar$d*exp(-ppar$b*tt)))
  mid <- (ppar$c+ppar$a)*(1-exp(-ppar$b*tt))
  right <- (ppar$a/(ppar$b-ppar$h))*(exp(-ppar$h*tt)-exp(-ppar$b*tt))
  return(left*(mid-right))
}

funcFRPZ <- function(ppar,tt){
  expr3 <- exp(-ppar$b * tt)
  expr5 <- 1 + ppar$d * expr3
  expr6 <- 1/expr5
  expr7 <- ppar$c + ppar$a
  expr11 <- ppar$a/(ppar$b - ppar$h)
  expr14 <- exp(-ppar$h * tt)
  expr17 <- expr7 * (1 - expr3) - expr11 * (expr14 - expr3)
  expr19 <- expr3 * ppar$b
  return(ppar$d * expr19/expr5^2 * expr17 + expr6 * 
           (expr7 * expr19 + expr11 * (expr14 * ppar$h - expr19)))
}

residFuncPZ <- function(ppar,observed,tt){
  observed - funcPZ(ppar,tt)
}

parInitPZ <- function(){
  return(list(a=50,b=0.1,c=1,d=1,h=0.01))
}

parLowerPZ <- function(){
  return(c(0,0,0,0,0))
}

parUpperPZ <- function(){
  return(c(Inf,Inf,Inf,Inf,Inf))
}

parGridPZ <- function(){
  veca <- 2^seq(from=0,to=14,length = 7)
  vecb <- 2^seq(from=-8,to=7,length = 7)
  vecc <- 2^seq(from=0,to=10,length = 3)
  vecd <- 2^seq(from=-3,to=3,length = 3)
  vech <- 2^seq(from=-5,to=5,length = 3)
  df <- list(a=veca,b=vecb,c=vecc,d=vecd,h=vech)
  return(df)
}

formulaPZ <- Mu ~ (1/(1+(d*exp(-b*t))))*(((c+a)*(1-exp(-b*t)))-((a/(b-h))*(exp(-h*t)-exp(-b*t))))


