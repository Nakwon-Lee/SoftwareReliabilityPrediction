# a(t) = a*exp(c*t) ( a:a(0), c:failure introducing rate ), b: rate of change of failure rate (estimated per failure failure rate)

funcYID1 <- function(ppar,tt){
  return(((ppar$a*ppar$b)/(ppar$b+ppar$c))*(exp(ppar$c*tt)-exp(ppar$b*tt)))
}

funcFRYID1 <- function(ppar,tt){
  expr3 <- ppar$a * (ppar$b/(ppar$b+ppar$c))
  expr5 <- exp(ppar$c*tt)
  expr7 <- exp(ppar$b*tt)
  return(expr3*((expr5*ppar$c) - (expr7*ppar$b)))
}

residFuncYID1 <- function(ppar,observed,tt){
  observed - funcYID1(ppar,tt)
}

parInitYID1 <- function(){
  # numrow <- nrow(df)
  # observed <- df$n[numrow]
  # ntvec <- vector()
  # for (i in 1:numrow){
  #   ntvec[i] <- df$n[i]/(i*8)
  # }
  # lm.out <- lm(df$n ~ ntvec)
  # initb <- abs(1/as.list(coef(lm.out))$ntvec)
  # print(list(a=observed,b=initb,c=0))
  return(list(a=50,b=0.1,c=0))
}

parLowerYID1 <- function(){
  return(c(0,0,0))
}

parUpperYID1 <- function(){
  return(c(Inf,Inf,Inf))
}

parGridYID1 <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=-8,to=7,length = 10)
  vecc <- 2^seq(from=-5,to=5,length = 10)
  df <- list(a=veca,b=vecb,c=vecc)
  return(df)
}

formulaYID1 <- Mu ~ ((a*b)/(b+c))*(exp(c*t)-exp(b*t))

