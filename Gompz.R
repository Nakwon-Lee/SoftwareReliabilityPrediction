#a: total failures, j: 1-failure rate, k: s-shaped parameter with a quite small value

funcGompz <- function(ppar,tt){
  return(ppar$a*(ppar$k^(ppar$j^(tt))))
}

funcFRGompz <- function(ppar,tt){
  expr1 <- ppar$j^tt
  expr2 <- ppar$k^expr1
  return(ppar$a * (expr2*(log(ppar$k) * (expr1 * log(ppar$j)))))
}

residFuncGompz <- function(ppar,observed,tt){
  observed - funcGompz(ppar,tt)
}

parInitGompz <- function(){
  return(list(a=100,j=0.9,k=0.0001))
}

parLowerGompz <- function(){
  return(c(0,0,0))
}

parUpperGompz <- function(){
  return(c(Inf,1,1))
}

parGridGompz <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecj <- 2^seq(from=-8,to=0,length = 10)
  veck <- 2^seq(from=-8,to=0,length = 10)
  df <- list(a=veca,j=vecj,k=veck)
  return(df)
}

formulaGompz <- Mu ~ a*(k^(j^t))