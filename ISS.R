#a: total failures, b(t)=b/(1+d*exp(-bt)) (b: failure rate, d: 0 = fixed failure rate and 1 = s-shaped rate start)

funcISS <- function(ppar,tt){
  return(ppar$a*((1-exp(-ppar$b*tt))/(1+(((1-ppar$d)/ppar$d)*exp(-ppar$b*tt)))))
}

funcFRISS <- function(ppar,tt){
  return((ppar$a*ppar$b*(ppar$d+1)*exp(ppar$b*tt))/((exp(ppar$b*tt)+ppar$d)^2))
}

residFuncISS <- function(ppar,observed,tt){
  observed - funcISS(ppar,tt)
}

parInitISS <- function(){
  return(list(a=100,b=0.1,d=1))
}

parLowerISS <- function(){
  return(c(0,0,0))
}

parUpperISS <- function(){
  return(c(Inf,Inf,Inf))
}

parGridISS <- function(){
  veca <- 2^seq(from=0,to=14,length = 10)
  vecb <- 2^seq(from=-8,to=7,length = 10)
  vecd <- 2^seq(from=-5,to=5,length = 10)
  df <- list(a=veca,b=vecb,d=vecd)
  return(df)
}

formulaISS <- Mu ~ a*((1-exp(-b*t))/(1+(((1-d)/d)*exp(-b*t))))