# paramEstimation <- function(df,cfg,presidFunc,pformula){
#   tt <- c(1:length(df))
#   obsrvd <- df
#   
#   targetfn <- function(Init,pOb,tt){
#     sum(presidFunc(Init,pOb,tt)^2)
#   }
#   
#   lsm.out <- NULL
#   lsm.out <- nlm(targetfn,cfg$parInit,obsrvd,tt)
#   tryCatch(
#     {
#       lsm.out <- nlm(targetfn,cfg$parInit,obsrvd,tt)
#     },
#     error = function(cond){
#       #message(cond)
#     }
#   )
#   
#   parm <- NULL
#   if(!is.null(lsm.out)){
#     parm <- lsm.out$estimate
#   }
#   return(parm)
# }

paramEstimation <- function(df,cfg,presidFunc,pformula){
  t <- c(1:length(df))

  observeddata <- df

  obtemp <- list(observeddata)
  tt <- list(t)
  obdata <- list()
  obdata["Mu"] <- obtemp
  obdata["t"] <- tt

  numtech <- 2

  nls.out <- NULL

  if(cfg$algo=="LM"){
    nls.out <- nls.lm(par=cfg$parInit,lower = cfg$parLower,
                      upper = cfg$parUpper,
                      fn=presidFunc,observed=observeddata,tt=t,
                      control=nls.lm.control(maxiter=cfg$miter))
  }else if(cfg$algo=="GN"){
    nls.out <- nls(formula = pformula,data = obdata,start = cfg$parInit,
                   lower = cfg$parLower, upper = cfg$parUpper,
                   control = nls.control(maxiter = cfg$miter))
  }else if (cfg$algo=="GP"){
    nls.out <- nls(formula = pformula,data = obdata,start = cfg$parInit,
                   algorithm = "plinear",
                   lower = cfg$parLower, upper = cfg$parUpper,
                   control = nls.control(maxiter = cfg$miter))
  }else if (cfg$algo=="PO"){
    nls.out <- nls(formula = pformula,data = obdata,start = cfg$parInit,
                   algorithm = "port",
                   lower = cfg$parLower, upper = cfg$parUpper,
                   control = nls.control(maxiter = cfg$miter))
  }else if (cfg$algo=="Best"){

    nls.out.list <- list()

    nls.out.list[1] <- list(NULL)
    nls.out.list[2] <- list(NULL)
    
    tryCatch(
      {
        nls.out.list[[1]] <- nls.lm(par=cfg$parInit,
                                    fn=presidFunc,observed=observeddata,tt=t,
                                    control=nls.lm.control(maxiter=cfg$miter,maxfev = cfg$miter*length(cfg$parInit)))
      },
      error = function(cond){
        #message(cond)
      }
    )
    
    if(!is.null(nls.out.list[[1]])){
      if(nls.out.list[[1]]$info==0 || 
         nls.out.list[[1]]$info==5 || 
         nls.out.list[[1]]$info==9 || 
         nls.out.list[[1]]$info==6 ||
         nls.out.list[[1]]$info==7 ||
         nls.out.list[[1]]$info==8){
        nls.out.list[1] <- list(NULL)
      }
    }
    
    tryCatch(
      {
        nls.out.list[[2]] <- nls(formula = pformula,data = obdata,start = cfg$parInit,
                                 control = nls.control(maxiter = cfg$miter))
      },
      error = function(cond){
        #message(cond)
      }
    )
    
    if(!is.null(nls.out.list[[2]])){
      if(!(nls.out.list[[2]]$convInfo$isConv)){
        nls.out.list[2] <- list(NULL)
      }
    }
    

    nls.out.sse.vec <- vector()

    for(i in 1:numtech){

      residlist <- NULL

      if(is.null(nls.out.list[[i]])){
        residlist <- observeddata
      }else{
        residlist <- residuals(nls.out.list[[i]])
      }
      nls.out.sse.vec[i] <- sum(residlist^2)
    }
    
    nls.out <- nls.out.list[[which.min(nls.out.sse.vec)]]
  }
  
  # thisinit <- bruteForceNls(residFunc = presidFunc,pgrid = cfg$parGrid,obsrvd = observeddata,pt = t)

  #validation
  
  retout <- NULL
  
  if (!is.null(nls.out)) {
    retout <- as.list(coef(nls.out))
  }

  # if (is.null(nls.out)) {
  #   
  #   retout <- thisinit
  #   
  # }else{
  # 
  #   validity <- TRUE
  #   parlist <- as.list(coef(nls.out))
  # 
  #   for (i in 1:length(cfg$parInit)){
  #     if (parlist[[i]] < cfg$parLower[i]){
  #       validity <- FALSE
  #       break
  #     }
  #     if (parlist[[i]] > cfg$parUpper[i]){
  #       validity <- FALSE
  #       break
  #     }
  #   }
  # 
  #   if (!validity) {
  #     
  #     retout <- thisinit
  #     
  #   }else{
  #     retout <- as.list(coef(nls.out))
  #   }
  # }
  
  return(retout)
}


bruteForceNls <- function(residFunc,pgrid,obsrvd,pt){
  griddf <- expand.grid(pgrid)
  ssevec <- vector()
  for (i in 1:nrow(griddf)) {
    ssevec[i] <- sum(residFunc(griddf[i,],obsrvd,pt)^2)
  }
  return(griddf[which.min(ssevec),])
}

