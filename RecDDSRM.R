RecDDSRM <- function(df,validationH,trainH,predictH){
  
  bestMSE <- Inf
  
  # grideps <- c((2^-5),(2^-4),(2^-3),(2^-2),(2^-1),
  #              (2^0),(2^1),(2^2),(2^3),(2^4),(2^5))
  
  grideps <- c((2^-5),(2^-4),(2^-3),(2^-2),(2^-1))
  gridcost <- c((2^-1),(2^0),(2^1),(2^2),(2^3),(2^4),
                (2^5),(2^6),(2^7),(2^8),(2^9),(2^10))
  gridK <- c((2^-16),(2^-15),(2^-14),(2^-13),(2^-12),
             (2^-11),(2^-10),(2^-9),(2^-8),(2^-7),(2^-6),
             (2^-5),(2^-4),(2^-3),(2^-2),(2^-1),(2^0))
  
  predname <- "P"
  
  fidxeps <- 1
  fidxcost <- 1
  fidxK <- 1
  
  fdim <- 1
  
  for(pd in 1:5){
    
    d <- 1
    
    if(validationH-pd < 2){
      d <- validationH-2
    }else{
      d <- pd
    }
    
    innamevec <- vector()
    for(i in 1:d){
      innamevec[i] <- paste0("D",as.character(i))
    }
    
    tmat <- matrix(nrow = (validationH-d),ncol = (d+1))
    
    for(i in 1:(validationH-d)){
      for(j in 1:d){
        tmat[i,j] <- df[(i-1+j)]
      }
      tmat[i,(d+1)] <- df[(i+d)]
    }
    
    valframe <- data.frame(tmat)
    
    colnames(valframe) <- c(innamevec,predname)
    
    tvars <- innamevec[1]
    if (d!=1){
      for(i in 2:length(innamevec)){
        tvars <- paste0(tvars,"+",innamevec[i])
      }
    }
    
    tformula <- as.formula(paste0("P ~ ",tvars))
    
    idxeps <- 1
    idxcost <- 1
    idxK <- 1
    
    locbestMSE <- Inf
    
    validationdf <- valframe[c(1:(validationH-d)),]
    
    for(i in 1:length(grideps)){
      for(j in 1:length(gridcost)){
        for(k in 1:length(gridK)){
          
          #print(paste0(grideps[i]," ",gridcost[j]," ",gridK[k]))
          
          tmodel <- NULL
          
          tryCatch(
            {
              tmodel <- svm(tformula,data = validationdf,
                            type = "eps-regression",kernel = "radial", 
                            epsilon = grideps[i], cost = gridcost[j],
                            gamma = gridK[k])
            },
            error = function(cond){
              #message(cond)
            }
          )
          
          if(!is.null(tmodel)){
            ttf <- predictionDDM(tmodel,validationdf,
                                 validationH,trainH,d,innamevec,predname)
            
            currMSE <- MSE2(validationH-d,trainH-d,df[c((d+1):trainH)],ttf[,predname])
            
            if(currMSE < locbestMSE){
              idxeps <- i
              idxcost <- j
              idxK <- k
              
              locbestMSE <- currMSE
            }
          }
          
        }
      }
    }
    
    if(locbestMSE < bestMSE){
      fdim <- d
      fidxeps <- idxeps
      fidxcost <- idxcost
      fidxK <- idxK
      
      bestMSE <- locbestMSE
    }
    
  }
  
  
  
  innvec <- vector()
  for(i in 1:fdim){
    innvec[i] <- paste0("D",as.character(i))
  }
  
  fmat <- matrix(nrow = (trainH-fdim),ncol = (fdim+1))
  
  for(i in 1:(trainH - fdim)){
    for(j in 1:fdim){
      fmat[i,j] <- df[(i-1+j)]
    }
    fmat[i,(fdim+1)] <- df[(i+fdim)]
  }
  
  trainframe <- data.frame(fmat)
  
  colnames(trainframe) <- c(innvec,predname)
  
  fvars <- innvec[1]
  if(fdim!=1){
    for(i in 2:length(innvec)){
      fvars <- paste0(fvars,"+",innvec[i])
    }
  }
  
  fformula <- as.formula(paste0("P ~ ",fvars))
  
  #eps, cost, K(gamma)
  parsvm <- list()
  
  parsvm$eps <- grideps[fidxeps]
  parsvm$cost <- gridcost[fidxcost]
  parsvm$K <- gridK[fidxK]
  
  fddmodel <- svm(fformula,data = trainframe,
                 type = "eps-regression",kernel = "radial", 
                 epsilon = parsvm$eps, cost = parsvm$cost, gamma = parsvm$K)
  
  # print(fdim)
  # print(paste0(trainH," ",predictH))
  # print(trainframe)
  
  trainframe <- predictionDDM(fddmodel,trainframe,
                              trainH,predictH,fdim,innvec,predname)
  
  # return estimated, fdim
  retlist <- list()
  retlist$est <- trainframe[c(1:(predictH-fdim)),predname]
  retlist$fdim <- fdim
  
  return(retlist)
}

predictionDDM <- function(model,ptdf,tH,pH,pd,innvec,predn){
  
  tdf <- ptdf
  
  tlist <- list()
  if(pd > 1){
    for(i in 2:length(innvec)){
      tlist[innvec[(i-1)]] <- tdf[(tH-pd),innvec[i]]
    }
  }
  tlist[innvec[pd]] <- tdf[(tH-pd),predn]
  tlist[predn] <- 0
  
  # print(pd)
  # print(innvec)
  # print(tlist)
  
  xdf <- rbind(tdf,tlist)
  
  for (i in 1:nrow(tdf)) {
    temp <- predict(model,tdf[i,])
    tdf[i,predn] <- temp[1]
  }
  
  for(i in 1:(pH-tH)){
    temp <- predict(model,xdf[(tH-pd+i),])
    xdf[(tH-pd+i),predn] <- temp[1]
    
    ttlist <- list()
    if(pd > 1){
      for(j in 2:length(innvec)){
        ttlist[innvec[(j-1)]] <- xdf[(tH-pd+i),innvec[j]]
      }
    }
    ttlist[innvec[pd]] <- xdf[(tH-pd+i),predn]
    ttlist[predn] <- 0
    
    xdf <- rbind(xdf,ttlist)
  }
  
  xdf <- xdf[-(pH-pd+1),]
  xdf <- xdf[-c(1:nrow(tdf)),]
  
  xdf <- rbind(tdf,xdf)
  
  return(xdf)
}