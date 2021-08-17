getDDMEstforSingle <- function(df,
                               models,
                               isFCTBF){
  numrows <- nrow(df)
  
  retlist <- list()
  
  dflearn <- NULL
  observed <- NULL
  
  if(isFCTBF=="FC"){
    dflearn <- df$n
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    dflearn <- vector()
    for(i in 1:length(tobsvd)){
      dflearn[i] <- sum(tobsvd[1:i])
    }
    observed <- vector()
    for(i in 1:length(tobsvd)){
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  trainH <- numrows
  
  retlist$trainH <- trainH
  
  estimated <- NULL
  fdim <- NULL
  
  dflearnsub <- dflearn[1:numrows]
  
  print(dflearnsub)
  
  for(j in 1:length(models)){
    if(models[j]=="SVR"){
      
      maxval <- max(dflearnsub)
      minval <- min(dflearnsub)
      
      denominator <- maxval-minval
      
      stopifnot(denominator>0)
      
      nmzdfsub <- 0.1+((0.9-0.1)*((dflearnsub-minval)/denominator))

      retddm <- RecDDSRM(nmzdfsub,floor(0.9*numrows),numrows,floor(numrows*(10/7)))
      
      estimated <- retddm$est
      fdim <- retddm$fdim
      
      estimated <- (((estimated-0.1)*denominator)/0.8)+minval
      
    }else if(models[j]=="ANN"){
      
      nmzdfsub <- (2*((dflearnsub-min(dflearnsub))/(max(dflearnsub)-min(dflearnsub))))-1
      
      retddm <- DDMRecANN(nmzdfsub,trainH,numrows)
      
      estimated <- retddm$est
      fdim <- retddm$fdim
      
      estimated <- (((estimated+1)*(max(dflearnsub)-min(dflearnsub)))/2)+min(dflearnsub)
      
    }else{
      stopifnot(FALSE)
    }
    
    print(estimated)
    
    stopifnot(length(estimated)==(floor(numrows*(10/7))-fdim))
    
    estimated <- ifelse(estimated<0,0,estimated)
    
    stopifnot(estimated>=0)
    
    estimated <- c(dflearn[1:fdim],estimated)
    
    frestimated <- NULL
    if(isFCTBF=="FC"){
      frestimated <- calcFCFR(estimated)
    }else if(isFCTBF=="TBF"){
      frestimated <- calcTBFFR(estimated)
    }
    
    stopifnot(length(estimated)==floor(numrows*(10/7)))
    
    retlist[[models[j]]] <- list()
    
    retlist[[models[j]]]$dim <- fdim
    retlist[[models[j]]]$EstElap <- estimated
    retlist[[models[j]]]$EstFr <- frestimated
  }
  
  return(retlist)
}

getDDMEstOrig <- function(df,model,trainh,predh){
  
  retlist <- list()
  retlist$trainH <- trainh
  
  dflearn <- df$n
  
  estimated <- NULL
  fdim <- NULL
  
  dflearnsub <- dflearn[1:trainh]
  
  if(model=="SVR"){
    
    maxval <- max(dflearnsub)
    minval <- min(dflearnsub)
    
    denominator <- maxval-minval
    
    stopifnot(denominator>0)
    
    nmzdfsub <- 0.1+((0.9-0.1)*((dflearnsub-minval)/denominator))
    
    retddm <- RecDDSRM(nmzdfsub,floor(0.9*trainh),trainh,predh)
    
    estimated <- retddm$est
    fdim <- retddm$fdim
    
    estimated <- (((estimated-0.1)*denominator)/0.8)+minval
    
  }else if(model=="ANN"){
    
    nmzdfsub <- (2*((dflearnsub-min(dflearnsub))/(max(dflearnsub)-min(dflearnsub))))-1
    
    retddm <- DDMRecANN(nmzdfsub,trainh,predh)
    
    estimated <- retddm$est
    fdim <- retddm$fdim
    
    estimated <- (((estimated+1)*(max(dflearnsub)-min(dflearnsub)))/2)+min(dflearnsub)
    
  }else{
    stopifnot(FALSE)
  }
  
  frestimated <- calcFCFR(estimated)
  
  stopifnot(length(estimated)==(predh-fdim))
  estimated <- ifelse(estimated<0,0,estimated)
  stopifnot(estimated>=0)
  estimated <- c(dflearn[1:fdim],estimated)
  
  retlist[[model]] <- list()
  retlist[[model]]$dim <- fdim
  retlist[[model]]$EstElap <- estimated
  retlist[[model]]$EstFr <- frestimated
  
  return(retlist)
}

getDDMEst <- function(df,
                      models,
                      totpert,
                      isFCTBF){
  
  numrows <- nrow(df)
  
  numrowssubs <- list()
  
  retlist <- list()
  
  dflearn <- NULL
  observed <- NULL
  
  if(isFCTBF=="FC"){
    dflearn <- df$n
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    dflearn <- vector()
    for(i in 1:length(tobsvd)){
      dflearn[i] <- sum(tobsvd[1:i])
    }
    observed <- vector()
    for(i in 1:length(tobsvd)){
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  for(i in 1:totpert){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    tempstart <- (((10-totpert)-1)+i)
    trainH <- floor((tempstart*0.1)*numrows)
    numrowssubs[[(length(numrowssubs)+1)]] <- trainH
    
    retlist[[currpert]] <- list()
    retlist[[currpert]]$trainH <- trainH
    
    estimated <- NULL
    fdim <- NULL
    
    dflearnsub <- dflearn[1:trainH]
    
    for(j in 1:length(models)){
      if(models[j]=="SVR"){
        
        maxval <- max(dflearnsub)
        minval <- min(dflearnsub)
        
        denominator <- maxval-minval
        
        stopifnot(denominator>0)
        
        nmzdfsub <- 0.1+((0.9-0.1)*((dflearnsub-minval)/denominator))
        
        if (i==1){
          retddm <- RecDDSRM(nmzdfsub,floor(numrows*(3/10)),numrowssubs[[i]],numrows)
        }else{
          retddm <- RecDDSRM(nmzdfsub,numrowssubs[[(i-1)]],numrowssubs[[i]],numrows)
        }
        
        estimated <- retddm$est
        fdim <- retddm$fdim
        
        estimated <- (((estimated-0.1)*denominator)/0.8)+minval
        
      }else if(models[j]=="ANN"){
        
        nmzdfsub <- (2*((dflearnsub-min(dflearnsub))/(max(dflearnsub)-min(dflearnsub))))-1
        
        retddm <- DDMRecANN(nmzdfsub,trainH,numrows)
        
        estimated <- retddm$est
        fdim <- retddm$fdim
        
        estimated <- (((estimated+1)*(max(dflearnsub)-min(dflearnsub)))/2)+min(dflearnsub)
        
      }else{
        stopifnot(FALSE)
      }
      
      stopifnot(length(estimated)==(numrows-fdim))
      
      # estimated <- ifelse(estimated<0,0,estimated)
      # 
      # stopifnot(estimated>=0)
      
      estimated <- c(dflearn[1:fdim],estimated)
      
      frestimated <- NULL
      if(isFCTBF=="FC"){
        frestimated <- calcFCFR(estimated)
      }else if(isFCTBF=="TBF"){
        frestimated <- calcTBFFR(estimated)
      }
      
      stopifnot(length(estimated)==numrows)
      
      retlist[[currpert]][[models[j]]] <- list()
      
      retlist[[currpert]][[models[j]]]$dim <- fdim
      retlist[[currpert]][[models[j]]]$EstElap <- estimated
      retlist[[currpert]][[models[j]]]$EstFr <- frestimated
    }
    
  }
  
  return(retlist)
}

calcFCFR <- function(pEstmtd){
  frest <- vector()
  frest[1] <- pEstmtd[1]
  for(i in 2:(length(pEstmtd)-1)){
    frest[i] <- ((pEstmtd[(i+1)] - pEstmtd[i]) + (pEstmtd[i] - pEstmtd[(i-1)]))/2 
  }
  frest[length(pEstmtd)] <- pEstmtd[length(pEstmtd)] - pEstmtd[(length(pEstmtd)-1)]
  
  ret <- ifelse(frest<0,0,frest)
  
  return(ret)
}

calcTBFFR <- function(pEstmtd){
  frest <- vector()
  frest[1] <- pEstmtd[1]
  for (i in 2:length(pEstmtd)) {
    frest[i] <- pEstmtd[i]-pEstmtd[(i-1)]
  }
  
  ret <- ifelse(frest<0,0,frest)
  
  return(ret)
}

getDDMEstfive <- function(df,
                          models,
                          totpert,
                          isFCTBF){
  numrows <- nrow(df)
  
  numrowssubs <- list()
  
  retlist <- list()
  
  dflearn <- NULL
  observed <- NULL
  
  if(isFCTBF=="FC"){
    dflearn <- df$n
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    dflearn <- vector()
    for(i in 1:length(tobsvd)){
      dflearn[i] <- sum(tobsvd[1:i])
    }
    observed <- vector()
    for(i in 1:length(tobsvd)){
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  if(numrows>=20){
    for(i in 1:totpert){
      currpert <- paste0("T",as.character((((10-totpert-1)+i)*10)+5))
      tempstart <- ((((10-totpert)-1)+i)*10)+5
      trainH <- floor((tempstart/100)*numrows)
      numrowssubs[[(length(numrowssubs)+1)]] <- trainH
      
      retlist[[currpert]] <- list()
      retlist[[currpert]]$trainH <- trainH
      
      dflearnsub <- dflearn[1:trainH]
      
      estimated <- NULL
      fdim <- NULL
      
      for(j in 1:length(models)){
        if(models[j]=="SVR"){
          
          maxval <- max(dflearnsub)
          minval <- min(dflearnsub)
          
          denominator <- maxval-minval
          
          stopifnot(denominator>0)
          
          nmzdfsub <- 0.1+((0.9-0.1)*((dflearnsub-minval)/denominator))
          
          if (i==1){
            retddm <- RecDDSRM(nmzdfsub,floor(numrows*(35/100)),numrowssubs[[i]],numrows)
          }else{
            retddm <- RecDDSRM(nmzdfsub,numrowssubs[[(i-1)]],numrowssubs[[i]],numrows)
          }
          
          estimated <- retddm$est
          fdim <- retddm$fdim
          
          estimated <- (((estimated-0.1)*denominator)/0.8)+minval
          
        }else if(models[j]=="ANN"){
          
          nmzdfsub <- (2*((dflearnsub-min(dflearnsub))/(max(dflearnsub)-min(dflearnsub))))-1
          
          retddm <- DDMRecANN(nmzdfsub,trainH,numrows)
          
          estimated <- retddm$est
          fdim <- retddm$fdim
          
          estimated <- (((estimated+1)*(max(dflearnsub)-min(dflearnsub)))/2)+min(dflearnsub)
          
        }else{
          stopifnot(FALSE)
        }
        
        stopifnot(length(estimated)==(numrows-fdim))
        
        # estimated <- ifelse(estimated<0,0,estimated)
        # 
        # stopifnot(estimated>=0)
        
        estimated <- c(dflearn[1:fdim],estimated)
        
        frestimated <- NULL
        if(isFCTBF=="FC"){
          frestimated <- calcFCFR(estimated)
        }else if(isFCTBF=="TBF"){
          frestimated <- calcTBFFR(estimated)
        }
        
        stopifnot(length(estimated)==numrows)
        
        retlist[[currpert]][[models[j]]] <- list()
        
        retlist[[currpert]][[models[j]]]$dim <- fdim
        retlist[[currpert]][[models[j]]]$EstElap <- estimated
        retlist[[currpert]][[models[j]]]$EstFr <- frestimated
      }
    }
  }
  
  return(retlist)
}


getDDMGOF <- function(pEst,
                      pEstfive,
                      df,
                      modelvec,
                      goffullvec,
                      crivec,
                      totpert,
                      isFCTBF){
  
  colnameleft <- c("Model","Pert")
  
  observed <- NULL
  
  if(isFCTBF=="FC"){
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    observed <- vector()
    for (i in 1:length(tobsvd)) {
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  leftmatrix <- matrix(nrow=(totpert*length(modelvec)),ncol=length(colnameleft))
  rightmatrix <- matrix(nrow=(totpert*length(modelvec)),ncol=(length(goffullvec)-1))
  lastmatrix <- matrix(nrow = (totpert*length(modelvec)),ncol = 1)
  
  for(i in 1:totpert){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    
    thisddm <- pEst[[currpert]]
    
    trainH <- thisddm$trainH
    
    for(j in 1:length(modelvec)){
      
      estimated <- NULL
      frestimated <- NULL
      dim <- NULL
      
      thismodelddm <- thisddm[[modelvec[j]]]
      estimated <- thismodelddm$EstElap
      estimated <- ifelse(estimated<0,0,estimated)
      frestimated <- thismodelddm$EstFr
      frestimated <- ifelse(frestimated<0,0,frestimated)
      dim <- thismodelddm$dim
      
      currmodel <- modelvec[j]
      
      currrowidx  <- j+((i-1)*length(modelvec))
      
      fittrainH <- (trainH-dim)
      fitobserved <- observed[-c(1:dim)]
      fitestimated <- estimated[-c(1:dim)]

      fitfrestimated <- NULL
      fitfrestimated <- frestimated[-c(1:dim)]
      
      leftmatrix[currrowidx,1] <- currmodel
      leftmatrix[currrowidx,2] <- currpert
      if(isFCTBF=="FC"){
        rightmatrix[currrowidx,1] <- FC.MSE(fittrainH,0,fitobserved,fitestimated)
        rightmatrix[currrowidx,2] <- FC.MAE(fittrainH,0,fitobserved,fitestimated)
        rightmatrix[currrowidx,3] <- FC.Rsquare(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,4] <- FC.Noise(fittrainH,fitfrestimated)
        rightmatrix[currrowidx,5] <- FC.Bias(1,fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,6] <- FC.Variation(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,7] <- FC.PRR(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,8] <- FC.WLSE(fittrainH,fitobserved,fitestimated,0.2)
        rightmatrix[currrowidx,9] <- FC.EP(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,10] <- FC.MEOP(1,fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,11] <- FC.TOUP(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,12] <- PRatio(length(estimated),trainH)
        rightmatrix[currrowidx,13] <- FC.EP(length(estimated),observed,estimated)
        rightmatrix[currrowidx,14] <- FC.MEOP(trainH,length(estimated),observed,estimated)
        rightmatrix[currrowidx,15] <- abs(FC.Bias(1,fittrainH,fitobserved,fitestimated))
        lastmatrix[currrowidx,1] <- GROUP(FC.Bias(trainH,length(estimated),observed,estimated))
      }else if(isFCTBF=="TBF"){
        rightmatrix[currrowidx,1] <- TBF.MSE(fittrainH,0,fitobserved,fitestimated)
        rightmatrix[currrowidx,2] <- TBF.MAE(fittrainH,0,fitobserved,fitestimated)
        rightmatrix[currrowidx,3] <- TBF.Rsquare(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,4] <- TBF.Noise(fittrainH,fitfrestimated)
        rightmatrix[currrowidx,5] <- TBF.Bias(1,fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,6] <- TBF.Variation(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,7] <- TBF.PRR(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,8] <- TBF.WLSE(fittrainH,fitobserved,fitestimated,0.2)
        rightmatrix[currrowidx,9] <- TBF.PE(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,10] <- TBF.AE(1,fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,11] <- TBF.TOUP(fittrainH,fitobserved,fitestimated)
        rightmatrix[currrowidx,12] <- PRatio(length(estimated),trainH)
        rightmatrix[currrowidx,13] <- TBF.PE(length(estimated),observed,estimated)
        rightmatrix[currrowidx,14] <- TBF.AE(trainH,length(estimated),observed,estimated)
        rightmatrix[currrowidx,15] <- abs(TBF.Bias(1,fittrainH,fitobserved,fitestimated))
        lastmatrix[currrowidx,1] <- GROUP(TBF.Bias(trainH,length(estimated),observed,estimated))
      }
    }
  }
  
  leftmatrix2 <- NULL
  rightmatrix2 <- NULL
  addrows <- ((totpert-1)*totpert)/2
  leftmatrix2 <- matrix(nrow=(addrows*length(modelvec)),ncol=length(colnameleft))
  rightmatrix2 <- matrix(nrow=(addrows*length(modelvec)),ncol=(length(goffullvec)-1))
  lastmatrix2 <- matrix(nrow=(addrows*length(modelvec)),ncol=1)
  
  tempidx <- 1
  
  for(i in 1:(totpert-1)){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    thisddm <- pEst[[currpert]]
    tempstart <- NULL
    tempstart <- (((10-totpert)-1)+i)
    trainH <- thisddm$trainH
    
    for(j in i:(totpert-1)){
      
      predH <- NULL
      tempend <- NULL
      tempend <- ((10-totpert)+j)
      predH <- floor((tempend*0.1)*length(observed))
      
      addpert <- paste0("P",as.character(tempstart),as.character(tempend))
      
      for(k in 1:length(modelvec)){
        
        estimated <- NULL
        frestimated <- NULL
        dim <- NULL
        
        thismodelddm <- thisddm[[modelvec[k]]]
        estimated <- thismodelddm$EstElap
        estimated <- ifelse(estimated<0,0,estimated)
        frestimated <- thismodelddm$EstFr
        frestimated <- ifelse(frestimated<0,0,frestimated)
        dim <- thismodelddm$dim
        
        currmodel <- modelvec[k]
        
        fitpredH <- predH-dim
        fittrainH <- (trainH-dim)
        fitobserved <- observed[-c(1:dim)]
        fitobserved <- fitobserved[1:fitpredH]
        fitestimated <- estimated[-c(1:dim)]
        fitestimated <- fitestimated[1:fitpredH]
        
        fitfrestimated <- NULL
        fitfrestimated <- frestimated[-c(1:dim)]
        fitfrestimated <- fitfrestimated[1:fitpredH]
        
        testimated <- estimated[1:predH]
        observedsub <- observed[1:predH]
        
        currrowidx <- tempidx
        
        leftmatrix2[currrowidx,1] <- currmodel
        leftmatrix2[currrowidx,2] <- addpert
        
        if(isFCTBF=="FC"){
          rightmatrix2[currrowidx,1] <- FC.MSE(fittrainH,0,fitobserved,fitestimated)
          rightmatrix2[currrowidx,2] <- FC.MAE(fittrainH,0,fitobserved,fitestimated)
          rightmatrix2[currrowidx,3] <- FC.Rsquare(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,4] <- FC.Noise(fittrainH,fitfrestimated)
          rightmatrix2[currrowidx,5] <- FC.Bias(1,fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,6] <- FC.Variation(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,7] <- FC.PRR(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,8] <- FC.WLSE(fittrainH,fitobserved,fitestimated,0.2)
          rightmatrix2[currrowidx,9] <- FC.EP(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,10] <- FC.MEOP(1,fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,11] <- FC.TOUP(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,12] <- PRatio(predH,trainH)
          rightmatrix2[currrowidx,13] <- FC.EP(predH,observedsub,testimated)
          rightmatrix2[currrowidx,14] <- FC.MEOP(trainH,predH,observedsub,testimated)
          rightmatrix2[currrowidx,15] <- abs(FC.Bias(1,fittrainH,fitobserved,fitestimated))
          lastmatrix2[currrowidx,1] <- GROUP(FC.Bias(trainH,predH,observedsub,testimated))
        }else if(isFCTBF=="TBF"){
          rightmatrix2[currrowidx,1] <- TBF.MSE(fittrainH,0,fitobserved,fitestimated)
          rightmatrix2[currrowidx,2] <- TBF.MAE(fittrainH,0,fitobserved,fitestimated)
          rightmatrix2[currrowidx,3] <- TBF.Rsquare(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,4] <- TBF.Noise(fittrainH,fitfrestimated)
          rightmatrix2[currrowidx,5] <- TBF.Bias(1,fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,6] <- TBF.Variation(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,7] <- TBF.PRR(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,8] <- TBF.WLSE(fittrainH,fitobserved,fitestimated,0.2)
          rightmatrix2[currrowidx,9] <- TBF.PE(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,10] <- TBF.AE(1,fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,11] <- TBF.TOUP(fittrainH,fitobserved,fitestimated)
          rightmatrix2[currrowidx,12] <- PRatio(predH,trainH)
          rightmatrix2[currrowidx,13] <- TBF.PE(predH,observedsub,testimated)
          rightmatrix2[currrowidx,14] <- TBF.AE(trainH,predH,observedsub,testimated)
          rightmatrix2[currrowidx,15] <- abs(TBF.Bias(1,fittrainH,fitobserved,fitestimated))
          lastmatrix2[currrowidx,1] <- GROUP(TBF.Bias(trainH,predH,observedsub,testimated))
        }
        tempidx <- tempidx + 1
      }
    }
  }
  
  leftmatrix3 <- NULL
  rightmatrix3 <- NULL
  lastmatrix3 <- NULL
  
  if(nrow(df)>=20){
    addrows <- (((totpert-1)*totpert)/2)+totpert
    leftmatrix3 <- matrix(nrow=(addrows*length(modelvec)),ncol=2)
    rightmatrix3 <- matrix(nrow=(addrows*length(modelvec)),ncol=(length(goffullvec)-1))
    lastmatrix3 <- matrix(nrow = (addrows*length(modelvec)),ncol = 1)
    
    tempidx <- 1
    
    for (i in 1:totpert) {
      currpert <- paste0("T",as.character((((10-totpert-1)+i)*10)+5))
      
      thisddm <- pEstfive[[currpert]]
      trainH <- thisddm$trainH
      
      for(j in 1:length(modelvec)){
        
        estimated <- NULL
        frestimated <- NULL
        dim <- NULL
        
        thismodelddm <- thisddm[[modelvec[j]]]
        estimated <- thismodelddm$EstElap
        estimated <- ifelse(estimated<0,0,estimated)
        frestimated <- thismodelddm$EstFr
        frestimated <- ifelse(frestimated<0,0,frestimated)
        dim <- thismodelddm$dim
        
        currmodel <- modelvec[j]
        
        currrowidx  <- tempidx
        
        fittrainH <- (trainH-dim)
        fitobserved <- observed[-c(1:dim)]
        fitestimated <- estimated[-c(1:dim)]
        
        fitfrestimated <- NULL
        fitfrestimated <- frestimated[-c(1:dim)]
        
        leftmatrix3[currrowidx,1] <- currmodel
        leftmatrix3[currrowidx,2] <- currpert
        if(isFCTBF=="FC"){
          rightmatrix3[currrowidx,1] <- FC.MSE(fittrainH,0,fitobserved,fitestimated)
          rightmatrix3[currrowidx,2] <- FC.MAE(fittrainH,0,fitobserved,fitestimated)
          rightmatrix3[currrowidx,3] <- FC.Rsquare(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,4] <- FC.Noise(fittrainH,fitfrestimated)
          rightmatrix3[currrowidx,5] <- FC.Bias(1,fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,6] <- FC.Variation(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,7] <- FC.PRR(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,8] <- FC.WLSE(fittrainH,fitobserved,fitestimated,0.2)
          rightmatrix3[currrowidx,9] <- FC.EP(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,10] <- FC.MEOP(1,fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,11] <- FC.TOUP(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,12] <- PRatio(length(estimated),trainH)
          rightmatrix3[currrowidx,13] <- FC.EP(length(estimated),observed,estimated)
          rightmatrix3[currrowidx,14] <- FC.MEOP(trainH,length(estimated),observed,estimated)
          rightmatrix3[currrowidx,15] <- abs(FC.Bias(1,fittrainH,fitobserved,fitestimated))
          lastmatrix3[currrowidx,1] <- GROUP(FC.Bias(trainH,length(estimated),observed,estimated))
        }else if(isFCTBF=="TBF"){
          rightmatrix3[currrowidx,1] <- TBF.MSE(fittrainH,0,fitobserved,fitestimated)
          rightmatrix3[currrowidx,2] <- TBF.MAE(fittrainH,0,fitobserved,fitestimated)
          rightmatrix3[currrowidx,3] <- TBF.Rsquare(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,4] <- TBF.Noise(fittrainH,fitfrestimated)
          rightmatrix3[currrowidx,5] <- TBF.Bias(1,fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,6] <- TBF.Variation(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,7] <- TBF.PRR(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,8] <- TBF.WLSE(fittrainH,fitobserved,fitestimated,0.2)
          rightmatrix3[currrowidx,9] <- TBF.PE(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,10] <- TBF.AE(1,fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,11] <- TBF.TOUP(fittrainH,fitobserved,fitestimated)
          rightmatrix3[currrowidx,12] <- PRatio(length(estimated),trainH)
          rightmatrix3[currrowidx,13] <- TBF.PE(length(estimated),observed,estimated)
          rightmatrix3[currrowidx,14] <- TBF.AE(trainH,length(estimated),observed,estimated)
          rightmatrix3[currrowidx,15] <- abs(TBF.Bias(1,fittrainH,fitobserved,fitestimated))
          lastmatrix3[currrowidx,1] <- GROUP(TBF.Bias(trainH,length(estimated),observed,estimated))
        }
        tempidx <- tempidx + 1
      }
    }
    
    for(i in 1:(totpert-1)){
      currpert <- paste0("T",as.character((((10-totpert-1)+i)*10)+5))
      thisddm <- pEstfive[[currpert]]
      tempstart <- NULL
      tempstart <- (((10-totpert-1)+i)*10)+5
      trainH <- thisddm$trainH
      
      for(j in i:(totpert-1)){
        predH <- NULL
        tempend <- NULL
        tempend <- (((10-totpert)+j)*10)+5
        predH <- floor((tempend/100)*length(estimated))
        
        addpert <- paste0("T",as.character(tempstart),as.character(tempend))
        
        for(k in 1:length(modelvec)){
          
          estimated <- NULL
          frestimated <- NULL
          dim <- NULL
          
          thismodelddm <- thisddm[[modelvec[k]]]
          estimated <- thismodelddm$EstElap
          estimated <- ifelse(estimated<0,0,estimated)
          frestimated <- thismodelddm$EstFr
          frestimated <- ifelse(frestimated<0,0,frestimated)
          dim <- thismodelddm$dim
          
          currmodel <- modelvec[k]
          
          fitpredH <- predH-dim
          fittrainH <- (trainH-dim)
          fitobserved <- observed[-c(1:dim)]
          fitobserved <- fitobserved[1:fitpredH]
          fitestimated <- estimated[-c(1:dim)]
          fitestimated <- fitestimated[1:fitpredH]
          
          fitfrestimated <- NULL
          fitfrestimated <- frestimated[-c(1:dim)]
          fitfrestimated <- fitfrestimated[1:fitpredH]
          
          testimated <- estimated[1:predH]
          observedsub <- observed[1:predH]
          
          currrowidx <- tempidx
          
          leftmatrix3[currrowidx,1] <- currmodel
          leftmatrix3[currrowidx,2] <- addpert
          
          if(isFCTBF=="FC"){
            rightmatrix3[currrowidx,1] <- FC.MSE(fittrainH,0,fitobserved,fitestimated)
            rightmatrix3[currrowidx,2] <- FC.MAE(fittrainH,0,fitobserved,fitestimated)
            rightmatrix3[currrowidx,3] <- FC.Rsquare(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,4] <- FC.Noise(fittrainH,fitfrestimated)
            rightmatrix3[currrowidx,5] <- FC.Bias(1,fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,6] <- FC.Variation(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,7] <- FC.PRR(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,8] <- FC.WLSE(fittrainH,fitobserved,fitestimated,0.2)
            rightmatrix3[currrowidx,9] <- FC.EP(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,10] <- FC.MEOP(1,fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,11] <- FC.TOUP(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,12] <- PRatio(predH,trainH)
            rightmatrix3[currrowidx,13] <- FC.EP(predH,observedsub,testimated)
            rightmatrix3[currrowidx,14] <- FC.MEOP(trainH,predH,observedsub,testimated)
            rightmatrix3[currrowidx,15] <- abs(FC.Bias(1,fittrainH,fitobserved,fitestimated))
            lastmatrix3[currrowidx,1] <- GROUP(FC.Bias(trainH,predH,observedsub,testimated))
          }else if(isFCTBF=="TBF"){
            rightmatrix3[currrowidx,1] <- TBF.MSE(fittrainH,0,fitobserved,fitestimated)
            rightmatrix3[currrowidx,2] <- TBF.MAE(fittrainH,0,fitobserved,fitestimated)
            rightmatrix3[currrowidx,3] <- TBF.Rsquare(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,4] <- TBF.Noise(fittrainH,fitfrestimated)
            rightmatrix3[currrowidx,5] <- TBF.Bias(1,fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,6] <- TBF.Variation(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,7] <- TBF.PRR(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,8] <- TBF.WLSE(fittrainH,fitobserved,fitestimated,0.2)
            rightmatrix3[currrowidx,9] <- TBF.PE(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,10] <- TBF.AE(1,fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,11] <- TBF.TOUP(fittrainH,fitobserved,fitestimated)
            rightmatrix3[currrowidx,12] <- PRatio(predH,trainH)
            rightmatrix3[currrowidx,13] <- TBF.PE(predH,observedsub,testimated)
            rightmatrix3[currrowidx,14] <- TBF.AE(trainH,predH,observedsub,testimated)
            rightmatrix3[currrowidx,15] <- abs(TBF.Bias(1,fittrainH,fitobserved,fitestimated))
            lastmatrix3[currrowidx,1] <- GROUP(TBF.Bias(trainH,predH,observedsub,testimated))
          }
          tempidx <- tempidx + 1
        }
      }
    }
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnameleft
  rightdf <- data.frame(rightmatrix)
  lastdf <- data.frame(lastmatrix)
  
  rightdf <- cbind(rightdf,lastdf)
  names(rightdf) <- goffullvec
  
  leftdf2 <- data.frame(leftmatrix2)
  names(leftdf2) <- colnameleft
  
  rightdf2 <- data.frame(rightmatrix2)
  lastdf2 <- data.frame(lastmatrix2)
  rightdf2 <- cbind(rightdf2,lastdf2)
  names(rightdf2) <- goffullvec
  
  tempaddddm <- NULL
  
  if(!is.null(leftmatrix3)){
    leftdf3 <- data.frame(leftmatrix3)
    names(leftdf3) <- colnameleft
    
    rightdf3 <- data.frame(rightmatrix3)
    lastdf3 <- data.frame(lastmatrix3)
    rightdf3 <- cbind(rightdf3,lastdf3)
    names(rightdf3) <- goffullvec
    
    tempaddddm <- na.omit(rbind(cbind(leftdf2,rightdf2),cbind(leftdf3,rightdf3)))
  }else{
    tempaddddm <- na.omit(cbind(leftdf2,rightdf2))
  }
  
  retinstance <- list()
  
  retinstance$ddm <- na.omit(cbind(leftdf,rightdf))
  
  retinstance$ddm$Group <- as.factor(retinstance$ddm$Group)
  
  retinstance$addddm <- tempaddddm
  
  retinstance$addddm$Group <- as.factor(retinstance$addddm$Group)
  
  return(retinstance)
}

getDDMGOFforaData <- function(pEst,df,
                              ptrainh,ppredh,
                              pmodel,pgoffts){
  
  leftmatrix <- matrix(nrow=1,ncol=1)
  rightmatrix <- matrix(nrow=1,ncol=length(pgoffts))
  lastmatrix <- matrix(nrow=1,ncol=2)
  
  observed <- df$n
  
  currmodel <- pmodel
  
  thismodelddm <- pEst[[currmodel]]
  estimated <- thismodelddm$EstElap
  estimated <- ifelse(estimated<0,0,estimated)
  frestimated <- thismodelddm$EstFr
  frestimated <- ifelse(frestimated<0,0,frestimated)
  dim <- thismodelddm$dim
  
  fittrainH <- (ptrainh-dim)
  fitpredH <- (ppredh-dim)
  fitobserved <- observed[-c(1:dim)]
  fitestimated <- estimated[-c(1:dim)]
  fitfrestimated <- frestimated[-c(1:dim)]
  
  leftmatrix[1,1] <- currmodel
    
  rightmatrix[1,1] <- FC.MSE(fittrainH,0,fitobserved,fitestimated)
  rightmatrix[1,2] <- FC.MAE(fittrainH,0,fitobserved,fitestimated)
  rightmatrix[1,3] <- FC.Rsquare(fittrainH,fitobserved,fitestimated)
  rightmatrix[1,4] <- FC.Noise(fittrainH,fitfrestimated)
  rightmatrix[1,5] <- FC.Bias(1,fittrainH,fitobserved,fitestimated)
  rightmatrix[1,6] <- FC.Variation(fittrainH,fitobserved,fitestimated)
  rightmatrix[1,7] <- FC.PRR(fittrainH,fitobserved,fitestimated)
  rightmatrix[1,8] <- FC.WLSE(fittrainH,fitobserved,fitestimated,0.2)
  rightmatrix[1,9] <- FC.EP(fittrainH,fitobserved,fitestimated)
  rightmatrix[1,10] <- FC.MEOP(1,fittrainH,fitobserved,fitestimated)
  
  lastmatrix[1,1] <- FC.EP(ppredh,observed,estimated)
  lastmatrix[1,2] <- FC.MEOP(ptrainh,ppredh,observed,estimated)
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- c('Model')
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- pgoffts
  lastdf <- data.frame(lastmatrix)
  names(lastdf) <- c('EP','MEOP')
  
  return(cbind(leftdf,rightdf,lastdf))
}