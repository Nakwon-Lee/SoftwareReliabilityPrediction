#| Pert | MetaKnowledge ~ |

getMETAInfo <- function(df,
                        metafullvec,
                        totpert,
                        isFCTBF){
  colnameleft <- c("Pert","Model")
  
  numrows <- nrow(df)
  
  numrowssubs <- vector()
  
  if(totpert!=0){
    for(i in 1:totpert){
      pert <- 0.1*((10-totpert-1)+i)
      tnumrows <- floor(numrows * pert)
      numrowssubs[i] <- tnumrows
    }
  }
  
  observed <- NULL
  elapobsvd <- NULL
  
  if(isFCTBF=="FC"){
    observed <- df$f
    elapobsvd <- df$n
  }else if(isFCTBF=="TBF"){
    observed <- df$f
    elapobsvd <- vector()
    for(i in 1:length(df$f)){
      elapobsvd[i] <- sum(df$f[1:i])
    }
  }
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  
  leftmatrix <- matrix(nrow=totpert,ncol=length(colnameleft))
  rightmatrix <- matrix(nrow=totpert,ncol=length(metafullvec))
  
  for(i in 1:totpert){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    
    currrowidx <- i
    
    observedsub <- observed[1:numrowssubs[i]]
    elapobsvdsub <- elapobsvd[1:numrowssubs[i]]
    
    leftmatrix[currrowidx,1] <- currpert
    leftmatrix[currrowidx,2] <- "None"
    rightmatrix[currrowidx,1] <- Variance(observedsub)
    rightmatrix[currrowidx,2] <- Inclination(c(1:numrowssubs[i]),observedsub)
    rightmatrix[currrowidx,3] <- AutoCorr(c(1:numrowssubs[i]),observedsub)
    rightmatrix[currrowidx,4] <- MetaNO(observedsub)
    rightmatrix[currrowidx,5] <- PRatio(numrows,numrowssubs[[i]])
    rightmatrix[currrowidx,6] <- NumP(numrowssubs[[i]])
    rightmatrix[currrowidx,7] <- LapFact(observedsub)
    rightmatrix[currrowidx,8] <- SubAddi(observedsub)
    rightmatrix[currrowidx,9] <- Variance(elapobsvdsub)
    rightmatrix[currrowidx,10] <- Inclination(c(1:numrowssubs[i]),elapobsvdsub)
    rightmatrix[currrowidx,11] <- AutoCorr(c(1:numrowssubs[i]),elapobsvdsub)
    rightmatrix[currrowidx,12] <- MetaNO(elapobsvdsub)
    rightmatrix[currrowidx,13] <- LapFact(elapobsvdsub)
    rightmatrix[currrowidx,14] <- SubAddi(elapobsvdsub)
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnameleft
  
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- metafullvec
  
  retinstance <- cbind(leftdf,rightdf)
  
  return(retinstance)
}

getMETAInfoAddi <- function(df,
                            metafullvec,
                            totpert,
                            isFCTBF){
  colnameleft <- c("Pert","Model")
  
  numrows <- nrow(df)
  
  pertlist <- list()
  numrowslist <- list()
  
  for(i in 1:(totpert-1)){
    trainH <- NULL
    tempstart <- NULL
    tempstart <- (((10-totpert)-1)+i)
    trainH <- floor((tempstart*0.1)*numrows)
    for(j in i:(totpert-1)){
      predH <- NULL
      tempend <- NULL
      tempend <- ((10-totpert)+j)
      predH <- floor((tempend*0.1)*numrows)
      
      pertlist[[(length(pertlist)+1)]] <- paste0("P",as.character(tempstart),as.character(tempend))
      temphorizon <- c(trainH,predH)
      stopifnot(trainH<predH)
      numrowslist[[(length(numrowslist)+1)]] <- temphorizon
    }
  }
  
  stopifnot(length(pertlist)==(((totpert-1)*totpert)/2))
  stopifnot(length(pertlist)==length(numrowslist))
  
  pertlistaddi <- list()
  numrowslistaddi <- list()
  
  if(numrows>=20){
    for (i in 1:totpert) {
      trainH <- NULL
      tempstart <- NULL
      tempstart <- ((((10-totpert)-1)+i)*10)+5
      trainH <- floor((tempstart/100)*numrows)
      
      predH <- NULL
      predH <- numrows
      
      pertlistaddi[[length(pertlistaddi)+1]] <- paste0("T",as.character(tempstart))
      temphorizon <- c(trainH,predH)
      stopifnot(trainH<predH)
      numrowslistaddi[[(length(numrowslistaddi)+1)]] <- temphorizon
    }
    
    for(i in 1:(totpert-1)){
      trainH <- NULL
      tempstart <- NULL
      tempstart <- ((((10-totpert)-1)+i)*10)+5
      trainH <- floor((tempstart/100)*numrows)
      for(j in i:(totpert-1)){
        predH <- NULL
        tempend <- NULL
        tempend <- (((10-totpert)+j)*10)+5
        predH <- floor((tempend/100)*numrows)
        
        pertlistaddi[[(length(pertlistaddi)+1)]] <- paste0("T",as.character(tempstart),as.character(tempend))
        temphorizon <- c(trainH,predH)
        stopifnot(trainH<predH)
        numrowslistaddi[[(length(numrowslistaddi)+1)]] <- temphorizon
      }
    }
  }
  
  if(length(pertlistaddi)>0){
    pertlist <- c(pertlist,pertlistaddi)
    numrowslist <- c(numrowslist,numrowslistaddi)
  }
  
  observed <- NULL
  elapobsvd <- NULL
  
  if(isFCTBF=="FC"){
    observed <- df$f
    elapobsvd <- df$n
  }else if(isFCTBF=="TBF"){
    observed <- df$f
    elapobsvd <- vector()
    for(i in 1:length(df$f)){
      elapobsvd[i] <- sum(df$f[1:i])
    }
  }
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  
  leftmatrix <- matrix(nrow=length(pertlist),ncol=length(colnameleft))
  rightmatrix <- matrix(nrow=length(pertlist),ncol=length(metafullvec))
  
  for(i in 1:length(pertlist)){
    currpert <- pertlist[[i]]
    
    currrowidx <- i
    
    observedsub <- observed[1:numrowslist[[i]][1]]
    elapobsvdsub <- elapobsvd[1:numrowslist[[i]][1]]
    
    leftmatrix[currrowidx,1] <- currpert
    leftmatrix[currrowidx,2] <- "None"
    rightmatrix[currrowidx,1] <- Variance(observedsub)
    rightmatrix[currrowidx,2] <- Inclination(c(1:numrowslist[[i]][1]),observedsub)
    rightmatrix[currrowidx,3] <- AutoCorr(c(1:numrowslist[[i]][1]),observedsub)
    rightmatrix[currrowidx,4] <- MetaNO(observedsub)
    rightmatrix[currrowidx,5] <- PRatio(numrowslist[[i]][2],numrowslist[[i]][1])
    rightmatrix[currrowidx,6] <- NumP(numrowslist[[i]][1])
    rightmatrix[currrowidx,7] <- LapFact(observedsub)
    rightmatrix[currrowidx,8] <- SubAddi(observedsub)
    rightmatrix[currrowidx,9] <- Variance(elapobsvdsub)
    rightmatrix[currrowidx,10] <- Inclination(c(1:numrowslist[[i]][1]),elapobsvdsub)
    rightmatrix[currrowidx,11] <- AutoCorr(c(1:numrowslist[[i]][1]),elapobsvdsub)
    rightmatrix[currrowidx,12] <- MetaNO(elapobsvdsub)
    rightmatrix[currrowidx,13] <- LapFact(elapobsvdsub)
    rightmatrix[currrowidx,14] <- SubAddi(elapobsvdsub)
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnameleft
  
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- metafullvec
  
  retinstance <- cbind(leftdf,rightdf)
  
  return(retinstance)
}

getMETADPInfo <- function(df,crivec,totpert){
  
  nump <- length(crivec)
  
  colnameleft <- c("Pert","Model")
  
  numrows <- nrow(df)
  
  numrowssubs <- list()
  
  if(totpert!=0){
    for(i in 1:totpert){
      pert <- 0.1*((10-totpert-1)+i)
      tnumrows <- floor(numrows * pert)
      numrowssubs[[i]] <- tnumrows
    }
  }
  
  dfsubs <- list()
  
  if(totpert!=0){
    for(i in 1:totpert){
      dfsubs[[i]] <- df[1:numrowssubs[[i]],]
    }
  }
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  
  leftmatrix <- matrix(nrow=totpert,ncol=2)
  rightmatrix <- matrix(nrow=totpert,ncol=(length(crivec)+1))
  
  for(i in 1:totpert){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    
    currrowidx <- i
    
    # min-max normalization
    dps <- dfsubs[[i]]$n
    
    normdps <- (dps-min(dps))/(max(dps)-min(dps))
    
    pointsI <- seq(from=1,to=numrowssubs[[i]],length=nump)
    
    leftmatrix[currrowidx,1] <- currpert
    leftmatrix[currrowidx,2] <- currmodel
    
    rightmatrix[currrowidx,1] <- 0
    for(j in 2:(length(pointsI)-1)){
      tfloor <- floor(pointsI[j])
      tceiling <- ceiling(pointsI[j])
      
      if((tceiling-tfloor)==0){
        rightmatrix[currrowidx,j] <- normdps[pointsI[j]]
      }else{
        tslope <- (normdps[tceiling]-normdps[tfloor])/(tceiling-tfloor)
        tyintercept <- (normdps[tceiling]-(tslope*tceiling))
        
        rightmatrix[currrowidx,j] <- (tslope*pointsI[j])+tyintercept
      }
    }
    rightmatrix[currrowidx,length(pointsI)] <- 1
    rightmatrix[currrowidx,(length(pointsI)+1)] <- PRatio(numrows,numrowssubs[[i]])
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnameleft
  
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- c(crivec,"PRatio")
  
  retinstance <- cbind(leftdf,rightdf)
  
  return(retinstance)
}

getMETADPInfoAddi <- function(df,totpert,crivec){
  
  nump <- length(crivec)
  
  colnameleft <- c("Pert","Model")
  
  currmodel <- "RecDDM"
  
  numrows <- nrow(df)
  
  pertlist <- list()
  numrowslist <- list()
  
  for(i in 1:(totpert-1)){
    trainH <- NULL
    tempstart <- NULL
    tempstart <- (((10-totpert)-1)+i)
    trainH <- floor((tempstart*0.1)*numrows)
    for(j in i:(totpert-1)){
      predH <- NULL
      tempend <- NULL
      tempend <- ((10-totpert)+j)
      predH <- floor((tempend*0.1)*numrows)
      
      pertlist[[(length(pertlist)+1)]] <- paste0("P",as.character(tempstart),as.character(tempend))
      temphorizon <- c(trainH,predH)
      stopifnot(trainH<predH)
      numrowslist[[(length(numrowslist)+1)]] <- temphorizon
    }
  }
  
  stopifnot(length(pertlist)==(((totpert-1)*totpert)/2))
  stopifnot(length(pertlist)==length(numrowslist))
  
  pertlistaddi <- list()
  numrowslistaddi <- list()
  
  if(numrows>=20){
    for (i in 1:totpert) {
      trainH <- NULL
      tempstart <- NULL
      tempstart <- ((((10-totpert)-1)+i)*10)+5
      trainH <- floor((tempstart/100)*numrows)
      
      predH <- NULL
      predH <- numrows
      
      pertlistaddi[[length(pertlistaddi)+1]] <- paste0("T",as.character(tempstart))
      temphorizon <- c(trainH,predH)
      stopifnot(trainH<predH)
      numrowslistaddi[[(length(numrowslistaddi)+1)]] <- temphorizon
    }
    
    for(i in 1:(totpert-1)){
      trainH <- NULL
      tempstart <- NULL
      tempstart <- ((((10-totpert)-1)+i)*10)+5
      trainH <- floor((tempstart/100)*numrows)
      for(j in i:(totpert-1)){
        predH <- NULL
        tempend <- NULL
        tempend <- (((10-totpert)+j)*10)+5
        predH <- floor((tempend/100)*numrows)
        
        pertlistaddi[[(length(pertlistaddi)+1)]] <- paste0("T",as.character(tempstart),as.character(tempend))
        temphorizon <- c(trainH,predH)
        stopifnot(trainH<predH)
        numrowslistaddi[[(length(numrowslistaddi)+1)]] <- temphorizon
      }
    }
  }
  
  if(length(pertlistaddi)>0){
    pertlist <- c(pertlist,pertlistaddi)
    numrowslist <- c(numrowslist,numrowslistaddi)
  }
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  
  leftmatrix <- matrix(nrow=length(pertlist),ncol=2)
  rightmatrix <- matrix(nrow=length(pertlist),ncol=(length(crivec)+1))
  
  for(i in 1:length(pertlist)){
    currpert <- pertlist[[i]]
    
    currrowidx <- i
    
    dfsub <- df[1:numrowslist[[i]][1],]
    
    # min-max normalization
    dps <- dfsub$n
    normdps <- (dps-min(dps))/(max(dps)-min(dps))
    
    pointsI <- seq(from=1,to=numrowslist[[i]][1],length=nump)
    
    leftmatrix[currrowidx,1] <- currpert
    leftmatrix[currrowidx,2] <- currmodel
    
    rightmatrix[currrowidx,1] <- 0
    for(j in 2:(length(pointsI)-1)){
      tfloor <- floor(pointsI[j])
      tceiling <- ceiling(pointsI[j])
      
      if((tceiling-tfloor)==0){
        rightmatrix[currrowidx,j] <- normdps[pointsI[j]]
      }else{
        tslope <- (normdps[tceiling]-normdps[tfloor])/(tceiling-tfloor)
        tyintercept <- (normdps[tceiling]-(tslope*tceiling))
        
        rightmatrix[currrowidx,j] <- (tslope*pointsI[j])+tyintercept
      }
    }
    rightmatrix[currrowidx,length(pointsI)] <- 1
    rightmatrix[currrowidx,(length(pointsI)+1)] <- PRatio(numrowslist[[i]][2],numrowslist[[i]][1])
  }
    
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnameleft
  
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- c(crivec,"PRatio")
  
  retinstance <- cbind(leftdf,rightdf)
  
  return(retinstance)
}

getMETAInfoOrig <- function(df,
                        metafts,
                        ptrainh){
  
  observed <- df$f
  elapobsvd <- df$n
  
  leftmatrix <- matrix(nrow=1,ncol=1)
  rightmatrix <- matrix(nrow=1,ncol=length(metafts))
  
  observedsub <- observed[1:ptrainh]
  elapobsvdsub <- elapobsvd[1:ptrainh]
  
  leftmatrix[1,1] <- 'None'
  
  rightmatrix[1,1] <- Variance(observedsub)
  rightmatrix[1,2] <- Inclination(c(1:ptrainh),observedsub)
  rightmatrix[1,3] <- AutoCorr(c(1:ptrainh),observedsub)
  rightmatrix[1,4] <- MetaNO(observedsub)
  rightmatrix[1,5] <- NumP(ptrainh)
  rightmatrix[1,6] <- LapFact(observedsub)
  rightmatrix[1,7] <- SubAddi(observedsub)
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- c('Model')
  
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- metafts
  
  retinstance <- cbind(leftdf,rightdf)
  
  return(retinstance)
}