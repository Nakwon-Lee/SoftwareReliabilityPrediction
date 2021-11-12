getDataPoints <- function(df,
                          totpert,
                          nump,isFCTBF){
  
  numrows <- nrow(df)
  
  numrowssubs <- vector()
  
  for(i in 1:totpert){
    pert <- 0.1*((10-totpert-1)+i)
    tnumrows <- floor(numrows * pert)
    numrowssubs[i] <- tnumrows
  }
  
  observed <- NULL
  
  if(isFCTBF=="FC"){
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- df$f
    observed <- vector()
    for(i in 1:length(tobsvd)){
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  retlist <- list()
  
  for(i in 1:totpert){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    
    retlist[[currpert]] <- list()
    
    observedsub <- observed[1:numrowssubs[i]]
    
    # min-max normalization
    
    normdps <- (observedsub-min(observedsub))/(max(observedsub)-min(observedsub))
    
    retlist[[currpert]]$raw <- normdps
    
    pointsI <- seq(from=1,to=numrowssubs[i],length=nump)
    
    tdpvec <- vector()
    
    tdpvec[1] <- 0
    for(j in 2:(length(pointsI)-1)){
      tfloor <- floor(pointsI[j])
      tceiling <- ceiling(pointsI[j])
      
      if((tceiling-tfloor)==0){
        tdpvec[j] <- normdps[pointsI[j]]
      }else{
        tslope <- (normdps[tceiling]-normdps[tfloor])/(tceiling-tfloor)
        tyintercept <- (normdps[tceiling]-(tslope*tceiling))
        
        tdpvec[j] <- (tslope*pointsI[j])+tyintercept
      }
    }
    tdpvec[length(pointsI)] <- 1
    
    retlist[[currpert]]$norm <- tdpvec
  }
  
  return(retlist)
}

getDataPointsSingle <- function(df,cpert,nump){
  
  numrows <- floor(nrow(df) * getPercentage(cpert))
  
  observed <- df$n[1:numrows]
  
  retlist <- list()
  
  normdps <- (observed-min(observed))/(max(observed)-min(observed))
  
  retlist$raw <- normdps
  
  pointsI <- seq(from=1,to=numrows,length=nump)
  
  tdpvec <- vector()
  
  tdpvec[1] <- 0
  for(j in 2:(length(pointsI)-1)){
    tfloor <- floor(pointsI[j])
    tceiling <- ceiling(pointsI[j])
    
    if((tceiling-tfloor)==0){
      tdpvec[j] <- normdps[pointsI[j]]
    }else{
      tslope <- (normdps[tceiling]-normdps[tfloor])/(tceiling-tfloor)
      tyintercept <- (normdps[tceiling]-(tslope*tceiling))
      
      tdpvec[j] <- (tslope*pointsI[j])+tyintercept
    }
  }
  tdpvec[length(pointsI)] <- 1
  
  stopifnot(!any(is.na(tdpvec)))
  
  retlist$norm <- tdpvec
  
  return(retlist)
}

getDataPointsaddi <- function(df,
                              totpert,
                              nump,isFCTBF){
  
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
  
  if(isFCTBF=="FC"){
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- df$f
    observed <- vector()
    for(i in 1:length(tobsvd)){
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  retlist <- list()
  
  for(i in 1:length(pertlist)){
    currpert <- pertlist[[i]]
    
    retlist[[currpert]] <- list()
    
    observedsub <- observed[1:numrowslist[[i]][1]]
    
    normdps <- (observedsub-min(observedsub))/(max(observedsub)-min(observedsub))
    
    retlist[[currpert]]$raw <- normdps
    
    pointsI <- seq(from=1,to=numrowslist[[i]][1],length=nump)
    
    tdpvec <- vector()
    
    tdpvec[1] <- 0
    for(j in 2:(length(pointsI)-1)){
      tfloor <- floor(pointsI[j])
      tceiling <- ceiling(pointsI[j])
      
      if((tceiling-tfloor)==0){
        tdpvec[j] <- normdps[pointsI[j]]
      }else{
        tslope <- (normdps[tceiling]-normdps[tfloor])/(tceiling-tfloor)
        tyintercept <- (normdps[tceiling]-(tslope*tceiling))
        
        tdpvec[j] <- (tslope*pointsI[j])+tyintercept
      }
    }
    tdpvec[length(pointsI)] <- 1
    
    retlist[[currpert]]$norm <- tdpvec
  }
  
  return(retlist)
}
