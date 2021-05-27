normalization <- function(df, vars){
  nvars <- paste0("N",vars)
  
  minmaxlist <- list()
  
  for(i in 1:length(vars)){
    
    minval <- min(df[,vars[i]])
    maxval <- max(df[,vars[i]])
    
    minmaxlist[[vars[i]]] <- c(minval,maxval)
    
    df[nvars[i]] <- (maxval - df[,vars[i]])/(maxval - minval)
  }
  
  retlist <- list()
  retlist$df <- df
  retlist$minmax <- minmaxlist
  
  return(retlist)
}

normtestset <- function(df, vars, minmax){
  nvars <- paste0("N",vars)
  for(i in 1:length(vars)){
    df[nvars[i]] <- (minmax[[vars[i]]][2] - df[,vars[i]])/(minmax[[vars[i]]][2] - minmax[[vars[i]]][1])
  }
  return(df)
}

powerSet <- function(pX){
  pick <- pX[1]
  
  retset <- list()
  
  if(length(pX)==1){
    retset <- list(c("em"),c("em",pick))
  }else{
    remainingset <- pX[2:length(pX)]
    
    nopickset <- powerSet(remainingset)
    
    pickset <- list()
    lenpickset <- length(pickset)
    
    for(i in 1:length(nopickset)){
      tvec <- c(pick,nopickset[[i]])
      pickset[[(lenpickset+i)]] <- tvec
    }
    
    retset <- append(retset,pickset)
    retset <- append(retset,nopickset)
  }
  
  return(retset)
}

gridSearch <- function(trainsetlist,Grd,thismetavec,thisnmetavec,totpert,target,key){
  
  bestGrd <- list()
  
  if(!key){
    bestGrd$eps <- sample(Grd$eps,1)
    bestGrd$cost <- sample(Grd$cost,1)
    bestGrd$K <- sample(Grd$K,1)
    
    return(bestGrd)
  }
  
  bestGrd$eps <- Grd$eps[1]
  bestGrd$cost <- Grd$cost[1]
  bestGrd$K <- Grd$K[1]
  
  bestfitval <- 1
  
  for(i in 1:length(Grd$eps)){
    for(j in 1:length(Grd$cost)){
      for(k in 1:length(Grd$K)){
        
        # constructing trainset and test set for cross-validation
        
        newnewlist <- trainsetlist
        
        numdatasets <- length(newnewlist)
        
        Xcross <- 6
        
        sublists <- list()
        
        for(l in 1:Xcross){
          sublists[[l]] <- list()
        }
        
        for(l in 1:numdatasets){
          curridxfords <- sample(1:length(newnewlist),1)
          
          thisidx <- (((l-1) %/% Xcross)+1)
          
          sublists[[thisidx]][[(length(sublists[[thisidx]])+1)]] <- newnewlist[[curridxfords]]
          
          newnewlist <- newnewlist[-c(curridxfords)]
        }
        
        fitval <- NULL
        
        fitvalvec <- vector()
        
        # for(l in 1:Xcross){
        for(l in 1:1){
          
          subtrainset <- NULL
          
          subtrainlist <- sublists[-c(l)]
          
          for(m in 1:length(subtrainlist)){
            for(n in 1:length(subtrainlist[[m]])){
              if(is.null(subtrainset)){
                subtrainset <- subtrainlist[[m]][[n]]
              }else{
                subtrainset <- rbind(subtrainset,subtrainlist[[m]][[n]])
              }
            }
          }
          
          temptrain <- NULL
          temptrain <- normalization(subtrainset,thismetavec)
          
          ntrainset <- NULL
          ntrainset <- temptrain$df
          ntrainminmax <- NULL
          ntrainminmax <- temptrain$minmax
          
          newtraindata <- NULL
          newtraindata <- ntrainset[,c(thisnmetavec,target)]
          
          currGrd <- list()
          currGrd$eps <- Grd$eps[i]
          currGrd$cost <- Grd$cost[j]
          currGrd$K <- Grd$K[k]
          
          rpartmeta <- NULL
          rpartmeta <- buildDTMetaLearning3(newtraindata,currGrd,thisnmetavec,target)
          
          testlist <- NULL
          
          testlist <- sublists[[l]]
          
          errorvec <- vector()
          
          for(m in 1:length(testlist)){
            
            ntestset <- NULL
            ntestset <- normtestset(testlist[[m]],thismetavec,ntrainminmax)
            
            predictedbestmeta <- NULL
            predictedbestmeta <- predictUsingDTMeta3(rpartmeta,totpert,ntestset)
            
            for(n in 1:totpert){
              currpert <- paste0("P",as.character((n+(10-totpert-1))))
              
              bestmodel <- NULL
              bestmodel <- predictedbestmeta[[currpert]][[1]]
              
              testcurr <- NULL
              testcurr <- testlist[[m]][testlist[[m]]$Pert==currpert,]
              
              if(nrow(testcurr)!=1){
                print(currpert)
                print(testcurr)
              }
              
              stopifnot(nrow(testcurr)==1)
              
              currerror <- abs((testcurr[1,target] - bestmodel))
              
              errorvec[(length(errorvec)+1)] <- currerror
            }
          }
          
          stopifnot(length(errorvec)==(length(testlist)*totpert))
          
          meanerrorofpred <- NULL
          meanerrorofpred <- mean(errorvec)
          
          fitvalvec[(length(fitvalvec)+1)] <- meanerrorofpred
        } #Xcross
        
        fitval <- mean(fitvalvec)
        
        # is the current fitval best?
        
        if(bestfitval > fitval){
          bestfitval <- fitval
          
          bestGrd <- list()
          bestGrd$eps <- Grd$eps[i]
          bestGrd$cost <- Grd$cost[j]
          bestGrd$K <- Grd$K[k]
        }
        
      } #k
    } #j
  } #i
  
  return(bestGrd)
}



MovingAverage <- function(pVals,order){
  
  retvals <- vector()
  
  for(i in 1:(order-1)){
    retvals[i] <- mean(pVals[1:i])
  }
  
  for(i in order:length(pVals)){
    retvals[i] <- mean(pVals[(i-(order-1)):i])
  }
  
  return(retvals)
}

RapairDDMEST <- function(ddmest,df,models,totpert,pRound,isFCTBF){
  numrows <- nrow(df)
  
  numrowssubs <- list()
  
  dflearn <- NULL
  
  if(isFCTBF=="FC"){
    dflearn <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    dflearn <- vector()
    for(i in 1:length(tobsvd)){
      dflearn[i] <- sum(tobsvd[1:i])
    }
  }
  
  retlist <- list()
  
  for(i in 1:totpert){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    tempstart <- (((10-totpert)-1)+i)
    trainH <- floor((tempstart*0.1)*numrows)
    numrowssubs[[(length(numrowssubs)+1)]] <- trainH
    
    retlist[[currpert]] <- list()
    
    estimated <- NULL
    fdim <- NULL
    
    dflearnsub <- dflearn[1:trainH]
    
    for(j in 1:length(models)){
      if(models[j]=="SVR"){
        
        maxval <- max(dflearnsub)
        minval <- min(dflearnsub)
        
        denominator <- maxval-minval
        
        estimated <- ddmest[[pRound]][[currpert]][[models[j]]]$Est
        fdim <- ddmest[[pRound]][[currpert]][[models[j]]]$dim
        
        estimatedr <- estimated[-c(1:fdim)]
        
        estimatedr <- (((estimatedr-0.1)*denominator)/0.8)+minval
        
        estimated <- c(estimated[1:fdim],estimatedr)
        
      }else if(models[j]=="ANN"){
        
        estimated <- ddmest[[pRound]][[currpert]][[models[j]]]$Est
        fdim <- ddmest[[pRound]][[currpert]][[models[j]]]$dim
        
        estimatedr <- estimated[-c(1:fdim)]
        
        estimatedr <- (((estimatedr+1)*(max(dflearnsub)-min(dflearnsub)))/2)+min(dflearnsub)
        
        estimated <- c(estimated[1:fdim],estimatedr)
        
      }else{
        stopifnot(FALSE)
      }
      
      stopifnot(length(estimated)==numrows)
      
      retlist[[currpert]][[models[j]]] <- estimated
    }
    
  }
  
  return(retlist)
}

RapairDDMESTfive <- function(ddmest,df,models,totpert,pRound,isFCTBF){
  numrows <- nrow(df)
  
  numrowssubs <- list()
  
  retlist <- list()
  
  dflearn <- NULL
  observed <- NULL
  
  if(isFCTBF=="FC"){
    dflearn <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    dflearn <- vector()
    for(i in 1:length(tobsvd)){
      dflearn[i] <- sum(tobsvd[1:i])
    }
  }
  
  if(numrows>=20){
    for(i in 1:totpert){
      currpert <- paste0("T",as.character((((10-totpert-1)+i)*10)+5))
      tempstart <- ((((10-totpert)-1)+i)*10)+5
      trainH <- floor((tempstart/100)*numrows)
      numrowssubs[[(length(numrowssubs)+1)]] <- trainH
      
      retlist[[currpert]] <- list()
      
      dflearnsub <- dflearn[1:trainH]
      
      estimated <- NULL
      fdim <- NULL
    
      for(j in 1:length(models)){
        if(models[j]=="SVR"){
          
          maxval <- max(dflearnsub)
          minval <- min(dflearnsub)
          
          denominator <- maxval-minval
          
          estimated <- ddmest[[pRound]][[currpert]][[models[j]]]$Est
          fdim <- ddmest[[pRound]][[currpert]][[models[j]]]$dim
          
          estimatedr <- estimated[-c(1:fdim)]
          
          estimatedr <- (((estimatedr-0.1)*denominator)/0.8)+minval
          
          estimated <- c(estimated[1:fdim],estimatedr)
          
        }else if(models[j]=="ANN"){
          
          estimated <- ddmest[[pRound]][[currpert]][[models[j]]]$Est
          fdim <- ddmest[[pRound]][[currpert]][[models[j]]]$dim
          
          estimatedr <- estimated[-c(1:fdim)]
          
          estimatedr <- (((estimatedr+1)*(max(dflearnsub)-min(dflearnsub)))/2)+min(dflearnsub)
          
          estimated <- c(estimated[1:fdim],estimatedr)
          
        }else{
          stopifnot(FALSE)
        }
        
        stopifnot(length(estimated)==numrows)
        
        retlist[[currpert]][[models[j]]] <- estimated
      }
      
    }
  }
  
  return(retlist)
}

normVal <- function(val,valvec){
  
  normval <- NULL
  
  avalvec <- c(val,valvec)
  
  if(length(avalvec)>1){
    infidx <- which(is.infinite(avalvec))
    
    noinfsub <- NULL
    if(length(infidx)!=0){
      noinfsub <- avalvec[-c(infidx)]
    }else{
      noinfsub <- avalvec
    }
    
    minval <- min(noinfsub)
    maxval <- max(noinfsub)
    
    denominator <- maxval-minval
    
    if(is.infinite(val)){
      normval <- -1
    }else{
      if(denominator>0){
        normval <- 0.1 + ((maxval - val)*0.9)/denominator
      }else{
        normval <- 1
      }
    }
    
  }else{
    normval <- 1
  }
  
  return(normval)
}
