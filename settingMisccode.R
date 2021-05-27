# setting MSE

settingMSE <- function(goflist,ngoflist,totpert,measures){
  
  numtestingcases <- length(goflist) * totpert
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)*2))
  
  for (i in 1:length(goflist)){
    
    cat(" ",i)
    
    picked <- i
    
    testset <- goflist[[picked]]
    
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      
      testcurr <- NULL
      testcurr <- testset[testset$Pert==currpert,]
      
      msemin <- NULL
      msemin <- max(testcurr$NMSE)
      
      bestrow <- testcurr[testcurr$NMSE==msemin,]
      
      ngofsub <- ngoflist[[picked]][ngoflist[[picked]]$Pert==currpert,]
      ngofsub <- ngofsub[ngofsub$Model==bestrow[1,"Model"],]
      
      csvmat[(((i-1)*totpert)+j),1] <- mean(bestrow[,measures[1]])
      csvmat[(((i-1)*totpert)+j),3] <- ngofsub[1,nmeasures[1]]
      csvmat[(((i-1)*totpert)+j),2] <- mean(bestrow[,measures[2]])
      csvmat[(((i-1)*totpert)+j),4] <- ngofsub[1,nmeasures[2]]
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures)
  
  print(" ")
  print("settingMSE")
  
  for (i in 1:length(measures)) {
    csvdf[measures[i]] <- as.numeric(csvdf[,measures[i]])
    csvdf[nmeasures[i]] <- as.numeric(csvdf[,nmeasures[i]])
  }
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
    prefix <- paste0(nmeasures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,nmeasures[i]])))
  }
  
  return(csvdf)
}

settingDDM <- function(ngofddmlist,totpert,measures){
  
  numtestingcases <- length(ngofddmlist) * totpert
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)))
  
  for (i in 1:length(ngofddmlist)){
    cat(" ",i)
    
    picked <- i
    
    testset <- ngofddmlist[[picked]]
    
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      
      testcurr <- NULL
      testcurr <- testset[testset$Pert==currpert,]
      
      for (k in 1:length(measures)) {
        csvmat[(((i-1)*totpert)+j),k] <- testcurr[testcurr$Model=="SVR",measures[k]]
      }
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures)
  
  print(" ")
  print("settingDDM")
  
  for (i in 1:length(measures)) {
    csvdf[measures[i]] <- as.numeric(csvdf[,measures[i]])
  }
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
  }
  
  return(csvdf)
}

settingDDMnoout <- function(ngofddmlist,totpert,measures){
  
  numtestingcases <- length(ngofddmlist) * totpert
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)))
  
  for (i in 1:length(ngofddmlist)){
    
    picked <- i
    
    testset <- ngofddmlist[[picked]]
    
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      
      testcurr <- NULL
      testcurr <- testset[testset$Pert==currpert,]
      
      for (k in 1:length(measures)) {
        csvmat[(((i-1)*totpert)+j),k] <- testcurr[testcurr$Model=="SVR",measures[k]]
      }
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures)
  
  # print(" ")
  # print("settingDDM")
  
  for (i in 1:length(measures)) {
    csvdf[measures[i]] <- as.numeric(csvdf[,measures[i]])
  }
  
  # for (i in 1:length(measures)) {
  #   prefix <- paste0(measures[i]," Avg.: ")
  #   print(c(prefix,median(csvdf[,measures[i]])))
  # }
  
  return(csvdf)
}



# settingBest

settingBest <- function(ngoflist,totpert,measures){
  
  numtestingcases <- length(ngoflist) * totpert
  
  nmeasures <- paste0("R",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)*3))
  
  for (i in 1:length(ngoflist)){
    cat(" ",i)
    
    picked <- i
    
    goftestset <- ngoflist[[picked]]
    
    for(j in 1:totpert){
      currpert <- paste0("P",as.character((j+(10-totpert-1))))
      
      testcurr <- goftestset[goftestset$Pert==currpert,c("Model",measures,nmeasures)]
      
      testcurr.ordered <- testcurr[order(testcurr[,measures[1]]),]
      
      csvmat[(((i-1)*totpert)+j),1] <- testcurr.ordered[1,measures[1]]
      csvmat[(((i-1)*totpert)+j),3] <- testcurr.ordered[1,nmeasures[1]]
      csvmat[(((i-1)*totpert)+j),5] <- testcurr.ordered[1,"Model"]
      
      testcurr.ordered <- testcurr[order(testcurr[,measures[2]]),]
      
      csvmat[(((i-1)*totpert)+j),2] <- testcurr.ordered[1,measures[2]]
      csvmat[(((i-1)*totpert)+j),4] <- testcurr.ordered[1,nmeasures[2]]
      csvmat[(((i-1)*totpert)+j),6] <- testcurr.ordered[1,"Model"]
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures,paste0("Best",measures))
  
  print(" ")
  print("settingBest")
  
  for (i in 1:length(measures)) {
    csvdf[measures[i]] <- as.numeric(csvdf[,measures[i]])
    csvdf[nmeasures[i]] <- as.numeric(csvdf[,nmeasures[i]])
  }
  
  for (i in 1:length(measures)) {
    csvdf[paste0("Best",measures[i])] <- as.factor(csvdf[,paste0("Best",measures[i])])
  }
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
    prefix <- paste0(nmeasures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,nmeasures[i]])))
  }
  
  return(csvdf)
}

settingSRGMfromDDM <- function(goflist,totpert,measures){
  
  numtestingcases <- length(goflist) * totpert
  
  nmeasures <- paste0("N",measures)
  
  ccnamecri <- paste0("F",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)*2))
  
  for (i in 1:length(goflist)){
    
    cat(" ",i)
    
    picked <- i
    
    testset <- goflist[[picked]]
    
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      
      testcurr <- NULL
      testcurr <- testset[testset$Pert==currpert,]
      
      for (k in 1:length(measures)) {
        testset.ordered <- testcurr[order(testcurr$FMEOP),]
        csvmat[(((i-1)*totpert)+j),k] <- testset.ordered[1,measures[k]]
        csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- testset.ordered[1,nmeasures[k]]
      }
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures)
  
  print(" ")
  print("settingSRGMfromDDM")
  
  for (i in 1:length(measures)) {
    csvdf[measures[i]] <- as.numeric(csvdf[,measures[i]])
    csvdf[nmeasures[i]] <- as.numeric(csvdf[,nmeasures[i]])
  }
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
    prefix <- paste0(nmeasures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,nmeasures[i]])))
  }
  
  return(csvdf)
}