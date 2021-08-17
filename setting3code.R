setting30LOO <- function(ngoflist,nmetalist,
                         metavec,nrmetavec,
                         totpert,measures){
  
  numtestingcases <-  (length(nmetalist) * totpert)
  
  bestname <- paste0("Best.",measures)
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases, ncol = ((length(measures)*3)+2))
  
  for (i in 1:length(nmetalist)){
    # for(i in 10:10){
    cat(" ",i)
    
    picked <- i
    
    trainsetlist <- nmetalist[-picked]
    trainset <- trainsetlist[[1]]
    for(j in 2:length(trainsetlist)){
      trainset <- rbind(trainset,trainsetlist[[j]])
    }
    testset <- nmetalist[[picked]]
    
    for(j in 1:totpert){
      currpert <- paste0("P",as.character((j+(10-totpert-1))))
      
      
      trainsetsub <- trainset[trainset$Pert==currpert,]
      testsetsub <- testset[testset$Pert==currpert,]
      
      
      temptrainsub <- NULL
      temptrainsub <- normalization(trainsetsub,nrmetavec)
      
      ntrainsetsub <- NULL
      ntrainsetsub <- temptrainsub$df
      ntrainminmax <- NULL
      ntrainminmax <- temptrainsub$minmax
      
      for (k in 1:length(bestname)) {
        ntrainsetsub[bestname[k]] <- as.character(ntrainsetsub[,bestname[k]])
        ntrainsetsub[bestname[k]] <- as.factor(ntrainsetsub[,bestname[k]])
      }
      
      ntestsetsub <- NULL
      ntestsetsub <- normtestset(testsetsub,nrmetavec,ntrainminmax)
      
      right <- NULL
      right <- metavec[1]
      for (k in 2:length(metavec)){
        right <- paste0(right,"+",metavec[k])
      }
      
      for (k in 1:length(measures)) {
        targetformula <- NULL
        targetformula <- as.formula(paste0(bestname[k],"~",right))
        
        # SMOTE here
        countdf <- as.data.frame(table(ntrainsetsub[,bestname[k]]))
        majornum <- max(countdf$Freq)
        countbest <- countdf[countdf$Freq==majornum,]
        majorclass <- countbest[1,"Var1"]
        
        countminor <- countdf[countdf$Freq!=majornum,]
        minorclassvec <- countminor[,"Var1"]

        newntrainsetsub <- NULL

        majorsub <- ntrainsetsub[ntrainsetsub[,bestname[k]]==majorclass,c(metavec,bestname[k])]
        majorsubnofac <- majorsub
        majorsubnofac[bestname[k]] <- as.character(majorsubnofac[,bestname[k]])
        newntrainsetsub <- majorsubnofac

        if(length(minorclassvec) > 0){
          for(l in 1:length(minorclassvec)){

            minorsub <- ntrainsetsub[ntrainsetsub[,bestname[k]]==minorclassvec[l],c(metavec,bestname[k])]

            augtrainset <- rbind(majorsub,minorsub)
            augtrainset[bestname[k]] <- as.character(augtrainset[,bestname[k]])
            augtrainset[bestname[k]] <- as.factor(augtrainset[,bestname[k]])

            augtrainset <- smotethis(augtrainset,targetformula,metavec,bestname[k],ubConf,pK)

            augtrainset[bestname[k]] <- as.character(augtrainset[,bestname[k]])

            newntrainsetsub <- rbind(newntrainsetsub,augtrainset[augtrainset[,bestname[k]]==minorclassvec[l],])
          }
        }

        newntrainsetsub[bestname[k]] <- as.factor(newntrainsetsub[,bestname[k]])
        # newntrainsetsub <- ntrainsetsub # no SMOTE
        # SMOTE end
        
        bestmodel <- NULL
        
        # trainset contains only one model (do not need model selection)
        if(length(minorclassvec) == 0){
          bestmodel <- majorclass
        }else{
          rpartEP <- setting3train(newntrainsetsub,metavec,bestname[k])
          selectionEPret <- setting3predict(rpartEP,metavec,ntestsetsub)
          
          # thismax <- NULL
          # thismax <- max(selectionEPret)
          # bestmodelidx <- NULL
          # bestmodelidx <- which(selectionEPret==thismax)
          # 
          # bestmodel <- colnames(selectionEPret)[bestmodelidx]
          
          bestmodel <- selectionEPret
        }
        
        testcurr <- NULL
        testcurr <- ngoflist[[i]][ngoflist[[i]]$Pert==currpert,c("Model",measures[k],nmeasures[k])]
        
        testbestrow <- NULL
        # testbestrow <- testcurr[testcurr$Model %in% c(bestmodel),c(measures[k])]
        testbestrow <- testcurr[testcurr$Model==bestmodel,]
        
        if(nrow(testbestrow)>0){
          csvmat[(((i-1)*totpert)+j),k] <- mean(testbestrow[,measures[k]])
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- mean(testbestrow[,nmeasures[k]])
        }else{
          csvmat[(((i-1)*totpert)+j),k] <- Inf
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- -1
        }
        
        fatorlev <- NULL
        fatorlev <- levels(nmetalist[[i]][,bestname[k]])
        csvmat[(((i-1)*totpert)+j),(k+(length(measures)*2))] <- fatorlev[nmetalist[[i]][nmetalist[[i]]$Pert==currpert,bestname[k]]]
        
      }
      
      csvmat[(((i-1)*totpert)+j),((length(measures)*3)+1)] <- i
      csvmat[(((i-1)*totpert)+j),((length(measures)*3)+2)] <- currpert
    }
  }
  csvdf <- NULL
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures,bestname,'Num','Pert')
  
  print(" ")
  print("setting30")
  
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

setting31LOO <- function(ngoflist,nmetalist,
                         metavec,nrmetavec,
                         totpert,measures){
  
  numtestingcases <-  length(nmetalist) * totpert
  
  bestname <- paste0("Best.",measures)
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases, ncol = length(measures)*3)
  
  for (i in 1:length(nmetalist)){
  # for(i in 10:10){
    cat(" ",i)
    
    picked <- i
    
    trainsetlist <- nmetalist[-picked]
    trainset <- trainsetlist[[1]]
    for(j in 2:length(trainsetlist)){
      trainset <- rbind(trainset,trainsetlist[[j]])
    }
    testset <- nmetalist[[picked]]
    
    for(j in 1:totpert){
      currpert <- paste0("P",as.character((j+(10-totpert-1))))
      
      
      trainsetsub <- trainset[trainset$Pert==currpert,]
      testsetsub <- testset[testset$Pert==currpert,]
      
      
      temptrainsub <- NULL
      temptrainsub <- normalization(trainsetsub,nrmetavec)
      
      ntrainsetsub <- NULL
      ntrainsetsub <- temptrainsub$df
      ntrainminmax <- NULL
      ntrainminmax <- temptrainsub$minmax
      
      for (k in 1:length(bestname)) {
        ntrainsetsub[bestname[k]] <- as.character(ntrainsetsub[,bestname[k]])
        ntrainsetsub[bestname[k]] <- as.factor(ntrainsetsub[,bestname[k]])
      }
      
      ntestsetsub <- NULL
      ntestsetsub <- normtestset(testsetsub,nrmetavec,ntrainminmax)
      
      right <- NULL
      right <- metavec[1]
      for (k in 2:length(metavec)){
        right <- paste0(right,"+",metavec[k])
      }
      
      for (k in 1:length(measures)) {
        targetformula <- NULL
        targetformula <- as.formula(paste0(bestname[k],"~",right))
        
        # SMOTE here
        countdf <- as.data.frame(table(ntrainsetsub[,bestname[k]]))
        majornum <- max(countdf$Freq)
        countbest <- countdf[countdf$Freq==majornum,]
        majorclass <- countbest[1,"Var1"]

        countminor <- countdf[countdf$Freq!=majornum,]
        minorclassvec <- countminor[,"Var1"]
        # 
        # newntrainsetsub <- NULL
        # 
        # majorsub <- ntrainsetsub[ntrainsetsub[,bestname[k]]==majorclass,c(metavec,bestname[k])]
        # majorsubnofac <- majorsub
        # majorsubnofac[bestname[k]] <- as.character(majorsubnofac[,bestname[k]])
        # newntrainsetsub <- majorsubnofac
        # 
        # if(length(minorclassvec) > 0){
        #   for(l in 1:length(minorclassvec)){
        #     
        #     minorsub <- ntrainsetsub[ntrainsetsub[,bestname[k]]==minorclassvec[l],c(metavec,bestname[k])]
        #     
        #     augtrainset <- rbind(majorsub,minorsub)
        #     augtrainset[bestname[k]] <- as.character(augtrainset[,bestname[k]])
        #     augtrainset[bestname[k]] <- as.factor(augtrainset[,bestname[k]])
        #     
        #     augtrainset <- smotethis(augtrainset,targetformula,metavec,bestname[k],ubConf,pK)
        #     
        #     augtrainset[bestname[k]] <- as.character(augtrainset[,bestname[k]])
        #     
        #     newntrainsetsub <- rbind(newntrainsetsub,augtrainset[augtrainset[,bestname[k]]==minorclassvec[l],])
        #   }
        # }
        # 
        # newntrainsetsub[bestname[k]] <- as.factor(newntrainsetsub[,bestname[k]])
        newntrainsetsub <- ntrainsetsub # no SMOTE
        # SMOTE end
        
        bestmodel <- NULL
        
        # trainset contains only one model (do not need model selection)
        if(length(minorclassvec) == 0){
          bestmodel <- majorclass
        }else{
          rpartEP <- setting3train(newntrainsetsub,metavec,bestname[k])
          selectionEPret <- setting3predict(rpartEP,metavec,ntestsetsub)
          
          # thismax <- NULL
          # thismax <- max(selectionEPret)
          # bestmodelidx <- NULL
          # bestmodelidx <- which(selectionEPret==thismax)
          # 
          # bestmodel <- colnames(selectionEPret)[bestmodelidx]
          
          bestmodel <- selectionEPret
        }
        
        testcurr <- NULL
        testcurr <- ngoflist[[i]][ngoflist[[i]]$Pert==currpert,c("Model",measures[k],nmeasures[k])]
        
        testbestrow <- NULL
        # testbestrow <- testcurr[testcurr$Model %in% c(bestmodel),c(measures[k])]
        testbestrow <- testcurr[testcurr$Model==bestmodel,]
        
        if(nrow(testbestrow)>0){
          csvmat[(((i-1)*totpert)+j),k] <- mean(testbestrow[,measures[k]])
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- mean(testbestrow[,nmeasures[k]])
        }else{
          csvmat[(((i-1)*totpert)+j),k] <- Inf
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- -1
        }
        
        fatorlev <- NULL
        fatorlev <- levels(nmetalist[[i]][,bestname[k]])
        csvmat[(((i-1)*totpert)+j),(k+(length(measures)*2))] <- fatorlev[nmetalist[[i]][nmetalist[[i]]$Pert==currpert,bestname[k]]]
        
      }
    }
  }
  csvdf <- NULL
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures,bestname)
  
  print(" ")
  print("setting31")
  
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




setting32LOO <- function(ngoflist,
                         nmetalist,addnmetalist,
                         metavec,nrmetavec,
                         totpert,measures){
  
  numtestingcases <-  length(nmetalist) * totpert
  
  bestname <- paste0("Best.",measures)
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases, ncol = length(measures)*3)
  
  for (i in 1:length(nmetalist)){
  # for(i in 10:10){
    cat(" ",i)
    
    picked <- i
    
    
    trainsetlist <- nmetalist[-picked]
    addtrainsetlist <- addnmetalist[-picked]
    
    trainset <- trainsetlist[[1]]
    addtrainset <- addtrainsetlist[[1]]
    
    for(j in 2:length(trainsetlist)){
      trainset <- rbind(trainset,trainsetlist[[j]])
      addtrainset <- rbind(addtrainset,addtrainsetlist[[j]])
    }
    
    # here!!
    trainset <- rbind(trainset,addtrainset)
    # here!!
    
    testset <- nmetalist[[picked]]
    
    temptrain <- NULL
    temptrain <- normalization(trainset,nrmetavec)
    
    for(j in 1:totpert){
      currpert <- paste0("P",as.character((j+(10-totpert-1))))
      
      testsetsub <- testset[testset$Pert==currpert,]
      
      ntrainset <- NULL
      ntrainset <- temptrain$df
      ntrainminmax <- NULL
      ntrainminmax <- temptrain$minmax
      
      for (k in 1:length(bestname)) {
        ntrainset[bestname[k]] <- as.character(ntrainset[,bestname[k]])
        ntrainset[bestname[k]] <- as.factor(ntrainset[,bestname[k]])
      }
      
      ntestsetsub <- NULL
      ntestsetsub <- normtestset(testsetsub,nrmetavec,ntrainminmax)
      
      right <- NULL
      right <- metavec[1]
      for (k in 2:length(metavec)){
        right <- paste0(right,"+",metavec[k])
      }
      
      for (k in 1:length(measures)) {
        targetformula <- NULL
        targetformula <- as.formula(paste0(bestname[k],"~",right))
        
        # SMOTE here
        countdf <- as.data.frame(table(ntrainset[,bestname[k]]))
        majornum <- max(countdf$Freq)
        countbest <- countdf[countdf$Freq==majornum,]
        majorclass <- countbest[1,"Var1"]

        countminor <- countdf[countdf$Freq!=majornum,]
        minorclassvec <- countminor[,"Var1"]
        # 
        # newntrainset <- NULL
        # 
        # majorsub <- ntrainset[ntrainset[,bestname[k]]==majorclass,c(metavec,bestname[k])]
        # majorsubnofac <- majorsub
        # majorsubnofac[bestname[k]] <- as.character(majorsubnofac[,bestname[k]])
        # newntrainset <- majorsubnofac
        # 
        # if(length(minorclassvec) > 0){
        #   for(l in 1:length(minorclassvec)){
        #     
        #     minorsub <- ntrainset[ntrainset[,bestname[k]]==minorclassvec[l],c(metavec,bestname[k])]
        #     
        #     augtrainset <- rbind(majorsub,minorsub)
        #     augtrainset[bestname[k]] <- as.character(augtrainset[,bestname[k]])
        #     augtrainset[bestname[k]] <- as.factor(augtrainset[,bestname[k]])
        #     
        #     augtrainset <- smotethis(augtrainset,targetformula,metavec,bestname[k],ubConf,pK)
        #     
        #     augtrainset[bestname[k]] <- as.character(augtrainset[,bestname[k]])
        #     
        #     newntrainset <- rbind(newntrainset,augtrainset[augtrainset[,bestname[k]]==minorclassvec[l],])
        #   }
        # }
        # 
        # newntrainset[bestname[k]] <- as.factor(newntrainset[,bestname[k]])
        newntrainset <- ntrainset # no SMOTE
        # SMOTE end
        
        bestmodel <- NULL
        
        if(length(minorclassvec) == 0){
          # trainset contains only one model (do not need model selection)
          bestmodel <- majorclass
        }else{
          rpartEP <- setting3train(newntrainset,metavec,bestname[k])
          selectionEPret <- setting3predict(rpartEP,metavec,ntestsetsub)
          
          # thismax <- NULL
          # thismax <- max(selectionEPret)
          # bestmodelidx <- NULL
          # bestmodelidx <- which(selectionEPret==thismax)
          # 
          # bestmodel <- colnames(selectionEPret)[bestmodelidx]
          
          bestmodel <- selectionEPret
        }
        
        testcurr <- NULL
        testcurr <- ngoflist[[i]][ngoflist[[i]]$Pert==currpert,c("Model",measures[k],nmeasures[k])]
        
        # print(bestmodel)
        # print(testcurr)
        
        testbestrow <- NULL
        # testbestrow <- testcurr[testcurr$Model %in% bestmodel,c(measures[k])]
        testbestrow <- testcurr[testcurr$Model==bestmodel,]
        
        
        if(nrow(testbestrow)>0){
          csvmat[(((i-1)*totpert)+j),k] <- mean(testbestrow[,measures[k]])
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- mean(testbestrow[,nmeasures[k]])
        }else{
          csvmat[(((i-1)*totpert)+j),k] <- Inf
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- -1
        }
        
        fatorlev <- NULL
        fatorlev <- levels(nmetalist[[i]][,bestname[k]])
        csvmat[(((i-1)*totpert)+j),(k+(length(measures)*2))] <- fatorlev[nmetalist[[i]][nmetalist[[i]]$Pert==currpert,bestname[k]]]
        
      }
    }
  }
  csvdf <- NULL
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures,bestname)
  
  print(" ")
  print("setting32")
  
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

setting33LOO <- function(ngoflist,
                         nmetalist,addnmetalist,
                         dplist,adddplist,
                         basicpert,addipert,
                         metavec,nrmetavec,
                         totpert,measures){
  
  numtestingcases <-  length(nmetalist) * totpert
  
  bestname <- paste0("Best.",measures)
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases, ncol = length(measures)*3)
  
  for (i in 1:length(nmetalist)){
  # for(i in 10:10){
    cat(" ",i)
    
    picked <- i
    
    
    trainsetlist <- nmetalist[-picked]
    addtrainsetlist <- addnmetalist[-picked]
    
    testset <- nmetalist[[picked]]
    
    dptrainlist <- dplist[-picked]
    adddptrainlist <- adddplist[-picked]
    
    dptest <- dplist[[picked]]
    
    for(j in 1:length(basicpert)){
      currpert <- basicpert[j]
      
      testsetsub <- testset[testset$Pert==currpert,]
      
      dptestsub <- dptest[[currpert]]$norm
      
      candidatelist <- list()
      for(k in 1:length(trainsetlist)){
        for (l in 1:length(basicpert)) {
          templist <- list()
          templist$dp <- dptrainlist[[k]][[basicpert[l]]]$norm
          templist$train <- trainsetlist[[k]][trainsetlist[[k]]$Pert==basicpert[l],]
          candidatelist[[(length(candidatelist)+1)]] <- templist
        }
        for (l in 1:length(addipert)) {
          if(!is.null(adddptrainlist[[k]][[addipert[l]]])){
            templist <- list()
            templist$dp <- adddptrainlist[[k]][[addipert[l]]]$norm
            templist$train <- addtrainsetlist[[k]][addtrainsetlist[[k]]$Pert==addipert[l],]
            candidatelist[[(length(candidatelist)+1)]] <- templist
          }
        }
      }
      
      pvalvec <- vector()
      for (k in 1:length(candidatelist)) {
        ksret <- ks.test(dptestsub,candidatelist[[k]]$dp)
        candidatelist[[k]]$pval <- ksret$p.value
        pvalvec[(length(pvalvec)+1)] <- ksret$p.value
      }
      
      candidatelist <- candidatelist[order(-pvalvec)]
      
      trainset <- candidatelist[[1]]$train
      for (k in 2:100) {
        trainset <- rbind(trainset,candidatelist[[k]]$train)
      }
      pval <- 0.05
      for (k in 101:length(candidatelist)) {
        if(candidatelist[[k]]$pval>=pval){
          trainset <- rbind(trainset,candidatelist[[k]]$train)
        }
      }

      temptrain <- NULL
      temptrain <- normalization(trainset,nrmetavec)
      
      ntrainset <- NULL
      ntrainset <- temptrain$df
      ntrainminmax <- NULL
      ntrainminmax <- temptrain$minmax
      
      for (k in 1:length(bestname)) {
        ntrainset[bestname[k]] <- as.character(ntrainset[,bestname[k]])
        ntrainset[bestname[k]] <- as.factor(ntrainset[,bestname[k]])
      }
      
      ntestsetsub <- NULL
      ntestsetsub <- normtestset(testsetsub,nrmetavec,ntrainminmax)
      
      right <- NULL
      right <- metavec[1]
      for (k in 2:length(metavec)){
        right <- paste0(right,"+",metavec[k])
      }
      
      for (k in 1:length(measures)) {
        targetformula <- NULL
        targetformula <- as.formula(paste0(bestname[k],"~",right))
        
        # SMOTE here
        countdf <- as.data.frame(table(ntrainset[,bestname[k]]))
        majornum <- max(countdf$Freq)

        countbest <- countdf[countdf$Freq==majornum,]
        majorclass <- countbest[1,"Var1"]

        countminor <- countdf[countdf$Freq!=majornum,]
        minorclassvec <- countminor[,"Var1"]
        # 
        # newntrainset <- NULL
        # 
        # majorsub <- ntrainset[ntrainset[,bestname[k]]==majorclass,c(metavec,bestname[k])]
        # majorsubnofac <- majorsub
        # majorsubnofac[bestname[k]] <- as.character(majorsubnofac[,bestname[k]])
        # newntrainset <- majorsubnofac
        # 
        # if(length(minorclassvec) > 0){
        #   for(l in 1:length(minorclassvec)){
        #     
        #     minorsub <- ntrainset[ntrainset[,bestname[k]]==minorclassvec[l],c(metavec,bestname[k])]
        #     
        #     augtrainset <- rbind(majorsub,minorsub)
        #     augtrainset[bestname[k]] <- as.character(augtrainset[,bestname[k]])
        #     augtrainset[bestname[k]] <- as.factor(augtrainset[,bestname[k]])
        #     
        #     augtrainset <- smotethis(augtrainset,targetformula,metavec,bestname[k],ubConf,pK)
        #     
        #     augtrainset[bestname[k]] <- as.character(augtrainset[,bestname[k]])
        #     
        #     newntrainset <- rbind(newntrainset,augtrainset[augtrainset[,bestname[k]]==minorclassvec[l],])
        #   }
        # }
        # 
        # newntrainset[bestname[k]] <- as.factor(newntrainset[,bestname[k]])
        newntrainset <- ntrainset # no SMOTE
        # SMOTE end
        
        bestmodel <- NULL
        
        # trainset contains only one model (do not need model selection)
        if(length(minorclassvec) == 0){
          bestmodel <- majorclass
        }else{
          rpartEP <- setting3train(newntrainset,metavec,bestname[k])
          selectionEPret <- setting3predict(rpartEP,metavec,ntestsetsub)
          
          # thismax <- NULL
          # thismax <- max(selectionEPret)
          # bestmodelidx <- NULL
          # bestmodelidx <- which(selectionEPret==thismax)
          # 
          # bestmodel <- colnames(selectionEPret)[bestmodelidx]
          
          bestmodel <- selectionEPret
        }
        
        testcurr <- NULL
        testcurr <- ngoflist[[i]][ngoflist[[i]]$Pert==currpert,c("Model",measures[k],nmeasures[k])]
        
        testbestrow <- NULL
        # testbestrow <- testcurr[testcurr$Model %in% bestmodel,c(measures[k])]
        testbestrow <- testcurr[testcurr$Model==bestmodel,]
        
        if(nrow(testbestrow)>0){
          csvmat[(((i-1)*totpert)+j),k] <- mean(testbestrow[,measures[k]])
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- mean(testbestrow[,nmeasures[k]])
        }else{
          csvmat[(((i-1)*totpert)+j),k] <- Inf
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- -1
        }
        
        fatorlev <- NULL
        fatorlev <- levels(nmetalist[[i]][,bestname[k]])
        csvmat[(((i-1)*totpert)+j),(k+(length(measures)*2))] <- fatorlev[nmetalist[[i]][nmetalist[[i]]$Pert==currpert,bestname[k]]]
        
      }
    }
  }
  
  csvdf <- NULL
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures,bestname)
  
  print(" ")
  print("setting33")
  
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

