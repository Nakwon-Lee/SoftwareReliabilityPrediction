
# setting1 Leave-one-out wrapper

setting1LOO <- function(goflist,ngoflist,
                        pEst,datalist,
                        gofregvec,gofclvec,
                        totpert,
                        targets,targetcl,
                        measures,isFCTBF){
  
  numtestingcases <-  length(goflist) * totpert
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases, ncol = ((length(measures)*2)+2))
  
  nmeasures <- paste0("N",measures)
  
  for (i in 1:length(goflist)){
  # for (i in 2:2){
    cat(" ",i)
    
    picked <- i
    
    trainsetlist <- goflist[-picked]
    
    trainset <- trainsetlist[[1]]
    for(j in 2:length(trainsetlist)){
      trainset <- rbind(trainset,trainsetlist[[j]])
    }
    testset <- goflist[[picked]]
    
    observed <- NULL
    if(isFCTBF=="FC"){
      observed <- datalist[[picked]]$n
    }else if(isFCTBF=="TBF"){
      tobsvd <- MovingAverage(datalist[[picked]]$f,7)
      observed <- vector()
      for(j in 1:length(tobsvd)){
        observed[j] <- sum(tobsvd[1:j])
      }
    }
    
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      
      trainsetsub <- trainset[trainset$Pert==currpert,c(gofregvec,targets)]
      testsetsub <- testset[testset$Pert==currpert,]
      
      if(nrow(testsetsub)==0){
        # there is no model for this test data
        
        csvmat[(((i-1)*totpert)+j),1] <- Inf
        csvmat[(((i-1)*totpert)+j),3] <- -1
        csvmat[(((i-1)*totpert)+j),2] <- Inf
        csvmat[(((i-1)*totpert)+j),4] <- -1
        
      }else if(nrow(testsetsub)==1){
        # there is only one appropriate model for this test data
        
        for (k in 1:length(measures)) {
          tep <- testsetsub[1,c(measures[k])]
          csvmat[(((i-1)*totpert)+j),k] <- tep
          csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- normVal(tep,ngoflist[[picked]][ngoflist[[picked]]$Pert==currpert,measures[k]])
        }
        
      }else{
        
        adtreemd <- setting1traincl(trainset[,c(gofclvec,targetcl)],gofclvec,targetcl)
        underover <- setting1predictcl(adtreemd,testsetsub)
        
        for (k in 1:length(measures)) {
          rpartEP <- setting1train(trainsetsub,gofregvec,targets[k])
          selectionEPret <- setting1predict(rpartEP,gofregvec,testsetsub)
          
          testcurr <- selectionEPret
          testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
          testcurr.ordered.top <- NULL
          if(nrow(testcurr.ordered)>3){
            testcurr.ordered.top <- testcurr.ordered[1:3,]
          }else{
            testcurr.ordered.top <- testcurr.ordered
          }
          
          pEst[[currpert]][[picked]][["dummy"]] <- list()
          pEst[[currpert]][[picked]][["dummy"]]$EstElap <- rep(0,times = length(observed))
          
          under <- list()
          under$Model <- "dummy"
          under$Weight <- 0
          over <- list()
          over$Model <- "dummy"
          over$Weight <- 0
          
          for(l in 1:nrow(testcurr.ordered.top)){
            bestmodel <- testcurr.ordered.top[l,"Model"]
            if(underover[underover$Model==bestmodel,"Pred"]=="U"){
              if(under$Weight < testcurr.ordered.top[l,"Pred"]){
                under$Model <- bestmodel
                under$Weight <- testcurr.ordered.top[l,"Pred"]
              }
            }else if(underover[underover$Model==bestmodel,"Pred"]=="O"){
              if(over$Weight < testcurr.ordered.top[l,"Pred"]){
                over$Model <- bestmodel
                over$Weight <- testcurr.ordered.top[l,"Pred"]
              }
            }
          }
          
          wavgEst <- NULL
          
          if(under$Model == "dummy" && over$Model == "dummy"){
            
            pEstsub <- pEst[[currpert]][[picked]]
            pEstsub.this <- NULL
            
            pEstsub.this <- pEstsub[[testcurr.ordered.top[1,"Model"]]]$EstElap
            
            wavgEst <- pEstsub.this
            
          }else{
            # prediction using under and over weights
            pEstsub <- pEst[[currpert]][[picked]]
            pEstsub.Under <- NULL
            pEstsub.Over <- NULL
            pEstsub.Under <- pEstsub[[under$Model]]$EstElap
            pEstsub.Over <- pEstsub[[over$Model]]$EstElap
            
            wavgEst <- ((pEstsub.Under * under$Weight) + (pEstsub.Over * over$Weight)) / (under$Weight + over$Weight)
            
          }
          
          # print(pEstsub[[testcurr.ordered.top[1,"Model"]]]$EstElap)
          # print(testcurr.ordered.top[1,"Model"])
          # print(testcurr.ordered.top[1,"Pred"])
          
          if(measures[k]=="EP"){
            tep <- FC.EP(length(observed),observed,wavgEst)
            csvmat[(((i-1)*totpert)+j),k] <- tep
            csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- tep
          }else if(measures[k]=="PE"){
            tpe <- TBF.PE(length(observed),observed,wavgEst)
            csvmat[(((i-1)*totpert)+j),k] <- tpe
            csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- tpe
          }else if(measures[k]=="MEOP"){
            tmeop <- FC.MEOP(pEstsub$trainH,length(observed),observed,wavgEst)
            csvmat[(((i-1)*totpert)+j),k] <- tmeop
            csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- tmeop
          }else if(measures[k]=="AE"){
            tae <- TBF.AE(pEstsub$trainH,length(observed),observed,wavgEst)
            csvmat[(((i-1)*totpert)+j),k] <- tae
            csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- tae
          }
        }
      }
      
      csvmat[(((i-1)*totpert)+j),((length(measures)*2)+1)] <- picked
      csvmat[(((i-1)*totpert)+j),((length(measures)*2)+2)] <- currpert
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures,'Num','Pert')
  
  print(" ")
  print("setting1")
  
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


setting111LOO <- function(goflist,
                        gofregvec,
                        totpert,
                        target,
                        measure,
                        fntrain,fnpredict){
  
  csvdf <- setting111LOOnoout(goflist = goflist,gofregvec = gofregvec,
                              totpert = totpert,target = target,
                              measure = measure,
                              fntrain = fntrain,fnpredict = fnpredict)
  
  print(" ")
  print("setting111")
  prefix <- paste0(measure," Median.: ")
  print(c(prefix,median(csvdf[,measure])))
  
  return(csvdf)
}

setting111LOOnoout <- function(goflist,
                          gofregvec,
                          totpert,
                          target,
                          measure,
                          fntrain,fnpredict){
  
  resultvec <- vector()
  selvec <- vector()
  setnumvec <- vector()
  pertvec <- vector()
  
  for (i in 1:length(goflist)){
    
    picked <- i
    
    trainsetlist <- goflist[-picked]
    
    trainset <- trainsetlist[[1]]
    for(j in 2:length(trainsetlist)){
      trainset <- rbind(trainset,trainsetlist[[j]])
    }
    testset <- goflist[[picked]]
    
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      
      setnumvec[(length(setnumvec)+1)] <- i
      pertvec[(length(pertvec)+1)] <- currpert
      
      trainsetsub <- trainset[trainset$Pert==currpert,c(gofregvec,target)]
      testsetsub <- testset[testset$Pert==currpert,]
      
      rpartEP <- fntrain(trainsetsub,gofregvec,target)
      selectionEPret <- fnpredict(rpartEP,gofregvec,testsetsub)
      
      testcurr <- selectionEPret
      testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
      
      resultvec[(length(resultvec)+1)] <- testcurr.ordered[1,measure]
      selvec[(length(selvec)+1)] <- testcurr.ordered[1,"Model"]
    }
  }
  
  csvdf <- data.frame(resultvec)
  names(csvdf) <- c(measure)
  csvdf$Num <- setnumvec
  csvdf$Pert <- pertvec
  
  return(csvdf)
}

setting111LOOPert <- function(goflist,
                              targetdf,
                              gofregvec,
                              pcri){
  
  traindf <- listrbinder(goflist)
  testdf <- targetdf
  
  rpartEP <- settingGLMtrain(traindf,gofregvec,pcri)
  
  selectionEPret <- settingBasicpredict(rpartEP,gofregvec,testdf)
  
  testcurr <- selectionEPret
  testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
  
  return(testcurr.ordered[1,'Model'])
}

execGOF111 <- function(goflist,
                       gofregvec,
                       totpert,
                       ptarget,
                       fntrain,fnpredict,
                       tgoflist){
  
  trainset <- goflist[[1]]
  for(j in 2:length(goflist)){
    trainset <- rbind(trainset,goflist[[j]])
  }
  
  models <- list()
  for (j in 1:totpert) {
    currpert <- paste0("P",as.character((10-totpert-1)+j))
    trainsetsub <- trainset[trainset$Pert==currpert,c(gofregvec,paste0("N",ptarget))]
    
    models[[j]] <- fntrain(trainsetsub,gofregvec,paste0("N",ptarget))
  }
  
  selectionEPret <- list()
  for (i in 1:length(tgoflist)) {
    testset <- tgoflist[[i]]
    selectionEPret[[i]] <- list()
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      testsetsub <- testset[testset$Pert==currpert,]
      tempret <- fnpredict(models[[j]],gofregvec,testsetsub)
      if(is.null(selectionEPret[[i]]$df)){
        selectionEPret[[i]]$df <- tempret
      }else{
        selectionEPret[[i]]$df <- rbind(selectionEPret[[i]]$df,tempret)
      }
      testcurr <- tempret
      testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
      if(is.null(selectionEPret[[i]]$sel)){
        selectionEPret[[i]]$sel <- testcurr.ordered[1,]
      }else{
        selectionEPret[[i]]$sel <- rbind(selectionEPret[[i]]$sel,testcurr.ordered[1,])
      }
    }
  }
  
  return(selectionEPret)
}

execGOF111Pert <- function(goflist,
                       gofregvec,
                       cpert,
                       ptarget,
                       fntrain,fnpredict,
                       tgoflist){
  
  currpert <- cpert
  
  trainset <- listrbinder(goflist)
  trainsetsub <- trainset[trainset$Pert==currpert,c(gofregvec,paste0("N",ptarget))]
  
  model <- fntrain(trainsetsub,gofregvec,paste0("N",ptarget))
  
  selectionEPret <- list()
  for (i in 1:length(tgoflist)) {
    testset <- tgoflist[[i]]
    selectionEPret[[i]] <- list()
    
    
    
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      testsetsub <- testset[testset$Pert==currpert,]
      tempret <- fnpredict(models[[j]],gofregvec,testsetsub)
      if(is.null(selectionEPret[[i]]$df)){
        selectionEPret[[i]]$df <- tempret
      }else{
        selectionEPret[[i]]$df <- rbind(selectionEPret[[i]]$df,tempret)
      }
      testcurr <- tempret
      testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
      if(is.null(selectionEPret[[i]]$sel)){
        selectionEPret[[i]]$sel <- testcurr.ordered[1,]
      }else{
        selectionEPret[[i]]$sel <- rbind(selectionEPret[[i]]$sel,testcurr.ordered[1,])
      }
    }
  }
  
  return(selectionEPret)
}

# setting111testdata <- function(goflist,ngoflist,
#                                gofregvec,
#                                totpert,
#                                targets,
#                                measures,isFCTBF,
#                                fntrain,fnpredict){
#   
#   numtestingcases <-  length(goflist) * totpert
#   
#   csvmat <- NULL
#   csvmat <- matrix(nrow = numtestingcases, ncol = (length(measures)*2))
#   
#   csvdf <- NULL
#   
#   trainsetlist <- goflist
#   
#   trainset <- trainsetlist[[1]]
#   for(j in 2:length(trainsetlist)){
#     trainset <- rbind(trainset,trainsetlist[[j]])
#   }
#   
#   for (i in 1:length(goflist)){
#     for (j in 1:totpert) {
#       currpert <- paste0("P",as.character((10-totpert-1)+j))
#     }
#   }
#   
#   for (k in 1:length(measures)) {
#     rpartEP <- fntrain(trainset,gofregvec,targets[k])
#     selectionEPret <- fnpredict(rpartEP,gofregvec,trainset)
#     
#     if(is.null(csvdf)){
#       csvdf <- data.frame(selectionEPret$Pred)
#     }
#   }
# }


# setting2 Leave-one-out wrapper

setting2LOO <- function(goflist,addgoflist,
                        ngoflist,
                        gofregvec,
                        totpert,
                        targets,
                        measures,isFCTBF,
                        fntrain,fnpredict){
  
  numtestingcases <-  length(goflist) * totpert
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)*2))
  
  for (i in 1:length(goflist)){
    
    cat(" ",i)
    
    picked <- i
    
    trainsetlist <- goflist[-picked]
    addtrainsetlist <- addgoflist[-picked]
    trainset <- trainsetlist[[1]]
    addtrainset <- addtrainsetlist[[1]]
    for(j in 2:length(trainsetlist)){
      trainset <- rbind(trainset,trainsetlist[[j]])
      addtrainset <- rbind(addtrainset,addtrainsetlist[[j]])
    }
    
    # here!!
    trainset <- rbind(trainset,addtrainset)
    # here!!
    testset <- goflist[[picked]]
    
    for (k in 1:length(measures)) {
      rpartEP <- fntrain(trainset,gofregvec,targets[k])
      
      for (j in 1:totpert) {
        
        currpert <- paste0("P",as.character((10-totpert-1)+j))
        
        testsetsub <- testset[testset$Pert==currpert,]
        
        selectionEPret <- fnpredict(rpartEP,gofregvec,testsetsub)
        
        testcurr <- selectionEPret
        testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
        
        ngofsub <- ngoflist[[picked]][ngoflist[[picked]]$Pert==currpert,]
        ngofsub <- ngofsub[ngofsub$Model==testcurr.ordered[1,"Model"],]
        
        csvmat[(((i-1)*totpert)+j),k] <- testcurr.ordered[1,measures[k]]
        csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- ngofsub[1,nmeasures[k]]
      }
      
    }
    
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures)
  
  print(" ")
  print("setting2")
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
    prefix <- paste0(nmeasures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,nmeasures[i]])))
  }
  
  return(csvdf)
}



setting13LOO <- function(goflist,addgoflist,
                         ngoflist,
                         dplist,adddplist,
                         gofregvec,
                         basicpert,addipert,
                         totpert,
                         targets,
                         measures,isFCTBF,
                         pvalthr,
                         fntrain,fnpredict){
  
  numtestingcases <-  length(goflist) * totpert
  
  nmeasures <- paste0("R",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)*2))
  
  for (i in 1:length(goflist)) {
    cat(" ",i)
    
    picked <- i
    
    trainsetlist <- goflist[-picked]
    addtrainsetlist <- addgoflist[-picked]
    
    testset <- goflist[[picked]]
    
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
      pval <- pvalthr
      for (k in 101:length(candidatelist)) {
        if(candidatelist[[k]]$pval>=pval){
          trainset <- rbind(trainset,candidatelist[[k]]$train)
        }
      }
      
      
      
      for (k in 1:length(measures)) {
        
        rpartEP <- fntrain(trainset,gofregvec,targets[k])
        
        selectionEPret <- fnpredict(rpartEP,gofregvec,testsetsub)
        
        testcurr <- selectionEPret
        testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
        
        ngofsub <- ngoflist[[picked]][ngoflist[[picked]]$Pert==currpert,]
        ngofsub <- ngofsub[ngofsub$Model==testcurr.ordered[1,"Model"],]
        
        csvmat[(((i-1)*totpert)+j),k] <- testcurr.ordered[1,measures[k]]
        csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- ngofsub[1,nmeasures[k]]
      }
      
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures)
  
  print(" ")
  print("setting13")
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
    prefix <- paste0(nmeasures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,nmeasures[i]])))
  }
  
  return(csvdf)
  
}



setting14LOO <- function(ngoflist,
                        gofregvec,
                        totpert,
                        targets,
                        measures,isFCTBF,
                        fntrain,fnpredict){
  
  numtestingcases <-  length(ngoflist) * totpert
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases, ncol = (length(measures)*2))
  
  for (i in 1:length(ngoflist)){
    cat(" ",i)
    
    picked <- i
    
    trainsetlist <- ngoflist[-picked]
    
    trainset <- trainsetlist[[1]]
    for(j in 2:length(trainsetlist)){
      trainset <- rbind(trainset,trainsetlist[[j]])
    }
    testset <- ngoflist[[picked]]
    
    for (j in 1:totpert) {
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      
      trainsetsub <- trainset[trainset$Pert==currpert,c(gofregvec,targets)]
      testsetsub <- testset[testset$Pert==currpert,]
      
      for (k in 1:length(measures)) {
        
        rpartEP <- fntrain(trainsetsub,gofregvec,targets[k])
        
        selectionEPret <- fnpredict(rpartEP,gofregvec,testsetsub)
        
        testcurr <- selectionEPret
        testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
        
        ngofsub <- testsetsub[testsetsub$Model==testcurr.ordered[1,"Model"],]
        
        csvmat[(((i-1)*totpert)+j),k] <- testcurr.ordered[1,measures[k]]
        csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- ngofsub[1,nmeasures[k]]
      }
      
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures)
  
  print(" ")
  print("setting14")
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
    prefix <- paste0(nmeasures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,nmeasures[i]])))
  }
  
  return(csvdf)
}



setting15LOO <- function(ngoflist,naddgoflist,
                        gofregvec,
                        totpert,
                        targets,
                        measures,isFCTBF,
                        fntrain,fnpredict){
  
  numtestingcases <-  length(ngoflist) * totpert
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)*2))
  
  for (i in 1:length(ngoflist)){
    
    cat(" ",i)
    
    picked <- i
    
    trainsetlist <- ngoflist[-picked]
    addtrainsetlist <- naddgoflist[-picked]
    trainset <- trainsetlist[[1]]
    addtrainset <- addtrainsetlist[[1]]
    for(j in 2:length(trainsetlist)){
      trainset <- rbind(trainset,trainsetlist[[j]])
      addtrainset <- rbind(addtrainset,addtrainsetlist[[j]])
    }
    
    # here!!
    trainset <- rbind(trainset,addtrainset)
    # here!!
    testset <- ngoflist[[picked]]
    
    
    for (k in 1:length(measures)) {
      
      rpartEP <- fntrain(trainset,gofregvec,targets[k])
      
      for (j in 1:totpert) {
        
        currpert <- paste0("P",as.character((10-totpert-1)+j))
        
        testsetsub <- testset[testset$Pert==currpert,]
        
        selectionEPret <- fnpredict(rpartEP,gofregvec,testsetsub)
        
        testcurr <- selectionEPret
        testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
        
        ngofsub <- testsetsub[testsetsub$Model==testcurr.ordered[1,"Model"],]
        
        csvmat[(((i-1)*totpert)+j),k] <- testcurr.ordered[1,measures[k]]
        csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- ngofsub[1,nmeasures[k]]
      }
      
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures)
  
  print(" ")
  print("setting15")
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
    prefix <- paste0(nmeasures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,nmeasures[i]])))
  }
  
  return(csvdf)
}

setting16LOO <- function(ngoflist,naddgoflist,
                         dplist,adddplist,
                         gofregvec,
                         basicpert,addipert,
                         totpert,
                         targets,
                         measures,isFCTBF,
                         pvalthr,
                         fntrain,fnpredict){
  
  
  numtestingcases <-  length(ngoflist) * totpert
  
  nmeasures <- paste0("N",measures)
  
  csvmat <- NULL
  csvmat <- matrix(nrow = numtestingcases,ncol = (length(measures)*3))
  
  
  for (i in 1:length(ngoflist)){
    # for (i in 22:22){
    
    cat(" ",i)
    
    picked <- i
    
    trainsetlist <- ngoflist[-picked]
    addtrainsetlist <- naddgoflist[-picked]
    
    testset <- ngoflist[[picked]]
    
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
      pval <- pvalthr
      for (k in 101:length(candidatelist)) {
        if(candidatelist[[k]]$pval>=pval){
          trainset <- rbind(trainset,candidatelist[[k]]$train)
        }
      }
      
      
      for (k in 1:length(measures)) {
        rpartEP <- fntrain(trainset,gofregvec,targets[k])
        
        selectionEPret <- fnpredict(rpartEP,gofregvec,testsetsub)
        
        testcurr <- selectionEPret
        testcurr.ordered <- testcurr[order(-testcurr$Pred,-testcurr$NMSE),]
        
        ngofsub <- testsetsub[testsetsub$Model==testcurr.ordered[1,"Model"],]
        
        csvmat[(((i-1)*totpert)+j),k] <- testcurr.ordered[1,measures[k]]
        csvmat[(((i-1)*totpert)+j),(k+length(measures))] <- ngofsub[1,nmeasures[k]]
        csvmat[(((i-1)*totpert)+j),(k+(length(measures)*2))] <- testcurr.ordered[1,"Model"]
      }
    }
  }
  
  csvdf <- data.frame(csvmat)
  
  names(csvdf) <- c(measures,nmeasures)
  
  for (i in 1:length(measures)) {
    csvdf[measures[i]] <- as.numeric(csvdf[,measures[i]])
    csvdf[nmeasures[i]] <- as.numeric(csvdf[,nmeasures[i]])
  }
  
  print(" ")
  print("setting16")
  
  for (i in 1:length(measures)) {
    prefix <- paste0(measures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,measures[i]])))
    prefix <- paste0(nmeasures[i]," Avg.: ")
    print(c(prefix,median(csvdf[,nmeasures[i]])))
  }
  
  return(csvdf)
}


