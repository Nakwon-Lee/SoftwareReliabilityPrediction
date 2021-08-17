settingDGS <- function(goflist,gofnlist,gofddmlist,
                       metadf,pcri,totpert,
                       goffts,metafts){
  
  resultvec <- vector()
  
  for (i in 1:length(goflist)) {
    cat(paste0(i," "))
    
    picked <- c(i)
    
    DGSretdf <- execDGS(goflist = goflist[-picked],gofnlist = gofnlist[-picked],
                        gofddmlist = gofddmlist[-picked],metadf = metadf[-picked],
                        pcri = pcri,totpert = totpert,
                        goffts = goffts,metafts = metafts,
                        psearch = "SEARCH",
                        tmetadf = metadf[picked],tgofddmlist = gofddmlist[picked],
                        tgoflist = goflist[picked])
    
    for (j in 1:length(picked)) {
      testirows <- c(1:totpert)
      testirows <- testirows + (totpert * (j-1))
      
      temdgsret <- DGSretdf[testirows,]
      
      for (k in 1:nrow(temdgsret)) {
        if (temdgsret[k,paste0("DGS.RET.",pcri)]=="DDM"){
          trows <- which(gofddmlist[[picked[j]]]$Model=="SVR" & gofddmlist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
          resultvec[(length(resultvec)+1)] <- gofddmlist[[picked[j]]][trows,][1,pcri]
        }else{
          trows <- which(goflist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
          tgofdf <- goflist[[picked[j]]][trows,]
          resultvec[(length(resultvec)+1)] <- tgofdf[tgofdf$Model==temdgsret[k,paste0("DGS.RET.",pcri)],][1,pcri]
        }
      }
    }
  }
    
  retdf <- data.frame(resultvec)
  colnames(retdf) <- c(pcri)
  
  print(c(pcri,median(retdf[,pcri])))
  
  return(retdf)
}

settingDGSsepPert <- function(goflist,gofnlist,gofddmlist,
                              metadf,pcri,cpert,
                              goffts,metafts){
  
  
  
  resultvec <- vector()
  
  for (i in 1:length(goflist)) {
    cat(paste0(i," "))
    
    picked <- c(i)
    
    DGSretdf <- execDGS(goflist = goflist[-picked],gofnlist = gofnlist[-picked],
                        gofddmlist = gofddmlist[-picked],metadf = metadf[-picked],
                        pcri = pcri,totpert = totpert,
                        goffts = goffts,metafts = metafts,
                        psearch = "SEARCH",
                        tmetadf = metadf[picked],tgofddmlist = gofddmlist[picked],
                        tgoflist = goflist[picked])
    
    for (j in 1:length(picked)) {
      testirows <- c(1:totpert)
      testirows <- testirows + (totpert * (j-1))
      
      temdgsret <- DGSretdf[testirows,]
      
      for (k in 1:nrow(temdgsret)) {
        if (temdgsret[k,paste0("DGS.RET.",pcri)]=="DDM"){
          trows <- which(gofddmlist[[picked[j]]]$Model=="SVR" & gofddmlist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
          resultvec[(length(resultvec)+1)] <- gofddmlist[[picked[j]]][trows,][1,pcri]
        }else{
          trows <- which(goflist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
          tgofdf <- goflist[[picked[j]]][trows,]
          resultvec[(length(resultvec)+1)] <- tgofdf[tgofdf$Model==temdgsret[k,paste0("DGS.RET.",pcri)],][1,pcri]
        }
      }
    }
  }
  
  retdf <- data.frame(resultvec)
  colnames(retdf) <- c(pcri)
  
  print(c(pcri,median(retdf[,pcri])))
  
  return(retdf)
}

# settingDGStest <- function(goflist,gofnlist,gofddmlist,
#                        metadf,pcri,totpert,
#                        goffts,metafts,
#                        gofretdf,ddmretdf){
#   
#   param <- list()
#   param$orig <- 1
#   param$margin <- 1
#   
#   resultvec <- vector()
#   DGSretdf <- NULL
#   
#   for (i in 1:length(goflist)) {
#     cat(paste0(i," "))
#     
#     picked <- c(i)
#     
#     tempret <- execDGS(goflist = goflist[-picked],gofnlist = gofnlist[-picked],
#                        gofddmlist = gofddmlist[-picked],metadf = metadf[-picked],
#                        pcri = pcri,totpert = totpert,
#                        goffts = goffts,metafts = metafts,
#                        psearch = "GIVEN",givenparam = param,
#                        tmetadf = metadf[picked],tgofddmlist = gofddmlist[picked],
#                        tgoflist = goflist[picked])
#     
#     if (is.null(DGSretdf)){
#       DGSretdf <- tempret
#     }else{
#       DGSretdf <- rbind(DGSretdf,tempret)
#     }
#   }
#   
#   resultvec <- ifelse(DGSretdf[,paste0("DGS.RET.",pcri)]=="DDM",ddmretdf[,pcri],gofretdf[,pcri])
#   
#   retdf <- data.frame(resultvec)
#   colnames(retdf) <- c(pcri)
#   
#   print(c(pcri,median(retdf[,pcri])))
#   
#   return(retdf)
# }

settingDGSsavedEvaluator <- function(goflist,gofnlist,gofddmlist,
                                     metadf,pcri,totpert,
                                     goffts,metafts,givenprefix){
  
  resultvec <- vector()
  
  for (i in 1:length(goflist)) {
    cat(paste0(i," "))
    
    picked <- c(i)
    
    DGSretdf <- execDGS(goflist = goflist[-picked],gofnlist = gofnlist[-picked],
                        gofddmlist = gofddmlist[-picked],metadf = metadf[-picked],
                        pcri = pcri,totpert = totpert,
                        goffts = goffts,metafts = metafts,
                        psearch = "SAVED",givenevalfile = paste0(givenprefix,".",i,".",pcri,".RData"),
                        tmetadf = metadf[picked],tgofddmlist = gofddmlist[picked],
                        tgoflist = goflist[picked])
    
    for (j in 1:length(picked)) {
      testirows <- c(1:totpert)
      testirows <- testirows + (totpert * (j-1))
      
      temdgsret <- DGSretdf[testirows,]
      
      for (k in 1:nrow(temdgsret)) {
        if (temdgsret[k,paste0("DGS.RET.",pcri)]=="DDM"){
          trows <- which(gofddmlist[[picked[j]]]$Model=="SVR" & gofddmlist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
          resultvec[(length(resultvec)+1)] <- gofddmlist[[picked[j]]][trows,][1,pcri]
        }else{
          trows <- which(goflist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
          tgofdf <- goflist[[picked[j]]][trows,]
          resultvec[(length(resultvec)+1)] <- tgofdf[tgofdf$Model==temdgsret[k,paste0("DGS.RET.",pcri)],][1,pcri]
        }
      }
    }
  }
  
  retdf <- data.frame(resultvec)
  colnames(retdf) <- c(pcri)
  
  print(c(pcri,median(retdf[,pcri])))
  
  return(retdf)
}

settingDGSfixedParam <- function(goflist,gofnlist,gofddmlist,
                                     metadf,pcri,totpert,
                                     goffts,metafts,
                                 pfiltvec=NULL,
                                 givenorimar=NULL,
                                 givenfts,
                                 givenfile = NULL,
                                 goffromddm){
  param <- NULL
  
  if(is.null(givenfile)){
    param <- list()
    param$orig <- 1
    param$margin <- 1
  }else{
    load(givenfile,envir = environment())
    param <- searchedEvaluator[[pcri]]$param
  }
  
  if(!is.null(givenorimar)){
    param <- givenorimar
  }
  
  print(param)
  print(givenfts)
  
  DGSres <- NULL
  resultvec <- vector()
  resultvec2 <- vector()
  
  for (i in 1:length(goflist)) {
    cat(paste0(i," "))

    picked <- c(i)

    DGSretdf <- execDGS(goflist = goflist[-picked],gofnlist = gofnlist[-picked],
                        gofddmlist = gofddmlist[-picked],metadf = metadf[-picked],
                        pcri = pcri,totpert = totpert,
                        goffts = goffts,metafts = metafts,
                        psearch = "GIVEN",givenparam = param,
                        tmetadf = metadf[picked],tgofddmlist = gofddmlist[picked],
                        tgoflist = goflist[picked],givenfts = givenfts)
    
    if(is.null(DGSres)){
      DGSres <- DGSretdf
    }else{
      DGSres <- rbind(DGSres,DGSretdf)
    }

    for (j in 1:length(picked)) {
      testirows <- c(1:totpert)
      testirows <- testirows + (totpert * (j-1))

      temdgsret <- DGSretdf[testirows,]
      
      tgofddm <- gofddmlist[[picked[j]]]
      
      tgoffromddm <- goffromddm[[picked[j]]]
      
      tgofdf <- goflist[[picked[j]]]

      for (k in 1:nrow(temdgsret)) {
        if (temdgsret[k,paste0("DGS.RET.",pcri)]=="DDM"){
          
          currgofddm <- tgofddm[tgofddm$Pert==temdgsret[k,"Pert"],]
          resultvec2[(length(resultvec2)+1)] <- currgofddm[currgofddm$Model=='SVR',][1,pcri]
          
          currsrgmddm <- tgoffromddm[tgoffromddm$Pert==temdgsret[k,"Pert"],]
          currsrgmddmodered <- currsrgmddm[order(currsrgmddm$MSE),]
          
          resultvec[(length(resultvec)+1)] <- currsrgmddmodered[1,pcri]
          
        }else{
          
          currgofdf <- tgofdf[tgofdf$Pert==temdgsret[k,"Pert"],]
          resultvec[(length(resultvec)+1)] <- currgofdf[currgofdf$Model==temdgsret[k,paste0("DGS.RET.",pcri)],][1,pcri]
          resultvec2[(length(resultvec2)+1)] <- currgofdf[currgofdf$Model==temdgsret[k,paste0("DGS.RET.",pcri)],][1,pcri]
        }
      }
    }
  }
  
  if(!is.null(pfiltvec)){
    filvec <- pfiltvec
    
    DGSres$FIL <- filvec
    
    filtnum <- nrow(na.omit(DGSres)) / nrow(DGSres)
    filtDGSres <- na.omit(DGSres)
    
    filtddmdf <- filtDGSres[which(filtDGSres$FIL=="DDM"),]
    filtgofdf <- filtDGSres[which(filtDGSres$FIL!="DDM"),]
    
    ddmacc <- nrow(filtddmdf[which(filtddmdf$FIL==filtddmdf[,paste0("DGS.RET.",pcri)]),]) / nrow(filtddmdf)
    gofacc <- nrow(filtgofdf[which(filtgofdf$FIL==filtgofdf[,paste0("DGS.RET.",pcri)]),]) / nrow(filtgofdf)
    
    retdf <- data.frame(resultvec)
    colnames(retdf) <- c(pcri)
    
    ret <- list()
    ret$val <- median(retdf[,pcri])
    ret$avg <- mean(retdf[,pcri])
    ret$acc <- (ddmacc + gofacc)/2
    ret$filtnum <- filtnum
    
    print(ret)
  }
  
  retdf <- data.frame(resultvec)
  colnames(retdf) <- c(pcri)
  
  print(c(pcri,median(retdf[,pcri])))
  
  print(c("DDM ",pcri,median(resultvec2)))
  
  return(retdf)
}

settingDGSfixedParamPert <- function(goflist,gofnlist,gofddmlist,
                                 metadf,pcri,totpert,
                                 goffts,metafts,
                                 pfiltvec=NULL,
                                 givenorimar=NULL,
                                 givenfts,
                                 givenfile = NULL){
  param <- NULL
  
  if(is.null(givenfile)){
    param <- list()
    param$orig <- 1
    param$margin <- 1
  }else{
    load(givenfile,envir = environment())
    param <- searchedEvaluator[[pcri]]$param
  }
  
  if(!is.null(givenorimar)){
    param <- givenorimar
  }
  
  print(param)
  print(givenfts)
  
  DGSres <- NULL
  resultvec <- vector()
  
  for (i in 1:length(goflist)) {
    cat(paste0(i," "))
    
    picked <- c(i)
    
    DGSretdf <- execDGS(goflist = goflist[-picked],gofnlist = gofnlist[-picked],
                        gofddmlist = gofddmlist[-picked],metadf = metadf[-picked],
                        pcri = pcri,totpert = totpert,
                        goffts = goffts,metafts = metafts,
                        psearch = "GIVEN",givenparam = param,
                        tmetadf = metadf[picked],tgofddmlist = gofddmlist[picked],
                        tgoflist = goflist[picked],givenfts = givenfts)
    
    if(is.null(DGSres)){
      DGSres <- DGSretdf
    }else{
      DGSres <- rbind(DGSres,DGSretdf)
    }
    
    for (j in 1:length(picked)) {
      testirows <- c(1:totpert)
      testirows <- testirows + (totpert * (j-1))
      
      temdgsret <- DGSretdf[testirows,]
      
      for (k in 1:nrow(temdgsret)) {
        if (temdgsret[k,paste0("DGS.RET.",pcri)]=="DDM"){
          trows <- which(gofddmlist[[picked[j]]]$Model=="SVR" & gofddmlist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
          resultvec[(length(resultvec)+1)] <- gofddmlist[[picked[j]]][trows,][1,pcri]
        }else{
          trows <- which(goflist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
          tgofdf <- goflist[[picked[j]]][trows,]
          resultvec[(length(resultvec)+1)] <- tgofdf[tgofdf$Model==temdgsret[k,paste0("DGS.RET.",pcri)],][1,pcri]
        }
      }
    }
  }
  
  if(!is.null(pfiltvec)){
    filvec <- pfiltvec
    
    DGSres$FIL <- filvec
    
    filtnum <- nrow(na.omit(DGSres)) / nrow(DGSres)
    filtDGSres <- na.omit(DGSres)
    
    filtddmdf <- filtDGSres[which(filtDGSres$FIL=="DDM"),]
    filtgofdf <- filtDGSres[which(filtDGSres$FIL!="DDM"),]
    
    ddmacc <- nrow(filtddmdf[which(filtddmdf$FIL==filtddmdf[,paste0("DGS.RET.",pcri)]),]) / nrow(filtddmdf)
    gofacc <- nrow(filtgofdf[which(filtgofdf$FIL==filtgofdf[,paste0("DGS.RET.",pcri)]),]) / nrow(filtgofdf)
    
    retdf <- data.frame(resultvec)
    colnames(retdf) <- c(pcri)
    
    ret <- list()
    ret$val <- median(retdf[,pcri])
    ret$avg <- mean(retdf[,pcri])
    ret$acc <- (ddmacc + gofacc)/2
    ret$filtnum <- filtnum
    
    print(ret)
  }
  
  retdf <- data.frame(resultvec)
  colnames(retdf) <- c(pcri)
  
  print(c(pcri,median(retdf[,pcri])))
  
  return(retdf)
}

settingDGSEvalGensep <- function(pkey,
                                 goflist,gofnlist,gofddmlist,
                                 metadf,pcri,totpert,
                                 goffts,metafts,givenprefix){
  
  cat(paste0(pkey," "))
  
  picked <- c(pkey)
  
  execDGSEvaluatorGen(goflist = goflist[-picked],gofnlist = gofnlist[-picked],
                      gofddmlist = gofddmlist[-picked],metadf = metadf[-picked],
                      pcri = pcri,totpert = totpert,
                      goffts = goffts,metafts = metafts,
                      givenevalfile = paste0(givenprefix,".",pkey,".",pcri,".RData"))
}

# settingDGSEvalNCV <- function(pgofldf,pddmodf,goflist,gofnlist,gofddmlist,metadf,crivec,totpert,goffts,metafts){
#   
#   rcrivec <- paste0("R",crivec)
#   
#   DGSpred <- list()
#   for(i in 1:length(crivec)){
#     DGSpred[[crivec[i]]] <- NULL
#   }
#   
#   # normalizing meta features for training and test data
#   metafulldf <- NULL
#   for (i in 1:length(metadf)) {
#     if(is.null(metafulldf)){
#       metafulldf <- metadf[[i]]
#     }else{
#       metafulldf <- rbind(metafulldf,metadf[[i]])
#     }
#   }
#   metafullndf <- normalization(metafulldf,metafts)$df
#   
#   fullfeat <- c(goffts,paste0("N",metafts))
#   
#   evalfullvec <- list()
#   for (i in 1:length(crivec)){
#     evalfullvec[[crivec[i]]] <- vector()
#   }
#   
#   for (i in 1:length(goflist)) {
#     
#     print(i)
#     
#     picked <- c(i)
#     testirows <- c(1:totpert)
#     testirows <- testirows + (totpert * (i-1))
#     testmetaftsdf <- metafullndf[testirows,]
#     testgofiftsdf <- gofddmlist[[i]][gofddmlist[[i]]$Model=="SVR",]
#     
#     ftsdf <- cbind(testgofiftsdf,testmetaftsdf)
#     
#     for (j in 1:length(crivec)){
#       ddmeval <- searchforEvaluatorNoCV(plist = goflist[-picked],pnlist = gofnlist[-picked],
#                                         pmetalist = metadf[-picked],regvars = goffts,ppert = totpert,
#                                         pcri = crivec[j],pddmres = pddmodf[-testirows,],
#                                         pgofddmlist = gofddmlist[-picked],goffeats = goffts,metafeats = metafts)
#       
#       print(paste0(crivec[j]," ",ddmeval$performance))
#       print(ddmeval$fts)
#       print(ddmeval$string)
#       
#       #classification
#       predretdf <- setting1predictcl(ddmeval$evaluator,ddmeval$fts,ftsdf)
#       
#       #regression
#       # predretdf <- setting1predict(ddmeval$evaluator,ddmeval$fts,ftsdf)
#       # predretdf$Pred <- ifelse(predretdf$Pred>0.5,"DDM","GOF")
#       
#       if(is.null(DGSpred[[crivec[j]]])){
#         DGSpred[[crivec[j]]] <- predretdf
#       }else{
#         DGSpred[[crivec[j]]] <- rbind(DGSpred[[crivec[j]]],predretdf)
#       }
#       
#       evalvec <- filteringLabeler(pddmodf[testirows,crivec[j]],pgofldf[testirows,crivec[j]],ddmeval$orimar)
#       evalfullvec[[crivec[j]]] <- c(evalfullvec[[crivec[j]]],evalvec)
#     }
#   }
#   
#   for (i in 1:length(crivec)) {
#     vvec <- which(!is.na(evalfullvec[[crivec[i]]]))
#     ttt <- length(which(evalfullvec[[crivec[i]]][vvec]==DGSpred[[crivec[i]]]$Pred[vvec]))
#     bbb <- length(vvec)
#     print(paste0(crivec[i]," top: ",ttt," bot: ",bbb))
#     print(paste0("evalres: ",(ttt/bbb)))
#   }
#   
#   DGSret <- list()
#   colnvec <- vector()
#   
#   for (i in 1:length(crivec)){
#     DGSret[[crivec[i]]] <- ifelse(DGSpred[[crivec[i]]]$Pred=="DDM",pddmodf[,crivec[i]],pgofldf[,crivec[i]])
#     colnvec[(length(colnvec)+1)] <- crivec[i]
#     DGSret[[rcrivec[i]]] <- ifelse(DGSpred[[crivec[i]]]$Pred=="DDM",pddmodf[,rcrivec[i]],pgofldf[,rcrivec[i]])
#     colnvec[(length(colnvec)+1)] <- rcrivec[i]
#   }
#   
#   retdf <- NULL
#   
#   for (i in 1:length(colnvec)){
#     if(is.null(retdf)){
#       retdf <- data.frame(DGSret[[colnvec[i]]])
#       colnames(retdf) <- c(colnvec[i])
#     }else{
#       retdf[colnvec[i]] <- DGSret[[colnvec[i]]]
#     }
#   }
#   
#   for (i in 1:length(crivec)) {
#     prefix <- paste0(crivec[i]," Avg.: ")
#     print(c(prefix,median(retdf[,crivec[i]])))
#     prefix <- paste0(rcrivec[i]," Avg.: ")
#     print(c(prefix,median(retdf[,rcrivec[i]])))
#   }
#   
#   return(retdf)
# }

execDGS <- function(pgofldf=NULL,pddmodf=NULL,ptgofdf = NULL,
                    goflist,
                    gofnlist,gofddmlist,metadf,
                    pcri,totpert,goffts,metafts,
                    givenevalfile=NULL,psearch=NULL,
                    givenparam=NULL,givenfts=NULL,
                    tmetadf,tgofddmlist,tgoflist){
  
  goflret <- pgofldf
  if (is.null(goflret)){
    goflret <- setting111LOOnoout(goflist = goflist, gofregvec = paste0("N",goffts),
                                  totpert = totpert,target = paste0("N",pcri),
                                  measure = pcri,
                                  fntrain = settingGLMtrain, fnpredict = settingBasicpredict)
  }
  
  ddmret <- pddmodf
  if (is.null(ddmret)){
    ddmret <- settingDDMnoout(ngofddmlist = gofddmlist,
                         totpert = totpert,measures = c(pcri))
  }
  
  metafulldf <- NULL
  for (i in 1:length(metadf)) {
    if(is.null(metafulldf)){
      metafulldf <- metadf[[i]]
    }else{
      metafulldf <- rbind(metafulldf,metadf[[i]])
    }
  }
  metafullndf <- normalization(metafulldf,metafts)$df
  
  goffulldf <- NULL
  for (i in 1:length(gofddmlist)) {
    if(is.null(goffulldf)){
      goffulldf <- gofddmlist[[i]]
    }else{
      goffulldf <- rbind(goffulldf,gofddmlist[[i]])
    }
  }
  goffulldf <- goffulldf[goffulldf$Model=="SVR",]
  
  fulldf <- cbind(metafullndf,goffulldf)
  fullfeat <- c(goffts,paste0("N",metafts))
  
  ddmeval <- NULL
  
  if(is.null(psearch)){
    ddmeval <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                               pfulldf = fulldf,pfullfeats = fullfeat,
                               pcri = pcri)
  }else{
    if(psearch=="SEARCH"){
      ddmeval <- searchforEvaluator(plist = goflist,pnlist = gofnlist,pmetalist = metadf,
                                    regvars = paste0("N",goffts),ppert = totpert,pcri = pcri,
                                    pgofddmlist = gofddmlist, goffeats = goffts,metafeats = metafts,
                                    pddmres = ddmret,pgofres = goflret,
                                    pfulldf = fulldf,pfullfeats = fullfeat)
    }else if(psearch=="SAVED"){
      load(givenevalfile) # load the ddmeval
    }else if(psearch=="GIVEN"){
      ddmeval <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                                 pfulldf = fulldf,pfullfeats = fullfeat,
                                 orimar = givenparam,pselfts = givenfts,
                                 pcri = pcri)
    }else if(psearch=="DIRECT"){
      xddm <- log2(ddmret[,pcri])
      ygof <- log2(goflret[,pcri])
      origvec <- sort(c(xddm,ygof))
      margvec <- sort(abs(ygof-xddm))
      
      param <- list()
      param$oidx <- ceiling(length(origvec)/4)
      param$orig <- origvec[param$oidx]
      param$midx <- ceiling(length(margvec)/4)
      param$margin <- margvec[param$midx]
      
      ddmeval <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                                 pfulldf = fulldf,pfullfeats = fullfeat,
                                 orimar = param,
                                 pcri = pcri)
    }
  }
  
  tmetafulldf <- NULL
  for (i in 1:length(tmetadf)) {
    if(is.null(tmetafulldf)){
      tmetafulldf <- tmetadf[[i]]
    }else{
      tmetafulldf <- rbind(tmetafulldf,tmetadf[[i]])
    }
  }
  tmetafulldf2 <- rbind(tmetafulldf,metafulldf)
  tmetafullndf <- normalization(tmetafulldf2,metafts)$df
  tmetafullndf <- tmetafullndf[1:nrow(tmetafulldf),]
  
  tgoffulldf <- NULL
  for (i in 1:length(tgofddmlist)) {
    if(is.null(tgoffulldf)){
      tgoffulldf <- tgofddmlist[[i]]
    }else{
      tgoffulldf <- rbind(tgoffulldf,tgofddmlist[[i]])
    }
  }
  tgoffulldf <- tgoffulldf[tgoffulldf$Model=="SVR",]
  
  stopifnot(nrow(tmetafullndf)==nrow(tgoffulldf))
  
  tfulldf <- cbind(tmetafullndf,tgoffulldf)
  
  #classification
  predretdf <- setting1predictcl(ddmeval$evaluator,ddmeval$fts,tfulldf)
  
  gofretlist <- ptgofdf
  if(is.null(gofretlist)){
    gofretlist <- execGOF111(goflist = goflist,gofregvec = paste0("N",goffts),totpert = totpert,
                             ptarget = pcri,fntrain = settingGLMtrain,fnpredict = settingBasicpredict,
                             tgoflist = tgoflist)
  }
  
  gofretfulldf <- NULL
  for (i in 1:length(gofretlist)) {
    if (is.null(gofretfulldf)){
      gofretfulldf <- gofretlist[[i]]$sel
    }else{
      gofretfulldf <- rbind(gofretfulldf,gofretlist[[i]]$sel)
    }
  }
  
  stopifnot(nrow(predretdf)==nrow(gofretfulldf))
  
  predretdf[paste0("DGS.RET.",pcri)] <- ifelse(predretdf$Pred=="DDM","DDM",gofretfulldf$Model)
  
  # for (i in 1:length(totpert)) {
  #   if (predretdf$Pred=="DDM"){
  #     fitSRGMfromDDMSingle(df=tdataf,pEst = tddmest,modelvec = psrgmvec,
  #                          pert = (i+3),isFCTBF = 'FC')
  #   }
  # }
  
  return(predretdf)
}

# N data, N+1 prediction
execDGSPert <- function(pgofldf,pddmodf,ptgofdf = NULL,
                    goflist,
                    gofnlist,gofddmlist,metadf,
                    pcri,cpert,goffts,metafts,
                    givenevalfile=NULL,psearch=NULL,
                    givenparam=NULL,givenfts=NULL,
                    tmetadf,tgofddmlist,tgoflist,
                    modelSelection,dataDriven){
  
  modelselres <- NULL
  for (i in 1:length(gofnlist)) {
    if(is.null(modelselres)){
      modelselres <- modelSelection(traininglist,targetlist)
    }else{
      modelselres <- rbind(modelselres,modelSelection(traininglist,targetlist))
    }
  }
  
  datadrivenEst <- list()
  for (i in 1:length(datadf)) {
    datadrivenEst[[i]] <- dataDriven(targetdf)
  }
  
  datadrivenres <- NULL
  for (i in 1:length(datadf)) {
    if(is.null(datadrivenres)){
      datadrivenres <- CalcGOFDDMlist(targetdf,datadrivenEst[[i]])
    }else{
      datadrivenres <- rbind(datadrivenres,CalcGOFDDMlist(targetdf,datadrivenEst[[i]]))
    }
  }
  
  
  
  goflret <- pgofldf
  
  ddmret <- pddmodf
  
  metafulldf <- listrbinder(metadf)
  metafulldf <- metafulldf[metafulldf$Pert==cpert,]
  metafullndf <- normalization(metafulldf,metafts)$df
  
  goffulldf <- listbinder(gofddmlist)
  goffulldf <- goffulldf[goffulldf$Pert==cpert,]
  goffulldf <- goffulldf[goffulldf$Model=="SVR",]
  
  fulldf <- cbind(metafullndf,goffulldf)
  fullfeat <- c(goffts,paste0("N",metafts))
  
  ddmeval <- NULL
  
  if(is.null(psearch)){
    ddmeval <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                               pfulldf = fulldf,pfullfeats = fullfeat,
                               pcri = pcri)
  }else{
    if(psearch=="SEARCH"){
      ddmeval <- searchforEvaluator(plist = goflist,pnlist = gofnlist,pmetalist = metadf,
                                    regvars = paste0("N",goffts),ppert = cpert,pcri = pcri,
                                    pgofddmlist = gofddmlist, goffeats = goffts,metafeats = metafts,
                                    pddmres = ddmret,pgofres = goflret,
                                    pfulldf = fulldf,pfullfeats = fullfeat)
    }else if(psearch=="GIVEN"){
      ddmeval <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                                 pfulldf = fulldf,pfullfeats = fullfeat,
                                 orimar = givenparam,pselfts = givenfts,
                                 pcri = pcri)
    }
  }
  
  tmetafulldf <- NULL
  for (i in 1:length(tmetadf)) {
    if(is.null(tmetafulldf)){
      tmetafulldf <- tmetadf[[i]]
    }else{
      tmetafulldf <- rbind(tmetafulldf,tmetadf[[i]])
    }
  }
  tmetafulldf2 <- rbind(tmetafulldf,metafulldf)
  tmetafullndf <- normalization(tmetafulldf2,metafts)$df
  tmetafullndf <- tmetafullndf[1:nrow(tmetafulldf),]
  
  tgoffulldf <- NULL
  for (i in 1:length(tgofddmlist)) {
    if(is.null(tgoffulldf)){
      tgoffulldf <- tgofddmlist[[i]]
    }else{
      tgoffulldf <- rbind(tgoffulldf,tgofddmlist[[i]])
    }
  }
  tgoffulldf <- tgoffulldf[tgoffulldf$Model=="SVR",]
  
  stopifnot(nrow(tmetafullndf)==nrow(tgoffulldf))
  
  tfulldf <- cbind(tmetafullndf,tgoffulldf)
  
  #classification
  predretdf <- setting1predictcl(ddmeval$evaluator,ddmeval$fts,tfulldf)
  
  gofretlist <- ptgofdf
  if(is.null(gofretlist)){
    gofretlist <- execGOF111(goflist = goflist,gofregvec = paste0("N",goffts),cpert = cpert,
                             ptarget = pcri,fntrain = settingGLMtrain,fnpredict = settingBasicpredict,
                             tgoflist = tgoflist)
  }
  
  gofretfulldf <- NULL
  for (i in 1:length(gofretlist)) {
    if (is.null(gofretfulldf)){
      gofretfulldf <- gofretlist[[i]]$sel
    }else{
      gofretfulldf <- rbind(gofretfulldf,gofretlist[[i]]$sel)
    }
  }
  
  stopifnot(nrow(predretdf)==nrow(gofretfulldf))
  
  predretdf[paste0("DGS.RET.",pcri)] <- ifelse(predretdf$Pred=="DDM","DDM",gofretfulldf$Model)
  
  return(predretdf)
}



execDGSEvaluatorGen <- function(pgofldf=NULL,pddmodf=NULL,goflist,
                                gofnlist,gofddmlist,metadf,
                                pcri,totpert,goffts,metafts,
                                givenevalfile){
  
  goflret <- pgofldf
  if (is.null(goflret)){
    goflret <- setting111LOOnoout(goflist = goflist, gofregvec = paste0("N",goffts),
                                  totpert = totpert,target = paste0("N",pcri),
                                  measure = pcri,
                                  fntrain = settingGLMtrain, fnpredict = settingBasicpredict)
  }
  
  ddmret <- pddmodf
  if (is.null(ddmret)){
    ddmret <- settingDDMnoout(ngofddmlist = gofddmlist,
                              totpert = totpert,measures = c(pcri))
  }
  
  metafulldf <- NULL
  for (i in 1:length(metadf)) {
    if(is.null(metafulldf)){
      metafulldf <- metadf[[i]]
    }else{
      metafulldf <- rbind(metafulldf,metadf[[i]])
    }
  }
  metafullndf <- normalization(metafulldf,metafts)$df
  
  goffulldf <- NULL
  for (i in 1:length(gofddmlist)) {
    if(is.null(goffulldf)){
      goffulldf <- gofddmlist[[i]]
    }else{
      goffulldf <- rbind(goffulldf,gofddmlist[[i]])
    }
  }
  goffulldf <- goffulldf[goffulldf$Model=="SVR",]
  
  fulldf <- cbind(metafullndf,goffulldf)
  fullfeat <- c(goffts,paste0("N",metafts))
  
  ddmeval <- searchforEvaluator(plist = goflist,pnlist = gofnlist,pmetalist = metadf,
                                regvars = paste0("N",goffts),ppert = totpert,pcri = pcri,
                                pgofddmlist = gofddmlist, goffeats = goffts,metafeats = metafts,
                                pddmres = ddmret,pgofres = goflret,
                                pfulldf = fulldf,pfullfeats = fullfeat)
  
  print(ddmeval$fts)
  
  save(ddmeval,file = givenevalfile)
}
