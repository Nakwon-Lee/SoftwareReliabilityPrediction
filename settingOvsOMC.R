modelgenAMETA <- function(pdflist){
  dppoints <- 30
  dpscales <- 30
  
  singleSRGMmodels <- c("GO","GG","Gompz","ISS","MD","MO","YID1","YID2",
                        "DSS","PNZ","PZ","PZI","Logi")
  
  load('ddmestsave.RData')
  stopifnot(exists('ddmestlist'))
  
  srgmestlist <- list()
  #ddmestlist <- list()
  predqlist <- list()
  historydplist <- list()
  for (i in 1:length(pdflist)) {
    trainhori <- floor((nrow(pdflist[[i]])*2)/3)
    predictionhori <- nrow(pdflist[[i]])
    historydplist[[i]] <- getDataPointsSingleLength(pdflist[[i]],trainhori,dppoints)
    srgmestlist[[i]] <- getSRGMEstOrig(dflist[[i]],singleSRGMmodels,trainhori,predictionhori)
    #ddmestlist[[i]] <- getDDMEstforSingle(df = dflist[[i]],isFCTBF = 'FC',models = c("SVR"))
    predqlist[[i]] <- getPredQOrigExt(dflist[[i]],srgmestlist[[i]],singleSRGMmodels,trainhori,predictionhori)
    predqlist[[i]] <- rbind(predqlist[[i]],getPredQOrigExt(dflist[[i]],ddmestlist[[i]],
                                                           c('SVR'),trainhori,predictionhori))
  }
  
  htempvec <- vector()
  for (i in 1:length(historydplist)) {
    htempvec <- c(htempvec,historydplist[[i]]$norm)
  }
  
  historydpmat <- matrix(htempvec,byrow = TRUE,ncol = dppoints)
  historydpmat <- round(historydpmat * dpscales)
  historymetadf <- data.frame(historydpmat)
  colnames(historymetadf) <- paste0('p',c(1:dppoints))
  
  for (i in 1:nrow(historymetadf)) {
    
    historycridf <- predqlist[[i]]
    
    modelstocompare <- historycridf$Model
    
    for (k in 1:length(modelstocompare)) {
      arow <- historycridf[which(historycridf$Model==modelstocompare[k]),]
      if (nrow(arow)==1){
        historymetadf[i,paste0('MEOP','.',modelstocompare[k])] <- arow[1,'MEOP']
      }else if(nrow(arow)>1){
        stop('more than one rows for a model')
      }else{
        historymetadf[i,paste0('MEOP','.',modelstocompare[k])] <- NA
      }
    }
  }
  
  modelstocompare2 <- vector()
  for (j in 1:length(modelstocompare)) {
    if(length(na.omit(historymetadf[,paste0('MEOP','.',modelstocompare[j])]))>=30){
      modelstocompare2[length(modelstocompare2)+1] <- modelstocompare[j]
    }
  }
  
  pmetafts <- paste0('p',c(1:dppoints))
  
  clmodellist <- list()
  
  modelpairmat <- combinations(length(modelstocompare2),2,modelstocompare2)
  for (i in 1:length(modelpairmat[,1])) {
    elem1 <- paste0('MEOP','.',modelpairmat[i,1])
    elem2 <- paste0('MEOP','.',modelpairmat[i,2])
    winvec <- ifelse(historymetadf[,elem1]<historymetadf[,elem2],0,1)
    dfforthispair <- historymetadf[,pmetafts]
    dfforthispair$Class <- winvec
    
    dfforthispair <- na.omit(dfforthispair)
    
    # Unix
    imagearr <- array(dim = c(nrow(dfforthispair),length(pmetafts),length(pmetafts)))
    for (j in 1:nrow(dfforthispair)) {
      imagearr[j,,] <- c(failuredataToImage(unlist(dfforthispair[j,pmetafts]),length(pmetafts)))
    }
    
    #Windows
    # imagearr <- array(dim = c(length(pmetafts),length(pmetafts)),nrow(dfforthispair))
    # for (j in 1:nrow(dfforthispair)) {
    #   imagearr[,,j] <- c(failuredataToImage(unlist(dfforthispair[j,pmetafts]),length(pmetafts)))
    # }
    
    labelvec <- dfforthispair$Class
    
    if (length(unique(labelvec))==1){
      if (unique(labelvec)[1]==0){
        clmodel <- modelpairmat[i,1]
      }else{
        clmodel <- modelpairmat[i,2]
      }
    }else{
      clmodel <- modelConstructionSVD21(imagearr,labelvec,pmetafts)
    }
    
    clmodellist[[paste0(modelpairmat[i,1],'.',modelpairmat[i,2])]] <- clmodel
  }
  
  clmodellist$modelpairmat <- modelpairmat
  
  return(clmodellist)
}

settingOvOMCPert <- function(pdflist,pgoflist,pddmlist,
                             cpert,pcri,pmodelvec,
                             pddmestlist,psrgmestlist){
  resultvec <- vector()
  binselvec <- vector()
  
  dppoints <- 30
  dpscales <- 30
  
  for(i in 1:length(pdflist)){
    picked <- c(i)
    
    targetdf <- pdflist[[picked]]
    
    tgofdf <- pgoflist[[picked]]
    tgofdf <- tgofdf[which(tgofdf$Pert==cpert),c('Model',pcri)]
    tddmdf <- pddmlist[[picked]]
    tddmdf <- tddmdf[which(tddmdf$Pert==cpert),c('Model',pcri)]
    tcridf <- rbind(tgofdf,tddmdf)
    
    modelstocompare <- tcridf$Model
    
    hori <- horiCalculator(cpert,nrow(targetdf))
    
    # histroymetalist <- list()
    # historydflist <- pdflist[-picked]
    # for (j in 1:length(historydflist)) {
    #   histroymetalist[[j]] <- getMETAInfoOrig(historydflist[[j]],pmetafts,hori$fit)
    #   stopifnot(nrow(histroymetalist[[j]])==1)
    #   histroymetalist[[j]]$Pert <- c(cpert)
    # }
    
    # histroymetalist <- pmetalist[-picked]
    # historymetadf <- listrbinder(histroymetalist)
    # historymetadf <- historymetadf[which(historymetadf$Pert==cpert),]
    
    historydflist <- pdflist[-picked]
    historydplist <- list()
    for (j in 1:length(historydflist)) {
      historydplist[[j]] <- getDataPointsSingle(historydflist[[j]],cpert,dppoints)
    }
    
    htempvec <- vector()
    for (j in 1:length(historydplist)) {
      htempvec <- c(htempvec,historydplist[[j]]$norm)
    }
    historydpmat <- matrix(htempvec,byrow = TRUE,ncol = dppoints)
    historydpmat <- round(historydpmat * dpscales)
    historymetadf <- data.frame(historydpmat)
    colnames(historymetadf) <- paste0('p',c(1:dppoints))
    
    historygoflist <- pgoflist[-picked]
    historyddmlist <- pddmlist[-picked]
    
    for (j in 1:nrow(historymetadf)) {
      
      # hgofdf <- historygoflist[[((j-1)%/%6)+1]]
      # hddmdf <- historyddmlist[[((j-1)%/%6)+1]]
      # hgofdf <- hgofdf[which(hgofdf$Pert==historymetadf[j,'Pert']),c('Model',pcri)]
      # hddmdf <- hddmdf[which(hddmdf$Pert==historymetadf[j,'Pert']),c('Model',pcri)]
      
      hgofdf <- historygoflist[[j]]
      hddmdf <- historyddmlist[[j]]
      hgofdf <- hgofdf[which(hgofdf$Pert==cpert),c('Model',pcri)]
      hddmdf <- hddmdf[which(hddmdf$Pert==cpert),c('Model',pcri)]
      
      historycridf <- rbind(hgofdf,hddmdf)
      
      for (k in 1:length(modelstocompare)) {
        arow <- historycridf[which(historycridf$Model==modelstocompare[k]),]
        if (nrow(arow)==1){
          historymetadf[j,paste0(pcri,'.',modelstocompare[k])] <- arow[1,pcri]
        }else if(nrow(arow)>1){
          stop('more than one rows for a model')
        }else{
          historymetadf[j,paste0(pcri,'.',modelstocompare[k])] <- NA
        }
      }
    }
    
    modelstocompare2 <- vector()
    for (j in 1:length(modelstocompare)) {
      if(length(na.omit(historymetadf[,paste0(pcri,'.',modelstocompare[j])]))>=30){
        modelstocompare2[length(modelstocompare2)+1] <- modelstocompare[j]
      }
    }
    
    # targetmetadf <- pmetalist[[picked]][which(pmetalist[[picked]]$Pert==cpert),]
    
    targetmetamat <- matrix(getDataPointsSingle(pdflist[[picked]],cpert,dppoints)$norm,byrow = TRUE,ncol = dppoints)
    targetmetamat <- round(targetmetamat * dpscales)
    targetmetadf <- data.frame(targetmetamat)
    colnames(targetmetadf) <- paste0('p',c(1:dppoints))
    
    retlist <- selectionOvOMC(pdflist[[picked]],paste0('p',c(1:dppoints)),modelstocompare2,
                                historymetadf,pcri,targetmetadf)
    
    idxvec <- vector()
    optvec <- vector()
    selvec <- vector()
    
    # print(modelstocompare2)
    
    for (j in 1:length(retlist$pairs[,1])) {
      elem1 <- retlist$pairs[j,1]
      elem2 <- retlist$pairs[j,2]
      idxvec[j] <- paste0(elem1,'.',elem2)
      
      if(tcridf[tcridf$Model==elem1,][1,pcri]<tcridf[tcridf$Model==elem2,][1,pcri]){
        optvec[j] <- elem1
      }else{
        optvec[j] <- elem2
      }
      selvec[j] <- retlist$binary[[idxvec[j]]]
    }
    
    selectionQ <- data.frame(idx=idxvec,opt=optvec,sel=selvec)
    
    # print(nrow(selectionQ[selectionQ$opt==selectionQ$sel,]))
    # print(combn(modelstocompare2,2))
    
    crival <- NA
    key <- 0
    while(is.na(crival)){
      key <- key + 1
      tempgofdf <- pgoflist[[picked]]
      tempgofdf <- tempgofdf[which(tempgofdf$Pert==cpert & tempgofdf$Model==names(retlist$selected)[key]),]
      crival <- tempgofdf[1,pcri]
      # # estimated <- predictionUsingSelectedModel(targetdf,names(retlist$selected)[key],hori)
      # estimated <- predictionUsingSelectedModelPre(names(retlist$selected)[key],pddmestlist[[1]][[picked]],psrgmestlist[[1]][[picked]])
      # if (pcri=='EP'){
      #   crival <- FC.EP(hori$pred,targetdf$n,estimated$EstElap)
      # }else if (pcri=='MEOP'){
      #   crival <- FC.MEOP(hori$fit,hori$pred,targetdf$n,estimated$EstElap)
      # }
    }
    
    resultvec[(length(resultvec)+1)] <- crival
    binselvec[(length(binselvec)+1)] <- nrow(selectionQ[selectionQ$opt==selectionQ$sel,])/ncol(combn(modelstocompare2,2))
  }
  
  retdf <- data.frame(resultvec)
  colnames(retdf) <- c(pcri)
  
  retdf$binsel <- binselvec
  
  print(retdf)
  print(c(pcri,median(retdf[,pcri])))
  
  return(retdf)
}

selectionOvOMC <- function(targetdf,pmetafts,pmodelvec,phmetadf,pcri,ptargetmetadf){
  
  clmodellist <- list()
  
  modelpairmat <- combinations(length(pmodelvec),2,pmodelvec)
  for (i in 1:length(modelpairmat[,1])) {
    elem1 <- paste0(pcri,'.',modelpairmat[i,1])
    elem2 <- paste0(pcri,'.',modelpairmat[i,2])
    # winvec <- ifelse(phmetadf[,elem1]<phmetadf[,elem2],modelpairmat[i,1],modelpairmat[i,2])
    winvec <- ifelse(phmetadf[,elem1]<phmetadf[,elem2],0,1)
    dfforthispair <- phmetadf[,pmetafts]
    # dfforthispair$Class <- as.factor(winvec)
    dfforthispair$Class <- winvec
    
    dfforthispair <- na.omit(dfforthispair)
    
    # Unix
    imagearr <- array(dim = c(nrow(dfforthispair),length(pmetafts),length(pmetafts)))
    for (j in 1:nrow(dfforthispair)) {
      imagearr[j,,] <- c(failuredataToImage(unlist(dfforthispair[j,pmetafts]),length(pmetafts)))
    }
    
    #Windows
    # imagearr <- array(dim = c(length(pmetafts),length(pmetafts)),nrow(dfforthispair))
    # for (j in 1:nrow(dfforthispair)) {
    #   imagearr[,,j] <- c(failuredataToImage(unlist(dfforthispair[j,pmetafts]),length(pmetafts)))
    # }
    
    labelvec <- dfforthispair$Class
    
    # print(phmetadf[,elem1])
    # print(phmetadf[,elem2])
    # print(paste0(modelpairmat[i,1],'.',modelpairmat[i,2]))
    # print(dfforthispair)
    
    # clmodel <- modelConstruction(dfforthispair,pmetafts,'Class')
    
    # print(imagearr[1,,])
    # print(dim(imagearr))
    # print(labelvec)
    # print(dim(labelvec))
    
    # clmodel <- modelConstructionImage(imagearr,labelvec,length(pmetafts))
    
    # clmodel <- modelConstructionSVD(dfforthispair,pmetafts)
    
    clmodel <- modelConstructionSVD21(imagearr,labelvec,pmetafts)
    
    clmodellist[[paste0(modelpairmat[i,1],'.',modelpairmat[i,2])]] <- clmodel
  }
  
  # Unix
  targetimgarr <- array(dim = c(1,length(pmetafts),length(pmetafts)))
  targetimgarr[1,,] <- failuredataToImage(unlist(ptargetmetadf[1,pmetafts]),length(pmetafts))
  
  # Windows
  # targetimgarr <- array(dim = c(length(pmetafts),length(pmetafts),1))
  # targetimgarr[,,1] <- failuredataToImage(unlist(ptargetmetadf[1,pmetafts]),length(pmetafts))
  
  # return(modelSelection(ptargetmetadf,clmodellist,pmodelvec,pcri,modelpairmat))
  # return(modelSelectionImage(ptargetmetadf,targetimgarr,clmodellist,pmodelvec,pcri,modelpairmat))
  # return(modelSelectionSVD(ptargetmetadf,clmodellist,pmodelvec,pcri,modelpairmat,pmetafts))
  return(modelSelectionSVD21(ptargetmetadf,targetimgarr,clmodellist,pmodelvec,pcri,modelpairmat,pmetafts))
}

modelConstruction <- function(pdataf,pfts,response){
  return(setting1traincl(pdataf,pfts,response))
}

modelConstructionImage <- function(pimgarr,plabvec,size){
  
  model <- keras_model_sequential()
  model %>%
    layer_flatten(input_shape = c(size, size)) %>%
    layer_dense(units = 64, activation = 'relu') %>%
    layer_dense(units = 2, activation = 'softmax')
  
  model %>% compile(
    optimizer = 'adam', 
    loss = 'sparse_categorical_crossentropy',
    metrics = c('accuracy')
  )
  
  model %>% fit(pimgarr, plabvec,
                validation_split = 0.2,
                epochs = 10, verbose = 0)
  
  return(model)
}

modelConstructionSVD <- function(pdatap, pfts){
  alphamatrics <- list()
  for (i in 1:2) {
    alphamatrics[[i]] <- t(data.matrix(pdatap[which(pdatap$Class==(i-1)),pfts]))
  }
  
  ret <- list()
  
  for (i in 1:2) {
    svdret <- svd(alphamatrics[[i]])
    
    ret[[i]] <- list()
    
    ret[[i]]$leftsingular <- svdret$u
    ret[[i]]$rightsingular <- svdret$v
    ret[[i]]$singularmatrix <- svdret$d
  }
  
  # validation on the training set (eliminating some left singular vectors)
  
  for (k in 1:2) {
    leftsingularqvec <- vector()
    for (i in 1:ncol(ret[[k]]$leftsingular)) {
      classaprojmean <- mean(t(ret[[k]]$leftsingular[,i]) %*% alphamatrics[[1]])
      classbprojmean <- mean(t(ret[[k]]$leftsingular[,i]) %*% alphamatrics[[2]])
      
      leftsingularqvec[i] <- abs(classaprojmean - classbprojmean)
    }
    
    ordervec <- order(leftsingularqvec,decreasing = TRUE)
    ret[[k]]$leftsingular <- ret[[k]]$leftsingular[,ordervec[1:max(1,floor(length(ordervec)*0.2))],drop=FALSE]
  }
  
  return(ret)
}

modelConstructionSVD1 <- function(pdatap, pfts){
  alphamatrics <- list()
  for (i in 1:2) {
    alphamatrics[[i]] <- t(data.matrix(pdatap[which(pdatap$Class==(i-1)),pfts]))
  }
  
  ret <- list()
  
  for (i in 1:2) {
    svdret <- svd(alphamatrics[[i]])
    
    ret[[i]] <- list()
    
    ret[[i]]$leftsingular <- svdret$u
    ret[[i]]$rightsingular <- svdret$v
    ret[[i]]$singularmatrix <- svdret$d
  }
  
  # validation on the training set (eliminating some left singular vectors)
  
  for (k in 1:2) {
    leftsingularqvec <- vector()
    for (i in 1:ncol(ret[[k]]$leftsingular)) {
      classaprojmean <- mean(t(ret[[k]]$leftsingular[,i]) %*% alphamatrics[[1]])
      classbprojmean <- mean(t(ret[[k]]$leftsingular[,i]) %*% alphamatrics[[2]])
      
      leftsingularqvec[i] <- abs(classaprojmean - classbprojmean)
    }
    
    ordervec <- order(leftsingularqvec,decreasing = TRUE)
    ret[[k]]$leftsingular <- ret[[k]]$leftsingular[,ordervec[1:max(1,floor(length(ordervec)*0.2))],drop=FALSE]
  }
  
  u <- cbind(ret[[1]]$leftsingular,ret[[2]]$leftsingular)
  
  trainingset <- data.frame(t(t(u) %*% t(data.matrix(pdatap[,pfts]))))
  colnames(trainingset) <- c(paste0('A',c(1:ncol(ret[[1]]$leftsingular))),paste0('B',c(1:ncol(ret[[2]]$leftsingular))))
  trainingset$Class <- pdatap$Class
  
  ret$model <- setting1traincl(trainingset,c(paste0('A',c(1:ncol(ret[[1]]$leftsingular))),paste0('B',c(1:ncol(ret[[2]]$leftsingular)))),'Class')
  
  return(ret)
}

modelConstructionSVD2 <- function(pimagearr, plabelvec, pfts){
  
  pdatap <- matrix(pimagearr,ncol = (length(pfts) * length(pfts)), nrow = dim(pimagearr)[1])
  pdatap <- t(pdatap)
  
  alphamatrics <- list()
  for (i in 1:2) {
    alphamatrics[[i]] <- pdatap[,which(plabelvec==(i-1))]
  }
  
  ret <- list()
  
  for (i in 1:2) {
    svdret <- svd(alphamatrics[[i]])
    
    ret[[i]] <- list()
    
    ret[[i]]$leftsingular <- svdret$u
    ret[[i]]$rightsingular <- svdret$v
    ret[[i]]$singularmatrix <- svdret$d
  }
  
  # validation on the training set (eliminating some left singular vectors)
  
  for (k in 1:2) {
    leftsingularqvec <- vector()
    for (i in 1:ncol(ret[[k]]$leftsingular)) {
      classaprojmean <- mean(t(ret[[k]]$leftsingular[,i]) %*% alphamatrics[[1]])
      classbprojmean <- mean(t(ret[[k]]$leftsingular[,i]) %*% alphamatrics[[2]])
      
      leftsingularqvec[i] <- abs(classaprojmean - classbprojmean)
    }
    
    ordervec <- order(leftsingularqvec,decreasing = TRUE)
    ret[[k]]$leftsingular <- ret[[k]]$leftsingular[,ordervec[1:max(1,floor(length(ordervec)*0.2))],drop=FALSE]
  }
  
  return(ret)
}

modelConstructionSVD21 <- function(pimagearr, plabelvec, pfts){
  
  pdatap <- matrix(pimagearr,ncol = (length(pfts) * length(pfts)), nrow = dim(pimagearr)[1])
  pdatap <- t(pdatap)
  
  alphamatrics <- list()
  for (i in 1:2) {
    alphamatrics[[i]] <- pdatap[,which(plabelvec==(i-1))]
  }
  
  ret <- list()
  
  for (i in 1:2) {
    
    svdret <- svd(alphamatrics[[i]])
    
    ret[[i]] <- list()
    
    ret[[i]]$leftsingular <- svdret$u
    ret[[i]]$rightsingular <- svdret$v
    ret[[i]]$singularmatrix <- svdret$d
  }
  
  # validation on the training set (eliminating some left singular vectors)
  
  for (k in 1:2) {
    leftsingularqvec <- vector()
    for (i in 1:ncol(ret[[k]]$leftsingular)) {
      classaprojmean <- mean(t(ret[[k]]$leftsingular[,i]) %*% alphamatrics[[1]])
      classbprojmean <- mean(t(ret[[k]]$leftsingular[,i]) %*% alphamatrics[[2]])
      
      leftsingularqvec[i] <- abs(classaprojmean - classbprojmean)
    }
    
    ordervec <- order(leftsingularqvec,decreasing = TRUE)
    ret[[k]]$leftsingular <- ret[[k]]$leftsingular[,ordervec[1:max(1,floor(length(ordervec)*0.2))],drop=FALSE]
  }
  
  u <- cbind(ret[[1]]$leftsingular,ret[[2]]$leftsingular)
  
  trainingset <- data.frame(t(t(u) %*% pdatap))
  colnames(trainingset) <- c(paste0('A',c(1:ncol(ret[[1]]$leftsingular))),paste0('B',c(1:ncol(ret[[2]]$leftsingular))))
  trainingset$Class <- plabelvec
  
  ret$model <- setting1traincl(trainingset,c(paste0('A',c(1:ncol(ret[[1]]$leftsingular))),paste0('B',c(1:ncol(ret[[2]]$leftsingular)))),'Class')
  
  return(ret)
}

modelConstructionSVD3 <- function(pimagearr, plabelvec, pfts){
  
  pdatap <- matrix(pimagearr,ncol = (length(pfts) * length(pfts)), nrow = dim(pimagearr)[1])
  pdatap <- t(pdatap)
  
  svdret <- svd(pdatap)
  
  
  
  for (i in 1:nrow(svdret$v)) {
    
  }
  
  alphamatrics <- list()
  for (i in 1:2) {
    alphamatrics[[i]] <- pdatap[,which(plabelvec==(i-1))]
  }
  
  ret <- list()
  
  for (i in 1:2) {
    svdret <- svd(alphamatrics[[i]])
    
    ret[[i]] <- list()
    
    ret[[i]]$leftsingular <- svdret$u
    ret[[i]]$rightsingular <- svdret$v
    ret[[i]]$singularmatrix <- svdret$d
  }
  
  return(ret)
}

modelSelection <- function(pdataf,pmodellist,pmodelvec,pcri,modelpairs){
  #evaluate all combinations
  
  retlist <- list()
  retlist$binary <- list()
  
  wincountvec <- rep(0,length(pmodelvec))
  for (i in 1:length(modelpairs[,1])) {
    elem1 <- modelpairs[i,1]
    elem2 <- modelpairs[i,2]
    
    idxname <- paste0(elem1,'.',elem2)
    
    predval <- predict(pmodellist[[idxname]],pdataf,type='class')
    
    if (predval == elem1){
      wincountvec[which(pmodelvec==elem1)[1]] <- wincountvec[which(pmodelvec==elem1)[1]] + 1
      retlist$binary[[idxname]] <- elem1
    }else{
      wincountvec[which(pmodelvec==elem2)[1]] <- wincountvec[which(pmodelvec==elem2)[1]] + 1
      retlist$binary[[idxname]]<- elem2
    }
  }
  
  names(wincountvec) <- pmodelvec
  wincountvec <-sort(wincountvec,decreasing = TRUE)
  
  selected <- wincountvec
  
  retlist$selected <- selected
  retlist$pairs <- modelpairs
  
  return(retlist)
}

modelSelectionImage <- function(pdataf,pimgarr,pmodellist,pmodelvec,pcri,modelpairs){
  retlist <- list()
  retlist$binary <- list()
  
  wincountvec <- rep(0,length(pmodelvec))
  for (i in 1:length(modelpairs[,1])) {
    elem1 <- modelpairs[i,1]
    elem2 <- modelpairs[i,2]
    
    idxname <- paste0(elem1,'.',elem2)
    
    predicted <- pmodellist[[idxname]] %>% predict(pimgarr)
    
    predval <- which.max(predicted[1,])
    
    if (predval == 1){
      wincountvec[which(pmodelvec==elem1)[1]] <- wincountvec[which(pmodelvec==elem1)[1]] + 1
      retlist$binary[[idxname]] <- elem1
    }else if (predval == 2){
      wincountvec[which(pmodelvec==elem2)[1]] <- wincountvec[which(pmodelvec==elem2)[1]] + 1
      retlist$binary[[idxname]]<- elem2
    }else{
      stop("no valid predicted class")
    }
  }
  
  names(wincountvec) <- pmodelvec
  wincountvec <-sort(wincountvec,decreasing = TRUE)
  
  selected <- wincountvec
  
  retlist$selected <- selected
  retlist$pairs <- modelpairs
  
  return(retlist)
}

modelSelectionSVD <- function(pdataf,pmodellist,pmodelvec,pcri,modelpairs,pmetafts){
  retlist <- list()
  retlist$binary <- list()
  
  wincountvec <- rep(0,length(pmodelvec))
  for (i in 1:length(modelpairs[,1])) {
    elem1 <- modelpairs[i,1]
    elem2 <- modelpairs[i,2]
    
    idxname <- paste0(elem1,'.',elem2)
    
    predval <- predictSVD(pmodellist[[idxname]],pmetafts,pdataf)
    
    if (predval == 1){
      wincountvec[which(pmodelvec==elem1)[1]] <- wincountvec[which(pmodelvec==elem1)[1]] + 1
      retlist$binary[[idxname]] <- elem1
    }else if (predval == 2){
      wincountvec[which(pmodelvec==elem2)[1]] <- wincountvec[which(pmodelvec==elem2)[1]] + 1
      retlist$binary[[idxname]]<- elem2
    }else{
      stop("no valid predicted class")
    }
  }
  
  names(wincountvec) <- pmodelvec
  wincountvec <-sort(wincountvec,decreasing = TRUE)
  
  selected <- wincountvec
  
  retlist$selected <- selected
  retlist$pairs <- modelpairs
  
  return(retlist)
}

predictSVD <- function(psvdret,pmetafts,pdataf){
  targetmat <- matrix(unlist(pdataf[1,pmetafts]),ncol = 1 ,nrow = length(pmetafts))
  residvec <- vector()
  for (i in 1:length(psvdret)) {
    # k <- min(10,ncol(psvdret[[i]]$leftsingular))
    u <- psvdret[[i]]$leftsingular
    residvec[i] <- norm((diag(nrow(targetmat)) - (u %*% t(u))) %*% targetmat)
  }
  return(which.min(residvec))
}


modelSelectionSVD1 <- function(pdataf,pmodellist,pmodelvec,pcri,modelpairs,pmetafts){
  retlist <- list()
  retlist$binary <- list()
  
  wincountvec <- rep(0,length(pmodelvec))
  for (i in 1:length(modelpairs[,1])) {
    elem1 <- modelpairs[i,1]
    elem2 <- modelpairs[i,2]
    
    idxname <- paste0(elem1,'.',elem2)
    
    predval <- predictSVD1(pmodellist[[idxname]],pmetafts,pdataf)
    
    if (predval == 1){
      wincountvec[which(pmodelvec==elem1)[1]] <- wincountvec[which(pmodelvec==elem1)[1]] + 1
      retlist$binary[[idxname]] <- elem1
    }else if (predval == 2){
      wincountvec[which(pmodelvec==elem2)[1]] <- wincountvec[which(pmodelvec==elem2)[1]] + 1
      retlist$binary[[idxname]]<- elem2
    }else{
      stop("no valid predicted class")
    }
  }
  
  names(wincountvec) <- pmodelvec
  wincountvec <-sort(wincountvec,decreasing = TRUE)
  
  selected <- wincountvec
  
  retlist$selected <- selected
  retlist$pairs <- modelpairs
  
  return(retlist)
}

predictSVD1 <- function(psvdret,pmetafts,pdataf){
  targetmat <- matrix(unlist(pdataf[1,pmetafts]),ncol = 1 ,nrow = length(pmetafts))
  
  u <- cbind(psvdret[[1]]$leftsingular,psvdret[[2]]$leftsingular)
  
  testset <- data.frame(t(t(u) %*% targetmat))
  colnames(testset) <- c(paste0('A',c(1:ncol(psvdret[[1]]$leftsingular))),paste0('B',c(1:ncol(psvdret[[2]]$leftsingular))))
  
  return(ifelse(predict(psvdret$model,testset)<0.5,1,2))
}

modelSelectionSVD2 <- function(pdataf,ptargetimgarr,pmodellist,pmodelvec,pcri,modelpairs,pmetafts){
  retlist <- list()
  retlist$binary <- list()
  
  wincountvec <- rep(0,length(pmodelvec))
  for (i in 1:length(modelpairs[,1])) {
    elem1 <- modelpairs[i,1]
    elem2 <- modelpairs[i,2]
    
    idxname <- paste0(elem1,'.',elem2)
    
    predval <- predictSVD2(pmodellist[[idxname]],pmetafts,ptargetimgarr)
    
    if (predval == 1){
      wincountvec[which(pmodelvec==elem1)[1]] <- wincountvec[which(pmodelvec==elem1)[1]] + 1
      retlist$binary[[idxname]] <- elem1
    }else if (predval == 2){
      wincountvec[which(pmodelvec==elem2)[1]] <- wincountvec[which(pmodelvec==elem2)[1]] + 1
      retlist$binary[[idxname]]<- elem2
    }else{
      stop("no valid predicted class")
    }
  }
  
  names(wincountvec) <- pmodelvec
  wincountvec <-sort(wincountvec,decreasing = TRUE)
  
  selected <- wincountvec
  
  retlist$selected <- selected
  retlist$pairs <- modelpairs
  
  return(retlist)
}

predictSVD2 <- function(psvdret,pmetafts,ptargetimgarr){
  targetmat <- matrix(c(ptargetimgarr[1,,]),ncol = 1 ,nrow = length(pmetafts) * length(pmetafts))
  residvec <- vector()
  for (i in 1:length(psvdret)) {
    # k <- min(5,ncol(psvdret[[i]]$leftsingular))
    u <- psvdret[[i]]$leftsingular
    
    residvec[i] <- binaryimageresid(targetmat,u)
    
    # residvec[i] <- norm((diag(nrow(targetmat)) - (u %*% t(u))) %*% targetmat)
  }
  return(which.min(residvec))
}

modelSelectionSVD21 <- function(pdataf,ptargetimgarr,pmodellist,pmodelvec,pcri,modelpairs,pmetafts){
  retlist <- list()
  retlist$binary <- list()
  
  wincountvec <- rep(0,length(pmodelvec))
  for (i in 1:length(modelpairs[,1])) {
    elem1 <- modelpairs[i,1]
    elem2 <- modelpairs[i,2]
    
    idxname <- paste0(elem1,'.',elem2)
    
    if (typeof(pmodellist[[idxname]])=='character'){
      if(pmodellist[[idxname]]==elem1){
        predval = 1
      }else{
        predval = 2
      }
    }else{
      predval <- predictSVD21(pmodellist[[idxname]],pmetafts,ptargetimgarr)
    }
    
    if (predval == 1){
      wincountvec[which(pmodelvec==elem1)[1]] <- wincountvec[which(pmodelvec==elem1)[1]] + 1
      retlist$binary[[idxname]] <- elem1
    }else if (predval == 2){
      wincountvec[which(pmodelvec==elem2)[1]] <- wincountvec[which(pmodelvec==elem2)[1]] + 1
      retlist$binary[[idxname]]<- elem2
    }else{
      stop("no valid predicted class")
    }
  }
  
  names(wincountvec) <- pmodelvec
  wincountvec <-sort(wincountvec,decreasing = TRUE)
  
  selected <- wincountvec
  
  retlist$selected <- selected
  retlist$pairs <- modelpairs
  
  return(retlist)
}

predictSVD21 <- function(psvdret,pmetafts,ptargetimgarr){
  
  targetmat <- matrix(c(ptargetimgarr[1,,]),ncol = 1 ,nrow = length(pmetafts) * length(pmetafts))
  
  u <- cbind(psvdret[[1]]$leftsingular,psvdret[[2]]$leftsingular)
  
  testset <- data.frame(t(t(u) %*% targetmat))
  colnames(testset) <- c(paste0('A',c(1:ncol(psvdret[[1]]$leftsingular))),paste0('B',c(1:ncol(psvdret[[2]]$leftsingular))))
  
  return(ifelse(predict(psvdret$model,testset)<0.5,1,2))
}

binaryimageresid <- function(ptargetmat,pu){
  x <- t(pu) %*% ptargetmat
  ux <- pu %*% x
  ux <- ifelse(ux<0.5,0,1)
  return(sum(abs(ptargetmat - ux)))
}

predictionUsingSelectedModel <- function(pdataf,selected,hori){
  return(getSRGMEstOrigSingle(selected,
                       pdataf$n,
                       hori$fit,
                       hori$pred)$estimated)
}

predictionUsingSelectedModelPre <- function(selected,ddmest,srgmest){
  ret <- NULL
  if(length(which(names(ddmest)==selected))>0){
    ret <- ddmest[[selected]]
  }else if(length(which(names(srgmest)==selected))>0){
    ret <- srgmest[[selected]]
  }else{
    stop('no matching model')
  }
  return(ret)
}

failuredataToImage <- function(pdp,size){
  stopifnot(length(pdp)==size)
  
  basevec <- vector()
  for (i in 1:size) {
    basevec <- c(basevec,rep(0,times = (size-pdp[i])),rep(1,times = pdp[i]))
  }
  
  return(matrix(data = basevec,nrow = size,ncol = size))
}
