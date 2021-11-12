settingOvOMCPert <- function(pdflist,pgoflist,pddmlist,pmetalist,
                             cpert,pcri,pmetafts,pmodelvec,
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
    winvec <- ifelse(phmetadf[,elem1]<phmetadf[,elem2],1,2)
    dfforthispair <- phmetadf[,pmetafts]
    # dfforthispair$Class <- as.factor(winvec)
    dfforthispair$Class <- winvec
    
    dfforthispair <- na.omit(dfforthispair)
    
    imagevec <- vector()
    for (i in 1:nrow(dfforthispair)) {
      imagevec <- c(imagevec,c(failuredataToImage(unlist(dfforthispair[i,pmetafts]),length(pmetafts))))
    }
    
    imagearr <- array(imagevec,dim = c(length(pmetafts),length(pmetafts),nrow(dfforthispair)))
    
    labelvec <- dfforthispair$Class
    
    # print(phmetadf[,elem1])
    # print(phmetadf[,elem2])
    # print(paste0(modelpairmat[i,1],'.',modelpairmat[i,2]))
    # print(dfforthispair)
    
    # clmodel <- modelConstruction(dfforthispair,pmetafts,'Class')
    
    clmodel <- modelConstructionImage(imagearr,labelvec,length(pmetafts))
    
    clmodellist[[paste0(modelpairmat[i,1],'.',modelpairmat[i,2])]] <- clmodel
  }
  
  targetimgarr <- failuredataToImage(unlist(ptargetmetadf[1,pmetafts]),length(pmetafts))
  
  # return(modelSelection(ptargetmetadf,clmodellist,pmodelvec,pcri,modelpairmat))
  return(modelSelectionImage(ptargetmetadf,targetimgarr,clmodellist,pmodelvec,pcri,modelpairmat))
}

modelConstruction <- function(pdataf,pfts,response){
  return(setting1traincl(pdataf,pfts,response))
}

modelConstructionImage <- function(pimgarr,plabvec,size){
  
  model <- keras_model_sequential()
  model %>%
    layer_flatten(input_shape = c(size, size)) %>%
    layer_dense(units = 128, activation = 'relu') %>%
    layer_dense(units = 10, activation = 'softmax')
  
  model %>% compile(
    optimizer = 'adam', 
    loss = 'sparse_categorical_crossentropy',
    metrics = c('accuracy')
  )
  
  model %>% fit(pimgarr, plabvec, epochs = 5, verbose = 2)
  
  return(model)
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
    
    predval <- pmodellist[[idxname]] %>% predict_classes(pimgarr)
    
    if (predval == 1){
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
