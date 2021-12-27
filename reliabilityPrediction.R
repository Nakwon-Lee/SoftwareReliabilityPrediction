reliabilityPrediction <- function(ddf,dmodels){
  singleSRGMmodels <- c("GO","GG","Gompz","ISS","MD","MO","YID1","YID2","DSS","PNZ","PZ","PZI","Logi")
  trainhori <- nrow(ddf)
  predictionhori <- floor(trainhori*1.5)
  singleestlist <- getSRGMEstOrig(ddf,singleSRGMmodels,trainhori,predictionhori)
  # singleddmest <- getDDMEstforSingle(ddf,c('SVR'),'FC')
  load('targetddmest.RData')
  singleestlist <- c(singleestlist,singleddmest)
  
  retlist <- singleestlist
  
  if ('GOFC' %in% names(dmodels)){
    goffts <- c("MSE","MAE","Rsquare","Noise","Bias","Variation","PRR","WLSE","CEP","CMEOP")
    gofdf <- getSRGMGOFOrigExt(ddf,singleestlist,singleSRGMmodels,goffts,trainhori,predictionhori)
    gofdf$Pred <- predict(dmodels$GOFC,gofdf)
    
    ordered.pred <- gofdf[order(gofdf$Pred),]
    
    retlist$GOFC <- singleestlist[[ordered.pred[1,'Model']]]
  }
  if ('META' %in% names(dmodels)){
    metafts <- c("Variance","Inclination","AutoCorr","MetaNO","NumP","LapFact","SubAddi")
    metadf <- getMETAInfoOrig(ddf,metafts,trainhori)
    retlist$META <- singleestlist[[predict(dmodels$META,metadf,type = 'class')]]
  }
  if ('AMETA' %in% names(dmodels)){
    dppoints <- 30
    dpscales <- 30
    pmetafts <- paste0('p',c(1:dppoints))
    targetmetamat <- matrix(getDataPointsSingleLength(ddf,trainhori,dppoints)$norm,byrow = TRUE,ncol = dppoints)
    targetmetamat <- round(targetmetamat * dpscales)
    targetmetadf <- data.frame(targetmetamat)
    colnames(targetmetadf) <- pmetafts
    
    targetimgarr <- array(dim = c(1,length(pmetafts),length(pmetafts)))
    targetimgarr[1,,] <- failuredataToImage(unlist(targetmetadf[1,pmetafts]),length(pmetafts))
    
    selretlist <- modelSelectionSVD21(targetmetadf,targetimgarr,dmodels$AMETA,c(singleSRGMmodels,'SVR'),
                               'MEOP',dmodels$AMETA$modelpairmat,pmetafts)
    
    crival <- NA
    key <- 0
    while(is.na(crival)){
      key <- key + 1
      estimated <- singleestlist[[names(selretlist$selected)[key]]]
      retlist$AMETA <- estimated
      crival <- FC.MEOP(1,trainhori,ddf$n,estimated$EstElap)
    }
  }
  
  return(retlist)
}