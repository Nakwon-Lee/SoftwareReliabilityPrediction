settingRESEDAPert <- function(pdflist,pgoffts,pmetafts,
                              pcri,pmodelvec,
                              pgaparam,cpert,
                              pddmestlist,
                              psrgmestlist,
                              pmodelforddm,
                              evalGen,preDict,
                              histdf=NULL,
                              evalfts){
  
  resultvec <- vector()
  selvec <- vector()
  
  for (i in 1:length(pdflist)) {
    cat(i)
    
    picked <- c(i)
    
    targetdf <- pdflist[[picked]]
    
    evaluator <- evalGen(pdflist = pdflist[-picked],pgoffts = pgoffts,
                                    pmetafts = pmetafts,
                                    ppert = cpert,
                                    pcri = pcri,
                                    pmodelvec = pmodelvec,
                                    pddmmodel = pmodelforddm,pgaparam = pgaparam,
                                    pddmestlist = pddmestlist[-picked],
                                    psrgmestlist = psrgmestlist[-picked],
                         evalfts=evalfts)
    
    hori <- horiCalculator(cpert,nrow(targetdf))
    
    predresult <- preDict(pevaluator = evaluator,ptargetdf = targetdf[1:hori$fit,],
                                  ppert = cpert,
                                  pmodelvec = pmodelvec,pdflist = pdflist[-picked],
                                  pcri = pcri,pgoffts = pgoffts,
                                  pmetafts = pmetafts,pddmmodel = pmodelforddm,
                                  pfith = hori$fit,ppredh = hori$pred,
                                  psrgmestlist = psrgmestlist[-picked],
                                  ptargetddmest = pddmestlist[[picked]],
                                  ptargetsrgmest = psrgmestlist[[picked]])
    estimated <- predresult$estimated
    
    if (pcri=='EP'){
      crival <- FC.EP(hori$pred,targetdf$n,estimated$EstElap)
    }else if (pcri=='MEOP'){
      crival <- FC.MEOP(hori$fit,hori$pred,targetdf$n,estimated$EstElap)
    }
    
    stopifnot(!is.na(crival))
    
    resultvec[(length(resultvec)+1)] <- crival
    selvec[(length(selvec)+1)] <- predresult$sel
  }
  
  if (!is.null(histdf)){
    histdfsub <- histdf[which(histdf$Pert==cpert),]
    
    stopifnot(nrow(histdfsub)==length(selvec))
    
    colsel <- paste0(pcri,'.SEL')
    
    totcount <- nrow(histdfsub)
    gofcount <- length(which(histdfsub[,colsel]=='GOF'))
    ddmcount <- length(which(histdfsub[,colsel]=='DDM'))
    
    crtcount <- length(which(histdfsub[,colsel]==selvec))
    crcntgof <- length(which(histdfsub[,colsel]==selvec & histdfsub[,colsel]=='GOF'))
    crcntddm <- length(which(histdfsub[,colsel]==selvec & histdfsub[,colsel]=='DDM'))

    print(c('acc:',(crtcount/totcount)))
    print(c('acc.gof:',(crcntgof/gofcount)))
    print(c('acc.ddm:',(crcntddm/ddmcount)))
  }
  
  retdf <- data.frame(resultvec)
  colnames(retdf) <- c(pcri)
  
  print(c(pcri,median(retdf[,pcri])))
  
  return(retdf)
}

settingRESEDAv2Pert <- function(pdflist,pddmmodel,ppert,pgoffts,pcri,pmodelvec,pddmestlist=NULL){
  
  datadrivenest <- pddmestlist
  if (is.null(datadrivenest)){
    datadrivenest <- dataDriven(pdflist,pddmmodel,ppert)
  }
  print('end datadriven')
  
  resultvec <- vector()
  
  for (i in 1:length(pdflist)) {
    
    print(i)
    
    targetdf <- pdflist[[i]]
    
    hori <- horiCalculator(ppert,nrow(targetdf))
    
    estimated <- predictionRESEDAv2(targetdf,pddmmodel,
                                    hori$fit,hori$pred,
                                    pmodelvec,pgoffts,
                                    ptargetddmest=datadrivenest[[i]])
    
    print('get RESEDA result')
    
    if (pcri=='EP'){
      crival <- FC.EP(hori$pred,targetdf$n,estimated$EstElap)
    }else if (pcri=='MEOP'){
      crival <- FC.MEOP(hori$fit,hori$pred,targetdf$n,estimated$EstElap)
    }
    
    stopifnot(!is.na(crival))
    
    resultvec[(length(resultvec)+1)] <- crival
  }
  
  retdf <- data.frame(resultvec)
  colnames(retdf) <- c(pcri)
  
  print(c(pcri,median(retdf[,pcri])))
  
  return(retdf)
}

# N data
evaluatorGenRESEDA <- function(pdflist,pgoffts,pmetafts,ppert,
                               pcri,pmodelvec,pddmmodel,
                               pgaparam,
                               pddmestlist=NULL,
                               psrgmestlist=NULL){
  
  cmodelsest <- psrgmestlist
  if (is.null(cmodelsest)){
    cmodelsest <- modelFitting(pdflist,pmodelvec,ppert)
  }
  print('end modelfit')
  
  cgofnlist <- calcGOF(pdflist,cmodelsest,pmodelvec,pgoffts,ppert)
  print('end calcGOF')
  
  mselectedlist <- modelSelection(cgofnlist,paste0('N',pgoffts),paste0('N',pcri))
  modelselres <- list()
  for (i in 1:length(mselectedlist)) {
    modelselres[[i]] <- cgofnlist[[i]][which(cgofnlist[[i]][,'Model']==mselectedlist[[i]]),]
  }
  print('end model selection')
  
  datadrivenest <- pddmestlist
  if (is.null(datadrivenest)){
    datadrivenest <- dataDriven(pdflist,pddmmodel,ppert)
  }
  print('end datadriven')
  
  datadrivenres <- calcGOFDDM(datadrivenest,pdflist,ppert,pgoffts)
  print('end calcGOFDDM')
  
  metares <- calcMeta(pdflist,pmetafts,ppert)
  print('end calcMeta')
  
  histdatalist <- list()
  for (i in 1:length(modelselres)) {
    tempdf <- cbind(datadrivenres[[i]],metares[[i]])
    tempdf[1,'SEL'] <- ifelse(modelselres[[i]][1,pcri]>=datadrivenres[[i]][1,pcri],'DDM','GOF')
    histdatalist[[i]] <- tempdf
  }
  print('end get hist data')
  
  searchedfeats <- searchForFeatures(histdatalist,c(pgoffts,pmetafts),pgaparam)
  print('end search feats')
  
  evaluator <- settingDGStraincl(listrbinder(histdatalist),searchedfeats,'SEL')
  
  return(evaluator)
}

evaluatorGenRESEDAv2 <- function(pdflist,pgoffts,pmetafts,ppert,
                               pcri,pmodelvec,pddmmodel,
                               pgaparam,
                               pddmestlist=NULL,
                               psrgmestlist=NULL,
                               evalfts){
  
  cmodelsest <- psrgmestlist
  if (is.null(cmodelsest)){
    cmodelsest <- modelFitting(pdflist,pmodelvec,ppert)
  }
  # print('end modelfit')
  
  cgofnlist <- calcGOF(pdflist,cmodelsest,pmodelvec,pgoffts,ppert)
  # print('end calcGOF')
  
  mselectedlist <- modelSelection(cgofnlist,paste0('N',pgoffts),paste0('N',pcri))
  modelselres <- list()
  for (i in 1:length(mselectedlist)) {
    modelselres[[i]] <- cgofnlist[[i]][which(cgofnlist[[i]][,'Model']==mselectedlist[[i]]),]
  }
  # print('end model selection')
  
  datadrivenest <- pddmestlist
  if (is.null(datadrivenest)){
    datadrivenest <- dataDriven(pdflist,pddmmodel,ppert)
  }
  # print('end datadriven')
  
  datadrivenres <- calcGOFDDM(datadrivenest,pdflist,ppert,pgoffts)
  # print('end calcGOFDDM')
  
  metares <- calcMeta(pdflist,pmetafts,ppert)
  # print('end calcMeta')
  
  histdatalist <- list()
  for (i in 1:length(modelselres)) {
    tempdf <- cbind(datadrivenres[[i]],metares[[i]])
    tempdf[1,'SEL'] <- ifelse(modelselres[[i]][1,pcri]>=datadrivenres[[i]][1,pcri],'DDM','GOF')
    tempdf[1,'DIFF'] <- modelselres[[i]][1,pcri] - datadrivenres[[i]][1,pcri]
    histdatalist[[i]] <- tempdf
  }
  # print('end get hist data')
  
  evaluator <- settingDGStraincl(listrbinder(histdatalist),evalfts,'SEL')
  # evaluator <- settingRESEDAtrainRg(listrbinder(histdatalist),evalfts,'DIFF')
  
  return(evaluator)
}

predictionRESEDA <- function(pevaluator,ptargetdf,ppert,
                             pdflist,pmodelvec,pcri,pgoffts,pmetafts,pddmmodel,
                             pfith,ppredh,
                             psrgmestlist=NULL,
                             ptargetddmest=NULL,
                             ptargetsrgmest=NULL){
  
  datadrivenest <- ptargetddmest
  if (is.null(datadrivenest)){
    datadrivenest <- getDDMEstOrig(ptargetdf,pddmmodel,pfith,ppredh)
  }
  # print('pred: end ddm est')
  
  datadrivenres <- getDDMGOFforaData(datadrivenest,ptargetdf,pfith,pfith,pddmmodel,pgoffts)
  # print('pred: end ddm gof')
  
  metares <- getMETAInfoOrig(ptargetdf,pmetafts,pfith)
  # print('pred: end meta')
  
  selorddm <- setting1predictcl(pevaluator,cbind(datadrivenres,metares))
  selectionres <- ifelse(selorddm[1,'Pred']=='GOF','GOF','DDM')
  # selorddm <- settingRESEDApredictRg(pevaluator,cbind(datadrivenres,metares))
  # selectionres <- ifelse(selorddm[1,'Pred']<0,'GOF','DDM')
  # print('pred: end selection')
  
  if (selectionres=='GOF'){
    
    modelest <- ptargetsrgmest
    if (is.null(modelest)){
      modelest <- getSRGMEstOrig(ptargetdf,pmodelvec,pfith,ppredh)
    }
    # print('pred gof: end srgm est')
    
    modelgof <- getSRGMGOFOrig(ptargetdf,modelest,pmodelvec,pgoffts,pfith,pfith)
    # print('pred gof: end srgm gof')
    
    cmodelsest <- psrgmestlist
    if (is.null(cmodelsest)){
      cmodelsest <- modelFitting(pdflist,pmodelvec,ppert)
    }
    # print('pred gof: end modelfitting')
    
    cgofnlist <- calcGOF(pdflist,cmodelsest,pmodelvec,pgoffts,ppert)
    # print('pred gof: end calcGOF')
    
    mselected <- setting111LOOPert(cgofnlist,modelgof,pgoffts,pcri)
    # print('pred gof: end model selection')
    
    return(list(estimated=modelest[[mselected]],sel='GOF'))
    
  }else if (selectionres=='DDM'){
    augtargetdf <- data.frame(n=c(ptargetdf$n[1:pfith],datadrivenest$SVR$EstElap[(pfith+1):ppredh]))
    modelest <- getSRGMEstOrig(augtargetdf,pmodelvec,ppredh,ppredh)
    # print('pred ddm: end srgmfromddmest')
    
    srgmgoffromddm <- getSRGMGOFOrig(augtargetdf,modelest,pmodelvec,pgoffts,ppredh,ppredh)
    # print('pred ddm: end srgmfromddmgof')
    
    srgmgoffromddm <- srgmgoffromddm[order(srgmgoffromddm$MSE),]
    
    return(list(estimated=modelest[[srgmgoffromddm[1,'Model']]],sel='DDM'))
  }
}

predictionRESEDAv2 <- function(ptargetdf,pddmmodel,pfith,ppredh,pmodelvec,pgoffts,ptargetddmest=NULL){
  
  datadrivenest <- ptargetddmest
  if (is.null(datadrivenest)){
    datadrivenest <- getDDMEstOrig(ptargetdf,pddmmodel,pfith,ppredh)
  }
  print('pred: end ddm est')
  
  augtargetdf <- data.frame(n=c(ptargetdf$n[1:pfith],datadrivenest[[pddmmodel]]$EstElap[(pfith+1):ppredh]))
  modelest <- getSRGMEstOrig(augtargetdf,pmodelvec,ppredh,ppredh)
  print('pred ddm: end srgmfromddmest')
  
  srgmgoffromddm <- getSRGMGOFOrig(augtargetdf,modelest,pmodelvec,pgoffts,ppredh,ppredh)
  print('pred ddm: end srgmfromddmgof')
  
  srgmgoffromddm <- srgmgoffromddm[order(srgmgoffromddm$MSE),]
  
  selectedmodel <- gofModelSelection(pgofdf =  srgmgoffromddm,
                                    pgoffts = pgoffts)
  
  return(modelest[[selectedmodel]])
}

modelFitting <- function(datadflist,modelvec,ppert){
  
  paramlist <- list()
  for (i in 1:length(datadflist)) {
    paramlist[[i]] <- list()
    paramlist[[i]]$df <- datadflist[[i]]
    paramlist[[i]]$modelvec <- modelvec
    paramlist[[i]]$hori <- horiCalculator(ppert,nrow(datadflist[[i]]))
  }
  if (Sys.info()['sysname']=='Windows'){
    return(lapply(paramlist,modelFittingApply))
  }else{
    return(mclapply(paramlist,modelFittingApply,mc.cores = detectCores()))
  }
}

modelFittingApply <- function(x){
  return(getSRGMEstOrig(df=x$df,
                        modelvec=x$modelvec,
                        trainh=x$hori$fit,predh=x$hori$pred))
}

calcGOF <- function(datadflist,pestlist,pmodelvec,pgofvec,ppert){
  stopifnot(length(datadflist)==length(pestlist))
  paramlist <- list()
  for (i in 1:length(datadflist)) {
    paramlist[[i]] <- list()
    paramlist[[i]]$datadf <- datadflist[[i]]
    paramlist[[i]]$pest <- pestlist[[i]]
    paramlist[[i]]$pmodelvec <- pmodelvec
    paramlist[[i]]$pgofvec <- pgofvec
    paramlist[[i]]$hori <- horiCalculator(ppert,nrow(datadflist[[i]]))
  }
  if (Sys.info()['sysname']=='Windows'){
    return(lapply(paramlist,calcGOFApply))
  }else{
    return(mclapply(paramlist,calcGOFApply,mc.cores = detectCores()))
  }
}

calcGOFApply <- function(x){
  return(getSRGMGOFOrig(df=x$datadf,pEst=x$pest,
                        modelvec = x$pmodelvec,
                        gofvec = x$pgofvec,
                        trainh = x$hori$fit,predh = x$hori$pred))
}

calcGOFExt <- function(datadflist,pestlist,pmodelvec,pgofvec,ppert){
  stopifnot(length(datadflist)==length(pestlist))
  paramlist <- list()
  for (i in 1:length(datadflist)) {
    paramlist[[i]] <- list()
    paramlist[[i]]$datadf <- datadflist[[i]]
    paramlist[[i]]$pest <- pestlist[[i]]
    paramlist[[i]]$pmodelvec <- pmodelvec
    paramlist[[i]]$pgofvec <- pgofvec
    paramlist[[i]]$hori <- horiCalculator(ppert,nrow(datadflist[[i]]))
  }
  if (Sys.info()['sysname']=='Windows'){
    return(lapply(paramlist,calcGOFExtApply))
  }else{
    return(mclapply(paramlist,calcGOFExtApply,mc.cores = detectCores()))
  }
}

calcGOFExtApply <- function(x){
  return(getSRGMGOFOrigExt(df=x$datadf,pEst=x$pest,
                        modelvec = x$pmodelvec,
                        gofvec = x$pgofvec,
                        trainh = x$hori$fit,predh = x$hori$pred))
}

calcGOFDDM <- function(pestlist,pdflist,ppert,pgofvec){
  stopifnot(length(pestlist)==length(pdflist))
  paramlist <- list()
  for (i in 1:length(pdflist)) {
    paramlist[[i]] <- list()
    paramlist[[i]]$pest <- pestlist[[i]]
    paramlist[[i]]$pdf <- pdflist[[i]]
    paramlist[[i]]$hori <- horiCalculator(ppert,nrow(pdflist[[i]]))
    paramlist[[i]]$model <- 'SVR'
    paramlist[[i]]$pgofvec <- pgofvec
  }
  if (Sys.info()['sysname']=='Windows'){
    return(lapply(paramlist,calcGOFDDMApply))
  }else{
    return(mclapply(paramlist,calcGOFDDMApply,mc.cores = detectCores()))
  }
}

calcGOFDDMApply <- function(x){
  return(getDDMGOFforaData(pEst = x$pest,df=x$pdf,
                           ptrainh = x$hori$fit,
                           ppredh = x$hori$pred,
                           pmodel = x$model,pgoffts = x$pgofvec))
}

calcGOFDDMExt <- function(pestlist,pdflist,ppert,pgofvec){
  stopifnot(length(pestlist)==length(pdflist))
  paramlist <- list()
  for (i in 1:length(pdflist)) {
    paramlist[[i]] <- list()
    paramlist[[i]]$pest <- pestlist[[i]]
    paramlist[[i]]$pdf <- pdflist[[i]]
    paramlist[[i]]$hori <- horiCalculator(ppert,nrow(pdflist[[i]]))
    paramlist[[i]]$model <- 'SVR'
    paramlist[[i]]$pgofvec <- pgofvec
  }
  if (Sys.info()['sysname']=='Windows'){
    return(lapply(paramlist,calcGOFDDMExtApply))
  }else{
    return(mclapply(paramlist,calcGOFDDMExtApply,mc.cores = detectCores()))
  }
}

calcGOFDDMExtApply <- function(x){
  return(getDDMGOFExtforaData(pEst = x$pest,df=x$pdf,
                           ptrainh = x$hori$fit,
                           ppredh = x$hori$pred,
                           pmodel = x$model,pgoffts = x$pgofvec))
}

calcMeta <- function(pdflist,pmetafts,ppert){
  paramlist <- list()
  for (i in 1:length(pdflist)) {
    paramlist[[i]] <- list()
    paramlist[[i]]$pdf <- pdflist[[i]]
    paramlist[[i]]$pmetafts <- pmetafts
    paramlist[[i]]$hori <- horiCalculator(ppert,nrow(pdflist[[i]]))
  }
  if (Sys.info()['sysname']=='Windows'){
    return(lapply(paramlist,calcMetaApply))
  }else{
    return(mclapply(paramlist,calcMetaApply,mc.cores = detectCores()))
  }
}

calcMetaApply <- function(x){
  return(getMETAInfoOrig(df=x$pdf,metafts = x$pmetafts,
                         ptrainh = x$hori$fit))
}

modelSelection <- function(pgoflist,pgoffts,pcri){
  paramlist <- list()
  for (i in 1:length(pgoflist)) {
    paramlist[[i]] <- list()
    paramlist[[i]]$pgoflist <- pgoflist[-i]
    paramlist[[i]]$ptdf <- pgoflist[[i]]
    paramlist[[i]]$pgoffts <- pgoffts
    paramlist[[i]]$pcri <- pcri
  }
  if (Sys.info()['sysname']=='Windows'){
    return(lapply(paramlist,modelSelectionApply))
  }else{
    return(mclapply(paramlist,modelSelectionApply,mc.cores = detectCores()))
  }
}

modelSelectionApply <- function(x){
  return(setting111LOOPert(goflist = x$pgoflist,
                           targetdf = x$ptdf,
                           gofregvec = x$pgoffts,
                           pcri = x$pcri))
}

gofModelSelection <- function(pselector=NULL,pgofdf,pgoffts){
  gofdf <- pgofdf
  if (is.null(pselector)){
    gofdf <- gofdf[order(gofdf$MSE),]
    return(gofdf[1,'Model'])
  }
}

dataDriven <- function(pdflist,pmodel,ppert){
  paramlist <- list()
  for (i in 1:length(pdflist)) {
    paramlist[[i]] <- list()
    paramlist[[i]]$pdf <- pdflist[[i]]
    paramlist[[i]]$pmodel <- pmodel
    paramlist[[i]]$hori <- horiCalculator(ppert,nrow(pdflist[[i]]))
  }
  if (Sys.info()['sysname']=='Windows'){
    return(lapply(paramlist,dataDrivenApply))
  }else{
    return(mclapply(paramlist,dataDrivenApply,mc.cores = detectCores()))
  }
}

dataDrivenApply <- function(x){
  return(getDDMEstOrig(df = x$pdf,model = x$pmodel,
                       trainh = x$hori$fit,
                       predh = x$hori$pred))
}

searchForFeatures <- function(phistlist,fullfeat,gaparam){
  
  ftslenbits <- length(fullfeat)
  
  S.Decoder <- function(x){
    
    gx <- x
    
    if(binary2decimal(gx)==0){
      gx <- decimal2binary(1,length = ftslenbits)
    }
    
    retfts <- fullfeat[c(which(gx==1))]
    
    return(retfts)
  }
  
  S.Evaluation <- function(phistlist,feats){
    
    resvec <- vector()
    
    for (i in 1:length(phistlist)) {
      traingdf <- listrbinder(phistlist[-i])
      testdf <- phistlist[[i]]
      
      classifier <- settingDGStraincl(traingdf,feats,'SEL')
      predretdf <- setting1predictcl(classifier,testdf)
      
      resvec[(length(resvec)+1)] <- ifelse(predretdf[1,'Pred']=='GOF','GOF','DDM') 
    }

    fulldf <- listrbinder(phistlist)
    fulldf['DGS.SEL'] <- factor(resvec,c('GOF','DDM'))
    
    totcount <- nrow(fulldf)
    gofcount <- length(which(fulldf[,'SEL']=='GOF'))
    ddmcount <- length(which(fulldf[,'SEL']=='DDM'))
    
    crtcount <- length(which(fulldf[,'SEL']==fulldf[,'DGS.SEL']))
    crcntgof <- length(which(fulldf[,'SEL']==fulldf[,'DGS.SEL'] & fulldf[,'SEL']=='GOF'))
    crcntddm <- length(which(fulldf[,'SEL']==fulldf[,'DGS.SEL'] & fulldf[,'SEL']=='DDM'))
    
    retlist <- list()
    retlist$acc <- crtcount/totcount
    retlist$acc.gof <- crcntgof/gofcount
    retlist$acc.ddm <- crcntddm/ddmcount
    
    return(retlist)
  }
  
  S.Fitness <- function(x){
    
    pfts <- S.Decoder(x)
    
    currcri <- S.Evaluation(phistlist = phistlist,feats = pfts)
    
    currfitness <- (currcri$acc + currcri$acc.gof + currcri$acc.ddm)/3
    
    return(currfitness)
  }
  
  print("GA for feats start!")
  
  lenbits <- ftslenbits
  
  GAret <- ga(type = 'binary',fitness = S.Fitness,nBits = lenbits,
              popSize = gaparam$popsize,maxiter = gaparam$miter,
              keepBest = TRUE,parallel = FALSE,monitor = gaMonitor)
  
  # print(GAret@bestSol[[length(GAret@bestSol)]])
  
  finalsol <- GAret@solution[1,]
  
  print(finalsol)
  
  retfts <- S.Decoder(finalsol)
  
  print(retfts)
  
  return(retfts)
}

horiCalculator <- function(ppert,totlen){
  
  ratio <- NULL
  if (ppert=='P4'){
    ratio <- 0.4
  }else if (ppert=='P5'){
    ratio <- 0.5
  }else if (ppert=='P6'){
    ratio <- 0.6
  }else if (ppert=='P7'){
    ratio <- 0.7
  }else if (ppert=='P8'){
    ratio <- 0.8
  }else if (ppert=='P9'){
    ratio <- 0.9
  }
  
  fithori <- floor(totlen * ratio)
  
  rethori <- list()
  rethori$fit <- fithori
  rethori$pred <- totlen
  
  return(rethori)
}