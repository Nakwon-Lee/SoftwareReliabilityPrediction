searchforEvaluator <- function(plist,pnlist,pmetalist,
                               regvars,ppert,pcri,pddmres=NULL,pgofres=NULL,
                               pfulldf,pfullfeats,pgofddmlist,
                               goffeats,metafeats){
  
  kmax <- 100
  
  goflret <- pgofres
  if(is.null(goflret)){
    goflret <- setting111LOOnoout(goflist = plist, gofregvec = regvars,
                                  totpert = ppert,targets = paste0("N",pcri),
                                  measures = pcri,isFCTBF = "FC",
                                  fntrain = settingGLMtrain, fnpredict = settingBasicpredict)
  }
  
  ddmret <- pddmres
  if(is.null(ddmret)){
    ddmret <- settingDDMnoout(pgofddmlist,ppert,pcri)
  }
  
  xddm <- log2(ddmret[,pcri])
  ygof <- log2(goflret[,pcri])
  origvec <- sort(c(xddm,ygof))
  margvec <- sort(abs(ygof-xddm))
  
  param <- list()
  param$oidx <- ceiling(length(origvec)/4)
  param$orig <- origvec[param$oidx]
  param$midx <- ceiling(length(margvec)/4)
  param$margin <- margvec[param$midx]
  
  gofsubretlist <- list()
  for (i in 1:length(plist)) {
    gofsubretlist[[i]] <- setting111LOOnoout(goflist = plist[-i], gofregvec = regvars,
                                          totpert = ppert,targets = paste0("N",pcri),
                                          measures = pcri,isFCTBF = "FC",
                                          fntrain = settingGLMtrain, fnpredict = settingBasicpredict)
  }
  
  ddmsubretlist <- list()
  for (i in 1:length(pgofddmlist)) {
    ddmsubretlist[[i]] <- settingDDMnoout(pgofddmlist[-i],ppert,pcri)
  }
  
  tgofsublist <- list()
  for (i in 1:length(plist)) {
    tgofsublist[[i]] <- execGOF111(goflist = plist[-i],gofregvec = regvars,totpert = ppert,
               ptarget = pcri,fntrain = settingGLMtrain,fnpredict = settingBasicpredict,
               tgoflist = plist[c(i)])
  }
  
  cricache <- list()
  currpar <- param
  currcri <- Evaluation(plist = plist,ppert = ppert,
                        pgofddmlist = pgofddmlist,pnlist = pnlist,pmetalist = pmetalist,
                        pcri = pcri,
                        gofretlist = gofsubretlist,ddmretlist = ddmsubretlist,tgofretlist = tgofsublist,
                        param = currpar,goffeats = goffeats,metafeats = metafeats,
                        pcache = cricache)
  cricache <- currcri$cache
  bestpar <- currpar
  bestcri <- currcri$val

  for (k in 1:kmax) {
    
    cat(".")
    
    preddf <- NULL
    temP <- temperature(pk = k,pkmax = kmax)
    newpar <- neighbour(param = currpar,pOvec = origvec,pMvec = margvec,bestpar)
    
    newcri <- Evaluation(plist = plist,ppert = ppert,
                         pgofddmlist = pgofddmlist,pnlist = pnlist,pmetalist = pmetalist,
                         pcri = pcri,
                         gofretlist = gofsubretlist,ddmretlist = ddmsubretlist,tgofretlist = tgofsublist,
                         param = newpar,goffeats = goffeats,metafeats = metafeats,
                         pcache = cricache)
    cricache <- newcri$cache
    
    # print(paste0(k," param: ",newcri))
    
    if(newcri$val < bestcri){
      bestcri <- newcri$val
      bestpar <- newpar
    }
    
    if (Acceptance(currcri$val,newcri$val,temP) >= (sample(x = c(1:1000),size = 1)/1000)){
      currpar <- newpar
      currcri <- newcri
    }
  }
  
  # end of the search: "bestpar" is the searched parameter for evaluator generation
  # after the search, we generate evaluator using the searched parameter (i.e., origin and margin)
  
  # metafulldf <- NULL
  # for (i in 1:length(pmetalist)) {
  #   if(is.null(metafulldf)){
  #     metafulldf <- pmetalist[[i]]
  #   }else{
  #     metafulldf <- rbind(metafulldf,pmetalist[[i]])
  #   }
  # }
  # metafullndf <- normalization(metafulldf,metafts)$df
  # 
  # goffulldf <- NULL
  # for (i in 1:length(pgofddmlist)) {
  #   if(is.null(goffulldf)){
  #     goffulldf <- pgofddmlist[[i]]
  #   }else{
  #     goffulldf <- rbind(goffulldf,pgofddmlist[[i]])
  #   }
  # }
  # goffulldf <- goffulldf[goffulldf$Model=="SVR",]
  # 
  # fulldf <- cbind(metafullndf,goffulldf)
  # fullfeat <- c(goffts,paste0("N",metafts))
  
  evalter <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                             pfulldf = pfulldf,pfullfeats = pfullfeats,
                             orimar = bestpar,pcri = pcri)
  
  evalter$orimar <- bestpar
  
  return(evalter)
}

# searchforEvaluatorNoCV <- function(plist,pnlist,pmetalist,regvars,ppert,pcri,pddmres,pgofddmlist,goffeats,metafeats){
#   
#   kmax <- 2
#   
#   goflret <- NULL
#   goflret <- setting111LOOnoout(goflist = plist, gofregvec = regvars,
#                                 totpert = ppert,targets = paste0("N",pcri),
#                                 measures = pcri,isFCTBF = "FC",
#                                 fntrain = settingGLMtrain, fnpredict = settingBasicpredict)
#   
#   ddmret <- pddmres
#   if(is.null(ddmret)){
#     ddmret <- settingDDM(pgofddmlist,ppert,pcri)
#   }
#   
#   xddm <- ddmret[,pcri]
#   ygof <- goflret[,pcri]
#   origvec <- sort(c(xddm,ygof))
#   margvec <- sort(abs(ygof-xddm))
#   
#   param <- list()
#   param$oidx <- ceiling(length(origvec)/2)
#   param$orig <- origvec[param$oidx]
#   param$midx <- ceiling(length(margvec)/2)
#   param$margin <- margvec[param$midx]
#   
#   # normalizing meta features for training and test data
#   metafulldf <- NULL
#   for (i in 1:length(pmetalist)) {
#     if(is.null(metafulldf)){
#       metafulldf <- pmetalist[[i]]
#     }else{
#       metafulldf <- rbind(metafulldf,pmetalist[[i]])
#     }
#   }
#   metafullndf <- normalization(metafulldf,metafeats)$df
#   
#   cricache <- list()
#   
#   currpar <- param
#   currcri <- EvaluationNoCV(plist = plist,ppert = ppert,metafullndf = metafullndf,
#                          pgofddmlist = pgofddmlist,pnlist = pnlist,pmetalist = pmetalist,
#                          regvars = regvars,pcri = pcri,ddmret = ddmret,goflret = goflret,
#                          param = currpar,goffeats = goffeats,metafeats = metafeats,
#                          pcache = cricache)
#   cricache <- currcri$cache
#   
#   bestpar <- currpar
#   bestcri <- currcri$val
#   
#   # print(paste0("init param: ",currcri$val))
#   
#   for (k in 1:kmax) {
#     
#     preddf <- NULL
#     temP <- temperature(pk = k,pkmax = kmax)
#     newpar <- neighbour(param = currpar,pOvec = origvec,pMvec = margvec)
#     
#     newcri <- EvaluationNoCV(plist = plist,ppert = ppert,metafullndf = metafullndf,
#                           pgofddmlist = pgofddmlist,pnlist = pnlist,pmetalist = pmetalist,
#                           regvars = regvars,pcri = pcri,ddmret = ddmret,goflret = goflret,
#                           param = newpar,goffeats = goffeats,metafeats = metafeats,
#                           pcache = cricache)
#     cricache <- newcri$cache
#     
#     # print(paste0(k," param: ",newcri$val))
#     
#     if(newcri$val < bestcri){
#       bestcri <- newcri$val
#       bestpar <- newpar
#     }
#     
#     if (Acceptance(currcri$val,newcri$val,temP) >= (sample(x = c(1:1000),size = 1)/1000)){
#       currpar <- newpar
#       currcri <- newcri
#     }
#     
#     # vvec <- which(!is.na(evalfullvec))
#     # evalfiltered <- evalfullvec[vvec]
#     # predfiltered <- preddf$Pred[vvec]
#     # correct.all <- length(which(evalfiltered==predfiltered))
#     # count.all <- length(vvec)
#     # correct.ddm <- length(which(evalfiltered==predfiltered & evalfiltered=="DDM"))
#     # count.ddm <- length(evalfiltered=="DDM")
#     # correct.gof <- length(which(evalfiltered==predfiltered & evalfiltered=="GOF"))
#     # count.gof <- length(evalfiltered=="GOF")
#     # 
#     # crit.acc <- (correct.all/count.all)*100
#     # crit.rat <- ((correct.ddm/count.ddm)*100-(correct.gof/count.gof)*100)^2
#     
#   }
#   
#   # end of the search: "bestpar" is the searched parameter for evaluator generation
#   # after the search, we generate evaluator using the searched parameter (i.e., origin and margin)
#   
#   evalter <- DDMevaluatorGen(plist = plist,pnlist = pnlist,pmetalist = pmetalist,
#                               regvars = regvars,ppert = ppert,pcri = pcri,
#                               pddmres = pddmres,pgofres = goflret,
#                               orimar = bestpar,pgofddmlist = pgofddmlist,
#                               goffeats = goffeats,metafeats = metafeats)
#   
#   evalter$orimar <- bestpar
#   
#   return(evalter)
# }


DDMevaluatorGen <- function(pddmres,pgofres,
                            pfulldf,pfullfeats,
                            orimar=NULL,pcri){
  
  fulldf <- pfulldf
  fullfeat <- pfullfeats
  
  param <- orimar
  if(is.null(param)){
    param <- list()
    param$orig <- 1
    param$margin <- 1
    fulldf[paste0("DGS.SEL.",pcri)] <- as.factor(filteringLabeler(pddmres[,pcri],pgofres[,pcri],param))
  }else{
    fulldf[paste0("DGS.SEL.",pcri)] <- as.factor(filteringLabeler(pddmres[,pcri],pgofres[,pcri],param))
  }
  
  traindf <- fulldf[,c(pfullfeats,paste0("DGS.SEL.",pcri))]
  traindfsub <- na.omit(traindf)
  
  traindfsub[paste0("DGS.SEL.BIN.",pcri)] <- ifelse(traindfsub[paste0("DGS.SEL.",pcri)]=="DDM",1,0)
  
  fullfeat <- featureSelection(ptrain = traindfsub,pfts = fullfeat,pcri = pcri)
  
  # SMOTE
  # traindfsub <- overSamplingBal(origdf = traindfsub,vars = fullfeat,target = paste0("DGS.SEL.",pcri))
  
  classifier <- setting1traincl(traindfsub,fullfeat,paste0("DGS.SEL.",pcri))
  
  testidf <- fulldf[,c(pfullfeats,paste0("DGS.SEL.",pcri))]
  predretdf <- setting1predictcl(classifier,fullfeat,testidf)
  fulldf[paste0("DGS.PRD.",pcri)] <- predretdf$Pred
  
  ccc <- length(which(fulldf[paste0("DGS.PRD.",pcri)]==fulldf[paste0("DGS.SEL.",pcri)]))
  dddvec <- fulldf[paste0("DGS.SEL.",pcri)][!is.na(fulldf[paste0("DGS.SEL.",pcri)])]
  ddd <- length(dddvec)
  mmm <- length(which(fulldf[paste0("DGS.PRD.",pcri)]==fulldf[paste0("DGS.SEL.",pcri)] & fulldf[paste0("DGS.PRD.",pcri)]=="DDM"))
  mmmall <- length(which(dddvec=="DDM"))
  ggg <- length(which(fulldf[paste0("DGS.PRD.",pcri)]==fulldf[paste0("DGS.SEL.",pcri)] & fulldf[paste0("DGS.PRD.",pcri)]=="GOF"))
  gggall <- length(which(dddvec=="GOF"))
  
  ret <- list()
  ret$evaluator <- classifier
  ret$performance <- ccc / ddd
  ret$string <- paste0(ccc," / ",ddd," DDM ",mmm," / ",mmmall,"(",mmm/mmmall,")"," GOF ",ggg," / ",gggall,"(",ggg/gggall,")")
  ret$fts <- fullfeat
  ret$param <- param

  return(ret)
}

filteringLabeler <- function(xax,yax,porimar){
  
  stopifnot(length(xax)==length(yax))
  
  selvec <- vector()
  
  xddm <- log2(xax)
  ysel <- log2(yax)
  
  for (j in 1:length(xax)){
    
    selvec[j] <- NA
    
    if(xddm[j] > porimar$orig){
      if(ysel[j] < (xddm[j] - porimar$margin)){
        selvec[j] <- "GOF"
      }
    }
    
    if(ysel[j] > porimar$orig){
      if(ysel[j] > (xddm[j] + porimar$margin)){
        selvec[j] <- "DDM"
      }
    }
  }
  
  return(selvec)
}

featureSelection <- function(ptrain,pfts,pcri){
  
  tform <- as.simple.formula(pfts,pcri)
  
  retfts <- FSelector::cfs(tform,ptrain)
  
  # retfts <- pfts
  # key <- TRUE
  # while(key){
  #   tmodel <- temptrain(ptrain,retfts,paste0("DGS.SEL.BIN.",pcri))
  #   retvif <- car::vif(tmodel)
  #   if (length(which(retvif>10))>0){
  #     vartorm <- names(retvif)[which.max(retvif)]
  #     idxtorm <- which(retfts==vartorm)
  #     retfts <- retfts[-idxtorm]
  #   }else{
  #     key <- FALSE
  #   }
  # }
  
  return(retfts)
}

overSamplingBal <- function(origdf,vars,target){
  
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  retdf <- DMwR::SMOTE(form = targetformula,data = origdf,perc.over = 5000,perc.under = 100,k = 5)
  
  return(retdf)
}

temperature <- function(pk,pkmax){
  
  tk <- pk/pkmax
  
  return(0.07/(tk+0.07))
}

neighbour <- function(param,pOvec,pMvec,bparam){
  
  neipar <- param
  
  keyvec <- c(1:10)
  
  rmidx <- vector()
  
  if(neipar$oidx==length(pOvec)){
    rmidx <- c(rmidx,c(1,5,6))
  }
  if(neipar$oidx==1){
    rmidx <- c(rmidx,c(2,7,8))
  }
  if(neipar$midx==length(pMvec)){
    rmidx <- c(rmidx,c(3,5,7))
  }
  if(neipar$midx==1){
    rmidx <- c(rmidx,c(4,6,8))
  }
  
  rmidx <- unique(rmidx)
  
  if(length(rmidx)>0){
    keyvec <- keyvec[-rmidx]
  }
  
  key <- sample(keyvec,size = 1)
  
  if(key==1){
    neipar$oidx <- neipar$oidx + 1 
  }else if(key==2){
    neipar$oidx <- neipar$oidx - 1
  }else if(key==3){
    neipar$midx <- neipar$midx + 1
  }else if(key==4){
    neipar$midx <- neipar$midx - 1
  }else if(key==5){
    neipar$oidx <- neipar$oidx + 1
    neipar$midx <- neipar$midx + 1
  }else if(key==6){
    neipar$oidx <- neipar$oidx + 1
    neipar$midx <- neipar$midx - 1
  }else if(key==7){
    neipar$oidx <- neipar$oidx - 1
    neipar$midx <- neipar$midx + 1
  }else if(key==8){
    neipar$oidx <- neipar$oidx - 1
    neipar$midx <- neipar$midx - 1
  }else if(key==9){
    neipar$oidx <- sample(c(1:pOvec))
    neipar$midx <- sample(c(1:pMvec))
  }else{#key==10
    neipar$oidx <- bparam$oidx
    neipar$midx <- bparam$midx
  }
  
  neipar$orig <- pOvec[neipar$oidx]
  neipar$margin <- pMvec[neipar$midx]
  
  return(neipar)
}

Acceptance <- function(pcfit,pnfit,ptmp){
  ret <- 0
  if (pnfit < pcfit){
    ret <- 1
  }else{
    ret <- exp(-(pnfit-pcfit)/ptmp)
  }
  return(ret)
}

Evaluation <- function(plist,ppert,pgofddmlist,
                       pnlist,pmetalist,pcri,
                       gofretlist,ddmretlist,tgofretlist,
                       param,goffeats,metafeats,
                       pcache){
  
  if (is.null(pcache[[paste0(param$oidx,".",param$midx)]])){
    resultvec <- vector()
    
    for (i in 1:length(plist)) {
      picked <- c(i)
      
      DGSretdf <- execDGS(pgofldf = gofretlist[[i]],pddmodf = ddmretlist[[i]],ptgofdf = tgofretlist[[i]],
                          goflist = plist[-picked],gofnlist = pnlist[-picked],
                          gofddmlist = pgofddmlist[-picked],metadf = pmetalist[-picked],
                          pcri = pcri,totpert = ppert,
                          goffts = goffeats,metafts = metafeats,
                          psearch = "GIVEN",givenparam = param,
                          tmetadf = pmetalist[picked],tgofddmlist = pgofddmlist[picked])
      
      for (j in 1:length(picked)) {
        testirows <- c(1:ppert)
        testirows <- testirows + (ppert * (j-1))
        
        temdgsret <- DGSretdf[testirows,]
        
        for (k in 1:nrow(temdgsret)) {
          if (temdgsret[k,paste0("DGS.RET.",pcri)]=="DDM"){
            trows <- which(pgofddmlist[[picked[j]]]$Model=="SVR" & pgofddmlist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
            resultvec[(length(resultvec)+1)] <- pgofddmlist[[picked[j]]][trows,][1,pcri]
          }else{
            trows <- which(plist[[picked[j]]]$Pert==temdgsret[k,"Pert"])
            tgofdf <- plist[[picked[j]]][trows,]
            resultvec[(length(resultvec)+1)] <- tgofdf[tgofdf$Model==temdgsret[k,paste0("DGS.RET.",pcri)],][1,pcri]
          }
        }
      }
    }
    
    retdf <- data.frame(resultvec)
    colnames(retdf) <- c(pcri)
    
    ret <- list()
    ret$val <- median(retdf[,pcri])
    pcache[[paste0(param$oidx,".",param$midx)]] <- median(retdf[,pcri])
    ret$cache <- pcache
    
    return(ret)
  }else{
    ret <- list()
    ret$val <- pcache[[paste0(param$oidx,".",param$midx)]]
    ret$cache <- pcache
    
    return(ret)
  }
}

EvaluationNoCV <- function(plist,ppert,metafullndf,pgofddmlist,pnlist,
                           pmetalist,regvars,pcri,ddmret,goflret,gofretlist = NULL,
                           param,goffeats,metafeats,
                           pcache){
  
  if (is.null(pcache[[paste0(param$oidx,".",param$midx)]])){
    teval <- DDMevaluatorGen(plist = plist,pnlist = pnlist,pmetalist = pmetalist,
                              regvars = regvars,ppert = ppert,pcri = pcri,
                              pddmres = ddmret,pgofres = goflret,orimar = param,pgofddmlist = pgofddmlist,
                              goffeats = goffeats,metafeats = metafeats)
    
    goffulldf <- NULL
    for (i in 1:length(pnlist)) {
      if(is.null(goffulldf)){
        goffulldf <- pnlist[[i]]
      }else{
        goffulldf <- rbind(goffulldf,pnlist[[i]])
      }
    }
    goffulldf <- goffulldf[goffulldf$Model=="SVR",]
    
    ftsdf <- cbind(goffulldf,metafullndf)
    
    #classification
    preddf <- setting1predictcl(teval$evaluator,teval$fts,ftsdf)
    DGSret <- ifelse(preddf$Pred=="DDM",ddmret[,pcri],goflret[,pcri])
    retdf <- NULL
    retdf <- data.frame(DGSret)
    colnames(retdf) <- pcri
    
    ret <- list()
    ret$val <- median(retdf[,pcri])
    pcache[[paste0(param$oidx,".",param$midx)]] <- median(retdf[,pcri])
    ret$cache <- pcache
    
    return(ret)
  }else{
    ret <- list()
    ret$val <- pcache[[paste0(param$oidx,".",param$midx)]]
    ret$cache <- pcache
    
    return(ret)
  }
}