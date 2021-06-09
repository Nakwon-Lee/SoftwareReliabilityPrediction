searchforEvaluator <- function(plist,pnlist,pmetalist,
                               regvars,ppert,pcri,pddmres=NULL,pgofres=NULL,
                               pfulldf=NULL,pfullfeats=NULL,pgofddmlist,
                               goffeats,metafeats,
                               pkmax = NULL){
  
  kmax <- pkmax
  if(is.null(kmax)){
    kmax <- 1000
  }
  
  print(kmax)
  
  goflret <- pgofres
  if(is.null(goflret)){
    goflret <- setting111LOOnoout(goflist = plist, gofregvec = regvars,
                                  totpert = ppert,target = paste0("N",pcri),
                                  measure = pcri,
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
  param$oidx <- 1
  param$orig <- origvec[param$oidx]
  param$midx <- 1
  param$margin <- margvec[param$midx]
  
  gofsubretlist <- list()
  for (i in 1:length(plist)) {
    gofsubretlist[[i]] <- setting111LOOnoout(goflist = plist[-i], gofregvec = regvars,
                                          totpert = ppert,target = paste0("N",pcri),
                                          measure = pcri,
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
    
    temP <- temperature(pk = k,pkmax = kmax)
    newpar <- neighbour(param = currpar,pOvec = origvec,pMvec = margvec,bparam = bestpar)
    
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
  
  fulldf <- pfulldf
  fullfeat <- pfullfeats
  
  if(is.null(fulldf)){
    metafulldf <- NULL
    for (i in 1:length(pmetalist)) {
      if(is.null(metafulldf)){
        metafulldf <- pmetalist[[i]]
      }else{
        metafulldf <- rbind(metafulldf,pmetalist[[i]])
      }
    }
    metafullndf <- normalization(metafulldf,metafeats)$df
    
    goffulldf <- NULL
    for (i in 1:length(pgofddmlist)) {
      if(is.null(goffulldf)){
        goffulldf <- pgofddmlist[[i]]
      }else{
        goffulldf <- rbind(goffulldf,pgofddmlist[[i]])
      }
    }
    goffulldf <- goffulldf[goffulldf$Model=="SVR",]
    
    fulldf <- cbind(metafullndf,goffulldf)
    fullfeat <- c(goffeats,paste0("N",metafeats))
  }
  
  evalter <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                             pfulldf = fulldf,pfullfeats = fullfeat,
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

searchforOriMarSA <- function(plist,pnlist,pmetalist,
                              regvars,ppert,pcri,pddmres=NULL,pgofres=NULL,
                              pfulldf=NULL,pfullfeats=NULL,pgofddmlist,
                              goffeats,metafeats,
                              pkmax = NULL,
                              pparam = NULL,gaparam = NULL,
                              pparallel='seq'){
  
  kmax <- pkmax
  if(is.null(kmax)){
    kmax <- 1000
  }
  
  goflret <- pgofres
  if(is.null(goflret)){
    goflret <- setting111LOOnoout(goflist = plist, gofregvec = regvars,
                                  totpert = ppert,target = paste0("N",pcri),
                                  measure = pcri,
                                  fntrain = settingGLMtrain, fnpredict = settingBasicpredict)
  }
  
  ddmret <- pddmres
  if(is.null(ddmret)){
    ddmret <- settingDDMnoout(pgofddmlist,ppert,pcri)
  }
  
  param <- pparam
  if(is.null(param)){
    xddm <- log2(ddmret[,pcri])
    ygof <- log2(goflret[,pcri])
    origvec <- sort(c(xddm,ygof))
    margvec <- sort(abs(ygof-xddm))
    
    param <- list()
    param$oidx <- 1
    param$orig <- origvec[param$oidx]
    param$midx <- 1
    param$margin <- margvec[param$midx]
  }
  
  gofsubretlist <- list()
  for (i in 1:length(plist)) {
    gofsubretlist[[i]] <- setting111LOOnoout(goflist = plist[-i], gofregvec = regvars,
                                             totpert = ppert,target = paste0("N",pcri),
                                             measure = pcri,
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
  
  fullfeat <- c(goffeats,paste0("N",metafeats))
  
  valrange <- 1:((2^length(fullfeat))-1)
  bits <- decimal2binary(max(valrange))
  lenbits <- length(bits)
  
  cricache <- list()
  
  Fitness <- function(x){
    
    if(binary2decimal(x)==0){
      x <- decimal2binary(1,length = lenbits)
    }
    
    currfts <- fullfeat[c(which(x==1))]
    
    print(currfts)
    print(cricache)
    
    currcri <- Evaluation(plist = plist,ppert = ppert,
                          pgofddmlist = pgofddmlist,pnlist = pnlist,pmetalist = pmetalist,
                          pcri = pcri,
                          gofretlist = gofsubretlist,ddmretlist = ddmsubretlist,tgofretlist = tgofsublist,
                          param = param,goffeats = goffeats,metafeats = metafeats,
                          pcache = cricache,pfts = currfts,pftnum = binary2decimal(x))
    
    print(currcri$cache)
    print(paste0('this fitness: ',-currcri$val))
    
    return(-currcri$val)
  }
  
  print('GA start!')
  
  if(pparallel=='seq'){
    if(is.null(gaparam)){
      GAret <- ga(type = 'binary',fitness = Fitness,nBits = lenbits,
                  popSize = 20,maxiter = 100,run = 100,keepBest = TRUE)
    }else{
      GAret <- ga(type = 'binary',fitness = Fitness,nBits = lenbits,
                  popSize = gaparam$popsize,maxiter = gaparam$miter,keepBest = TRUE)
    }
  }else if(pparallel=='par'){
    if(is.null(gaparam)){
      GAret <- ga(type = 'binary',fitness = Fitness,nBits = lenbits,
                  popSize = 20,maxiter = 100,run = 100,keepBest = TRUE,parallel = TRUE)
    }else{
      GAret <- ga(type = 'binary',fitness = Fitness,nBits = lenbits,
                  popSize = gaparam$popsize,maxiter = gaparam$miter,keepBest = TRUE,parallel = TRUE)
    }
  }
  
  print(GAret@bestSol[[length(GAret@bestSol)]])
  
  finalsol <- GAret@solution
  
  print(finalsol)
  
  selfts <- fullfeat[c(which(finalsol[1,]==1))]
  
  print(selfts)
  
  fulldf <- pfulldf
  fullfeat <- pfullfeats
  
  if(is.null(fulldf)){
    metafulldf <- NULL
    for (i in 1:length(pmetalist)) {
      if(is.null(metafulldf)){
        metafulldf <- pmetalist[[i]]
      }else{
        metafulldf <- rbind(metafulldf,pmetalist[[i]])
      }
    }
    metafullndf <- normalization(metafulldf,metafeats)$df
    
    goffulldf <- NULL
    for (i in 1:length(pgofddmlist)) {
      if(is.null(goffulldf)){
        goffulldf <- pgofddmlist[[i]]
      }else{
        goffulldf <- rbind(goffulldf,pgofddmlist[[i]])
      }
    }
    goffulldf <- goffulldf[goffulldf$Model=="SVR",]
    
    fulldf <- cbind(metafullndf,goffulldf)
    fullfeat <- c(goffeats,paste0("N",metafeats))
  }
  
  evalter <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                             pfulldf = fulldf,pfullfeats = fullfeat,
                             pselfts = selfts,
                             orimar = param,pcri = pcri)
  
  evalter$orimar <- param
  
  return(evalter)
}

searchforOriMarGA <- function(plist,pnlist,pmetalist,
                              ppert,pcri,pddmres,pgofres,
                              pgofddmlist,
                              pgofretlist,pddmretlist,tgofretlist,
                              goffeats,metafeats,
                              gaparam = NULL,
                              pparallel='seq'){
  
  xddm <- log2(pddmres[,pcri])
  ygof <- log2(pgofres[,pcri])
  origvec <- sort(c(xddm,ygof))
  margvec <- sort(abs(ygof-xddm))
  
  ovecrange <- c(1:floor(length(origvec)/2))
  mvecrange <- c(1:floor((length(margvec)*3)/4))
  
  print(paste0('Ori max: ',max(ovecrange)))
  print(paste0('Mag max: ',max(mvecrange)))
  
  ovcbits <- decimal2binary((max(ovecrange)-1))
  ovclenbits <- length(ovcbits)
  
  mvcbits <- decimal2binary((max(mvecrange)-1))
  mvclenbits <- length(mvcbits)
  
  Decoder1 <- function(x){
    
    gx <- x
    
    ori <- min((max(ovecrange)-1),binary2decimal(gx[1:ovclenbits]))+1
    
    marlen <- c(1:mvclenbits) + ovclenbits
    mar <- min((max(mvecrange)-1),binary2decimal(gx[c(marlen)]))+1
    
    retp <- list()
    retp$oidx <- ori
    retp$orig <- origvec[retp$oidx]
    retp$midx <- mar
    retp$margin <- margvec[retp$midx]
    
    return(retp)
  }
  
  cricache <- list()
  
  Fitness1 <- function(x){
    
    ppar <- Decoder1(x)
    
    filtvec <- filteringLabeler(pddmres[,pcri],pgofres[,pcri],ppar)
    
    currcri <- Evaluation(plist = plist,ppert = ppert,
                          pgofddmlist = pgofddmlist,pnlist = pnlist,pmetalist = pmetalist,
                          pcri = pcri,pfilvec = filtvec,
                          gofretlist = pgofretlist,ddmretlist = pddmretlist,tgofretlist = tgofretlist,
                          param = ppar,goffeats = goffeats,metafeats = metafeats,
                          pcache = cricache)
    
    currfitness <- currcri$val
    
    print(paste0('this fitness: ',currfitness))
    
    return(-currfitness)
  }
  
  print('GA start!')
  
  lenbits <- ovclenbits + mvclenbits
  
  if(pparallel=='seq'){
    if(is.null(gaparam)){
      GAret <- ga(type = 'binary',fitness = Fitness1,nBits = lenbits,
                  popSize = 20,maxiter = 100,run = 100,
                  keepBest = TRUE,monitor = gaMonitor)
    }else{
      GAret <- ga(type = 'binary',fitness = Fitness1,nBits = lenbits,
                  popSize = gaparam$popsize,maxiter = gaparam$miter,
                  keepBest = TRUE,monitor = gaMonitor)
    }
  }else if(pparallel=='par'){
    if(is.null(gaparam)){
      GAret <- ga(type = 'binary',fitness = Fitness1,nBits = lenbits,
                  popSize = 20,maxiter = 100,run = 100,
                  keepBest = TRUE,parallel = TRUE,monitor = gaMonitor)
    }else{
      GAret <- ga(type = 'binary',fitness = Fitness1,nBits = lenbits,
                  popSize = gaparam$popsize,maxiter = gaparam$miter,
                  keepBest = TRUE,parallel = TRUE,monitor = gaMonitor)
    }
  }
  
  print(GAret@bestSol[[length(GAret@bestSol)]])
  
  orimarsol <- GAret@solution[1,]
  
  print(orimarsol)
  
  finalorimar <- Decoder1(orimarsol)
  
  print(finalorimar)
  
  return(finalorimar)
}

searchforOriMarExhaustive <- function(plist,pnlist,pmetalist,
                                      ppert,pcri,pddmres,pgofres,
                                      pgofddmlist,
                                      pgofretlist,pddmretlist,tgofretlist,
                                      goffeats,metafeats,
                                      gaparam = NULL,
                                      pparallel='seq'){
  
  xddm <- log2(pddmres[,pcri])
  ygof <- log2(pgofres[,pcri])
  origvec <- sort(c(xddm,ygof))
  margvec <- sort(abs(ygof-xddm))
  
  ovecrange <- c(1:floor(length(origvec)/2))
  mvecrange <- c(1:floor((length(margvec)*3)/4))
  
  print(paste0('Ori max: ',max(ovecrange)))
  print(paste0('Mag max: ',max(mvecrange)))
  
  pairsvec <- list()
  for (i in ovecrange) {
    for (j in mvecrange) {
      pairsvec[[(length(pairsvec)+1)]] <- c(i,j)
    }
  }
  
  cricache <- list()
  
  Fitness <- function(x){
    
    ori <- x[1]
    mar <- x[2]
    
    retp <- list()
    retp$oidx <- ori
    retp$orig <- origvec[retp$oidx]
    retp$midx <- mar
    retp$margin <- margvec[retp$midx]
    
    filtvec <- filteringLabeler(pddmres[,pcri],pgofres[,pcri],retp)
    
    currcri <- Evaluation(plist = plist,ppert = ppert,
                          pgofddmlist = pgofddmlist,pnlist = pnlist,pmetalist = pmetalist,
                          pcri = pcri,pfilvec = filtvec,
                          gofretlist = pgofretlist,ddmretlist = pddmretlist,tgofretlist = tgofretlist,
                          param = retp,goffeats = goffeats,metafeats = metafeats,
                          pcache = cricache)
    
    currfitness <- currcri$val
    
    print(paste0('this fitness: ',currfitness))
    
    return(currfitness)
  }
  
  results <- unlist(mclapply(pairsvec,Fitness,mc.cores = detectCores()))
  thepair <- pairsvec[[which.min(results)]]
  
  retpar <- list()
  retpar$oidx <- thepair[1]
  retpar$orig <- origvec[retpar$oidx]
  retpar$midx <- thepair[2]
  retpar$margin <- margvec[retpar$midx]
  
  print(retpar)
  
  return(retpar)
}

searchforFtsGA <- function(plist,pnlist,pmetalist,
                                ppert,pcri,pddmres,pgofres,
                           pgofretlist,pddmretlist,tgofretlist,
                                pgofddmlist,
                                goffeats,metafeats,
                                pparam,gaparam2 = NULL,
                                pparallel='seq'){
  
  param <- pparam
  
  filtvecwpar <- filteringLabeler(pddmres[,pcri],pgofres[,pcri],param)
  
  fullfeat <- c(goffeats,paste0("N",metafeats))
  
  ftslenbits <- length(fullfeat)
  
  Decoder2 <- function(x){
    
    gx <- x
    
    if(binary2decimal(gx)==0){
      gx <- decimal2binary(1,length = ftslenbits)
    }
    
    retfts <- fullfeat[c(which(gx==1))]
    
    return(retfts)
  }
  
  cricache <- list()
  
  Fitness2 <- function(x){
    
    pfts <- Decoder2(x)
    
    currcri <- Evaluation(plist = plist,ppert = ppert,
                          pgofddmlist = pgofddmlist,pnlist = pnlist,pmetalist = pmetalist,
                          pcri = pcri,pfilvec = filtvecwpar,
                          gofretlist = pgofretlist,ddmretlist = pddmretlist,tgofretlist = tgofretlist,
                          param = param,goffeats = goffeats,metafeats = metafeats,
                          pcache = cricache,pfts = pfts,pftnum = binary2decimal(x))
    
    currfitness <- currcri$acc
    
    print(paste0('this fitness: ',currfitness))
    
    return(currfitness)
  }
  
  print("GA for feats start!")
  
  lenbits <- ftslenbits
  
  if(pparallel=='seq'){
    if(is.null(gaparam2)){
      GAret <- ga(type = 'binary',fitness = Fitness2,nBits = lenbits,
                  popSize = 20,maxiter = 100,run = 100,
                  keepBest = TRUE,monitor = gaMonitor)
    }else{
      GAret <- ga(type = 'binary',fitness = Fitness2,nBits = lenbits,
                  popSize = gaparam2$popsize,maxiter = gaparam2$miter,
                  keepBest = TRUE,monitor = gaMonitor)
    }
  }else if(pparallel=='par'){
    if(is.null(gaparam2)){
      GAret <- ga(type = 'binary',fitness = Fitness2,nBits = lenbits,
                  popSize = 20,maxiter = 100,run = 100,
                  keepBest = TRUE,parallel = TRUE,monitor = gaMonitor)
    }else{
      GAret <- ga(type = 'binary',fitness = Fitness2,nBits = lenbits,
                  popSize = gaparam2$popsize,maxiter = gaparam2$miter,
                  keepBest = TRUE,parallel = TRUE,monitor = gaMonitor)
    }
  }
  
  print(GAret@bestSol[[length(GAret@bestSol)]])
  
  finalsol <- GAret@solution[1,]
  
  print(finalsol)
  
  retfts <- Decoder2(finalsol)
  
  print(retfts)
  
  return(retfts)
}

searchforALLGASep <- function(plist,pnlist,pmetalist,
                           regvars,ppert,pcri,pddmres=NULL,pgofres=NULL,
                           pfulldf=NULL,pfullfeats=NULL,pgofddmlist,
                           goffeats,metafeats,
                           porimar = NULL, pfts = NULL,
                           gaparam = NULL,gaparam2 = NULL,
                           pparallel='seq',
                           searchforom,searchforft){
  
  goflret <- pgofres
  if(is.null(goflret)){
    goflret <- setting111LOOnoout(goflist = plist, gofregvec = regvars,
                                  totpert = ppert,target = paste0("N",pcri),
                                  measure = pcri,
                                  fntrain = settingGLMtrain, fnpredict = settingBasicpredict)
  }
  
  ddmret <- pddmres
  if(is.null(ddmret)){
    ddmret <- settingDDMnoout(pgofddmlist,ppert,pcri)
  }
  
  gofsubretlist <- NULL
  load(file = paste0('gofsubretlist.',pcri,'.RData'))
  if(is.null(gofsubretlist)){
    gofsubretlist <- list()
    for (i in 1:length(plist)) {
      gofsubretlist[[i]] <- setting111LOOnoout(goflist = plist[-i], gofregvec = regvars,
                                               totpert = ppert,target = paste0("N",pcri),
                                               measure = pcri,
                                               fntrain = settingGLMtrain, fnpredict = settingBasicpredict)
    }
    save(gofsubretlist,file = paste0('gofsubretlist.',pcri,'.RData'))
  }
  
  ddmsubretlist <- NULL
  load(file = paste0('ddmsubretlist.',pcri,'.RData'))
  if(is.null(ddmsubretlist)){
    ddmsubretlist <- list()
    for (i in 1:length(pgofddmlist)) {
      ddmsubretlist[[i]] <- settingDDMnoout(pgofddmlist[-i],ppert,pcri)
    }
    save(ddmsubretlist,file = paste0('ddmsubretlist.',pcri,'.RData'))
  }
  
  tgofsublist <- NULL
  load(file = paste0('tgofsublist.',pcri,'.RData'))
  if(is.null(tgofsublist)){
    tgofsublist <- list()
    for (i in 1:length(plist)) {
      tgofsublist[[i]] <- execGOF111(goflist = plist[-i],gofregvec = regvars,totpert = ppert,
                                     ptarget = pcri,fntrain = settingGLMtrain,fnpredict = settingBasicpredict,
                                     tgoflist = plist[c(i)])
    }
    save(tgofsublist,file = paste0('tgofsublist.',pcri,'.RData'))
  }
  
  corimar <- porimar
  if(is.null(corimar)){
    corimar <- searchforom(plist = plist,pnlist = pnlist,pmetalist = pmetalist,
                                  ppert = ppert,pcri = pcri,pddmres = ddmret,pgofres = goflret,
                                  pgofddmlist = pgofddmlist,
                                  pgofretlist = gofsubretlist,pddmretlist = ddmsubretlist,tgofretlist = tgofsublist,
                                  goffeats = goffeats,metafeats = metafeats,
                                  gaparam = gaparam,
                                  pparallel = pparallel)
  }
  save(corimar,file = paste0('tmp.orimar.',pcri,'.RData'))
  
  cfts <- pfts
  if(is.null(cfts)){
    cfts <- searchforft(plist = plist,pnlist = pnlist,pmetalist = pmetalist,
                           ppert = ppert,pcri = pcri,pddmres = ddmret,pgofres = goflret,
                           pgofretlist = gofsubretlist,pddmretlist = ddmsubretlist,tgofretlist = tgofsublist,
                           pgofddmlist = pgofddmlist,
                           goffeats = goffeats,metafeats = metafeats,
                           pparam = corimar,gaparam2 = gaparam2,
                           pparallel = pparallel)
  }
  save(cfts,file = paste0('tmp.fts.',pcri,'.RData'))
  
  finalpar <- corimar
  finalpar$fts <- cfts

  selfts <- finalpar$fts
  
  fulldf <- pfulldf
  fullfeat <- pfullfeats
  
  if(is.null(fulldf)){
    metafulldf <- NULL
    for (i in 1:length(pmetalist)) {
      if(is.null(metafulldf)){
        metafulldf <- pmetalist[[i]]
      }else{
        metafulldf <- rbind(metafulldf,pmetalist[[i]])
      }
    }
    metafullndf <- normalization(metafulldf,metafeats)$df
    
    goffulldf <- NULL
    for (i in 1:length(pgofddmlist)) {
      if(is.null(goffulldf)){
        goffulldf <- pgofddmlist[[i]]
      }else{
        goffulldf <- rbind(goffulldf,pgofddmlist[[i]])
      }
    }
    goffulldf <- goffulldf[goffulldf$Model=="SVR",]
    
    fulldf <- cbind(metafullndf,goffulldf)
    fullfeat <- c(goffeats,paste0("N",metafeats))
  }
  
  evalter <- DDMevaluatorGen(pddmres = ddmret,pgofres = goflret,
                             pfulldf = fulldf,pfullfeats = fullfeat,
                             pselfts = selfts,
                             orimar = finalpar,pcri = pcri)
  
  evalter$orimar <- finalpar
  
  return(evalter)
}


DDMevaluatorGen <- function(pddmres,pgofres,
                            pfulldf,pfullfeats,
                            pselfts=NULL,
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
  
  traindf <- fulldf[,c(fullfeat,paste0("DGS.SEL.",pcri))]
  traindfsub <- na.omit(traindf)
  
  traindfsub[paste0("DGS.SEL.BIN.",pcri)] <- ifelse(traindfsub[,paste0("DGS.SEL.",pcri)]=="DDM",1,0)
  
  if(is.null(pselfts)){
    fullfeat <- featureSelection(ptrain = traindfsub,pfts = fullfeat,pcri = pcri)
  }else{
    fullfeat <- pselfts
  }
  
  # SMOTE
  # traindfsub <- overSamplingBal(origdf = traindfsub,vars = fullfeat,target = paste0("DGS.SEL.",pcri))
  
  classifier <- settingDGStraincl(traindfsub,fullfeat,paste0("DGS.SEL.",pcri))
  
  predretdf <- setting1predictcl(classifier,fullfeat,fulldf)
  
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
  
  # tablenames <- names(table(ptrain[,paste0("DGS.SEL.",pcri)]))
  # 
  # retfts <- NULL
  # if(length(tablenames)==2){
  #   tform <- as.simple.formula(pfts,paste0("DGS.SEL.",pcri))
  #   retfts <- FSelector::cfs(tform,ptrain)
  # }else{
  #   retfts <- pfts
  #   key <- TRUE
  #   while(key){
  #     tmodel <- temptrain(ptrain,retfts,paste0("DGS.SEL.BIN.",pcri))
  #     retvif <- regclass::VIF(tmodel)
  #     if (length(which(retvif>10))>0){
  #       vartorm <- names(retvif)[which.max(retvif)]
  #       idxtorm <- which(retfts==vartorm)
  #       retfts <- retfts[-idxtorm]
  #     }else{
  #       key <- FALSE
  #     }
  #   }
  # }
  
  retfts <- pfts
  key <- TRUE
  while(key){
    tmodel <- temptrain(ptrain,retfts,paste0("DGS.SEL.BIN.",pcri))
    retvif <- regclass::VIF(tmodel)
    if (length(which(retvif>10))>0){
      vartorm <- names(retvif)[which.max(retvif)]
      idxtorm <- which(retfts==vartorm)
      retfts <- retfts[-idxtorm]
    }else{
      key <- FALSE
    }
  }
  
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
  
  ovecrange <- c(1:floor(length(pOvec)/2))
  mvecrange <- c(1:floor((length(pMvec)*3)/4))
  
  keyvec <- c(1:10)
  
  rmidx <- vector()
  
  if(neipar$oidx==ovecrange[length(ovecrange)]){
    rmidx <- c(rmidx,c(1,5,6))
  }
  if(neipar$oidx==1){
    rmidx <- c(rmidx,c(2,7,8))
  }
  if(neipar$midx==mvecrange[length(mvecrange)]){
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
    neipar$oidx <- sample(ovecrange,size = 1)
    neipar$midx <- sample(mvecrange,size = 1)
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
                       pfilvec,
                       gofretlist,ddmretlist,tgofretlist,
                       param,goffeats,metafeats,
                       pcache,pfts = NULL,pftnum=NULL){
  
  isHit <- FALSE
  
  if(is.null(pfts)){
    isHit <- any(paste0(param$oidx,".",param$midx) %in% names(pcache))
  }else{
    isHit <- any(paste0(param$oidx,".",param$midx,'.',pftnum) %in% names(pcache))
  }
  
  if (isHit){
    ret <- list()
    if(is.null(pfts)){
      ret$val <- pcache[[paste0(param$oidx,".",param$midx)]]
    }else{
      ret$val <- pcache[[paste0(param$oidx,".",param$midx,'.',pftnum)]]
    }
    ret$cache <- pcache
    
    return(ret)
  }else{
    resultvec <- vector()
    DGSres <- NULL
    
    for (i in 1:length(plist)) {
      picked <- c(i)
      
      DGSretdf <- execDGS(pgofldf = gofretlist[[i]],pddmodf = ddmretlist[[i]],ptgofdf = tgofretlist[[i]],
                          goflist = plist[-picked],gofnlist = pnlist[-picked],
                          gofddmlist = pgofddmlist[-picked],metadf = pmetalist[-picked],
                          pcri = pcri,totpert = ppert,
                          goffts = goffeats,metafts = metafeats,
                          psearch = "GIVEN",givenparam = param,givenfts = pfts,
                          tmetadf = pmetalist[picked],tgofddmlist = pgofddmlist[picked],
                          tgoflist = plist[picked])
      
      if(is.null(DGSres)){
        DGSres <- DGSretdf
      }else{
        DGSres <- rbind(DGSres,DGSretdf)
      }
      
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
    
    filvec <- pfilvec
    
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
    if(is.null(pfts)){
      pcache[[paste0(param$oidx,".",param$midx)]] <- median(retdf[,pcri])
    }else{
      pcache[[paste0(param$oidx,".",param$midx,'.',pftnum)]] <- median(retdf[,pcri])
    }
    ret$avg <- mean(retdf[,pcri])
    ret$cache <- pcache
    ret$acc <- (ddmacc + gofacc)/2
    ret$filtnum <- filtnum
    
    return(ret)
  }
}