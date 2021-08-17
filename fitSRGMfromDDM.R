fitSRGMfromDDM <- function(df,pEst,
                           modelvec,
                           totpert,
                           isFCTBF){
    
  maxiter <- 1000
  
  numrows <- nrow(df)
  
  numrowssubs <- list()
  
  if(totpert!=0){
    for(i in 1:totpert){
      pert <- 0.1*((10-totpert-1)+i)
      tnumrows <- floor(numrows * pert)
      numrowssubs[[i]] <- tnumrows
    }
  }
  
  dflearn <- NULL
  
  if(isFCTBF=="FC"){
    dflearn <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    dflearn <- 1/tobsvd
  }
  
  retlist <- list()
  
  for(j in 1:totpert){
    
    currpert <- paste0("P",as.character((10-totpert-1)+j))
    
    retlist[[currpert]] <- list()
    retlist[[currpert]]$trainH <- numrowssubs[[j]]
    
    print(currpert)
    
    for(i in 1:length(modelvec)){
      
      currmodel <- modelvec[i] #current model name
      
      retparam <- NULL
      estimated <- NULL
      tbfestimated <- NULL
      frestimated <- NULL
      
      numparam <- 0
      
      dflearnsub <- NULL
      if(isFCTBF=="FC"){
        ddmest <- pEst[[currpert]]$SVR$EstElap
        dflearnsub <- c(dflearn[1:numrowssubs[[j]]],ddmest[numrowssubs[[j]]:numrows])
      }else{
        ddmest <- pEst[[currpert]]$SVR$EstFr
        dflearnsub <- 1/ddmest
      }
      
      pcfg <- list()
      
      pcfg$miter <- maxiter
      pcfg$algo <- "Best"
      
      funcForm <- NULL
      funcFRForm <- NULL
      residForm <- NULL
      formulaThis <- NULL
      parmEstThis <- NULL
      
      if (modelvec[i]=="GO") {
        pcfg$parInit <- parInitGO()
        pcfg$parLower <- parLowerGO()
        pcfg$parUpper <- parUpperGO()
        pcfg$parGrid <- parGridGO()
        
        funcForm <- funcGO
        funcFRForm <- funcFRGO
        residForm <- residFuncGO
        formulaThis <- formulaGO
        
        numparam <- 2
        
      }else if(modelvec[i]=="GG"){
        pcfg$parInit <- parInitGG()
        pcfg$parLower <- parLowerGG()
        pcfg$parUpper <- parUpperGG()
        pcfg$parGrid <- parGridGG()
        
        funcForm <- funcGG
        funcFRForm <- funcFRGG
        residForm <- residFuncGG
        formulaThis <- formulaGG
        
        numparam <- 3
        
      }else if(modelvec[i]=="Gompz"){
        pcfg$parInit <- parInitGompz()
        pcfg$parLower <- parLowerGompz()
        pcfg$parUpper <- parUpperGompz()
        pcfg$parGrid <- parGridGompz()
        
        funcForm <- funcGompz
        funcFRForm <- funcFRGompz
        residForm <- residFuncGompz
        formulaThis <- formulaGompz
        
        numparam <- 3
        
      }else if(modelvec[i]=="ISS"){
        pcfg$parInit <- parInitISS()
        pcfg$parLower <- parLowerISS()
        pcfg$parUpper <- parUpperISS()
        pcfg$parGrid <- parGridISS()
        
        funcForm <- funcISS
        funcFRForm <- funcFRISS
        residForm <- residFuncISS
        formulaThis <- formulaISS
        
        numparam <- 3
        
      }else if(modelvec[i]=="MD"){
        pcfg$parInit <- parInitMD()
        pcfg$parLower <- parLowerMD()
        pcfg$parUpper <- parUpperMD()
        pcfg$parGrid <- parGridMD()
        
        funcForm <- funcMD
        funcFRForm <- funcFRMD
        residForm <- residFuncMD
        formulaThis <- formulaMD
        
        numparam <- 3
        
      }else if(modelvec[i]=="MO"){
        pcfg$parInit <- parInitMO()
        pcfg$parLower <- parLowerMO()
        pcfg$parUpper <- parUpperMO()
        pcfg$parGrid <- parGridMO()
        
        funcForm <- funcMO
        funcFRForm <- funcFRMO
        residForm <- residFuncMO
        formulaThis <- formulaMO
        
        numparam <- 2
        
      }else if(modelvec[i]=="YID1"){
        pcfg$parInit <- parInitYID1()
        pcfg$parLower <- parLowerYID1()
        pcfg$parUpper <- parUpperYID1()
        pcfg$parGrid <- parGridYID1()
        
        funcForm <- funcYID1
        funcFRForm <- funcFRYID1
        residForm <- residFuncYID1
        formulaThis <- formulaYID1
        
        numparam <- 3
        
      }else if(modelvec[i]=="YID2"){
        pcfg$parInit <- parInitYID2()
        pcfg$parLower <- parLowerYID2()
        pcfg$parUpper <- parUpperYID2()
        pcfg$parGrid <- parGridYID2()
        
        funcForm <- funcYID2
        funcFRForm <- funcFRYID2
        residForm <- residFuncYID2
        formulaThis <- formulaYID2
        
        numparam <- 3
        
      }else if(modelvec[i]=="DSS"){
        pcfg$parInit <- parInitDSS()
        pcfg$parLower <- parLowerDSS()
        pcfg$parUpper <- parUpperDSS()
        pcfg$parGrid <- parGridDSS()
        
        funcForm <- funcDSS
        funcFRForm <- funcFRDSS
        residForm <- residFuncDSS
        formulaThis <- formulaDSS
        
        numparam <- 2
        
      }else if(modelvec[i]=="PNZ"){
        pcfg$parInit <- parInitPNZ()
        pcfg$parLower <- parLowerPNZ()
        pcfg$parUpper <- parUpperPNZ()
        pcfg$parGrid <- parGridPNZ()
        
        funcForm <- funcPNZ
        funcFRForm <- funcFRPNZ
        residForm <- residFuncPNZ
        formulaThis <- formulaPNZ
        
        numparam <- 4
        
      }else if(modelvec[i]=="PZ"){
        pcfg$parInit <- parInitPZ()
        pcfg$parLower <- parLowerPZ()
        pcfg$parUpper <- parUpperPZ()
        pcfg$parGrid <- parGridPZ()
        
        funcForm <- funcPZ
        funcFRForm <- funcFRPZ
        residForm <- residFuncPZ
        formulaThis <- formulaPZ
        
        numparam <- 5
        
      }else if(modelvec[i]=="PZI"){
        pcfg$parInit <- parInitPZI()
        pcfg$parLower <- parLowerPZI()
        pcfg$parUpper <- parUpperPZI()
        pcfg$parGrid <- parGridPZI()
        
        funcForm <- funcPZI
        funcFRForm <- funcFRPZI
        residForm <- residFuncPZI
        formulaThis <- formulaPZI
        
        numparam <- 3
        
      }else if(modelvec[i]=="Logi"){
        pcfg$parInit <- parInitLogi()
        pcfg$parLower <- parLowerLogi()
        pcfg$parUpper <- parUpperLogi()
        pcfg$parGrid <- parGridLogi()
        
        funcForm <- funcLogi
        funcFRForm <- funcFRLogi
        residForm <- residFuncLogi
        formulaThis <- formulaLogi
        
        numparam <- 3
      }else if(modelvec[i]=="JM"){
        pcfg$parInit <- parInitJM()
        pcfg$parLower <- parLowerJM()
        pcfg$parUpper <- parUpperJM()
        pcfg$parGrid <- parGridJM()
        
        funcForm <- funcJM
        residForm <- residFuncJM
        formulaThis <- formulaJM
        parmEstThis <- parmEstJM
        
        numparam <- 2
      }else if(modelvec[i]=="GEO"){
        pcfg$parInit <- parInitGEO()
        pcfg$parLower <- parLowerGEO()
        pcfg$parUpper <- parUpperGEO()
        pcfg$parGrid <- parGridGEO()
        
        funcForm <- funcGEO
        residForm <- residFuncGEO
        formulaThis <- formulaGEO
        parmEstThis <- parmEstGEO
        
        numparam <- 2
      }else if(modelvec[i]=="MB"){
        pcfg$parInit <- parInitMB()
        pcfg$parLower <- parLowerMB()
        pcfg$parUpper <- parUpperMB()
        pcfg$parGrid <- parGridMB()
        
        funcForm <- funcMB
        residForm <- residFuncMB
        formulaThis <- formulaMB
        
        numparam <- 2
      }else if(modelvec[i]=="LV"){
        pcfg$parInit <- parInitLV()
        pcfg$parLower <- parLowerLV()
        pcfg$parUpper <- parUpperLV()
        pcfg$parGrid <- parGridLV()
        
        funcForm <- funcLV
        residForm <- residFuncLV
        formulaThis <- formulaLV
        
        numparam <- 3
      }else{
        stopifnot(FALSE)
      }
      
      parmlist <- NULL
      
      if(isFCTBF=="FC"){
        retparam <- paramEstimation(dflearnsub,pcfg,residForm,formulaThis)
        if(!is.null(retparam)){
          parmlist <- retparam
        }
      }else if(isFCTBF=="TBF"){
        retparam <- paramEstimation(dflearnsub,pcfg,residForm,formulaThis)
        if(!is.null(retparam)){
          parmlist <- retparam
        }
      }
      
      if(is.null(parmlist)){
        if(isFCTBF=="FC"){
          estimated <- rep(NA,times = numrows)
          frestimated <- rep(NA,times = numrows)
        }else if(isFCTBF=="TBF"){
          estimated <- rep(NA,times = numrows)
          frestimated <- rep(NA,times = numrows)
          tbfestimated <- rep(NA,times = numrows)
        }
      }else{
        if(isFCTBF=="FC"){
          estimated <- funcForm(parmlist,df$i)
          frestimated <- funcFRForm(parmlist,df$i)
          
          estimated <- ifelse(estimated<0,0,estimated)
          frestimated <- ifelse(frestimated<0,0,frestimated)
          
        }else if(isFCTBF=="TBF"){
          frestimated <- funcForm(parmlist,df$i)
          frestimated <- ifelse(frestimated<0,0,frestimated)
          
          if(modelvec[i]=="JM"){
            if(parmlist$a<length(frestimated)){
              for(k in ceiling(parmlist$a):length(frestimated)){
                frestimated[k] <- 0
              }
            }
          }
          if(modelvec[i]=="MB"){
            if(parmlist$a<length(frestimated)){
              for(k in ceiling(parmlist$a):length(frestimated)){
                frestimated[k] <- 0
              }
            }
          }
          
          tbfestimated <- vector()
          estimated <- vector()
          for(k in 1:length(frestimated)){
            if(frestimated[k]>0){
              tbfestimated[k] <- 1/frestimated[k]
              estimated[k] <- sum(tbfestimated[1:k])
            }else{
              tbfestimated[k] <- Inf
              estimated[k] <- Inf
            }
          }
          
          estimated <- ifelse(estimated<0,0,estimated)
          tbfestimated <- ifelse(tbfestimated<0,0,tbfestimated)
        }
      }
      
      retlist[[currpert]][[modelvec[i]]] <- list()
      
      if(isFCTBF=="FC"){
        retlist[[currpert]][[modelvec[i]]]$NumParam <- numparam
        retlist[[currpert]][[modelvec[i]]]$EstElap <- estimated
        retlist[[currpert]][[modelvec[i]]]$EstFr <- frestimated
      }else if(isFCTBF=="TBF"){
        retlist[[currpert]][[modelvec[i]]]$NumParam <- numparam
        retlist[[currpert]][[modelvec[i]]]$EstElap <- estimated
        retlist[[currpert]][[modelvec[i]]]$EstFr <- tbfestimated
      }
      
    }
  }
  
  return(retlist)
}

fitSRGMfromDDMSingle <- function(df,pEst,
                           modelvec,
                           trainh,
                           predh){
  
  maxiter <- 1000
  
  dflearn <- df$n
  
  retlist <- list()
    
  for(i in 1:length(modelvec)){
    
    currmodel <- modelvec[i] #current model name
    
    retparam <- NULL
    estimated <- NULL
    tbfestimated <- NULL
    frestimated <- NULL
    
    numparam <- 0
    
    dflearnsub <- NULL
    if(isFCTBF=="FC"){
      ddmest <- pEst$SVR$EstElap
      dflearnsub <- c(dflearn[1:numrows],ddmest[(numrows+1):trainH])
    }else{
      ddmest <- pEst$SVR$EstFr
      dflearnsub <- 1/ddmest
    }
    
    pcfg <- list()
    
    pcfg$miter <- maxiter
    pcfg$algo <- "Best"
    
    funcForm <- NULL
    funcFRForm <- NULL
    residForm <- NULL
    formulaThis <- NULL
    parmEstThis <- NULL
    
    if (modelvec[i]=="GO") {
      pcfg$parInit <- parInitGO()
      pcfg$parLower <- parLowerGO()
      pcfg$parUpper <- parUpperGO()
      pcfg$parGrid <- parGridGO()
      
      funcForm <- funcGO
      funcFRForm <- funcFRGO
      residForm <- residFuncGO
      formulaThis <- formulaGO
      
      numparam <- 2
      
    }else if(modelvec[i]=="GG"){
      pcfg$parInit <- parInitGG()
      pcfg$parLower <- parLowerGG()
      pcfg$parUpper <- parUpperGG()
      pcfg$parGrid <- parGridGG()
      
      funcForm <- funcGG
      funcFRForm <- funcFRGG
      residForm <- residFuncGG
      formulaThis <- formulaGG
      
      numparam <- 3
      
    }else if(modelvec[i]=="Gompz"){
      pcfg$parInit <- parInitGompz()
      pcfg$parLower <- parLowerGompz()
      pcfg$parUpper <- parUpperGompz()
      pcfg$parGrid <- parGridGompz()
      
      funcForm <- funcGompz
      funcFRForm <- funcFRGompz
      residForm <- residFuncGompz
      formulaThis <- formulaGompz
      
      numparam <- 3
      
    }else if(modelvec[i]=="ISS"){
      pcfg$parInit <- parInitISS()
      pcfg$parLower <- parLowerISS()
      pcfg$parUpper <- parUpperISS()
      pcfg$parGrid <- parGridISS()
      
      funcForm <- funcISS
      funcFRForm <- funcFRISS
      residForm <- residFuncISS
      formulaThis <- formulaISS
      
      numparam <- 3
      
    }else if(modelvec[i]=="MD"){
      pcfg$parInit <- parInitMD()
      pcfg$parLower <- parLowerMD()
      pcfg$parUpper <- parUpperMD()
      pcfg$parGrid <- parGridMD()
      
      funcForm <- funcMD
      funcFRForm <- funcFRMD
      residForm <- residFuncMD
      formulaThis <- formulaMD
      
      numparam <- 3
      
    }else if(modelvec[i]=="MO"){
      pcfg$parInit <- parInitMO()
      pcfg$parLower <- parLowerMO()
      pcfg$parUpper <- parUpperMO()
      pcfg$parGrid <- parGridMO()
      
      funcForm <- funcMO
      funcFRForm <- funcFRMO
      residForm <- residFuncMO
      formulaThis <- formulaMO
      
      numparam <- 2
      
    }else if(modelvec[i]=="YID1"){
      pcfg$parInit <- parInitYID1()
      pcfg$parLower <- parLowerYID1()
      pcfg$parUpper <- parUpperYID1()
      pcfg$parGrid <- parGridYID1()
      
      funcForm <- funcYID1
      funcFRForm <- funcFRYID1
      residForm <- residFuncYID1
      formulaThis <- formulaYID1
      
      numparam <- 3
      
    }else if(modelvec[i]=="YID2"){
      pcfg$parInit <- parInitYID2()
      pcfg$parLower <- parLowerYID2()
      pcfg$parUpper <- parUpperYID2()
      pcfg$parGrid <- parGridYID2()
      
      funcForm <- funcYID2
      funcFRForm <- funcFRYID2
      residForm <- residFuncYID2
      formulaThis <- formulaYID2
      
      numparam <- 3
      
    }else if(modelvec[i]=="DSS"){
      pcfg$parInit <- parInitDSS()
      pcfg$parLower <- parLowerDSS()
      pcfg$parUpper <- parUpperDSS()
      pcfg$parGrid <- parGridDSS()
      
      funcForm <- funcDSS
      funcFRForm <- funcFRDSS
      residForm <- residFuncDSS
      formulaThis <- formulaDSS
      
      numparam <- 2
      
    }else if(modelvec[i]=="PNZ"){
      pcfg$parInit <- parInitPNZ()
      pcfg$parLower <- parLowerPNZ()
      pcfg$parUpper <- parUpperPNZ()
      pcfg$parGrid <- parGridPNZ()
      
      funcForm <- funcPNZ
      funcFRForm <- funcFRPNZ
      residForm <- residFuncPNZ
      formulaThis <- formulaPNZ
      
      numparam <- 4
      
    }else if(modelvec[i]=="PZ"){
      pcfg$parInit <- parInitPZ()
      pcfg$parLower <- parLowerPZ()
      pcfg$parUpper <- parUpperPZ()
      pcfg$parGrid <- parGridPZ()
      
      funcForm <- funcPZ
      funcFRForm <- funcFRPZ
      residForm <- residFuncPZ
      formulaThis <- formulaPZ
      
      numparam <- 5
      
    }else if(modelvec[i]=="PZI"){
      pcfg$parInit <- parInitPZI()
      pcfg$parLower <- parLowerPZI()
      pcfg$parUpper <- parUpperPZI()
      pcfg$parGrid <- parGridPZI()
      
      funcForm <- funcPZI
      funcFRForm <- funcFRPZI
      residForm <- residFuncPZI
      formulaThis <- formulaPZI
      
      numparam <- 3
      
    }else if(modelvec[i]=="Logi"){
      pcfg$parInit <- parInitLogi()
      pcfg$parLower <- parLowerLogi()
      pcfg$parUpper <- parUpperLogi()
      pcfg$parGrid <- parGridLogi()
      
      funcForm <- funcLogi
      funcFRForm <- funcFRLogi
      residForm <- residFuncLogi
      formulaThis <- formulaLogi
      
      numparam <- 3
    }else if(modelvec[i]=="JM"){
      pcfg$parInit <- parInitJM()
      pcfg$parLower <- parLowerJM()
      pcfg$parUpper <- parUpperJM()
      pcfg$parGrid <- parGridJM()
      
      funcForm <- funcJM
      residForm <- residFuncJM
      formulaThis <- formulaJM
      parmEstThis <- parmEstJM
      
      numparam <- 2
    }else if(modelvec[i]=="GEO"){
      pcfg$parInit <- parInitGEO()
      pcfg$parLower <- parLowerGEO()
      pcfg$parUpper <- parUpperGEO()
      pcfg$parGrid <- parGridGEO()
      
      funcForm <- funcGEO
      residForm <- residFuncGEO
      formulaThis <- formulaGEO
      parmEstThis <- parmEstGEO
      
      numparam <- 2
    }else if(modelvec[i]=="MB"){
      pcfg$parInit <- parInitMB()
      pcfg$parLower <- parLowerMB()
      pcfg$parUpper <- parUpperMB()
      pcfg$parGrid <- parGridMB()
      
      funcForm <- funcMB
      residForm <- residFuncMB
      formulaThis <- formulaMB
      
      numparam <- 2
    }else if(modelvec[i]=="LV"){
      pcfg$parInit <- parInitLV()
      pcfg$parLower <- parLowerLV()
      pcfg$parUpper <- parUpperLV()
      pcfg$parGrid <- parGridLV()
      
      funcForm <- funcLV
      residForm <- residFuncLV
      formulaThis <- formulaLV
      
      numparam <- 3
    }else{
      stopifnot(FALSE)
    }
    
    parmlist <- NULL
    
    if(isFCTBF=="FC"){
      retparam <- paramEstimation(dflearnsub,pcfg,residForm,formulaThis)
      if(!is.null(retparam)){
        parmlist <- retparam
      }
    }else if(isFCTBF=="TBF"){
      retparam <- paramEstimation(dflearnsub,pcfg,residForm,formulaThis)
      if(!is.null(retparam)){
        parmlist <- retparam
      }
    }
    
    if(is.null(parmlist)){
      if(isFCTBF=="FC"){
        estimated <- rep(NA,times = trainH)
        frestimated <- rep(NA,times = trainH)
      }else if(isFCTBF=="TBF"){
        estimated <- rep(NA,times = trainH)
        frestimated <- rep(NA,times = trainH)
        tbfestimated <- rep(NA,times = trainH)
      }
    }else{
      if(isFCTBF=="FC"){
        estimated <- funcForm(parmlist,c(1:trainH))
        frestimated <- funcFRForm(parmlist,c(1:trainH))
        
        estimated <- ifelse(estimated<0,0,estimated)
        frestimated <- ifelse(frestimated<0,0,frestimated)
        
      }else if(isFCTBF=="TBF"){
        frestimated <- funcForm(parmlist,c(1:trainH))
        frestimated <- ifelse(frestimated<0,0,frestimated)
        
        if(modelvec[i]=="JM"){
          if(parmlist$a<length(frestimated)){
            for(k in ceiling(parmlist$a):length(frestimated)){
              frestimated[k] <- 0
            }
          }
        }
        if(modelvec[i]=="MB"){
          if(parmlist$a<length(frestimated)){
            for(k in ceiling(parmlist$a):length(frestimated)){
              frestimated[k] <- 0
            }
          }
        }
        
        tbfestimated <- vector()
        estimated <- vector()
        for(k in 1:length(frestimated)){
          if(frestimated[k]>0){
            tbfestimated[k] <- 1/frestimated[k]
            estimated[k] <- sum(tbfestimated[1:k])
          }else{
            tbfestimated[k] <- Inf
            estimated[k] <- Inf
          }
        }
        
        estimated <- ifelse(estimated<0,0,estimated)
        tbfestimated <- ifelse(tbfestimated<0,0,tbfestimated)
      }
    }
    
    retlist[[modelvec[i]]] <- list()
    
    if(isFCTBF=="FC"){
      retlist[[modelvec[i]]]$NumParam <- numparam
      retlist[[modelvec[i]]]$EstElap <- estimated
      retlist[[modelvec[i]]]$EstFr <- frestimated
    }else if(isFCTBF=="TBF"){
      retlist[[modelvec[i]]]$NumParam <- numparam
      retlist[[modelvec[i]]]$EstElap <- estimated
      retlist[[modelvec[i]]]$EstFr <- tbfestimated
    }
    
  }
  
  return(retlist)
}


getSRGMfromDDMGOF <- function(df,pEst,pDDMEst,
                     modelvec,
                     crivec,
                     colnmleft,
                     totpert,
                     isFCTBF){
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  
  leftmatrix <- matrix(nrow=(length(modelvec)*totpert),ncol=length(colnmleft))
  rightmatrix <- matrix(nrow=(length(modelvec)*totpert),ncol=(length(crivec)+1))
  
  ccriname <- paste0("F",crivec)
  
  observed <- NULL
  tbfobserved <- NULL
  frobserved <- NULL
  
  if(isFCTBF=="FC"){
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    frobserved <- 1/tobsvd
    tbfobserved <- tobsvd
    observed <- vector()
    for(i in 1:length(tobsvd)){
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  numrows <- nrow(df)
  
  for(i in 1:length(modelvec)){
    
    currmodel <- modelvec[i] #current model name
    
    for(j in 1:totpert){
      
      currpert <- paste0("P",as.character((10-totpert-1)+j))
      
      gofobsrvd <- pDDMEst[[currpert]]$SVR$EstElap
      
      estimated <- NULL
      tbfestimated <- NULL
      frestimated <- NULL
      numparam <- 0
      
      pEstsub <- pEst[[currpert]][[currmodel]]
      trainH <- pEst[[currpert]]$trainH
      estimated <- pEstsub$EstElap
      if(isFCTBF=="FC"){
        frestimated <- pEstsub$EstFr
      }else if(isFCTBF=="TBF"){
        tbfestimated <- pEstsub$EstFr
      }
      numparam <- pEstsub$NumParam
      
      gofobsrvd <- c(observed[1:trainH],gofobsrvd[(trainH+1):numrows])

      stopifnot(length(observed)==length(gofobsrvd))
      
      if(is.null(numparam)){
        print(currmodel)
        print(currpert)
      }
      
      currrowidx  <- j+((i-1)*totpert)
      
      leftmatrix[currrowidx,1] <- currmodel
      leftmatrix[currrowidx,2] <- currpert
      if(isFCTBF=="FC"){
        # rightmatrix[currrowidx,1] <- FC.MSE((numrows-(trainH+1)),1,gofobsrvd[(trainH+1):numrows],estimated[(trainH+1):numrows])
        rightmatrix[currrowidx,1] <- FC.MSE(numrows,1,gofobsrvd,estimated)
        rightmatrix[currrowidx,2] <- FC.EP(numrows,observed,estimated)
        rightmatrix[currrowidx,3] <- FC.MEOP(trainH,numrows,observed,estimated)
      }else if(isFCTBF=="TBF"){
        rightmatrix[currrowidx,1] <- TBF.PE(numrows,gofobsrvd,estimated)
        rightmatrix[currrowidx,2] <- TBF.AE(trainH,numrows,gofobsrvd,estimated)
        rightmatrix[currrowidx,3] <- TBF.PE(numrows,observed,estimated)
        rightmatrix[currrowidx,4] <- TBF.AE(trainH,numrows,observed,estimated)
      }
    }
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnmleft
  
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- c('MSE',crivec)
  
  tnamevec <- c('MSE',crivec)
  
  for (i in 1:length(tnamevec)) {
    rightdf[tnamevec[i]] <- as.numeric(rightdf[,tnamevec[i]]) 
  }
  
  retinstance <- cbind(leftdf,rightdf)
  
  retinstance <- na.omit(retinstance)
  
  nmzretinstance <- NULL
  
  for(i in 1:totpert){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    subdf <- retinstance[retinstance$Pert==currpert,]
    for(j in 1:length(crivec)){
      
      valvec  <- subdf[,crivec[j]]
      
      if(length(valvec)>1){
        
        infidx <- which(is.infinite(valvec))
        
        noinfsub <- NULL
        if(length(infidx)!=0){
          noinfsub <- valvec[-c(infidx)]
        }else{
          noinfsub <- valvec
        }
        
        minval <- min(noinfsub)
        maxval <- max(noinfsub)
        
        denominator <- maxval-minval
        
        nmzvec <- vector()
        for (k in 1:length(valvec)) {
          if(is.infinite(valvec[k])){
            nmzvec[k] <- 0
          }else{
            if(denominator==0){
              nmzvec[k] <- 1
            }else{
              nmzvec[k] <- 0.1 + ((maxval - valvec[k])*0.9)/denominator
            }
          }
        }
        subdf[paste0("N",crivec[j])] <- nmzvec
      }else if(length(valvec)==1){
        subdf[paste0("N",crivec[j])] <- c(1)
      }
      
    }
    
    if(is.null(nmzretinstance)){
      nmzretinstance <- subdf
    }else{
      nmzretinstance <- rbind(nmzretinstance,subdf)
    }
  }
  
  return(nmzretinstance)
}

getSRGMfromDDMGOFSingle <- function(df,pEst,pDDMEst,
                              modelvec,
                              crivec,
                              colnmleft,
                              isFCTBF){
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  
  leftmatrix <- matrix(nrow=length(modelvec),ncol=length(colnmleft))
  rightmatrix <- matrix(nrow=length(modelvec),ncol=1)
  
  observed <- NULL
  tbfobserved <- NULL
  frobserved <- NULL
  
  if(isFCTBF=="FC"){
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    frobserved <- 1/tobsvd
    tbfobserved <- tobsvd
    observed <- vector()
    for(i in 1:length(tobsvd)){
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  numrows <- nrow(df)
  
  trainH <- pEst$trainH
  
  for(i in 1:length(modelvec)){
    
    currmodel <- modelvec[i] #current model name
      
    gofobsrvd <- pDDMEst$SVR$EstElap
    
    estimated <- NULL
    tbfestimated <- NULL
    frestimated <- NULL
    numparam <- 0
    
    pEstsub <- pEst[[currmodel]]
    estimated <- pEstsub$EstElap
    
    if(isFCTBF=="FC"){
      frestimated <- pEstsub$EstFr
    }else if(isFCTBF=="TBF"){
      tbfestimated <- pEstsub$EstFr
    }
    numparam <- pEstsub$NumParam
    
    gofobsrvd <- c(observed[1:trainH],gofobsrvd[(trainH+1):numrows])
    
    currrowidx  <- i
    
    leftmatrix[currrowidx,1] <- currmodel
    leftmatrix[currrowidx,2] <- "None"
    if(isFCTBF=="FC"){
      rightmatrix[currrowidx,1] <- FC.MSE(numrows,1,gofobsrvd,estimated)
    }else if(isFCTBF=="TBF"){
      rightmatrix[currrowidx,1] <- TBF.PE(trainH,gofobsrvd,estimated)
      rightmatrix[currrowidx,2] <- TBF.AE((numrows+1),trainH,gofobsrvd,estimated)
    }
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnmleft
  
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- c('MSE')
  
  tnamevec <- c('MSE')
  
  for (i in 1:length(tnamevec)) {
    rightdf[tnamevec[i]] <- as.numeric(rightdf[,tnamevec[i]]) 
  }
  
  retinstance <- cbind(leftdf,rightdf)
  
  retinstance <- na.omit(retinstance)
  
  return(retinstance)
}