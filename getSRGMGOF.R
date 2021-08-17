getSRGMEst <- function(df,
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
    
    for(i in 1:length(modelvec)){
      
      currmodel <- modelvec[i] #current model name

      retparam <- NULL
      estimated <- NULL
      tbfestimated <- NULL
      frestimated <- NULL
      
      numparam <- 0
      
      dflearnsub <- dflearn[1:numrowssubs[[j]]]
      
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








getSRGMEstfive <- function(df,
                       modelvec,
                       totpert,
                       isFCTBF){
  
  maxiter <- 1000
  pc <- 0.2
  
  numrows <- nrow(df)
  
  retlist <- list()
  
  if(numrows>=20){
    
    pertlist <- list()
    numrowslist <- list()
    
    for (i in 1:totpert) {
      trainH <- NULL
      tempstart <- NULL
      tempstart <- ((((10-totpert)-1)+i)*10)+5
      trainH <- floor((tempstart/100)*numrows)
      
      predH <- NULL
      predH <- numrows
      
      pertlist[[(length(pertlist)+1)]] <- paste0("T",as.character(tempstart))
      temphorizon <- c(trainH,predH)
      stopifnot(trainH<predH)
      numrowslist[[(length(numrowslist)+1)]] <- temphorizon
    }
  
    dflearn <- NULL
    
    if(isFCTBF=="FC"){
      dflearn <- df$n
    }else if(isFCTBF=="TBF"){
      tobsvd <- MovingAverage(df$f,7)
      dflearn <- 1/tobsvd
    }
    
    for(j in 1:length(pertlist)){
      
      currpert <- pertlist[[j]]
      
      retlist[[currpert]] <- list()
      retlist[[currpert]]$trainH <- numrowslist[[j]][1]
    
      for(i in 1:length(modelvec)){
      
        currmodel <- modelvec[i]
  
        retparam <- NULL
        estimated <- NULL
        frestimated <- NULL
        tbfestimated <- NULL
        
        numparam <- 0
        
        dflearnsub <- dflearn[1:numrowslist[[j]][1]]
        
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
  }
  
  return(retlist)
}

getSRGMEstOrig <- function(df,
                       modelvec,
                       trainh,
                       predh){
  
  maxiter <- 1000
  
  dflearn <- df$n
  
  retlist <- list()
  retlist$trainH <- trainh
  
  dflearnsub <- dflearn[1:trainh]

  for(i in 1:length(modelvec)){
    
    currmodel <- modelvec[i] #current model name
    
    retparam <- NULL
    estimated <- NULL
    frestimated <- NULL
    
    numparam <- 0
    
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
    
    retparam <- paramEstimation(dflearnsub,pcfg,residForm,formulaThis)
    if(!is.null(retparam)){
      parmlist <- retparam
    }
    
    if(is.null(parmlist)){
      estimated <- rep(NA,times = predh)
      frestimated <- rep(NA,times = predh)
    }else{
      estimated <- funcForm(parmlist,c(1:predh))
      frestimated <- funcFRForm(parmlist,c(1:predh))
      estimated <- ifelse(estimated<0,0,estimated)
      frestimated <- ifelse(frestimated<0,0,frestimated)
    }
    
    retlist[[modelvec[i]]] <- list()
    retlist[[modelvec[i]]]$NumParam <- numparam
    retlist[[modelvec[i]]]$EstElap <- estimated
    retlist[[modelvec[i]]]$EstFr <- frestimated
  }
  
  return(retlist)
}




#getSRGMGOF



getSRGMGOF<-function(df,pEst,
                     modelvec,
                     goffullvec,nrgofvec,
                     crivec,
                     colnmleft,
                     totpert,
                     isFCTBF){
  
  pc <- 0.2
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  lastmatrix <- NULL
  
  leftmatrix <- matrix(nrow=(length(modelvec)*totpert),ncol=length(colnmleft))
  rightmatrix <- matrix(nrow=(length(modelvec)*totpert),ncol=(length(goffullvec)-1))
  lastmatrix <- matrix(nrow=(length(modelvec)*totpert),ncol=1)
  
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
      
      if(is.null(numparam)){
        print(currmodel)
        print(currpert)
      }
      
      currrowidx  <- j+((i-1)*totpert)
      
      leftmatrix[currrowidx,1] <- currmodel
      leftmatrix[currrowidx,2] <- currpert
      if(isFCTBF=="FC"){
        rightmatrix[currrowidx,1] <- FC.MSE(trainH,numparam,observed,estimated)
        rightmatrix[currrowidx,2] <- FC.MAE(trainH,numparam,observed,estimated)
        rightmatrix[currrowidx,3] <- FC.Rsquare(trainH,observed,estimated)
        rightmatrix[currrowidx,4] <- FC.Noise(trainH,frestimated)
        rightmatrix[currrowidx,5] <- FC.Bias(1,trainH,observed,estimated)
        rightmatrix[currrowidx,6] <- FC.Variation(trainH,observed,estimated)
        rightmatrix[currrowidx,7] <- FC.PRR(trainH,observed,estimated)
        rightmatrix[currrowidx,8] <- FC.WLSE(trainH,observed,estimated,pc)
        rightmatrix[currrowidx,9] <- FC.EP(trainH,observed,estimated)
        rightmatrix[currrowidx,10] <- FC.MEOP(1,trainH,observed,estimated)
        rightmatrix[currrowidx,11] <- FC.TOUP(trainH,observed,estimated)
        rightmatrix[currrowidx,12] <- PRatio(numrows,trainH)
        rightmatrix[currrowidx,13] <- FC.EP(numrows,observed,estimated)
        rightmatrix[currrowidx,14] <- FC.MEOP(trainH,numrows,observed,estimated)
        rightmatrix[currrowidx,15] <- abs(FC.Bias(1,trainH,observed,estimated))
        lastmatrix[currrowidx,1] <- GROUP(FC.Bias(trainH,numrows,observed,estimated))
      }else if(isFCTBF=="TBF"){
        rightmatrix[currrowidx,1] <- TBF.MSE(trainH,numparam,observed,estimated)
        rightmatrix[currrowidx,2] <- TBF.MAE(trainH,numparam,observed,estimated)
        rightmatrix[currrowidx,3] <- TBF.Rsquare(trainH,observed,estimated)
        rightmatrix[currrowidx,4] <- TBF.Noise(trainH,tbfestimated)
        rightmatrix[currrowidx,5] <- TBF.Bias(1,trainH,observed,estimated)
        rightmatrix[currrowidx,6] <- TBF.Variation(trainH,observed,estimated)
        rightmatrix[currrowidx,7] <- TBF.PRR(trainH,observed,estimated)
        rightmatrix[currrowidx,8] <- TBF.WLSE(trainH,observed,estimated,pc)
        rightmatrix[currrowidx,9] <- TBF.PE(trainH,observed,estimated)
        rightmatrix[currrowidx,10] <- TBF.AE(1,trainH,observed,estimated)
        rightmatrix[currrowidx,11] <- TBF.TOUP(trainH,observed,estimated)
        rightmatrix[currrowidx,12] <- PRatio(numrows,trainH)
        rightmatrix[currrowidx,13] <- TBF.PE(numrows,observed,estimated)
        rightmatrix[currrowidx,14] <- TBF.AE(trainH,numrows,observed,estimated)
        rightmatrix[currrowidx,15] <- abs(TBF.Bias(1,trainH,observed,estimated))
        lastmatrix[currrowidx,1] <- GROUP(TBF.Bias(trainH,numrows,observed,estimated))
      }
    }
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnmleft
  
  rightdf <- data.frame(rightmatrix)
  
  lastdf <- data.frame(lastmatrix)
  
  rightdf <- cbind(rightdf,lastdf)
  names(rightdf) <- goffullvec
  
  retinstance <- cbind(leftdf,rightdf)
  
  retinstance <- na.omit(retinstance)
  
  retinstance$Group <- as.factor(retinstance$Group)
  
  nmzretinstance <- NULL
  
  for(i in 1:totpert){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    subdf <- retinstance[retinstance$Pert==currpert,]
    for(j in 1:length(nrgofvec)){
      
      valvec  <- subdf[,nrgofvec[j]]
      
      if(length(valvec)>1){
        
        infidx <- which(is.infinite(valvec))
        
        noinfsub <- NULL
        if(length(infidx)>0){
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
            if(nrgofvec[j]=="Rsquare"){
              nmzvec[k] <- 1
            }else{
              nmzvec[k] <- 0
            }
          }else{
            if(denominator==0){
              if(nrgofvec[j]=="Rsquare"){
                nmzvec[k] <- 0
              }else{
                nmzvec[k] <- 1
              }
            }else{
              if(nrgofvec[j]=="Rsquare"){
                nmzvec[k] <- 0 + ((maxval - valvec[k])*0.9)/denominator
              }else{
                nmzvec[k] <- 0.1 + ((maxval - valvec[k])*0.9)/denominator
              }
            }
          }
        }
        subdf[paste0("N",nrgofvec[j])] <- nmzvec
        
      }else if(length(valvec)==1){
        if(nrgofvec[j]=="Rsquare"){
          subdf[paste0("N",nrgofvec[j])] <- c(0)
        }else{
          subdf[paste0("N",nrgofvec[j])] <- c(1)
        }
      }
    }
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
        minval <- ifelse(minval<0.000001,0.000001,minval)
        
        ratvec <- vector()
        for (k in 1:length(valvec)) {
          if(is.infinite(valvec[k])){
            ratvec[k] <- Inf
          }else{
            if(minval==0){
              ratvec[k] <- valvec[k]
            }else{
              ratvec[k] <- (valvec[k] - minval) / minval
            }
          }
        }
        subdf[paste0("R",crivec[j])] <- ratvec
      }else if(length(valvec)==1){
        subdf[paste0("R",crivec[j])] <- c(0)
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







getSRGMGOFaddi<-function(df,pEst,pEstfive,
                         modelvec,
                         goffullvec,nrgofvec,
                         crivec,
                         colnmleft,
                         totpert,
                         isFCTBF){
  
  pc <- 0.2
  
  numrows <- nrow(df)
  
  pertlist <- list()
  
  observed <- NULL
  tbfobserved <- NULL
  frobserved <- NULL
  
  if(isFCTBF=="FC"){
    observed <- df$n
  }else if(isFCTBF=="TBF"){
    tobsvd <- MovingAverage(df$f,7)
    tbfobserved <- tobsvd
    frobserved <- 1/tobsvd
    observed <- vector()
    for(i in 1:length(tobsvd)){
      observed[i] <- sum(tobsvd[1:i])
    }
  }
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  addrows <- ((totpert-1)*totpert)/2
  leftmatrix <- matrix(nrow=(addrows*length(modelvec)),ncol=length(colnmleft))
  rightmatrix <- matrix(nrow=(addrows*length(modelvec)),ncol=(length(goffullvec)-1))
  lastmatrix <- matrix(nrow=(addrows*length(modelvec)),ncol=1)
  
  tempidx <- 1
  
  for(i in 1:(totpert-1)){
    currpert <- paste0("P",as.character((10-totpert-1)+i))
    pEstsub <- pEst[[currpert]]
    tempstart <- NULL
    tempstart <- (((10-totpert)-1)+i)
    trainH <- pEstsub$trainH
    
    for(j in i:(totpert-1)){
      
      predH <- NULL
      tempend <- NULL
      tempend <- ((10-totpert)+j)
      predH <- floor((tempend*0.1)*length(observed))
      
      addpert <- paste0("P",as.character(tempstart),as.character(tempend))
      
      pertlist[[(length(pertlist)+1)]] <- addpert
      
      for(k in 1:length(modelvec)){
        
        currmodel <- modelvec[k]
        
        estimated <- NULL
        frestimated <- NULL
        tbfestimated <- NULL
        numparam <- 0
        
        pEstsubsub <- pEstsub[[currmodel]]
        
        estimated <- pEstsubsub$EstElap
        if(isFCTBF=="FC"){
          frestimated <- pEstsubsub$EstFr
        }else if(isFCTBF=="TBF"){
          tbfestimated <- pEstsubsub$EstFr
        }
        numparam <- pEstsubsub$NumParam
        
        currrowidx <- tempidx
        
        leftmatrix[currrowidx,1] <- currmodel
        leftmatrix[currrowidx,2] <- addpert
        if(isFCTBF=="FC"){
          rightmatrix[currrowidx,1] <- FC.MSE(trainH,numparam,observed,estimated)
          rightmatrix[currrowidx,2] <- FC.MAE(trainH,numparam,observed,estimated)
          rightmatrix[currrowidx,3] <- FC.Rsquare(trainH,observed,estimated)
          rightmatrix[currrowidx,4] <- FC.Noise(trainH,frestimated)
          rightmatrix[currrowidx,5] <- FC.Bias(1,trainH,observed,estimated)
          rightmatrix[currrowidx,6] <- FC.Variation(trainH,observed,estimated)
          rightmatrix[currrowidx,7] <- FC.PRR(trainH,observed,estimated)
          rightmatrix[currrowidx,8] <- FC.WLSE(trainH,observed,estimated,pc)
          rightmatrix[currrowidx,9] <- FC.EP(trainH,observed,estimated)
          rightmatrix[currrowidx,10] <- FC.MEOP(1,trainH,observed,estimated)
          rightmatrix[currrowidx,11] <- FC.TOUP(trainH,observed,estimated)
          rightmatrix[currrowidx,12] <- PRatio(predH,trainH)
          rightmatrix[currrowidx,13] <- FC.EP(predH,observed,estimated)
          rightmatrix[currrowidx,14] <- FC.MEOP(trainH,predH,observed,estimated)
          rightmatrix[currrowidx,15] <- abs(FC.Bias(1,trainH,observed,estimated))
          lastmatrix[currrowidx,1] <- GROUP(FC.Bias(trainH,predH,observed,estimated))
        }else if(isFCTBF=="TBF"){
          rightmatrix[currrowidx,1] <- TBF.MSE(trainH,numparam,observed,estimated)
          rightmatrix[currrowidx,2] <- TBF.MAE(trainH,numparam,observed,estimated)
          rightmatrix[currrowidx,3] <- TBF.Rsquare(trainH,observed,estimated)
          rightmatrix[currrowidx,4] <- TBF.Noise(trainH,tbfestimated)
          rightmatrix[currrowidx,5] <- TBF.Bias(1,trainH,observed,estimated)
          rightmatrix[currrowidx,6] <- TBF.Variation(trainH,observed,estimated)
          rightmatrix[currrowidx,7] <- TBF.PRR(trainH,observed,estimated)
          rightmatrix[currrowidx,8] <- TBF.WLSE(trainH,observed,estimated,pc)
          rightmatrix[currrowidx,9] <- TBF.PE(trainH,observed,estimated)
          rightmatrix[currrowidx,10] <- TBF.AE(1,trainH,observed,estimated)
          rightmatrix[currrowidx,11] <- TBF.TOUP(trainH,observed,estimated)
          rightmatrix[currrowidx,12] <- PRatio(predH,trainH)
          rightmatrix[currrowidx,13] <- TBF.PE(predH,observed,estimated)
          rightmatrix[currrowidx,14] <- TBF.AE(trainH,predH,observed,estimated)
          rightmatrix[currrowidx,15] <- abs(TBF.Bias(1,trainH,observed,estimated))
          lastmatrix[currrowidx,1] <- GROUP(TBF.Bias(trainH,predH,observed,estimated))
        }
        
        tempidx <- tempidx + 1
      }
    }
  }
  
  leftmatrix2 <- NULL
  rightmatrix2 <- NULL
  lastmatrix2 <- NULL
  
  if(nrow(df)>=20){
    addrows <- (((totpert-1)*totpert)/2)+totpert
    leftmatrix2 <- matrix(nrow=(addrows*length(modelvec)),ncol=2)
    rightmatrix2 <- matrix(nrow=(addrows*length(modelvec)),ncol=(length(goffullvec)-1))
    lastmatrix2 <- matrix(nrow = (addrows*length(modelvec)),ncol = 1)
    
    tempidx <- 1
    
    for (i in 1:totpert) {
      currpert <- paste0("T",as.character((((10-totpert-1)+i)*10)+5))
      pertlist[[(length(pertlist)+1)]] <- currpert
      
      pEstsub <- pEstfive[[currpert]]
      trainH <- pEstsub$trainH
      
      for(j in 1:length(modelvec)){
        
        currmodel <- modelvec[j]
        
        estimated <- NULL
        frestimated <- NULL
        tbfestimated <- NULL
        numparam <- 0
        
        pEstsubsub <- pEstsub[[currmodel]]
        
        estimated <- pEstsubsub$EstElap
        if(isFCTBF=="FC"){
          frestimated <- pEstsubsub$EstFr
        }else if(isFCTBF=="TBF"){
          tbfestimated <- pEstsubsub$EstFr
        }
        numparam <- pEstsubsub$NumParam
        
        currrowidx <- tempidx
        
        leftmatrix2[currrowidx,1] <- currmodel
        leftmatrix2[currrowidx,2] <- currpert
        if(isFCTBF=="FC"){
          rightmatrix2[currrowidx,1] <- FC.MSE(trainH,numparam,observed,estimated)
          rightmatrix2[currrowidx,2] <- FC.MAE(trainH,numparam,observed,estimated)
          rightmatrix2[currrowidx,3] <- FC.Rsquare(trainH,observed,estimated)
          rightmatrix2[currrowidx,4] <- FC.Noise(trainH,frestimated)
          rightmatrix2[currrowidx,5] <- FC.Bias(1,trainH,observed,estimated)
          rightmatrix2[currrowidx,6] <- FC.Variation(trainH,observed,estimated)
          rightmatrix2[currrowidx,7] <- FC.PRR(trainH,observed,estimated)
          rightmatrix2[currrowidx,8] <- FC.WLSE(trainH,observed,estimated,pc)
          rightmatrix2[currrowidx,9] <- FC.EP(trainH,observed,estimated)
          rightmatrix2[currrowidx,10] <- FC.MEOP(1,trainH,observed,estimated)
          rightmatrix2[currrowidx,11] <- FC.TOUP(trainH,observed,estimated)
          rightmatrix2[currrowidx,12] <- PRatio(length(estimated),trainH)
          rightmatrix2[currrowidx,13] <- FC.EP(length(estimated),observed,estimated)
          rightmatrix2[currrowidx,14] <- FC.MEOP(trainH,length(estimated),observed,estimated)
          rightmatrix2[currrowidx,15] <- abs(FC.Bias(1,trainH,observed,estimated))
          lastmatrix2[currrowidx,1] <- GROUP(FC.Bias(trainH,length(estimated),observed,estimated))
        }else if(isFCTBF=="TBF"){
          rightmatrix2[currrowidx,1] <- TBF.MSE(trainH,numparam,observed,estimated)
          rightmatrix2[currrowidx,2] <- TBF.MAE(trainH,numparam,observed,estimated)
          rightmatrix2[currrowidx,3] <- TBF.Rsquare(trainH,observed,estimated)
          rightmatrix2[currrowidx,4] <- TBF.Noise(trainH,tbfestimated)
          rightmatrix2[currrowidx,5] <- TBF.Bias(1,trainH,observed,estimated)
          rightmatrix2[currrowidx,6] <- TBF.Variation(trainH,observed,estimated)
          rightmatrix2[currrowidx,7] <- TBF.PRR(trainH,observed,estimated)
          rightmatrix2[currrowidx,8] <- TBF.WLSE(trainH,observed,estimated,pc)
          rightmatrix2[currrowidx,9] <- TBF.PE(trainH,observed,estimated)
          rightmatrix2[currrowidx,10] <- TBF.AE(1,trainH,observed,estimated)
          rightmatrix2[currrowidx,11] <- TBF.TOUP(trainH,observed,estimated)
          rightmatrix2[currrowidx,12] <- PRatio(length(estimated),trainH)
          rightmatrix2[currrowidx,13] <- TBF.PE(length(estimated),observed,estimated)
          rightmatrix2[currrowidx,14] <- TBF.AE(trainH,length(estimated),observed,estimated)
          rightmatrix2[currrowidx,15] <- abs(TBF.Bias(1,trainH,observed,estimated))
          lastmatrix2[currrowidx,1] <- GROUP(TBF.Bias(trainH,length(estimated),observed,estimated))
        }
        tempidx <- tempidx + 1
      }
    }
    
    for(i in 1:(totpert-1)){
      currpert <- paste0("T",as.character((((10-totpert-1)+i)*10)+5))
      pEstsub <- pEstfive[[currpert]]
      tempstart <- NULL
      tempstart <- (((10-totpert-1)+i)*10)+5
      trainH <- pEstsub$trainH
      
      for(j in i:(totpert-1)){
        predH <- NULL
        tempend <- NULL
        tempend <- (((10-totpert)+j)*10)+5
        predH <- floor((tempend/100)*length(estimated))
        
        addpert <- paste0("T",as.character(tempstart),as.character(tempend))
        pertlist[[(length(pertlist)+1)]] <- addpert
        
        for(k in 1:length(modelvec)){
          
          currmodel <- modelvec[k]
          
          estimated <- NULL
          frestimated <- NULL
          tbfestimated <- NULL
          numparam <- 0
          
          pEstsubsub <- pEstsub[[currmodel]]
          
          estimated <- pEstsubsub$EstElap
          if(isFCTBF=="FC"){
            frestimated <- pEstsubsub$EstFr
          }else if(isFCTBF=="TBF"){
            tbfestimated <- pEstsubsub$EstFr
          }
          numparam <- pEstsubsub$NumParam
          
          currrowidx <- tempidx
          
          leftmatrix2[currrowidx,1] <- currmodel
          leftmatrix2[currrowidx,2] <- addpert
          if(isFCTBF=="FC"){
            rightmatrix2[currrowidx,1] <- FC.MSE(trainH,numparam,observed,estimated)
            rightmatrix2[currrowidx,2] <- FC.MAE(trainH,numparam,observed,estimated)
            rightmatrix2[currrowidx,3] <- FC.Rsquare(trainH,observed,estimated)
            rightmatrix2[currrowidx,4] <- FC.Noise(trainH,frestimated)
            rightmatrix2[currrowidx,5] <- FC.Bias(1,trainH,observed,estimated)
            rightmatrix2[currrowidx,6] <- FC.Variation(trainH,observed,estimated)
            rightmatrix2[currrowidx,7] <- FC.PRR(trainH,observed,estimated)
            rightmatrix2[currrowidx,8] <- FC.WLSE(trainH,observed,estimated,pc)
            rightmatrix2[currrowidx,9] <- FC.EP(trainH,observed,estimated)
            rightmatrix2[currrowidx,10] <- FC.MEOP(1,trainH,observed,estimated)
            rightmatrix2[currrowidx,11] <- FC.TOUP(trainH,observed,estimated)
            rightmatrix2[currrowidx,12] <- PRatio(predH,trainH)
            rightmatrix2[currrowidx,13] <- FC.EP(predH,observed,estimated)
            rightmatrix2[currrowidx,14] <- FC.MEOP(trainH,predH,observed,estimated)
            rightmatrix2[currrowidx,15] <- abs(FC.Bias(1,trainH,observed,estimated))
            lastmatrix2[currrowidx,1] <- GROUP(FC.Bias(trainH,predH,observed,estimated))
          }else if(isFCTBF=="TBF"){
            rightmatrix2[currrowidx,1] <- TBF.MSE(trainH,numparam,observed,estimated)
            rightmatrix2[currrowidx,2] <- TBF.MAE(trainH,numparam,observed,estimated)
            rightmatrix2[currrowidx,3] <- TBF.Rsquare(trainH,observed,estimated)
            rightmatrix2[currrowidx,4] <- TBF.Noise(trainH,tbfestimated)
            rightmatrix2[currrowidx,5] <- TBF.Bias(1,trainH,observed,estimated)
            rightmatrix2[currrowidx,6] <- TBF.Variation(trainH,observed,estimated)
            rightmatrix2[currrowidx,7] <- TBF.PRR(trainH,observed,estimated)
            rightmatrix2[currrowidx,8] <- TBF.WLSE(trainH,observed,estimated,pc)
            rightmatrix2[currrowidx,9] <- TBF.PE(trainH,observed,estimated)
            rightmatrix2[currrowidx,10] <- TBF.AE(1,trainH,observed,estimated)
            rightmatrix2[currrowidx,11] <- TBF.TOUP(trainH,observed,estimated)
            rightmatrix2[currrowidx,12] <- PRatio(predH,trainH)
            rightmatrix2[currrowidx,13] <- TBF.PE(predH,observed,estimated)
            rightmatrix2[currrowidx,14] <- TBF.AE(trainH,predH,observed,estimated)
            rightmatrix2[currrowidx,15] <- abs(TBF.Bias(1,trainH,observed,estimated))
            lastmatrix2[currrowidx,1] <- GROUP(TBF.Bias(trainH,predH,observed,estimated))
          }
          tempidx <- tempidx + 1
          
        }
      }
    }
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- colnmleft
  
  rightdf <- data.frame(rightmatrix)
  lastdf <- data.frame(lastmatrix)
  rightdf <- cbind(rightdf,lastdf)
  
  names(rightdf) <- goffullvec
  
  retinstance <- NULL
  if(!is.null(leftmatrix2)){
    
    leftdf2 <- data.frame(leftmatrix2)
    names(leftdf2) <- colnmleft
  
    rightdf2 <- data.frame(rightmatrix2)
    lastdf2 <- data.frame(lastmatrix2)
    rightdf2 <- cbind(rightdf2,lastdf2)
    names(rightdf2) <- goffullvec
    
    retinstance <- na.omit(rbind(cbind(leftdf,rightdf),cbind(leftdf2,rightdf2)))
  }else{
    retinstance <- na.omit(cbind(leftdf,rightdf))
  }
  
  retinstance$Group <- as.factor(retinstance$Group)
  
  nmzretinstance <- NULL
  
  for(i in 1:length(pertlist)){
    currpert <- pertlist[i]
    subdf <- retinstance[retinstance$Pert==currpert,]
    
    for(j in 1:length(nrgofvec)){
      
      valvec  <- subdf[,nrgofvec[j]]
      
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
            if(nrgofvec[j]=="Rsquare"){
              nmzvec[k] <- 1
            }else{
              nmzvec[k] <- 0
            }
          }else{
            if(denominator==0){
              if(nrgofvec[j]=="Rsquare"){
                nmzvec[k] <- 0
              }else{
                nmzvec[k] <- 1
              }
            }else{
              if(nrgofvec[j]=="Rsquare"){
                nmzvec[k] <- 0 + ((maxval - valvec[k])*0.9)/denominator
              }else{
                nmzvec[k] <- 0.1 + ((maxval - valvec[k])*0.9)/denominator
              }
            }
          }
        }
        subdf[paste0("N",nrgofvec[j])] <- nmzvec
        
      }else if(length(valvec)==1){
        if(nrgofvec[j]=="Rsquare"){
          subdf[paste0("N",nrgofvec[j])] <- c(0)
        }else{
          subdf[paste0("N",nrgofvec[j])] <- c(1)
        }
      }
      
    }
    
    for(j in 1:length(crivec)){
      
      valvec  <- subdf[,crivec[j]]
      
      if(length(valvec)>1){
        
        infidx <- which(valvec==Inf)
        
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
        minval <- ifelse(minval<0.000001,0.000001,minval)
        
        ratvec <- vector()
        for (k in 1:length(valvec)) {
          if(is.infinite(valvec[k])){
            ratvec[k] <- Inf
          }else{
            if(minval==0){
              ratvec[k] <- valvec[k]
            }else{
              ratvec[k] <- (valvec[k] - minval) / minval
            }
          }
        }
        subdf[paste0("R",crivec[j])] <- ratvec
      }else if(length(valvec)==1){
        subdf[paste0("R",crivec[j])] <- c(0)
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

getSRGMGOFOrig <-function(df,pEst,
                     modelvec,
                     gofvec,
                     trainh,predh){
  
  pc <- 0.2
  
  leftmatrix <- NULL
  rightmatrix <- NULL
  lastmatrix <- NULL
  
  leftmatrix <- matrix(nrow=length(modelvec),ncol=1)
  rightmatrix <- matrix(nrow=length(modelvec),ncol=length(gofvec))
  lastmatrix <- matrix(nrow=length(modelvec),ncol=2)
  
  observed <- NULL
  observed <- df$n
  
  for(i in 1:length(modelvec)){
    
    currmodel <- modelvec[i] #current model name
      
    estimated <- NULL
    frestimated <- NULL
    numparam <- 0
      
    pEstsub <- pEst[[currmodel]]
    estimated <- pEstsub$EstElap
    frestimated <- pEstsub$EstFr

    numparam <- pEstsub$NumParam
    
    currrowidx  <- i
    
    leftmatrix[currrowidx,1] <- currmodel
    
    rightmatrix[currrowidx,1] <- FC.MSE(trainh,numparam,observed,estimated)
    rightmatrix[currrowidx,2] <- FC.MAE(trainh,numparam,observed,estimated)
    rightmatrix[currrowidx,3] <- FC.Rsquare(trainh,observed,estimated)
    rightmatrix[currrowidx,4] <- FC.Noise(trainh,frestimated)
    rightmatrix[currrowidx,5] <- FC.Bias(1,trainh,observed,estimated)
    rightmatrix[currrowidx,6] <- FC.Variation(trainh,observed,estimated)
    rightmatrix[currrowidx,7] <- FC.PRR(trainh,observed,estimated)
    rightmatrix[currrowidx,8] <- FC.WLSE(trainh,observed,estimated,pc)
    rightmatrix[currrowidx,9] <- FC.EP(trainh,observed,estimated)
    rightmatrix[currrowidx,10] <- FC.MEOP(1,trainh,observed,estimated)
    
    lastmatrix[currrowidx,1] <- FC.EP(predh,observed,estimated)
    lastmatrix[currrowidx,2] <- FC.MEOP(trainh,predh,observed,estimated)
  }
  
  leftdf <- data.frame(leftmatrix)
  names(leftdf) <- c('Model')
  
  rightdf <- data.frame(rightmatrix)
  names(rightdf) <- gofvec
  
  lastdf <- data.frame(lastmatrix)
  names(lastdf) <- c('EP','MEOP')
  
  retinstance <- cbind(leftdf,rightdf,lastdf)
  
  retinstance <- na.omit(retinstance)
  
  for(j in 1:length(gofvec)){
    
    valvec  <- retinstance[,gofvec[j]]
    
    if(length(valvec)>1){
      
      infidx <- which(is.infinite(valvec))
      
      noinfsub <- NULL
      if(length(infidx)>0){
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
          if(gofvec[j]=="Rsquare"){
            nmzvec[k] <- 1
          }else{
            nmzvec[k] <- 0
          }
        }else{
          if(denominator==0){
            if(gofvec[j]=="Rsquare"){
              nmzvec[k] <- 0
            }else{
              nmzvec[k] <- 1
            }
          }else{
            if(gofvec[j]=="Rsquare"){
              nmzvec[k] <- 0 + ((maxval - valvec[k])*0.9)/denominator
            }else{
              nmzvec[k] <- 0.1 + ((maxval - valvec[k])*0.9)/denominator
            }
          }
        }
      }
      retinstance[paste0("N",gofvec[j])] <- nmzvec
      
    }else if(length(valvec)==1){
      if(gofvec[j]=="Rsquare"){
        retinstance[paste0("N",gofvec[j])] <- c(0)
      }else{
        retinstance[paste0("N",gofvec[j])] <- c(1)
      }
    }
  }
  
  crivec <- c('EP','MEOP')
    
  for(j in 1:length(crivec)){
    
    valvec  <- retinstance[,crivec[j]]
    
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
      retinstance[paste0("N",crivec[j])] <- nmzvec
    }else if(length(valvec)==1){
      retinstance[paste0("N",crivec[j])] <- c(1)
    }
  }
  
  return(retinstance)
}
