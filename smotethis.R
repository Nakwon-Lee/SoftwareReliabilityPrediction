smotethis <- function(df,pformula,pnamevec,ptarget,pConf,pK){
  
  newtraindata <- NULL
  newtraindata <- df[,c(pnamevec,ptarget)]
  
  classtab <- array(table(newtraindata[,ptarget]))
  pknn <- min(classtab)
  if(pknn>5){
    pknn <- 5
  }
  
  countdf <- as.data.frame(table(newtraindata[,ptarget]))
  stopifnot(nrow(countdf)==2)
  majornum <- max(countdf$Freq)
  minornum <- countdf[countdf$Freq!=majornum,"Freq"]
  minorclass <- countdf[countdf$Freq==minornum,"Var1"]
  
  stopifnot(length(minornum)!=0)
  
  mmratio <- 600
  
  if(pknn > 0 && minornum > 1){
    
    # ttt <- NULL
    # ttt <- ubRacing(pformula,newtraindata, "randomForest", positive = 1, ubConf = pConf, ntree = 5)
    
    # rrr <- NULL
    # rrr <- ubBalance(X = newtraindata[,c(pnamevec)],Y = newtraindata[,c(ptarget)], type = ttt$best, percOver = 300, percUnder = 150)
    # newtraindata <- cbind(rrr$X,rrr$Y)
    # names(newtraindata) <- c(pnamevec,ptarget)
    
    # st SMOTE
    temptrain <- performanceEstimation::smote(pformula,newtraindata,perc.over = mmratio, k = pknn, perc.under = 100)
    # ed SMOTE

    # ttt <- NULL
    # ttt <- ADAS(newtraindata[,c(pnamevec)],newtraindata[,ptarget],K = pknn)
    # newtraindata <- ttt$data
    # names(newtraindata) <- c(pnamevec,ptarget)
    # newtraindata[,ptarget] <- as.factor(newtraindata[,ptarget])
    
    temptrain <- na.omit(temptrain)
    tempminor <- temptrain[temptrain[,ptarget]==minorclass,]
    tempminornew <- tempminor[(minornum+1):nrow(tempminor),]
    if(nrow(tempminornew) > (majornum-minornum)){
      tempminornew <- tempminornew[sample(nrow(tempminornew),(majornum-minornum)),]
    }
    
    newtraindata <- rbind(tempminor[1:minornum,],tempminornew)
  }
  
  return(newtraindata)
}