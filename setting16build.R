settingXGBtrain <- function(ptrain,vars,target){
  
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  ptraindt <- ptrain[,c(vars)]
  plabel <- ptrain[,c(target)]

  rpartmodel <- xgboost(data = as.matrix(ptraindt),
                        label = as.matrix(plabel), verbose = 0,
                        max.depth = 7, eta = 1, nthread = 2,
                        nrounds = 100, objective = "reg:logistic")
  
  # LogiReg <- make_Weka_classifier("weka/classifiers/functions/Logistic")
  # rpartmodel <- LogiReg(targetformula,data=ptrain)
  
  # rpartmodel <- glm(targetformula,data=ptrain,family = gaussian(link = "log"))
  # rpartmodel <- gam(targetformula,family = gaussian(link = "identity"),data=ptrain)
  
  # rpartmodel <- rlm(targetformula,data=ptrain)
  
  # REPTree <- make_Weka_classifier("weka/classifiers/trees/REPTree")
  # rpartmodel <- REPTree(targetformula,data=ptrain)
  
  # rpartmodel <- pcr(targetformula,data=ptrain,validation = "CV")
  
  return(rpartmodel)
}

settingXGBpredict <- function(fittedtree,pvec,df){
  retbest <- NULL
  for(i in 1:nrow(df)){
    temp <- predict(fittedtree,as.matrix(df[i,c(pvec)]))
    df[i,"Pred"] <- temp[1]
  }
  retbest <- df
  return(retbest)
}

settingBasicpredict <- function(fittedtree,pvec,df){
  retbest <- NULL
  for(i in 1:nrow(df)){
    temp <- predict(fittedtree,df[i,c(pvec)])
    df[i,"Pred"] <- temp[1]
  }
  retbest <- df
  return(retbest)
}

settingGLMtrain <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  targetformula <- as.formula(paste0(target,"~",right))
  rpartmodel <- glm(targetformula,data=ptrain)
  return(rpartmodel)
}

settingGAMtrain <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  targetformula <- as.formula(paste0(target,"~",right))
  rpartmodel <- gam(targetformula,family = gaussian(link = "inverse"),data=ptrain)
  return(rpartmodel)
}

settingRLMtrain <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  targetformula <- as.formula(paste0(target,"~",right))
  rpartmodel <- rlm(targetformula,data=ptrain)
  return(rpartmodel)
}

settingPCRtrain <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  targetformula <- as.formula(paste0(target,"~",right))
  rpartmodel <- pcr(targetformula,data=ptrain,validation = "CV")
  return(rpartmodel)
}

settingPLStrain <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  targetformula <- as.formula(paste0(target,"~",right))
  rpartmodel <- plsr(targetformula,data=ptrain,validation = "CV")
  return(rpartmodel)
}

settingPLS2train <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  targetformula <- as.formula(paste0(target,"~",right))
  rpartmodel <- plsr(targetformula,data=ptrain,method = "widekernelpls",validation = "LOO")
  return(rpartmodel)
}

settingRDFtrain <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  targetformula <- as.formula(paste0(target,"~",right))
  rpartmodel <- randomForest(targetformula,data=ptrain)
  return(rpartmodel)
}
