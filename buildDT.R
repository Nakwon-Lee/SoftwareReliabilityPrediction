setting1traincl <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  # ADTree <- make_Weka_classifier("weka/classifiers/trees/ADTree")
  # adtmodel <- ADTree(targetformula,data=ptrain)
  
  # adtmodel <- best.svm(targetformula,data=ptrain)
  
  adtmodel <- rpart(targetformula,data=ptrain)
  
  # adtmodel <- best.randomForest(targetformula,data=ptrain)
  
  # adtmodel <- lm(targetformula,data=ptrain)
  
  # REPTree <- make_Weka_classifier("weka/classifiers/trees/REPTree")
  # adtmodel <- REPTree(targetformula,data=ptrain)
  
  return(adtmodel)
}

setting1trglm <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  adtmodel <- glm(targetformula,data=ptrain)
  
  return(adtmodel)
}

temptrain <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  adtmodel <- lm(targetformula,data=ptrain)
  
  return(adtmodel)
}

setting1predictcl <- function(fittedtree,pvec,df){
  retbest <- NULL
  
  for(i in 1:nrow(df)){
    temp <- predict(fittedtree,df[i,],type="class")
    df[i,"Pred"] <- temp[1]
  }
  
  retbest <- df
  
  return(retbest)
}

setting1train <- function(ptrain,vars,target){
  
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  REPTree <- make_Weka_classifier("weka/classifiers/trees/REPTree")
  rpartmodel <- REPTree(targetformula,data=ptrain)
  
  return(rpartmodel)
}

setting1predict <- function(fittedtree,pvec,df){
  
  retbest <- NULL
  
  for(i in 1:nrow(df)){
    temp <- predict(fittedtree,df[i,pvec])
    df[i,"Pred"] <- temp[1]
  }
  
  retbest <- df
  
  return(retbest)
}

setting17train <- function(ptrain,vars,target){
  
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

  # rpartmodel <- glm(targetformula,data=ptrain,family = gaussian(link = "identity"))
  # rpartmodel <- gam(targetformula,family = gaussian(link = "identity"),data=ptrain)
  
  # rpartmodel <- rlm(targetformula,data=ptrain)
  
  # REPTree <- make_Weka_classifier("weka/classifiers/trees/REPTree")
  # rpartmodel <- REPTree(targetformula,data=ptrain)
  
  # rpartmodel <- pcr(targetformula,data=ptrain,validation = "CV")
  
  return(rpartmodel)
}

setting17predict <- function(fittedtree,pvec,df){
  
  retbest <- NULL
  
  for(i in 1:nrow(df)){
    temp <- predict(fittedtree,as.matrix(df[i,c(pvec)]))
    # temp <- predict(fittedtree,df[i,c(pvec)])
    df[i,"Pred"] <- temp[1]
  }
  
  retbest <- df
  
  return(retbest)
}

setting3train <- function(ptrain,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  rpartmodel <- J48(targetformula, data = ptrain)
  
  # rpartmodel <- rpart(targetformula, data = ptrain, method = "class",
  #                     control = rpart.control(minsplit = 10, cp=0.01))
  
  return(rpartmodel)
}

setting3predict <- function(fittedtree,vars,df){
  
  stopifnot(nrow(df)==1)
  
  predictedlist <- NULL
  predictedlist <- predict(fittedtree,df[,vars])
  
  return(predictedlist)
}

setting44train <- function(ptrain,vars,target){
  
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  tttrain <- ptrain[c(vars,target)]
  tttrain[target] <- ifelse(tttrain[target]=="1",1,0)

  rpartmodel <- neuralnet(targetformula,data=ptrain,
                          stepmax = 1e+7,hidden = 3,
                          threshold = 0.01)
  
  # rpartmodel <- rpart(targetformula, data = ptrain, method = "class",
  #                     control = rpart.control(minsplit = 10, cp=0.01))
  
  return(rpartmodel)
}

setting44predict <- function(fittedtree, vars, df){
  stopifnot(nrow(df)==1)
  
  predictedlist <- NULL
  predictedlist <- predict(fittedtree,df[1,vars])
  
  return(predictedlist)
}

predictNN <- function(object, newdata, rep = 1, all.units = FALSE) {
  
  weights <- object$weights[[rep]]
  num_hidden_layers <- length(weights) - 1
  
  # Init prediction with data, subset if necessary
  if (ncol(newdata) == length(object$model.list$variables)) {
    pred <- as.matrix(newdata)
  } else {
    pred <- as.matrix(newdata[, object$model.list$variables])
  }
  
  
  # Init units if requested
  if (all.units) {
    units <- list(pred)
  }
  
  # Hidden layers
  if (num_hidden_layers > 0) {
    for (i in 1:num_hidden_layers) {
      pred <- object$act.fct(cbind(1, pred) %*% weights[[i]])
      
      # Save unit outputs if requested
      if (all.units) {
        units <- append(units, list(pred))
      }
    }
  }
  
  # Output layer: Only apply activation function if non-linear output
  pred <- cbind(1, pred) %*% weights[[num_hidden_layers + 1]]
  if (!object$linear.output) {
    pred <- object$act.fct(pred)
  }
  
  # Save unit outputs if requested
  if (all.units) {
    units <- append(units, list(pred))
  } 
  
  # Return result
  if (all.units) {
    units
  } else {
    pred
  }
}

setting5train <- function(ptrain,Grd,vars,target){
  right <- vars[1]
  if(length(vars)>1){
    for (i in 2:length(vars)){
      right <- paste0(right,"+",vars[i])
    }
  }
  
  targetformula <- as.formula(paste0(target,"~",right))
  
  rpartmodel <- glm(targetformula,data=ptrain)
  
  # rpartmodel <- svm(targetformula,data = ptrain,
  #               type = "eps-regression",kernel = "radial",
  #               epsilon = Grd$eps, cost = Grd$cost,
  #               gamma = Grd$K)
  
  # rpartmodel <- randomForest(targetformula, data=ptrain)
  
  # rpartmodel <- rpart(targetformula,data=ptrain,method = "anova",
  #                     control=rpart.control(minsplit = 1,minbucket = 1,
  #                                           cp=0.001))
  
  # ptraindt <- ptrain[,c(vars)]
  # plabel <- ptrain[,c(target)]
  # 
  # rpartmodel <- xgboost(data = as.matrix(ptraindt),
  #                       label = as.matrix(plabel), verbose = 0,
  #                       max.depth = 20, eta = 1, nthread = 2,
  #                       nrounds = 1000, objective = "reg:squarederror")
  
  return(rpartmodel)
}

setting5predict <- function(fittedtree, vars, df){
  
  stopifnot(nrow(df)==1)
  
  subdf <- df[,vars]
  
  temp <- predict(fittedtree,subdf) # RandForest, SVM, GLM, DT
  # predictedlist <- predict(fittedtree,as.matrix(df[,c(vars)])) # XGB
  
  df[1,"Pred"] <- temp[1]
  
  # retbest <- round(predictedlist[[1]]) # GLM
  # retbest <- predictedlist # DT, RandForest, SVM, XGB
  
  return(df)
}
