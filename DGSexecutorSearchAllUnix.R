rm(list = ls())
ls()

library(smotefamily)
# library(DMwR)
library(minpack.lm)
library(rpart)
library(nls2)
library(e1071)
library(mlr)
library(neuralnet)
library(randomForest)
library(xgboost)
library(unbalanced)
library(RWeka)
library(RWekajars)
library(MASS)
library(pls)
library(mgcv)
library(regclass)
library(FSelector)
library(GA)
library(parallel)

setwd('~/git/SoftwareReliabilityPrediction')

load("envs/envQ1.data.RData")
load("envs/envQ1.vars.RData")
load('envs/envQ1.gofres.RData')
load('envs/envQ1.ddmres.RData')

source('setting1code.R')
source('setting16build.R')
source('settingMisccode.R')
source('settingDGS.R')
source('EvaluatorGen.R')
source('buildDT.R')
source('Misc.R')

arguments <- commandArgs(trailingOnly = TRUE)

gaparam <- list()
gaparam$popsize <- as.integer(arguments[1])
gaparam$miter <- as.integer(arguments[2])

searchedAll <- list()
  
if(arguments[3]=='P'){
  print('Parallel')
  searchedAll[[arguments[4]]] <- searchforALLGA(plist = vFCDGoflist,pnlist = vFCDGofNlist,pmetalist = vFCDMetalist,
                                                   regvars = paste0('N',c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP")),
                                                   ppert = vtotpert,pcri = arguments[4],pddmres = setting7retFCD,pgofres = setting111retFCD,
                                                   pgofddmlist = vFCDGofDDMlist,
                                                   goffeats = c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP"),
                                                   metafeats = vNRMetaVec,
                                                   gaparam = gaparam,pparallel = 'par')
}else{
  print('Sequential')
  searchedAll[[arguments[4]]] <- searchforALLGA(plist = vFCDGoflist,pnlist = vFCDGofNlist,pmetalist = vFCDMetalist,
                                                   regvars = paste0('N',c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP")),
                                                   ppert = vtotpert,pcri = arguments[4],pddmres = setting7retFCD,pgofres = setting111retFCD,
                                                   pgofddmlist = vFCDGofDDMlist,
                                                   goffeats = c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP"),
                                                   metafeats = vNRMetaVec,
                                                   gaparam = gaparam)
}
save(searchedAll,file = paste0('envQ1aug.evalAll.GA.',arguments[1],'.',arguments[2],'.',arguments[4],'.RData'))
