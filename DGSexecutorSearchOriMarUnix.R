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

searchedEvaluator <- searchforEvaluator(plist = vFCDGoflist,pnlist = vFCDGofNlist,pmetalist = vFCDMetalist,
                                        regvars = paste0('N',c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP")),
                                        ppert = vtotpert,pcri = vCriVec[1],pddmres = setting7retFCD,pgofres = setting111retFCD,
                                        pgofddmlist = vFCDGofDDMlist,
                                        goffeats = c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP"),
                                        metafeats = vNRMetaVec)
save(searchedEvaluator,file = paste0('envQ1aug.evalFull.',vCriVec[1],'.RData'))

searchedEvaluator <- searchforEvaluator(plist = vFCDGoflist,pnlist = vFCDGofNlist,pmetalist = vFCDMetalist,
                                        regvars = paste0('N',c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP")),
                                        ppert = vtotpert,pcri = vCriVec[2],pddmres = setting7retFCD,pgofres = setting111retFCD,
                                        pgofddmlist = vFCDGofDDMlist,
                                        goffeats = c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP"),
                                        metafeats = vNRMetaVec)
save(searchedEvaluator,file = paste0('envQ1aug.evalFull.',vCriVec[2],'.RData'))

