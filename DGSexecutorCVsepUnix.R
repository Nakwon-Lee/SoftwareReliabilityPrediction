rm(list = ls())
ls()

library(smotefamily)
library(DMwR)
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
library(car)
library(FSelector)

setwd('~/git/SR.SRGMR')

load("envs/envQ1.data.RData")
load("envs/envQ1.vars.RData")

source('setting1code.R')
source('setting16build.R')
source('settingMisccode.R')
source('settingDGS.R')
source('EvaluatorGen.R')
source('buildDT.R')
source('Misc.R')

Msted <- c(1:4)
  
thisprefix <- "evaluator.CV/NoSMOTE.100.CFS/envQ1aug"

for (Mkey in Msted) {
  settingDGSEvalGensep(pkey = Mkey,
                       goflist = vFCDGoflist,gofnlist = vFCDGofNlist,gofddmlist = vFCDGofDDMlist,
                       metadf = vFCDMetalist,pcri = vCriVec[1],totpert = vtotpert,
                       goffts = c("MSE","MAE","Rsquare","Noise","Variation","PRR","WLSE","CEP","CMEOP","Bias2"),
                       metafts = vNRMetaVec,
                       givenprefix = thisprefix)
}

for (Mkey in Msted) {
  settingDGSEvalGensep(pkey = Mkey,
                       goflist = vFCDGoflist,gofnlist = vFCDGofNlist,gofddmlist = vFCDGofDDMlist,
                       metadf = vFCDMetalist,pcri = vCriVec[2],totpert = vtotpert,
                       goffts = c("MSE","MAE","Rsquare","Noise","Variation","PRR","WLSE","CEP","CMEOP","Bias2"),
                       metafts = vNRMetaVec,
                       givenprefix = thisprefix)
}