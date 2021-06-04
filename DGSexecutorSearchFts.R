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

setwd('D:/SR.SRGMR')

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

load('evaluator.CV/NoSMOTE.SAv2.VIF.FULL/envQ1aug.evalFull.1000.RData')

gaparam <- list()
gaparam$popsize <- 20
gaparam$miter <- 100

searchedFeats <- list()

for (i in 1:length(vCriVec)) {
  searchedFeats[[vCriVec[i]]] <- searchforFeatures(plist = vFCDGoflist,pnlist = vFCDGofNlist,pmetalist = vFCDMetalist,
                                                        regvars = paste0('N',c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP")),
                                                        ppert = vtotpert,pcri = vCriVec[i],pddmres = setting7retFCD,pgofres = setting111retFCD,
                                                        pgofddmlist = vFCDGofDDMlist,
                                                        goffeats = c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP"),
                                                        metafeats = vNRMetaVec,
                                                   pparam = searchedEvaluator[[vCriVec[i]]]$param,
                                                   gaparam = gaparam,pparallel = 'par')
  save(searchedFeats,file = paste0('envQ1aug.evalFeats.GA.',gaparam$popsize,'.',gaparam$miter,'.WIN.RData'))
}
