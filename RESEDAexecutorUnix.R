rm(list = ls())
ls()

library(smotefamily)
library(performanceEstimation)
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
library(foreach)
library(doParallel)

setwd('~/git/SoftwareReliabilityPrediction')

source('setting1code.R')
source('buildDT.R')
source('settingMisccode.R')
source('setting16build.R')
source('setting3code.R')
source('smotethis.R')
source('Misc.R')
source('fitSRGMfromDDM.R')
source('GO.R')
source('GG.R')
source('Gompz.R')
source('ISS.R')
source('MD.R')
source('MO.R')
source('YID1.R')
source('YID2.R')
source('DSS.R')
source('PNZ.R')
source('PZ.R')
source('PZI.R')
source('Logi.R')
source('JM.R')
source('GEO.R')
source('MB.R')
source('LV.R')
source('EvaluatorGen.R')
source('settingRESEDAPert.R')
source('getSRGMGOF.R')
source('getDDMGOF.R')
source('parameterEstimation.R')
source('Criteria.R')
source('getMETAInfo.R')
source('MetaKnow.R')

load('envs/data.fcddflist.RData')
load('envs/data.fcdsrgmestlist.RData')
load('envs/data.fcdddmestlist.RData')

arguments <- commandArgs(trailingOnly = TRUE)

gaparam <- list()
gaparam$popsize <- as.integer(arguments[1])
gaparam$miter <- as.integer(arguments[2])

print(gaparam)

ppertvec <- c('P4','P5','P6','P7','P8','P9')

settingDGSretFCD <- list()
for (i in 1:length(ppertvec)) {
  settingDGSretFCD[[ppertvec[i]]] <- settingRESEDAPert(pdflist = vFCDlist,
                                                            pgoffts = c("MSE","MAE","Rsquare","Noise","Bias2","Variation","PRR","WLSE","CEP","CMEOP"),
                                                            pmetafts = c('Variance','Inclination','AutoCorr','MetaNO','NumP','LapFact','SubAddi'),
                                                            pcri = arguments[3],
                                                       pmodelvec = c("GO","GG","Gompz","ISS","MD","MO","YID1","YID2","DSS","PNZ","PZ","PZI","Logi"),
                                                            pgaparam = gaparam,cpert = ppertvec[i],
                                                       pddmestlist = tddmestlist[[ppertvec[i]]],
                                                       psrgmestlist = tsrgmestlist[[ppertvec[i]]],
                                                       pmodelforddm = 'SVR',
                                                       evalGen = evaluatorGenRESEDA,preDict = predictionRESEDA)
}

save(settingDGSretFCD,file = paste0('RESEDAres',arguments[3],'.RData'))
