rm(list = ls())
ls()

library(smotefamily)
library(performanceEstimation)
library(minpack.lm)
library(rpart)
library(nls2)
library(e1071)
library(unbalanced)
library(MASS)
library(pls)
library(mgcv)
library(regclass)
library(GA)
library(parallel)
library(foreach)
library(doParallel)
library(gtools)
library(tensorflow)
library(keras)
library(reticulate)

setwd('~/git/SoftwareReliabilityPrediction')

load('envs/data.fcdsrgmestlist.RData')
load('envs/data.fcdddmestlist.RData')
load('vdatalist.RData')
load('vgofddmlist.RData')
load('vgoflist.RData')

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
source('getDataPoints.R')
source('MetaKnow.R')
source('settingOvsOMC.R')

ppertvec <- c('P7','P8','P9')

settingOvOretFCD <- list()
for (i in 1:length(ppertvec)) {
  settingOvOretFCD[[ppertvec[i]]] <- settingOvOMCPert(vFCDlist,vFCDGoflist,vFCDGofDDMlist,
                                                      ppertvec[i],'EP',
                                                      c("GO","GG","Gompz","ISS","MD","MO","YID1","YID2","DSS","PNZ","PZ","PZI","Logi","SVR"),
                                                      tddmestlist,tsrgmestlist)
}

for (i in 1:length(ppertvec)) {
  settingOvOretFCD[[ppertvec[i]]] <- cbind(settingOvOretFCD[[ppertvec[i]]],settingOvOMCPert(vFCDlist,vFCDGoflist,vFCDGofDDMlist,
                                                                                            ppertvec[i],'MEOP',
                                                                                            c("GO","GG","Gompz","ISS","MD","MO","YID1","YID2","DSS","PNZ","PZ","PZI","Logi","SVR"),
                                                                                            tddmestlist,tsrgmestlist))
}