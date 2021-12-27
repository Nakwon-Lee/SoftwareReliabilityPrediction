rm(list = ls())
ls()

library(nls2)
library(minpack.lm)
library(rpart)
library(e1071)
library(performanceEstimation)
library(gtools)

# library(smotefamily)
# library(mlr)
# library(neuralnet)
# library(randomForest)
# library(xgboost)
# library(unbalanced)
# library(RWeka)
# library(RWekajars)
# library(MASS)
# library(pls)
# library(mgcv)
# library(regclass)
# library(FSelector)
# library(GA)
# library(parallel)
# library(foreach)
# library(doParallel)

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
source('Misc.R')
source('parameterEstimation.R')
source('Criteria.R')
source('getSRGMGOF.R')
source('getMETAInfo.R')
source('setting1code.R')
source('setting3code.R')
source('settingOvsOMC.R')
source('buildDT.R')
source('MetaKnow.R')
source('getDDMGOF.R')
source('RecDDSRM.R')
source('smotethis.R')
source('getDataPoints.R')


# cli execution start

pargs <- commandArgs()
setwd(pargs[1])         # pargs[1]: path where the program (R files) installed
load(pargs[2])          # pargs[2]: path where the RData file for learning-based models installed
                        # the file includes "pmodels"
pmodelmode <- pargs[3]  # pargs[3]: the name of the model to learn 
                        # {GOFC,META,HYC,HYR,AMETA}
load(pargs[4])          # pargs[4]: the name of RData 
                        # file containing the list "dflist" that have dataframes 
                        # of each historical failure data
# cli execution end

# manual execution start
# setwd('D:/SR.SRGMR/')
# load('selectionmodels.RData')
# pmodelmode <- 'AMETA'
# load('historicalfailuredata.RData')
# manual execution end

stopifnot(exists('dflist'))

if(!exists('pmodels')){
  pmodels <- list()
}

if(pmodelmode=="GOFC"){
  pmodels$GOFC <- modelgenGOFC(dflist)
}else if(pmodelmode=="META"){
  pmodels$META <- modelgenMETA(dflist)
}else if(pmodelmode=="HYC"){
  pmodels$HYC <- modelgenHYC(dflist)
}else if(pmodelmode=="HYR"){
  pmodels$HYR <- modelgenHYR(dflist)
}else if(pmodelmode=="AMETA"){
  pmodels$AMETA <- modelgenAMETA(dflist)
}else{
  stop("Invalid model name for learning")
}

save(pmodels,file = 'selectionmodels.RData')
