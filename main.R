rm(list = ls())
ls()

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
source('getDDMGOF.R')
source('getMETAInfo.R')
source('reliabilityPrediction.R')
source('RecDDSRM.R')
source('MetaKnow.R')
source('getDataPoints.R')
source('settingOvsOMC.R')

# cli execution start
pargs <- commandArgs()
setwd(pargs[1])   # pargs[1]: path where the program (R files) installed
load(pargs[2])    # pargs[2]: path where the RData file for learning-based models installed
                  # the file includes "pmodels"
fdata <- read.csv(pargs[3],header=TRUE) # pargs[3]: the path of the csv file s
                                        # for target failure data for reliability
                                        # prediction
# cli execution end

# manual execution start
# setwd('D:/SR.SRGMR/')
# load('selectionmodels.RData')
# fdata <- read.csv('D:/SR.ReliabilityData/Data/KCC3.csv',header=TRUE)
# manual execution end

stopifnot(exists('pmodels'))
resultlist <- reliabilityPrediction(fdata,pmodels)