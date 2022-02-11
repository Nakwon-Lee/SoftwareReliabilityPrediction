Software reliability prediction tool

---command-line instruction for reliability prediction---

$ Rscript main.R PATH_PROGRAM PATH_MODELS PATH_TARGET
PATH_PROGRAM = path where the software reliability prediction tool is installed
PATH_MODELS = the RData file name storing the object "pmodels" that includes
			  selection models produced by learning-based techniques
PATH_TARGET = the csv file name for the target failure data which follows the
              format: i is the unit time; n is the number of cumulative failures
			  at each unit time
			  ex)
			  | i | n | f |
			  | 1 | 0 | 0 |
			  | 2 | 1 | 1 |
			  | 3 | 3 | 2 |
			  | . | . | . |
			  | . | . | . |
			  | 20| 10| 0 |

ex)
Rscript main.R ~/SoftwareReliabilityTool selectionmodels.RData exampletarget.csv
If the tool is installed in SoftwareReliabilityTool directoty of the home directoty
			  
---command-line instruction for generating a selection model for a learning-based
technique---

$ Rscript modelGen.R PATH_PROGRAM PATH_MODELS MODEL_NAME PATH_HISTORICALDATA
MODEL_NAME = the name of the model to learn {GOFC,META,AMETA}
             GOFC = the goodness-of-fit based model selection
			 META = the meta-data based model selection
			 AMETA = the feature learning based model selection

PATH_HISTORICALDATA = the RData file name storing the object "dflist" which 
                      contains the historical failure data sets
				
---pre-learned selection models---
We provide the pre-learned selection models in the file "selectionmodels.RData"
which is learned from 36 existing failure data we have.

---API list---
getSRGMGOF::getSRGMEstOrig(df,modelvec,trainh,predh)

get results of estimated and predicted failure counts and cumulative number of
failures for a given failure data set using single srgm models

-input-
df: the dataframe of a given failure data set
modelvec: a vector of single srgm model names for prediction
          currently supporting 13 models
		  {"GO","GG","Gompz","ISS","MD","MO","YID1",
		  "YID2","DSS","PNZ","PZ","PZI","Logi"}
trainh: unit time for model fitting
predh: the last unit time for the prediction

-value-
$trainH: the unit time for model fitting
a list of predicted results named by each model name of the modelvec
for each element of the list:
$(Model)$NumParam: number of parameters of the srgm model
$(Model)$EstElap: the estimated and predicted cumulative number of failures from
          the unit time 1 to the unit time predh
$(Model)$EstFr: the estimated and predicted failure counts at each unit time from
          the unit time 1 to the unit time predh

		  
getDDMGOF::getDDMEstOrig(df,model,trainh,predh)

get results of estimated and predicted failure counts and cumulative number of
failures for a given failure data set using a single ddm model

-input-
df: the dataframe of a given failure data set
model: a single ddm model name for prediction
          currently supporting 2 models
		  {"SVR","ANN"}
trainh: unit time for model learning
predh: the last unit time for the prediction

-value-
$trainH: the unit time for model fitting
a list of predicted results named by each model name of the modelvec
for each element of the list:
$(Model)$dim: the window size of the previous data points for model learning
$(Model)$EstElap: the estimated and predicted cumulative number of failures from
          the unit time 1 to the unit time predh
$(Model)$EstFr: the estimated and predicted failure counts at each unit time from
          the unit time 1 to the unit time predh


getMETAInfo::getMETAInfoOrig(df,ptrainh)

get meta-knowledge metric values of a given failure data set
currently supporting seven metrics
{"Variance","Inclination","AutoCorr","MetaNO","NumP","LapFact","SubAddi"}

-input-
df: the dataframe of a given failure data set
ptrainh: the last unit time for meta-knowledge metric calculation

-value-
a dataframe of one object with 8 variables {"Model","Variance","Inclination",
"AutoCorr","MetaNO","NumP","LapFact","SubAddi"} where Model has dummy value
Each variable indicates each metric value of the given failure data set.


getDataPoints::getDataPointsSingleLength(df,trainh,nump)

get normalized (between 0 to 1) cumulative number of failures
with fixed unit time

-input-
df: the dataframe of a given failure data set
trainh: the last unit time of the observed data
nump: the fixed unit time for the data generation

-value-
$raw: the normalized (between 0 to 1) cumulative number of failures
      from the unit time 1 to the unit time trainh
$norm: the normalized (between 0 to 1) cumulative number of failures
       of the fixed unit time (nump)
	   

