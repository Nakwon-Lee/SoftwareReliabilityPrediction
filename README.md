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
			  | i | n |
			  | 1 | 0 |
			  | 2 | 1 |
			  | . | . |
			  | . | . |
			  | . | . |
			  | 20| 10|

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
which is learned from 36 existing failure data we have