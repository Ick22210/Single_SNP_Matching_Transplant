#need to change this to be where the RealDataAnalysisCode.R file is located
fileLocation = paste0("/home/vlynn/Paper1_CodeToSendToCollabs/")
source(paste0(fileLocation,"RealDataAnalysisCode_Updated.R"))

#set working directory to where the data is located
#this will need to be changed for your data
setwd(fileLocation)

ex_datafile = "Example_Datafile.csv"
ex_pairIds = "Example_Pair_Ids.txt"
ex_covariates = "Example_Covariates.txt"
ex_outfile =  "ExampleOut"
runDonorRecipScoreMethod_par(datafile = ex_datafile, matchedPairs = ex_pairIds, 
                         covFile = ex_covariates, pvalueThreshold = 1, pvalueSig = 0.05,
                         numCores = 4, outFile = ex_outfile)