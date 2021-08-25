The main code file is RealDataAnalysisCode.R. 
You shouldn't need to open this file, but it contains the main program for running the code, so it will need to be sourced.  
Code requires R packages survival and epicalc.
I ran my example using R version 3.6.2.

The main function is called "runDonorRecipScoreMethod()". It has four needed inputs. 

1) Main datafile. This can be in any file format as long as  the data is formatted like a matrix. 
Number of rows should equal 2*number of D/R pairs (ie all Donors and Recips from the pairs), 
number of columns should correspond to SNPs. Columns should have values between 0-2 for additive model. 
Rownames should be the Ids of the Donor or Recipient the genotype information corresponds to in order to do D/R matching.

2) Id File. A text file with two columns. First column is Donor ID, second is matched Recipient ID. 
This file must be in .txt format or the code won't read it in correctly.

3) Covariate File. A text file with at least 3 columns. First column is Recipient Ids, 
second column is a measure of time to event (integer value), third column is an event indicator of event (1) vs censored (0). 
Any other columns will correspond to covariates that you wish to include in association analyses.

4) Name of output file. What you want to call the outputted text file. 
Output file gives Bonferroni corrected p-value threshold, and two tables (similar to Table 3 from the manuscript). 
First table is significant results of joint analyses (p <= 0.05). 
Columns are "SNP rsid", "Score used in analysis", "Joint HR for Score with 95% CI", "Joint HR for Recipient Genotype with 95% CI", 
"Joint P value". 

Second table gives all results of univariates analyses done after significant joint testing. 
Columns are "SNP rsid", "Score or Genotype used in analysis", "HR and 95% CI", "P value".

I have included sample files of all the inputs (Example_Datafile.csv, Example_Pair_Ids.txt, Example_Covariates.txt)
and an example of the output you get after running the code (ExampleOut.txt).

In addition, the RMarkdown file (RealDataAnalysis_Example.Rmd) has the example R code for how I obtained the output.
I will paste that here for simplicity, or you can look at the Rmd or the RealDataAnalysis_Example.pdf file for the code.

#source the RealDataAnalysisCode.R file 
#need to change this to be where the RealDataAnalysisCode.R file is located
fileLocation = paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Project_One_Work/CodeToSendToCollabs/")
source(paste0(fileLocation, "RealDataAnalysisCode.R"))
ex_datafile = "Example_Datafile.csv"
ex_pairIds = "Example_Pair_Ids.txt"
ex_covariates = "Example_Covariates.txt"
ex_outfile = "ExampleOut.txt"
runDonorRecipScoreMethod(datafile = ex_datafile,
matchedPairs = ex_pairIds,
covFile = ex_covariates,
outFile = ex_outfile)

