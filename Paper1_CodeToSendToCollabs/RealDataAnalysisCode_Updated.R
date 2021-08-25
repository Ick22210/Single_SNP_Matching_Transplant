###########################################
## Code to run analyses on real data
##
## Written by: Victoria Arthur
## Last edited: 2/2/2021
###########################################
#Input data format:
# Data should be in matrix-type form, any file format is fine (text file, csv, etc) 
# with number of rows = 2*number of D/R pairs (ie all Donors and Recips from the pairs)
# and num of columns should correspond to SNPs, 
# columns should have values between 0-2 for additive model
# Should also have rownames that are the Ids of the Donor and Recips for matching
##
# Also need a file of Ids, first column is Recipient ID, second column is matched Donor ID
# This file should be in ".txt" format
##
# Additional covariates should be in ".txt" format as well
# first column: Recipient ID (rownames)
# second column: Outcome time_to_event
# third column: Outcome event_indicator 
# remaining columns: all other covariates

#function to run the analyses
runDonorRecipScoreMethod = function(datafile = "genotypeData.csv", matchedPairs = "pairedIDs.txt", covFile = "covariates.txt", pvalueThreshold = 0.05, pvalueSig = 0.05, outFile = "out.txt"){
  #Inputs:
  # datafile is the name of the file containing the genotype data
  # matchedPairs is a text file containing the data on which donor id is matched to which recipient id
  # covFile is a text file containing additional covariate information
  # p-value threshold is the highest p-value that will be output for analyses, can change depending on what you want to output
  # p-value Sig is the significance level you want to reach in the joint analyses so that the SNP gets run through the univariate analyses
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  suppressMessages(require(survival))
  suppressMessages(require(lmtest))
  
  #read in the data
  #determine what file type
  fileName = read.table(text = datafile, sep = ".", as.is = TRUE)$V2
  if(fileName == "csv"){ #if it's a csv file, read in with read.csv
    geneticData = read.csv(datafile, header = TRUE, row.names = 1)
  } else { #otherwise should be good with a read.table
    geneticData = read.table(datafile, header = TRUE, row.names = 1)
  }
  
  ###############################################
  #Determine D/R pairs and match
  ###############################################  
  #read in the file with the donor/recipient pair information
  pairings = read.table(matchedPairs)
  #make sure column names are specified
  colnames(pairings) = c("R_Id", "D_Id")
  #assign each pair a number
  pairings$PairNumber = seq(1:length(pairings$R_Id))
  #split the original data into donor and recipient subsets
  donorData = geneticData[rownames(geneticData) %in% pairings$D_Id,]
  recipData = geneticData[rownames(geneticData) %in% pairings$R_Id,]
  #add an id column based on row names for merging
  donorData$D_Id = rownames(donorData)
  recipData$R_Id = rownames(recipData)
  #now merge the genotype data with the ids so that the genotype data has the pair numbers
  donorGenotypes_numbered = merge(donorData, pairings, by.x = "D_Id", by.y = "D_Id")
  recipGenotypes_numbered = merge(recipData, pairings, by.x = "R_Id", by.y = "R_Id")
  #make everything a data frame 
  donorGenotypes_numbered.df = as.data.frame(donorGenotypes_numbered)
  recipGenotypes_numbered.df = as.data.frame(recipGenotypes_numbered)
  #subset to only include Pairs with both D and R info
  donorGenotypes_numbered.df.subset = donorGenotypes_numbered.df[(donorGenotypes_numbered.df$PairNumber %in% recipGenotypes_numbered.df$PairNumber),]
  recipGenotypes_numbered.df.subset = recipGenotypes_numbered.df[(recipGenotypes_numbered.df$PairNumber %in% donorGenotypes_numbered.df$PairNumber),]
  #order by Pair number
  donorGenotypes_numbered.df.ordered = donorGenotypes_numbered.df.subset[order(donorGenotypes_numbered.df.subset$PairNumber),]
  recipGenotypes_numbered.df.ordered = recipGenotypes_numbered.df.subset[order(recipGenotypes_numbered.df.subset$PairNumber),]
  #make rownames into R_Id or D_Id
  rownames(donorGenotypes_numbered.df.ordered) = donorGenotypes_numbered.df.ordered$D_Id
  rownames(recipGenotypes_numbered.df.ordered) = recipGenotypes_numbered.df.ordered$R_Id
  #remove non-SNP columns
  donorGenotypes_numbered.df.ordered.snps = donorGenotypes_numbered.df.ordered[,2:(ncol(donorGenotypes_numbered.df.ordered)-2)]
  recipGenotypes_numbered.df.ordered.snps = recipGenotypes_numbered.df.ordered[,2:(ncol(recipGenotypes_numbered.df.ordered)-2)]
  
  ###############################################
  #Calculate all the Scores
  ###############################################  
  #calculate the difference between the two subjects
  diffScore = abs(donorGenotypes_numbered.df.ordered.snps - recipGenotypes_numbered.df.ordered.snps)
  
  #######################
  ## IBS Score
  ####################### 
  #if diff = 0, score = 0
  #if diff = 1, score = 1
  #if diff = 2, score = 2
  IBSMismatch = diffScore
  IBSMismatch = as.matrix(IBSMismatch)
  #add rownames to be R Ids
  rownames(IBSMismatch) = rownames(recipGenotypes_numbered.df.ordered.snps)
  #######################
  ## Incomp Score
  #######################
  #if diff = 0, score = 0
  #otherwise score = 1
  incomp = diffScore
  incomp = as.matrix(incomp)
  incomp[diffScore == 0] = 0
  incomp[diffScore != 0] = 1
  #add rownames to be R Ids
  rownames(incomp) = rownames(recipGenotypes_numbered.df.ordered.snps)
  #######################
  ## AMS
  #######################
  #mismatch if D has allele not in R
  #sum across both alleles in genotype
  #Score is either 0, 1, or 2
  alloMismatch = matrix(0, nrow = nrow(diffScore), ncol = ncol(diffScore)) #make default value 0
  alloMismatch[(donorGenotypes_numbered.df.ordered.snps == 0) & (recipGenotypes_numbered.df.ordered.snps == 2)] = 2 #Donor AA, Recip aa
  alloMismatch[(donorGenotypes_numbered.df.ordered.snps == 2) & (recipGenotypes_numbered.df.ordered.snps == 0)] = 2 #Donor aa, Recip AA
  alloMismatch[(donorGenotypes_numbered.df.ordered.snps == 1) & (recipGenotypes_numbered.df.ordered.snps == 2)] = 1 #Donor Aa, Recip aa
  alloMismatch[(donorGenotypes_numbered.df.ordered.snps == 1) & (recipGenotypes_numbered.df.ordered.snps == 0)] = 1 #Donor Aa, Recip AA
  alloMismatch[is.na(donorGenotypes_numbered.df.ordered.snps) | is.na(recipGenotypes_numbered.df.ordered.snps)] = NA #make sure NAs are preserved
  #match row and column names from the original data set
  #row names should be R Ids
  rownames(alloMismatch) = rownames(recipGenotypes_numbered.df.ordered.snps)
  colnames(alloMismatch) = colnames(donorGenotypes_numbered.df.ordered.snps)
  #######################
  ## Binary  Mismatch
  #######################
  #mismatch if D has allele not in R
  #Score is either 0 or 1
  binMismatch = matrix(0, nrow = nrow(diffScore), ncol = ncol(diffScore)) #make default value 0
  binMismatch[(donorGenotypes_numbered.df.ordered.snps == 1) & (recipGenotypes_numbered.df.ordered.snps == 2)] = 1 #Donor Aa, Recip aa
  binMismatch[(donorGenotypes_numbered.df.ordered.snps == 0) & (recipGenotypes_numbered.df.ordered.snps == 2)] = 1 #Donor AA, Recip aa
  binMismatch[(donorGenotypes_numbered.df.ordered.snps == 2) & (recipGenotypes_numbered.df.ordered.snps == 0)] = 1 #Donor aa, Recip AA
  binMismatch[(donorGenotypes_numbered.df.ordered.snps == 1) & (recipGenotypes_numbered.df.ordered.snps == 0)] = 1 #Donor Aa, Recip AA
  binMismatch[is.na(donorGenotypes_numbered.df.ordered.snps) | is.na(recipGenotypes_numbered.df.ordered.snps)] = NA #make sure NAs are preserved
  #match row and column names from the original data set
  #row names should be R Ids
  rownames(binMismatch) = rownames(recipGenotypes_numbered.df.ordered.snps)
  colnames(binMismatch) = colnames(donorGenotypes_numbered.df.ordered.snps)

  ###############################################
  # Run joint analyses
  ###############################################
  #need to combine the score data with the recipient data
  #add _recip to recipient genos to make col names different
  colnames(recipGenotypes_numbered.df.ordered.snps) = paste0(colnames(recipGenotypes_numbered.df.ordered.snps),"_recip")
  IBSMismatch_RGenos = merge(IBSMismatch, recipGenotypes_numbered.df.ordered.snps, by = "row.names")
  incomp_RGenos = merge(incomp, recipGenotypes_numbered.df.ordered.snps, by = "row.names")
  alloMismatch_RGenos = merge(alloMismatch, recipGenotypes_numbered.df.ordered.snps, by = "row.names")
  binMismatch_RGenos = merge(binMismatch, recipGenotypes_numbered.df.ordered.snps, by = "row.names")
  #fix the row names
  rownames(IBSMismatch_RGenos) = IBSMismatch_RGenos$Row.names
  rownames(incomp_RGenos) = incomp_RGenos$Row.names
  rownames(alloMismatch_RGenos) = alloMismatch_RGenos$Row.names
  rownames(binMismatch_RGenos) = binMismatch_RGenos$Row.names
  #remove the Row.names column
  IBSMismatch_RGenos = IBSMismatch_RGenos[,-1]
  incomp_RGenos = incomp_RGenos[,-1]
  alloMismatch_RGenos = alloMismatch_RGenos[,-1]
  binMismatch_RGenos = binMismatch_RGenos[,-1]
  #read in the covariate file
  covariates = read.table(covFile, header = TRUE)
  #set row names to be Ids from 1st row
  rownames(covariates) = covariates[,1]
  #remove the ids column
  covariates = covariates[,-1]
  #pull the col names for time to event and event indicator
  timeToEvent = colnames(covariates)[[1]]
  eventIndicator = colnames(covariates)[[2]]
  survFormulaInfo = paste0("Surv(", timeToEvent,", ", eventIndicator, ")")
  #determine number of covariates
  numCovs = ncol(covariates) - 2
  #pull additional covariate names
  covNames = colnames(covariates)[3:ncol(covariates)]
  #covariates to add to survival model
  allCovariates = covNames[[1]]
  for(jj in 2:length(covNames)){
    allCovariates = paste0(allCovariates, " + ", covNames[[jj]])
  }
  #merge the covariates with the new scores (By R_Id)
  covariates_wIBSMismatch_RGenos = merge(covariates, IBSMismatch_RGenos, by.x = 'row.names', by.y = 'row.names')
  covariates_wIncomp_RGenos = merge(covariates, incomp_RGenos, by.x = 'row.names', by.y = 'row.names')
  covariates_wAMS_RGenos = merge(covariates, alloMismatch_RGenos, by.x = 'row.names', by.y = 'row.names')
  covariates_wBinMM_RGenos = merge(covariates, binMismatch_RGenos, by.x = 'row.names', by.y = 'row.names')
  covariates_wRGenos = merge(covariates, recipGenotypes_numbered.df.ordered.snps, by.x = 'row.names', by.y = 'row.names')
  #pull list of SNPs and Recip SNPs
  SNPs = colnames(IBSMismatch)
  RSNPs = colnames(recipGenotypes_numbered.df.ordered.snps)
  #store significant SNPs
  sigSNPsIBS = c()
  sigSNPsIncomp = c()
  sigSNPsAMS = c()
  sigSNPsBinMM = c()
  IBSTable = matrix(0,nrow = 0, ncol = 6)
  IncompTable = matrix(0,nrow = 0, ncol = 5)
  AMSTable = matrix(0,nrow = 0, ncol = 6)
  BinMMTable = matrix(0,nrow = 0, ncol = 5)
  colnames(IBSTable) = colnames(AMSTable) =  c("SNP", "Jt Score HR (95% CI)", "Jt R Geno HR (95% CI)", "Jt P-value", "Freq of Score = 1", "Freq of Score = 2")
  colnames(IncompTable) = colnames(BinMMTable) = c("SNP", "Jt Score HR (95% CI)", "Jt R Geno HR (95% CI)", "Jt P-value", "Freq of Score = 1")
  univarTable = matrix(0,nrow = 0, ncol = 4)
  colnames(univarTable) = c("SNP", "Score/Genotype", "HR (95% CI)", "P-value")
  
  #calc Bonferroni p-value
  Bonferroni = 0.05/length(SNPs)
  
  for(ii in 1:length(SNPs)){
    ###########################################
    ##### IBS Mismatch Score
    ###########################################
    IBSSNPs = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ",SNPs[ii]," + ",RSNPs[ii])), data = na.omit(covariates_wIBSMismatch_RGenos))
    IBSSNPsNull = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates)), data = na.omit(covariates_wIBSMismatch_RGenos)) 
    #perform joint test
    jttestIBS = lrtest(IBSSNPs, IBSSNPsNull)
    
    if(jttestIBS$`Pr(>Chisq)`[2] <= pvalueThreshold){ #if joint p-value passes the threshold
      #pull joint p-value
      jtp = round(jttestIBS$`Pr(>Chisq)`[2], digits = 6)
      #pull HRs and 95% CIs
      summaryStats = summary(IBSSNPs)
      IBSSNPInfo = round(summaryStats$conf.int[SNPs[[ii]],], digits = 2)
      IBSHR = paste0(IBSSNPInfo[1]," (",IBSSNPInfo[3],", ",IBSSNPInfo[4],")")
      RSNPInfo = round(summaryStats$conf.int[RSNPs[[ii]],], digits = 2)
      RSNPHR = paste0(RSNPInfo[1]," (",RSNPInfo[3],", ",RSNPInfo[4],")")
      #Calc freq of score values
      freq1 = sum(na.omit(covariates_wIBSMismatch_RGenos[,colnames(covariates_wIBSMismatch_RGenos) == SNPs[ii]])==1)/length(na.omit(covariates_wIBSMismatch_RGenos[,colnames(covariates_wIBSMismatch_RGenos) == SNPs[ii]]))
      freq2 = sum(na.omit(covariates_wIBSMismatch_RGenos[,colnames(covariates_wIBSMismatch_RGenos) == SNPs[ii]])==2)/length(na.omit(covariates_wIBSMismatch_RGenos[,colnames(covariates_wIBSMismatch_RGenos) == SNPs[ii]]))
      #add to joint table
      tableRow = cbind(SNPs[[ii]], IBSHR, RSNPHR, jtp, round(freq1, digits = 2), round(freq2, digits = 2))
      IBSTable = rbind(IBSTable, tableRow)
    }
    #add SNP to sigSNPs if it passes the significance value
    if(jttestIBS$`Pr(>Chisq)`[2] <= pvalueSig){
      sigSNPsIBS = c(sigSNPsIBS, ii) #add SNP number to sigSNPs
    }
    ###########################################
    ##### Incompatibility Score
    ###########################################
    IncompSNPs = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ",SNPs[ii]," + ",RSNPs[ii])), data = na.omit(covariates_wIncomp_RGenos))
    IncompSNPsNull = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates)), data = na.omit(covariates_wIncomp_RGenos)) 
    #perform joint test
    jttestIncomp = lrtest(IncompSNPs, IncompSNPsNull)

    if(jttestIncomp$`Pr(>Chisq)`[2] <= pvalueThreshold){ #if joint p-value is significant
      #pull joint p-value
      jtp = round(jttestIncomp$`Pr(>Chisq)`[2], digits = 6)
      #pull HRs and 95% CIs
      summaryStats = summary(IncompSNPs)
      IncompSNPInfo = round(summaryStats$conf.int[SNPs[[ii]],], digits = 2)
      IncompHR = paste0(IncompSNPInfo[1]," (",IncompSNPInfo[3],", ",IncompSNPInfo[4],")")
      RSNPInfo = round(summaryStats$conf.int[RSNPs[[ii]],], digits = 2)
      RSNPHR = paste0(RSNPInfo[1]," (",RSNPInfo[3],", ",RSNPInfo[4],")")
      #Calc freq of score values
      freq1 = sum(na.omit(covariates_wIncomp_RGenos[,colnames(covariates_wIncomp_RGenos) == SNPs[ii]])==1)/length(na.omit(covariates_wIncomp_RGenos[,colnames(covariates_wIncomp_RGenos) == SNPs[ii]]))
      #add to joint table
      tableRow = cbind(SNPs[[ii]], IncompHR, RSNPHR, jtp, round(freq1,digits=2))
      IncompTable = rbind(IncompTable, tableRow)
    }
    #add SNP to sigSNPs if it passes the significance value
    if(jttestIncomp$`Pr(>Chisq)`[2] <= pvalueSig){
      sigSNPsIncomp = c(sigSNPsIncomp, ii) #add SNP number to sigSNPs
    }
    ###########################################
    ##### Allogenomics Mismatch Score
    ###########################################
    AMSSNPs = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ",SNPs[ii]," + ",RSNPs[ii])), data = na.omit(covariates_wAMS_RGenos))
    AMSSNPsNull = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates)), data = na.omit(covariates_wAMS_RGenos)) 
    #perform joint test
    jttestAMS = lrtest(AMSSNPs, AMSSNPsNull)

    if(jttestAMS$`Pr(>Chisq)`[2] <= pvalueThreshold){ #if joint p-value is significant
      #pull joint p-value
      jtp = round(jttestAMS$`Pr(>Chisq)`[2], digits = 6)
      #pull HRs and 95% CIs
      summaryStats = summary(AMSSNPs)
      AMSSNPInfo = round(summaryStats$conf.int[SNPs[[ii]],], digits = 2)
      AMSHR = paste0(AMSSNPInfo[1]," (",AMSSNPInfo[3],", ",AMSSNPInfo[4],")")
      RSNPInfo = round(summaryStats$conf.int[RSNPs[[ii]],], digits = 2)
      RSNPHR = paste0(RSNPInfo[1]," (",RSNPInfo[3],", ",RSNPInfo[4],")")
      #Calc freq of score values
      freq1 = sum(na.omit(covariates_wAMS_RGenos[,colnames(covariates_wAMS_RGenos) == SNPs[ii]])==1)/length(na.omit(covariates_wAMS_RGenos[,colnames(covariates_wAMS_RGenos) == SNPs[ii]]))
      freq2 = sum(na.omit(covariates_wAMS_RGenos[,colnames(covariates_wAMS_RGenos) == SNPs[ii]])==2)/length(na.omit(covariates_wAMS_RGenos[,colnames(covariates_wAMS_RGenos) == SNPs[ii]]))
      #add to joint table
      tableRow = cbind(SNPs[[ii]], AMSHR, RSNPHR, jtp, round(freq1,digits=2), round(freq2,digits = 2))
      AMSTable = rbind(AMSTable, tableRow)
    }
    #add SNP to sigSNPs if it passes the significance value
    if(jttestAMS$`Pr(>Chisq)`[2] <= pvalueSig){
      sigSNPsAMS = c(sigSNPsAMS, ii) #add SNP number to sigSNPs
    }
    ###########################################
    ##### Binary Mismatch Score
    ###########################################
    BinMMSNPs = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ",SNPs[ii]," + ",RSNPs[ii])), data = na.omit(covariates_wBinMM_RGenos))
    BinMMSNPsNull = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates)), data = na.omit(covariates_wBinMM_RGenos)) 
    #perform joint test
    jttestBinMM = lrtest(BinMMSNPs, BinMMSNPsNull)

    if(jttestBinMM$`Pr(>Chisq)`[2] <= pvalueThreshold){ #if joint p-value is significant
      #pull joint p-value
      jtp = round(jttestBinMM$`Pr(>Chisq)`[2], digits = 6)
      #pull HRs and 95% CIs
      summaryStats = summary(BinMMSNPs)
      BinMMSNPInfo = round(summaryStats$conf.int[SNPs[[ii]],], digits = 2)
      BinMMHR = paste0(BinMMSNPInfo[1]," (",BinMMSNPInfo[3],", ",BinMMSNPInfo[4],")")
      RSNPInfo = round(summaryStats$conf.int[RSNPs[[ii]],], digits = 2)
      RSNPHR = paste0(RSNPInfo[1]," (",RSNPInfo[3],", ",RSNPInfo[4],")")
      #Calc freq of score values
      freq1 = sum(na.omit(covariates_wBinMM_RGenos[,colnames(covariates_wBinMM_RGenos) == SNPs[ii]])==1)/length(na.omit(covariates_wBinMM_RGenos[,colnames(covariates_wBinMM_RGenos) == SNPs[ii]]))
      #add to joint table
      tableRow = cbind(SNPs[[ii]], BinMMHR, RSNPHR, jtp, round(freq1,digits = 2))
      BinMMTable = rbind(BinMMTable, tableRow)
    }
    #add SNP to sigSNPs if it passes the significance value
    if(jttestBinMM$`Pr(>Chisq)`[2] <= pvalueSig){
      sigSNPsBinMM = c(sigSNPsBinMM, ii) #add SNP number to sigSNPs
    }
  }
  #write tables out for joint analyses
  write.csv(as.data.frame(IBSTable), file = paste0(outFile,"_IBS.csv"))
  write.csv(as.data.frame(IncompTable), file = paste0(outFile,"_Incomp.csv"))
  write.csv(as.data.frame(AMSTable), file = paste0(outFile,"_AMS.csv"))
  write.csv(as.data.frame(BinMMTable), file = paste0(outFile,"_BinMM.csv"))

  #If joint analyses are significant, run separate univariate analyses
  for(jj in sigSNPsIBS){
    ###########################################
    ##### IBS Mismatch Score
    ###########################################
    IBSModel = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ", SNPs[jj], sep = "")), data = na.omit(covariates_wIBSMismatch_RGenos))
    #pull HR and 95% CI
    summaryStats = summary(IBSModel)
    IBSSNPInfo = round(summaryStats$conf.int[SNPs[[jj]],], digits = 2)
    IBSHR = paste0(IBSSNPInfo[1]," (",IBSSNPInfo[3],", ",IBSSNPInfo[4],")")
    #pull p value
    coefficients = coef(summaryStats)
    pval = round(coefficients[SNPs[[jj]],5], digits = 6)
    #add to joint table
    tableRow = cbind(SNPs[[jj]], "IBS Mismatch Score", IBSHR, pval)
    univarTable = rbind(univarTable, tableRow)
  }
  for(jj in sigSNPsIncomp){
    ###########################################
    ##### Incompatibility Score
    ###########################################
    IncompModel = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ", SNPs[jj], sep = "")), data = na.omit(covariates_wIncomp_RGenos))
    #pull HR and 95% CI
    summaryStats = summary(IncompModel)
    IncompSNPInfo = round(summaryStats$conf.int[SNPs[[jj]],], digits = 2)
    IncompHR = paste0(IncompSNPInfo[1]," (",IncompSNPInfo[3],", ",IncompSNPInfo[4],")")
    #pull p value
    coefficients = coef(summaryStats)
    pval = round(coefficients[SNPs[[jj]],5], digits = 6)
    #add to joint table
    tableRow = cbind(SNPs[[jj]], "Incompatibility Score", IncompHR, pval)
    univarTable = rbind(univarTable, tableRow)
  }
  for(jj in sigSNPsAMS){
    ###########################################
    ##### Allogenomics Mismatch Score
    ###########################################
    AMSModel = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ", SNPs[jj], sep = "")), data = na.omit(covariates_wAMS_RGenos))
    #pull HR and 95% CI
    summaryStats = summary(AMSModel)
    AMSSNPInfo = round(summaryStats$conf.int[SNPs[[jj]],], digits = 2)
    AMSHR = paste0(AMSSNPInfo[1]," (",AMSSNPInfo[3],", ",AMSSNPInfo[4],")")
    #pull p value
    coefficients = coef(summaryStats)
    pval = round(coefficients[SNPs[[jj]],5], digits = 6)
    #add to joint table
    tableRow = cbind(SNPs[[jj]], "AMS Score", AMSHR, pval)
    univarTable = rbind(univarTable, tableRow)
  }
  for(jj in sigSNPsBinMM){
    ###########################################
    ##### Allogenomics Mismatch Score
    ###########################################
    BinMMModel = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ", SNPs[jj], sep = "")), data = na.omit(covariates_wBinMM_RGenos))
    #pull HR and 95% CI
    summaryStats = summary(BinMMModel)
    BinMMSNPInfo = round(summaryStats$conf.int[SNPs[[jj]],], digits = 2)
    BinMMHR = paste0(BinMMSNPInfo[1]," (",BinMMSNPInfo[3],", ",BinMMSNPInfo[4],")")
    #pull p value
    coefficients = coef(summaryStats)
    pval = round(coefficients[SNPs[[jj]],5], digits = 6)
    #add to joint table
    tableRow = cbind(SNPs[[jj]], "Binary Mismatch Score", BinMMHR, pval)
    univarTable = rbind(univarTable, tableRow)
  }
  #we only want R geno once, even if multiple Scores are significant,
  #so determine that there are significant SNPs
  if(length(sigSNPsAMS)>0 | length(sigSNPsIncomp)>0 |length(sigSNPsBinMM)>0 |length(sigSNPsIBS)>0){
    for(jj in unique(c(sigSNPsAMS, sigSNPsBinMM, sigSNPsIBS, sigSNPsIncomp))){
      #score
      RGenoModel = coxph(formula = as.formula(paste0(survFormulaInfo, " ~ ", allCovariates," + ", RSNPs[jj], sep = "")), data = na.omit(covariates_wRGenos))
      #pull HR and 95% CI
      summaryStats = summary(RGenoModel)
      RGenoSNPInfo = round(summaryStats$conf.int[RSNPs[[jj]],], digits = 2)
      RGenoHR = paste0(RGenoSNPInfo[1]," (",RGenoSNPInfo[3],", ",RGenoSNPInfo[4],")")
      #pull p value
      coefficients = coef(summaryStats)
      pval = round(coefficients[RSNPs[[jj]],5], digits = 6)
      #add to joint table
      tableRow = cbind(SNPs[[jj]], "Recipient Genotype", RGenoHR, pval)
      univarTable = rbind(univarTable, tableRow)
    }
  }

  #write out the univariate analysis results
  write.csv(as.data.frame(univarTable), file = paste0(outFile,"_Univariate_Combined.csv"))

  options(warn = oldw)  
}
