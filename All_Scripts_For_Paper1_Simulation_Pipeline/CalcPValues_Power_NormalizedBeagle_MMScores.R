##########################################
### Calculate p-values (Power Analysis)
###########################################

#to obtain the arguments from the bash files (chr, ss, rel)
args <- commandArgs()

chr <- args[6]
numSamples <- args[7]
simNum <- args[8]
gene <- args[9]
RAND <- args[10]

numSamples = as.integer(numSamples)
RAND = as.integer(RAND)

set.seed(RAND)
#setwd to location of plink data
path = paste0("/home/vlynn/Simulating_With_Haps/",gene,"_Results_",numSamples,"Pairs")
setwd(path)

########################
# Single SNP Analysis
######################## 
#read in data
MM_SingleSNPs = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_LancetMismatch_Unsummed_Simulation",simNum,".csv"))
MM_SingleSNPs = MM_SingleSNPs[,2:ncol(MM_SingleSNPs)] #dim numSamples x nSNPs
MM_SingleSNPs.mat = as.matrix(MM_SingleSNPs)

BinMM_SingleSNPs = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_BinaryMismatch_Unsummed_Simulation",simNum,".csv"))
BinMM_SingleSNPs = BinMM_SingleSNPs[,2:ncol(BinMM_SingleSNPs)]
BinMM_SingleSNPs.mat = as.matrix(BinMM_SingleSNPs)

#need to assign AR status to each subject
#incidence is between 20-40% according to recent paper (double check with Brendan)
#Focus on 15 or 30% for Power analysis since 40 is a bit large 
prob1 = 0.15
prob2 = 0.3

#define effect size
#maybe start smaller with larger step?
effect_size_range = c(0.37)

#define numSNPS
numSNPs = ncol(MM_SingleSNPs)

#binary phenotype generated using rbinom,
#with p = Beta_0 + Beta_1*Score
#need to alter Beta_0 such that prevalence of phenotype = prob1 or prob2

#define new matrix for phenotypes
P_MM_15incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
P_MM_30incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
P_BinMM_15incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
P_BinMM_30incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)

for(ii in 1:numSNPs){
  #calculate Score*Beta_1
  #dim is numSamples x 1 (for each SNP)
  MM_ScoreBeta1 = MM_SingleSNPs.mat[,ii] * effect_size_range
  BinMM_ScoreBeta1 = BinMM_SingleSNPs.mat[,ii] * effect_size_range
  #beta0
  #15% phenotype prevalence
  BinMM_Beta0_15 = MM_Beta0_15 = -1.87
  BinMM_Beta0_30 = MM_Beta0_30 = -0.95
  
  #calc linear predictors
  linPred_MM_15 = MM_Beta0_15 + MM_ScoreBeta1
  linPred_BinMM_15 = BinMM_Beta0_15 + BinMM_ScoreBeta1
  
  linPred_MM_30 = MM_Beta0_30 + MM_ScoreBeta1
  linPred_BinMM_30 = BinMM_Beta0_30 + BinMM_ScoreBeta1
  
  #calculate p values for binomials
  p_MM_15 = exp(linPred_MM_15)/(1 + exp(linPred_MM_15))
  p_BinMM_15 = exp(linPred_BinMM_15)/(1 + exp(linPred_BinMM_15))
  
  #calculate p values for binomials
  p_MM_30 = exp(linPred_MM_30)/(1 + exp(linPred_MM_30))
  p_BinMM_30 = exp(linPred_BinMM_30)/(1 + exp(linPred_BinMM_30))
  
  #use binomial dist to get phenotypes
  ##
  P_MM_15incidence[,ii] = rbinom(numSamples,1,p_MM_15)
  P_MM_30incidence[,ii] = rbinom(numSamples,1,p_MM_30)
  ##
  P_BinMM_15incidence[,ii] = rbinom(numSamples,1,p_BinMM_15)
  P_BinMM_30incidence[,ii] = rbinom(numSamples,1,p_BinMM_30)
  
  P_MM_15incidence = as.data.frame(P_MM_15incidence)
  P_MM_30incidence = as.data.frame(P_MM_30incidence)
  
  P_BinMM_15incidence = as.data.frame(P_BinMM_15incidence)
  P_BinMM_30incidence = as.data.frame(P_BinMM_30incidence)
  
}
#need fake SNP names for the single SNP analysis, probably
snpNamesMM = list()
for(ii in 1:numSNPs){snpNamesMM[ii] = paste0("rs",ii,"_MM")}
snpNamesMM = as.matrix(unlist(snpNamesMM), nrow = 1)

#want to predict AR using scores, start simple with logistic regression
# Single SNP Analysis
#make list to store all betas and p vals
SS_BetaswSEAll_15Percent = list()
SS_PvalsAll_15Percent = list()
SS_BetaswSEAll_30Percent = list()
SS_PvalsAll_30Percent = list()
for(ii in 1:numSNPs){
  #True phenotype is generated with Single SNP Mismatch
  rejPred15PercentMMSNPMMSNP = glm(P_MM_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
  rejPred30PercentMMSNPMMSNP = glm(P_MM_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_MM_30incidence)
  rejPred15PercentMMSNPBinMMSNP = glm(P_MM_15incidence[,ii] ~ BinMM_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
  rejPred30PercentMMSNPBinMMSNP = glm(P_MM_30incidence[,ii] ~ BinMM_SingleSNPs.mat[,ii], family=binomial, data = P_MM_30incidence)
  
  #True phenotype is generated with Single SNP Binary Mismatch
  rejPred15PercentBinMMSNPMMSNP = glm(P_BinMM_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_BinMM_15incidence)
  rejPred30PercentBinMMSNPMMSNP = glm(P_BinMM_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_BinMM_30incidence)
  rejPred15PercentBinMMSNPBinMMSNP = glm(P_BinMM_15incidence[,ii] ~ BinMM_SingleSNPs.mat[,ii], family=binomial, data = P_BinMM_15incidence)
  rejPred30PercentBinMMSNPBinMMSNP = glm(P_BinMM_30incidence[,ii] ~ BinMM_SingleSNPs.mat[,ii], family=binomial, data = P_BinMM_30incidence)
  
  #need to pull estimate and SE values for later
  rejPred15Percent_MMSNP_BetaswSE = c(summary(rejPred15PercentMMSNPMMSNP)$coefficients[2,1:2], summary(rejPred15PercentMMSNPBinMMSNP)$coefficients[2,1:2])
  rejPred30Percent_MMSNP_BetaswSE = c(summary(rejPred30PercentMMSNPMMSNP)$coefficients[2,1:2], summary(rejPred30PercentMMSNPBinMMSNP)$coefficients[2,1:2])
  
  rejPred15Percent_BinMMSNP_BetaswSE = c(summary(rejPred15PercentBinMMSNPMMSNP)$coefficients[2,1:2], summary(rejPred15PercentBinMMSNPBinMMSNP)$coefficients[2,1:2])
  rejPred30Percent_BinMMSNP_BetaswSE = c(summary(rejPred30PercentBinMMSNPMMSNP)$coefficients[2,1:2], summary(rejPred30PercentBinMMSNPBinMMSNP)$coefficients[2,1:2])
    
  #stack Betas and SE values
  BetaswSEAll_15Percent = rbind(rejPred15Percent_MMSNP_BetaswSE, rejPred15Percent_BinMMSNP_BetaswSE)
  BetaswSEAll_30Percent = rbind(rejPred30Percent_MMSNP_BetaswSE, rejPred30Percent_BinMMSNP_BetaswSE)
    
  SS_BetaswSEAll_15Percent[[ii]] = BetaswSEAll_15Percent
  SS_BetaswSEAll_30Percent[[ii]] = BetaswSEAll_30Percent  
  
  #need to pull p-values from summary tables
  rejPred15Percent_MMSNP_Pvals = c(summary(rejPred15PercentMMSNPMMSNP)$coefficients[2,4], summary(rejPred15PercentMMSNPBinMMSNP)$coefficients[2,4])
  rejPred30Percent_MMSNP_Pvals = c(summary(rejPred30PercentMMSNPMMSNP)$coefficients[2,4], summary(rejPred30PercentMMSNPBinMMSNP)$coefficients[2,4])
  
  rejPred15Percent_BinMMSNP_Pvals = c(summary(rejPred15PercentBinMMSNPMMSNP)$coefficients[2,4], summary(rejPred15PercentBinMMSNPBinMMSNP)$coefficients[2,4])
  rejPred30Percent_BinMMSNP_Pvals = c(summary(rejPred30PercentBinMMSNPMMSNP)$coefficients[2,4], summary(rejPred30PercentBinMMSNPBinMMSNP)$coefficients[2,4])
  
  #stack p-values
  PvalsAll_15Percent = rbind(rejPred15Percent_MMSNP_Pvals, rejPred15Percent_BinMMSNP_Pvals)
  PvalsAll_30Percent = rbind(rejPred30Percent_MMSNP_Pvals, rejPred30Percent_BinMMSNP_Pvals)
  
  SS_PvalsAll_15Percent[[ii]] = PvalsAll_15Percent #10 columns x 5 rows for each effect size
  SS_PvalsAll_30Percent[[ii]] = PvalsAll_30Percent #10 columns x 5 rows for each effect size
}
#unlist and make into matrix
SS_BetaswSEAll_15Percent.mat = do.call(rbind, SS_BetaswSEAll_15Percent)
SS_BetaswSEAll_30Percent.mat = do.call(rbind, SS_BetaswSEAll_30Percent)
SS_PvalsAll_15Percent.mat = do.call(rbind, SS_PvalsAll_15Percent)
SS_PvalsAll_30Percent.mat = do.call(rbind, SS_PvalsAll_30Percent)

#output Betas
write.csv(SS_BetaswSEAll_15Percent.mat, paste0("BetasWSE_SingleSNPs_MMScores_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,"_Simulation",simNum,".csv"))
write.csv(SS_BetaswSEAll_30Percent.mat, paste0("BetasWSE_SingleSNPs_MMScores_30PercentAR_BinaryPheno_EffectSize_",effect_size_range,"_Simulation",simNum,".csv"))

#output p-values
write.csv(SS_PvalsAll_15Percent.mat, paste0("PValues_SingleSNPs_MMScores_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,"_Simulation",simNum,".csv"))
write.csv(SS_PvalsAll_30Percent.mat, paste0("PValues_SingleSNPs_MMScores_30PercentAR_BinaryPheno_EffectSize_",effect_size_range,"_Simulation",simNum,".csv"))

