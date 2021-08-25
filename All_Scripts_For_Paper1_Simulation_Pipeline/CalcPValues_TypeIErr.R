########################################
### Calculate p-values (Type I Error)
########################################

#added mismatch score (8/30/19)
#changed phenotype generation to mimic power analysis, but with effect size 0
#changed phenotype generation back to binomial (11/04/2019)

#to obtain the arguments from the bash files (chr, ss, rel)
args <- commandArgs()

chr <- args[6]
numSamples <- args[7]
simNum <- args[8]
gene <- args[9]

numSamples = as.integer(numSamples)

#setwd to location of plink data
path = paste0("/home/vlynn/Simulating_With_Haps/",gene,"_Results_",numSamples,"Pairs")
setwd(path)

########################
# Single SNP Analysis
######################## 
#read in data
D_Genos = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_DonorGenotypes_Simulation",simNum,".csv"))
D_Genos = D_Genos[,2:ncol(D_Genos)] #dim numSamples x nSNPs
D_Genos.mat = as.matrix(D_Genos)
R_Genos = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_RecipientGenotypes_Simulation",simNum,".csv"))
R_Genos = R_Genos[,2:ncol(R_Genos)] #dim numSamples x nSNPs
R_Genos.mat = as.matrix(R_Genos)

IBS_SingleSNPs = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_IBSScores_Unsummed_Simulation",simNum,".csv"))
IBS_SingleSNPs = IBS_SingleSNPs[,2:ncol(IBS_SingleSNPs)] #dim numSamples x nSNPs
IBS_SingleSNPs.mat = as.matrix(IBS_SingleSNPs)
Incomp_SingleSNPs = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_IncompScores_Unsummed_Simulation",simNum,".csv"))
Incomp_SingleSNPs = Incomp_SingleSNPs[,2:ncol(Incomp_SingleSNPs)] #dim numSamples x nSNPs
Incomp_SingleSNPs.mat = as.matrix(Incomp_SingleSNPs)
MM_SingleSNPs = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_LancetMismatch_Unsummed_Simulation",simNum,".csv"))
MM_SingleSNPs = MM_SingleSNPs[,2:ncol(MM_SingleSNPs)] #dim numSamples x nSNPs
MM_SingleSNPs.mat = as.matrix(MM_SingleSNPs)
BeagleIBD_SingleSNPs = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_BeagleIBDProbabilities_Simulation",simNum,".csv"))
BeagleIBD_SingleSNPs.t = t(BeagleIBD_SingleSNPs)
BeagleIBD_SingleSNPs = BeagleIBD_SingleSNPs.t[2:nrow(BeagleIBD_SingleSNPs.t),] #dim numSamples x nSNPs
#make Beagle IBD NAs into 0s
BeagleIBD_SingleSNPs[is.na(BeagleIBD_SingleSNPs)] = 0

#check the beagle
#read in the donor geno subset file
DgenoSubset = read.csv(paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_DonorGenotypesRemovedSNPs_Simulation",simNum,".csv"))
BeagleIBD_SingleSNPs.subset = DgenoSubset$x

if(length(BeagleIBD_SingleSNPs.subset) > 0){
  BeagleIBD_SingleSNPs = BeagleIBD_SingleSNPs[,-BeagleIBD_SingleSNPs.subset]
}else{
  BeagleIBD_SingleSNPs = BeagleIBD_SingleSNPs
}
BeagleIBD_SingleSNPs.mat = as.matrix(BeagleIBD_SingleSNPs) #dim numSamples x nSNPs

#normalize values
BeagleIBD_SingleSNPs.mat.norm = scale(BeagleIBD_SingleSNPs.mat)

#need to assign AR status to each subject
#incidence is between 20-40% according to recent paper (double check with Brendan)
# do same as power analysis
prob1 = 0.15
prob2 = 0.3

#for Type I error, H0 is true, so AR should not be related to scores
#thus assign AR status using random binomial with prob of AR
#define numSNPS
numSNPs = ncol(D_Genos)

rejStatus15Percent = rbinom(numSamples, 1, prob1)
rejStatus30Percent = rbinom(numSamples, 1, prob2)

#merge everything
singleSNPsProb15 = cbind(R_Genos.mat, IBS_SingleSNPs.mat, Incomp_SingleSNPs.mat, MM_SingleSNPs.mat, BeagleIBD_SingleSNPs.mat.norm, rejStatus15Percent)
singleSNPsProb30 = cbind(R_Genos.mat, IBS_SingleSNPs.mat, Incomp_SingleSNPs.mat, MM_SingleSNPs.mat, BeagleIBD_SingleSNPs.mat.norm, rejStatus30Percent)

#need fake SNP names for the single SNP analysis, probably
snpNamesDonor = list()
snpNamesRecip = list()
snpNamesIBS = list()
snpNamesIncomp = list()
snpNamesMM = list()
snpNamesBeagle = list()

for(ii in 1:ncol(D_Genos)){snpNamesDonor[ii] = paste0("rs",ii,"_D")}
for(ii in 1:ncol(R_Genos)){snpNamesRecip[ii] = paste0("rs",ii,"_R")}
for(ii in 1:ncol(IBS_SingleSNPs)){snpNamesIBS[ii] = paste0("rs",ii,"_IBS")}
for(ii in 1:ncol(Incomp_SingleSNPs)){snpNamesIncomp[ii] = paste0("rs",ii,"_Incomp")}
for(ii in 1:ncol(MM_SingleSNPs)){snpNamesMM[ii] = paste0("rs",ii,"_MM")}
for(ii in 1:ncol(BeagleIBD_SingleSNPs)){snpNamesBeagle[ii] = paste0("rs",ii,"_Beagle")}

snpNamesRecip = as.matrix(unlist(snpNamesRecip), nrow = 1)
snpNamesIBS = as.matrix(unlist(snpNamesIBS), nrow = 1)
snpNamesIncomp = as.matrix(unlist(snpNamesIncomp), nrow = 1)
snpNamesMM = as.matrix(unlist(snpNamesMM),nrow = 1)
snpNamesBeagle = as.matrix(unlist(snpNamesBeagle), nrow = 1)

colnames(singleSNPsProb15) = c(snpNamesRecip, snpNamesIBS, snpNamesIncomp, snpNamesMM, snpNamesBeagle,"Phenotype")
colnames(singleSNPsProb30) = c(snpNamesRecip, snpNamesIBS, snpNamesIncomp, snpNamesMM, snpNamesBeagle,"Phenotype")

#want to predict AR using scores, start simple with logistic regression
# Single SNP Analysis
#make list to store all betas and p vals
SS_BetaswSEAll_15Percent = list()
SS_PvalsAll_15Percent = list()
SS_BetaswSEAll_30Percent = list()
SS_PvalsAll_30Percent = list()

for(ii in 1:numSNPs){
  #True phenotype is generated with null model
  rejPred15PercentRGenoRGeno = glm(Phenotype ~ R_Genos.mat[,ii], family=binomial, data = as.data.frame(singleSNPsProb15))
  rejPred30PercentRGenoRGeno = glm(Phenotype ~ R_Genos.mat[,ii], family=binomial, data = as.data.frame(singleSNPsProb30))
  
  rejPred15PercentIBSSNPIBSSNP = glm(Phenotype ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = as.data.frame(singleSNPsProb15))
  rejPred30PercentIBSSNPIBSSNP = glm(Phenotype ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = as.data.frame(singleSNPsProb30))
  
  rejPred15PercentIncompSNPIncompSNP = glm(Phenotype ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = as.data.frame(singleSNPsProb15))
  rejPred30PercentIncompSNPIncompSNP = glm(Phenotype ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = as.data.frame(singleSNPsProb30))
  
  rejPred15PercentMMSNPMMSNP = glm(Phenotype ~ MM_SingleSNPs.mat[,ii], family=binomial, data = as.data.frame(singleSNPsProb15))
  rejPred30PercentMMSNPMMSNP = glm(Phenotype ~ MM_SingleSNPs.mat[,ii], family=binomial, data = as.data.frame(singleSNPsProb30))
  
  rejPred15PercentBeagleSNPBeagleSNP = glm(Phenotype ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = as.data.frame(singleSNPsProb15))
  rejPred30PercentBeagleSNPBeagleSNP = glm(Phenotype ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = as.data.frame(singleSNPsProb30))
  
  #need to pull estimate and SE values for later
  rejPred15Percent_BetaswSE = c(summary(rejPred15PercentRGenoRGeno)$coefficients[2,1:2], summary(rejPred15PercentIBSSNPIBSSNP)$coefficients[2,1:2], summary(rejPred15PercentIncompSNPIncompSNP)$coefficients[2,1:2], summary(rejPred15PercentMMSNPMMSNP)$coefficients[2,1:2],summary(rejPred15PercentBeagleSNPBeagleSNP)$coefficients[2,1:2])
  rejPred30Percent_BetaswSE = c(summary(rejPred30PercentRGenoRGeno)$coefficients[2,1:2], summary(rejPred30PercentIBSSNPIBSSNP)$coefficients[2,1:2], summary(rejPred30PercentIncompSNPIncompSNP)$coefficients[2,1:2], summary(rejPred30PercentMMSNPMMSNP)$coefficients[2,1:2],summary(rejPred30PercentBeagleSNPBeagleSNP)$coefficients[2,1:2])

  #stack Betas and SE values
  SS_BetaswSEAll_15Percent[[ii]] = rejPred15Percent_BetaswSE
  SS_BetaswSEAll_30Percent[[ii]] = rejPred30Percent_BetaswSE
  
  #need to pull p-values from summary tables
  rejPred15Percent_Pvals = c(summary(rejPred15PercentRGenoRGeno)$coefficients[2,4], summary(rejPred15PercentIBSSNPIBSSNP)$coefficients[2,4], summary(rejPred15PercentIncompSNPIncompSNP)$coefficients[2,4], summary(rejPred15PercentMMSNPMMSNP)$coefficients[2,4], summary(rejPred15PercentBeagleSNPBeagleSNP)$coefficients[2,4])
  rejPred30Percent_Pvals = c(summary(rejPred30PercentRGenoRGeno)$coefficients[2,4], summary(rejPred30PercentIBSSNPIBSSNP)$coefficients[2,4], summary(rejPred30PercentIncompSNPIncompSNP)$coefficients[2,4], summary(rejPred30PercentMMSNPMMSNP)$coefficients[2,4], summary(rejPred30PercentBeagleSNPBeagleSNP)$coefficients[2,4])
  
  #stack p-values
  SS_PvalsAll_15Percent[[ii]] = rejPred15Percent_Pvals 
  SS_PvalsAll_30Percent[[ii]] = rejPred30Percent_Pvals 
}
#unlist and make into matrix
SS_BetaswSEAll_15Percent.mat = do.call(rbind, SS_BetaswSEAll_15Percent)
SS_BetaswSEAll_30Percent.mat = do.call(rbind, SS_BetaswSEAll_30Percent)
SS_PvalsAll_15Percent.mat = do.call(rbind, SS_PvalsAll_15Percent)
SS_PvalsAll_30Percent.mat = do.call(rbind, SS_PvalsAll_30Percent)

#output Betas and SE values
write.csv(SS_BetaswSEAll_15Percent.mat, paste0("BetasWSE_",gene,"_15PercentAR_BinomOutcome_SingleSNPs_Simulation",simNum,".csv"))
write.csv(SS_BetaswSEAll_30Percent.mat, paste0("BetasWSE_",gene,"_30PercentAR_BinomOutcome_SingleSNPs_Simulation",simNum,".csv"))
#output p-values
write.csv(SS_PvalsAll_15Percent.mat, paste0("PValues_",gene,"_15PercentAR_BinomOutcome_SingleSNPs_Simulation",simNum,".csv"))
write.csv(SS_PvalsAll_30Percent.mat, paste0("PValues_",gene,"_30PercentAR_BinomOutcome_SingleSNPs_Simulation",simNum,".csv"))
