#####
#TNFAIP8L3

effect_size_range = 0.37
numSims = 5000

#TNFAIP8L3 - SS = 1000
#read in MAFs
minorAlleleFreqsTNFAIP8L3 = c(0.355, 0.214, 0.181, 0.392, 0.214, 0.104, 0.210, 0.086, 0.198, 0.352, 0.300, 0.142, 0.063, 0.088, 0.169, 0.143, 0.169, 0.291, 0.139, 0.288, 0.140, 0.066, 0.367, 0.291, 0.174, 0.469, 0.366, 0.154, 0.274, 0.458, 0.293, 0.067, 0.162, 0.053, 0.095, 0.474, 0.169, 0.313, 0.298, 0.410, 0.114, 0.108, 0.330, 0.074, 0.390, 0.362, 0.079, 0.076)

path = paste0("/home/vlynn/Simulating_With_Haps/TNFAIP8L3_MargPower_Betas")
setwd(path)

#create matrix to hold betas
SSFiles15 = matrix(NA, nrow = numSims, ncol = 1)

#place filenames in the matrix
for(jj in 1:numSims){
  SSFiles15[jj,] = paste0("BetasWSE_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,"_Simulation",jj,".csv")
}

#read in all files
#dim is 9 x numSims
allFiles15 = sapply(SSFiles15,read.csv)

#make into a matrix of dimension (numSNPs*5) x (9*numSims)
allFiles15.mat = matrix(unlist(allFiles15), nrow = length(allFiles15[[1]]))

#separate based on true model
RGenos_TNFAIP8L3_1000_15AR = allFiles15.mat[which(allFiles15.mat[,1] == 5),]
IBS_TNFAIP8L3_1000_15AR = allFiles15.mat[which(allFiles15.mat[,1] == 2),]
Incomp_TNFAIP8L3_1000_15AR = allFiles15.mat[which(allFiles15.mat[,1] == 3),]
MM_TNFAIP8L3_1000_15AR = allFiles15.mat[which(allFiles15.mat[,1] == 4),]
BeagleIBD_TNFAIP8L3_1000_15AR = allFiles15.mat[which(allFiles15.mat[,1] == 1),]

#want to pull only 1 snp per MAF
#add MAF to matrices
RGenos_TNFAIP8L3_1000_15AR_MAF = cbind(minorAlleleFreqsTNFAIP8L3, RGenos_TNFAIP8L3_1000_15AR)
IBS_TNFAIP8L3_1000_15AR_MAF = cbind(minorAlleleFreqsTNFAIP8L3, IBS_TNFAIP8L3_1000_15AR)
Incomp_TNFAIP8L3_1000_15AR_MAF = cbind(minorAlleleFreqsTNFAIP8L3, Incomp_TNFAIP8L3_1000_15AR)
MM_TNFAIP8L3_1000_15AR_MAF = cbind(minorAlleleFreqsTNFAIP8L3, MM_TNFAIP8L3_1000_15AR)
BeagleIBD_TNFAIP8L3_1000_15AR_MAF = cbind(minorAlleleFreqsTNFAIP8L3, BeagleIBD_TNFAIP8L3_1000_15AR)

#pull first instance of each MAF
maf005 = which(RGenos_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.05 & RGenos_TNFAIP8L3_1000_15AR_MAF[,1] < 0.056)
maf010 = which(RGenos_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.09 & RGenos_TNFAIP8L3_1000_15AR_MAF[,1] < 0.112)
maf10subset = maf010[seq(1,length(maf010),3)]
maf020 = which(RGenos_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.19 & RGenos_TNFAIP8L3_1000_15AR_MAF[,1] < 0.203)
maf030 = which(RGenos_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.299 & RGenos_TNFAIP8L3_1000_15AR_MAF[,1] < 0.301)
maf040 = which(RGenos_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.393 & RGenos_TNFAIP8L3_1000_15AR_MAF[,1] < 0.42)
maf050 = which(RGenos_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.47 & RGenos_TNFAIP8L3_1000_15AR_MAF[,1] < 0.5)
allMAF = c(maf005, maf10subset, maf020, maf030, maf040, maf050)

#subset to only one of each MAF
RGenos_TNFAIP8L3_1000_15AR_MAFSubset = RGenos_TNFAIP8L3_1000_15AR_MAF[allMAF,]
IBS_TNFAIP8L3_1000_15AR_MAFSubset = IBS_TNFAIP8L3_1000_15AR_MAF[allMAF,]
Incomp_TNFAIP8L3_1000_15AR_MAFSubset = Incomp_TNFAIP8L3_1000_15AR_MAF[allMAF,]
MM_TNFAIP8L3_1000_15AR_MAFSubset = MM_TNFAIP8L3_1000_15AR_MAF[allMAF,]
BeagleIBD_TNFAIP8L3_1000_15AR_MAFSubset = BeagleIBD_TNFAIP8L3_1000_15AR_MAF[allMAF,]

#remove MAF column
RGenos_TNFAIP8L3_1000_15AR_Subset = RGenos_TNFAIP8L3_1000_15AR_MAFSubset[,2:ncol(RGenos_TNFAIP8L3_1000_15AR_MAFSubset)]
IBS_TNFAIP8L3_1000_15AR_Subset = IBS_TNFAIP8L3_1000_15AR_MAFSubset[,2:ncol(IBS_TNFAIP8L3_1000_15AR_MAFSubset)]
Incomp_TNFAIP8L3_1000_15AR_Subset = Incomp_TNFAIP8L3_1000_15AR_MAFSubset[,2:ncol(Incomp_TNFAIP8L3_1000_15AR_MAFSubset)]
MM_TNFAIP8L3_1000_15AR_Subset = MM_TNFAIP8L3_1000_15AR_MAFSubset[,2:ncol(MM_TNFAIP8L3_1000_15AR_MAFSubset)]
BeagleIBD_TNFAIP8L3_1000_15AR_Subset = BeagleIBD_TNFAIP8L3_1000_15AR_MAFSubset[,2:ncol(BeagleIBD_TNFAIP8L3_1000_15AR_MAFSubset)]

#reformat matrices to nrow=6*nSims, ncol = 11
startVal = seq(1,ncol(RGenos_TNFAIP8L3_1000_15AR_Subset),11)
endVal = seq(11,ncol(RGenos_TNFAIP8L3_1000_15AR_Subset),11)

RGenos_All_Matrices = list()
for(ii in 1:length(startVal)){
  RGenos_All_Matrices[[ii]] = RGenos_TNFAIP8L3_1000_15AR_Subset[,startVal[[ii]]:endVal[[ii]]]
}
IBS_All_Matrices = list()
for(ii in 1:length(startVal)){
  IBS_All_Matrices[[ii]] = IBS_TNFAIP8L3_1000_15AR_Subset[,startVal[[ii]]:endVal[[ii]]]
}
Incomp_All_Matrices = list()
for(ii in 1:length(startVal)){
  Incomp_All_Matrices[[ii]] = Incomp_TNFAIP8L3_1000_15AR_Subset[,startVal[[ii]]:endVal[[ii]]]
}
MM_All_Matrices = list()
for(ii in 1:length(startVal)){
  MM_All_Matrices[[ii]] = MM_TNFAIP8L3_1000_15AR_Subset[,startVal[[ii]]:endVal[[ii]]]
}
BeagleIBD_All_Matrices = list()
for(ii in 1:length(startVal)){
  BeagleIBD_All_Matrices[[ii]] = BeagleIBD_TNFAIP8L3_1000_15AR_Subset[,startVal[[ii]]:endVal[[ii]]]
}

RGenos_TNFAIP8L3_1000_15AR_Subset_Reformat = do.call(rbind,RGenos_All_Matrices)
IBS_TNFAIP8L3_1000_15AR_Subset_Reformat = do.call(rbind,IBS_All_Matrices)
Incomp_TNFAIP8L3_1000_15AR_Subset_Reformat = do.call(rbind,Incomp_All_Matrices)
MM_TNFAIP8L3_1000_15AR_Subset_Reformat = do.call(rbind,MM_All_Matrices)
BeagleIBD_TNFAIP8L3_1000_15AR_Subset_Reformat = do.call(rbind,BeagleIBD_All_Matrices)

#add the MAF column back in
pullMAF = minorAlleleFreqsTNFAIP8L3[allMAF]
MAF_rep = rep(pullMAF,numSims)

RGenos_TNFAIP8L3_1000_15AR_Subset_ReformatMAF = cbind(MAF_rep,RGenos_TNFAIP8L3_1000_15AR_Subset_Reformat[,2:ncol(RGenos_TNFAIP8L3_1000_15AR_Subset_Reformat)])
IBS_TNFAIP8L3_1000_15AR_Subset_ReformatMAF = cbind(MAF_rep,IBS_TNFAIP8L3_1000_15AR_Subset_Reformat[,2:ncol(IBS_TNFAIP8L3_1000_15AR_Subset_Reformat)])
Incomp_TNFAIP8L3_1000_15AR_Subset_ReformatMAF = cbind(MAF_rep,Incomp_TNFAIP8L3_1000_15AR_Subset_Reformat[,2:ncol(Incomp_TNFAIP8L3_1000_15AR_Subset_Reformat)])
MM_TNFAIP8L3_1000_15AR_Subset_ReformatMAF = cbind(MAF_rep,MM_TNFAIP8L3_1000_15AR_Subset_Reformat[,2:ncol(MM_TNFAIP8L3_1000_15AR_Subset_Reformat)])
BeagleIBD_TNFAIP8L3_1000_15AR_Subset_ReformatMAF = cbind(MAF_rep,BeagleIBD_TNFAIP8L3_1000_15AR_Subset_Reformat[,2:ncol(BeagleIBD_TNFAIP8L3_1000_15AR_Subset_Reformat)])

#pull estimates only (no SE needed)
estNums = seq(2,11,2)

RGenos_ests = RGenos_TNFAIP8L3_1000_15AR_Subset_ReformatMAF[,c(1,estNums)]
IBS_ests = IBS_TNFAIP8L3_1000_15AR_Subset_ReformatMAF[,c(1,estNums)]
Incomp_ests = Incomp_TNFAIP8L3_1000_15AR_Subset_ReformatMAF[,c(1,estNums)]
MM_ests = MM_TNFAIP8L3_1000_15AR_Subset_ReformatMAF[,c(1,estNums)]
BeagleIBD_ests = BeagleIBD_TNFAIP8L3_1000_15AR_Subset_ReformatMAF[,c(1,estNums)]

#get column averages based on rows with the same MAF only
RGenos_AveBeta_MAF05 = colMeans(RGenos_ests[MAF_rep == 0.053,])
RGenos_AveBeta_MAF10 = colMeans(RGenos_ests[MAF_rep == 0.104,])
RGenos_AveBeta_MAF20 = colMeans(RGenos_ests[MAF_rep == 0.198,])
RGenos_AveBeta_MAF30 = colMeans(RGenos_ests[MAF_rep == 0.300,])
RGenos_AveBeta_MAF40 = colMeans(RGenos_ests[MAF_rep == 0.410,])
RGenos_AveBeta_MAF50 = colMeans(RGenos_ests[MAF_rep == 0.474,])

RGenos_SEBeta_MAF05 = apply(RGenos_ests[MAF_rep == 0.053,],2,sd)
RGenos_SEBeta_MAF10 = apply(RGenos_ests[MAF_rep == 0.104,],2,sd)
RGenos_SEBeta_MAF20 = apply(RGenos_ests[MAF_rep == 0.198,],2,sd)
RGenos_SEBeta_MAF30 = apply(RGenos_ests[MAF_rep == 0.300,],2,sd)
RGenos_SEBeta_MAF40 = apply(RGenos_ests[MAF_rep == 0.410,],2,sd)
RGenos_SEBeta_MAF50 = apply(RGenos_ests[MAF_rep == 0.474,],2,sd)

IBS_AveBeta_MAF05 = colMeans(IBS_ests[MAF_rep == 0.053,])
IBS_AveBeta_MAF10 = colMeans(IBS_ests[MAF_rep == 0.104,])
IBS_AveBeta_MAF20 = colMeans(IBS_ests[MAF_rep == 0.198,])
IBS_AveBeta_MAF30 = colMeans(IBS_ests[MAF_rep == 0.300,])
IBS_AveBeta_MAF40 = colMeans(IBS_ests[MAF_rep == 0.410,])
IBS_AveBeta_MAF50 = colMeans(IBS_ests[MAF_rep == 0.474,])

IBS_SEBeta_MAF05 = apply(IBS_ests[MAF_rep == 0.053,],2,sd)
IBS_SEBeta_MAF10 = apply(IBS_ests[MAF_rep == 0.104,],2,sd)
IBS_SEBeta_MAF20 = apply(IBS_ests[MAF_rep == 0.198,],2,sd)
IBS_SEBeta_MAF30 = apply(IBS_ests[MAF_rep == 0.300,],2,sd)
IBS_SEBeta_MAF40 = apply(IBS_ests[MAF_rep == 0.410,],2,sd)
IBS_SEBeta_MAF50 = apply(IBS_ests[MAF_rep == 0.474,],2,sd)

Incomp_AveBeta_MAF05 = colMeans(Incomp_ests[MAF_rep == 0.053,])
Incomp_AveBeta_MAF10 = colMeans(Incomp_ests[MAF_rep == 0.104,])
Incomp_AveBeta_MAF20 = colMeans(Incomp_ests[MAF_rep == 0.198,])
Incomp_AveBeta_MAF30 = colMeans(Incomp_ests[MAF_rep == 0.300,])
Incomp_AveBeta_MAF40 = colMeans(Incomp_ests[MAF_rep == 0.410,])
Incomp_AveBeta_MAF50 = colMeans(Incomp_ests[MAF_rep == 0.474,])

Incomp_SEBeta_MAF05 = apply(Incomp_ests[MAF_rep == 0.053,],2,sd)
Incomp_SEBeta_MAF10 = apply(Incomp_ests[MAF_rep == 0.104,],2,sd)
Incomp_SEBeta_MAF20 = apply(Incomp_ests[MAF_rep == 0.198,],2,sd)
Incomp_SEBeta_MAF30 = apply(Incomp_ests[MAF_rep == 0.300,],2,sd)
Incomp_SEBeta_MAF40 = apply(Incomp_ests[MAF_rep == 0.410,],2,sd)
Incomp_SEBeta_MAF50 = apply(Incomp_ests[MAF_rep == 0.474,],2,sd)

MM_AveBeta_MAF05 = colMeans(MM_ests[MAF_rep == 0.053,])
MM_AveBeta_MAF10 = colMeans(MM_ests[MAF_rep == 0.104,])
MM_AveBeta_MAF20 = colMeans(MM_ests[MAF_rep == 0.198,])
MM_AveBeta_MAF30 = colMeans(MM_ests[MAF_rep == 0.300,])
MM_AveBeta_MAF40 = colMeans(MM_ests[MAF_rep == 0.410,])
MM_AveBeta_MAF50 = colMeans(MM_ests[MAF_rep == 0.474,])

MM_SEBeta_MAF05 = apply(MM_ests[MAF_rep == 0.053,],2,sd)
MM_SEBeta_MAF10 = apply(MM_ests[MAF_rep == 0.104,],2,sd)
MM_SEBeta_MAF20 = apply(MM_ests[MAF_rep == 0.198,],2,sd)
MM_SEBeta_MAF30 = apply(MM_ests[MAF_rep == 0.300,],2,sd)
MM_SEBeta_MAF40 = apply(MM_ests[MAF_rep == 0.410,],2,sd)
MM_SEBeta_MAF50 = apply(MM_ests[MAF_rep == 0.474,],2,sd)

BeagleIBD_AveBeta_MAF05 = colMeans(BeagleIBD_ests[MAF_rep == 0.053,])
BeagleIBD_AveBeta_MAF10 = colMeans(BeagleIBD_ests[MAF_rep == 0.104,])
BeagleIBD_AveBeta_MAF20 = colMeans(BeagleIBD_ests[MAF_rep == 0.198,])
BeagleIBD_AveBeta_MAF30 = colMeans(BeagleIBD_ests[MAF_rep == 0.300,])
BeagleIBD_AveBeta_MAF40 = colMeans(BeagleIBD_ests[MAF_rep == 0.410,])
BeagleIBD_AveBeta_MAF50 = colMeans(BeagleIBD_ests[MAF_rep == 0.474,])

BeagleIBD_SEBeta_MAF05 = apply(BeagleIBD_ests[MAF_rep == 0.053,],2,sd)
BeagleIBD_SEBeta_MAF10 = apply(BeagleIBD_ests[MAF_rep == 0.104,],2,sd)
BeagleIBD_SEBeta_MAF20 = apply(BeagleIBD_ests[MAF_rep == 0.198,],2,sd)
BeagleIBD_SEBeta_MAF30 = apply(BeagleIBD_ests[MAF_rep == 0.300,],2,sd)
BeagleIBD_SEBeta_MAF40 = apply(BeagleIBD_ests[MAF_rep == 0.410,],2,sd)
BeagleIBD_SEBeta_MAF50 = apply(BeagleIBD_ests[MAF_rep == 0.474,],2,sd)

#calculate 95% CIs for the Beta values
#Lower Bounds
RGenos_LB_MAF05 = RGenos_AveBeta_MAF05 - 1.96*RGenos_SEBeta_MAF05
RGenos_LB_MAF10 = RGenos_AveBeta_MAF10 - 1.96*RGenos_SEBeta_MAF10
RGenos_LB_MAF20 = RGenos_AveBeta_MAF20 - 1.96*RGenos_SEBeta_MAF20
RGenos_LB_MAF30 = RGenos_AveBeta_MAF30 - 1.96*RGenos_SEBeta_MAF30
RGenos_LB_MAF40 = RGenos_AveBeta_MAF40 - 1.96*RGenos_SEBeta_MAF40
RGenos_LB_MAF50 = RGenos_AveBeta_MAF50 - 1.96*RGenos_SEBeta_MAF50

IBS_LB_MAF05 = IBS_AveBeta_MAF05 - 1.96*IBS_SEBeta_MAF05
IBS_LB_MAF10 = IBS_AveBeta_MAF10 - 1.96*IBS_SEBeta_MAF10
IBS_LB_MAF20 = IBS_AveBeta_MAF20 - 1.96*IBS_SEBeta_MAF20
IBS_LB_MAF30 = IBS_AveBeta_MAF30 - 1.96*IBS_SEBeta_MAF30
IBS_LB_MAF40 = IBS_AveBeta_MAF40 - 1.96*IBS_SEBeta_MAF40
IBS_LB_MAF50 = IBS_AveBeta_MAF50 - 1.96*IBS_SEBeta_MAF50

Incomp_LB_MAF05 = Incomp_AveBeta_MAF05 - 1.96*Incomp_SEBeta_MAF05
Incomp_LB_MAF10 = Incomp_AveBeta_MAF10 - 1.96*Incomp_SEBeta_MAF10
Incomp_LB_MAF20 = Incomp_AveBeta_MAF20 - 1.96*Incomp_SEBeta_MAF20
Incomp_LB_MAF30 = Incomp_AveBeta_MAF30 - 1.96*Incomp_SEBeta_MAF30
Incomp_LB_MAF40 = Incomp_AveBeta_MAF40 - 1.96*Incomp_SEBeta_MAF40
Incomp_LB_MAF50 = Incomp_AveBeta_MAF50 - 1.96*Incomp_SEBeta_MAF50

MM_LB_MAF05 = MM_AveBeta_MAF05 - 1.96*MM_SEBeta_MAF05
MM_LB_MAF10 = MM_AveBeta_MAF10 - 1.96*MM_SEBeta_MAF10
MM_LB_MAF20 = MM_AveBeta_MAF20 - 1.96*MM_SEBeta_MAF20
MM_LB_MAF30 = MM_AveBeta_MAF30 - 1.96*MM_SEBeta_MAF30
MM_LB_MAF40 = MM_AveBeta_MAF40 - 1.96*MM_SEBeta_MAF40
MM_LB_MAF50 = MM_AveBeta_MAF50 - 1.96*MM_SEBeta_MAF50

BeagleIBD_LB_MAF05 = BeagleIBD_AveBeta_MAF05 - 1.96*BeagleIBD_SEBeta_MAF05
BeagleIBD_LB_MAF10 = BeagleIBD_AveBeta_MAF10 - 1.96*BeagleIBD_SEBeta_MAF10
BeagleIBD_LB_MAF20 = BeagleIBD_AveBeta_MAF20 - 1.96*BeagleIBD_SEBeta_MAF20
BeagleIBD_LB_MAF30 = BeagleIBD_AveBeta_MAF30 - 1.96*BeagleIBD_SEBeta_MAF30
BeagleIBD_LB_MAF40 = BeagleIBD_AveBeta_MAF40 - 1.96*BeagleIBD_SEBeta_MAF40
BeagleIBD_LB_MAF50 = BeagleIBD_AveBeta_MAF50 - 1.96*BeagleIBD_SEBeta_MAF50

#Upper Bounds
RGenos_UB_MAF05 = RGenos_AveBeta_MAF05 + 1.96*RGenos_SEBeta_MAF05
RGenos_UB_MAF10 = RGenos_AveBeta_MAF10 + 1.96*RGenos_SEBeta_MAF10
RGenos_UB_MAF20 = RGenos_AveBeta_MAF20 + 1.96*RGenos_SEBeta_MAF20
RGenos_UB_MAF30 = RGenos_AveBeta_MAF30 + 1.96*RGenos_SEBeta_MAF30
RGenos_UB_MAF40 = RGenos_AveBeta_MAF40 + 1.96*RGenos_SEBeta_MAF40
RGenos_UB_MAF50 = RGenos_AveBeta_MAF50 + 1.96*RGenos_SEBeta_MAF50

IBS_UB_MAF05 = IBS_AveBeta_MAF05 + 1.96*IBS_SEBeta_MAF05
IBS_UB_MAF10 = IBS_AveBeta_MAF10 + 1.96*IBS_SEBeta_MAF10
IBS_UB_MAF20 = IBS_AveBeta_MAF20 + 1.96*IBS_SEBeta_MAF20
IBS_UB_MAF30 = IBS_AveBeta_MAF30 + 1.96*IBS_SEBeta_MAF30
IBS_UB_MAF40 = IBS_AveBeta_MAF40 + 1.96*IBS_SEBeta_MAF40
IBS_UB_MAF50 = IBS_AveBeta_MAF50 + 1.96*IBS_SEBeta_MAF50

Incomp_UB_MAF05 = Incomp_AveBeta_MAF05 + 1.96*Incomp_SEBeta_MAF05
Incomp_UB_MAF10 = Incomp_AveBeta_MAF10 + 1.96*Incomp_SEBeta_MAF10
Incomp_UB_MAF20 = Incomp_AveBeta_MAF20 + 1.96*Incomp_SEBeta_MAF20
Incomp_UB_MAF30 = Incomp_AveBeta_MAF30 + 1.96*Incomp_SEBeta_MAF30
Incomp_UB_MAF40 = Incomp_AveBeta_MAF40 + 1.96*Incomp_SEBeta_MAF40
Incomp_UB_MAF50 = Incomp_AveBeta_MAF50 + 1.96*Incomp_SEBeta_MAF50

MM_UB_MAF05 = MM_AveBeta_MAF05 + 1.96*MM_SEBeta_MAF05
MM_UB_MAF10 = MM_AveBeta_MAF10 + 1.96*MM_SEBeta_MAF10
MM_UB_MAF20 = MM_AveBeta_MAF20 + 1.96*MM_SEBeta_MAF20
MM_UB_MAF30 = MM_AveBeta_MAF30 + 1.96*MM_SEBeta_MAF30
MM_UB_MAF40 = MM_AveBeta_MAF40 + 1.96*MM_SEBeta_MAF40
MM_UB_MAF50 = MM_AveBeta_MAF50 + 1.96*MM_SEBeta_MAF50

BeagleIBD_UB_MAF05 = BeagleIBD_AveBeta_MAF05 + 1.96*BeagleIBD_SEBeta_MAF05
BeagleIBD_UB_MAF10 = BeagleIBD_AveBeta_MAF10 + 1.96*BeagleIBD_SEBeta_MAF10
BeagleIBD_UB_MAF20 = BeagleIBD_AveBeta_MAF20 + 1.96*BeagleIBD_SEBeta_MAF20
BeagleIBD_UB_MAF30 = BeagleIBD_AveBeta_MAF30 + 1.96*BeagleIBD_SEBeta_MAF30
BeagleIBD_UB_MAF40 = BeagleIBD_AveBeta_MAF40 + 1.96*BeagleIBD_SEBeta_MAF40
BeagleIBD_UB_MAF50 = BeagleIBD_AveBeta_MAF50 + 1.96*BeagleIBD_SEBeta_MAF50

#bind mean, lb, ub based for all into single file (for each true model?)
RGenos_CIs = rbind(RGenos_AveBeta_MAF05,RGenos_LB_MAF05,RGenos_UB_MAF05,
                  RGenos_AveBeta_MAF10,RGenos_LB_MAF10,RGenos_UB_MAF10,
                  RGenos_AveBeta_MAF20,RGenos_LB_MAF20,RGenos_UB_MAF20,
                  RGenos_AveBeta_MAF30,RGenos_LB_MAF30,RGenos_UB_MAF30,
                  RGenos_AveBeta_MAF40,RGenos_LB_MAF40,RGenos_UB_MAF40,
                  RGenos_AveBeta_MAF50,RGenos_LB_MAF50,RGenos_UB_MAF50)
write.csv(RGenos_CIs, paste0("RGenos_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,".csv"))

IBS_CIs = rbind(IBS_AveBeta_MAF05,IBS_LB_MAF05,IBS_UB_MAF05,
                   IBS_AveBeta_MAF10,IBS_LB_MAF10,IBS_UB_MAF10,
                   IBS_AveBeta_MAF20,IBS_LB_MAF20,IBS_UB_MAF20,
                   IBS_AveBeta_MAF30,IBS_LB_MAF30,IBS_UB_MAF30,
                   IBS_AveBeta_MAF40,IBS_LB_MAF40,IBS_UB_MAF40,
                   IBS_AveBeta_MAF50,IBS_LB_MAF50,IBS_UB_MAF50)

write.csv(IBS_CIs, paste0("IBS_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,".csv"))

Incomp_CIs = rbind(Incomp_AveBeta_MAF05,Incomp_LB_MAF05,Incomp_UB_MAF05,
                   Incomp_AveBeta_MAF10,Incomp_LB_MAF10,Incomp_UB_MAF10,
                   Incomp_AveBeta_MAF20,Incomp_LB_MAF20,Incomp_UB_MAF20,
                   Incomp_AveBeta_MAF30,Incomp_LB_MAF30,Incomp_UB_MAF30,
                   Incomp_AveBeta_MAF40,Incomp_LB_MAF40,Incomp_UB_MAF40,
                   Incomp_AveBeta_MAF50,Incomp_LB_MAF50,Incomp_UB_MAF50)

write.csv(Incomp_CIs, paste0("Incomp_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,".csv"))

MM_CIs = rbind(MM_AveBeta_MAF05,MM_LB_MAF05,MM_UB_MAF05,
                   MM_AveBeta_MAF10,MM_LB_MAF10,MM_UB_MAF10,
                   MM_AveBeta_MAF20,MM_LB_MAF20,MM_UB_MAF20,
                   MM_AveBeta_MAF30,MM_LB_MAF30,MM_UB_MAF30,
                   MM_AveBeta_MAF40,MM_LB_MAF40,MM_UB_MAF40,
                   MM_AveBeta_MAF50,MM_LB_MAF50,MM_UB_MAF50)

write.csv(MM_CIs, paste0("MM_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,".csv"))

BeagleIBD_CIs = rbind(BeagleIBD_AveBeta_MAF05,BeagleIBD_LB_MAF05,BeagleIBD_UB_MAF05,
                   BeagleIBD_AveBeta_MAF10,BeagleIBD_LB_MAF10,BeagleIBD_UB_MAF10,
                   BeagleIBD_AveBeta_MAF20,BeagleIBD_LB_MAF20,BeagleIBD_UB_MAF20,
                   BeagleIBD_AveBeta_MAF30,BeagleIBD_LB_MAF30,BeagleIBD_UB_MAF30,
                   BeagleIBD_AveBeta_MAF40,BeagleIBD_LB_MAF40,BeagleIBD_UB_MAF40,
                   BeagleIBD_AveBeta_MAF50,BeagleIBD_LB_MAF50,BeagleIBD_UB_MAF50)

write.csv(BeagleIBD_CIs, paste0("BeagleIBD_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,".csv"))
