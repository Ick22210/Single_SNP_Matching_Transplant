#####
#TNFAIP8L3

effect_size_range = 0.37
numSims = 5000

#TNFAIP8L3 - SS = 1000
#read in MAFs
minorAlleleFreqsTNFAIP8L3 = c(0.355, 0.214, 0.181, 0.392, 0.214, 0.104, 0.210, 0.086, 0.198, 0.352, 0.300, 0.142, 0.063, 0.088, 0.169, 0.143, 0.169, 0.291, 0.139, 0.288, 0.140, 0.066, 0.367, 0.291, 0.174, 0.469, 0.366, 0.154, 0.274, 0.458, 0.293, 0.067, 0.162, 0.053, 0.095, 0.474, 0.169, 0.313, 0.298, 0.410, 0.114, 0.108, 0.330, 0.074, 0.390, 0.362, 0.079, 0.076)

path = paste0("/home/vlynn/Simulating_With_Haps/TNFAIP8L3_Results_1000Pairs")
setwd(path)

#create matrix to hold betas
SSFiles15 = matrix(NA, nrow = numSims, ncol = 1)

#place filenames in the matrix
for(jj in 1:numSims){
  SSFiles15[jj,] = paste0("BetasWSE_SingleSNPs_MMScores_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,"_Simulation",jj,".csv")
}

#read in all files
#dim is 9 x numSims
allFiles15 = sapply(SSFiles15,read.csv)

#make into a matrix of dimension (numSNPs*5) x (9*numSims)
allFiles15.mat = matrix(unlist(allFiles15), nrow = length(allFiles15[[1]]))

#separate based on true model
MM_TNFAIP8L3_1000_15AR = allFiles15.mat[which(allFiles15.mat[,1] == 2),]
BinMM_TNFAIP8L3_1000_15AR = allFiles15.mat[which(allFiles15.mat[,1] == 1),]

#want to pull only 1 snp per MAF
#add MAF to matrices
MM_TNFAIP8L3_1000_15AR_MAF = cbind(minorAlleleFreqsTNFAIP8L3, MM_TNFAIP8L3_1000_15AR)
BinMM_TNFAIP8L3_1000_15AR_MAF = cbind(minorAlleleFreqsTNFAIP8L3, BinMM_TNFAIP8L3_1000_15AR)

#pull first instance of each MAF
maf005 = which(MM_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.05 & MM_TNFAIP8L3_1000_15AR_MAF[,1] < 0.056)
maf010 = which(MM_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.09 & MM_TNFAIP8L3_1000_15AR_MAF[,1] < 0.112)
maf10subset = maf010[seq(1,length(maf010),3)]
maf020 = which(MM_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.19 & MM_TNFAIP8L3_1000_15AR_MAF[,1] < 0.203)
maf030 = which(MM_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.299 & MM_TNFAIP8L3_1000_15AR_MAF[,1] < 0.301)
maf040 = which(MM_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.393 & MM_TNFAIP8L3_1000_15AR_MAF[,1] < 0.42)
maf050 = which(MM_TNFAIP8L3_1000_15AR_MAF[,1] >= 0.47 & MM_TNFAIP8L3_1000_15AR_MAF[,1] < 0.5)
allMAF = c(maf005, maf10subset, maf020, maf030, maf040, maf050)

#subset to only one of each MAF
MM_TNFAIP8L3_1000_15AR_MAFSubset = MM_TNFAIP8L3_1000_15AR_MAF[allMAF,]
BinMM_TNFAIP8L3_1000_15AR_MAFSubset = BinMM_TNFAIP8L3_1000_15AR_MAF[allMAF,]

#remove MAF column
MM_TNFAIP8L3_1000_15AR_Subset = MM_TNFAIP8L3_1000_15AR_MAFSubset[,2:ncol(MM_TNFAIP8L3_1000_15AR_MAFSubset)]
BinMM_TNFAIP8L3_1000_15AR_Subset = BinMM_TNFAIP8L3_1000_15AR_MAFSubset[,2:ncol(BinMM_TNFAIP8L3_1000_15AR_MAFSubset)]

#reformat matrices to nrow=6*nSims, ncol = 11
startVal = seq(1,ncol(MM_TNFAIP8L3_1000_15AR_Subset),5)
endVal = seq(5,ncol(MM_TNFAIP8L3_1000_15AR_Subset),5)

MM_All_Matrices = list()
BinMM_All_Matrices = list()

for(ii in 1:length(startVal)){
  MM_All_Matrices[[ii]] = MM_TNFAIP8L3_1000_15AR_Subset[,startVal[[ii]]:endVal[[ii]]]
}
for(ii in 1:length(startVal)){
  BinMM_All_Matrices[[ii]] = BinMM_TNFAIP8L3_1000_15AR_Subset[,startVal[[ii]]:endVal[[ii]]]
}

MM_TNFAIP8L3_1000_15AR_Subset_Reformat = do.call(rbind,MM_All_Matrices)
BinMM_TNFAIP8L3_1000_15AR_Subset_Reformat = do.call(rbind,BinMM_All_Matrices)

#add the MAF column back in
pullMAF = minorAlleleFreqsTNFAIP8L3[allMAF]
MAF_rep = rep(pullMAF,numSims)

MM_TNFAIP8L3_1000_15AR_Subset_ReformatMAF = cbind(MAF_rep,MM_TNFAIP8L3_1000_15AR_Subset_Reformat[,2:ncol(MM_TNFAIP8L3_1000_15AR_Subset_Reformat)])

BinMM_TNFAIP8L3_1000_15AR_Subset_ReformatMAF = cbind(MAF_rep,BinMM_TNFAIP8L3_1000_15AR_Subset_Reformat[,2:ncol(BinMM_TNFAIP8L3_1000_15AR_Subset_Reformat)])

#pull estimates only (no SE needed)
estNums = seq(2,4,2)

MM_ests = MM_TNFAIP8L3_1000_15AR_Subset_ReformatMAF[,c(1,estNums)]
BinMM_ests = BinMM_TNFAIP8L3_1000_15AR_Subset_ReformatMAF[,c(1,estNums)]

#get column averages based on rows with the same MAF only
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

BinMM_AveBeta_MAF05 = colMeans(BinMM_ests[MAF_rep == 0.053,])
BinMM_AveBeta_MAF10 = colMeans(BinMM_ests[MAF_rep == 0.104,])
BinMM_AveBeta_MAF20 = colMeans(BinMM_ests[MAF_rep == 0.198,])
BinMM_AveBeta_MAF30 = colMeans(BinMM_ests[MAF_rep == 0.300,])
BinMM_AveBeta_MAF40 = colMeans(BinMM_ests[MAF_rep == 0.410,])
BinMM_AveBeta_MAF50 = colMeans(BinMM_ests[MAF_rep == 0.474,])

BinMM_SEBeta_MAF05 = apply(BinMM_ests[MAF_rep == 0.053,],2,sd)
BinMM_SEBeta_MAF10 = apply(BinMM_ests[MAF_rep == 0.104,],2,sd)
BinMM_SEBeta_MAF20 = apply(BinMM_ests[MAF_rep == 0.198,],2,sd)
BinMM_SEBeta_MAF30 = apply(BinMM_ests[MAF_rep == 0.300,],2,sd)
BinMM_SEBeta_MAF40 = apply(BinMM_ests[MAF_rep == 0.410,],2,sd)
BinMM_SEBeta_MAF50 = apply(BinMM_ests[MAF_rep == 0.474,],2,sd)


#calculate 95% CIs for the Beta values
#Lower Bounds
MM_LB_MAF05 = MM_AveBeta_MAF05 - 1.96*MM_SEBeta_MAF05
MM_LB_MAF10 = MM_AveBeta_MAF10 - 1.96*MM_SEBeta_MAF10
MM_LB_MAF20 = MM_AveBeta_MAF20 - 1.96*MM_SEBeta_MAF20
MM_LB_MAF30 = MM_AveBeta_MAF30 - 1.96*MM_SEBeta_MAF30
MM_LB_MAF40 = MM_AveBeta_MAF40 - 1.96*MM_SEBeta_MAF40
MM_LB_MAF50 = MM_AveBeta_MAF50 - 1.96*MM_SEBeta_MAF50

#Upper Bounds
MM_UB_MAF05 = MM_AveBeta_MAF05 + 1.96*MM_SEBeta_MAF05
MM_UB_MAF10 = MM_AveBeta_MAF10 + 1.96*MM_SEBeta_MAF10
MM_UB_MAF20 = MM_AveBeta_MAF20 + 1.96*MM_SEBeta_MAF20
MM_UB_MAF30 = MM_AveBeta_MAF30 + 1.96*MM_SEBeta_MAF30
MM_UB_MAF40 = MM_AveBeta_MAF40 + 1.96*MM_SEBeta_MAF40
MM_UB_MAF50 = MM_AveBeta_MAF50 + 1.96*MM_SEBeta_MAF50

BinMM_LB_MAF05 = BinMM_AveBeta_MAF05 - 1.96*BinMM_SEBeta_MAF05
BinMM_LB_MAF10 = BinMM_AveBeta_MAF10 - 1.96*BinMM_SEBeta_MAF10
BinMM_LB_MAF20 = BinMM_AveBeta_MAF20 - 1.96*BinMM_SEBeta_MAF20
BinMM_LB_MAF30 = BinMM_AveBeta_MAF30 - 1.96*BinMM_SEBeta_MAF30
BinMM_LB_MAF40 = BinMM_AveBeta_MAF40 - 1.96*BinMM_SEBeta_MAF40
BinMM_LB_MAF50 = BinMM_AveBeta_MAF50 - 1.96*BinMM_SEBeta_MAF50

#Upper Bounds
BinMM_UB_MAF05 = BinMM_AveBeta_MAF05 + 1.96*BinMM_SEBeta_MAF05
BinMM_UB_MAF10 = BinMM_AveBeta_MAF10 + 1.96*BinMM_SEBeta_MAF10
BinMM_UB_MAF20 = BinMM_AveBeta_MAF20 + 1.96*BinMM_SEBeta_MAF20
BinMM_UB_MAF30 = BinMM_AveBeta_MAF30 + 1.96*BinMM_SEBeta_MAF30
BinMM_UB_MAF40 = BinMM_AveBeta_MAF40 + 1.96*BinMM_SEBeta_MAF40
BinMM_UB_MAF50 = BinMM_AveBeta_MAF50 + 1.96*BinMM_SEBeta_MAF50

path = paste0("/home/vlynn/Simulating_With_Haps/TNFAIP8L3_MargPower_Betas")
setwd(path)

#bind mean, lb, ub based for all into single file (for each true model?)
MM_CIs = rbind(MM_AveBeta_MAF05,MM_LB_MAF05,MM_UB_MAF05,
                   MM_AveBeta_MAF10,MM_LB_MAF10,MM_UB_MAF10,
                   MM_AveBeta_MAF20,MM_LB_MAF20,MM_UB_MAF20,
                   MM_AveBeta_MAF30,MM_LB_MAF30,MM_UB_MAF30,
                   MM_AveBeta_MAF40,MM_LB_MAF40,MM_UB_MAF40,
                   MM_AveBeta_MAF50,MM_LB_MAF50,MM_UB_MAF50)

write.csv(MM_CIs, paste0("MM_MismatchOnly_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,".csv"))

BinMM_CIs = rbind(BinMM_AveBeta_MAF05,BinMM_LB_MAF05,BinMM_UB_MAF05,
                   BinMM_AveBeta_MAF10,BinMM_LB_MAF10,BinMM_UB_MAF10,
                   BinMM_AveBeta_MAF20,BinMM_LB_MAF20,BinMM_UB_MAF20,
                   BinMM_AveBeta_MAF30,BinMM_LB_MAF30,BinMM_UB_MAF30,
                   BinMM_AveBeta_MAF40,BinMM_LB_MAF40,BinMM_UB_MAF40,
                   BinMM_AveBeta_MAF50,BinMM_LB_MAF50,BinMM_UB_MAF50)

write.csv(BinMM_CIs, paste0("BinMM_MismatchOnly_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range,".csv"))
