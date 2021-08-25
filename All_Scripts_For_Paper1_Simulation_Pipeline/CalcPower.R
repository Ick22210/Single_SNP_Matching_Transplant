#######################
# Calc Power
#######################

calcPower = function(x,numSims){
  sum(x <= 0.05)/numSims 
}

#to obtain the arguments from the bash files
args <- commandArgs()

chr <- args[6]
numSamples <- args[7]
numSims <- args[8]
gene <- args[9]

numSamples = as.integer(numSamples)
numSims = as.integer(numSims)

#setwd to location of p-value data
path = paste0("/home/vlynn/Simulating_With_Haps/",gene,"_Results_",numSamples,"Pairs")
setwd(path)

########################
# Single SNP Analysis
########################
#read in data:
#loop through effect sizes
effect_size_range = c(-0.37, -0.30, -0.22, -0.14, 0.14, 0.22, 0.30, 0.37)
length_effect_size = length(effect_size_range)

for(ii in 1:length_effect_size){
  #read in files for effect size ii
  SSFiles15 = matrix(NA, nrow = numSims, ncol = 1)
  
  for(jj in 1:numSims){
      SSFiles15[jj,1] = paste0("PValues_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range[ii],"_Simulation",jj,".csv")
  }
  
  #read in all files
  #dim is 6 x numSims
  allFiles15 = sapply(SSFiles15,read.csv)
  
  #make into a matrix of dimension (numSNPs*5) x (9*numSims)
  allFiles15.mat = matrix(unlist(allFiles15), ncol= (6*numSims), nrow = length(allFiles15[[1]]))
  
  #Order is
  #RGeno, IBSSNP, IncompSNP, MMSNP, BeagleSNP for columns
  #label in first column is whatever was used to actually make the phenotypes
  
  #pull out the rows that were used to make each of the models
  #DGenosSeq = seq(2,length(allFiles15),11)
  RGenosSeq = seq(2,length(allFiles15),6)
  IBSSeq = seq(3,length(allFiles15),6)
  IncompSeq = seq(4,length(allFiles15),6)
  MMSeq = seq(5,length(allFiles15),6)
  BeagleIBDSeq = seq(6,length(allFiles15),6)
  
  #allFiles15.DGenos = allFiles15.mat[,c(1,DGenosSeq)]
  allFiles15.RGenos = allFiles15.mat[,c(1,RGenosSeq)]
  allFiles15.IBS = allFiles15.mat[,c(1,IBSSeq)]
  allFiles15.Incomp = allFiles15.mat[,c(1,IncompSeq)]
  allFiles15.MM = allFiles15.mat[,c(1,MMSeq)]
  allFiles15.BeagleIBD = allFiles15.mat[,c(1,BeagleIBDSeq)]
  
  #separate based on SNP #
  total_length = nrow(allFiles15.RGenos)
  numSNPs = total_length/5
  
  SNPsStartSeq = seq(1,total_length,5)
  SNPsEndSeq = seq(5,total_length,5)
  
  #rows for these are in the order: DGenos, RGenos, IBS, Incomp, Beagle IBD
  #these are the values that were used to make the phenotypes
  #allFiles15.DGenos.SNPs = list()
  allFiles15.RGenos.SNPs = list()
  allFiles15.IBS.SNPs = list()
  allFiles15.Incomp.SNPs = list()
  allFiles15.MM.SNPs = list()
  allFiles15.BeagleIBD.SNPs = list()
  
  for(jj in 1:numSNPs){
    #allFiles15.DGenos.SNPs[[jj]] = allFiles15.DGenos[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles15.RGenos.SNPs[[jj]] = allFiles15.RGenos[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles15.IBS.SNPs[[jj]] = allFiles15.IBS[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles15.Incomp.SNPs[[jj]] = allFiles15.Incomp[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles15.MM.SNPs[[jj]] = allFiles15.MM[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles15.BeagleIBD.SNPs[[jj]] = allFiles15.BeagleIBD[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
}
  
  #then for each row, determine how many p-values are <= 0.05
  #Power is looking at Pr(rejecting H0|HA is true), so we need those sims that reject H0
  mat_power = matrix(NA, nrow=1, ncol = 5)
  mat_power.list = rep(list(mat_power),numSNPs)
  
  #allFiles15.DGenos.SNP.Power = mat_power.list
  allFiles15.RGenos.SNP.Power = mat_power.list
  allFiles15.IBS.SNP.Power = mat_power.list
  allFiles15.Incomp.SNP.Power = mat_power.list
  allFiles15.MM.SNP.Power = mat_power.list
  allFiles15.BeagleIBD.SNP.Power = mat_power.list
  
  for(jj in 1:numSNPs){
      #allFiles15.DGenos.SNP.Power[[jj]] = apply(allFiles15.DGenos.SNPs[[jj]], 1, calcPower, numSims = numSims)
      allFiles15.RGenos.SNP.Power[[jj]] = apply(allFiles15.RGenos.SNPs[[jj]], 1, calcPower, numSims = numSims)
      allFiles15.IBS.SNP.Power[[jj]] = apply(allFiles15.IBS.SNPs[[jj]], 1, calcPower, numSims = numSims)
      allFiles15.Incomp.SNP.Power[[jj]] = apply(allFiles15.Incomp.SNPs[[jj]], 1, calcPower, numSims = numSims)
      allFiles15.MM.SNP.Power[[jj]] = apply(allFiles15.MM.SNPs[[jj]], 1, calcPower, numSims = numSims)
      allFiles15.BeagleIBD.SNP.Power[[jj]] = apply(allFiles15.BeagleIBD.SNPs[[jj]], 1, calcPower, numSims = numSims)
    }
  
  #make everything a matrix again 
  #dim numSNPs x 5 (RGeno, IBS, Incomp, Beagle)
  #allFiles15.DGenos.Power.mat = matrix(unlist(allFiles15.DGenos.SNP.Power), nrow = numSNPs, ncol = 6, byrow = TRUE)
  allFiles15.RGenos.Power.mat = matrix(unlist(allFiles15.RGenos.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  allFiles15.IBS.Power.mat = matrix(unlist(allFiles15.IBS.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  allFiles15.Incomp.Power.mat = matrix(unlist(allFiles15.Incomp.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  allFiles15.MM.Power.mat = matrix(unlist(allFiles15.MM.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  allFiles15.BeagleIBD.Power.mat = matrix(unlist(allFiles15.BeagleIBD.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  
  #write files out
  #write.csv(allFiles15.DGenos.Power.mat, paste0("DGenotypes_Power_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles15.RGenos.Power.mat, paste0("RGenotypes_Power_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles15.IBS.Power.mat, paste0("IBSScoresSS_Power_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles15.Incomp.Power.mat, paste0("IncompScoresSS_Power_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles15.MM.Power.mat, paste0("MMScoresSS_Power_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles15.BeagleIBD.Power.mat, paste0("BeagleIBDProbs_Power_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
}

###########################################################################################################3
#Outcome Incidence 30%

for(ii in 1:length_effect_size){
  #read in files for effect size ii
  SSFiles30 = matrix(NA, nrow = numSims, ncol = 1)
  
  for(jj in 1:numSims){
    SSFiles30[jj,1] = paste0("PValues_SingleSNPs_30PercentAR_BinaryPheno_EffectSize_",effect_size_range[ii],"_Simulation",jj,".csv")
  }
  
  #read in all files
  #dim is 9 x numSims
  allFiles30 = sapply(SSFiles30,read.csv)
  
  #make into a matrix of dimension (numSNPs*5) x (9*numSims)
  allFiles30.mat = matrix(unlist(allFiles30), ncol= (6*numSims), nrow = length(allFiles30[[1]]))
  
  #Order is
  #RGeno, IBSSNP, IncompSNP, MMSNP, BeagleSNP for columns
  #label in first column is whatever was used to actually make the phenotypes
  
  #pull out the rows that were used to make each of the models
  #DGenosSeq = seq(2,length(allFiles30),11)
  RGenosSeq = seq(2,length(allFiles30),6)
  IBSSeq = seq(3,length(allFiles30),6)
  IncompSeq = seq(4,length(allFiles30),6)
  MMSeq = seq(5,length(allFiles30),6)
  BeagleIBDSeq = seq(6,length(allFiles30),6)
  
  #allFiles30.DGenos = allFiles30.mat[,c(1,DGenosSeq)]
  allFiles30.RGenos = allFiles30.mat[,c(1,RGenosSeq)]
  allFiles30.IBS = allFiles30.mat[,c(1,IBSSeq)]
  allFiles30.Incomp = allFiles30.mat[,c(1,IncompSeq)]
  allFiles30.MM = allFiles30.mat[,c(1,MMSeq)]
  allFiles30.BeagleIBD = allFiles30.mat[,c(1,BeagleIBDSeq)]
  
  #rows for these are in the order: DGenos, RGenos, IBS, Incomp, Beagle IBD
  #these are the values that were used to make the phenotypes
  #allFiles30.DGenos.SNPs = list()
  allFiles30.RGenos.SNPs = list()
  allFiles30.IBS.SNPs = list()
  allFiles30.Incomp.SNPs = list()
  allFiles30.MM.SNPs = list()
  allFiles30.BeagleIBD.SNPs = list()
  
  for(jj in 1:numSNPs){
    #allFiles30.DGenos.SNPs[[jj]] = allFiles30.DGenos[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles30.RGenos.SNPs[[jj]] = allFiles30.RGenos[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles30.IBS.SNPs[[jj]] = allFiles30.IBS[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles30.Incomp.SNPs[[jj]] = allFiles30.Incomp[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles30.MM.SNPs[[jj]] = allFiles30.MM[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
    allFiles30.BeagleIBD.SNPs[[jj]] = allFiles30.BeagleIBD[SNPsStartSeq[jj]:SNPsEndSeq[jj],]
  }
  
  #then for each row, determine how many p-values are <= 0.05
  #Power is looking at Pr(rejecting H0|HA is true), so we need those sims that reject H0
  mat_power = matrix(NA, nrow=1, ncol = 5)
  mat_power.list = rep(list(mat_power),numSNPs)
  
  #allFiles30.DGenos.SNP.Power = mat_power.list
  allFiles30.RGenos.SNP.Power = mat_power.list
  allFiles30.IBS.SNP.Power = mat_power.list
  allFiles30.Incomp.SNP.Power = mat_power.list
  allFiles30.MM.SNP.Power = mat_power.list
  allFiles30.BeagleIBD.SNP.Power = mat_power.list
  
  for(jj in 1:numSNPs){
    #allFiles30.DGenos.SNP.Power[[jj]] = apply(allFiles30.DGenos.SNPs[[jj]], 1, calcPower, numSims = numSims)
    allFiles30.RGenos.SNP.Power[[jj]] = apply(allFiles30.RGenos.SNPs[[jj]], 1, calcPower, numSims = numSims)
    allFiles30.IBS.SNP.Power[[jj]] = apply(allFiles30.IBS.SNPs[[jj]], 1, calcPower, numSims = numSims)
    allFiles30.Incomp.SNP.Power[[jj]] = apply(allFiles30.Incomp.SNPs[[jj]], 1, calcPower, numSims = numSims)
    allFiles30.MM.SNP.Power[[jj]] = apply(allFiles30.MM.SNPs[[jj]], 1, calcPower, numSims = numSims)
    allFiles30.BeagleIBD.SNP.Power[[jj]] = apply(allFiles30.BeagleIBD.SNPs[[jj]], 1, calcPower, numSims = numSims)
  }
  
  #make everything a matrix again 
  #dim numSNPs x 5 (DGeno, RGeno, IBS, Incomp, Beagle)
  #allFiles30.DGenos.Power.mat = matrix(unlist(allFiles30.DGenos.SNP.Power), nrow = numSNPs, ncol = 6, byrow = TRUE)
  allFiles30.RGenos.Power.mat = matrix(unlist(allFiles30.RGenos.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  allFiles30.IBS.Power.mat = matrix(unlist(allFiles30.IBS.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  allFiles30.Incomp.Power.mat = matrix(unlist(allFiles30.Incomp.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  allFiles30.MM.Power.mat = matrix(unlist(allFiles30.MM.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  allFiles30.BeagleIBD.Power.mat = matrix(unlist(allFiles30.BeagleIBD.SNP.Power), nrow = numSNPs, ncol=5, byrow = TRUE)
  
  #write files out
  #write.csv(allFiles30.DGenos.Power.mat, paste0("DGenotypes_Power_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles30.RGenos.Power.mat, paste0("RGenotypes_Power_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles30.IBS.Power.mat, paste0("IBSScoresSS_Power_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles30.Incomp.Power.mat, paste0("IncompScoresSS_Power_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles30.MM.Power.mat, paste0("MMScoresSS_Power_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles30.BeagleIBD.Power.mat, paste0("BeagleIBDProbs_Power_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
}