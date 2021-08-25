#######################
# Calc Joint Power
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
    SSFiles15[jj,] = paste0("PValues_ScoreAndJoint_RGeno_SingleSNPs_15PercentAR_EffectSize_",effect_size_range[ii],"_Simulation",jj,".csv")
  }
  
  #read in all files
  #dim is 9 x numSims
  allFiles15 = sapply(SSFiles15,read.csv)
  
  #make into a matrix of dimension (numSNPs*5) x (9*numSims)
  allFiles15.mat = matrix(unlist(allFiles15), ncol= (4*numSims), nrow = length(allFiles15[[1]]))
  
  #Order is
  #IBSSNP, IncompSNP, Mismatch, BeagleSNP for rows

  #pull out the rows that were used to make each of the models
  IBSSeq = seq(1,nrow(allFiles15.mat),4)
  IncompSeq = seq(2,nrow(allFiles15.mat),4)
  MMSeq = seq(3,nrow(allFiles15.mat),4)
  BeagleIBDSeq = seq(4,nrow(allFiles15.mat),4)
  
  allFiles15.IBS = allFiles15.mat[IBSSeq,]
  allFiles15.Incomp = allFiles15.mat[IncompSeq,]
  allFiles15.MM = allFiles15.mat[MMSeq,]
  allFiles15.BeagleIBD = allFiles15.mat[BeagleIBDSeq,]
  
  #separate based on SNP #
  numSNPs = nrow(allFiles15.IBS)
  labelCols = seq(1,ncol(allFiles15.mat),4)
  scoreCols = seq(2,ncol(allFiles15.mat),4)
  rGenoCols = seq(3,ncol(allFiles15.mat),4)
  jointCols = seq(4,ncol(allFiles15.mat),4)
  
  scoreLabelRGeno = c(labelCols, rGenoCols, scoreCols)
  jointLabelRGeno = c(labelCols, jointCols, rGenoCols)
  jointLabelScore = c(labelCols, jointCols, scoreCols)
  
  #rows for these are in the order: IBS, Incomp, MM, Beagle IBD
  #these are the values that were used to make the phenotypes
  allFiles15.IBS.SNPs.ScorePvals = list()
  allFiles15.Incomp.SNPs.ScorePvals = list()
  allFiles15.MM.SNPs.ScorePvals = list() 
  allFiles15.BeagleIBD.SNPs.ScorePvals = list()
  
  allFiles15.IBS.SNPs.RGenoPvals = list()
  allFiles15.Incomp.SNPs.RGenoPvals = list()
  allFiles15.MM.SNPs.RGenoPvals = list() 
  allFiles15.BeagleIBD.SNPs.RGenoPvals = list()
  
  allFiles15.IBS.SNPs.JtPvals = list()
  allFiles15.Incomp.SNPs.JtPvals = list()
  allFiles15.MM.SNPs.JtPvals = list() 
  allFiles15.BeagleIBD.SNPs.JtPvals = list()
  
  for(jj in 1:numSNPs){
    allFiles15.IBS.SNPs.ScorePvals[[jj]] = allFiles15.IBS[jj,-jointLabelRGeno]
    allFiles15.Incomp.SNPs.ScorePvals[[jj]] = allFiles15.Incomp[jj,-jointLabelRGeno]
    allFiles15.MM.SNPs.ScorePvals[[jj]] = allFiles15.MM[jj,-jointLabelRGeno]
    allFiles15.BeagleIBD.SNPs.ScorePvals[[jj]] = allFiles15.BeagleIBD[jj,-jointLabelRGeno]
    
    allFiles15.IBS.SNPs.RGenoPvals[[jj]] = allFiles15.IBS[jj,-jointLabelScore]
    allFiles15.Incomp.SNPs.RGenoPvals[[jj]] = allFiles15.Incomp[jj,-jointLabelScore]
    allFiles15.MM.SNPs.RGenoPvals[[jj]] = allFiles15.MM[jj,-jointLabelScore]
    allFiles15.BeagleIBD.SNPs.RGenoPvals[[jj]] = allFiles15.BeagleIBD[jj,-jointLabelScore]
    
    allFiles15.IBS.SNPs.JtPvals[[jj]] = allFiles15.IBS[jj,-scoreLabelRGeno]
    allFiles15.Incomp.SNPs.JtPvals[[jj]] = allFiles15.Incomp[jj,-scoreLabelRGeno]
    allFiles15.MM.SNPs.JtPvals[[jj]] = allFiles15.MM[jj,-scoreLabelRGeno]
    allFiles15.BeagleIBD.SNPs.JtPvals[[jj]] = allFiles15.BeagleIBD[jj,-scoreLabelRGeno]
  }
  
  #unlist and have 1 column for each SNP
  allFiles15.IBS.SNPs.ScorePvals.mat = matrix(unlist(allFiles15.IBS.SNPs.ScorePvals), ncol = numSNPs)
  allFiles15.Incomp.SNPs.ScorePvals.mat = matrix(unlist(allFiles15.Incomp.SNPs.ScorePvals), ncol = numSNPs)
  allFiles15.MM.SNPs.ScorePvals.mat = matrix(unlist(allFiles15.MM.SNPs.ScorePvals), ncol = numSNPs)
  allFiles15.BeagleIBD.SNPs.ScorePvals.mat = matrix(unlist(allFiles15.BeagleIBD.SNPs.ScorePvals), ncol = numSNPs)
  
  allFiles15.IBS.SNPs.RGenoPvals.mat = matrix(unlist(allFiles15.IBS.SNPs.RGenoPvals), ncol = numSNPs)
  allFiles15.Incomp.SNPs.RGenoPvals.mat = matrix(unlist(allFiles15.Incomp.SNPs.RGenoPvals), ncol = numSNPs)
  allFiles15.MM.SNPs.RGenoPvals.mat = matrix(unlist(allFiles15.MM.SNPs.RGenoPvals), ncol = numSNPs)
  allFiles15.BeagleIBD.SNPs.RGenoPvals.mat = matrix(unlist(allFiles15.BeagleIBD.SNPs.RGenoPvals), ncol = numSNPs)
  
  allFiles15.IBS.SNPs.JtPvals.mat = matrix(unlist(allFiles15.IBS.SNPs.JtPvals), ncol = numSNPs)
  allFiles15.Incomp.SNPs.JtPvals.mat = matrix(unlist(allFiles15.Incomp.SNPs.JtPvals), ncol = numSNPs)
  allFiles15.MM.SNPs.JtPvals.mat = matrix(unlist(allFiles15.MM.SNPs.JtPvals), ncol = numSNPs)
  allFiles15.BeagleIBD.SNPs.JtPvals.mat = matrix(unlist(allFiles15.BeagleIBD.SNPs.JtPvals), ncol = numSNPs)
  
  #then for each column, determine how many p-values are <= 0.05
  #Power is looking at Pr(rejecting H0|HA is true), so we need those sims that reject H0
  mat_power = matrix(NA, nrow=1, ncol = numSNPs)

  allFiles15.IBS.SNP.Score.Power = mat_power
  allFiles15.Incomp.SNP.Score.Power = mat_power
  allFiles15.MM.SNP.Score.Power = mat_power
  allFiles15.BeagleIBD.SNP.Score.Power = mat_power
  
  allFiles15.IBS.SNP.RGeno.Power = mat_power
  allFiles15.Incomp.SNP.RGeno.Power = mat_power
  allFiles15.MM.SNP.RGeno.Power = mat_power
  allFiles15.BeagleIBD.SNP.RGeno.Power = mat_power
  
  allFiles15.IBS.SNP.Joint.Power = mat_power
  allFiles15.Incomp.SNP.Joint.Power = mat_power
  allFiles15.MM.SNP.Joint.Power = mat_power
  allFiles15.BeagleIBD.SNP.Joint.Power = mat_power
  
  allFiles15.IBS.SNP.Score.Power = apply(allFiles15.IBS.SNPs.ScorePvals.mat, 2, calcPower, numSims = numSims)
  allFiles15.Incomp.SNP.Score.Power = apply(allFiles15.Incomp.SNPs.ScorePvals.mat, 2,  calcPower, numSims = numSims)
  allFiles15.MM.SNP.Score.Power = apply(allFiles15.MM.SNPs.ScorePvals.mat, 2, calcPower, numSims = numSims)
  allFiles15.BeagleIBD.SNP.Score.Power = apply(allFiles15.BeagleIBD.SNPs.ScorePvals.mat, 2, calcPower, numSims = numSims)

  allFiles15.IBS.SNP.RGeno.Power = apply(allFiles15.IBS.SNPs.RGenoPvals.mat, 2, calcPower, numSims = numSims)
  allFiles15.Incomp.SNP.RGeno.Power = apply(allFiles15.Incomp.SNPs.RGenoPvals.mat, 2,  calcPower, numSims = numSims)
  allFiles15.MM.SNP.RGeno.Power = apply(allFiles15.MM.SNPs.RGenoPvals.mat, 2, calcPower, numSims = numSims)
  allFiles15.BeagleIBD.SNP.RGeno.Power = apply(allFiles15.BeagleIBD.SNPs.RGenoPvals.mat, 2, calcPower, numSims = numSims)
  
  allFiles15.IBS.SNP.Joint.Power = apply(allFiles15.IBS.SNPs.JtPvals.mat, 2, calcPower, numSims = numSims)
  allFiles15.Incomp.SNP.Joint.Power = apply(allFiles15.Incomp.SNPs.JtPvals.mat, 2, calcPower, numSims = numSims)
  allFiles15.MM.SNP.Joint.Power = apply(allFiles15.MM.SNPs.JtPvals.mat, 2, calcPower, numSims = numSims)
  allFiles15.BeagleIBD.SNP.Joint.Power = apply(allFiles15.BeagleIBD.SNPs.JtPvals.mat, 2, calcPower, numSims = numSims)
  
  #merge matrices? 
  allFiles15.IBS.SNP.All.Power = rbind(allFiles15.IBS.SNP.Score.Power, allFiles15.IBS.SNP.RGeno.Power, allFiles15.IBS.SNP.Joint.Power)
  allFiles15.Incomp.SNP.All.Power = rbind(allFiles15.Incomp.SNP.Score.Power, allFiles15.Incomp.SNP.RGeno.Power, allFiles15.Incomp.SNP.Joint.Power)
  allFiles15.MM.SNP.All.Power = rbind(allFiles15.MM.SNP.Score.Power, allFiles15.MM.SNP.RGeno.Power, allFiles15.MM.SNP.Joint.Power)
  allFiles15.BeagleIBD.All.Power = rbind(allFiles15.BeagleIBD.SNP.Score.Power, allFiles15.BeagleIBD.SNP.RGeno.Power, allFiles15.BeagleIBD.SNP.Joint.Power)
  
  #write files out
  write.csv(allFiles15.IBS.SNP.All.Power, paste0("IBSScoresSS_ScoreAndJointPower_RGeno_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles15.Incomp.SNP.All.Power, paste0("IncompScoresSS_ScoreAndJointPower_RGeno_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles15.MM.SNP.All.Power, paste0("MMScoresSS_ScoreAndJointPower_RGeno_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles15.BeagleIBD.All.Power, paste0("BeagleIBDProbs_ScoreAndJointPower_RGeno_Chr",chr,"_",numSamples,"Pairs_15PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
}

###########################################################################################################3
#AR Incidence 30%

for(ii in 1:length_effect_size){
  #read in files for effect size ii
  SSFiles30 = matrix(NA, nrow = numSims, ncol = 1)
  
  for(jj in 1:numSims){
    SSFiles30[jj,] = paste0("PValues_ScoreAndJoint_RGeno_SingleSNPs_30PercentAR_EffectSize_",effect_size_range[ii],"_Simulation",jj,".csv")
  }
  
  #read in all files
  #dim is 9 x numSims
  allFiles30 = sapply(SSFiles30,read.csv)
  
  #make into a matrix of dimension (numSNPs*5) x (9*numSims)
  allFiles30.mat = matrix(unlist(allFiles30), ncol= (4*numSims), nrow = length(allFiles30[[1]]))
  
  #Order is
  #IBSSNP, IncompSNP, Mismatchm BeagleSNP for rows
  
  #pull out the rows that were used to make each of the models
  IBSSeq = seq(1,nrow(allFiles30.mat),4)
  IncompSeq = seq(2,nrow(allFiles30.mat),4)
  MMSeq = seq(3,nrow(allFiles30.mat),4)
  BeagleIBDSeq = seq(4,nrow(allFiles30.mat),4)
  
  allFiles30.IBS = allFiles30.mat[IBSSeq,]
  allFiles30.Incomp = allFiles30.mat[IncompSeq,]
  allFiles30.MM = allFiles30.mat[MMSeq,]
  allFiles30.BeagleIBD = allFiles30.mat[BeagleIBDSeq,]
  
  #separate based on SNP #
  numSNPs = nrow(allFiles30.IBS)
  labelCols = seq(1,ncol(allFiles30.mat),4)
  scoreCols = seq(2,ncol(allFiles30.mat),4)
  rGenoCols = seq(3,ncol(allFiles30.mat),4)
  jointCols = seq(4,ncol(allFiles30.mat),4)
  
  scoreLabelRGeno = c(labelCols, rGenoCols, scoreCols)
  jointLabelRGeno = c(labelCols, jointCols, rGenoCols)
  jointLabelScore = c(labelCols, jointCols, scoreCols)
  
  #rows for these are in the order: IBS, Incomp, MM, Beagle IBD
  #these are the values that were used to make the phenotypes
  allFiles30.IBS.SNPs.ScorePvals = list()
  allFiles30.Incomp.SNPs.ScorePvals = list()
  allFiles30.MM.SNPs.ScorePvals = list() 
  allFiles30.BeagleIBD.SNPs.ScorePvals = list()
  
  allFiles30.IBS.SNPs.RGenoPvals = list()
  allFiles30.Incomp.SNPs.RGenoPvals = list()
  allFiles30.MM.SNPs.RGenoPvals = list() 
  allFiles30.BeagleIBD.SNPs.RGenoPvals = list()
  
  allFiles30.IBS.SNPs.JtPvals = list()
  allFiles30.Incomp.SNPs.JtPvals = list()
  allFiles30.MM.SNPs.JtPvals = list() 
  allFiles30.BeagleIBD.SNPs.JtPvals = list()
  
  for(jj in 1:numSNPs){
    allFiles30.IBS.SNPs.ScorePvals[[jj]] = allFiles30.IBS[jj,-jointLabelRGeno]
    allFiles30.Incomp.SNPs.ScorePvals[[jj]] = allFiles30.Incomp[jj,-jointLabelRGeno]
    allFiles30.MM.SNPs.ScorePvals[[jj]] = allFiles30.MM[jj,-jointLabelRGeno]
    allFiles30.BeagleIBD.SNPs.ScorePvals[[jj]] = allFiles30.BeagleIBD[jj,-jointLabelRGeno]
    
    allFiles30.IBS.SNPs.RGenoPvals[[jj]] = allFiles30.IBS[jj,-jointLabelScore]
    allFiles30.Incomp.SNPs.RGenoPvals[[jj]] = allFiles30.Incomp[jj,-jointLabelScore]
    allFiles30.MM.SNPs.RGenoPvals[[jj]] = allFiles30.MM[jj,-jointLabelScore]
    allFiles30.BeagleIBD.SNPs.RGenoPvals[[jj]] = allFiles30.BeagleIBD[jj,-jointLabelScore]
    
    allFiles30.IBS.SNPs.JtPvals[[jj]] = allFiles30.IBS[jj,-scoreLabelRGeno]
    allFiles30.Incomp.SNPs.JtPvals[[jj]] = allFiles30.Incomp[jj,-scoreLabelRGeno]
    allFiles30.MM.SNPs.JtPvals[[jj]] = allFiles30.MM[jj,-scoreLabelRGeno]
    allFiles30.BeagleIBD.SNPs.JtPvals[[jj]] = allFiles30.BeagleIBD[jj,-scoreLabelRGeno]
  }
  
  #unlist and have 1 column for each SNP
  allFiles30.IBS.SNPs.ScorePvals.mat = matrix(unlist(allFiles30.IBS.SNPs.ScorePvals), ncol = numSNPs)
  allFiles30.Incomp.SNPs.ScorePvals.mat = matrix(unlist(allFiles30.Incomp.SNPs.ScorePvals), ncol = numSNPs)
  allFiles30.MM.SNPs.ScorePvals.mat = matrix(unlist(allFiles30.MM.SNPs.ScorePvals), ncol = numSNPs)
  allFiles30.BeagleIBD.SNPs.ScorePvals.mat = matrix(unlist(allFiles30.BeagleIBD.SNPs.ScorePvals), ncol = numSNPs)
  
  allFiles30.IBS.SNPs.RGenoPvals.mat = matrix(unlist(allFiles30.IBS.SNPs.RGenoPvals), ncol = numSNPs)
  allFiles30.Incomp.SNPs.RGenoPvals.mat = matrix(unlist(allFiles30.Incomp.SNPs.RGenoPvals), ncol = numSNPs)
  allFiles30.MM.SNPs.RGenoPvals.mat = matrix(unlist(allFiles30.MM.SNPs.RGenoPvals), ncol = numSNPs)
  allFiles30.BeagleIBD.SNPs.RGenoPvals.mat = matrix(unlist(allFiles30.BeagleIBD.SNPs.RGenoPvals), ncol = numSNPs)
  
  allFiles30.IBS.SNPs.JtPvals.mat = matrix(unlist(allFiles30.IBS.SNPs.JtPvals), ncol = numSNPs)
  allFiles30.Incomp.SNPs.JtPvals.mat = matrix(unlist(allFiles30.Incomp.SNPs.JtPvals), ncol = numSNPs)
  allFiles30.MM.SNPs.JtPvals.mat = matrix(unlist(allFiles30.MM.SNPs.JtPvals), ncol = numSNPs)
  allFiles30.BeagleIBD.SNPs.JtPvals.mat = matrix(unlist(allFiles30.BeagleIBD.SNPs.JtPvals), ncol = numSNPs)
  
  #then for each row, determine how many p-values are <= 0.05
  #Power is looking at Pr(rejecting H0|HA is true), so we need those sims that reject H0
  mat_power = matrix(NA, nrow=1, ncol = numSNPs)
  
  allFiles30.IBS.SNP.Score.Power = mat_power
  allFiles30.Incomp.SNP.Score.Power = mat_power
  allFiles30.MM.SNP.Score.Power = mat_power
  allFiles30.BeagleIBD.SNP.Score.Power = mat_power
  
  allFiles30.IBS.SNP.RGeno.Power = mat_power
  allFiles30.Incomp.SNP.RGeno.Power = mat_power
  allFiles30.MM.SNP.RGeno.Power = mat_power
  allFiles30.BeagleIBD.SNP.RGeno.Power = mat_power
  
  allFiles30.IBS.SNP.Joint.Power = mat_power
  allFiles30.Incomp.SNP.Joint.Power = mat_power
  allFiles30.MM.SNP.Joint.Power = mat_power
  allFiles30.BeagleIBD.SNP.Joint.Power = mat_power
  
  allFiles30.IBS.SNP.Score.Power = apply(allFiles30.IBS.SNPs.ScorePvals.mat, 2, calcPower, numSims = numSims)
  allFiles30.Incomp.SNP.Score.Power = apply(allFiles30.Incomp.SNPs.ScorePvals.mat, 2,  calcPower, numSims = numSims)
  allFiles30.MM.SNP.Score.Power = apply(allFiles30.MM.SNPs.ScorePvals.mat, 2, calcPower, numSims = numSims)
  allFiles30.BeagleIBD.SNP.Score.Power = apply(allFiles30.BeagleIBD.SNPs.ScorePvals.mat, 2, calcPower, numSims = numSims)
  
  allFiles30.IBS.SNP.RGeno.Power = apply(allFiles30.IBS.SNPs.RGenoPvals.mat, 2, calcPower, numSims = numSims)
  allFiles30.Incomp.SNP.RGeno.Power = apply(allFiles30.Incomp.SNPs.RGenoPvals.mat, 2,  calcPower, numSims = numSims)
  allFiles30.MM.SNP.RGeno.Power = apply(allFiles30.MM.SNPs.RGenoPvals.mat, 2, calcPower, numSims = numSims)
  allFiles30.BeagleIBD.SNP.RGeno.Power = apply(allFiles30.BeagleIBD.SNPs.RGenoPvals.mat, 2, calcPower, numSims = numSims)
  
  allFiles30.IBS.SNP.Joint.Power = apply(allFiles30.IBS.SNPs.JtPvals.mat, 2, calcPower, numSims = numSims)
  allFiles30.Incomp.SNP.Joint.Power = apply(allFiles30.Incomp.SNPs.JtPvals.mat, 2, calcPower, numSims = numSims)
  allFiles30.MM.SNP.Joint.Power = apply(allFiles30.MM.SNPs.JtPvals.mat, 2, calcPower, numSims = numSims)
  allFiles30.BeagleIBD.SNP.Joint.Power = apply(allFiles30.BeagleIBD.SNPs.JtPvals.mat, 2, calcPower, numSims = numSims)
  
  #merge matrices? 
  allFiles30.IBS.SNP.All.Power = rbind(allFiles30.IBS.SNP.Score.Power, allFiles30.IBS.SNP.RGeno.Power, allFiles30.IBS.SNP.Joint.Power)
  allFiles30.Incomp.SNP.All.Power = rbind(allFiles30.Incomp.SNP.Score.Power, allFiles30.Incomp.SNP.RGeno.Power, allFiles30.Incomp.SNP.Joint.Power)
  allFiles30.MM.SNP.All.Power = rbind(allFiles30.MM.SNP.Score.Power, allFiles30.MM.SNP.RGeno.Power, allFiles30.MM.SNP.Joint.Power)
  allFiles30.BeagleIBD.All.Power = rbind(allFiles30.BeagleIBD.SNP.Score.Power, allFiles30.BeagleIBD.SNP.RGeno.Power, allFiles30.BeagleIBD.SNP.Joint.Power)
  
  #write files out
  write.csv(allFiles30.IBS.SNP.All.Power, paste0("IBSScoresSS_ScoreAndJointPower_RGeno_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles30.Incomp.SNP.All.Power, paste0("IncompScoresSS_ScoreAndJointPower_RGeno_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles30.MM.SNP.All.Power, paste0("MMScoresSS_ScoreAndJointPower_RGeno_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
  write.csv(allFiles30.BeagleIBD.All.Power, paste0("BeagleIBDProbs_ScoreAndJointPower_RGeno_Chr",chr,"_",numSamples,"Pairs_30PercentAR_",gene,"_ES_",effect_size_range[ii],".csv"))
}