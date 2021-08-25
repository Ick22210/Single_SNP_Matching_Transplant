#######################
# Calc Type I Error
#######################

calcTypeIErr = function(x,numSims){
 sum(x <= 0.05)/numSims 
}
calcTypeIErrBonferroni = function(x,numSims,numSNPs){
  sum(x <= 0.05/numSNPs)/numSims 
}

#to obtain the arguments from the bash files (chr, ss)
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
##AR 15%
#read in data:
SSFiles = list()
for(ii in 1:numSims){
    SSFiles[ii] = paste0("PValues_",gene,"_15PercentAR_BinomOutcome_SingleSNPs_Simulation",ii,".csv")
}

AllFiles = list()
for(ii in 1:numSims){
  AllFiles[[ii]] = read.csv(SSFiles[[ii]])
}

everyFile = do.call(cbind, AllFiles)

#figure out where each p-value segment is
#DGenosSeq = seq(2,7*numSims,7)
RGenosSeq = seq(2,6*numSims,6)
IBSSeq = seq(3,6*numSims,6)
IncompSeq = seq(4,6*numSims,6)
MismatchSeq = seq(5,6*numSims,6)
BeagleIBDSeq = seq(6,6*numSims,6)

#pull each of the 5 p-value scenarios out
#D Genos
#everyFile.DGenos = everyFile[,DGenosSeq]
#R Genos
everyFile.RGenos = everyFile[,RGenosSeq]
#IBS Score
everyFile.IBS = everyFile[,IBSSeq]
#Incomp Score
everyFile.Incomp = everyFile[,IncompSeq]
#Mismatch Score
everyFile.Mismatch = everyFile[,MismatchSeq]
#BeagleIBD Score
everyFile.BeagleIBD = everyFile[,BeagleIBDSeq]

numSNPs = nrow(everyFile.RGenos)

#then for each row, determine how many p-values are <= 0.05
#D Genos - not corrected
#DGenos.AR15.TIE = apply(everyFile.DGenos, 1, calcTypeIErr, numSims = numSims)

#R Genos - not corrected
RGenos.AR15.TIE = apply(everyFile.RGenos, 1, calcTypeIErr, numSims = numSims)

#IBS Score - not corrected
IBS.AR15.TIE = apply(everyFile.IBS, 1, calcTypeIErr, numSims = numSims)

#Incomp Score
Incomp.AR15.TIE = apply(everyFile.Incomp, 1, calcTypeIErr, numSims = numSims)

#Mismatch Score
Mismatch.AR15.TIE = apply(everyFile.Mismatch, 1, calcTypeIErr, numSims = numSims)

#BeagleIBD Score
BeagleIBD.AR15.TIE = apply(everyFile.BeagleIBD, 1, calcTypeIErr, numSims = numSims)

###Bonf corrected
#D Genos - Bonf corrected
#DGenos.AR15.TIEBonf  = apply(everyFile.DGenos, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#R Genos - Bonf corrected
RGenos.AR15.TIEBonf  = apply(everyFile.RGenos, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#IBS Score - Bonf corrected
IBS.AR15.TIEBonf  = apply(everyFile.IBS, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#Incomp Score - Bonf corrected
Incomp.AR15.TIEBonf  = apply(everyFile.Incomp, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#Mismatch Score - Bonf corrected
Mismatch.AR15.TIEBonf = apply(everyFile.Mismatch, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#BeagleIBD Score - Bonf Corrected
BeagleIBD.AR15.TIEBonf  = apply(everyFile.BeagleIBD, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

##AR 30%
#read in data:
SSFiles = list()
for(ii in 1:numSims){
  SSFiles[ii] = paste0("PValues_",gene,"_30PercentAR_BinomOutcome_SingleSNPs_Simulation",ii,".csv")
}

AllFiles = list()
for(ii in 1:numSims){
  AllFiles[[ii]] = read.csv(SSFiles[[ii]])
}

everyFile = do.call(cbind, AllFiles)

#pull each of the 5 p-value scenarios out
#D Genos
#everyFile.DGenos = everyFile[,DGenosSeq]
#R Genos
everyFile.RGenos = everyFile[,RGenosSeq]
#IBS Score
everyFile.IBS = everyFile[,IBSSeq]
#Incomp Score
everyFile.Incomp = everyFile[,IncompSeq]
#Mismatch Score
everyFile.Mismatch = everyFile[,MismatchSeq]
#BeagleIBD Score
everyFile.BeagleIBD = everyFile[,BeagleIBDSeq]

#numSNPs = nrow(everyFile.DGenos)

#then for each row, determine how many p-values are <= 0.05
#D Genos - not corrected
#DGenos.AR30.TIE = apply(everyFile.DGenos, 1, calcTypeIErr, numSims = numSims)

#R Genos - not corrected
RGenos.AR30.TIE = apply(everyFile.RGenos, 1, calcTypeIErr, numSims = numSims)

#IBS Score - not corrected
IBS.AR30.TIE = apply(everyFile.IBS, 1, calcTypeIErr, numSims = numSims)

#Incomp Score
Incomp.AR30.TIE = apply(everyFile.Incomp, 1, calcTypeIErr, numSims = numSims)

#Mismatch Score
Mismatch.AR30.TIE = apply(everyFile.Mismatch, 1, calcTypeIErr, numSims = numSims)

#BeagleIBD Score
BeagleIBD.AR30.TIE = apply(everyFile.BeagleIBD, 1, calcTypeIErr, numSims = numSims)

###Bonf corrected
#D Genos - Bonf corrected
#DGenos.AR30.TIEBonf  = apply(everyFile.DGenos, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#R Genos - Bonf corrected
RGenos.AR30.TIEBonf  = apply(everyFile.RGenos, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#IBS Score - Bonf corrected
IBS.AR30.TIEBonf  = apply(everyFile.IBS, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#Incomp Score - Bonf corrected
Incomp.AR30.TIEBonf  = apply(everyFile.Incomp, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#Mismatch Score - Bonf corrected
Mismatch.AR30.TIEBonf = apply(everyFile.Mismatch, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#BeagleIBD Score - Bonf Corrected
BeagleIBD.AR30.TIEBonf  = apply(everyFile.BeagleIBD, 1, calcTypeIErrBonferroni, numSims = numSims, numSNPs = numSNPs)

#this gives 10 lists of all the different options
#DGenosTIE.mat = rbind(DGenos.AR15.TIE, DGenos.AR30.TIE)
RGenosTIE.mat = rbind(RGenos.AR15.TIE, RGenos.AR30.TIE)
IBSTIE.mat = rbind(IBS.AR15.TIE, IBS.AR30.TIE)
IncompTIE.mat = rbind(Incomp.AR15.TIE, Incomp.AR30.TIE)
MismatchTIE.mat = rbind(Mismatch.AR15.TIE, Mismatch.AR30.TIE)
BeagleIBDTIE.mat = rbind(BeagleIBD.AR15.TIE, BeagleIBD.AR30.TIE)

#DGenosTIEBonf.mat = rbind(DGenos.AR15.TIEBonf, DGenos.AR30.TIEBonf)
RGenosTIEBonf.mat = rbind(RGenos.AR15.TIEBonf, RGenos.AR30.TIEBonf)
IBSTIEBonf.mat = rbind(IBS.AR15.TIEBonf, IBS.AR30.TIEBonf)
IncompTIEBonf.mat = rbind(Incomp.AR15.TIEBonf, Incomp.AR30.TIEBonf)
MismatchTIEBonf.mat = rbind(Mismatch.AR15.TIEBonf, Mismatch.AR30.TIEBonf)
BeagleIBDTIEBonf.mat = rbind(BeagleIBD.AR15.TIEBonf, BeagleIBD.AR30.TIEBonf)

#write files out
#write.csv(DGenosTIE.mat, paste0("DGenotypes_TypeIError_AR15AndAR30_Uncorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(RGenosTIE.mat, paste0("RGenotypes_TypeIError_AR15AndAR30_Uncorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(IBSTIE.mat, paste0("IBSScoresSS_TypeIError_AR15AndAR30_Uncorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(IncompTIE.mat, paste0("IncompScoresSS_TypeIError_AR15AndAR30_Uncorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(MismatchTIE.mat, paste0("MismatchScoresSS_TypeIError_AR15AndAR30_Uncorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(BeagleIBDTIE.mat, paste0("BeagleIBDProbs_TypeIError_AR15AndAR30_Uncorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))

#write.csv(DGenosTIEBonf.mat, paste0("DGenotypes_TypeIError_AR15AndAR30_BonfCorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(RGenosTIEBonf.mat, paste0("RGenotypes_TypeIError_AR15AndAR30_BonfCorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(IBSTIEBonf.mat, paste0("IBSScoresSS_TypeIError_AR15AndAR30_BonfCorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(IncompTIEBonf.mat, paste0("IncompScoresSS_TypeIError_AR15AndAR30_BonfCorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(MismatchTIEBonf.mat, paste0("MismatchScoresSS_TypeIError_AR15AndAR30_Bonfcorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))
write.csv(BeagleIBDTIEBonf.mat, paste0("BeagleIBDProbs_TypeIError_AR15AndAR30_BonfCorrected_Chr",chr,"_",numSamples,"Pairs",gene,".csv"))