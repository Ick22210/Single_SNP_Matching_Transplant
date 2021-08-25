#added in the Lancet mismatch score and changed IBS to mismatch score (8/30/19)
#Edited AMS to be 0,1,or 2 to check if that changes results
#load ARTP2 package
library(ARTP2)
library(dplyr)

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

#make lists of names of plink data
bedFile = paste0("Sim_",simNum,"_Var_Chr",chr,"_",numSamples,"Pairs_",gene,".bed")

bimFile = paste0("Sim_",simNum,"_Var_Chr",chr,"_",numSamples,"Pairs_",gene,".bim")

famFile = paste0("Sim_",simNum,"_Var_Chr",chr,"_",numSamples,"Pairs_",gene,".fam")

#read in plink binary files
plinkFile = read.bed(bed = bedFile, bim = bimFile, fam = famFile)

#match pairs are even and odd column of each file
#define variables for number of columns and rows of the df
nrowPlink = nrow(plinkFile) #this is number of generated subjects (2*SS)
ncolPlink = ncol(plinkFile) #this is number of SNPs

##Want to save the genotypes for Ds and Rs
#Ds are odd numbered rows, Rs are even numbered rows
donorGenotypes = plinkFile[seq(1,nrowPlink,2),]
DGenosMat = matrix(unlist(donorGenotypes), ncol = ncolPlink, byrow = F)
#determine if there are any non-variable SNPs and remove them
DGenosMat.subset = vector()
for(ii in 1:ncolPlink){
  if(length(unique(DGenosMat[,ii])) <= 1){
    DGenosMat.subset = c(DGenosMat.subset, ii)
  } 
}
if(length(DGenosMat.subset) > 0){
  DGenosMat = DGenosMat[,-DGenosMat.subset]
}else {
  DGenosMat = DGenosMat
}
write.csv(DGenosMat.subset, paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_DonorGenotypesRemovedSNPs_Simulation",simNum,".csv"))
write.csv(DGenosMat, paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_DonorGenotypes_Simulation",simNum,".csv"))

recipientGenotypes = plinkFile[seq(2,nrowPlink,2),]
RGenosMat = matrix(unlist(recipientGenotypes), ncol = ncolPlink, byrow = F)
#determine if there are any non-variable SNPs and remove them
RGenosMat.subset = vector()
for(ii in 1:ncolPlink){
  if(length(unique(RGenosMat[,ii])) <= 1){
    RGenosMat.subset = c(RGenosMat.subset, ii)
  } 
}
if(length(RGenosMat.subset) > 0){
  RGenosMat = RGenosMat[,-RGenosMat.subset]
}else {
  RGenosMat = RGenosMat
}
write.csv(RGenosMat, paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_RecipientGenotypes_Simulation",simNum,".csv"))

#redefine number of SNPs
ncolPlink = ncol(DGenosMat)

#########################
# Matches and Scores
#########################
#calculate the difference between the two subjects
diffsPlink = abs(DGenosMat - RGenosMat)

#need to convert difference to score
#if diff = 0, score = 0
#if diff = 1, score is unchanged
#if diff = 2, score = 2
scoresPlink = diffsPlink

#need to keep the unsummed scores also
RIBSScoresMat = matrix(unlist(scoresPlink), ncol = ncolPlink, byrow = F)
write.csv(RIBSScoresMat, paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_IBSScores_Unsummed_Simulation",simNum,".csv"))

#######################
## Incomp Score
#######################
#initialize a list of empty dfs with same number of columns as original
plink_mismatch = diffsPlink
for(ii in 1:ncol(diffsPlink)){
  plink_mismatch[diffsPlink[,ii] == 0,ii] = 0
}
for(ii in 1:ncol(diffsPlink)){
  plink_mismatch[diffsPlink[,ii] != 0,ii] = 1
}

#need to keep unsummed scores also
RIncompScoresMat = matrix(unlist(plink_mismatch), ncol = ncolPlink, byrow = F)
write.csv(RIncompScoresMat, paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_IncompScores_Unsummed_Simulation",simNum,".csv"))

#########################
# Matches and Scores
#########################
#mismatch if D has allele not in R
#if R =2, AND D != 2, then mismatch
#if R = 0, and D != 0, then mismatch
mismatchPlink = matrix(0, nrow = nrow(diffsPlink), ncol = ncol(diffsPlink))

mismatchPlink[(DGenosMat == 0) & (RGenosMat == 2)] = 2
mismatchPlink[(DGenosMat == 2) & (RGenosMat == 0)] = 2
mismatchPlink[(DGenosMat == 1) & (RGenosMat == 2)] = 1
mismatchPlink[(DGenosMat == 1) & (RGenosMat == 0)] = 1
mismatchPlink[is.na(DGenosMat) | is.na(RGenosMat)] = NA

rownames(mismatchPlink) = rownames(DGenosMat)
colnames(mismatchPlink) = colnames(DGenosMat)

#need to keep the unsummed scores also
write.csv(mismatchPlink, paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_LancetMismatch_Unsummed_Simulation",simNum,".csv"))
