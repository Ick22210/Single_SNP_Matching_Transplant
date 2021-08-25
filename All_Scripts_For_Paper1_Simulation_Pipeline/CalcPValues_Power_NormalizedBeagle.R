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
#Focus on 15 or 30% for Power analysis since 40 is a bit large 
prob1 = 0.15
prob2 = 0.3

#define effect size
#maybe start smaller with larger step?
effect_size_range = c(-0.37, -0.30, -0.22, -0.14, 0.14, 0.22, 0.30, 0.37)

#define numSNPS
numSNPs = ncol(D_Genos)

#binary phenotype generated using rbinom,
#with p = Beta_0 + Beta_1*Score
#need to alter Beta_0 such that prevalence of phenotype = prob1 or prob2

for(aa in 1:length(effect_size_range)){
  #define new matrix for phenotypes
  P_RGeno_15incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_RGeno_30incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_IBS_15incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_IBS_30incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_Incomp_15incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_Incomp_30incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_MM_15incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_MM_30incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_BeagleIBD_15incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  P_BeagleIBD_30incidence = matrix(NA,nrow = numSamples,ncol = numSNPs)
  
  for(ii in 1:numSNPs){
    #calculate Score*Beta_1
    #dim is numSamples x 1 (for each SNP)
    RGenos_ScoreBeta1 = R_Genos.mat[,ii] * effect_size_range[aa]
    IBS_ScoreBeta1 = IBS_SingleSNPs.mat[,ii] * effect_size_range[aa]
    Incomp_ScoreBeta1 = Incomp_SingleSNPs.mat[,ii] * effect_size_range[aa]
    MM_ScoreBeta1 = MM_SingleSNPs.mat[,ii] * effect_size_range[aa]
    BeagleIBD_ScoreBeta1 = BeagleIBD_SingleSNPs.mat.norm[,ii] * effect_size_range[aa]
    
    #beta0
    if(aa == 1){ #done
      #15% phenotype prevalence
      RGenos_Beta0_15 = -1
      IBS_Beta0_15 = MM_Beta0_15 = -1.7
      Incomp_Beta0_15 = -1.65
      BeagleIBD_Beta0_15 = -1.75
      #30% phenotype prevalence
      RGenos_Beta0_30 = -0.15
      IBS_Beta0_30 = Incomp_Beta0_30 = MM_Beta0_30 = -0.8
      BeagleIBD_Beta0_30 = -0.85
    }else if(aa == 2){ #done
      #15% phenotype prevalence
      RGenos_Beta0_15 = -1.3
      IBS_Beta0_15 = -1.55
      Incomp_Beta0_15 = -1.6
      MM_Beta0_15 = -1.65
      BeagleIBD_Beta0_15 = -1.75
      #30% phenotype prevalence
      RGenos_Beta0_30 = -0.4
      IBS_Beta0_30 = -0.65
      Incomp_Beta0_30 = -0.7
      MM_Beta0_30 = -0.8
      BeagleIBD_Beta0_30 = -0.85
    }else if(aa == 3){ #done
      #15% phenotype prevalence
      RGenos_Beta0_15 = -1.4
      IBS_Beta0_15 = Incomp_Beta0_15 = -1.6
      MM_Beta0_15 = -1.65
      BeagleIBD_Beta0_15 = -1.75
      #30% phenotype prevalence
      RGenos_Beta0_30 = -0.5
      IBS_Beta0_30 = Incomp_Beta0_30 = -0.7
      MM_Beta0_30 = -0.75
      BeagleIBD_Beta0_30 = -0.85
    }else if(aa == 4){ #done
      #15% phenotype prevalence
      RGenos_Beta0_15 = -1.5
      IBS_Beta0_15 = Incomp_Beta0_15 = -1.65
      MM_Beta0_15 = -1.675
      BeagleIBD_Beta0_15 = -1.72
      #30% phenotype prevalence
      RGenos_Beta0_30 = -0.65
      IBS_Beta0_30 = Incomp_Beta0_30 = -0.75
      MM_Beta0_30 = -0.8   
      BeagleIBD_Beta0_30 = -0.85
    }else if(aa == 5){ #done
      #15% phenotype prevalence
      RGenos_Beta0_15 = -1.95
      IBS_Beta0_15 = Incomp_Beta0_15 = -1.8
      MM_Beta0_15 = -1.78
      BeagleIBD_Beta0_15 = -1.72
      #30% phenotype prevalence
      RGenos_Beta0_30 = -1.05
      IBS_Beta0_30 = Incomp_Beta0_30 = MM_Beta0_30 = -0.9
      BeagleIBD_Beta0_30 = -0.85
    }else if(aa == 6){ #done
      #15% phenotype prevalence
      RGenos_Beta0_15 = -2.05
      IBS_Beta0_15 = Incomp_Beta0_15 = -1.85
      MM_Beta0_15 = -1.8
      BeagleIBD_Beta0_15 = -1.72
      #30% phenotype prevalence
      RGenos_Beta0_30 = -1.2
      IBS_Beta0_30 = Incomp_Beta0_30 = -0.95
      MM_Beta0_30 = -1.0    
      BeagleIBD_Beta0_30 = -0.85
    }else if(aa == 7){ #done
      #15% phenotype prevalence
      RGenos_Beta0_15 = -2.175
      IBS_Beta0_15 = -1.925
      Incomp_Beta0_15 = -1.9
      MM_Beta0_15 = -1.825
      BeagleIBD_Beta0_15 = -1.75
      #30% phenotype prevalence
      RGenos_Beta0_30 = -1.3
      IBS_Beta0_30 = -1.0
      Incomp_Beta0_30 = -1.0
      MM_Beta0_30 = -0.925    
      BeagleIBD_Beta0_30 = -0.85
    }else if(aa == 8){ #done
      #15% phenotype prevalence
      RGenos_Beta0_15 = -2.3
      IBS_Beta0_15 = -2.0
      Incomp_Beta0_15 = -1.95
      MM_Beta0_15 = -1.87
      BeagleIBD_Beta0_15 = -1.75
      #30% phenotype prevalence
      RGenos_Beta0_30 = -1.425
      IBS_Beta0_30 = -1.1
      Incomp_Beta0_30 = -1.05
      MM_Beta0_30 = -0.95
      BeagleIBD_Beta0_30 = -0.85
    }
    
    #calc linear predictors
    linPred_RGenos_15 = RGenos_Beta0_15 + RGenos_ScoreBeta1
    linPred_IBS_15 = IBS_Beta0_15 + IBS_ScoreBeta1
    linPred_Incomp_15 = Incomp_Beta0_15 + Incomp_ScoreBeta1
    linPred_MM_15 = MM_Beta0_15 + MM_ScoreBeta1
    linPred_BeagleIBD_15 = BeagleIBD_Beta0_15 + BeagleIBD_ScoreBeta1
    
    linPred_Rgenos_30 = RGenos_Beta0_30 + RGenos_ScoreBeta1
    linPred_IBS_30 = IBS_Beta0_30 + IBS_ScoreBeta1
    linPred_Incomp_30 = Incomp_Beta0_30 + Incomp_ScoreBeta1
    linPred_MM_30 = MM_Beta0_30 + MM_ScoreBeta1
    linPred_BeagleIBD_30 = BeagleIBD_Beta0_30 + BeagleIBD_ScoreBeta1
    
    #calculate p values for binomials
    p_Rgenos_15 = exp(linPred_RGenos_15)/(1 + exp(linPred_RGenos_15))
    p_IBS_15 = exp(linPred_IBS_15)/(1 + exp(linPred_IBS_15))
    p_Incomp_15 = exp(linPred_Incomp_15)/(1 + exp(linPred_Incomp_15))
    p_MM_15 = exp(linPred_MM_15)/(1 + exp(linPred_MM_15))
    p_BeagleIBD_15 = exp(linPred_BeagleIBD_15)/(1 + exp(linPred_BeagleIBD_15))
    #calculate p values for binomials
    p_Rgenos_30 = exp(linPred_Rgenos_30)/(1 + exp(linPred_Rgenos_30))
    p_IBS_30 = exp(linPred_IBS_30)/(1 + exp(linPred_IBS_30))
    p_Incomp_30 = exp(linPred_Incomp_30)/(1 + exp(linPred_Incomp_30))
    p_MM_30 = exp(linPred_MM_30)/(1 + exp(linPred_MM_30))
    p_BeagleIBD_30 = exp(linPred_BeagleIBD_30)/(1 + exp(linPred_BeagleIBD_30))
    
    #use binomial dist to get phenotypes
    P_RGeno_15incidence[,ii] = rbinom(numSamples,1,p_Rgenos_15)
    P_RGeno_30incidence[,ii] = rbinom(numSamples,1,p_Rgenos_30)
    ##
    P_IBS_15incidence[,ii] = rbinom(numSamples,1,p_IBS_15)
    P_IBS_30incidence[,ii] = rbinom(numSamples,1,p_IBS_30)
    ##
    P_Incomp_15incidence[,ii] = rbinom(numSamples,1,p_Incomp_15)
    P_Incomp_30incidence[,ii] = rbinom(numSamples,1,p_Incomp_30)
    ##
    P_MM_15incidence[,ii] = rbinom(numSamples,1,p_MM_15)
    P_MM_30incidence[,ii] = rbinom(numSamples,1,p_MM_30)
    ##
    P_BeagleIBD_15incidence[,ii] = rbinom(numSamples,1,p_BeagleIBD_15)
    P_BeagleIBD_30incidence[,ii] = rbinom(numSamples,1,p_BeagleIBD_30)
    
    P_RGeno_15incidence = as.data.frame(P_RGeno_15incidence)
    P_RGeno_30incidence = as.data.frame(P_RGeno_30incidence)
    P_IBS_15incidence = as.data.frame(P_IBS_15incidence)
    P_IBS_30incidence = as.data.frame(P_IBS_30incidence)
    P_Incomp_15incidence = as.data.frame(P_Incomp_15incidence)
    P_Incomp_30incidence = as.data.frame(P_Incomp_30incidence)
    P_MM_15incidence = as.data.frame(P_MM_15incidence)
    P_MM_30incidence = as.data.frame(P_MM_30incidence)
    P_BeagleIBD_15incidence = as.data.frame(P_BeagleIBD_15incidence)
    P_BeagleIBD_30incidence = as.data.frame(P_BeagleIBD_30incidence)
  }
  #need fake SNP names for the single SNP analysis, probably
  snpNamesRecip = list()
  snpNamesIBS = list()
  snpNamesIncomp = list()
  snpNamesMM = list()
  snpNamesBeagle = list()
  
  for(ii in 1:numSNPs){snpNamesRecip[ii] = paste0("rs",ii,"_R")}
  for(ii in 1:numSNPs){snpNamesIBS[ii] = paste0("rs",ii,"_IBS")}
  for(ii in 1:numSNPs){snpNamesIncomp[ii] = paste0("rs",ii,"_Incomp")}
  for(ii in 1:numSNPs){snpNamesMM[ii] = paste0("rs",ii,"_MM")}
  for(ii in 1:numSNPs){snpNamesBeagle[ii] = paste0("rs",ii,"_Beagle")}
  
  snpNamesRecip = as.matrix(unlist(snpNamesRecip), nrow = 1)
  snpNamesIBS = as.matrix(unlist(snpNamesIBS), nrow = 1)
  snpNamesIncomp = as.matrix(unlist(snpNamesIncomp), nrow = 1)
  snpNamesMM = as.matrix(unlist(snpNamesMM), nrow = 1)
  snpNamesBeagle = as.matrix(unlist(snpNamesBeagle), nrow = 1)
  
  #want to predict AR using scores, start simple with logistic regression
  # Single SNP Analysis
  #make list to store all betas and p vals
  SS_BetaswSEAll_15Percent = list()
  SS_PvalsAll_15Percent = list()
  SS_BetaswSEAll_30Percent = list()
  SS_PvalsAll_30Percent = list()
  for(ii in 1:numSNPs){
    # #True phenotype is generated with D Genotype
    # rejPred15PercentDGenoDGeno = glm(P_DGeno_15incidence[[aa]][,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_DGeno_15incidence)
    # rejPred30PercentDGenoDGeno = glm(P_DGeno_30incidence[[aa]][,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_DGeno_30incidence)
    # rejPred15PercentDGenoRGeno = glm(P_DGeno_15incidence[[aa]][,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_DGeno_15incidence)
    # rejPred30PercentDGenoRGeno = glm(P_DGeno_30incidence[[aa]][,ii] ~ R_Genos[,ii], family=binomial, data = P_DGeno_30incidence)
    # rejPred15PercentDGenoIBSSNP = glm(P_DGeno_15incidence[[aa]][,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_DGeno_15incidence)
    # rejPred30PercentDGenoIBSSNP = glm(P_DGeno_30incidence[[aa]][,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_DGeno_30incidence)
    # rejPred15PercentDGenoIncompSNP = glm(P_DGeno_15incidence[[aa]][,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_DGeno_15incidence)
    # rejPred30PercentDGenoIncompSNP = glm(P_DGeno_30incidence[[aa]][,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_DGeno_30incidence)
    # rejPred15PercentDGenoMMSNP = glm(P_DGeno_15incidence[[aa]][,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_DGeno_15incidence)
    # rejPred30PercentDGenoMMSNP = glm(P_DGeno_30incidence[[aa]][,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_DGeno_30incidence)
    # rejPred15PercentDGenoBeagleSNP = glm(P_DGeno_15incidence[[aa]][,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_DGeno_15incidence)
    # rejPred30PercentDGenoBeagleSNP = glm(P_DGeno_30incidence[[aa]][,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_DGeno_30incidence)
 
    #True phenotype is generated with R Genotype
    # rejPred15PercentRGenoDGeno = glm(P_RGeno_15incidence[[aa]][,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    # rejPred30PercentRGenoDGeno = glm(P_RGeno_30incidence[[aa]][,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    rejPred15PercentRGenoRGeno = glm(P_RGeno_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoRGeno = glm(P_RGeno_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    rejPred15PercentRGenoIBSSNP = glm(P_RGeno_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoIBSSNP = glm(P_RGeno_30incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    rejPred15PercentRGenoIncompSNP = glm(P_RGeno_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoIncompSNP = glm(P_RGeno_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    rejPred15PercentRGenoMMSNP = glm(P_RGeno_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoMMSNP = glm(P_RGeno_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    rejPred15PercentRGenoBeagleSNP = glm(P_RGeno_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoBeagleSNP = glm(P_RGeno_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_RGeno_30incidence)
    
    #True phenotype is generated with Single SNP IBS
    # rejPred15PercentIBSSNPDGeno = glm(P_IBS_15incidence[,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_IBS_15incidence)
    # rejPred30PercentIBSSNPDGeno = glm(P_IBS_30incidence[,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_IBS_30incidence)
    rejPred15PercentIBSSNPRGeno = glm(P_IBS_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPRGeno = glm(P_IBS_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_IBS_30incidence)
    rejPred15PercentIBSSNPIBSSNP = glm(P_IBS_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPIBSSNP = glm(P_IBS_30incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_IBS_30incidence)
    rejPred15PercentIBSSNPIncompSNP = glm(P_IBS_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPIncompSNP = glm(P_IBS_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_IBS_30incidence)
    rejPred15PercentIBSSNPMMSNP = glm(P_IBS_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPMMSNP = glm(P_IBS_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_IBS_30incidence)
    rejPred15PercentIBSSNPBeagleSNP = glm(P_IBS_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPBeagleSNP = glm(P_IBS_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_IBS_30incidence)
    
    #True phenotype is generated with Single SNP Incompatibility
    # rejPred15PercentIncompSNPDGeno = glm(P_Incomp_15incidence[,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    # rejPred30PercentIncompSNPDGeno = glm(P_Incomp_30incidence[,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred15PercentIncompSNPRGeno = glm(P_Incomp_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPRGeno = glm(P_Incomp_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred15PercentIncompSNPIBSSNP = glm(P_Incomp_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPIBSSNP = glm(P_Incomp_30incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred15PercentIncompSNPIncompSNP = glm(P_Incomp_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPIncompSNP = glm(P_Incomp_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred15PercentIncompSNPMMSNP = glm(P_Incomp_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPMMSNP = glm(P_Incomp_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred15PercentIncompSNPBeagleSNP = glm(P_Incomp_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPBeagleSNP = glm(P_Incomp_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_Incomp_15incidence)
    
    #True phenotype is generated with Single SNP Mismatch
    # rejPred15PercentMMSNPDGeno = glm(P_MM_15incidence[,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_MM_15incidence)
    # rejPred30PercentMMSNPDGeno = glm(P_MM_30incidence[,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred15PercentMMSNPRGeno = glm(P_MM_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPRGeno = glm(P_MM_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred15PercentMMSNPIBSSNP = glm(P_MM_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPIBSSNP = glm(P_MM_30incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred15PercentMMSNPIncompSNP = glm(P_MM_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPIncompSNP = glm(P_MM_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred15PercentMMSNPMMSNP = glm(P_MM_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPMMSNP = glm(P_MM_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred15PercentMMSNPBeagleSNP = glm(P_MM_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPBeagleSNP = glm(P_MM_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_MM_15incidence)
    
    #True phenotype is generated with Single SNP Beagle IBD
    # rejPred15PercentBeagleSNPDGeno = glm(P_BeagleIBD_15incidence[,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    # rejPred30PercentBeagleSNPDGeno = glm(P_BeagleIBD_30incidence[,ii] ~ D_Genos.mat[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    rejPred15PercentBeagleSNPRGeno = glm(P_BeagleIBD_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPRGeno = glm(P_BeagleIBD_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    rejPred15PercentBeagleSNPIBSSNP = glm(P_BeagleIBD_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPIBSSNP = glm(P_BeagleIBD_30incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    rejPred15PercentBeagleSNPIncompSNP = glm(P_BeagleIBD_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPIncompSNP = glm(P_BeagleIBD_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    rejPred15PercentBeagleSNPMMSNP = glm(P_BeagleIBD_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPMMSNP = glm(P_BeagleIBD_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    rejPred15PercentBeagleSNPBeagleSNP = glm(P_BeagleIBD_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPBeagleSNP = glm(P_BeagleIBD_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    
    #need to pull estimate and SE values for later
    # rejPred15Percent_DGeno_BetaswSE = c(summary(rejPred15PercentDGenoDGeno)$coefficients[2,1:2], summary(rejPred15PercentDGenoRGeno)$coefficients[2,1:2], summary(rejPred15PercentDGenoIBSSNP)$coefficients[2,1:2], summary(rejPred15PercentDGenoIncompSNP)$coefficients[2,1:2], summary(rejPred15PercentDGenoMMSNP)$coefficients[2,1:2],summary(rejPred15PercentDGenoBeagleSNP)$coefficients[2,1:2],summary(rejPred15PercentDGenoIBS_Score)$coefficients[2,1:2],summary(rejPred15PercentDGenoIncomp_Score)$coefficients[2,1:2],summary(rejPred15PercentDGenoMM_Score)$coefficients[2,1:2],summary(rejPred15PercentDGenoBeagle_Score)$coefficients[2,1:2])
    rejPred15Percent_RGeno_BetaswSE = c(summary(rejPred15PercentRGenoRGeno)$coefficients[2,1:2], summary(rejPred15PercentRGenoIBSSNP)$coefficients[2,1:2], summary(rejPred15PercentRGenoIncompSNP)$coefficients[2,1:2], summary(rejPred15PercentRGenoMMSNP)$coefficients[2,1:2], summary(rejPred15PercentRGenoBeagleSNP)$coefficients[2,1:2])
    rejPred15Percent_IBSSNP_BetaswSE = c(summary(rejPred15PercentIBSSNPRGeno)$coefficients[2,1:2], summary(rejPred15PercentIBSSNPIBSSNP)$coefficients[2,1:2], summary(rejPred15PercentIBSSNPIncompSNP)$coefficients[2,1:2], summary(rejPred15PercentIBSSNPMMSNP)$coefficients[2,1:2], summary(rejPred15PercentIBSSNPBeagleSNP)$coefficients[2,1:2])
    rejPred15Percent_IncompSNP_BetaswSE = c(summary(rejPred15PercentIncompSNPRGeno)$coefficients[2,1:2], summary(rejPred15PercentIncompSNPIBSSNP)$coefficients[2,1:2], summary(rejPred15PercentIncompSNPIncompSNP)$coefficients[2,1:2], summary(rejPred15PercentIncompSNPMMSNP)$coefficients[2,1:2], summary(rejPred15PercentIncompSNPBeagleSNP)$coefficients[2,1:2])
    rejPred15Percent_MMSNP_BetaswSE = c(summary(rejPred15PercentMMSNPRGeno)$coefficients[2,1:2], summary(rejPred15PercentMMSNPIBSSNP)$coefficients[2,1:2], summary(rejPred15PercentMMSNPIncompSNP)$coefficients[2,1:2], summary(rejPred15PercentMMSNPMMSNP)$coefficients[2,1:2], summary(rejPred15PercentMMSNPBeagleSNP)$coefficients[2,1:2])
    rejPred15Percent_BeagleSNP_BetaswSE = c(summary(rejPred15PercentBeagleSNPRGeno)$coefficients[2,1:2], summary(rejPred15PercentBeagleSNPIBSSNP)$coefficients[2,1:2], summary(rejPred15PercentBeagleSNPIncompSNP)$coefficients[2,1:2], summary(rejPred15PercentBeagleSNPMMSNP)$coefficients[2,1:2], summary(rejPred15PercentBeagleSNPBeagleSNP)$coefficients[2,1:2])
    
    # rejPred30Percent_DGeno_BetaswSE = c(summary(rejPred30PercentDGenoDGeno)$coefficients[2,1:2], summary(rejPred30PercentDGenoRGeno)$coefficients[2,1:2], summary(rejPred30PercentDGenoIBSSNP)$coefficients[2,1:2], summary(rejPred30PercentDGenoIncompSNP)$coefficients[2,1:2], summary(rejPred30PercentDGenoMMSNP)$coefficients[2,1:2],summary(rejPred30PercentDGenoBeagleSNP)$coefficients[2,1:2],summary(rejPred30PercentDGenoIBS_Score)$coefficients[2,1:2],summary(rejPred30PercentDGenoIncomp_Score)$coefficients[2,1:2],summary(rejPred30PercentDGenoMM_Score)$coefficients[2,1:2],summary(rejPred30PercentDGenoBeagle_Score)$coefficients[2,1:2])
    rejPred30Percent_RGeno_BetaswSE = c(summary(rejPred30PercentRGenoRGeno)$coefficients[2,1:2], summary(rejPred30PercentRGenoIBSSNP)$coefficients[2,1:2], summary(rejPred30PercentRGenoIncompSNP)$coefficients[2,1:2], summary(rejPred30PercentRGenoMMSNP)$coefficients[2,1:2], summary(rejPred30PercentRGenoBeagleSNP)$coefficients[2,1:2])
    rejPred30Percent_IBSSNP_BetaswSE = c(summary(rejPred30PercentIBSSNPRGeno)$coefficients[2,1:2], summary(rejPred30PercentIBSSNPIBSSNP)$coefficients[2,1:2], summary(rejPred30PercentIBSSNPIncompSNP)$coefficients[2,1:2], summary(rejPred30PercentIBSSNPMMSNP)$coefficients[2,1:2], summary(rejPred30PercentIBSSNPBeagleSNP)$coefficients[2,1:2])
    rejPred30Percent_IncompSNP_BetaswSE = c(summary(rejPred30PercentIncompSNPRGeno)$coefficients[2,1:2], summary(rejPred30PercentIncompSNPIBSSNP)$coefficients[2,1:2], summary(rejPred30PercentIncompSNPIncompSNP)$coefficients[2,1:2], summary(rejPred30PercentIncompSNPMMSNP)$coefficients[2,1:2], summary(rejPred30PercentIncompSNPBeagleSNP)$coefficients[2,1:2])
    rejPred30Percent_MMSNP_BetaswSE = c(summary(rejPred30PercentMMSNPRGeno)$coefficients[2,1:2], summary(rejPred30PercentMMSNPIBSSNP)$coefficients[2,1:2], summary(rejPred30PercentMMSNPIncompSNP)$coefficients[2,1:2], summary(rejPred30PercentMMSNPMMSNP)$coefficients[2,1:2], summary(rejPred30PercentMMSNPBeagleSNP)$coefficients[2,1:2])
    rejPred30Percent_BeagleSNP_BetaswSE = c(summary(rejPred30PercentBeagleSNPRGeno)$coefficients[2,1:2], summary(rejPred30PercentBeagleSNPIBSSNP)$coefficients[2,1:2], summary(rejPred30PercentBeagleSNPIncompSNP)$coefficients[2,1:2], summary(rejPred30PercentBeagleSNPMMSNP)$coefficients[2,1:2], summary(rejPred30PercentBeagleSNPBeagleSNP)$coefficients[2,1:2])
    
    #stack Betas and SE values
    BetaswSEAll_15Percent = rbind(rejPred15Percent_RGeno_BetaswSE, rejPred15Percent_IBSSNP_BetaswSE, rejPred15Percent_IncompSNP_BetaswSE, rejPred15Percent_MMSNP_BetaswSE, rejPred15Percent_BeagleSNP_BetaswSE)
    BetaswSEAll_30Percent = rbind(rejPred30Percent_RGeno_BetaswSE, rejPred30Percent_IBSSNP_BetaswSE, rejPred30Percent_IncompSNP_BetaswSE, rejPred30Percent_MMSNP_BetaswSE, rejPred30Percent_BeagleSNP_BetaswSE)
    
    SS_BetaswSEAll_15Percent[[ii]] = BetaswSEAll_15Percent
    SS_BetaswSEAll_30Percent[[ii]] = BetaswSEAll_30Percent
    
    #need to pull p-values from summary tables
    # rejPred15Percent_DGeno_Pvals = c(summary(rejPred15PercentDGenoDGeno)$coefficients[2,4], summary(rejPred15PercentDGenoRGeno)$coefficients[2,4], summary(rejPred15PercentDGenoIBSSNP)$coefficients[2,4], summary(rejPred15PercentDGenoIncompSNP)$coefficients[2,4], summary(rejPred15PercentDGenoMMSNP)$coefficients[2,4], summary(rejPred15PercentDGenoBeagleSNP)$coefficients[2,4])
    rejPred15Percent_RGeno_Pvals = c(summary(rejPred15PercentRGenoRGeno)$coefficients[2,4], summary(rejPred15PercentRGenoIBSSNP)$coefficients[2,4], summary(rejPred15PercentRGenoIncompSNP)$coefficients[2,4], summary(rejPred15PercentRGenoMMSNP)$coefficients[2,4], summary(rejPred15PercentRGenoBeagleSNP)$coefficients[2,4])
    rejPred15Percent_IBSSNP_Pvals = c(summary(rejPred15PercentIBSSNPRGeno)$coefficients[2,4], summary(rejPred15PercentIBSSNPIBSSNP)$coefficients[2,4], summary(rejPred15PercentIBSSNPIncompSNP)$coefficients[2,4],summary(rejPred15PercentIBSSNPMMSNP)$coefficients[2,4], summary(rejPred15PercentIBSSNPBeagleSNP)$coefficients[2,4])
    rejPred15Percent_IncompSNP_Pvals = c(summary(rejPred15PercentIncompSNPRGeno)$coefficients[2,4], summary(rejPred15PercentIncompSNPIBSSNP)$coefficients[2,4], summary(rejPred15PercentIncompSNPIncompSNP)$coefficients[2,4], summary(rejPred15PercentIncompSNPMMSNP)$coefficients[2,4],summary(rejPred15PercentIncompSNPBeagleSNP)$coefficients[2,4])
    rejPred15Percent_MMSNP_Pvals = c(summary(rejPred15PercentMMSNPRGeno)$coefficients[2,4], summary(rejPred15PercentMMSNPIBSSNP)$coefficients[2,4], summary(rejPred15PercentMMSNPIncompSNP)$coefficients[2,4], summary(rejPred15PercentMMSNPMMSNP)$coefficients[2,4],summary(rejPred15PercentMMSNPBeagleSNP)$coefficients[2,4])
    rejPred15Percent_BeagleSNP_Pvals = c(summary(rejPred15PercentBeagleSNPRGeno)$coefficients[2,4], summary(rejPred15PercentBeagleSNPIBSSNP)$coefficients[2,4], summary(rejPred15PercentBeagleSNPIncompSNP)$coefficients[2,4], summary(rejPred15PercentBeagleSNPMMSNP)$coefficients[2,4], summary(rejPred15PercentBeagleSNPBeagleSNP)$coefficients[2,4])
    
    # rejPred30Percent_DGeno_Pvals = c(summary(rejPred30PercentDGenoDGeno)$coefficients[2,4], summary(rejPred30PercentDGenoRGeno)$coefficients[2,4], summary(rejPred30PercentDGenoIBSSNP)$coefficients[2,4], summary(rejPred30PercentDGenoIncompSNP)$coefficients[2,4], summary(rejPred30PercentDGenoMMSNP)$coefficients[2,4], summary(rejPred30PercentDGenoBeagleSNP)$coefficients[2,4])
    rejPred30Percent_RGeno_Pvals = c(summary(rejPred30PercentRGenoRGeno)$coefficients[2,4], summary(rejPred30PercentRGenoIBSSNP)$coefficients[2,4], summary(rejPred30PercentRGenoIncompSNP)$coefficients[2,4], summary(rejPred30PercentRGenoMMSNP)$coefficients[2,4], summary(rejPred30PercentRGenoBeagleSNP)$coefficients[2,4])
    rejPred30Percent_IBSSNP_Pvals = c(summary(rejPred30PercentIBSSNPRGeno)$coefficients[2,4], summary(rejPred30PercentIBSSNPIBSSNP)$coefficients[2,4], summary(rejPred30PercentIBSSNPIncompSNP)$coefficients[2,4],summary(rejPred30PercentIBSSNPMMSNP)$coefficients[2,4], summary(rejPred30PercentIBSSNPBeagleSNP)$coefficients[2,4])
    rejPred30Percent_IncompSNP_Pvals = c(summary(rejPred30PercentIncompSNPRGeno)$coefficients[2,4], summary(rejPred30PercentIncompSNPIBSSNP)$coefficients[2,4], summary(rejPred30PercentIncompSNPIncompSNP)$coefficients[2,4], summary(rejPred30PercentIncompSNPMMSNP)$coefficients[2,4],summary(rejPred30PercentIncompSNPBeagleSNP)$coefficients[2,4])
    rejPred30Percent_MMSNP_Pvals = c(summary(rejPred30PercentMMSNPRGeno)$coefficients[2,4], summary(rejPred30PercentMMSNPIBSSNP)$coefficients[2,4], summary(rejPred30PercentMMSNPIncompSNP)$coefficients[2,4], summary(rejPred30PercentMMSNPMMSNP)$coefficients[2,4],summary(rejPred30PercentMMSNPBeagleSNP)$coefficients[2,4])
    rejPred30Percent_BeagleSNP_Pvals = c(summary(rejPred30PercentBeagleSNPRGeno)$coefficients[2,4], summary(rejPred30PercentBeagleSNPIBSSNP)$coefficients[2,4], summary(rejPred30PercentBeagleSNPIncompSNP)$coefficients[2,4], summary(rejPred30PercentBeagleSNPMMSNP)$coefficients[2,4], summary(rejPred30PercentBeagleSNPBeagleSNP)$coefficients[2,4])
    
    #stack p-values
    PvalsAll_15Percent = rbind(rejPred15Percent_RGeno_Pvals, rejPred15Percent_IBSSNP_Pvals, rejPred15Percent_IncompSNP_Pvals, rejPred15Percent_MMSNP_Pvals, rejPred15Percent_BeagleSNP_Pvals)
    PvalsAll_30Percent = rbind(rejPred30Percent_RGeno_Pvals, rejPred30Percent_IBSSNP_Pvals, rejPred30Percent_IncompSNP_Pvals, rejPred30Percent_MMSNP_Pvals, rejPred30Percent_BeagleSNP_Pvals)
    
    SS_PvalsAll_15Percent[[ii]] = PvalsAll_15Percent #10 columns x 5 rows for each effect size
    SS_PvalsAll_30Percent[[ii]] = PvalsAll_30Percent #10 columns x 5 rows for each effect size
  }
  #unlist and make into matrix
  SS_BetaswSEAll_15Percent.mat = do.call(rbind, SS_BetaswSEAll_15Percent)
  SS_BetaswSEAll_30Percent.mat = do.call(rbind, SS_BetaswSEAll_30Percent)
  SS_PvalsAll_15Percent.mat = do.call(rbind, SS_PvalsAll_15Percent)
  SS_PvalsAll_30Percent.mat = do.call(rbind, SS_PvalsAll_30Percent)
  
  #output Betas and SE values
  write.csv(SS_BetaswSEAll_15Percent.mat, paste0("BetasWSE_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_BetaswSEAll_30Percent.mat, paste0("BetasWSE_SingleSNPs_30PercentAR_BinaryPheno_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  
  #output p-values
  write.csv(SS_PvalsAll_15Percent.mat, paste0("PValues_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_PvalsAll_30Percent.mat, paste0("PValues_SingleSNPs_30PercentAR_BinaryPheno_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  
}

