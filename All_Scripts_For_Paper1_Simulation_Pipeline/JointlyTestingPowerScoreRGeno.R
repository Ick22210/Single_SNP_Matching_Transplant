##########################################
### Calculate p-values (Power Analysis)
### Jointly testing with Score/RGeno
#Running with AMS as 0,1,or 2 (1/24/20)
###########################################

#load car library
suppressMessages(library(epicalc))

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
#Jinbo said focus on 15 or 30% for Power analysis since 40 is a bit large 
prob1 = 0.15
prob2 = 0.3

#define effect size
#maybe start smaller with larger step?
effect_size_range = c(-0.37, -0.30, -0.22, -0.14, 0.14, 0.22, 0.30, 0.37)

#define numSNPS
numSNPs = ncol(R_Genos)

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
  #want to predict AR using scores, start simple with logistic regression
  # Single SNP Analysis
  #make list to store all betas and p vals
  SS_BetaswSEAll_15Percent = list()
  SS_BetaswSEAll_Rgenos_15Percent = list()
  SS_PvalsAll_15Percent = list()
  SS_BetaswSEAll_30Percent = list()
  SS_BetaswSEAll_Rgenos_30Percent = list()
  SS_PvalsAll_30Percent = list()
  
  for(ii in 1:numSNPs){
    #True phenotype is generated with Single SNP IBS
    ##Code for true model (only correct for Score)
    rejPred15PercentIBSSNPIBSSNP = glm(P_IBS_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPIBSSNP = glm(P_IBS_30incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_IBS_30incidence)
    #model only using R geno to fit
    rejPred15PercentIBSSNPRGeno = glm(P_IBS_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPRGeno = glm(P_IBS_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_IBS_30incidence)
    ##Code for model that also corrects for R Genotype
    rejPred15PercentIBSSNPIBSSNPRGeno = glm(P_IBS_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii] + R_Genos.mat[,ii], family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPIBSSNPRGeno = glm(P_IBS_30incidence[,ii] ~ IBS_SingleSNPs.mat[,ii] + R_Genos.mat[,ii], family=binomial, data = P_IBS_30incidence)
    ##Jointly test significance of IBS Score and R Genos
    #make null models
    rejPred15PercentIBSSNPIBSSNPRGenoNull = glm(P_IBS_15incidence[,ii] ~ 1, family=binomial, data = P_IBS_15incidence)
    rejPred30PercentIBSSNPIBSSNPRGenoNull = glm(P_IBS_30incidence[,ii] ~ 1, family=binomial, data = P_IBS_30incidence)
    ##tells to test against the null model (i.e., only the intercept is not set to 0)
    lhIBS15 = lrtest(rejPred15PercentIBSSNPIBSSNPRGeno, rejPred15PercentIBSSNPIBSSNPRGenoNull)
    lhIBS30 = lrtest(rejPred30PercentIBSSNPIBSSNPRGeno, rejPred30PercentIBSSNPIBSSNPRGenoNull)
    #pull the p-values
    lhIBS15pval = lhIBS15$p.value
    lhIBS30pval = lhIBS30$p.value
    
    #True phenotype is generated with Single SNP Incompatibility
    ##Code for true model (only correct for Score)
    rejPred15PercentIncompSNPIncompSNP = glm(P_Incomp_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPIncompSNP = glm(P_Incomp_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    #model only using R geno to fit
    rejPred15PercentIncompSNPRGeno = glm(P_Incomp_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPRGeno = glm(P_Incomp_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    ##Code for model that also corrects for R Genotype
    rejPred15PercentIncompSNPIncompSNPRGeno = glm(P_Incomp_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii]+ R_Genos.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPIncompSNPRGeno = glm(P_Incomp_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii]+ R_Genos.mat[,ii], family=binomial, data = P_Incomp_15incidence)
    ##Jointly test significance of IBS Score and R Genos
    #make null models
    rejPred15PercentIncompSNPIncompSNPRGenoNull = glm(P_Incomp_15incidence[,ii] ~ 1, family=binomial, data = P_Incomp_15incidence)
    rejPred30PercentIncompSNPIncompSNPRGenoNull = glm(P_Incomp_30incidence[,ii] ~ 1, family=binomial, data = P_Incomp_15incidence)
    ##tells to test against the null model (i.e., only the intercept is not set to 0)
    lhIncomp15 = lrtest(rejPred15PercentIncompSNPIncompSNPRGeno, rejPred15PercentIncompSNPIncompSNPRGenoNull)
    lhIncomp30 = lrtest(rejPred30PercentIncompSNPIncompSNPRGeno, rejPred30PercentIncompSNPIncompSNPRGenoNull)
    #pull the p-values
    lhIncomp15pval = lhIncomp15$p.value
    lhIncomp30pval = lhIncomp30$p.value
    
    #True phenotype is generated with Single SNP Mismatch
    ##Code for true model (only correct for Score)
    rejPred15PercentMMSNPMMSNP = glm(P_MM_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPMMSNP = glm(P_MM_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_MM_15incidence)
    #model only using R geno to fit
    rejPred15PercentMMSNPRGeno = glm(P_MM_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPRGeno = glm(P_MM_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_MM_15incidence)
    ##Code for model that also corrects for R Genotype
    rejPred15PercentMMSNPMMSNPRGeno = glm(P_MM_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii] + R_Genos.mat[,ii], family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPMMSNPRGeno = glm(P_MM_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii] + R_Genos.mat[,ii], family=binomial, data = P_MM_15incidence)
    ##Jointly test significance of IBS Score and R Genos
    #make null models
    rejPred15PercentMMSNPMMSNPRGenoNull = glm(P_MM_15incidence[,ii] ~ 1, family=binomial, data = P_MM_15incidence)
    rejPred30PercentMMSNPMMSNPRGenoNull = glm(P_MM_30incidence[,ii] ~ 1, family=binomial, data = P_MM_15incidence)
    ##tells to test against the null model (i.e., only the intercept is not set to 0)
    lhMM15 = lrtest(rejPred15PercentMMSNPMMSNPRGeno, rejPred15PercentMMSNPMMSNPRGenoNull)
    lhMM30 = lrtest(rejPred30PercentMMSNPMMSNPRGeno, rejPred30PercentMMSNPMMSNPRGenoNull)
    #pull the p-values
    lhMM15pval = lhMM15$p.value
    lhMM30pval = lhMM30$p.value
    
    #True phenotype is generated with Single SNP Beagle IBD
    ##Code for true model (only correct for Score)
    rejPred15PercentBeagleSNPBeagleSNP = glm(P_BeagleIBD_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPBeagleSNP = glm(P_BeagleIBD_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    #model only using R geno to fit
    rejPred15PercentBeagleSNPRGeno = glm(P_BeagleIBD_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPRGeno = glm(P_BeagleIBD_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    ##Code for model that also corrects for R Genotype
    rejPred15PercentBeagleSNPBeagleSNPRGeno = glm(P_BeagleIBD_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii] + R_Genos.mat[,ii], family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPBeagleSNPRGeno = glm(P_BeagleIBD_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii] + R_Genos.mat[,ii], family=binomial, data = P_BeagleIBD_30incidence)
    ##Jointly test significance of IBS Score and R Genos
    #make null models
    rejPred15PercentBeagleSNPBeagleSNPRGenoNull = glm(P_BeagleIBD_15incidence[,ii] ~ 1, family=binomial, data = P_BeagleIBD_15incidence)
    rejPred30PercentBeagleSNPBeagleSNPRGenoNull = glm(P_BeagleIBD_30incidence[,ii] ~ 1, family=binomial, data = P_BeagleIBD_30incidence)
    ##tells to test against the null model (i.e., only the intercept is not set to 0)
    lhBeagleIBD15 = lrtest(rejPred15PercentBeagleSNPBeagleSNPRGeno, rejPred15PercentBeagleSNPBeagleSNPRGenoNull)
    lhBeagleIBD30 = lrtest(rejPred30PercentBeagleSNPBeagleSNPRGeno, rejPred30PercentBeagleSNPBeagleSNPRGenoNull)
    #pull the p-values
    lhBeagleIBD15pval = lhBeagleIBD15$p.value
    lhBeagleIBD30pval = lhBeagleIBD30$p.value
    
    #need to pull estimate for later
    #correct models
    rejPred15Percent_IBSSNP_BetaswSE = c(summary(rejPred15PercentIBSSNPIBSSNP)$coefficients[2,1], summary(rejPred15PercentIBSSNPRGeno)$coefficients[2,1])
    rejPred15Percent_IncompSNP_BetaswSE = c(summary(rejPred15PercentIncompSNPIncompSNP)$coefficients[2,1], summary(rejPred15PercentIncompSNPRGeno)$coefficients[2,1])
    rejPred15Percent_MMSNP_BetaswSE = c(summary(rejPred15PercentMMSNPMMSNP)$coefficients[2,1], summary(rejPred15PercentMMSNPRGeno)$coefficients[2,1])
    rejPred15Percent_BeagleSNP_BetaswSE = c(summary(rejPred15PercentBeagleSNPBeagleSNP)$coefficients[2,1], summary(rejPred15PercentBeagleSNPRGeno)$coefficients[2,1])
    #models with correction for R geno
    rejPred15Percent_IBSSNPRGenos_BetaswSE = c(summary(rejPred15PercentIBSSNPIBSSNPRGeno)$coefficients[2,1],summary(rejPred15PercentIBSSNPIBSSNPRGeno)$coefficients[3,1])
    rejPred15Percent_IncompSNPRGenos_BetaswSE = c(summary(rejPred15PercentIncompSNPIncompSNPRGeno)$coefficients[2,1],summary(rejPred15PercentIncompSNPIncompSNPRGeno)$coefficients[3,1])
    rejPred15Percent_MMSNPRGenos_BetaswSE = c(summary(rejPred15PercentMMSNPMMSNPRGeno)$coefficients[2,1],summary(rejPred15PercentMMSNPMMSNPRGeno)$coefficients[3,1])
    rejPred15Percent_BeagleSNPRGenos_BetaswSE = c(summary(rejPred15PercentBeagleSNPBeagleSNPRGeno)$coefficients[2,1],summary(rejPred15PercentBeagleSNPBeagleSNPRGeno)$coefficients[3,1])
    
    #30 AR
    rejPred30Percent_IBSSNP_BetaswSE = c(summary(rejPred30PercentIBSSNPIBSSNP)$coefficients[2,1], summary(rejPred30PercentIBSSNPRGeno)$coefficients[2,1])
    rejPred30Percent_IncompSNP_BetaswSE = c(summary(rejPred30PercentIncompSNPIncompSNP)$coefficients[2,1], summary(rejPred30PercentIncompSNPRGeno)$coefficients[2,1])
    rejPred30Percent_MMSNP_BetaswSE = c(summary(rejPred30PercentMMSNPMMSNP)$coefficients[2,1], summary(rejPred30PercentMMSNPRGeno)$coefficients[2,1])
    rejPred30Percent_BeagleSNP_BetaswSE = c(summary(rejPred30PercentBeagleSNPBeagleSNP)$coefficients[2,1], summary(rejPred30PercentBeagleSNPRGeno)$coefficients[2,1])
    #models with correction for R geno
    rejPred30Percent_IBSSNPRGenos_BetaswSE = c(summary(rejPred30PercentIBSSNPIBSSNPRGeno)$coefficients[2,1],summary(rejPred30PercentIBSSNPIBSSNPRGeno)$coefficients[3,1])
    rejPred30Percent_IncompSNPRGenos_BetaswSE = c(summary(rejPred30PercentIncompSNPIncompSNPRGeno)$coefficients[2,1],summary(rejPred30PercentIncompSNPIncompSNPRGeno)$coefficients[3,1])
    rejPred30Percent_MMSNPRGenos_BetaswSE = c(summary(rejPred30PercentMMSNPMMSNPRGeno)$coefficients[2,1],summary(rejPred30PercentMMSNPMMSNPRGeno)$coefficients[3,1])
    rejPred30Percent_BeagleSNPRGenos_BetaswSE = c(summary(rejPred30PercentBeagleSNPBeagleSNPRGeno)$coefficients[2,1],summary(rejPred30PercentBeagleSNPBeagleSNPRGeno)$coefficients[3,1])
    
    #stack Betas and SE values
    BetaswSEAll_15Percent = rbind(rejPred15Percent_IBSSNP_BetaswSE, rejPred15Percent_IncompSNP_BetaswSE, rejPred15Percent_MMSNP_BetaswSE, rejPred15Percent_BeagleSNP_BetaswSE)
    BetaswSEAll_Rgenos_15Percent = rbind(rejPred15Percent_IBSSNPRGenos_BetaswSE, rejPred15Percent_IncompSNPRGenos_BetaswSE, rejPred15Percent_MMSNPRGenos_BetaswSE, rejPred15Percent_BeagleSNPRGenos_BetaswSE)
    BetaswSEAll_30Percent = rbind(rejPred30Percent_IBSSNP_BetaswSE, rejPred30Percent_IncompSNP_BetaswSE, rejPred30Percent_MMSNP_BetaswSE, rejPred30Percent_BeagleSNP_BetaswSE)
    BetaswSEAll_Rgenos_30Percent = rbind(rejPred30Percent_IBSSNPRGenos_BetaswSE, rejPred30Percent_IncompSNPRGenos_BetaswSE, rejPred30Percent_MMSNPRGenos_BetaswSE, rejPred30Percent_BeagleSNPRGenos_BetaswSE)
    
    SS_BetaswSEAll_15Percent[[ii]] = BetaswSEAll_15Percent
    SS_BetaswSEAll_Rgenos_15Percent[[ii]] = BetaswSEAll_Rgenos_15Percent
    SS_BetaswSEAll_30Percent[[ii]] = BetaswSEAll_30Percent
    SS_BetaswSEAll_Rgenos_30Percent[[ii]] = BetaswSEAll_Rgenos_30Percent
    
    #need to pull p-values from summary tables
    rejPred15Percent_IBSSNP_Pvals = c(summary(rejPred15PercentIBSSNPIBSSNP)$coefficients[2,4], summary(rejPred15PercentIBSSNPRGeno)$coefficients[2,4], lhIBS15pval)
    rejPred15Percent_IncompSNP_Pvals = c(summary(rejPred15PercentIncompSNPIncompSNP)$coefficients[2,4], summary(rejPred15PercentIncompSNPRGeno)$coefficients[2,4], lhIncomp15pval)
    rejPred15Percent_MMSNP_Pvals = c(summary(rejPred15PercentMMSNPMMSNP)$coefficients[2,4], summary(rejPred15PercentMMSNPRGeno)$coefficients[2,4], lhMM15pval)
    rejPred15Percent_BeagleSNP_Pvals = c(summary(rejPred15PercentBeagleSNPBeagleSNP)$coefficients[2,4], summary(rejPred15PercentBeagleSNPRGeno)$coefficients[2,4], lhBeagleIBD15pval)
    
    rejPred30Percent_IBSSNP_Pvals = c(summary(rejPred30PercentIBSSNPIBSSNP)$coefficients[2,4], summary(rejPred30PercentIBSSNPRGeno)$coefficients[2,4], lhIBS30pval)
    rejPred30Percent_IncompSNP_Pvals = c(summary(rejPred30PercentIncompSNPIncompSNP)$coefficients[2,4], summary(rejPred30PercentIncompSNPRGeno)$coefficients[2,4], lhIncomp30pval)
    rejPred30Percent_MMSNP_Pvals = c(summary(rejPred30PercentMMSNPMMSNP)$coefficients[2,4], summary(rejPred30PercentMMSNPRGeno)$coefficients[2,4], lhMM30pval)
    rejPred30Percent_BeagleSNP_Pvals = c(summary(rejPred30PercentBeagleSNPBeagleSNP)$coefficients[2,4], summary(rejPred30PercentBeagleSNPRGeno)$coefficients[2,4], lhBeagleIBD30pval)
    
    #stack p-values
    PvalsAll_15Percent = rbind(rejPred15Percent_IBSSNP_Pvals, rejPred15Percent_IncompSNP_Pvals, rejPred15Percent_MMSNP_Pvals, rejPred15Percent_BeagleSNP_Pvals)
    PvalsAll_30Percent = rbind(rejPred30Percent_IBSSNP_Pvals, rejPred30Percent_IncompSNP_Pvals, rejPred30Percent_MMSNP_Pvals, rejPred30Percent_BeagleSNP_Pvals)
    
    SS_PvalsAll_15Percent[[ii]] = PvalsAll_15Percent #10 columns x 5 rows for each effect size
    SS_PvalsAll_30Percent[[ii]] = PvalsAll_30Percent #10 columns x 5 rows for each effect size
  }
  #unlist and make into matrix
  SS_BetaswSEAll_15Percent.mat = do.call(rbind, SS_BetaswSEAll_15Percent)
  SS_BetaswSEAll_Rgenos_15Percent.mat = do.call(rbind, SS_BetaswSEAll_Rgenos_15Percent)
  SS_BetaswSEAll_30Percent.mat = do.call(rbind, SS_BetaswSEAll_30Percent)
  SS_BetaswSEAll_Rgenos_30Percent.mat = do.call(rbind, SS_BetaswSEAll_Rgenos_30Percent)
  SS_PvalsAll_15Percent.mat = do.call(rbind, SS_PvalsAll_15Percent)
  SS_PvalsAll_30Percent.mat = do.call(rbind, SS_PvalsAll_30Percent)
  
  #output Betas and SE values
  write.csv(SS_BetaswSEAll_15Percent.mat, paste0("True_BetasWSE_SingleSNPs_15PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_BetaswSEAll_Rgenos_15Percent.mat, paste0("FittedWRGeno_BetasWSE_SingleSNPs_15PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_BetaswSEAll_30Percent.mat, paste0("True_BetasWSE_SingleSNPs_30PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_BetaswSEAll_Rgenos_30Percent.mat, paste0("FittedWRGeno_BetasWSE_SingleSNPs_30PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  
  #output p-values
  write.csv(SS_PvalsAll_15Percent.mat, paste0("PValues_ScoreAndJoint_SingleSNPs_15PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_PvalsAll_30Percent.mat, paste0("PValues_ScoreAndJoint_SingleSNPs_30PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  
}

