##########################################
### Calculate p-values (Power Analysis)
### Jointly testing with Score/RGeno
### R Geno is true model
###########################################

#load epicalc library
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
    #True phenotype is generated with R geno
    # true model (only using R geno to fit)
    rejPred15PercentRGenoRGeno = glm(P_RGeno_15incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoRGeno = glm(P_RGeno_30incidence[,ii] ~ R_Genos.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    ## correct for IBS Score
    rejPred15PercentRGenoIBSSNP = glm(P_RGeno_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoIBSSNP = glm(P_RGeno_30incidence[,ii] ~ IBS_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    ## correct for Incomp Score
    rejPred15PercentRGenoIncompSNP = glm(P_RGeno_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoIncompSNP = glm(P_RGeno_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    ## correct for MM Score
    rejPred15PercentRGenoMMSNP = glm(P_RGeno_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoMMSNP = glm(P_RGeno_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    ## correct for Beagle Score
    rejPred15PercentRGenoBeagleSNP = glm(P_RGeno_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoBeagleSNP = glm(P_RGeno_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii], family=binomial, data = P_RGeno_30incidence)
    
    ###Joint Model
    ##Code for model that corrects for IBS Score and R Genotype
    rejPred15PercentRGenoIBSSNPRGeno = glm(P_RGeno_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii] + R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoIBSSNPRGeno = glm(P_RGeno_15incidence[,ii] ~ IBS_SingleSNPs.mat[,ii] + R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    ##Jointly test significance of IBS Score and R Genos
    ##make null models
    rejPred15PercentRGenoIBSSNPRGenoNull = glm(P_RGeno_15incidence[,ii] ~ 1, family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoIBSSNPRGenoNull = glm(P_RGeno_15incidence[,ii] ~ 1, family=binomial, data = P_RGeno_15incidence)
    ##tells to test against the null model (i.e., only the intercept is not set to 0)
    lhIBS15 = lrtest(rejPred15PercentRGenoIBSSNPRGeno, rejPred15PercentRGenoIBSSNPRGenoNull)
    lhIBS30 = lrtest(rejPred30PercentRGenoIBSSNPRGeno, rejPred30PercentRGenoIBSSNPRGenoNull)
    #pull the p-values
    lhIBS15pval = lhIBS15$p.value
    lhIBS30pval = lhIBS30$p.value
    
    ##Code for model that corrects for Incomp Score and R Genotype
    rejPred15PercentRGenoIncompSNPRGeno = glm(P_RGeno_15incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii]+ R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoIncompSNPRGeno = glm(P_RGeno_30incidence[,ii] ~ Incomp_SingleSNPs.mat[,ii]+ R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    ##Jointly test significance of Incomp Score and R Genos
    ##make null models
    rejPred15PercentRGenoIncompSNPRGenoNull = glm(P_RGeno_15incidence[,ii] ~ 1, family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoIncompSNPRGenoNull = glm(P_RGeno_30incidence[,ii] ~ 1, family=binomial, data = P_RGeno_15incidence)
    ##tells to test against the null model (i.e., only the intercept is not set to 0)
    lhIncomp15 = lrtest(rejPred15PercentRGenoIncompSNPRGeno, rejPred15PercentRGenoIncompSNPRGenoNull)
    lhIncomp30 = lrtest(rejPred30PercentRGenoIncompSNPRGeno, rejPred30PercentRGenoIncompSNPRGenoNull)
    #pull the p-values
    lhIncomp15pval = lhIncomp15$p.value
    lhIncomp30pval = lhIncomp30$p.value
    
    ##Code for model that also corrects for R Genotype
    rejPred15PercentRGenoMMSNPRGeno = glm(P_RGeno_15incidence[,ii] ~ MM_SingleSNPs.mat[,ii] + R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoMMSNPRGeno = glm(P_RGeno_30incidence[,ii] ~ MM_SingleSNPs.mat[,ii] + R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    ##Jointly test significance of IBS Score and R Genos
    ##make null models
    rejPred15PercentRGenoMMSNPRGenoNull = glm(P_RGeno_15incidence[,ii] ~ 1, family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoMMSNPRGenoNull = glm(P_RGeno_30incidence[,ii] ~ 1, family=binomial, data = P_RGeno_15incidence)
    ##tells to test against the null model (i.e., only the intercept is not set to 0)
    lhMM15 = lrtest(rejPred15PercentRGenoMMSNPRGeno, rejPred15PercentRGenoMMSNPRGenoNull)
    lhMM30 = lrtest(rejPred30PercentRGenoMMSNPRGeno, rejPred30PercentRGenoMMSNPRGenoNull)
    #pull the p-values
    lhMM15pval = lhMM15$p.value
    lhMM30pval = lhMM30$p.value
    
    ##Code for model that also corrects for R Genotype
    rejPred15PercentRGenoBeagleSNPRGeno = glm(P_RGeno_15incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii] + R_Genos.mat[,ii], family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoBeagleSNPRGeno = glm(P_RGeno_30incidence[,ii] ~ BeagleIBD_SingleSNPs.mat.norm[,ii] + R_Genos.mat[,ii], family=binomial, data = P_RGeno_30incidence)
    ##Jointly test significance of IBS Score and R Genos
    ##make null models
    rejPred15PercentRGenoBeagleSNPRGenoNull = glm(P_RGeno_15incidence[,ii] ~ 1, family=binomial, data = P_RGeno_15incidence)
    rejPred30PercentRGenoBeagleSNPRGenoNull = glm(P_RGeno_30incidence[,ii] ~ 1, family=binomial, data = P_RGeno_30incidence)
    ##tells to test against the null model (i.e., only the intercept is not set to 0)
    lhBeagleIBD15 = lrtest(rejPred15PercentRGenoBeagleSNPRGeno, rejPred15PercentRGenoBeagleSNPRGenoNull)
    lhBeagleIBD30 = lrtest(rejPred30PercentRGenoBeagleSNPRGeno, rejPred30PercentRGenoBeagleSNPRGenoNull)
    #pull the p-values
    lhBeagleIBD15pval = lhBeagleIBD15$p.value
    lhBeagleIBD30pval = lhBeagleIBD30$p.value
    
    #need to pull estimates for later
    #correct model
    rejPred15Percent_RGenoIBSSNP_BetaswSE = c(summary(rejPred15PercentRGenoIBSSNP)$coefficients[2,1], summary(rejPred15PercentRGenoRGeno)$coefficients[2,1])
    rejPred15Percent_RGenoIncompSNP_BetaswSE = c(summary(rejPred15PercentRGenoIncompSNP)$coefficients[2,1], summary(rejPred15PercentRGenoRGeno)$coefficients[2,1])
    rejPred15Percent_RGenoMMSNP_BetaswSE = c(summary(rejPred15PercentRGenoMMSNP)$coefficients[2,1], summary(rejPred15PercentRGenoRGeno)$coefficients[2,1])
    rejPred15Percent_RGenoBeagleSNP_BetaswSE = c(summary(rejPred15PercentRGenoBeagleSNP)$coefficients[2,1], summary(rejPred15PercentRGenoRGeno)$coefficients[2,1])
    #models with correction for R geno
    rejPred15Percent_IBSSNPRGenos_BetaswSE = c(summary(rejPred15PercentRGenoIBSSNPRGeno)$coefficients[2,1],summary(rejPred15PercentRGenoIBSSNPRGeno)$coefficients[3,1])
    rejPred15Percent_IncompSNPRGenos_BetaswSE = c(summary(rejPred15PercentRGenoIncompSNPRGeno)$coefficients[2,1],summary(rejPred15PercentRGenoIncompSNPRGeno)$coefficients[3,1])
    rejPred15Percent_MMSNPRGenos_BetaswSE = c(summary(rejPred15PercentRGenoMMSNPRGeno)$coefficients[2,1],summary(rejPred15PercentRGenoMMSNPRGeno)$coefficients[3,1])
    rejPred15Percent_BeagleSNPRGenos_BetaswSE = c(summary(rejPred15PercentRGenoBeagleSNPRGeno)$coefficients[2,1],summary(rejPred15PercentRGenoBeagleSNPRGeno)$coefficients[3,1])
    
    #30 AR
    rejPred30Percent_IBSSNP_BetaswSE = c(summary(rejPred30PercentRGenoIBSSNP)$coefficients[2,1], summary(rejPred30PercentRGenoRGeno)$coefficients[2,1])
    rejPred30Percent_IncompSNP_BetaswSE = c(summary(rejPred30PercentRGenoIncompSNP)$coefficients[2,1], summary(rejPred30PercentRGenoRGeno)$coefficients[2,1])
    rejPred30Percent_MMSNP_BetaswSE = c(summary(rejPred30PercentRGenoMMSNP)$coefficients[2,1], summary(rejPred30PercentRGenoRGeno)$coefficients[2,1])
    rejPred30Percent_BeagleSNP_BetaswSE = c(summary(rejPred30PercentRGenoBeagleSNP)$coefficients[2,1], summary(rejPred30PercentRGenoRGeno)$coefficients[2,1])
    #models with correction for R geno
    rejPred30Percent_IBSSNPRGenos_BetaswSE = c(summary(rejPred30PercentRGenoIBSSNPRGeno)$coefficients[2,1],summary(rejPred30PercentRGenoIBSSNPRGeno)$coefficients[3,1])
    rejPred30Percent_IncompSNPRGenos_BetaswSE = c(summary(rejPred30PercentRGenoIncompSNPRGeno)$coefficients[2,1],summary(rejPred30PercentRGenoIncompSNPRGeno)$coefficients[3,1])
    rejPred30Percent_MMSNPRGenos_BetaswSE = c(summary(rejPred30PercentRGenoMMSNPRGeno)$coefficients[2,1],summary(rejPred30PercentRGenoMMSNPRGeno)$coefficients[3,1])
    rejPred30Percent_BeagleSNPRGenos_BetaswSE = c(summary(rejPred30PercentRGenoBeagleSNPRGeno)$coefficients[2,1],summary(rejPred30PercentRGenoBeagleSNPRGeno)$coefficients[3,1])
    
    #stack Betas and SE values
    BetaswSEAll_15Percent = rbind(rejPred15Percent_RGenoIBSSNP_BetaswSE, rejPred15Percent_RGenoIncompSNP_BetaswSE, rejPred15Percent_RGenoMMSNP_BetaswSE, rejPred15Percent_RGenoBeagleSNP_BetaswSE)
    BetaswSEAll_Rgenos_15Percent = rbind(rejPred15Percent_IBSSNPRGenos_BetaswSE, rejPred15Percent_IncompSNPRGenos_BetaswSE, rejPred15Percent_MMSNPRGenos_BetaswSE, rejPred15Percent_BeagleSNPRGenos_BetaswSE)
    BetaswSEAll_30Percent = rbind(rejPred30Percent_IBSSNP_BetaswSE, rejPred30Percent_IncompSNP_BetaswSE, rejPred30Percent_MMSNP_BetaswSE, rejPred30Percent_BeagleSNP_BetaswSE)
    BetaswSEAll_Rgenos_30Percent = rbind(rejPred30Percent_IBSSNPRGenos_BetaswSE, rejPred30Percent_IncompSNPRGenos_BetaswSE, rejPred30Percent_MMSNPRGenos_BetaswSE, rejPred30Percent_BeagleSNPRGenos_BetaswSE)
    
    SS_BetaswSEAll_15Percent[[ii]] = BetaswSEAll_15Percent
    SS_BetaswSEAll_Rgenos_15Percent[[ii]] = BetaswSEAll_Rgenos_15Percent
    SS_BetaswSEAll_30Percent[[ii]] = BetaswSEAll_30Percent
    SS_BetaswSEAll_Rgenos_30Percent[[ii]] = BetaswSEAll_Rgenos_30Percent
    
    #need to pull p-values from summary tables
    rejPred15Percent_IBSSNP_Pvals = c(summary(rejPred15PercentRGenoIBSSNP)$coefficients[2,4], summary(rejPred15PercentRGenoRGeno)$coefficients[2,4], lhIBS15pval)
    rejPred15Percent_IncompSNP_Pvals = c(summary(rejPred15PercentRGenoIncompSNP)$coefficients[2,4], summary(rejPred15PercentRGenoRGeno)$coefficients[2,4], lhIncomp15pval)
    rejPred15Percent_MMSNP_Pvals = c(summary(rejPred15PercentRGenoMMSNP)$coefficients[2,4], summary(rejPred15PercentRGenoRGeno)$coefficients[2,4], lhMM15pval)
    rejPred15Percent_BeagleSNP_Pvals = c(summary(rejPred15PercentRGenoBeagleSNP)$coefficients[2,4], summary(rejPred15PercentRGenoRGeno)$coefficients[2,4], lhBeagleIBD15pval)
    
    rejPred30Percent_IBSSNP_Pvals = c(summary(rejPred30PercentRGenoIBSSNP)$coefficients[2,4], summary(rejPred30PercentRGenoRGeno)$coefficients[2,4], lhIBS30pval)
    rejPred30Percent_IncompSNP_Pvals = c(summary(rejPred30PercentRGenoIncompSNP)$coefficients[2,4], summary(rejPred30PercentRGenoRGeno)$coefficients[2,4], lhIncomp30pval)
    rejPred30Percent_MMSNP_Pvals = c(summary(rejPred30PercentRGenoMMSNP)$coefficients[2,4], summary(rejPred30PercentRGenoRGeno)$coefficients[2,4], lhMM30pval)
    rejPred30Percent_BeagleSNP_Pvals = c(summary(rejPred30PercentRGenoBeagleSNP)$coefficients[2,4], summary(rejPred30PercentRGenoRGeno)$coefficients[2,4], lhBeagleIBD30pval)
    
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
  write.csv(SS_BetaswSEAll_15Percent.mat, paste0("RGenoTrue_BetasWSE_SingleSNPs_15PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_BetaswSEAll_Rgenos_15Percent.mat, paste0("FittedWScore_BetasWSE_SingleSNPs_15PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_BetaswSEAll_30Percent.mat, paste0("RGenoTrue_BetasWSE_SingleSNPs_30PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_BetaswSEAll_Rgenos_30Percent.mat, paste0("FittedWScore_BetasWSE_SingleSNPs_30PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  
  #output p-values
  write.csv(SS_PvalsAll_15Percent.mat, paste0("PValues_ScoreAndJoint_RGeno_SingleSNPs_15PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  write.csv(SS_PvalsAll_30Percent.mat, paste0("PValues_ScoreAndJoint_RGeno_SingleSNPs_30PercentAR_EffectSize_",effect_size_range[aa],"_Simulation",simNum,".csv"))
  
}
