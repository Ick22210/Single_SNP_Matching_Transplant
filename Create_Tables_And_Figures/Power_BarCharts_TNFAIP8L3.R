library(RColorBrewer)
library(ggplot2)

#effect_size_range = 0.22
#effect_size_range = 0.3
effect_size_range = 0.37

#TNFAIP8L3 - SS = 1000
#read in MAFs
minorAlleleFreqsTNFAIP8L3 = c(0.355, 0.214, 0.181, 0.392, 0.214, 0.104, 0.210, 0.086, 0.198, 0.352, 0.300, 0.142, 0.063, 0.088, 0.169, 0.143, 0.169, 0.291, 0.139, 0.288, 0.140, 0.066, 0.367, 0.291, 0.174, 0.469, 0.366, 0.154, 0.274, 0.458, 0.293, 0.067, 0.162, 0.053, 0.095, 0.474, 0.169, 0.313, 0.298, 0.410, 0.114, 0.108, 0.330, 0.074, 0.390, 0.362, 0.079, 0.076)

#set working directory to where data is
setwd("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Haplotypes/Gene_Files/Power_Haplotypes/TNFAIP8L3")

#read in files
RGenos_TNFAIP8L3_1000_15AR = read.csv(paste0("RGenotypes_Power_Chr15_1000Pairs_15PercentAR_TNFAIP8L3_ES_",effect_size_range,".csv"))

IBS_TNFAIP8L3_1000_15AR = read.csv(paste0("IBSScoresSS_Power_Chr15_1000Pairs_15PercentAR_TNFAIP8L3_ES_",effect_size_range,".csv"))

Incomp_TNFAIP8L3_1000_15AR = read.csv(paste0("IncompScoresSS_Power_Chr15_1000Pairs_15PercentAR_TNFAIP8L3_ES_",effect_size_range,".csv"))

MM_TNFAIP8L3_1000_15AR = read.csv(paste0("MMScoresSS_Power_Chr15_1000Pairs_15PercentAR_TNFAIP8L3_ES_",effect_size_range,".csv"))

BeagleIBD_TNFAIP8L3_1000_15AR = read.csv(paste0("BeagleIBDProbs_Power_Chr15_1000Pairs_15PercentAR_TNFAIP8L3_ES_",effect_size_range,".csv"))

#remove 1st column
RGenos_TNFAIP8L3_1000_15AR = RGenos_TNFAIP8L3_1000_15AR[,2:ncol(RGenos_TNFAIP8L3_1000_15AR)]
IBS_TNFAIP8L3_1000_15AR = IBS_TNFAIP8L3_1000_15AR[,2:ncol(IBS_TNFAIP8L3_1000_15AR)]
Incomp_TNFAIP8L3_1000_15AR = Incomp_TNFAIP8L3_1000_15AR[,2:ncol(Incomp_TNFAIP8L3_1000_15AR)]
MM_TNFAIP8L3_1000_15AR = MM_TNFAIP8L3_1000_15AR[,2:ncol(MM_TNFAIP8L3_1000_15AR)]
BeagleIBD_TNFAIP8L3_1000_15AR = BeagleIBD_TNFAIP8L3_1000_15AR[,2:ncol(BeagleIBD_TNFAIP8L3_1000_15AR)]

#add column names
colnames(RGenos_TNFAIP8L3_1000_15AR) = colnames(IBS_TNFAIP8L3_1000_15AR) = colnames(Incomp_TNFAIP8L3_1000_15AR) = colnames(MM_TNFAIP8L3_1000_15AR) = colnames(BeagleIBD_TNFAIP8L3_1000_15AR) = c("RGenos, AR 15%", "IBS, AR 15%", "Incomp, AR 15%",  "Mismatch, AR 15%","BeagleIBD, AR 15%")

#add column with model phenotype variable
ModelPheno_Column = c(rep(c("R Genotype"), nrow(RGenos_TNFAIP8L3_1000_15AR)), rep(c("IBS Mismatch Score"), nrow(RGenos_TNFAIP8L3_1000_15AR)), rep(c("Incompatibility Score"), nrow(RGenos_TNFAIP8L3_1000_15AR)), rep(c("Allogenomics Mismatch Score"), nrow(RGenos_TNFAIP8L3_1000_15AR)), rep(c("IBD Similarity Score"), nrow(RGenos_TNFAIP8L3_1000_15AR)))

#add MAFs for SNPs
minorAlleleFreqsTNFAIP8L3.rep = rep(minorAlleleFreqsTNFAIP8L3, length(effect_size_range))
minorAlleleFreqsTNFAIP8L3.rep.all = rep(minorAlleleFreqsTNFAIP8L3.rep, 5)

#bind everything together
All_TNFAIP8L3_1000_15AR = rbind(RGenos_TNFAIP8L3_1000_15AR, IBS_TNFAIP8L3_1000_15AR, Incomp_TNFAIP8L3_1000_15AR, MM_TNFAIP8L3_1000_15AR, BeagleIBD_TNFAIP8L3_1000_15AR)

All_TNFAIP8L3_1000_15AR_all = cbind(minorAlleleFreqsTNFAIP8L3.rep.all, All_TNFAIP8L3_1000_15AR, ModelPheno_Column)

#separate based on MAF
maf005 = which(All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all >= 0.05 & All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all < 0.056)

maf010 = which(All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all >= 0.09 & All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all < 0.112)

maf10subset = maf010[seq(1,length(maf010),3)]

maf020 = which(All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all >= 0.19 & All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all < 0.203)

maf030 = which(All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all >= 0.299 & All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all < 0.301)

maf040 = which(All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all >= 0.393 & All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all < 0.42)

maf050 = which(All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all >= 0.47 & All_TNFAIP8L3_1000_15AR_all$minorAlleleFreqsTNFAIP8L3.rep.all < 0.5)

allMAF = c(maf005, maf10subset, maf020, maf030, maf040, maf050)

All_TNFAIP8L3_1000_15AR_allSubset = All_TNFAIP8L3_1000_15AR_all[allMAF,]

TNFAIP8L3Colors = brewer.pal(5, "Spectral")

#plot for single snp and gene based
ggplot(All_TNFAIP8L3_1000_15AR_allSubset, aes(fill=ModelPheno_Column, 
        y=`RGenos, AR 15%`, x=minorAlleleFreqsTNFAIP8L3.rep.all)) + 
        geom_bar(position="dodge", color="black", stat="identity") + ggtitle("Estimated Power - True Pheno R Geno, 15% Outcome, True OR 1.45") +
        xlab("SNP MAF") + ylab("Estimated Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Model Covariate") + scale_fill_brewer(palette="Spectral") #+  geom_hline(yintercept=0.65, color = "red") +  geom_hline(yintercept=0.85, color = "blue")

ggplot(All_TNFAIP8L3_1000_15AR_allSubset, aes(fill=ModelPheno_Column, 
        y=`IBS, AR 15%`, x=minorAlleleFreqsTNFAIP8L3.rep.all)) + 
        geom_bar(position="dodge", color="black", stat="identity") + ggtitle("Estimated Power - True Pheno IBS Single SNP, 15% Outcome, True OR 1.45") +
        xlab("SNP MAF") + ylab("Estimated Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Model Covariate") + scale_fill_brewer(palette="Spectral") #+  geom_hline(yintercept=0.65, color = "red") +  geom_hline(yintercept=0.85, color = "blue")

ggplot(All_TNFAIP8L3_1000_15AR_allSubset, aes(fill=ModelPheno_Column, 
        y=`Incomp, AR 15%`, x=minorAlleleFreqsTNFAIP8L3.rep.all)) + 
        geom_bar(position="dodge", color="black", stat="identity") + ggtitle("Estimated Power - True Pheno Incomp. Single SNP, 15% AR, True OR 1.45") +
        xlab("SNP MAF") + ylab("Estimated Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Model Covariate") + scale_fill_brewer(palette="Spectral")# +  geom_hline(yintercept=0.65, color = "red") +  geom_hline(yintercept=0.85, color = "blue")

ggplot(All_TNFAIP8L3_1000_15AR_allSubset, aes(fill=ModelPheno_Column, 
        y=`Mismatch, AR 15%`, x=minorAlleleFreqsTNFAIP8L3.rep.all)) + 
        geom_bar(position="dodge", color="black", stat="identity") + ggtitle("Estimated Power - True Pheno Mismatch Single SNP, 15% AR, True OR 1.45") +
        xlab("SNP MAF") + ylab("Estimated Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Model Covariate") + scale_fill_brewer(palette="Spectral") #+  geom_hline(yintercept=0.65, color = "red") +  geom_hline(yintercept=0.85, color = "blue")

ggplot(All_TNFAIP8L3_1000_15AR_allSubset, aes(fill=ModelPheno_Column, 
        y=`BeagleIBD, AR 15%`, x=minorAlleleFreqsTNFAIP8L3.rep.all)) + 
        geom_bar(position="dodge", color="black", stat="identity") + ggtitle("Estimated Power - True Pheno IBD Probability Single SNP, 15% AR, True OR 1.25") +
        xlab("SNP MAF") + ylab("Estimated Power") + coord_cartesian(ylim = c(0, 1)) + labs(fill = "Model Covariate") + scale_fill_brewer(palette="Spectral") #+  geom_hline(yintercept=0.65, color = "red") +  geom_hline(yintercept=0.85, color = "blue")

