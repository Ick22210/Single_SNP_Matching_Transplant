#read in files
Rgenos = read.csv("RGenos_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_0.37.csv")
MM = read.csv("MM_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_0.37.csv")
Incomp = read.csv("Incomp_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_0.37.csv")
IBS = read.csv("IBS_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_0.37.csv")
BeagleIBD = read.csv("BeagleIBD_CIs_SingleSNPs_15PercentAR_BinaryPheno_EffectSize_0.37.csv")

colnames(Rgenos) = colnames(MM) = colnames(Incomp) = colnames(IBS) = colnames(BeagleIBD) = c("Name", "MAF", "R Geno", "IBS Mismatch", "Incomp", "AMS", "IBD Similarity")

RGenos_OR = exp(Rgenos[,3:ncol(Rgenos)])
MM_OR = exp(MM[,3:ncol(MM)])
Incomp_OR = exp(Incomp[,3:ncol(Incomp)])
IBS_OR = exp(IBS[,3:ncol(IBS)])
BeagleIBD_OR = exp(BeagleIBD[,3:ncol(BeagleIBD)])

#pull the averages only
seq = seq(1,nrow(Rgenos),3)

RGenos_OR_Ave = RGenos_OR[seq,]
MM_OR_Ave = MM_OR[seq,]
Incomp_OR_Ave = Incomp_OR[seq,] 
IBS_OR_Ave = IBS_OR[seq,]
BeagleIBD_OR_Ave = BeagleIBD_OR[seq,]

MAF = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)

RGenos_OR_Ave_MAF = cbind(MAF,RGenos_OR_Ave)
MM_OR_Ave_MAF = cbind(MAF,MM_OR_Ave)
Incomp_OR_Ave_MAF = cbind(MAF,Incomp_OR_Ave)
IBS_OR_Ave_MAF = cbind(MAF,IBS_OR_Ave)
BeagleIBD_OR_Ave_MAF = cbind(MAF,BeagleIBD_OR_Ave)

#bind all and export
All = rbind(RGenos_OR_Ave_MAF, IBS_OR_Ave_MAF, BeagleIBD_OR_Ave_MAF, Incomp_OR_Ave_MAF, MM_OR_Ave_MAF)
write.csv(All, "All_ORs_ByMAF.csv")
