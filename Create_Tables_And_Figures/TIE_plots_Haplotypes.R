SS = 1000
GENE = "TNFAIP8L3"
CHR = 15

library(RColorBrewer)
######################################################################################################################################################
setwd(paste0("C:/Users/ickme/Box Sync/V_Arthur_Dissertation_Work/Haplotypes/Gene_Files/TIE_Haplotypes/",GENE,"/SS_",SS,"_NSims_5000"))

#dim is 3 x numSNPs
# DGenosTIE = read.csv(paste0("DGenotypes_TypeIError_AR15AndAR30_Uncorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
RGenosTIE = read.csv(paste0("RGenotypes_TypeIError_AR15AndAR30_Uncorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
IBSTIE = read.csv(paste0("IBSScoresSS_TypeIError_AR15AndAR30_Uncorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
IncompTIE = read.csv(paste0("IncompScoresSS_TypeIError_AR15AndAR30_Uncorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
MismatchTIE = read.csv(paste0("MismatchScoresSS_TypeIError_AR15AndAR30_Uncorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
BeagleIBDTIE = read.csv(paste0("BeagleIBDProbs_TypeIError_AR15AndAR30_Uncorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))

#remove first column
# DGenosTIE = DGenosTIE[,2:ncol(DGenosTIE)]
RGenosTIE = RGenosTIE[,2:ncol(RGenosTIE)]
IBSTIE = IBSTIE[,2:ncol(IBSTIE)]
IncompTIE = IncompTIE[,2:ncol(IncompTIE)]
MismatchTIE = MismatchTIE[,2:ncol(MismatchTIE)]
BeagleIBDTIE = BeagleIBDTIE[,2:ncol(BeagleIBDTIE)]

# DGenosTIEBonf = read.csv(paste0("DGenotypes_TypeIError_AR15AndAR30_BonfCorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
RGenosTIEBonf = read.csv(paste0("RGenotypes_TypeIError_AR15AndAR30_BonfCorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
IBSTIEBonf = read.csv(paste0("IBSScoresSS_TypeIError_AR15AndAR30_BonfCorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
IncompTIEBonf = read.csv(paste0("IncompScoresSS_TypeIError_AR15AndAR30_BonfCorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
MismatchTIEBonf = read.csv(paste0("MismatchScoresSS_TypeIError_AR15AndAR30_BonfCorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))
BeagleIBDTIEBonf = read.csv(paste0("BeagleIBDProbs_TypeIError_AR15AndAR30_BonfCorrected_Chr",CHR,"_",SS,"Pairs",GENE,".csv"))

#remove first column
# DGenosTIEBonf = DGenosTIEBonf[,2:ncol(DGenosTIEBonf)]
RGenosTIEBonf = RGenosTIEBonf[,2:ncol(RGenosTIEBonf)]
IBSTIEBonf = IBSTIEBonf[,2:ncol(IBSTIEBonf)]
IncompTIEBonf = IncompTIEBonf[,2:ncol(IncompTIEBonf)]
MismatchTIEBof = MismatchTIEBonf[,2:ncol(MismatchTIEBonf)]
BeagleIBDTIEBonf = BeagleIBDTIEBonf[,2:ncol(BeagleIBDTIEBonf)]

#################################
###Plot MAF vs TIE?
minorAlleleFreqsTNFAIP8L3 = c(0.355, 0.214, 0.181, 0.392, 0.214, 0.104, 0.210, 0.086, 0.198, 0.352, 0.300, 0.142, 0.063, 0.088, 0.169, 0.143, 0.169,
                              0.291, 0.139, 0.288, 0.140, 0.066, 0.367, 0.291, 0.174, 0.469, 0.366, 0.154, 0.274, 0.458, 0.293, 0.067, 0.162, 0.053, 0.095,
                              0.474, 0.169, 0.313, 0.298, 0.410, 0.114, 0.108, 0.330, 0.074, 0.390, 0.362, 0.079, 0.076)
######################################################################################################################
#subset to fewer MAFs
MAF = minorAlleleFreqsTNFAIP8L3

MAFsubset = MAF[match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]

# DGenosTIEsubset = DGenosTIE[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
RGenosTIEsubset = RGenosTIE[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
IBSTIEsubset = IBSTIE[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
IncompTIEsubset = IncompTIE[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
MismatchTIEsubset = MismatchTIE[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
BeagleIBDTIEsubset = BeagleIBDTIE[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]

# DGenosTIEBonfsubset = DGenosTIEBonf[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
RGenosTIEBonfsubset = RGenosTIEBonf[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
IBSTIEBonfsubset = IBSTIEBonf[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
IncompTIEBonfsubset = IncompTIEBonf[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
MismatchTIEBonfsubset = MismatchTIEBof[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]
BeagleIBDTIEBonfsubset = BeagleIBDTIEBonf[,match(unique(round(MAF,digits = 2)), round(MAF,digits = 2))]

########################################
numSNPs = length(MAF)
bonfLine = 0.05/numSNPs
evenseq = seq(2,12,2)
GeneColors = brewer.pal(12,"Paired")
GeneColors = GeneColors[evenseq]

########################################
#plot MAF vs TIE for each SNP (uncorrected)
par(mfrow=c(2,2))

#Uncorrected
# plot(MAFsubset,DGenosTIEsubset[1,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[1],pch=15,xlab="MAF",ylab = "Estimated Type I Error", main="Type I Error - 15% Incidence")
plot(MAFsubset,RGenosTIEsubset[1,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[2],pch=17,xlab="MAF",ylab = "Estimated Type I Error", main="Type I Error - 15% Incidence")
points(MAFsubset,IBSTIEsubset[1,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[3],pch=18)
points(MAFsubset,IncompTIEsubset[1,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[4],pch=19)
points(MAFsubset,MismatchTIEsubset[1,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[5],pch=12)
points(MAFsubset,BeagleIBDTIEsubset[1,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[6],pch=8)
abline(h=0.05)

#Uncorrected
# plot(MAFsubset,DGenosTIEsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[1],pch=15,xlab="MAF",ylab = "Estimated Type I Error", main="Type I Error - 30% Incidence")
plot(MAFsubset,RGenosTIEsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[2],pch=17,xlab="MAF",ylab = "Estimated Type I Error", main="Type I Error - 30% Incidence")
points(MAFsubset,IBSTIEsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[3],pch=18)
points(MAFsubset,IncompTIEsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[4],pch=19)
points(MAFsubset,MismatchTIEsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[5],pch=12)
points(MAFsubset,BeagleIBDTIEsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.1),col=GeneColors[6],pch=8)
abline(h=0.05)

################################################
#plot MAF vs TIE for each SNP (Bonf corrected)

#All methods - Bonf - 15% incidence
# plot(MAFsubset,DGenosTIEBonfsubset[1,],xlim = c(0.025,0.5), ylim = c(0,0.002), col=GeneColors[1],pch=15,xlab="MAF",ylab = "Est Type I Error - Bonf. Corr.", main="Type I Error - 15% Incidence")
plot(MAFsubset,RGenosTIEBonfsubset[1,],xlim = c(0.025,0.5), ylim = c(0,0.002),col=GeneColors[2],pch=17,xlab="MAF",ylab = "Est Type I Error - Bonf. Corr.", main="Type I Error - 15% Incidence")
points(MAFsubset,IBSTIEBonfsubset[1,],xlim = c(0.025,0.5), ylim = c(0,0.002),col=GeneColors[3],pch=18)
points(MAFsubset,IncompTIEBonfsubset[1,],xlim = c(0.025,0.5), ylim = c(0,0.002),col=GeneColors[4],pch=19)
points(MAFsubset,MismatchTIEBonfsubset[1,],xlim = c(0.025,0.5),ylim = c(0,0.002),col=GeneColors[5],pch=12)
points(MAFsubset,BeagleIBDTIEBonfsubset[1,],xlim = c(0.025,0.5), ylim = c(0,0.002),col=GeneColors[6],pch=8)
abline(h=bonfLine)

#All methods - Bonf - 30% incidence
# plot(MAFsubset,DGenosTIEBonfsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.002),col=GeneColors[1],pch=15,xlab="MAF",ylab = "Est Type I Error - Bonf. Corr.", main="Type I Error - 30% Incidence")
plot(MAFsubset,RGenosTIEBonfsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.002),col=GeneColors[2],pch=17,xlab="MAF",ylab = "Est Type I Error - Bonf. Corr.", main="Type I Error - 30% Incidence")
points(MAFsubset,IBSTIEBonfsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.002),col=GeneColors[3],pch=18)
points(MAFsubset,IncompTIEBonfsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.002),col=GeneColors[4],pch=19)
points(MAFsubset,MismatchTIEBonfsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.002),col=GeneColors[5],pch=12)
points(MAFsubset,BeagleIBDTIEBonfsubset[2,],xlim = c(0.025,0.5),ylim = c(0,0.002),col=GeneColors[6],pch=8)
abline(h=bonfLine)

par(mfrow=c(1,1))

#plot legend
plot(MAFsubset,RGenosTIEBonfsubset[1,],xlim = c(0.05,0.5),ylim = c(0,0.15), type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",xpd=FALSE,legend=c("R Genotype", "IBS Mismatch Score", "Incompatibility Score", "Allogenomics Mismatch Score", "IBD Similarity Score"),  col=GeneColors[2:6], pch=c(17,18,19,12,8), bty="n")
