############################################
## Sampling Haplotypes from HapMap3 Data ###
############################################

#to obtain the arguments from the bash files (chr, ss)
args <- commandArgs()

chr <- args[6]
i <- args[7]
numPairs <- args[8]
RAND <- args[9]
Gene <- args[10]

numPairs = as.integer(numPairs)
RAND = as.integer(RAND)

#setwd to location of plink data
path = paste0("/home/vlynn/Simulating_With_Haps/Gene_Files")
setwd(path)

#read in marker information from Haploview
markerInfo=read.table(paste0("All_Marker_",Gene,".txt"), header = F) 

#remove SNPs with MAF less than 0.05 and with percent genotyped over 75 (less than 25% missing)
#and with HW p values less than 0.001
snpSubset=which((markerInfo[,7]>=75)&(markerInfo[,10]>=0.05)&(markerInfo[,6]>0.001))  

#subset the data for the VCF file
markerSubset = markerInfo[snpSubset,]

#read in the haplotype frequency file
hapFreqFile=read.table(paste0("HapFreqAllMarkers_",Gene,"_spaced.txt"), header = F)

###Minor allele frequencies
allMinorAlleles=substring(markerInfo[,11],3,3) #this lists the minor alleles
#these are the minor alleles as A,G,T,C values for VCF files
minorAlleleSubset.letters = allMinorAlleles[snpSubset]
allMinorAlleles[allMinorAlleles=="A"]=1 #recodes the MAs as numbers, so that they match the haplotype freq file
allMinorAlleles[allMinorAlleles=="C"]=2
allMinorAlleles[allMinorAlleles=="G"]=3
allMinorAlleles[allMinorAlleles=="T"]=4
minorAlleleSubset=as.numeric(allMinorAlleles)[snpSubset] #this lists the minor alleles as numbers for all SNPs that pass filtering

#this lists the major alleles
allMajorAlleles = substring(markerInfo[,11],1,1) 
#this lists the major alleles as numbers for all SNPs that pass filtering
majorAlleleSubset.letters = allMajorAlleles[snpSubset] 

#vectorizes the haplotype frequencies and combines them
haplotypeFreq.vec=as.vector(t(hapFreqFile[,1:(ncol(hapFreqFile)-3)]))

#checks if haplotypeFreq.vec values are equal to MA, basically like, yes minor allele = 1, no = 0
#dimension is number of haplotypes possible x number of SNPs in haplotype
ma.haplotypeFreq.vec = matrix(as.numeric(haplotypeFreq.vec==rep(minorAlleleSubset,nrow(hapFreqFile))),nrow=nrow(hapFreqFile),ncol=length(minorAlleleSubset),byrow=T)

#divides haplotype frequencies by total sum (which is slightly less than 1) to get accurate freq values
correctedHapFreqs = hapFreqFile[,ncol(hapFreqFile)-1] /sum(hapFreqFile[,ncol(hapFreqFile)-1])

#not sure what KK stands for...
numDonors = numRecip = numPairs
KK = numDonors + numRecip

#should put in a random seed here to make sure sample differs each time
set.seed(RAND)

#pulling indices of the haplotype frequencies based on prob of haplotype, pulls KK of them
idHap1 = sample(1:length(correctedHapFreqs),KK,replace = T,prob = correctedHapFreqs) 
#pulling indices of the haplotype frequencies based on prob of haplotype, pulls KK of them
idHap2 = sample(1:length(correctedHapFreqs),KK,replace = T,prob = correctedHapFreqs) 
#this is adding the 2 haplotypes to make full genotype (I think I use this?)
#file is length of haplotype (number of SNPs) x KK in dimension
sampledGenotype = ma.haplotypeFreq.vec[idHap1, ] + ma.haplotypeFreq.vec[idHap2, ] 

############################################################################################
##In order to make Beagle Files, need to make these genotypes into VCF files
#Then VCF -> plink -> beagle format

#create vcf template that will get filled in by the needed information, but contains all the regular info
vcfTemplate = matrix(NA, nrow=length(minorAlleleSubset), ncol=(9 + KK))
#nrow is number of SNPs, ncols is  9 + number of subjects generated

#first col should be chrom number
vcfTemplate[,1] = chr

#second column is position, need to fill in later

#third column is id, which is just a . for msprime vcfs, 6 (quality) and 8 (info) also have .s
vcfTemplate[,3] = "."
vcfTemplate[,6] = "."
vcfTemplate[,8] = "."

#4th column is reference allele, and 5th is alternative allele, msprime vcf has A and T
vcfTemplate[,4] = majorAlleleSubset.letters
vcfTemplate[,5] = minorAlleleSubset.letters

#7th column is filter, msprime vcf has PASS
vcfTemplate[,7] = "PASS"

#9th column is format, msprime vcf has GT for genotype
vcfTemplate[,9] = "GT"

#add the vcf template to the final list
finalVCFsVar = vcfTemplate

###################################################
#now work on adding the positions to the vcfs
finalVCFsVar[,2] = markerSubset$V3
###################################################
#adding genotypes to the final vcfs

#make new lists to put output in
hapsSubsetVarMerged = matrix(paste(t(ma.haplotypeFreq.vec[idHap1,]), t(ma.haplotypeFreq.vec[idHap2,]), sep="|"), nrow=length(minorAlleleSubset), ncol=(KK))

#add to the final vcfs file
finalVCFsVar[,10:(KK+9)] = hapsSubsetVarMerged

##########################################################
#now need to export VCF file, sans comments, as separate files
write.table(finalVCFsVar, file=paste0("/home/vlynn/Simulating_With_Haps/",Gene,"_Results_",numPairs,"Pairs/Sim_",i,"_Var_Chr",chr,"_",numDonors,"Pairs_",Gene,".vcf"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
