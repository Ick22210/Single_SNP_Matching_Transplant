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

#make list of names of ibd files
varibdFile = paste0("Sim_",simNum,"_Var_Chr",chr,"_",numSamples,"_BeagleIBDOutput_Sample.Sim_",simNum,"_Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_Output_Beagle_Sample.chr-",chr,".dat.ibd")

#read in Beagle IBD files
varbeagleIBD = read.table(varibdFile, skip = 2)

#need to pull what was column 2
varseqColII = seq(2,ncol(varbeagleIBD),4)
varcolII = varbeagleIBD[,varseqColII]

#need to save the items in the column for later
varbeagleIBDProbsMat = matrix(unlist(varcolII), ncol = ncol(varcolII), byrow = F)
write.csv(varbeagleIBDProbsMat, paste0("Var_Chr",chr,"_",numSamples,"Pairs_",gene,"_BeagleIBDProbabilities_Simulation",simNum,".csv"))

