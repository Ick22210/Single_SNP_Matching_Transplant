#!/bin/bash

CHR=${1?Error: no chr given}
SS=${2?Error: no sample size given}
NSIMS=${3?Error: no number of simulations given}
RAND=${4?Error: no random number given}
GENE=${5?Error: no gene name given}
START=${6?Error: no starting value given}

export CHR
export SS
export NSIMS
export GENE
export START

for ((j=START;j<=NSIMS;j++));
do		
	export RAND	
	export j 

	#pull haplotypes from hapmap3 data using RAND
	#sh PullHaplotypesForSimulations.sh
	
	#increment random number for next sim
	#RAND=$((RAND + (2*SS)))
	
	#add headers to the created vcf files
	#sh AddVCFHeaders.sh
	
	#convert VCF to plink file format
	#sh ConvertVCFToPlink.sh	

	#calc binary mismatch score with R
	#sh CalcBinMM.sh
	
	#calc p-values for Power (in R)
	sh CalcPValues_PowerAnalysis_MMScores.sh
	
done