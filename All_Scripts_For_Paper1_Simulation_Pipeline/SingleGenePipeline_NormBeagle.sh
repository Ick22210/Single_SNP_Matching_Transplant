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

#create directory to hold files that should not be erased
mkdir /home/vlynn/Simulating_With_Haps/${GENE}_Results_${SS}Pairs

for ((j=START;j<=NSIMS;j++));
do		
	export RAND	
	export j 
	
	#pull haplotypes from hapmap3 data using RAND
	sh PullHaplotypesForSimulations.sh
	
	#increment random number for next sim
	RAND=$((RAND + (2*SS)))
	
	#add headers to the created vcf files
	sh AddVCFHeaders.sh
	
	#convert VCF to plink file format
	sh ConvertVCFToPlink.sh
	
	#convert plink to beagle
	sh ConvertPlinkToBeagle.sh
	
	#run beagle ibd to get IBD probabilities
	sh RunBeagleIBD.sh &

	#calc ibs scores based on PLINK files with R
	sh CalcIBSAndIncompScoresVariableOnly.sh &
	
	wait
	
	#process the Beagle Output with R
	sh ReformatBeagleIBDScores.sh
	
	#calc p-values for type I error (in R)
	sh CalcPValues_TypeIErr.sh 
	
	#calc p-values for Power (in R)
	sh CalcPValues_PowerAnalysis.sh
	
	# clean up files
	find /home/vlynn/Simulating_With_Haps/${GENE}_Results_${SS}Pairs/ -type f -not -name '*.csv' -print0 | xargs -0 rm --
	
done
