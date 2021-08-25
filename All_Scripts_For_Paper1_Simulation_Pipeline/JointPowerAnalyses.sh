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
	
	#calc p-values for Power (in R)
	sh Run_JointlyTestingPowerScoreRGeno.sh
	
done