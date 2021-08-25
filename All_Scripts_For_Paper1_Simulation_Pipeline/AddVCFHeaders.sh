#!/bin/bash

path="/home/vlynn/Simulating_With_Haps/${GENE}_Results_${SS}Pairs"
cd $path

cat /home/vlynn/Simulating_With_Haps/All_Scripts_For_Paper1_Simulation_Pipeline/VCFheaderfiles/VCFheader_${SS}.txt Sim_${j}_Var_Chr${CHR}_${SS}Pairs_${GENE}.vcf > Sim_${j}_file3.txt; mv Sim_${j}_file3.txt Sim_${j}_Var_Chr${CHR}_${SS}Pairs_${GENE}.vcf

	