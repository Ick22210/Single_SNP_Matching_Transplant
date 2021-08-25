#!/bin/bash

path="/home/vlynn/Simulating_With_Haps/${GENE}_Results_${SS}Pairs"
cd $path

/home/vlynn/Simulating_With_Haps/All_Scripts_For_Paper1_Simulation_Pipeline/Plink_19/plink --vcf Sim_${j}_Var_Chr${CHR}_${SS}Pairs_${GENE}.vcf --make-bed --out Sim_${j}_Var_Chr${CHR}_${SS}Pairs_${GENE}
