#!/bin/bash

path="/home/vlynn/Simulating_With_Haps/${GENE}_Results_${SS}Pairs"
cd $path 

/home/vlynn/plink_linux_x86_64_20190304/plink --bfile Sim_${j}_Var_Chr${CHR}_${SS}Pairs_${GENE} --recode beagle --out Sim_${j}_Var_Chr${CHR}_${SS}Pairs_${GENE}_Output_Beagle_Sample
