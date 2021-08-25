#!/bin/bash

path="/home/vlynn/Simulating_With_Haps/${GENE}_Results_${SS}Pairs"
cd $path 

/home/vlynn/jdk-12.0.1/bin/java -jar /home/vlynn/beagle.jar \
unphased=Sim_${j}_Var_Chr${CHR}_${SS}Pairs_${GENE}_Output_Beagle_Sample.chr-${CHR}.dat \
missing=0 \
markers=Sim_${j}_Var_Chr${CHR}_${SS}Pairs_${GENE}_Output_Beagle_Sample.chr-${CHR}.map \
ibdpairs=/home/vlynn/Simulating_With_Haps/Related_Pairs_Files/RelatedPairings_${SS}.txt \
out=Sim_${j}_Var_Chr${CHR}_${SS}_BeagleIBDOutput_Sample
