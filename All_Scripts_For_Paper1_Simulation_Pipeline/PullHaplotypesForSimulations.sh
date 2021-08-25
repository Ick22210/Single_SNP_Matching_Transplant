#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh 
fi
   
module load R/3.6.1
Rscript SamplingHaplotypesFromHapMap3.R $CHR $j $SS $RAND $GENE --nosave
