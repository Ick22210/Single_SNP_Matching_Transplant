#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
module load R/3.5.3
Rscript CalcBinMM.R $CHR $SS $j $GENE --no-save
