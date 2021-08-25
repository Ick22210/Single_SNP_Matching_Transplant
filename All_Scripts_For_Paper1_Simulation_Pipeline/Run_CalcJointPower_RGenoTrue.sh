#!/bin/bash
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
module load R/3.5.3
CHR=${1?Error: no chr given}
SS=${2?Error: no sample size given}
NSIMS=${3?Error: no number of simulations given}
GENE=${4?Error: no gene name given}

export CHR
export SS
export NSIMS
export GENE

Rscript CalcJointPower_RGenoTrue.R $CHR $SS $NSIMS $GENE --no-save