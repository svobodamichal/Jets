#!/bin/bash

#NOT USED

#source ../set_paths.sh
ANALYSISDIR="/gpfs01/star/pwg/svomich/Jets/"
echo "ADir: $ANALYSISDIR"
STARLIB_VER=pro
starver $STARLIB_VER

TRG="MB"
export RMTYPE="BG_sp" #deltapT distribution: BG_sp|inplane|outplane
export V2CORR=0 # correct delta pT for event plane bias
#export CENTRAL=0
export CENTRAL=1

if [ $V2CORR -eq 1 ]; then
	SUFF="_v2"
else
	SUFF="_normal"
fi

PTLEADCUTS="2"

if [ $CENTRAL -eq 1 ]; then
#CENTSUFF="_central"
CENTSUFF="central"
#PTLEADCUTS="5 6 7"
else
#CENTSUFF="_peripheral"
CENTSUFF="peripheral"
#PTLEADCUTS="4 5 6 7"
fi

#export V2PATH="$ANALYSISDIR/macros/EP_corrections" #path to correction histograms
export V2PATH="$ANALYSISDIR/StRoot/macros/EP_corrections" #path to correction histograms
#export PATH_TO_DELTA_PT_HISTOGRAMS="$ANALYSISDIR/out/$TRG/embedding${CENTSUFF}_main"
export PATH_TO_DELTA_PT_HISTOGRAMS="$ANALYSISDIR/responseM/$CENTSUFF/"

#mkdir -p $PATH_TO_DELTA_PT_HISTOGRAMS/rmatrix$SUFF
mkdir -p $PATH_TO_DELTA_PT_HISTOGRAMS/rmatrix_$SUFF

for RPARAM in 0.2 0.3 0.4 #0.5
do
	export RPARAM
for PTLEAD in `echo $PTLEADCUTS`
do
   export PTLEAD
	root4star -l buildResponseM.C -q -b
	#root4star -l '/gpfs01/star/pwg/licenrob/jets/StRoot/macros/loadSharedHFLibraries.C' 'buildResponseM.C' -q -b
done
done
