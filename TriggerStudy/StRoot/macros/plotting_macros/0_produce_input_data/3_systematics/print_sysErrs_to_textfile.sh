#!/bin/bash
#this macro creates textfiles in ./sysErr subdirectory with values of relative errors of individual systematic effects
SUFFIX="GPC2" #"newBinning2" #"eta" # 
TOYMODEL=0
PPBASE="PYTHIA" #PYTHIA | STAR
PPBIAS=0 #1: use pythia with pTlead cut as the pp baseline for RAA | 0: use PYTHIA with no pTlead cut
#TEST=1 # 0:chi2 1:Camberra 2: Kolmogorov-Smirnov 
RATIO="sys"
PTLEAD_DENOM=5.0

for BININGCH in 1 #0 4 #binning set
do

	for SYSTYPE in "unfold" "corr"
	do
	mkdir -p ./sysErr/${SUFFIX}/${SYSTYPE}/bin${BININGCH}/
	done
	
	for SYSTEM in peri cent #pp # 
do
	if [ $SYSTEM == cent ]; then #central collisions
		PTLEADCUTS="5" # 6 7"
		PTLEADCUTS_DENOM="5"
		CENTDIR="central"
	elif [ $SYSTEM == peri ]; then #peripheral collisions
		PTLEADCUTS="5" #"4 5 6 7"
		PTLEADCUTS_DENOM="5"
		CENTDIR="peripheral"
	else #pp data
		PTLEADCUTS="5" #"4 5 6" 
		PTLEADCUTS_DENOM="4"
		CENTDIR="pp"
	fi

	for R in 0.2 0.3 0.4 #0.5
	do
	for PTLEAD in `echo $PTLEADCUTS`
	do
	for UNC_T in correlated unfolding #normalization
	do
		root -q -l -b systematics.C\($UNC_T,$RATIO,$SYSTEM,$R,$PTLEAD,$PTLEAD_DENOM,\"$SUFFIX\",$BININGCH,\"normal\",$PPBASE,$PPBIAS,$TOYMODEL,0,1\)
	done #uncert
	done #pTlead
	done #R
done #centrality
done #bining
