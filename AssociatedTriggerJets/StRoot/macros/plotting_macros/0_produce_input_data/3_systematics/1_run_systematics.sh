#!/bin/bash
#this macro produces root files containing TGraphs of final results with several kinds of systematic errors
SUFFIX="GPC2" #"newBinning2" #"eta" # 
TOYMODEL=1
PPBASE="PYTHIA" #PYTHIA | STAR
#TEST=1 # 0:chi2 1:Camberra 2: Kolmogorov-Smirnov 

for PPBIAS in 0 #1: use pythia with pTlead cut as the pp baseline for RAA | 0: use PYTHIA with no pTlead cut
do
if [ $PPBIAS -eq 0 ]; then
	SUFFIX2=${SUFFIX}
else
	SUFFIX2="${SUFFIX}_ppBiased"
fi

#if [ $TEST -eq 0 ]; then
#	TESTDIR="chi2"
#elif [ $TEST -eq 1 ]; then
#	TESTDIR="Camb"
#else
#	TESTDIR="KS"
#fi


for BININGCH in 1 #0 4 #1 2 3 4 #2 #binning set
do
	for SYSTEM in cent #peri #pp # 
do
	TOYDIR=""
	if [ $SYSTEM == cent ]; then #central collisions
		PTLEADCUTS="5 6 7"
		PTLEADCUTS_DENOM="5 6"
		RATIOS="RAA RR plpl RCP"
		if [ $TOYMODEL -eq 1 ]; then
		RATIOS="closure RAA RR plpl"
		TOYDIR="toymodel/"
		fi
		CENTDIR="central"
	elif [ $SYSTEM == peri ]; then #peripheral collisions
		PTLEADCUTS="4 5 6 7"
		PTLEADCUTS_DENOM="4 5 6"
		RATIOS="RAA RR plpl"
		CENTDIR="peripheral"
	else #pp data
		PTLEADCUTS="5" #"4 5 6" 
		PTLEADCUTS_DENOM="4"
		RATIOS="RAA RR plpl"
		CENTDIR="pp"
	fi

	OUTDIR="../../../plotting_out/systematics/${TOYDIR}${CENTDIR}/${SUFFIX2}/bining$BININGCH"
	mkdir -p $OUTDIR
	echo "WARNING: Removing files in $OUTDIR"
	rm $OUTDIR/*.root
	#rm $OUTDIR/*_pT6.root


	for RAT_T in `echo $RATIOS` 
	do
	for R in 0.2 0.3 0.4 #0.5
	do
	for PTLEAD in `echo $PTLEADCUTS`
	do
	for PTLEAD_DENOM in `echo $PTLEADCUTS_DENOM`
	do
	for UNC_T in unfolding #correlated #normalization
	do
		root -q -l -b systematics.C\($UNC_T,$RAT_T,$SYSTEM,$R,$PTLEAD,$PTLEAD_DENOM,\"$SUFFIX\",$BININGCH,\"normal\",$PPBASE,$PPBIAS,$TOYMODEL,1\)
	done #uncert
	done #pTlead_denom
	done #pTlead
	done #R
	done #ratio
done #centrality
done #bining
done #pythia bias
