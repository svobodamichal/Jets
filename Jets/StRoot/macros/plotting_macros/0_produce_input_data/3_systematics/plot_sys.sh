#!/bin/bash
#this macro creates figures of systematic errors for several systems, Rs,pTls,...

SUFFIX="GPC2" #name of the subdirectory containing input data
TOYMODEL=0
PTLEAD_DENOM=5.0 #pTlead cut of denumerator in spectra ratio with two different pTlead cuts

	for BININGCH in 1 #0 1 4 #1 2 3 4 #2 #binning set
	do
	
	for SYSTEM in  cent #pp # 
	do
	if [ $SYSTEM == cent ]; then #central collisions
		PTLEADCUTS="5" # 6 7"
	elif [ $SYSTEM == peri ]; then #peripheral collisions
		PTLEADCUTS="5" #4 5 6 7"
	else #pp data
		PTLEADCUTS="5" #"4 5 6" 
	fi

	for R in 0.2 0.3 0.4 #0.5
	do
	
	for PTLEAD in `echo $PTLEADCUTS`
	do
	for SYSTYPE in "unfolding" #"correlated"
	do
		root -q -l -b systematics.C\($SYSTYPE,sys,$SYSTEM,$R,$PTLEAD,$PTLEAD_DENOM,\"$SUFFIX\",$BININGCH,\"normal\","PYTHIA",0,$TOYMODEL,0\)
	done
	done #pTlead
	done #R
	done #centrality
	done #bining
