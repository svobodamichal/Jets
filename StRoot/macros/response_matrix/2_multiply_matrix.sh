#!/bin/bash
#export JETTYPE="pythia"
source ../set_paths.sh
starver $STARVLIB_VER

export REVERSE=0 #0: BGxDete 1: DetexBG


for CENTRAL in 0 1 #central or peripheral collisions
do

if [ $CENTRAL -eq 1 ]; then
	CSUFFIX="_central" #centrality suffix
	PTLEADCUTS="5 6 7"
elif [ $CENTRAL -eq 0 ]; then
	CSUFFIX="_peripheral"
	PTLEADCUTS="4 6 5 7"
else
	CSUFFIX="_pp"
	PTLEADCUTS=1 #"0 3" #"4 6 5"
fi

for SYSSUF in  "_main" #"_nfit14" "_nfit18" #"_main" "_pythiaU" "_pythiaUG" "_RRho02" "_pythiaG"  "_RRho04" "_nrem-1"  #suffix of the systematic effect
do
	DATA_PATH="$ANALYSISDIR/out/MB/embedding${CSUFFIX}${SYSSUF}"
	TOY_PATH="$TOYMODELDIR/DataOut/pythia/jetonly"

	if [ $SYSSUF == "_main" ]; then
		TSUFFIX_ARR="_normal" # _v2 _pp _g _u _m5 _p5" #suffix of the detector effect
	else
		TSUFFIX_ARR="_normal"
	fi


	for TSUFFIX_OUT in `echo $TSUFFIX_ARR`
	do
		#export TSUFFIX_OUT #output directory suffix
		TSUFFIX1="_normal" #1st input dir suffix
		TSUFFIX2=$TSUFFIX_OUT #2nd input dir suffix
	
		if [ $SYSSUF == "_global" ]; then
			TSUFFIX2="_global"
		elif [ $SYSSUF == "_nfit12" ]; then
			TSUFFIX2="_nfit12"
		elif [ $SYSSUF == "_nfit20" ]; then
			TSUFFIX2="_nfit20"
		fi
	
		if [ $TSUFFIX_OUT == "_v2"  ]; then
			TSUFFIX1="_v2"
		fi

		#echo "$TSUFFIX1 X $TSUFFIX2 -> $TSUFFIX_OUT"
		
		export DIR_IN1="$DATA_PATH/rmatrix$TSUFFIX1"
		export DIR_OUT="$DATA_PATH/rmatrix$TSUFFIX_OUT"

		mkdir -p $DIR_OUT

		for RPARAM in 0.2 0.3 0.4 
		do
			export RPARAM
			export DIR_IN2="$TOY_PATH/pyEmb_R${RPARAM}$CSUFFIX$TSUFFIX2"
			for PTTHRESH in `echo $PTLEADCUTS` 
		   do
			   export PTTHRESH
			   root -b -l multiply_matrix.C -q
done
done
done
done
done
