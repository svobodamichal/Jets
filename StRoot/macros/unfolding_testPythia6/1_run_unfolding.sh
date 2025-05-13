#!/bin/bash
source ../set_paths.sh
starver $STARLIB_VER

#SCRIPT_NAME=`basename -- $0`
#_usage() {
#    echo "Usage: ${SCRIPT_NAME} RMATRIX_TYPE (BG_sp | BG_dete | dete )" 
#    exit 1
#}

#export RMATRIX_TYPE=$1 #BG_sp BG_pyt dete BG_dete - correction for BG (using single particle / pythia jet), detector effects, BG+detector effects 
#export RMATRIX_TYPE=BGD #BG_pyt dete BG_dete - correction for BG (using single particle / pythia jet), detector effects, BG+detector effects 
export RMATRIX_TYPE=emb #RM from embedding
echo "RMATRIX_TYPE: $RMATRIX_TYPE"
#check arguments
#[ -n "$RMATRIX_TYPE" ] || _usage

#prior_type=(flat pythiadete pythia powlaw3 powlaw45 powlaw5 powlaw55 levy levy_alex)
prior_type=(flat flat pythia powlaw4 powlaw45 powlaw5 powlaw55 tsalis_1 tsalis_2 tsalis_3 tsalis_4 tsalis_5 tsalis_6 tsalis_7 tsalis_8 tsalis_9)
#TEST NEW PRIORS
#prior_type=(flat mtsalis_1 mtsalis_2 mtsalis_3 mtsalis_4 mtsalis_5 mtsalis_6 gauspol expol gammapol_1 gammapol_2 gammapol_3 gammapol_4 gammapol_5 gauspol_1 gauspol_2 gauspol_3 gauspol_4 gauspol_5)

BASEPATH="$ANALYSISDIR"
#export SVD=0 # SVD unfolding instead of Bayes
#for SVD in 0 1
for SVD in 0 1
do
export SVD # SVD unfolding instead of Bayes
export SMOOTH=0 #smooth unfolded solutions in between iterations #originally 0
#export NBINS=200
export NBINS=VAR #variable binning, used only for output directory name
export NITER=8 #number of iterations
#export NITER=40 #number of iterations
# FOR PRIOR DISTRIBUTION
export PTCUT=0.2 #min track pT cut
export SECONDUNFOLD=0 #unfold already unfolded results (eg. using a different RM)
	export INPUTITER=4 #if unfolding already unfolded results, which iteration to unfold
#export PTCUTOFF=0 #from witch pT to start with unfolding
#EFFICORR=1 # do efficiency correction
EFFICORR=0 # do efficiency correction
#SUFF2="_GPC2" #output dir suffix
SUFF2="_TEST" #output dir suffix

#for CENTRAL in 0 1 #1:central|0:peripheral|2:pp collisions
for CENTRAL in 1 0 #1:central|0:peripheral|2:pp collisions test
do
#USE2DHISTO=1 # is the histogram with the measured distribution a 2D histogram?
USE2DHISTO=0 # is the histogram with the measured distribution a 2D histogram?

TRG="HT2" #HT2, MB
#TRG="MB"
export TRG

if [ $CENTRAL -eq 1 ]; then
	#SUFFIX="_central"
	SUFFIX="central"
	#PTLEADCUTS="5 6 7"
	PTLEAD_MIN=0
	PTLEAD_MAX=9
	#PTLEAD_MAX=7 #test only
elif [ $CENTRAL -eq 0 ]; then
	#SUFFIX="_peripheral"
	SUFFIX="peripheral"  
	#PTLEADCUTS="4 5 6 7"
	#PTLEAD_MIN=4
	PTLEAD_MIN=0
	PTLEAD_MAX=9
	#PTLEAD_MAX=5 #test only
else #p+p
	SUFFIX="_pp" 
	#PTLEADCUTS=3 #"0 1 3 4 5 6" #"4 5 6"
	PTLEAD_MIN=3
	PTLEAD_MAX=5
	USE2DHISTO=0
	TRG="MBHT"
fi
export SUFFIX
export USE2DHISTO
export PTLEAD_MIN #start with this pTlead cut
export PTLEAD_MAX #finish with this pTlead cut

#WRKDIR="$BASEPATH/out/$TRG"
WRKDIR="$ANALYSISDIR"

if [ $RMATRIX_TYPE == "BG_sp" ]; then
	export EFFICORR=0 
else
	export EFFICORR
fi
if [ $EFFICORR -eq 0 ]; then
EFFSUF=""
else
EFFSUF="_eff"
fi

if [ $SVD -eq 0 ]; then
UTYPE="Bayes"
else
UTYPE="SVD"
fi

for SYSSUF in "_main" #"_pythia" 
do

   if [ $SYSSUF == "_main" ]; then
      TSUFFIX_ARR="_normal" # _pp _g _u _m5 _p5 _v2"
   else
      TSUFFIX_ARR="_normal"
   fi

for TSUFF in `echo $TSUFFIX_ARR`
do
export TSUFF
SYS="closure" #SYSTEMATIC UNCERTAINTY VARIATION
export SYS


for RPARAM in 0.2 0.3 0.4 #0.3
#for RPARAM in 0.2 #0.4 #0.3 test
do
export RPARAM
#export DATA_PATH="$WRKDIR/inclusive${SUFFIX}${SYSSUF}"
export DATA_PATH="$WRKDIR/raw_spectra_to_unfold/${TRG}"
#export PRIOR_PATH="$WRKDIR/prior"
export PRIOR_PATH="$WRKDIR/responseM_${TRG}/prior"
#export RMATRIX_PATH="$WRKDIR/embedding${SUFFIX}${SYSSUF}/rmatrix${TSUFF}"
export RMATRIX_PATH="$WRKDIR/responseM_${TRG}/Pythia6"
#export EPSILON_PATH="$TOYMODELDIR/DataOut/pythia/jetonly/pyEmb_R${RPARAM}_${SUFFIX}${SYSSUF2}${TSUFF}/"
export EPSILON_PATH="$WRKDIR/matchingeffi"

for BININGCH in 6 #7 - equidistant #-test binning #5 #-full jets  #6 #0 1 4 #1 2 3 4 #choice of bining arrays 0: nu=nm, 1: nu<nm, 2: nu<nm
do
export BININGCH 
#for PRIOR in 15 #2 4 5 6 7 8 9 10 11 12 13 14 15 #0: truth, 1: flat, 2: biased pythia, 3: pT^(-3), 4:pT^(-4) 5: pT^(-5) 6:pT^(-6) 7: levy I 8: levy II
#do
#	export PRIOR
#	if [ $SECONDUNFOLD -eq 0 ]; then
		#OUT_DIR="../out_test/unfolding/"
		#OUT_DIR="${DATA_PATH}/Unfolded_R${RPARAM}_${UTYPE}_${NBINS}bins_bining${BININGCH}_${RMATRIX_TYPE}${SUFF2}${TSUFF}" #/${prior_type[$PRIOR]}"
		#OUT_DIR="$ANALYSISDIR/out_test/unfolding/${UTYPE}/${SUFFIX}"
		OUT_DIR="$ANALYSISDIR/Unfolded_${TRG}${EFFSUF}_${SYS}/unfolding/${UTYPE}/${SUFFIX}"
#	else
#		OUT_DIR=$DATA_PATH"/Unfolded_R${RPARAM}_${UTYPE}_${NBINS}bins_U2_initer${INPUTITER}/"${prior_type[$PRIOR]}
#	fi
   echo "creating directory: $OUT_DIR"
	#we will create a subdirectory for each prior
	#for PRIOR in 2 4 5 6 7 8 9 10 11 12 13 14 15 #0: truth, 1: flat, 2: biased pythia, 3: pT^(-3), 4:pT^(-4) 5: pT^(-5) 6:pT^(-6) 7: levy I 8: levy II
	#for PRIOR in {1..} # TEST NEW PRIORS
	#do
		mkdir -p "$OUT_DIR/Pythia6"
	#done
	export OUT_DIR

#	for PTLEAD in `echo $PTLEADCUTS`
#	do
#		export PTLEAD
	  root -l -b -q run_unfolding.C 2>/dev/null
#	done #pT threshold
#done #prior
done #bining
done #R
done #systematic set 1
done #systematic set 2
done #centrality
done #unfolding type
