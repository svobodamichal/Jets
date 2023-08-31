#!/bin/bash
source set_paths.sh
SCRIPT_NAME=`basename -- $0`
_usage() {
    echo "Usage: ${SCRIPT_NAME} RMATRIX_TYPE (BG_sp | BG_pyt | dete | BGD)"
    exit 1
}

export RMATRIX_TYPE=$1 #BG_sp BG_pyt dete BG_dete - correction for BG (using single particle / pythia jet), detector effects, BG+detector effects
echo "TYPE: $RMATRIX_TYPE"
#check arguments
[ -n "$RMATRIX_TYPE" ] || _usage


TRG="MB_QA"
BASEDIR=$ANALYSISDIR
LOGDIR="$BASEDIR/submitter/log"
ERRDIR="$BASEDIR/submitter/err"
export WORKDIR="$BASEDIR/macros/response_matrix"
export NEVENTS=250000000 #number of generated events

for CENTRAL in 1 0
do
if [ $CENTRAL -eq 1 ]; then
CENTSUFF="_central"
PTLEADCUTS="5 6 7"
elif [ $CENTRAL -eq 0 ]; then
CENTSUFF="_peripheral"
PTLEADCUTS="3 4 5 6 7"
else #p+p
CENTSUFF="_pp"
PTLEADCUTS="3 4" #"0 1 2 3 4 5 6 7"
fi

#initial pT-distribution shapes for the unfolding
#0: from input histogram, 1: flat, 2: biased pythia, 3: pT^(-3), 4:pT^(-4.5) 5: pT^(-5) 6:pT^(-5.5) 7-15: Tsalis function parametrizations

export PRIOR_START=2 #starting prior
export PRIOR_STOP=15 #last prior
export PRIOR_SKIP1=3 #skip this prior
export PRIOR_SKIP2=3 #skip this prior
export PRIOR_SKIP3=3 #skip this prior


#PRIORS="2 4 5 6 7 8 9 10 11 12 13 14 15" 
#for PRIOR in `echo $PRIORS`
#do
#	export PRIOR

#datasets (for systematic studies)
for SYSSUF in "_main" #"_nfit14" "_nfit18" #"_pythiaU" "_pythiaUG" "_RRho04" "_nrem-1" "_pythiaG" "_RRho02" #"_main" ## "_global" 
do

if [ $SYSSUF == "_main" ]; then
	#embedding sets (for systematic studies)
	TSUFFIX_ARR="_normal" #_pp _g _u _m5 _p5 _v2"  
else
   TSUFFIX_ARR="_normal"
fi


for TSUFF in `echo $TSUFFIX_ARR`
do
for RPARAM in 0.2 0.3 0.4
do
export RPARAM
export RM_PATH="$ANALYSISDIR/out/${TRG}/embedding${CENTSUFF}${SYSSUF}/rmatrix${TSUFF}"
export PYEMB_PATH="$TOYMODELDIR/DataOut/pythia/jetonly/pyEmb_R${RPARAM}${CENTSUFF}${TSUFF}"
export PRIOR_PATH="$ANALYSISDIR/out/prior"

for PTTHRESH in `echo $PTLEADCUTS`
do
      export PTTHRESH
		mkdir -p ${RM_PATH}


#Prepare job submission
if [ ! -e tmp ]; then
	mkdir -p tmp
fi

TEMPLATE_NAME="buildRMROO_${CENTSUFF}_${RMATRIX_TYPE}_R${RPARAM}_pTlead${PTTHRESH}${SYSSUF}_${TSUFF}_prior${PRIOR_START}.xml"
JOBNAME="RMROO_"
NFILES=1 #how many files merge into one batch
FILES_PER_HOUR="0.1"

#===========================
#create submission xml file
#===========================
echo "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" > tmp/$TEMPLATE_NAME
echo "<job maxFilesPerProcess=\"$NFILES\" simulateSubmission = \"false\" filesPerHour = \"$FILES_PER_HOUR\">" >> tmp/$TEMPLATE_NAME
echo "<command>" >> tmp/$TEMPLATE_NAME
echo "setenv RMATRIX_TYPE $RMATRIX_TYPE" >> tmp/$TEMPLATE_NAME
echo "setenv WORKDIR $WORKDIR" >> tmp/$TEMPLATE_NAME
echo "setenv NEVENTS $NEVENTS" >> tmp/$TEMPLATE_NAME
echo "setenv PRIOR_START $PRIOR_START" >> tmp/$TEMPLATE_NAME
echo "setenv PRIOR_STOP $PRIOR_STOP" >> tmp/$TEMPLATE_NAME
echo "setenv PRIOR_SKIP1 $PRIOR_SKIP1" >> tmp/$TEMPLATE_NAME
echo "setenv PRIOR_SKIP2 $PRIOR_SKIP2" >> tmp/$TEMPLATE_NAME
echo "setenv PRIOR_SKIP3 $PRIOR_SKIP3" >> tmp/$TEMPLATE_NAME
echo "setenv RPARAM $RPARAM" >> tmp/$TEMPLATE_NAME
echo "setenv RM_PATH $RM_PATH" >> tmp/$TEMPLATE_NAME
echo "setenv PYEMB_PATH $PYEMB_PATH" >> tmp/$TEMPLATE_NAME
echo "setenv PRIOR_PATH $PRIOR_PATH" >> tmp/$TEMPLATE_NAME
echo "setenv PTTHRESH $PTTHRESH" >> tmp/$TEMPLATE_NAME
echo "  starver $STARLIB_VER" >> tmp/$TEMPLATE_NAME
echo "  source $ANALYSISDIR/set_paths.csh" >> tmp/$TEMPLATE_NAME
echo "  cd $ANALYSISDIR/macros/response_matrix" >> tmp/$TEMPLATE_NAME
echo "  pwd" >> tmp/$TEMPLATE_NAME
echo "  root4star -q -b -l buildResponseROO.C" >> tmp/$TEMPLATE_NAME
echo "  </command>" >> tmp/$TEMPLATE_NAME
echo "  <stdout URL=\"file:$LOGDIR/$JOBNAME\$JOBID.log\"/>" >> tmp/$TEMPLATE_NAME
echo "  <stderr URL=\"file:$ERRDIR/$JOBNAME\$JOBID.err\"/>" >> tmp/$TEMPLATE_NAME
echo "  <output fromScratch=\"*.root\" toURL=\"file:$RM_PATH/\"/>" >> tmp/$TEMPLATE_NAME
echo "  <SandBox>" >> tmp/$TEMPLATE_NAME
echo "  <Package>" >> tmp/$TEMPLATE_NAME
echo "  <File>file:$ANALYSISDIR/macros/response_matrix/buildResponseROO.C</File>" >> tmp/$TEMPLATE_NAME
echo "  </Package>" >> tmp/$TEMPLATE_NAME
echo "  </SandBox>" >> tmp/$TEMPLATE_NAME
echo "</job>" >> tmp/$TEMPLATE_NAME

#let's submit
	cd tmp
	star-submit $TEMPLATE_NAME 
	cd ..

done #pTtresh
done #R
done #TYPE SUFFIX
done #SYS SUFFIX
#done #prior
done #centrality
