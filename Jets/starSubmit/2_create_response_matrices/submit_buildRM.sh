#!/bin/bash
source set_paths.sh
TRG="HT2"
BASEDIR=$ANALYSISDIR
LOGDIR="$BASEDIR/starSubmit/2_create_response_matrices/log"
ERRDIR="$BASEDIR/starSubmit/2_create_response_matrices/err"
export WORKDIR="$BASEDIR/StRoot/macros/response_matrix"
export RMTYPE="BG_sp" #deltapT distribution: BG_sp|inplane|outplane
export V2CORR=0 # correct delta pT for event plane bias
export V2PATH="$BASEDIR/StRoot/macros/EP_corrections" #path to correction histograms
SCRIPTNAME=run_buildRM.csh

TIMESTAMP=$(date +%s)

#create response matrices for different systematic effects
for SYSSUF in "main" #"RRho02" "RRho04" "pythiaUG" "nrem-1"  # #"main" "nfit12" "nfit20"  "global" #
do
#centrality class
for CENTRAL in 0 1
do
export CENTRAL
if [ $CENTRAL -eq 0 ]; then
	CENTSUFF="peripheral" 
	PTLEADCUTS="4 5 6 7 9"
	#PTLEADCUTS="9"  
elif [ $CENTRAL -eq 1 ]; then
	CENTSUFF="central" 
	PTLEADCUTS="5 6 7 9"
	#PTLEADCUTS="9"  
else
	CENTSUFF="_pp" 
	PTLEADCUTS="0 1 2 3 4 5 6 7" 
fi

export PATH_TO_DELTA_PT_HISTOGRAMS="$ANALYSISDIR/responseM_${TRG}/${CENTSUFF}"

for RPARAM in 0.2 0.3 0.4 #0.5
#for RPARAM in 0.2 0.4 #0.5
do
export RPARAM
for PTLEAD in `echo $PTLEADCUTS`
do
	export PTLEAD
   JOBNAME="buildRM_R${RPARAM}_pTl${PTLEAD}"

mkdir -p ${PATH_TO_DELTA_PT_HISTOGRAMS}/rmatrix_normal
#echo "creating dir ${PATH_TO_DELTA_PT_HISTOGRAMS}/rmatrix_normal"
if [ $V2CORR -eq 1 ]; then
	mkdir -p ${PATH_TO_DELTA_PT_HISTOGRAMS}/rmatrix_v2
fi


#Prepare job submission
if [ ! -e tmp ]; then
	mkdir -p tmp
fi

TEMPLATE_NAME="BuildRM_${CENTRAL}_${PTLEAD}_${RPARAM}_${SYSSUF}_${TIMESTAMP}.xml"
NFILES=1 #how many files merge into one batch
FILES_PER_HOUR=5

#===========================
#create submission xml file
#===========================
echo "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" > tmp/$TEMPLATE_NAME
echo "<job maxFilesPerProcess=\"$NFILES\" simulateSubmission = \"false\" filesPerHour = \"$FILES_PER_HOUR\">" >> tmp/$TEMPLATE_NAME
echo "<command>" >> tmp/$TEMPLATE_NAME
echo "setenv WORKDIR $WORKDIR" >> tmp/$TEMPLATE_NAME
echo "setenv RMTYPE $RMTYPE" >> tmp/$TEMPLATE_NAME
echo "setenv V2CORR $V2CORR" >> tmp/$TEMPLATE_NAME
echo "setenv V2PATH $V2PATH" >> tmp/$TEMPLATE_NAME
echo "setenv CENTRAL $CENTRAL" >> tmp/$TEMPLATE_NAME
echo "setenv PATH_TO_DELTA_PT_HISTOGRAMS $PATH_TO_DELTA_PT_HISTOGRAMS" >> tmp/$TEMPLATE_NAME
echo "setenv RPARAM $RPARAM" >> tmp/$TEMPLATE_NAME
echo "setenv PTLEAD $PTLEAD"  >> tmp/$TEMPLATE_NAME
echo "  starver $STARLIB_VER" >> tmp/$TEMPLATE_NAME
#echo "  source $ANALYSISDIR/set_paths.csh" >> tmp/$TEMPLATE_NAME
echo "  cd $ANALYSISDIR/" >> tmp/$TEMPLATE_NAME
echo "  pwd" >> tmp/$TEMPLATE_NAME
echo "  root4star -l 'StRoot/macros/loadSharedHFLibraries.C' 'StRoot/macros/response_matrix/buildResponseM.C' -q -b" >> tmp/$TEMPLATE_NAME
echo "  </command>" >> tmp/$TEMPLATE_NAME
echo "  <stdout URL=\"file:$LOGDIR/\$JOBID.log\"/>" >> tmp/$TEMPLATE_NAME
echo "  <stderr URL=\"file:$ERRDIR/\$JOBID.err\"/>" >> tmp/$TEMPLATE_NAME
echo "  <SandBox>" >> tmp/$TEMPLATE_NAME
echo "  <Package>" >> tmp/$TEMPLATE_NAME
echo "  <File>file:$ANALYSISDIR/StRoot/macros/response_matrix/buildResponseM.C</File>" >> tmp/$TEMPLATE_NAME
echo "  </Package>" >> tmp/$TEMPLATE_NAME
echo "  </SandBox>" >> tmp/$TEMPLATE_NAME
echo "</job>" >> tmp/$TEMPLATE_NAME

#let's submit
	cd tmp
	star-submit $TEMPLATE_NAME 
	cd ..



done #pTlead
done #R
done #centrality
done #systematic error type
