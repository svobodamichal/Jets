#!/bin/bash
STARLIB_VER="SL19c"
BASEPATH="/star/u/rusnak/JET_analysis_run11/STARJet"
export STARJETBASEDIR="$BASEPATH" 
#export ANALYSISDIR="$BASEPATH/analysis_$STARLIB_VER"
export ANALYSISDIR="/gpfs01/star/pwg/svomich/Jets"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ANALYSISDIR/lib"

#export ROOUNFOLD="$BASEPATH/software_$STARLIB_VER/RooUnfold/v-trunk-custom"
export ROOUNFOLD="$ANALYSISDIR/RooUnfold/v-trunk-custom"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ROOUNFOLD"

export FASTJETDIR="$BASEPATH/software_$STARLIB_VER/fastjet3"
export PATH="$PATH:$FASTJETDIR/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$FASTJETDIR/lib"

export TOYMODELDIR="/star/u/rusnak/JET_analysis_run11/toymodel"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TOYMODELDIR/Production"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TOYMODELDIR/Analysis"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TOYMODELDIR/Unfolding"
