#!/bin/bash
source ../set_paths.sh
starver $STARLIB_VER
echo "starver $STARLIB_VER"

TRG="MB"
export RPARAM=0.3
export PTTHRESH=9
export RMATRIX_TYPE="BGD" #deltapT distribution: BG_sp|inplane|outplane
export RM_PATH="$ANALYSISDIR/out/${TRG}/embedding_central_main/rmatrix_normal"
export PYEMB_PATH="$TOYMODELDIR/DataOut/pythia/jetonly/pyEmb_R${RPARAM}_central_normal"
export PRIOR_PATH="$ANALYSISDIR/out/prior"
export PRIOR_START=2
export PRIOR_STOP=15
export PRIOR_SKIP1=3
export PRIOR_SKIP2=3
export PRIOR_SKIP3=3
export NEVENTS=100000


	root4star -l buildResponseROO.C -q -b
