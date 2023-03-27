#!/bin/bash
EVO="GPC2" #input data set ID label
EVO_TOY="GPC2" #input data set ID label for Toymodel
EVO_SYS_CORR="GPC2" #ID label for data set with systematic (correlated) uncertainties
LABEL=""
EXT="eps"
BINING=1
PTLEAD_MAIN=5.0
PTLEAD_DENUMC=5.0 #pTlead cut value for the denumerator of spectra ratio with two different pTlead cuts, central collisions
PTLEAD_DENUMP=5.0 #pTlead cut value for the denumerator of spectra ratio with two different pTlead cuts, peripheral collisions
FIGURE_PATH="../plotting_out/obr/$EVO"

IS4PAPER=1 #are these plots for paper or for thesis/AN
#if [ $IS4PAPER -eq 0 ]; then
#LABEL="THIS THESIS"
#fi


#which plots do we want to draw?
DO_RAW_R=0 #raw spectrum, R-dependency
DO_RAW_PTL=0 #raw spectrum, pTlead dependency
DO_PYTHIAvsSTAR=0 #uncorrected STAR data vs smeared PYTHIA
DO_STARvsTOY=0 #STAR data vs TOYMODEL
DO_AREA=0 #jet area
DO_DPT=0 #delta-pT
DO_RM=0 #Response Matrix
DO_EFFI=0 #jet reconstruction efficiency
DO_JER=0 #jet energy resolution
DO_TOY_CLOSURE=0 #toymodel closure test
DO_PYTHIAPP=1 #PYTHIA 6 pp-baseline plots
DO_FINAL=0 #final plots (spectra, RAA, RCP, ratios R1/R2, pTlead1/pTlead2)
DO_FINAL_BIASED=0 #final plots, biased pp baseline in RAA
DO_COMPARE=0 #comparison plots with theory

#scripts to be used
SCRIPTDIR1=1_side_results
SCRIPT1=plot_whatever
SCRIPTDIR2=4_main_results
SCRIPT2=plot_FINAL
SCRIPTDIR3="RAA_RCP_theory_comparison"
SCRIPT3=RAA_RCP
SCRIPTDIR4="RR_ratio_theory_comparison"
SCRIPT4=R2R4
SCRIPTDIR5="3_additional_studies/JER/"
SCRIPT5=plot_JER
SCRIPTDIR6="3_additional_studies/BGdete_on_PYTHIA"
SCRIPT6=plot
SCRIPTDIR7="2_pythia_baseline"
SCRIPT7a=plotptleadratio2
SCRIPT7b=plotPhenixNew2

#create output dirs
mkdir -p $FIGURE_PATH/results/comparison
mkdir -p $FIGURE_PATH/results/central
mkdir -p $FIGURE_PATH/results/toymodel/central
mkdir -p $FIGURE_PATH/results/peripheral
mkdir -p $FIGURE_PATH/uncorrected
mkdir -p $FIGURE_PATH/area
mkdir -p $FIGURE_PATH/dpT
mkdir -p $FIGURE_PATH/effi
mkdir -p $FIGURE_PATH/pp
mkdir -p $FIGURE_PATH/RM
mkdir -p $FIGURE_PATH/toy


PTLEAD_CENT="5.0 6.0 7.0"
PTLEAD_PERI="4.0 5.0 6.0 7.0"

cd $SCRIPTDIR1
echo "cd $SCRIPTDIR1"
pwd
#jet area and rho
if [ $DO_AREA -eq 1 ]; then
for PYTHIAEMB in 0 1 #plot jet area for SP or PYTHIA embedded jets?
do	
for CENTRALITY in "cent" "peri"
do
root -l -b <<EOF
.L $SCRIPT1.C
plot_area_and_rho($CENTRALITY, 0, 0, "$EVO","$LABEL",$IS4PAPER,$PYTHIAEMB,"$EXT")
.q
EOF
done
done
fi

#STAR vs Toymodel
if [ $DO_STARvsTOY -eq 1 ]; then
root -l -b <<EOF
.L $SCRIPT1.C
plot_measured_vs_toymodel(0.0, cent,"$EVO","$EVO_TOY", "$EXT")
.q
EOF
fi

for SYS in "peri" "cent"
do


if [ $SYS == "peri" ]; then
PTLEAD_ARR=$PTLEAD_PERI
else
PTLEAD_ARR=$PTLEAD_CENT
fi

#raw spectrum - R dependence
if [ $DO_RAW_R -eq 1 ]; then
root -l -b <<EOF	
.L $SCRIPT1.C
plot_measured_R(0.0, $SYS, "$EVO" ,"$LABEL", "$EXT")
.q
EOF
fi


for R in 0.2 0.3 0.4
do
#raw spectrum - pTlead dependence
if [ $DO_RAW_PTL -eq 1 ]; then
root -l -b <<EOF
.L $SCRIPT1.C
plot_measured_pTlead($R, $SYS, "$EVO", "$LABEL", "$EXT")
.q
EOF
fi

#delta-pT histograms
if [ $DO_DPT -eq 1 ]; then
root -l -b <<EOF
.L $SCRIPT1.C
plot_dpT($R, $PTLEAD_MAIN, 0, $SYS,"$EVO","$LABEL","$EXT")
.q
EOF
fi

#Response Matrix and its Y projection
if [ $DO_RM -eq 1 ]; then
root -l -b <<EOF
.L $SCRIPT1.C
plot_RM($R, $PTLEAD_MAIN, 0, $SYS, "$EVO", "$LABEL", "$EXT")
.q
EOF
fi


#Response Matrix and its Y projection
if [ $DO_EFFI -eq 1 ]; then
root -l -b <<EOF
.L $SCRIPT1.C
plot_epsilon($R, $PTLEAD_MAIN, $SYS, "$EVO","$LABEL","$EXT")
.q
EOF
fi

done #R loop
done #SYS

cd ..
#echo "cd .."
#pwd

#===================
#final results
#===================
cd $SCRIPTDIR2
echo "cd $SCRIPTDIR2"
pwd

if [ $DO_TOY_CLOSURE -eq 1 ]; then
root -l -b <<EOF
.L $SCRIPT2.C
$SCRIPT2($PTLEAD_MAIN, $PTLEAD_DENUMC, cent, "$EVO", "$EVO_SYS_CORR",$BINING, 1, 1, 0, "$LABEL","$EXT")
.q
EOF
fi #toymodel closure test

for SYS in "peri" "cent"
do

PTLEAD_ARR=$PTLEAD_CENT
PTLEAD_DENUM=$PTLEAD_DENUMC
if [ $SYS == "peri" ]; then
PTLEAD_DENUM=$PTLEAD_DENUMP
PTLEAD_ARR=$PTLEAD_PERI
fi

for PTLEAD in `echo $PTLEAD_ARR`
do
if [ $DO_FINAL -eq 1 ]; then
echo "plotting $SYS-results for pTlead $PTLEAD and pTlead_denum $PTLEAD_DENUM"
root -l -b <<EOF
.L $SCRIPT2.C
$SCRIPT2($PTLEAD, $PTLEAD_DENUM, $SYS, "$EVO", "$EVO_SYS_CORR", $BINING, 0, 0, 0, "$LABEL","$EXT")
.q
EOF
fi #final results with unbiased pythia

#final results, biased pythia
if [ $DO_FINAL_BIASED -eq 1 ]; then
root -l -b <<EOF
.L $SCRIPT2.C
$SCRIPT2($PTLEAD, $PTLEAD_DENUM, $SYS, "$EVO", "$EVO_SYS_CORR", $BINING, 0, 0, 1, "$LABEL","$EXT")
.q
EOF
fi #final results with biased pythia
done #PTLEAD
done #SYS

if [ $DO_COMPARE -eq 1 ]; then
for CENTRAL in cent #peri
do
#if [ CENTRAL == "cent" ]; then
PTLEAD_ARR=$PTLEAD_CENT
#else
#PTLEAD_ARR=$PTLEAD_PERI
#fi
for PTLEAD in `echo $PTLEAD_ARR`
do
#comparison with theory - RAA, RCP
cd $SCRIPTDIR3
echo "cd $SCRIPTDIR3"
pwd
root -l -b <<EOF
.L $SCRIPT3.C
$SCRIPT3($CENTRAL, $PTLEAD, $BINING,"$EVO", "$LABEL",0,theory,"$EXT")
.q
EOF
root -l -b <<EOF
.L $SCRIPT3.C
$SCRIPT3($CENTRAL, $PTLEAD, $BINING,"$EVO", "$LABEL",0,hadrons,"$EXT")
.q
EOF

#comparison with theory - RAA, RCP, biased PYTHIA
root -l -b <<EOF
.L $SCRIPT3.C
$SCRIPT3($CENTRAL, $PTLEAD, $BINING,"$EVO", "$LABEL",1,theory,"$EXT")
.q
EOF
root -l -b <<EOF
.L $SCRIPT3.C
$SCRIPT3($CENTRAL, $PTLEAD, $BINING,"$EVO", "$LABEL",1,hadrons,"$EXT")
.q
EOF

cd ..
#echo "cd .."
#pwd
done #PTLEAD
done #CENTRAL

for PTLEAD in `echo $PTLEAD_CENT`
do
#comparison with theory - R1/R2 ratios
cd $SCRIPTDIR4
echo "cd $SCRIPTDIR4"
pwd
root -l -b <<EOF
.L $SCRIPT4.C
$SCRIPT4($PTLEAD,$BINING,"$EVO","$LABEL","$EXT")
.q
EOF
cd ..
#echo "cd .."
#pwd
done #PTLEAD
fi #comparison plots with theory

cd ..

#Jet energy resolution
if [ $DO_JER -eq 1 ]; then
cd $SCRIPTDIR5
echo "cd $SCRIPTDIR5"
pwd
for PTLEAD in `echo $PTLEAD_CENT`
do
for R in 0.2 0.3 0.4
do
for QUARK in "u" "g"
do
root -l -b <<EOF
.L $SCRIPT5.C
plot_JER_AllinOne($R,$PTLEAD,"$QUARK","$EVO","$EXT")
.q
EOF
done
done
done
cd ..
cd ..
#echo "cd .."
#pwd
fi #JER


#magnitude of the detector and delta-pT effects
if [ $DO_PYTHIAvsSTAR -eq 1 ]; then
cd $SCRIPTDIR6
echo "cd $SCRIPTDIR6"
pwd
root -l -b <<EOF
.L $SCRIPT6.C
$SCRIPT6(0.3,$PTLEAD_MAIN,"$EVO","$EXT")
.q
EOF
cd ..
cd ..
fi #smeared Pythia

#PYTHIA6 pp baseline
if [ $DO_PYTHIAPP -eq 1 ]; then
cd $SCRIPTDIR7
echo "cd $SCRIPTDIR7"
pwd
root -l -b <<EOF
.L $SCRIPT7a.C
$SCRIPT7a("$EVO", "$EXT")
.q
EOF

root -l -b <<EOF
.L $SCRIPT7b.C
$SCRIPT7b("$EVO", "$EXT")
.q
EOF
cd ..
fi #PYTHIA6 pp baseline

