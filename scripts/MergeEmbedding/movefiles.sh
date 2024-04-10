#!/bin/bash
#moves files into separate directories and merges them


merge_files() { 
    echo " "
    echo "Merging root files..."
    setup 64b          # star system is a 32bit system, to use ROOT6 we need to switch to 64bit
    setup root 6.20.08 # setup the root version which can merge files in parallel
 
    # beware not to use the same folder for your new analysis and old one, or delete your old root files
    hadd -f -k -j ${outputFile} ${inputDir}/*.root # -f force overwrite(current outputFile), -k skip corrupt files, -j parallelize
 
}

PTHARD=(3 5 7 9 11 15 20 25 30 40 50 -1)

year=2024
files=(${year}*)

for (( i=0; i<${#PTHARD[@]}-1 ; i+=1 ))
do
PTMIN="${PTHARD[i]}"
PTMAX="${PTHARD[i+1]}"
dir="${files[i]}"
inputDir="Pythia6_${PTMIN}_${PTMAX}"
mv ${dir} ${inputDir}
outputFile="Pythia6_${PTMIN}_${PTMAX}"
./merge.sh ${inputDir} ${outputFile} &
#merge_files &

done #pThard bins
