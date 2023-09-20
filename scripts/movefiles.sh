#!/bin/bash
#moves files into separate directories and merges them
PTHARD=(3 5 7 9 11 15 20 25 30 40 50 -1)

for (( i=0; i<${#PTHARD[@]}-1 ; i+=1 ))
do
PTMIN="${PTHARD[i]}"
PTMAX="${PTHARD[i+1]}"
./merge.sh Pythia6_${PTMIN}_${PTMAX} Pythia6_${PTMIN}_${PTMAX} &

done #pThard bins
