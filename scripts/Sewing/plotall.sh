#!/bin/bash

rm pythia6_all.root

PTHARD=(3 5 7 9 11 15 20 25 30 40 50 -1)

for (( i=0; i<${#PTHARD[@]}-1 ; i+=1 ))
do
PTMIN="${PTHARD[i]}"
PTMAX="${PTHARD[i+1]}"

root 'plothisto.cxx("'${PTMIN}'_'${PTMAX}'")' -b -q

done #pThard bins

hadd pythia6_all.root pythia6_normalized_*.root 
