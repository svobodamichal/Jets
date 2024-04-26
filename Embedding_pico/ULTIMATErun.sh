#!/bin/bash

pthat=(3 5 7 9 11 15 20 25 30 40 50 -1)
#PYTHIA6 cross sections [mb] for pT-hat bins (3 5 7 9 11 15 20 25 30 40 50 -1) to be used as weights
weights=(1.616e+0 1.355e-01 2.288e-02 5.524e-03 2.203e-03 3.437e-04 4.681e-05 8.532e-06 2.178e-06 1.198e-07 6.939e-09)
#pthat=(50 -1)
#weights=(6.939e-09)


for (( i=0; i<${#pthat[@]}-1 ; i+=1 ))
do
pthatmin="${pthat[i]}"
pthatmax="${pthat[i+1]}"
xweight="${weights[i]}"

echo "$pthatmin $pthatmax $xweight"
list="filelists/pythia6picoDsts_${pthatmin}_${pthatmax}.list"

#pthatmin="25"
#pthatmax="30"
#xweight="8.532e-06"

#echo "$pthatmin $pthatmax $xweight"
#list="test.list"

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

echo executing submitPicoHFJetMaker_filelist.csh f0r $list
csh starSubmit/submitPicoHFJetMaker_filelist.csh $path $list $pthatmin $pthatmax $xweight &
sleep 60s

done
