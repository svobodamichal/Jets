#!/bin/tcsh

starpro

set pthat=(3 5 7 9 11 15 20 25 30 40 50 -1)
#PYTHIA6 cross sections [mb] for pT-hat bins (3 5 7 9 11 15 20 25 30 40 50 -1) to be used as weights
set weights=(1.616e+0 1.355e-01 2.288e-02 5.524e-03 2.203e-03 3.437e-04 4.681e-05 8.532e-06 2.178e-06 1.198e-07 6.939e-09)

set pthatmin = 50
set pthatmax = -1
set xweight = 6.939e-09
root4star -l -b -q 'StRoot/macros/loadSharedHFLibraries.C' 'StRoot/macros/runPicoHFJetMaker.C("test.list","output_test",0,"BadRunList_14.list","PicoDst",'${pthatmin}', '${pthatmax}', '$xweight')'
#root -l -b -q 'StRoot/macros/loadSharedHFLibraries.C' 'StRoot/macros/runPicoHFJetMaker.C++("test.list","output_test",0,"BadRunList_14.list","HotTowerList_14.list","PicoDst")'
#gdb --quiet --args root4star -l -b -q 'StRoot/macros/loadSharedHFLibraries.C' 'StRoot/macros/runPicoHFJetMaker.C("test.list","output_test",0,"BadRunList_14.list","PicoDst",'${pthatmin}', '${pthatmax}', '$xweight')'
