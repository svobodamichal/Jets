#!/bin/tcsh

#starpro

#root4star -l -b -q 'runPicoDpmAnaMaker.C("test_pico.list","output_test",0,"BadRunList_MB.list","picoHFtree","root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_2016/ReversedFullField/P16ij/2016",0)'

#root4star -l -b -q 'runPicoDpmAnaMaker.C("test_pico.list","output_test",0,"BadRunList_MB.list","picoHFtree","/gpfs01/star/pwg/robotmon/myDpmAnalysis/",0)'

#root4star -l -b -q 'runPicoHFJetMaker.C("root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu_200_production_2016/ReversedFullField/P16ij/2016/038/17038050/st_physics_17038050_raw_1000002.picoDst.root","output_test",0,"PicoDst")'

#root4star -l -b -q 'runPicoHFJetMaker.C("test.list","output_test",0,"BadRunList_14.list","PicoDst")'
#root -l -b -q 'StRoot/macros/loadSharedHFLibraries.C' 'jobs/2019-04-23_07-19/StRoot/macros/runPicoHFJetMaker.C++("test.list","output_test",0,"BadRunList_14.list","PicoDst")'
#root -l -b -q 'StRoot/macros/loadSharedHFLibraries.C' 'StRoot/macros/runPicoHFJetMaker.C++("test.list","output_test",0,"BadRunList_14.list","HotTowerList_14.list","PicoDst")'
gdb --quiet --args root4star -l -b -q 'StRoot/macros/loadSharedHFLibraries.C' 'StRoot/macros/runPicoHFJetMaker.C("test.list","output_test",0,"BadRunList_14.list","PicoDst")'
#root4star -l -b -q 'StRoot/macros/loadSharedHFLibraries.C' 'StRoot/macros/runPicoHFJetMaker.C("test.list","output_test",0,"BadRunList_14.list","PicoDst")'
