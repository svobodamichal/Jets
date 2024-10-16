#!/bin/bash

#copylist
#cp /global/homes/k/kvapil/StPicoDpmMakerSL16d/picoLists/picoList_all.list ./

#divide list
#list=${1:-"pico_list_all.list"}
list=${1:-"pico_low_14.list"}
#list=${1:-"pico_mid_14.list"}
#list=${1:-"pico_high_14.list"}

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

echo executing submitPicoHFJetMaker.csh f0r $list
csh starSubmit/submitPicoHFJetMaker.csh $path $list




