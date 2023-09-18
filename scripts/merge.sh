#!/bin/bash


# baseDir=2015-08-28_16-27

baseDir=$1
output_name=${2:-output}

# ----------------------------------------------------

mergeFolder () {

    inSplitList=inlist.${mergeLevel}.list

    let mergeLevel=mergeLevel+1
    outSplitListFolder=splitLists.$mergeLevel

    # -- split
    mkdir -p $outSplitListFolder
    split -l 10 -d -a 3 $inSplitList ${outSplitListFolder}/inlist.sublist.${mergeLevel}. #5 files
    
    for ii in `find  ${outSplitListFolder} -name "inlist.sublist.${mergeLevel}.*"` ; do 
	sub=`echo $ii | cut -d '.' -f 5`
	newFile=merged_output_${mergeLevel}_${sub}.root
	echo "merge  $ii -> $newFile"
	hadd $newFile `cat $ii`  
    done

    nFiles=`ls merged_output_${mergeLevel}_*.root | wc -l`
    if [ $nFiles -gt 1 ] ; then
	ls merged_output_${mergeLevel}_*.root > inlist.${mergeLevel}.list
	mergeFolder
    else
	mv merged_output_${mergeLevel}_000.root ${output_name}.root
	mkdir -p attic
	mv merged_output*.root attic
	mv splitLists* attic
	mv *.list attic
    fi
}

# ----------------------------------------------------
#version for more than 10k files in a folder
mergeFolder10k () {

    inSplitList=inlist.${mergeLevel}.list

    let mergeLevel=mergeLevel+1
    outSplitListFolder=splitLists.$mergeLevel

    # -- split
    mkdir -p $outSplitListFolder
    split -l 10 -d -a 4 $inSplitList ${outSplitListFolder}/inlist.sublist.${mergeLevel}. #5 files
    
    for ii in `find  ${outSplitListFolder} -name "inlist.sublist.${mergeLevel}.*"` ; do 
	sub=`echo $ii | cut -d '.' -f 5`
	newFile=merged_output_${mergeLevel}_${sub}.root
	echo "merge  $ii -> $newFile"
	hadd $newFile `cat $ii`  
    done

    nFiles=`ls merged_output_${mergeLevel}_*.root | wc -l`
    if [ $nFiles -gt 1 ] ; then
	ls merged_output_${mergeLevel}_*.root > inlist.${mergeLevel}.list
	mergeFolder10k
    else
	mv merged_output_${mergeLevel}_0000.root ${output_name}.root
	mkdir -p attic
	mv merged_output*.root attic
	mv splitLists* attic
	mv *.list attic
    fi
}

# ----------------------------------------------------

if [ -d ${baseDir}/merge ] ; then
    echo "files in ${baseDir} already merged"
		exit
fi

mkdir -p ${baseDir}/merge

inList=${baseDir}/merge/inlist.0.list
    
if [ -f ${inList} ] ; then
    rm ${inList}
fi

touch ${inList}

current=`pwd`
for ii in `find ${current}/${baseDir} -name "*.root" | sort` ; do
    size=`stat -c %s $ii`
    if [ $size -gt 5000 ] ; then
	echo "$ii" >> ${inList}
	((nFilesToMerge++))
    fi
done

echo "merging ${nFilesToMerge} files"


pushd $baseDir/merge > /dev/null
mergeLevel=0
if [ $nFilesToMerge -gt 10000 ] ; then
	mergeFolder10k
else 
	mergeFolder
fi

popd > /dev/null
