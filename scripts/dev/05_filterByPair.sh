#!/bin/bash

# script to filter fastq files and keep only a single read per sequence, preferably the forward read
# usage: bash 05_filterByPair.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set trimming directory
filterRegionOut=$outputsPath"/filteredByRegionLength"

# make a new directory for analysis
filterReadsOut=$outputsPath"/filteredByPair"
mkdir $filterReadsOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $filterReadsOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $filterReadsOut

# name tmp data directory
tmpDir=$filterRegionOut"/tmp"

# make tmp data directory
mkdir $tmpDir

# status message
echo "Beginning analysis..."

# filter to keep only a single read per sequence, preferably the forward read
for f1 in $filterRegionOut"/"*_pForward\.flt\.fq; do
	# status message
	echo "Processing $f1"
	# get sequence run name
	newName=$(basename $f1 | sed "s/_L001_p.*\.flt\.fq//g")
	# set inputs files
	forwardFile=$f1
	reverseFile=$filterRegionOut"/"$newName"_L001_pReverse.flt.fq"
	reverseNames=$tmpDir"/"$newName"_L001_pReverse.names.flt.fq"
	# retrieve read names from headers
	cat $reverseFile | awk 'NR % 4 == 1' | cut -d" " -f1 > $reverseNames
	# name output file
	outFile=$filterReadsOut"/"$newName".flt.fq"
	# copy the forward read data
	cat $forwardFile  > $outFile
	# loop over each reverse read
	while read line; do
		# status message
		echo "Processing $line"
		# check if current read is already in the forward read data
		if ! grep -q $line $outFile; then # doesn't exsist
			# add the current reverse read
			cat $reverseFile | grep -A1 $line >> $outFile
		fi
	done < $reverseNames
	# status message
	echo "$forwardFile processed!"
done

# clean up
rm -r $tmpDir

# status message
echo "Analysis complete!"
