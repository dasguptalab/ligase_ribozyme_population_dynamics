#!/bin/bash

# script to cluster sequences using clustalo
# usage: bash 05_formatReads.sh 

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_local.txt" | tr -d " " | sed "s/outputs://g")

# set trimming directory
filterOut=$outputsPath"/filteredByRegion"

# make a new directory for analysis
formatOut=$outputsPath"/formatted"
mkdir $formatOut
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $formatOut directory already exsists... please remove before proceeding."
	exit 1
fi

# move to the new directory
cd $formatOut

# status message
echo "Beginning analysis..."

# combine all forward and reverse reads
for f1 in $filterOut"/"*; do
	# get sequence run name
	newName=$(basename $f1 | sed "s/_L001_p.*\.flt\.fq//g")
	# re-format sequences
	cat $f1 | cut -d" " -f1 | awk 'NR%4==1 || NR%4==2' | sed "s/^@.*/>$newName/g" >> $formatOut"/combined.fmt.fasta"
done

# status message
echo "Analysis complete!"
