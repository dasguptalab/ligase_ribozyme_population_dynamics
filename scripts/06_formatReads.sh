#!/bin/bash

# script to cluster sequences using clustalo
# usage: bash 06_formatReads.sh

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputs/inputPaths_HPC.txt" | tr -d " " | sed "s/outputs://g")

# set trimming directory
inputsPath=$outputsPath"/filtered_combined"

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

# re-format read headers for each run
for f1 in $inputsPath"/"*; do
	# create sequence run name
	#newName=$(basename $f1 | sed "s/_L001_p.*\.flt\.fq//g")
	# get file name
	fileName=$(basename $f1 | sed "s/\.fa//g")
	# status message
	echo "Processing $fileName ..."
	# re-format read headers
	#cat $f1 | cut -d" " -f1 | awk 'NR%4==1 || NR%4==2' | sed "s/^@.*/>$newName/g" >> $formatOut"/combined.fmt.fasta"
	cat $f1 | cut -d" " -f1 | sed "s/^@/>/g" >> $formatOut"/"$fileName".fmt.fa"
done

# status message
echo "Analysis complete!"
